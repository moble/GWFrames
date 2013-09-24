/* $Id: DocumentationServlet.java 570 2010-05-20 14:23:11Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.webapp;

import uk.ac.ed.ph.snuggletex.InputError;
import uk.ac.ed.ph.snuggletex.SerializationMethod;
import uk.ac.ed.ph.snuggletex.SnuggleEngine;
import uk.ac.ed.ph.snuggletex.SnuggleInput;
import uk.ac.ed.ph.snuggletex.SnuggleLogicException;
import uk.ac.ed.ph.snuggletex.SnuggleSession;
import uk.ac.ed.ph.snuggletex.WebPageOutputOptions;
import uk.ac.ed.ph.snuggletex.WebPageOutputOptionsTemplates;
import uk.ac.ed.ph.snuggletex.DOMOutputOptions.ErrorOutputOptions;
import uk.ac.ed.ph.snuggletex.WebPageOutputOptions.WebPageType;
import uk.ac.ed.ph.snuggletex.definitions.W3CConstants;
import uk.ac.ed.ph.snuggletex.internal.util.IOUtilities;
import uk.ac.ed.ph.snuggletex.jeuclid.JEuclidUtilities;
import uk.ac.ed.ph.snuggletex.jeuclid.SimpleMathMLImageSavingCallback;
import uk.ac.ed.ph.snuggletex.utilities.MessageFormatter;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.xml.transform.Transformer;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Rather trivial (and badly factored) servlet that serves up the documentation pages for
 * SnuggleTeX. Main points:
 * <ul>
 *   <li>
 *     Documentation is written in SnuggleTeX format as <tt>.tex</tt> files under
 *     {@link #TEX_SOURCE_BASE_RESOURCE}
 *   </li>
 *   <li>
 *     Servlet responds to requests under a particular base URL as defined in <tt>web.xml</tt>.
 *     The extra path info is mapped to a "resource path".
 *   </li>
 *   <li>
 *     File extension determines which web output to serve.
 *   </li>
 *   <li>
 *     Servlet stores generated versions of documentation in a temp directory. If not found,
 *     it tries to recreate from scratch using the <tt>.tex</tt> file.
 *   </li>
 *   <li>
 *     Limited caching is available. If off, we recreate each documentation resource each time,
 *     otherwise we create once and keep forever.
 *   </li>
 * </ul>
 *
 * @author  David McKain
 * @version $Revision: 570 $
 */
public final class DocumentationServlet extends BaseServlet {
    
    private static final long serialVersionUID = -4091098938512353014L;

    static final Logger logger = LoggerFactory.getLogger(DocumentationServlet.class);
    
    /** <tt>init-param</tt> controlling whether we are caching or not */
    private static final String CACHING_PARAM = "caching";
    
    /** Location of XSLT for formatting the resulting web pages */
    private static final String FORMAT_OUTPUT_XSLT_URI = "classpath:/format-output.xsl";
    
    /** Location of macros TeX file (relative to webapp) */
    private static final String MACROS_RESOURCE_LOCATION = "/WEB-INF/macros.tex";
    
    /** Base Location for <tt>.tex</tt> source files (relative to webapp) */
    private static final String TEX_SOURCE_BASE_RESOURCE = "/WEB-INF/docs";
    
    /** Maps supported extensions to Content Types, for all extensions we support here */
    private final Map<String, String> extensionToContentTypeMap;

    /** Maps extension to {@link WebPageType}, for those generated that way */
    private final Map<String, WebPageType> extensionToWebPageTypeMap;

    /** Whether to cache result or not */
    private boolean caching;
    
    /** Directory in which Files created and cached by this servlet will be stored. */
    private File baseDirectory;
    
    @Override
    public void init() throws ServletException {
        /* Set up base directory */
        if (baseDirectory==null) {
            try {
                baseDirectory = File.createTempFile("snuggletex-", ".dir");
            }
            catch (IOException e) {
                throw new ServletException("Could not create initial tempfile for storing documentation", e);
            }
            if (!baseDirectory.delete()) {
                throw new ServletException("Could not delete tempfile at " + baseDirectory
                        + " for re-creating as a directory");
            }
            if (!baseDirectory.mkdir()) {
                throw new ServletException("Could not create base directory at " + baseDirectory
                        + " for storing documentation");
            }
            logger.info("Set base directory for documentation as {}", baseDirectory);
        }
        
        /* Check whether caching is turned on or not */
        caching = "true".equals(getServletConfig().getInitParameter(CACHING_PARAM));
    }
    
    public DocumentationServlet() {
        /* Define all supported content types */
        this.extensionToContentTypeMap = new HashMap<String, String>();
        extensionToContentTypeMap.put("xhtml", "application/xhtml+xml");
        extensionToContentTypeMap.put("xml", "application/xhtml+xml");
        extensionToContentTypeMap.put("cxml", "application/xhtml+xml");
        extensionToContentTypeMap.put("htm", "text/html");
        extensionToContentTypeMap.put("html", "text/html");
        extensionToContentTypeMap.put("tex", "application/x-latex");
        extensionToContentTypeMap.put("png", "image/png"); /* (Required for JEuclid MathML images) */
        
        /* Determines which extensions to use for standard web page outputs */
        this.extensionToWebPageTypeMap = new HashMap<String, WebPageType>();
        extensionToWebPageTypeMap.put("xhtml", WebPageType.MOZILLA);
        extensionToWebPageTypeMap.put("xml", WebPageType.UNIVERSAL_STYLESHEET);
        extensionToWebPageTypeMap.put("cxml", WebPageType.CROSS_BROWSER_XHTML);
        extensionToWebPageTypeMap.put("htm", WebPageType.MATHPLAYER_HTML);
        extensionToWebPageTypeMap.put("html", WebPageType.PROCESSED_HTML);
    }
    
    @Override
    protected void doGet(HttpServletRequest request, HttpServletResponse response)
            throws ServletException, IOException {
        /* Determine which resource to serve up, and how it's going to be served */
        /* Work out what is to be served using the file extension to determine what to do */
        String resourcePath = request.getPathInfo();
        int lastDotPosition = resourcePath.lastIndexOf(".");
        if (lastDotPosition==-1) {
            response.sendError(HttpServletResponse.SC_NOT_FOUND, "Could not locate '.' in resourcePath " + resourcePath);
            return;
        }
        String resourceBaseName = resourcePath.substring(0, lastDotPosition);
        String extension = resourcePath.substring(lastDotPosition+1);
        File resourceFile = mapResourcePath(resourcePath);
        
        /* Make sure content type is known */
        String contentType = extensionToContentTypeMap.get(extension);
        if (contentType==null) {
            response.sendError(HttpServletResponse.SC_NOT_FOUND, "Unsupported documentation resource extension " + extension);
            return;
        }
        
        /* Decide whether to (re)generate resource or not. This will happen if:
         * 
         * 1. It's a "proper" documentation page (i.e. not PNG).
         * 2. Cached resource doesn't already exist (or caching is turned off)
         */
        if (!"png".equals(extension) && (!resourceFile.exists() || !caching)) {
            resourceFile = generateResource(request, resourcePath, resourceBaseName, extension);
            if (resourceFile==null) {
                response.sendError(HttpServletResponse.SC_NOT_FOUND, "Documentation page does not exist");
                return;
            }
        }
        
        /* Serve up */
        response.setContentType(contentType);
        response.setContentLength((int) resourceFile.length());
        IOUtilities.transfer(new FileInputStream(resourceFile), response.getOutputStream(), true, false);
    }
    
    /**
     * Maps a requested "resource path" to a file under our {@link #baseDirectory} "cache".
     */
    private File mapResourcePath(final String resourcePath) {
        return new File(baseDirectory + File.separator + resourcePath.replace("/", File.separator));
    }
    
    /**
     * Generates the documentation resource at the given path.
     * 
     * @return resulting File, or null if the source TeX file for this resource couldn't be
     *   located or if the file extension couldn't be understood.
     */
    private File generateResource(final HttpServletRequest request, final String resourcePath,
            final String resourceBaseName, final String extension)
            throws ServletException, IOException {
        logger.info("Generating Resource at {}", resourcePath);
        
        /* Locate TeX source */
        String texSourceName = resourceBaseName + ".tex";
        String texSourceResourcePath = TEX_SOURCE_BASE_RESOURCE + texSourceName;
        InputStream texSourceStream = getServletContext().getResourceAsStream(texSourceResourcePath);
        if (texSourceStream==null) {
            logger.info("Could not locate source resource at {}", texSourceResourcePath);
            return null;
        }
        
        /* Decide what to do based on extension */
        File resultFile = null;
        if (extension.equals("tex")) {
            /* Just copy TeX resource over */
            resultFile = IOUtilities.ensureFileCreated(mapResourcePath(resourcePath));
            IOUtilities.transfer(texSourceStream, new FileOutputStream(resultFile));
        }
        else if (extension.equals("png")) {
            throw new SnuggleLogicException("PNG generation is done as a side-effect so shouldn't have been requested directly");
        }
        else {
            /* Use SnuggleTeX standard web page builder (possibly with JEuclid stuff) */
            WebPageType webPageType = extensionToWebPageTypeMap.get(extension);
            if (webPageType==null) {
                logger.info("Resource extension {} not understood", extension);
                return null;
            } 
            String imageOutputDirectortyResourcePath = resourceBaseName;
            String imageOutputBaseUrl = request.getContextPath() + request.getServletPath() + resourceBaseName;
            resultFile = generateSnuggledFile(request, texSourceStream, texSourceResourcePath,
                    webPageType, resourcePath, imageOutputDirectortyResourcePath, imageOutputBaseUrl);
        }
        return resultFile;
    }
    
    private File generateSnuggledFile(final HttpServletRequest request,
            final InputStream texSourceStream, final String texSourceResourcePath,
            final WebPageType webPageType, final String outputResourcePath,
            final String imageOutputDirectoryResourcePath, final String imageOutputBaseURL)
            throws ServletException, IOException {
        /* Parse macros.tex and source resource */
        InputStream macrosResource = ensureReadResource(MACROS_RESOURCE_LOCATION);
        SnuggleEngine engine = createSnuggleEngine();
        SnuggleSession session = engine.createSession();
        session.parseInput(new SnuggleInput(macrosResource, "Web resource at " + MACROS_RESOURCE_LOCATION));
        session.parseInput(new SnuggleInput(texSourceStream, "Web resource at " + texSourceResourcePath));
        
        /* Work out basic SnuggleTeX options */
        WebPageOutputOptions options = WebPageOutputOptionsTemplates.createWebPageOptions(webPageType);
        options.setErrorOutputOptions(ErrorOutputOptions.XHTML);
        options.setMathVariantMapping(true);
        options.setAddingMathSourceAnnotations(true);
        options.setIndenting(true);
        options.setIncludingStyleElement(false);
        if (webPageType==WebPageType.PROCESSED_HTML) {
            /* Create folder for storing MathML images. */
            File imageOutputDirectory = IOUtilities.ensureDirectoryCreated(mapResourcePath(imageOutputDirectoryResourcePath));
            ImageSavingCallback callback = new ImageSavingCallback(imageOutputDirectory, imageOutputBaseURL);
            
            /* We'll actually generate XHTML 1.0 Strict here so that I can be
             * anally retentive about having valid documentation! */
            options.setDoctypePublic(W3CConstants.XHTML_10_STRICT_PUBLIC_IDENTIFIER);
            options.setDoctypeSystem(W3CConstants.XHTML_10_STRICT_SYSTEM_IDENTIFIER);
            options.setSerializationMethod(SerializationMethod.XHTML);
            
            /* Configure JEuclid post-processor, with down-conversion */
            JEuclidUtilities.setupJEuclidPostProcessors(options, true, callback);
        }
        else if (webPageType==WebPageType.UNIVERSAL_STYLESHEET) {
            /* Point to our own version of the USS if required */
            options.setClientSideXSLTStylesheetURLs(request.getContextPath() + "/includes/pmathml.xsl");
        }
        
        /* Set up stylesheet to format the output */
        Transformer stylesheet = getStylesheet(request, FORMAT_OUTPUT_XSLT_URI);
        stylesheet.setParameter("page-type", webPageType!=WebPageType.PROCESSED_HTML ? webPageType.name() : null);
        options.setStylesheets(stylesheet);
        
        /* Generate output file */
        File outputFile = IOUtilities.ensureFileCreated(mapResourcePath(outputResourcePath));
        boolean success = session.writeWebPage(options, new FileOutputStream(outputFile));
        
        /* Log any errors or failures */
        List<InputError> errors = session.getErrors();
        if (!errors.isEmpty()) {
            logger.warn("Errors occurred generating resource {}" ,outputResourcePath);
            for (InputError error : errors) {
                logger.warn("Error: " + MessageFormatter.formatErrorAsString(error));
            }
        }
        if (!success) {
            logger.warn("Failed to generate resulting resource at " + outputFile);
        }
        return outputFile;
    }
    
    /**
     * Implementation of {@link SimpleMathMLImageSavingCallback} that stores images in the
     * given output directory with the given base URL, using a very simple naming scheme.
     */
    private class ImageSavingCallback extends SimpleMathMLImageSavingCallback {
        
        private final File imageOutputDirectory;
        private final String imageOutputBaseURL;
        
        public ImageSavingCallback(final File imageOutputDirectory, final String imageOutputBaseURL) {
            this.imageOutputDirectory = imageOutputDirectory;
            this.imageOutputBaseURL = imageOutputBaseURL;
        }
        
        @Override
        public File getImageOutputFile(int mathmlCounter) {
            return new File(imageOutputDirectory, getImageName(mathmlCounter));
        }
        
        @Override
        public OutputStream getImageOutputStream(int mathmlCounter) {
            /* Not used here */
            return null;
        }
        
        @Override
        public String getImageURL(int mathmlCounter) {
            return imageOutputBaseURL + "/" + getImageName(mathmlCounter);
        }
 
        private String getImageName(int mathmlCounter) {
            /* We just append the page base name with the image number, which will be unique
             * so is fine for us.
             */
            return "mathml-" + mathmlCounter + ".png";
        }
        
        public void imageSavingFailed(Object imageOutputObject, int mathmlCounter, String contentType,
                Throwable exception) {
            logger.warn("Image saving failed", exception);
        }
    }
}
