/* $Id: BaseServlet.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.webapp;

import static uk.ac.ed.ph.snuggletex.utilities.MathMLUtilities.isMathMLElement;

import uk.ac.ed.ph.snuggletex.DownConvertingPostProcessor;
import uk.ac.ed.ph.snuggletex.SerializationSpecifier;
import uk.ac.ed.ph.snuggletex.SnuggleConstants;
import uk.ac.ed.ph.snuggletex.SnuggleEngine;
import uk.ac.ed.ph.snuggletex.WebPageOutputOptions;
import uk.ac.ed.ph.snuggletex.WebPageOutputOptionsTemplates;
import uk.ac.ed.ph.snuggletex.WebPageOutputOptions.WebPageType;
import uk.ac.ed.ph.snuggletex.definitions.W3CConstants;
import uk.ac.ed.ph.snuggletex.upconversion.UpConversionUtilities;
import uk.ac.ed.ph.snuggletex.utilities.ClassPathURIResolver;
import uk.ac.ed.ph.snuggletex.utilities.SerializationOptions;
import uk.ac.ed.ph.snuggletex.utilities.StylesheetCache;
import uk.ac.ed.ph.snuggletex.utilities.StylesheetManager;

import java.io.InputStream;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.xml.transform.Templates;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerConfigurationException;

import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 * Trivial base class for servlets in the demo webapp
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
abstract class BaseServlet extends HttpServlet {
    
    private static final long serialVersionUID = -2577813908466694931L;
    
    public String ensureGetContextInitParam(String propertyName) throws ServletException {
        String result = getServletContext().getInitParameter(propertyName);
        if (result==null) {
            throw new ServletException("Context init-param " + propertyName + " is not set");
        }
        return result;
    }

    /**
     * Helper that reads in a resource from the webapp hierarchy, throwing a {@link ServletException}
     * if the resource could not be found.
     * 
     * @param resourcePathInsideWebpp path of Resource to load, relative to base of webapp.
     * @return resulting {@link InputStream}, which will not be null
     * @throws ServletException
     */
    protected InputStream ensureReadResource(String resourcePathInsideWebpp) throws ServletException {
        InputStream result = getServletContext().getResourceAsStream(resourcePathInsideWebpp);
        if (result==null) {
            throw new ServletException("Could not read in required web resource at " + resourcePathInsideWebpp);
        }
        return result;
    }
    
    protected StylesheetManager getStylesheetManager() {
        return (StylesheetManager) getServletContext().getAttribute(ContextInitialiser.STYLESHEET_MANAGER_ATTRIBUTE_NAME);
    }
    
    protected StylesheetCache getStylesheetCache() {
        return getStylesheetManager().getStylesheetCache();
    }
    
    protected SnuggleEngine createSnuggleEngine() {
        return new SnuggleEngine(getStylesheetManager());
    }
    
    /**
     * Compiles the XSLT stylesheet at the given location within the webapp,
     * using {@link ClassPathURIResolver} to locate the stylesheet and anything it wants to import.
     * <p>
     * It also sets some core parameters based on certain properties set for the webapp.
     * 
     * @param request Request being processed (so we can pass the context path to the XSLT)
     * @param classPathUri location of XSLT to compile.
     * 
     * @return resulting {@link Templates} representing the compiled stylesheet.
     * @throws ServletException if XSLT could not be found or could not be compiled.
     */
    protected Transformer getStylesheet(HttpServletRequest request, String classPathUri) throws ServletException {
        Transformer result;
        try {
            result = getStylesheetManager().getStylesheet(classPathUri).newTransformer();
        }
        catch (TransformerConfigurationException e) {
            throw new ServletException("Could not create Transformer from Templates", e);
        }
        result.setParameter("context-path", request.getContextPath());
        result.setParameter("snuggletex-version", ensureGetContextInitParam(ContextInitialiser.SNUGGLETEX_VERSION_PROPERTY_NAME));
        result.setParameter("maven-site-url", ensureGetContextInitParam(ContextInitialiser.MAVEN_SITE_URL_PROPERTY_NAME));
        return result;
    }
    
    protected SerializationSpecifier createMathMLSourceSerializationOptions() {
        SerializationSpecifier result = new SerializationOptions();
        result.setIndenting(true);
        result.setUsingNamedEntities(true);
        return result;
    }

    /**
     * Convenience method which picks the most appropriate MathML-based {@link WebPageType}
     * for the current UserAgent, returning {@link WebPageType#PROCESSED_HTML}
     * if the UserAgent does not appear to support MathML as this is the only sensible
     * option in that case.
     * 
     * @param request
     */
    protected WebPageType chooseBestWebPageType(final HttpServletRequest request) {
        String userAgent = request.getHeader("User-Agent");
        WebPageType result = WebPageType.PROCESSED_HTML;
        if (userAgent!=null) {
            if (userAgent.contains("MathPlayer ")) {
                result  = WebPageType.MATHPLAYER_HTML;
            }
            else if (userAgent.contains("Gecko/")) {
                result  = WebPageType.MOZILLA;
            }
        }
        return result;
    }
    
    /**
     * Convenience method to choose the most appropriate base {@link WebPageOutputOptions}
     * for the current UserAgent, using {@link #chooseBestWebPageType(HttpServletRequest)}
     * to determine the underlying {@link WebPageType}.
     * 
     * @param request
     */
    protected WebPageOutputOptions chooseBestBaseWebPageOutputOptions(final HttpServletRequest request) {
        /* Set common options */
        WebPageType webPageType = chooseBestWebPageType(request);
        WebPageOutputOptions result = WebPageOutputOptionsTemplates.createWebPageOptions(webPageType);
        result.setMathVariantMapping(true);
        result.setAddingMathSourceAnnotations(true);
        result.setIndenting(true);
        result.setIncludingStyleElement(false);

        /* Then additional suitable options */
        if (webPageType==WebPageType.PROCESSED_HTML) {
            result.setDoctypePublic(W3CConstants.XHTML_10_STRICT_PUBLIC_IDENTIFIER);
            result.setDoctypeSystem(W3CConstants.XHTML_10_STRICT_SYSTEM_IDENTIFIER);
            
            /* If browser can't handle MathML, we'll add post-processors to down-convert
             * simple expressions to XHTML + CSS and replace the remaining MathML islands
             * with dynamically generated images. */
            result.setDOMPostProcessors(
                    new DownConvertingPostProcessor(),
                    new MathMLToImageLinkPostProcessor(request.getContextPath())
            );
        }
        else if (webPageType==WebPageType.MATHPLAYER_HTML) {
            result.setDoctypePublic(W3CConstants.XHTML_10_STRICT_PUBLIC_IDENTIFIER);
            result.setDoctypeSystem(W3CConstants.XHTML_10_STRICT_SYSTEM_IDENTIFIER);
        }
        else {
            result.setDoctypePublic(W3CConstants.XHTML_11_MATHML_20_PUBLIC_IDENTIFIER);
            result.setDoctypeSystem(W3CConstants.XHTML_11_MATHML_20_SYSTEM_IDENTIFIER);
        }
        return result;
    }
    
    protected boolean isMathMLCapable(final HttpServletRequest request) {
        String userAgent = request.getHeader("User-Agent");
        return userAgent!=null && (userAgent.contains("MathPlayer ") || userAgent.contains("Gecko/"));
    }
    
    protected boolean isInternetExplorer(final HttpServletRequest request) {
        String userAgent = request.getHeader("User-Agent");
        return userAgent!=null && userAgent.contains("MSIE");
    }

    /**
     * Untangles the given {@link NodeList} to find an expected single MathML element,
     * and possibly some whitespace and possibly any number of <c:upconversion-options/>
     * elements.
     */
    protected Element extractMathMLElement(final NodeList resultNodeList, final boolean allowUpConversionOptionsElements) {
        /* Make sure there is exactly one MathML element, and any number of upconversion options
         * specifiers.
         */
        Element result = null;
        for (int i=0, size=resultNodeList.getLength(); i<size; i++) {
            Node node = resultNodeList.item(i);
            if (isMathMLElement(node)) {
                if (result!=null) {
                    return null;
                }
                result = (Element) node;
            }
            else if (node.getNodeType()==Node.TEXT_NODE && node.getNodeValue().trim().length()==0) {
                /* Ignore whitespace */
                continue;
            }
            else if (allowUpConversionOptionsElements
                    && node.getNodeType()==Node.ELEMENT_NODE
                    && SnuggleConstants.SNUGGLETEX_NAMESPACE.equals(node.getNamespaceURI())
                    && UpConversionUtilities.UPCONVERSION_OPTIONS_XML_LOCAL_NAME.equals(node.getLocalName())) {
                /* Allow <s:upconversion-options/> */
                continue;
            }
            else {
                /* Anything else is an not allowed */
                return null;
            }
        }
        return result;
    }
}