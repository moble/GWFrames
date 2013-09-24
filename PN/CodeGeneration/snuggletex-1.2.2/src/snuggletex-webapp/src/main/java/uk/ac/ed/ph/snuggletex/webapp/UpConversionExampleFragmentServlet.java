/* $Id:FullLaTeXInputDemoServlet.java 158 2008-07-31 10:48:14Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.webapp;

import static uk.ac.ed.ph.snuggletex.utilities.MathMLUtilities.extractAnnotationString;
import static uk.ac.ed.ph.snuggletex.utilities.MathMLUtilities.isolateAnnotationXML;
import static uk.ac.ed.ph.snuggletex.utilities.MathMLUtilities.isolateFirstSemanticsBranch;
import static uk.ac.ed.ph.snuggletex.utilities.MathMLUtilities.serializeDocument;
import static uk.ac.ed.ph.snuggletex.utilities.MathMLUtilities.serializeElement;

import uk.ac.ed.ph.snuggletex.DOMOutputOptions;
import uk.ac.ed.ph.snuggletex.InputError;
import uk.ac.ed.ph.snuggletex.SerializationSpecifier;
import uk.ac.ed.ph.snuggletex.SnuggleEngine;
import uk.ac.ed.ph.snuggletex.SnuggleInput;
import uk.ac.ed.ph.snuggletex.SnuggleSession;
import uk.ac.ed.ph.snuggletex.WebPageOutputOptions;
import uk.ac.ed.ph.snuggletex.WebPageOutputOptionsTemplates;
import uk.ac.ed.ph.snuggletex.DOMOutputOptions.ErrorOutputOptions;
import uk.ac.ed.ph.snuggletex.SnuggleSession.EndOutputAction;
import uk.ac.ed.ph.snuggletex.WebPageOutputOptions.WebPageType;
import uk.ac.ed.ph.snuggletex.definitions.W3CConstants;
import uk.ac.ed.ph.snuggletex.internal.util.XMLUtilities;
import uk.ac.ed.ph.snuggletex.upconversion.MathMLUpConverter;
import uk.ac.ed.ph.snuggletex.upconversion.UpConvertingPostProcessor;
import uk.ac.ed.ph.snuggletex.upconversion.internal.UpConversionPackageDefinitions;
import uk.ac.ed.ph.snuggletex.utilities.MessageFormatter;

import java.io.IOException;
import java.net.URLDecoder;
import java.util.ArrayList;
import java.util.List;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.xml.transform.Transformer;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

/**
 * Variant of {@link UpConversionDemoServlet} that demonstrates a page fragment
 * containing similar results. This is used for showing examples inside DHTML dialog
 * boxes.
 * 
 * @author  David McKain
 * @version $Revision:158 $
 */
public final class UpConversionExampleFragmentServlet extends BaseServlet {
    
    private static final long serialVersionUID = 4376587500238353176L;
    
    /** Logger so that we can log what users are trying out to allow us to improve things */
    private static Logger logger = LoggerFactory.getLogger(UpConversionExampleFragmentServlet.class);
    
    /** Location of XSLT controlling page layout */
    private static final String DISPLAY_XSLT_LOCATION = "classpath:/upconversion-example-fragment.xsl";
    
    /** Default assumptions to use here */
    private static final String DEFAULT_UPCONVERSION_OPTIONS =
        "\\assumeSymbol{e}{exponentialNumber}\n"
        + "\\assumeSymbol{f}{function}\n"
        + "\\assumeSymbol{f_n}{function}\n"
        + "\\assumeSymbol{g}{function}\n"
        + "\\assumeSymbol{i}{imaginaryNumber}\n"
        + "\\assumeSymbol{\\pi}{constantPi}\n"
        + "\\assumeSymbol{\\gamma}{eulerGamma}";
    
    @Override
    protected void doGet(HttpServletRequest request, HttpServletResponse response)
            throws ServletException, IOException {
        doRequest(request, response);
    }
    
    private void doRequest(HttpServletRequest request, HttpServletResponse response)
            throws ServletException, IOException {
        /* Read in input LaTeX, which must be provided */
        String rawInputLaTeX = request.getParameter("input");
        if (rawInputLaTeX==null) {
            response.sendError(HttpServletResponse.SC_BAD_REQUEST, "No input provided");
            return;
        }
        String inputLaTeX = rawInputLaTeX.replaceAll("\\s+", " ");
        
        /* I've noticed that some search engines appear to be double-encoding the query string,
         * so I'll detect this now and decode if required.
         */
        if (inputLaTeX.contains("%")) {
            inputLaTeX = URLDecoder.decode(inputLaTeX, "UTF-8");
        }
        
        /* Parse the LaTeX */
        SnuggleEngine engine = createSnuggleEngine();
        engine.addPackage(UpConversionPackageDefinitions.getPackage());
        SnuggleSession session = engine.createSession();
        if (inputLaTeX.endsWith("$") || inputLaTeX.endsWith("\\]") || inputLaTeX.endsWith("\\)")) {
            /* Author has explicitly ended Math mode, so is probably doing some custom assumptions */
            session.parseInput(new SnuggleInput(inputLaTeX, "Query Input"));
        }
        else {
            /* Parse whole thing in Math mode using default assumptions */
            session.parseInput(new SnuggleInput(DEFAULT_UPCONVERSION_OPTIONS, "Default Assumptions Input"));
            session.parseInput(new SnuggleInput("\\[ " + inputLaTeX + " \\]", "Query Input"));
        }
        
        /* Create raw DOM, without any up-conversion for the time being. I've done this
         * so that we can show how much the PMathML hopefully improves after up-conversion!
         */
        Document resultDocument = XMLUtilities.createNSAwareDocumentBuilder().newDocument();
        Element resultRoot = resultDocument.createElement("root");
        resultDocument.appendChild(resultRoot);
        DOMOutputOptions domOptions = new DOMOutputOptions();
        domOptions.setMathVariantMapping(true);
        domOptions.setAddingMathSourceAnnotations(true);
        domOptions.setErrorOutputOptions(ErrorOutputOptions.NO_OUTPUT);
        session.buildDOMSubtree(resultRoot, domOptions);
        
        /* See if parsing succeeded and generated a single <math/> element. We'll only continue
         * up-converting if this happened.
         */
        NodeList resultNodeList = resultRoot.getChildNodes();
        List<InputError> errors = session.getErrors();
        
        Element mathMLElement = null;
        String parallelMathML = null;
        String pMathMLInitial = null;
        String pMathMLUpConverted = null;
        String cMathML = null;
        String maximaInput = null;
        List<Element> parsingErrors = null;
        UpConvertingPostProcessor upConvertingPostProcessor = new UpConvertingPostProcessor();
        SerializationSpecifier sourceSerializationOptions = createMathMLSourceSerializationOptions();
        boolean badInput = false;
        
        /* Check for any errors and that the shape of the result is as expected */
        if (!errors.isEmpty()) {
            /* Input error occurred */
            parsingErrors = new ArrayList<Element>();
            for (InputError error : errors) {
                parsingErrors.add(MessageFormatter.formatErrorAsXML(resultDocument, error, true));
            }
        }
        else {
            /* Make sure there is exactly one MathML element, and any number of upconversion options
             * specifiers.
             */
            mathMLElement = extractMathMLElement(resultNodeList, true);
            if (mathMLElement==null) {
                badInput = true;
            }
        }
        
        /* Happy path! */
        if (mathMLElement!=null) {
            /* We found exactly one MathML element out of raw results */
            pMathMLInitial = serializeElement(mathMLElement, sourceSerializationOptions);
            
            /* Do up-conversion and extract wreckage */
            MathMLUpConverter upConverter = new MathMLUpConverter(getStylesheetCache());
            Document upConvertedMathDocument = upConverter.upConvertSnuggleTeXMathML(mathMLElement.getOwnerDocument(), upConvertingPostProcessor.getUpconversionOptions());
            
            mathMLElement = extractMathMLElement(upConvertedMathDocument.getDocumentElement().getChildNodes(), false);
            parallelMathML = serializeElement(mathMLElement, sourceSerializationOptions);
            
            pMathMLUpConverted = serializeDocument(isolateFirstSemanticsBranch(mathMLElement), sourceSerializationOptions);
            Document cMathMLDocument = isolateAnnotationXML(mathMLElement, MathMLUpConverter.CONTENT_MATHML_ANNOTATION_NAME);
            cMathML = cMathMLDocument!=null ? serializeDocument(cMathMLDocument, sourceSerializationOptions) : null;
            maximaInput = extractAnnotationString(mathMLElement, MathMLUpConverter.MAXIMA_ANNOTATION_NAME);
        }
        
        /* Only log failures, as this would normally be caused by bad documentation but could
         * also be due to a clever clogs user! */
        if (!errors.isEmpty()) {
            logger.error("Input: {}", inputLaTeX);
            logger.error("Final MathML: {}", parallelMathML);
            logger.error("Error count: {}", errors.size());
            for (InputError error : errors) {
                logger.error("Error: " + MessageFormatter.formatErrorAsString(error));
            }
        }
        
        /* We'll cheat slightly and bootstrap off the SnuggleTeX web page generation process,
         * even though most of the interesting page content is going to be fed in as stylesheet
         * parameters.
         * 
         * These fragments are always going to generate plain old XHTML, so we'll use the
         * PROCESSED_HTML output for this.
         */
        WebPageOutputOptions webOptions = WebPageOutputOptionsTemplates.createWebPageOptions(WebPageType.PROCESSED_HTML);
        webOptions.setDoctypePublic(W3CConstants.XHTML_10_STRICT_PUBLIC_IDENTIFIER);
        webOptions.setDoctypeSystem(W3CConstants.XHTML_10_STRICT_SYSTEM_IDENTIFIER);
        webOptions.setMathVariantMapping(true);
        webOptions.setAddingMathSourceAnnotations(true);
        webOptions.setIndenting(true);
        webOptions.setIncludingStyleElement(false);
        
        /* Create XSLT to generate the resulting page */
        Transformer viewStylesheet = getStylesheet(request, DISPLAY_XSLT_LOCATION);
        viewStylesheet.setParameter("latex-input", inputLaTeX);
        viewStylesheet.setParameter("is-bad-input", Boolean.valueOf(badInput));
        viewStylesheet.setParameter("parsing-errors", parsingErrors);
        viewStylesheet.setParameter("parallel-mathml", parallelMathML);
        viewStylesheet.setParameter("pmathml-initial", pMathMLInitial);
        viewStylesheet.setParameter("pmathml-upconverted", pMathMLUpConverted);
        viewStylesheet.setParameter("cmathml", cMathML);
        viewStylesheet.setParameter("maxima-input", maximaInput);
        webOptions.setStylesheets(viewStylesheet);
        
        /* Generate and serve the resulting web page */
        try {
            session.writeWebPage(webOptions, response, response.getOutputStream(), EndOutputAction.FLUSH);
        }
        catch (Exception e) {
            throw new ServletException("Unexpected Exception", e);
        }
    }
}