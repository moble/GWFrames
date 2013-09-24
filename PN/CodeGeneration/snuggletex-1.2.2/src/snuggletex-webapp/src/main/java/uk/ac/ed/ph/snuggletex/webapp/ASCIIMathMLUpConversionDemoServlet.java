/* $Id:FullLaTeXInputDemoServlet.java 158 2008-07-31 10:48:14Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.webapp;

import uk.ac.ed.ph.snuggletex.SerializationSpecifier;
import uk.ac.ed.ph.snuggletex.SnuggleEngine;
import uk.ac.ed.ph.snuggletex.SnuggleInput;
import uk.ac.ed.ph.snuggletex.SnuggleSession;
import uk.ac.ed.ph.snuggletex.WebPageOutputOptions;
import uk.ac.ed.ph.snuggletex.WebPageOutputOptionsTemplates;
import uk.ac.ed.ph.snuggletex.WebPageOutputOptions.WebPageType;
import uk.ac.ed.ph.snuggletex.upconversion.MathMLUpConverter;
import uk.ac.ed.ph.snuggletex.utilities.MathMLUtilities;

import java.io.IOException;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.xml.transform.Transformer;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * Servlet demonstrating the up-conversion process on MathML generated
 * dynamically in the browser by ASCIIMathML.
 * 
 * @author  David McKain
 * @version $Revision:158 $
 */
public final class ASCIIMathMLUpConversionDemoServlet extends BaseServlet {
    
    private static final long serialVersionUID = 2261754980279697343L;

    /** Logger so that we can log what users are trying out to allow us to improve things */
    private static Logger logger = LoggerFactory.getLogger(ASCIIMathMLUpConversionDemoServlet.class);
    
    /** Location of XSLT controlling page layout */
    private static final String DISPLAY_XSLT_LOCATION = "classpath:/asciimathml-upconversion-demo.xsl";
    
    /** Generates initial input form with some demo JavaScript. */
    @Override
    protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException,
            IOException {
        generatePage(request, response, true, null, null, null, null, null, null);
    }
    
    /** Handles the posted raw input & PMathML extracted from ASCIIMathML. */
    @Override
    protected void doPost(HttpServletRequest request, HttpServletResponse response)
            throws ServletException, IOException {
        /* Get the raw ASCIIMathML input and Presentation MathML created by the ASCIIMathML
         * JavaScript code.
         */
        String asciiMathInput = request.getParameter("asciiMathInput");
        String asciiMathOutput = request.getParameter("asciiMathML");
        if (asciiMathInput==null || asciiMathOutput==null) {
            logger.warn("Could not extract data from ASCIIMath: asciiMathInput={}, asciiMathOutput={}", asciiMathInput, asciiMathOutput);
            response.sendError(HttpServletResponse.SC_BAD_REQUEST, "Could not extract data passed by ASCIIMathML");
            return;
        }
        asciiMathOutput = asciiMathOutput.trim();
        
        /* Do up-conversion and extract wreckage */
        MathMLUpConverter upConverter = new MathMLUpConverter(getStylesheetManager());
        SerializationSpecifier sourceSerializationOptions = createMathMLSourceSerializationOptions();
        Document upConvertedMathDocument = upConverter.upConvertASCIIMathML(asciiMathOutput, null);
        Element mathElement = upConvertedMathDocument.getDocumentElement(); /* NB: Document is <math/> here */
        String parallelMathML = MathMLUtilities.serializeElement(mathElement, sourceSerializationOptions);
        String pMathMLUpConverted = MathMLUtilities.serializeDocument(MathMLUtilities.isolateFirstSemanticsBranch(mathElement), sourceSerializationOptions);
        Document cMathMLDocument = MathMLUtilities.isolateAnnotationXML(mathElement, MathMLUpConverter.CONTENT_MATHML_ANNOTATION_NAME);
        String cMathML = cMathMLDocument!=null ? MathMLUtilities.serializeDocument(cMathMLDocument, sourceSerializationOptions) : null;
        String maximaInput = MathMLUtilities.extractAnnotationString(mathElement, MathMLUpConverter.MAXIMA_ANNOTATION_NAME);
        
        logger.info("ASCIIMathML Input: {}", asciiMathInput);
        logger.info("Raw ASCIIMathML Output: {}", asciiMathOutput);
        logger.info("Final parallel MathML: {}", parallelMathML);
        
        generatePage(request, response, false, asciiMathInput, mathElement,
                parallelMathML, pMathMLUpConverted, cMathML, maximaInput);
    }
    
    private void generatePage(HttpServletRequest request, HttpServletResponse response,
            boolean isNewForm, String asciiMathInput, Element parallelMathMLElement,
            String parallelMathML, String pMathMLUpConverted,
            String cMathML, String maximaInput)
            throws IOException, ServletException {
        /* Generate final page using the same process as other demos, which is a bit
         * cheaty here but saves rewriting code. In this case, we'll pass an empty
         * input to SnuggleTeX and ignore the output it gives!
         */
        WebPageOutputOptions options = WebPageOutputOptionsTemplates.createWebPageOptions(WebPageType.CROSS_BROWSER_XHTML);
        options.setIndenting(true);
        options.setIncludingStyleElement(false);
        
        SnuggleEngine engine = createSnuggleEngine();
        SnuggleSession session = engine.createSession();
        session.parseInput(new SnuggleInput("", "Dummy Input"));
        
        /* Create XSLT to generate the resulting page */
        Transformer viewStylesheet = getStylesheet(request, DISPLAY_XSLT_LOCATION);
        viewStylesheet.setParameter("is-mathml-capable", isMathMLCapable(request));
        viewStylesheet.setParameter("is-internet-explorer", isInternetExplorer(request));
        viewStylesheet.setParameter("is-new-form", Boolean.valueOf(isNewForm));
        if (!isNewForm) {
            viewStylesheet.setParameter("ascii-math-input", asciiMathInput);
            viewStylesheet.setParameter("parallel-mathml-element", parallelMathMLElement);
            viewStylesheet.setParameter("parallel-mathml", parallelMathML);
            viewStylesheet.setParameter("pmathml-upconverted", pMathMLUpConverted);
            viewStylesheet.setParameter("cmathml", cMathML);
            viewStylesheet.setParameter("maxima-input", maximaInput);
        }
        options.setStylesheets(viewStylesheet);
        
        /* Generate and serve the resulting web page */
        try {
            session.writeWebPage(options, response, response.getOutputStream());
        }
        catch (Exception e) {
            throw new ServletException("Unexpected Exception", e);
        }
    }
}