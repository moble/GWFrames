/* $Id: MathMLUpConverter.java 560 2010-05-19 09:25:55Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.upconversion;

import uk.ac.ed.ph.snuggletex.DOMOutputOptions;
import uk.ac.ed.ph.snuggletex.SnuggleConstants;
import uk.ac.ed.ph.snuggletex.SnuggleLogicException;
import uk.ac.ed.ph.snuggletex.SnuggleRuntimeException;
import uk.ac.ed.ph.snuggletex.SnuggleSession;
import uk.ac.ed.ph.snuggletex.internal.util.XMLUtilities;
import uk.ac.ed.ph.snuggletex.utilities.SaxonTransformerFactoryChooser;
import uk.ac.ed.ph.snuggletex.utilities.SimpleStylesheetCache;
import uk.ac.ed.ph.snuggletex.utilities.StylesheetCache;
import uk.ac.ed.ph.snuggletex.utilities.StylesheetManager;
import uk.ac.ed.ph.snuggletex.utilities.TransformerFactoryChooser;

import java.io.IOException;
import java.io.StringReader;

import javax.xml.transform.Templates;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;
import javax.xml.transform.dom.DOMResult;
import javax.xml.transform.dom.DOMSource;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

/**
 * Standalone utility class for "up-converting" MathML Documents created by either SnuggleTeX
 * or ASCIIMathML as follows:
 * <ul>
 *   <li>Presentation MathML is "enhanced" with extra semantics</li>
 *   <li>Conversion to Content MathML is optionally attempted (assuming PMathML follows certain conventions)</li>
 *   <li>Conversion to Maxima input is optionally attempted (assuming PMathML follows certain conventions)</li>
 * </ul>
 * This can be invoked within the normal SnuggleTeX parsing process using by adding a
 * {@link UpConvertingPostProcessor} to the List returned by
 * {@link DOMOutputOptions#getDOMPostProcessors()}.
 * 
 * <h2>Usage Notes</h2>
 * 
 * <ul>
 *   <li>
 *     If the up-conversion process fails because the input cannot be understood, you can find
 *     out what went wrong by passing the resulting {@link Document} to
 *     {@link UpConversionUtilities#extractUpConversionFailures(Document)}, to get a List
 *     of {@link UpConversionFailure} Objects that can be inspected further.
 *   </li>
 *   <li>An implementation of this class is thread-safe</li>
 *   <li>
 *     If you use XSLT in your own application, consider using the constructor that
 *     takes a {@link StylesheetCache} to integrate with your own caching mechanism.
 *   </li> 
 *     If you don't use XSLT in your own application, you can use the default constructor
 *     that creates a cache used by the instance of this class. In this case, you'll want
 *     to ensure your instance has "application-wide" context to maximise performance.
 *   </li>
 * </ul>
 * 
 * @since 1.1.0
 * 
 * @author  David McKain
 * @version $Revision: 560 $
 */
public class MathMLUpConverter {

    /* (Names of the various annotations. These are also defined in the XSLT so, if changing, we
     * must ensure that things are manually kept in sync. I could pass all of these as parameters
     * but it seems like overkill here!)
     */
    
    public static final String SNUGGLETEX_ANNOTATION_NAME = "SnuggleTeX";
    public static final String LATEX_ANNOTATION_NAME = "LaTeX";
    public static final String CONTENT_MATHML_ANNOTATION_NAME = "MathML-Content";
    public static final String CONTENT_FAILURES_ANNOTATION_NAME = "MathML-Content-upconversion-failures";
    public static final String MAXIMA_ANNOTATION_NAME = "Maxima";
    public static final String MAXIMA_FAILURES_ANNOTATION_NAME = "Maxima-upconversion-failures";
    public static final String ASCIIMATH_INPUT_ANNOTATION_NAME = "ASCIIMathInput";

    /** "Base" location for the XSLT stylesheets used here */
    private static final String UPCONVERTER_BASE_LOCATION = "classpath:/uk/ac/ed/ph/snuggletex/upconversion";
    
    /** Location of the initial XSLT for fixing up ASCIIMathML */
    private static final String ASCIIMATH_FIXER_XSL_LOCATION = UPCONVERTER_BASE_LOCATION + "/asciimathml-fixer.xsl";
    
    /** Location of the main up-converting XSLT */
    private static final String UPCONVERTER_XSL_LOCATION = UPCONVERTER_BASE_LOCATION + "/snuggletex-upconverter.xsl";
    
    /** Helper to manage the {@link StylesheetCache} */
    private final StylesheetManager stylesheetManager;
    
    /**
     * Creates a new up-converter using a simple internal cache, hard-coded to use Saxon for
     * XSLT 2.0 support.
     * <p>
     * Use this constructor if you don't use XSLT yourself. In this case, you'll want your
     * instance of this class to be reused as much as possible to get the benefits of caching.
     */
    public MathMLUpConverter() {
        this(SaxonTransformerFactoryChooser.getInstance(), new SimpleStylesheetCache());
    }
    
    /**
     * Creates a new up-converter using the given {@link StylesheetCache} to cache internal XSLT
     * stylesheets, hard-coded to use Saxon for XSLT 2.0 support.
     * <p>
     * Use this constructor if you do your own XSLT caching that you want to integrate in, or
     * if the default doesn't do what you want.
     */
    public MathMLUpConverter(StylesheetCache stylesheetCache) {
        this(SaxonTransformerFactoryChooser.getInstance(), stylesheetCache);
    }
    
    /**
     * Creates a new up-converter using the given {@link TransformerFactoryChooser} and
     * {@link StylesheetCache} to cache internal XSLT stylesheets.
     * <p>
     * Use this constructor if you do your own XSLT caching that you want to integrate in, or
     * if the default doesn't do what you want and/or if you want to override the choice of 
     * Saxon for XSLT 2.0 support (which is currently a bad idea!)
     * 
     * @param transformerFactoryChooser {@link TransformerFactoryChooser} to use, which must not
     *   be null.
     * @param stylesheetCache {@link StylesheetCache} to use, which may be null to avoid caching.
     */
    public MathMLUpConverter(final TransformerFactoryChooser transformerFactoryChooser, StylesheetCache stylesheetCache) {
        this(new StylesheetManager(transformerFactoryChooser, stylesheetCache));
    }
    
    /**
     * Creates a new up-converter using the given {@link StylesheetManager} for managing XSLT
     * stylesheets and stuff.
     */
    public MathMLUpConverter(StylesheetManager stylesheetManager) {
        this.stylesheetManager = stylesheetManager;
    }
    
    //----------------------------------------------------------------------
    
    /**
     * Up-converts the SnuggleTeX XHTML output Document (assumed to contain islands of MathML),
     * creating a new result {@link Document}.
     * 
     * @param document DOM {@link Document}, assumed to have been generated by
     *   {@link SnuggleSession#buildDOMSubtree()} or similar.
     * @param upConversionOptions {@link UpConversionOptions} to use, or null to
     *   use defaults.
     */
    public Document upConvertSnuggleTeXMathML(final Document document, final UpConversionOptions upConversionOptions) {
        Document resultDocument = XMLUtilities.createNSAwareDocumentBuilder().newDocument();
        try {
            /* Create required XSLT */
            Templates upconverterStylesheet = stylesheetManager.getStylesheet(UPCONVERTER_XSL_LOCATION, true);
            Transformer upconverter = upconverterStylesheet.newTransformer();
            
            /* Pass UpConversionOptions, using default as required */
            upconverter.setParameter("{" + SnuggleConstants.SNUGGLETEX_NAMESPACE + "}global-upconversion-options",
                    createUpConversionOptionsElement(upConversionOptions));
            
            /* Do the transform */
            upconverter.transform(new DOMSource(document), new DOMResult(resultDocument));
        }
        catch (TransformerConfigurationException e) {
            throw new SnuggleRuntimeException("Could not instantiate Transformer from Templates", e);
        }
        catch (TransformerException e) {
            throw new SnuggleLogicException("Up-conversion XSLT transform failed", e);
        }
        return resultDocument;
    }
    
    /**
     * DEVELOPER NOTE: I have included the default value for this parameter within the XSLT as well
     * so that people can reuse it easier in a standalone fashion. Make sure you update the XSLT
     * if any of the defaults change here.
     */
    private Element createUpConversionOptionsElement(final UpConversionOptions upConversionOptions) {
        Document document = XMLUtilities.createNSAwareDocumentBuilder().newDocument();
        Element root = (Element) document.appendChild(document.createElementNS(SnuggleConstants.SNUGGLETEX_NAMESPACE, "root"));
        UpConversionUtilities.appendUpConversionOptionsElement(document, root, upConversionOptions, true);
        return (Element) root.getFirstChild();
    }

    /**
     * Up-converts the ASCIIMathML output represented by the given {@link Document}. This
     * would normally consist of a single MathML element, but can be arbitrary XML containing
     * MathML islands.
     * 
     * @param asciiMathMLDocument DOM {@link Document}, assumed to have been generated by
     *   ASCIIMathML
     * @param upConversionOptions {@link UpConversionOptions} to use, or null to
     *   use defaults.
     */
    public Document upConvertASCIIMathML(final Document asciiMathMLDocument, final UpConversionOptions upConversionOptions) {
        /* First of all we convert the ASCIIMathML into something equivalent to SnuggleTeX output */
        Document fixedDocument = XMLUtilities.createNSAwareDocumentBuilder().newDocument();
        try {
            Templates fixerStylesheet = stylesheetManager.getStylesheet(ASCIIMATH_FIXER_XSL_LOCATION, true);
            fixerStylesheet.newTransformer().transform(new DOMSource(asciiMathMLDocument), new DOMResult(fixedDocument));
        }
        catch (TransformerConfigurationException e) {
            throw new SnuggleRuntimeException("Could not instantiate Transformer from Templates", e);
        }
        catch (TransformerException e) {
            throw new SnuggleLogicException("ASCIIMathML fix-up XSLT transform failed", e);
        }
        /* Then do the normal SnuggleTeX up-conversion */
        return upConvertSnuggleTeXMathML(fixedDocument, upConversionOptions);
    }

    /**
     * Convenience version of {@link #upConvertASCIIMathML(Document, UpConversionOptions)} that
     * takes the raw ASCIIMathML as a String. You would normally extract this from ASCIIMathML
     * using some of our JavaScript helper functions.
     * 
     * @param rawASCIIMathML ASCIIMathML output string
     * @param upConversionOptions {@link UpConversionOptions} to use, or null to
     *   use defaults.
     */
    public Document upConvertASCIIMathML(final String rawASCIIMathML, final UpConversionOptions upConversionOptions) {
        Document inputDocument;
        try {
            inputDocument = XMLUtilities.createNSAwareDocumentBuilder().parse(new InputSource(new StringReader(rawASCIIMathML)));
        }
        catch (SAXException e) {
            throw new SnuggleRuntimeException("Could not parse XML generated by ASCIIMathML: " + rawASCIIMathML, e);
        }
        catch (IOException e) {
            throw new SnuggleRuntimeException("Unexpected Exception reading ASCIIMathML XML String data", e);
        }
        return upConvertASCIIMathML(inputDocument, upConversionOptions);
    }
}
