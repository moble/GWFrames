/* $Id: AbstractGoodXMLTest.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex;

import static org.easymock.EasyMock.createStrictControl;

import uk.ac.ed.ph.snuggletex.internal.util.XMLUtilities;
import uk.ac.ed.ph.snuggletex.testutil.EasyMockContentHandler;

import java.io.StringReader;
import java.io.StringWriter;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.sax.SAXResult;
import javax.xml.transform.stream.StreamResult;

import org.easymock.IMocksControl;
import org.w3c.dom.Document;
import org.xml.sax.InputSource;
import org.xml.sax.XMLReader;

/**
 * Base for tests which take LaTeX input, parse it and compare the resulting XML
 * against the specified output.
 * 
 * @author  David McKain
 * @version $Revision: 525 $
 */
public abstract class AbstractGoodXMLTest extends AbstractGoodTest {
    
    private static final Logger log = Logger.getLogger(AbstractGoodXMLTest.class.getName());
    
    protected final String expectedXML;
    
    /**
     * Subclasses should call this to initialise {@link #inputLaTeX} and {@link #expectedXML}
     * as appropriate.
     */
    protected AbstractGoodXMLTest(String inputLaTeX, String expectedXML) {
        super(inputLaTeX);
        this.expectedXML = expectedXML;
    }
    
    protected abstract void fixupDocument(Document document);
    
    protected boolean showTokensOnFailure() {
        return true;
    }
    
    protected void verifyResultDocument(TransformerFactory transformerFactory, Document resultDocument) throws Throwable {
        /* Create mock to handle the SAX streams */
        IMocksControl control = createStrictControl();
        EasyMockContentHandler saxControl = new EasyMockContentHandler(control);
        
        /* Fire expected output at handler */
        SAXParserFactory parserFactory = SAXParserFactory.newInstance();
        parserFactory.setNamespaceAware(true);
        SAXParser parser = parserFactory.newSAXParser();
        XMLReader reader = parser.getXMLReader();
        InputSource inputSource = new InputSource(new StringReader(expectedXML));
        reader.setContentHandler(saxControl);
        reader.parse(inputSource);
        
        /* Now replay and fire actual resulting XML to mock as SAX stream */
        control.replay();
        saxControl.replay();
        
        /* Finally verify everything */
        Transformer serializer = transformerFactory.newTransformer();
        serializer.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION, "yes");
        serializer.transform(new DOMSource(resultDocument), new SAXResult(saxControl));
        control.verify();
        saxControl.verify();
    }
    
    public void runTest() throws Throwable {
        String output=null;
        try {
            Document resultDocument = runSnuggleProcessSuccessfully();
    
            /* Let subclass fudge up the resulting document if required */
            fixupDocument(resultDocument);
            
            /* Serialize the output */
            StringWriter outputWriter = new StringWriter();
            TransformerFactory transformerFactory = XMLUtilities.createJAXPTransformerFactory();
            Transformer serializer = transformerFactory.newTransformer();
            serializer.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION, "yes");
            serializer.transform(new DOMSource(resultDocument), new StreamResult(outputWriter));
            output = outputWriter.toString();
            
            /* Verify result */
            verifyResultDocument(transformerFactory, resultDocument);
        }
        catch (Throwable e) {
            log.severe("Input was: " + inputLaTeX);
            if (showTokensOnFailure()) {
                if (rawDump!=null) {
                    log.severe("Raw dump was: " + rawDump);
                }
                if (fixedDump!=null) {
                    log.severe("Fixed dump was: " + fixedDump);
                }
            }
            if (output!=null) {
                log.severe("Expected output: " + expectedXML);
                log.severe("Actual output:   " + output);
            }
            log.log(Level.SEVERE, "Error was: ", e);
            throw e;
        }
    }
}