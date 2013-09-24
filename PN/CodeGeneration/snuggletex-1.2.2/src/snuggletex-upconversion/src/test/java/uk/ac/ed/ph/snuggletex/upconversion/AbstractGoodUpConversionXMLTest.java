/* $Id:MathTests.java 179 2008-08-01 13:41:24Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.upconversion;

import uk.ac.ed.ph.snuggletex.AbstractGoodMathTest;
import uk.ac.ed.ph.snuggletex.AbstractGoodXMLTest;
import uk.ac.ed.ph.snuggletex.SnuggleEngine;
import uk.ac.ed.ph.snuggletex.SnuggleSession;
import uk.ac.ed.ph.snuggletex.definitions.W3CConstants;
import uk.ac.ed.ph.snuggletex.upconversion.internal.UpConversionPackageDefinitions;

import java.util.List;

import javax.xml.transform.TransformerFactory;

import junit.framework.Assert;

import org.w3c.dom.Document;

/**
 * Base class for up-conversion tests. This is pretty much the same as {@link AbstractGoodMathTest}
 * but only wraps the input in "$...$" delimiters if they don't already end with one of these.
 * This allows assumptions to be set up in advance.
 * 
 * @author  David McKain
 * @version $Revision:179 $
 */
public abstract class AbstractGoodUpConversionXMLTest extends AbstractGoodXMLTest {
    
    public AbstractGoodUpConversionXMLTest(final String inputFragment, final String expectedMathMLContent) {
        super(inputFragment.endsWith("$") ? inputFragment : "$" + inputFragment + "$",
            "<math xmlns='" + W3CConstants.MATHML_NAMESPACE + "'>"
            + expectedMathMLContent.replaceAll("(?m)^\\s+", "").replaceAll("(?m)\\s+$", "").replace("\n", "")
            + "</math>");
    }
    
    @Override
    protected void fixupDocument(Document document) {
        AbstractGoodMathTest.extractMathElement(document);
    }
    

    /**
     * Overridden to cope with up-conversion failures, checking them against the given error code.
     */
    @Override
    protected void verifyResultDocument(TransformerFactory transformerFactory, Document resultDocument) throws Throwable {
        List<UpConversionFailure> upConversionFailures = UpConversionUtilities.extractUpConversionFailures(resultDocument);
        String result;
        if (upConversionFailures.isEmpty()) {
            /* Should have succeeded, so verify as normal */
            super.verifyResultDocument(transformerFactory, resultDocument);
        }
        else {
            /* Make sure we get the correct error code(s) */
            result = expectedXML.replaceAll("<.+?>", ""); /* (Yes, it's not really XML in this case!) */
            if (result.charAt(0)!='!') {
                Assert.fail("Did not expect up-conversion errors!");
            }
            String[] expectedErrorCodes = result.substring(1).split(",\\s*");
            Assert.assertEquals(expectedErrorCodes.length, upConversionFailures.size());
            for (int i=0; i<expectedErrorCodes.length; i++) {
                Assert.assertEquals(expectedErrorCodes[i], upConversionFailures.get(i).getErrorCode().toString());
            }
        }
    }
    
    @Override
    protected SnuggleSession createSnuggleSession() {
        SnuggleEngine engine = new SnuggleEngine();
        engine.addPackage(UpConversionPackageDefinitions.getPackage());
        
        return engine.createSession();
    }
    
    @Override
    protected boolean showTokensOnFailure() {
        return false;
    }
}
