/* $Id:MathTests.java 179 2008-08-01 13:41:24Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.upconversion;

import uk.ac.ed.ph.snuggletex.AbstractGoodMathTest;
import uk.ac.ed.ph.snuggletex.AbstractGoodTest;
import uk.ac.ed.ph.snuggletex.DOMOutputOptions;
import uk.ac.ed.ph.snuggletex.MathTests;
import uk.ac.ed.ph.snuggletex.SnuggleEngine;
import uk.ac.ed.ph.snuggletex.SnuggleSession;
import uk.ac.ed.ph.snuggletex.testutil.TestFileHelper;
import uk.ac.ed.ph.snuggletex.upconversion.internal.UpConversionPackageDefinitions;
import uk.ac.ed.ph.snuggletex.utilities.MathMLUtilities;

import java.util.Collection;
import java.util.List;
import java.util.logging.Logger;

import junit.framework.Assert;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * Same idea as {@link MathTests}, but tests the initial up-conversion to more
 * semantic Presentation MathML.
 * 
 * @author  David McKain
 * @version $Revision:179 $
 */
@RunWith(Parameterized.class)
public class MathUpConversionToMaximaTests extends AbstractGoodTest {
    
    public static final String TEST_RESOURCE_NAME = "math-upconversion-maxima-tests.txt";
    
    private static final Logger log = Logger.getLogger(MathUpConversionToMaximaTests.class.getName());
    
    private final String expectedMaxima;
    private final UpConvertingPostProcessor upconverter;
    
    public MathUpConversionToMaximaTests(final String inputFragment, final String expectedMaxima) {
        super(inputFragment.endsWith("$") ? inputFragment : "$" + inputFragment + "$");
        this.expectedMaxima = expectedMaxima;
        
        /* Set up up-converter */
        UpConversionOptions upConversionOptions = new UpConversionOptions();
        upConversionOptions.setSpecifiedOption(UpConversionOptionDefinitions.DO_CONTENT_MATHML_NAME, "true");
        upConversionOptions.setSpecifiedOption(UpConversionOptionDefinitions.DO_MAXIMA_NAME, "true");
        upconverter = new UpConvertingPostProcessor(upConversionOptions);
    }
    
    @Parameters
    public static Collection<String[]> data() throws Exception {
        return TestFileHelper.readAndParseSingleLineInputTestResource(TEST_RESOURCE_NAME);
    }
    
    @Test
    public void runTest() throws Throwable {
        String maximaAnnotation = null;
        Element mathmlElement = null;
        String errorCodes = null;
        try {
            /* Run usual process, expecting success */
            Document resultDocument = runSnuggleProcessSuccessfully();
            
            /* Result should be of the correct form */
            mathmlElement = AbstractGoodMathTest.extractMathElement(resultDocument);
            
            /* Get any up-conversion errors, if found */
            List<UpConversionFailure> upConversionFailures = UpConversionUtilities.extractUpConversionFailures(resultDocument);
            StringBuilder errorCodeBuilder = new StringBuilder();
            for (UpConversionFailure error : upConversionFailures) {
                if (errorCodeBuilder.length()!=0) {
                    errorCodeBuilder.append(", ");
                }
                errorCodeBuilder.append(error.getErrorCode());
            }
            errorCodes = errorCodeBuilder.toString();
            if (upConversionFailures.isEmpty()) {
                /* Should have succeeded... */
                /* Extract Maxima annotation */
                maximaAnnotation = MathMLUtilities.extractAnnotationString(mathmlElement, MathMLUpConverter.MAXIMA_ANNOTATION_NAME);
                
                /* Compare with expected */
                Assert.assertEquals(expectedMaxima, maximaAnnotation);
            }
            else {
                /* Make sure we get the correct error code(s) */
                if (expectedMaxima.charAt(0)!='!') {
                    Assert.fail("Did not expect up-conversion errors!");
                }
                String[] expectedErrorCodes = expectedMaxima.substring(1).split(",\\s*");
                Assert.assertEquals(expectedErrorCodes.length, upConversionFailures.size());
                for (int i=0; i<expectedErrorCodes.length; i++) {
                    Assert.assertEquals(expectedErrorCodes[i], upConversionFailures.get(i).getErrorCode().toString());
                }
            }
        }
        catch (Throwable e) {
            log.severe("Input was: " + inputLaTeX);
            log.severe("Resulting MathML was " + (mathmlElement!=null ? MathMLUtilities.serializeElement(mathmlElement) : null));
            log.severe("Resulting Maxima annotation was: " + maximaAnnotation);
            log.severe("Resulting Error codes were:      " + errorCodes);
            log.severe("Expected result would have been: " + expectedMaxima);
            throw e;
        }
    }

    @Override
    protected DOMOutputOptions createDOMOutputOptions() {
        DOMOutputOptions result = super.createDOMOutputOptions();
        result.setDOMPostProcessors(upconverter);
        return result;
    }
    
    @Override
    protected SnuggleSession createSnuggleSession() {
        SnuggleEngine engine = new SnuggleEngine();
        engine.addPackage(UpConversionPackageDefinitions.getPackage());
        
        return engine.createSession();
    }
}
