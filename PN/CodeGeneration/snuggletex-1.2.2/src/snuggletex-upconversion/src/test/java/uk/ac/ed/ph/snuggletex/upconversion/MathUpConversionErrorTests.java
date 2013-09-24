/* $Id:MathTests.java 179 2008-08-01 13:41:24Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.upconversion;

import uk.ac.ed.ph.snuggletex.AbstractErrorTest;
import uk.ac.ed.ph.snuggletex.MathTests;
import uk.ac.ed.ph.snuggletex.SnuggleEngine;
import uk.ac.ed.ph.snuggletex.SnuggleSession;
import uk.ac.ed.ph.snuggletex.testutil.TestFileHelper;
import uk.ac.ed.ph.snuggletex.upconversion.internal.UpConversionPackageDefinitions;

import java.util.Collection;

import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

/**
 * Same idea as {@link MathTests}, but tests the up-conversion to Content
 * MathML.
 * 
 * @author  David McKain
 * @version $Revision:179 $
 */
@RunWith(Parameterized.class)
public class MathUpConversionErrorTests extends AbstractErrorTest {
    
    public static final String TEST_RESOURCE_NAME = "math-upconversion-error-tests.txt";
    
    @Parameters
    public static Collection<String[]> data() throws Exception {
        return TestFileHelper.readAndParseSingleLineInputTestResource(TEST_RESOURCE_NAME);
    }

    public MathUpConversionErrorTests(String inputLaTeX, String expectedErrorCodeStrings) {
        super(inputLaTeX, expectedErrorCodeStrings);
    }
    
    @Override
    protected SnuggleSession createSnuggleSession() {
        SnuggleEngine engine = new SnuggleEngine();
        engine.addPackage(UpConversionPackageDefinitions.getPackage());
        
        return engine.createSession();
    }
}
