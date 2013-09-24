/* $Id:MathTests.java 179 2008-08-01 13:41:24Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex;

import uk.ac.ed.ph.snuggletex.testutil.TestFileHelper;

import java.util.Collection;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

/**
 * Set of simple maths-based tests that read their data in from <tt>{@link #TEST_RESOURCE_NAME}</tt>.
 * The input is a single line of LaTeX which will be put into <tt>$...$</tt> and parsed
 * then compared with the given multi-line XML. The enclosing <tt>math</tt> element is
 * automatically added to the XML for convenience. See the sample file for examples.
 * 
 * @author  David McKain
 * @version $Revision:179 $
 */
@RunWith(Parameterized.class)
public class MathTests extends AbstractGoodMathTest {
    
    public static final String TEST_RESOURCE_NAME = "math-tests.txt";
    
    @Parameters
    public static Collection<String[]> data() throws Exception {
        return TestFileHelper.readAndParseSingleLineInputTestResource(TEST_RESOURCE_NAME);
    }
    
    public MathTests(final String inputLaTeXMaths, final String expectedMathMLContent) {
        super(inputLaTeXMaths, expectedMathMLContent);
    }
    
    @Override
    @Test
    public void runTest() throws Throwable {
        super.runTest();
    }
}
