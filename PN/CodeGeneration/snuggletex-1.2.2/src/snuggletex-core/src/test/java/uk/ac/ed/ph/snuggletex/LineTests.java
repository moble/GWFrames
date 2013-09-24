/* $Id:LineTests.java 179 2008-08-01 13:41:24Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex;

import uk.ac.ed.ph.snuggletex.definitions.W3CConstants;
import uk.ac.ed.ph.snuggletex.testutil.TestFileHelper;

import java.util.Collection;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;
import org.w3c.dom.Document;

/**
 * Set of simple tests that read in data from <tt>{@link #TEST_RESOURCE_NAME}</tt>.
 * These all take one input line, parse and compare with XML, which
 * may be specified on multiple input lines for convenience.
 * 
 * @author  David McKain
 * @version $Revision:179 $
 */
@RunWith(Parameterized.class)
public class LineTests extends AbstractGoodXMLTest {
    
    public static final String TEST_RESOURCE_NAME = "line-tests.txt";
    
    @Parameters
    public static Collection<String[]> data() throws Exception {
        return TestFileHelper.readAndParseSingleLineInputTestResource(TEST_RESOURCE_NAME);
    }
    
    public LineTests(final String inputLaTeX, final String expectedXML) {
        super(inputLaTeX,
                "<body xmlns='" + W3CConstants.XHTML_NAMESPACE + "'>"
                + expectedXML.replaceAll("(?m)^ +", "").replaceAll("(?m) +$", "").replace("\n", "")
                + "</body>"
        );
    }
    
    @Override
    protected void fixupDocument(Document document) {
        /* Nothing to do */
    }
    
    @Override
    @Test
    public void runTest() throws Throwable {
        super.runTest();
    }
}
