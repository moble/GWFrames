/* $Id: TestFileHelper.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.testutil;

import uk.ac.ed.ph.snuggletex.SnuggleRuntimeException;
import uk.ac.ed.ph.snuggletex.internal.util.IOUtilities;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collection;

/**
 * Trivial helper class to read in test file data
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
public final class TestFileHelper {

    /**
     * Reads in the given "single line" test file, assuming it is of the format
     * <pre>
     * single line input
     * 1 or more lines of output
     * ==== (divider token, at least 4 characters)
     * ...
     * </pre>
     * 
     * Returns a List of [input,output] pairs.
     * 
     * @throws Exception
     */
    public static Collection<String[]> readAndParseSingleLineInputTestResource(String resourceName) throws Exception {
        String testData = ensureGetResource(resourceName);
        testData = testData.replaceAll("(?m)^#.*$(\\s+)(^|$)", "");
        String[] testItems = testData.split("(?m)\\s*^={4,}\\s*");
        Collection<String[]> result = new ArrayList<String[]>(testItems.length);
        for (String testItem : testItems) {
            result.add(testItem.split("\n+", 2));
        }
        return result;
    }
    
    /**
     * Reads in the given "multiple line" test file, assuming it is of the format
     * <pre>
     * 1 or more lines of output input
     * ---- (divider token, at least 4 characters)
     * 1 or more lines of output
     * ==== (divider token, at least 4 characters)
     * ...
     * </pre>
     * 
     * Returns a List of [input,output] pairs.
     * 
     * @throws Exception
     */
    public static Collection<String[]> readAndParseMultiLineInputTestResource(String resourceName) throws Exception {
        String testData = ensureGetResource(resourceName);
        testData = testData.replaceAll("(?m)^#.*$(\\s+)(^|$)", "");
        String[] testItems = testData.split("(?m)\\s*^={4,}\\s*");
        Collection<String[]> result = new ArrayList<String[]>(testItems.length);
        for (String testItem : testItems) {
            result.add(testItem.split("(?m)\\s*-{4,}\\s*", 2));
        }
        return result;
    }
    
    private static String ensureGetResource(String resourceName) throws IOException {
        InputStream resourceStream = TestFileHelper.class.getClassLoader().getResourceAsStream(resourceName);
        if (resourceStream==null) {
            throw new SnuggleRuntimeException("Could not load Resource '" + resourceName
                    + "' via ClassLoader - check the ClassPath!");
        }
        return IOUtilities.readUnicodeStream(resourceStream);
    }
}
