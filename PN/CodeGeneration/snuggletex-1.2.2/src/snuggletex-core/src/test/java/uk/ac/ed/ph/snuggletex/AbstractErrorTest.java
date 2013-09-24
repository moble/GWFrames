/* $Id: AbstractErrorTest.java 564 2010-05-20 10:35:51Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex;

import uk.ac.ed.ph.snuggletex.utilities.MessageFormatter;

import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import junit.framework.Assert;

import org.junit.Test;

/**
 * Base for tests that expect specific errors to be recorded.
 *
 * @author  David McKain
 * @version $Revision: 564 $
 */
public abstract class AbstractErrorTest {
    
    private static final Logger logger = Logger.getLogger(AbstractErrorTest.class.getName());
    
    private final String inputLaTeX;
    private final String expectedErrorCodeStrings;
    
    protected AbstractErrorTest(final String inputLaTeX, final String expectedErrorCodeStrings) {
        this.inputLaTeX = inputLaTeX;
        this.expectedErrorCodeStrings = expectedErrorCodeStrings;
    }
    
    protected SnuggleSession createSnuggleSession() {
        return new SnuggleEngine().createSession();
    }
    
    @Test
    public void runTest() throws Throwable {
        /* Convert %n in input LaTeX to a newline (cheap hack for simple multi-line inputs!) */
        String inputLaTeXLines = inputLaTeX.replace("%n", "\n");
        
        /* Parse document and build XML */
        SnuggleSession session = createSnuggleSession();
        try {
            session.parseInput(new SnuggleInput(inputLaTeXLines));
            session.buildXMLString();
        }
        catch (Exception e) {
            logger.log(Level.WARNING, "Parsing failed unexpectedly on input " + inputLaTeX, e);
            throw e;
        }
        /* Extract list of error codes */
        List<InputError> errors = session.getErrors();
        String[] actualErrorCodes = new String[errors.size()];
        for (int i=0; i<actualErrorCodes.length; i++) {
            actualErrorCodes[i] = errors.get(i).getErrorCode().getName();
        }
        
        /* We allow a single wildcard (...) at the end of the expectedErrorCodeStrings to denote
         * 'any other stuff afterwards'.
         */
        String expectedErrorCodeData = expectedErrorCodeStrings.trim();
        boolean isWildcard = false;
        if (expectedErrorCodeData.endsWith("...")) {
            isWildcard = true;
            expectedErrorCodeData = expectedErrorCodeData.substring(0, expectedErrorCodeData.length()-3);
        }
        String[] expectedErrorCodes = expectedErrorCodeData.split(",\\s*");
        
        /* FINISH ME OFF! */
        
        /* Now check things */
        if (!compare(expectedErrorCodes, actualErrorCodes, isWildcard)) {
            StringBuilder messageBuilder = new StringBuilder("Test failed!\nInput was: ").append(inputLaTeX)
                .append("\nExpected error codes were: ").append(Arrays.toString(expectedErrorCodes))
                .append("\nActual error codes were:   ").append(Arrays.toString(actualErrorCodes));
            for (int i=0; i<errors.size(); i++) {
                messageBuilder.append("\nError ").append(i+1)
                    .append(" was: ")
                    .append(MessageFormatter.formatErrorAsString(errors.get(i)));
            }
            logger.warning(messageBuilder.toString());
            Assert.fail("Failed on " + inputLaTeX);
        }
    }
    
    private boolean compare(String[] expected, String[] actual, boolean expectedWildcard) {
        if (!expectedWildcard) {
            return Arrays.equals(expected, actual);
        }
        /* Wildcard case */
        if (actual.length < expected.length) {
            return false;
        }
        for (int i=0; i<expected.length; i++) {
            if (!expected[i].equals(actual[i])) {
                return false;
            }
        }
        return true;
    }
}
