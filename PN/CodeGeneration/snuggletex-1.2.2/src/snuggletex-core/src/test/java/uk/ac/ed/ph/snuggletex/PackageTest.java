/* $Id:DefinitionMapTest.java 179 2008-08-01 13:41:24Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex;

import uk.ac.ed.ph.snuggletex.SnuggleInput;
import uk.ac.ed.ph.snuggletex.SnuggleEngine;
import uk.ac.ed.ph.snuggletex.SnuggleSession;
import uk.ac.ed.ph.snuggletex.definitions.BuiltinCommand;
import uk.ac.ed.ph.snuggletex.definitions.Globals;
import uk.ac.ed.ph.snuggletex.definitions.TextFlowContext;
import uk.ac.ed.ph.snuggletex.dombuilding.DoNothingHandler;
import uk.ac.ed.ph.snuggletex.tokens.FlowToken;

import java.util.List;

import junit.framework.Assert;

import org.junit.Test;

/**
 * Tests the {@link SnugglePackage} class.
 *
 * @author  David McKain
 * @version $Revision:179 $
 */
public class PackageTest {
    
    @Test
    public void testCustomDefinition() throws Exception {
        /* We'll create a custom built-in command called \\bob that simply does nothing of interest */
        SnugglePackage pkg = new SnugglePackage("test");
        BuiltinCommand bobCommand = pkg.addSimpleCommand("bob", Globals.TEXT_MODE_ONLY,
                new DoNothingHandler(), TextFlowContext.ALLOW_INLINE);
        
        SnuggleEngine engine = new SnuggleEngine();
        engine.addPackage(pkg);
        
        SnuggleSession session = engine.createSession();
        session.parseInput(new SnuggleInput("\\bob"));

        /* Verify that we got exactly one command token for '\\bob' */
        List<FlowToken> tokens = session.getParsedTokens();
        Assert.assertEquals(1, tokens.size());
        Assert.assertTrue(tokens.get(0).isCommand(bobCommand));
    }
}
