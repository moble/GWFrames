/* $Id: CommandArgumentsTest.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex;

import static junit.framework.Assert.assertEquals;

import uk.ac.ed.ph.snuggletex.definitions.BuiltinCommand;
import uk.ac.ed.ph.snuggletex.definitions.CoreErrorCode;
import uk.ac.ed.ph.snuggletex.definitions.Globals;
import uk.ac.ed.ph.snuggletex.definitions.TextFlowContext;
import uk.ac.ed.ph.snuggletex.definitions.UserDefinedCommand;
import uk.ac.ed.ph.snuggletex.dombuilding.DoNothingHandler;
import uk.ac.ed.ph.snuggletex.tokens.CommandToken;

import org.junit.Before;
import org.junit.Test;

/**
 * Tests the extraction of arguments on {@link BuiltinCommand}s. (There is no need to do
 * this here for {@link UserDefinedCommand}}s - they are more easily tested using {@link LineTests}
 * or something similar.
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
public final class CommandArgumentsTest {
    
    private SnuggleSession session;
    
    @Before
    public void setup() {
        /* We'll create a custom built-in command called \\bob that simply does nothing of interest */
        SnugglePackage pkg = new SnugglePackage("temp");
        pkg.addComplexCommandSameArgMode("bob", true, 1,
                Globals.TEXT_MODE_ONLY, new DoNothingHandler(), TextFlowContext.ALLOW_INLINE);
        pkg.addComplexCommandSameArgMode("bobopt", true, 0,
                Globals.TEXT_MODE_ONLY, new DoNothingHandler(), TextFlowContext.ALLOW_INLINE);
        
        SnuggleEngine engine = new SnuggleEngine();
        engine.addPackage(pkg);
        
        session = engine.createSession();
    }
    
    @Test
    public void testOptionalArgumentOnBuiltin() throws Exception {
        String input = "\\bob[a]{thing}";
        session.parseInput(new SnuggleInput(input));
        
        CommandToken commandToken = (CommandToken) session.getParsedTokens().get(0);
        assertEquals("a", commandToken.getOptionalArgument().getSlice().extract().toString());
    }
    
    @Test
    public void testUnfinishedOptionalArgumentOnBuiltin() throws Exception {
        String input = "\\bobopt[a";
        session.parseInput(new SnuggleInput(input));
        
        CommandToken commandToken = (CommandToken) session.getParsedTokens().get(0);
        assertEquals("a", commandToken.getOptionalArgument().getSlice().extract().toString());
    }
    
    @Test
    public void testRequiredArgumentOnBuiltin() throws Exception {
        String input = "\\bob[a]{thing}";
        session.parseInput(new SnuggleInput(input));
        
        CommandToken commandToken = (CommandToken) session.getParsedTokens().get(0);
        assertEquals("thing", commandToken.getArguments()[0].getSlice().extract().toString());
    }
    
    @Test
    public void testUnfinishedRequiredArgumentOnBuiltin() throws Exception {
        String input = "\\bob{thing";
        session.parseInput(new SnuggleInput(input));
        
        CommandToken commandToken = (CommandToken) session.getParsedTokens().get(0);
        assertEquals("thing", commandToken.getArguments()[0].getSlice().extract().toString());
        assertEquals(1, session.getErrors().size());
        assertEquals(CoreErrorCode.TTEG00, session.getErrors().get(0).getErrorCode());
    }
}
