/* $Id: UnassumeSymbolHandler.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.upconversion.internal;

import uk.ac.ed.ph.snuggletex.internal.DOMBuilder;
import uk.ac.ed.ph.snuggletex.internal.SnuggleParseException;
import uk.ac.ed.ph.snuggletex.internal.DOMBuilder.OutputContext;
import uk.ac.ed.ph.snuggletex.tokens.CommandToken;
import uk.ac.ed.ph.snuggletex.upconversion.IllegalUpconversionOptionException;
import uk.ac.ed.ph.snuggletex.upconversion.UpConversionOptions;
import uk.ac.ed.ph.snuggletex.upconversion.UpConversionUtilities;

import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

/**
 * Handler for the <tt>\\unassumeSymbol{target}</tt> command.
 * 
 * @since 1.2.0
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
public final class UnassumeSymbolHandler extends UpConversionHandlerBase {
    
    public void handleCommand(DOMBuilder builder, Element parentElement, CommandToken token)
            throws SnuggleParseException {
        /* First argument is the symbol that this assumption will apply to, which is
         * a blob of MathML. This will be "checked" in the XSLT to make sure we can handle
         * this correctly.
         */
        builder.pushOutputContext(OutputContext.MATHML_INLINE);
        NodeList assumptionTargetRaw = builder.extractNodeListValue(token.getArguments()[0]);
        builder.popOutputContext();
        Element assumptionTarget = ensureLegalSymbolTarget(builder, parentElement, token, assumptionTargetRaw);
        if (assumptionTarget==null) {
            return;
        }
        
        /* Make the change */
        UpConversionOptions options = ensureGetAuthorUpconversionOptions(builder);
        try {
            options.unassumeSymbol(assumptionTarget);
        }
        catch (IllegalUpconversionOptionException e) {
            builder.appendOrThrowError(parentElement, token, e.getErrorCode(), (Object[]) e.getArguments());
            return;
 
        }
        
        /* Now output all current assumptions for the XSLT to use */
        UpConversionUtilities.appendUpConversionOptionsElement(builder.getDocument(), parentElement, options, false);
    }
}
