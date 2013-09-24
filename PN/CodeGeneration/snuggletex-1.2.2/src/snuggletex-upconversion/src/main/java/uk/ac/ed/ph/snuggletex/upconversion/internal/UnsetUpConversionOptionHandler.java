/* $Id: UnsetUpConversionOptionHandler.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.upconversion.internal;

import uk.ac.ed.ph.snuggletex.internal.DOMBuilder;
import uk.ac.ed.ph.snuggletex.internal.SnuggleParseException;
import uk.ac.ed.ph.snuggletex.tokens.CommandToken;
import uk.ac.ed.ph.snuggletex.upconversion.IllegalUpconversionOptionException;
import uk.ac.ed.ph.snuggletex.upconversion.UpConversionOptions;
import uk.ac.ed.ph.snuggletex.upconversion.UpConversionUtilities;

import org.w3c.dom.Element;

/**
 * Handler for the <tt>\\unsetUpConversionOption{name}</tt> command.
 * 
 * @since 1.2.0
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
public final class UnsetUpConversionOptionHandler extends UpConversionHandlerBase {
    
    public void handleCommand(DOMBuilder builder, Element parentElement, CommandToken token)
            throws SnuggleParseException {
        /* Get property name */
        String propertyName = builder.extractStringValue(token.getArguments()[0]);
        
        /* Make the change */
        UpConversionOptions options = ensureGetAuthorUpconversionOptions(builder);
        try {
            options.clearSpecifiedOption(propertyName);
        }
        catch (IllegalUpconversionOptionException e) {
            builder.appendOrThrowError(parentElement, token, e.getErrorCode(), (Object[]) e.getArguments());
            return;
        }
        
        /* Now output all current assumptions for the XSLT to use */
        UpConversionUtilities.appendUpConversionOptionsElement(builder.getDocument(), parentElement, options, false);
    }
}
