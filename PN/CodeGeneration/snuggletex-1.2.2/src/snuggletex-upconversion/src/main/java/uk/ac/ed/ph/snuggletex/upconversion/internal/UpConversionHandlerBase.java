/* $Id: UpConversionHandlerBase.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.upconversion.internal;

import static uk.ac.ed.ph.snuggletex.upconversion.UpConversionOptionDefinitions.OPTIONS_VARIABLE_NAME;
import static uk.ac.ed.ph.snuggletex.upconversion.UpConversionOptionDefinitions.UPCONVERSION_VARIABLE_NAMESPACE;

import uk.ac.ed.ph.snuggletex.dombuilding.CommandHandler;
import uk.ac.ed.ph.snuggletex.internal.DOMBuilder;
import uk.ac.ed.ph.snuggletex.internal.SnuggleParseException;
import uk.ac.ed.ph.snuggletex.internal.VariableManager;
import uk.ac.ed.ph.snuggletex.tokens.CommandToken;
import uk.ac.ed.ph.snuggletex.upconversion.UpConversionErrorCode;
import uk.ac.ed.ph.snuggletex.upconversion.UpConversionOptions;
import uk.ac.ed.ph.snuggletex.utilities.MathMLUtilities;

import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

/**
 * Base handler for the various assume-based commands.
 * 
 * @since 1.2.0
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
abstract class UpConversionHandlerBase implements CommandHandler {
    
    protected UpConversionOptions ensureGetAuthorUpconversionOptions(DOMBuilder builder) {
        VariableManager variableManager = builder.getVariableManager();
        UpConversionOptions options = (UpConversionOptions) variableManager.getVariable(UPCONVERSION_VARIABLE_NAMESPACE, OPTIONS_VARIABLE_NAME);
        if (options==null) {
            options = new UpConversionOptions();
            variableManager.setVariable(UPCONVERSION_VARIABLE_NAMESPACE, OPTIONS_VARIABLE_NAME, options);
        }
        return options;
    }
    
    protected Element ensureLegalSymbolTarget(DOMBuilder builder, Element parentElement, CommandToken token,
            NodeList targetNodeList)
            throws SnuggleParseException {
        /* Target must be a single Element Node... */
        if (targetNodeList.getLength()!=1 && !MathMLUtilities.isMathMLElement(targetNodeList.item(0))) {
            /* Error: unsupported symbol construct */
            builder.appendOrThrowError(parentElement, token, UpConversionErrorCode.UAESY1);
            return null;
        }
        /* ...and either a <mn>, <mi> or <msub> with similar content */
        Element targetElement = (Element) targetNodeList.item(0);
        String localName = targetElement.getLocalName();
        if ("msub".equals(localName)) {
            if (ensureLegalSymbolTarget(builder, parentElement, token, targetElement.getChildNodes())==null) {
                return null;
            }
        }
        else if (!("mn".equals(localName) || "mi".equals(localName))) {
            /* Error: unsupported symbol construct */
            builder.appendOrThrowError(parentElement, token, UpConversionErrorCode.UAESY1);
            return null;
        }
        return targetElement;
    }
}
