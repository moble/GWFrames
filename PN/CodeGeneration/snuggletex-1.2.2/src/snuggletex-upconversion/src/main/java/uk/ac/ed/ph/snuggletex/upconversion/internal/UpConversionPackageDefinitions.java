/* $Id: UpConversionPackageDefinitions.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.upconversion.internal;

import static uk.ac.ed.ph.snuggletex.definitions.Globals.TEXT_MODE_ONLY;
import static uk.ac.ed.ph.snuggletex.definitions.LaTeXMode.LR;
import static uk.ac.ed.ph.snuggletex.definitions.LaTeXMode.MATH;
import static uk.ac.ed.ph.snuggletex.definitions.TextFlowContext.IGNORE;

import uk.ac.ed.ph.snuggletex.SnugglePackage;
import uk.ac.ed.ph.snuggletex.SnuggleRuntimeException;
import uk.ac.ed.ph.snuggletex.definitions.LaTeXMode;
import uk.ac.ed.ph.snuggletex.upconversion.UpConversionErrorCode;

import java.util.MissingResourceException;
import java.util.ResourceBundle;

/**
 * This defines the {@link SnugglePackage} providing up-conversion functionality.
 * 
 * @since 1.2.0
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
public final class UpConversionPackageDefinitions {
    
    /** Location of {@link ResourceBundle} providing error messages for this bundle */
    public static final String ERROR_MESSAGES_PROPERTIES_BASENAME = "uk/ac/ed/ph/snuggletex/upconversion/error-messages";
    
    private static final SnugglePackage upConversionPackage;
    
    public static final SnugglePackage getPackage() {
        return upConversionPackage;
    }
    
    static {
        upConversionPackage = new SnugglePackage("UpConversion");
        
        /* Set up error messages for this package */
        upConversionPackage.addErrorCodes(UpConversionErrorCode.values());
        try {
            upConversionPackage.setErrorMessageBundle(ResourceBundle.getBundle(ERROR_MESSAGES_PROPERTIES_BASENAME));
        }
        catch (MissingResourceException e) {
            throw new SnuggleRuntimeException(e);
        }
        
        /* Special commands for configuring the up-conversion process. */
        upConversionPackage.addComplexCommand("setUpConversionOption", false, 2, TEXT_MODE_ONLY, new LaTeXMode[] { LR, LR }, new SetUpConversionOptionHandler(), IGNORE);
        upConversionPackage.addComplexCommand("unsetUpConversionOption", false, 1, TEXT_MODE_ONLY, new LaTeXMode[] { LR }, new UnsetUpConversionOptionHandler(), IGNORE);
        upConversionPackage.addComplexCommand("assumeSymbol", false, 2, TEXT_MODE_ONLY, new LaTeXMode[] { MATH, LR }, new AssumeSymbolHandler(), IGNORE);
        upConversionPackage.addComplexCommand("unassumeSymbol", false, 1, TEXT_MODE_ONLY, new LaTeXMode[] { MATH }, new UnassumeSymbolHandler(), IGNORE);
     }
}
