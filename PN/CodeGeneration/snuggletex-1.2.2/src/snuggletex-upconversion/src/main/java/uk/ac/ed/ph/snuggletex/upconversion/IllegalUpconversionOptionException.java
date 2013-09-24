/* $Id: IllegalUpconversionOptionException.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.upconversion;

import uk.ac.ed.ph.snuggletex.utilities.MessageFormatter;

/**
 * This Exception is thrown when an attempt is made to set an illegal assertion,
 * option or value using the {@link UpConversionOptions} API.
 * <p>
 * These contain the same underlying {@link UpConversionErrorCode}s as the errors that would
 * be obtained if setting the same options using the equivalent LaTeX macros.
 * 
 * @since 1.2.0
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
public final class IllegalUpconversionOptionException extends IllegalArgumentException {
    
    private static final long serialVersionUID = 784979608771641983L;

    /** Underlying error code */
    private final UpConversionErrorCode errorCode;
    
    /** Additional arguments, used to formulate a full error message */
    private final String[] arguments;

    public IllegalUpconversionOptionException(final UpConversionErrorCode errorCode, final String... arguments) {
        super(MessageFormatter.getErrorMessage(errorCode, (Object[]) arguments));
        this.errorCode = errorCode;
        this.arguments = arguments;
    }

    public UpConversionErrorCode getErrorCode() {
        return errorCode;
    }
    
    public String[] getArguments() {
        return arguments;
    }
}
