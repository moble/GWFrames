/* $Id: UpConversionFailure.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.upconversion;

import uk.ac.ed.ph.snuggletex.ErrorCode;
import uk.ac.ed.ph.snuggletex.InputError;

import java.io.Serializable;
import java.util.Arrays;

/**
 * This is the equivalent of an {@link InputError} for the up-conversion process,
 * indicating that a particular island of Presentation MathML created by SnuggleTeX
 * could not be up-converted to either Content MathML or Maxima input form. This
 * will be due to one of the following reasons:
 * <ul>
 *   <li>
 *     The MathML falls outside the scope of what is supported by the up-conversion
 *     process. (I.e. too complex)
 *   </li>
 *   <li>
 *     The MathML is deemed to make no sense. (This is, of course, subjective as it
 *     assumes the MathML is entered in a particular common form so one might reasonably
 *     argue that this is just a special case of the last reason!)
 *   </li>
 * </ul>
 * 
 * See {@link UpConversionUtilities} for methods to help with this, including a method
 * for converting these into error message Strings.
 * 
 * @since 1.1.0
 * 
 * @see UpConversionUtilities
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
public final class UpConversionFailure implements Serializable {
    
    private static final long serialVersionUID = -7888377912814662998L;
    
    /** Error code */
    private final ErrorCode errorCode;
    
    /** Any additional arguments about the error. These are interpolated into error messages */
    private final Object[] arguments;
    
    /** 
     * XPath String corresponding to the first error Element. This is given relative to the
     * "top" of the original MathML (normally a <math/> Element but might be a <semantics/> or
     * <annotation-xml/> and assumes that MathML namespace is the default.)
     * 
     * @since 1.2.0
     */
    private final String xPath;
    
    /** Context within the MathML that the error occurred, as a serialized DOM subtree */
    private final String context;
    
    public UpConversionFailure(final ErrorCode errorCode, final String xPath, final String context, final Object[] arguments) {
        this.errorCode = errorCode;
        this.xPath = xPath;
        this.context = context;
        this.arguments = arguments;
    }

    
    public ErrorCode getErrorCode() {
        return errorCode;
    }

    public Object[] getArguments() {
        return arguments;
    }
    
    public String getXPath() {
        return xPath;
    }

    public String getContext() {
        return context;
    }
    
    @Override
    public String toString() {
        return getClass().getName()
            + "(errorCode=" + errorCode
            + ",arguments=" + Arrays.toString(arguments)
            + ",xPath=" + xPath
            + ",context=" + context
            + ")";
    }
}