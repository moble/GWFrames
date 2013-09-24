/* $Id:ErrorCode.java 179 2008-08-01 13:41:24Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.upconversion;

import uk.ac.ed.ph.snuggletex.ErrorCode;
import uk.ac.ed.ph.snuggletex.ErrorGroup;

/**
 * Enumerates the various types of client-induced errors that can arise when using the SnuggleTeX
 * Up-Conversion process.
 * 
 * @since 1.2.0
 * 
 * @author  David McKain
 * @version $Revision:179 $
 */
public enum UpConversionErrorCode implements ErrorCode {
    
    /* ================ Failures in the Up-conversion Process, via extension (U) =============== */
    
    /* Presentation to Content MathML Failures */
    UCFG00(UpConversionErrorGroup.UCF),
    UCFG01(UpConversionErrorGroup.UCF),
    UCFG02(UpConversionErrorGroup.UCF),
    UCFG03(UpConversionErrorGroup.UCF),
    UCFOP0(UpConversionErrorGroup.UCF),
    UCFOP1(UpConversionErrorGroup.UCF),
    UCFOP2(UpConversionErrorGroup.UCF),
    UCFOP3(UpConversionErrorGroup.UCF),
    UCFOP4(UpConversionErrorGroup.UCF),
    UCFOP5(UpConversionErrorGroup.UCF),
    UCFFN0(UpConversionErrorGroup.UCF),
    UCFFN1(UpConversionErrorGroup.UCF),
    UCFFN2(UpConversionErrorGroup.UCF),
    UCFFX0(UpConversionErrorGroup.UCF),
    UCFFX1(UpConversionErrorGroup.UCF),
    UCFFX2(UpConversionErrorGroup.UCF),
    
    /* Content MathML to Maxima Failures */
    UMFG00(UpConversionErrorGroup.UMF),
    UMFG01(UpConversionErrorGroup.UMF),
    UMFG02(UpConversionErrorGroup.UMF),
    UMFG03(UpConversionErrorGroup.UMF),
    UMFG04(UpConversionErrorGroup.UMF),
    UMFFX0(UpConversionErrorGroup.UMF),
    UMFOP0(UpConversionErrorGroup.UMF),
    
    /* Options and Assumption errors */
    UAEOP0(UpConversionErrorGroup.UAE),
    UAEOP1(UpConversionErrorGroup.UAE),
    UAEOP2(UpConversionErrorGroup.UAE),
    UAESY0(UpConversionErrorGroup.UAE),
    UAESY1(UpConversionErrorGroup.UAE),
    UAESY2(UpConversionErrorGroup.UAE),
    
    ;

    private final UpConversionErrorGroup errorGroup;
    
    private UpConversionErrorCode(UpConversionErrorGroup errorGroup) {
        this.errorGroup = errorGroup;
    }
    
    public String getName() {
        return name();
    }
    
    public ErrorGroup getErrorGroup() {
        return errorGroup;
    }
}
