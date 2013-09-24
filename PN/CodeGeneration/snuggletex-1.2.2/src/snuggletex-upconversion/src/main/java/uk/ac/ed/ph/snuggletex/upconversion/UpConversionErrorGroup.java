/* $Id:ErrorCode.java 179 2008-08-01 13:41:24Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.upconversion;

import uk.ac.ed.ph.snuggletex.ErrorGroup;
import uk.ac.ed.ph.snuggletex.SnugglePackage;
import uk.ac.ed.ph.snuggletex.upconversion.internal.UpConversionPackageDefinitions;

/**
 * Enumerates the various {@link ErrorGroup}s for the up-conversion module.
 * 
 * @since 1.2.0
 * 
 * @author  David McKain
 * @version $Revision:179 $
 */
public enum UpConversionErrorGroup implements ErrorGroup {
    
    UCF(), /* Presentation to Content MathML Failures */
    UMF(), /* Content MathML to Maxima Failures */
    UAE(), /* Options and Assumption errors */
    
    ;
    
    public String getName() {
        return name();
    }

    public SnugglePackage getPackage() {
        return UpConversionPackageDefinitions.getPackage();
    }
}
