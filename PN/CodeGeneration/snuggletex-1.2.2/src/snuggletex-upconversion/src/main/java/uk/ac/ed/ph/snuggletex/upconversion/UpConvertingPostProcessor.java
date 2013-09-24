/* $Id: UpConvertingPostProcessor.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.upconversion;

import uk.ac.ed.ph.snuggletex.DOMOutputOptions;
import uk.ac.ed.ph.snuggletex.DOMPostProcessor;
import uk.ac.ed.ph.snuggletex.utilities.StylesheetManager;

import org.w3c.dom.Document;

/**
 * Implementation of {@link DOMPostProcessor} that bootstraps into {@link MathMLUpConverter},
 * providing the functionality offered by it.
 * 
 * @since 1.1.0
 * 
 * @see MathMLUpConverter
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
public final class UpConvertingPostProcessor implements DOMPostProcessor {
    
    private UpConversionOptions upconversionOptions;
    
    public UpConvertingPostProcessor() {
        this(null);
    }
    
    public UpConvertingPostProcessor(UpConversionOptions upconversionOptions) {
        this.upconversionOptions = upconversionOptions;
    }

    public UpConversionOptions getUpconversionOptions() {
        return upconversionOptions;
    }
    
    public void setUpconversionOptions(UpConversionOptions upconversionOptions) {
        this.upconversionOptions = upconversionOptions;
    }
    
    public Document postProcessDOM(Document workDocument, final DOMOutputOptions unused,
            StylesheetManager stylesheetManager) {
        return new MathMLUpConverter(stylesheetManager)
            .upConvertSnuggleTeXMathML(workDocument, upconversionOptions);
    }
}
