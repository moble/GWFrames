/* $Id: JEuclidUtilities.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.jeuclid;

import uk.ac.ed.ph.snuggletex.DOMOutputOptions;
import uk.ac.ed.ph.snuggletex.DOMPostProcessor;
import uk.ac.ed.ph.snuggletex.DownConvertingPostProcessor;
import uk.ac.ed.ph.snuggletex.WebPageOutputOptions;
import uk.ac.ed.ph.snuggletex.WebPageOutputOptionsTemplates;
import uk.ac.ed.ph.snuggletex.WebPageOutputOptions.WebPageType;

/**
 * Some utility methods for using the JEuclid-based "MathML to Image" conversion
 * functionality.
 *
 * @see WebPageOutputOptionsTemplates
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
public final class JEuclidUtilities {
    
    /**
     * Takes an existing {@link DOMOutputOptions} and configures it to convert
     * MathML to images using JEuclid, optionally down-converting simple expressions to XHTML
     * beforehand.
     * <p>
     * This works by <strong>replacing</strong> any existing {@link DOMPostProcessor}s. If you
     * want to use other {@link DOMPostProcessor}s, then you will have to work out whether they
     * fit in with this process and configure things manually.
     * 
     * @param options existing {@link DOMOutputOptions} Object
     * @param downConvertFirst
     * @param callback
     */
    public static void setupJEuclidPostProcessors(DOMOutputOptions options,
            boolean downConvertFirst, MathMLImageSavingCallback callback) {
        if (downConvertFirst) {
            options.setDOMPostProcessors(
                    new DownConvertingPostProcessor(),
                    new JEuclidMathMLPostProcessor(callback)
            );
        }
        else {
            options.setDOMPostProcessors(new JEuclidMathMLPostProcessor(callback));
        }
    }
    
    /**
     * Creates a new {@link WebPageOutputOptions} suitably configured for converting MathML
     * to images, with optional down-conversion.
     * <p>
     * Note that the resulting {@link WebPageType} of this will be {@link WebPageType#PROCESSED_HTML}.
     * 
     * @param downConvertFirst
     * @param callback
     */
    public static WebPageOutputOptions createWebPageOptions(boolean downConvertFirst, MathMLImageSavingCallback callback) {
        WebPageOutputOptions options = WebPageOutputOptionsTemplates.createWebPageOptions(WebPageType.PROCESSED_HTML);
        setupJEuclidPostProcessors(options, downConvertFirst, callback);
        return options;
    }
}