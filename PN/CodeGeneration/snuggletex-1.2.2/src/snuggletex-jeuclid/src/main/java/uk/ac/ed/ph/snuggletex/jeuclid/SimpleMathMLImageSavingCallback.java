/* $Id: SimpleMathMLImageSavingCallback.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.jeuclid;

import java.io.File;
import java.io.OutputStream;

import net.sourceforge.jeuclid.MutableLayoutContext;
import net.sourceforge.jeuclid.context.LayoutContextImpl;
import net.sourceforge.jeuclid.context.Parameter;

/**
 * Partial convenience implementation of {@link MathMLImageSavingCallback} that
 * assumes that the same Content Type will be used to produce each MathML image and
 * restricts the number of configurable features somewhat.
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
public abstract class SimpleMathMLImageSavingCallback implements MathMLImageSavingCallback {
    
    public static String DEFAULT_FONT_SIZE = "16.0";
    public static String DEFAULT_CONTENT_TYPE = "image/png";
    public static boolean DEFAULT_ANTI_ALIASING = true;
    
    private String fontSize;
    private boolean antiAliasing;
    private String imageContentType;
    
    private final LayoutContextImpl layoutContext;
    
    public SimpleMathMLImageSavingCallback() {
        this.layoutContext = new LayoutContextImpl(LayoutContextImpl.getDefaultLayoutContext());
        setFontSize(DEFAULT_FONT_SIZE);
        setImageContentType(DEFAULT_CONTENT_TYPE);
        setAntiAliasing(DEFAULT_ANTI_ALIASING);
    }
    
    public String getFontSize() {
        return fontSize;
    }

    public void setFontSize(String fontSize) {
        this.fontSize = fontSize;
        this.layoutContext.setParameter(Parameter.MATHSIZE, Float.valueOf(fontSize));
    }


    public boolean isAntiAliasing() {
        return antiAliasing;
    }

    public void setAntiAliasing(boolean antiAliasing) {
        this.antiAliasing = antiAliasing;
        this.layoutContext.setParameter(Parameter.ANTIALIAS, Boolean.valueOf(antiAliasing));
    }

    
    public String getImageContentType() {
        return imageContentType;
    }
    
    public void setImageContentType(String imageType) {
        this.imageContentType = imageType;
    }
    
    //----------------------------------------------------
    
    public final String getImageContentType(int mathmlCounter) {
        return imageContentType;
    }
    
    public final MutableLayoutContext getLayoutContext(int mathmlCounter) {
        return layoutContext;
    }
    
    //----------------------------------------------------
    
    /** @see MathMLImageSavingCallback#getImageOutputFile(int) */
    public abstract File getImageOutputFile(int mathmlCounter);
    
    /** @see MathMLImageSavingCallback#getImageOutputStream(int) */
    public abstract OutputStream getImageOutputStream(int mathmlCounter);
    
    /** @see MathMLImageSavingCallback#getImageURL(int) */
    public abstract String getImageURL(int mathmlCounter);

    /**
     * (A "do nothing" implementation for convenience.)
     * 
     * @see MathMLImageSavingCallback#imageSavingSucceeded(Object, int, String)
     */
    public void imageSavingSucceeded(Object imageFileOrOutputStream, int mathmlCounter,
            String contentType) {
        /* (Do nothing by default) */
    }
}
