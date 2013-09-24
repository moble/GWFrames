/* $Id: MathMLImageSavingCallback.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.jeuclid;

import java.io.File;
import java.io.OutputStream;

import net.sourceforge.jeuclid.MutableLayoutContext;

import org.w3c.dom.Document;

/**
 * Trivial callback interface used by {@link JEuclidMathMLPostProcessor} to determine
 * where to save each generating MathML image and the URL to use in the resulting
 * HTML as it traverses a given DOM {@link Document}.
 * <p>
 * You should implement this to do whatever is appropriate for your needs.
 * <p>
 * This can be used independently of SnuggleTeX.
 * 
 * @see JEuclidMathMLPostProcessor
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
public interface MathMLImageSavingCallback {
    
    /**
     * Implement to return the desired MIME type for the resulting image.
     * <p>
     * Consult the JEuclid documentation for what is supported; you may need to add additional
     * libraries for certain MIME types.
     * <p>
     * SnuggleTeX contains everything required to support PNG. Note that conversion to GIF
     * does not generally work well so should be avoided...!
     * 
     * @param mathmlCounter identifies the position of the image within the document being processed,
     *   which can be used to ensure unique file names.
     */
    String getImageContentType(int mathmlCounter);
    
    /**
     * Implement if you want the resulting image to be written to a {@link File} of your choice.
     * Return null if you want the image to be written to the result of
     * {@link #getImageOutputStream(int)} instead.
     * 
     * @param mathmlCounter identifies the position of the image within the document being processed,
     *   which can be used to ensure unique file names.
     */
    File getImageOutputFile(int mathmlCounter);
    
    /**
     * Implement this if you want the resulting image to be written to an {@link OutputStream} of
     * your choice. The method {@link #getImageOutputFile(int)} is checked first, so must return
     * null in this case.
     * 
     * @param mathmlCounter identifies the position of the image within the document being processed,
     *   which can be used to ensure unique file names.
     */
    OutputStream getImageOutputStream(int mathmlCounter);
    
    /**
     * Implement to return the URL String that will be put into the <tt>img</tt> <tt>src</tt>
     * attribute to refer to the image.
     * 
     * @param mathmlCounter identifies the position of the image within the document being processed,
     *   which can be used to ensure unique file names.
     */
    String getImageURL(int mathmlCounter);

    /**
     * Implement to fill in the JEuclid {@link MutableLayoutContext} specifying how you want to
     * render this image.
     * 
     * @param mathmlCounter identifies the position of the image within the document being processed,
     *   which can be used to ensure unique file names.
     */
    MutableLayoutContext getLayoutContext(int mathmlCounter);
    
    /**
     * Called back once a MathML image has been saved successfully. Implementors can do anything they
     * need to do at this point.
     * 
     * @param imageFileOrOutputStream {@link File} or {@link OutputStream} written
     * @param mathmlCounter identifies the position of the image within the document being processed,
     *   which can be used to ensure unique file names.
     * @param contentType content type of the saved image File
     */
    void imageSavingSucceeded(Object imageFileOrOutputStream, int mathmlCounter, String contentType);
    
    /**
     * Called back if MathML image could not be saved for some reason.
     * 
     * @param imageFileOrOutputStream {@link File} or {@link OutputStream} that would have
     *   been written if successful
     * @param mathmlCounter identifies the position of the image within the document being processed,
     *   which can be used to ensure unique file names.
     * @param contentType content type of the saved image File
     * @param exception cause of the failure, which may be null
     */
    void imageSavingFailed(Object imageFileOrOutputStream, int mathmlCounter, String contentType,
            Throwable exception);
}