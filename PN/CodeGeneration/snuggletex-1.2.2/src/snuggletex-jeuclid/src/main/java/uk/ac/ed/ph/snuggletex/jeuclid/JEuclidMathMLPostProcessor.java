/* $Id: JEuclidMathMLPostProcessor.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.jeuclid;

import uk.ac.ed.ph.snuggletex.definitions.W3CConstants;
import uk.ac.ed.ph.snuggletex.utilities.MathMLPostProcessor;
import uk.ac.ed.ph.snuggletex.utilities.SnuggleUtilities;

import java.awt.Dimension;
import java.io.File;
import java.io.OutputStream;

import net.sourceforge.jeuclid.MutableLayoutContext;
import net.sourceforge.jeuclid.converter.Converter;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

/**
 * Extension of {@link MathMLPostProcessor} that uses JEuclid to convert MathML elements into
 * XHTML + image replacements.
 * 
 * @see MathMLImageSavingCallback
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
public final class JEuclidMathMLPostProcessor extends MathMLPostProcessor {
    
    private final MathMLImageSavingCallback imageSavingCallback;
    
    public JEuclidMathMLPostProcessor(final MathMLImageSavingCallback callback) {
        this.imageSavingCallback = callback;
    }
    
    @Override
    protected void handleMathMLIsland(final Element inputMathIsland, Document outputDocument,
            Node outputParentNode, final int mathmlCounter) {
        /* We use JEuclid to create image rendition of Node, adding some appropriate XHTML to the
         * outputDocument instead of the original MathML.
         * 
         * First we determine whether we're saving to File (first choice) or OutputStream.
         */
        File imageOutputFile = imageSavingCallback.getImageOutputFile(mathmlCounter);
        Object imageOutputObject = imageOutputFile;
        OutputStream imageOutputStream = null;
        if (imageOutputFile==null) {
            imageOutputStream = imageSavingCallback.getImageOutputStream(mathmlCounter);
            imageOutputObject = imageOutputStream;
            if (imageOutputStream==null) {
                throw new IllegalArgumentException("Both getImageOutputFile() and getImageOutputStream() returned null");
            }
        }
        String contentType = imageSavingCallback.getImageContentType(mathmlCounter);
        MutableLayoutContext layoutContext = imageSavingCallback.getLayoutContext(mathmlCounter);
        Dimension imageDimension = null;
        try {
            /* Call up JEuclid, saving to either File or Stream as determine above */
            Converter converter = Converter.getInstance();
            if (imageOutputFile!=null) {
                imageDimension = converter.convert(inputMathIsland, imageOutputFile, contentType, layoutContext);
            }
            else {
                imageDimension = converter.convert(inputMathIsland, imageOutputStream, contentType, layoutContext);
            }
            if (imageDimension!=null) {
                /* Success - inform callback */
                imageSavingCallback.imageSavingSucceeded(imageOutputObject, mathmlCounter, contentType);
            }
            else {
                /* Failed */
                imageDimension = new Dimension(0, 0);
                imageSavingCallback.imageSavingFailed(imageOutputObject, mathmlCounter, contentType, null);
            }
        }
        catch (Exception e) {
            imageDimension = new Dimension(0, 0);
            imageSavingCallback.imageSavingFailed(imageOutputObject, mathmlCounter, contentType, e);
        }

        /* Next we extract the SnuggleTeX annotation within the MathML element, if applicable, which contains the
         * original LaTeX input for this math region. This is used to create an "alt" attribute.
         */
        String snuggleTeXEncoding = SnuggleUtilities.extractSnuggleTeXAnnotation(inputMathIsland);
        if (snuggleTeXEncoding!=null) {
            snuggleTeXEncoding = snuggleTeXEncoding
                .replaceAll("%\\s+", "") /* Strip LaTeX comments */
                .replaceAll("\\s+", " "); /* Normalise whitespace */
        }
        
        /* Next we add <div> or <span> to the output Document instead of the input <math/> element */
        boolean isBlock = inputMathIsland.getAttribute("display").equals("block");
        Element divOrSpan = outputDocument.createElementNS(W3CConstants.XHTML_NAMESPACE, isBlock ? "div" : "span");
        divOrSpan.setAttribute("class", "mathml-math");
        outputParentNode.appendChild(divOrSpan);
        
        /* Then put an <img/> inside the <div> or <span> */
        Element imgElement = outputDocument.createElementNS(W3CConstants.XHTML_NAMESPACE, "img");
        imgElement.setAttribute("src", imageSavingCallback.getImageURL(mathmlCounter));
        imgElement.setAttribute("width", Integer.toString(imageDimension.width));
        imgElement.setAttribute("height", Integer.toString(imageDimension.height));
        if (snuggleTeXEncoding!=null) {
            imgElement.setAttribute("alt", snuggleTeXEncoding);
        }
        divOrSpan.appendChild(imgElement);
    }
}