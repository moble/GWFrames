/* $Id: MathMLToImageLinkPostProcessor.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.webapp;

import uk.ac.ed.ph.snuggletex.SnuggleConstants;
import uk.ac.ed.ph.snuggletex.SnuggleRuntimeException;
import uk.ac.ed.ph.snuggletex.definitions.W3CConstants;
import uk.ac.ed.ph.snuggletex.utilities.MathMLPostProcessor;
import uk.ac.ed.ph.snuggletex.utilities.MathMLUtilities;

import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

/**
 * Simple {@link MathMLPostProcessor} that replaces any MathML islands
 * with HTML links calling back on {@link MathInputToImageServlet} to
 * render them as images.
 * 
 * TODO: This might be quite a useful utility so may want to
 * be generalised and moved out of the demo webapp.
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
public final class MathMLToImageLinkPostProcessor extends MathMLPostProcessor {
    
    private final String contextPath;
    
    public MathMLToImageLinkPostProcessor(final String contextPath) {
        this.contextPath = contextPath;
    }
    
    @Override
    protected void handleMathMLIsland(Element inputMathIsland, Document outputDocument,
            Node outputParentNode, int mathmlCounter) {
        /* Extract SnuggleTeX annotation */
        String snuggleInput = MathMLUtilities.extractAnnotationString(inputMathIsland, SnuggleConstants.SNUGGLETEX_MATHML_SOURCE_ANNOTATION_ENCODING);
        if (snuggleInput==null) {
            throw new SnuggleRuntimeException("Expected to find SnuggleTeX annotation inside MathML - fix this process!");
        }
        /* Replace with appropriate image servlet link */
        boolean isBlock = "block".equals(inputMathIsland.getAttribute("display"));
        Element replacement = outputDocument.createElementNS(W3CConstants.XHTML_NAMESPACE, isBlock ? "div" : "span");
        replacement.setAttribute("class", "mathml-math");
        
        Element image = outputDocument.createElementNS(W3CConstants.XHTML_NAMESPACE, "img");
        try {
            image.setAttribute("src", contextPath + "/MathInputToImage.png?input="
                    + URLEncoder.encode(snuggleInput, "UTF-8"));
        }
        catch (UnsupportedEncodingException e) {
            throw new SnuggleRuntimeException("Unexpected Exception", e);
        }
        replacement.appendChild(image);
        outputParentNode.appendChild(replacement);
    }

}
