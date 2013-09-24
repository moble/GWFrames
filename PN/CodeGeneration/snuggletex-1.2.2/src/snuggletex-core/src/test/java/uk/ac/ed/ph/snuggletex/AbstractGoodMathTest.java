/* $Id:MathTests.java 179 2008-08-01 13:41:24Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex;

import uk.ac.ed.ph.snuggletex.definitions.W3CConstants;
import uk.ac.ed.ph.snuggletex.utilities.MathMLUtilities;

import java.util.logging.Logger;

import junit.framework.Assert;
import junit.framework.AssertionFailedError;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 * Base class for tests which read in LaTeX and parse it in Math mode, expecting to get a single
 * <math/> element in the resulting DOM.
 * 
 * @author  David McKain
 * @version $Revision:179 $
 */
public abstract class AbstractGoodMathTest extends AbstractGoodXMLTest {
    
    private static final Logger log = Logger.getLogger(AbstractGoodMathTest.class.getName());
    
    public AbstractGoodMathTest(final String inputLaTeXMaths, final String expectedMathMLContent) {
        super("$" + inputLaTeXMaths + "$",
            "<math xmlns='" + W3CConstants.MATHML_NAMESPACE + "'>"
            + expectedMathMLContent.replaceAll("(?m)^\\s+", "").replaceAll("(?m)\\s+$", "").replace("\n", "")
            + "</math>");
    }
    
    @Override
    protected void fixupDocument(Document document) {
        extractMathElement(document);
    }
    
    public static Element extractMathElement(Document document) {
        /* Should only have 1 child of doc root (<body/>) element here, which should be <math/>.
         * We'll make that the new root Node */
        try {
            Node rootElement = document.getChildNodes().item(0);
            
            Element firstMathElement = null;
            NodeList childNodes = rootElement.getChildNodes();
            for (int i=0, length=childNodes.getLength(); i<length; i++) {
                Node childNode = childNodes.item(i);
                if (MathMLUtilities.isMathMLElement(childNode, "math")) {
                    if (firstMathElement!=null) {
                        Assert.fail("Found more than one <math/> children");
                    }
                    firstMathElement = (Element) childNode;
                }
                else if (childNode.getNodeType()==Node.ELEMENT_NODE) {
                    Assert.fail("Found unexpected element under root");
                }
                else if (childNode.getNodeType()==Node.TEXT_NODE && childNode.getNodeValue().matches("\\S")) {
                    Assert.fail("Found non-whitespace text Node");
                }
            }
            if (firstMathElement==null) {
                Assert.fail("No <math/> child found");
            }
            document.removeChild(rootElement);
            document.appendChild(firstMathElement);
            
            return firstMathElement;
        }
        catch (AssertionFailedError e) {
            log.severe("Resulting DOM Document did not have expected structure. Got:\n"
                    + MathMLUtilities.serializeDocument(document));
            throw e;
        }
    }
}
