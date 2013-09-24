/* $Id: UpConversionUtilities.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.upconversion;

import static uk.ac.ed.ph.snuggletex.SnuggleConstants.SNUGGLETEX_NAMESPACE;
import static uk.ac.ed.ph.snuggletex.upconversion.UpConversionOptionDefinitions.OPTION_DEFINITIONS;

import uk.ac.ed.ph.snuggletex.ErrorCode;
import uk.ac.ed.ph.snuggletex.SnuggleConstants;
import uk.ac.ed.ph.snuggletex.SnuggleLogicException;
import uk.ac.ed.ph.snuggletex.internal.util.ConstraintUtilities;
import uk.ac.ed.ph.snuggletex.internal.util.XMLUtilities;
import uk.ac.ed.ph.snuggletex.upconversion.UpConversionOptionDefinitions.OptionValueDefinition;
import uk.ac.ed.ph.snuggletex.utilities.MathMLUtilities;
import uk.ac.ed.ph.snuggletex.utilities.MessageFormatter;

import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 * Some static utility methods for the up-conversion process, including methods for extracting
 * errors.
 * 
 * @since 1.1.0
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
public final class UpConversionUtilities {
    
    /** 
     * Local name of XML container element when writing out options using
     * {@link #appendUpConversionOptionsElement(Document, Element, UpConversionOptions, boolean)}
     */
    public static final String UPCONVERSION_OPTIONS_XML_LOCAL_NAME = "upconversion-options";
    
    /**
     * Returns a full error message for the given {@link UpConversionFailure}, using
     * the SnuggleTeX {@link MessageFormatter} class to do the hard work.
     * 
     * @param failure
     */
    public static String getErrorMessage(UpConversionFailure failure) {
        return MessageFormatter.getErrorMessage(failure.getErrorCode(), failure.getArguments());
    }
    
    /**
     * Returns a full error message for the given {@link UpConversionErrorCode} String and arguments,
     * using the SnuggleTeX {@link MessageFormatter} class to do the hard work.
     */
    public static String getErrorMessage(String upConversionErrorCodeAsString, Object... arguments) {
        /* Look up ErrorCode in enum */
        ErrorCode errorCode;
        try {
            errorCode = UpConversionErrorCode.valueOf(upConversionErrorCodeAsString);
        }
        catch (IllegalArgumentException e) {
            throw new SnuggleLogicException("Could not look up UpConversionErrorCode " + upConversionErrorCodeAsString);
        }
        return MessageFormatter.getErrorMessage(errorCode, arguments);
    }
    
    /**
     * Populates and returns an {@link UpConversionFailure} for the given {@link Element},
     * which is assumed to be a <tt>s:fail</tt> element having the appropriate structure.
     * 
     * @param sFailElement
     */
    public static UpConversionFailure extractUpConversionFailure(final Element sFailElement) {
        ConstraintUtilities.ensureNotNull(sFailElement, "sFailElement");
        if (!(SNUGGLETEX_NAMESPACE.equals(sFailElement.getNamespaceURI()) && sFailElement.getLocalName().equals("fail"))) {
            throw new IllegalArgumentException("Element is not an <s:fail/> element");
        }
        /* Yay! This is one of our <s:fail/> elements. First extract ErrorCode */
        String codeAttribute = sFailElement.getAttribute("code");
        ErrorCode errorCode;
        try {
            errorCode = UpConversionErrorCode.valueOf(codeAttribute);
        }
        catch (IllegalArgumentException e) {
            throw new SnuggleLogicException("Error code '" + codeAttribute + "' not defined");
        }
        /* Now get arguments and context */
        NodeList childNodes = sFailElement.getChildNodes();
        Node child;
        List<String> arguments = new ArrayList<String>();
        String xPath = null;
        String context = null;
        for (int i=0, size=childNodes.getLength(); i<size; i++) {
            child = childNodes.item(i);
            if (child.getNodeType()==Node.ELEMENT_NODE) {
                if (!SNUGGLETEX_NAMESPACE.equals(child.getNamespaceURI())) {
                    throw new SnuggleLogicException("Didn't expect child of <s:fail/> in namespace " + child.getNamespaceURI());
                }
                if (child.getLocalName().equals("arg")) {
                    arguments.add(XMLUtilities.extractTextElementValue((Element) child));
                }
                else if (child.getLocalName().equals("xpath")) {
                    if (context!=null) {
                        throw new SnuggleLogicException("Did not expect more than 1 <s:xpath/> element inside <s:fail/>");
                    }
                    xPath = XMLUtilities.extractTextElementValue((Element) child);
                }
                else if (child.getLocalName().equals("context")) {
                    if (context!=null) {
                        throw new SnuggleLogicException("Did not expect more than 1 <s:context/> element inside <s:fail/>");
                    }
                    context = MathMLUtilities.serializeElement((Element) child);
                }
                else {
                    throw new SnuggleLogicException("Didn't expect child of <s:fail/> with local name " + child.getLocalName());
                }
            }
        }
        if (context==null) {
            throw new SnuggleLogicException("No <s:context/> element found inside <s:fail/>");
        }
        return new UpConversionFailure(errorCode, xPath, context, arguments.toArray(new String[arguments.size()]));
    }
    
    public static List<UpConversionFailure> extractUpConversionFailures(Document upConvertedDocument) {
        return extractUpConversionFailures(upConvertedDocument.getDocumentElement());
    }
    
    public static List<UpConversionFailure> extractUpConversionFailures(Element startSearchElement) {
        List<UpConversionFailure> result = new ArrayList<UpConversionFailure>();
        walkDOM(startSearchElement, result);
        return result;
    }

    /**
     * Checks to see whether the given element is of the form:
     * 
     * <![CDATA[
     * <s:fail code="...">
     *   <s:arg>...</s:arg>
     *   ...
     *   <s:context>...</s:context>
     * </s:fail>
     * ]]>
     * 
     * If so, creates an {@link UpConversionFailure} for this and adds to the result. If not,
     * traverses downwards recursively.
     * 
     * @param searchElement
     * @param resultBuilder
     */
    private static void walkDOM(Element searchElement, List<UpConversionFailure> resultBuilder) {
        if (SNUGGLETEX_NAMESPACE.equals(searchElement.getNamespaceURI()) && searchElement.getLocalName().equals("fail")) {
            resultBuilder.add(extractUpConversionFailure(searchElement));
        }
        else {
            /* Descend into child elements */
            NodeList childNodes = searchElement.getChildNodes();
            Node child;
            for (int i=0, size=childNodes.getLength(); i<size; i++) {
                child = childNodes.item(i);
                if (child.getNodeType()==Node.ELEMENT_NODE) {
                    walkDOM((Element) child, resultBuilder);
                }
            }
        }
    }
    
    //-------------------------------------------------------
    
    /**
     * Writes out the given {@link UpConversionOptions} Object as a child Element of the
     * given element, using a simple ad-hoc XML format.
     * <p>
     * This is used internally to pass Java-specified options to the stylesheets that perform
     * the up-conversion process and may be useful in other situations as well.
     * 
     * @since 1.2.0
     * 
     * @param document Document to append to
     * @param containerElement Element that will be the parent of the resulting DOM fragment
     * @param options {@link UpConversionOptions} to write, may be null.
     * @param applyDefaults if true, then default values for each option are written out
     *   when not explicitly specified or if {@link UpConversionOptions} is null.
     */
    public static void appendUpConversionOptionsElement(Document document, Element containerElement, 
            final UpConversionOptions options, final boolean applyDefaults) {
        Element optionsContainer = appendSnuggleElement(document, containerElement, UpConversionUtilities.UPCONVERSION_OPTIONS_XML_LOCAL_NAME);

        if (options==null) {
            if (applyDefaults) {
                /* Write out default option values */
                for (Entry<String, OptionValueDefinition> entry : OPTION_DEFINITIONS.entrySet()) {
                    String name = entry.getKey();
                    String defaultValue = entry.getValue().getDefaultValue();
                    Element optionElement = appendSnuggleElement(document, optionsContainer, "option");
                    optionElement.setAttribute("name", name);
                    optionElement.setAttribute("value", defaultValue);
                }
                /* (There are no default symbol assumptions) */
            }
            else {
                /* (Nothing to do) */
            }
        }
        else {
            /* Write out specified options, substituting defaults if requested */
            for (String name : OPTION_DEFINITIONS.keySet()) {
                if (applyDefaults || options.isOptionSpecified(name)) {
                    String value = options.getEffectiveOptionValue(name, applyDefaults);
                    Element optionElement = appendSnuggleElement(document, optionsContainer, "option");
                    optionElement.setAttribute("name", name);
                    optionElement.setAttribute("value", value);
                }
            }
            /* Write out specified Symbol assumptions */
            for (ElementWrapper elementWrapper : options.getAssumedSymbols()) {
                String assumptionType = options.getSymbolAssumptionType(elementWrapper);
                
                Element assumeElement = appendSnuggleElement(document, optionsContainer, "symbol");
                assumeElement.setAttribute("assume", assumptionType);
                Node assumptionTargetCopy = elementWrapper.getSymbolElement().cloneNode(true);
                assumeElement.appendChild(assumptionTargetCopy);
            }
        }
    }
    
    private static Element appendSnuggleElement(Document document, Element parentElement, String localName) {
        return (Element) parentElement.appendChild(document.createElementNS(SnuggleConstants.SNUGGLETEX_NAMESPACE,
                "s:" + localName));
    }
}
