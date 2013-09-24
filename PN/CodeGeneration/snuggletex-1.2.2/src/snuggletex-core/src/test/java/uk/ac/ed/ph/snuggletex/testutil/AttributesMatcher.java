/* $Id: AttributesMatcher.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.testutil;

import org.easymock.EasyMock;
import org.easymock.IArgumentMatcher;
import org.xml.sax.Attributes;
import org.xml.sax.helpers.AttributesImpl;

/**
 * EasyMock argument matcher that compares SAX {@link Attributes} Objects.
 * 
 * @author  David McKain
 * @version $Revision: 525 $
 */
public final class AttributesMatcher implements IArgumentMatcher {
    
    private final AttributesImpl expectedAttrs;
    
    public AttributesMatcher(Attributes expected) {
        this.expectedAttrs = new AttributesImpl(expected);
    }
    
    public boolean matches(Object actual) {
        if (!(actual instanceof Attributes)) {
            return false;
        }
        Attributes actualAttrs = (Attributes) actual;
        
        /* Go through each attribute, comparing with an expected one. 
         * We'll ignore 'xmlns' attributes.
         * 
         * We succeed if there is a one-to-one correspondence between the non-xmlns attributes
         * in the list of expected attributes with those in the actual attributes.
         * 
         * Since there's no order here, we need to compare both ways
         */
        int expectedIndex,actualIndex;
        int expectedSize = expectedAttrs.getLength();
        int actualSize = actualAttrs.getLength();
        int expectedMatchCount = 0;
        int actualMatchCount = 0;
        int expectedNonXmlnsCount = 0;
        int actualNonXmlnsCount = 0;
        for (expectedIndex=0; expectedIndex<expectedSize; expectedIndex++) {
            if (expectedAttrs.getQName(expectedIndex).startsWith("xmlns")) {
                continue;
            }
            expectedNonXmlnsCount++;
            for (actualIndex=0; actualIndex<actualSize; actualIndex++) {
                if (actualAttrs.getQName(actualIndex).equals(expectedAttrs.getQName(expectedIndex))
                        && actualAttrs.getLocalName(actualIndex).equals(expectedAttrs.getLocalName(expectedIndex))
                        && actualAttrs.getType(actualIndex).equals(expectedAttrs.getType(expectedIndex))
                        && actualAttrs.getURI(actualIndex).equals(expectedAttrs.getURI(expectedIndex))
                        && actualAttrs.getValue(actualIndex).equals(expectedAttrs.getValue(expectedIndex))) {
                    expectedMatchCount++;
                }
            }
        }
        for (actualIndex=0; actualIndex<actualSize; actualIndex++) {
            if (actualAttrs.getQName(actualIndex).startsWith("xmlns")) {
                continue;
            }
            actualNonXmlnsCount++;
            for (expectedIndex=0; expectedIndex<expectedSize; expectedIndex++) {
                if (actualAttrs.getQName(actualIndex).equals(expectedAttrs.getQName(expectedIndex))
                        && actualAttrs.getLocalName(actualIndex).equals(expectedAttrs.getLocalName(expectedIndex))
                        && actualAttrs.getType(actualIndex).equals(expectedAttrs.getType(expectedIndex))
                        && actualAttrs.getURI(actualIndex).equals(expectedAttrs.getURI(expectedIndex))
                        && actualAttrs.getValue(actualIndex).equals(expectedAttrs.getValue(expectedIndex))) {
                    actualMatchCount++;
                }
            }
        }
        return expectedMatchCount==actualMatchCount
            && expectedMatchCount==expectedNonXmlnsCount
            && expectedNonXmlnsCount==actualNonXmlnsCount;
    }
    
    public void appendTo(StringBuffer buffer) {
        buffer.append("Attributes(");
        for (int i=0; i<expectedAttrs.getLength(); i++) {
            if (i>0) {
                buffer.append(",");
            }
            buffer.append(expectedAttrs.getQName(i))
                .append("=>")
                .append(expectedAttrs.getValue(i));
        }
        buffer.append(")");
    }
    
    public static <T extends Attributes> T eqAtts(T atts) {
        EasyMock.reportMatcher(new AttributesMatcher(atts));
        return atts;
    }
}