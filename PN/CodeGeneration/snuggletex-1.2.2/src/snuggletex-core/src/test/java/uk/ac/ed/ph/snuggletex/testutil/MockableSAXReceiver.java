/* $Id: MockableSAXReceiver.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.testutil;

import org.xml.sax.Attributes;
import org.xml.sax.ContentHandler;

/**
 * Altered version of the standard SAX {@link ContentHandler} interface that simplifies
 * the interface for handling text Nodes by assuming that incoming text (and whitespace)
 * is coalesced.
 * <p>
 * This makes construction of mock Objects easier.
 * 
 * @see EasyMockContentHandler
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
interface MockableSAXReceiver {

    void startDocument();
    void skippedEntity(String name);
    void processingInstruction(String target, String data);
    void startElement(String uri, String localName, String qName, Attributes atts);
    void coalescedText(String text);
    void coalescedIgnorableWhitespace(String text);
    void endElement(String uri, String localName, String qName);
    void endDocument();

}
