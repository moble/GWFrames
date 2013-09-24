/* $Id: EasyMockContentHandler.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.testutil;

import static org.easymock.EasyMock.eq;

import org.easymock.IMocksControl;
import org.xml.sax.Attributes;
import org.xml.sax.ContentHandler;
import org.xml.sax.Locator;

/**
 * This SAX {@link ContentHandler} makes it easy to record and verify SAX events.
 * 
 * <h2>Usage</h2>
 * 
 * <ul>
 *   <li>Create an instance of this handler using an existing {@link IMocksControl}</li>
 *   <li>Set it to receive SAX events that you want to "record"</li>
 *   <li>Call {@link #replay()} when you do the same with your {@link IMocksControl}</li>
 *   <li>Fire the real SAX events at this handler</li>
 *   <li>Call {@link #verify()} when you do the same with your {@link IMocksControl}</li>
 * </ul>
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
public final class EasyMockContentHandler implements ContentHandler {
    
    private final MockableSAXReceiver target;
    private final StringBuilder charactersBuilder;
    private final StringBuilder whitespaceBuilder;
    private boolean replayed;
    
    public EasyMockContentHandler(IMocksControl control) {
        this.target = control.createMock(MockableSAXReceiver.class);
        this.charactersBuilder = new StringBuilder();
        this.whitespaceBuilder = new StringBuilder();
        this.replayed = false;
    }
    
    //----------------------------------------------
    
    public void reset() {
        ensureState();
        this.replayed = false;
    }
    
    public void replay() {
        ensureState();
        this.replayed = true;
    }
    
    public void verify() {
        ensureState();
    }
    
    private void ensureState() {
        if (charactersBuilder.length()>0 || whitespaceBuilder.length()>0) {
            throw new IllegalStateException("Unfinished character or whitespace data");
        }        
    }
    
    private void stopBuildingText() {
        if (charactersBuilder.length()>0) {
            target.coalescedText(charactersBuilder.toString());
            charactersBuilder.setLength(0);
        }
        if (whitespaceBuilder.length()>0) {
            target.coalescedIgnorableWhitespace(whitespaceBuilder.toString());
            whitespaceBuilder.setLength(0);
        }
    }
    
    public void startDocument() {
        ensureState();
        target.startDocument();
    }

    public void startElement(String uri, String localName, String qName, Attributes atts) {
        stopBuildingText();
        if (!replayed) {
            target.startElement(eq(uri), eq(localName), eq(qName), AttributesMatcher.eqAtts(atts));
        }
        else {
            target.startElement(uri, localName, qName, atts);
        }
    }

    public void startPrefixMapping(String prefix, String uri) {
        /* Ignore this */
    }

    public void characters(char[] ch, int start, int length) {
        charactersBuilder.append(ch, start, length);
    }

    public void endDocument() {
        ensureState();
        target.endDocument();
    }

    public void endElement(String uri, String localName, String qName) {
        stopBuildingText();
        target.endElement(uri, localName, qName);
    }

    public void endPrefixMapping(String prefix) {
        /* Ignore this */
    }

    public void ignorableWhitespace(char[] ch, int start, int length) {
        whitespaceBuilder.append(ch, start, length);
    }

    public void processingInstruction(String piTarget, String data) {
        stopBuildingText();
        target.processingInstruction(piTarget, data);
    }

    public void setDocumentLocator(Locator locator) {
        /* Ignore this */
    }

    public void skippedEntity(String name) {
        target.skippedEntity(name);
    }
}
