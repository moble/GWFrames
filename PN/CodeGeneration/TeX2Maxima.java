import uk.ac.ed.ph.snuggletex.SnuggleInput;
import uk.ac.ed.ph.snuggletex.SnuggleEngine;
import uk.ac.ed.ph.snuggletex.SnuggleSession;
import uk.ac.ed.ph.snuggletex.XMLStringOutputOptions;
import uk.ac.ed.ph.snuggletex.WebPageOutputOptions;
import uk.ac.ed.ph.snuggletex.upconversion.UpConvertingPostProcessor;
import uk.ac.ed.ph.snuggletex.upconversion.internal.UpConversionPackageDefinitions;
import uk.ac.ed.ph.snuggletex.utilities.MathMLUtilities;

import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.w3c.dom.Element;

import java.io.IOException;

//// These are needed for the Java 6 version of readFile
import java.io.File;
import java.io.FileInputStream;
import java.nio.channels.FileChannel;
import java.nio.MappedByteBuffer;
import java.nio.charset.Charset;

//// These are needed for the Java 7 version of readFile
// import java.nio.charset.Charset;
// import java.nio.ByteBuffer;
// import java.nio.file.Paths;

// // These are for testing
// import org.w3c.dom.Node;
// import org.w3c.dom.traversal.DocumentTraversal;
// import org.w3c.dom.traversal.NodeFilter;
// import org.w3c.dom.traversal.NodeIterator;

public class TeX2Maxima {
    
    // This is works in Java 6 and greater, which allows support for Mac OS X 10.6
    private static String readFile(String path) throws IOException {
	FileInputStream stream = new FileInputStream(new File(path));
	try {
	    FileChannel fc = stream.getChannel();
	    MappedByteBuffer bb = fc.map(FileChannel.MapMode.READ_ONLY, 0, fc.size());
	    /* Instead of using default, pass in a decoder. */
	    return Charset.defaultCharset().decode(bb).toString();
	}
	finally {
	    stream.close();
	}
    }
    
    //// This only works with Java 7, which is not available on Mac OS X 10.6
    // static String readFile(String path, Charset encoding) throws IOException {
    // 	byte[] encoded = Files.readAllBytes(Paths.get(path));
    // 	return encoding.decode(ByteBuffer.wrap(encoded)).toString();
    // }
    
    public static void main(String[] args) throws IOException {
        /* Create vanilla SnuggleEngine and new SnuggleSession */
        SnuggleEngine engine = new SnuggleEngine();
	engine.addPackage(UpConversionPackageDefinitions.getPackage());
        SnuggleSession session = engine.createSession();
        
        /* Parse some very basic Math Mode input */
	// String inputString = new String(readFile("Test.tex"));
	// System.out.println("inputString = \"\"\"" + inputString + "\"\"\"");
        // SnuggleInput input = new SnuggleInput(inputString);
	File file = new File("Test2.tex");
        SnuggleInput input = new SnuggleInput(file);
        session.parseInput(input);
	UpConvertingPostProcessor upConverter = new UpConvertingPostProcessor();
	
        XMLStringOutputOptions xmlStringOutputOptions = new XMLStringOutputOptions();
        xmlStringOutputOptions.addDOMPostProcessors(upConverter);
        xmlStringOutputOptions.setIndenting(true);
        xmlStringOutputOptions.setUsingNamedEntities(true);
        
	NodeList entries = session.buildDOMSubtree(xmlStringOutputOptions);
	int length = entries.getLength();
	for(int i = 0; i < length; i = i+1) {
	    try {
		System.out.print("(Element) entry " + i + " : " + ((Element) entries.item(i)) + "\n" );
		System.out.println("'" + MathMLUtilities.extractAnnotationString((Element) entries.item(i), "Maxima") + "'");
	    } catch (ClassCastException e) {
		System.out.print("ClassCastException: entry " + i + " : " + entries.item(i) + "\n" );
	    } catch (IllegalArgumentException e) {
		System.out.print("IllegalArgumentException: entry " + i + " : " + entries.item(i) + "\n" );
	    }
	}
	
        WebPageOutputOptions webPageOutputOptions = new WebPageOutputOptions();
        webPageOutputOptions.addDOMPostProcessors(upConverter);
        webPageOutputOptions.setIndenting(true);
        webPageOutputOptions.setUsingNamedEntities(true);
	Document doc = session.createWebPage(webPageOutputOptions);
	NodeList mathElements = doc.getElementsByTagName("math");
	length = mathElements.getLength();
	System.out.print("And now...\n\n");
	for(int i = 0; i < length; i = i+1) {
	    System.out.print("(Element) entry " + i + " : " + ((Element) mathElements.item(i)) + "\n" );
	    System.out.println("'" + MathMLUtilities.extractAnnotationString((Element) mathElements.item(i), "Maxima") + "'");
	}
        
    }
}
