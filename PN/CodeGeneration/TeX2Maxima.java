import uk.ac.ed.ph.snuggletex.SnugglePackage;
import uk.ac.ed.ph.snuggletex.SnuggleInput;
import uk.ac.ed.ph.snuggletex.SnuggleEngine;
import uk.ac.ed.ph.snuggletex.SnuggleSession;
import uk.ac.ed.ph.snuggletex.XMLStringOutputOptions;
import uk.ac.ed.ph.snuggletex.WebPageOutputOptions;
import uk.ac.ed.ph.snuggletex.upconversion.UpConvertingPostProcessor;
import uk.ac.ed.ph.snuggletex.upconversion.internal.UpConversionPackageDefinitions;
import uk.ac.ed.ph.snuggletex.utilities.MathMLUtilities;

import static uk.ac.ed.ph.snuggletex.definitions.Globals.MATH_MODE_ONLY;
import static uk.ac.ed.ph.snuggletex.definitions.Globals.PARA_MODE_ONLY;
import static uk.ac.ed.ph.snuggletex.definitions.Globals.TEXT_MODE_ONLY;
import static uk.ac.ed.ph.snuggletex.definitions.LaTeXMode.LR;
import static uk.ac.ed.ph.snuggletex.definitions.LaTeXMode.MATH;
import static uk.ac.ed.ph.snuggletex.definitions.LaTeXMode.PARAGRAPH;
import static uk.ac.ed.ph.snuggletex.definitions.LaTeXMode.VERBATIM;
import static uk.ac.ed.ph.snuggletex.definitions.TextFlowContext.ALLOW_INLINE;
import static uk.ac.ed.ph.snuggletex.definitions.TextFlowContext.IGNORE;
import static uk.ac.ed.ph.snuggletex.definitions.TextFlowContext.START_NEW_XHTML_BLOCK;
import uk.ac.ed.ph.snuggletex.dombuilding.MathEnvironmentHandler;


import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.w3c.dom.Element;

import java.io.IOException;
import java.io.File;

public class TeX2Maxima {
    
    public static void main(String[] args) throws IOException {
	// Basic object
	SnuggleEngine engine = new SnuggleEngine();
	
	// Add ability to convert to Maxima
	engine.addPackage(UpConversionPackageDefinitions.getPackage());
	
	// Inexplicably, snuggletex does not handle equation
	// environments, say we need to add it explicitly.  I found
	// the exemplar for the following code in
	// CorePackageDefinitions.java.
	SnugglePackage localPackage = new SnugglePackage("localPackage");
	localPackage.addEnvironment("equation", TEXT_MODE_ONLY, MATH, null, new MathEnvironmentHandler(), ALLOW_INLINE);
	localPackage.addEnvironment("equation*", TEXT_MODE_ONLY, MATH, null, new MathEnvironmentHandler(), ALLOW_INLINE);
	engine.addPackage(localPackage);
	
	// Loop over all files on the command line
	for(int i_arg=0; i_arg<args.length; ++i_arg) {
	    SnuggleSession session = engine.createSession();
	    System.out.println("\nProcessing "+args[i_arg]+"\n=====================");
	    File file = new File(args[i_arg]);
	    SnuggleInput input = new SnuggleInput(file);
	    session.parseInput(input);
	    
	    // Make sure we use the Maxima converter
	    WebPageOutputOptions webPageOutputOptions = new WebPageOutputOptions();
	    UpConvertingPostProcessor upConverter = new UpConvertingPostProcessor();
	    webPageOutputOptions.addDOMPostProcessors(upConverter);
	    webPageOutputOptions.setIndenting(true);
	    webPageOutputOptions.setUsingNamedEntities(true);
	    
	    // Create the DOM document and search for <math> elements
	    Document doc = session.createWebPage(webPageOutputOptions);
	    NodeList mathElements = doc.getElementsByTagName("math");
	    int length = mathElements.getLength();
	    for(int i = 0; i < length; ++i) {
		String result = MathMLUtilities.extractAnnotationString((Element) mathElements.item(i), "Maxima");
		if(result!=null) {
		    System.out.println("'" + result + "'");
		} else {
		    System.out.println("Failed to parse '" + mathElements.item(i) + "'");
		}
	    }
	    System.out.println("\n");
	}
    }
}
