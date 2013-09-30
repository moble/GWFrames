import uk.ac.ed.ph.snuggletex.SnugglePackage;
import uk.ac.ed.ph.snuggletex.SnuggleInput;
import uk.ac.ed.ph.snuggletex.SnuggleEngine;
import uk.ac.ed.ph.snuggletex.SnuggleSession;
import uk.ac.ed.ph.snuggletex.XMLStringOutputOptions;
import uk.ac.ed.ph.snuggletex.WebPageOutputOptions;
import uk.ac.ed.ph.snuggletex.upconversion.UpConvertingPostProcessor;
import uk.ac.ed.ph.snuggletex.upconversion.internal.UpConversionPackageDefinitions;
import uk.ac.ed.ph.snuggletex.utilities.MathMLUtilities;

import static uk.ac.ed.ph.snuggletex.definitions.Globals.ALL_MODES;
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
import uk.ac.ed.ph.snuggletex.dombuilding.EqnArrayHandler;
import uk.ac.ed.ph.snuggletex.dombuilding.DoNothingHandler;
import uk.ac.ed.ph.snuggletex.semantics.TabularInterpretation;

import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.w3c.dom.Element;
import org.w3c.dom.Attr;

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
        TabularInterpretation tabularInterpretation = new TabularInterpretation();
	DoNothingHandler doNothingHandler = new DoNothingHandler();
	localPackage.addEnvironment("equation",  TEXT_MODE_ONLY, MATH, null, new MathEnvironmentHandler(), ALLOW_INLINE);
	localPackage.addEnvironment("equation*", TEXT_MODE_ONLY, MATH, null, new MathEnvironmentHandler(), ALLOW_INLINE);
        localPackage.addEnvironment("align",  PARA_MODE_ONLY, MATH, tabularInterpretation, new EqnArrayHandler(), START_NEW_XHTML_BLOCK);
        localPackage.addEnvironment("align*", PARA_MODE_ONLY, MATH, tabularInterpretation, new EqnArrayHandler(), START_NEW_XHTML_BLOCK);
	localPackage.addComplexCommandSameArgMode("label", false, 1, ALL_MODES, doNothingHandler, IGNORE);
	engine.addPackage(localPackage);
	
	// Loop over all files on the command line
	for(int i_arg=0; i_arg<args.length; ++i_arg) {
	    SnuggleSession session = engine.createSession();
	    String processingString = "Processing "+args[i_arg];
	    System.out.println(processingString+"\n"+ new String(new char[processingString.length()]).replace("\0", "="));
	    File file = new File(args[i_arg]);
	    SnuggleInput input = new SnuggleInput(file);
	    session.parseInput(input);
	    
	    // Make sure we use the Maxima converter
	    WebPageOutputOptions webPageOutputOptions = new WebPageOutputOptions();
	    UpConvertingPostProcessor upConverter = new UpConvertingPostProcessor();
	    webPageOutputOptions.addDOMPostProcessors(upConverter);
	    webPageOutputOptions.setIndenting(true);
	    webPageOutputOptions.setUsingNamedEntities(true);
	    webPageOutputOptions.setAddingMathSourceAnnotations(true);
	    
	    // Create the DOM document and search for <math> elements
	    Document doc = session.createWebPage(webPageOutputOptions);
	    NodeList mathElements = doc.getElementsByTagName("math");
	    int length = mathElements.getLength();
	    for(int i = 0; i < length; ++i) {
		String result = MathMLUtilities.extractAnnotationString((Element) mathElements.item(i), "Maxima");
		if(result!=null) {
		    System.out.println("'" + result + "'");
		} else {
		    try {
			result = MathMLUtilities.extractAnnotationString((Element) mathElements.item(i), "LaTeX");
			try {
			    String error;
			    try {
				NodeList annotationList = MathMLUtilities.extractAnnotationXML((Element) mathElements.item(i),
											       "MathML-Content-upconversion-failures");
				error = ((Element) annotationList.item(0)).getAttribute("message");
			    } catch (Exception e) {
				NodeList annotationList = MathMLUtilities.extractAnnotationXML((Element) mathElements.item(i),
											       "Maxima-upconversion-failures");
				error = ((Element) annotationList.item(0)).getAttribute("message");
			    }
			    System.err.println("Error: '" + error + "' for input\n       " + result);
			} catch (Exception e) {
			    System.err.println("Error: Unkown failure on input\n       '" + result);
			}
		    } catch (Exception e) {
			System.err.println("Don't even know what failed to parse...");
		    }
		}
	    }
	    System.out.println("");
	    
	    // // For debugging, output XML tree
	    // XMLStringOutputOptions options = new XMLStringOutputOptions();
	    // // options.setSerializationMethod(SerializationMethod.XHTML);
	    // options.addDOMPostProcessors(upConverter);
	    // options.setIndenting(true);
	    // options.setEncoding("UTF-8");
	    // options.setAddingMathSourceAnnotations(true);
	    // options.setUsingNamedEntities(true); /* (Only used if caller has an XSLT 2.0 processor) */
	    // System.out.println(session.buildXMLString(options)+"\n");
	}
    }
}
