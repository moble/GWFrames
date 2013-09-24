/* $Id: AssumptionExample.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.upconversion.samples;

import uk.ac.ed.ph.snuggletex.SnuggleEngine;
import uk.ac.ed.ph.snuggletex.SnuggleInput;
import uk.ac.ed.ph.snuggletex.SnuggleSession;
import uk.ac.ed.ph.snuggletex.XMLStringOutputOptions;
import uk.ac.ed.ph.snuggletex.upconversion.UpConversionOptionDefinitions;
import uk.ac.ed.ph.snuggletex.upconversion.UpConversionOptions;
import uk.ac.ed.ph.snuggletex.upconversion.UpConvertingPostProcessor;
import uk.ac.ed.ph.snuggletex.upconversion.IllegalUpconversionOptionException;
import uk.ac.ed.ph.snuggletex.upconversion.internal.UpConversionPackageDefinitions;

import java.io.IOException;

/**
 * Basic example of up-converting some simple LaTeX input to Content MathML and Maxima forms.
 * 
 * <h2>Running Notes</h2>
 * 
 * You will need the following in your ClassPath:
 * 
 * <ul>
 *   <li>snuggletex-core.jar</li>
 *   <li>snuggletex-upconversion.jar</li>
 *   <li>saxon9.jar, saxon9-dom.jar</li> (These are required as the conversion process uses XSLT 2.0)
 * </ul>
 * 
 * @since 1.2.0
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
public final class AssumptionExample {
    
    public static void main(String[] args) throws IOException, IllegalUpconversionOptionException {
        /* We will up-convert this LaTeX input */
        String input = "\\setUpConversionOption{doContentMathML}{true} \\assumeSymbol{f}{function} $f(x)+\\{1,2\\}$";
        
        /* Set up SnuggleEngine, remembering to register package providing up-conversion support */
        SnuggleEngine engine = new SnuggleEngine();
        engine.addPackage(UpConversionPackageDefinitions.getPackage());
        
        /* Create session in usual way */
        SnuggleSession session = engine.createSession();
        
        /* Parse input. I won't bother checking it here */
        session.parseInput(new SnuggleInput(input));
        
        /* Create suitable default options */
        UpConversionOptions defaultOptions = new UpConversionOptions();
        defaultOptions.setSpecifiedOption(UpConversionOptionDefinitions.DO_CONTENT_MATHML_NAME, "true");
        defaultOptions.setSpecifiedOption(UpConversionOptionDefinitions.DO_MAXIMA_NAME, "true");
        defaultOptions.setSpecifiedOption(UpConversionOptionDefinitions.ADD_OPTIONS_ANNOTATION_NAME, "true");
        
        /* Create UpConvertingPostProcessor to perform the up-conversion. */
        UpConvertingPostProcessor upConverter = new UpConvertingPostProcessor(defaultOptions);
        
        /* Create simple XML String output */
        XMLStringOutputOptions xmlStringOutputOptions = new XMLStringOutputOptions();
        xmlStringOutputOptions.addDOMPostProcessors(upConverter);
        xmlStringOutputOptions.setIndenting(true);
        xmlStringOutputOptions.setUsingNamedEntities(true);
        
        /* Do the up-conversion process */
        String result = session.buildXMLString(xmlStringOutputOptions);
        System.out.println("Up-Conversion process generated: " + result);
        System.out.println("Errors: " + session.getErrors());
    }
}
