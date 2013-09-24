/* $Id: ErrorCodeDocumentBuilder.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.webapp;

import static uk.ac.ed.ph.snuggletex.utilities.SnuggleUtilities.quoteTextForInput;

import uk.ac.ed.ph.snuggletex.ErrorCode;
import uk.ac.ed.ph.snuggletex.ErrorGroup;
import uk.ac.ed.ph.snuggletex.SnugglePackage;
import uk.ac.ed.ph.snuggletex.definitions.CorePackageDefinitions;
import uk.ac.ed.ph.snuggletex.internal.util.IOUtilities;
import uk.ac.ed.ph.snuggletex.upconversion.internal.UpConversionPackageDefinitions;
import uk.ac.ed.ph.snuggletex.utilities.MessageFormatter;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Builds a SnuggleTeX documentation page showing all of the error codes. This is called
 * as part of the webapp build process.
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
public final class ErrorCodeDocumentBuilder {
    
    private final File outputFile;
    
    public ErrorCodeDocumentBuilder(final File outputFile) {
        this.outputFile = outputFile;
    }
    
    public void run() throws IOException {
        IOUtilities.ensureFileCreated(outputFile);
        PrintWriter outputWriter = new PrintWriter(outputFile);
        outputWriter.println("\\pageId{errorCodes}");
        outputWriter.println("\n(In the tables below, \\{0\\} et\\ al\\ are placeholders for details "
                + "specific to each error instance that are substituted in when formatting error messages)");
        
        doSnugglePackage(outputWriter, CorePackageDefinitions.getPackage());
        doSnugglePackage(outputWriter, UpConversionPackageDefinitions.getPackage());
        
        outputWriter.close();
    }
    
    public void doSnugglePackage(PrintWriter outputWriter, SnugglePackage snugglePackage) {
        outputWriter.println("\n\\subsection*{Package: " + quoteTextForInput(snugglePackage.getName()) + "}\n");
        for (ErrorGroup errorGroup : snugglePackage.getErrorGroups()) {
            doErrorGroup(outputWriter, errorGroup);
        }
    }
    
    public void doErrorGroup(PrintWriter outputWriter, ErrorGroup errorGroup) {
        outputWriter.println("\n\\subsubsection*{Error Group: "
                + quoteTextForInput(MessageFormatter.formatErrorGroup(errorGroup))
                + "}");
        outputWriter.println("\\begin{tabular}{|c|l|}\n\\hline");
        for (ErrorCode errorCode : errorGroup.getPackage().getErrorCodes(errorGroup)) {
            String codeName = errorCode.getName();
            outputWriter.println("\\anchor{"
                    + codeName
                    + "}\\verb|"
                    + codeName
                    + "| & "
                    + quoteTextForInput(errorCode.getErrorGroup().getPackage().getErrorMessageBundle().getString(codeName))
                    + " \\\\ \\hline");
        }
        outputWriter.println("\\end{tabular}");
    }

    
    public static void main(String[] args) throws IOException {
        if (args.length!=1) {
            System.err.println("Please supply an output file");
            System.exit(1);
        }
        new ErrorCodeDocumentBuilder(new File(args[0])).run();
    }
}
