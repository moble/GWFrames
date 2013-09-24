/*
 * This provides some basic code for managing the ASCIIMathML input widget
 * used in the ASCIIMathML input demo in SnuggleTeX.
 *
 * The general ideas may be useful in other scenarios, so feel free to use
 * and/or build on this as is necessary.
 *
 * NOTE:
 *
 * This code uses the lovely jQuery library, but avoids using the
 * $(...) function just in case your code also uses some other library
 * like prototype that defines its own $(...) function.
 * (In this case, you will still want to read:
 *
 * http://docs.jquery.com/Using_jQuery_with_Other_Libraries
 *
 * to make sure you do whatver is necessary to make sure that both
 * libraries co-exist correctly.)
 *
 * Requirements:
 *
 * ASCIIMathML.js
 * jquery.js
 *
 * Author: David McKain
 *
 * $Id:web.xml 158 2008-07-31 10:48:14Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh
 * All Rights Reserved
 */

/************************************************************/

/* (Reset certain defaults chosen by ASCIIMathML) */
var mathcolor = "";
var mathfontfamily = "";

/**
 * Simple hash that will keep track of the current value of each
 * ASCIIMathML input box, keyed on its ID. This is used by
 * updatePreviewIfChanged() to determine whether the MathML preview
 * should be updated or not.
 */
var inputTextByIdMap = {};

/**
 * Hacked version of AMdisplay() from ASCIIMathMLeditor.js that allows
 * us to specify which element to display the resulting MathML
 * in and where the raw input is going to come from.
 *
 * @param {String} mathModeInput ASCIIMathML input string
 * @param {String} previewElementId ID of the XHTML element that will contain the
 *   resulting MathML preview, replacing any existing child Nodes.
 */
function updatePreview(mathModeInput, previewElementId, previewSourceId) {
    /* Escape use of backquote symbol to prevent exiting math mode */
    mathModeInput = mathModeInput.replace(/`/, "\\`");
    var input = "` " + mathModeInput + " `";

    /* Do ASCIIMathML processing. (This is based on AMdisplay() from ASCIIMathMLeditor.js) */
    var outnode = document.getElementById(previewElementId);
    var newnode = AMcreateElementXHTML("div");
    newnode.setAttribute("id", previewElementId);
    outnode.parentNode.replaceChild(newnode, outnode);
    outnode = document.getElementById(previewElementId);
    var n = outnode.childNodes.length;
    for (var i=0; i<n; i++) {
        outnode.removeChild(outnode.firstChild);
    }
    outnode.appendChild(document.createComment(input + "``"));
    AMprocessNode(outnode, true);

    /* Now maybe update preview source box */
    if (previewSourceId!=null) {
        var source = extractMathML(previewElementId);
        jQuery("#" + previewSourceId).text(source);
    }
}

/**
 * Checks the content of the <input/> element having the given inputBoxId,
 * and calls {@link #updatePreview} if its contents have changed since the
 * last call to this.
 *
 * @param {String} inputBoxId ID of the ASCIIMath entry <input/> element
 * @param {String} previewElementId ID of the XHTML element that will contain the
 *   resulting MathML preview, replacing any existing child Nodes.
 */
function updatePreviewIfChanged(inputBoxId, previewElementId, previewSourceId) {
    var inputSelector = jQuery("#" + inputBoxId);
    var newValue = inputSelector.get(0).value;
    var oldValue = inputTextByIdMap[inputBoxId];
    if (oldValue==null || newValue!=oldValue) {
        updatePreview(newValue, previewElementId, previewSourceId);
    }
    inputTextByIdMap[inputBoxId] = newValue;
}

/**
 * Extracts the source MathML contained within the ASCIIMath preview element
 * having the given ID
 *
 * @param {String} previewElementId ID of the XHTML parent element
 *   containing the MathML to be extracted.
 */
function extractMathML(previewElementId) {
    var previewElement = document.getElementById(previewElementId);
    var mathNode = previewElement.getElementsByTagName("math")[0];
    return AMnode2string(mathNode, "").substring(1); // (First char is newline)
}

/* Fixed up version of the function of the same name in ASCIIMathMLeditor.js,
 * with the following changes:
 *
 * * Used '\n' for line breaks
 * * Attribute values are escape correctly
 */
function AMnode2string(inNode, indent) {
    var i, str = "";
    if (inNode.nodeType == 1) {
        var name = inNode.nodeName.toLowerCase(); // (IE fix)
        str = "\n" + indent + "<" + name;
        for (i=0; i < inNode.attributes.length; i++) {
            var attrValue = inNode.attributes[i].nodeValue;
            if (attrValue!="italic" &&
                    attrValue!="" &&  //stop junk attributes
                    attrValue!="inherit" && // (mostly IE)
                    attrValue!=undefined) {
                str += " " + inNode.attributes[i].nodeName
                    + "=\"" + AMescapeValue(inNode.attributes[i].nodeValue) + "\"";
            }
        }
        if (name == "math") str += " xmlns=\"http://www.w3.org/1998/Math/MathML\"";
        str += ">";
        for(i=0; i<inNode.childNodes.length; i++) {
            str += AMnode2string(inNode.childNodes[i], indent+"  ");
        }
        if (name != "mo" && name != "mi" && name != "mn") str += "\n"+indent;
        str += "</" + name + ">";
    }
    else if (inNode.nodeType == 3) {
        str += AMescapeValue(inNode.nodeValue);
    }
    return str;
}

function AMescapeValue(value) {
    var str = "";
    for (i=0; i<value.length; i++) {
        if (value.charCodeAt(i)<32 || value.charCodeAt(i)>126) str += "&#"+value.charCodeAt(i)+";";
        else if (value.charAt(i)=="<") str += "&lt;";
        else if (value.charAt(i)==">") str += "&gt;";
        else if (value.charAt(i)=="&") str += "&amp;";
        else str += value.charAt(i);
    }
    return str;
}

/**
 * Sets up an ASCIIMathML input using the elements provided, binding
 * the appropriate event handlers to make everything work correctly.
 *
 * @param {String} inputBoxId ID of the ASCIIMath entry <input/> element
 * @param {String} previewElementId ID of the XHTML element that will be
 *   used to hold the resulting MathML preview. Note that all of its child
 *   Nodes will be removed.
 * @param {String} mathmlResultControlId ID of the hidden <input/> field that will
 *   hold the resulting MathML on submission.
 * @param {String} previewSourceId optional ID of the XHTML element that will show
 *   the generated MathML source as it is built. This may be null to suppress this
 *   behaviour.
 */
function setupASCIIMathMLInputWidget(inputBoxId, previewElementId, mathmlResultControlId, previewSourceId) {
    /* Set up submit handler for the form */
    jQuery("#" + inputBoxId).closest("form").bind("submit", function(evt) {
        var mathmlResultControl = document.getElementById(mathmlResultControlId);
        var mathml = extractMathML(previewElementId, mathmlResultControlId);
        mathmlResultControl.value = mathml;
        return true;
    });
    var inputSelector = jQuery("#" + inputBoxId);
    var initialInput = inputSelector.get(0).value;

    /* Set up initial preview */
    updatePreview(initialInput, previewElementId, previewSourceId);

    /* Set up handler to update preview when required */
    inputSelector.bind("change keyup keydown", function() {
        updatePreviewIfChanged(inputBoxId, previewElementId, previewSourceId);
    });

    /* TODO: Do we want to set up a timer as well? If so, we probably want
     * one to be global to a page, rather than each interaction.
     */
}

/**
 * Registers a new ASCIIMathML input widget. This calls {@link #setupASCIIMathMLInputWidget}
 * once the document has finished loading to bind everything together correctly.
 *
 * @param {String} inputBoxId ID of the ASCIIMath entry <input/> element
 * @param {String} previewElementId ID of the XHTML element that will be
 *   used to hold the resulting MathML preview. Note that all of its child
 *   Nodes will be removed.
 * @param {String} mathmlResultControlId ID of the hidden <input/> field that will
 *   hold the resulting MathML on submission.
 * @param {String} previewSourceId optional ID of the XHTML element that will show
 *   the generated MathML source as it is built. This may be null to suppress this
 *   behaviour.
 */
function registerASCIIMathMLInputWidget(inputBoxId, previewElementId, mathmlResultControlId, previewSourceId) {
    jQuery(document).ready(function() {
        setupASCIIMathMLInputWidget(inputBoxId, previewElementId, mathmlResultControlId, previewSourceId);
    });
}
