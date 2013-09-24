/* Random JavaScript for the SnuggleTeX webapp
 *
 * Prerequisites:
 *
 * jquery.js
 *
 * $Id: webapp.js 527 2010-01-08 11:50:20Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */

/* ============================================================ */

var popup = null;

jQuery(document).ready(function() {
    /* Initialise dialog box showing example material */
    popup = jQuery("#popup");
    popup.dialog({
      autoOpen: false,
      title: 'Example',
      width: 600,
      height: 400
    });
    /* Attach handlers to dialog popups links for example links */
    jQuery(".dialog").bind("click", function(event) {
        var latexInput = this.getAttribute('title');
        popup.load(this.href + " .exampleResult", null, function() {
            jQuery(".exampleResult").tabs();
            popup.dialog('option', 'title', 'Example: ' + latexInput);
            popup.dialog('open');
        });
        return false;
    });
});
