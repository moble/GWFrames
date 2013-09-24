<?xml version="1.0"?>
<!--

$Id: full-latex-input-demo.xsl 575 2010-05-21 10:40:58Z davemckain $

Overrides format-output.xsl to add in functionality for
the full LaTeX input demo.

Copyright (c) 2010, The University of Edinburgh.
All Rights Reserved

-->
<xsl:stylesheet version="2.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:xs="http://www.w3.org/2001/XMLSchema"
  xmlns:m="http://www.w3.org/1998/Math/MathML"
  xmlns="http://www.w3.org/1999/xhtml"
  xpath-default-namespace="http://www.w3.org/1999/xhtml"
  exclude-result-prefixes="m xs">

  <xsl:import href="demo-utilities.xsl"/>
  <xsl:import href="format-output.xsl"/>

  <!-- Override page ID -->
  <xsl:variable name="pageId" select="'fullLaTeXInputDemo'" as="xs:string"/>

  <!-- LaTeX input - this will be put into a textarea -->
  <xsl:param name="latex-input" as="xs:string" required="yes"/>

  <xsl:template match="body" mode="make-content">
    <!-- Do input form -->
    <h3>Input</h3>
    <p>
      This demo lets you enter a chunk of (text mode) LaTeX for SnuggleTeX
      to convert into XHTML and MathML. You can include mathematics within your
      LaTeX in the usual way.
    </p>
    <p>
      Simply enter some LaTeX into the box below and hit <tt>Go!</tt> to see
      the results.
    </p>
    <p>
      (If you just want to convert a single mathematical
      formula then you will find the
      <a href="{$context-path}/MathInputDemo">Simple Math Input Demo</a>
      more useful in this case.)
    </p>
    <form method="post" class="input" action="{$context-path}/FullLaTeXInputDemo">
      <fieldset><!-- NB: This additional fieldset container - as well as the div - fixes the 'width:100%' bug for textareas in IE -->
        <div class="inputBox">
          <textarea id="inputBox" name="input" rows="20" cols="120">
            <xsl:value-of select="$latex-input"/>
          </textarea>
          <input type="submit" value="Go!" />
          <input type="button" value="Clear Form" onclick="document.getElementById('inputBox').value=''" />
        </div>
      </fieldset>
    </form>

    <!-- Output -->
    <h3>Output </h3>
    <xsl:if test="not($is-mathml-capable)">
      <xsl:call-template name="maybe-make-mathml-legacy-output-warning"/>
    </xsl:if>
    <div class="result">
      <xsl:copy-of select="node()"/>
    </div>
  </xsl:template>

</xsl:stylesheet>
