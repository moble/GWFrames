<?xml version="1.0"?>
<!--

$Id: upconversion-demo.xsl 527 2010-01-08 11:50:20Z davemckain $

Overrides format-output.xsl to add in functionality for
demonstrating LaTeX -> Presentation MathML -> Content MathML -> Maxima
up-conversion process.

Copyright (c) 2010, The University of Edinburgh.
All Rights Reserved

-->
<xsl:stylesheet version="2.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:xs="http://www.w3.org/2001/XMLSchema"
  xmlns:m="http://www.w3.org/1998/Math/MathML"
  xmlns:s="http://www.ph.ed.ac.uk/snuggletex"
  xmlns="http://www.w3.org/1999/xhtml"
  xpath-default-namespace="http://www.w3.org/1999/xhtml"
  exclude-result-prefixes="xs m s">

  <!-- Import basic formatting stylesheet -->
  <xsl:import href="format-output.xsl"/>
  <xsl:import href="demo-utilities.xsl"/>

  <xsl:param name="latex-input" as="xs:string" required="yes"/>
  <xsl:param name="upconversion-options" as="xs:string" required="yes"/>
  <xsl:param name="is-bad-input" as="xs:boolean" required="yes"/>
  <xsl:param name="parsing-errors" as="element(s:error)*"/>
  <xsl:param name="parallel-mathml" as="xs:string?"/>
  <xsl:param name="pmathml-initial" as="xs:string?"/>
  <xsl:param name="pmathml-upconverted" as="xs:string?"/>
  <xsl:param name="cmathml" as="xs:string?"/>
  <xsl:param name="maxima-input" as="xs:string?"/>

  <!-- Override page ID -->
  <xsl:variable name="pageId" select="'upConversionDemo'" as="xs:string"/>

  <xsl:template match="body" mode="make-content">
    <!-- Do input form -->
    <h3>Input</h3>
    <p>
      This demo lets you try out the experimental
      <a href="{s:fix-href('docs://upconversion')}">Semantic Enrichment</a>
      process, converting simple LaTeX inputs to semantically richer Presentation MathML,
      Content MathML and Maxima input syntax.
    </p>
    <p>
      Simply enter a LaTeX math mode expression
      into the box below and hit <tt>Go!</tt> to see the resulting outputs.
    </p>
    <p>
      You can also click the <tt>Show Options</tt> button if you want to see
      or mess around with the various configurable options and assumptions
      that SnuggleTeX will use when processing your input.
    </p>
    <form method="post" class="input" action="{$context-path}/UpConversionDemo">
      <div class="inputBox">
        LaTeX Math Mode Input:
        <input id="input" name="input" type="text" value="{$latex-input}"/>
        <input type="submit" value="Go!" />
        <input id="showOptions" onclick="jQuery('#options').show(); jQuery('#showOptions').attr('disabled', true)" type="button" value="Show Options" />
        <input type="button" value="Clear" onclick="document.getElementById('input').value=''" />
      </div>
      <div id="options" style="display:none">
        <h3>Options and Assumptions</h3>
        <div class="optionsBox">
          <textarea id="upConversionOptions" name="upConversionOptions" rows="8"><xsl:value-of select="$upconversion-options"/></textarea>
        </div>
      </div>
    </form>
    <xsl:choose>
      <xsl:when test="$is-bad-input">
        <!-- Bad input -->
        <xsl:apply-templates select="." mode="handle-bad-input"/>
      </xsl:when>
      <xsl:when test="exists($parsing-errors)">
        <!-- SnuggleTeX Parsing Error(s) -->
        <xsl:apply-templates select="." mode="handle-failed-input"/>
      </xsl:when>
      <xsl:otherwise>
        <!-- Successful Parsing -->
        <xsl:apply-templates select="." mode="handle-successful-input"/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <!--
  This template shows the details when the input has been
  successfully passed through SnuggleTeX, though with the
  possibility that up-conversion has not been entirely successful.
  -->
  <xsl:template match="body" mode="handle-successful-input">
    <h3>MathML Rendering</h3>
    <xsl:call-template name="maybe-make-mathml-legacy-output-warning"/>
    <div class="result">
      <!-- (Strip out any SnuggleTeX elements from the result. These might appear in legacy
      outputs as they're intended for the up-converter but won't have gotten that far) -->
      <xsl:copy-of select="node()[not(self::s:*)]"/>
    </div>

    <h3>Raw Presentation MathML</h3>
    <p>
      This is the Presentation MathML normally produced by SnuggleTeX:
    </p>
    <pre class="result">
      <xsl:value-of select="$pmathml-initial"/>
    </pre>

    <h3>Enhanced Presentation MathML</h3>
    <p>
      This shows the result of attempting to infer semantics from the raw
      Presentation MathML, as described in
      <a href="documentation/pmathml-enhancement.html">Presentation MathML Enhancement</a>:
    </p>
    <pre class="result">
      <xsl:value-of select="$pmathml-upconverted"/>
    </pre>

    <h3>Content MathML</h3>
    <xsl:variable name="content-failures" as="element(s:fail)*" select="m:math/m:semantics/m:annotation-xml[@encoding='MathML-Content-upconversion-failures']/*"/>
    <xsl:choose>
      <xsl:when test="not(exists($cmathml)) and not(exists($content-failures))">
        <p>(You requested not to obtain Content MathML.)</p>
      </xsl:when>
      <xsl:otherwise>
        <p>
          This shows the result of an attempted
          <a href="documentation/content-mathml.html">conversion to Content MathML</a>:
        </p>
        <xsl:choose>
          <xsl:when test="exists($content-failures)">
            <p>
              The conversion from Presentation MathML to Content MathML was not successful
              for this input.
            </p>
            <xsl:call-template name="format-upconversion-failures">
              <xsl:with-param name="failures" select="$content-failures"/>
            </xsl:call-template>
          </xsl:when>
          <xsl:otherwise>
            <pre class="result">
              <xsl:value-of select="$cmathml"/>
            </pre>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:otherwise>
    </xsl:choose>

    <h3>Maxima Input Syntax</h3>
    <xsl:variable name="maxima-failures" as="element(s:fail)*" select="m:math/m:semantics/m:annotation-xml[@encoding='Maxima-upconversion-failures']/*"/>
    <xsl:choose>
      <xsl:when test="not(exists($maxima-input)) and not(exists($maxima-failures))">
        <p>(You requested not to obtain Maxima input.)</p>
      </xsl:when>
      <xsl:otherwise>
        <p>
          This shows the result of an attempted
          <a href="documentation/maxima-input.html">conversion to Maxima Input syntax</a>:
        </p>
        <xsl:choose>
          <xsl:when test="exists($content-failures)">
            <p>
              Conversion to Maxima Input is reliant on the conversion to Content MathML
              being successful, which was not the case here.
            </p>
          </xsl:when>
          <xsl:when test="exists($maxima-failures)">
            <p>
              The conversion from Content MathML to Maxima Input was not successful for
              this input.
            </p>
            <xsl:call-template name="format-upconversion-failures">
              <xsl:with-param name="failures" select="$maxima-failures"/>
            </xsl:call-template>
          </xsl:when>
          <xsl:otherwise>
            <pre class="result">
              <xsl:value-of select="$maxima-input"/>
            </pre>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:otherwise>
    </xsl:choose>

    <h3>MathML Parallel Markup</h3>
    <p>
      This shows the enhanced Presentation MathML with other forms encapsulated
      as annotations:
    </p>
    <pre class="result">
      <xsl:value-of select="$parallel-mathml"/>
    </pre>
  </xsl:template>

  <!-- Show SnuggleTeX failure details. -->
  <xsl:template match="body" mode="handle-failed-input">
    <h3>Result: LaTeX Errors in Input</h3>
    <p>
      The input you entered contained LaTeX errors as follows:
    </p>
    <xsl:call-template name="format-parsing-errors">
      <xsl:with-param name="parsing-errors" select="$parsing-errors"/>
    </xsl:call-template>
  </xsl:template>

  <!-- Show bad input details -->
  <xsl:template match="body" mode="handle-bad-input">
    <h3>Result: Bad Input</h3>
    <p>
      Sorry, your input was not successfully parsed as Math Mode input.
    </p>
  </xsl:template>
</xsl:stylesheet>
