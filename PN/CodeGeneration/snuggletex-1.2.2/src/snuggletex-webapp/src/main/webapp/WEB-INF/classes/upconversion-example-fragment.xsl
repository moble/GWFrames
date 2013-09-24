<?xml version="1.0"?>
<!--

$Id: upconversion-example-fragment.xsl 527 2010-01-08 11:50:20Z davemckain $

Formats the mini pages containing up-conversion example results.
(These are normally loaded using AJAX with most of the content thrown away,
but we'll generated full pages here anyway as that can be useful.)

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

  <xsl:import href="demo-utilities.xsl"/>

  <xsl:param name="is-mathml-capable" as="xs:boolean" select="false()"/>
  <xsl:param name="is-internet-explorer" as="xs:boolean" select="false()"/>
  <xsl:param name="latex-input" as="xs:string" required="yes"/>
  <xsl:param name="is-bad-input" as="xs:boolean" required="yes"/>
  <xsl:param name="parsing-errors" as="element(s:error)*"/>
  <xsl:param name="parallel-mathml" as="xs:string?"/>
  <xsl:param name="pmathml-initial" as="xs:string?"/>
  <xsl:param name="pmathml-upconverted" as="xs:string?"/>
  <xsl:param name="cmathml" as="xs:string?"/>
  <xsl:param name="maxima-input" as="xs:string?"/>

  <xsl:template match="html">
    <xsl:copy>
      <xsl:copy-of select="@*"/>
      <head>
        <link rel="stylesheet" href="{$context-path}/includes/core.css" />
        <link rel="stylesheet" href="{$context-path}/includes/webapp.css" />
        <link rel="stylesheet" href="{$context-path}/includes/jquery-ui-1.7.2.custom.css" />
        <script type="text/javascript" src="{$context-path}/includes/jquery.js"></script>
        <script type="text/javascript" src="{$context-path}/includes/jquery-ui-1.7.2.custom.js"></script>
        <script type="text/javascript" src="{$context-path}/includes/webapp.js"></script>
        <script type="text/javascript">
          jQuery(document).ready(function() {
              jQuery(".exampleResult").tabs();
          });
        </script>
      </head>
      <body>
        <div class="exampleResult">
          <xsl:choose>
            <xsl:when test="$is-bad-input">
              <!-- Bad input -->
              <xsl:call-template name="handle-successful-input"/>
            </xsl:when>
            <xsl:when test="exists($parsing-errors)">
              <!-- SnuggleTeX Parsing Error(s) -->
              <xsl:call-template name="handle-failed-input"/>
            </xsl:when>
            <xsl:otherwise>
              <!-- Successful Parsing -->
              <xsl:call-template name="handle-successful-input"/>
            </xsl:otherwise>
          </xsl:choose>
        </div>
      </body>
    </xsl:copy>
  </xsl:template>

  <xsl:template name="handle-successful-input">
    <xsl:variable name="mathml" as="element(m:math)" select="//m:math[1]"/>
    <xsl:variable name="content-failures" as="element(s:fail)*" select="$mathml/m:semantics/m:annotation-xml[@encoding='MathML-Content-upconversion-failures']/*"/>
    <xsl:variable name="maxima-failures" as="element(s:fail)*" select="$mathml/m:semantics/m:annotation-xml[@encoding='Maxima-upconversion-failures']/*"/>
    <!-- (We generate JQuery UI tabs microformat) -->
    <ul>
      <li><a href="#ex-parallel">Parallel Markup</a></li>
      <li><a href="#ex-rpmathml">Raw PMathML</a></li>
      <li><a href="#ex-epmathml">Enhanced PMathML</a></li>
      <li><a href="#ex-cmathml">Content MathML</a></li>
      <li><a href="#ex-maxima">Maxima</a></li>
    </ul>
    <div id="ex-parallel">
      <pre>
        <xsl:value-of select="$parallel-mathml"/>
      </pre>
    </div>
    <div id="ex-rpmathml">
      <pre>
        <xsl:value-of select="$pmathml-initial"/>
      </pre>
    </div>
    <div id="ex-epmathml">
      <pre>
        <xsl:value-of select="$pmathml-upconverted"/>
      </pre>
    </div>
    <div id="ex-cmathml">
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
          <pre>
            <xsl:value-of select="$cmathml"/>
          </pre>
        </xsl:otherwise>
      </xsl:choose>
    </div>
    <div id="ex-maxima">
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
          <pre>
            <xsl:value-of select="$maxima-input"/>
          </pre>
        </xsl:otherwise>
      </xsl:choose>
    </div>
  </xsl:template>

  <!-- Show SnuggleTeX failure details. -->
  <xsl:template name="handle-failed-input">
    <h3>Result: LaTeX Errors in Input</h3>
    <p>
      The input LaTeX contained errors as follows:
    </p>
    <xsl:call-template name="format-parsing-errors">
      <xsl:with-param name="parsing-errors" select="$parsing-errors"/>
    </xsl:call-template>
  </xsl:template>

  <!-- Show bad input details -->
  <xsl:template name="handle-bad-input">
    <h3>Result: Bad Input</h3>
    <p>
      The input LaTeX was not successfully parsed as Math Mode input.
    </p>
  </xsl:template>

</xsl:stylesheet>
