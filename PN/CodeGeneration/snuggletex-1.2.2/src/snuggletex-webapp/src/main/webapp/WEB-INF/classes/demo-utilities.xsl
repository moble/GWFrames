<?xml version="1.0"?>
<!--

$Id: demo-utilities.xsl 575 2010-05-21 10:40:58Z davemckain $

Some handy utility templates used by the various demos.

Copyright (c) 2010, The University of Edinburgh.
All Rights Reserved

-->
<xsl:stylesheet version="2.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:xs="http://www.w3.org/2001/XMLSchema"
  xmlns:m="http://www.w3.org/1998/Math/MathML"
  xmlns:s="http://www.ph.ed.ac.uk/snuggletex"
  xmlns:mu="ext://uk.ac.ed.ph.snuggletex.utilities.MathMLUtilities"
  xmlns="http://www.w3.org/1999/xhtml"
  xpath-default-namespace="http://www.w3.org/1999/xhtml"
  exclude-result-prefixes="xs m s mu">

  <xsl:import href="base.xsl"/>

  <xsl:param name="is-mathml-capable" as="xs:boolean" required="yes"/>
  <xsl:param name="is-internet-explorer" as="xs:boolean" required="yes"/>

  <xsl:template name="maybe-make-mathml-legacy-output-warning">
    <xsl:if test="not($is-mathml-capable)">
      <div class="warning">
        <p>
          Your browser does not support MathML so we have converted any MathML in this
          page to either XHTML + CSS or images, depending on how complex the MathML is.
        </p>
        <p>
          (Note: this statement will not necessarily be true if you are viewing
          a cached or proxied version of this page, e.g. within the Google cache or
          something similar. In this case, try viewing the page directly if you can.)
        </p>
        <xsl:if test="$is-internet-explorer">
          <p>
            You are using Internet Explorer, so you might want to consider installing
            the MathPlayer plugin, which will add in support for MathML.
          </p>
        </xsl:if>
      </div>
    </xsl:if>
  </xsl:template>

  <xsl:template name="format-parsing-errors">
    <xsl:param name="parsing-errors" as="element(s:error)+"/>
    <table class="failures">
      <thead>
        <tr>
          <th>Error Code</th>
          <th>Message</th>
        </tr>
      </thead>
      <tbody>
        <xsl:for-each select="$parsing-errors">
          <tr>
            <td>
              <xsl:call-template name="make-error-code-link">
                <xsl:with-param name="error-code" select="@code"/>
              </xsl:call-template>
            </td>
            <td>
              <pre>
                <xsl:value-of select="."/>
              </pre>
            </td>
          </tr>
        </xsl:for-each>
      </tbody>
    </table>
  </xsl:template>

  <xsl:template name="format-upconversion-failures">
    <xsl:param name="failures" as="element(s:fail)+"/>
    <table class="failures">
      <thead>
        <tr>
          <th>Failure Code</th>
          <th>Message</th>
          <th>XPath</th>
          <th>Context</th>
        </tr>
      </thead>
      <tbody>
        <xsl:for-each select="$failures">
          <tr>
            <td>
              <xsl:call-template name="make-error-code-link">
                <xsl:with-param name="error-code" select="@code"/>
              </xsl:call-template>
            </td>
            <td><xsl:value-of select="@message"/></td>
            <td><xsl:value-of select="s:xpath"/></td>
            <td>
              <pre>
                <!--
                We'll strip off the enclosing <s:fail/>, which also conveniently
                removes namespace declarations.
                -->
                <xsl:value-of select="replace(
                  replace(
                    replace(mu:serializeElement(s:context), '^&lt;s:context.+?&gt; ?[\n\r]+', '', 's'),
                    '[\n\r]+ ?&lt;/s:context&gt;$', '', 's'),
                    '^   ', '', 'm') (: de-indent 1 level :)
                  "/>
              </pre>
            </td>
          </tr>
        </xsl:for-each>
      </tbody>
    </table>
  </xsl:template>

  <xsl:template name="make-error-code-link">
    <xsl:param name="error-code" as="xs:string"/>
    <a href="{$context-path}/documentation/error-codes.html#{$error-code}">
      <xsl:value-of select="$error-code"/>
    </a>
  </xsl:template>

</xsl:stylesheet>


