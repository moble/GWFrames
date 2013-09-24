<!--

$Id: snuggletex-utilities.xsl 525 2010-01-05 14:07:36Z davemckain $

Some common definitions for the up-conversion process.

NOTE: This uses an extension function to convert raw error codes
to full messages, which is convenient when inspecting the raw
XML output.

Copyright (c) 2010, The University of Edinburgh
All Rights Reserved

-->
<xsl:stylesheet version="2.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:xs="http://www.w3.org/2001/XMLSchema"
  xmlns:s="http://www.ph.ed.ac.uk/snuggletex"
  xmlns:util="ext://uk.ac.ed.ph.snuggletex.upconversion.UpConversionUtilities"
  xmlns:m="http://www.w3.org/1998/Math/MathML"
  xmlns="http://www.w3.org/1998/Math/MathML"
  exclude-result-prefixes="xs s m util">

  <!-- Creates an "up-conversion failure" element <s:fail/> with the given arguments -->
  <xsl:function name="s:make-error" as="element(s:fail)">
    <xsl:param name="code" as="xs:string"/>
    <xsl:param name="context" as="element()+"/>
    <xsl:param name="arguments" as="xs:string*"/>
    <s:fail code="{$code}" message="{util:getErrorMessage($code, $arguments)}">
      <xsl:for-each select="$arguments">
        <s:arg><xsl:value-of select="."/></s:arg>
      </xsl:for-each>
      <s:xpath>
        <xsl:value-of select="s:make-math-xpath($context[1])"/>
      </s:xpath>
      <s:context>
        <xsl:copy-of select="$context"/>
      </s:context>
    </s:fail>
  </xsl:function>

  <!--
  Makes an XPath expression from the given (MathML) element, rooted at
  the enclosing <math/> container.

  This should generate a valid XPath, provided the default namespace is
  set correctly.
  -->
  <xsl:function name="s:make-math-xpath" as="xs:string">
    <xsl:param name="element" as="element()"/>
    <xsl:sequence select="s:build-math-xpath($element, '')"/>
  </xsl:function>

  <xsl:function name="s:build-math-xpath" as="xs:string">
    <xsl:param name="element" as="element()"/>
    <xsl:param name="xpath-tail" as="xs:string"/>
    <xsl:variable name="path-element" as="xs:string" select="concat(
      local-name($element),
      '[', 1 + count($element/preceding-sibling::*[local-name()=local-name($element)]), ']',
      (if ($xpath-tail != '') then '/' else ''),
      $xpath-tail)"/>
    <xsl:variable name="math-parent" select="$element/parent::m:*" as="element()?"/>
    <xsl:sequence select="if (exists($math-parent))
        then s:build-math-xpath($math-parent, $path-element)
        else $path-element"/>
  </xsl:function>

</xsl:stylesheet>
