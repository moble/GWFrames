<!--

$Id: pmathml-utilities.xsl 525 2010-01-05 14:07:36Z davemckain $

This stylesheet contains some utility templates for handling Presentation
MathML.

I've factored this out as we use some of these in both the Presentation
MathML enrichment stylesheet and the PMathML -> CMathML stylesheet.

Copyright (c) 2010, The University of Edinburgh
All Rights Reserved

-->
<xsl:stylesheet version="2.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:xs="http://www.w3.org/2001/XMLSchema"
  xmlns:s="http://www.ph.ed.ac.uk/snuggletex"
  xmlns:m="http://www.w3.org/1998/Math/MathML"
  xmlns="http://www.w3.org/1998/Math/MathML"
  exclude-result-prefixes="xs m s"
  xpath-default-namespace="http://www.w3.org/1998/Math/MathML">

  <!--
  Helper to wrap the given elements inside an <mrow/> if that is
  deemed necessary, which is quite common in certain MathML
  constructs.
  -->
  <xsl:template name="s:maybe-wrap-in-mrow" as="element()">
    <xsl:param name="elements" as="element()*" required="yes"/>
    <xsl:choose>
      <xsl:when test="count($elements)=1">
        <xsl:copy-of select="$elements"/>
      </xsl:when>
      <xsl:otherwise>
        <mrow>
          <xsl:copy-of select="$elements"/>
        </mrow>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <!--
  Handy utility function for checking whether an element is an <mo/>
  -->
  <xsl:function name="s:is-operator" as="xs:boolean">
    <xsl:param name="element" as="element()"/>
    <xsl:sequence select="boolean($element[self::mo])"/>
  </xsl:function>

  <xsl:function name="s:is-power" as="xs:boolean">
    <xsl:param name="element" as="element()"/>
    <xsl:sequence select="boolean($element[self::msup or self::msubsup])"/>
  </xsl:function>

  <xsl:function name="s:get-power-base" as="element()?">
    <xsl:param name="element" as="element()"/>
    <xsl:choose>
      <xsl:when test="$element[self::msup]">
        <xsl:sequence select="$element/*[1]"/>
      </xsl:when>
      <xsl:when test="$element[self::msubsup]">
        <msub>
          <xsl:copy-of select="$element/*[position()!=3]"/>
        </msub>
      </xsl:when>
    </xsl:choose>
  </xsl:function>

  <xsl:function name="s:get-power-exponent" as="element()?">
    <xsl:param name="element" as="element()"/>
    <xsl:choose>
      <xsl:when test="$element[self::msup]">
        <xsl:sequence select="$element/*[2]"/>
      </xsl:when>
      <xsl:when test="$element[self::msubsup]">
        <xsl:sequence select="$element/*[3]"/>
      </xsl:when>
    </xsl:choose>
  </xsl:function>

</xsl:stylesheet>
