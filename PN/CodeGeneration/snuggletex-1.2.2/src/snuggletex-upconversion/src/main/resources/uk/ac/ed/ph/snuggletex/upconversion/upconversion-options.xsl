<!--

$Id: upconversion-options.xsl 525 2010-01-05 14:07:36Z davemckain $

This module contains templates for managing up-conversion options.

Copyright (c) 2010, The University of Edinburgh
All Rights Reserved

-->
<xsl:stylesheet version="2.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:xs="http://www.w3.org/2001/XMLSchema"
  xmlns:s="http://www.ph.ed.ac.uk/snuggletex"
  xmlns="http://www.w3.org/1998/Math/MathML"
  exclude-result-prefixes="xs s"
  xpath-default-namespace="http://www.w3.org/1998/Math/MathML">

  <xsl:function name="s:get-boolean-option" as="xs:boolean">
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)"/>
    <xsl:param name="name" as="xs:string"/>
    <xsl:sequence select="boolean($upconversion-options/s:option[@name=$name]/@value='true')"/>
  </xsl:function>

  <xsl:function name="s:get-upconversion-option" as="xs:string?">
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)"/>
    <xsl:param name="name" as="xs:string"/>
    <xsl:sequence select="$upconversion-options/s:option[@name=$name]/@value"/>
  </xsl:function>

  <xsl:function name="s:get-symbol-assumption" as="element(s:symbol)?">
    <xsl:param name="element" as="element()"/>
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)"/>
    <xsl:sequence select="$upconversion-options/s:symbol[deep-equal(*, $element)]"/>
  </xsl:function>

  <xsl:function name="s:is-assumed-symbol" as="xs:boolean">
    <xsl:param name="element" as="element()"/>
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)"/>
    <xsl:param name="assume" as="xs:string"/>
    <xsl:sequence select="exists($upconversion-options/s:symbol[@assume=$assume and deep-equal($element, *)])"/>
  </xsl:function>

  <xsl:function name="s:is-assumed-function" as="xs:boolean">
    <xsl:param name="element" as="element()"/>
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)"/>
    <xsl:sequence select="s:is-assumed-symbol($element, $upconversion-options, 'function')"/>
  </xsl:function>

</xsl:stylesheet>
