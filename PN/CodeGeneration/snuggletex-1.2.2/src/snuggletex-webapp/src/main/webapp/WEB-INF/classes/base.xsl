<?xml version="1.0"?>
<!--

$Id: base.xsl 525 2010-01-05 14:07:36Z davemckain $

Very trivial base stylesheet for the various types of pages/fragments
generated within the SnuggleTeX webapp.

Copyright (c) 2010, The University of Edinburgh.
All Rights Reserved

-->
<xsl:stylesheet version="2.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns="http://www.w3.org/1999/xhtml"
  xmlns:xs="http://www.w3.org/2001/XMLSchema"
  xmlns:m="http://www.w3.org/1998/Math/MathML"
  xmlns:s="http://www.ph.ed.ac.uk/snuggletex"
  xpath-default-namespace="http://www.w3.org/1999/xhtml"
  exclude-result-prefixes="s m xs">

  <!-- SnuggleTeX version string -->
  <xsl:param name="snuggletex-version" as="xs:string" required="yes"/>

  <!-- URL for the SnuggleTeX Maven site, used to generate documentation links -->
  <xsl:param name="maven-site-url" as="xs:string" required="yes"/>

  <!-- Context path for this webapp, used to generate internal links -->
  <xsl:param name="context-path" as="xs:string" required="yes"/>

</xsl:stylesheet>
