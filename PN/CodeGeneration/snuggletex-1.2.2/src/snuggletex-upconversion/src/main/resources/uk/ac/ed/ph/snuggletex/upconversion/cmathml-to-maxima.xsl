<!--

$Id: cmathml-to-maxima.xsl 576 2010-05-21 11:23:19Z davemckain $

This stylesheet is intended to be pipelined after the P->C stylesheet
and converts the subset of Content MathML produced by that stylesheet
into Maxima input.

IMPORTANT NOTE: This stylesheet is NOT intended to be applied to more general
Content MathML elements as it assumes certain post-conditions of the earlier
conversion to Content MathML, making life very easy here.

Copyright (c) 2010, The University of Edinburgh
All Rights Reserved

TODO: Handle the lack of support for log to base 10 (or indeed other bases)

-->
<xsl:stylesheet version="2.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:xs="http://www.w3.org/2001/XMLSchema"
  xmlns:s="http://www.ph.ed.ac.uk/snuggletex"
  xmlns:local="http://www.ph.ed.ac.uk/snuggletex/cmathml-to-maxima"
  xmlns="http://www.w3.org/1998/Math/MathML"
  exclude-result-prefixes="xs s local"
  xpath-default-namespace="http://www.w3.org/1998/Math/MathML">

  <xsl:import href="snuggletex-utilities.xsl"/>

  <!-- ************************************************************ -->

  <!-- Entry Point -->
  <xsl:template name="s:cmathml-to-maxima" as="node()*">
    <xsl:param name="elements" as="item()*"/>
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)"/>
    <xsl:apply-templates select="$elements" mode="cmathml-to-maxima">
      <xsl:with-param name="upconversion-options" select="$upconversion-options" tunnel="yes"/>
    </xsl:apply-templates>
  </xsl:template>

  <!-- ************************************************************ -->

  <!-- Supported non-alphanumeric identifiers, mapping Unicode character to Maxima input -->
  <xsl:variable name="local:identifier-dictionary" as="element(ci)+">
    <ci maxima-input="%ell">&#x2113;</ci>
    <ci maxima-input="%alpha">&#x3b1;</ci>
    <ci maxima-input="%beta">&#x3b2;</ci>
    <ci maxima-input="%gamma">&#x3b3;</ci>
    <ci maxima-input="%delta">&#x3b4;</ci>
    <ci maxima-input="%epsilon">&#x3f5;</ci>
    <ci maxima-input="%zeta">&#x3b6;</ci>
    <ci maxima-input="%eta">&#x3b7;</ci>
    <ci maxima-input="%theta">&#x3b8;</ci>
    <ci maxima-input="%iota">&#x3b9;</ci>
    <ci maxima-input="%kappa">&#x3ba;</ci>
    <ci maxima-input="%lambda">&#x3bb;</ci>
    <ci maxima-input="%mu">&#x3bc;</ci>
    <ci maxima-input="%nu">&#x3bd;</ci>
    <ci maxima-input="%xi">&#x3be;</ci>
    <ci maxima-input="%pi">&#x3c0;</ci>
    <ci maxima-input="%rho">&#x3c1;</ci>
    <ci maxima-input="%sigma">&#x3c3;</ci>
    <ci maxima-input="%tau">&#x3c4;</ci>
    <ci maxima-input="%upsilon">&#x3c5;</ci>
    <ci maxima-input="%phi">&#x3c6;</ci>
    <ci maxima-input="%chi">&#x3c7;</ci>
    <ci maxima-input="%psi">&#x3c8;</ci>
    <ci maxima-input="%omega">&#x3c9;</ci>
    <ci maxima-input="%Gamma">&#x393;</ci>
    <ci maxima-input="%Delta">&#x394;</ci>
    <ci maxima-input="%Theta">&#x398;</ci>
    <ci maxima-input="%Lambda">&#x39b;</ci>
    <ci maxima-input="%Xi">&#x39e;</ci>
    <ci maxima-input="%Pi">&#x3a0;</ci>
    <ci maxima-input="%Sigma">&#x3a3;</ci>
    <ci maxima-input="%Upsilon">&#x3a5;</ci>
    <ci maxima-input="%Phi">&#x3a6;</ci>
    <ci maxima-input="%Psi">&#x3a8;</ci>
    <ci maxima-input="%Omega">&#x3a9;</ci>
  </xsl:variable>

  <!-- Supported functions, named after CMathML elements -->
  <xsl:variable name="local:supported-functions" as="element()+">
    <!-- The resulting Maxima function name is encoded within an input Content MathML element -->
    <sin maxima-function="sin"/>
    <cos maxima-function="cos"/>
    <tan maxima-function="tan"/>
    <sec maxima-function="sec"/>
    <csc maxima-function="csc"/>
    <cot maxima-function="cot"/>
    <arcsin maxima-function="asin"/>
    <arccos maxima-function="acos"/>
    <arctan maxima-function="atan"/>
    <arcsec maxima-function="asec"/>
    <arccsc maxima-function="acsc"/>
    <arccot maxima-function="acot"/>
    <sinh maxima-function="sinh"/>
    <cosh maxima-function="cosh"/>
    <tanh maxima-function="tanh"/>
    <sech maxima-function="sech"/>
    <csch maxima-function="csch"/>
    <coth maxima-function="coth"/>
    <arcsinh maxima-function="asinh"/>
    <arccosh maxima-function="acosh"/>
    <arctanh maxima-function="atanh"/>
    <arcsech maxima-function="asech"/>
    <arccsch maxima-function="acsch"/>
    <arccoth maxima-function="acoth"/>
    <exp maxima-function="exp"/>
    <ln maxima-function="log"/>
    <determinant maxima-function="determinant"/>
    <real maxima-function="realpart"/>
    <imaginary maxima-function="imagpart"/>
    <gcd maxima-function="gcd" require-nary="true"/>
    <lcm maxima-function="lcm" require-nary="true"/>
    <max maxima-function="max" require-nary="true"/>
    <min maxima-function="min" require-nary="true"/>
  </xsl:variable>

  <xsl:function name="local:is-supported-function" as="xs:boolean">
    <xsl:param name="element" as="element()"/>
    <xsl:sequence select="boolean($local:supported-functions[local-name()=$element/local-name()])"/>
  </xsl:function>

  <xsl:function name="local:get-supported-function" as="element()?">
    <xsl:param name="element" as="element()"/>
    <xsl:sequence select="$local:supported-functions[local-name()=$element/local-name()]"/>
  </xsl:function>

  <!-- Supported prefix/infix/postfix operators -->
  <xsl:variable name="local:operators" as="element()+">
    <eq maxima-unapplied-operator="=" maxima-nary-infix-operator=" = "/>
    <neq maxima-unapplied-operator="#" maxima-nary-infix-operator=" # "/>
    <lt maxima-unapplied-operator="&lt;" maxima-nary-infix-operator=" &lt; "/>
    <gt maxima-unapplied-operator="&gt;" maxima-nary-infix-operator=" &gt; "/>
    <leq maxima-unapplied-operator="&lt;=" maxima-nary-infix-operator=" &lt;= "/>
    <geq maxima-unapplied-operator="&gt;=" maxima-nary-infix-operator=" &gt;= "/>
    <plus maxima-unapplied-operator="+" maxima-nary-infix-operator=" + " maxima-unary-prefix-operator="+"/>
    <minus maxima-unapplied-operator="-" maxima-nary-infix-operator=" - " maxima-unary-prefix-operator="-"/>
    <times maxima-unapplied-operator="*" maxima-nary-infix-operator=" * "/>
    <divide maxima-unapplied-operator="/" maxima-nary-infix-operator=" / "/>
    <power maxima-unapplied-operator="^" maxima-nary-infix-operator="^"/>
    <factorial maxima-unapplied-operator="!" maxima-unary-postfix-operator="!"/>
    <not maxima-unapplied-operator="not" maxima-unary-function="not"/>
    <and maxima-unapplied-operator="and" maxima-nary-infix-operator=" and "/>
    <or maxima-unapplied-operator="or" maxima-nary-infix-operator=" or "/>
    <union maxima-unapplied-operator="union" maxima-nary-function="union"/>
    <intersect maxima-unapplied-operator="intersection" maxima-nary-function="intersection"/>
    <setdiff maxima-unapplied-operator="setdifference" maxima-nary-function="setdifference"/>
  </xsl:variable>

  <xsl:function name="local:is-operator" as="xs:boolean">
    <xsl:param name="element" as="element()"/>
    <xsl:sequence select="boolean($local:operators[local-name()=$element/local-name()])"/>
  </xsl:function>

  <xsl:function name="local:get-operator" as="element()?">
    <xsl:param name="element" as="element()"/>
    <xsl:sequence select="$local:operators[local-name()=$element/local-name()]"/>
  </xsl:function>

  <!-- ************************************************************ -->
  <!-- "Functional" helpers  -->
  <!-- (Recall that these may return either xs:string or an <s:fail/> element) -->

  <xsl:function name="local:to-maxima" as="node()*">
    <xsl:param name="element" as="element()"/>
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)"/>
    <xsl:apply-templates select="$element" mode="cmathml-to-maxima">
      <xsl:with-param name="upconversion-options" select="$upconversion-options" tunnel="yes"/>
    </xsl:apply-templates>
  </xsl:function>

  <xsl:function name="local:to-maxima-map" as="node()*">
    <xsl:param name="elements" as="element()*"/>
    <xsl:param name="joiner" as="xs:string"/>
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)"/>
    <xsl:for-each select="$elements">
      <xsl:copy-of select="local:to-maxima(., $upconversion-options)"/>
      <xsl:if test="position() != last()">
        <xsl:value-of select="$joiner"/>
      </xsl:if>
    </xsl:for-each>
  </xsl:function>

  <!-- ************************************************************ -->

  <!--
  Unapplied infix operator

  Example:

  <plus/>

  Output:

  operator("+")
  -->
  <xsl:template match="*[local:is-operator(.)]" mode="cmathml-to-maxima" as="text()">
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)" tunnel="yes"/>
    <xsl:variable name="operator" select="local:get-operator(.)" as="element()"/>
    <xsl:value-of select="local:unapply-operator($operator, false(), $upconversion-options)"/>
  </xsl:template>

  <xsl:function name="local:unapply-operator" as="xs:string">
    <xsl:param name="operator" as="element()"/>
    <xsl:param name="negate" as="xs:boolean"/>
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)"/>
    <xsl:variable name="maxima-raw-name" as="xs:string?" select="$operator/@maxima-unapplied-operator"/>
    <xsl:choose>
      <xsl:when test="$maxima-raw-name">
        <xsl:copy-of select="local:make-unapplied-operator(if ($negate)
          then concat('not', $maxima-raw-name) else $maxima-raw-name,
          $upconversion-options)"/>
      </xsl:when>
      <xsl:otherwise>
        <xsl:message terminate="yes">
          Operator <xsl:value-of select="$operator/local-name()"/> cannot be
          used in an unapplied context.
        </xsl:message>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:function>

  <xsl:function name="local:make-unapplied-operator" as="xs:string">
    <xsl:param name="operator" as="xs:string"/>
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)"/>
    <xsl:variable name="maxima-operator-function" select="s:get-upconversion-option($upconversion-options, 'maximaOperatorFunction')" as="xs:string"/>
    <xsl:value-of select="concat($maxima-operator-function, '(&quot;', $operator, '&quot;)')"/>
  </xsl:function>

  <!--
  Unapplied and negated infix operator

  Example:

  <apply>
    <not/>
    <plus/>
  </apply>

  For better or worse, this would output:

  operator("not+")
  -->
  <xsl:template match="apply[count(*)=2 and *[1][self::not] and local:is-operator(*[2])]" mode="cmathml-to-maxima" as="text()">
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)" tunnel="yes"/>
    <xsl:variable name="operator" as="element()" select="local:get-operator(*[2])"/>
    <xsl:value-of select="local:unapply-operator($operator, true(), $upconversion-options)"/>
  </xsl:template>

  <!--
  Half-arsed Applied infix operator

  Example:

  <apply>
    <plus/>
  </apply>

  This one doesn't really make any sense; we'll pretend it's an unapplied
  operator for the time being and output:

  operator("+")
  -->
  <xsl:template match="apply[count(*)=1 and local:is-operator(*[1])]" mode="cmathml-to-maxima" as="text()">
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)" tunnel="yes"/>
    <xsl:variable name="operator" as="element()" select="local:get-operator(*[1])"/>
    <xsl:value-of select="local:unapply-operator($operator, false(), $upconversion-options)"/>
  </xsl:template>

  <!--
  Applied infix operator.

  Example:

  <apply>
    <plus/>
    <ci>x</ci>
    <cn>5</cn>
  </apply>
  -->
  <xsl:template match="apply[count(*)&gt;1 and local:is-operator(*[1]) and not(local:is-operator(*[2]))]" mode="cmathml-to-maxima" as="node()+">
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)" tunnel="yes"/>
    <xsl:variable name="operator" as="element()" select="local:get-operator(*[1])"/>
    <xsl:variable name="operands" as="element()+" select="*[position()!=1]"/>
    <xsl:choose>
      <xsl:when test="count($operands)=1">
        <!-- Unary case -->
        <xsl:choose>
          <xsl:when test="$operator/@maxima-unary-prefix-operator">
            <!-- Prefix operator, e.g. '-' in unary context -->
            <xsl:text>(</xsl:text>
            <xsl:value-of select="$operator/@maxima-unary-prefix-operator"/>
            <xsl:copy-of select="local:to-maxima($operands[1], $upconversion-options)"/>
            <xsl:text>)</xsl:text>
          </xsl:when>
          <xsl:when test="$operator/@maxima-unary-function">
            <!-- Function operator e.g. not() -->
            <xsl:value-of select="$operator/@maxima-unary-function"/>
            <xsl:text>(</xsl:text>
            <xsl:copy-of select="local:to-maxima($operands[1], $upconversion-options)"/>
            <xsl:text>)</xsl:text>
          </xsl:when>
          <xsl:when test="$operator/@maxima-unary-postfix-operator">
            <!-- Postfix operator -->
            <xsl:text>(</xsl:text>
            <xsl:copy-of select="local:to-maxima($operands[1], $upconversion-options)"/>
            <xsl:value-of select="$operator/@maxima-unary-postfix-operator"/>
            <xsl:text>)</xsl:text>
          </xsl:when>
          <xsl:otherwise>
            <xsl:copy-of select="s:make-error('UMFOP0', ., ($operator))"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:otherwise>
        <!-- nary case (NOTE: Earlier stylesheet will ensure binary when required) -->
        <xsl:choose>
          <xsl:when test="$operator/@maxima-nary-infix-operator">
            <!-- Operator like '+' -->
            <xsl:text>(</xsl:text>
            <xsl:copy-of select="local:to-maxima-map($operands, $operator/@maxima-nary-infix-operator, $upconversion-options)"/>
            <xsl:text>)</xsl:text>
          </xsl:when>
          <xsl:when test="$operator/@maxima-nary-function">
            <!-- Function operator, e.g. union() -->
            <xsl:value-of select="$operator/@maxima-nary-function"/>
            <xsl:text>(</xsl:text>
            <xsl:for-each select="$operands">
              <xsl:copy-of select="local:to-maxima(., $upconversion-options)"/>
              <xsl:if test="position()!=last()">
                <xsl:text>, </xsl:text>
              </xsl:if>
            </xsl:for-each>
            <xsl:text>)</xsl:text>
          </xsl:when>
          <xsl:otherwise>
            <xsl:message terminate="yes">
              Operator <xsl:value-of select="$operator/local-name()"/> cannot
              be used in an n-ary context.
            </xsl:message>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <!--
  Square Root

  Example:

  <apply>
    <root/>
    <ci>x</ci>
  </apply>

  -->
  <xsl:template match="apply[*[1][self::root] and not(degree) and count(*)=2]" mode="cmathml-to-maxima" as="node()+">
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)" tunnel="yes"/>
    <xsl:variable name="operand" as="element()" select="*[2]"/>
    <xsl:text>sqrt(</xsl:text>
    <xsl:copy-of select="local:to-maxima($operand, $upconversion-options)"/>
    <xsl:text>)</xsl:text>
  </xsl:template>

  <!--
  nth Root

  Example:

  <apply>
    <root/>
    <degree>
      <ci>n</ci>
    </degree>
    <ci>x</ci>
  </apply>
  -->
  <xsl:template match="apply[*[1][self::root] and degree and count(*)=3]" mode="cmathml-to-maxima" as="node()+">
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)" tunnel="yes"/>
    <xsl:variable name="operand" as="element()" select="*[3]"/>
    <xsl:text>((</xsl:text>
    <xsl:copy-of select="local:to-maxima($operand, $upconversion-options)"/>
    <xsl:text>)^(1/</xsl:text>
    <xsl:copy-of select="local:to-maxima(degree/*, $upconversion-options)"/>
    <xsl:text>))</xsl:text>
  </xsl:template>

  <!--
  Unapplied function

  Example:

  <sin/>
  -->
  <xsl:template match="*[local:is-supported-function(.)]" mode="cmathml-to-maxima" as="text()">
    <xsl:value-of select="local:get-supported-function(.)/@maxima-function"/>
  </xsl:template>

  <!--
  Applied Function, allowing n-ary application when supported.

  Example:

  <apply>
    <sin/>
    <ci>x</ci>
  </apply>

  -->
  <xsl:template match="apply[count(*) &gt;= 2 and *[1][local:is-supported-function(.)]]" mode="cmathml-to-maxima" as="node()+">
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)" tunnel="yes"/>
    <xsl:variable name="function" as="element()" select="local:get-supported-function(*[1])"/>
    <xsl:variable name="arguments" as="element()+" select="*[position()!=1]"/>
    <xsl:choose>
      <xsl:when test="count($arguments)=1 and $function/@require-nary='true'">
        <!-- Fail: function cannot be used in a unary context -->
        <xsl:copy-of select="s:make-error('UMFFX0', ., ($function/@maxima-function))"/>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="$function/@maxima-function"/>
        <xsl:text>(</xsl:text>
        <xsl:copy-of select="local:to-maxima-map($arguments, ', ', $upconversion-options)"/>
        <xsl:text>)</xsl:text>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <!--
  Half-arsed Applied Function

  Example:

  <apply>
    <sin/>
  </apply>

  This one doesn't really make any sense; we'll pretend it's an unapplied
  function for the time being

  -->
  <xsl:template match="apply[count(*)=1 and *[1][local:is-supported-function(.)]]" mode="cmathml-to-maxima" as="text()">
    <xsl:value-of select="local:get-supported-function(*[1])/@maxima-function"/>
  </xsl:template>

  <!--
  Power of a function. For example:

  <apply>
    <apply>
      <power/>
      <sin/>
      <cn>2</cn>
    </apply>
    <ci>x</ci>
  </apply>
  -->
  <xsl:template match="apply[count(*) &gt;= 2 and *[1][self::apply and *[1][self::power] and local:is-supported-function(*[2]) and *[3][self::cn]]]" mode="cmathml-to-maxima" as="node()+">
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)" tunnel="yes"/>
    <xsl:variable name="function" as="element()" select="local:get-supported-function(*[1]/*[2])"/>
    <xsl:variable name="power" as="element()" select="*[1]/*[3]"/>
    <xsl:variable name="arguments" as="element()+" select="*[position()!=1]"/>
    <xsl:choose>
      <xsl:when test="count($arguments)=1 and $function/@require-nary='true'">
        <!-- Fail: function must be used in n-ary context -->
        <xsl:copy-of select="s:make-error('UMFFX0', ., ($function/@maxima-function))"/>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="$function/@maxima-function"/>
        <xsl:text>(</xsl:text>
        <xsl:copy-of select="local:to-maxima-map($arguments, ', ', $upconversion-options)"/>
        <xsl:text>)^</xsl:text>
        <xsl:copy-of select="local:to-maxima($power, $upconversion-options)"/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template match="apply[count(*)=1 and *[1][self::apply and *[1][self::power] and local:is-supported-function(*[2]) and *[3][self::cn]]]" mode="cmathml-to-maxima">
    <xsl:message terminate="yes">
      Power of function <xsl:value-of select="local:get-supported-function(*[1]/*[2])/@maxima-function"/>
      was expected to take at least one argument
    </xsl:message>
  </xsl:template>

  <!--
  Custom function application.

  <apply>
    <ci type="function">f</ci>
    <ci>x</ci>
  </apply>
  -->
  <xsl:template match="apply[count(*) &gt;= 2 and *[1][self::ci and @type='function']]" mode="cmathml-to-maxima" as="node()+">
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)" tunnel="yes"/>
    <xsl:variable name="function" as="element(ci)" select="*[1]"/>
    <xsl:variable name="arguments" as="element()+" select="*[position()!=1]"/>
    <xsl:apply-templates select="$function" mode="cmathml-to-maxima"/>
    <xsl:text>(</xsl:text>
    <xsl:copy-of select="local:to-maxima-map($arguments, ', ', $upconversion-options)"/>
    <xsl:text>)</xsl:text>
  </xsl:template>

  <!--
  Inverse of custom function application.

  <apply>
    <apply>
      <inverse/>
      <ci type="function">f</ci>
    </apply>
    <ci>x</ci>
  </apply>
  -->
  <xsl:template match="apply[count(*) &gt;= 2 and *[1][self::apply and *[1][self::inverse] and *[2][self::ci and @type='function']]]" mode="cmathml-to-maxima" as="node()+">
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)" tunnel="yes"/>
    <xsl:variable name="function" as="element(ci)" select="*[1]/*[2]"/>
    <xsl:variable name="arguments" as="element()+" select="*[position()!=1]"/>
    <xsl:value-of select="s:get-upconversion-option($upconversion-options, 'maximaInverseFunction')"/>
    <xsl:text>[</xsl:text>
    <xsl:apply-templates select="$function" mode="cmathml-to-maxima"/>
    <xsl:text>](</xsl:text>
    <xsl:copy-of select="local:to-maxima-map($arguments, ', ', $upconversion-options)"/>
    <xsl:text>)</xsl:text>
  </xsl:template>

  <!--
  Power of custom function application.

  <apply>
    <apply>
      <power/>
      <ci type="function">f</ci>
      <cn>5</cn>
    </apply>
    <ci>x</ci>
  </apply>
  -->
  <xsl:template match="apply[count(*) &gt;= 2 and *[1][self::apply and *[1][self::power] and *[2][self::ci and @type='function']]]" mode="cmathml-to-maxima" as="node()+">
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)" tunnel="yes"/>
    <xsl:variable name="function" as="element(ci)" select="*[1]/*[2]"/>
    <xsl:variable name="power" as="element()" select="*[1]/*[3]"/>
    <xsl:variable name="arguments" as="element()+" select="*[position()!=1]"/>
    <xsl:apply-templates select="$function" mode="cmathml-to-maxima"/>
    <xsl:text>(</xsl:text>
    <xsl:copy-of select="local:to-maxima-map($arguments, ', ', $upconversion-options)"/>
    <xsl:text>)^</xsl:text>
    <xsl:copy-of select="local:to-maxima($power, $upconversion-options)"/>
  </xsl:template>

  <!--
  Unapplied Inverse of custom function

  <apply>
    <inverse/>
    <ci type="function">f</ci>
  </apply>
  -->
  <xsl:template match="apply[count(*)=2 and *[1][self::inverse] and *[2][self::ci and @type='function']]" mode="cmathml-to-maxima" as="node()+">
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)" tunnel="yes"/>
    <xsl:variable name="function" as="element(ci)" select="*[2]"/>
    <xsl:value-of select="s:get-upconversion-option($upconversion-options, 'maximaInverseFunction')"/>
    <xsl:text>[</xsl:text>
    <xsl:apply-templates select="$function" mode="cmathml-to-maxima"/>
    <xsl:text>]</xsl:text>
  </xsl:template>

  <!--
  Unapplied Power of custom function

  <apply>
    <power/>
    <ci type="function">f</ci>
    <cn>5</cn>
  </apply>

  NOTE: The end result here isn't very good as the best we can do is construct an anonymous function but
  we don't know its arity.
  -->
  <xsl:template match="apply[count(*)=3 and *[1][self::power] and *[2][self::ci and @type='function']]" mode="cmathml-to-maxima" as="node()+" priority="5">
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)" tunnel="yes"/>
    <xsl:variable name="function" as="element(ci)" select="*[2]"/>
    <xsl:variable name="power" as="element()" select="*[3]"/>
    <xsl:variable name="arguments" as="element()+" select="*[position()!=1]"/>
    <xsl:text>lambda([x], </xsl:text>
    <xsl:apply-templates select="$function" mode="cmathml-to-maxima"/>
    <xsl:text>(x)^</xsl:text>
    <xsl:copy-of select="local:to-maxima($power, $upconversion-options)"/>
    <xsl:text>)</xsl:text>
  </xsl:template>

  <!-- ************************************************************ -->

  <!-- Maxima doesn't actually support intervals! -->
  <xsl:template match="interval" mode="cmathml-to-maxima" as="element(s:fail)">
    <xsl:copy-of select="s:make-error('UMFG02', ., ())"/>
  </xsl:template>

  <xsl:template match="set" mode="cmathml-to-maxima" as="node()+">
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)" tunnel="yes"/>
    <xsl:text>{</xsl:text>
    <xsl:copy-of select="local:to-maxima-map(*, ', ', $upconversion-options)"/>
    <xsl:text>}</xsl:text>
  </xsl:template>

  <xsl:template match="list|vector" mode="cmathml-to-maxima" as="node()+">
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)" tunnel="yes"/>
    <xsl:text>[</xsl:text>
    <xsl:copy-of select="local:to-maxima-map(*, ', ', $upconversion-options)"/>
    <xsl:text>]</xsl:text>
  </xsl:template>

  <!-- ************************************************************ -->

  <!--
  Handles specially marked-up units created via the special \units{.} macro.

  This will have produced CMathML of the following form:

  <semantics definitionURL="http://www.ph.ed.ac.uk/snuggletex/units">
    <csymbol>kg</csymbol>
  </semantics>
  -->
  <xsl:template match="semantics[@definitionURL='http://www.ph.ed.ac.uk/snuggletex/units']" mode="cmathml-to-maxima" as="node()*">
    <xsl:apply-templates mode="cmathml-to-maxima"/>
  </xsl:template>

  <xsl:template match="semantics[@definitionURL='http://www.ph.ed.ac.uk/snuggletex/units']/csymbol" mode="cmathml-to-maxima" as="text()">
    <xsl:param name="upconversion-options" as="element(s:upconversion-options)" tunnel="yes"/>
    <xsl:variable name="maxima-units-function" select="s:get-upconversion-option($upconversion-options, 'maximaUnitsFunction')" as="xs:string"/>
    <xsl:value-of select="concat($maxima-units-function, '(&quot;', ., '&quot;)')"/>
  </xsl:template>

  <!-- ************************************************************ -->

  <xsl:template match="emptyset" mode="cmathml-to-maxima" as="text()">
    <xsl:text>{}</xsl:text>
  </xsl:template>

  <xsl:template match="infinity" mode="cmathml-to-maxima" as="text()">
    <!-- NB: This represents real positive infinity only! -->
    <xsl:text>inf</xsl:text>
  </xsl:template>

  <!-- ************************************************************ -->

  <xsl:template match="exponentiale" mode="cmathml-to-maxima" as="text()">
    <xsl:text>%e</xsl:text>
  </xsl:template>

  <xsl:template match="imaginaryi" mode="cmathml-to-maxima" as="text()">
    <xsl:text>%i</xsl:text>
  </xsl:template>

  <xsl:template match="pi" mode="cmathml-to-maxima" as="text()">
    <xsl:text>%pi</xsl:text>
  </xsl:template>

  <xsl:template match="eulergamma" mode="cmathml-to-maxima" as="text()">
    <xsl:text>%gamma</xsl:text>
  </xsl:template>

  <!-- ************************************************************ -->

  <!--
  Helper function to map a flattened <ci/> element name into a suitable
  Maxima input, coping with non-alphanumeric characters as required
  -->
  <xsl:function name="local:map-identifier" as="node()">
    <xsl:param name="element" as="element(ci)"/>
    <xsl:param name="flattened" as="xs:string"/>
    <xsl:variable name="name" select="normalize-space($flattened)" as="xs:string"/>
    <xsl:choose>
      <xsl:when test="matches($name, '^[a-zA-Z_][a-zA-Z0-9_]*$')">
        <!-- Safe to map to a Maxima variable of the same name -->
        <xsl:value-of select="$name"/>
      </xsl:when>
      <xsl:otherwise>
        <!-- Use the identifier dictionary to map it to some Maxima input -->
        <xsl:variable name="maxima-input" as="xs:string?"
          select="$local:identifier-dictionary[.=$name]/@maxima-input"/>
        <xsl:choose>
          <xsl:when test="exists($maxima-input)">
            <xsl:value-of select="$maxima-input"/>
          </xsl:when>
          <xsl:otherwise>
            <!-- Fail: no suitable Maxima input form for identifier -->
            <xsl:copy-of select="s:make-error('UMFG03', $element, ($name))"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:function>

  <!-- Map simple identifiers and unapplied functions over as-is -->
  <xsl:template match="ci[count(node())=1 and text()]" mode="cmathml-to-maxima" as="node()">
    <xsl:copy-of select="local:map-identifier(., string(.))"/>
  </xsl:template>

  <!-- Map subscripts in a reasonable (but rather limited) way -->
  <xsl:template match="ci[count(node())=1 and msub]" mode="cmathml-to-maxima" as="node()+">
    <!-- Drill down into the subscripts (which are PMathML!) -->
    <xsl:apply-templates select="msub" mode="ci-subscripted">
      <xsl:with-param name="ci" select="."/>
    </xsl:apply-templates>
  </xsl:template>

  <xsl:template match="msub[*[1][self::mi]]" mode="ci-subscripted" as="node()+">
    <xsl:param name="ci" as="element(ci)"/>
    <xsl:copy-of select="local:map-identifier($ci, string(*[1]))"/>
    <xsl:text>[</xsl:text>
    <xsl:apply-templates select="*[2]" mode="ci-subscripted">
      <xsl:with-param name="ci" select="$ci"/>
    </xsl:apply-templates>
    <xsl:text>]</xsl:text>
  </xsl:template>

  <xsl:template match="mn" mode="ci-subscripted" as="node()+">
    <xsl:param name="ci" as="element(ci)"/>
    <xsl:copy-of select="node()"/>
  </xsl:template>

  <xsl:template match="mi" mode="ci-subscripted" as="node()">
    <xsl:param name="ci" as="element(ci)"/>
    <xsl:copy-of select="local:map-identifier($ci, string(.))"/>
  </xsl:template>

  <xsl:template match="mfenced" mode="ci-subscripted" as="node()+">
    <xsl:param name="ci" as="element(ci)"/>
    <xsl:for-each select="*">
      <xsl:apply-templates select="." mode="ci-subscripted">
        <xsl:with-param name="ci" select="$ci"/>
      </xsl:apply-templates>
      <xsl:if test="position()!=last()">
        <xsl:text>,</xsl:text>
      </xsl:if>
    </xsl:for-each>
  </xsl:template>

  <xsl:template match="*" mode="ci-subscripted" as="node()">
    <xsl:param name="ci" as="element(ci)"/>
    <!-- Fail: cannot create subscripted variable -->
    <xsl:copy-of select="s:make-error('UMFG04', $ci, ())"/>
  </xsl:template>

  <!-- Don't know what to do in other cases -->
  <xsl:template match="ci" mode="cmathml-to-maxima">
    <xsl:message terminate="yes">
      Did not expect <ci/> element with content <xsl:copy-of select="node()"/>
    </xsl:message>
  </xsl:template>

  <!-- ************************************************************ -->

  <xsl:template match="cn" mode="cmathml-to-maxima" as="text()">
    <xsl:value-of select="if (starts-with(., '-'))
        then concat('(', string(.), ')')
        else string(.)"/>
  </xsl:template>

  <!-- ************************************************************ -->

  <!-- Catch-all for the <apply/> cases we can't or won't handle here -->
  <xsl:template match="apply" mode="cmathml-to-maxima" as="element(s:fail)">
    <xsl:copy-of select="s:make-error('UMFG01', ., (local-name(*[1])))"/>
  </xsl:template>

  <!-- Default catch-all for everything else -->
  <xsl:template match="*" mode="cmathml-to-maxima" as="element(s:fail)">
    <xsl:copy-of select="s:make-error('UMFG00', ., (local-name(.)))"/>
  </xsl:template>

</xsl:stylesheet>
