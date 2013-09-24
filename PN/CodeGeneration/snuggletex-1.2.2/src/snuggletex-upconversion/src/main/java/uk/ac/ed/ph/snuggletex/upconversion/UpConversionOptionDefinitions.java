/* $Id: UpConversionOptionDefinitions.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.upconversion;

import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Set;

/**
 * Defines the various options and values for the up-conversion process.
 * 
 * @since 1.2.0
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
public final class UpConversionOptionDefinitions {
    
    public static final String UPCONVERSION_VARIABLE_NAMESPACE = "upconversion";
    public static final String OPTIONS_VARIABLE_NAME = "options";
    
    public static final String DO_CONTENT_MATHML_NAME = "doContentMathML";
    public static final String DO_MAXIMA_NAME = "doMaxima";
    public static final String ADD_OPTIONS_ANNOTATION_NAME = "addOptionsAnnotation";
    
    public static final String ROUND_BRACKET_HANDLING = "roundBracketHandling";
    public static final String ROUND_FENCE_HANDLING = "roundFenceHandling";
    public static final String SQUARE_FENCE_HANDLING = "squareFenceHandling";
    public static final String CURLY_FENCE_HANDLING = "curlyFenceHandling";
    public static final String EMPTY_FENCE_HANDLING = "emptyFenceHandling";
    
    public static final String MAXIMA_OPERATOR_FUNCTION_NAME = "maximaOperatorFunction";
    public static final String MAXIMA_UNITS_FUNCTION_NAME = "maximaUnitsFunction";
    public static final String MAXIMA_INVERSE_FUNCTION_NAME = "maximaInverseFunction";
    
    public static final Set<String> SYMBOL_ASSUMPTION_TYPES = makeHashSet(new String[] {
            "function",
            "imaginaryNumber",
            "exponentialNumber",
            "constantPi",
            "eulerGamma"
    });
    
    /** Defines the value space for booleans */
    public static final Set<String> BOOLEAN_VALUES = makeHashSet(new String[] {
            "true",
            "false"
    });
    
    /** Defines the values space for handling brackets */
    public static final Set<String> BRACKET_HANDLING_VALUES = makeHashSet(new String[] {
            "grouping",  /* (Treats brackets as grouping only with no special meaning) */
            "error", /* (Causes up-conversion to fail, instead of using default behaviour) */
            "list",
            "set",
            "vector",
    });
    
    /** Defines all of the allowed options for controlling the up-conversion process. */
    public static final LinkedHashMap<String, OptionValueDefinition> OPTION_DEFINITIONS = makeDefinitionMap(new Object[] {
            /* General processing control */
            DO_CONTENT_MATHML_NAME, BOOLEAN_VALUES, "true",
            DO_MAXIMA_NAME, BOOLEAN_VALUES, "true",
            ADD_OPTIONS_ANNOTATION_NAME, BOOLEAN_VALUES, "false",
            
            /* Handling of brackets of various types */
            ROUND_BRACKET_HANDLING, BRACKET_HANDLING_VALUES, "grouping",
            ROUND_FENCE_HANDLING, BRACKET_HANDLING_VALUES, "vector",
            SQUARE_FENCE_HANDLING, BRACKET_HANDLING_VALUES, "list",
            CURLY_FENCE_HANDLING, BRACKET_HANDLING_VALUES, "set",
            EMPTY_FENCE_HANDLING, BRACKET_HANDLING_VALUES, "list",
            
            /* Maxima output control */
            MAXIMA_INVERSE_FUNCTION_NAME, null, "inverse",
            MAXIMA_OPERATOR_FUNCTION_NAME, null, "operator",
            MAXIMA_UNITS_FUNCTION_NAME, null, "units",
    });
    
    /**
     * Defines the default and allowed values for a particular Option.
     *
     * @author  David McKain
     * @version $Revision: 525 $
     */
    public static class OptionValueDefinition {
        
        /** Set of allowed values, or null to indicate "anything allowed" */
        private final Set<String> valueSpace;
        
        /** Default value to use during processing if nothing has been specified */
        private final String defaultValue;
        
        public OptionValueDefinition(final Set<String> valueSpace, final String defaultValue) {
            this.valueSpace = valueSpace;
            this.defaultValue = defaultValue;
        }

        public Set<String> getValueSpace() {
            return valueSpace;
        }

        public String getDefaultValue() {
            return defaultValue;
        }
    }
    
    @SuppressWarnings("unchecked")
    private static LinkedHashMap<String, OptionValueDefinition> makeDefinitionMap(Object[] inputs) {
        LinkedHashMap<String, OptionValueDefinition> result = new LinkedHashMap<String, OptionValueDefinition>();
        for (int i=0; i<inputs.length; ) {
            String propertyName = (String) inputs[i++];
            Set<String> valueSpace = (Set<String>) inputs[i++];
            String defaultValue = (String) inputs[i++];
            result.put(propertyName, new OptionValueDefinition(valueSpace, defaultValue));
        }
        return result;
    }
    
    private static Set<String> makeHashSet(String[] inputs) {
        Set<String> result = new HashSet<String>();
        for (String input : inputs) {
            result.add(input);
        }
        return result;
    }
}
