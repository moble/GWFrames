/* $Id: UpConversionOptions.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.upconversion;

import static uk.ac.ed.ph.snuggletex.upconversion.UpConversionOptionDefinitions.OPTION_DEFINITIONS;
import static uk.ac.ed.ph.snuggletex.upconversion.UpConversionOptionDefinitions.SYMBOL_ASSUMPTION_TYPES;

import uk.ac.ed.ph.snuggletex.upconversion.UpConversionOptionDefinitions.OptionValueDefinition;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.w3c.dom.Element;

/**
 * Encapsulates the various options and assumptions that can be used when running the
 * up-conversion process.
 * <p>
 * (These can be set explicitly via this Java API, or using custom LaTeX macros as well.)
 *
 * @since 1.2.0
 * 
 * @author  David McKain
 * @version $Revision: 525 $
 */
public final class UpConversionOptions {
    
    private final Map<String, String> specifiedOptionMap;
    
    private final Map<ElementWrapper, String> symbolAssumptions;
    
    public UpConversionOptions() {
        this.specifiedOptionMap = new HashMap<String, String>();
        this.symbolAssumptions = new HashMap<ElementWrapper, String>();
    }
    
    //-------------------------------------------------------------
    
    public boolean isOptionSpecified(final String name) {
        return specifiedOptionMap.containsKey(name);
    }
    
    public Set<String> getSpecifiedOptionNames() {
        return specifiedOptionMap.keySet();
    }
    
    public String getSpecifiedOptionValue(final String name) {
        if (!OPTION_DEFINITIONS.containsKey(name)) {
            throw new IllegalUpconversionOptionException(UpConversionErrorCode.UAEOP0, name);
        }
        return specifiedOptionMap.get(name);
    }
    
    public String getEffectiveOptionValue(final String name, final boolean applyDefault) {
        if (!OPTION_DEFINITIONS.containsKey(name)) {
            throw new IllegalUpconversionOptionException(UpConversionErrorCode.UAEOP0, name);
        }
        return specifiedOptionMap.containsKey(name) ? specifiedOptionMap.get(name) : (applyDefault ? OPTION_DEFINITIONS.get(name).getDefaultValue() : null);
    }
    
    /**
     * @throws IllegalUpconversionOptionException
     */
    public void setSpecifiedOption(final String name, final String value) {
        OptionValueDefinition valueDefinition = OPTION_DEFINITIONS.get(name);
        if (valueDefinition==null) {
            throw new IllegalUpconversionOptionException(UpConversionErrorCode.UAEOP0, name);
        }
        Set<String> valueSpace = valueDefinition.getValueSpace();
        if (valueSpace!=null && !valueSpace.contains(value)) {
            throw new IllegalUpconversionOptionException(UpConversionErrorCode.UAEOP1, name, value);
        }
        specifiedOptionMap.put(name, value);
    }
    
    /**
     * @throws IllegalUpconversionOptionException
     */
    public void clearSpecifiedOption(final String name) {
        if (specifiedOptionMap.containsKey(name)) {
            specifiedOptionMap.remove(name);
        }
        else {
            throw new IllegalUpconversionOptionException(UpConversionErrorCode.UAEOP2, name);
        }
    }
    
    //-------------------------------------------------------------
    
    public Set<ElementWrapper> getAssumedSymbols() {
        return symbolAssumptions.keySet();
    }
    
    public String getSymbolAssumptionType(Element element) {
        return getSymbolAssumptionType(new ElementWrapper(element));
    }
    
    public String getSymbolAssumptionType(ElementWrapper elementWrapper) {
        return symbolAssumptions.get(elementWrapper);
    }
    
    /**
     * @throws IllegalUpconversionOptionException
     */
    public void assumeSymbol(Element element, String assumptionType) {
        if (!SYMBOL_ASSUMPTION_TYPES.contains(assumptionType)) {
            throw new IllegalUpconversionOptionException(UpConversionErrorCode.UAESY0, assumptionType);
        }
        symbolAssumptions.put(new ElementWrapper(element), assumptionType);
    }
    
    /**
     * @throws IllegalUpconversionOptionException
     */
    public void unassumeSymbol(Element element) {
        ElementWrapper wrapper = new ElementWrapper(element);
        if (symbolAssumptions.containsKey(wrapper)) {
            symbolAssumptions.remove(wrapper);
        }
        else {
            throw new IllegalUpconversionOptionException(UpConversionErrorCode.UAESY2);
        }
    }
}
