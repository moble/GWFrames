/* $Id: ContextInitialiser.java 556 2010-04-30 12:47:23Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.webapp;

import uk.ac.ed.ph.snuggletex.utilities.SaxonTransformerFactoryChooser;
import uk.ac.ed.ph.snuggletex.utilities.SimpleStylesheetCache;
import uk.ac.ed.ph.snuggletex.utilities.StylesheetManager;

import javax.servlet.ServletContext;
import javax.servlet.ServletContextEvent;
import javax.servlet.ServletContextListener;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Fairly typical {@link ServletContextListener} that just sets up a few shared resources
 * and sticks them in the {@link ServletContext} for access by servlets.
 *
 * @author  David McKain
 * @version $Revision: 556 $
 */
public final class ContextInitialiser implements ServletContextListener {
    
    private static final Logger logger = LoggerFactory.getLogger(ContextInitialiser.class);
    
    public static final String SNUGGLETEX_VERSION_PROPERTY_NAME = "snuggletex.version";
    public static final String MAVEN_SITE_URL_PROPERTY_NAME = "maven.site.url";
    
    public static final String STYLESHEET_MANAGER_ATTRIBUTE_NAME = "stylesheetManager";
    
    public void contextInitialized(ServletContextEvent servletContextEvent) {
        ServletContext servletContext = servletContextEvent.getServletContext();
        
        /* Create and store StylesheetManager, hard-coded to use Saxon with caching
         * turned on as we're loading via the ClassPath so there's no point trying to turn
         * it off. */
        StylesheetManager stylesheetManager = new StylesheetManager();
        stylesheetManager.setStylesheetCache(new SimpleStylesheetCache());
        stylesheetManager.setTransformerFactoryChooser(SaxonTransformerFactoryChooser.getInstance());
        servletContext.setAttribute(STYLESHEET_MANAGER_ATTRIBUTE_NAME, stylesheetManager);
        logger.info("Context initialised");
    }
    
    public void contextDestroyed(ServletContextEvent servletContextEvent) {
        logger.info("Context destroyed");
    }
}