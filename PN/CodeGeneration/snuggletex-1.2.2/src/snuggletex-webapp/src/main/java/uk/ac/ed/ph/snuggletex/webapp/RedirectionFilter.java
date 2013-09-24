/* $Id: RedirectionFilter.java 525 2010-01-05 14:07:36Z davemckain $
 *
 * Copyright (c) 2010, The University of Edinburgh.
 * All Rights Reserved
 */
package uk.ac.ed.ph.snuggletex.webapp;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.servlet.Filter;
import javax.servlet.FilterChain;
import javax.servlet.FilterConfig;
import javax.servlet.ServletException;
import javax.servlet.ServletRequest;
import javax.servlet.ServletResponse;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Simple little {@link Filter} that allows us to apply redirections (specified in web.xml)
 * to cater for the case where pages move about or get better names.
 *
 * @author  David McKain
 * @version $Revision: 525 $
 */
public final class RedirectionFilter implements Filter {
    
    private static final Logger logger = LoggerFactory.getLogger(RedirectionFilter.class);
    
    private final LinkedHashMap<Pattern, String> redirectionMap = new LinkedHashMap<Pattern, String>();
    
    public void init(FilterConfig filterConfig) {
        String redirectionData = filterConfig.getInitParameter("redirections");
        if (redirectionData!=null) {
            redirectionMap.clear();
            String[] fields = redirectionData.split("\\s+");
            Pattern pattern;
            String target;
            for (int i=0; i<fields.length; ) {
                pattern = Pattern.compile(fields[i++]);
                target = fields[i++];
                logger.info("Registering redirect " + pattern + " => " + target);
                redirectionMap.put(pattern, target);
            }
        }
    }
    
    public void doFilter(ServletRequest request, ServletResponse response, FilterChain filterChain)
            throws IOException, ServletException {
        HttpServletRequest httpRequest = (HttpServletRequest) request;
        
        /* See if within-context URL matches one of our patterns. If so, redirect */
        String requestUrl = WebUtilities.getWithinContextRequestUrl(httpRequest);
        for (Entry<Pattern, String> mapEntry : redirectionMap.entrySet()) {
            Pattern matchPattern = mapEntry.getKey();
            String replacement = mapEntry.getValue();
            Matcher matcher = matchPattern.matcher(requestUrl);
            if (matcher.matches()) {
                HttpServletResponse httpResponse = (HttpServletResponse) response;
                httpResponse.sendRedirect(httpRequest.getContextPath() + matcher.replaceFirst(replacement));
                return;
            }
        }
        /* If still here, then let request proceed as normal */
        filterChain.doFilter(request, response);
    }
    
    public void destroy() {
        /* (Nothing to do) */
    }
}
