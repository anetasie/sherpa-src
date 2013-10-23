Build the Sherpa stand-alone documentation

Edit the main RST file 

     source/index.rst


Compile the HTML, JavaScript, CSS with

% Make html

Built documentation appears in 

      build


To keep track of individual tarball downloads using Google Analytics, I edit the
HTML file to explicitly include JavaScript "onClick" functions for analytics

     index.html


e.g.

     <a class="reference external" href="http://cxc.harvard.edu/contrib/sherpa/sherpa-4.3.0.tar.gz" onClick="javascript: _gaq.push(['_trackPageview', '/sherpa-4.3.0.tar.gz']);">sherpa-4.3.0.tar.gz</a>


     <a class="reference external" href="http://cxc.harvard.edu/contrib/sherpa/sherpa-4.3.0-py2.6-python.org-macosx10.6-intel.dmg" onClick="javascript: _gaq.push(['_trackPageview', '/sherpa-4.3.0-py2.6-python.org-macosx10.6-intel.dmg']);">sherpa-4.3.0-py2.6-python.org-macosx10.6-intel.dmg</a>


     <a class="reference external" href="http://cxc.harvard.edu/contrib/sherpa/sherpa-4.3.0-py2.7-python.org-macosx10.6-intel.dmg" onClick="javascript: _gaq.push(['_trackPageview', '/sherpa-4.3.0-py2.7-python.org-macosx10.6-intel.dmg']);">sherpa-4.3.0-py2.7-python.org-macosx10.6-intel.dmg</a>

     
     <a class="reference external" href="http://cxc.harvard.edu/contrib/sherpa/NOTES-4.3.0.txt" onClick="javascript: _gaq.push(['_trackPageview', '/NOTES-4.3.0.txt']);">NOTES</a>


The google analytics ID tag for username 'ascdsdev' is

    UA-27024519-1


Join new UNIX group

     cxcweb_contrib


Copy the contents of 'build' over to the webserver location

     
     /proj/web-cxc-dmz/htdocs/contrib/sherpa
