  Building the Python distribution installer

  Python-2.7.2

  These scripts build all of the standard library dependencies
  for Python from source.

    db-4.7.25
    ncurses-5.9
    openssl-0.9.8e
    bzip2-1.0.6
    readline-6.2
    sqlite-autoconf-3070701
    tcl8.5.10
    tk8.5.10
    zlib-1.2.5

  The scripts also build Python and the Python dependencies
  from source required for Sherpa.

    Python-2.7.2   
    numpy-1.6.1
    pyfits-3.0
    ipython-0.11
    matplotlib-1.1.0
      freetype-2.4.6
      libpng-1.2.46
      
  The scripts also build the follow Python package managing
  modules from source

    setuptools-0.6c11
    distribute-0.6.21
    enstaller-4.4.1-1


  Optional Python modules built from source include

    APLpy-0.9.6
    ATpy-0.9.5.1
    asciitable-0.7.1
    pywcs-1.10-4.7
    pyregion-1.0.1
      pyparsing-1.5.6
    vo-0.7.2
    

  The final product is a shell script installer that will place
  the complete Python distribution anywhere.


  * Building Linux

  % ./make_driver.sh --build=Linux64 --version=2.7.2

  After about 18 minutes...

  Builds Python-2.7.2-Linux64.sh

  * Installing on Linux

  % ./Python-2.7.2-Linux64.sh

  Follow the terminal instructions.


  Enjoy!
