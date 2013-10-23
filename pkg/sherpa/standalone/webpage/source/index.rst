.. Sherpa documentation master file, created by
   sphinx-quickstart on Wed Nov  4 16:22:54 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Intro
=====


Sherpa is a modeling and fitting application for Python. It contains a
powerful language for combining simple models into complex expressions
that can be fit to the data using a variety of statistics and
optimization methods. It is easily extensible to include user models,
statistics and optimization methods.


What can you do with Sherpa?
----------------------------

 * Fit 1D (multiple) data including: spectra, surface brightness
   profiles, light curves, general ASCII arrays 

 * Fit 2D images/surfaces in Poisson/Gaussian regime

 * Build complex model expressions

 * Import and use your own models

 * Use appropriate statistics for modeling Poisson or Gaussian data

 * Import new statistics, with priors if required by analysis

 * Visualize a parameter space with simulations or using 1D/2D cuts of
   the parameter space 

 * Calculate confidence levels on the best fit model parameters

 * Choose a robust optimization method for the fit: Levenberg-Marquardt [lm]_,
   Nelder-Mead Simplex [nm]_ or Monte Carlo/Differential Evolution [mc]_.


Documentation
=============

For detailed documentation see:

http://cxc.cfa.harvard.edu/sherpa


New Features in Version 4.4.1 Patch
===================================

 * Support for XSPEC 12.7.1 has been added.  Interfaces to the
   following additive models has been provided:

   * xseplogpar
   * xsgadem
   * xsvgadem
   * xslogpar
   * xsoptxagn
   * xsoptxagnf
   * xspexmon

 * See the release notes NOTES_

    
Features In Version 4.4
=======================

 * 24 new optical models for astronomical data, in module
   sherpa.astro.optical.

 * A new template model, for comparing 1-D data to libraries of
   precalculated models.

 * The pyBLoCXS module is now integrated into Sherpa. pyBLoCXS is a
   sophisticated Markov chain Monte Carlo (MCMC) based algorithm
   designed to carry out Bayesian Low-Count X-ray Spectral (BLoCXS)
   analysis in the Sherpa environment.

 * There is a new PPP statistic function, plot_pvalue(), to compute
    the likelihood ratios of the null and the alternative models.

 * Many significant improvements to functions ``set_full_model()`` and
   ``set_full_bkg_model()``, allowing users to explicitly define
   instrument responses and convolutions that are applied to specified
   model components. Instrument responses, PSF models and table models
   may be used in model expressions defined with these functions.

 * Many bug fixes.

 * See a complete list of features in the release notes NOTES4.4.0_.


Dependencies
============

The new Sherpa can be downloaded, built, installed, and used independently of
CIAO.  I/O support and plotting can be supplemented using PyFITS and
matplotlib.  The current (Sherpa 4.4.1) source tarball_ has the following
prerequisites:

.. _tarball: http://cxc.cfa.harvard.edu/contrib/sherpa/index.html#download


Compilers:
----------

 * GCC compilers with gfortran, versions 3.4.5 or later

 * Sun Studio with f95 WS 10 or later

 * Xcode with gfortran 4.2.3 or later (disk images of gfortran for Mac can be
   found at http://r.research.att.com/tools/)

    - Mac OSX 10.7 Snow Leopard, Xcode 3.2 or later

    - Mac OSX 10.6 Snow Leopard, Xcode 3.2 or later

    - Mac OSX 10.5 Leopard, Xcode 3.1, 3.1.1, 3.1.2, or 3.1.3 (Intel processors
      preferred)


Required:
---------

 :Python 2.6 or 2.7 [py]_: Python is an interpreted, interactive,
  object-oriented, extensible programming language.


 :NumPy 1.5.1 or later [np]_: Powerful N-dimensional array object for Python.


 :FFTW 3.3 or later [fftw]_: Fast fourier transform library speeds up Sherpa
  PSF model convolution. Be sure to build with "--enable-float" also.


Optional:
---------

 :CIAO dmRegion Library 4.4 [reg]_: The Sherpa region module accepts CIAO region
  filtering syntax parsed by the dmRegion library.  Regions are two dimensional
  filters that can be used to include or exclude data points in a 2D data set.


 :DS9 5.6 or later [ds9]_: DS9 is the SAO imager used for astronomical imaging
  and data visualization.


 :IPython 0.10 [ipy]_: IPython provides a nicer interpreter interface with readline
  support compared to regular Python.


 :matplotlib 1.0.0 or later preferred [mpl]_: Sherpa can use matplotlib for line and
  contour plotting and image visualization.  For smooth behavior in Sherpa, set
  the configuration switch "interactive=True" in ~/.matplotlib/matplotlibrc.


 :PyFITS 2.3.1 or later preferred [pyfits]_: Sherpa uses PyFITS as a FITS file reader
  backend.


 :WCS 3.7.0 or later [wcs]_:  The World Coordinate System library provides the
  transforms necessary to support the physical and world coordinate systems
  found in FITS file headers.


 :XPA 2.1.9 or later [xpa]_: Sherpa uses XPA applications (e.g. xpaset, xpaget)
  to communicate with DS9.


 :XSPEC 12.7.0 [xspec]_: The XSPEC spectral models are available in Sherpa by
  linking to certain dynamic libraries found in an XSPEC installation.



Download
========

Current version **4.4.1** June 27, 2012

 ==========================  ==================================================
 source tarball              sherpa-4.4.1.tar.gz_
 Release Notes               NOTES_
 ==========================  ==================================================

For previous versions see downloads_

.. _downloads: http://cxc.cfa.harvard.edu/contrib/sherpa/downloads

.. _sherpa-4.4.1.tar.gz: http://cxc.cfa.harvard.edu/contrib/sherpa/sherpa-4.4.1.tar.gz

.. _NOTES: http://cxc.cfa.harvard.edu/contrib/sherpa/NOTES-4.4.1.txt

.. _NOTES4.4.0: http://cxc.cfa.harvard.edu/contrib/sherpa/NOTES-4.4.0.txt


Build and Install
=================

To build and install the package, do the following ::

  $ tar xzf sherpa-4.4.1.tar.gz
  $ cd sherpa-4.4.1
  $ python setup.py [config-vars] install --prefix=<dest-dir>


config-vars is an optional list of arguments in the format var=value that 
specifies where to find prerequisites required at build time. The following 
variables can be set

===================  ================================================================================
Variable 	     Default
===================  ================================================================================
fftw_library_dir     /usr/lib (only needed if in alternate location, e.g. /usr/lib64)
fftw_include_dir     /usr/include (only needed if in alternate location)
wcs_library_dir      None (if not given, World Coordinate System, WCS, module is not built)
wcs_include_dir      None (if not given, World Coordinate System, WCS, module is not built)
reg_library_dir      None (if not given, Region 2D filtering module is not built)
reg_include_dir      None (if not given, Region 2D filtering module is not built)
fortran_library_dir  None (may be needed to find libgfortran.so or libf95.a depending on how XSPEC was built)
fortran_lib          None ('gfortran', or 'f95' depending on fcompiler )
xspec_library_dir    None (if not given, XSPEC module is not built)
cfitsio_lib 	     cfitsio  (or 'cfitsio_3.25' depending on version found in XSPEC installation)
cfitsio_library_dir  <xspec_library_dir>
===================  ================================================================================

**NOTE:** The version of gfortran used to compile XSPEC and Sherpa must be
identical.  This issue with gfortran has been known to cause compiler errors.


For example, to use the FFTW in /soft/fftw and the XSPEC library in
/opt/local/headas/lib, the command to install Sherpa would be ::

  $ python setup.py \
     fftw_library_dir=/soft/fftw/lib \
     fftw_include_dir=/soft/fftw/include \
     xspec_library_dir=/opt/local/headas/lib \
     cfitsio_library_dir=/opt/local/headas/lib \
     install ...

The setup.py script distributed with Sherpa uses the standard Python distutils
package. For more information on using it, see Installing Python Modules.


Platforms
---------

Notes about building Sherpa from source on various platforms

 * **Mac OSX 10.6 Snow Leopard**

   On Mac 10.6, there are various different flavors of Python.

   - The Python 2.6 that comes with the installation of OSX 10.6 is typically called
     *Apple* Python, located here ::

        /System/Library/Frameworks/Python.framework/Versions/Current

     This version is built as a three-way binary: i386, x86_64, and ppc.


   - The Python 2.6 disk image from www.python.org for OSX 10.6 is typically
     called *python.org* Python.  It installs here ::

        /Library/Frameworks/Python.framework/Versions/Current  

     This version is built as a two-way binary: i386 and ppc.  Using this disk
     image of Python 2.6 can be problematic when building Sherpa from source on
     10.6 because the default Mach-O architecture is x86_64.

   *Apple* Python or a Python installation built from source are the recommended
   options for building 64-bit Sherpa from source on Mac 10.6.  The version of
   *Apple* Python is usually out-of-date, so a source build of Python + NumPy is
   preferred.  Users are always free to investigate building *fat* binaries.  It
   can be done!


 * **Mac OSX 10.5 Leopard**
 
   On Mac 10.5, there are again different flavors of Python.

   - The Python 2.5 that comes with the installation of OSX 10.5 is typically
     called *Apple* Python, located here ::

       /System/Library/Frameworks/Python.framework/Versions/Current

     This version is built as a two-way binary: i386 and ppc.


   - The Python 2.6 disk image from www.python.org for OSX 10.5 is typically called
     *python.org* Python.  It installs here ::

       /Library/Frameworks/Python.framework/Versions/Current  

     This version is built as a two-way binary: i386 and ppc.  

   The default Mach-O architecture is i386, so a *python.org* Python
   installation or a Python installation from source are the recommended
   options.


.. * **Fedora**

.. * **Ubuntu**

.. * **CentOS**

.. * **Red Hat**


Configuration Files
===================

Sherpa comes with a configuration file "sherpa.rc". This file is located in
<dest-dir>/lib/python2.X/site-packages/sherpa/sherpa.rc.  For a stand-alone
installation, this file should be copied from source in to the home directory as
"~/.sherpa.rc". Be sure to indicate the IO and Plotting back-ends as "pyfits"
and "pylab" depending on configuration.

matplotlib comes with a configuration file "matplotlibrc". For smooth behavior
with Sherpa, be sure to indicate "interactive=True" in
~/.matplotlib/matplotlibrc.


Environment Variables
=====================

If Sherpa is installed in a non-standard location, be sure to update the 
following variables accordingly. ::

    PYTHONPATH

    HEADAS (for XSPEC spectral tables) e.g. export HEADAS=<xspec_install>/spectral


References
==========

.. [ds9] http://hea-www.harvard.edu/RD/ds9

.. [fftw] http://www.fftw.org

.. [ipy] http://ipython.scipy.org

.. [lm] Lecture Notes in Mathematics 630: Numerical Analysis, G.A. Watson (Ed.), Springer-Verlag: Berlin, 1978, pp. 105-116

.. [mc] Differential Evolution: A Simple and Efficient Adaptive Scheme for Global Optimization over Continuous Spaces, J. Global Optimization 11, Storn, R. and Price, K., 1997, pp. 341-359 http://www.icsi.berkeley.edu/~storn/code.html

.. [mpl] Hunter, JD (2007).  *Matplotlib: A 2D graphics environment.*
   Computing in Science and Engineering. 9:
   90-95. http://matplotlib.sourceforge.net.

.. [nm] Computer Journal, J.A. Nelder and R. Mead, 1965, vol 7, pp. 308-313.  Jeffrey C. Lagarias, James A. Reeds, Margaret H. Wright, Paul E. Wright "Convergence Properties of the Nelder-Mead Simplex Algorithm in Low Dimensions", SIAM Journal on Optimization,Vol. 9, No. 1 (1998), pp. 112-147. http://citeseer.ist.psu.edu/3996.html.  Wright, M. H. (1996) "Direct Search Methods: Once Scorned, Now Respectable" in Numerical Analysis 1995 (Proceedings of the 1995 Dundee Biennial Conference in Numerical Analysis) (D.F. Griffiths and G.A. Watson, eds.), 191-208, Addison Wesley Longman, Harlow, United Kingdom. http://citeseer.ist.psu.edu/155516.html

.. [np] http://numpy.scipy.org

.. [py] http://www.python.org

.. [pyfits] http://www.stsci.edu/resources/software_hardware/pyfits

.. [reg] http://cxc.harvard.edu/ciao

.. [wcs] http://tdc-www.harvard.edu/wcstools

.. [xpa] http://hea-www.harvard.edu/saord/xpa

.. [xspec] http://heasarc.gsfc.nasa.gov/docs/xanadu/xspec

.. toctree::
   :maxdepth: 2
