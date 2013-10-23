# 
#  Copyright (C) 2007  Smithsonian Astrophysical Observatory
#
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along
#  with this program; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#


#
# Makefile for Sherpa Python package
#

MK_TOP = ../..
include $(MK_TOP)/Makefile.master
include $(MK_TOP)/include/Makefile.scidev

CP_F = cp -f
CP_PF = cp -pf
RM_F = /bin/rm -f
RM_RF = /bin/rm -rf

SL_HELP_DIR  = sherpa/docs
SL_HELP_FILE = sherpa.hlp
DOC_DIR = $(SHARE_DIR)/slsh/local-packages/help

# Take F95 variable, which defines full path to Fortran compiler
# Get directory path, substitute trailing "bin/" with "lib"
# Assumes F95 from Makefile.scidev defines full path to compiler

FORTRAN_LIB=$(subst -l,,$(word 2,$(F95LIB)))
FORTRAN_LIB_DIR=$(subst -L,,$(firstword $(F95LIB)))

SETUP_VARS = xspec_library_dir=$(XSPEC12)/base/lib \
	cfitsio_lib=$(CFITSIO_LIB_NAME) cfitsio_library_dir=$(CFITSIO)/lib\
	reg_library_dir=$(LIB_DIR) reg_include_dir=$(INC_DIR) \
	fftw_library_dir=$(FFTW)/lib fftw_include_dir=$(FFTW)/include \
	wcs_library_dir=$(WCS_LIB) wcs_include_dir=$(WCS_INC) \
	fortran_library_dir=$(FORTRAN_LIB_DIR) fortran_lib=$(FORTRAN_LIB)

# We really do need to export the compilers and compiler flags to the
# shell so that the Python setup.py script picks up the compilers and
# flags defined in Makefile.scidev.  Don't remove the export from
# SHERPA_SETUP_PY unless you have figured out another way to do this.

# As of CIAO 4.1, certain compiler flags (such -m32) may be appended
# to the compiler targets CC, CXX etc.  We need to peel those flags
# out and add them to CFLAGS instead.  The method below works
# for both GNU and Solaris make.  SMD 12/04/2008

SHERPA_SETUP_PY = \
	CC_COMP="`echo $(CC) | awk '{print $$1}'`" ; \
	CC_EXTRA_FLAGS="`echo $(CC) | awk '{ for (i = 2; i <= NF; i++) print $$i }'`" ; \
	CXX_COMP="`echo $(CCC) | awk '{print $$1}'`" ; \
	CXX_EXTRA_FLAGS="`echo $(CCC) | awk '{ for (i = 2; i <= NF; i++) print $$i}'`" ; \
	F77_COMP="`echo $(F77) | awk '{print $$1}'`" ; \
	F77_EXTRA_FLAGS="`echo $(F77) | awk '{ for (i = 2; i <= NF; i++) print $$i }'`" ; \
	F90_COMP="`echo $(F95) | awk '{print $$1}'`" ; \
	F90_EXTRA_FLAGS="`echo $(F95) | awk '{ for (i = 2; i <= NF; i++) print $$i }'`" ; \
	CC=$$CC_COMP ; CXX=$$CXX_COMP ; \
	F77=$$F77_COMP ; F90=$$F90_COMP ; \
	CFLAGS="$(CFLAGS:-Xc=) $$CC_EXTRA_FLAGS" ; \
	CXXFLAGS="$(CCFLAGS) $$CXX_EXTRA_FLAGS"; \
	F77FLAGS="$(FFLAGS) $$F77_EXTRA_FLAGS" ; \
	F90FLAGS="$(FFLAGS) $$F90_EXTRA_FLAGS" ; \
	DFLAGS="$(DFLAGS)" ; OCFLAGS="$(OCFLAGS)" ; \
	FOPT="$(OFLAGS)" ; FDEBUG="$(DFLAGS)"; \
	PYTHONPATH="$(PYTHON)/lib/python$(PYTHONVER)/site-packages:$(NUMPYMOD)" ; \
	if [ `$(UNAME)` = "Darwin" ] ; then \
	   DYLD_FRAMEWORK_PATH="$(PYTHON)/Library/Frameworks" ; \
	   DYLD_LIBRARY_PATH="$(PYTHONLDLIBS):$(DYLD_LIBRARY_PATH)" ; \
	   export DYLD_LIBRARY_PATH DYLD_FRAMEWORK_PATH ; \
	else \
	   LD_LIBRARY_PATH="$(PYTHONLDLIBS):$(LD_LIBRARY_PATH)" ; \
	   export LD_LIBRARY_PATH ; \
	fi ; \
	export CC CXX F77 F90 FDEBUG FOPT CFLAGS CXXFLAGS F77FLAGS F90FLAGS ; \
	export DFLAGS OCFLAGS PYTHONPATH ; \
	if [ "$$CC_EXTRA_FLAGS" != "" ] ; then \
	   LDFLAGS="$$CC_EXTRA_FLAGS $(SHARE_FLAGS)" ; \
	   export LDFLAGS ; \
	fi ; \
	cd sherpa ; $(PYTHONEXE) setup.py $(SETUP_VARS)

# Wrapper script to call IPython with Sherpa profile
SHERPA_SCRIPT = sherpa
SHERPA_SCRIPT_PY = sherpa_py
SHERPA_SCRIPT_SL = sherpa_slsh.sl

all: build

build:
	@echo ""
	@echo "   /-------------------------------------------------------\ "
	@echo "   |                   Building Sherpa                     | "
	@echo "   \-------------------------------------------------------/ "
	@echo ""
	$(SHERPA_SETUP_PY) build
	cd tests; $(MAKE) $(MKMACROS)
	cd doc; $(MAKE) $(MKMACROS)
	cd scripts; $(MAKE) $(MKMACROS)

install:
	@echo "   /-------------------------------------------------------\ "
	@echo "   |                  Installing Sherpa                    | "
	@echo "   \-------------------------------------------------------/ "
	@echo ""
##	$(SHERPA_SETUP_PY) install --install-lib=$(PYTHONMOD)
	$(SHERPA_SETUP_PY) install --prefix=$(ASCDS_BLD_INSTALL)
	if [ -n "$(SHARE_DIR)" -a -d "$(SHARE_DIR)" ]; then \
	  $(CP_PF) sherpa/slang/sherpa.sl $(SHARE_DIR)/slsh/local-packages; \
	  $(CP_PF) sherpa/sherpa/$(SHERPA_SCRIPT_SL) $(SHARE_DIR)/slsh/local-packages; \
	  chmod 644 $(SHARE_DIR)/slsh/local-packages/sherpa.sl; \
	  chmod 644 $(SHARE_DIR)/slsh/local-packages/$(SHERPA_SCRIPT_SL); \
	fi;
	if [ -n "$(BIN_DIR)" -a -d "$(BIN_DIR)" ]; then \
	  $(CP_PF) sherpa/sherpa/$(SHERPA_SCRIPT) $(BIN_DIR); \
	  $(CP_PF) sherpa/sherpa/$(SHERPA_SCRIPT_PY) $(BIN_DIR); \
	  chmod 755 $(BIN_DIR)/$(SHERPA_SCRIPT); \
	  chmod 755 $(BIN_DIR)/$(SHERPA_SCRIPT_PY); \
	  if [ `echo $(BIN_DIR) | grep -c "/binexe" | tr -d '\n'` -eq 1 -a -f "$(MK_TOP)/config/runciao" ]; then \
	    $(CP_PF) $(MK_TOP)/config/runciao $(ASCDS_BLD_INSTALL)/bin/$(SHERPA_SCRIPT); \
	    $(CP_PF) $(MK_TOP)/config/runciao $(ASCDS_BLD_INSTALL)/bin/$(SHERPA_SCRIPT_PY); \
	    chmod 755 $(ASCDS_BLD_INSTALL)/bin/$(SHERPA_SCRIPT); \
	    chmod 755 $(ASCDS_BLD_INSTALL)/bin/$(SHERPA_SCRIPT_PY); \
	  fi; \
	fi;
	if [ -n "$(DOC_DIR)" -a -d "$(DOC_DIR)" ]; then \
	  $(CP_PF) $(SL_HELP_DIR)/$(SL_HELP_FILE) $(DOC_DIR); \
	  chmod 644 $(DOC_DIR)/$(SL_HELP_FILE); \
	fi;
	cd tests; $(MAKE) $(MKMACROS) install
	cd doc; $(MAKE) $(MKMACROS) install
	cd scripts; $(MAKE) $(MKMACROS) install

uninstall:
	@echo "   /-------------------------------------------------------\ "
	@echo "   |                 Uninstalling Sherpa                   | "
	@echo "   \-------------------------------------------------------/ "
	@echo ""
	if [ -n "$(PYTHONMOD)" -a -d "$(PYTHONMOD)" \
	     -a -d "$(PYTHONMOD)/sherpa" ]; then \
	  $(RM_RF) $(PYTHONMOD)/sherpa; \
	fi;
	if [ -n "$(SHARE_DIR)" -a -d "$(SHARE_DIR)" ]; then \
	  $(RM_F) $(SHARE_DIR)/slsh/local-packages/sherpa.sl; \
	  $(RM_F) $(SHARE_DIR)/slsh/local-packages/$(SHERPA_SCRIPT_SL); \
	fi;
	if [ -n "$(BIN_DIR)" -a -d "$(BIN_DIR)" ]; then \
	  $(RM_F) $(BIN_DIR)/$(SHERPA_SCRIPT); \
	  $(RM_F) $(BIN_DIR)/$(SHERPA_SCRIPT_PY); \
	  if [ `echo $(BIN_DIR) | grep -c "/binexe" | tr -d '\n'` -eq 1 ]; then \
	    $(RM_F) $(ASCDS_BLD_INSTALL)/bin/$(SHERPA_SCRIPT); \
	    $(RM_F) $(ASCDS_BLD_INSTALL)/bin/$(SHERPA_SCRIPT_PY); \
	  fi; \
	fi;
	if [ -n "$(DOC_DIR)" -a -d "$(DOC_DIR)" ]; then \
	  $(RM_F) $(DOC_DIR)/$(SL_HELP_FILE); \
	fi;
	cd tests; $(MAKE) $(MKMACROS) uninstall
	cd doc; $(MAKE) $(MKMACROS) uninstall
	cd scripts; $(MAKE) $(MKMACROS) uninstall

dist:
	@echo "   /-------------------------------------------------------\ "
	@echo "   |          Creating Sherpa source distribution          | "
	@echo "   \-------------------------------------------------------/ "
	@echo ""
	$(SHERPA_SETUP_PY) sdist

clean:
	@echo "   /-------------------------------------------------------\ "
	@echo "   |           Cleaning Sherpa build directories           | "
	@echo "   \-------------------------------------------------------/ "
	@echo ""
	$(RM_RF) sherpa/build sherpa/dist sherpa/MANIFEST
	cd tests; $(MAKE) $(MKMACROS) clean
	cd doc; $(MAKE) $(MKMACROS) clean
	cd scripts; $(MAKE) $(MKMACROS) clean
