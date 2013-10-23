#-----------------------------------------------------------------------
# Makefile for SciAlg Sherpa scripts 
#
# To make 
#	% make all
# 
# To install type:
#	% make install
#
#-----------------------------------------------------------------------
MK_TOP = ../../..
include $(MK_TOP)/Makefile.master
include $(MK_TOP)/include/Makefile.scidev

PYTHONSCRIPTS = \
	lev3_f2cratio.py \
	lev3_iss.py \
	lev3_mho.py \
	lev3_specfit.py \
	modelflux \
	srcextent

PAR_FILES = \
	lev3_specfit.par \
	lev3_f2cratio.par \
	lev3_mho.par \
	lev3_iss.par \
	modelflux.par \
	srcextent.par

EXEC	= \
	lev3_specfit.py \
	lev3_f2cratio.py \
	lev3_mho.py \
	lev3_iss.py \
	modelflux \
	srcextent


#PIPES	= 


#INTERPRETED_FILES = 


#CAL_FILES = 


MAKECLEAN =

MAKEALL = $(PYTHONSCRIPTS)

include $(MK_TOP)/Makefile.all

#-----------------------------------------------------------------------
# 			MAKEFILE DEPENDENCIES	
#-----------------------------------------------------------------------

install: announce1
	cd tests; $(MAKE) $(MKMACROS) install


announce1:
	@echo "   /---------------------------------------------------------\ "
	@echo "   |                scialg sherpa scripts                    | "
	@echo "   \---------------------------------------------------------/ "

#
