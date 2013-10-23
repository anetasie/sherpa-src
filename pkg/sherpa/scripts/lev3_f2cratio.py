#! /usr/bin/env python

# 
#  Copyright (C) 2009  Smithsonian Astrophysical Observatory
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



#---------------------------------------------------------------
#
# Level 3 pipeline script
# Start in March, 2007
#
# calculate the ratio of model flux to model counts
#
# source model - 1: absorption * powlaw  
# source model - 2: absorption * blackbody
# all model parameters are fixed
# Input:  par file, incl PI/PHA source data file
# output: par file, incl ratio_model1 and ratio_model2
#
#----------------------------------------------------------------
#


#---------------------------------------------------
#
# import Python modules
#
#--------------------------------------------------
import sys
import os
import os.path
import logging
import numpy
import paramio
import sherpa.astro.ui as sherpa
import sherpa.astro.xspec as xspec


class SpecRatio(object):

    def __init__(self, pfname):
        #
        # Define absorption, powerlaw, and blackbody models
        #--------------------------------------------------
        self.p1   = xspec.XSpowerlaw("p1")
        self.abs1 = xspec.XSphabs("abs1")
        self.b1   = xspec.XSbbody("b1")

        self.d_model_desc = { "p1" : "power-law",
                              "abs1" : "absorption H column",
                              "b1" : "blackbody"
                              }

        self.d_model_names = { "p1.PhoIndex" : "index",
                               "p1.norm" : "alpha",
                               "abs1.nH" : "nH",
                               "b1.norm" : "norm",
                               "b1.kT"   : "kT"
                               }

        self.d_model_ids = { "p1" : "pl",
                             "b1" : "bb",
                             "abs1" : "abs"
                             }

        self.pfname         = pfname
        #---------------------------------------------------
        # parameter struct
        #---------------------------------------------------
        self.params = {
            "spectrum" : None,  # Name of spectrum file
            "outfile" : None,   # Output w/ source characteristics, flux etc.
            "e_band" : None,    # energy range (keV)
            "indx" : None,      # photon index of power law model
            "norm" : None,      # amplitude of power law model
            "nh" : None,        # column density of H absorption model
            "nh_gal" : None,    # column density of H absorption model for fit X
            "kt" : None,        # temperature of blackbody model (keV)
            "verbose" : None
            }

    #
    # Class Functions 
    #
    def quit(self):
        """ a graceful exit from sherpa """
        sys.exit()

    def get_verbose(self):
        return self.params["verbose"]

    def set_verbose(self, verb):
        log = logging.getLogger("sherpa")

        # Conversion of CIAO verbosity levels to Python logging levels
        convert = { 0 : 50, 1 : 40, 2 : 30, 3 : 20, 4 : 10, 5 : 0 }
        log.setLevel(convert[verb])
        del log

    def get_parameters(self):
        """ Get parameters from param file """

        params = {}
        try:
            pfile = paramio.paramopen(self.pfname, "rw", sys.argv)
        except:
            print "Error: could not open parameter file '%s'" % self.pfname
            raise

        params["specfile"]= paramio.pgetstr(pfile,"specfile")
        params["outfile"] = paramio.pgetstr(pfile,"outfile")
        self.e_band       = paramio.pgetstr(pfile,"e_band")
        self.d_indx       = paramio.pgetd(pfile,"indx")
        self.d_norm       = paramio.pgetd(pfile,"norm")
        # Convert nH from 10^20/cm^2 to 10^22/cm^2
        # Because measured value coming from colden is in units
        # of 10^20/cm^2, but XSPEC absorption needs units of 10^22/cm^2
        self.d_nh         = paramio.pgetd(pfile,"nH") / 100.
        self.d_nh_gal     = paramio.pgetd(pfile,"nH_gal") / 100.
        self.d_kt         = paramio.pgetd(pfile,"kT")
        params["verbose"] = paramio.pgeti(pfile,"verbose")
        paramio.paramclose(pfile)

        # print "%g\t%g\n" % (self.d_nh, self.d_nh_gal)

        if self.e_band == "" or "," not in self.e_band:
            raise IOError("could not understand e_band input of %s" % self.e_band)

        self.e_band  = [float(val) for val in self.e_band.split(",")]
        self.e_band  = numpy.asarray(self.e_band, float)

        # Remove outfile if it exists
        if os.path.isfile(params["outfile"]):
            os.remove(params["outfile"])

        # Open outfile
        try:
            self.d_fout = file(params["outfile"], "a+")
        except:
            if params["outfile"] != "":
                os.system("touch " + params["outfile"])
            print "Error: could not open output file '%s'" % params["outfile"]
            raise

        return params


    def init_powlaw_model(self):
        """ initialize powerlaw model """
        src_model = self.abs1 * self.p1

        self.abs1.nh     = self.d_nh_gal
        self.p1.phoindex = self.d_indx
        self.p1.norm     = self.d_norm
        
        sherpa.freeze(src_model)

        return src_model


    def init_bbody_model(self):
        """ initialize bbody model """ 
        src_model = self.abs1 * self.b1

        self.abs1.nH = self.d_nh
        self.b1.kT   = self.d_kt
        self.b1.norm = 1.0

        sherpa.freeze(src_model)

        return src_model


    def write_params(self, params):
        """ write out to file the parameter information """
        model = sherpa.get_source()
        lhs = (type(model.lhs).__name__.lower() + '[' +
               self.d_model_ids[model.lhs.name] + ']')
        rhs = (type(model.rhs).__name__.lower() + '[' +
               self.d_model_ids[model.rhs.name] + ']')

        self.d_fout.write("# Source Model: %s * %s\n" %
                          (lhs, rhs))

        for par in params:
            id = getattr(par, "modelname", None)
            if id is None:
                id  = par.name
            format = "%s,r,h,%s,,,\"%s %s, %s \"\n"
            arg    = "INDEF"
            if not numpy.isnan(par.val):
                format = "%s,r,h,%5.1f,,,\"%s %s, %s \"\n"
                arg    = par.val
                if par.name == "nH":
                    arg = arg*100.

            if par.name == "nH":
                par.units = "10^20 atoms / cm^2"

            self.d_fout.write(format % (self.d_model_names[par.fullname], arg, 
                                        self.d_model_desc[par.modelname],
                                        par.name, par.units))


    def write_ratio(self, ratio):
        """ write out the f2c ratio to file """
        model = sherpa.get_source()
        name = model.rhs.name
        self.d_fout.write("f2c_%s,r,h,%.2e,,,\"f/c ratio (e_band = %.1f-%.1f keV), %s\"\n" %
                          (self.d_model_ids[name], ratio["val"], self.e_band[0],
                           self.e_band[1], ratio["units"]))


    def get_name(self, dataset):
        """ get filename from dataset """
        if dataset is not None:
            return dataset.name
        
        return "none"


    def write_name(self, descr, arg):
        """ write arguments to file """
        self.d_fout.write("# %s:\n#  %s\n" % (descr, arg))


    def write_data(self):
        """ write out filenames from data set to file """
        
        sname = self.get_name(sherpa.get_data())
        aname = self.get_name(sherpa.get_arf())
        rname = self.get_name(sherpa.get_rmf())
        
        self.d_fout.write("#\n")
        self.write_name("Source Data",   os.path.abspath(sname))
        self.write_name("Instrument ARF",os.path.abspath(aname))
        self.write_name("Instrument RMF",os.path.abspath(rname))
        self.d_fout.write("#\n")


    def calc_f2c_ratio(self):
        """ calculate flux/counts ratio """

        ratio = {}

        # calulate flux for the e-range on the dataset
        flux = sherpa.calc_energy_flux(lo=self.e_band[0],
                                       hi=self.e_band[1])
        
        # update the value with model counts
        counts = sherpa.calc_model_sum(lo=self.e_band[0],
                                       hi=self.e_band[1]);
        
        if counts == 0.0:
            raise("Error: counts of source model %s is zero" %
                  sherpa.get_source().name)
        
        ratio["val"] = flux/counts
        ratio["units"] = "ergs/count";
        
        return ratio;


    def do_ratio(self, model):
        """ run flux/counts ratio thread for current model """
        
        ratio = self.calc_f2c_ratio()
        self.write_params(model.pars)
        self.write_ratio(ratio)


    def run_lev3_ratio(self):
        """ run flux/counts ratio threads for bbody and powerlaw """

        params = self.get_parameters()

        self.set_verbose(params["verbose"])

        sherpa.load_pha(params["specfile"])

        self.write_data()

        model = self.init_powlaw_model()
        sherpa.set_source(model)
        self.do_ratio(model)

        model = self.init_bbody_model()
        sherpa.set_source(model)
        self.do_ratio(self.b1)

        self.d_fout.close()


#------------------------------------
# main routine
#        
#-----------------------------------
# get parameter file name, the root of  the executable
#

if __name__ == "__main__":

    # retrieve the script's name
    name = __file__.lstrip("./").rstrip(".py") + ".par"
 
    spec = SpecRatio("lev3_f2cratio.par")
    spec.run_lev3_ratio()
    spec.quit()

# --- end of the main routine call ----
