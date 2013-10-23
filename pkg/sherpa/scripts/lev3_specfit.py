#!/usr/bin/env python

# 
#  Copyright (C) 2008  Smithsonian Astrophysical Observatory
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
# February, 2008
#
# calculate spectrals and corresponding fluxes of  stack models
#
# source model - 1: absorption * powlaw  
# source model - 2: absorption * blackbody
# 
# Input:  par file, incl PHA source data file
# output: par file, incl model parameter fits and fluxes
#
#----------------------------------------------------------------


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
import pychips as chips
import sherpa.astro.ui as sherpa
import sherpa.astro.xspec as xspec
from itertools import izip

def _check_if_nan(arg, format="%.4e"):
    val = "INDEF"
    if not (numpy.isnan(arg) or numpy.isinf(arg)):
        val = format % arg
    return val

class SpecFit(object):

    def __init__(self, pfname):
        #
        # Define absorption, powerlaw, and blackbody models
        #--------------------------------------------------
        self.p1   = xspec.XSpowerlaw("p1")
        self.abs1 = xspec.XSphabs("abs1")
        self.b1   = xspec.XSbbody("b1")

        self.d_model_desc = { "p1" : "power-law",
                              "abs1" : "absorption",
                              "b1" : "blackbody"
                              }

        self.d_model_names = { "PhoIndex" : "PhoIndx",
                               "norm" : "norm",
                               "nH" : "nH",
                               "kT"   : "kT"
                               }
        
        self.d_model_ids = { "p1" : "pl",
                             "abs1" : "abs",
                             "b1" : "bb"
                             }

        self.d_netcts_thresh = 150  # threshold for spectral fits 0.5-7.0 keV
        self.d_mingrpcts   = 16     # grouping number per bin, default

        self.d_nH          = 0.03   # H absorption column, default
        self.d_phoindex    = 1.7    # Photon index
        self.d_kT          = 1.0    # Temperature

        self.max_rstat     = 1000.0 # Maximum threshold reduced stat for proj
        self.max_remin     = 25     # Maximum reminimizations for proj 

        # Setting the Parameters Hard Limits
        #

        self.d_nH_min=1.e-7
        self.d_nH_max=1.e5
        self.d_PhoIndx_min=-10
        self.d_PhoIndx_max=10
        self.d_kT_min=0.0001
        self.d_kT_max=100.

        #
        # controlling parameters
        #
        self.d_bfilename    = None
        self.pfname         = pfname

        self.d_fout   = None
        self.d_eband  = None #0-1 energy band for fitting
        self.d_frange = None #2-3 energy band for flux calc
        self.flux_units = "ergs/cm**2/s"
        self.units = "energy"

        #---------------------------------------------------
        # parameter struct
        #---------------------------------------------------
        self.params = {
            "phafile" : None,       # Name of PHA file
            "bkgfile" : None,       # Name of bkg file
            "outfile" : None,       # Output source characteristics, flux etc.
            "erange_fit" : None,    # energy range (keV) for fitting
            "erange_flux" : None,   # energy range (keV) for flux integration
            "verbose" : None
            }

        sherpa.set_stat("chi2xspecvar")

    #
    # Class Functions 
    #
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

        params["phafile"]  = paramio.pgetstr(pfile,"phafile")
        params["bkgfile"]  = paramio.pgetstr(pfile,"bkgfile")
        params["outfile"]  = paramio.pgetstr(pfile,"outfile")

        if paramio.pgetstr(pfile,"nH") != 'INDEF':
            # Convert nH from 10^20/cm^2 to 10^22/cm^2
            # Because measured value coming from colden is in units
            # of 10^20/cm^2, but XSPEC absorption needs units of 10^22/cm^2
            self.d_nH      = paramio.pgetd(  pfile,"nH") / 100.

        self.d_mingrpcts   = paramio.pgeti(  pfile,"mincts")
        fit_range          = paramio.pgetstr(pfile,"e_band")
        flux_range         = paramio.pgetstr(pfile,"e_range")
        params["verbose"]  = paramio.pgeti(  pfile,"verbose")
        paramio.paramclose(pfile)

        self.d_bfilename = params["bkgfile"]
        
        if fit_range == "" or "," not in fit_range:
            raise "Error: could not understand e_band input of %s" % fit_range

        if flux_range == "" or "," not in flux_range:
            raise "Error: could not understand e_range input of %s" % flux_range

        params["erange_fit"]  = fit_range
        params["erange_flux"] = flux_range

        self.d_eband  = [float(val) for val in fit_range.split(",")]
        self.d_frange = [float(val) for val in flux_range.split(",")]

        self.d_eband  = numpy.asarray(self.d_eband, float)
        self.d_frange = numpy.asarray(self.d_frange, float)

        # Remove outfile if it exists
        #() = remove(params.outfile)
        if os.path.isfile(params["outfile"]):
            os.remove(params["outfile"])

        try:
            self.d_fout = file(params["outfile"], "a+")
        except:
            if params["outfile"] != "":
                os.system("touch " + params["outfile"])
            print "Error: could not open output file '%s'" % params["outfile"]
            raise

        return params


    def subtract_bkg(self, bkg):
        """ subtract bkg dataset """

        verb = self.get_verbose()
        srcf = sherpa.get_data().name

        if sherpa.get_bkg() is None:
            if bkg == "none":
                return

            try:
                sherpa.load_bkg(bkg)
            except:
                print (("Warning: Cannot load bkg file, '%s'.\n" % bkg) +
                       "proceeding without background data.\n")
                return

        if verb >= 3:
            print("--- bkg -> %s" % bkg)

        # Apply src grouping to bkg
        sherpa.get_bkg().grouping = sherpa.get_data().grouping
        sherpa.get_bkg().group()

        src_dep = sherpa.get_data().to_fit()[0]
        bkg_dep = sherpa.get_bkg().to_fit()[0]

        sherpa.subtract()

        # Includes subtracted bkg in calculation
        net_dep = sherpa.get_counts(filter=True)
        scts = src_dep[(src_dep > 0.0)].sum() # Only sum positive counts
        bcts = bkg_dep[(bkg_dep > 0.0)].sum()
        ncts = net_dep[(net_dep > 0.0)].sum()

        if verb >= 3:
            print("src_sum = %s" % scts)
            print("bkg_sum = %s" % bcts)
            print("net_cts = %s" % ncts)

        if ncts <= 0:
            sherpa.unsubtract()
            print ("Warning: zero or negative net counts: %s.\n" % ncts +
                   "Check if the src PHA, '%s' matches \n the bkg, '%s'.\n" %
                   (srcf, sherpa.get_bkg().name) +
                   "I will proceed without the background data.\n")


    def write_stat(self, res):
        """ write the statistic to file """

        statname   = sherpa.get_stat_name()
        dof        = res.dof
        stat       = sherpa.calc_stat()
        rstat      = numpy.divide(stat,dof)
        model      = sherpa.get_source() #sherpa.get_model()
        model_name = "(%s %s %s)" %(self.d_model_ids[model.lhs.name],
                                    "*",
                                    self.d_model_ids[model.rhs.name])
        id         = model.rhs.name

        self.d_fout.write("#------------\n") #clear and readble

        format  = "dof_%s,r,h,%s,,,\"degree of freedom for %s fit\"\n"
        dof = _check_if_nan(dof, "%d")
        self.d_fout.write( format  % (self.d_model_ids[id], dof, model_name))

        # write stat info 
        #
        format = "stat_%s,r,h,%s,,,\"%s statistic for %s fit\"\n"
        stat = _check_if_nan(stat)
        self.d_fout.write( format % (self.d_model_ids[id], stat, statname,
                                     model_name))

        format = "rstat_%s,r,h,%s,,,\"%s reduced statistic for %s fit\"\n"
        rstat = _check_if_nan(rstat)
        self.d_fout.write( format % (self.d_model_ids[id], rstat, statname,
                                     model_name))


    def write_fit(self, errors=None):
        """ write out to file the fit information """
        model = sherpa.get_source() #sherpa.get_model()

        lhs = (type(model.lhs).__name__.lower() + '[' +
               self.d_model_ids[model.lhs.name] + ']')
        rhs = (type(model.rhs).__name__.lower() + '[' +
               self.d_model_ids[model.rhs.name] + ']')

        self.d_fout.write("#\n#------------ Source Model: %s * %s\n" %
                          (lhs, rhs))

        if self.get_verbose() >= 3:
            print model

        if errors is None:
            for par in model.pars:
                id  = par.modelname

                format = "%s,r,h,%s,,,\"%s %s, %s\"\n"
                val = _check_if_nan(par.val)
                if par.val == par.hard_min or par.val == par.hard_max:
                    val = "INDEF"

                param_name = (self.d_model_names[par.name] +
                              "_" + self.d_model_ids[id])
                if par.name == "nH":
                    # convert nh from 10^22 to 10^20
                    val = _check_if_nan(par.val*100.)
                    if par.val == par.hard_min or par.val == par.hard_max:
                        val = "INDEF"
                    par.units="10^20 atoms / cm^2"
                    param_name = (self.d_model_names[par.name] +
                                  "_" + self.d_model_ids[id] + "_" +
                                  self.d_model_ids[model.rhs.name])
 
                self.d_fout.write(format % (param_name, val,
                                            self.d_model_desc[id],
                                            par.name, par.units))
                format = "%s_lo,r,h,%s,,,\"lower confidence limit on %s %s\"\n"
                self.d_fout.write(format % (param_name, "INDEF",
                                            self.d_model_desc[id],
                                            par.name))
                format = "%s_hi,r,h,%s,,,\"upper confidence limit on %s %s\"\n"
                self.d_fout.write(format % (param_name, "INDEF",
                                            self.d_model_desc[id],
                                            par.name))
            return

        for par in model.pars:
           if par.fullname in errors.parnames:
               id  = par.modelname
               idx = list(errors.parnames).index(par.fullname)

               param_name = (self.d_model_names[par.name] +
                             "_" + self.d_model_ids[id])
               if par.name == "nH":
                   param_name = (self.d_model_names[par.name] +
                                 "_" + self.d_model_ids[id] + "_" +
                                 self.d_model_ids[model.rhs.name])

               format = "%s,r,h,%s,,,\"%s %s, %s\"\n"
               val = _check_if_nan(par.val)
               if par.name == "nH":
                   # convert nh from 10^22 to 10^20
                   val = _check_if_nan(par.val*100.)
                   par.units="10^20 atoms / cm^2"
               if par.val == par.hard_min or par.val == par.hard_max:
                   val = "INDEF"

               self.d_fout.write(format % (param_name, val, 
                                           self.d_model_desc[id],
                                           par.name, par.units))

               # assume confidence limit hit hard minimum in proj
               lo = errors.parmins[idx]
               format = "%s_lo,r,h,%s,,,\"lower confidence limit on %s %s\"\n"
               arg    = "INDEF"

               if lo is not None:
                   arg = "%.4e" % (par.val - numpy.abs(lo))
                   if par.name == "nH":
                       # convert nh from 10^22 to 10^20
                       arg = "%.4e" % ((par.val - numpy.abs(lo))*100.)
               if (numpy.isnan(par.val) or numpy.isinf(par.val) or
                   par.val <= par.min):
                   arg = "INDEF"

               self.d_fout.write(format % (param_name, arg,
                                           self.d_model_desc[id],
                                           par.name))

               # assume confidence limit hit hard maximum in proj
               hi = errors.parmaxes[idx]
               format = "%s_hi,r,h,%s,,,\"upper confidence limit on %s %s\"\n"
               arg    = "INDEF"

               if hi is not None:
                   arg = "%.4e" % (par.val + numpy.abs(hi))
                   if par.name == "nH":
                       # convert nh from 10^22 to 10^20
                       arg = "%.4e" % ((par.val + numpy.abs(hi))*100.)
               if (numpy.isnan(par.val) or numpy.isinf(par.val) or
                   par.val >= par.max):
                   arg = "INDEF"

               self.d_fout.write(format % (param_name, arg,
                                           self.d_model_desc[id],
                                           par.name))
 

    def write_flux(self, flux, flux_lo, flux_hi):
        """ Write out flux and the error """

        id = sherpa.get_source().rhs.name #sherpa.get_model().rhs.name

        val = flux
        lo  = flux - numpy.abs(flux_lo-flux)
        hi  = flux + numpy.abs(flux_hi-flux)

        self.d_fout.write("#------------\n") #clear and readble
        format = "flux_%s,r,h,%s,,,\"Flux (e_band = %.1f-%.1f keV), %s\"\n"
        self.d_fout.write(format % (self.d_model_ids[id],
                                    _check_if_nan(val,"%.2e"), self.d_frange[0],
                                    self.d_frange[1], self.flux_units))

        format = "flux_%s_lo,r,h,%s,,,\"lower confidence limit on flux\"\n"
        self.d_fout.write(format % (self.d_model_ids[id], 
                                    _check_if_nan(lo,"%.2e")))

        format = "flux_%s_hi,r,h,%s,,,\"upper confidence limit on flux\"\n"
        self.d_fout.write(format % (self.d_model_ids[id], 
                                    _check_if_nan(hi, "%.2e")))


    def get_name(self, dataset):
        """ get filename from dataset """
        if dataset is not None:
            return dataset.name

        return "none"


    def write_name(self, descr, arg):
        """ write arguments to file """
        self.d_fout.write("# %-18s - %s\n" % (descr, arg))


    def write_info(self, data):
        """ Write data/bkg/grouping, etc. of the fittings """ 
        self.d_fout.write("#\n")

        sname = self.get_name(sherpa.get_data())
        aname = self.get_name(sherpa.get_arf())
        rname = self.get_name(sherpa.get_rmf())

        bname = self.get_name(sherpa.get_bkg())

        if self.get_verbose() >= 5:
            print(os.path.basename(sname))
            print(os.path.dirname(sname))
            print(os.path.abspath(sname))


        self.write_name("Source Data",   os.path.abspath(sname))
        self.write_name("Instrument ARF",os.path.abspath(aname))
        self.write_name("Instrument RMF",os.path.abspath(rname))
        self.write_name("Background",    os.path.abspath(bname))

        if sherpa.get_data().subtracted:
            self.write_name("bkg subtracted", "yes")

        if self.d_mingrpcts > 0:
            self.d_fout.write("# %-18s - %d\n" %("min grp counts",
                                                 self.d_mingrpcts))

        self.d_fout.write("# %-18s - noticed channels %s\n" %
                          ("filters", sherpa.get_data().get_noticed_expr()))

        self.d_fout.write("# %-18s - %s\n" % ("", self.filt_expr))


    def get_fluxproj(self, res):
        """ Calculate Flux and the projections (upper and lower values) """

        # Save fitted parameter values
        old_vals = sherpa.get_model().thawedpars

        # Calculate flux for the e-range on the dataset and the current
        # stack of src/bkg models
        flux = sherpa.calc_energy_flux(lo=self.d_frange[0],
                                       hi=self.d_frange[1])

        #set low components of fitted parameters of whole src/bkg stack
        for name, val, min in izip(res.parnames, res.parvals, res.parmins):
            par = eval("self." + name)
            if min is None:
                val = par.hard_min
            else:
                val = val - numpy.abs(min)
            sherpa.set_par(par, val, min=val)

        # set low components and calc the corresponding flux 
        flux_lo = sherpa.calc_energy_flux(lo=self.d_frange[0],
                                          hi=self.d_frange[1])
        
        for name, val, max in izip(res.parnames, res.parvals, res.parmaxes):
            par = eval("self." + name)
            if max is None:
                val = par.hard_max
            else:
                val = val + numpy.abs(max)
            sherpa.set_par(par, val, max=val)

        # set low components and calc the corresponding flux 
        flux_hi = sherpa.calc_energy_flux(lo=self.d_frange[0],
                                          hi=self.d_frange[1]) 

        # replace fitted parameter values in model
        sherpa.get_model().thawedpars = old_vals

        return flux, flux_lo, flux_hi


    def guess_norm(self, param):
        counts = sherpa.get_data().to_guess()[0]
        pos = 0
        max = counts.max()
        min = counts.min()

        if((max > 0.0 and min >= 0.0) or
           (max > 0.0 and min < 0.0 and (abs(min) <= max ))):
            pos = counts.argmax()
            if min != 0.0 or max != 0.0:
                param.max = max*1000
                param.min = max/1000
            param.val = max

        elif((max > 0.0 and min < 0.0 and abs(min) > max ) or
             (max == 0.0 and min < 0.0 ) or ( max < 0.0 )):
            pos = counts.argmin()
            if min != 0.0 or max != 0.0:
                param.max = min/1000
                param.min = min*1000
            param.val = min


    def init_powlaw_model(self):
        """ initialize powerlaw model """
        src_model = self.abs1 * self.p1

        self.abs1.nh = self.d_nH
        self.p1.phoindex = self.d_phoindex
        self.guess_norm(self.p1.norm)

        # param limits
        self.abs1.nh.min = self.d_nH_min
        self.abs1.nh.max = self.d_nH_max
        self.p1.phoindex.min = self.d_PhoIndx_min
        self.p1.phoindex.max = self.d_PhoIndx_max

        return src_model


    def init_bbody_model(self):
        """ initialize bbody model """ 

        # reinitialize absorb to be used again.
        self.abs1 = xspec.XSphabs("abs1")
        
        src_model = self.abs1 * self.b1
        
        self.abs1.nH = self.d_nH
        self.b1.kT = self.d_kT
        self.guess_norm(self.b1.norm)

        # Setting the param limits here
        self.abs1.nh.min = self.d_nH_min
        self.abs1.nh.max = self.d_nH_max
        self.b1.kT.min = self.d_kT_min
        self.b1.kT.max = self.d_kT_max

        return src_model


    def do_fit(self, src_model):
        """ do fitting  """

        sherpa.set_source(src_model)

        if self.get_verbose() >= 3:
            print(" ---")
            print sherpa.get_model()

        try:
            sherpa.fit()
            fit_results = sherpa.get_fit_results()
        except:
            print ("Error: fitting source model failed '%s' with %s" % 
                   (src_model, fit_results.message))
            raise
            
        if self.get_verbose() >= 3:
            print fit_results.format()

        # Get errors on best-fit parameters
        try:
            if hasattr(sherpa.get_proj(), "max_rstat"):
                sherpa.get_proj().max_rstat = self.max_rstat
            sherpa.get_proj().maxfits = self.max_remin
            sherpa.proj()
            proj_results = sherpa.get_proj_results()
        except:
            print ("Error: running projection failed")
            self.write_fit()
            self.write_stat(fit_results)
            flux = sherpa.calc_energy_flux(lo=self.d_frange[0],
                                           hi=self.d_frange[1])
            self.write_flux(flux, numpy.nan, numpy.nan)
            return

        if self.get_verbose() >= 3:
            print proj_results.format()

        # write the fits and flux outputs
        self.write_fit(proj_results)
        self.write_stat(fit_results)
        
        # Calculate flux and the proj
        flux, flux_lo, flux_hi = self.get_fluxproj(proj_results)
        self.write_flux(flux, flux_lo, flux_hi)


    def set_src_filter(self, threshold):
        """ setup filter """
        cnts = sherpa.get_counts(filter=True) # subtracted,grouped src counts
        sherpa.get_data().mask = (cnts >= threshold)


    def run_lev3_specfit(self):
        """ spectral thread """

        # Get parameters
        self.params = self.get_parameters()
        
        self.set_verbose(self.params["verbose"])

        # Read data and set up instrument models, RMF and ARF
        # Note that RMF, ARF files are always specified in
        # the header of the PHA file.  So, when PHA file is read
        # in here, the responses are also automatically set.

        sherpa.load_pha(self.params["phafile"])  # do auto bkg

        sherpa.plot_data()
        chips.log_scale()
        path    = os.path.dirname(os.path.abspath(self.params["outfile"]))
        outfile = os.path.basename(self.params["outfile"])
        if '.' in outfile:
            file_sub_strs = outfile.split('.')
            dev_null = file_sub_strs.pop()
            outfile = '.'.join(file_sub_strs)
        outfile = os.path.join(path,outfile)
        chips.set_preference("export.clobber", "true")
        chips.print_window(outfile + '_spectrum',["format","png"])

        # 
        #
        # do group,  based on grouping criteria, before subtract
        #

        dep = sherpa.get_counts(filter=True)
        totcnts = dep[(dep > 0.0)].sum()

        if self.get_verbose() >= 3:
            print "--- total SRC data counts: %s" % totcnts

        if (self.d_mingrpcts != 0 and totcnts < self.d_mingrpcts):
            print ("Error: can't proceed due to src counts, %d," % totcnts +
                   " under the minimum group# (%d)\n" %  self.d_mingrpcts)
            return

        sherpa.group_counts(self.d_mingrpcts)

        if self.get_verbose() >= 3:
            print("--- grouping counts per bin: %s" % self.d_mingrpcts)

        #------------------------------------------
        # subtract background data (if not already loaded)
        #------------------------------------------
        self.subtract_bkg(self.params["bkgfile"])

        #
        # set filter to zero if the data is less than 0.0 (negative)
        #
        self.set_src_filter(0.0)

        #
        # set filter on energy bands
        #

        self.filt_expr = "ignore %s :%s, %s:" % (self.units, self.d_eband[0],
                                                 self.d_eband[1])

        sherpa.ignore(None,self.d_eband[0])
        sherpa.ignore(self.d_eband[1],None)

        dep = sherpa.get_data().to_fit()[0]
        netcnts = dep[(dep > 0.0)].sum()

        if self.get_verbose() > 1:
            print ("Net Counts in the filter: %d" % netcnts)

        if netcnts < self.d_netcts_thresh:
            print ("Error: can't proceed due to the src counts, %d," % netcnts +
                   " under the current threshold (%d)\n" % self.d_netcts_thresh)
            return

        if self.get_verbose() > 1:
            print("--- %s" % self.filt_expr)

        self.write_info(self.params["phafile"])

        src_model = self.init_powlaw_model()
        self.do_fit(src_model)

        sherpa.plot_fit_delchi()
        chips.current_plot("plot1")
        chips.log_scale()

        chips.set_preference("export.clobber", "true")
        chips.print_window(outfile + '_powlaw',["format","png"])

        src_model = self.init_bbody_model()
        self.do_fit(src_model)

        sherpa.plot_fit_delchi()
        chips.current_plot("plot1")
        chips.log_scale()

        chips.print_window(outfile + '_bbody',["format","png"])

        #Close the output parameter file
        self.d_fout.close()

#------------------------------------
# main routine
#        
#-----------------------------------
# get parameter file name, the root of  the executable
#

if __name__ == "__main__":

    # retrieve the script's name
    #name = __file__.lstrip("./").rstrip(".py") + ".par"

    spec = SpecFit("lev3_specfit.par" ); #name)

    # FIXME: what does set_ignore_all do?
    #() = set_ignore_all()
    
    spec.run_lev3_specfit()

# --- end of the main routine call ----
