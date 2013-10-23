#!/usr/bin/env python

# 
#  Copyright (C) 2008  Smithsonian Astrophysical Observatory
#
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
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


import os
import os.path
import sys
import numpy
import logging
import paramio

# Quiet numpy chatter
_ = numpy.seterr(all='ignore')

# Use smallest possible float in comparisons
tinyval = numpy.float(numpy.finfo(numpy.float32).tiny)

# Assume the shape of the intrinsic source is a circle.
# Therefore, a_rss and its error have to be corrected by this
# factor.
circle_cfactor = 1.0 / numpy.sqrt(2.0)

def _check_if_nan(arg, format="%.4e"):
    val = "INDEF"
    if not (numpy.isnan(arg) or numpy.isinf(arg)):
        val = format % arg
    return val


def _check_if_indef(arg):
    if str(arg).startswith("INDEF"):
        return arg
    return numpy.float(arg)

def read_parameter_file(filename):
    try:
        pfile = paramio.paramopen(filename, "rw", sys.argv)
    except:
        print "Error: could not open parameter file '%s'" % filename
        raise
    
    srcfile = paramio.pgetstr(pfile,"srcfile")
    psffile = paramio.pgetstr(pfile,"psffile")
    outfile = paramio.pgetstr(pfile,"outfile")
    flag    = paramio.pgetd(pfile,"extent_flag_threshold")

    return srcfile, psffile, outfile, flag


def get_parameters(filename):
    """ Get parameters from param file """

    params = {}
    try:
        pfile = paramio.paramopen(filename, "rw") #, sys.argv)
    except:
        print "Error: could not open parameter file '%s'" % filename
        raise

    sigma_1     = paramio.pgetstr(pfile,"mho_sigma_1")
    sigma_1_err = paramio.pgetstr(pfile,"mho_sigma_1_err")
    sigma_2     = paramio.pgetstr(pfile,"mho_sigma_2")
    sigma_2_err = paramio.pgetstr(pfile,"mho_sigma_2_err")

    paramio.paramclose(pfile)

    params["sigma_1"]     = _check_if_indef( sigma_1 )
    params["sigma_1_err"] = _check_if_indef( sigma_1_err )
    params["sigma_2"]     = _check_if_indef( sigma_2 )
    params["sigma_2_err"] = _check_if_indef( sigma_2_err )

    return params


def write_parameters(outfile, a_rss, a_rss_err, extended):
    """ Write parameters out to file """

    # Remove outfile if it exists
    if os.path.isfile(outfile):
        os.remove(outfile)

    try:
        fp = file(outfile, "a+")
    except:
        if outfile != "":
            os.system("touch " + outfile)
        print "Error: could not open output file '%s'" % outfile
        raise

    fp.write("a_rss,r,h,%s,,,\"%s\"\n" %
             (_check_if_nan(a_rss), "intrinsic source size [pixels]"))

    fp.write("a_rss_err,r,h,%s,,,\"%s\"\n" %
             (_check_if_nan(a_rss_err), "intrinsic source size error [pixels]"))

    fp.write("extended,i,h,%i,,,\"%s\"\n" %
             (extended, "extended source flag"))

    fp.write("mode,s,h,\"ql\",,,\n")
    fp.close()


def calc_a_rss(src, psf):
    var = (numpy.square(src["sigma_1"]) +
           numpy.square(src["sigma_2"]) -
           numpy.square(psf["sigma_1"]) -
           numpy.square(psf["sigma_2"]))
    var = numpy.array([0,var]).max()
    return numpy.sqrt(var) * circle_cfactor 


def calc_a_rss_err(a_rss, src, psf):
    # revert back to a_rss without sqrt(2)
    # the factor is replaced at the end
    a_rss /= circle_cfactor
    if a_rss < tinyval:
        a_rss = numpy.sqrt(numpy.square(psf["sigma_1"]) +
                           numpy.square(psf["sigma_2"]))
    return (1.0 / a_rss ) * circle_cfactor * numpy.sqrt(
        numpy.square(src["sigma_1"]) * numpy.square(src["sigma_1_err"]) +
        numpy.square(src["sigma_2"]) * numpy.square(src["sigma_2_err"]) +
        numpy.square(psf["sigma_1"]) * numpy.square(psf["sigma_1_err"]) +
        numpy.square(psf["sigma_2"]) * numpy.square(psf["sigma_2_err"]))


def calc_extended(a_rss, a_rss_err, psf, extent_flag_threshold):
    psf_size = circle_cfactor * numpy.sqrt(
        (numpy.square(psf["sigma_1"]) +
         numpy.square(psf["sigma_2"])))

    psf_size_err = circle_cfactor * numpy.sqrt(
        (numpy.square(psf["sigma_1_err"]) +
         numpy.square(psf["sigma_2_err"])))

    ##return ((a_rss - a_rss_err) > (psf_size + psf_size_err))
    ## both a_rss and a_rss_err already account for the PSF 
    ## size and its associated uncertainty
    return (a_rss > (a_rss_err * extent_flag_threshold))


def run_iss(outfile, src, psf, flag):

    # Assume worst, PSF and SRC sigmas and errs are INDEF
    a_rss = numpy.nan
    a_rss_err = numpy.nan
    extended = 0

    # Check if sizes are INDEFs -- if any are, exit immediately
    if (src["sigma_1"] == "INDEF" or
        src["sigma_2"] == "INDEF" or
        psf["sigma_1"] == "INDEF" or
        psf["sigma_2"] == "INDEF"):
        write_parameters(outfile, a_rss, a_rss_err, extended)
        return

    # Can at least calculate instrinsic source size, so do that
    # Report it even if you can't calculate an error.
    a_rss = calc_a_rss(src, psf)
    if _check_if_nan(a_rss) == "INDEF":
        write_parameters(outfile, a_rss, a_rss_err, extended)
        return

    # Don't calc errors or extended flag if INDEFS in errors
    if (src["sigma_1_err"] == "INDEF" or
        src["sigma_2_err"] == "INDEF" or
        psf["sigma_1_err"] == "INDEF" or
        psf["sigma_2_err"] == "INDEF"):
        pass
    else:
        a_rss_err = calc_a_rss_err(a_rss, src, psf)
        extended  = calc_extended(a_rss, a_rss_err, psf, flag)

    write_parameters(outfile, a_rss, a_rss_err, extended)


if __name__ == "__main__":
    srcfile, psffile, outfile, flag = read_parameter_file("lev3_iss.par")

    src = get_parameters(srcfile)
    psf = get_parameters(psffile)

    run_iss(outfile, src, psf, flag)
