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


import os
import os.path
import sys
import numpy
import logging
import paramio
import pycrates
from exceptions import ValueError
from sherpa.stats import Stat
from sherpa.data import BaseData
from sherpa.utils import NoNewAttributesAfterInit
from sherpa.astro.ui import *

# Quiet numpy chatter
_ = numpy.seterr(all='ignore')

def _check_if_nan(arg, format="%.4e"):
    val = "INDEF"
    if not (numpy.isnan(arg) or numpy.isinf(arg)):
        val = format % arg
    return val

class MHOData(DataIMG):
    def eval_model_to_fit(self, modelfunc):
        return [None]
    
    def to_fit(self, staterrfunc=None):
        return ([1],[1],[1])


class MHOStruct(BaseData):
    def __init__(self, x, y, x0, y0, x_0, x_1, i, image):
        BaseData.__init__(self)


class MHOStat(Stat):
    def __init__(self, calc_stat, name='mho_stat'):
        self.calc_stat = calc_stat;
        Stat.__init__(self, name)


#  Author: John C. Houck <houck@space.mit.edu>
#  Written: July 2007
#  Python port: Feb 2008    Brian Refsdal <brefsdal@cfa.harvard.edu>
#
# This algorithm is intended to generate an elliptical source region
# given an approximate source position (logical X,Y) and a crude estimate
# for the source size, such as may be produced by a source detection
# algorithm.
#
# Interface:
#
#  mho = MHO(x0, y0, crude_src_size, img_size, bin_factor)
#
#      x0, y0 = Approximate sky X, Y coordinates for a detected source.
#
#      crude_src_size = Crude estimate of the source size [pixels].
#               For best results, the estimate should be within a
#               factor ~3 of the actual source size.  For point sources
#               and moderately extended sources, use the approximate
#               PSF sigma at (x0, y0) [e.g. radius enclosing 40% of
#               the flux].
#
#      img_size = Size [pixels] of the region centered on (x0,y0) to be 
#                 examined by the algorithm.  If img_size=0, the algorithm
#                 will choose a value.
#
#      bin_factor = Integer factor by which the extracted image should be
#                   binned.  If bin_factor=0, the algorithm will choose a value.
#
#
#  x0, y0, pars[], errs[] = mho.mho_find_source_extent(events, shape)
#
#       Source region parameters
#
#                x0, y0  = sky X,Y coordinates of the region center [pixels]
#                pars[]  = (sigma_1, sigma_2, phi)
#                sigma_* = Gaussian sigma values [pixels]
#                phi     = Ellipse position angle [radians]
#                          The +X axis is phi=0.
#                          The +Y axis is phi=PI/2.
#
#      evtfile = lev3 events file of sky X, Y coordinates of detected events
#
#      shape   = estimated shape of source [True  : 2-D gaussian,
#                                           False : uniform circular disk]
#
#
#  Note:  'pixel' means ACIS pixel = 0.492 arcsec.
#
# Mho_Max_Size_Factor sets the largest size scale and the maximum axis
# ratio for source regions.  If the value is "too large", things
# will run a lot slower (big images) and confusion problems
# in crowded fields will probably get worse.
Mho_Max_Size_Factor = 3.0

def shift(a, size):
    a_size = len(a)
    if size > 0:
        return numpy.concatenate((a[-(a_size-size):], a[:size]))
    elif size < 0:
        return numpy.concatenate((a[size:], a[:a_size+size]))
    
    return a


def extract_image(tx, ty, x0, y0, xsize, ysize, bin_factor):
    hwx = xsize/2; hwy = ysize/2

    i = (numpy.abs(tx-x0) < hwx) & (numpy.abs(ty-y0) < hwy)

    nx = int(numpy.round(xsize/bin_factor))
    ny = int(numpy.round(ysize/bin_factor))

    # odd dimensions work better, assuming that the
    # image is centered on the source of interest

    if (nx/2)*2 == nx: nx += 1
    if (ny/2)*2 == ny: ny += 1

    x = numpy.linspace(x0-hwx, x0+hwx, nx)
    y = numpy.linspace(y0-hwy, y0+hwy, ny)
    image = histogram2d(ty[i], tx[i], y, x)

    f = MHOStruct( x, y, x0, y0, None, None, None, image)

    f.x -= x0; f.y -= y0
    # grid will give pixel center coordinates
    # except in over-flow bin
    dx = f.x[1] - f.x[0]
    dy = f.y[1] - f.y[0]
    f.x += 0.5*dx; f.y += 0.5*dy

    (f.x_0, f.x_1) = numpy.meshgrid(f.x, f.y)
    f.i = f.image.nonzero()

    return f


def sign_change(y):
    sy = numpy.sign(y)
    sy_i = numpy.arange(sy.size)
    i = sy_i[(sy != shift(sy,1))]
    if len(i) == 0:
        return i
    n = len(y)
    return i[(i < (n-2))]


def deriv(y):
    return (y - shift(y,1))


def elliptical_mh(x, y, pars):
    phi = pars[2]

    lhs = x*numpy.cos(phi) + y*numpy.sin(phi)
    rhs = -x*numpy.sin(phi) + y*numpy.cos(phi)
    r2  = (numpy.square(numpy.divide(lhs,pars[0])) +
           numpy.square(numpy.divide(rhs,pars[1])))

    return (2. - r2) * numpy.exp(-r2/2.)


def randomize(f, val):
    if val != 0:
        return (val * (1.0 + f * numpy.random.standard_normal()))
    return (f * numpy.random.standard_normal())


def get_parameters(filename):
    """ Get parameters from param file """

    params = {}
    try:
        pfile = paramio.paramopen(filename, "rw", sys.argv)
    except:
        print "Error: could not open parameter file '%s'" % filename
        raise

    params["evtfile"]  = paramio.pgetstr(pfile,"evtfile")
    out_filename       = paramio.pgetstr(pfile,"outfile")
    params["x0"]       = paramio.pgetd(  pfile,"x0")
    params["y0"]       = paramio.pgetd(  pfile,"y0")
    params["srcsize"]  = paramio.pgetd(  pfile,"srcsize")
    params["imgsize"]  = paramio.pgetd(  pfile,"imgsize")
    params["binfactor"]= paramio.pgeti(  pfile,"binfactor")
    src_shape          = paramio.pgetstr(pfile,"shape")
    params["mincounts"]= paramio.pgeti(  pfile,"mincounts")
    params["minthresh"]= paramio.pgeti(  pfile,"mincounts.p_min")
    params["sigmafactor"]= paramio.pgetd(pfile,"sigmafactor")
    params["verbose"]  = paramio.pgeti(  pfile,"verbose")
    paramio.paramclose(pfile)

    params["shape"]=False
    if src_shape.startswith("gaussian"):
        params["shape"]=True

    # Remove outfile if it exists
    if os.path.isfile(out_filename):
        os.remove(out_filename)

    try:
        params["outfp"] = file(out_filename, "a+")
    except:
        if out_filename != "":
            os.system("touch " + out_filename)
        print "Error: could not open output file '%s'" % out_filename
        raise

    return params


def write_parameters(fp, x0, y0, sigma1, sigma2, phi, sigma1_err, sigma2_err):

    phi *= 180/numpy.pi # switch over to degrees

    # The L3 standard, the position angle zero point
    # should lie on the sky +Y axis (North DEC).
    if (0 <= phi) and (phi < 90):
        phi += 90
    elif (90 <= phi) and (phi < 180):
        phi -= 90

    fp.write("mho_x0,r,h,%s,,,\"%s\"\n" %
             (_check_if_nan(x0), "center about x0 [pixels]"))

    fp.write("mho_y0,r,h,%s,,,\"%s\"\n" %
             (_check_if_nan(y0), "center about y0 [pixels]"))

    fp.write("mho_sigma_1,r,h,%s,,,\"%s\"\n" %
             (_check_if_nan(sigma1),
              "first gaussian sigma [pixels]"))

    fp.write("mho_sigma_1_err,r,h,%s,,,\"%s\"\n" %
             (_check_if_nan(sigma1_err),
              "first gaussian sigma uncertainty [pixels]"))

    fp.write("mho_sigma_2,r,h,%s,,,\"%s\"\n" %
             (_check_if_nan(sigma2), "second gaussian sigma [pixels]"))

    fp.write("mho_sigma_2_err,r,h,%s,,,\"%s\"\n" %
             (_check_if_nan(sigma2_err),
              "second gaussian sigma uncertainty [pixels]"))

    fp.write("mho_phi,r,h,%s,,,\"%s\"\n" %
             (_check_if_nan(phi),"ellipse position angle [degrees]"))

    fp.write("mode,s,h,\"ql\",,,\n")
    fp.close()


# Conversion of CIAO verbosity levels to Python logging levels
verbosity_conversion = { 0 : 50, 1 : 40, 2 : 30, 3 : 20, 4 : 10, 5 : 0 }

def get_verbose():
    log = logging.getLogger("sherpa")
    level = log.level
    del log
    return level


def set_verbose(verb):
    log = logging.getLogger("sherpa")
    log.setLevel(verbosity_conversion[verb])
    del log


class MHO(NoNewAttributesAfterInit):

    def __init__(self, x0, y0, crude_src_size, img_size=0, bin_factor=0,
                 min_counts=15, min_thresh=4):
        self.f = None
        self.x0 = x0
        self.y0 = y0
        self.crude_src_size = crude_src_size
        self.img_size = img_size
        self.bin_factor = bin_factor
        self.min_counts = min_counts
        self.min_thresh = min_thresh


    def correlation_integral(self, x, y, pars):
        if len(self.f.i) == 0:
            return 0.0

        dx = self.f.x[1] - self.f.x[0]
        dy = self.f.y[1] - self.f.y[0]

        w = elliptical_mh(self.f.x_0[self.f.i]-x,
                          self.f.x_1[self.f.i]-y, pars)

        return (numpy.sum(w * self.f.image[self.f.i]) * dx * dy)


    def axes_constraint_model(self, *args):
        return self.axes_constraint_stat(*args)[0]


    def axes_constraint_stat(self, *args):

        p = [par.val for par in get_model().pars]

        # If the source's physical location is known, the correlation
        # maximum will occur in that pixel, so we don't need to compute
        # the other pixels.
        dims = self.f.image.shape
        y = self.f.y[dims[0]/2]
        x = self.f.x[dims[1]/2]

        W_max = self.correlation_integral (x, y, p)

        return [numpy.divide(-W_max , numpy.sqrt(p[0]*p[1])),
                numpy.array([])]


    def randomize_near_value(self, f, list):

        if list is None:
            list = get_model().pars

        for par in list:
            if par.frozen:
                continue

            val = randomize(f, par.val)
            while(val < par.min or par.max <= val):
                val = randomize(f, par.val)
                
            par.val = val


    def compute_stat(self, a):
        get_model().thawedpars = [a, a, 0]
        return -calc_stat()


    def correlation_scale(self):

        _as = self.f.x[(self.f.x > 1)]
        w = numpy.asarray([self.compute_stat(i) for i in _as], dtype=float)

        w1 = deriv(w)
        i1 = sign_change(w1)

        if len(i1) == 0:
            return _as[0]

        w2 = deriv(w1)
        i2 = sign_change(w2)

        ilist = numpy.concatenate([i1, i2])
        ilist.sort()
        i = ilist[0]

        if w2[i] > 0:
            if w[0] > w[i]:
                return _as[0]
            else:
                i = ilist[1]

        return _as[i]


    def mho_optimize_axes(self):

        a_max_factor = Mho_Max_Size_Factor

        if not (self.f.image > 0).any():
            return [-1.0, -1.0, 0.0]

        load_user_model(self.axes_constraint_model, "constraint")
        add_user_pars("constraint", ["a1", "a2", "phi"])
        set_source("constraint")

        # overload the fitting callback function with our own
        #get_stat().calc_stat = axes_constraint_stat
        set_stat(MHOStat(self.axes_constraint_stat))

        # correlation_scale almost always gives the best
        # initial estimate for a.  Is there a situation
        # where the input a_guess is consistently better?

        a_scale = None #a_guess;
        if (a_scale is None) or (a_scale < self.crude_src_size):
            a_scale = self.correlation_scale()
            if a_scale < self.crude_src_size:
                a_scale = self.crude_src_size

        # Prefer a1 >= a2.
        # 0.6 is a common Chandra PSF axis ratio.
        #
        # a_max is set conservatively low to minimize trouble in
        # crowded fields.  For that reason, there's some motivation to
        # set a_min generously low -- the PSF size provides a
        # natural lower bound to filter on and, if the correlation
        # scale is relatively large, the generous lower bound provides
        # wiggle room to fit highly elliptical source regions.
        phi_start= 0.0
        a_min    = a_scale * 0.2
        a2_start = a_scale * 0.8
        a1_start = a_scale * 1.2
        a_max    = a_scale * a_max_factor

        if a_max <= a1_start:
            raise ValueError("Mho_Max_Size_Factor=%s is too small" %
                                   Mho_Max_Size_Factor)

        set_par(get_model().a1, a1_start, a_min, a_max, 0)
        set_par(get_model().a2, a2_start, a_min, a_max, 0)
        set_par(get_model().phi, phi_start, -1.5*numpy.pi, 1.5*numpy.pi, 0)

        set_method("neldermead")
        set_method_opt("ftol", 1.e-4)

        num_loops = 0
        while(True):
            if get_verbose() < 30:
                print get_model()

            # adjust phi
            freeze(get_model())
            thaw(get_model().phi)
            fit()

            if get_verbose() < 30:
                print get_fit_results().format()

            # adjust all params
            thaw(get_model())
            fit()
            
            if get_verbose() < 30:
                print get_fit_results().format()

            # avoid quitting when axes are at their range limits
            pars = get_model().thawedpars
            pars = numpy.asarray(pars[:2])
            pegged = (numpy.divide(numpy.abs(1 - pars.max()),a_max) < 0.01
                      or numpy.divide(numpy.abs(1 - pars.min()),a_min) < 0.01)
            if get_fit_results().succeeded and pegged == 0:
                break

            num_loops += 1
            if num_loops >= 5:
                break

            # re-set the parameters only if we're going to try
            # another fit
            get_model().thawedpars = [a1_start, a2_start,
                                      numpy.pi*numpy.random.uniform()]

            self.randomize_near_value(0.05, None)

        # Always report a1 >= a2
        pars = numpy.asarray(get_model().thawedpars)
        if pars[0] < pars[1]:
            pars = pars[[1,0,2]]
            pars[2] += numpy.pi/2.

        # Always report 0 <= phi <= PI
        pars[2] = numpy.mod(pars[2],numpy.pi)
        if pars[2] < 0:
            pars[2] += numpy.pi

        get_model().thawedpars = list(pars)
        return pars


    def position_size_stat(self, *args):
        p = [par.val for par in get_model().pars]
        pars = numpy.array([p[2], p[2], 0])
        return [numpy.divide(-self.correlation_integral(p[0], p[1], pars),p[2]),
                numpy.array([])]


    def position_size_model(self, *args):
        return self.position_size_stat(*args)[0]


    def improve_source_position(self):
        
        if not (self.f.image > 0).any():
            return
                
        nx = len(self.f.x)
        ny = len(self.f.y)
        ox = nx/4
        oy = ny/4


        a_min = 0.5 * self.crude_src_size
        a_max = 0.5 * numpy.array([self.f.x[nx-ox-1],
                                   self.f.y[ny-oy-1]]).min()
        a = numpy.array([self.crude_src_size,
                         numpy.sqrt(a_min * a_max)]).min()

        if a_min > a_max:
            #raise ValueError("a_min[%g] is greater than a_max[%g]" % (a_min, a_max))
            a_max,a_min = a_min,a_max

        # Look for a better source position
        # near the middle of the sub-image
        load_user_model(self.position_size_model, "constraint")
        add_user_pars("constraint",
                      ["x", "y", "a"],
                      [0.0,0.0,a],
                      [self.f.x[ox],self.f.y[oy],a_min],
                      [self.f.x[nx-ox-1],self.f.y[ny-oy-1],a_max],
                      parfrozen=[False,False,False])
        set_source("constraint")

        # overload the calc_stat function with our own
        set_stat(MHOStat(self.position_size_stat))

        initial = calc_stat()

        set_method("neldermead")
        set_method_opt("ftol", 1.e-4)

        if get_verbose() < 30:
            print get_model()
        
        fit()

        if get_verbose() < 30:
            print get_fit_results().format()

        final = calc_stat()
        if final < initial:
            pars = get_model().thawedpars
            self.x0 = self.f.x0 + pars[0]
            self.y0 = self.f.y0 + pars[1]
            self.a  = pars[2]
            return True
        return False


    def calc_errors(self, counts, pars):
        """ calculate uncertainties  """
        if counts == 0:
            return numpy.array([numpy.inf]*2)

        sig_errs = 0.3 * numpy.sqrt(numpy.divide(15.,counts)) * pars[:2]
        return sig_errs


    def counts_error(self, evtfile):
        print ("Error: less than %s counts detected in %s" %
               (self.min_thresh, evtfile))
        return [self.x0, self.y0,
                numpy.array([numpy.nan]*3),
                numpy.array([numpy.nan]*2)]


    def mho_find_source_extent(self, evtfile, source_shape=True):
        adjustable_bin_factor = False

        tbl = pycrates.TABLECrate(evtfile)
        if tbl.get_column("X") is None:
            return self.counts_error(evtfile)

        X = pycrates.get_colvals(tbl, "X")
        Y = pycrates.get_colvals(tbl, "Y")

        if X.size < self.min_thresh:
            return self.counts_error(evtfile)

        if self.img_size <= 0:
            # img_size should be _much_ larger than the max allowed
            # elliptical region size
            max_major_axis = 2 * (Mho_Max_Size_Factor * self.crude_src_size)
            self.img_size = 10 * max_major_axis

        if self.bin_factor <= 0:
            adjustable_bin_factor = True
            if self.img_size < 128:
                self.bin_factor = 1
            else:
                self.bin_factor = int(numpy.round(self.img_size/128.))

        self.f = extract_image(X,Y,self.x0, self.y0, self.img_size,
                               self.img_size, self.bin_factor)

        set_data(MHOData("",None,None,None))

        if adjustable_bin_factor and self.bin_factor > 1:
            # If the source region looks bright, don't use coarse binning.
            while True:
                dims = self.f.image.shape
                ny = dims[0]
                nx = dims[1]
                slice = self.f.image[ny/2-ny/8:ny/2+ny/8+1,
                                     nx/2-nx/8:nx/2+nx/8+1]
                if slice.max() < 100 or self.bin_factor == 1:
                    break
                self.bin_factor /= 2
                self.f = extract_image(X,Y,self.x0, self.y0, self.img_size,
                                       self.img_size, self.bin_factor)

        if self.improve_source_position():
            self.f = extract_image(X,Y,self.x0, self.y0, self.img_size,
                                   self.img_size, self.bin_factor)

        pars = self.mho_optimize_axes()

        if source_shape:
            # source = 2-D Gaussian => sigma_i = a_i / sqrt(3)
            pars[[0,1]] /= numpy.sqrt(3)

        else:
            # source = uniform circular disk =>  R = a * sqrt(2)
            pars[[0,1]] *= numpy.sqrt(2)

        # determine number of counts inside source ellipse
        reg = "[sky=ellipse(%s,%s,%s,%s,%s)]" % (self.x0, self.y0,
                                                 pars[0], pars[1], pars[2])

        tbl = pycrates.TABLECrate(evtfile + reg)
        if tbl.get_column("X") is None:
            counts = 0
        else:
            counts = tbl.get_column("X").get_values().size

        if counts < self.min_counts:
            print ("Error: less than %s counts detected in " % self.min_counts +
                   "detected in ellipse(%s,%s,%s,%s,%s)" %
                   ( self.x0, self.y0, pars[0], pars[1], pars[2]))
            return [self.x0, self.y0, pars, [numpy.nan]*2]

        errs = self.calc_errors(counts, pars)

        return [self.x0, self.y0, pars, errs]


if __name__ == "__main__":
    params = get_parameters("lev3_mho.par")

    set_verbose(params["verbose"])

    # modify the guess by wavdetect to desired sigma
    params["srcsize"] *= params["sigmafactor"]

    mho = MHO(params["x0"], params["y0"], params["srcsize"], params["imgsize"],
              params["binfactor"], params["mincounts"], params["minthresh"])

    x0, y0, pars, errs = mho.mho_find_source_extent(params["evtfile"],
                                                    params["shape"])
    write_parameters(params["outfp"], x0, y0, pars[0], pars[1], pars[2],
                     errs[0], errs[1])
