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


#!/usr/bin/env python
import numpy
import numpy.random
from sherpa.estmethods import *
from sherpa.astro.utils import calc_energy_flux
from itertools import izip
from sherpa.utils import divide_run_parallel
from sherpa.utils.err import EstErr

import logging
warning = logging.getLogger("sherpa").warning


__all__ = ['get_sample_uncorr', 'get_sample_corr',
           'sample_param_uncorr', 'sample_param_corr',
           'calc_flux', 'sample_flux']

def get_sample_uncorr(val, sigma, num=1, dist=numpy.random.normal):
    return dist(val, sigma, num)


def get_sample_corr(vals, cov, num=1, dist=numpy.random.multivariate_normal):
    return dist(vals, cov, num)


def sample_param_uncorr(fit, num=1, dist=numpy.random.normal):
    samples=[]
    oldestmethod = fit.estmethod

    fit.estmethod = Covariance()
    try:
        r = fit.est_errors()
    finally:
        fit.estmethod = oldestmethod

    thawedpars = [par for par in fit.model.pars if not par.frozen]
    for par, val, lo, hi in izip(thawedpars, r.parvals, r.parmins, r.parmaxes):
        sigma = None
        if lo is not None and hi is not None:
            sigma = numpy.abs(lo)
        else:
            warning("Covariance failed for '%s', trying Projection..." %
                    par.fullname)
            fit.estmethod = Projection()
            try:
                t = fit.est_errors(parlist = (par,))
                if t.parmins[0] is not None and t.parmaxes[0] is not None:
                    sigma = numpy.abs(t.parmins[0])
                else:
                    warning('1 sigma bounds for parameter ' +
                            par.fullname +
                            ' could not be found, using soft limit minimum')
                    sigma = numpy.abs(par.min)
            finally:
                fit.estmethod = oldestmethod
        samples.append(get_sample_uncorr(val, sigma, num, dist))
    samples = numpy.asarray(samples).transpose()
    return samples


def sample_param_corr(fit, num=1, dist=numpy.random.multivariate_normal):
    oldestmethod = fit.estmethod
    fit.estmethod = Covariance()
    try:
        r = fit.est_errors()
    finally:
        fit.estmethod = oldestmethod

    cov = r.extra_output
    if cov is None:
        raise EstErr('nocov')

    vals = fit.model.thawedpars
    samples = get_sample_corr(vals, cov, num, dist)

    return samples


# def calc_flux(fit, data, src, samples, method=calc_energy_flux,
#               lo=None, hi=None):
#     #fluxes=numpy.array([])
#     fluxes = []

#     old_model_vals  = fit.model.thawedpars
#     try:
#         # FIXME add parallel support here
#         for sample in samples:
#             fit.model.thawedpars = sample
#             flux = method(data, src, lo, hi)
#             #flux = [ method(data, src, lo, hi) ]
#             #flux.extend(sample)
#             #fluxes = numpy.append(fluxes, [flux])
#             fluxes.append(numpy.concatenate([[flux],fit.model.thawedpars]))
#     finally:
#         fit.model.thawedpars = old_model_vals

#     return numpy.asarray(fluxes)


def calc_flux(fit, data, src, samples, method=calc_energy_flux,
              lo=None, hi=None):

    def eval(slice, *args, **kwargs):
        fluxes = []
        old_model_vals  = fit.model.thawedpars
        try:
            for sample in slice:
                fit.model.thawedpars = sample
                flux = method(data, src, lo, hi)
                fluxes.append(numpy.concatenate([[flux],fit.model.thawedpars]))
        finally:
            fit.model.thawedpars = old_model_vals

        return numpy.asarray(fluxes)

    return divide_run_parallel(eval, samples)


def sample_flux(fit, data, src, method=calc_energy_flux, correlated=False, num=1,
                lo=None, hi=None):
    samples=None
    if correlated:
        samples = sample_param_corr(fit, num, numpy.random.multivariate_normal)
    else:
        samples = sample_param_uncorr(fit, num, numpy.random.normal)
    return calc_flux(fit, data, src, samples, method, lo, hi)
