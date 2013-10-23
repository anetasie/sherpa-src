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

"""
Classes for plotting, analysis of astronomical data sets
"""

from sherpa.astro.data import DataPHA
from sherpa.plot import *
from sherpa.astro.utils import compile_energy_grid, bounds_check
from sherpa.utils.err import PlotErr, IOErr
from sherpa.utils import parse_expr, dataspace1d, histogram1d, filter_bins
from numpy import iterable, set_printoptions, array2string, asarray
from itertools import izip

__all__ = ('SourcePlot','ARFPlot', 'BkgDataPlot', 'BkgModelPlot', 'BkgFitPlot',
           'BkgSourcePlot', 'BkgDelchiPlot', 'BkgResidPlot', 'BkgRatioPlot',
           'BkgChisqrPlot', 'OrderPlot', 'ModelHistogram', 'BkgModelHistogram')


class ModelHistogram(Histogram):
    "Derived class for creating 1D PHA model histogram plots"
    histo_prefs = backend.get_model_histo_defaults()

    def __init__(self):
        self.xlo = None
        self.xhi = None
        self.y  = None
        self.xlabel = None
        self.ylabel = None
        self.title = 'Model'
        Histogram.__init__(self)

    def __str__(self):
        set_printoptions(precision=4, threshold=6)

        xlo = self.xlo
        if self.xlo is not None:
            xlo = array2string(self.xlo)

        xhi = self.xhi
        if self.xhi is not None:
            xhi = array2string(self.xhi)

        y = self.y
        if self.y is not None:
            y = array2string(self.y)
        
        return (('xlo    = %s\n' +
                 'xhi    = %s\n' +
                 'y      = %s\n' +
                 'xlabel = %s\n' +
                 'ylabel = %s\n' +
                 'title  = %s\n' +
                 'histo_prefs = %s') %
                ( xlo,
                  xhi,
                  y,
                  self.xlabel,
                  self.ylabel,
                  self.title,
                  self.histo_prefs))

    def prepare(self, data, model, stat=None):

        old_filter = parse_expr(data.get_filter())
        old_group = data.grouped

        try:
            if old_group:
                data.ungroup()
                for interval in old_filter:
                    data.notice(*interval)

            (self.xlo, self.y, yerr, xerr,
             self.xlabel, self.ylabel) = data.to_plot(yfunc=model)
            self.y = self.y[1]

            if data.units != 'channel':
                elo, ehi = data._get_ebins(group=False)
                self.xlo = data.apply_filter(elo, data._min)
                self.xhi = data.apply_filter(ehi, data._max)
                if data.units == 'wavelength':
                    self.xlo = data._hc/self.xlo
                    self.xhi = data._hc/self.xhi
            else:
                self.xhi = self.xlo + 1.

        finally:
            if old_group:
                data.ignore()
                data.group()
                for interval in old_filter:
                    data.notice(*interval)


    def plot(self, overplot=False, clearwindow=True):
        Histogram.plot(self, self.xlo, self.xhi, self.y, title=self.title,
                       xlabel=self.xlabel, ylabel=self.ylabel,
                       overplot=overplot, clearwindow=clearwindow)


class SourcePlot(ModelHistogram):
    "Derived class for creating plots of the unconvolved source model"

    def __init__(self):
        self.flux = None
        self.units = None
        self.mask  = None
        self.title = 'Source'
        ModelHistogram.__init__(self)

    def __str__(self):
        return (('xlo    = %s\n' +
                 'xhi    = %s\n' +
                 'flux   = %s\n' +
                 'y      = %s\n' +
                 'xlabel = %s\n' +
                 'ylabel = %s\n' +
                 'units  = %s\n' +
                 'title  = %s\n' +
                 'histo_prefs = %s') %
                ( self.xlo,
                  self.xhi,
                  self.flux,
                  self.y,
                  self.xlabel,
                  self.ylabel,
                  self.units,
                  self.title,
                  self.histo_prefs))

    def prepare(self, data, src, lo=None, hi=None):
        # Note: src is source model before folding
        if not isinstance(data, DataPHA):
            raise IOErr('notpha', data.name)

        lo, hi = bounds_check(lo, hi)

        self.units = data.units
        self.xlabel = data.get_xlabel()
        self.title  = 'Source Model of %s' % data.name

        self.xlo, self.xhi = data.get_indep(filter=False)

        self.mask = filter_bins( (lo,), (hi,), (self.xlo,) )

        self.y = src(self.xlo, self.xhi)

        self.ylabel = 'Photons/sec/cm^2/keV'
        if data.units == "wavelength":
            (self.xlo, self.xhi) = (self.xhi, self.xlo)
            self.ylabel = 'Photons/sec/cm^2/Angstrom'

        # calculate photon flux per bin, photons/s/cm2/<quantity>
        self.flux = self.y / abs(self.xhi - self.xlo)


    def plot(self, overplot=False, clearwindow=True):
        xlo = self.xlo
        xhi = self.xhi
        flux= self.flux

        if self.mask is not None:
            xlo = self.xlo[self.mask]
            xhi = self.xhi[self.mask]
            flux= self.flux[self.mask]

        Histogram.plot(self, xlo, xhi, flux, title=self.title,
                       xlabel=self.xlabel, ylabel=self.ylabel,
                       overplot=overplot, clearwindow=clearwindow)


class ARFPlot(ModelHistogram):
    "Derived class for creating plots of ancillary response"

    def __init__(self):
        ModelHistogram.__init__(self)

    def prepare(self, arf, data=None):
        self.xlo = arf.energ_lo
        self.xhi = arf.energ_hi
        self.y = arf.specresp

        self.title = arf.name
        self.xlabel = arf.get_xlabel()
        self.ylabel = arf.get_ylabel()

        if data is not None:
            if not isinstance(data, DataPHA):
                raise PlotErr('notpha', data.name)
            if data.units == "wavelength":
                self.xlabel = 'Wavelength (Angstrom)'
                self.xlo = data._hc/self.xlo
                self.xhi = data._hc/self.xhi


class BkgDataPlot(DataPlot):
    "Derived class for creating plots of background counts"
    def __init__(self):
        DataPlot.__init__(self)


class BkgModelPlot(ModelPlot):
    "Derived class for creating plots of background model"
    def __init__(self):
        ModelPlot.__init__(self)
        self.title = 'Background Model Contribution'

class BkgFitPlot(FitPlot):
    "Derived class for creating plots of background counts with fitted model"
    def __init__(self):
        FitPlot.__init__(self)

class BkgDelchiPlot(DelchiPlot):
    "Derived class for creating background plots of 1D delchi chi ((data-model)/error)"
    def __init__(self):
        DelchiPlot.__init__(self)

class BkgResidPlot(ResidPlot):
    "Derived class for creating background plots of 1D residual (data-model)"
    def __init__(self):
        ResidPlot.__init__(self)

    def prepare(self, data, model, stat):
        ResidPlot.prepare(self, data, model, stat)
        self.title = 'Residuals of %s - Bkg Model' % data.name

class BkgRatioPlot(RatioPlot):
    "Derived class for creating background plots of 1D ratio (data:model)"
    def __init__(self):
        RatioPlot.__init__(self)

    def prepare(self, data, model, stat):
        RatioPlot.prepare(self, data, model, stat)
        self.title = 'Ratio of %s : Bkg Model' % data.name

class BkgChisqrPlot(ChisqrPlot):
    "Derived class for creating background plots of 1D chi**2 ((data-model)/error)**2"
    def __init__(self):
        ChisqrPlot.__init__(self)

class BkgSourcePlot(SourcePlot):
    "Derived class for plotting the background unconvolved source model"
    def __init__(self):
        SourcePlot.__init__(self)


class OrderPlot(ModelHistogram):
    """
    Derived class for creating plots of the convolved source model using 
    selected multiple responses
    """

    def __init__(self):
        self.orders=None
        self.colors=None
        self.use_default_colors=True
        ModelHistogram.__init__(self)

    def __str__(self):
        set_printoptions(precision=4, threshold=6)

        xlo = self.xlo
        if self.xlo is not None:
            xlo = array2string(asarray(self.xlo))

        xhi = self.xhi
        if self.xhi is not None:
            xhi = array2string(asarray(self.xhi))

        y = self.y
        if self.y is not None:
            y = array2string(asarray(self.y))

        return (('xlo    = %s\n' +
                 'xhi    = %s\n' +
                 'y      = %s\n' +
                 'xlabel = %s\n' +
                 'ylabel = %s\n' +
                 'title  = %s\n' +
                 'histo_prefs = %s') %
                ( xlo,
                  xhi,
                  y,
                  self.xlabel,
                  self.ylabel,
                  self.title,
                  self.histo_prefs))

    def prepare(self, data, model, orders=None, colors=None):
        self.orders = data.response_ids

        if orders is not None:
            if iterable(orders):
                self.orders = list(orders)
            else:
                self.orders = [orders]

        if colors is not None:
            self.use_default_colors=False
            if iterable(colors):
                self.colors = list(colors)
            else:
                self.colors = [colors]
        else:
            self.colors=[]
            top_color = '0xffffff'
            bot_color = '0x0000bf'
            num = len(self.orders)
            jump = (int(top_color, 16) - int(bot_color,16))/(num+1)
            for order in self.orders:
                self.colors.append(top_color)
                top_color = hex(int(top_color,16)-jump)

        if not self.use_default_colors and len(colors) != len(orders):
            raise PlotErr('ordercolors', len(orders), len(colors))

        old_filter = parse_expr(data.get_filter())
        old_group = data.grouped

        try:
            if old_group:
                data.ungroup()
                for interval in old_filter:
                    data.notice(*interval)

            self.xlo=[]
            self.xhi=[]
            self.y=[]
            (xlo, y, yerr,xerr,
             self.xlabel, self.ylabel) = data.to_plot(model)
            y = y[1]
            if data.units != 'channel':
                elo, ehi = data._get_ebins(group=False)
                xlo = data.apply_filter(elo, data._min)
                xhi = data.apply_filter(ehi, data._max)
                if data.units == 'wavelength':
                    xlo = data._hc/xlo
                    xhi = data._hc/xhi
            else:
                xhi = xlo + 1.

            for order in self.orders:
                self.xlo.append(xlo)
                self.xhi.append(xhi)
                if len(data.response_ids) > 2:
                    if order < 1 or order > len(model.rhs.orders):
                        raise PlotErr('notorder', order)
                    y = data.apply_filter(model.rhs.orders[order-1])
                    y = data._fix_y_units(y,True)
                    if data.exposure:
                        y = data.exposure * y
                self.y.append(y)

        finally:
            if old_group:
                data.ignore()
                data.group()
                for interval in old_filter:
                    data.notice(*interval)

        self.title = 'Model Orders %s' % str(self.orders)

        if len(self.xlo) != len(self.y):
            raise PlotErr("orderarrfail")


    def plot(self, overplot=False, clearwindow=True):
        default_color = self.histo_prefs['linecolor']
        count = 0
        for xlo, xhi, y, color in izip(self.xlo, self.xhi, self.y, self.colors):
            if count != 0:
                overplot=True
                self.histo_prefs['linecolor']=color
            Histogram.plot(self, xlo, xhi, y, title=self.title,
                           xlabel=self.xlabel, ylabel=self.ylabel,
                           overplot=overplot, clearwindow=clearwindow)
            count += 1

        self.histo_prefs['linecolor'] = default_color

class BkgModelHistogram(ModelHistogram):
    "Derived class for creating 1D background PHA model histogram plots"

    def __init__(self):
        ModelHistogram.__init__(self)


class FluxHistogram(ModelHistogram):
    "Derived class for creating 1D flux distribution plots"
    def __init__(self):
        self.modelvals=None
        ModelHistogram.__init__(self)

    def __str__(self):
        set_printoptions(precision=4, threshold=6)

        vals = self.modelvals
        if self.modelvals is not None:
            vals = array2string(asarray(self.modelvals))

        return '\n'.join(['modelvals = %s' % vals,
                          ModelHistogram.__str__(self)])


    def prepare(self, fluxes, bins):
        y = asarray(fluxes[:,0])
        self.modelvals = asarray(fluxes[:,1:])
        self.xlo, self.xhi = dataspace1d(y.min(), y.max(), numbins=bins+1)[:2]
        y = histogram1d(y, self.xlo, self.xhi)
	self.y = y/float(y.max())


class EnergyFluxHistogram(FluxHistogram):
    "Derived class for creating 1D energy flux distribution plots"

    def __init__(self):
        FluxHistogram.__init__(self)
        self.title = "Energy flux distribution"
        self.xlabel = "energy flux [ ergs cm^{-2} sec^{-1} ]"
        self.ylabel = "frequency"


class PhotonFluxHistogram(FluxHistogram):
    "Derived class for creating 1D photon flux distribution plots"

    def __init__(self):
        FluxHistogram.__init__(self)
        self.title = "Photon flux distribution"
        self.xlabel = "photon flux [ photons cm^{-2} sec^{-1} ]"
        self.ylabel = "frequency"
