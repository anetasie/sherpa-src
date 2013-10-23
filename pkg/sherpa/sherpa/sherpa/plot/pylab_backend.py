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

import numpy
import pylab
from sherpa.utils import get_keyword_defaults
from sherpa.utils.err import NotImplementedErr


__all__ = ('clear_window','point','plot','histo','contour','set_subplot','init',
           'get_split_plot_defaults', 'get_plot_defaults', 'begin', 'end',
           'get_data_plot_defaults', 'get_model_plot_defaults', 'exceptions',
           'get_fit_plot_defaults', 'get_resid_plot_defaults',
           'get_ratio_plot_defaults', 'get_contour_defaults',
           'get_data_contour_defaults', 'get_model_contour_defaults',
           'get_fit_contour_defaults', 'get_resid_contour_defaults',
           'get_ratio_contour_defaults','get_confid_plot_defaults',
           'get_confid_contour_defaults', 'set_window_redraw', 'set_jointplot',
           'get_model_histo_defaults', 'get_histo_defaults')

def init():
    pass

def begin():
    pass

def end():
    set_window_redraw(True)
    if pylab.isinteractive():
        pylab.draw()

def exceptions():
    pass

def _choose(test, iftrue, iffalse=None):
    if test:
        return iftrue
    return iffalse


_errorbar_defaults = get_keyword_defaults(pylab.Axes.errorbar)


def clear_window():
    pylab.clf()

def set_window_redraw(redraw):
    if redraw:
        pylab.draw()


def point(x, y, overplot=True, clearwindow=False,
          symbol=None,
          color=None):

    if overplot:
        axes = pylab.gca()
    else:
        if clearwindow:
            clear_window()
        axes = pylab.gca()
        
    if color is None:
        str = '%s'%(symbol)
    else:
        str = '%s%s'%(color,symbol)

    point = axes.plot(numpy.array([x]),numpy.array([y]),str)[0]
    
    #pylab.draw()


def histo(xlo, xhi, y, yerr=None, title=None, xlabel=None, ylabel=None,
          overplot=False, clearwindow=True,
          yerrorbars=False,
          ecolor=_errorbar_defaults['ecolor'],
          capsize=_errorbar_defaults['capsize'],
          barsabove=_errorbar_defaults['barsabove'],
          xlog=False,
          ylog=False,
          linestyle='steps',
          linecolor=None,
          color=None,
          marker='None',
          markerfacecolor=None,
          markersize=None):

    plot(xlo, y, yerr, None, title, xlabel, ylabel, overplot, clearwindow,
         False, yerrorbars, ecolor, capsize, barsabove, xlog, ylog, linestyle,
         linecolor, color, marker, markerfacecolor, markersize, False, False)


def plot(x, y, yerr=None, xerr=None, title=None, xlabel=None, ylabel=None,
         overplot=False, clearwindow=True,
         xerrorbars=False,
         yerrorbars=False,
         ecolor=_errorbar_defaults['ecolor'],
         capsize=_errorbar_defaults['capsize'],
         barsabove=_errorbar_defaults['barsabove'],
         xlog=False,
         ylog=False,
         linestyle='steps',
         linecolor=None,
         color=None,
         marker='None',
         markerfacecolor=None,
         markersize=None,
         xaxis=False,
         ratioline=False):

    if overplot:
        axes = pylab.gca()
    else:
        if clearwindow:
            clear_window()
        axes = pylab.gca()

        xscale = _choose(xlog, 'log', 'linear')
        yscale = _choose(ylog, 'log', 'linear')
        axes.set_xscale(xscale)
        axes.set_yscale(yscale)

        if title:
            axes.set_title(title)
        if xlabel:
            axes.set_xlabel(xlabel)
        if ylabel:
            axes.set_ylabel(ylabel)

    # Even if we're doing an error bar plot, we do a normal plot first so
    # that we can take advantage of the default color cycling
    line = axes.plot(x, y)[0]

    if xerrorbars or yerrorbars:
        line.set_visible(False)
        if color is None:
            color = line.get_color()
        if markerfacecolor is None:
            markerfacecolor = color
        if xerr is not None:
            xerr = xerr / 2.
        xerr = _choose(xerrorbars, xerr)
        yerr = _choose(yerrorbars, yerr)
        line = axes.errorbar(x, y, yerr, xerr, ecolor=ecolor, capsize=capsize,
                             barsabove=barsabove, color=color,
                             markerfacecolor=markerfacecolor)[0]

    for var in ('linestyle', 'color', 'marker', 'markerfacecolor',
                'markersize'):
        val = locals()[var]
        if val is not None:
            getattr(line, 'set_' + var)(val)

    if xaxis:
        axes.axhspan(ymin=0, ymax=0, xmin=0, xmax=1)
        
    if ratioline:
        axes.axhspan(ymin=1, ymax=1, xmin=0, xmax=1)

    #pylab.draw()
   
def contour(x0, x1, y, levels=None, title=None, xlabel=None, ylabel=None,
            overcontour=False, clearwindow=True,
            linewidths=None,
            colors=None):
    
    if overcontour:
        axes = pylab.gca()
    else:
        if clearwindow:
            clear_window()
        axes = pylab.gca()

        if title:
            axes.set_title(title)
        if xlabel:
            axes.set_xlabel(xlabel)
        if ylabel:
            axes.set_ylabel(ylabel)

    x0 = numpy.unique1d(x0)
    x1 = numpy.unique1d(x1)
    y  = numpy.asarray(y)

    if x0.size * x1.size != y.size:
        raise NotImplementedErr('contourgrids')

    y = y.reshape(x1.size, x0.size)

    if levels is None:
        line = axes.contour(x0, x1, y, colors=colors, linewidths=linewidths)
    else:
        line = axes.contour(x0, x1, y, levels, colors=colors,
                            linewidths=linewidths)

    #pylab.draw()


def set_subplot(row, col, nrows, ncols, clearaxes=True,
                left=None,
                right=None,
                bottom=None,
                top=None,
                wspace=0.3,
                hspace=0.4):

    pylab.subplots_adjust(left=left, right=right, bottom=bottom, top=top,
                          wspace=wspace, hspace=hspace)

    num = row*ncols + col + 1

    # As of numpy 0.9.8, these need to be cast to int to prevent errors
    # in matplotlib
    nrows = int(nrows)
    ncols = int(ncols)
    num   = int(num)
    
    pylab.subplot(nrows, ncols, num)

    if clearaxes:
        pylab.cla()

set_jointplot = set_subplot

def get_split_plot_defaults():
    return get_keyword_defaults(set_subplot, 1)

    
def get_plot_defaults():
    return get_keyword_defaults(plot, 7)


def get_point_defaults():
    return get_keyword_defaults(point, 2)

def get_histo_defaults():
    return get_keyword_defaults(histo, 6)

def get_confid_point_defaults():
    d = get_point_defaults()
    d['symbol']='+'
    return d


def get_data_plot_defaults():
    d = get_plot_defaults()

    d['yerrorbars'] = True
    d['linestyle'] = 'None'
    d['marker'] = '.'
    return d

def get_model_histo_defaults():
    d = get_histo_defaults()
    return d

def get_model_plot_defaults():
    d = get_plot_defaults()

    d['linestyle'] = '-'
    d['marker'] = 'None'

    return d


def get_fit_plot_defaults():
    return {}


def get_resid_plot_defaults():
    d = get_data_plot_defaults()
    d['xerrorbars'] = True
    d['capsize'] = 0
    #d['marker'] = '_'
    d['xaxis'] = True
    return d


def get_ratio_plot_defaults():
    d = get_data_plot_defaults()
    d['xerrorbars'] = True
    d['capsize'] = 0
    #d['marker'] = '_'
    d['ratioline'] = True
    return d


def get_confid_plot_defaults():
    d = get_plot_defaults()

    d['linestyle'] = '-'
    d['marker'] = 'None'
    return d


def get_contour_defaults():
    return get_keyword_defaults(contour, 6)


get_data_contour_defaults = get_contour_defaults
get_model_contour_defaults = get_contour_defaults


def get_fit_contour_defaults():
    return {}


get_confid_contour_defaults = get_data_contour_defaults
get_resid_contour_defaults = get_data_contour_defaults
get_ratio_contour_defaults = get_data_contour_defaults
