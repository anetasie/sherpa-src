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
Objects and utilities used by multiple Sherpa subpackages
"""

import operator
import inspect
from itertools import izip
from types import FunctionType as function
from types import MethodType as instancemethod
import string
import numpy
import numpy.random
import numpytest
import numpy.fft
from sherpa.utils._utils import *
from sherpa.utils._psf import convolve, extract_kernel, normalize, \
    set_origin, clear_kernel_fft, pad_bounding_box

_multi=False
_ncpus=1

try:
    import multiprocessing
    _multi=True
    _ncpus = multiprocessing.cpu_count()
except:
    pass


__all__ = ('calc_ftest', 'calc_mlr', 'erf', 'igamc', 'lgam', 'sao_fcmp', 
           'NoNewAttributesAfterInit', 'convolve', 'extract_kernel',
           'SherpaTestCase', 'needs_data', 'SherpaTest',
           'bool_cast', 'export_method', 'get_keyword_names',
           'get_keyword_defaults', 'get_func_usage', 'get_num_args',
           'print_fields', 'calc_total_error', 'poisson_noise', 'rebin',
           'clear_kernel_fft', 'numpy_convolve', 'normalize', 'set_origin',
           'is_binary_file', 'guess_amplitude', 'guess_amplitude2d', 'guess_fwhm',
           'guess_reference', 'guess_position', 'guess_bounds', 'get_midpoint',
           'get_peak', 'get_valley', 'get_fwhm', 'guess_radius','histogram1d',
           'guess_amplitude_at_ref', 'get_position', 'histogram2d',
           'sum_intervals', 'pad_bounding_box', 'interpolate', 'create_expr',
           'dataspace1d', 'dataspace2d', 'parse_expr', 'gamma', 'igam',
           'incbet', 'divide_run_parallel', 'param_apply_limits',
           '_guess_ampl_scale', 'sao_arange', 'filter_bins',
           'bisection', 'demuller', 'new_muller', 'apache_muller', 'zeroin')

_guess_ampl_scale = 1.e+3

###############################################################################
#
# Types
#
###############################################################################



# Default numeric types (these match the typedefs in extension.hh)
SherpaInt = numpy.intp
SherpaUInt = numpy.uintp
SherpaFloat = numpy.float_


class NoNewAttributesAfterInit(object):
    """

    Prevents attribute deletion and setting of new attributes after
    __init__ has been called.  Derived classes must call
    NoNewAttributesAfterInit.__init__ after all other initialization.

    """

    __initialized = False  # Use name mangling

    def __init__(self):
        self.__initialized = True

    def __delattr__(self, name):
        if self.__initialized and hasattr(self, name):
            raise AttributeError(("'%s' object attribute '%s' cannot be " +
                                  "deleted") % (type(self).__name__, name))
        object.__delattr__(self, name)

    def __setattr__(self, name, val):
        if self.__initialized and (not hasattr(self, name)):
            raise AttributeError("'%s' object has no attribute '%s'" %
                                 (type(self).__name__, name))

        if self.__initialized and hasattr(self, name):
            if callable(getattr(self, name)) and not callable(val):
                raise AttributeError(("'%s' object attribute '%s' cannot be " +
                                      "replaced with a non-callable attribute")%
                                     (type(self).__name__, name))
            elif not callable(getattr(self, name)) and callable(val):
                raise AttributeError(("'%s' object attribute '%s' cannot be " +
                                      "replaced with a callable attribute") %
                                     (type(self).__name__, name))

        object.__setattr__(self, name, val)



###############################################################################
#
# Testing stuff
#
###############################################################################



class SherpaTestCase(numpytest.NumpyTestCase):
    "Base class for Sherpa unit tests"

    datadir = None

    def assertEqualWithinTol(self, first, second, tol=1e-7, msg=None):
        """

        Test first and second for floating-point equality to within
        the given tolerance.  If first and second are arrays, then
        they will be tested for equality on an element-by-element
        basis.

        """

        self.assert_(not numpy.any(sao_fcmp(first, second, tol)), msg)

    def assertNotEqualWithinTol(self, first, second, tol=1e-7, msg=None):
        """

        Test first and second for floating-point inequality to within
        the given tolerance.  If first and second are arrays, then
        they will be tested for inequality on an element-by-element
        basis.

        """

        self.assert_(numpy.all(sao_fcmp(first, second, tol)), msg)


def needs_data(meth):
    """

    Decorator for test_* methods of SherpaTestCase subclasses that
    indicates that the corresponding test requires external data
    (i.e. data not distributed with Sherpa itself).  Its effect is to
    make the test a no-op if SherpaTestCase.datadir is None.

    """

    def new_meth(self):
        if SherpaTestCase.datadir is None:
            return
        meth(self)
    new_meth.__name__ = meth.__name__
    return new_meth


class SherpaTest(numpytest.NumpyTest):
    "Sherpa test suite manager"

    def test(self, level=1, verbosity=1, datadir=None):
        old_datadir = SherpaTestCase.datadir
        SherpaTestCase.datadir = datadir

        try:
            numpytest.NumpyTest.test(self, level, verbosity)
        finally:
            SherpaTestCase.datadir = old_datadir



###############################################################################
#
# Utilities
#
###############################################################################


# at what precisions do we assume equality in energy grids?
eps = numpy.finfo(numpy.float32).eps

def filter_bins(mins, maxes, axislist):
    mask=None

    for lo,hi,axis in izip(mins,maxes,axislist):

        if (lo is None) and (hi is None):
            continue

        if lo is None:
            #axismask = axis <= hi
            axismask = (sao_fcmp(hi,axis,eps) >= 0)

        elif hi is None:
            #axismask = axis >= lo
            axismask = (sao_fcmp(lo,axis,eps) <= 0)

        else:
            #axismask = (lo <= axis) & (axis <= hi)
            axismask = (sao_fcmp(lo,axis,eps) <= 0) & (sao_fcmp(hi,axis,eps) >= 0)

        if mask is None:
            mask = axismask
        else:
            mask &= axismask

    return mask


def bool_cast(val):
    "Converts a string (true|False|on|OFF|etc...)  to a boolean value"

    if type(val) in (tuple, list, numpy.ndarray):
        return numpy.asarray([bool_cast(item) for item in val], bool)

    elif type(val) == str:
        # since built in bool() only returns false for empty strings
        vlo = val.lower()
        if vlo in ('false','off','no','0','f','n'):
            return False

        elif vlo in ('true','on','yes','1','t','y'):
            return True

        raise TypeError("unknown boolean value: '%s'" % str(val))

    else:
        # use built in bool cast
        return bool(val)


def export_method(meth, name=None, modname=None):
    """

    Given a bound instance method, return a simple function that wraps
    it.  The only difference between the interface of the original
    method and the generated function is that the latter doesn't
    include 'self' in its argument list.  This means that when the
    wrapper function is called with an incorrect number of arguments,
    the error message does not include 'self' in argument counts.  The
    only reason to generate such a wrapper is to hide from a user the
    fact that they're using an instance method rather than a simple
    function.  If meth is not an instance method, it is returned
    unchanged.

    If name is None, the generated function will have the same name as
    the input method.  Otherwise, name must be a string containing the
    desired name of the new method.  If modname is not None, it must
    be a string and will be used as the module name for the generated
    function.  Note that the caller is responsible for assigning the
    returned function to an appropriate name in the calling scope.

    """
    
    if type(meth) is not instancemethod:
        return meth

    if name is None:
        name = meth.func_name

    if name == meth.func_name:
        old_name = '_old_' + name
    else:
        old_name = meth.func_name

    # Make an argument list string, removing 'self'
    argspec = inspect.getargspec(meth)
    argspec[0].pop(0)
    argspec = inspect.formatargspec(argspec[0], argspec[1], argspec[2])

    # Create a wrapper function with no default arguments
    g = {old_name: meth}
    if modname is not None:
        g['__name__'] = modname
    fdef = 'def %s%s:  return %s%s' % (name, argspec, old_name, argspec)
    exec fdef in g

    # Create another new function from the one we just made, this time
    # adding the default arguments and doc string from the original method
    new_meth = g[name]
    new_meth = function(new_meth.func_code, new_meth.func_globals,
                        new_meth.func_name, meth.func_defaults,
                        new_meth.func_closure)
    new_meth.func_doc = meth.func_doc

    return new_meth


def get_keyword_names(func, skip=0):
    """

    Return a list containing the names of func's keyword arguments,
    skipping the first skip keywords.

    """

    argspec = inspect.getargspec(func)
    if argspec[3] is None:
        return []
    first = len(argspec[0]) - len(argspec[3])
    return argspec[0][first+skip:]


def get_keyword_defaults(func, skip=0):
    """

    Return a dictionary containing the default values of func's keyword
    arguments, skipping the first skip keywords.

    """

    argspec = inspect.getargspec(func)
    if argspec[3] is None:
        return {}
    first = len(argspec[0]) - len(argspec[3])
    return dict(izip(argspec[0][first+skip:], argspec[3][skip:]))


def get_func_usage(func):
    """
    
    Returns a string describing the function signature, similar to pydoc
    but flexible enough to pipe to slsh
    
    """
    argspec = inspect.getargspec(func)

    #Note: tot not used
    tot, num_args, num_kargs = get_num_args(func)

    msg = 'Usage: %s(' % func.__name__

    for arg in argspec[0]:
        if num_args > 0:
            msg = '%s %s,' % (msg, arg)
            num_args = num_args - 1
        
        else:
            msg = '%s [ %s,' % (msg, arg)

    if argspec[1] is not None:
        msg = '%s [ *%s,' % (msg, argspec[1])
        num_kargs = num_kargs + 1

    if argspec[2] is not None:
        msg = '%s [ **%s,' % (msg, argspec[2])
        num_kargs = num_kargs + 1

    msg = msg.strip(',')

    for i in xrange(num_kargs):
        msg = '%s]' % msg
    
    return '%s )' % msg


def get_num_args(func):
    """

    Return a tuple of the number of arguments.
    ( total number of args, number of non-keyword args, number of keyword args)
    
    """

    argspec = inspect.getargspec(func)
    num_args = 0
    num_kargs = 0
    
    if len(argspec[0]) != 0:
        num_args = len(argspec[0])
    
    if argspec[3] is not None:
        num_kargs = len(argspec[3])

    return (num_args, (num_args - num_kargs), num_kargs)


def print_fields(names, vals, converters={}):
    """

    Given a list of strings names and mapping vals, where names is a
    subset of vals.keys(), return a listing of name/value pairs
    printed one per line in the format '<name> = <value>'.  If a value
    is a NumPy array, print it in the format
    '<data type name>[<array size>]'.  Otherwise, use str(value).

    """

    width = max(len(n) for n in names)
    fmt = '%%-%ds = %%s' % width
    lines = []
    for n in names:
        v = vals[n]

        conv = None
        if converters:
            for t in converters:
                if isinstance(v, t):
                    conv = converters[t]
                    break

        if conv is not None:
            v = conv(v)
        elif isinstance(v, numpy.ndarray):
            v = '%s[%d]' % (numpy.typeNA[v.dtype.type], v.size)
        else:
            v = str(v)
        lines.append(fmt % (n, v))
    return '\n'.join(lines)


def create_expr(vals, mask=None, format='%s', delim='-'):
    """
    collapse a list of channels into an expression using hyphens
    and commas to indicate filtered intervals.
    """
    expr=[]

    if len(vals) == 0:
        return ''
    elif len(vals) == 1:
        return format % vals[0]

    diffs = numpy.apply_along_axis(numpy.diff, 0, vals)
    if mask is not None:
        index = numpy.arange(len(mask))
        diffs = numpy.apply_along_axis(numpy.diff, 0, index[mask])

    for ii, delta in enumerate(diffs):
        if ii == 0:
            expr.append(format % vals[ii])
            if delta != 1 or len(diffs)==1:
                expr.append(',')
            continue
        if delta == 1:
            if expr[-1] == ',':
                expr.append(format % vals[ii])
            if expr[-1] != delim:
                expr.append(delim)
        else:
            if not expr[-1] in (',',delim):
                expr.append(',')
            expr.append(format % vals[ii])
            expr.append(',')
    if expr[-1] in (',',delim):
        expr.append(format % vals[-1])

    return ''.join(expr)


def parse_expr(expr):
    """
    parse a filter expression into numerical components for notice/ignore
    e.g. ':2,4:5,7:8,10:'
    """
    res = []
    if expr is None or str(expr).strip() == '':
        res.append((None,None))
        return res
    vals = str(expr).strip().split(',')
    for val in vals:
        lo=None; hi=None
        interval = val.strip().split(':')
        if len(interval) == 1:
            lo = interval[0]
            if lo == '':
                lo=None
            hi = lo
        elif len(interval) > 1:
            lo=interval[0]
            hi=interval[1]
            if lo == '':
                lo=None
            if hi == '':
                hi=None
        else:
            raise TypeError("interval syntax requires a tuple, 'lo:hi'")
        if lo is not None:
            try:
                lo = float(lo)
            except ValueError:
                raise TypeError("Invalid lower bound '%s'" % str(lo))
        if hi is not None:
            try:
                hi = float(hi)
            except ValueError:
                raise TypeError("Invalid upper bound '%s'" % str(hi))
        res.append((lo, hi))
    return res


def calc_total_error(staterror=None, syserror=None):
    """

    Add statistical and systematic errors in quadrature.  The
    arguments should be numpy arrays or None.  If both arguments are
    arrays, a new array containing the element-wise square root of the
    sum of their squares is returned.  If one argument is None, the
    other is returned unaltered.  If both arguments are None, None is
    returned.

    """

    if (staterror is None) and (syserror is None):
	error = None
    elif (staterror is not None) and (syserror is None):
	error = staterror
    elif (staterror is None) and (syserror is not None):
	error = syserror
    else:
	error = numpy.sqrt(staterror*staterror + syserror*syserror)
    return error


def poisson_noise(x):
    """

    Return a random value from a Poisson distribution with mean x.  If
    x is an array, return an array of such values.  Uses
    numpy.random.poisson to generate the values.

    """

    x = numpy.asarray(x)

    # Using numpy.where() and indexing doesn't work with 0-d arrays, so
    # handle them separately
    if x.shape == ():
        x = SherpaFloat(x)
        if x <= 0.0:
            x = 0.0
        else:
            x = numpy.random.poisson(x)
        return SherpaFloat(x)

    x_out = numpy.zeros(x.shape, SherpaFloat)
    good = numpy.where(x > 0.0)
    x_out[good] = numpy.random.poisson(x[good])

    return x_out

def _convolve( a, b ):
    if len(a) != len(b):
        raise TypeError("Input arrays are not equal in length, a: %s b: %s" %
                      (len(a), len(b)))

    imag = numpy.fft.fft(a)*numpy.fft.fft(b)
    return numpy.asarray( numpy.fft.ifft(imag), dtype=SherpaFloat )

def numpy_convolve( a, b ):

    if a.ndim > 1 or b.ndim > 1:
        raise TypeError("numpy_convolution is 1D only")

    c = numpy.concatenate((a, numpy.zeros(len(b))))
    d = numpy.concatenate((b, numpy.zeros(len(a))))

    if len(a) > len(b):
        return _convolve(c,d)[:len(a)]
    
    return _convolve(c,d)[:len(b)]


def dataspace1d(start, stop, step=1, numbins=None):
    """
    Populates an integrated grid

    if numbins is None (default) -> numpy.arange(start,stop,step)

    if numbins is not None -> numpy.linspace(start, stop, numbins)

    """
    if start >= stop:
        raise TypeError("input should be start < stop, found start=%s stop=%s" %
                      (start, stop))

    if numbins is None:
        if step <= 0:
            raise TypeError("input should be step > 0, found step=%s" % step)

        if step >= (stop-start):
            raise TypeError("input has produced less than 2 bins, found start=%s stop=%s step=%s" % (start, stop, step))

    #xx = numpy.arange(start, stop, step, dtype=float)
    xx = sao_arange(start, stop, step)
    if numbins is not None:
        if numbins <= 1:
            raise TypeError("input should be numbins > 1, found numbins=%s" % numbins)
            
        xx = numpy.linspace(start, stop, numbins+1)

    xlo = numpy.array(xx[:-1])
    xhi = numpy.array(xx[1:])
    y = numpy.zeros(len(xlo), dtype=float)

    return xlo, xhi, y


def dataspace2d(dim):
    """
    Populates a blank image dataset
    """
    if not numpy.iterable(dim):
        raise TypeError("dim must be an array of dimensions")

    if len(dim) < 2:
        raise TypeError("dimensions for dataspace2d must be > 1")

    if dim[0] < 1 or dim[1] < 1:
        raise TypeError("dimensions should be > 0, found dim0 %s dim1 %s"
                      % (dim[0], dim[1]))

    x0 = numpy.arange(dim[0], dtype=numpy.float)+1.
    x1 = numpy.arange(dim[1], dtype=numpy.float)+1.

    x0, x1 = numpy.meshgrid(x0, x1)
    shape = tuple(x0.shape)
    x0 = x0.ravel()
    x1 = x1.ravel()
    y = numpy.zeros(numpy.prod(dim))

    return x0, x1, y, shape


def histogram1d( x, x_lo, x_hi ):
    x_lo = numpy.asarray(x_lo)
    x_hi = numpy.asarray(x_hi)

    x_lo.sort()
    x_hi.sort()

    return hist1d( numpy.asarray(x), x_lo, x_hi )


def histogram2d( x, y, x_grid, y_grid ):
    x_grid = numpy.asarray(x_grid)
    y_grid = numpy.asarray(y_grid)

    x_grid.sort()
    y_grid.sort()

    vals = hist2d( numpy.asarray(x), numpy.asarray(y), x_grid, y_grid )
    return vals.reshape( (len(x_grid),len(y_grid)) )


def interpolate( x, xp, yp ):
    n = len(xp)
    if n != len(yp):
        raise TypeError("input arrays are not of equal length")

    xp = numpy.asarray(xp, dtype=float)
    xp.sort()
    n1 = numpy.searchsorted(xp, x)

    if n1 == 0:
        return yp[n1]

    n0 = n1 - 1
    if n1 == n:
        if n0 == 0:
            return yp[n0]
        n1 = n0 - 1

    x0 = xp[n0]
    x1 = xp[n1]

    if x == x1:
        return yp[n0]
    if x1 == x0:
        return yp[n0]

    return yp[n0] + (yp[n1] - yp[n0]) / (x1 - x0) * ( x - x0 )


def is_binary_file( filename ):
    """
    boolean determining if the file 'filename' is binary

    """
    fd = open( filename, 'r')
    lines = fd.readlines(1024)

    if len(lines) == 0:
        return False

    # If a non-printable character is found in first 1024 --> binary
    for line in lines:
        for char in line:
            if char not in string.printable:
                return True

    return False


def get_midpoint(a):
    #return numpy.abs(a.max() - a.min())/2. + a.min()
    return numpy.abs(a.max() + a.min())/2.


def get_peak(y, x, xhi=None):
    return x[y.argmax()]


def get_valley(y, x, xhi=None):
    return x[y.argmin()]

def get_fwhm(y, x, xhi=None):
    half_max_val = y.max() / 2.0
    x_max = x[y.argmax()]
    for ii, val in enumerate(y[:y.argmax()]):
        if val >= half_max_val:
            return 2.0*numpy.abs(x[ii] - x_max)
    return x_max

def guess_fwhm(y, x, xhi=None, scale=1000):
    fwhm = get_fwhm(y, x, xhi)
    return {'val': fwhm, 'min': fwhm/scale, 'max': fwhm*scale }


def param_apply_limits(param_limits, par, limits=True, values=True):
    """

    apply the dictionary of guess values to parameter, also, save the
    defaults for rollback.

    """
    if limits and values:
        default_val = par.val
        par.set(param_limits['val'], param_limits['min'], param_limits['max'], 
                default_min = par.min, default_max = par.max)

        # set original value as default outside the property interface
        par._default_val = default_val

        # set guessed flag to enable user-defined limit reset
        par._guessed=True

    elif values:
        default_val = par.val
        par.set(param_limits['val'])

        # set original value as default outside the property interface
        par._default_val = default_val

    else:
        par.set(min = param_limits['min'], max = param_limits['max'],
                default_min = par.min, default_max = par.max)

        # set guessed flag to enable user-defined limit reset
        par._guessed=True


def get_amplitude_position(arr, mean=False):
    """
    Guess model amplitude and position of array
    returns (xval, xmin, xmax, xpos)

    """

    pos = xpos = xmin = xmax = xval = 0
    max = arr.max()
    min = arr.min()
    if((max > 0.0 and min >= 0.0) or
       (max > 0.0 and min < 0.0 and (abs(min) <= max ))):
        xpos = arr.argmax()
        if mean:
            xpos = numpy.where(arr==max)

        xmax = max*_guess_ampl_scale
        xmin = max/_guess_ampl_scale
        xval = max

    elif((max > 0.0 and min < 0.0 and abs(min) > max ) or
         (max == 0.0 and min < 0.0 ) or ( max < 0.0 )):
        xpos = arr.argmin()
        if mean:
            xpos = numpy.where(arr==min)

        xmax = min/_guess_ampl_scale
        xmin = min*_guess_ampl_scale
        xval = min

    return (xval, xmin, xmax, xpos)


def guess_amplitude(y, x, xhi=None):
    """
    Guess model parameter amplitude (val, min, max)

    """

    val, ymin, ymax, pos = get_amplitude_position(y)

    min = None; max = None
    if ymin != 0.0 or ymax != 0.0:
        min = ymin
        max = ymax

    if xhi is not None:
        binsize = numpy.abs(xhi[pos] - x[pos])
        min /= binsize
        max /= binsize
        val /= binsize

    return { 'val' : val, 'min' : min, 'max' : max }



def guess_amplitude_at_ref(r, y, x, xhi=None):
    """
    Guess model parameter amplitude (val, min, max)

    only valid for beta1d, bpl1d, powlaw1d

    """

    t = 1.0
    if x[1] > x[0] and r < x[0]:
        t = numpy.abs(y[0] + y[1])/2.0
    elif x[1] > x[0] and r > x[-1]:
        t = numpy.abs(y[-1] + y[-2])/2.0
    elif x[1] < x[0] and r > x[0]:
        t = numpy.abs( y[0] + y[1])/2.0
    elif x[1] < x[0] and r < x[-1]:
        t = numpy.abs( y[-1] + y[-2])/2.0
    else:
        for i in range(len(x)-1):
            if ( ( r >= x[i] and r < x[i+1] ) or ( r >= x[i+1] and r < x[i] ) ):
                t = numpy.abs(y[i] + y[i+1])/2.0
                break

    if t == 0.0:
        totband = 0.0
        dv = 0.0
        i = 1
        for j in range(len(x)):
            dv = x[i] - x[i-1]
            t += y[i]*dv
            totband += dv
            i += 1
        t /= totband;

    return { 'val':t, 'min':t/_guess_ampl_scale, 'max':t*_guess_ampl_scale }


def guess_amplitude2d(y, x0lo, x0hi, x1lo=None, x1hi=None):
    """
    Guess 2D model parameter amplitude (val, min, max)

    """

    limits = guess_amplitude(y, x0lo)

    if x1lo is not None and x1hi is not None:
        binsize = numpy.abs((x0hi[0]-x0lo[0])*(x1hi[0]-x1lo[0]))
        limits['min'] /= binsize
        limits['max'] /= binsize
        limits['val'] /= binsize

    return limits


def guess_reference(pmin, pmax, x, xhi=None):
    """
    Guess model parameter reference (val, min, max)

    """

    min = x.min()
    max = x.max()

    if min >= 1: pmin = 1
    if max <= 1: pmax = 1

    val = 0.0
    if min < 1.0 and max > 1.0:
        val = 1.0
    else:
        refval = numpy.floor((min+max)/2.)
        if refval < pmin or refval > pmax: 
            refval=(min+max)/2.
        val = refval

    return { 'val':val, 'min':None, 'max':None }


def get_position(y, x, xhi=None):
    """
    Get 1D model parameter positions pos (val, min, max)

    """
    xval, xmin, xmax, xpos = get_amplitude_position(y, mean=True)

    val = numpy.mean(x[xpos])
    min = x.min()
    max = x.max()
    if xhi is not None:
        max = xhi.max()

    return { 'val':val, 'min':min, 'max':max }


def guess_position(y, x0lo, x0hi, x1lo=None, x1hi=None):
    """
    Guess 2D model parameter positions xpos, ypos ({val0, min0, max0},
                                                   {val1, min1, max1})

    """
    if x1lo is None and x1hi is None:
        x0, x1 = x0lo, x0hi
    else:
        x0, x1 = x0lo, x1lo

    #pos = int(y.argmax())
    # return the average of location of brightest pixels
    pos = numpy.where(y==y.max())

    return tuple(({ 'val':numpy.mean(x0[pos]), 'min':x0.min(), 'max':x0.max() },
                  { 'val':numpy.mean(x1[pos]), 'min':x1.min(), 'max':x1.max() }
                  ))


def guess_bounds(x, xhi=True):
    """
    Guess model parameters xlo, xhi (val, min, max)

    """
    min = x.min()
    max = x.max()
    lo = (min+(max-min)/2.)
    if xhi:
        lo = (min+(max-min)/3.)
        hi = (min+2.*(max-min)/3.)
        return tuple(( { 'val':lo, 'min':min, 'max':max },
                       { 'val':hi, 'min':min, 'max':max } ))

    return { 'val':lo, 'min':min, 'max':max }


def guess_radius(x0lo, x0hi, x1lo=None, x1hi=None):
    """
    Guess 2D model parameter radius (val, min, max)

    """
    if x1lo is None and x1hi is None:
        x0, x1 = x0lo, x0hi
    else:
        x0, x1 = x0lo, x1lo

    delta = numpy.apply_along_axis(numpy.diff, 0, x0)[0]
    rad = numpy.abs(10*delta)

    return { 'val':rad, 'min':rad/_guess_ampl_scale, 'max':rad*_guess_ampl_scale }


def divide_run_parallel(func, array, *args, **kwargs):

    if not numpy.iterable(array):
        raise TypeError('No input array(s) found')

    array = numpy.asarray(array)

    if not _multi or _ncpus < 2 or len(array) < _ncpus:
        return func(array, *args, **kwargs)

    def worker(pid, queue, err, func, slice, args, kwargs):
        try:
            result = func(slice, *args, **kwargs)
        except Exception, e:
            # Serialize just the class and args instead of instance
            # Pickle has some crazy interpretation of the Exception
            # interface.  This causes problems for classes that 
            # inherit from SherpaErr
            err.put( (e.__class__,e.args) )
            return

        queue.put( (pid, result) )


    tasks = []
    die = (lambda tasks : [task.terminate() for task in tasks
                           if task.exitcode is None])

    ncpus = _ncpus
    manager = multiprocessing.Manager()
    q   = manager.Queue()
    err = manager.Queue()

    # partition the array space according to num cpus
    start = 0
    nelem = len(array)
    chunk = int(nelem / ncpus)
    resid = int(nelem % ncpus)
    stop  = chunk

    # process ids keep the pieces in order for final result array
    pid = 0

    # less the last CPU, partition the array space for each process
    for ii in xrange(ncpus-1):
        tasks.append(multiprocessing.Process(target=worker,
                             args=(pid, q, err, func, array[start:stop],
                                   args, kwargs)))
        start += chunk
        stop += chunk
        pid += 1


    # last CPU includes residual array elements (mod ncpus)
    tasks.append(multiprocessing.Process(target=worker,
                         args=(pid, q, err, func, array[start:stop+resid],
                               args, kwargs)))

    for task in tasks:
        task.start()

    for task in tasks:
        task.join()

    if not err.empty():
        die(tasks)
        # Obtain the class and args from the child process
        # create a new instance with args and then clobber
        # them to ensure this instance matches the one from
        # the child.  This works-around the constructor of
        # exceptions that inherit from SherpaErr.
        cls, args = err.get()
        ex = cls(*args)
        ex.args = args
        raise ex

    # use hash table to quickly sort pieces of result arrays by pid
    results = {}
    while not q.empty():
        id, array = q.get()
        results[id] = array

    return numpy.concatenate([results[key] for key in range(ncpus)])


############################### Root of all evil ##############################

def printf(format, *args):
    """Format args with the first argument as format string, and write.
    Return the last arg, or format itself if there are no args."""
    sys.stdout.write(str(format) % args)
    return if_(args, args[-1], format)

def caller(n=1):
    """Return the name of the calling function n levels up in the frame stack.
    >>> caller(0)
    'caller'
    >>> def f():
    ...     return caller()
    >>> f()
    'f'
    """
    import inspect
    return  inspect.getouterframes(inspect.currentframe())[n][3]

def func_wrapper( func, history ):
    '''A function wrapper to count the number of times beingg called'''
    nfev = [0]
    def func_counter_history( x, *args ):
        nfev[0] += 1
        y = func( x, *args )
        history[ 0 ].append( x )
        history[ 1 ].append( y )
        return y
    return nfev, func_counter_history

def is_in( arg, seq ):
    for x in seq:
        if arg is x: return True
    return False

def is_iterable( x ):
    return isinstance(x, (list, tuple, numpy.ndarray)) or numpy.iterable( x )
    
def is_sequence( start, mid, end ):
    return (start < mid) and (mid < end)

def Knuth_close( x, y, tol, myop=operator.__or__ ):
    """ The following text was taken verbatim from:
        
    http://www.boost.org/doc/libs/1_35_0/libs/test/doc/components/test_tools/floating_point_comparison.html#Introduction

    In most cases it is unreasonable to use an operator==(...)
    for a floating-point values equality check. The simple solution
    like abs(f1-f2) <= e does not work for very small or very big values.
    This floating-point comparison algorithm is based on the more
    confident solution presented by D. E. Knuth in 'The art of computer
    programming (vol II)'. For a given floating point values u and v and
    a tolerance e:
    
    | u - v | <= e * |u| and | u - v | <= e * |v|                    (1)
    defines a "very close with tolerance e" relationship between u and v
    
    | u - v | <= e * |u| or   | u - v | <= e * |v|                   (2)
    defines a "close enough with tolerance e" relationship between
    u and v. Both relationships are commutative but are not transitive.
    The relationship defined by inequations (1) is stronger that the
    relationship defined by inequations (2) (i.e. (1) => (2) )."""
    
    diff = abs( x - y )
    if 0.0 == x or 0.0 == y:
        return diff <= tol
    return myop( diff <= tol * abs( x ), diff <= tol * abs( y ) )

def safe_div( num, denom ):
    import sys
    dbl_max = sys.float_info.max
    dbl_min = sys.float_info.min

    # avoid overflow
    if denom < 1 and num > denom * dbl_max:
        return dbl_max

    # avoid underflow
    if 0.0 == num or denom > 1 and num < denom * dbl_min:
        return 0
    
    return num / denom

def Knuth_boost_close( x, y, tol, myop=operator.__or__ ):
    """ The following text was taken verbatim from:
        
    http://www.boost.org/doc/libs/1_35_0/libs/test/doc/components/test_tools/floating_point_comparison.html#Introduction

    In most cases it is unreasonable to use an operator==(...)
    for a floating-point values equality check. The simple solution
    like abs(f1-f2) <= e does not work for very small or very big values.
    This floating-point comparison algorithm is based on the more
    confident solution presented by D. E. Knuth in 'The art of computer
    programming (vol II)'. For a given floating point values u and v and
    a tolerance e:
    
    | u - v | <= e * |u| and | u - v | <= e * |v|                    (1)
    defines a "very close with tolerance e" relationship between u and v
    
    | u - v | <= e * |u| or   | u - v | <= e * |v|                   (2)
    defines a "close enough with tolerance e" relationship between
    u and v. Both relationships are commutative but are not transitive.
    The relationship defined by inequations (1) is stronger that the
    relationship defined by inequations (2) (i.e. (1) => (2) ).
    Because of the multiplication in the right side of inequations,
    that could cause an unwanted underflow condition, the implementation
    is using modified version of the inequations (1) and (2) where all
    underflow, overflow conditions could be guarded safely:
    
    | u - v | / |u| <= e and | u - v | / |v| <= e          	     (1`)
    | u - v | / |u| <= e or  | u - v | / |v| <= e                    (2`)"""
    
    diff = abs( x - y )
    if 0.0 == x or 0.0 == y:
        return diff <= tol
    diff_x = safe_div( diff, x )
    diff_y = safe_div( diff, y )
    return myop( diff_x <= tol, diff_y <= tol )

def list_to_open_interval( arg ):
    if False == numpy.iterable( arg ):
        return arg
    str = '(%e, %e)' % (arg[0],arg[1])
    return str

def mysgn( arg ):
    if arg == 0.0:
        return 0
    elif arg < 0.0:
        return -1
    else:
        return 1

class OutOfBoundErr:
    pass

class QuadEquaRealRoot:
    """ solve for the real roots of the quadratic equation:
    a * x^2 + b * x + c = 0"""
    
    def __call__( self, a, b, c ):
        
        if 0.0 == a:
            #
            # 0 * x^2 + b * x + c = 0
            #

            if 0.0 != b:
                #
                # 0 * x^2 + b * x + c = 0
                # the folowing still works even if c == 0
                #
                answer = - c / b
                return [ answer, answer ]
            
            else:
                
                #
                # 0 * x^2 + 0 * x + c = 0
                #
                # a == 0, b == 0, so if c == 0 then all numbers work so
                # returning nan is not right. However if c != 0 then no
                # roots exist.
                #
                return [ None, None ]

        elif 0.0 == b:
            
            #
            # a * x^2 + 0 * x + c = 0
            #
            if 0.0 == c:
                
                # a * x^2 + 0 * x + 0 = 0
                return [ 0.0, 0.0 ]
            else:
                
                # a * x^2 + 0 * x + c = 0
                if mysgn( a ) == mysgn( c ):
                    return [ None, None ]
                answer = numpy.sqrt( c / a )
                return [ -answer, answer ]

        elif 0.0 == c:
            
            #
            # a * x^2 + b * x + 0 = 0
            #
            return [ 0.0, - b / a ]

        else:

            discriminant = b * b - 4.0 * a * c
            print 'disc=', discriminant
            sqrt_disc = numpy.sqrt( discriminant )
            t = - ( b + mysgn( b ) * sqrt_disc ) / 2.0
            return [ c / t, t / a ]

            
def bisection( fcn, xa, xb, fa=None, fb=None, args=(), maxfev=48, tol=1.0e-6 ):

    history = [ [], [] ]
    nfev, myfcn = func_wrapper( fcn, history )

    try:

        if fa is None:
            fa = myfcn( xa, *args )
        if abs( fa ) <= tol:
            return [ [xa, fa], [ [xa, fa], [xa, fa] ], nfev[0] ]
        
        if fb is None:
            fb = myfcn( xb, *args )
        if abs( fb ) <= tol:
            return [ [xb, fb], [ [xb, fb], [xb, fb] ], nfev[0] ]

        if mysgn( fa ) == mysgn( fb ):
            sys.stderr.write( __name__ + ': ' + fcn.__name__ +
                              ' fa * fb < 0 is not met\n' )
            return [ [None, None], [ [None, None], [None, None] ], nfev[0] ]

        while nfev[0] < maxfev:

            if abs( fa ) > tol and abs( fb ) > tol:

                xc = ( xa + xb ) / 2.0
                fc = myfcn( xc, *args )

                if abs( xa - xb ) < min( tol * abs( xb ), tol / 10.0 ):
                    return [ [xc, fc], [ [xa, fa], [xb, fb] ], nfev[0] ]

                if mysgn( fa ) !=  mysgn( fc ):
                    xb = xc; fb = fc                
                else:
                    xa = xc; fa = fc
                    
            else:
                if abs( fa ) <= tol:
                    return [ [xa, fa], [ [xa, fa], [xb, fb] ], nfev[0] ]
                else:
                    return [ [xb, fb], [ [xa, fa], [xb, fb] ], nfev[0] ]

                
        xc = ( xa + xb ) / 2.0
        fc = myfcn( xc, *args )
        return [ [xc, fc], [ [xa, fa], [xb, fb] ], nfev[0] ]

    except OutOfBoundErr:
        return [ [None, None], [ [xa, fa], [xb, fb] ], nfev[0] ]

def quad_coef( x, f ):
    """
    p( x ) = f( xc ) + A ( x - xc ) + B ( x - xc ) ( x - xb )
           = f( xc ) + A ( x - xc ) + B ( ( x - xc ) ( x - xc ) +
                                            ( x - xc ) ( xc - xb ) )
           = f( xc ) + ( A + B ( xc - xb ) ) ( x - xc ) + B ( x - xc )^2
    
           = f( xc ) + C ( x - xc ) + B ( x - xc )^2 ; C = A + B ( xc - xb )
                    
           = f( xc ) + C x - C xc + B ( x^2  - 2 x xc + xc^2 )
                    
           = B x^2 + ( C - 2 * B xc ) x + f( xc ) - C xc  + B xc^2
    
           = B x^2 + ( C - 2 * B x[2] ) x + f[ 2 ] + x[2] * ( B x[ 2 ] - C )"""
    

    [B,C] = transformed_quad_coef( x, f )
    B_x2 = B * x[ 2 ]
    return [ B, C - 2 * B_x2, f[ 2 ] + x[ 2 ] * ( B_x2 - C ) ]
     
def transformed_quad_coef( x, f ):
    """
    p( x ) = f( xc ) + A ( x - xc ) + B ( x - xc ) ( x - xb )
    
       where A and B are the divided differences:
        
                                   f( xc ) - f( xb )
                               A = -----------------
                                        xc - xb
    
    
                               1     ( f( xc ) - f( xb )   f( xb ) - f( xa ) )
                        B = -------  ( ----------------- - ----------------- )
                            xc - xa  (    xc - xb               xb - xa      )
    
        p( x ) = f( xc ) + A ( x - xc ) + B ( x - xc ) ( x - xb )
               = f( xc ) + A ( x - xc ) + B ( ( x - xc ) ( x - xc ) +
                                            ( x - xc ) ( xc - xb ) )
               = f( xc ) + ( A + B ( xc - xb ) ) ( x - xc ) + B ( x - xc )^2
    
               = f( xc ) + C ( x - xc ) + B ( x - xc )^2
               
        where  C = A + B ( xc - xb )
    
        The root of p( x ), using the quadratic formula:
    
                            1  (                                   )
                  x - xc = --- ( - C +/- sqrt( C^2 - 4 f( xc ) B ) )
                           2 B (                                   )
                           
        Rationalize the numerator to avoid subtractive cancellation:
    
                                     2 f( xc )
                  x - xc = -------------------------------
                           C +/- sqrt( C^2 - 4 f( xc ) B )
    
        The sign should be chosen to maximize the denominator.  Therefore,
        the next point in the iteration is:
        
    
                                       2 f( xc )
                  x = xc - --------------------------------------
                           C + sgn( C ) sqrt( C^2 - 4 f( xc ) B )
                           
                       {    -1,  x < 0
        where sgn(x) = {
                       {     1,  x >= 0"""

    xa = x[ 0 ]; xb = x[ 1 ]; xc = x[ 2 ]
    fa = f[ 0 ]; fb = f[ 1 ]; fc = f[ 2 ]

    xc_xb = xc - xb
    fc_fb = fc - fb
    A = fc_fb / xc_xb
    fb_fa = fb - fa
    xb_xa = xb - xa
    xc_xa = xc - xa
    B = ( A - fb_fa / xb_xa ) / xc_xa
    C = A + B * xc_xb
    return [ B, C ]

def demuller( fcn, xa, xb, xc, fa=None, fb=None, fc=None, args=(), 
              maxfev=32, tol=1.0e-6 ):
    """
    p( x ) = f( xc ) + A ( x - xc ) + B ( x - xc ) ( x - xb )
    
    The general case:
     
                                     2 f( x )
                                           n
               x   = x  -  ----------------------------------------
                n+1   n    C  + sgn( C  ) sqrt( C^2 - 4 f( x  ) B )
                            n         n          n          n    n
    
                           1     ( f( x  ) - f( x   )   f( x   ) - f( x   )  )
                                 (     n         n-1        n-1        n-2   )
                   B  = -------  ( ------------------ - -------------------  )
                    n   x - x    (    x - x                 x   - x          )
                         n   n-2 (     n   n-1               n-1   n-2       )
                         
                        
                        f( x  ) - f( x   )
                            n         n-1
                   A  = -----------------
                    n       x  - x
                             n    n-1
    
                   C  = A  + B ( x - x   )
                    n    n    n   n   n-1
    
    
    The convergence rate for Muller's method can be shown to be
    the real root of the cubic x - x^3, that is,
    p = (a + 4 / a + 1) / 3 where a = (19 + 3√33)^1/3.
    In other words: O(h^p) where p ≈ 1.839286755.
    
    >>> import math
    >>> a = pow( 19.0 + 3 * math.sqrt( 33 ), 1/3.0 )
    >>> p = ( a + 4 / a + 1 ) / 3.0
    >>> print p
    1.83928675521"""
    

    history = [ [], [] ]
    nfev, myfcn = func_wrapper( fcn, history )

    try:
     
        if fa is None:
            fa = myfcn( xa, *args )
        if abs( fa ) <= tol:
            return [ [xa, fa], [ [xa, fa], [xa, fa] ], nfev[0] ]

        if fb is None:
            fb = myfcn( xb, *args )
        if abs( fb ) <= tol:
            return [ [xb, fb], [ [xb, fb], [xb, fb] ], nfev[0] ]

        if fc is None:
            fc = myfcn( xc, *args )
        if abs( fc ) <= tol:
            return [ [xc, fc], [ [xc, fc], [xc, fc] ], nfev[0] ]

        while nfev[0] < maxfev:
            
            [B, C] = transformed_quad_coef( [xa, xb, xc], [fa, fb, fc] )
            
            discriminant = max( C * C - 4.0 * fc * B, 0.0 )

            xd = xc - 2.0 * fc / ( C + mysgn( C ) * numpy.sqrt( discriminant ) )

            fd = myfcn( xd, *args )
            #print 'fd(%e)=%e' % (xd, fd)
            if abs( fd ) <= tol:
                return [ [xd, fd], [ [None, None], [None, None] ], nfev[0] ]

            xa = xb
            fa = fb
            xb = xc
            fb = fc
            xc = xd
            fc = fd

        #print 'demuller(): maxfev exceeded'
        return [ [xd, fd], [ [None, None], [None, None] ], nfev[0] ]

    except ZeroDivisionError:
        
        #print 'demuller(): fixme ZeroDivisionError'
        #for x, y in izip( history[0], history[1] ):
        #    print 'f(%e)=%e' % ( x, y )
        return [ [xd, fd], [ [None, None], [None, None] ], nfev[0] ]


def new_muller( fcn, xa, xb, fa=None, fb=None, args=(), maxfev=32, tol=1.e-6 ):

    def regula_falsi( x0, x1, f0, f1 ):
        if f0 < f1:
            xl = x0; fl = f0
            xh = x1; fh = f1
        else:
            xl = x1; fl = f1
            xh = x0; fh = f0
        x = xl + ( xh - xl ) * fl / ( fl - fh )
        if is_sequence( x0, x, x1 ):
            return x
        else:
            return ( x0 + x1 ) / 2.0
        
    history = [ [], [] ]
    nfev, myfcn = func_wrapper( fcn, history )

    try:

        if fa is None:
            fa = myfcn( xa, *args )
        if abs( fa ) <= tol:
            return [ [xa, fa], [ [xa, fa], [xa, fa] ], nfev[0] ]

        if fb is None:
            fb = myfcn( xb, *args )
        if abs( fb ) <= tol:
            return [ [xb, fb], [ [xb, fb], [xb, fb] ], nfev[0] ]

        if mysgn( fa ) == mysgn( fb ):
            sys.stderr.write( __name__ + ': ' + fcn.__name__ +
                              ' fa * fb < 0 is not met\n' )
            return [ [None, None], [ [None, None], [None, None] ], nfev[0] ]

        while nfev[0] < maxfev:

            xc = ( xa + xb ) / 2.0
            fc = myfcn( xc, *args )
            if abs( fc ) <= tol:
                return [ [xc, fc], [ [xa, fa], [xb, fb] ], nfev[0] ]

            tran = transformed_quad_coef( [xa, xb, xc], [fa, fb, fc] )
            B = tran[ 0 ]
            C = tran[ 1 ]
            
            discriminant = max( C * C - 4.0 * fc * B, 0.0 )

            xd = xc - 2.0 * fc / ( C + mysgn( C ) * numpy.sqrt( discriminant ) )

            fd = myfcn( xd, *args )
            #print 'fd(%e)=%e' % (xd, fd)
            if abs( fd ) <= tol:
                return [ [xd, fd], [ [xa, fa], [xb, fb] ], nfev[0] ]

            if mysgn( fa ) != mysgn( fc ):
                xb = xc; fb = fc
                continue

            if mysgn( fd ) != mysgn( fc ) and xc < xd:
                xa = xc; fa = fc
                xb = xd; fb = fd
                continue

            if mysgn( fb ) != mysgn( fd ):
                xa = xd; fa = fd
                continue

            if mysgn( fa ) != mysgn( fd ):
                xb = xd; fb = fd
                continue

            if mysgn( fc ) != mysgn( fd ) and xd < xc:
                xa = xd; fa = fd
                xb = xc; fb = fc
                continue

            if mysgn( fc ) != mysgn( fd ):
                xa = xc; fa = fc
                continue
            
        #print 'new_muller(): maxfev exceeded'
        return [ [xd, fd], [ [xa, fa], [xb, fb] ], nfev[0] ]

    except ZeroDivisionError, OutOfBoundErr:
        
        #print 'new_muller(): fixme ZeroDivisionError'
        #for x, y in izip( history[0], history[1] ):
        #    print 'f(%e)=%e' % ( x, y )
        return [ [xd, fd], [ [xa, fa], [xb, fb] ], nfev[0] ]

#
# /*
#  * Licensed to the Apache Software Foundation (ASF) under one or more
#  * contributor license agreements.  See the NOTICE file distributed with
#  * this work for additional information regarding copyright ownership.
#  * The ASF licenses this file to You under the Apache License, Version 2.0
#  * (the "License"); you may not use this file except in compliance with
#  * the License.  You may obtain a copy of the License at
#  *
#  *      http://www.apache.org/licenses/LICENSE-2.0
#  *
#  * Unless required by applicable law or agreed to in writing, software
#  * distributed under the License is distributed on an "AS IS" BASIS,
#  * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  * See the License for the specific language governing permissions and
#  * limitations under the License.
#  */
# 
def apache_muller( fcn, xa, xb, fa=None, fb=None, args=(), maxfev=32,
                   tol=1.0e-6 ):

    history = [ [], [] ]
    nfev, myfcn = func_wrapper( fcn, history )

    try:

        if fa is None:
            fa = myfcn( xa, *args )
        if abs( fa ) <= tol:
            return [ [xa, fa], [ [xa,fa], [xa,fa] ], nfev[0] ]

        if fb is None:
            fb = myfcn( xb, *args )
        if abs( fb ) <= tol:
            return [ [xb, fb], [ [xb,fb], [xb,fb] ], nfev[0] ]

        if mysgn( fa ) == mysgn( fb ):
            sys.stderr.write( __name__ + ': ' + fcn.__name__ +
                              ' fa * fb < 0 is not met\n' )
            return [ [None, None], [ [None, None], [None, None] ], nfev[0] ]

        xc = ( xa + xb ) / 2.0
        fc = myfcn( xc, *args )
        #print 'MullerBound() fc(%.14e)=%.14e' % (xc,fc)
        if abs( fc ) <= tol:
            return [ [ xc, fc ], [ [xc, fc], [xc, fc] ], nfev[0] ]

        xbest = xa; fbest = fa
        if abs( fb ) < abs( fa ):
            xbest = xb; fbest = fb
        if abs( fc ) < abs( fbest ):
            xbest = xc; fbest = fc
            
        oldx = 1.0e128
        while nfev[0] < maxfev:

            tran = transformed_quad_coef( [xa, xb, xc], [fa, fb, fc] )
            B = tran[ 0 ]
            C = tran[ 1 ]
            discriminant = max( C * C - 4.0 * fc * B, 0.0 )
            den = mysgn( C ) * numpy.sqrt( discriminant )
            xplus = xc - 2.0 * fc / ( C + den )
            if C != den:
                xminus = xc - 2.0 * fc / ( C - den )
            else:
                xminus = 1.0e128

            if is_sequence(xa, xplus, xb):
                x = xplus
            else:
                x = xminus;

            #print 'xa=', xa, '\tx=', x, '\txb=', xb, '\txc=', xc

            #fubar = quad_coef( [xa,xb,xc], [fa,fb,fc] )
            #quad = QuadEquaRealRoot( )
            #print quad( fubar[0], fubar[1], fubar[2] )
            #print

            #sanity check
            if False == is_sequence( xa, x, xb ):
                x = ( xa + xb ) / 2.0

            y = myfcn( x, *args );
            #print 'MullerBound() y(%.14e)=%.14e' % (x,y)
            if abs( y ) < abs( fbest ):
                xbest = x; fbest = y
            tolerance = min( tol * abs( x ), tol )
            if abs( y ) <= tol or abs( x - oldx ) <= tolerance:
                return [ [x,y], [ [xa, fa], [xb, fb] ], nfev[0] ]

            mybisect = (x < xc and (xc - xa) > 0.95 * (xb - xa)) or \
                       (x > xc and (xb - xc) > 0.95 * (xb - xa)) or \
                       (x == xc)

            if False == mybisect:
                if x > xc:
                    xa = xc
                    fa = fc
                if x < xc:
                    xb = xc
                    fb = fc
                xc = x; fc = y;
                oldx = x
            else:
                xmid = ( xa + xb ) / 2.0
                fmid = myfcn( xmid, *args )
                if abs( fmid ) < abs( fbest ):
                    xbest = xmid; fbest = fmid
                #print 'MullerBound() fmid(%.14e)=%.14e' % (xmid,fmid)
                if abs( fmid ) <= tol:
                    return [ [xmid, fmid], [ [xa, fa], [xb, fb] ], nfev[0] ]
                if mysgn( fa ) + mysgn( fmid ) == 0:
                    xb = xmid
                    fb = fmid
                else:
                    xa = xmid
                    fa = fmid
                xc = ( xa + xb ) / 2.0
                fc = myfcn( xc, *args )
                if abs( fc ) < abs( fbest ):
                    xbest = xc; fbest = fc                
                #print 'MullerBound() fc(%.14e)=%.14e' % (xc,fc)
                if abs( fc ) <= tol:
                    return [ [xc, fc], [ [xa, fa], [xb, fb] ], nfev[0] ]
                oldx = 1.0e128

        #
        # maxfev has exceeded, return the minimum so far
        #
        return [ [ xbest, fbest ], [ [xa, fa], [xb, fb] ], nfev[0] ]

    #
    # Something drastic has happened
    #
    except ZeroDivisionError, OutOfBoundErr:
        
        return [ [xbest, fbest], [ [xa, fa], [xb, fb] ], nfev[0] ]


def zeroin( fcn, xa, xb, fa=None, fb=None, args=(), maxfev=32, tol=1.0e-2 ):
    """
    *
    *			    Brent's root finder
    *	       obtains a zero of a function of one variable
    *
    * Synopsis
    *	Zeroin returns an approximate location for the root with accuracy
    *	4*DBL_EPSILON*abs(x) + tol
    *
    * Algorithm
    *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
    *	computations. M., Mir, 1980, p.180 of the Russian edition
    *
    * The function makes use of a bisection procedure combined with
    * a linear or quadratic inverse interpolation.
    * At each step the code operates three abscissae - a, b, and c:
    *	b - the last and the best approximation to the root
    *	a - the last but one approximation
    *	c - the last but one or even an earlier approximation such that
    *		1) |f(b)| <= |f(c)|
    *		2) f(b) and f(c) have opposite signs, i.e. b and c encompass
    *		   the root
    * Given these abscissae, the code computes two new approximations,
    * one by the  bisection procedure and the other one from interpolation
    * (if a,b, and c are all different the quadratic interpolation is used,
    * linear otherwise). If the approximation obtained by the interpolation
    * looks reasonable (i.e. falls within the current interval [b,c], not
    * too close to the end points of the interval), the point is accepted
    * as a new approximation to the root. Otherwise, the result of the
    * bissection is used.
    ************************************************************************"""

    history = [ [], [] ]
    nfev, myfcn = func_wrapper( fcn, history )

    try:
        
        if fa is None:
            fa = myfcn( xa, *args )
            if abs( fa ) <= tol:
                return [ [xa, fa], [ [xa, fa], [xb, fb] ], nfev[0] ]

        if fb is None:
            fb = myfcn( xb, *args )
        if abs( fb ) <= tol:
            return [ [xb, fb], [ [xa, fa], [xb, fb] ], nfev[0] ]    

        if mysgn( fa ) == mysgn( fb ):
            sys.stderr.write( __name__ + ': ' + fcn.__name__ +
                              ' fa * fb < 0 is not met\n' )
            return [ [None, None], [ [None, None], [None, None] ], nfev[0] ]
    
        xc = xa
        fc = fa
        DBL_EPSILON = numpy.float_(numpy.finfo(numpy.float32).eps)
        while nfev[0] < maxfev:
        
            prev_step = xb - xa

            if abs( fc ) < abs( fb ):
                xa = xb; fa = fb
                xb = xc; fb = fc
                xc = xa; fc = fa

            tol_act = 2.0 * DBL_EPSILON * abs( xb ) + tol / 2.0
            new_step = ( xc - xb ) / 2.0

            if abs( fb ) <= tol:
                return [ [xb, fb], [ [xa, fa], [xb, fb] ], nfev[0] ]

            if abs( new_step ) <= tol_act:
                if mysgn( fb ) != mysgn( fa ):
                    tmp = apache_muller( fcn, xa, xb, fa, fb, args=args,
                                         maxfev=maxfev-nfev[0], tol=tol )
                    tmp[ -1 ] += nfev[0]
                    return tmp
                elif mysgn( fb ) != mysgn( fc ):
                    tmp = apache_muller( fcn, xb, xc, fb, fc, args=args,
                                         maxfev=maxfev-nfev[0], tol=tol )
                    tmp[ -1 ] += nfev[0]
                    return tmp                
                else:
                    return [ [xb, fb], [ [xa, fa], [xb, fb] ], nfev[0] ]
            
            if abs( prev_step ) >= tol_act and abs( fa ) > abs( fb ):

                cb = xc - xb
                if xa == xc:
                    t1 = fb / fa
                    p = cb * t1
                    q = 1.0 - t1
                else:
                    t1 = fb / fc
                    t2 = fb / fa
                    q = fa / fc
                    p = t2 * ( cb * q * (q-t1) - (xb-xa) * (t1-1.0) )
                    q = (q-1.0) * (t1-1.0) * (t2-1.0)

                if p > 0:
                    q = -q
                else:
                    p = -p
                
                if 2*p < (1.5*cb*q-abs(tol_act*q)) and 2*p < abs(prev_step*q):
                    new_step = p / q

            if abs( new_step ) < tol_act:
                if new_step > 0:
                    new_step = tol_act
                else:
                    new_step = -tol_act
            xa = xb
            fa = fb
            xb += new_step
            fb = myfcn( xb, *args )
            #print 'fa(%f)=%f\tfb(%f)=%f\tfc(%f)=%f' % (xa,fa,xb,fb,xc,fc)
            if fb > 0 and fc > 0 or fb < 0 and fc < 0:
                xc = xa
                fc = fa

        return [ [xb, fb], [ [xa, fa], [xc, fc] ], nfev[0] ]

    except ZeroDivisionError, OutOfBoundErr:
        return [ [xb, fb], [ [xa, fa], [xc, fc] ], nfev[0] ]


############################### Root of all evil ##############################

