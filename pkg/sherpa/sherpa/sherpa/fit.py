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

from itertools import izip
import logging
import os
import signal
from numpy import power, arange, array, abs, iterable, sqrt, where, \
     ones_like, isnan, isinf
from sherpa.utils import NoNewAttributesAfterInit, print_fields, erf, igamc, \
    bool_cast, is_in, is_iterable, list_to_open_interval
from sherpa.utils.err import FitErr, EstErr
from sherpa.data import DataSimulFit
from sherpa.estmethods import Covariance, EstNewMin
from sherpa.models import SimulFitModel, Parameter
from sherpa.optmethods import LevMar, NelderMead
from sherpa.stats import Chi2, Chi2Gehrels, Cash, CStat, Chi2ModVar, LeastSq, \
    Likelihood


warning = logging.getLogger(__name__).warning
info = logging.getLogger(__name__).info

__all__ = ('FitResults', 'ErrorEstResults', 'Fit')


class FitResults(NoNewAttributesAfterInit):

    _fields = ('datasets', 'methodname', 'statname', 'succeeded',
               'parnames', 'parvals', 'statval', 'istatval',
               'dstatval', 'numpoints', 'dof', 'qval', 'rstat', 'message',
               'nfev')

    def __init__(self, fit, results, init_stat):
        _vals   = fit.data.eval_model_to_fit(fit.model)
        _dof    = len(_vals) - len(tuple(results[1]))
        _qval   = None
        _rstat  = None
        _covarerr = results[4].get('covarerr')
        if (isinstance(fit.stat, (CStat,Chi2)) and
            not isinstance(fit.stat, LeastSq)):
            if results[2] >= 0.0:
                _qval = igamc(_dof/2., results[2]/2.)
            _rstat = results[2]/_dof

        self.succeeded    = results[0]
        self.parnames     = tuple(p.fullname for p in fit.model.pars
                                  if not p.frozen)
        self.parvals      = tuple(results[1])
        self.istatval     = init_stat
        self.statval      = results[2]
        self.dstatval     = abs(results[2]-init_stat)
        self.numpoints    = len(_vals)
        self.dof          = _dof
        self.qval         = _qval
        self.rstat        = _rstat
        self.message      = results[3]
        if _covarerr is not None:
            self.covarerr      = tuple(_covarerr)
        else:
            self.covarerr      = None
        self.nfev         = results[4].get('nfev')
        self.extra_output = results[4]
        self.modelvals    = _vals
        self.methodname   = type(fit.method).__name__.lower()
        self.statname     = type(fit.stat).__name__.lower()
        self.datasets     = None # To be filled by calling function
        NoNewAttributesAfterInit.__init__(self)

    def __nonzero__(self):
        return self.succeeded

    def __repr__(self):
        return '<Fit results instance>'

    def __str__(self):
        return print_fields(self._fields, vars(self))

    def format(self):
        s = ''
        if (self.datasets != None):
            if len(self.datasets) == 1:
                s = 'Dataset               = %s\n' % str(self.datasets[0])
            else:
                s = 'Datasets              = %s\n' % str(self.datasets).strip("()")
        s += 'Method                = %s\n' % self.methodname
        s += 'Statistic             = %s\n' % self.statname 
        s += 'Initial fit statistic = %g\n' % self.istatval
        s += 'Final fit statistic   = %g' % self.statval
        if self.nfev is not None:
            s += ' at function evaluation %d' % self.nfev

        s += '\nData points           = %g' % self.numpoints
        s += '\nDegrees of freedom    = %g' % self.dof
        
        if self.qval is not None:
            s += '\nProbability [Q-value] = %g' % self.qval
        if self.rstat is not None:
            s += '\nReduced statistic     = %g' % self.rstat
        s += '\nChange in statistic   = %g' % self.dstatval

        if self.covarerr is None:
            for name, val in izip(self.parnames, self.parvals):
                s += '\n   %-12s   %-12g' % (name, val)
        else:
            for name, val, covarerr in izip(self.parnames, self.parvals, self.covarerr):
                s += '\n   %-12s   %-12g +/- %-12g' % (name, val, covarerr)

        return s

class ErrorEstResults(NoNewAttributesAfterInit):

    _fields = ('datasets', 'methodname', 'fitname', 'statname', 'sigma',
               'percent', 'parnames', 'parvals', 'parmins', 'parmaxes',
               'nfits')

    def __init__(self, fit, results, parlist = []):
        if (parlist == []):
            parlist = [p for p in fit.model.pars if not p.frozen]

        from sherpa.estmethods import est_success, est_hardmin, est_hardmax, est_hardminmax
        
        warning_hmin      = "hard minimum hit for parameter "
        warning_hmax      = "hard maximum hit for parameter "
        self.datasets     = None # To be set by calling function
        self.methodname   = type(fit.estmethod).__name__.lower()
        self.fitname      = type(fit.method).__name__.lower()
        self.statname     = type(fit.stat).__name__.lower()
        self.sigma        = fit.estmethod.sigma
        self.percent      = erf(self.sigma / sqrt(2.0)) * 100.0
        self.parnames     = tuple(p.fullname for p in parlist if not p.frozen)
        self.parvals      = tuple(p.val for p in parlist if not p.frozen)
        self.parmins      = ()
        self.parmaxes     = ()
        self.nfits        = 0
        success           = True
        for i in range(len(parlist)):
            if (results[2][i] != est_success):
                success = False
            if (results[2][i] == est_hardmin or
                results[2][i] == est_hardminmax):
                self.parmins  = self.parmins + (None,)
                warning(warning_hmin + self.parnames[i])
            else:
                self.parmins  = self.parmins + (results[0][i],)

            if (results[2][i] == est_hardmax or
                results[2][i] == est_hardminmax):
                self.parmaxes  = self.parmaxes + (None,)
                warning(warning_hmax + self.parnames[i])
            else:
                self.parmaxes  = self.parmaxes + (results[1][i],)

        self.nfits = results[3]
        self.extra_output = results[4]

        NoNewAttributesAfterInit.__init__(self)

    def __repr__(self):
        return '<%s results instance>' % self.methodname

    def __str__(self):
        return print_fields(self._fields, vars(self))

    def format(self):
        s = ""
        if (self.datasets != None):
            if len(self.datasets) == 1:
                s = 'Dataset               = %s\n' % str(self.datasets[0])
            else:
                s = 'Datasets              = %s\n' % str(self.datasets).strip("()")
        s += 'Confidence Method     = %s\n' % self.methodname
        s += 'Fitting Method        = %s\n' % self.fitname
        s += 'Statistic             = %s\n' % self.statname 
        
        s += "%s %g-sigma (%2g%%) bounds:" % (self.methodname, self.sigma,
                                              self.percent)

        def myformat( hfmt, str, lowstr, lownum, highstr, highnum ):
            str += hfmt % ('Param', 'Best-Fit', 'Lower Bound', 'Upper Bound')
            str += hfmt % ('-'*5, '-'*8, '-'*11, '-'*11)

            for name, val, lower, upper in izip(self.parnames, self.parvals,
                                                self.parmins, self.parmaxes):

                str += '\n   %-12s %12g ' % (name, val)
                if is_iterable( lower ):
                    str += ' '
                    str += list_to_open_interval( lower )
                elif (lower is None):
                    str += lowstr % '-----'
                else:
                    str += lownum % lower
                if is_iterable( upper ):
                    str += '  '
                    str += list_to_open_interval( upper )
                elif (upper is None):
                    str += highstr % '-----'
                else:
                    str += highnum % upper

            return str

        low = map( is_iterable, self.parmins )
        high = map( is_iterable, self.parmaxes )
        in_low = is_in( True, low )
        in_high = is_in( True, high )
        mymethod = self.methodname == 'confidence'

        lowstr = '%12s '
        lownum = '%12g '
        highstr = '%12s'
        highnum = '%12g'

        if True == in_low and True == in_high and mymethod:
            hfmt = '\n   %-12s %12s %29s %29s'
            lowstr = '%29s '
            lownum = '%29g '
            highstr = '%30s'
            highnum = '%30g'
        elif True == in_low and False == in_high and mymethod:
            hfmt = '\n   %-12s %12s %29s %12s'
            lowstr = '%29s '
            lownum = '%29g '
            highstr = '%13s'
            highnum = '%13g'
        elif False == in_low and True == in_high and mymethod:
            hfmt = '\n   %-12s %12s %12s %29s'
            highstr = '%29s'
            highnum = '%29g'            
        else:
            hfmt = '\n   %-12s %12s %12s %12s'        

        return myformat( hfmt, s, lowstr, lownum, highstr, highnum )
        


class Fit(NoNewAttributesAfterInit):

    def __init__(self, data, model, stat=None, method=None, estmethod=None):
        self.data = data
        self.model = model

        if stat is None:
            stat = Chi2Gehrels()
        if method is None:
            method = LevMar()
        if estmethod is None:
            estmethod = Covariance()

        self.stat = stat
        self.method = method
        self.estmethod = estmethod
        # Confidence limit code freezes one parameter
        # at a time.  Keep a record here of which one
        # that is, in case an exception is raised and
        # this parameter needs to be thawed in the
        # exception handler.
        self.thaw_indices = ()
        iter = 0
        for current_par in self.model.pars:
            if current_par.frozen is True:
                pass
            else:
                self.thaw_indices = self.thaw_indices + (iter,)
            iter = iter + 1
        self.current_frozen = -1

        # The number of times that reminimization has occurred
        # during an attempt to compute confidence limits.  If
        # that number equals self.estmethod.maxfits, cease all
        # further attempt to reminimize.
        self.refits = 0
        self._file = None
        self._nfev = 0
        NoNewAttributesAfterInit.__init__(self)
    
    def __str__(self):
	return (('data      = %s\n' +
                 'model     = %s\n' +
                 'stat      = %s\n' +
                 'method    = %s\n' +
                 'estmethod = %s') %
                (self.data.name,
                 self.model.name,
                 type(self.stat).__name__,
                 type(self.method).__name__,
                 type(self.estmethod).__name__))


    def guess(self, **kwargs):
        """
        kwargs = { 'limits' : True, 'values' : False }
        """
        self.model.guess(*self.data.to_guess(), **kwargs)


    def calc_stat(self):
        dep, staterror, syserror = self.data.to_fit(self.stat.calc_staterror)
        model = self.data.eval_model_to_fit(self.model)
        return self.stat.calc_stat(dep, model, staterror, syserror)[0]

    def calc_chisqr(self):
        if not isinstance(self.stat, Chi2):
            return None
        
        dep, staterror, syserror = self.data.to_fit(self.stat.calc_staterror)
        model = self.data.eval_model_to_fit(self.model)
        stat = self.stat.calc_stat(dep, model, staterror, syserror)[1]
        return stat*stat

    # SIGINT (i.e., typing ctrl-C) can dump the user to the Unix prompt,
    # when signal is sent from G95 compiled code.  What we want is to
    # get to the Sherpa prompt instead.  Typically the user only thinks
    # to interrupt during long fits or projection, so look for SIGINT
    # here, and if it happens, raise the KeyboardInterrupt exception
    # instead of aborting.
    def _sig_handler(self, signum, frame):
        raise KeyboardInterrupt()
    
    def _get_callback(self, outfile=None, clobber=False):
        if len(self.model.thawedpars) == 0:
            #raise FitError('model has no thawed parameters')
            raise FitErr( 'nothawedpar' )

        signal.signal(signal.SIGINT, self._sig_handler)
        
        dep, staterror, syserror = self.data.to_fit(self.stat.calc_staterror)

        self._nfev = 0
        if outfile is not None:
            if os.path.isfile(outfile) and not clobber:
                #raise FitError("'%s' exists, and clobber==False" % outfile)
                raise FitErr( 'noclobererr', outfile )
            self._file = file(outfile, 'w')
            names = ['# nfev statistic']
            names.extend(['%s' % par.fullname for par in self.model.pars
                          if not par.frozen])
            print >> self._file, ' '.join(names)

        def cb(pars):
            # We need to store the new parameter values in order to support
            # linked parameters

            self.model.thawedpars = pars
            model = self.data.eval_model_to_fit(self.model)
            stat = self.stat.calc_stat(dep, model, staterror, syserror)

            if self._file is not None:
                vals = ['%5e %5e' % (self._nfev, stat[0])]
                vals.extend(['%5e' % val for val in self.model.thawedpars])
                print >> self._file, ' '.join(vals)

            self._nfev+=1
            return stat

        return cb

    def fit(self, outfile=None, clobber=False):
        dep, staterror, syserror = self.data.to_fit(self.stat.calc_staterror)
        if not iterable(dep) or len(dep) == 0:
            #raise FitError('no noticed bins found in data set')
            raise FitErr( 'nobins' )

        if ((iterable(staterror) and 0.0 in staterror) and
            isinstance(self.stat, Chi2) and
            type(self.stat) != Chi2 and
            type(self.stat) != Chi2ModVar):
            #raise FitError('zeros found in uncertainties, consider using' +
            #               ' calculated uncertainties')
            raise FitErr( 'binhas0' )

        if (getattr(self.data, 'subtracted', False) and
            isinstance(self.stat, Likelihood) ):
            #raise FitError('%s statistics cannot be used with background'
            #               % self.stat.name + ' subtracted data')
            raise FitErr( 'statnotforbackgsub', self.stat.name )


        init_stat = self.calc_stat()
        output = self.method.fit(self._get_callback(outfile, clobber),
                                 self.model.thawedpars,
                                 self.model.thawedparmins,
                                 self.model.thawedparmaxes)
        # LevMar always calculate chisquare, so call calc_stat
        # just in case statistics is something other then chisquare
        self.model.thawedpars = output[1]
        if self._file is not None:
            self._file.close()
            self._file=None

        tmp = list(output)
        tmp[2] = self.calc_stat()
        output = tuple(tmp)
        # end of the gymnastics 'cause one cannot write to a tuple
        return FitResults(self, output, init_stat)

    def simulfit(self, *others):
        if len(others) == 0:
            return self.fit()

        fits = (self,) + others
        d = DataSimulFit('simulfit data', tuple(f.data for f in fits))
        m = SimulFitModel('simulfit model', tuple(f.model for f in fits))

        f = Fit(d, m, self.stat, self.method)
        return f.fit()

    def est_errors(self, methoddict=None, parlist=None):
        # Define functions to freeze and thaw a parameter before
        # we call fit function -- projection can call fit several
        # times, for each parameter -- that parameter must be frozen
        # while the others freely vary.        
        def freeze_par(pars, parmins, parmaxes, i):
            # Freeze the indicated parameter; return
            # its place in the list of all parameters,
            # and the current values of the parameters,
            # and the hard mins amd maxs of the parameters
            self.model.pars[self.thaw_indices[i]].val = pars[i]
            self.model.pars[self.thaw_indices[i]].frozen = True
            self.current_frozen = self.thaw_indices[i]
            keep_pars = ones_like(pars)
            keep_pars[i] = 0
            current_pars = pars[where(keep_pars)]
            current_parmins = parmins[where(keep_pars)]
            current_parmaxes = parmaxes[where(keep_pars)]
            return (current_pars, current_parmins, current_parmaxes)

        def thaw_par(i):
            if (i < 0):
                pass
            else:
                self.model.pars[self.thaw_indices[i]].frozen = False
                self.current_frozen = -1

        # confidence needs to know which parameter it is working on.
        def get_par_name( ii ):
            return self.model.pars[self.thaw_indices[ii]].fullname
        
        # Call from a parameter estimation method, to report
        # that limits for a given parameter have been found
        def report_progress(i, lower, upper):
            if (i < 0):
                pass
            else:
                name = self.model.pars[self.thaw_indices[i]].fullname
                if isnan(lower) or isinf(lower):
                    info("%s \tlower bound: -----" % name)
                else:
                    info("%s \tlower bound: %g" % (name, lower))
                if isnan(upper) or isinf(upper):
                    info("%s \tupper bound: -----" % name)
                else:
                    info("%s \tupper bound: %g" % (name, upper))


        # If starting fit statistic is chi-squared or C-stat,
        # can calculate reduced fit statistic -- if it is
        # more than 3, don't bother calling method to estimate
        # parameter limits.

        if (type(self.stat) is LeastSq):
            #raise FitError('cannot estimate confidence limits with ' +
            #               type(self.stat).__name__)
            raise EstErr( 'noerr4least2', type(self.stat).__name__)

        
        if (type(self.stat) is not Cash):
            dep, staterror, syserror = self.data.to_fit(
                self.stat.calc_staterror)

            if not iterable(dep) or len(dep) == 0:
                #raise FitError('no noticed bins found in data set')
                raise FitErr( 'nobins' )

            # For chi-squared and C-stat, reduced statistic is
            # statistic value divided by number of degrees of
            # freedom.

            # Degress of freedom are number of data bins included
            # in fit, minus the number of thawed parameters.
            dof = len(dep) - len(self.model.thawedpars)
            if (dof < 1):
                #raise FitError('degrees of freedom are zero or lower')
                raise EstErr( 'nodegfreedom' )
            
            if (hasattr(self.estmethod, "max_rstat") and
                (self.calc_stat() / dof) > self.estmethod.max_rstat):
                #raise FitError('reduced statistic larger than ' +
                #               str(self.estmethod.max_rstat))
                raise EstErr( 'rstat>max', str(self.estmethod.max_rstat) )

        # If statistic is chi-squared, change fitting method to
        # Levenberg-Marquardt; else, switch to NelderMead.  (We
        # will do fitting during projection, and therefore don't
        # want to use LM with a stat other than chi-squared).

        # If current method is not LM or NM, warn it is not a good
        # method for estimating parameter limits.
        if (type(self.estmethod) is not Covariance and
            type(self.method) is not NelderMead and
            type(self.method) is not LevMar):
            warning(self.method.name + " is inappropriate for confidence " +
                    "limit estimation")
        
        oldmethod = self.method
        if (hasattr(self.estmethod, "fast") and
            bool_cast(self.estmethod.fast) is True and
            methoddict is not None):
            if (isinstance(self.stat, Likelihood) ):
                if (type(self.method) is not NelderMead):
                    self.method = methoddict['neldermead']
                    warning("Setting optimization to " + self.method.name
                            + " for confidence limit search")
            else:
                if (type(self.method) is not LevMar):
                    self.method = methoddict['levmar']
                    warning("Setting optimization to " + self.method.name
                            + " for confidence limit search")

        # Now, set up before we call the confidence limit function
        # Keep track of starting values, will need to set parameters
        # back to starting values when we are done.
        startpars = self.model.thawedpars
        startsoftmins = self.model.thawedparmins
        startsoftmaxs = self.model.thawedparmaxes
        starthardmins = self.model.thawedparhardmins
        starthardmaxs = self.model.thawedparhardmaxes

        # If restricted to soft_limits, only send soft limits to
        # method, and do not reset model limits
        if (bool_cast(self.estmethod.soft_limits) is True):
            starthardmins = self.model.thawedparmins
            starthardmaxs = self.model.thawedparmaxes
        else:
            self.model.thawedparmins = starthardmins
            self.model.thawedparmaxes = starthardmaxs
        
        self.current_frozen = -1

        # parnums is the list of indices of the thawed parameters
        # we want to visit.  For example, if there are three thawed
        # parameters, and we want to derive limits for only the first
        # and third, then parnums = [0,2].  We construct the list by
        # comparing each parameter in parlist to the thawed model
        # parameters.  (In the default case, when parlist is None,
        # that means get limits for all thawed parameters, so parnums
        # is [0, ... , numpars - 1], if the number of thawed parameters
        # is numpars.)
        parnums = []
        if parlist is not None:
            allpars = [p for p in self.model.pars if not p.frozen]
            for p in parlist:
                count = 0
                match = False
                for par in allpars:
                    if p is par:
                        parnums.append(count)
                        match = True
                    count = count + 1
                if (match == False):
                    raise EstErr('noparameter', p.fullname)
            parnums = array(parnums)
        else:
            parlist = [p for p in self.model.pars if not p.frozen]
            parnums = arange(len(startpars))
            
        # If we are here, we are ready to try to derive confidence limits.
        # General rule:  if failure because a hard limit was hit, find
        # out which parameter it was so we can tell the user.
        # If a new minimum statistic was found, start over, with parameter
        # values that yielded new lower statistic as the new starting point.
        output = None
        results = None
        oldremin = -1.0
        if (hasattr(self.estmethod, "remin")):
            oldremin = self.estmethod.remin
        try:
            output = self.estmethod.compute(self._get_callback(),
                                            self.method.fit,
                                            self.model.thawedpars,
                                            startsoftmins,
                                            startsoftmaxs,
                                            starthardmins,
                                            starthardmaxs,
                                            parnums,
                                            freeze_par, thaw_par,
                                            report_progress, get_par_name)
        except EstNewMin, e:
            # If maximum number of refits has occurred, don't
            # try to reminimize again.
            if (hasattr(self.estmethod, "maxfits") and
                not(self.refits < self.estmethod.maxfits-1)):
                self.refits = 0
                thaw_par(self.current_frozen)
                self.model.thawedpars = startpars
                self.model.thawedparmins = startsoftmins
                self.model.thawedparmaxes = startsoftmaxs
                self.method = oldmethod
                if (hasattr(self.estmethod, "remin")):
                    self.estmethod.remin = -1.0
                warning("Maximum number of reminimizations reached")
            
            # First report results of new fit, then call
            # compute limits for those new best-fit parameters
            for p in parlist:
                p.frozen = False
            self.current_frozen = -1

            if e.args != ():
                self.model.thawedpars = e.args[0]

            self.model.thawedparmins = startsoftmins
            self.model.thawedparmaxes = startsoftmaxs
            results = self.fit()
            self.refits = self.refits + 1
            warning("New minimum statistic found while computing confidence limits")
            warning("New best-fit parameters:\n" + results.format())

            # Now, recompute errors for new best-fit parameters
            results = self.est_errors(methoddict, parlist)
            self.model.thawedparmins = startsoftmins
            self.model.thawedparmaxes = startsoftmaxs
            self.method = oldmethod
            if (hasattr(self.estmethod, "remin")):
                self.estmethod.remin = oldremin
            return results
        except:
            for p in parlist:
                p.frozen = False
            self.current_frozen = -1
            self.model.thawedpars = startpars
            self.model.thawedparmins = startsoftmins
            self.model.thawedparmaxes = startsoftmaxs
            self.method = oldmethod
            if (hasattr(self.estmethod, "remin")):
                self.estmethod.remin = oldremin
            raise

        for p in parlist:
            p.frozen = False
        self.current_frozen = -1
        self.model.thawedpars = startpars
        self.model.thawedparmins = startsoftmins
        self.model.thawedparmaxes = startsoftmaxs
        results = ErrorEstResults(self, output, parlist)
        self.method = oldmethod
        if (hasattr(self.estmethod, "remin")):
            self.estmethod.remin = oldremin
        
        return results
