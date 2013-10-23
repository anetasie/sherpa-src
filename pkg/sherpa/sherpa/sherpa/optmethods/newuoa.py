import numpy
import _newuoa


# I know the following util stuff lives in optfct.py but I doubt if this
# file will be around for long, so this is a quick and dirty way of doing
# it for now.
EPSILON = numpy.float_(numpy.finfo(numpy.float32).eps)

def _check_args(x0, xmin, xmax):
    x = numpy.array(x0, numpy.float_)  # Make a copy
    xmin = numpy.asarray(xmin, numpy.float_)
    xmax = numpy.asarray(xmax, numpy.float_)    
    if (x.shape != xmin.shape) or (x.shape != xmax.shape):
        raise TypeError('input array sizes do not match')
    _move_within_limits(x, xmin, xmax)
    return x, xmin, xmax

def _move_within_limits(x, xmin, xmax):
    below = numpy.flatnonzero(x < xmin)
    if below.size > 0:
        x[below] = xmin[below]

    above = numpy.flatnonzero(x > xmax)
    if above.size > 0:
        x[above] = xmax[above]

def _outside_limits(x, xmin, xmax):
    return (numpy.any(x < xmin) or numpy.any(x > xmax))

def my_is_nan( x ):
    fubar = filter( lambda xx: xx != xx or xx is numpy.nan or numpy.isnan(xx) and numpy.isfinite(xx), x )
    if len( fubar ) > 0:
        return True
    else:
        return False

def newuoa( fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None, verbose=-1,
            rhobeg=None, rhoend=None ):
    """    
C
C     This subroutine seeks the least value of a function of many variables,
C     by a trust region method that forms quadratic models by interpolation.
C     There can be some freedom in the interpolation conditions, which is
C     taken up by minimizing the Frobenius norm of the change to the second
C     derivative of the quadratic model, beginning with a zero matrix. The
C     arguments of the subroutine are as follows.
C
C     N must be set to the number of variables and must be at least two.
C     NPT is the number of interpolation conditions. Its value must be in the
C       interval [N+2,(N+1)(N+2)/2].
C     Initial values of the variables must be set in X(1),X(2),...,X(N). They
C       will be changed to the values that give the least calculated F.
C     RHOBEG and RHOEND must be set to the initial and final values of a trust
C       region radius, so both must be positive with RHOEND<=RHOBEG. Typically
C       RHOBEG should be about one tenth of the greatest expected change to a
C       variable, and RHOEND should indicate the accuracy that is required in
C       the final values of the variables.
C     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
C       amount of printing. Specifically, there is no output if IPRINT=0 and
C       there is output only at the return if IPRINT=1. Otherwise, each new
C       value of RHO is printed, with the best vector of variables so far and
C       the corresponding value of the objective function. Further, each new
C       value of F with its variables are output if IPRINT=3.
C     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
C     The array W will be used for working space. Its length must be at least
C     (NPT+13)*(NPT+N)+3*N*(N+3)/2.
C
C     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
C     the value of the objective function for the variables X(1),X(2),...,X(N).
C
C     Partition the working space array, so that different parts of it can be
C     treated separately by the subroutine that performs the main calculation.
C
"""
    x, xmin, xmax = _check_args( x0, xmin, xmax )

    if rhobeg is None:
        rhobeg = 10.0
    if rhoend is None:
        rhoend = ftol
    if maxfev is None:
        maxfev = 1024 * len(x)

    orig_fcn = fcn
    def fcn(x_new):
        if _outside_limits(x_new, xmin, xmax) or my_is_nan(x_new):        
            return 1.0e128
        return orig_fcn(x_new)

    ierr = 0
    if len(x) >= 2:
        nfev, fval = _newuoa.newuoa( x, rhobeg, rhoend, verbose,
                                     maxfev, fcn )
    else:
        fval = fcn(x)
        nfev = 1
        ierr = 2
        

    if nfev >= maxfev: ierr = 1

    key = {
        0: (True, 'successful termination'),
        1: (False,
            ('number of function evaluations has exceeded maxfev=%d' % maxfev)),
        2: (False,
            ('the number of free par (%d) must be >= 2' % len(x)))
        } 

    status, msg = key.get(ierr, (False, 'unknown status flag (%d)' % ierr))
    rv = (status, x, fval)
    rv += (msg, {'info': ierr, 'nfev': nfev})
    return rv


