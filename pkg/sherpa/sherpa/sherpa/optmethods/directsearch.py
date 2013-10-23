import numpy

class OptError( Exception ):
    """Base class for exceptions in this module."""
    pass

class InputError( OptError ):
    """Exception raised for errors in the input."""
    def __init__( self, low, x, high ):
        self.low = low
        self.xpar = x
        self.high = high
    def __str__(self):
        msg = 'The input arrays: low (%s), x (%s) and high (%s) sizes must match' % ( self.low, self.x, self.high )
        return msg

class MaxfevError( OptError ):
    """Exception raised when the maximum number of function evaluation is exceeded"""
    def __init__( self, maxfev ):
        self.maxfev = maxfev
    def __str__( self ):
        msg = 'number of function evaluations has exceeded maxfev=%d' % self.maxfev
        return msg
    
class OutOfBounds( OptError ):
    """Exception raised for errors in the input (low <= x <= high)"""
    def __init__( self, low, x, high ):
        self.low = low
        self.x = x
        self.high = high
    def __str__( self ):
        msg = 'The following condition must be valid: low (%s) <= x (%s) <= high (%s)' % ( self.low, self.x, self.high )
        return msg


class Opt( object ):
    
    def __init__( self, fcn, low, x, high ):

        self.low, self.xpar, self.high = self.check_args( low, x, high )
        
        self.nfev = 0
        self.fcn = fcn
        self.FUNC_MAX = numpy.float_(numpy.finfo(numpy.float_).max) / 1.0e8

        if self.outside_limits( x ):
            raise OutOfBounds( low, x, high )

        
    def check_args( self, xmin, x, xmax ):
        x = numpy.array( x, numpy.float_ )          # Make a copy
        xmin = numpy.asarray( xmin, numpy.float_ )  # Make a copy
        xmax = numpy.asarray( xmax, numpy.float_  ) # Make a copy
        if ( x.shape != xmin.shape ) or ( x.shape != xmax.shape ):
            raise InputError( xmin, x, xmax )
        return xmin, x, xmax
    
    def outside_limits( self, x ):
        return ( numpy.any( x < self.low ) or numpy.any( x > self.high ) )

class DirectSearch( Opt ):
    
    def __init__( self, fcn, x, xmin, xmax, maxfev, step ):
        Opt.__init__( self, fcn, xmin, x, xmax )
        self.maxfev = maxfev
        if step is None:
            step = 0.4*numpy.ones(self.xpar.shape, numpy.float_,
                                  numpy.isfortran(self.xpar))
        if ( self.xpar.shape != step.shape ):
            raise InputError( xmin, self.xpar, step ) # !!!!!!!fix me!!!!!!!!
        self.reflection_coef  = 1.0          # rho
        self.expansion_coef   = 2.0          # chi
        self.contraction_coef = 0.5          # gamma
        self.shrink_coef      = 0.5          # sigma
        self.npar = len( x )
        self.simplex = self.init_simplex( self.xpar, step )

    def calc_centroid( self ):
        centroid = numpy.zeros((self.npar+1,), dtype=numpy.float_ )
        for vertex in range( self.npar ):
            for jj in range( self.npar ):
                centroid[ jj ] += self.simplex[ vertex ][ jj ] / self.npar
        return centroid

    def check_convergence( self, tol ):
        num = 2.0 * abs( self.simplex[0][-1] - self.simplex[-1][-1] )
        denom = abs( self.simplex[0][-1] ) + abs( self.simplex[-1][-1] ) + 1.0
        if ( num / denom > tol ):
            return False
        func_vals = self.get_func_vals( )
        if numpy.std( func_vals ) > tol:
            return False
        return True
    
    def contract_in_out( self, verbose, centroid, reflection_pt ):
        
        if self.simplex[-2][-1] <= reflection_pt[-1] and reflection_pt[-1] < self.simplex[-1][-1]:
            
            # Perform outside contraction
            rho_chi = self.reflection_coef * self.contraction_coef
            outside_contraction_pt = self.move_vertex( centroid, rho_chi )

            if outside_contraction_pt[-1] <= reflection_pt[-1]:
                self.simplex[-1] = outside_contraction_pt
                if verbose:
                    print '\taccept outside contraction point'
                return False
            else:
                return True
            
        elif reflection_pt[-1] >= self.simplex[-1][-1]:

            # Perform an inside contraction
            inside_contraction_pt = self.move_vertex( centroid,
                                                      - self.contraction_coef )
                    
            if inside_contraction_pt[-1] < self.simplex[-1][-1]:
                self.simplex[-1] = inside_contraction_pt
                if verbose:
                    print '\taccept inside contraction point'
                return False
            else:
                return True

        else:
            print 'something is wrong with contract_in_out'
            return True
        
    def eval_user_func( self, x ):
        if self.outside_limits( x ):
            return self.FUNC_MAX
        else:
            if None != self.fcn:
                self.nfev += 1
                if ( self.nfev < self.maxfev ):
                    return self.fcn( x )
                else:
                    raise MaxfevError( self.maxfev )
        
    def get_func_vals( self ):
        func_vals = []
        for vertex in self.simplex:
            func_vals.append( vertex[-1] )
        return func_vals

    def get_result( self, ierr ):
        
        x = self.simplex[0][:-1]
        fval = self.simplex[0][-1]
        key = {
            0: (True, 'Optimization terminated successfully'),
            1: (False, 'improper input parameters'),
            2: (False, 'improper values for x, xmin or xmax'),
            3: (False,
                'number of function evaluations has exceeded %d' % self.maxfev)
            }
        status, msg = key.get( ierr,
                               (False, 'unknown status flag (%d)' % ierr) )
            
        rv = (status, x, fval)
        rv += (msg, {'info': status, 'nfev': self.nfev})
        return rv

    def order_simplex( self ):
        func_vals = self.get_func_vals( )
        # return a list of indices goes from largest to smallest
        index_from_largest_to_smallest = numpy.argsort( func_vals )
        # re-arrange so self.simplex goes from smallest to largest fct val
        self.simplex = numpy.take( self.simplex,
                                   index_from_largest_to_smallest, 0 )

    def init_simplex( self, xguess, step, initsimplex=0 ):
        npars = len(xguess)
        simplex = numpy.zeros( (npars+1,npars+1), dtype=numpy.float_ )
        for ii in range( npars + 1 ):
            simplex[ ii ][:-1] = numpy.array( xguess, copy=True )
        if 0 == initsimplex:
            for ii in range( npars ):
                simplex[ ii + 1 ][ ii ] += step[ ii ]
        else:
            delta = 0.025
            for ii in range( npars ):
                if 0.0 == simplex[ ii + 1 ][ ii ]:
                    simplex[ ii + 1 ][ ii ] = delta
                else:
                    simplex[ ii + 1 ][ ii ] *= (1.0+delta)
        for ii in range( npars + 1 ):
            simplex[ ii ][-1] = self.eval_user_func( simplex[ ii ][:-1] )
        return simplex

    def move_vertex( self, centroid, coef ):
        vertex = ( 1.0 + coef ) * centroid - coef * self.simplex[-1]
        vertex[-1] = self.eval_user_func( vertex[ : -1 ] )
        return vertex

    def shrink( self ):
        npars_plus_1 = self.npar + 1
        for j in range( 1, npars_plus_1 ):
            self.simplex[j] = self.simplex[0] + self.shrink_coef * ( self.simplex[j] - self.simplex[0] )
            self.simplex[j][-1] = self.eval_user_func( self.simplex[j][:-1] )

if __name__ == '__main__':

    def tst_OutOfBounds( num ):
        try:
            x = num * [ 1.0 ]
            low = num * [ -10.0 ]
            high = num * [ 10.0 ]
            ok = Opt( None, low, x, high )
            notok = Opt( None, high, x,  low )
        except OutOfBounds, o:
            print o
        except OptError:
            print 'caught base class'

    def tst_arg( num ):
        try:
            x = num
            low = num * [ -10.0 ]
            high = num * [ 10.0 ]
            notok = Opt( None, x, high, low )
        except InputError, i:
            print i
        except OptError:
            print 'caught base class'

    num = 3
    tst_OutOfBounds( num )
    tst_arg( num )
