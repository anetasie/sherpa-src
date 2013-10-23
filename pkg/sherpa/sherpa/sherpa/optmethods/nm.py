from copy import copy
import numpy
import directsearch

EPSILON = numpy.float_(numpy.finfo(numpy.float32).eps)

class NelderMeadClassic( directsearch.DirectSearch ):

    def __init__( self, fcn, xmin, x, xmax, maxfev, step=None ):
        directsearch.DirectSearch.__init__( self, fcn, xmin, x, xmax, maxfev,
                                            step )
        
    def __call__( self, ftol=EPSILON, verbose=0, iter=0 ):
            
        try :
            
            while ( self.nfev < self.maxfev ):

                self.order_simplex( )

                if verbose:
                    print 'f%s=%f' % (self.simplex[0][:-1],self.simplex[0][-1])

                if self.check_convergence( ftol ):
                    break

                centroid = self.calc_centroid( )
                    
                reflection_pt = (1.0 + self.reflection_coef) * centroid - self.reflection_coef * self.simplex[-1]
                reflection_pt[-1] = self.eval_user_func( reflection_pt[:-1] )

                if reflection_pt[-1] <= self.simplex[0][-1]:

                    expansion_pt = self.expansion_coef * reflection_pt + (1.0-self.expansion_coef) * centroid
                    expansion_pt[-1] = self.eval_user_func( expansion_pt[:-1] )
                    if expansion_pt[-1] < self.simplex[0][-1]:
                        self.simplex[-1] = copy( expansion_pt )
                    else:
                        self.simplex[-1] = copy( reflection_pt )

                elif reflection_pt[-1] >= self.simplex[-2][-1]:
                    
                    if reflection_pt[-1] < self.simplex[-1][-1]:
                        self.simplex[-1] = copy( reflection_pt )
                    contraction_pt = self.contraction_coef * self.simplex[-1] + (1.0 - self.contraction_coef ) * centroid
                    contraction_pt[-1] = self.eval_user_func( contraction_pt[:-1] )
                    if contraction_pt[-1] < self.simplex[-1][-1]:
                        self.simplex[-1] = copy( contraction_pt )
                    else:
                        for ii in range( 1, self.npar + 1):
                            self.simplex[ ii ] = self.shrink_coef * ( self.simplex[ ii ] + self.simplex[0] )
                            self.simplex[ ii ][-1] = self.eval_user_func( self.simplex[ii][:-1] )
                else:
                    self.simplex[-1] = copy( reflection_pt )
                    
            if iter < 0:
                self.xpar = self.simplex[0][:-1]
                self.simplex = self.init_simplex( self.xpar, step=None,
                                                  initsimplex=1 )
                iter += 1
                self.__call__( ftol=ftol, verbose=verbose, iter=iter )
            return self.get_result( 0 )

        except directsearch.InputError:
            return self.get_result( 1 )            
        except directsearch.OutOfBounds:
            return self.get_result( 2 )
        except directsearch.MaxfevError:
            return self.get_result( 3 )

        
class NelderMeadSimplex( directsearch.DirectSearch ):

    def __init__( self, fcn, xmin, x, xmax, maxfev, step=None ):
        directsearch.DirectSearch.__init__( self, fcn, xmin, x, xmax, maxfev,
                                            step )

    def __call__( self, ftol=EPSILON, verbose=0, iter=0 ):

        try :
            
            rho_chi = self.reflection_coef * self.expansion_coef

            while ( self.nfev < self.maxfev ):

                self.order_simplex( )

                if verbose:
                    print 'f%s=%f' % (self.simplex[0][:-1],smallest_funcval)

                if self.check_convergence( ftol ):
                    break

                centroid = self.calc_centroid( )
            
                smallest_funcval = self.simplex[0][-1]
                secondlargest_funcval = self.simplex[-2][-1]
                largest_vertex = self.simplex[-1]
                largest_funcval = largest_vertex[-1]

                reflection_pt = self.move_vertex( centroid, self.reflection_coef )

                if smallest_funcval <= reflection_pt[-1] and \
                       reflection_pt[-1] < secondlargest_funcval:

                    self.simplex[-1] = reflection_pt[:]
                    if verbose:
                        print '\taccept reflection point'

                elif reflection_pt[-1] < smallest_funcval:

                    # calculate the expansion point
                    expansion_pt = self.move_vertex( centroid, rho_chi )
            
                    if expansion_pt[-1] < reflection_pt[-1]:
                        self.simplex[-1] = expansion_pt
                        if verbose:
                            print '\taccept expansion point'
                    else:
                        self.simplex[-1] = reflection_pt
                        if verbose:
                            print '\taccept reflection point'

                else: 
                
                    shrinkme = self.contract_in_out( verbose, centroid, reflection_pt )
                    if shrinkme:
                        self.shrink()
                        if verbose:
                            print '\tshrink'

            if iter < 0:
                self.xpar = self.simplex[0][:-1]
                self.simplex = self.init_simplex( self.xpar, step=None,
                                                  initsimplex=1 )
                iter += 1
                self.__call__( ftol=ftol, verbose=verbose, iter=iter )
            return self.get_result( 0 )

        except directsearch.InputError:
            return self.get_result( 1 )            
        except directsearch.OutOfBounds:
            return self.get_result( 2 )
        except directsearch.MaxfevError:
            return self.get_result( 3 )

def rosen(x):  # The Rosenbrock function
    x = numpy.asarray(x)
    return numpy.sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0,axis=0)

if __name__ == '__main__':
    npar = 5
    x0 = npar * [-1.2,1.0]
    xmin = npar * [-10,-10.]
    xmax = npar * [10,10]
    maxfev = npar * 1024
    nms = NelderMeadSimplex( rosen, x0, xmin, xmax, maxfev )
    print nms( )
    nmc = NelderMeadClassic( rosen, x0, xmin, xmax, maxfev )
    print nmc( )
