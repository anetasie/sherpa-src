#ifndef NelderMead_hh
#define NelderMead_hh

// 
//  Copyright (C) 2007  Smithsonian Astrophysical Observatory
//
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with this program; if not, write to the Free Software Foundation, Inc.,
//  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//


//
// Nelder, J.A. and Mead, R., "A Simplex Method for Function Minimization",
// Computer Journal, Vol. 7, Issue 4 (1965), 308-313.
//
// The implementation is based on the following two papers (the description
// of the algorithm is identical in the two papers):
//
// Jeffrey C. Lagarias, James A. Reeds, Margaret H. Wright, Paul E. Wright
// "Convergence Properties of the Nelder-Mead Simplex Algorithm in Low
// Dimensions", SIAM Journal on Optimization,Vol. 9, No. 1 (1998),
// pages 112-147.
// http://citeseer.ist.psu.edu/3996.html
//
// Wright, M. H. (1996) "Direct Search Methods: Once Scorned, Now Respectable"
// in Numerical Analysis 1995 (Proceedings of the 1995 Dundee Biennial
// Conference in Numerical Analysis) (D.F. Griffiths and G.A. Watson, eds.),
// 191-208, Addison Wesley Longman, Harlow, United Kingdom.
// http://citeseer.ist.psu.edu/155516.html
//
// Note: Some (ie most) of the comments within the code were taken directly
// from the M. Wright paper titled: "Direct Search Methods: Once Scorned,
// Now Respectable".  Note, the notations 
//   f , f , f    from the paper translate to following code
//    1   n   n+1
//
//  f      ==>  simplex( index[Smallest], npar )
//   1
//  f      ==>  simplex( index[NextLargest], npar )
//   n
//  f      ==>  simplex( index[Largest], npar )
//   n+1
//
// Sep 2006 Original version written by D. T. Nguyen
// Jan 2008 Removed the dependency of an obsolete multi-dimensional
//          array with one of my own. Modifid code to be called from
//          python. Did some work to avoid the 'false convergence'
//          on a few problems sent in by the users.
//

#include "DirectSearch.hh"
#include "PyWrapper.hh"

namespace sherpa {

  template< typename Func, typename Data >
  class NelderMead : public sherpa::DirectSearch< Func, Data > {

  public:
    
    NelderMead( int numpar, double* par, const double* lo, const double* hi,
		Func func, Data xtra, double shrinkcoef=0.5 )
      : DirectSearch< Func, Data >( numpar, par, lo, hi, func, xtra ),
	shrink_coef( shrinkcoef ), centroid( numpar + 1 ),
	reflection( numpar + 1 ), expansion( numpar + 1 ),
	contraction( numpar + 1 ) {

      check_shrink_coefficients( );

    }

    // de
    NelderMead( int numpar, double* par, const double* lo, const double* hi,
		Func func, Data xtra, int mfct, double shrinkcoef=0.5 )
      : DirectSearch< Func, Data >( numpar, par, lo, hi, func, xtra ),
	shrink_coef( shrinkcoef ), centroid( numpar + 1 ),
	reflection( numpar + 1 ), expansion( numpar + 1 ),
	contraction( numpar + 1 ) {

      check_shrink_coefficients( );

    }
    /*
    NelderMead( int numpar, double* par, const double* lo, const double* hi,
		Func func, Data xtra, double tol, int maxnfev, int& nfev,
		double& fstat, int& err,
		int mfct=0, double* fvec=0, double shrinkcoef=0.5 ) :
      DirectSearch< Func, Data >( numpar, par, lo, hi, func, xtra ), 
      shrink_coef( shrinkcoef ), centroid( numpar + 1 ),
      reflection( numpar + 1 ), expansion( numpar + 1 ),
      contraction( numpar + 1 ) {

      check_shrink_coefficients( );

      int nm_verbose=0, init_simplex=0;
      const int tmp[]={ 0, 1 };
      std::vector< int > final_simplex( tmp, tmp + sizeof(tmp) / sizeof(int) );
      std::vector< double > step( NPAR, 1.2 );
      
      err = this->operator( )( par, nm_verbose, init_simplex, final_simplex,
			       tol, &step[0], maxnfev, nfev, fstat );

    }
    */
    int minimize( double* model_par, double tol, int maxnfev, int& nfev,
		  double& fmin ) {

      int verbose=0, init_simplex=0;
      const int tmp[]={ 0, 1 };
      std::vector< int > final_simplex( tmp, tmp + sizeof(tmp) / sizeof(int) );
      std::vector< double > step( NPAR, 1.2 );

      return this->operator( )( model_par, verbose, init_simplex,
				final_simplex, tol, &step[0], maxnfev, nfev,
				fmin );
    }
    // de

    int operator( )( double* model_par, int verbose, int initsimplex,
		     std::vector< int >& finalsimplex,
		     double tolerance, const double* step,
		     int maxnfev, int& nfev, double& fmin,
		     double* g=NULL, double* h=NULL )  {

      try {

	nfev = 0;

	std::vector< double > funcvals( NPAR + 1 );
	int err_status=EXIT_SUCCESS;
	double tol_sqr = tolerance * tolerance;

	DirectSearch<Func,Data>::init_simplex( initsimplex, model_par, step );
	err_status =
	  DirectSearch<Func,Data>::eval_init_simplex( maxnfev, model_par,
						      nfev );
	if ( EXIT_SUCCESS != err_status )
	  return err_status;

	for ( ; nfev < maxnfev; ) {

	  // 1. Order the npar + 1 vertices to satisfy....
	  SIMPLEX.sort( );
	  fmin = SIMPLEX[ 0 ][ NPAR ];

	  // Need to have order before deciding to quit,
	  // cause the order_index[Smallest] must be known.
	  if ( SIMPLEX.check_convergence( tolerance, tol_sqr, funcvals,
					  finalsimplex[0] ) )
	    break;

	  if ( verbose ) {
	    fprintf( stdout, "nfev=%d\tfmin=%.10e\n", nfev, fmin );
	    if ( verbose > 2 )
	      SIMPLEX.print_simplex( );
	  }

	  calculate_centroid( );
		
	  // 2. Reflect.
	  reflect( verbose, maxnfev, nfev, err_status );
	  if ( EXIT_SUCCESS != err_status )
	    break;

	  //
	  // If f  <= f  < f  ,
	  //     1     r    n
	  if ( SIMPLEX[ 0 ][ NPAR ] <= reflection[ NPAR ] &&
	       reflection[ NPAR ] < SIMPLEX[ NPAR-1 ][ NPAR ] ) {

	    // accept the reflected point x and terminate iteration.
	    SIMPLEX.copy_row( reflection, NPAR );
	    if ( verbose > 1 )
	      fprintf( stdout, "\t\taccept reflection point.\n" );

	  } else if ( reflection[ NPAR ] < SIMPLEX[ 0 ][ NPAR ] ) {

	    // 3. Expand. If f  < f  ,
	    //                r    1  
	    expand( verbose, maxnfev, nfev, err_status );
	    if ( EXIT_SUCCESS != err_status )
	      break;

	  } else {

	    //
	    // 4. Contract. If f  >= f  perform a contraction between centroid
	    //                  r     n
	    // 
	    // and the better of x    and x .
	    //                    n+1      r
	    //
	    bool shrinkme = true;
	    if ( reflection[ NPAR ] >= SIMPLEX[ NPAR-1 ][ NPAR ] ) {

	      shrinkme = contract( verbose, maxnfev, nfev, err_status );
	      if ( EXIT_SUCCESS != err_status )
		break;

	    }

	    if ( shrinkme ) {
	      // 5. Perform a shrink step.
	      shrink( verbose, maxnfev, nfev, err_status );
	      if ( EXIT_SUCCESS != err_status )
		break;
	    }

	  }

	}                                        // for ( ; nfev < maxnfev; ) {
	
	for ( int ii = 0; ii < NPAR; ++ii )
	  model_par[ ii ] = SIMPLEX[ 0 ][ ii ];
  
	fmin = SIMPLEX[ 0 ][ NPAR ];

	std::vector< int >::iterator current = finalsimplex.begin() + 1;
	std::vector< int >::iterator the_end = finalsimplex.end();

	if ( current != the_end && nfev < maxnfev ) {
      
	  int mynfev=0;
	  //double old_fmin = fmin;

	  std::vector< int > myfinalsimplex( current, the_end );
	  err_status = this->operator( )( model_par, verbose, initsimplex,
					  myfinalsimplex, tolerance, step,
					  maxnfev - nfev, mynfev,
					  fmin, g, h );

	  nfev += mynfev;
	  if ( EXIT_SUCCESS != err_status )
	    return err_status;

	  /*
	    printf( "operator(iter=%d): mynfev=%d,\tnfev=%d\tfmin=%.20e\n",
	    iter, mynfev, nfev, fmin );
	    printf("\told_fmin=%.20e\n\tnew_fmin=%.20e\n", old_fmin, fmin );
	  */

	}

	if ( NULL != g )
	  for ( int ii = 0; ii <= NPAR; ++ii )
	    for ( int jj = 0; jj < NPAR; ++jj )
	      g[ ii * NPAR + jj ] = SIMPLEX[ ii ][ jj ];

	if ( NULL != h )
	  for ( int ii = 0; ii <= NPAR; ++ii )
	    h[ ii ] = SIMPLEX[ ii ][ NPAR ];

	return err_status;


      } catch( std::runtime_error& re ) {
	/*
	if ( verbose )
	  std::cerr << re.what( ) << '\n';
	*/
	throw re;
      } catch( std::exception& e ) {
	/*
	if ( verbose )
	  std::cerr << e.what( ) << '\n';
	*/
	throw e;
      }

    }                                                            // operator( )


  private:

    const double shrink_coef;

    std::vector<double> centroid, reflection, expansion, contraction;

    void calculate_centroid( )  {

      for ( int ii = 0; ii < NPAR; ++ii ) {
	centroid[ ii ] = 0.0;
	for ( int jj = 0; jj < NPAR; ++jj ) {
	  // The last vextex is to be avoided so ; jj <= NPAR; would be wrong!
	  centroid[ ii ] += SIMPLEX[ jj ][ ii ];
	}
	centroid[ ii ] /= double( NPAR );

      }
      
    }                                                     // calculate_centroid


    //
    // reflection_coef > 0, expansion_coef > 1,
    // expansion_coef > reflection_coef, 0 < contraction_coef < 1,
    // and 0 < shrinkage_coef < 1.
    //
    void check_shrink_coefficients( ) const {
      
      if ( shrink_coef <= 0.0 || shrink_coef >= 1.0 )
	throw std::runtime_error( "The shrink coefficient must be "
				  "within (0,1)" );
      
    }                                              // check_shrink_coefficients

    // 4.
    // return of true ==> perform step 5 (shrink) otherwise do not shrink 
    bool contract( int verbose, int maxnfev, int& nfev, int& err_status ) {

      sherpa::Array2d<double>::Row worst_vertex = SIMPLEX[ NPAR ];
      
      if ( SIMPLEX[ NPAR-1 ][ NPAR ] <= reflection[ NPAR ] &&
	   reflection[ NPAR ] < SIMPLEX[ NPAR ][ NPAR ] ) {

	//
	// 4a. Outside.  If f  <= f  < f    (i.e., x  is stricly better then
	//                   n     r    n+1         r
	// x    ), perform an outside contraction: calculate
	//  n+1
	//      _                 _                       _
	// x  = x + gamma * ( x - x ) = ( 1 + rho gamma ) x  - rho gamma x
	//  c                  r                                          n+1
	//
	//                                                                (2.6)
	//
	// and evaluate f  = f( x  ).
	//               c       c

	double rho_chi = DirectSearch<Func,Data>::reflection_coef *
	  DirectSearch<Func,Data>::contraction_coef;
	err_status = move_vertex( maxnfev, rho_chi, centroid, &worst_vertex[0],
		       contraction, nfev );
	if ( verbose > 1 )
	  fprintf( stdout, "\tOutside contraction\n" );
	if ( EXIT_SUCCESS != err_status )
	  return false;
    
	if ( contraction[ NPAR ] <= reflection[ NPAR ] ) {

	  // If f  <= f  , accept x  
	  //     c     r           c
	  SIMPLEX.copy_row( contraction, NPAR );
	  if ( verbose > 1 )
	    fprintf( stdout, "\t\taccept contraction point.\n" );
	  // and terminate the iteration;
	  return false;

	} else
	  // otherwise, go to step 5 (perform a shrink).
	  return true;

      } else if ( reflection[ NPAR ] >= SIMPLEX[ NPAR ][ NPAR ] ) {

	//
	// 4b. Inside. If f  >= f   , perform an inside contraction: calculate
	//                 r     n+1
	//       _           _                         _
	// x   = x - gamma ( x - x   ) = ( 1 - gamma ) x + gamma x        (2.7)
	//  cc                    n+1                             n+1
	//
	// and evaluate f   = f( x   ).
	//               cc       cc
	//

	err_status =
	  move_vertex( maxnfev, -DirectSearch<Func,Data>::contraction_coef,
		       centroid, &worst_vertex[0], contraction, nfev );
	if ( verbose > 1 )
	  fprintf( stdout, "\tInside contraction\n" );
	if ( EXIT_SUCCESS != err_status )
	  return false;

	if ( contraction[ NPAR ] < SIMPLEX[ NPAR ][ NPAR ] ) {

	  //
	  // If f   < f   , accept x  
	  //     cc    n+1          cc
	  //
	  SIMPLEX.copy_row( contraction, NPAR );
	  // and terminate the iteration;
	  if ( verbose > 1 )
	    fprintf( stdout, "\t\taccept contraction point.\n" );
	  return false;

	} else
	  // otherwise, go to step 5 (perform a shrink).
	  return true;

      } else {

	throw std::runtime_error( "ERROR: Unknown contract case\n" );

      }

      return false;

    }                                                               // contract



    // 3.
    void expand( int verbose, int maxnfev, int& nfev, int& err_status ) {

      if ( verbose > 1 )
	fprintf( stdout, "\tExpand\n" );

      //
      // calculate the expansion point x  :
      //                                e
      //      _               _                      _
      // x  = x + chi * ( x - x ) =  ( 1 + rho chi ) x  - rho chi x       (2.5)
      //  e                r                                       n+1
      //
      // and evaluate f  = f( x )
      //               e       e
      //
      double rho_chi = DirectSearch<Func,Data>::reflection_coef *
	DirectSearch<Func,Data>::expansion_coef;
      sherpa::Array2d<double>::Row worst_vertex = SIMPLEX[ NPAR ];
      err_status = move_vertex( maxnfev, rho_chi, centroid, &worst_vertex[0],
				expansion, nfev );
      if ( EXIT_SUCCESS != err_status )
	return;
  
      if ( expansion[ NPAR ] < reflection[ NPAR ] ) {

	//
	// If f  < f  , accept x  and terminate the iteration;
	//     e    r           e
	SIMPLEX.copy_row( expansion, NPAR );
	if ( verbose > 1 )
	  fprintf( stdout, "\t\taccept expansion point.\n" );

      } else {

	//
	// otherwise, (if f  >= f  ), accept x  and terminate the iteration.
	//                 e     r            r
	SIMPLEX.copy_row( reflection, NPAR );
	if ( verbose > 1 )
	  fprintf( stdout, "\t\taccept reflection point.\n" );

      }

      return;

    }                                                                // expand


    //
    // Move vertex to the position:
    //                                _
    //           x     = ( 1 + coef ) x - coef x
    //            new                           n+1
    //       _
    // where x and x    are the centroid and the worst vertex of the simplex,
    //              n+1
    //
    // respectively. Then evaluate the function at the new position.
    //
    //      _           _                        _
    // x  = x + rho * ( x - x   ) =  ( 1 + rho ) x  - rho x               (2.4)
    //  r                    n+1                           n+1
    //      _               _                      _
    // x  = x + chi * ( x - x ) =  ( 1 + rho chi ) x  - rho chi x         (2.5)
    //  e                r                                       n+1
    //      _                 _                       _
    // x  = x + gamma * ( x - x ) = ( 1 + rho gamma ) x  - rho gamma x    (2.6)
    //  c                  r                                          n+1
    //       _           _                         _
    // x   = x - gamma ( x - x   ) = ( 1 - gamma ) x + gamma x            (2.7)
    //  cc                    n+1                             n+1
    //
    int move_vertex( int maxnfev, double coef, 
		     const std::vector< double >& xbar,
		     const double* largest_vertex,
		     std::vector< double >& new_vertex,
		     int& nfev )  {

      //                                _
      //           x     = ( 1 + coef ) x - coef x
      //            new                           n+1
      //       _
      // where x and x    are the centroid and the worst vertex of the simplex,
      //              n+1
      //
      // respectively. Then evaluate the function at the new position
      //
      int err = EXIT_SUCCESS;
      double coef_plus_1 = 1.0 + coef;
      for ( int ii = 0; ii < NPAR; ++ii )
	new_vertex[ ii ] = coef_plus_1 * xbar[ ii ] -
	  coef * largest_vertex[ ii ];

      eval_user_func( maxnfev, &new_vertex[0], new_vertex[ NPAR ], nfev, err );

      return err;

    }                                                            // move_vertex


    // 1.
    void reflect( int verbose, int maxnfev, int& nfev, int& err_status ) {

      if ( verbose > 1 )
	fprintf( stdout, "\tReflect\n" );

      //
      // Compute the reflection point x  from
      //                               r
      //      _           _                        _
      // x  = x + rho * ( x - x   ) =  ( 1 + rho ) x  - rho x             (2.4)
      //  r                    n+1                           n+1
      //       _   
      // where x =  sum( x , i=1..n) / n   is the centroid of the n best 
      //                  i                points (all vertices except
      //                                   for x   )
      //                                        n+1
      //
      // Evaluate f  = f( x  )
      //           r       r
      //
      sherpa::Array2d<double>::Row worst_vertex = SIMPLEX[ NPAR ];
      err_status =
	move_vertex( maxnfev, DirectSearch<Func,Data>::reflection_coef,
		     centroid, &worst_vertex[0], reflection, nfev );

    }                                                                // reflect


    int search( double* model_par, int verbose, int initsimplex,
		int finalsimplex, double tolerance, const double* step,
		int maxnfev, int& nfev, double& fmin )  {


      try {

	std::vector< double > funcvals( NPAR + 1 );

	int err_status = EXIT_SUCCESS;
	double tol_sqr = tolerance * tolerance;

	nfev = 0;

	DirectSearch<Func,Data>::init_simplex( initsimplex, step );
	err_status =
	  DirectSearch<Func,Data>::eval_init_simplex( maxnfev, nfev );
	if ( EXIT_SUCCESS != err_status )
	  return err_status;

	for ( ; nfev < maxnfev; ) {

	  // 1. Order the NPAR + 1 vertices to satisfy....
	  SIMPLEX.sort( );
	  fmin = SIMPLEX[ 0 ][ NPAR ];

	  // Need to have order before deciding to quit,
	  // cause the order_index[Smallest] must be known.
	  if ( SIMPLEX.check_convergence( tolerance, tol_sqr, funcvals,
					  finalsimplex ) )
	    break;

	  if ( verbose ) {
	    fprintf( stdout, "nfev=%d\tfmin=%.10e\n", nfev, fmin );
	    if ( verbose > 2 )
	      SIMPLEX.print_simplex( );
	  }

	  calculate_centroid( NPAR );

	  sherpa::Array2d<double>::Row worst_vertex =
	    SIMPLEX[ NPAR ];

	  if ( verbose > 1 )
	    fprintf( stdout, "\tReflect\n" );
	  err_status =
	    move_vertex( maxnfev, DirectSearch<Func,Data>::reflection_coef,
			 centroid, &worst_vertex[0], reflection, nfev );
	  if ( EXIT_SUCCESS != err_status )
	    break;

	  if ( reflection[ NPAR ] < SIMPLEX[ 0 ][ NPAR ] ) {

	    //
	    // reflection or expand
	    //

	    if ( verbose > 1 )
	      fprintf( stdout, "\tExpand\n" );

	    double rho_chi =
	      DirectSearch<Func,Data>::reflection_coef *
	      DirectSearch<Func,Data>::expansion_coef;
	    err_status =
	      move_vertex( maxnfev, rho_chi, centroid, &worst_vertex[0],
			   expansion, nfev );
	    if ( EXIT_SUCCESS != err_status )
	      break;

	    if ( expansion[ NPAR ] < reflection[ NPAR ] ) {
	      if ( verbose > 1 )
		fprintf( stdout, "\t\taccept expansion point.\n" );
	      SIMPLEX.copy_row( expansion, NPAR );
	    } else {
	      if ( verbose > 1 )
		fprintf( stdout, "\t\taccept reflection point.\n" );
	      SIMPLEX.copy_row( reflection, NPAR );
	    }

	  } else {

	    //
	    // contract or shrink
	    //

	    int contract_shrink_me = false;

	    for ( int ii = 0; ii <= NPAR; ++ii )
	      if ( ii != NPAR )
		if ( reflection[ NPAR ] < SIMPLEX[ ii ][ NPAR ] ) {
		  SIMPLEX.copy_row( reflection, NPAR );
		  contract_shrink_me = true;
		}

	    if ( false == contract_shrink_me ) {

	      if ( reflection[ NPAR ] <= SIMPLEX[ NPAR ][ NPAR ] )
		SIMPLEX.copy_row( reflection, NPAR );

	      if ( verbose > 1 )
		fprintf( stdout, "\tInside contraction\n" );

	      // inside contraction
	      err_status =
		move_vertex( maxnfev,
			     -DirectSearch<Func,Data>::contraction_coef,
			     centroid, &worst_vertex[0], contraction, nfev );
	      if ( EXIT_SUCCESS != err_status )
		break;

	      if ( contraction[ NPAR ] > SIMPLEX[ NPAR ][ NPAR ] ) {
		shrink( verbose, maxnfev, nfev, err_status );
		if ( EXIT_SUCCESS != err_status )
		  break;

	      } else {

		if ( verbose > 1 )
		  fprintf( stdout, "\t\taccept contraction point.\n" );

		SIMPLEX.copy_row( contraction, NPAR );

	      }

	    }

	  }      

	}                                        // for ( ; nfev < maxnfev; ) {
	
	for ( int ii = 0; ii < NPAR; ++ii )
	  model_par[ ii ] = SIMPLEX[ 0 ][ ii ];
    
	fmin = SIMPLEX[ 0 ][ NPAR ];

	return err_status;

      } catch( std::runtime_error& re ) {
	throw re;
      } catch( std::exception& e ) {
	throw e;
      }

    }                                                                 // search

    // 5.
    void shrink( int verbose, int maxnfev, int& nfev, int& err_status ) {

      if ( verbose > 1 )
	fprintf( stdout, "\tShrink\n" );

      for ( int ii = 1; ii <= NPAR; ++ii ) {
	
	//
	// Evaluate f at the n points v  = x  + sigma * ( x  - x  ),
	//                             i    1              i    1
	// i = 2, ..., n+1.  The (unordered) vertices of the simplex at
	// the next iteration consist of x , v , ..., v
	//                                1   2        n+1
	for ( int jj = 0; jj < NPAR; ++jj )
	  SIMPLEX[ ii ][ jj ] = shrink_coef * SIMPLEX[ ii ][ jj ] +
	    ( 1.0 - shrink_coef ) * SIMPLEX[ 0 ][ jj ];
	sherpa::Array2d<double>::Row rowii = SIMPLEX[ ii ];
	eval_user_func( maxnfev, &rowii[0], SIMPLEX[ ii ][ NPAR ], nfev,
			err_status );
	if ( nfev >= maxnfev )
	  err_status = OptErr::MaxFev;
	if ( EXIT_SUCCESS != err_status )
	  return;

      }                          // for ( int ii = 0; ii <= this->npar; ++ii )

    }                                                            // void shrink

  };

}                                                          // namespace sherpa

#endif
