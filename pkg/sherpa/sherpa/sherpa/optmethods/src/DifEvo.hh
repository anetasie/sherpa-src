#ifndef DifEvo_hh
#define DifEvo_hh

//
// An implementation of the Differential Evolution for continous
// function optimization, an algorithm by Kenneth Price and Rainer Storn.
// See: http://www.icsi.berkeley.edu/~storn/code.html
//
//  /***************************************************************
//  **                                                            **
//  **        D I F F E R E N T I A L     E V O L U T I O N       **
//  **                                                            **
//  ** Program: de.c                                              **
//  ** Version: 3.6                                               **
//  **                                                            **
//  ** Authors: Dr. Rainer Storn                                  **
//  **          c/o ICSI, 1947 Center Street, Suite 600           **
//  **          Berkeley, CA 94707                                **
//  **          Tel.:   510-642-4274 (extension 192)              **
//  **          Fax.:   510-643-7684                              **
//  **          E-mail: storn@icsi.berkeley.edu                   **
//  **          WWW: http://http.icsi.berkeley.edu/~storn/        **
//  **          on leave from                                     **
//  **          Siemens AG, ZFE T SN 2, Otto-Hahn Ring 6          **
//  **          D-81739 Muenchen, Germany                         **
//  **          Tel:    636-40502                                 **
//  **          Fax:    636-44577                                 **
//  **          E-mail: rainer.storn@zfe.siemens.de               **
//  **                                                            **
//  **          Kenneth Price                                     **
//  **          836 Owl Circle                                    **
//  **          Vacaville, CA 95687                               **
//  **          E-mail: kprice@solano.community.net               ** 
//  **                                                            **
//  ** This program implements some variants of Differential      **
//  ** Evolution (DE) as described in part in the techreport      **
//  ** tr-95-012.ps of ICSI. You can get this report either via   **
//  ** ftp.icsi.berkeley.edu/pub/techreports/1995/tr-95-012.ps.Z  **
//  ** or via WWW: http://http.icsi.berkeley.edu/~storn/litera.html*
//  ** A more extended version of tr-95-012.ps is submitted for   **
//  ** publication in the Journal Evolutionary Computation.       ** 
//  **                                                            **
//  ** You may use this program for any purpose, give it to any   **
//  ** person or change it according to your needs as long as you **
//  ** are referring to Rainer Storn and Ken Price as the origi-  **
//  ** nators of the the DE idea.                                 **
//  ** If you have questions concerning DE feel free to contact   **
//  ** us. We also will be happy to know about your experiences   **
//  ** with DE and your suggestions of improvement.               **
//  **                                                            **
//  ***************************************************************/
//
// Differential Evolution Solver Class
// Based on algorithms developed by Dr. Rainer Storn & Kenneth Price
// Written By: Lester E. Godwin
//             PushCorp, Inc.
//             Dallas, Texas
//             972-840-0208 x102
//             godwin@pushcorp.com
// Created: 6/8/98
// Last Modified: 6/8/98
// Revision: 1.0
//

#include "RanOpt.hh"
#include "Simplex.hh"

namespace sherpa {

  template < typename Func, typename Data, typename Algo >
  class DifEvo : public sherpa::RanOpt< Func, Data > {

    typedef void (DifEvo<Func,Data,Algo>::*StrategyFuncPtr)( int,
							    sherpa::Simplex&,
							    const double* model_par,
							    std::vector<double>& );

  public:

    enum Strategy { Best1Exp, Rand1Exp, RandToBest1Exp, Best2Exp, Rand2Exp,
		    Best1Bin, Rand1Bin, RandToBest1Bin, Best2Bin, Rand2Bin };

    DifEvo( int numpar, double* par, const double* lo, const double* hi,
	    Func func, Data xtra, double xprob, double scale, int strategy,
	    int seed, int mfct=0 )
      : RanOpt< Func, Data >( numpar, par, lo, hi, func, xtra, seed ), 
	cross_over_probability( xprob ), scale_factor( scale ),
	strategy_func_ptr( 0 ), local_opt( numpar, par, lo, hi, func, xtra,
					   mfct ) {
      switch ( strategy  ) {
      case Best1Exp:
	// strategy DE0 not in the paper
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::best1exp;
	break;
      case Rand1Exp:
	// strategy DE1 in the techreport
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::rand1exp;
	break;
      case RandToBest1Exp:
	// similiar to DE2 but generally better
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::randtobest1exp;
	break;
      case Best2Exp:
	// is another powerful strategy worth trying
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::best2exp;
	break;
      case Rand2Exp:
	// seems to be a robust optimizer for many functions
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::rand2exp;
	break;
      case Best1Bin:
	// Essentially same strategies but BINOMIAL CROSSOVER
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::best1bin;
	break;
      case Rand1Bin:
	// Essentially same strategies but BINOMIAL CROSSOVER
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::rand1bin;
	break;
      case RandToBest1Bin:
	// Essentially same strategies but BINOMIAL CROSSOVER
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::randtobest1bin;
	break;
      case Best2Bin:
	// Essentially same strategies but BINOMIAL CROSSOVER
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::best2bin;
	break;
      case Rand2Bin:
	// Essentially same strategies but BINOMIAL CROSSOVER
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::rand2bin;
	break;
      default:
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::best1exp;
	break;
      }

      //
      // Choice of strategy
      // We have tried to come up with a sensible naming-convention: DE/x/y/z
      // DE :  stands for Differential Evolution
      // x  :  a string which denotes the vector to be perturbed
      // y  :  number of difference vectors taken for perturbation of x
      // z  :  crossover method (exp = exponential, bin = binomial)
      //
      // There are some simple rules which are worth following:
      // 1)  F is usually between 0.5 and 1 (in rare cases > 1)
      // 2)  CR is between 0 and 1 with 0., 0.3, 0.7 and 1. being worth to be 
      //     tried first
      // 3)  To start off NP = 10*D is a reasonable choice. Increase NP if 
      //     misconvergence happens.                                       
      // 4)  If you increase NP, F usually has to be decreased
      // 5)  When the DE/best... schemes fail DE/rand... usually works and
      //     vice versa
      //
   
      return;

    }
	
    int operator( )( double* model_par, int verbose, int maxnfev, double tol, 
		     int population_size, int& nfev, double& fstat ) {

      nfev = 0;
      int ierr = EXIT_SUCCESS;
      fstat = std::numeric_limits< double >::max( );
      population_size = std::abs( population_size );

      //
      // For each row of the 2d-array population and children:
      // the columns [ 0, npar - 1 ] contain the parameters, and
      // the column npar contains the function values:
      // (*usrfunc)( population(ii,0), ... population(ii,npar-1) ) =
      //                                                     population(ii,npar);
      //
      // The array shall have dimension: population( population_size, npar + 1 )
      //
      // Will use the class Simplex since it has all the the infrastructure
      // that is needed to check for convergence although it is not a simplex
      // in the classic sense.
      //
      sherpa::Simplex population( population_size, NPAR + 1 );
      for ( int ii = 0; ii < population_size; ++ii ) {
	for ( int jj = 0; jj < NPAR; ++jj )
	  population[ ii ][ jj ] =
	    random_number_ee( Opt<Func,Data>::lo_bound[ jj ],
			      Opt<Func,Data>::hi_bound[ jj ] );
	population[ ii ][ NPAR ] =
	  std::numeric_limits< double >::max( );
      }

      //
      // allocate an extra element to store the function value
      //
      std::vector< double > trial_solution( NPAR + 1 );
      std::vector< double > fct_vals( population_size );
      double tol_sqr = tol * tol;
      int simplex_tst = 0;

      ierr = local_opt.minimize( model_par, tol, maxnfev - nfev, nfev, fstat );
      if ( EXIT_SUCCESS != ierr )
	return ierr;

      for ( ; nfev < maxnfev; ) {

	for ( int candidate=0; candidate < population_size && nfev < maxnfev;
	      ++candidate ) {

	  population.copy_row( candidate, trial_solution );
	  (this->*strategy_func_ptr)( candidate, population, model_par,
				      trial_solution );

	  ierr = local_opt.eval_user_func( maxnfev, &trial_solution[0],
					   trial_solution[ NPAR ],
					   nfev, ierr );

	  if ( trial_solution[ NPAR ] <
	       population[ candidate ][ NPAR ] ) {
	    population.copy_row( trial_solution, candidate );

	    if ( trial_solution[ NPAR ] < fstat ) {

	      ierr =
		local_opt.minimize( &trial_solution[ 0 ], tol, maxnfev - nfev,
				    nfev, trial_solution[ NPAR ] );
	      if ( EXIT_SUCCESS != ierr )
		return ierr;

	      fstat = trial_solution[ NPAR ];
	      update_par( trial_solution, model_par );

	    }  // if ( trial_solution[ NPAR ] < fstat ) {

	    population.sort( );
	    if ( population.check_convergence( tol, tol_sqr, fct_vals,
					       simplex_tst ) )
	      return EXIT_SUCCESS;

	  } // if ( trial_solution[ NPAR ] < population( ...

	}  // for ( int candidate=0; candidate < population_size && 

      }                                            // for ( ; nfev < maxnfev; )
      return ierr;

    }


  private:

    double cross_over_probability, scale_factor;
    StrategyFuncPtr strategy_func_ptr;
    Algo local_opt;

    //
    // EXPONENTIAL CROSSOVER
    //
    // DE/best/1/exp
    // Our oldest strategy but still not bad. However, we have found several
    // optimization problems where misconvergence occurs.
    //    
    void best1exp( int candidate, sherpa::Simplex& population,
		   const double* model_par,
		   std::vector< double >& trial_solution )  {

      int r1, r2;
      select_samples( candidate, population.get_nrows( ), &r1, &r2 );

      int n = RanOpt<Func,Data>::random_number( 0, NPAR - 1 );
      for ( int ii = 0;
	    RanOpt<Func,Data>::random_number( ) < cross_over_probability &&
	      ii < NPAR; ++ii ) {
	trial_solution[ n ] = model_par[ n ] +
	  scale_factor * ( population[ r1 ][ n ] - population[ r2 ][ n ] );
	n = ( n + 1 ) % NPAR;
      }

      return;
    }

    //
    // DE/rand/1/exp
    // This is one of my favourite strategies. It works especially well when
    // the "bestit[]"-schemes experience misconvergence. Try e.g.
    // F=0.7 and CR=0.5 as a first guess.
    //
    void rand1exp( int candidate, sherpa::Simplex& population,
		   const double* model_par,
		   std::vector< double >& trial_solution ) {

      int r1, r2, r3;
      select_samples( candidate, population.get_nrows(), &r1, &r2, &r3 );
      
      int n = random_number( 0, NPAR - 1 );
      for ( int ii = 0;
	    RanOpt<Func,Data>::random_number( ) < cross_over_probability &&
	      ii < NPAR; ++ii ) {
	trial_solution[ n ] = population[ r1 ][ n ] +
	  + scale_factor * ( population[ r2 ][ n ] - population[ r3 ][ n ] );
	n = (n + 1) % NPAR;
      }
      
      return;
      
    }

    //
    // DE/rand-to-best/1/exp
    // This strategy seems to be one of the best strategies. Try F=0.85 and
    // CR=1. If you get misconvergence try to increase NP. If this doesn't
    // help you should play around with all three control variables.
    //
    void randtobest1exp( int candidate, sherpa::Simplex& population,
			 const double* model_par,
			 std::vector< double >& trial_solution ) {

      int r1, r2;  
      select_samples( candidate, population.get_nrows(), &r1, &r2 );

      int n = random_number( 0, NPAR - 1 );
      for ( int ii = 0;
	    RanOpt<Func,Data>::random_number( ) < cross_over_probability &&
	      ii < NPAR; ++ii ) {
	trial_solution[n] += scale_factor * ( model_par[ n ] -
					      trial_solution[ n ] ) +
	  scale_factor * ( population[ r1 ][ n ] - population[ r2 ][ n ] );
	n = (n + 1) % NPAR;
      }

      return;

    }

    void best2exp( int candidate, sherpa::Simplex& population,
		   const double* model_par,
		   std::vector< double >& trial_solution ) {

      int r1, r2, r3, r4;
      select_samples( candidate, population.get_nrows(), &r1, &r2, &r3, &r4 );

      int n = random_number( 0, NPAR - 1 );
      for ( int ii = 0;
	    RanOpt<Func,Data>::random_number( ) < cross_over_probability &&
	      ii < NPAR; ++ii ) {
	trial_solution[n] = model_par[ n ] + 
	  scale_factor * ( population[ r1 ][ n ] + population[ r2 ][ n ] -
			   - population[ r3 ][ n ] - population[ r4 ][ n ] );
	n = (n + 1) % NPAR;
      }

      return;

    }

    void rand2exp( int candidate, sherpa::Simplex& population,
		   const double* model_par,
		   std::vector< double >& trial_solution ) {

      int r1, r2, r3, r4, r5;
      select_samples( candidate, population.get_nrows(),
		      &r1, &r2, &r3, &r4, &r5 );

      int n = random_number( 0, NPAR - 1 );
      for ( int ii = 0;
	    RanOpt<Func,Data>::random_number( ) < cross_over_probability && 
	      ii < NPAR; ++ii ) {
	trial_solution[n] = population[ r1 ][ n ] +
	  scale_factor * ( population[ r2 ][ n ] + population[ r3 ][ n ] -
			   population[ r4 ][ n ] - population[ r5 ][ n ] );
	n = (n + 1) % NPAR;
      }

      return;

    }

    void best1bin( int candidate, sherpa::Simplex& population,
		   const double* model_par,
		   std::vector< double >& trial_solution ) {

      int r1, r2;
      select_samples( candidate, population.get_nrows(), &r1, &r2 );

      int n = random_number( 0, NPAR - 1 );  
      for ( int ii = 0; ii < NPAR; ++ii ) {
	if ( RanOpt<Func,Data>::random_number( ) < cross_over_probability ||
	     NPAR - 1 == ii )
	  trial_solution[n] = model_par[ n ] +
	    scale_factor * ( population[ r1 ][ n ] - population[ r2 ][ n ] );
	n = (n + 1) % NPAR;
      }

      return;

    }

    void rand1bin( int candidate, sherpa::Simplex& population,
		   const double* model_par,
		   std::vector< double >& trial_solution ) {

      int r1, r2, r3;
      select_samples( candidate, population.get_nrows(), &r1, &r2, &r3 );

      int n = random_number( 0, NPAR - 1 );  
      for ( int ii = 0; ii < NPAR; ++ii ) {
	if ( RanOpt<Func,Data>::random_number( ) < cross_over_probability ||
	     NPAR - 1 == ii )
	  trial_solution[n] = population[ r1 ][ n ] +
	    scale_factor * ( population[ r2 ][ n ] - population[ r3 ][ n ] );
	n = (n + 1) % NPAR;
      }

      return;

    }

    void randtobest1bin( int candidate, sherpa::Simplex& population,
			 const double* model_par,
			 std::vector< double >& trial_solution ) {

      int r1, r2;
      select_samples( candidate, population.get_nrows(), &r1, &r2 );

      int n = random_number( 0, NPAR - 1 );  
      for ( int ii = 0; ii < NPAR; ++ii ) {
	if ( RanOpt<Func,Data>::random_number( ) < cross_over_probability || NPAR - 1 == ii )
	  trial_solution[n] += scale_factor * (model_par[ n ] - trial_solution[ n ] ) +
	    scale_factor * ( population[ r1 ][ n ] - population[ r2 ][ n ] );
	n = (n + 1) % NPAR;
      }
  
      return;

    }

    void best2bin( int candidate, sherpa::Simplex& population,
		   const double* model_par,
		   std::vector< double >& trial_solution ) {

      int r1, r2, r3, r4;
      select_samples( candidate, population.get_nrows(), &r1, &r2, &r3, &r4 );

      int n = random_number( 0, NPAR - 1 );
      for ( int ii = 0; ii < NPAR; ++ii ) {
	if ( RanOpt<Func,Data>::random_number( ) < cross_over_probability || NPAR - 1 == ii )
	  trial_solution[n] = model_par[ n ] +
	    scale_factor * ( population[ r1 ][ n ] + population[ r2 ][ n ] -
			     population[ r3 ][ n ] - population[ r4 ][ n ] );
	n = (n + 1) % NPAR;
      }

      return;

    }

    void rand2bin( int candidate, sherpa::Simplex& population,
		   const double* model_par,
		   std::vector< double >& trial_solution ) {

      int r1, r2, r3, r4, r5;
      select_samples( candidate, population.get_nrows(), &r1, &r2, &r3, &r4, &r5 );

      int n = random_number( 0, NPAR - 1 );
      for ( int ii = 0; ii < NPAR; ++ii ) {
	// perform NPAR binomial trials
	if ( RanOpt<Func,Data>::random_number( ) < cross_over_probability || NPAR - 1 == ii )
	  trial_solution[n] = population[ r1 ][ n ] + 
	    scale_factor * ( population[ r2 ][ n ] + population[ r3 ][ n ] -
			     population[ r4 ][ n ] - population[ r5 ][ n ] );
	n = (n + 1) % NPAR;
      }

      return;

    }


    void select_samples( int candidate, int npop, int* r1, int* r2=0,
			 int* r3=0, int* r4=0, int* r5=0 ) {
      if ( r1 ) {
	do {
	  *r1 = RanOpt<Func,Data>::random_number( 0, npop - 1 );
	} while (*r1 == candidate);
      }
  
      if ( r2 ) {
	do {
	  *r2 = RanOpt<Func,Data>::random_number( 0, npop - 1 );
	} while ( (*r2 == candidate) || (*r2 == *r1) );
      }
  
      if ( r3 ) {
	do {
	  *r3 = RanOpt<Func,Data>::random_number( 0, npop - 1 );
	}
	while ( (*r3 == candidate) || (*r3 == *r2) || (*r3 == *r1) );
      }
  
      if ( r4 ) {
	do {
	  *r4 = RanOpt<Func,Data>::random_number( 0, npop - 1 );
	} while ( (*r4 == candidate) || (*r4 == *r3) || (*r4 == *r2) ||
		  (*r4 == *r1) );
      }
  
      if ( r5 ) {
	do {
	  *r5 = RanOpt<Func,Data>::random_number( 0, npop - 1 );
	} while ( (*r5 == candidate) || (*r5 == *r4) || (*r5 == *r3) ||
		  (*r5 == *r2) || (*r5 == *r1) );
      }  
  
      return;

    }


    void update_par( const std::vector< double >& from, double* to ) const {

      for ( int ii = 0; ii < NPAR; ++ii )
	to[ ii ] = from[ ii ];
      
    }

  };                                                            // class DifEvo

}                                                                  // namespace

#endif
