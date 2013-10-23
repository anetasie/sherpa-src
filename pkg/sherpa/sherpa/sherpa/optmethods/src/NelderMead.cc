#ifdef testNelderMead

#include <iostream>
#include <stdexcept>

#include "NelderMead.hh"
#include "fcmp.h"
#include "tstoptfct.hh"
#include "sherpa/functor.hh"

template <typename Real>
void print_pars( const char* name, int nfev, Real stat, Real answer,
		 int n, const std::vector< Real >& x,
		 Real tol=
		 1.0e4*std::sqrt( std::numeric_limits< Real >::epsilon() ) ) {

 
  std::cout << "NelderMead_" << name << '\t';
  if ( 0 == _sao_fcmp( stat, answer, std::sqrt(tol) ) )
    std::cout << nfev << '\t';
  else
    std::cout << -nfev << '\t';
  std::cout << answer << '\t';
  std::cout << stat << '\t';
  std::cout << x[0];
  for ( int ii = 1; ii < n; ++ii )
    std::cout << ',' << x[ii];
  std::cout << '\n';

}

template< typename Init, typename Fct >
void justdoit( Init init, Fct fct, int npars, 
	       std::vector< double >& pars, std::vector< double >& lo,
	       std::vector< double >& hi, std::vector< double >& step,
	       double tol, const char* header ) {

  std::vector< int > finalsimplex;

  try {

    for ( int iii = 0; iii < 2; ++iii ) {

      finalsimplex.push_back( iii );
      
      int mfcts;
      double answer;
      
      init( npars, mfcts, answer, &pars[0], &lo[0], &hi[0] );

      sherpa::NelderMead< Fct, void* > nm( npars, &pars[0], &lo[0], &hi[0],
					   fct, NULL );

      int verbose=0, initsimplex=0, maxnfev=npars*npars*1024, nfev;
      double fmin;
      nm( &pars[0], verbose, initsimplex, finalsimplex, tol, &step[0], maxnfev, nfev, fmin, NULL, NULL );
      
      print_pars( header, nfev, fmin, answer, npars, pars );
      
    }
    
  } catch( const sherpa::OptErr& oe ) {
    
    std::cerr << oe << '\n';
    
  }

  return;

}

using namespace sherpa;

void tstuncopt( int npars, double tol ) {

  int size = npars * npars * 4;
  std::vector< double > pars( size ), step( size ), lo( size ),
    hi( size );
  
  for ( int ii = 0; ii < npars; ++ii )
    step[ ii ] = 1.2;
  
    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::Rosenbrock<double,void*> );

      justdoit( fct_ptr( tstoptfct::RosenbrockInit<double> ),
		fct, npars, pars, lo, hi, step, tol, "Rosenbrock" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::FreudensteinRoth<double,void*> );

      justdoit( fct_ptr( tstoptfct::FreudensteinRothInit<double> ),
		fct, npars, pars, lo, hi, step, tol, "FreudensteinRoth" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::PowellBadlyScaled<double,void*> );

      justdoit( fct_ptr( tstoptfct::PowellBadlyScaledInit<double> ),
		fct, npars, pars, lo, hi, step, tol, "PowellBadlyScaled" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::BrownBadlyScaled<double,void*> );

      justdoit( fct_ptr( tstoptfct::BrownBadlyScaledInit<double> ),
		fct, 2, pars, lo, hi, step, tol, "BrownBadlyScaled" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::Beale<double,void*> );

      justdoit( fct_ptr( tstoptfct::BealeInit<double> ),
		fct, npars, pars, lo, hi, step, tol, "Beale" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::JennrichSampson<double,void*> );

      justdoit( fct_ptr( tstoptfct::JennrichSampsonInit<double> ),
		fct, npars, pars, lo, hi, step, tol, "JennrichSampson" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::HelicalValley<double,void*> );

      justdoit( fct_ptr( tstoptfct::HelicalValleyInit<double> ),
		fct, 3, pars, lo, hi, step, tol, "HelicalValley" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::Bard<double,void*> );

      justdoit( fct_ptr( tstoptfct::BardInit<double> ),
		fct, 3, pars, lo, hi, step, tol, "Bard" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::Gaussian<double,void*> );

      justdoit( fct_ptr( tstoptfct::GaussianInit<double> ),
		fct, 3, pars, lo, hi, step, tol, "Gaussian" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::Meyer<double,void*> );

      justdoit( fct_ptr( tstoptfct::MeyerInit<double> ),
		fct, 3, pars, lo, hi, step, tol, "Meyer" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::GulfResearchDevelopment<double,void*> );

      justdoit( fct_ptr( tstoptfct::GulfResearchDevelopmentInit<double> ),
		fct, 3, pars, lo, hi, step, tol, "GulfResearchDevelopment" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::Box3d<double,void*> );

      justdoit( fct_ptr( tstoptfct::Box3dInit<double> ),
		fct, 3, pars, lo, hi, step, tol, "Box3d" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::PowellSingular<double,void*> );

      justdoit( fct_ptr( tstoptfct::PowellSingularInit<double> ),
		fct, 4, pars, lo, hi, step, tol, "PowellSingular" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::Wood<double,void*> );

      justdoit( fct_ptr( tstoptfct::WoodInit<double> ),
		fct, 4, pars, lo, hi, step, tol, "Wood" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::KowalikOsborne<double,void*> );

      justdoit( fct_ptr( tstoptfct::KowalikOsborneInit<double> ),
		fct, 4, pars, lo, hi, step, tol, "KowalikOsborne" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::BrownDennis<double,void*> );

      justdoit( fct_ptr( tstoptfct::BrownDennisInit<double> ),
		fct, 4, pars, lo, hi, step, tol, "BrownDennis" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::Osborne1<double,void*> );

      justdoit( fct_ptr( tstoptfct::Osborne1Init<double> ),
		fct, 5, pars, lo, hi, step, tol, "Osborne1" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::Biggs<double,void*> );

      justdoit( fct_ptr( tstoptfct::BiggsInit<double> ),
		fct, 6, pars, lo, hi, step, tol, "Biggs" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::Osborne2<double,void*> );

      justdoit( fct_ptr( tstoptfct::Osborne2Init<double> ),
		fct, 11, pars, lo, hi, step, tol, "Osborne2" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::Watson<double,void*> );

      justdoit( fct_ptr( tstoptfct::WatsonInit<double> ),
		fct, 6, pars, lo, hi, step, tol, "Watson" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::PenaltyI<double,void*> );

      justdoit( fct_ptr( tstoptfct::PenaltyIInit<double> ),
		fct, 4, pars, lo, hi, step, tol, "PenaltyI" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::PenaltyII<double,void*> );

      justdoit( fct_ptr( tstoptfct::PenaltyIIInit<double> ),
		fct, 4, pars, lo, hi, step, tol, "PenaltyII" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::VariablyDimensioned<double,void*> );

      justdoit( fct_ptr( tstoptfct::VariablyDimensionedInit<double> ),
		fct, npars, pars, lo, hi, step, tol, "VariablyDimensioned" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::Trigonometric<double,void*> );

      justdoit( fct_ptr( tstoptfct::TrigonometricInit<double> ),
		fct, npars, pars, lo, hi, step, tol, "Trigonometric" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::BrownAlmostLinear<double,void*> );

      justdoit( fct_ptr( tstoptfct::BrownAlmostLinearInit<double> ),
		fct, npars, pars, lo, hi, step, tol, "BrownAlmostLinear" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::DiscreteBoundary<double,void*> );

      justdoit( fct_ptr( tstoptfct::DiscreteBoundaryInit<double> ),
		fct, npars, pars, lo, hi, step, tol, "DiscreteBoundary" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::DiscreteIntegral<double,void*> );

      justdoit( fct_ptr( tstoptfct::DiscreteIntegralInit<double> ),
		fct, npars, pars, lo, hi, step, tol, "DiscreteIntegral" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::BroydenTridiagonal<double,void*> );

      justdoit( fct_ptr( tstoptfct::BroydenTridiagonalInit<double> ),
		fct, npars, pars, lo, hi, step, tol, "BroydenTridiagonal" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::BroydenBanded<double,void*> );

      justdoit( fct_ptr( tstoptfct::BroydenBandedInit<double> ),
		fct, npars, pars, lo, hi, step, tol, "BroydenBanded" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::LinearFullRank<double,void*> );

      justdoit( fct_ptr( tstoptfct::LinearFullRankInit<double> ),
		fct, npars, pars, lo, hi, step, tol, "LinearFullRank" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::LinearFullRank1<double,void*> );

      justdoit( fct_ptr( tstoptfct::LinearFullRank1Init<double> ),
		fct, npars, pars, lo, hi, step, tol, "LinearFullRank1" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::LinearFullRank0cols0rows<double,void*> );

      justdoit( fct_ptr( tstoptfct::LinearFullRank0cols0rowsInit<double> ),
		fct, npars, pars, lo, hi, step, tol, "LinearFullRank0cols0rows" );
    }

    {
      FctPtr< void, int, double*, double&, int&, void* >
	fct( tstoptfct::Chebyquad<double,void*> );

      justdoit( fct_ptr( tstoptfct::ChebyquadInit<double> ),
		fct, 9, pars, lo, hi, step, tol, "Chebyquad" );
    }

}

void tstglobal( int npars, double tol ) {

  int size = npars * npars * 4;
  std::vector< double > pars( size ), step( size ), lo( size ),
    hi( size );

  for ( int ii = 0; ii < npars; ++ii )
    step[ ii ] = 0.4;

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::McCormick<double,void*> );

    justdoit( fct_ptr( tstoptfct::McCormickInit<double> ),
	      fct, 2, pars, lo, hi, step, tol, "McCormick" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::BoxBetts<double,void*> );

    justdoit( fct_ptr( tstoptfct::BoxBettsInit<double> ),
	      fct, 3, pars, lo, hi, step, tol, "BoxBetts" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Paviani<double,void*> );

    justdoit( fct_ptr( tstoptfct::PavianiInit<double> ),
	      fct, 10, pars, lo, hi, step, tol, "Paviani" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::GoldsteinPrice<double,void*> );

    justdoit( fct_ptr( tstoptfct::GoldsteinPriceInit<double> ),
	      fct, 2, pars, lo, hi, step, tol, "GoldsteinPrice" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Shekel5<double,void*> );

    justdoit( fct_ptr( tstoptfct::Shekel5Init<double> ),
	      fct, 4, pars, lo, hi, step, tol, "Shekel5" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Shekel7<double,void*> );

    justdoit( fct_ptr( tstoptfct::Shekel7Init<double> ),
	      fct, 4, pars, lo, hi, step, tol, "Shekel7" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Shekel10<double,void*> );

    justdoit( fct_ptr( tstoptfct::Shekel10Init<double> ),
	      fct, 4, pars, lo, hi, step, tol, "Shekel10" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Levy<double,void*> );

    justdoit( fct_ptr( tstoptfct::LevyInit<double> ),
	      fct, 4, pars, lo, hi, step, tol, "Levy4" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Levy<double,void*> );

    justdoit( fct_ptr( tstoptfct::LevyInit<double> ),
	      fct, 5, pars, lo, hi, step, tol, "Levy5" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Levy<double,void*> );

    justdoit( fct_ptr( tstoptfct::LevyInit<double> ),
	      fct, 6, pars, lo, hi, step, tol, "Levy6" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Levy<double,void*> );

    justdoit( fct_ptr( tstoptfct::LevyInit<double> ),
	      fct, 7, pars, lo, hi, step, tol, "Levy7" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Griewank<double,void*> );

    justdoit( fct_ptr( tstoptfct::GriewankInit<double> ),
	      fct, 2, pars, lo, hi, step, tol, "Griewank" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::SixHumpCamel<double,void*> );

    justdoit( fct_ptr( tstoptfct::SixHumpCamelInit<double> ),
	      fct, 2, pars, lo, hi, step, tol, "SixHumpCamel" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Branin<double,void*> );

    justdoit( fct_ptr( tstoptfct::BraninInit<double> ),
	      fct, 2, pars, lo, hi, step, tol,
	      "Branin" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Shubert<double,void*> );

    justdoit( fct_ptr( tstoptfct::ShubertInit<double> ),
	      fct, 2, pars, lo, hi, step, tol, "Shubert" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Hansen<double,void*> );

    justdoit( fct_ptr( tstoptfct::HansenInit<double> ),
	      fct, 2, pars, lo, hi, step, tol,
	      "Hansen" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Cola<double,void*> );

    justdoit( fct_ptr( tstoptfct::ColaInit<double> ),
	      fct, 17, pars, lo, hi, step, tol, "Cola" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Ackley<double,void*> );

    justdoit( fct_ptr( tstoptfct::AckleyInit<double> ),
	      fct, 2, pars, lo, hi, step, tol, "Ackley" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Bohachevsky1<double,void*> );

    justdoit( fct_ptr( tstoptfct::Bohachevsky1Init<double> ),
	      fct, 2, pars, lo, hi, step, tol, "Bohachevsky1" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Bohachevsky2<double,void*> );

    justdoit( fct_ptr( tstoptfct::Bohachevsky2Init<double> ),
	      fct, 2, pars, lo, hi, step, tol, "Bohachevsky2" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Bohachevsky3<double,void*> );

    justdoit( fct_ptr( tstoptfct::Bohachevsky3Init<double> ),
	      fct, 2, pars, lo, hi, step, tol, "Bohachevsky3" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::DixonPrice<double,void*> );

    justdoit( fct_ptr( tstoptfct::DixonPriceInit<double> ),
	      fct, 25, pars, lo, hi, step, tol, "DixonPrice" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Easom<double,void*> );

    justdoit( fct_ptr( tstoptfct::EasomInit<double> ),
	      fct, 2, pars, lo, hi, step, tol, "Easom" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Rastrigin<double,void*> );

    justdoit( fct_ptr( tstoptfct::RastriginInit<double> ),
	      fct, 2, pars, lo, hi, step, tol, "Rastrigin" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Michalewicz<double,void*> );

    justdoit( fct_ptr( tstoptfct::MichalewiczInit<double> ),
	      fct, 2, pars, lo, hi, step, tol, "Michalewicz2" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Michalewicz<double,void*> );

    justdoit( fct_ptr( tstoptfct::MichalewiczInit<double> ),
	      fct, 5, pars, lo, hi, step, tol, "Michalewicz5" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Michalewicz<double,void*> );

    justdoit( fct_ptr( tstoptfct::MichalewiczInit<double> ),
	      fct, 10, pars, lo, hi, step, tol, "Michalewicz10" );
  }

  return;

}

int main( int argc, char* argv[] ) {

  try {

    int c, uncopt = 1, globalopt = 1;
    while ( --argc > 0 && (*++argv)[ 0 ] == '-' )
      while ( c = *++argv[ 0 ] )
	switch( c ) {
	case 'u':
	  uncopt = 0;
	  break;
	case 'g':
	  globalopt = 0;
	  break;
	default:
	  fprintf( stderr, "%s: illegal option '%c'\n", argv[ 0 ], c );
	  fprintf( stderr, "Usage %s [ -g ] [ -u ] [ npars ]\n", argv[ 0 ] );
	  return EXIT_FAILURE;
      }


    int npars=6;
    if ( argc == 1 )
      npars = atoi( *argv );
    
    if ( npars % 2 || npars < 2 ) {
      printf( "The minimum value for the free parameter must be an even "
	      "and it is greater then 2\n" );
      return EXIT_FAILURE;
    }

    double tol = 1.0e-8;

    std::cout << "#:tol=" << tol << '\n';
    std::cout << "# A negative value for the nfev signifies that the "
      "optimization method did not converge\n#\n";
    std::cout << "name\tnfev\tanswer\tstat\tpars\nS\tN\tN\tN\tN\n";

    if ( uncopt )
      tstuncopt( npars, tol );
    if ( globalopt )
      tstglobal( npars, tol );

    return EXIT_SUCCESS;

  } catch( std::exception& e ) {

    std::cerr << e.what( ) << '\n';
    return EXIT_FAILURE;

  }

}

/*
gcc -g -Wall -pedantic -ansi -c -O3 ../../utils/src/gsl/fcmp.c
g++ -g -Wall -pedantic -ansi -c -O3 -I../../include/ -I../../utils/src/gsl Simplex.cc
g++ -g -Wall -pedantic -ansi -O3 -I../../include/ -I../tests -I../../utils/src/gsl -DtestNelderMead NelderMead.cc Simplex.o fcmp.o
valgrind --tool=memcheck --leak-check=yes --show-reachable=yes a.out
==25763== Memcheck, a memory error detector.
==25763== Copyright (C) 2002-2006, and GNU GPL'd, by Julian Seward et al.
==25763== Using LibVEX rev 1658, a library for dynamic binary translation.
==25763== Copyright (C) 2004-2006, and GNU GPL'd, by OpenWorks LLP.
==25763== Using valgrind-3.2.1, a dynamic binary instrumentation framework.
==25763== Copyright (C) 2000-2006, and GNU GPL'd, by Julian Seward et al.
==25763== For more details, rerun with: -v
==25763==
#:tol=1e-08
# A negative value for the nfev signifies that the optimization method did not converge
#
name	nfev	answer	stat	pars
S	N	N	N	N
NelderMead_Rosenbrock	1197	0	0.00435378	0.963276,0.927824,1.04977,1.10213,1.02286,1.04643
NelderMead_Rosenbrock	1654	0	5.84053e-09	1.00004,1.00009,0.99997,0.999943,1.00001,1.00001
NelderMead_FreudensteinRoth	-530	0	146.953	11.4139,-0.896739,11.414,-0.89672,11.4117,-0.896875
NelderMead_FreudensteinRoth	-735	0	146.953	11.4128,-0.896788,11.4132,-0.896803,11.4129,-0.896808
NelderMead_PowellBadlyScaled	2435	0	1.04824e-08	1.0905e-05,9.17007,8.45513e-06,11.8271,1.11313e-05,8.98371
NelderMead_PowellBadlyScaled	2799	0	1.04787e-08	1.08905e-05,9.18229,8.4445e-06,11.842,1.11157e-05,8.9963
NelderMead_BrownBadlyScaled	-215	0	1748.78	999965,2.44628e-05
NelderMead_BrownBadlyScaled	421	0	9.54873e-09	1e+06,1.99991e-06
NelderMead_Beale	430	0	1.55834e-07	2.99999,0.500019,2.9995,0.499856,3.00077,0.500199
NelderMead_Beale	661	0	2.87446e-08	2.99976,0.49995,3.00005,0.500008,3.0001,0.5
NelderMead_JennrichSampson	661	373.086	373.087	0.257817,0.257827,0.257839,0.257817,0.257812,0.257847
NelderMead_JennrichSampson	915	373.086	373.087	0.257821,0.257829,0.257843,0.257806,0.257818,0.257828
NelderMead_HelicalValley	153	0	2.41172e-09	0.999998,-2.33949e-05,-3.92257e-05
NelderMead_HelicalValley	285	0	1.51587e-09	0.999997,1.48411e-05,2.40208e-05
NelderMead_Bard	167	0.00821487	0.00821488	0.0824145,1.1331,2.34363
NelderMead_Bard	288	0.00821487	0.00821488	0.0824095,1.13299,2.34375
NelderMead_Gaussian	111	1.12793e-08	1.14499e-08	0.39895,1,1.81238e-05
NelderMead_Gaussian	229	1.12793e-08	1.14499e-08	0.39895,1,1.81238e-05
NelderMead_Meyer	-30	87.9458	1.139e+07	0.0597053,4000.76,249.722
NelderMead_Meyer	2525	87.9458	87.9459	0.00560961,6181.35,345.224
NelderMead_GulfResearchDevelopment	488	0	6.03041e-11	58.0921,24.2541,1.53308
NelderMead_GulfResearchDevelopment	1356	0	8.68083e-15	50.0873,24.9913,1.50039
NelderMead_Box3d	226	0	1.64492e-10	1.00002,9.99958,0.999978
NelderMead_Box3d	320	0	1.64492e-10	1.00002,9.99958,0.999978
NelderMead_PowellSingular	388	0	1.66094e-12	0.000922767,-9.22902e-05,0.000441954,0.000441754
NelderMead_PowellSingular	576	0	3.73418e-15	4.80257e-06,-4.8136e-07,5.67939e-05,5.68199e-05
NelderMead_Wood	711	0	3.06143e-08	0.999946,0.999899,1.00004,1.0001
NelderMead_Wood	875	0	4.01606e-09	1.00003,1.00005,0.999977,0.999953
NelderMead_KowalikOsborne	278	0.000307505	0.000307506	0.192809,0.191312,0.123089,0.136076
NelderMead_KowalikOsborne	441	0.000307505	0.000307506	0.192806,0.191257,0.123032,0.136052
NelderMead_BrownDennis	234	85822.2	85822.2	-11.5944,13.2035,-0.403895,0.236169
NelderMead_BrownDennis	355	85822.2	85822.2	-11.5946,13.2038,-0.40314,0.23692
NelderMead_Osborne1	460	5.46489e-05	5.50277e-05	0.376257,2.04397,-1.57345,0.0130751,0.0217191
NelderMead_Osborne1	659	5.46489e-05	5.50277e-05	0.376257,2.04397,-1.57345,0.0130751,0.0217191
NelderMead_Biggs	716	0	4.3845e-33	1.00254,10,1,5,4,3
NelderMead_Biggs	1419	0	4.3845e-33	1.00254,10,1,5,4,3
NelderMead_Osborne2	-534	0.0401377	0.447557	1.14053,-0.0813793,0.388371,0.52296,0.256453,0.21697,5,7,2,4.5,5.5
NelderMead_Osborne2	-745	0.0401377	0.447557	1.14053,-0.0813988,0.388347,0.523085,0.256463,0.216905,5,7,2,4.5,5.5
NelderMead_Watson	627	0.00228767	0.00228767	-0.0157235,1.01243,-0.232966,1.2604,-1.51373,0.993008
NelderMead_Watson	894	0.00228767	0.00228767	-0.0157231,1.01243,-0.232989,1.26046,-1.51378,0.993028
NelderMead_PenaltyI	-1150	9.37629e-06	2.24998e-05	0.250536,0.249133,0.249885,0.250474
NelderMead_PenaltyI	-1364	9.37629e-06	2.24998e-05	0.250006,0.249975,0.249997,0.250053
NelderMead_PenaltyII	1855	9.37629e-06	9.50843e-06	0.199983,0.203119,0.575966,0.229747
NelderMead_PenaltyII	3536	9.37629e-06	9.37784e-06	0.199999,0.203599,0.46345,0.534856
NelderMead_VariablyDimensioned	303	0	3.475e-08	0.999994,0.999897,0.999961,1.00012,0.999966,0.999989
NelderMead_VariablyDimensioned	538	0	9.81666e-09	0.999977,0.999949,0.999976,1.00006,0.999977,1
NelderMead_Trigonometric	433	0	1.40747e-08	0.00705824,0.00696962,0.0067329,-0.117601,0.00653846,0.0064437
NelderMead_Trigonometric	676	0	2.53731e-09	0.00703642,0.00688127,0.00676746,-0.117661,0.00649293,0.00640494
NelderMead_BrownAlmostLinear	-852	1	5.10492e-10	1,0.999991,0.999979,1,1,1.00003
NelderMead_BrownAlmostLinear	-1089	1	5.10492e-10	1,0.999991,0.999979,1,1,1.00003
NelderMead_DiscreteBoundary	289	0	1.49323e-09	-0.0898696,-0.167835,-0.1536,-0.13246,-0.102081,-0.0593097
NelderMead_DiscreteBoundary	542	0	5.51653e-10	-0.0898787,-0.167848,-0.153594,-0.132449,-0.102034,-0.0592801
NelderMead_DiscreteIntegral	266	0	1.76342e-09	-0.065457,-0.118143,-0.154624,-0.169972,-0.157286,-0.106014
NelderMead_DiscreteIntegral	531	0	1.76342e-09	-0.065457,-0.118143,-0.154624,-0.169972,-0.157286,-0.106014
NelderMead_BroydenTridiagonal	230	0	1.04327e-07	-0.576006,-0.695934,-0.680257,-0.643005,-0.55646,-0.36606
NelderMead_BroydenTridiagonal	494	0	1.83315e-08	-0.576063,-0.695951,-0.680244,-0.642993,-0.556434,-0.366046
NelderMead_BroydenBanded	235	0	9.4919e-08	-0.42827,-0.476627,-0.519657,-0.558064,-0.593413,-0.593428
NelderMead_BroydenBanded	505	0	2.29414e-08	-0.428294,-0.476602,-0.519652,-0.558061,-0.593425,-0.593428
NelderMead_LinearFullRank	346	0	9.80237e-09	-0.999944,-1.00002,-0.999936,-0.999982,-1.00004,-1.00001
NelderMead_LinearFullRank	592	0	6.5498e-09	-1.00006,-0.99997,-1.00004,-1.00001,-1.00003,-1
NelderMead_LinearFullRank1	199	1.15385	1.15385	2.97371,2.20459,2.13584,0.159209,0.599986,-2.86608
NelderMead_LinearFullRank1	481	1.15385	1.15385	3.07835,2.28814,2.20579,0.177665,0.583544,-2.94494
NelderMead_LinearFullRank0cols0rows	189	2.66667	2.66667	2.14804,1.95635,0.261106,0.864515,-1.56416,2.61238
NelderMead_LinearFullRank0cols0rows	454	2.66667	2.66667	2.54627,1.95813,0.262525,0.864232,-1.56548,3.05755
NelderMead_Chebyquad	-670	0	0.0162113	0.0348735,0.161667,0.22691,0.423012,0.422982,0.611581,0.7,0.8,0.9
NelderMead_Chebyquad	-871	0	0.0162113	0.0348735,0.161667,0.22691,0.423012,0.422982,0.611581,0.7,0.8,0.9
NelderMead_McCormick	72	-1.91	-1.91322	-0.547215,-1.54724
NelderMead_McCormick	125	-1.91	-1.91322	-0.547215,-1.54724
NelderMead_BoxBetts	106	0	5.59281e-10	1.00007,9.99952,0.999948
NelderMead_BoxBetts	182	0	1.97261e-10	0.999998,10,1.00001
NelderMead_Paviani	-509	-45.7	-13.097	9.20636,9.20674,9.20642,9.20651,9.20586,9.20575,5,5,5,5
NelderMead_Paviani	-694	-45.7	-13.097	9.20626,9.20605,9.2062,9.20624,9.20614,9.20624,5,5,5,5
NelderMead_GoldsteinPrice	86	3	3	-4.52456e-05,-1.00002
NelderMead_GoldsteinPrice	157	3	3	8.1296e-06,-0.999997
NelderMead_Shekel5	-99	-10.1532	-5.10076	7.99928,7.99961,7.9991,7.99939
NelderMead_Shekel5	-239	-10.1532	-5.10077	7.99959,7.99961,7.99961,7.99965
NelderMead_Shekel7	-98	-10.4029	-5.12881	7.9994,7.99943,7.99961,7.99903
NelderMead_Shekel7	-237	-10.4029	-5.12882	7.99951,7.99965,7.99949,7.99962
NelderMead_Shekel10	-92	-10.5364	-5.17563	7.99989,7.99876,7.99943,7.99934
NelderMead_Shekel10	-235	-10.5364	-5.17565	7.99945,7.99938,7.99939,7.99946
NelderMead_Levy4	-728	-21.502	10.1306	0.0112008,0.342276,2.6177,6.999
NelderMead_Levy4	-939	-21.502	5.99894	0.999985,0.999799,1.00005,6.9979
NelderMead_Levy5	-142	-11.504	38.8617	3.96346,3.99613,3.99617,3.99623,3.9995
NelderMead_Levy5	-338	-11.504	38.8617	3.96377,3.99615,3.99623,3.99624,3.99947
NelderMead_Levy6	-169	-11.504	47.8505	3.96403,3.99616,3.99602,3.99616,3.99607,3.9995
NelderMead_Levy6	-399	-11.504	47.8504	3.96381,3.99616,3.99624,3.99623,3.99624,3.99945
NelderMead_Levy7	-174	-11.504	56.8395	3.96373,3.99634,3.99594,3.99626,3.99612,3.99613,4
NelderMead_Levy7	-401	-11.504	56.8394	3.96381,3.99615,3.99624,3.99624,3.99623,3.99625,4
NelderMead_Griewank	-44	0	4.91141	100.482,-97.6443
NelderMead_Griewank	-89	0	4.91141	100.481,-97.646
NelderMead_SixHumpCamel	87	-1.03	-1.03163	0.0898309,-0.71266
NelderMead_SixHumpCamel	145	-1.03	-1.03163	0.0898309,-0.71266
NelderMead_Branin	-68	0	0.618436	3.73049,-1.23055
NelderMead_Branin	-123	0	0.618436	3.73055,-1.23069
NelderMead_Shubert	-48	-24.06	-7.21602	6.85993,6.85985
NelderMead_Shubert	-109	-24.06	-7.21603	6.86013,6.86012
NelderMead_Hansen	-52	-176.54	-13.2572	3.34968,3.77213
NelderMead_Hansen	-116	-176.54	-13.2573	3.3496,3.77232
NelderMead_Cola	-1296	12.8154	164.519	1.48672,1.22147,2.32856,-0.607661,2.42371,2.33567,0,0,0,0,0,0,0,0,0,0,0
NelderMead_Cola	-1870	12.8154	164.519	1.47302,1.22133,2.32858,-0.607775,2.42368,2.33566,0,0,0,0,0,0,0,0,0,0,0
NelderMead_Ackley	-40	0	19.3325	16.9991,16.9992
NelderMead_Ackley	-96	0	19.3325	16.9988,16.9988
NelderMead_Bohachevsky1	122	0	1.25887e-08	-2.78979e-05,-6.54913e-06
NelderMead_Bohachevsky1	180	0	9.97341e-09	8.17281e-06,1.63857e-05
NelderMead_Bohachevsky2	122	0	2.88746e-08	3.37759e-05,-2.20893e-05
NelderMead_Bohachevsky2	183	0	2.00189e-09	1.14095e-05,-2.31143e-06
NelderMead_Bohachevsky3	139	0	1.1725e-09	2.18609e-05,-1.2337e-05
NelderMead_Bohachevsky3	199	0	1.1725e-09	2.18609e-05,-1.2337e-05
NelderMead_DixonPrice	-1180	0	30619.6	1.16587,0.79024,0.694572,0.72751,0.88738,1.30278,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5
NelderMead_DixonPrice	-1335	0	30619.6	1.16587,0.79024,0.694572,0.72751,0.88738,1.30278,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5
NelderMead_Easom	-31	-1	-0	25,25
NelderMead_Easom	-62	-1	-0	25,25
NelderMead_Rastrigin	-53	0	17.9092	2.98486,2.98493
NelderMead_Rastrigin	-115	0	17.9092	2.98484,2.98485
NelderMead_Michalewicz2	56	-1.8013	-1.8013	2.20297,1.57073
NelderMead_Michalewicz2	115	-1.8013	-1.8013	2.20293,1.57079
NelderMead_Michalewicz5	228	-4.68766	-4.68766	2.20305,1.57074,1.28517,1.92306,1.72046
NelderMead_Michalewicz5	411	-4.68766	-4.68766	2.2029,1.57082,1.285,1.92306,1.72046
NelderMead_Michalewicz10	-1263	-9.66015	-6.68961	2.20293,1.57075,1.28509,1.92312,1.7205,1.57082,1.5708,1.5708,1.5708,1.5708
NelderMead_Michalewicz10	-1481	-9.66015	-6.68961	2.2029,1.5708,1.285,1.92306,1.72047,1.57079,1.5708,1.5708,1.5708,1.5708

==25763==
==25763== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 15 from 1)
==25763== malloc/free: in use at exit: 0 bytes in 0 blocks.
==25763== malloc/free: 43,692 allocs, 43,692 frees, 4,385,292 bytes allocated.
==25763== For counts of detected errors, rerun with: -v
==25763== All heap blocks were freed -- no leaks are possible.
*/
#endif
