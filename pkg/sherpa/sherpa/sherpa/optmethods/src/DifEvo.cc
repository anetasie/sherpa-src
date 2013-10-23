#ifdef testDifEvo

#include "DifEvo.hh"

#include "fcmp.h"
#include "tstoptfct.hh"
#include "sherpa/functor.hh"

#include "minpack/levmar.hh"
#include "NelderMead.hh"

using namespace sherpa;

template <typename Real>
void print_pars( const char* prefix, const char* name, int nfev, Real stat,
		 Real answer, int n, const std::vector< Real >& x,
		 Real tol=
		 1.0e4*std::sqrt( std::numeric_limits< Real >::epsilon() ) ) {

  
  std::cout << prefix << name << '\t';
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

#define NPOP 16
#define MAXNFEV 64

template< typename Init, typename Fct >
void justdoit( Init init, Fct fct, int npars, 
	       std::vector< double >& pars, std::vector< double >& lo,
	       std::vector< double >& hi, double tol, const char* header ) {

  int mfcts, seed=1357;
  double answer;
  double xprob = 0.9, scale=1.0;
  int begin_strategy = 0, end_strategy = 10;

  for ( int strategy = begin_strategy; strategy < end_strategy; ++ strategy ) {

    init( npars, mfcts, answer, &pars[0], &lo[0], &hi[0] );

    sherpa::DifEvo< Fct, void*,
      sherpa::NelderMead< Fct, void* > > de_nm( npars, &pars[0], &lo[0],
						&hi[0], fct, NULL, xprob,
						scale, strategy, seed );

    sherpa::DifEvo< Fct, void*,
      sherpa::OptFunc< Fct, void* > > de( npars, &pars[0], &lo[0], &hi[0], fct,
					  NULL, xprob, scale, strategy, seed );

    int nfev, verbose = 0, population_size = NPOP * npars;
    int maxnfev = MAXNFEV * npars * population_size;
    double fmin;
    de_nm( &pars[0], verbose, maxnfev, tol, population_size, nfev, fmin );
    print_pars( "DifEvo_nm_", header, nfev, fmin, answer, npars, pars );

    de( &pars[0], verbose, maxnfev, tol, population_size, nfev, fmin );
    print_pars( "DifEvo_", header, nfev, fmin, answer, npars, pars );

  }

}

template< typename Init, typename Fct >
void justdoitlm( Init init, Fct fct, int npars, 
		 std::vector< double >& pars, std::vector< double >& lo,
		 std::vector< double >& hi, double tol, const char* header ) {

  int mfcts, seed=1357;
  double answer;
  double xprob = 0.9, scale=1.0;
  int begin_strategy = 0, end_strategy = 10;

  for ( int strategy = begin_strategy; strategy < end_strategy; ++ strategy ) {

    init( npars, mfcts, answer, &pars[0], &lo[0], &hi[0] );

    sherpa::DifEvo< Fct, void*,
      minpack::LevMar< Fct, void* > > de( npars, &pars[0], &lo[0], &hi[0],
					  fct, NULL, xprob, scale, strategy,
					  seed, mfcts );

    int nfev, verbose = 0, population_size = NPOP * npars;
    int maxnfev = MAXNFEV * npars * population_size;
    double fmin;
    de( &pars[0], verbose, maxnfev, tol, population_size, nfev, fmin );
    
    print_pars( "DifEvo_lm_", header, nfev, fmin, answer, npars, pars );

  }

}

void tstuncopt( int npars, double tol ) {

  const int verbose=0, size=npars*32;
  std::vector< double > par( size, 0 ), lo( size, -1.0e2 ), hi( size, 1.0e2 );

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Rosenbrock<double,void*> );

    justdoit( fct_ptr( tstoptfct::RosenbrockInit<double> ), fct, npars, par,
	      lo, hi, tol, "Rosenbrock" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::FreudensteinRoth<double,void*> );

    justdoit( fct_ptr( tstoptfct::FreudensteinRothInit<double> ),
	      fct, npars, par, lo, hi, tol, "FreudensteinRoth" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::PowellBadlyScaled<double,void*> );

    justdoit( fct_ptr( tstoptfct::PowellBadlyScaledInit<double> ),
	      fct, npars, par, lo, hi, tol, "PowellBadlyScaled" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::BrownBadlyScaled<double,void*> );

    justdoit( fct_ptr( tstoptfct::BrownBadlyScaledInit<double> ),
	      fct, 2, par, lo, hi, tol, "BrownBadlyScaled" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Beale<double,void*> );

    justdoit( fct_ptr( tstoptfct::BealeInit<double> ),
	      fct, npars, par, lo, hi, tol, "Beale" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::JennrichSampson<double,void*> );

    justdoit( fct_ptr( tstoptfct::JennrichSampsonInit<double> ),
	      fct, npars, par, lo, hi, tol,
	      "JennrichSampson" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::HelicalValley<double,void*> );

    justdoit( fct_ptr( tstoptfct::HelicalValleyInit<double> ),
	      fct, 3, par, lo, hi, tol, "HelicalValley" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Bard<double,void*> );

    justdoit( fct_ptr( tstoptfct::BardInit<double> ),
	      fct, 3, par, lo, hi, tol, "Bard" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Gaussian<double,void*> );

    justdoit( fct_ptr( tstoptfct::GaussianInit<double> ),
	      fct, 3, par, lo, hi, tol, "Gaussian" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Meyer<double,void*> );

    justdoit( fct_ptr( tstoptfct::MeyerInit<double> ),
	      fct, 3, par, lo, hi, tol, "Meyer" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::GulfResearchDevelopment<double,void*> );

    justdoit( fct_ptr( tstoptfct::GulfResearchDevelopmentInit<double> ),
	      fct, 3, par, lo, hi, tol, "GulfResearchDevelopment" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Box3d<double,void*> );

    justdoit( fct_ptr( tstoptfct::Box3dInit<double> ),
	      fct, 3, par, lo, hi, tol, "Box3d" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::PowellSingular<double,void*> );

    justdoit( fct_ptr( tstoptfct::PowellSingularInit<double> ),
	      fct, 4, par, lo, hi, tol, "PowellSingular" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Wood<double,void*> );

    justdoit( fct_ptr( tstoptfct::WoodInit<double> ),
	      fct, 4, par, lo, hi, tol, "Wood" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::KowalikOsborne<double,void*> );

    justdoit( fct_ptr( tstoptfct::KowalikOsborneInit<double> ),
	      fct, 4, par, lo, hi, tol, "KowalikOsborne" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::BrownDennis<double,void*> );

    justdoit( fct_ptr( tstoptfct::BrownDennisInit<double> ),
	      fct, 4, par, lo, hi, tol, "BrownDennis" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Osborne1<double,void*> );

    justdoit( fct_ptr( tstoptfct::Osborne1Init<double> ),
	      fct, 5, par, lo, hi, tol, "Osborne1" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Biggs<double,void*> );

    justdoit( fct_ptr( tstoptfct::BiggsInit<double> ),
	      fct, 6, par, lo, hi, tol, "Biggs" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Osborne2<double,void*> );

    justdoit( fct_ptr( tstoptfct::Osborne2Init<double> ),
	      fct, 11, par, lo, hi, tol, "Osborne2" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Watson<double,void*> );

    justdoit( fct_ptr( tstoptfct::WatsonInit<double> ),
	      fct, 6, par, lo, hi, tol, "Watson" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::PenaltyI<double,void*> );

    justdoit( fct_ptr( tstoptfct::PenaltyIInit<double> ),
	      fct, 4, par, lo, hi, tol, "PenaltyI" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::PenaltyII<double,void*> );

    justdoit( fct_ptr( tstoptfct::PenaltyIIInit<double> ),
	      fct, 4, par, lo, hi, tol, "PenaltyII" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::VariablyDimensioned<double,void*> );

    justdoit( fct_ptr( tstoptfct::VariablyDimensionedInit<double> ),
	      fct, npars, par, lo, hi, tol, "VariablyDimensioned" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Trigonometric<double,void*> );

    justdoit( fct_ptr( tstoptfct::TrigonometricInit<double> ),
	      fct, npars, par, lo, hi, tol, "Trigonometric" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::BrownAlmostLinear<double,void*> );

    justdoit( fct_ptr( tstoptfct::BrownAlmostLinearInit<double> ),
	      fct, npars, par, lo, hi, tol, "BrownAlmostLinear" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::DiscreteBoundary<double,void*> );

    justdoit( fct_ptr( tstoptfct::DiscreteBoundaryInit<double> ),
	      fct, npars, par, lo, hi, tol, "DiscreteBoundary" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::DiscreteIntegral<double,void*> );

    justdoit( fct_ptr( tstoptfct::DiscreteIntegralInit<double> ),
	      fct, npars, par, lo, hi, tol, "DiscreteIntegral" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::BroydenTridiagonal<double,void*> );

    justdoit( fct_ptr( tstoptfct::BroydenTridiagonalInit<double> ),
	      fct, npars, par, lo, hi, tol, "BroydenTridiagonal" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::BroydenBanded<double,void*> );

    justdoit( fct_ptr( tstoptfct::BroydenBandedInit<double> ),
	      fct, npars, par, lo, hi, tol, "BroydenBanded" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::LinearFullRank<double,void*> );

    justdoit( fct_ptr( tstoptfct::LinearFullRankInit<double> ),
	      fct, npars, par, lo, hi, tol, "LinearFullRank" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::LinearFullRank1<double,void*> );

    justdoit( fct_ptr( tstoptfct::LinearFullRank1Init<double> ),
	      fct, npars, par, lo, hi, tol, "LinearFullRank1" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::LinearFullRank0cols0rows<double,void*> );

    justdoit( fct_ptr( tstoptfct::LinearFullRank0cols0rowsInit<double> ),
	      fct, npars, par, lo, hi, tol, "LinearFullRank0cols0rows" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Chebyquad<double,void*> );

    justdoit( fct_ptr( tstoptfct::ChebyquadInit<double> ),
	      fct, 9, par, lo, hi, tol, "Chebyquad" );
  }

  return;

}

void tstuncoptlm( int npars, double tol ) {

  const int verbose=0, size=npars*32;
  std::vector< double > par( size, 0 ), lo( size, -1.0e2 ), hi( size, 1.0e2 );

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Rosenbrock<double,void*> );
    
    justdoitlm( sherpa::fct_ptr( tstoptfct::RosenbrockInit<double> ),
		fct, npars, par, lo, hi, tol, "Rosenbrock" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::FreudensteinRoth<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::FreudensteinRothInit<double> ),
		fct, npars, par, lo, hi, tol, "FreudensteinRoth" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::PowellBadlyScaled<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::PowellBadlyScaledInit<double> ),
		fct, npars, par, lo, hi, tol, "PowellBadlyScaled" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::BrownBadlyScaled<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::BrownBadlyScaledInit<double> ),
	      fct, 2, par, lo, hi, tol, "BrownBadlyScaled" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Beale<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::BealeInit<double> ),
	      fct, npars, par, lo, hi, tol, "Beale" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::JennrichSampson<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::JennrichSampsonInit<double> ),
	      fct, npars, par, lo, hi, tol, "JennrichSampson" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::HelicalValley<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::HelicalValleyInit<double> ),
	      fct, 3*npars, par, lo, hi, tol, "HelicalValley" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Bard<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::BardInit<double> ),
	      fct, 3*npars, par, lo, hi, tol, "Bard" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Gaussian<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::GaussianInit<double> ),
	      fct, 3, par, lo, hi, tol, "Gaussian" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Meyer<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::MeyerInit<double> ),
	      fct, 3, par, lo, hi, tol, "Meyer" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::GulfResearchDevelopment<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::GulfResearchDevelopmentInit<double> ),
	      fct, 3, par, lo, hi, tol, "GulfResearchDevelopment" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Box3d<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::Box3dInit<double> ),
	      fct, 3, par, lo, hi, tol, "Box3d" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::PowellSingular<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::PowellSingularInit<double> ),
	      fct, 4*npars, par, lo, hi, tol, "PowellSingular" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Wood<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::WoodInit<double> ),
	      fct, 4, par, lo, hi, tol, "Wood" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::KowalikOsborne<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::KowalikOsborneInit<double> ),
	      fct, 4, par, lo, hi, tol, "KowalikOsborne" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::BrownDennis<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::BrownDennisInit<double> ),
	      fct, 4, par, lo, hi, tol, "BrownDennis" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Osborne1<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::Osborne1Init<double> ),
	      fct, 5, par, lo, hi, tol, "Osborne1" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Biggs<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::BiggsInit<double> ),
	      fct, 6, par, lo, hi, tol, "Biggs" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Osborne2<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::Osborne2Init<double> ),
	      fct, 11, par, lo, hi, tol, "Osborne2" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Watson<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::WatsonInit<double> ),
	      fct, 6, par, lo, hi, tol, "Watson" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::PenaltyI<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::PenaltyIInit<double> ),
	      fct, 4, par, lo, hi, tol, "PenaltyI" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::PenaltyII<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::PenaltyIIInit<double> ),
	      fct, 4, par, lo, hi, tol, "PenaltyII" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::VariablyDimensioned<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::VariablyDimensionedInit<double> ),
	      fct, npars, par, lo, hi, tol, "VariablyDimensioned" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Trigonometric<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::TrigonometricInit<double> ),
	      fct, npars, par, lo, hi, tol, "Trigonometric" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::BrownAlmostLinear<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::BrownAlmostLinearInit<double> ),
	      fct, npars, par, lo, hi, tol, "BrownAlmostLinear" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::DiscreteBoundary<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::DiscreteBoundaryInit<double> ),
	      fct, npars, par, lo, hi, tol, "DiscreteBoundary" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::DiscreteIntegral<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::DiscreteIntegralInit<double> ),
	      fct, npars, par, lo, hi, tol, "DiscreteIntegral" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::BroydenTridiagonal<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::BroydenTridiagonalInit<double> ),
	      fct, npars, par, lo, hi, tol, "BroydenTridiagonal" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::BroydenBanded<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::BroydenBandedInit<double> ),
	      fct, npars, par, lo, hi, tol, "BroydenBanded" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::LinearFullRank<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::LinearFullRankInit<double> ),
	      fct, npars, par, lo, hi, tol, "LinearFullRank" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::LinearFullRank1<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::LinearFullRank1Init<double> ),
	      fct, npars, par, lo, hi, tol, "LinearFullRank1" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::LinearFullRank0cols0rows<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::LinearFullRank0cols0rowsInit<double> ),
	      fct, npars, par, lo, hi, tol, "LinearFullRank0cols0rows" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Chebyquad<double,void*> );

    justdoitlm( sherpa::fct_ptr( tstoptfct::ChebyquadInit<double> ),
	      fct, 9, par, lo, hi, tol, "Chebyquad" );
  }

  return;

}

void tstglobal( int npars, double tol ) {

  const int verbose = 0, size=npars*32;
  std::vector< double > par( size, 0 ), lo( size, -1.0e2 ), hi( size, 1.0e2 );

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::McCormick<double,void*> );

    justdoit( fct_ptr( tstoptfct::McCormickInit<double> ),
	      fct, 2, par, lo, hi, tol, "McCormick" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::BoxBetts<double,void*> );

    justdoit( fct_ptr( tstoptfct::BoxBettsInit<double> ),
	      fct, 3, par, lo, hi, tol, "BoxBetts" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Paviani<double,void*> );

    justdoit( fct_ptr( tstoptfct::PavianiInit<double> ),
	      fct, 10, par, lo, hi, tol, "Paviani" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::GoldsteinPrice<double,void*> );

    justdoit( fct_ptr( tstoptfct::GoldsteinPriceInit<double> ),
	      fct, 2, par, lo, hi, tol, "GoldsteinPrice" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Shekel5<double,void*> );

    justdoit( fct_ptr( tstoptfct::Shekel5Init<double> ),
	      fct, 4, par, lo, hi, tol, "Shekel5" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Shekel7<double,void*> );

    justdoit( fct_ptr( tstoptfct::Shekel7Init<double> ),
	      fct, 4, par, lo, hi, tol, "Shekel7" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Shekel10<double,void*> );

    justdoit( fct_ptr( tstoptfct::Shekel10Init<double> ),
	      fct, 4, par, lo, hi, tol, "Shekel10" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Levy<double,void*> );

    justdoit( fct_ptr( tstoptfct::LevyInit<double> ),
	      fct, 4, par, lo, hi, tol, "Levy4" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Levy<double,void*> );

    justdoit( fct_ptr( tstoptfct::LevyInit<double> ),
	      fct, 5, par, lo, hi, tol, "Levy5" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Levy<double,void*> );

    justdoit( fct_ptr( tstoptfct::LevyInit<double> ),
	      fct, 6, par, lo, hi, tol, "Levy6" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Levy<double,void*> );

    justdoit( fct_ptr( tstoptfct::LevyInit<double> ),
	      fct, 7, par, lo, hi, tol, "Levy7" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Griewank<double,void*> );

    justdoit( fct_ptr( tstoptfct::GriewankInit<double> ),
	      fct, 2, par, lo, hi, tol, "Griewank" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::SixHumpCamel<double,void*> );

    justdoit( fct_ptr( tstoptfct::SixHumpCamelInit<double> ),
	      fct, 2, par, lo, hi, tol, "SixHumpCamel" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Branin<double,void*> );

    justdoit( fct_ptr( tstoptfct::BraninInit<double> ),
	      fct, 2, par, lo, hi, tol,
	      "Branin" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Shubert<double,void*> );

    justdoit( fct_ptr( tstoptfct::ShubertInit<double> ),
	      fct, 2, par, lo, hi, tol, "Shubert" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Hansen<double,void*> );

    justdoit( fct_ptr( tstoptfct::HansenInit<double> ),
	      fct, 2, par, lo, hi, tol,
	      "Hansen" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Cola<double,void*> );

    justdoit( fct_ptr( tstoptfct::ColaInit<double> ),
	      fct, 17, par, lo, hi, tol, "Cola" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Ackley<double,void*> );

    justdoit( fct_ptr( tstoptfct::AckleyInit<double> ),
	      fct, 2, par, lo, hi, tol, "Ackley" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Bohachevsky1<double,void*> );

    justdoit( fct_ptr( tstoptfct::Bohachevsky1Init<double> ),
	      fct, 2, par, lo, hi, tol, "Bohachevsky1" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Bohachevsky2<double,void*> );

    justdoit( fct_ptr( tstoptfct::Bohachevsky2Init<double> ),
	      fct, 2, par, lo, hi, tol, "Bohachevsky2" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Bohachevsky3<double,void*> );

    justdoit( fct_ptr( tstoptfct::Bohachevsky3Init<double> ),
	      fct, 2, par, lo, hi, tol, "Bohachevsky3" );
  }


  // {
  //   FctPtr< void, int, double*, double&, int&, void* >
  //     fct( tstoptfct::DixonPrice<double,void*> );
  // 
  //   justdoit( fct_ptr( tstoptfct::DixonPriceInit<double> ),
  // 	      fct, 25, par, lo, hi, tol, "DixonPrice" );
  // }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Easom<double,void*> );

    justdoit( fct_ptr( tstoptfct::EasomInit<double> ),
	      fct, 2, par, lo, hi, tol, "Easom" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Rastrigin<double,void*> );

    justdoit( fct_ptr( tstoptfct::RastriginInit<double> ),
	      fct, 2, par, lo, hi, tol, "Rastrigin" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Michalewicz<double,void*> );

    justdoit( fct_ptr( tstoptfct::MichalewiczInit<double> ),
	      fct, 2, par, lo, hi, tol, "Michalewicz2" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Michalewicz<double,void*> );

    justdoit( fct_ptr( tstoptfct::MichalewiczInit<double> ),
	      fct, 5, par, lo, hi, tol, "Michalewicz5" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Michalewicz<double,void*> );

    justdoit( fct_ptr( tstoptfct::MichalewiczInit<double> ),
	      fct, 10, par, lo, hi, tol, "Michalewicz10" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::McCormick<double,void*> );

    justdoit( fct_ptr( tstoptfct::McCormickInit<double> ),
	      fct, 2, par, lo, hi, tol, "McCormick" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::BoxBetts<double,void*> );

    justdoit( fct_ptr( tstoptfct::BoxBettsInit<double> ),
	      fct, 3, par, lo, hi, tol, "BoxBetts" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Paviani<double,void*> );

    justdoit( fct_ptr( tstoptfct::PavianiInit<double> ),
	      fct, 10, par, lo, hi, tol, "Paviani" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::GoldsteinPrice<double,void*> );

    justdoit( fct_ptr( tstoptfct::GoldsteinPriceInit<double> ),
	      fct, 2, par, lo, hi, tol, "GoldsteinPrice" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Shekel5<double,void*> );

    justdoit( fct_ptr( tstoptfct::Shekel5Init<double> ),
	      fct, 4, par, lo, hi, tol, "Shekel5" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Shekel7<double,void*> );

    justdoit( fct_ptr( tstoptfct::Shekel7Init<double> ),
	      fct, 4, par, lo, hi, tol, "Shekel7" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Shekel10<double,void*> );

    justdoit( fct_ptr( tstoptfct::Shekel10Init<double> ),
	      fct, 4, par, lo, hi, tol, "Shekel10" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Griewank<double,void*> );

    justdoit( fct_ptr( tstoptfct::GriewankInit<double> ),
	      fct, 2, par, lo, hi, tol, "Griewank" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Ackley<double,void*> );

    justdoit( fct_ptr( tstoptfct::AckleyInit<double> ),
	      fct, 2, par, lo, hi, tol, "Ackley" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Bohachevsky1<double,void*> );

    justdoit( fct_ptr( tstoptfct::Bohachevsky1Init<double> ),
	      fct, 2, par, lo, hi, tol, "Bohachevsky1" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Bohachevsky2<double,void*> );

    justdoit( fct_ptr( tstoptfct::Bohachevsky2Init<double> ),
	      fct, 2, par, lo, hi, tol, "Bohachevsky2" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Bohachevsky3<double,void*> );

    justdoit( fct_ptr( tstoptfct::Bohachevsky3Init<double> ),
	      fct, 2, par, lo, hi, tol, "Bohachevsky3" );
  }


  // {
  //   FctPtr< void, int, double*, double&, int&, void* >
  //     fct( tstoptfct::DixonPrice<double,void*> );
  // 
  //  justdoit( fct_ptr( tstoptfct::DixonPriceInit<double> ),
  // 	      fct, 25, par, lo, hi, tol, "DixonPrice" );
  // }


  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Easom<double,void*> );

    justdoit( fct_ptr( tstoptfct::EasomInit<double> ),
	      fct, 2, par, lo, hi, tol, "Easom" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Rastrigin<double,void*> );

    justdoit( fct_ptr( tstoptfct::RastriginInit<double> ),
	      fct, 2, par, lo, hi, tol, "Rastrigin" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Michalewicz<double,void*> );

    justdoit( fct_ptr( tstoptfct::MichalewiczInit<double> ),
	      fct, 2, par, lo, hi, tol, "Michalewicz2" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Michalewicz<double,void*> );

    justdoit( fct_ptr( tstoptfct::MichalewiczInit<double> ),
	      fct, 5, par, lo, hi, tol, "Michalewicz5" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Michalewicz<double,void*> );

    justdoit( fct_ptr( tstoptfct::MichalewiczInit<double> ),
	      fct, 10, par, lo, hi, tol, "Michalewicz10" );
  }

  return;
}

int main( int argc, char* argv[] ) {

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

  int npars=2;
  if ( argc == 1 ) {
    npars = atoi( *argv );
    std::cout << "#\n# npars = " << npars << "\n#\n";
  }

  double tol=1.0e-6;
  std::cout << "#:tol=" << tol << '\n';
  std::cout << "#\n# A negative value for the nfev signifies that the "
    "optimization method did not converge\n#\n";
  std::cout << "name\tnfev\tanswer\tstat\tpars\nS\tN\tN\tN\tN\n";

  if ( uncopt ) {
    tstuncopt( npars, tol );
    tstuncoptlm( npars, tol );
  }
  if ( globalopt )
    tstglobal( npars, tol );

  return 0;

 }
/*
gcc -g -Wall -pedantic -ansi -c -O3 mt19937ar.c
gcc -g -Wall -pedantic -ansi -c -O3 ../../utils/src/gsl/fcmp.c
g++ -g -Wall -pedantic -ansi -c -O3 -I../../include/ -I../../utils/src/gsl Simplex.cc
g++ -g -Wall -pedantic -ansi -O3 -I. -I../../include/ -I../tests -I../../utils/src/gsl -DtestDifEvo DifEvo.cc Simplex.o fcmp.o mt19937ar.o
valgrind --tool=memcheck --leak-check=yes --show-reachable=yes a.out
==11191== Memcheck, a memory error detector.
==11191== Copyright (C) 2002-2006, and GNU GPL'd, by Julian Seward et al.
==11191== Using LibVEX rev 1658, a library for dynamic binary translation.
==11191== Copyright (C) 2004-2006, and GNU GPL'd, by OpenWorks LLP.
==11191== Using valgrind-3.2.1, a dynamic binary instrumentation framework.
==11191== Copyright (C) 2000-2006, and GNU GPL'd, by Julian Seward et al.
==11191== For more details, rerun with: -v
==11191== 
#
# A negative value for the nfev signifies that the optimization method did not converge
#
name	nfev	answer	stat	pars
S	N	N	N	N
DifEvo_McCormick	221	-1.91	-1.91322	-0.547614,-1.54759
DifEvo_McCormick	889	-1.91	-1.91322	-0.547614,-1.54759
DifEvo_McCormick	221	-1.91	-1.91322	-0.547614,-1.54759
DifEvo_McCormick	8192	-1.91	-1.91322	-0.547614,-1.54759
DifEvo_McCormick	1498	-1.91	-1.91322	-0.547404,-1.54723
DifEvo_McCormick	310	-1.91	-1.91322	-0.547614,-1.54759
DifEvo_McCormick	976	-1.91	-1.91322	-0.547013,-1.5473
DifEvo_McCormick	310	-1.91	-1.91322	-0.547614,-1.54759
DifEvo_McCormick	1129	-1.91	-1.91322	-0.547215,-1.54717
DifEvo_McCormick	1076	-1.91	-1.91322	-0.547614,-1.54759
DifEvo_BoxBetts	438	0	1.46276e-06	1.00262,9.96803,0.99759
DifEvo_BoxBetts	2588	0	6.92735e-09	0.999957,10.0018,1.00002
DifEvo_BoxBetts	438	0	1.46276e-06	1.00262,9.96803,0.99759
DifEvo_BoxBetts	12288	0	1.46276e-06	1.00262,9.96803,0.99759
DifEvo_BoxBetts	3919	0	2.16333e-09	0.999869,10.0014,1.00007
DifEvo_BoxBetts	568	0	1.46276e-06	1.00262,9.96803,0.99759
DifEvo_BoxBetts	2086	0	1.47751e-08	1.00046,9.99673,0.999703
DifEvo_BoxBetts	568	0	1.46276e-06	1.00262,9.96803,0.99759
DifEvo_BoxBetts	3093	0	1.1011e-08	0.999602,10.0026,1.00023
DifEvo_BoxBetts	3276	0	1.623e-08	0.999733,10.0009,1.0002
DifEvo_Paviani	15677	-45.7	-45.7784	9.34977,9.34997,9.35003,9.35181,9.3508,9.35024,9.3502,9.34811,9.34999,9.34946
DifEvo_Paviani	30721	-45.7	-45.7784	9.34977,9.34997,9.35003,9.35181,9.3508,9.35024,9.3502,9.34811,9.34999,9.34946
DifEvo_Paviani	15677	-45.7	-45.7784	9.34977,9.34997,9.35003,9.35181,9.3508,9.35024,9.3502,9.34811,9.34999,9.34946
DifEvo_Paviani	40960	-45.7	-45.7784	9.34977,9.34997,9.35003,9.35181,9.3508,9.35024,9.3502,9.34811,9.34999,9.34946
DifEvo_Paviani	40960	-45.7	-45.7784	9.34977,9.34997,9.35003,9.35181,9.3508,9.35024,9.3502,9.34811,9.34999,9.34946
DifEvo_Paviani	1890	-45.7	-45.7784	9.34977,9.34997,9.35003,9.35181,9.3508,9.35024,9.3502,9.34811,9.34999,9.34946
DifEvo_Paviani	40960	-45.7	-45.7784	9.34977,9.34997,9.35003,9.35181,9.3508,9.35024,9.3502,9.34811,9.34999,9.34946
DifEvo_Paviani	1890	-45.7	-45.7784	9.34977,9.34997,9.35003,9.35181,9.3508,9.35024,9.3502,9.34811,9.34999,9.34946
DifEvo_Paviani	40960	-45.7	-45.7784	9.34977,9.34997,9.35003,9.35181,9.3508,9.35024,9.3502,9.34811,9.34999,9.34946
DifEvo_Paviani	40960	-45.7	-45.7784	9.34977,9.34997,9.35003,9.35181,9.3508,9.35024,9.3502,9.34811,9.34999,9.34946
DifEvo_GoldsteinPrice	348	3	3	-4.53033e-05,-1.00002
DifEvo_GoldsteinPrice	927	3	3	-4.53033e-05,-1.00002
DifEvo_GoldsteinPrice	552	3	3	-4.53033e-05,-1.00002
DifEvo_GoldsteinPrice	8192	3	3	-4.53033e-05,-1.00002
DifEvo_GoldsteinPrice	1422	3	3	-4.53033e-05,-1.00002
DifEvo_GoldsteinPrice	372	3	3	-4.53033e-05,-1.00002
DifEvo_GoldsteinPrice	643	3	3	-4.53033e-05,-1.00002
DifEvo_GoldsteinPrice	270	3	3	-4.53033e-05,-1.00002
DifEvo_GoldsteinPrice	1094	3	3	-4.53033e-05,-1.00002
DifEvo_GoldsteinPrice	1213	3	3	-4.53033e-05,-1.00002
DifEvo_Shekel5	-1554	-10.1532	-5.10077	7.99949,7.99915,7.99954,7.99974
DifEvo_Shekel5	5313	-10.1532	-10.1532	3.9999,4.00035,4.00003,3.99989
DifEvo_Shekel5	-1554	-10.1532	-5.10077	7.99949,7.99915,7.99954,7.99974
DifEvo_Shekel5	-16384	-10.1532	-5.10077	7.99949,7.99915,7.99954,7.99974
DifEvo_Shekel5	9679	-10.1532	-10.1532	3.99995,4.00001,4.0002,4.00004
DifEvo_Shekel5	-691	-10.1532	-5.10077	7.99949,7.99915,7.99954,7.99974
DifEvo_Shekel5	4374	-10.1532	-10.1532	4.00032,4.00038,3.99997,4.00006
DifEvo_Shekel5	-691	-10.1532	-5.10077	7.99949,7.99915,7.99954,7.99974
DifEvo_Shekel5	7148	-10.1532	-10.1532	3.99984,4.00006,4.00011,4.00002
DifEvo_Shekel5	10167	-10.1532	-10.1532	3.99995,4.00025,4.00003,3.99991
DifEvo_Shekel7	-1901	-10.4029	-5.12882	7.99958,7.99983,7.99964,7.99975
DifEvo_Shekel7	5363	-10.4029	-10.4029	4.00052,4.00085,3.99944,3.99984
DifEvo_Shekel7	-1901	-10.4029	-5.12882	7.99958,7.99983,7.99964,7.99975
DifEvo_Shekel7	-16384	-10.4029	-5.12882	7.99958,7.99983,7.99964,7.99975
DifEvo_Shekel7	9135	-10.4029	-10.4029	4.00067,4.0005,3.99972,3.99951
DifEvo_Shekel7	-758	-10.4029	-5.12882	7.99958,7.99983,7.99964,7.99975
DifEvo_Shekel7	5179	-10.4029	-10.4029	4.00061,4.00057,3.99936,3.99978
DifEvo_Shekel7	-758	-10.4029	-5.12882	7.99958,7.99983,7.99964,7.99975
DifEvo_Shekel7	7584	-10.4029	-10.4029	4.00083,4.00088,3.99956,3.9994
DifEvo_Shekel7	12969	-10.4029	-10.4029	4.00059,4.00058,3.99944,3.99981
DifEvo_Shekel10	1694	-10.5364	-10.5364	4.00072,4.00045,3.99961,3.99969
DifEvo_Shekel10	5343	-10.5364	-10.5364	4.00057,4.00065,3.99974,3.99953
DifEvo_Shekel10	1694	-10.5364	-10.5364	4.00072,4.00045,3.99961,3.99969
DifEvo_Shekel10	-16384	-10.5364	-5.17564	7.99938,7.99918,7.99936,7.9992
DifEvo_Shekel10	12078	-10.5364	-10.5364	4.00066,4.00068,3.99952,3.99941
DifEvo_Shekel10	1100	-10.5364	-10.5364	4.00087,4.00084,3.99971,3.99933
DifEvo_Shekel10	4757	-10.5364	-10.5364	4.00071,4.00049,3.99967,3.99927
DifEvo_Shekel10	1100	-10.5364	-10.5364	4.00087,4.00084,3.99971,3.99933
DifEvo_Shekel10	8631	-10.5364	-10.5364	4.00125,4.00077,3.99957,3.9997
DifEvo_Shekel10	11297	-10.5364	-10.5364	4.00069,4.00051,3.99982,3.99979
DifEvo_Griewank	-567	0	0.0394593	12.559,-0.000449444
DifEvo_Griewank	1783	0	0.0073965	3.13953,-4.43728
DifEvo_Griewank	-627	0	0.0394589	12.5596,-0.000449444
DifEvo_Griewank	3284	0	6.68743e-09	4.9327e-05,0.000147848
DifEvo_Griewank	4406	0	1.3728e-08	5.27899e-05,-0.000222005
DifEvo_Griewank	1109	0	0.00739627	-3.14043,-4.43767
DifEvo_Griewank	1842	0	1.92659e-08	0.000116045,-0.000223757
DifEvo_Griewank	1109	0	0.00739627	-3.14043,-4.43767
DifEvo_Griewank	3070	0	2.56402e-08	-0.000225947,-2.01272e-05
DifEvo_Griewank	3778	0	8.50242e-10	-9.0804e-06,-5.68572e-05
DifEvo_SixHumpCamel	221	-1.03	-1.03163	0.0895984,-0.71269
DifEvo_SixHumpCamel	1381	-1.03	-1.03163	-0.0897643,0.712763
DifEvo_SixHumpCamel	221	-1.03	-1.03163	0.0895984,-0.71269
DifEvo_SixHumpCamel	8192	-1.03	-1.03163	0.0895984,-0.71269
DifEvo_SixHumpCamel	2521	-1.03	-1.03163	0.089959,-0.712517
DifEvo_SixHumpCamel	175	-1.03	-1.03163	0.0895984,-0.71269
DifEvo_SixHumpCamel	733	-1.03	-1.03163	0.0895984,-0.71269
DifEvo_SixHumpCamel	175	-1.03	-1.03163	0.0895984,-0.71269
DifEvo_SixHumpCamel	2017	-1.03	-1.03163	-0.0897924,0.712623
DifEvo_SixHumpCamel	3148	-1.03	-1.03163	0.0898212,-0.712675
DifEvo_Shubert	-407	-24.06	-21.5259	-0.491651,4.55775
DifEvo_Shubert	1110	-24.06	-24.0625	-0.491243,-0.491156
DifEvo_Shubert	-407	-24.06	-21.5259	-0.491651,4.55775
DifEvo_Shubert	8192	-24.06	-24.0625	-0.491334,-6.77432
DifEvo_Shubert	2326	-24.06	-24.0625	-6.77463,5.7921
DifEvo_Shubert	-328	-24.06	-21.5259	-0.491481,4.55745
DifEvo_Shubert	1732	-24.06	-24.0625	-0.491288,5.79178
DifEvo_Shubert	-407	-24.06	-21.5259	-0.491481,4.55745
DifEvo_Shubert	1890	-24.06	-24.0625	5.79176,-0.491392
DifEvo_Shubert	1662	-24.06	-24.0625	-6.77457,-6.77461
DifEvo_Hansen	678	-176.54	-176.542	4.97654,-1.4253
DifEvo_Hansen	2109	-176.54	-176.542	4.97671,-7.70823
DifEvo_Hansen	678	-176.54	-176.542	4.97654,-1.4253
DifEvo_Hansen	8192	-176.54	-176.542	-1.30678,4.85818
DifEvo_Hansen	4114	-176.54	-176.542	4.97646,-1.42513
DifEvo_Hansen	467	-176.54	-176.542	4.97637,-1.42542
DifEvo_Hansen	755	-176.54	-176.542	4.97629,-1.42526
DifEvo_Hansen	467	-176.54	-176.542	4.97637,-1.42542
DifEvo_Hansen	4232	-176.54	-176.542	-1.30671,-7.70831
DifEvo_Hansen	4226	-176.54	-176.542	-7.58989,-1.42513
DifEvo_Ackley	792	0	2.39448e-06	-3.12293e-07,7.86865e-07
DifEvo_Ackley	1379	0	4.43822e-06	1.43402e-06,6.36968e-07
DifEvo_Ackley	792	0	2.39448e-06	-3.12293e-07,7.86865e-07
DifEvo_Ackley	1931	0	2.04394e-06	-6.41251e-07,3.33164e-07
DifEvo_Ackley	2290	0	3.24083e-06	3.19246e-07,1.10042e-06
DifEvo_Ackley	573	0	2.72339e-06	-9.20229e-07,2.83322e-07
DifEvo_Ackley	909	0	1.38965e-05	2.78327e-06,-4.04848e-06
DifEvo_Ackley	573	0	2.72339e-06	-9.20229e-07,2.83322e-07
DifEvo_Ackley	1831	0	1.09558e-06	1.40256e-07,3.6106e-07
DifEvo_Ackley	1737	0	5.07638e-06	2.29979e-07,-1.77995e-06
DifEvo_Bohachevsky1	284	0	8.90415e-07	-0.000165047,-0.000122046
DifEvo_Bohachevsky1	1407	0	8.90415e-07	-0.000165047,-0.000122046
DifEvo_Bohachevsky1	284	0	8.90415e-07	-0.000165047,-0.000122046
DifEvo_Bohachevsky1	1815	0	3.5449e-07	5.96988e-05,9.50558e-05
DifEvo_Bohachevsky1	1858	0	8.90415e-07	-0.000165047,-0.000122046
DifEvo_Bohachevsky1	286	0	8.90415e-07	-0.000165047,-0.000122046
DifEvo_Bohachevsky1	1244	0	8.90415e-07	-0.000165047,-0.000122046
DifEvo_Bohachevsky1	286	0	8.90415e-07	-0.000165047,-0.000122046
DifEvo_Bohachevsky1	1777	0	8.90415e-07	-0.000165047,-0.000122046
DifEvo_Bohachevsky1	2420	0	1.20855e-08	2.41911e-05,-1.05007e-05
DifEvo_Bohachevsky2	284	0	8.27817e-07	-0.000233765,-4.18875e-05
DifEvo_Bohachevsky2	1355	0	2.85184e-07	8.5379e-05,-8.38889e-05
DifEvo_Bohachevsky2	284	0	8.27817e-07	-0.000233765,-4.18875e-05
DifEvo_Bohachevsky2	1685	0	8.27817e-07	-0.000233765,-4.18875e-05
DifEvo_Bohachevsky2	2274	0	2.44601e-08	-9.54606e-07,3.08501e-05
DifEvo_Bohachevsky2	396	0	8.27817e-07	-0.000233765,-4.18875e-05
DifEvo_Bohachevsky2	1292	0	8.27817e-07	-0.000233765,-4.18875e-05
DifEvo_Bohachevsky2	740	0	8.27817e-07	-0.000233765,-4.18875e-05
DifEvo_Bohachevsky2	1458	0	8.27817e-07	-0.000233765,-4.18875e-05
DifEvo_Bohachevsky2	1845	0	5.39408e-07	-1.68607e-05,-0.000144363
DifEvo_Bohachevsky3	284	0	3.83822e-07	-0.000199647,2.95373e-05
DifEvo_Bohachevsky3	1548	0	1.43438e-08	-4.07221e-05,4.88246e-05
DifEvo_Bohachevsky3	1506	0	3.83822e-07	-0.000199647,2.95373e-05
DifEvo_Bohachevsky3	1967	0	1.43935e-07	0.00020853,-0.000190635
DifEvo_Bohachevsky3	2575	0	9.87139e-08	0.00020864,-0.000124541
DifEvo_Bohachevsky3	556	0	3.83822e-07	-0.000199647,2.95373e-05
DifEvo_Bohachevsky3	1466	0	1.73513e-07	0.000286824,-0.000183214
DifEvo_Bohachevsky3	462	0	3.83822e-07	-0.000199647,2.95373e-05
DifEvo_Bohachevsky3	2047	0	6.10196e-08	-0.000131821,0.000122747
DifEvo_Bohachevsky3	2358	0	3.54981e-08	-6.07921e-05,7.50413e-05
DifEvo_Easom	-781	-1	-8.10999e-05	4.9764,1.30555
DifEvo_Easom	1154	-1	-0.999999	3.14203,3.14212
DifEvo_Easom	-991	-1	-8.10999e-05	4.9764,1.30555
DifEvo_Easom	8192	-1	-1	3.1412,3.14122
DifEvo_Easom	3131	-1	-1	3.14202,3.14123
DifEvo_Easom	-843	-1	-8.11007e-05	4.97884,1.30639
DifEvo_Easom	1094	-1	-1	3.14166,3.1419
DifEvo_Easom	-925	-1	-8.11007e-05	4.97884,1.30639
DifEvo_Easom	1652	-1	-0.999999	3.14123,3.14213
DifEvo_Easom	2570	-1	-1	3.14188,3.14181
DifEvo_Rastrigin	774	0	1.79606e-07	2.85961e-05,9.35771e-06
DifEvo_Rastrigin	1404	0	4.53968e-07	-4.69884e-05,-8.96251e-06
DifEvo_Rastrigin	774	0	1.79606e-07	2.85961e-05,9.35771e-06
DifEvo_Rastrigin	1815	0	1.09824e-06	-7.06252e-05,-2.34045e-05
DifEvo_Rastrigin	2627	0	2.83246e-07	-1.76469e-05,3.3411e-05
DifEvo_Rastrigin	-338	0	0.99496	9.19171e-06,0.994903
DifEvo_Rastrigin	1201	0	5.30971e-07	-4.61495e-05,-2.33794e-05
DifEvo_Rastrigin	-338	0	0.99496	9.19171e-06,0.994903
DifEvo_Rastrigin	1970	0	3.3732e-07	3.19813e-05,2.60282e-05
DifEvo_Rastrigin	1934	0	6.96922e-07	4.99665e-06,-5.90583e-05
DifEvo_Michalewicz2	264	-1.8013	-1.8013	2.20311,1.57076
DifEvo_Michalewicz2	627	-1.8013	-1.8013	2.20311,1.57076
DifEvo_Michalewicz2	264	-1.8013	-1.8013	2.20311,1.57076
DifEvo_Michalewicz2	8192	-1.8013	-1.8013	2.20311,1.57076
DifEvo_Michalewicz2	1006	-1.8013	-1.8013	2.20311,1.57076
DifEvo_Michalewicz2	471	-1.8013	-1.8013	2.20311,1.57076
DifEvo_Michalewicz2	712	-1.8013	-1.8013	2.20311,1.57076
DifEvo_Michalewicz2	471	-1.8013	-1.8013	2.20311,1.57076
DifEvo_Michalewicz2	767	-1.8013	-1.8013	2.20311,1.57076
DifEvo_Michalewicz2	805	-1.8013	-1.8013	2.20311,1.57076
DifEvo_Michalewicz5	3861	-4.68766	-4.68765	2.20298,1.57089,1.28517,1.92313,1.72036
DifEvo_Michalewicz5	9434	-4.68766	-4.68766	2.20301,1.571,1.28495,1.92308,1.72051
DifEvo_Michalewicz5	3861	-4.68766	-4.68765	2.20298,1.57089,1.28517,1.92313,1.72036
DifEvo_Michalewicz5	-20480	-4.68766	-4.53765	2.20279,1.57071,1.28519,1.92312,0.996674
DifEvo_Michalewicz5	18738	-4.68766	-4.68765	2.20317,1.57092,1.28517,1.92312,1.72048
DifEvo_Michalewicz5	-1406	-4.68766	-4.53765	2.20279,1.57071,1.28519,1.92312,0.996674
DifEvo_Michalewicz5	9708	-4.68766	-4.68765	2.20278,1.57051,1.285,1.92308,1.72048
DifEvo_Michalewicz5	-1406	-4.68766	-4.53765	2.20279,1.57071,1.28519,1.92312,0.996674
DifEvo_Michalewicz5	19626	-4.68766	-4.68766	2.20254,1.57073,1.28505,1.92307,1.7205
DifEvo_Michalewicz5	-12415	-4.68766	-4.53765	2.20279,1.57071,1.28519,1.92312,0.996674
DifEvo_Michalewicz10	28315	-9.66015	-9.66014	2.20289,1.5708,1.28497,1.92312,1.72035,1.57081,1.45433,1.75614,1.65571,1.5708
DifEvo_Michalewicz10	39942	-9.66015	-9.5765	2.20304,1.57104,1.28492,1.92313,1.72055,1.57085,1.45431,1.36052,1.28277,1.85849
DifEvo_Michalewicz10	28315	-9.66015	-9.66014	2.20289,1.5708,1.28497,1.92312,1.72035,1.57081,1.45433,1.75614,1.65571,1.5708
DifEvo_Michalewicz10	-40960	-9.66015	-8.63908	2.20304,1.57092,1.95495,1.92304,1.72042,1.57076,1.45442,1.75616,1.65574,1.21697
DifEvo_Michalewicz10	-23687	-9.66015	-8.63908	2.20304,1.57092,1.95495,1.92304,1.72042,1.57076,1.45442,1.75616,1.65574,1.21697
DifEvo_Michalewicz10	11619	-9.66015	-9.59815	2.20369,1.57043,1.2846,1.9231,1.7204,1.57078,1.45442,1.75604,1.65567,1.21699
DifEvo_Michalewicz10	-37887	-9.66015	-9.18803	2.20284,1.57051,2.21942,1.92308,1.72044,1.57082,2.22106,1.75613,1.95893,1.85853
DifEvo_Michalewicz10	11619	-9.66015	-9.59815	2.20369,1.57043,1.2846,1.9231,1.7204,1.57078,1.45442,1.75604,1.65567,1.21699
DifEvo_Michalewicz10	-40960	-9.66015	-8.63908	2.20304,1.57092,1.95495,1.92304,1.72042,1.57076,1.45442,1.75616,1.65574,1.21697
DifEvo_Michalewicz10	-40960	-9.66015	-8.63908	2.20304,1.57092,1.95495,1.92304,1.72042,1.57076,1.45442,1.75616,1.65574,1.21697
==11191== 
==11191== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 15 from 1)
==11191== malloc/free: in use at exit: 0 bytes in 0 blocks.
==11191== malloc/free: 8,766 allocs, 8,766 frees, 582,556 bytes allocated.
==11191== For counts of detected errors, rerun with: -v
==11191== All heap blocks were freed -- no leaks are possible.

*/
#endif
