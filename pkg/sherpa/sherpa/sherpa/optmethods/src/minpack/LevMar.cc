#ifdef testLevMar

#include "levmar.hh"

using namespace minpack;


#include "fcmp.h"
#include "tstoptfct.hh"
#include "sherpa/functor.hh"

template< typename Real >
void print_pars( const char* name, int nfev, Real fval, Real answer,
                int n, const std::vector< Real >& x,
		 const std::vector< Real >& err,
                Real tol=1.0e4*std::sqrt( std::numeric_limits< Real >::epsilon()
 ),
                const char* prefix="lmdif_" ) {


  std::cout << prefix << name << '\t';
  if ( 0 == _sao_fcmp( fval, answer, std::sqrt(tol) ) )
    std::cout << nfev << '\t';
  else
    std::cout << -nfev << '\t';
  std::cout << answer << '\t';
  std::cout << fval << '\t';
  std::cout << x[0];
  for ( int ii = 1; ii < n; ++ii )
    std::cout << ',' << x[ii];
  std::cout << '\t';
  std::cout << err[0];
  for ( int ii = 1; ii < n; ++ii )
    std::cout << ',' << err[ii];
  std::cout << '\n';

}

template< typename Init, typename Fct >
void justdoit( Init init, Fct fct, int npars, 
	       std::vector< double >& pars, std::vector< double >& lo,
	       std::vector< double >& hi,std::vector<double> covarerr,
	       double tol, const char* header ) {

  try {

      int mfcts;
      double answer;
      init( npars, mfcts, answer, &pars[0], &lo[0], &hi[0] );

      minpack::LevMar< Fct, void* > lm( npars, &pars[0], &lo[0], &hi[0],
					fct, NULL, mfcts );

      int maxnfev=128*npars, nfev;
      double fmin=0.0;

      int nprint = 0;
      double epsfcn = 1.0e-8, factor=100.0;

      lm( &pars[0], tol, tol, tol, maxnfev, epsfcn, factor,
	  nprint, nfev, fmin, &covarerr[0] );

      // lm.minimize( &pars[0], tol, maxnfev, nfev, fmin );
      
      print_pars( header, nfev, fmin, answer, npars, pars, covarerr );
      
  } catch( const sherpa::OptErr& oe ) {
    
    std::cerr << oe << '\n';
    
  }

  return;

}

int main( int argc, char* argv[] ) {

  int npars=16;
  if ( argc == 2 )
    npars = atoi( argv[1] );

  if ( npars % 2 || npars < 2 ) {
    printf( "The minimum value for the free parameter must be an even "
            "and it is greater then 2\n" );
    return EXIT_FAILURE;
  }

  std::cout << "#\n# sizeof(double) = " << sizeof( double ) << "\n";
  std::cout << "# A negative value for the nfev signifies that the "
    "optimization method did not converge\n#\n";
  std::cout << "name\tnfev\tanswer\tfval\tpars\terr\nS\tN\tN\tN\tN\tN\n";

  double tol = std::sqrt( std::numeric_limits< double >::epsilon() );

  // The factor of 8 for the parameter is because the Chebyquad has 9
  // free parameters and npars has to be at least 2.
  // The factor of 32 for fvec is because the number of functions may
  // be greater then the number of free parameters.
  std::vector< double > pars( 8 * npars ), lo( 8 * npars ), hi( 8 * npars ),
    covarerr( 8 * npars );

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Rosenbrock<double,void*> );
    
    justdoit( sherpa::fct_ptr( tstoptfct::RosenbrockInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "Rosenbrock" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::FreudensteinRoth<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::FreudensteinRothInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "FreudensteinRoth" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::PowellBadlyScaled<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::PowellBadlyScaledInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "PowellBadlyScaled" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::BrownBadlyScaled<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::BrownBadlyScaledInit<double> ),
	      fct, 2, pars, lo, hi, covarerr, tol, "BrownBadlyScaled" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Beale<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::BealeInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "Beale" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::JennrichSampson<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::JennrichSampsonInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "JennrichSampson" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::HelicalValley<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::HelicalValleyInit<double> ),
	      fct, 3*npars, pars, lo, hi, covarerr, tol, "HelicalValley" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Bard<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::BardInit<double> ),
	      fct, 3*npars, pars, lo, hi, covarerr, tol, "Bard" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Gaussian<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::GaussianInit<double> ),
	      fct, 3, pars, lo, hi, covarerr, tol, "Gaussian" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Meyer<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::MeyerInit<double> ),
	      fct, 3, pars, lo, hi, covarerr, tol, "Meyer" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::GulfResearchDevelopment<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::GulfResearchDevelopmentInit<double> ),
	      fct, 3, pars, lo, hi, covarerr, tol, "GulfResearchDevelopment" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Box3d<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::Box3dInit<double> ),
	      fct, 3, pars, lo, hi, covarerr, tol, "Box3d" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::PowellSingular<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::PowellSingularInit<double> ),
	      fct, 4*npars, pars, lo, hi, covarerr, tol, "PowellSingular" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Wood<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::WoodInit<double> ),
	      fct, 4, pars, lo, hi, covarerr, tol, "Wood" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::KowalikOsborne<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::KowalikOsborneInit<double> ),
	      fct, 4, pars, lo, hi, covarerr, tol, "KowalikOsborne" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::BrownDennis<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::BrownDennisInit<double> ),
	      fct, 4, pars, lo, hi, covarerr, tol, "BrownDennis" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Osborne1<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::Osborne1Init<double> ),
	      fct, 5, pars, lo, hi, covarerr, tol, "Osborne1" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Biggs<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::BiggsInit<double> ),
	      fct, 6, pars, lo, hi, covarerr, tol, "Biggs" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Osborne2<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::Osborne2Init<double> ),
	      fct, 11, pars, lo, hi, covarerr, tol, "Osborne2" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Watson<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::WatsonInit<double> ),
	      fct, 6, pars, lo, hi, covarerr, tol, "Watson" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::PenaltyI<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::PenaltyIInit<double> ),
	      fct, 4, pars, lo, hi, covarerr, tol, "PenaltyI" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::PenaltyII<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::PenaltyIIInit<double> ),
	      fct, 4, pars, lo, hi, covarerr, tol, "PenaltyII" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::VariablyDimensioned<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::VariablyDimensionedInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "VariablyDimensioned" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Trigonometric<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::TrigonometricInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "Trigonometric" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::BrownAlmostLinear<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::BrownAlmostLinearInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "BrownAlmostLinear" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::DiscreteBoundary<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::DiscreteBoundaryInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "DiscreteBoundary" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::DiscreteIntegral<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::DiscreteIntegralInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "DiscreteIntegral" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::BroydenTridiagonal<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::BroydenTridiagonalInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "BroydenTridiagonal" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::BroydenBanded<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::BroydenBandedInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "BroydenBanded" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::LinearFullRank<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::LinearFullRankInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "LinearFullRank" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::LinearFullRank1<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::LinearFullRank1Init<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "LinearFullRank1" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::LinearFullRank0cols0rows<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::LinearFullRank0cols0rowsInit<double> ),
	      fct, npars, pars, lo, hi, covarerr, tol, "LinearFullRank0cols0rows" );
  }

  {
    sherpa::FctPtr< void, int, int, double*, double*, int&, void* >
      fct( tstoptfct::Chebyquad<double,void*> );

    justdoit( sherpa::fct_ptr( tstoptfct::ChebyquadInit<double> ),
	      fct, 9, pars, lo, hi, covarerr, tol, "Chebyquad" );
  }

  return 0;

}

#endif
