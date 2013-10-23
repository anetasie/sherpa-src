#ifndef Opt_hh
#define Opt_hh

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


#include<cstdlib>
#include <iostream>
#include <limits>
#include <vector>
#include <stdexcept>

namespace sherpa {

  class OptErr {

    friend std::ostream& operator << ( std::ostream& os, const OptErr& opte ) {
      return opte.print( os ); }

  public:

    enum Err { Input=1, OutOfBounds, MaxFev, UsrFunc };

    OptErr( OptErr::Err e ) : err( e ) { }

    Err err;

  private:

    std::ostream& print( std::ostream& os ) const {

      char* msg[] = {
	"Input Err",
	"Par is out of bound",
	"Max number of function evaluation",
	"User Function error"
      };

      os << msg[ err ];

      return os;

    }

  };

  template< typename Func, typename Data >
  class Opt {

    friend std::ostream& operator << ( std::ostream& os, const Opt& opt ) {
      return opt.print( os ); }

  public:

    virtual ~Opt( ) { }
    
    Opt( int numpar, double* par, const double* lo, const double* hi,
	 Func func, Data data ) 
      : npar( numpar ), lo_bound( lo, lo + npar ), hi_bound( hi, hi + npar ),
	usr_func( func ), usr_data( data ) {

      if ( true == are_pars_outside_limits( par ) )
	throw sherpa::OptErr( OptErr::OutOfBounds );
      
    }

    std::ostream& print( std::ostream& os, double* par, double fstat=0,
			 int nfev=-1 ) const {

      if ( nfev > 0 )
	os << "nfev=" << nfev << ": ";
      print_par( os, par );
      if ( nfev > 0 )
	os << "=" << fstat << '\n';

      return os;

    }

    std::ostream& print_par( std::ostream& os, double* par ) const {

      os << "f( " << par[0];
      for ( int ii = 1; ii < npar; ++ii )
	os << ", " << par[ ii ];
      os << ")";
      return os;

    }

    virtual int eval_user_func( int maxnfev, double* mypar, double& fval,
				int& nfev, int& ierr ) {
      std::cerr << "Please define me!" << std::endl;
      return EXIT_FAILURE;
    }

  protected:

    const int npar;
    const std::vector< double > lo_bound;
    const std::vector< double > hi_bound;
    Func usr_func;
    Data usr_data;

    // If lo <= xpar <= hi return false, else return true
    bool are_pars_outside_limits( const double* pars ) const {
      for ( int ii = 0; ii < npar; ++ii )
	if ( pars[ ii ] < lo_bound[ ii ] || pars[ ii ] > hi_bound[ ii ] )
	  return true;
      return false;

    }
      
    double calc_standard_deviation_square( int num, const double* ptr ) const {

      //
      // The standard deviation algorithm is due to Donald E. Knuth (1998).
      // The art of Computer Programming, volume2: Seminumerical Algorithms,
      // 3rd edn., p 232
      //
      double mean = 0.0, stddev = 0.0;
      for ( int ii = 0; ii < num; ++ii ) {
	double delta = ptr[ ii ] - mean;
	mean += delta / double( ii + 1 );
	stddev += delta * ( ptr[ ii ] - mean );
      }

      if ( 1 != num )
	stddev /= double( num - 1);

      return stddev;

    }                                         // calc_standard_deviation_square

  };                                                               // class Opt

  template< typename Func, typename Data >
  class OptFunc : public sherpa::Opt< Func, Data > {

  public:

    OptFunc( int numpar, double* par, const double* lo, const double* hi,
	     Func func, Data data, int mfcts=0 ) : 
      sherpa::Opt<Func,Data>( numpar, par, lo, hi, func, data ) { }

    int eval_user_func( int maxnfev, double* mypar, double& fval,
			int& nfev, int& ierr ) {
      
      if ( true ==
	   Opt<Func,Data>::are_pars_outside_limits( mypar ) ) {
	fval = std::numeric_limits< double >::max( );
	return EXIT_SUCCESS;
      }

      ++nfev;
      
      Opt<Func,Data>::usr_func( Opt<Func,Data>::npar, mypar, fval, ierr,
				Opt<Func,Data>::usr_data );
      // this->print( std::cout, mypar, fval, Opt<Func,Data>::nfev );

      if ( EXIT_SUCCESS != ierr )
	return ierr;
      if ( nfev >= maxnfev ) {
	ierr = sherpa::OptErr::MaxFev;
      }

      return ierr;

    }                                                         // eval_user_func

    virtual int minimize( double* model_par, double tol, int maxnfev,
			  int& nfev, double& fmin ) {
      int ierr = EXIT_SUCCESS;
      return eval_user_func( maxnfev, model_par, fmin, nfev, ierr );
    }

  };                                                           // class OptFunc

#define NPAR  sherpa::Opt<Func,Data>::npar

}                                                           // namespace sherpa

#endif
