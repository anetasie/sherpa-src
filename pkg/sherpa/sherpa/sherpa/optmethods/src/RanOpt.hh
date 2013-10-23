#ifndef RanOpt_hh
#define RanOpt_hh

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


#include "mt19937ar.h"

#include "Opt.hh"

namespace sherpa {

  //
  // The base class for the global optimization class.
  //
  template< typename Func, typename Data >
  class RanOpt : public sherpa::Opt< Func, Data > {

  public:

    RanOpt( int numpar, double* par, const double* lo, const double* hi,
	    Func func, Data data, int seed )
      : sherpa::Opt< Func, Data >( numpar, par, lo, hi, func, data ) {
      
	init_genrand( seed );
	
    }

  protected:

    // return a pseudo random number in the [ 0, 1 ] real-interval
    double random_number( ) const { return genrand_real1( ); }

    // return a pseudo random number in the ( 0, 1 ) real-interval
    double random_number_ee( ) const { return genrand_real3( ); }

    // return a pseudo random number in the [ low, high ] integer-interval
    int random_number( int low, int high ) const {
      int answer = low + int( (high - low + 1) * random_number( ) );
      return answer;
    }

    // return a pseudo random number in the ( lo, high ) real-interval
    double random_number_ee( double low, double high ) const { 
      return low + ( high - low ) * random_number_ee( );
    }

  };                         // class RanOpt : public sherpa::Opt< Func, Data >

}                                                                  // namespace

#endif
