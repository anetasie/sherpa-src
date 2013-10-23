#ifndef Simplex_hh
#define Simplex_hh

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


#include <sherpa/Array2d.hh>

namespace sherpa {

  class Simplex : public sherpa::Array2d< double > {

  public:

    Simplex( int r=0, int c=0 ) : sherpa::Array2d< double >( r, c ),
				  npar( c - 1 ) { }

    static double calc_standard_deviation_square( int num, const double* ptr );

    bool check_convergence( double tolerance, double tol_sqr,
			    std::vector< double >& fctvals, 
			    int finalsimplex=0 );

    void print_simplex( ) const;

    void print_vertex( int vertex ) const;

    void sort( );

  private:


    const int npar;

    bool are_fct_vals_close_enough( double tolerance ) const;

    bool is_max_length_small_enough( double tol ) const;

    bool is_stddev_small_enough( std::vector< double >& fctvals,
				 double tolerance, double tol_sqr );

  };                                                           // class Simplex

}                                                                  // namespace

#endif
