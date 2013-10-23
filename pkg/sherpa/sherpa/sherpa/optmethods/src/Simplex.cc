#include <cmath>

#include <fcmp.h>

#include "Simplex.hh"

using namespace sherpa;

bool Simplex::are_fct_vals_close_enough( double tolerance ) const {

  if ( 0 == _sao_fcmp( this->operator()( 0, npar ), 
		       this->operator()( nrows-1, npar ),
		       tolerance ) )
    return true;

  return false;

}

double Simplex::calc_standard_deviation_square( int num, const double* ptr ) {

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

}

bool Simplex::check_convergence( double tolerance, double tol_sqr,
				 std::vector< double >& fctvals,
				 int finalsimplex ) {

  switch( finalsimplex ) {
  case 0:
    if ( false == is_max_length_small_enough( tolerance ) )
      return false;
    return true;
  case 2:
    {
      if ( false == is_max_length_small_enough( tolerance ) )
	return false;
      bool stddev = is_stddev_small_enough( fctvals, tolerance, tol_sqr );
      bool fctval = are_fct_vals_close_enough( tolerance );
      return stddev && fctval;
    }
    break;
  default:
    {
      if ( false == is_max_length_small_enough( tolerance ) )
	return false;
      bool stddev = is_stddev_small_enough( fctvals, tolerance, tol_sqr );
      bool fctval = are_fct_vals_close_enough( tolerance );
      return stddev || fctval;
    }
  }

  return false;

}                                                          // check_convergence

bool Simplex::is_max_length_small_enough( double tol )  const {

  const int index_smallest = 0;

  double maxof_x_i_minus_x_min = -1.0; // norm is always a positive number.
  for ( int ii = 0; ii <= npar; ++ii ) {
    double tmp = 0.0;
    if ( ii != index_smallest )
      for ( int jj = 0; jj < npar; ++jj )
	tmp += ( this->operator( )( ii, jj ) -
		 this->operator( )( index_smallest, jj ) ) *
	  ( this->operator( )( ii, jj ) -
	    this->operator( )( index_smallest, jj ) );

    maxof_x_i_minus_x_min = std::max( maxof_x_i_minus_x_min, tmp );
  }
  double norm_min = 0.0;
  for ( int ii = 0; ii < npar; ++ii )
    norm_min +=
      this->operator( )( index_smallest, ii ) *
      this->operator( )( index_smallest, ii );
  norm_min = norm_min > 1.0 ? norm_min : 1.0;
  if ( maxof_x_i_minus_x_min <= tol * norm_min )
    return true;

  return false;

}                                                 // is_max_length_small_enough

bool Simplex::is_stddev_small_enough( std::vector< double >& fctvals,
				      double tolerance, double tol_sqr ) {

  this->copy_col( npar, fctvals );
  double std_dev_sqr =
    calc_standard_deviation_square( npar + 1 , &fctvals[0] );
  if ( _sao_fcmp( std_dev_sqr, tol_sqr, tolerance ) <= 0 )
    return true;

  return false;

}                                                     // is_stddev_small_enough

void Simplex::print_simplex( ) const {

  for ( int ii = 0; ii <= npar; ++ii )
    print_vertex( ii );

}                                                              // print_simplex

void Simplex::print_vertex( int vertex ) const {

  fprintf( stdout, "\tf" );
  if ( 0 == vertex )
    fprintf( stdout, "'" );
  fprintf( stdout, "( %.10e", this->operator( )( vertex, 0 ) );
  for ( int ii = 1; ii < npar; ++ii )
    fprintf( stdout, ", %.10e", this->operator( )( vertex, ii ) );
  fprintf( stdout, " ) = %.10e\n", this->operator( )( vertex, npar ) );

}                                                               // print_vertex

void Simplex::sort( ) {

  // key should be stored as a member of the class to minimize de/alloc
  std::vector< double > key( ncols );
  // key should be stored as a member of the class to minimize de/alloc

  const int fvalindex = ncols - 1;

  for ( int jj = 1; jj < nrows; ++jj ) {

    for ( int ii = 0; ii < ncols; ++ii )
      key[ ii ] = this->operator( )( jj, ii );

    int ii = jj;
    for ( ; ii > 0 &&
	    this->operator( )( ii - 1, fvalindex ) > key[ fvalindex ];
	  --ii ) {

      for ( int kk = 0; kk < ncols; ++kk )
	this->operator( )( ii, kk ) =
	  this->operator( )( ii - 1, kk );

    }

    for ( int kk = 0; kk < ncols; ++kk )
      this->operator( )( ii, kk ) = key[ kk ];

  }

}                                                                       // sort
