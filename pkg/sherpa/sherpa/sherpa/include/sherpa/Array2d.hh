#ifndef Array2d_hh
#define Array2d_hh

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


#include <iostream>
#include <stdexcept>
#include <vector>

//
// A simple 2d array class written for NelderMead/MultiDirSearch.
//
// For heavy duty numerical computations, the user is advised to
// check out the following three packages:
//
//    http://math.nist.gov/tnt/index.html
//    http://www.oonumerics.org/blitz/
//    http://boost.org/libs/multi_array/doc/user.html
//
//
namespace sherpa {

  template < typename Real >
  class Array2d {

    friend std::ostream& operator << ( std::ostream& os,
				       const Array2d< Real >& a ) {
      return a.print( os );
    }

  public:

    class Row {

      Array2d& myarray2d;
      int const row;

    public:

      Row( Array2d& array2d, int r ) : myarray2d( array2d ), row( r ) { }
      Real& operator [ ] ( int col ) { return myarray2d( row, col ); }
      const Real& operator [ ] ( int col ) const {
	return myarray2d( row, col ); }

    };

    Array2d( int r=0, int c=0 ) : nrows( r ), ncols( c ), myvector( r * c ),
				  row_index( r ) {
      generate_row_index( r, c );
    }

    Real& operator ( ) ( int r, int c ) { return myvector[ r * ncols + c ]; }
    const Real& operator ( ) ( int r, int c ) const {
      // return myvector[ r * ncols + c ];
      return myvector[ row_index[ r ] + c ]; }

    Row operator [ ] (int row) { return Row( *this, row ); }
    const Row operator [ ] (int row) const { return Row( *this, row ); }

    // Deep copy
    void copy( const Array2d< Real >& arg ) {

      if ( nrows != arg.get_nrows( ) || ncols != arg.get_ncols( ) )
	resizeme( arg.get_nrows( ), arg.get_ncols( ) );

      for ( int ii = 0; ii < nrows; ++ii )
	for ( int jj = 0; jj < ncols; ++jj )
	  this->operator( )( ii, jj ) = arg( ii, jj );

    }

    // copy the c-th col of the Array2d to the pointer 'cpto'
    void copy_col( int c, Real* cpto ) {

      if ( c < 0 || c >= ncols )
	throw std::runtime_error( "index out of bounds" );
      for ( int ii = 0; ii < nrows; ++ii )
	cpto[ ii ] = this->operator[]( ii )[ c ];

    }
    void copy_col( int c, std::vector< Real >& cpto ) {
      copy_col( c, &cpto[ 0 ] );
    }

    // copy from the pointer 'cpfrom' to the c-th col of the array
    void copy_col( const Real* cpfrom, int c ) {

      if ( c < 0 || c >= ncols )
	throw std::runtime_error( "index out of bounds" );
      for ( int ii = 0; ii < nrows; ++ii )
	this->operator[]( ii )[ c ] = cpfrom[ ii ];

    }
    void copy_col( const std::vector< Real >& cpfrom, int c ) {
      copy_col( &cpfrom[ 0 ], c );
    }

    // copy the r-th row of the Array2d to the pointer 'cpto'
    void copy_row( int r, Real* cpto ) {

      if ( r < 0 || r >= nrows )
	throw std::runtime_error( "index out of bounds" );
      for ( int ii = 0; ii < ncols; ++ii )
	cpto[ ii ] = this->operator[]( r )[ ii ];

    }
    void copy_row( int r, std::vector< Real >& cpto ) {
      copy_row( r, &cpto[ 0 ] );
    }

    // Copy from the pointer 'cpfrom' to the r-th row of the array
    void copy_row( const Real* cpfrom, int r ) {

      if ( r < 0 || r >= nrows )
	throw std::runtime_error( "index out of bounds" );
      for ( int ii = 0; ii < ncols; ++ii )
	this->operator[]( r )[ ii ] = cpfrom[ ii ];

    }
    void copy_row( const std::vector< Real >& cpfrom, int r ) {
      copy_row( &cpfrom[ 0 ], r );
    }

    int get_ncols( void ) const { return ncols; }
    int get_nrows( void ) const { return nrows; }

    std::ostream& print( std::ostream& os ) const {

      for ( int ii = 0; ii < nrows; ++ii ) {
	os << this->operator( )( ii, 0 );
	for ( int jj = 1; jj < ncols; ++jj )
	  os << ' ' << this->operator( )( ii, jj );
	if ( nrows - 1 != ii )
	  os << '\n';
      }
      return os;
    }

  protected:

    int nrows;
    int ncols;
    std::vector< Real > myvector;
    std::vector< int > row_index;

    void resizeme( int nrow, int ncol ) {
      nrows = nrow;
      ncols = ncol;
      myvector.resize( nrows * ncols );
      generate_row_index( nrow, ncol );
    }

    void generate_row_index( int nrow, int ncol ) {
      row_index.resize( nrow );
      for ( int ii = 0; ii < nrow; ++ii )
	row_index[ ii ] = ii * ncol;
    }

  private:
    // purposedly declare but do not define
    Array2d& operator = (Array2d const&);
    Array2d( Array2d const& );

  };

}

#endif

#ifdef testArray2d

#include <iostream>

#include "ArrayNd.hh"
#include "StopWatch.hh"

template < typename Type >
void initLaplacian( Type& uu, int r, int c ) {

  for ( int ii = 0; ii < r; ++ii ) {
    for ( int jj = 0; jj < c; ++jj )
      uu[ ii ][ jj ] = ii * ii + jj * jj;
  }

}

template < typename Type >
void timeme( Type& uu, Type& laplacian, int r, int c, const char* header ) {

  StopWatch stopwatch( header );

  initLaplacian( uu, r, c );

  for( int kk=0; kk < 10; kk++ ) {
    for( int ii = 1; ii < r - 1; ii++ ) {
      for( int jj = 1; jj < c - 1; jj++ ) {
        laplacian[ ii ][ jj ] = - uu[ ii - 1 ][ jj ] - uu[ ii ][ jj - 1 ] +
          4.0 * uu[ ii ][ jj ] - uu[ ii ][ jj + 1 ] - uu[ ii + 1 ][ jj ];
      }
    }
  }

  double sum = 0.0;
  for ( int ii = 1; ii < r - 1; ++ii )
    for ( int jj = 1; jj < c - 1; jj++ )
      sum += laplacian[ ii ][ jj ];

  fprintf( stdout, "%.12f\t", sum );

}

template < typename Type >
Type*** alloc3d( int nx, int ny, int nz ) {

  Type ***ptr = new int** [ nx ];

  for( int x = 0; x < nx; ++x ) {
    ptr[ x ] = new Type* [ ny ];

    for( int y = 0; y < ny; ++ y )
      ptr[ x ][ y ] = new Type[ nz ];
  }

  return ptr;

}

template < typename Type >
void del3d( Type*** ptr, int nx, int ny ) {
  for( int x = 0; x < nx; ++x ) {
    for( int y = 0; y < ny; ++y )
      delete [] ptr[ x ][ y ], ptr[ x ][ y ] = 0;
    delete [] ptr[ x ], ptr[ x ] = 0;
  }
  delete [] ptr;
  ptr = 0;
}

template < typename Type >
Type** alloc2d( int r, int c ) {
  Type** ptr = new Type*[ r ];
  for ( int ii = 0; ii < r; ++ii )
    ptr[ ii ] = new Type[ c ];
  return ptr;
}

template < typename Type >
void del2d( Type** ptr, int r ) {

  for ( int ii = 0; ii < r; ++ii )
    delete [] ptr[ ii ];
  delete [] ptr;

}

void timeclassic( int r, int c ) {

  double** uu = alloc2d< double >( r, c );
  double** laplacian = alloc2d< double >( r, c );

  timeme( uu, laplacian, r, c, "classic" );

  del2d( laplacian, r );
  del2d( uu, r );

}

void timeparen( int r, int c ) {

  StopWatch stopwatch( "sherpa::Array2d using operator(i,j)" );

  sherpa::Array2d< double > uu( r, c ), laplacian( r, c );

  for ( int ii = 0; ii < r; ++ii ) {
    for ( int jj = 0; jj < c; ++jj )
      uu( ii, jj ) = ii * ii + jj * jj;
  }

  for( int kk=0; kk < 10; kk++ ) {
    for( int ii = 1; ii < r - 1; ii++ ) {
      for( int jj = 1; jj < c - 1; jj++ ) {
        laplacian( ii, jj ) = - uu( ii - 1, jj ) - uu( ii, jj - 1 ) +
          4.0 * uu( ii, jj ) - uu( ii, jj + 1 ) - uu( ii + 1, jj );
      }
    }
  }

  double sum = 0.0;
  for ( int ii = 1; ii < r - 1; ++ii )
    for ( int jj = 1; jj < c - 1; jj++ )
      sum += laplacian( ii, jj );

  fprintf( stdout, "%.12f\t", sum );

}

void timebracket( int r, int c ) {

  sherpa::Array2d< double > uu( r, c ), laplacian( r, c );

  for ( int ii = 0; ii < r; ++ii ) {
    for ( int jj = 0; jj < c; ++jj )
      uu[ ii ][ jj ] = ii * ii + jj * jj;
  }

  timeme( uu, laplacian, r, c, "sherpa::Array2d[][]" );

}

void timemulti_array( int r, int c ) {

  std::vector<size_t> dims( 2, r );
  multi_array< double, 2 > laplacian( dims ), uu( dims );

  timeme( uu, laplacian, r, c, "Multi_Array[][]" );

}

void timeArray( int r, int c ) {

  std::vector<size_t> dims( 2, r );
  Array< double, 2 > laplacian( dims ), uu( dims );

  timeme( uu, laplacian, r, c, "Array[][]" );

}



void timeBavestrelliArray( int r, int c ) {

  bavestrelli::Array<double, 2> laplacian( bavestrelli::ArraySizes(r)(c) ),
    uu(  bavestrelli::ArraySizes(r)(c) );

  timeme( uu, laplacian, r, c, "bavestrelli::Array[][]" );

}

template< typename A, typename B >
void cmpme( A& a, B& b, int r, int c ) {

  for ( int ii = 0; ii < r; ++ii )
    for ( int jj = 0; jj < c; ++jj )
      if ( a[ ii ][ jj ] != b[ ii ][ jj ] )
	printf( "a[%d][%d] = %g vs b[%d][%d] = %g\n",
		ii, jj, a[ ii ][ jj ], ii, jj, b[ ii ][ jj ] );

}


/*
void testme( int r, int c ) {

  sherpa::Array2d< int > foo( r, c );

  for ( int ii = 0; ii < r; ++ii )
    for ( int jj = 0; jj < c; ++jj )
      foo[ ii ][ jj ] = ii + jj;

  std::cout << foo << '\n';

  sherpa::Array2d< int >::Row bar = foo[ r / 2 ];
  int* ptr = &bar[ 0 ];
  std::cout << ptr[ 0 ];
  for ( int jj = 1; jj < c; ++jj )
    std::cout << '\t' << ptr[ jj ];
  std::cout << '\n';

}
*/

int main( int argc, char* argv[] ) {

  const int num = 4000;

  int row = num;
  int col = num;
  if ( argc == 3 ) {
    row = atoi( argv[1] );
    col = atoi( argv[2] );
  }


  timeclassic( row, col );
  timeparen( row, col );
  timebracket( row, col );
  //timemulti_array( row, col );
  //timeArray( row, col );
  timeBavestrelliArray( row, col );

  return 0;

}

#endif

//
//cp Array2d.hh tmp.cc; g++ -Wall -ansi -pedantic -O3 -DtestArray2d -DNDEBUG tmp.cc; a.out
//
