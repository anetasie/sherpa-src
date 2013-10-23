// 
//  Copyright (C) 2009  Smithsonian Astrophysical Observatory
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

#ifndef __sherpa_astro_xspec_extension_hh__
#define __sherpa_astro_xspec_extension_hh__

#include <sherpa/extension.hh>
#include <sherpa/constants.hh>
#include <vector>
#include <sstream>
#include <iostream>

namespace sherpa { namespace astro { namespace xspec {
 

  typedef sherpa::Array< float, NPY_FLOAT > FloatArray;
  typedef float FloatArrayType;


  template <npy_intp NumPars, bool HasNorm,
	    void (*XSpecFunc)( float* ear, int* ne, float* param, int* ifl, 
			       float* photar, float* photer )>
  PyObject* xspecmodelfct( PyObject* self, PyObject* args )
  {

#ifdef INIT_XSPEC
    if ( EXIT_SUCCESS != INIT_XSPEC() )
      return NULL;
#endif
    
    FloatArray pars;
    DoubleArray xlo;
    DoubleArray xhi;
    DoubleArray *x;

    if ( !PyArg_ParseTuple( args, (char*)"O&O&|O&",
			    (converter)convert_to_contig_array< FloatArray >,
			    &pars,
			    (converter)convert_to_contig_array< DoubleArray >,
			    &xlo,
			    (converter)convert_to_contig_array< DoubleArray >,
			    &xhi ) )
      return NULL;
    
    npy_intp npars = pars.get_size();
    
    if ( NumPars != npars ) {
      std::ostringstream err;
      err << "expected " << NumPars << " parameters, got " << npars;
      PyErr_SetString( PyExc_TypeError, err.str().c_str() );
      return NULL;
    }
    
    int nelem = int( xlo.get_size() );

    if ( nelem < 2 ) {
      std::ostringstream err;
      err << "input array must have at least 2 elements, found " << nelem;
      PyErr_SetString( PyExc_TypeError, err.str().c_str() );
      return NULL;
    }

    int ifl = 0;

    double hc = (sherpa::constants::c_ang<SherpaFloat>() *
		 sherpa::constants::h_kev<SherpaFloat>());
    bool is_wave = (xlo[0] > xlo[nelem-1]) ? true : false;

    // The XSPEC functions expect the input array to be of length ne+1
    int near = nelem;
    if( xhi )
      near++;

    std::vector<FloatArrayType> ear(near);

    for( int ii = 0; ii < nelem; ii++ ) {
      if( is_wave ) {

	// wave analysis swaps edges, e.g. wave_hi <--> energy_lo
	// if xhi is available use it 
	x = (xhi) ? &xhi : &xlo;

	if ( 0.0 == (*x)[ii] ) {
	  PyErr_SetString( PyExc_ValueError,
	             (char*)"XSPEC model evaluation failed, division by zero" );	  return NULL;
	}
	ear[ ii ] = ( FloatArrayType ) (hc / (*x)[ ii ]);
      }
      else
	ear[ ii ] = ( FloatArrayType ) xlo[ ii ];
    }

    if( xhi ) {
    
      if( is_wave ) {

	// wave analysis swaps edges, e.g. wave_lo <--> energy_hi
	// use xlo

	if ( 0.0 == xlo[ xlo.get_size() - 1 ] ) {
	  PyErr_SetString( PyExc_ValueError,
	             (char*)"XSPEC model evaluation failed, division by zero" );	  return NULL;
	}
	ear[ near - 1 ] = ( FloatArrayType ) (hc / xlo[ xlo.get_size() - 1 ]);
      }
      else
	ear[ near - 1 ] = ( FloatArrayType ) xhi[ xhi.get_size() - 1 ];
      
    }
    else
      nelem--;
    
    FloatArray result;
    if ( EXIT_SUCCESS != result.create( xlo.get_ndim(), xlo.get_dims() ) )
      return NULL;

    // Even though the XSPEC model function is Fortran, it could call
    // C++ functions, so swallow exceptions here

    try {

      XSpecFunc( &ear[0], &nelem, &pars[0], &ifl, &result[0], NULL );

    } catch(...) {

      PyErr_SetString( PyExc_ValueError,
		       (char*)"XSPEC model evaluation failed" );
      return NULL;

    }

    // Apply normalization if required
    if ( HasNorm )
      for ( int ii = 0; ii < nelem; ii++ )
	result[ii] *= pars[NumPars - 1];

    // The XSPEC functions expect the output array to be of length ne
    // (one less than the input array), so set the last element to
    // zero to avoid having random garbage in it
    if( !xhi )
      result[ result.get_size() - 1 ] = 0.0;

    return result.return_new_ref();

  }


  template <npy_intp NumPars, bool HasNorm,
	    void (*XSpecFunc)( const double* energy, int nFlux,
			       const double* params, int spectrumNumber,
			       double* flux, double* fluxError,
			       const char* initStr )>
  PyObject* xspecmodelfct_C( PyObject* self, PyObject* args )
  {

#ifdef INIT_XSPEC
    if ( EXIT_SUCCESS != INIT_XSPEC() )
      return NULL;
#endif
    
    DoubleArray pars;
    DoubleArray xlo;
    DoubleArray xhi;
    DoubleArray *x;

    if ( !PyArg_ParseTuple( args, (char*)"O&O&|O&",
			    (converter)convert_to_contig_array< DoubleArray >,
			    &pars,
			    (converter)convert_to_contig_array< DoubleArray >,
			    &xlo,
			    (converter)convert_to_contig_array< DoubleArray >,
			    &xhi ) )
      return NULL;
    
    npy_intp npars = pars.get_size();
    
    if ( NumPars != npars ) {
      std::ostringstream err;
      err << "expected " << NumPars << " parameters, got " << npars;
      PyErr_SetString( PyExc_TypeError, err.str().c_str() );
      return NULL;
    }
    
    int nelem = int( xlo.get_size() );

    if ( nelem < 2 ) {
      std::ostringstream err;
      err << "input array must have at least 2 elements, found " << nelem;
      PyErr_SetString( PyExc_TypeError, err.str().c_str() );
      return NULL;
    }

    double hc = (sherpa::constants::c_ang<SherpaFloat>() *
		 sherpa::constants::h_kev<SherpaFloat>());
    bool is_wave = (xlo[0] > xlo[nelem-1]) ? true : false;

    // The XSPEC functions expect the input array to be of length nFlux+1
    int near = nelem;
    if( xhi )
      near++;
    
    std::vector<SherpaFloat> ear(near);

    for( int ii = 0; ii < nelem; ii++ ) {
      if( is_wave ) {

	// wave analysis swaps edges, e.g. wave_hi <--> energy_lo
	// if xhi is available use it 
	x = (xhi) ? &xhi : &xlo;

	if ( 0.0 == (*x)[ii] ) {
	  PyErr_SetString( PyExc_ValueError,
		     (char*)"XSPEC model evaluation failed, division by zero" );	  return NULL;
	}
	ear[ ii ] = ( SherpaFloat ) (hc / (*x)[ ii ]);
      }
      else
	ear[ ii ] = ( SherpaFloat ) xlo[ ii ];
    }

    if( xhi ) {

      if( is_wave ) {

	// wave analysis swaps edges, e.g. wave_lo <--> energy_hi
	// use xlo 

	if ( 0.0 == xlo[ xlo.get_size() - 1 ] ) {
	  PyErr_SetString( PyExc_ValueError,
	             (char*)"XSPEC model evaluation failed, division by zero" );	  return NULL;
	}
	ear[ near - 1 ] = ( SherpaFloat ) (hc / xlo[ xlo.get_size() - 1 ]);
      }
      else
	ear[ near - 1 ] = ( SherpaFloat ) xhi[ xhi.get_size() - 1 ];
      
    }
    else
      nelem--;
    
    DoubleArray result;
    if ( EXIT_SUCCESS != result.create( xlo.get_ndim(), xlo.get_dims() ) )
      return NULL;

    // The XSPEC functions require fluxError to be non-NULL, so we create
    // it but discard it after the computation is done    
    DoubleArray error;
    if ( EXIT_SUCCESS != error.create( xlo.get_ndim(), xlo.get_dims() ) )
      return NULL;

    // Swallow C++ exceptions

    try {

      XSpecFunc( &ear[0], nelem, &pars[0], 0, &result[0], &error[0], NULL );

    } catch(...) {

      PyErr_SetString( PyExc_ValueError,
		       (char*)"XSPEC model evaluation failed" );
      return NULL;

    }

    // Apply normalization if required
    if ( HasNorm )
      for ( int ii = 0; ii < nelem; ii++ )
	result[ii] *= pars[NumPars - 1];

    // The XSPEC functions expect the output array to be of length nFlux
    // (one less than the input array), so set the last element to
    // zero to avoid having random garbage in it
    if( !xhi )
      result[ result.get_size() - 1 ] = 0.0;

    return result.return_new_ref();

  }


} } } /* namespace xspec, namespace astro, namespace sherpa */


#define _XSPECFCTSPEC(name, npars, has_norm) \
  FCTSPEC(name, (sherpa::astro::xspec::xspecmodelfct< npars, has_norm, \
                                                      name##_ >))

#define XSPECMODELFCT(name, npars)  _XSPECFCTSPEC(name, npars, false)
#define XSPECMODELFCT_NORM(name, npars)  _XSPECFCTSPEC(name, npars, true)

#define XSPECMODELFCT_C(name, npars) \
  FCTSPEC(name, (sherpa::astro::xspec::xspecmodelfct_C< npars, false, name >))

#define XSPECMODELFCT_C_NORM(name, npars) \
  FCTSPEC(name, (sherpa::astro::xspec::xspecmodelfct_C< npars, true, name >))


#endif /* __sherpa_astro_xspec_extension_hh__ */
