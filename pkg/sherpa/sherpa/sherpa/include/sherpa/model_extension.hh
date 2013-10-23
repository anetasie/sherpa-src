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

#ifndef __sherpa_model_extension_hh__
#define __sherpa_model_extension_hh__

#include <sherpa/extension.hh>
#include <sherpa/integration.hh>
#include <sstream>
#include <iostream>


namespace sherpa { namespace models {

  template <int (*PtFunc)( const DoubleArray& p, double x, double& val )>
  double integrand_model1d( double x, void* params )
  {

    DoubleArray& p = *( static_cast< DoubleArray* >( params ) );
    double val = 0.0;

    // FIXME: do something with function return value
    PtFunc( p, x, val );

    return val;

  }


  template <int (*PtFunc)( const DoubleArray& p, double x, double& val )>
  int integrated_model1d( const DoubleArray& p, double xlo, double xhi,
			  double &val )
  {

    // FIXME: make these user-settable function args!
    double epsabs = 1.0e-14;
    double epsrel = 0.0;
    unsigned int maxeval = 10000;

    double abserr = 0.0;

    return integrate_1d( (integrand_1d)(integrand_model1d< PtFunc >),
			 (void*)&p, xlo, xhi, maxeval,
			 epsabs, epsrel, val, abserr );

  }


  template <int (*PtFunc)( const DoubleArray& p, double x0, double x1,
			   double& val )>
  double integrand_model2d( unsigned int ndim, const double* x, void* params )
  {

    DoubleArray& p = *( static_cast< DoubleArray* >( params ) );
    double val = 0.0;

    // FIXME: throw error if ndim != 2

    // FIXME: do something with function return value
    PtFunc( p, x[0], x[1], val );

    return val;

  }


  template <int (*PtFunc)( const DoubleArray& p, double x0, double x1,
			   double& val )>
  int integrated_model2d( const DoubleArray& p, double x0lo, double x0hi,
			  double x1lo, double x1hi, double &val )
  {

    // FIXME: make these user-settable function args!
    double epsabs = 1.0e-14;
    double epsrel = 0.0;
    unsigned int maxeval = 100000;

    double xlo[2];
    double xhi[2];
    xlo[0] = x0lo;
    xlo[1] = x1lo;
    xhi[0] = x0hi;
    xhi[1] = x1hi;

    double abserr = 0.0;

    return integrate_Nd( (integrand_Nd)(integrand_model2d< PtFunc >),
			 (void*)&p, 2, xlo, xhi, maxeval,
			 epsabs, epsrel, val, abserr );

  }


  template <typename ArrayType,
	    typename DataType,
	    npy_intp NumPars,
	    int (*PtFunc)( const ArrayType& p, DataType x, DataType& val ),
	    int (*IntFunc)( const ArrayType& p, DataType xlo, DataType xhi,
			    DataType& val )>
  PyObject* modelfct1d( PyObject* self, PyObject* args, PyObject *kwds)
  {

    ArrayType pars;
    ArrayType xlo;
    ArrayType xhi;
    int integrate = 1;

    static char *kwlist[] = {(char*)"pars",(char*)"xlo",(char*)"xhi",(char*)"integrate", NULL};

    if ( !PyArg_ParseTupleAndKeywords(args, kwds, (char*)"O&O&|O&i", kwlist,
			   (converter)convert_to_array< ArrayType >, &pars,
			   (converter)convert_to_array< ArrayType >, &xlo,
			   (converter)convert_to_array< ArrayType >, &xhi,
			   &integrate) )
      return NULL;
    
    npy_intp npars = pars.get_size();

    if ( NumPars != npars ) {
      std::ostringstream err;
      err << "expected " << NumPars << " parameters, got " << npars;
      PyErr_SetString( PyExc_TypeError, err.str().c_str() );
      return NULL;
    }

    npy_intp nelem = xlo.get_size();

    if ( xhi && ( xhi.get_size() != nelem ) ) {
      std::ostringstream err;
      err << "1D model evaluation input array sizes do not match, "
	  << "xlo: " << nelem << " vs xhi: " << xhi.get_size();
      PyErr_SetString( PyExc_TypeError, err.str().c_str() );
      return NULL;
    }

    ArrayType result;
    if ( EXIT_SUCCESS != result.create( xlo.get_ndim(), xlo.get_dims() ) )
      return NULL;


    if ( !(xhi && integrate) ) {

      for ( npy_intp ii = 0; ii < nelem; ii++ )
	if ( EXIT_SUCCESS != PtFunc( pars, xlo[ii], result[ii] ) ) {
	  PyErr_SetString( PyExc_ValueError,
			   (char*)"model evaluation failed" );
	  return NULL;
	}

    } else {

      for ( npy_intp ii = 0; ii < nelem; ii++ )
	if ( EXIT_SUCCESS != IntFunc( pars, xlo[ii], xhi[ii], result[ii] ) ) {
	  PyErr_SetString( PyExc_ValueError,
			   (char*)"model evaluation failed" );
	  return NULL;
	}

    }


    return result.return_new_ref();

  }


  template <typename ArrayType,
	    typename DataType,
	    npy_intp NumPars,
	    int (*PtFunc)( const ArrayType& p, DataType x0, DataType x1,
			   DataType& val ),
	    int (*IntFunc)( const ArrayType& p, DataType x0lo, DataType x0hi,
			    DataType x1lo, DataType x1hi, DataType& val )>
  PyObject* modelfct2d( PyObject* self, PyObject* args, PyObject *kwds )
  {

    ArrayType pars;
    ArrayType x0lo;
    ArrayType x1lo;
    ArrayType x0hi;
    ArrayType x1hi;

    int integrate = 1;
    static char *kwlist[] = {(char*)"pars", (char*)"x0lo", (char*)"x1lo",
			     (char*)"x0hi", (char*)"x1hi", (char*)"integrate", NULL};
    if ( !PyArg_ParseTupleAndKeywords( args, kwds, (char*)"O&O&O&|O&O&i", kwlist,
			    (converter)convert_to_array< ArrayType >, &pars,
			    (converter)convert_to_array< ArrayType >, &x0lo,
			    (converter)convert_to_array< ArrayType >, &x1lo,
			    (converter)convert_to_array< ArrayType >, &x0hi,
			    (converter)convert_to_array< ArrayType >, &x1hi,
			    &integrate) )
      return NULL;

    npy_intp npars = pars.get_size();

    if ( NumPars != npars ) {
      std::ostringstream err;
      err << "expected " << NumPars << " parameters, got " << npars;
      PyErr_SetString( PyExc_TypeError, err.str().c_str() );      
      return NULL;
    }

    if ( x0hi && !x1hi )  {
      PyErr_SetString( PyExc_TypeError, (char*)"expected 3 or 5 arguments, got 4");
      return NULL;
    }

    npy_intp nelem = x0lo.get_size();

    if ( ( x1lo.get_size() != nelem ) || 
	 ( x0hi &&
	   ( ( x0hi.get_size() != nelem ) ||
	     ( x1hi.get_size() != nelem ) ) ) ) {
      PyErr_SetString( PyExc_TypeError,
		       (char*)"2D model evaluation input array sizes do not match" );
      return NULL;
    }

    ArrayType result;
    if ( EXIT_SUCCESS != result.create( x0lo.get_ndim(), x0lo.get_dims() ) )
      return NULL;

    if ( !(x0hi && integrate) ) {

      for ( npy_intp ii = 0; ii < nelem; ii++ )
	if ( EXIT_SUCCESS != PtFunc( pars, x0lo[ii], x1lo[ii], result[ii] ) ) {
	  PyErr_SetString( PyExc_ValueError,
			   (char*)"model evaluation failed" );
	  return NULL;
	}

    } else {

      for ( npy_intp ii = 0; ii < nelem; ii++ )
	if ( EXIT_SUCCESS != IntFunc( pars, x0lo[ii], x0hi[ii], x1lo[ii],
				      x1hi[ii], result[ii] ) ) {
	  PyErr_SetString( PyExc_ValueError,
			   (char*)"model evaluation failed" );
	  return NULL;
	}

    }

    return result.return_new_ref();

  }


}  }  /* namespace models, namespace sherpa */


#define SHERPAMODELMOD(name, fctlist) \
PyMODINIT_FUNC \
init##name(void) \
{ \
  import_array(); \
  if ( -1 == import_integration() ) \
    return; \
  Py_InitModule( (char*)#name, fctlist );	\
}


// Allow this to be customized on a per-file basis

#define MODSPEC(name, func) \
  { (char*)#name, (PyCFunction)((PyCFunctionWithKeywords)func), METH_VARARGS|METH_KEYWORDS, NULL }

#ifndef _MODELFCTPTR
#define _MODELFCTPTR(name) \
  sherpa::models::name< SherpaFloat, SherpaFloatArray >
#endif

#define _MODELFCTSPEC(name, ftype, npars) \
  MODSPEC(name, (sherpa::models::ftype< SherpaFloatArray, SherpaFloat, npars, \
                                        _MODELFCTPTR(name##_point), \
                                        _MODELFCTPTR(name##_integrated) >))

#define _MODELFCTSPEC_NOINT(name, ftype, intftype, npars) \
  MODSPEC(name, \
          (sherpa::models::ftype< SherpaFloatArray, SherpaFloat, npars, \
                                  _MODELFCTPTR(name##_point), \
                                  sherpa::models::intftype \
                                    < _MODELFCTPTR(name##_point) > >))

#define MODELFCT1D(name, npars)		_MODELFCTSPEC(name, modelfct1d, npars)
#define MODELFCT2D(name, npars)		_MODELFCTSPEC(name, modelfct2d, npars)
#define MODELFCT1D_NOINT(name, npars) \
  _MODELFCTSPEC_NOINT(name, modelfct1d, integrated_model1d, npars)
#define MODELFCT2D_NOINT(name, npars) \
  _MODELFCTSPEC_NOINT(name, modelfct2d, integrated_model2d, npars)

#endif /* __sherpa_model_extension_hh__ */
