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

#include <cmath>
#include <cfloat>
#include <vector>
#include <sstream>
#include <iostream>
#include <sherpa/extension.hh>
#include "tcd/tcd.h"


extern "C" {
  void init_psf();
}

typedef struct {
  tcdDComplex *kernel_fft;
  long* newAxes;
} tcdCData;

typedef struct {
  // Note that there is no semicolon after the PyObject_HEAD macro;
  // one is included in the macro definition.
  PyObject_HEAD
  tcdCData* tcdCVal;
} tcdPyData;

// New datatype dealloc method
static void tcdPyData_dealloc(tcdPyData* data) {

  if( data->tcdCVal ) {

    if( data->tcdCVal->kernel_fft )
      tcdFreeTransformD( &data->tcdCVal->kernel_fft );

    if( data->tcdCVal->newAxes )
      free(data->tcdCVal->newAxes);
    
    free( data->tcdCVal );
  }    
  PyObject_Del(data);
}

// New datatype initialization
static PyTypeObject tcdPyData_Type = {
  PyObject_HEAD_INIT(NULL)
  0,                              // ob_size
  (char*)"tcdPyData",                    // tp_name
  sizeof(tcdPyData),              // tp_basicsize
  0, 		                  // tp_itemsize
  (destructor)tcdPyData_dealloc,  // tp_dealloc
};

// New datatype method to create new region
static tcdPyData *tcdPyData_new(tcdCData* tcdCVal)
{
  tcdPyData *mydata = PyObject_New(tcdPyData, &tcdPyData_Type);
  mydata->tcdCVal = tcdCVal;  
  return mydata;
}

PyObject* tcdKernel = NULL;

typedef sherpa::Array< long, NPY_LONG > LongArray;

static void _normalize_kernel( double* res, long len ) {

  long ii;
  double res_sum = 0.0;

  for( ii = 0; ii < len; ii++ )
    res_sum += res[ii];


  if( std::fabs(res_sum - 1) > DBL_EPSILON )
    for( ii = 0; ii < len; ii++ )
      res[ii] /= res_sum;
  
}


static int _extract_kernel( const double* psf, const long* dims,
			    const long ndims, const double* center,
			    long* newdims, DoubleArray& res, int rad,
			    double* xmin, double* xmax, double* width,
			    double& psffrac) {

  long size=1;
  std::vector<long> nsize(ndims);
  std::vector<long> xlo(ndims);
  std::vector<long> xhi(ndims);
  
  if( NULL == psf )
    return EXIT_FAILURE;
  
  for(long ii = 0 ; ii < ndims ; ii++ ) {

    //If the "subregion" of the PSF to be extracted is actually
    // bigger than the PSF, then reset the szs array to be of the 
    // size of the PSF.
    if( newdims[ii] > dims[ii] )
      newdims[ii] = dims[ii];
    
    //
    // To try to eliminate algorithmic problems of 2D "centroid" PSFs and
    // 1D radial profile "one-sided" PSFs, I've implemented wrap-around in
    // extraction.  (See coord, below.)
    //
    double diff = width[ii];
    double mid = center[ii];
    //double mid = (xmin[ii]+xmax[ii])/2.0;
    if( rad == 1) mid = 0.0; 

    // Need to revisit--is size defined in image units, or number of bins?
    // (Moot for logical coord images.)
    // Answer:  number of bins.  Hence we take bin widths into
    // account below.
    double xlotmp = mid - ( newdims[ii] * diff / 2.0 );
    double xhitmp = mid + ( newdims[ii] * diff / 2.0 );
    double dlo,dhi;
    // this helps to ensure proper centering of radial profile PSF
    if ( rad == 1 ) {
      dlo = ( xmin[ii] - fabs( xlotmp - xmin[ii] ) ) / diff;
      dhi = ( xmin[ii] + fabs( xhitmp - xmin[ii] ) ) / diff;
    } else {
      // FIXME: xmin appears to be in index coords??
      dlo = ( xlotmp - (xmin[ii]-1) ) / diff;
      dhi = ( xhitmp - (xmin[ii]-1) ) / diff;
    }

    xlo[ii] = (long) dlo;
    xhi[ii] = (long) dhi;
    //nsize[ii] = xhi[ii] - xlo[ii];
    //size *= nsize[ii];
    nsize[ii] = newdims[ii];
    size *= newdims[ii];
  }

  npy_intp npydims[1];
  npydims[0] = size;
  if ( EXIT_SUCCESS != res.create( 1, npydims ) )
    return EXIT_FAILURE;  
  
  // Extract y- and x-arrays from old PSF, centered on the
  // center position determined above.
  std::vector<long> coord(ndims);
  switch( ndims ) {
    
  case 1:
    
    for (long ii = 0 ; ii < nsize[0] ; ii++ ) {
      coord[0] = xlo[0] + ii;
      while ( coord[0] < 0 ) coord[0] += dims[0];
      while ( coord[0] >= dims[0] ) coord[0] -= dims[0];
      res[ ii ] = psf[ coord[0] ];
    }
    break;
      
  case 2:
    
    for (long jj = 0 ; jj < nsize[1] ; jj++ ) {
      coord[1] = xlo[1] + jj;
      while ( coord[1] < 0 ) coord[1] += dims[1];
      while ( coord[1] >= dims[1] ) coord[1] -= dims[1];
      for (long ii = 0 ; ii < nsize[0] ; ii++ ) {    
	coord[0] = xlo[0] + ii;
	while ( coord[0] < 0 ) coord[0] += dims[0];
        while ( coord[0] >= dims[0] ) coord[0] -= dims[0];
	
	long newidx = jj * nsize[0] + ii;
	// Find indices relevent for new bounds wrt old psf size (dims[0])
	long oldidx = coord[1] * dims[0] + coord[0];
	
	res[ newidx  ] = psf[ oldidx ];
      }
    }
    
    break;
    
  default:
    return EXIT_FAILURE;
    
  } // end switch

  psffrac = 0.0;
  for(long ii = 0; ii < size; ii++ )
    psffrac += res[ii];
  
  return EXIT_SUCCESS;  
}

static void _set_origin( long* dims_kern, long* kOrigin, long max, long nAxes )
{
  // if more than one pixel qualifies as brightest, such as const2D
  // use the middle of subkernel -- assumes the user provided center at
  // time of kernel extraction, so that should be middle of subkernel.
  if( max == -1 ) {
    for( int ii = 0; ii < nAxes; ii++ )
      if( dims_kern[ii] % 2 == 0 )
	kOrigin[ii] = dims_kern[ii]/2-1;
      else
	kOrigin[ii] = dims_kern[ii]/2;
    return;
  }
  
  if ( nAxes == 1 ) {
    // Just in case, go back to the center ymaxpos is outside
    // the PSF data space for some reason.
    kOrigin[0] = ( max < 0 || max > dims_kern[0]-1 ) ? dims_kern[0]/2 : max;
    return;
  }

  else if( nAxes == 2 ) {
    // Any bin number we get back is really the y-position * xsize +
    // x-position; therefore we must convert the bin number to x- and
    // y-positions.
    
    // Can do this because it drops remainder
    long ybin = max / dims_kern[0];
    long xbin = max - (ybin * dims_kern[0]);
    
    // Just in case, go back to the center if xbin, ybin turn out
    // to be outside the PSF data space for some reason.
    kOrigin[0] = ( xbin < 0 || xbin > dims_kern[0]-1 ) ? dims_kern[0]/2 : xbin;
    kOrigin[1] = ( ybin < 0 || ybin > dims_kern[1]-1 ) ? dims_kern[1]/2 : ybin;
    return;
  }
}

static int _pad_data(const int dim, double* res, double* src,
		     long* padSize, long* lAxes)
{
  // 10 Nov 1999 pef
  // 18 Jun 2007 br - updated PSF to CIAO4 Beta
  //
  // Pad the data, but not the kernel (TCD takes care of that).
  //

  if ( dim == 1 ) {
    for ( int i = 0 ; i < lAxes[0] ; i++ )
        res[i] = src[i];
    
    return EXIT_SUCCESS;
  }

  else if ( dim == 2 ) {
    for ( int i = 0 ; i < padSize[1] ; i++ ) {
      for ( int j = 0 ; j < padSize[0] ; j++ ) {
        int newIndex = i*padSize[0] + j;
        if ( i < lAxes[1] && j < lAxes[0] ) {
          int oldIndex = i*lAxes[0] + j;
          res[newIndex] = src[oldIndex];
        }
      }
    }
    
    return EXIT_SUCCESS;
  }

  return EXIT_FAILURE;
}


static int _unpad_data(const int dim, double* res, double* output,
		       long* padSize, long* lAxes)
{
  // 10 Nov 1999 pef
  //
  // Unpack the convolved model from the output array.
  //
  if ( dim == 1 ) {
    for ( int i = 0 ; i < lAxes[0] ; i++ ) {
      res[i] = output[i];
    }

    return EXIT_SUCCESS;
  }

  else if ( dim == 2 ) {
    for ( int i = 0 ; i < lAxes[1] ; i++ ) {
      for ( int j = 0 ; j < lAxes[0] ; j++ ) {
        int oldIndex = i*lAxes[0] + j;
        int newIndex = i*padSize[0] + j;
	res[oldIndex] = output[newIndex];
      }
    }

    return EXIT_SUCCESS;
  }

  return EXIT_FAILURE;
}

static int _convolve( double* source, double* kernel, long* dims_src,
		      long* dims_kern, long* dOrigin, long* kOrigin,
		      const long nAxes, double*& output ) {
  
  // TCDlib variables
  tcdDComplex *fftData = NULL, *fftKern = NULL, *ffttemp = NULL;
  tcdCData *tcdVals = NULL;
  long  *newAxes = NULL, *newTemp = NULL;
  
  if( tcdKernel && ( ( tcdPyData* ) tcdKernel )->tcdCVal ) {
    ffttemp = ( tcdDComplex* ) ( ( tcdPyData* ) (tcdKernel ) )->tcdCVal->kernel_fft;
    newTemp = ( long* )        ( ( tcdPyData* ) (tcdKernel ) )->tcdCVal->newAxes;
    if( ffttemp ) fftKern = ffttemp;
    if( newTemp ) newAxes = newTemp;
  }
  
  //if ( tcdSUCCESS != tcdInitConvolveOutputD( nAxes, dims_src, &output ) )
  //  return EXIT_FAILURE;
  
  if( fftKern && newAxes ) {
    if( tcdSUCCESS != tcdFFTConvolveD( tcdCONVOLVE, tcdDOUBLE, source,
  				       nAxes, dims_src, dOrigin, tcdDOUBLE,
				       NULL, dims_kern, kOrigin, &output,
				       &newAxes, &fftData, &fftKern ) )
      return EXIT_FAILURE;
  }
  else {
  if( tcdSUCCESS != tcdFFTConvolveD( tcdCONVOLVE, tcdDOUBLE, source,
				     nAxes, dims_src, dOrigin, tcdDOUBLE,
				     kernel, dims_kern, kOrigin, &output,
				     &newAxes, &fftData, &fftKern ) )
    return EXIT_FAILURE;
  }
  
  if ( tcdSUCCESS != tcdFreeTransformD( &fftData ) )
    return EXIT_FAILURE;

  //   fftKern is freed in tcdKernel dealloc
  
  if( !tcdKernel ) {
    tcdVals = ( tcdCData* ) malloc( sizeof( tcdCData ) );
    tcdVals->newAxes = newAxes;
    tcdVals->kernel_fft = fftKern;
    tcdKernel = ( PyObject* ) tcdPyData_new( tcdVals );
  }
  
  return EXIT_SUCCESS;
}

static PyObject* extract_kernel( PyObject* self, PyObject* args )
{
  
  DoubleArray kernel;
  DoubleArray res;
  DoubleArray xlo;
  DoubleArray xhi;
  DoubleArray widths;
  LongArray dims_kern;
  LongArray dims_new;
  DoubleArray center;
  int radial;
  
  
  if ( !PyArg_ParseTuple( args, (char*)"O&O&O&O&O&O&O&i",
			  CONVERTME( DoubleArray ),
			  &kernel,			 
			  CONVERTME( LongArray ),
			  &dims_kern,
			  CONVERTME( LongArray ),
			  &dims_new,
			  CONVERTME( DoubleArray ),
			  &center,
			  CONVERTME( DoubleArray ),
			  &xlo,
			  CONVERTME( DoubleArray ),
			  &xhi,
			  CONVERTME( DoubleArray ),
			  &widths,
			  &radial) )
    return NULL;
  
  double frac = 0.0;
  const long nAxes = (long) dims_kern.get_size();

  
  if( dims_new.get_size() != nAxes ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "dims_new: " << dims_new.get_size()
	<< " vs dims_kern: " << nAxes;
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  if( nAxes != center.get_size() ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "dims_kern: " << nAxes
	<< " vs center: " << center.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  } 

  if( nAxes != xlo.get_size() ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "dims_kern: " << nAxes
	<< " vs xlo: " << xlo.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  if( nAxes != xhi.get_size() ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "dims_kern: " << nAxes
	<< " vs xhi: " << xhi.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  if( nAxes != widths.get_size() ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "dims_kern: " << nAxes
	<< " vs widths: " << widths.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  if( EXIT_SUCCESS != _extract_kernel( &kernel[0], &dims_kern[0], nAxes,
				       &center[0], &dims_new[0], res,
				       radial, &xlo[0], &xhi[0], &widths[0],
				       frac) ) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"kernel extraction failed" );
    return NULL;
  }
  
  return Py_BuildValue( (char*)"(NNd)",
			res.return_new_ref(),
			dims_new.return_new_ref(),
			frac);
}

static PyObject* convolve( PyObject* self, PyObject* args )
{
  
  DoubleArray source;
  DoubleArray kernel;
  LongArray dims_src;
  LongArray dims_kern;
  LongArray center;
  
  if ( !PyArg_ParseTuple( args, (char*)"O&O&O&O&O&",
			  CONVERTME( DoubleArray ),
			  &source,
			  CONVERTME( DoubleArray ),
			  &kernel,
			  CONVERTME( LongArray ),
			  &dims_src,
			  CONVERTME( LongArray ),
			  &dims_kern,
			  CONVERTME( LongArray ),
			  &center) )
    return NULL;
  
  if( dims_src.get_size() != dims_kern.get_size() ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "dims_src: " << dims_src.get_size()
	<< " vs dims_kern: " << dims_kern.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  if( dims_kern.get_size() != center.get_size() ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "dims_kern: " << dims_kern.get_size()
	<< " vs center: " << center.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  } 
  
  long kern_size = 1, src_size = 1;
  // src and kernel dims should be equal length
  for( npy_intp ii = 0; ii < dims_kern.get_size(); ii++ ) {
    kern_size *= dims_kern[ii];
    src_size *= dims_src[ii];
  }
  
  if( source.get_size() != src_size ) {
    std::ostringstream err;
    err << "input array size do not match dimensions, "
	<< "source size: " << source.get_size()
	<< " vs source dim: " << src_size;
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }
  
  if( kernel.get_size() != kern_size ) {
    std::ostringstream err;
    err << "input array size do not match dimensions, "
	<< "kernel size: " << kernel.get_size()
	<< " vs kernel dim: " << kern_size;
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }
  
  const long nAxes = (long) dims_kern.get_size();
  double* output = NULL;
  std::vector<long> dOrigin(nAxes, 0);

  if( EXIT_SUCCESS != _convolve( &source[0], &kernel[0], &dims_src[0],
				 &dims_kern[0], &dOrigin[0], &center[0],
				 nAxes, output ) ) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"tcd convolution failed" );
    return NULL;
  } 

  DoubleArray result;
  npy_intp cdims[1];

  cdims[0] = 1;
  for( long ii = 0; ii < nAxes; ii++ )
    cdims[0] *= (npy_intp) ( ( tcdPyData* ) (tcdKernel ) )->tcdCVal->newAxes[ii];
  // newAxes is freed in tcdKernel dealloc

  // FIXME: why doesn't python free this memory when it goes out of scope???
  //if ( EXIT_SUCCESS != result.create( 1, cdims, output ) ) {
  if ( EXIT_SUCCESS != result.create( 1, cdims ) ) {
    return NULL;
  }

  for( long ii = 0; ii < cdims[0]; ii++ )
    result[ ii ] = output[ ii ];

  if(output) free(output);  

  return result.return_new_ref();
}


static PyObject* normalize( PyObject* self, PyObject* args )
{
  
  DoubleArray data;
  
  if ( !PyArg_ParseTuple( args, (char*)"O&", CONVERTME(DoubleArray), &data) )
    return NULL;
  
  _normalize_kernel( &data[0], (long) data.get_size() );
  
  return data.return_new_ref();
}


static PyObject* set_origin( PyObject* self, PyObject* args )
{
  
  LongArray dims;
  LongArray res;
  long max=-1;
  
  if ( !PyArg_ParseTuple( args, (char*)"O&|l", CONVERTME(LongArray), &dims,
			  &max) )
    return NULL;
  
  if( EXIT_SUCCESS != res.zeros( dims.get_ndim(), dims.get_dims() ) )
    return NULL;      
  
  _set_origin( &dims[0], &res[0], max, (long) dims.get_size() );
  
  return res.return_new_ref();
}

static int _get_bin_number(double d, DoubleArray& x,
			   DoubleArray& xhi, int size) {
  int i;
  if ( x && xhi ) {
    if ( d < x[0] ) return 0;
    for ( i = 0 ; i < size ; i++ ) {
      if ( d >= x[i] && d < xhi[i] ) return i;
    }
    return size-1;
  } else if ( x ) {
    // Assumes sorted values.
    double metric = DBL_MAX;
    for ( i = 0 ; i < size ; i++ ) {
      if ( fabs(d-x[i]) < metric ) {
        metric = fabs(d-x[i]);
      } else {
        if ( i > 0 ) return i-1;
        return i;
      }
    }
    return size-1;
  }
  return 0;
}

static PyObject* get_bin_number( PyObject* self, PyObject* args )
{
  
  DoubleArray xlo;
  DoubleArray xhi;
  double val;
  int res;
  
  if ( !PyArg_ParseTuple( args, (char*)"dO&|O&",
			  &val,
			  CONVERTME(DoubleArray), &xlo,
			  CONVERTME(DoubleArray), &xhi) )
    return NULL;
  
  if( xlo && xhi ) {
    if( xlo.get_size() != xhi.get_size() ) {
      std::ostringstream err;
      err << "input array sizes do not match, "
	  << "xlo: " << xlo.get_size() << " vs xhi: " << xhi.get_size();
      PyErr_SetString( PyExc_TypeError, err.str().c_str() );
      return NULL;
    }
  }
  res = _get_bin_number(val, xlo, xhi, (int)xlo.get_size());
  
  return Py_BuildValue( (char*)"i", res);
}

static long _pad(const long length)
{
  // 05 Nov 1999 pef
  //
  // Determines the proper padding factor given length, which should be
  // lAxes + ksize/2.  Proper means prime factorization in terms of
  // powers of two, three, and five.
  //
  // Fixed 31 Jan 2002 (some factors were missing. -pef

  long padding[238] =
      {   2,   3,   4,   5,   6,   8,   9,  10,  12,  15,  16,  18,  20,
         24,  25,  27,  30,  32,  36,  40,  45,  48,  50,  54,  60,  64,
         72,  75,  80,  81,  90,  96, 100, 108, 120, 125, 128, 135, 144, 
        150, 160, 162, 180, 192, 200, 216, 225, 240, 243, 250, 256, 270, 
        288, 300, 320, 324, 360, 375, 384, 400, 405, 432, 450, 480, 486, 
        500, 512, 540, 576, 600, 625, 640, 648, 675, 720, 729, 750, 768, 
        800, 810, 864, 900, 960, 972,1000,1024,1080,1125,1152,1200,1215,
       1250,1280,1296,1350,1440,1458,1500,1536,1600,1620,1728,1800,1875,
       1920,1944,2000,2025,2048,2160,2187,2250,2304,2400,2430,2500,2560,
       2592,2700,2880,2916,3000,3072,3125,3200,3240,3375,3456,3600,3645,
       3750,3840,3888,4000,4050,4096,4320,4374,4500,4608,4800,4860,5000,
       5120,5184,5400,5625,5760,5832,6000,6075,6144,6250,6400,6480,6561,
       6750,6912,7200,7290,7500,7680,7776,8000,8100,8192,8640,8748,9000,
       9216,9375,9600,9720,10000,10125,10240,10368,10800,10935,11250,
       11520,11664,12000,12150,12288,12500,12800,12960,13122,13500,13824,
       14400,14580,15000,15360,15552,15625,16000,16200,16384,16875,17280,
       17496,18000,18225,18432,18750,19200,19440,19683,20000,20250,20480,
       20736,21600,21870,22500,23040,23328,24000,24300,24576,25000,25600,
       25920,26244,27000,27648,28125,28800,29160,30000,30375,30720,31104,
       31250,32000,32400};

  for ( int i = 0; i < 238 ; i++ ) {
    if ( padding[i] >= length ) {
      return padding[i];
    }
  }

  return padding[237];
}


static PyObject* get_padsize( PyObject* self, PyObject* args )
{
  
  long size;  
  if ( !PyArg_ParseTuple( args, (char*)"l", &size) )
    return NULL;
 
  return Py_BuildValue( (char*)"l", _pad( size ) );
}

static PyObject* pad_data( PyObject* self, PyObject* args )
{
  
  DoubleArray kernel;
  LongArray shape;
  LongArray padshape;
  DoubleArray res;
  
  if ( !PyArg_ParseTuple( args, (char*)"O&O&O&",
			  CONVERTME(DoubleArray), &kernel,
			  CONVERTME(LongArray), &shape,
			  CONVERTME(LongArray), &padshape) )
    return NULL;

  if ( shape.get_size() != padshape.get_size() ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "shape: " << shape.get_size()
	<< " vs padshape: " << padshape.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }
  
  long size=1, padsize=1;
  for( npy_intp ii = 0; ii < shape.get_size(); ii++ ) {
    size *= shape[ii];

    if( padshape[ii] < shape[ii] ) {
      std::ostringstream err;
      err << "pad size is smaller than data shape, "
	  << "padshape[" << ii << "]: " << padshape[ii]
	  << " < shape[" << ii << "]: " << shape[ii];
      PyErr_SetString( PyExc_TypeError, err.str().c_str() );
      return NULL;
    }
    
    padsize *= padshape[ii];
  }

  if ( kernel.get_size() != size ) {
    std::ostringstream err;
    err << "input array size do not match dimensions, "
	<< "kernel size: " << kernel.get_size()
	<< " vs kernel dim: " << size;
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  } 
  
  npy_intp dims[1];
  dims[0] = padsize;

  // Create a zeros array to padding size
  if( EXIT_SUCCESS != res.zeros( kernel.get_ndim(), dims ) )
    return NULL;

  if( EXIT_SUCCESS != _pad_data( (int)shape.get_size(), &res[0], &kernel[0],
				 &padshape[0], &shape[0] ) ) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"padding kernel failed - dimension unsupported" );
    return NULL;
  } 
  
  return res.return_new_ref();
}

static PyObject* unpad_data( PyObject* self, PyObject* args )
{
  DoubleArray kernel;
  LongArray shape;
  LongArray padshape;
  DoubleArray res;
  
  if ( !PyArg_ParseTuple( args, (char*)"O&O&O&",
			  CONVERTME(DoubleArray), &kernel,
			  CONVERTME(LongArray), &padshape,
			  CONVERTME(LongArray), &shape) )
    return NULL;

  if ( shape.get_size() != padshape.get_size() ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "shape: " << shape.get_size()
	<< " vs padshape: " << padshape.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }
  
  long size=1, padsize=1;
  for( npy_intp ii = 0; ii < shape.get_size(); ii++ ) {
    size *= shape[ii];

    if( padshape[ii] < shape[ii] ) {
      std::ostringstream err;
      err << "pad size is smaller than data shape, "
	  << "padshape[" << ii << "]: " << padshape[ii]
	  << " < shape[" << ii << "]: " << shape[ii];
      PyErr_SetString( PyExc_TypeError, err.str().c_str() );
      return NULL;
    }
    
    padsize *= padshape[ii];
  }

  if ( kernel.get_size() != padsize ) {
    std::ostringstream err;
    err << "input array size do not match dimensions, "
	<< "kernel size: " << kernel.get_size()
	<< " vs kernel dim: " << size;
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  } 
  
  npy_intp dims[1];
  dims[0] = size;

  // Create a zeros array to padding size
  if( EXIT_SUCCESS != res.create( kernel.get_ndim(), dims ) )
    return NULL;

  if( EXIT_SUCCESS != _unpad_data( (int)shape.get_size(), &res[0], &kernel[0],
				 &padshape[0], &shape[0] ) ) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"unpadding kernel failed-dimension unsupported" );
    return NULL;
  } 
  
  return res.return_new_ref();
}

static PyObject* clear_kernel_fft( PyObject* self, PyObject* args )
{

  if( tcdKernel ) {
    tcdPyData_dealloc( ( tcdPyData* ) tcdKernel );
    tcdKernel = NULL;
  }

  // Py_None MUST be incremented before being returned!!
  Py_INCREF(Py_None);
  return Py_None;
}




static PyObject* pad_bounding_box( PyObject* self, PyObject* args )
{
  
  DoubleArray kernel;
  DoubleArray res;
  IntArray mask;
   
  if ( !PyArg_ParseTuple( args, (char*)"O&O&",
			  CONVERTME(DoubleArray), &kernel,
			  CONVERTME(IntArray), &mask) )
    return NULL;

  if ( kernel.get_size() > mask.get_size() ) {
    std::ostringstream err;
    err << "kernel size: " << kernel.get_size()
	<< " is > than mask size: " << mask.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  // Create a zeros array to padding size
  if( EXIT_SUCCESS != res.zeros( mask.get_ndim(), mask.get_dims() ) )
    return NULL;
  
  int size = mask.get_size();
  int jj = 0;
  for( int ii = 0; ii < size; ++ii )
    if( mask[ ii ] )
      res[ ii ] = kernel[ jj++ ];
  
  return res.return_new_ref();
}

static PyMethodDef PsfFcts[] = {

  FCTSPEC(extract_kernel, extract_kernel),

  // PSF convolution
  FCTSPEC(convolve, convolve),

  FCTSPEC(normalize, normalize),

  FCTSPEC(set_origin, set_origin),

  FCTSPEC(get_bin_number, get_bin_number),

  FCTSPEC(get_padsize, get_padsize),

  FCTSPEC(pad_data, pad_data),

  FCTSPEC(unpad_data, unpad_data),

  FCTSPEC(clear_kernel_fft, clear_kernel_fft),
  
  FCTSPEC( pad_bounding_box, pad_bounding_box ),
  
  { NULL, NULL, 0, NULL }

};

//SHERPAMOD(_psf, PsfFcts)

// Initialize the module
PyMODINIT_FUNC init_psf(void) {

  PyObject* m;

  if( PyType_Ready(&tcdPyData_Type) < 0 )
    return;
  
  import_array();  // Must be present for NumPy.
  
  m = Py_InitModule3( (char*)"_psf", PsfFcts, NULL);

  if( m == NULL )
    return;

  Py_INCREF(&tcdPyData_Type);
  PyModule_AddObject(m, (char*)"tcdKernel", (PyObject*)&tcdPyData_Type);
}
