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

#define _INTEGRATIONMODULE
#include "sherpa/integration.hh"
#include <iostream>
#include <Python.h>
#include "gsl_errno.h"
#include "gsl_integration.h"
#include "adapt_integrate.h"


namespace sherpa { namespace integration {

  static int gsl_int_flag = 1;
  
  int integrate_1d( integrand_1d fct, void* params,
		    double xlo, double xhi,
		    unsigned int maxeval, double epsabs, double epsrel,
		    double& result, double& abserr )
  {

    if ( NULL == fct )
      return EXIT_FAILURE;

    int retval;
    gsl_function func;
    func.function = fct;
    func.params = params;

    size_t neval = size_t( maxeval );

    gsl_set_error_handler_off ();
 
    retval = gsl_integration_qng( &func, xlo, xhi, epsabs, epsrel, &result,
				  &abserr, &neval );
    
    if (retval != EXIT_SUCCESS) {
      if( gsl_int_flag ) {
	std::cerr << "WARNING: Gauss-Kronrod integration has failed, "
		  << "using less precise means"
		  << std::endl;
	gsl_int_flag = 0;
      }      
      result = 0.5 * ( xhi - xlo ) * ( fct( xlo, params) + fct( xhi, params) );

      return EXIT_SUCCESS;
    }
    
    return EXIT_SUCCESS;

  }


  int integrate_Nd( integrand_Nd fct, void* params,
                    unsigned int ndim, const double* xlo, const double* xhi,
                    unsigned int maxeval, double epsabs, double epsrel,
                    double& result, double& abserr )
  {

    if ( NULL == fct || NULL == xlo || NULL == xhi )
      return EXIT_FAILURE;

    if ( 0 != adapt_integrate( fct, params, ndim, xlo, xhi, maxeval,
			       epsabs, epsrel, &result, &abserr ) )
      return EXIT_FAILURE;

    return EXIT_SUCCESS;

  }


}  }  /* namespace integration, namespace sherpa */


PyMODINIT_FUNC
initintegration(void)
{

  static void *Integration_API[2];

  PyObject *m;
  PyObject *api_cobject;

  if ( NULL == ( m = Py_InitModule( (char*)"integration", NULL ) ) )
    return;

  Integration_API[0] = (void*)sherpa::integration::integrate_1d;
  Integration_API[1] = (void*)sherpa::integration::integrate_Nd;

  if ( NULL == ( api_cobject = PyCObject_FromVoidPtr( (void*)Integration_API,
						      NULL) ) )
    return;

  // Since the actual data is static, we can let PyModule_AddObject()
  // steal the reference
  PyModule_AddObject( m, (char*)"_C_API", api_cobject );

}
