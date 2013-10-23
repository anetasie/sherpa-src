#ifdef testNelderMead

// 
//  Copyright (C) 2007  Smithsonian Astrophysical Observatory
//
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
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


#include "NelderMead.hh"

#include "tests/tstopt.hh"

void tstnm( Init init, Fct fct, int npar, std::vector<double>& par,
	    std::vector<double>& lo, std::vector<double>& hi,
	    double tol, const char* fct_name, int npop, int maxfev,
	    double c1, double c2 ) {

  try {

    char header[64];

    //
    // you may think you are clever by eliminating the following overhead
    // and simply use the vector par, but believe me it this is necessary
    //
    std::vector<double> mypar( npar, 0.0 );

    std::vector<double> step( npar * npar * 4 );
    for ( int ii = 0; ii < npar; ++ii )
      step[ ii ] = 1.2;

    std::vector< int > finalsimplex;

    for ( int ii = 0; ii < 2; ++ii ) {

      finalsimplex.push_back( ii );

      for ( int initsimplex = 0; initsimplex < 2; ++initsimplex ) {
	
	int mfcts;
	double answer;
	
	init( npar, mfcts, answer, &par[0], &lo[0], &hi[0] );
	for ( int jj = 0; jj < npar; ++jj )
	  mypar[ jj ] = par[ jj ];
	sherpa::NelderMead< Fct, void* > nm( fct, NULL );
	
	int verbose=0, maxnfev=npar*npar*maxfev, nfev;
	double fmin;
	nm( verbose, maxnfev, tol, npar, initsimplex, finalsimplex, lo, hi,
	    step, mypar, nfev, fmin );
	
	sprintf( header, "NelderMead_%d_", initsimplex );
	print_pars( header, fct_name, nfev, fmin, answer, npar, mypar );

      }
	
    }
    
  } catch( const sherpa::OptErr& oe ) {
    
    std::cerr << oe << '\n';
    
  }

  return;

}  

int main( int argc, char* argv[] ) {

  try {

    int c, uncopt = 1, globalopt = 1;
    while ( --argc > 0 && (*++argv)[ 0 ] == '-' )
      while ( c = *++argv[ 0 ] )
	switch( c ) {
	case 'u':
	  uncopt = 0;
	  break;
	case 'g':
	  globalopt = 0;
	  break;
	default:
	  fprintf( stderr, "%s: illegal option '%c'\n", argv[ 0 ], c );
	  fprintf( stderr, "Usage %s [ -g ] [ -u ] [ npar ]\n", argv[ 0 ] );
	  return EXIT_FAILURE;
      }


    int npar=6;
    if ( argc == 1 )
      npar = atoi( *argv );
    
    if ( npar % 2 || npar < 2 ) {
      printf( "The minimum value for the free parameter must be an even "
	      "and it is greater then 2\n" );
      return EXIT_FAILURE;
    }

    double tol = 1.0e-8;
    std::cout << "#\n#:npar = " << npar << "\n";
    std::cout << "#:tol=" << tol << '\n';
    std::cout << "# A negative value for the nfev signifies that the "
      "optimization method did not converge\n#\n";
    std::cout << "name\tnfev\tanswer\tstat\tpar\nS\tN\tN\tN\tN\n";

    int npop=0, maxfev=1024;
    double c1=0.0, c2=0.0;
    if ( uncopt )
      tst_unc_opt( npar, tol, tstnm, npop, maxfev, c1, c2 );

    if ( globalopt )
      tst_global( npar, tol, tstnm, npop, maxfev, c1, c2 );

    return EXIT_SUCCESS;

  } catch( std::exception& e ) {

    std::cerr << e.what( ) << '\n';
    return EXIT_FAILURE;

  }

}

/*
==1929== Memcheck, a memory error detector
==1929== Copyright (C) 2002-2009, and GNU GPL'd, by Julian Seward et al.
==1929== Using Valgrind-3.5.0 and LibVEX; rerun with -h for copyright info
==1929== Command: tstnm
==1929==
#
#:npar = 6
#:tol=1e-08
# A negative value for the nfev signifies that the optimization method did not converge
#
name	nfev	answer	stat	par
S	N	N	N	N
NelderMead_0_Rosenbrock	1197	0	0.00435378	1.02286,1.04643,1.04977,1.10213,0.963276,0.927824
NelderMead_1_Rosenbrock	-1658	0	0.0506586	0.902054,0.813548,1.20214,1.44569,1.0131,1.02633
NelderMead_0_Rosenbrock	1654	0	5.84053e-09	1.00001,1.00001,0.99997,0.999943,1.00004,1.00009
NelderMead_1_Rosenbrock	2152	0	3.72013e-09	1.00005,1.00009,0.999968,0.999934,1.00001,1.00002
NelderMead_0_FreudensteinRoth	-530	0	146.953	11.4117,-0.896875,11.414,-0.89672,11.4139,-0.896739
NelderMead_1_FreudensteinRoth	-448	0	146.953	11.4135,-0.896757,11.4121,-0.896834,11.4133,-0.896775
NelderMead_0_FreudensteinRoth	-735	0	146.953	11.4129,-0.896808,11.4132,-0.896803,11.4128,-0.896788
NelderMead_1_FreudensteinRoth	-626	0	146.953	11.4128,-0.896819,11.4124,-0.89682,11.4124,-0.89681
NelderMead_0_PowellBadlyScaled	1429	0	8.03264e-08	8.38095e-06,11.9321,1.26161e-05,7.92642,9.29973e-06,10.7529
NelderMead_1_PowellBadlyScaled	1468	0	2.57977e-06	4.17305e-06,23.9559,1.48524e-05,6.73218,1.49442e-05,6.69038
NelderMead_0_PowellBadlyScaled	2551	0	7.73439e-08	8.83658e-06,11.3166,1.26077e-05,7.93166,9.42547e-06,10.6096
NelderMead_1_PowellBadlyScaled	1801	0	2.43226e-06	4.17384e-06,23.9585,1.48476e-05,6.73462,1.49367e-05,6.69459
NelderMead_0_BrownBadlyScaled	-716	0	2.97791e+12	1491.49,-99.9975,1351.81,6785.79,8263.52,-0.372906
NelderMead_1_BrownBadlyScaled	-383	0	4.90426e+11	1.09578e+06,479512,896591,86645,517175,5.58367e-05
NelderMead_0_BrownBadlyScaled	-971	0	2.96981e+12	1773.84,903.886,2757.35,-99.9968,10628.5,0.0279322
NelderMead_1_BrownBadlyScaled	-749	0	3.12746e+11	1.07107e+06,436191,1.14219e+06,250217,813975,2.58843e-05
NelderMead_0_Beale	372	0	1.89415e-08	3.00017,0.500055,3.00018,0.500035,3.00007,0.500026
NelderMead_1_Beale	454	0	2.08803e-08	2.99988,0.499969,3.00028,0.500084,3.00004,0.500014
NelderMead_0_Beale	587	0	1.2406e-08	2.99998,0.500006,3.00002,0.500014,2.99988,0.499986
NelderMead_1_Beale	676	0	1.01467e-08	3.00002,0.500012,2.99994,0.499967,3.00004,0.500007
NelderMead_0_JennrichSampson	475	373.086	373.087	0.257842,0.25782,0.257841,0.25782,0.257822,0.257825
NelderMead_1_JennrichSampson	453	373.086	373.087	0.257852,0.257803,0.257842,0.257805,0.257846,0.257806
NelderMead_0_JennrichSampson	756	373.086	373.087	0.257819,0.257829,0.257837,0.257817,0.257832,0.257821
NelderMead_1_JennrichSampson	701	373.086	373.087	0.257816,0.257833,0.25783,0.257825,0.257829,0.25783
NelderMead_0_HelicalValley	153	0	2.41172e-09	0.999998,-2.33949e-05,-3.92257e-05
NelderMead_1_HelicalValley	191	0	5.87823e-09	0.999996,-2.84244e-05,-4.01583e-05
NelderMead_0_HelicalValley	285	0	1.51587e-09	0.999997,1.48411e-05,2.40208e-05
NelderMead_1_HelicalValley	309	0	4.0028e-09	1.00001,-2.60728e-06,-1.22224e-06
NelderMead_0_Bard	167	0.00821487	0.00821488	0.0824145,1.1331,2.34363
NelderMead_1_Bard	192	0.00821487	0.00821488	0.0824133,1.13308,2.34366
NelderMead_0_Bard	288	0.00821487	0.00821488	0.0824095,1.13299,2.34375
NelderMead_1_Bard	313	0.00821487	0.00821488	0.0824133,1.13308,2.34366
NelderMead_0_Gaussian	111	1.12793e-08	1.14499e-08	0.39895,1,1.81238e-05
NelderMead_1_Gaussian	112	1.12793e-08	1.13965e-08	0.39895,0.999984,5.31944e-06
NelderMead_0_Gaussian	229	1.12793e-08	1.14499e-08	0.39895,1,1.81238e-05
NelderMead_1_Gaussian	221	1.12793e-08	1.13965e-08	0.39895,0.999984,5.31944e-06
NelderMead_0_Meyer	-30	87.9458	1.139e+07	0.0597053,4000.76,249.722
NelderMead_1_Meyer	-21	87.9458	1.80421e+08	0.071666,4000.34,250.232
NelderMead_0_Meyer	-1624	87.9458	3092.66	0.00983416,5723.15,329.538
NelderMead_1_Meyer	-270	87.9458	111643	0.101549,4002.82,263.992
NelderMead_0_GulfResearchDevelopment	488	0	6.03061e-11	58.0921,24.2541,1.53308
NelderMead_1_GulfResearchDevelopment	83	0	8.9791e-05	5.52721,3.20528,0.76279
NelderMead_0_GulfResearchDevelopment	1367	0	2.47267e-17	50.004,24.9996,1.50002
NelderMead_1_GulfResearchDevelopment	539	0	3.96308e-10	74.6596,23.0122,1.58805
NelderMead_0_Box3d	226	0	1.64492e-10	1.00002,9.99958,0.999978
NelderMead_1_Box3d	231	0	5.03599e-10	0.999952,10.0006,1.00004
NelderMead_0_Box3d	320	0	1.64492e-10	1.00002,9.99958,0.999978
NelderMead_1_Box3d	332	0	2.68907e-10	0.999999,10.0004,1.00001
NelderMead_0_PowellSingular	388	0	1.66094e-12	0.000922767,-9.22902e-05,0.000441954,0.000441754
NelderMead_1_PowellSingular	374	0	1.30643e-15	-6.31659e-05,6.31858e-06,1.79201e-05,1.79104e-05
NelderMead_0_PowellSingular	576	0	3.73417e-15	4.80257e-06,-4.8136e-07,5.67939e-05,5.68199e-05
NelderMead_1_PowellSingular	547	0	1.30643e-15	-6.31659e-05,6.31858e-06,1.79201e-05,1.79104e-05
NelderMead_0_Wood	711	0	3.06143e-08	0.999946,0.999899,1.00004,1.0001
NelderMead_1_Wood	813	0	7.06812e-09	0.999968,0.999937,1.00004,1.00008
NelderMead_0_Wood	875	0	4.01606e-09	1.00003,1.00005,0.999977,0.999953
NelderMead_1_Wood	976	0	5.36986e-09	1.00002,1.00005,0.999986,0.999972
NelderMead_0_KowalikOsborne	278	0.000307505	0.000307506	0.192809,0.191312,0.123089,0.136076
NelderMead_1_KowalikOsborne	310	0.000307505	0.000307506	0.192808,0.191298,0.123076,0.136069
NelderMead_0_KowalikOsborne	441	0.000307505	0.000307506	0.192806,0.191257,0.123032,0.136052
NelderMead_1_KowalikOsborne	484	0.000307505	0.000307506	0.192808,0.191298,0.123076,0.136069
NelderMead_0_BrownDennis	234	85822.2	85822.2	-11.5944,13.2035,-0.403895,0.236169
NelderMead_1_BrownDennis	224	85822.2	85822.2	-11.5949,13.2038,-0.403349,0.23675
NelderMead_0_BrownDennis	355	85822.2	85822.2	-11.5946,13.2038,-0.40314,0.23692
NelderMead_1_BrownDennis	340	85822.2	85822.2	-11.5949,13.2038,-0.403349,0.23675
NelderMead_0_Osborne1	460	5.46489e-05	5.50277e-05	0.376257,2.04397,-1.57345,0.0130751,0.0217191
NelderMead_1_Osborne1	-516	5.46489e-05	6.87987e-05	0.370176,1.54872,-1.07349,0.0118644,0.0243448
NelderMead_0_Osborne1	659	5.46489e-05	5.50277e-05	0.376257,2.04397,-1.57345,0.0130751,0.0217191
NelderMead_1_Osborne1	1567	5.46489e-05	5.46489e-05	0.37541,1.93582,-1.46466,0.0128675,0.0221228
NelderMead_0_Biggs	626	0	0	1,10,1,5,4,3
NelderMead_1_Biggs	249	0	0	1,10,1,5,4,3
NelderMead_0_Biggs	1252	0	0	1,10,1,5,4,3
NelderMead_1_Biggs	498	0	0	1,10,1,5,4,3
NelderMead_0_Osborne2	-2853	0.0401377	0.0443746	1.28814,0.40824,0.633314,0.579866,0.710138,1.00532,1.29459,5.15732,2.42452,4.58027,5.69192
NelderMead_1_Osborne2	-2426	0.0401377	0.0845621	1.29114,0.386615,0.62108,0.408116,0.669402,1.40689,0.698202,13.1647,2.3395,4.75763,5.70398
NelderMead_0_Osborne2	3870	0.0401377	0.0401683	1.30997,0.431442,0.633621,0.599272,0.753864,0.905773,1.36485,4.82508,2.3988,4.56888,5.67538
NelderMead_1_Osborne2	5172	0.0401377	0.0402209	1.30699,0.432518,0.634409,0.598327,0.753971,0.901382,1.36738,4.8245,2.39833,4.57019,5.67577
NelderMead_0_Watson	627	0.00228767	0.00228767	-0.0157235,1.01243,-0.232966,1.2604,-1.51373,0.993008
NelderMead_1_Watson	550	0.00228767	0.00228767	-0.0157257,1.01244,-0.233019,1.2605,-1.51379,0.993011
NelderMead_0_Watson	894	0.00228767	0.00228767	-0.0157231,1.01243,-0.232989,1.26046,-1.51378,0.993028
NelderMead_1_Watson	795	0.00228767	0.00228767	-0.0157257,1.01244,-0.233019,1.2605,-1.51379,0.993011
NelderMead_0_PenaltyI	-1150	9.37629e-06	2.24998e-05	0.250536,0.249133,0.249885,0.250474
NelderMead_1_PenaltyI	-1397	9.37629e-06	2.24998e-05	0.249972,0.250038,0.250028,0.249992
NelderMead_0_PenaltyI	-1364	9.37629e-06	2.24998e-05	0.250006,0.249975,0.249997,0.250053
NelderMead_1_PenaltyI	-1549	9.37629e-06	2.24998e-05	0.249991,0.249985,0.250005,0.250048
NelderMead_0_PenaltyII	1531	9.37629e-06	9.51886e-06	0.199972,0.235897,0.5592,0.218371
NelderMead_1_PenaltyII	149	9.37629e-06	9.45514e-06	0.200003,0.272851,0.346956,0.613102
NelderMead_0_PenaltyII	3185	9.37629e-06	9.37663e-06	0.199997,0.195042,0.484442,0.506477
NelderMead_1_PenaltyII	1239	9.37629e-06	9.39556e-06	0.199999,0.25432,0.460263,0.471467
NelderMead_0_VariablyDimensioned	303	0	3.475e-08	0.999994,0.999897,0.999961,1.00012,0.999966,0.999989
NelderMead_1_VariablyDimensioned	383	0	1.35619e-08	1.00003,0.999937,0.999987,0.999974,1.00002,1.00004
NelderMead_0_VariablyDimensioned	538	0	9.81666e-09	0.999977,0.999949,0.999976,1.00006,0.999977,1
NelderMead_1_VariablyDimensioned	605	0	1.02168e-08	0.999974,1.00007,0.999989,0.999965,1.00004,0.999974
NelderMead_0_Trigonometric	433	0	1.40747e-08	0.00705824,0.00696962,0.0067329,-0.117601,0.00653846,0.0064437
NelderMead_1_Trigonometric	386	0	1.57647e-09	0.00704857,0.00690729,0.00675235,-0.117636,0.00650421,0.00641273
NelderMead_0_Trigonometric	676	0	2.53731e-09	0.00703642,0.00688127,0.00676746,-0.117661,0.00649293,0.00640494
NelderMead_1_Trigonometric	626	0	1.45373e-09	0.00700695,0.00685928,0.00676118,-0.117613,0.00651907,0.00642281
NelderMead_0_BrownAlmostLinear	-337	1	5.10516e-09	1.00007,1.00006,1.00004,1.00005,1.00009,0.99965
NelderMead_1_BrownAlmostLinear	-357	1	3.4561e-09	0.999962,0.999984,1.00001,0.999977,1.00002,1.00004
NelderMead_0_BrownAlmostLinear	-574	1	7.08455e-10	0.999987,0.99998,0.999983,0.999985,0.999973,1.00011
NelderMead_1_BrownAlmostLinear	-571	1	2.67574e-09	1.00002,1.00003,1.00002,1.00003,0.999974,0.999917
NelderMead_0_DiscreteBoundary	289	0	1.49323e-09	-0.0898696,-0.167835,-0.1536,-0.13246,-0.102081,-0.0593097
NelderMead_1_DiscreteBoundary	262	0	1.66255e-09	-0.0899031,-0.16788,-0.15363,-0.132474,-0.102067,-0.0593257
NelderMead_0_DiscreteBoundary	542	0	5.51653e-10	-0.0898787,-0.167848,-0.153594,-0.132449,-0.102034,-0.0592801
NelderMead_1_DiscreteBoundary	530	0	3.40351e-10	-0.0898959,-0.167877,-0.153619,-0.132464,-0.102071,-0.0593007
NelderMead_0_DiscreteIntegral	266	0	1.76342e-09	-0.065457,-0.118143,-0.154624,-0.169972,-0.157286,-0.106014
NelderMead_1_DiscreteIntegral	237	0	1.4493e-09	-0.0654704,-0.118147,-0.154633,-0.17001,-0.157249,-0.106044
NelderMead_0_DiscreteIntegral	531	0	1.76342e-09	-0.065457,-0.118143,-0.154624,-0.169972,-0.157286,-0.106014
NelderMead_1_DiscreteIntegral	495	0	1.4493e-09	-0.0654704,-0.118147,-0.154633,-0.17001,-0.157249,-0.106044
NelderMead_0_BroydenTridiagonal	224	0	8.70649e-08	-0.576047,-0.695921,-0.680251,-0.642928,-0.556413,-0.366021
NelderMead_1_BroydenTridiagonal	217	0	5.43608e-08	-0.576064,-0.695959,-0.680252,-0.643012,-0.556464,-0.366023
NelderMead_0_BroydenTridiagonal	484	0	1.48415e-08	-0.576061,-0.695932,-0.680253,-0.64297,-0.55643,-0.366035
NelderMead_1_BroydenTridiagonal	459	0	4.40255e-08	-0.576072,-0.695946,-0.68028,-0.642972,-0.556422,-0.366046
NelderMead_0_BroydenBanded	242	0	1.41433e-07	-0.42828,-0.476565,-0.519642,-0.558101,-0.59347,-0.593434
NelderMead_1_BroydenBanded	269	0	1.24035e-07	-0.428273,-0.47657,-0.519687,-0.558068,-0.593434,-0.593409
NelderMead_0_BroydenBanded	518	0	8.57026e-09	-0.428298,-0.476594,-0.519662,-0.558081,-0.593434,-0.593431
NelderMead_1_BroydenBanded	532	0	2.11593e-08	-0.428289,-0.476593,-0.519663,-0.55807,-0.593443,-0.593423
NelderMead_0_LinearFullRank	249	0	1.18546e-08	-1.00003,-0.999933,-1,-1.00003,-0.999946,-0.999949
NelderMead_1_LinearFullRank	259	0	1.35489e-08	-0.99996,-1.00001,-1.00006,-1,-1.00008,-0.999957
NelderMead_0_LinearFullRank	494	0	1.18546e-08	-1.00003,-0.999933,-1,-1.00003,-0.999946,-0.999949
NelderMead_1_LinearFullRank	506	0	1.06097e-08	-0.999987,-1.00006,-0.999943,-0.999964,-0.999956,-1.00001
NelderMead_0_LinearFullRank1	199	1.15385	1.15385	2.97371,2.20459,2.13584,0.159209,0.599986,-2.86608
NelderMead_1_LinearFullRank1	213	1.15385	1.15385	2.15604,0.91343,1.2103,-0.975731,-0.292156,-0.336558
NelderMead_0_LinearFullRank1	481	1.15385	1.15385	3.07835,2.28814,2.20579,0.177665,0.583544,-2.94494
NelderMead_1_LinearFullRank1	522	1.15385	1.15385	2.24657,0.983871,1.26376,-0.970922,-0.313712,-0.387097
NelderMead_0_LinearFullRank0cols0rows	189	2.66667	2.66667	2.14804,1.95635,0.261106,0.864515,-1.56416,1.92177
NelderMead_1_LinearFullRank0cols0rows	196	2.66667	2.66667	1.56493,0.703205,0.818582,-0.603392,-0.223035,1.69523
NelderMead_0_LinearFullRank0cols0rows	458	2.66667	2.66667	2.54626,1.95813,0.262525,0.864232,-1.56549,2.27683
NelderMead_1_LinearFullRank0cols0rows	465	2.66667	2.66667	1.65756,0.765659,0.86782,-0.600641,-0.279777,1.78379
NelderMead_0_Chebyquad	898	0	2.42098e-05	0.0450991,0.215181,0.222438,0.429553,0.486155,0.593495,0.760509,0.805557,0.956184
NelderMead_1_Chebyquad	2226	0	7.24599e-07	0.0441004,0.23679,0.19815,0.414455,0.582454,0.50183,0.799581,0.765167,0.955723
NelderMead_0_Chebyquad	2046	0	1.50681e-09	0.0442028,0.199493,0.235614,0.416064,0.499987,0.583977,0.764389,0.800515,0.955792
NelderMead_1_Chebyquad	2617	0	4.6099e-07	0.0440289,0.236508,0.198236,0.414554,0.582641,0.501506,0.799873,0.764807,0.955755
NelderMead_0_McCormick	66	-1.91	-1.91322	-0.547127,-1.54715
NelderMead_1_McCormick	68	-1.91	-1.91322	-0.547201,-1.54716
NelderMead_0_McCormick	124	-1.91	-1.91322	-0.547214,-1.54716
NelderMead_1_McCormick	129	-1.91	-1.91322	-0.547201,-1.54716
NelderMead_0_BoxBetts	112	0	1.46974e-10	1.00004,9.9997,0.999971
NelderMead_1_BoxBetts	119	0	1.8843e-10	1.00003,9.99999,0.999983
NelderMead_0_BoxBetts	183	0	1.46974e-10	1.00004,9.9997,0.999971
NelderMead_1_BoxBetts	191	0	1.8843e-10	1.00003,9.99999,0.999983
NelderMead_0_Paviani	916	-45.7	-45.7785	9.35053,9.35058,9.35035,9.35037,9.34932,9.35131,9.34967,9.35021,9.35147,9.35032
NelderMead_1_Paviani	1216	-45.7	-45.7784	9.34889,9.35141,9.34997,9.34849,9.35104,9.35041,9.3501,9.35006,9.3512,9.35108
NelderMead_0_Paviani	1299	-45.7	-45.7785	9.35046,9.35051,9.3504,9.3503,9.34999,9.35023,9.35017,9.35022,9.35029,9.35034
NelderMead_1_Paviani	1679	-45.7	-45.7785	9.35019,9.35034,9.35017,9.35011,9.35015,9.35013,9.35021,9.35031,9.35031,9.35012
NelderMead_0_GoldsteinPrice	63	3	3	5.56372e-06,-1.00002
NelderMead_1_GoldsteinPrice	72	3	3	-2.83725e-05,-0.999983
NelderMead_0_GoldsteinPrice	134	3	3	-4.76258e-06,-0.999994
NelderMead_1_GoldsteinPrice	145	3	3	8.28096e-06,-0.999998
NelderMead_0_Shekel5	-115	-10.1532	-5.10076	7.99917,7.99921,7.99968,7.99955
NelderMead_1_Shekel5	-120	-10.1532	-5.10076	8,8,8,8
NelderMead_0_Shekel5	-275	-10.1532	-5.10077	7.99961,7.99965,7.99954,7.99961
NelderMead_1_Shekel5	-282	-10.1532	-5.10077	7.9996,7.99964,7.99961,7.99963
NelderMead_0_Shekel7	-117	-10.4029	-5.12881	7.99939,8,7.99911,7.99939
NelderMead_1_Shekel7	-121	-10.4029	-5.12881	7.99921,7.99973,7.9999,7.99899
NelderMead_0_Shekel7	-278	-10.4029	-5.12882	7.99951,7.9996,7.99952,7.9996
NelderMead_1_Shekel7	-288	-10.4029	-5.12882	7.99955,7.99964,7.99948,7.99962
NelderMead_0_Shekel10	-116	-10.5364	-5.17563	7.99911,8.00003,7.99955,7.99904
NelderMead_1_Shekel10	-124	-10.5364	-5.17563	7.99957,7.9988,7.99962,7.99955
NelderMead_0_Shekel10	-270	-10.5364	-5.17565	7.9995,7.99944,7.99945,7.99946
NelderMead_1_Shekel10	-291	-10.5364	-5.17565	7.99946,7.99949,7.99945,7.99942
NelderMead_0_Levy4	-528	-21.502	8.95668	2.64777,1.33177,1.30012,6.99795
NelderMead_1_Levy4	-756	-21.502	15.0168	0.00938728,-1.63593,1.99843,6.99607
NelderMead_0_Levy4	-667	-21.502	8.95668	2.64772,1.33182,1.30044,6.99796
NelderMead_1_Levy4	-897	-21.502	15.0141	0.0122165,-1.63428,1.99839,6.99846
NelderMead_0_Levy5	-100	-11.504	38.8617	3.96389,3.99626,3.99616,3.99627,3.9994
NelderMead_1_Levy5	-198	-11.504	32.0094	3.30626,4.32626,3.66394,3.3297,3.99914
NelderMead_0_Levy5	-233	-11.504	38.8617	3.9638,3.99613,3.99623,3.99625,3.99946
NelderMead_1_Levy5	-633	-11.504	2.54461	1.65918,1.32488,0.999992,1.00007,2.99363
NelderMead_0_Levy6	-99	-11.504	47.8506	3.96384,3.99586,3.99612,3.99607,3.99638,3.99953
NelderMead_1_Levy6	-228	-11.504	47.8505	3.96362,3.99621,3.9963,3.99611,3.99621,3.99951
NelderMead_0_Levy6	-267	-11.504	47.8504	3.96389,3.99615,3.99624,3.99622,3.99625,3.99946
NelderMead_1_Levy6	-478	-11.504	47.8504	3.96392,3.99615,3.99624,3.99622,3.99623,3.99945
NelderMead_0_Levy7	-115	-11.504	56.8393	3.96409,3.99599,3.99608,3.99603,3.99633,3.99604,3.99966
NelderMead_1_Levy7	-245	-11.504	56.8392	3.96439,3.99618,3.9963,3.99618,3.99616,3.99622,3.99946
NelderMead_0_Levy7	-356	-11.504	56.8392	3.9639,3.99615,3.99625,3.99624,3.99622,3.99625,3.99946
NelderMead_1_Levy7	-524	-11.504	56.8392	3.96385,3.99615,3.99623,3.99624,3.99625,3.99624,3.99946
NelderMead_0_Griewank	-40	0	4.91142	100.48,-97.6503
NelderMead_1_Griewank	-36	0	4.91141	100.483,-97.646
NelderMead_0_Griewank	-87	0	4.91141	100.481,-97.6459
NelderMead_1_Griewank	-85	0	4.91141	100.48,-97.6456
NelderMead_0_SixHumpCamel	148	-1.03	-1.03163	0.0898672,-0.712641
NelderMead_1_SixHumpCamel	89	-1.03	-1.03163	0.0898139,-0.712667
NelderMead_0_SixHumpCamel	209	-1.03	-1.03163	0.0898672,-0.712641
NelderMead_1_SixHumpCamel	153	-1.03	-1.03163	0.0898139,-0.712667
NelderMead_0_Branin	58	0.397889	0.397887	9.42484,2.47488
NelderMead_1_Branin	68	0.397889	0.397887	9.4249,2.47518
NelderMead_0_Branin	114	0.397889	0.397887	9.42476,2.47496
NelderMead_1_Branin	121	0.397889	0.397887	-3.14158,12.275
NelderMead_0_Shubert	-59	-24.06	-14.6909	8.82759,5.79172
NelderMead_1_Shubert	-51	-24.06	-7.21602	6.86033,6.86023
NelderMead_0_Shubert	128	-24.06	-24.0625	5.79182,5.7918
NelderMead_1_Shubert	-119	-24.06	-7.21603	6.86014,6.86014
NelderMead_0_Hansen	60	-176.54	-176.542	4.9766,4.85798
NelderMead_1_Hansen	58	-176.54	-176.542	4.97637,4.85802
NelderMead_0_Hansen	128	-176.54	-176.542	4.97649,4.85806
NelderMead_1_Hansen	123	-176.54	-176.542	4.97646,4.85804
NelderMead_0_Cola	-291	12.8154	256.203	1.93387,-0.57903,0.373911,0.0739108,0.0739108,0.0739108,0.0739108,0.0739108,0.0739108,0.673911,-0.0382,0.0739108,0.0739108,0.0739108,0.0739108,0.0739108,0.0739108
NelderMead_1_Cola	-270	12.8154	269.106	2.12345,-0.032006,0.248757,0.160369,0.160369,0.160369,0.160369,0.160369,0.160369,0.690699,0.127338,0.160369,0.160369,0.160369,0.160369,0.160369,0.160369
NelderMead_0_Cola	-773	12.8154	169.125	1.80169,-1.37731,0.670804,0.370804,0.370804,0.370804,0.370804,0.370804,0.370804,0.970804,-1.97723,0.370804,0.370804,0.370804,0.370804,0.370804,0.370804
NelderMead_1_Cola	-750	12.8154	197.195	2.17655,-1.25365,0.566106,0.477718,0.477718,0.477718,0.477718,0.477718,0.477718,1.00805,0.444687,0.477718,0.477718,0.477718,0.477718,0.477718,0.477718
NelderMead_0_Ackley	-47	0	19.3325	16.9991,16.999
NelderMead_1_Ackley	-42	0	19.3325	16.9997,16.9983
NelderMead_0_Ackley	-107	0	19.3325	16.9987,16.9988
NelderMead_1_Ackley	-99	0	19.3325	16.9988,16.9988
NelderMead_0_Bohachevsky1	107	0	2.80925e-08	3.68721e-05,1.60197e-05
NelderMead_1_Bohachevsky1	107	0	1.44232e-08	-2.52107e-05,1.25853e-05
NelderMead_0_Bohachevsky1	172	0	2.09614e-08	-3.76885e-05,-4.28061e-06
NelderMead_1_Bohachevsky1	171	0	7.07406e-09	1.04185e-05,-1.28198e-05
NelderMead_0_Bohachevsky2	113	0	1.54786e-08	-1.30841e-05,-2.25193e-05
NelderMead_1_Bohachevsky2	116	0	1.84323e-08	3.537e-05,-4.46684e-06
NelderMead_0_Bohachevsky2	180	0	8.38959e-09	1.91356e-05,-1.10643e-05
NelderMead_1_Bohachevsky2	182	0	4.57689e-09	-1.74279e-05,2.96785e-06
NelderMead_0_Bohachevsky3	113	0	2.6074e-09	1.94913e-05,-2.19285e-05
NelderMead_1_Bohachevsky3	116	0	7.57306e-09	-7.13698e-06,2.19882e-05
NelderMead_0_Bohachevsky3	177	0	2.6074e-09	1.94913e-05,-2.19285e-05
NelderMead_1_Bohachevsky3	179	0	2.94437e-09	-2.99086e-06,-8.60461e-06
NelderMead_0_Easom	-39	-1	-0	25,25
NelderMead_1_Easom	-39	-1	-0	25,25
NelderMead_0_Easom	-78	-1	-0	25,25
NelderMead_1_Easom	-78	-1	-0	25,25
NelderMead_0_Rastrigin	-68	0	7.95966	1.98983,1.98988
NelderMead_1_Rastrigin	-63	0	7.95967	1.99003,1.98996
NelderMead_0_Rastrigin	-139	0	7.95966	1.9899,1.9899
NelderMead_1_Rastrigin	-128	0	7.95966	1.98993,1.98994
NelderMead_0_Michalewicz2	57	-1.8013	-1.8013	2.20294,1.57083
NelderMead_1_Michalewicz2	57	-1.8013	-1.8013	2.20286,1.57074
NelderMead_0_Michalewicz2	125	-1.8013	-1.8013	2.20292,1.5708
NelderMead_1_Michalewicz2	125	-1.8013	-1.8013	2.20291,1.57081
NelderMead_0_Michalewicz5	-265	-4.68766	-4.37489	2.20315,1.57071,2.21932,1.92308,0.996771
NelderMead_1_Michalewicz5	-581	-4.68766	-4.21132	2.20294,1.57064,1.2851,2.48209,0.996649
NelderMead_0_Michalewicz5	-465	-4.68766	-4.3749	2.20291,1.57081,2.21934,1.92306,0.996673
NelderMead_1_Michalewicz5	-769	-4.68766	-4.21132	2.20296,1.57068,1.28499,2.48202,0.996671
NelderMead_0_Michalewicz10	-757	-9.66015	-7.54152	2.20317,1.57105,1.28507,1.38683,1.7204,1.57082,2.2211,1.58688,1.65572,1.5708
NelderMead_1_Michalewicz10	-1309	-9.66015	-7.38385	2.24268,1.4873,1.58613,1.12782,1.71622,1.56826,1.45042,1.63713,1.65397,1.57083
NelderMead_0_Michalewicz10	-1491	-9.66015	-8.63743	2.2029,1.5708,1.285,1.92305,1.72048,1.57079,1.87724,1.9251,1.65572,1.57079
NelderMead_1_Michalewicz10	-3389	-9.66015	-9.36542	2.14938,1.51456,1.29891,1.1225,1.70948,1.56669,1.44722,1.75381,1.65608,1.57376
==1929==
==1929== HEAP SUMMARY:
==1929==     in use at exit: 0 bytes in 0 blocks
==1929==   total heap usage: 103,087 allocs, 103,087 frees, 15,091,548 bytes allocated
==1929==
==1929== All heap blocks were freed -- no leaks are possible
==1929==
==1929== For counts of detected and suppressed errors, rerun with: -v
==1929== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 4 from 4)
*/

#endif
