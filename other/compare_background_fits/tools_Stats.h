
//Includes
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <complex>
#include <algorithm>

#include "TGraph.h"
#include "TMath.h"

using namespace std;

// compute the log-likelihood 
double ReturnLogL_highN( double N_bin, double U_exp ) {
    double LogL; // -2 ln( L )
	if ( N_bin > 0. ) {
		LogL = -2.* ( N_bin * log( U_exp / N_bin ) - U_exp + N_bin );
	}
	else {
		LogL = -2.* ( -U_exp );
	}
	return LogL;
}

/*
	Functions for computing S_up
	Old method
*/

double lnGamma( double z, double x, int n=100 )
{
	if ( z>0 && x>0 && n>2 ) {
		double mult1_1;
		if (x>100.) {
			double step1_1 = x/100.;
			mult1_1 = step1_1 * log(exp(-1.*100.));
		}
		else {
			mult1_1 = log(exp(-1.*x));
		}

		double mult1_2;
		if (z>10.) {
			double step = z/10.;
			mult1_2 = step*log(pow(x,10.));
		}
		else {
			mult1_2 = log(pow(x,z));
		}

		//cout<<"log(mult1_1) : "<<mult1_1<<endl;
		//cout<<"log(mult1_2) : "<<mult1_2<<endl;

		double mult2 = 0.;

		for (int i=n; i>0; i--) {

			mult2 = (double)i*( (double)i - z ) / ( x + 2.*(double)i+1.-z - mult2 );
			//cout<<"at i:"<<i<<" mult2 : "<<mult2<<endl;
		}

		mult2 = 1. / ( x +1.-z - mult2 );
		mult2 = log( mult2 ); // change to log value

		//return mult1 * mult2;
		return mult1_1 + mult1_2 + mult2; // add log values which is multiple
	}
	else return 0.;
}

double Alpha_nb_ln( double s_up, double nb, int n=100 )
{
	if ( s_up>0 && nb>0 ) {
		double gamma_part = lnGamma(1.+nb, s_up+nb, n) - lnGamma(1.+nb, nb, n);
 		gamma_part = exp( gamma_part );
		return 1. - gamma_part;
	}
	else return 0.;
}


double GetS_up ( double ExpEvts, double &alpha_out, double alpha_cut, int n=100 ) {
	double s_up = -0.01;
	alpha_out = 0.;
	while ( alpha_out < alpha_cut ) {
		s_up += 0.01;
		//alpha = Alpha_nb_ln( s_up, ExpEvts, n );
		alpha_out = Alpha_nb_ln( s_up, ExpEvts, n );
	}
	return s_up;
}

/*
	Functions for computing S_up
	New Method, Using ROOT (because we have to anyway...)
*/

double Mathematica_Gamma( double a, double z ) {
	// Mathematic incomplete gamma function output
	return TMath::Gamma(a) * ( 1. - TMath::Gamma(a,z) );
}

double Get_TMath_Alpha( double s_up, double nb ) {
	// return simplified 1-alpha value from Mathematica
	return 1. - ( Mathematica_Gamma( 1.+nb, s_up+nb ) / Mathematica_Gamma( 1.+nb, nb ) );
}

double GetS_up_TMath ( double ExpEvts, double &alpha_out, double alpha_cut ) {
	double s_up = 0.;
	alpha_out = 0.;
	while ( alpha_out < alpha_cut ) {
		s_up += 0.01;
		alpha_out = Get_TMath_Alpha( s_up, ExpEvts );
	}
	return s_up;
}
