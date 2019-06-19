////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	v2_final_plots.cxx 
////	A23 diffuse, make plots of the final cut parameter space
////
////	Nov 2018
////////////////////////////////////////////////////////////////////////////////

//C++
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <sys/stat.h>

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TChain.h"

// analysis custom
#include "tools_Cuts.h"
#include "tools_Stats.h"
#include "tools_CommandLine.h"
#include "tools_outputObjects.h"

double func(double *x, double *p){
	double num = x[0];
	double signal = p[0];
	double background = p[1];
	double term1 = TMath::Exp(-(signal+background));
	double term2 = TMath::Power(signal + background,num);
	double term3 = TMath::Gamma(num+1.);
	return term1*term2/term3;
}

double total_livetime=1142.46; //total livetime of the experiment in days
// double live_frac[5]={
// 	0.15674,
// 	0.12477,
// 	0.08275,
// 	0.38428,
// 	0.25146
// };

double live_frac[5]={
	0.7,
	0.3,
	0.08275,
	0.38428,
	0.25146
};

double getSignal(double Sbin1, double Sbin2){
	Sbin1*=live_frac[0];
	Sbin2*=live_frac[1]*2.; //because Brian was dumb and threw half as many neutrinos in config 2 and config 3
	return Sbin1 + Sbin2;
}

void splitSup(double Suptotal, double(&fracs)[5]){
	for(int i=0; i<5; i++){
		fracs[i]=Suptotal*live_frac[i];
	}
}

/*
	Takes as arguemnts the livetime fraction
*/
double getProductProb(double(&SupFracs)[5],  double(&totalIntegrals)[5], TF1 *theFuncs[5]){
	double ratios[5] = {0.};
	for(int i=0; i<2; i++){
		// cout<<"Sup frac is "<<SupFracs[i]<<endl;
		ratios[i] = theFuncs[i]->Integral(SupFracs[i],50.)/(totalIntegrals[i]);
	}
	double totalProb=1.;
	for(int i=0; i<2; i++){
		totalProb*=ratios[i];
	}
	// cout<<"          Prob 1 "<<ratios[0]<<" x prob 2 "<<ratios[1]<<" = "<<totalProb<<endl;
	return totalProb;
}



using namespace std;

int main(int argc, char **argv)
{
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	stringstream ss;
	gStyle->SetOptStat(11);

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

	if(argc<4){
		cout<< "Usage\n" << argv[0] << " <1-station> <2-pol> <3-slope>"<<endl;;
		return -1;
	}

	int station = atoi(argv[1]);
	int pol_select=atoi(argv[2]);
	double slope=double(atof(argv[3]));


	if(station!=2 && station!=3){
		printf("No good! You asked for station %d, but this code only works for stations 2 and 3 \n",station);
		return -1;
	}

	char inputfilename[500];
	sprintf(inputfilename,"A%d_optimize_pol%d_slope%d_save.root",station,pol_select,int(slope));
	TFile *InputFIle = TFile::Open(inputfilename,"READ");
	TTree *InputTree = (TTree*) InputFIle->Get("OutputTree");

	double intercept[30];
	double backgrounds[30];
	double signal[30];
	double pVal;
	int fitConverge;

	InputTree->SetBranchAddress("intercept",&intercept);
	InputTree->SetBranchAddress("backgrounds",&backgrounds);
	InputTree->SetBranchAddress("signal",&signal);
	InputTree->SetBranchAddress("fitConverge",&fitConverge);
	InputTree->SetBranchAddress("pVal",&pVal);

	// loop over configs, and get the info stored into these larger "mega" arrays, so we don't have to keep asking root to get stuff when we do the grid search
	// c++/the computer has more than enough memory to hold all these numbers at once

	double backgrounds_all[5][30];
	double signal_all[5][30];
	double pVal_all[5];
	int fitConverge_all[5];

	bool allConverge=true;
	bool allGood=true;

	for(int config=0; config<5; config++){

		InputTree->GetEntry(config);
		for(int bin=0; bin<30; bin++){
			backgrounds_all[config][bin] = backgrounds[bin];
			signal_all[config][bin] = signal[bin];
		}
		pVal_all[config] = pVal;
		fitConverge_all[config] = fitConverge;
		if(fitConverge!=1){
			allConverge=false;
		}
		if(pVal<0.05 || pVal>0.99){
			allGood=false;
		}

	}
	// now, we need to check for fit convergence and good p-values in all configurations
	// if not, we stop
	allConverge=true;
	if(!allConverge){
		printf("Not all fits converged! abandon ship! Convergence was: %d, %d, %d, %d, %d \n",fitConverge_all[0], fitConverge_all[1], fitConverge_all[2], fitConverge_all[3], fitConverge_all[4]);
		return -1;
	}
	allGood=true;
	if(!allGood){
		printf("Not all p-vals are reasonable! abandon ship! Values were: %.4f, %.4f, %.4f, %.4f, %.4f \n",pVal_all[0], pVal_all[1], pVal_all[2], pVal_all[3], pVal_all[4]);
		return -1;
	}

	stringstream ss1;

	TF1 *funcs[5];
	for(int i=0; i<5; i++){
		ss1.str("");
		ss1<<"func"<<i;
		funcs[i] = new TF1(ss1.str().c_str(),func,0,50,2);
		funcs[i]->SetParameters(0,0.); //always initiliaze signal to zero
	}

	// f.SetParameters(0.,0.7);
	// double total_integral = f.Integral(0,30);
	// double partial_integral = f.Integral(2,30);
	// cout<<"Probability is "<<partial_integral/total_integral<<endl;

	double total_integrals[5];

	int best_bins[5] = {0};
	double best_S_over_Sup = -10000.;

	// if all converge and all are good, start grid search (ugh);
	int numBins = 30;
	for(int bin1=0; bin1<numBins; bin1++){ // loop over config 1 bins
		
		double back1 = backgrounds_all[0][bin1];
		if(back1>1.) continue; // no background this large in single bin will ever be viable as part of final answer
		funcs[0]->SetParameter(1,back1); // set the mean of the background Poisson
		total_integrals[0] = funcs[0]->Integral(0,50.);
		// printf("For bin1 %d, at intercept %.2f we find %.6f num backgrounds, and total integral of %.6f \n",bin1,intercept[bin1],back1, total_integrals[0]);
		
		for(int bin2=0; bin2<numBins+9; bin2++){ // loop over config 2 bins
			
			double back2 = backgrounds_all[1][bin2];
			if(back2>1.) continue; // no background this large in single bin will ever be viable as part of final answer
			funcs[1]->SetParameter(1,back2); // set the mean of the background Poisson
			total_integrals[1] = funcs[1]->Integral(0,50.);
			// printf("     For bin2 %d, at intercept %.2f, we find %.6f num backgrounds , and total integral of %.6f \n",bin2,intercept[bin2],back2, total_integrals[1]);

			// adaptive step size algorithm for locating the ideal sup
			// this should be much faster than a brute fore scan

			double possibleSup = 100.;
			double fractionalSup[5];
			splitSup(possibleSup,fractionalSup);
			double goal = 0.1;
			double tolerance = 0.001;
			double thisProb = getProductProb(fractionalSup,total_integrals,funcs);
			double diff =  thisProb - goal;

			double step = possibleSup/2.;
			int counter = 0.;
			int maxNumTries=100;

			while(
				abs(diff) > tolerance
				&&
				counter<maxNumTries
			){
				// printf("                    For Try %d, prob is %.8f, which has a diff of %.8f\n",counter,thisProb,diff);
				if(thisProb>goal){
					// probability is too high! make sup smaller
					possibleSup+=step;
					step*=0.5;
				}
				else{
					// probability is too small! make sup bigger
					possibleSup-=step;
					step*=0.5;
				}
				splitSup(possibleSup,fractionalSup); // with our adjusted Sup, re-divvy up the Sup
				thisProb = getProductProb(fractionalSup,total_integrals,funcs);
				diff = thisProb - goal; // and check again how far off we are
				counter++; //don't forget to increment counter!
			}
			// cout<<"               I think I found the optimum with an Sup of "<<possibleSup<<" in only "<<counter<<" counts"<<endl;
			double totalSignal = getSignal(signal_all[0][bin1],signal_all[1][bin2]); // this function livetime adjusts the flux and fixes a numerical mistake brian made
			double this_S_over_Sup = totalSignal/possibleSup;
			if(this_S_over_Sup>best_S_over_Sup){
				// printf(GREEN"Found new optimum! Old is %.4f and new is S/Sup = %.4f/%.4f = %.4f \n"RESET,best_S_over_Sup, totalSignal, possibleSup, this_S_over_Sup);
				best_S_over_Sup = this_S_over_Sup;
				best_bins[0]=bin1;
				best_bins[1]=bin2;
				// printf(GREEN"     New optimum has intercepts %.2f, %.2f\n"RESET, intercept[bin1], intercept[bin2]);
				// printf("     New optimum has intercepts %.2f, %.2f, %.2f, %.2f, %.2f \n", intercept[bin1], intercept[bin2], intercept[bin3], intercept[bin4], intercept[bin5]);
			}


			// for(int bin3=0; bin3<numBins; bin3++){ // loop over config 3 bins
			// 	double back3 = backgrounds_all[2][bin3];
			// 	if(back3>1.) continue; // no background this large in single bin will ever be viable as part of final answer

			// 	for(int bin4=0; bin4<numBins; bin4++){ // loop over config 4 bins
			// 		double back4 = backgrounds_all[3][bin3];
			// 		if(back4>1.) continue; // no background this large in single bin will ever be viable as part of final answer
					
			// 		for(int bin5=0; bin5<numBins; bin5++){ // loop over config 5 bins
			// 			double back5 = backgrounds_all[4][bin3];
			// 			if(back5>1.) continue; // no background this large in single bin will ever be viable as part of final answer
					
			// 		} // loop over intercepts in config 5
			// 	} // loop over intercepts in config 4
			// } // loop over intercepts in config 3
		} // loop over intercepts in config 2
	} // loop over intercepts in config 1

	printf("I think I have the global S/Sup optimum of %.4f, at intercepts %.2f and %.2f \n",best_S_over_Sup,intercept[best_bins[0]],intercept[best_bins[1]]);

}