////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////    This time exhaustive search grid w/ S/Sup, where Sup is Sup for sum of backgrounds
////
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

double total_livetime=1142.46; //total livetime of the experiment in days
double live_frac[5]={
	0.15674,
	0.12477,
	0.08275,
	0.38428,
	0.25146
};

// for testing
// double live_frac[5]={
// 	0.21,
// 	0.17,
// 	0.11,
// 	0.51,
// 	0.25146
// };

double getSignal(double Sbin1, double Sbin2, double Sbin3, double Sbin4, double Sbin5){
	Sbin1*=live_frac[0];
	Sbin2*=live_frac[1]*2.; //because Brian was dumb and threw half as many neutrinos in config 2 and config 3
	Sbin3*=live_frac[2]*2.; //because Brian was dumb and threw half as many neutrinos in config 2 and config 3
	Sbin4*=live_frac[3];
	Sbin5*=live_frac[4];
	return Sbin1 + Sbin2 + Sbin3 + Sbin4 + Sbin5;
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
	sprintf(inputfilename,"/fs/scratch/PAS0654/ara/optimize/A%d_optimize_pol%d_slope%d.root",station,pol_select,int(slope));
	TFile *InputFIle = TFile::Open(inputfilename,"READ");
	TTree *InputTree = (TTree*) InputFIle->Get("OutputTree");

	double intercept[200];
	double backgrounds[200];
	double signal[200];
	double pVal;
	int fitConverge;

	InputTree->SetBranchAddress("intercept",&intercept);
	InputTree->SetBranchAddress("backgrounds",&backgrounds);
	InputTree->SetBranchAddress("signal",&signal);
	InputTree->SetBranchAddress("fitConverge",&fitConverge);
	InputTree->SetBranchAddress("pVal",&pVal);

	// loop over configs, and get the info stored into these larger "mega" arrays, so we don't have to keep asking root to get stuff when we do the grid search
	// c++/the computer has more than enough memory to hold all these numbers at once

	double backgrounds_all[5][200];
	double signal_all[5][200];
	double pVal_all[5];
	int fitConverge_all[5];

	bool allConverge=true;
	bool allGood=true;

	for(int config=0; config<5; config++){
		InputTree->GetEntry(config);
		for(int bin=0; bin<200; bin++){
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
	// allConverge=true;
	// if(!allConverge){
	// 	printf("Not all fits converged! abandon ship! Convergence was: %d, %d, %d, %d, %d \n",fitConverge_all[0], fitConverge_all[1], fitConverge_all[2], fitConverge_all[3], fitConverge_all[4]);
	// 	return -1;
	// }
	// allGood=true;
	// if(!allGood){
	// 	printf("Not all p-vals are reasonable! abandon ship! Values were: %.4f, %.4f, %.4f, %.4f, %.4f \n",pVal_all[0], pVal_all[1], pVal_all[2], pVal_all[3], pVal_all[4]);
	// 	return -1;
	// }

	int best_bins[5] = {0};
	double best_S_over_Sup = -10000.;
	int numBins = 80;

	if(allConverge && allGood){
		// if all converge and all are good, start grid search (ugh);
		for(int bin1=0; bin1<numBins; bin1++){ // loop over config 1 bins
		  	cout<<"Slope "<<slope<<" on bin1 bin "<<bin1<<" which is intercept "<<intercept[bin1];
		        //printf("Slope %.2f, on bin1 bin %d which is %.2f \n", slope, bin1, intercept[bin1]);
			double back1 = backgrounds_all[0][bin1];
			if(back1>1. || back1<1e-7) continue; // no background this large/small in single bin will ever be viable as part of final answer
			
			for(int bin2=0; bin2<numBins; bin2++){ // loop over config 2 bins
				double back2 = backgrounds_all[1][bin2];
				if(back2>1. || back2<1e-7) continue; // no background this large/small in single bin will ever be viable as part of final answer

				for(int bin3=0; bin3<numBins; bin3++){ // loop over config 3 bins
					double back3 = backgrounds_all[2][bin3];
					if(back3>1. || back3<1e-7) continue; // no background this large/small in single bin will ever be viable as part of final answer

						for(int bin4=0; bin4<numBins; bin4++){ // loop over config 4 bins
							double back4 = backgrounds_all[3][bin4];
							if(back4>1.|| back4<1e-7) continue; // no background this large/small in single bin will ever be viable as part of final answer

								for(int bin5=0; bin5<numBins; bin5++){ // loop over config 5 bins
									double back5 = backgrounds_all[4][bin5];
									if(back5>1.|| back5<1e-7) continue; // no background this large/small in single bin will ever be viable as part of final answer
								
									double achieved_alpha;
									double this_Sup = GetS_up(back1+back2+back3+back4+back5,achieved_alpha,0.9);
									// printf("For back %.2f and %.2f, Sup is %.2f \n",back1,back2,this_Sup );

									double totalSignal = getSignal(signal_all[0][bin1],signal_all[1][bin2],signal_all[2][bin3],signal_all[3][bin4],signal_all[4][bin5]); // this function livetime adjusts the flux and fixes a numerical mistake brian made
									
									double this_S_over_Sup = totalSignal/this_Sup;

									// printf("               Now, S/Sup = %.2f/%.2f = %.2f\n", totalSignal, possibleSup, this_S_over_Sup);
									if(this_S_over_Sup>best_S_over_Sup){
										// printf(GREEN"     Found new optimum! Old is %.4f and new is S/Sup = %.4f/%.4f = %.4f \n"RESET,best_S_over_Sup, totalSignal, possibleSup, this_S_over_Sup);
										best_S_over_Sup = this_S_over_Sup;
										best_bins[0]=bin1;
										best_bins[1]=bin2;
										best_bins[2]=bin3;
										best_bins[3]=bin4;
										best_bins[4]=bin5;
										// printf(GREEN"          New optimum has intercepts %.2f, %.2f\n"RESET, intercept[bin1], intercept[bin2]);
										// printf("     New optimum has intercepts %.2f, %.2f, %.2f, %.2f, %.2f \n", intercept[bin1], intercept[bin2], intercept[bin3], intercept[bin4], intercept[bin5]);
									}

								} // loop over intercepts in config 5
						}
				} // loop over intercepts in config 3
			} // loop over intercepts in config 2
		} // loop over intercepts in config 1
	}

	printf("Slope %.2f has global S/Sup optimum of %.4f, at intercepts %.2f and %.2f and %.2f and %.2f and %.2f \n",
														slope,
														best_S_over_Sup,
														intercept[best_bins[0]],
														intercept[best_bins[1]],
														intercept[best_bins[2]],
														intercept[best_bins[3]],
														intercept[best_bins[4]]);

	printf("Slope %.2f has optimum bins %d, %d, %d, %d, %d \n", 
												slope,
												best_bins[0], 
												best_bins[1], 
												best_bins[2], 
												best_bins[3], 
												best_bins[4]);

	printf("Slope %.2f has backgrounds are %.7f, %.7f, %.7f, %.47f, %.7f \n",
																	slope,
																	backgrounds_all[0][best_bins[0]],
																	backgrounds_all[1][best_bins[1]],
																	backgrounds_all[2][best_bins[2]],
																	backgrounds_all[3][best_bins[3]],
																	backgrounds_all[4][best_bins[4]]);
	
	// write this to file so I can study it later
	char outfile_name[200];
	sprintf(outfile_name,"/users/PAS0654/osu0673/A23_analysis_new2/optimize/A%d_optimize_pol%d.txt",station,pol_select);
	FILE *fout = fopen(outfile_name, "a");
	fprintf(fout,"%4.1f, %4.2f, %2.2f, %2.2f, %2.2f, %2.2f, %2.2f, %2d, %2d, %2d, %2d, %2d, %1.4f, %1.4f, %1.4f, %1.4f, %1.4f \n",
				slope,
				best_S_over_Sup,
				intercept[best_bins[0]],
				intercept[best_bins[1]],
				intercept[best_bins[2]],
				intercept[best_bins[3]],
				intercept[best_bins[4]],
				best_bins[0], 
				best_bins[1], 
				best_bins[2], 
				best_bins[3], 
				best_bins[4],
				backgrounds_all[0][best_bins[0]],
				backgrounds_all[1][best_bins[1]],
				backgrounds_all[2][best_bins[2]],
				backgrounds_all[3][best_bins[3]],
				backgrounds_all[4][best_bins[4]]
	);	
	fclose(fout);
}
