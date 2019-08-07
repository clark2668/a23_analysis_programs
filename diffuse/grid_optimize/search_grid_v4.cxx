////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////    This time exhaustive search grid w/ S/Sup, where Sup is Sup for sum of backgrounds
////	This time start our grid scan in the neighborhood of the right answer to save time
////	And, to boot, we will use proper background ucertainties
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
#include "tools_Stats.h"

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


// compute a number of backgrounds given the fit results, their errors, and their correlation
double getRandBack(TRandom *rand, double cut, double par0, double par0Err, double par1, double par1Err, double correlation){
	double rand_slope = rand->Gaus();
	double new_slope = par0 
						+ (par0Err*rand_slope);

	double rand_amp = rand->Gaus();
	double new_amp = 
					(par1)
					+ (par1Err * rand_slope * correlation)
					+ (par1Err * rand_amp * sqrt(1-pow(correlation,2.)));

	double new_back = 10.*(1./(new_slope*0.1)) * (-exp(new_slope*cut + new_amp)); // the 0.1 in the binWidthIntercept, which I know is 0.1....
	return new_back;
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
	double fitPar[2];
	double fitParErr[2];
	double covariance;
	double correlation;

	InputTree->SetBranchAddress("intercept",&intercept);
	InputTree->SetBranchAddress("backgrounds",&backgrounds);
	InputTree->SetBranchAddress("signal",&signal);
	InputTree->SetBranchAddress("fitConverge",&fitConverge);
	InputTree->SetBranchAddress("pVal",&pVal);
	InputTree->SetBranchAddress("fitPar0",&fitPar[0]);
	InputTree->SetBranchAddress("fitPar1",&fitPar[1]);
	InputTree->SetBranchAddress("fitParErr0",&fitParErr[0]);
	InputTree->SetBranchAddress("fitParErr1",&fitParErr[1]);
	InputTree->SetBranchAddress("covariance",&covariance);
	InputTree->SetBranchAddress("correlation",&correlation);

	// loop over configs, and get the info stored into these larger "mega" arrays, so we don't have to keep asking root to get stuff when we do the grid search
	// c++/the computer has more than enough memory to hold all these numbers at once

	cout<<"About to load lots of information, sit tight..."<<endl;

	double backgrounds_all[5][200];
	double signal_all[5][200];
	double pVal_all[5];
	int fitConverge_all[5];
	double fitPar_all[5][2];
	double fitParErr_all[5][2];
	double covariance_all[5];
	double correlation_all[5];

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
		fitPar_all[config][0] = fitPar[0];
		fitPar_all[config][1] = fitPar[1];
		fitParErr_all[config][0] = fitParErr[0];
		fitParErr_all[config][1] = fitParErr[1];
		covariance_all[config] = covariance;
		correlation_all[config] = covariance/(fitParErr[0]*fitParErr[1]); //build this by hand, because I'm an idiot...
		// correlation_all[config] = correlation; // at some point in the future, you can risk using this
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

	// just to quickly test our new procedure

	TRandom3 *random = new TRandom3(0);
	// TRandom *random = new TRandomMixMax();

	bool doTest = false;
	if(doTest){

		TH1D *hBackgrounds[5];
		for(int i=0; i<5; i++){
			hBackgrounds[i] = new TH1D("","",100,0,0.35);
		}
		TH1D *hBackgroundTot = new TH1D("","",100,0,0.35);
		TH1D *hSup = new TH1D("","",100,0,5);

		// this should make Sup higher for same choice of cuts, because we are allowing upward fluctuations
		int numExperiments = 10000;
		double average_Sup=0.;
		for(int i=0; i<numExperiments; i++){

			int bin1 = 21;
			int bin2 = 21;
			int bin3 = 21;
			int bin4 = 21;
			int bin5 = 21;

			double back1_test = getRandBack(random, intercept[bin1], fitPar_all[0][0], fitParErr_all[0][0], fitPar_all[0][1], fitParErr_all[0][1], correlation_all[0]);
			double back2_test = getRandBack(random, intercept[bin2], fitPar_all[1][0], fitParErr_all[1][0], fitPar_all[1][1], fitParErr_all[1][1], correlation_all[1]);
			double back3_test = getRandBack(random, intercept[bin3], fitPar_all[2][0], fitParErr_all[2][0], fitPar_all[2][1], fitParErr_all[2][1], correlation_all[2]);
			double back4_test = getRandBack(random, intercept[bin4], fitPar_all[3][0], fitParErr_all[3][0], fitPar_all[3][1], fitParErr_all[3][1], correlation_all[3]);
			double back5_test = getRandBack(random, intercept[bin5], fitPar_all[4][0], fitParErr_all[4][0], fitPar_all[4][1], fitParErr_all[4][1], correlation_all[4]);

			hBackgrounds[0]->Fill(back1_test);
			hBackgrounds[1]->Fill(back2_test);
			hBackgrounds[2]->Fill(back3_test);
			hBackgrounds[3]->Fill(back4_test);
			hBackgrounds[4]->Fill(back5_test);

			// sum backgrounds
			double total_back = back1_test + back2_test + back3_test + back4_test + back5_test;
			// printf("Back estimates are %.3f, %.3f, %.3f, %.3f, %.3f \n", back1_test, back2_test, back3_test, back4_test, back5_test );

			// cout<<"On experiment"<<i<<" with total background "<<total_back<<endl;

			hBackgroundTot->Fill(total_back);

			// double this_Sup;
			// // don't waste time (and not a small amount of time, mind you) actually computing Sup
			// // if we know the answer is going to be 2.3 anyway because the background is so far below zero
			// if(total_back<1e-5){
			// 	this_Sup=2.3;
			// }
			// else{
			// 	double achieved_alpha;
			// 	// this_Sup = GetS_up(total_back,achieved_alpha,0.9);
			// 	this_Sup = GetS_up_TMath(total_back,achieved_alpha,0.9);
			// }
			
			double achieved_alpha;
			double this_Sup = GetS_up_TMath(total_back,achieved_alpha,0.9);

			average_Sup+=this_Sup;
			hSup->Fill(this_Sup);
		} // loop pseudo experiments
		average_Sup/=double(numExperiments);
		cout<<"Average Sup is "<<average_Sup<<endl;

		gStyle->SetOptStat(111111);
		TCanvas *cVerify = new TCanvas("","",5*850,2*850);
		cVerify->Divide(5,2);
		for(int i=0; i<5; i++){
			cVerify->cd(i+1);
			hBackgrounds[i]->Draw("");
		}
		cVerify->cd(6);
			hBackgroundTot->Draw("");
		cVerify->cd(7);
			hSup->Draw("");
			hSup->SetLineWidth(2);
			gPad->SetLogy();
		cVerify->SaveAs("validate_new_grid_search.png");

		// cout<<"By the way, our random seed was "<<random->GetSeed()<<endl;

		return 0;
	}

	int best_bins[5] = {0};
	double best_S_over_Sup = -10000.;
	int numBins = 200;

	int start_bins[5] = {0};

	// let's work up a reasonable "first guess" for each, so we know the range to scan and save time
	// normally, we'd start by just choosing the best S/Sup we can find for each individual configuration, so let's do that

	if(allConverge && allGood){
		for(int config=0; config<5; config++){
			double testS[numBins];
			double testSup[numBins];
			double testSoverSup[numBins];
			double localBestSoverSup=-10;
			int best_bin;
			for(int bin=0; bin<numBins; bin++){
				double thisBack  = backgrounds_all[config][bin];
				if(thisBack<1e-6 || thisBack>1.){
					continue;
				}
				double thisSup;
				if(thisBack<1e-4){
					thisSup=2.31;
				}
				else{
					double achieved_alpha;
					thisSup = GetS_up(thisBack,achieved_alpha,0.9);
				}
				testS[bin] = signal_all[config][bin];
				testSup[bin]=thisSup;
				testSoverSup[bin]=testS[bin]/testSup[bin];
				if(testSoverSup[bin]>localBestSoverSup){
					best_bin=bin;
					localBestSoverSup=testSoverSup[bin];
				}
				// printf("For intercept %.2f, background is %.5f and Sup is %.5f and S is %.2f and S/Sup is %.3f \n", intercept[bin],thisBack, thisSup, testS[bin], testSoverSup[bin]);
			}
			// printf("For config %d, global best S/Sup is %.3f at bin %d and intercept %.2f \n",config,localBestSoverSup,best_bin, intercept[best_bin]);
			start_bins[config]=best_bin;
		}
		for(int config=0; config<5; config++){
			printf("Start bin for config %d is bin %d and intercept %.2f \n", config,start_bins[config],intercept[start_bins[config]]);
		}
	}

	if(allConverge && allGood){
		// if all converge and all are good, start grid search (ugh);
		// for(int bin1=0; bin1<numBins; bin1++){ // loop over config 1 bins
		for(int bin1=start_bins[0]+0; bin1<start_bins[0]+15; bin1++){ // loop over config 1 bins
		  	cout<<"Slope "<<slope<<" on bin1 bin "<<bin1<<" which is intercept "<<intercept[bin1]<<endl;
		  	// printf("Slope %.2f, on bin1 bin %d which is %.2f \n", slope, bin1, intercept[bin1]);
			double back1 = backgrounds_all[0][bin1];
			if(back1>1. || back1<1e-6) continue; // no background this large/small in single bin will ever be viable as part of final answer
			
			for(int bin2=start_bins[1]+0; bin2<start_bins[1]+15; bin2++){ // loop over config 2 bins
				double back2 = backgrounds_all[1][bin2];
				if(back2>1. || back2<1e-6) continue; // no background this large/small in single bin will ever be viable as part of final answer

				for(int bin3=start_bins[2]+0; bin3<start_bins[2]+15; bin3++){ // loop over config 3 bins
					double back3 = backgrounds_all[2][bin3];
					if(back3>1. || back3<1e-6) continue; // no background this large/small in single bin will ever be viable as part of final answer

					for(int bin4=start_bins[3]+0; bin4<start_bins[3]+15; bin4++){ // loop over config 4 bins
						double back4 = backgrounds_all[3][bin4];
						if(back4>1. || back4<1e-6) continue; // no background this large/small in single bin will ever be viable as part of final answer

						for(int bin5=start_bins[4]+0; bin5<start_bins[4]+15; bin5++){ // loop over config 5 bins
							// cout<<"     Slope "<<slope<<" on bin5 bin "<<bin5<<" which is intercept "<<intercept[bin5]<<endl;
							double back5 = backgrounds_all[4][bin5];
							if(back5>1. || back5<1e-6) continue; // no background this large/small in single bin will ever be viable as part of final answer

							// we needing das pseudo-experiments!
							// which allows the backgrounds to fluctuate with their fit uncertainties
							// so we can accurately sample the average S/Sup, instead of assuming we only get the average fit result all the time
							// this should make Sup higher for same choice of cuts, because we are allowing upward fluctuations
							int numExperiments = 7500;
							double average_Sup=0.;
							for(int i=0; i<numExperiments; i++){

								double back1_test = getRandBack(random, intercept[bin1], fitPar_all[0][0], fitParErr_all[0][0], fitPar_all[0][1], fitParErr_all[0][1], correlation_all[0]);
								double back2_test = getRandBack(random, intercept[bin2], fitPar_all[1][0], fitParErr_all[1][0], fitPar_all[1][1], fitParErr_all[1][1], correlation_all[1]);
								double back3_test = getRandBack(random, intercept[bin3], fitPar_all[2][0], fitParErr_all[2][0], fitPar_all[2][1], fitParErr_all[2][1], correlation_all[2]);
								double back4_test = getRandBack(random, intercept[bin4], fitPar_all[3][0], fitParErr_all[3][0], fitPar_all[3][1], fitParErr_all[3][1], correlation_all[3]);
								double back5_test = getRandBack(random, intercept[bin5], fitPar_all[4][0], fitParErr_all[4][0], fitPar_all[4][1], fitParErr_all[4][1], correlation_all[4]);

								// cout<<"back 1 is "<<back1_test<<endl;
								// cout<<"back 2 is "<<back2_test<<endl;
								// cout<<"back 3 is "<<back3_test<<endl;
								// cout<<"back 4 is "<<back4_test<<endl;
								// cout<<"back 5 is "<<back5_test<<endl;

								// cout<<"Fit par for 3 is "<<fitPar_all[2][0]<<", "<<fitParErr_all[2][0]<<" , "<<fitPar_all[2][1]<<" , "<<fitParErr_all[2][1]<<" , "<<correlation_all[2]<<endl;
								// cout<<"Fit par for 2 is "<<fitPar_all[1][0]<<", "<<fitParErr_all[1][0]<<" , "<<fitPar_all[1][1]<<" , "<<fitParErr_all[1][1]<<" , "<<correlation_all[1]<<endl;

								// sum backgrounds
								double total_back = back1_test + back2_test + back3_test + back4_test + back5_test;									
								double this_Sup;
								// // don't waste time (and not a small amount of time, mind you) actually computing Sup
								// // if we know the answer is going to be 2.3 anyway because the background is so far below zero
								if(total_back<3e-3){
									this_Sup=2.31;
								}
								else{
									double achieved_alpha;
									// this_Sup = GetS_up(total_back,achieved_alpha,0.9);
									this_Sup = GetS_up_TMath(total_back,achieved_alpha,0.9);
								}

								// compute Sup with much faster analytic technique; don't fail me now Eugene...
								// double achieved_alpha;
								// double this_Sup = GetS_up_TMath(total_back,achieved_alpha,0.9);
								average_Sup+=this_Sup;
								// cout<<"total back is "<<total_back<<endl;
								// cout<<"This sup is "<<this_Sup<<endl;
							} // loop pseudo experiments
							average_Sup/=double(numExperiments);
							// cout<<"          This Average Sup is "<<average_Sup<<endl;

							double total_signal = getSignal(signal_all[0][bin1],signal_all[1][bin2],signal_all[2][bin3],signal_all[3][bin4],signal_all[4][bin5]); // this function livetime adjusts the flux and fixes a numerical mistake brian made
							double this_S_over_Sup = total_signal/average_Sup;
							// cout<<"          This S is "<<total_signal<<endl;
							// cout<<"          This S over Sup is "<<this_S_over_Sup<<endl;
							printf("%d : %d : %d : %d : %d : %4.2f \n", bin1, bin2, bin3, bin4, bin5, this_S_over_Sup);

							if(this_S_over_Sup>best_S_over_Sup){
								best_S_over_Sup = this_S_over_Sup;
								best_bins[0]=bin1;
								best_bins[1]=bin2;
								best_bins[2]=bin3;
								best_bins[3]=bin4;
								best_bins[4]=bin5;
								printf(" New optimum has S/Sup of %.3f and intercepts %.2f, %.2f, %.2f, %.2f, %.2f at bins %d, %d, %d, %d, %d \n", 
											best_S_over_Sup, 
											intercept[bin1], intercept[bin2], intercept[bin3], intercept[bin4], intercept[bin5], 
											bin1, bin2, bin3, bin4, bin5);
							}
						} // loop over intercepts in config 5
					} // loop over intercepts in config 4
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

	printf("Slope %.2f has backgrounds %.7f, %.7f, %.7f, %.7f, %.7f \n",
																	slope,
																	backgrounds_all[0][best_bins[0]],
																	backgrounds_all[1][best_bins[1]],
																	backgrounds_all[2][best_bins[2]],
																	backgrounds_all[3][best_bins[3]],
																	backgrounds_all[4][best_bins[4]]);
	
	// write this to file so I can study it later, and dont' forget the random number seed!
	char outfile_name[200];
	sprintf(outfile_name,"/users/PAS0654/osu0673/A23_analysis_new2/optimize/A%d_optimize_pol%d_new.txt",station,pol_select);
	FILE *fout = fopen(outfile_name, "a");
	fprintf(fout,"%5.1f, %4.2f, %2.2f, %2.2f, %2.2f, %2.2f, %2.2f, %2d, %2d, %2d, %2d, %2d, %1.5f, %1.5f, %1.5f, %1.5f, %1.5f, %d \n",
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
				backgrounds_all[4][best_bins[4]],
				random->GetSeed()
	);	
	fclose(fout);
}