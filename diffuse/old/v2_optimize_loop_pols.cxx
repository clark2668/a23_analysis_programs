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

	if(argc<6){
		cout<< "Usage\n" << argv[0] << " <isSim> <station> <config> <year_or_energy> <ValForCuts filename>"<<endl;;
		return -1;
	}

	int isSim = atoi(argv[1]);
	int station = atoi(argv[2]);
	int config = atoi(argv[3]);
	double year_or_energy = double(atof(argv[4]));
	double slope=double(atof(argv[5]));


	if(station!=2 && station!=3){
		printf("No good! You asked for station %d, but this code only works for stations 2 and 3 \n",station);
		return -1;
	}

	vector<int> BadRunList=BuildBadRunList(station);


	/*
		Now, we also need a way to deal with the rotate cross correlation cut
	*/

	int numTotal=0;

	// just for V for a hot minute

	// let's see if we can do it with arrays....
	int numSNRbins=200;
	double numEventsPassed[2][numSNRbins];
	double numEventsPassed_diff[2][numSNRbins];

	double dRcut=0.1; // SNR bin incremets (bins in 0.1 SNR units)
	double Rcutmin=0.;
	double intercept[numSNRbins];
	for(int pol=0; pol<2; pol++){
		for(int bin=0; bin<numSNRbins; bin++){
			numEventsPassed[pol][bin]=0.;
			numEventsPassed_diff[pol][bin]=0.;
			intercept[bin] = Rcutmin + double(bin)*dRcut;
		}
	}
	
	// need to be able to make the final 2D distribution
	double max=0.05;
	TH2D *h2SNRvsCorr[2]; // SNR on Y axis, Corr on X axis, like in the TB
	h2SNRvsCorr[0]=new TH2D("","V Data",100,0,max,300,0,30);
	h2SNRvsCorr[1]=new TH2D("","H Data",100,0,max,300,0,30);

	// first, we go through the events once to construct the rotated cut
	for(int file_num=6; file_num<argc; file_num++){

		string chRun = "run";
		string file = string(argv[file_num]);
		size_t foundRun = file.find(chRun);
		string strRunNum = file.substr(foundRun+4,4);
		int runNum = atoi(strRunNum.c_str());

		TFile *inputFile = TFile::Open(argv[file_num]);
		if(!inputFile){
			cout<<"Can't open val for cuts file!"<<endl;
			return -1;
		}
		printf("Fitting loop data File %d: run %d \n", file_num, runNum);

		TTree *trees[4];
		trees[0] = (TTree*) inputFile->Get("VTree");
		trees[1] = (TTree*) inputFile->Get("HTree");
		trees[2] = (TTree*) inputFile->Get("AllTree");

		double corr_val[2];
		double snr_val[2];
		int WFRMS[2];
		double frac_of_power_notched_V[8];
		double frac_of_power_notched_H[8];
		int Refilt[2];

		trees[0]->SetBranchAddress("corr_val_V",&corr_val[0]);
		trees[0]->SetBranchAddress("snr_val_V",&snr_val[0]);
		trees[0]->SetBranchAddress("wfrms_val_V",&WFRMS[0]);
		trees[0]->SetBranchAddress("Refilt_V",&Refilt[0]);
		trees[1]->SetBranchAddress("corr_val_H",&corr_val[1]);
		trees[1]->SetBranchAddress("snr_val_H",&snr_val[1]);
		trees[1]->SetBranchAddress("wfrms_val_H",&WFRMS[1]);
		trees[1]->SetBranchAddress("Refilt_H",&Refilt[1]);

		int isCal;
		int isSoft;
		int isShort;
		int isCW;
		int isNewBox;
		int isSurf[2];
		int isBadEvent;
		double weight;
		int isSurfEvent_top[2];
		int unixTime;
		int isFirstFiveEvent;
		int hasBadSpareChanIssue;

		trees[2]->SetBranchAddress("cal",&isCal);
		trees[2]->SetBranchAddress("soft",&isSoft);
		trees[2]->SetBranchAddress("short",&isShort);
		trees[2]->SetBranchAddress("CW",&isCW);
		trees[2]->SetBranchAddress("box",&isNewBox);
		trees[2]->SetBranchAddress("surf_V",&isSurf[0]);
		trees[2]->SetBranchAddress("surf_H",&isSurf[1]);
		trees[2]->SetBranchAddress("bad",&isBadEvent);
		trees[2]->SetBranchAddress("weight",&weight);
		trees[2]->SetBranchAddress("surf_top_V",&isSurfEvent_top[0]);
		trees[2]->SetBranchAddress("surf_top_H",&isSurfEvent_top[1]);
		trees[2]->SetBranchAddress("unixTime",&unixTime);
		trees[2]->SetBranchAddress("isFirstFiveEvent",&isFirstFiveEvent);
		trees[2]->SetBranchAddress("hasBadSpareChanIssue",&hasBadSpareChanIssue);

		stringstream ss;
		for(int i=0; i<8; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			trees[0]->SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_V[i]);
		}
		for(int i=8; i<16; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			trees[1]->SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_H[i-8]);
		}
		
		int numEntries = trees[0]->GetEntries();
		numTotal+=numEntries;

		if(isBadRun(station,runNum,BadRunList)){
			continue;
		}

		//now to loop over events
		// numEntries=10000;
		for(int event=0; event<numEntries; event++){

			trees[0]->GetEvent(event);
			trees[1]->GetEvent(event);
			trees[2]->GetEvent(event);

			if( (isSoft || isBadEvent || hasBadSpareChanIssue || isFirstFiveEvent || isShort || isCal))
				continue;

			if(isBadLivetime(station,unixTime) && !isSim){
				continue;
			}
				
			for(int pol=0; pol<2; pol++){
				if(!WFRMS[pol] && !isNewBox && !isSurf[0] && !isSurf[1] && !isSurfEvent_top[pol]){
					bool failsCWPowerCut=false;
					if(Refilt[pol] && !WFRMS[pol]){
						vector<double> frac;
						for(int i=0; i<8; i++){
							if(pol==0) frac.push_back(frac_of_power_notched_V[i]);
							else if(pol==1) frac.push_back(frac_of_power_notched_H[i]);
						}
						sort(frac.begin(), frac.end(), std::greater<double>());
						if(frac[2]>0.06){
							failsCWPowerCut=true;
						}
					} //refiltered?
					if(!failsCWPowerCut){
						h2SNRvsCorr[pol]->Fill(corr_val[pol],snr_val[pol],weight);
	
						for(int bin=0; bin<numSNRbins; bin++){ //check all potential intercepts
							double this_y_val = (slope * corr_val[pol] ) + intercept[bin]; // compute the SNR to pass at this intercept
							// printf("For bin %d, with intercept %.2f, SNR to pass is %.2f \n", bin, intercept[bin], this_y_val);
							if(snr_val[pol]>=this_y_val){ // does this event pass?
								// printf("     This event has SNR %.2f, so it passes!\n",snr_val[pol]);
								// printf("          Old number of events in this bin %5f \n", numEventsPassed[bin]);
								numEventsPassed[pol][bin]+=1.;
								// printf(".         New number of events in this bin %5f \n", numEventsPassed[bin]);
							}
						}
					}
				}
			}
		} // loop over events
		inputFile->Close();
		delete inputFile;
	} // loop over files

	// okay, now save out the 2D histogram
	TCanvas *cSNRvsCorr = new TCanvas("","",2.1*850, 850);
	cSNRvsCorr->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		cSNRvsCorr->cd(pol+1);
		h2SNRvsCorr[pol]->Draw("colz");
		h2SNRvsCorr[pol]->GetYaxis()->SetTitle("3rd Highest VPeak/RMS");
		h2SNRvsCorr[pol]->GetXaxis()->SetTitle("Peak Corr");
		gPad->SetLogz();
	}
	char title[300];
	sprintf(title, "%s/optimize/A%d_config%d_%dEvents_SNRvsCorr.png",plotPath,station,config,numTotal);
	cSNRvsCorr->SaveAs(title);
	delete cSNRvsCorr;

	// for(int bin=0; bin<numSNRbins; bin++){
	// 	printf("Bin %d at cut %.2f has %5f events passing \n",bin,intercept[bin], numEventsPassed[bin]);
	// }

	// now to get the differential distribution up and running
	for(int pol=0; pol<2; pol++){
		for(int bin=0; bin<numSNRbins-1; bin++){
			numEventsPassed_diff[pol][bin] = numEventsPassed[pol][bin] - numEventsPassed[pol][bin+1];
			// printf("Pol %d Bin %d at cut %.2f has %5.f events passing, and next bin has %5f, so diff is %.5f \n", pol, bin, intercept[pol][bin],numEventsPassed[pol][bin],numEventsPassed[pol][bin+1],numEventsPassed_diff[pol][bin]);
		}
	}

	TH1D *hEventsVsSNR[2];
	TH1D *hEventsVsSNR[0] = new TH1D("DiffDistroV","DiffDistroV",numSNRbins,0,20);
	TH1D *hEventsVsSNR[1] = new TH1D("DiffDistroH","DiffDistroH",numSNRbins,0,20);

	int max_bin[2];
	double fit_start_bin[2];
	double start_of_fit[2];
	double end_of_fit[2];
	int last_filled_bin[2];

	for(int pol=0; pol<2; pol++){
		for(int bin=0; bin<numSNRbins; bin++){
			hEventsVsSNR[pol]->SetBinContent(bin+1,numEventsPassed_diff[pol][bin]);
		}
		max_bin[pol] = hEventsVsSNR[pol]->GetMaximumBin();
		fit_start_bin[pol] = max_bin[pol]+14;
		start_of_fit[pol] = hEventsVsSNR[pol]->GetBinCenter(fit_start_bin[pol]);
		last_filled_bin[pol] = hEventsVsSNR[pol]->FindLastBinAbove(0.,1);
		last_filled_bin[pol]+=5; // go up some more bins just to make sure the fit is over
		end_of_fit[pol] = hEventsVsSNR[pol]->GetBinCenter(last_filled_bin[pol]);
		// printf("Pol %d Start of fit is %.2f and end of fit is %.2f \n", pol, start_of_fit[pol], end_of_fit[pol]);
	}


	// now we exponential fit
	char equation[150];
	sprintf(equation,"exp([0]*x+[1])");
	char equation_name[2][150];
	TFit *fit[2];
	int status[2];
	double fitParams[2][2];
	for(int pol=0; pol<2; pol++){
		sprintf(equation_name,"ExpFit%d",pol);
		fit[pol] = new TF1(equation_name[pol],equation,start_of_fit[pol],end_of_fit[pol]);
		status[pol] = hEventsVsSNR[pol]->Fit(equation_name,"LL,R");
		fitParams[pol][0] = fit[pol]->GetParameter(0);
		fitParams[pol][1] = fit[pol]->GetParameter(1);
		printf("Pol %d Fit Parameters are %.2f and %.2f \n", fitParams[pol][0], fitParams[pol][1]);
	}

	double binWidthIntercept = hEventsVsSNR->GetBinWidth(1);
	double leftBoundary = start_of_fit - binWidthIntercept/2.;
	double rightBoundary = end_of_fit + binWidthIntercept/2.;
	int numBinsThisFit = (rightBoundary - leftBoundary)/binWidthIntercept + 1; // how many bins do we need in our histogram to actually do the fitting
	printf("Number of bins in this fit is %d \n", numBinsThisFit);
	char this_fit_title[400];
	sprintf(this_fit_title,"Fit_slope_%.2f ", slope);
	TH1D *hNumObserved = new TH1D(this_fit_title,"",numBinsThisFit,leftBoundary,rightBoundary);
	for(int bin=0; bin<numBinsThisFit; bin++){
		double originalContent = hEventsVsSNR->GetBinContent(bin+fit_start_bin);
		hNumObserved->SetBinContent(bin+1,originalContent);
	}

	/*
		Now we must prepare our *expectation* for the number of events
		in a bin given that we now have our exponential model
	*/

	double numExpected[numBinsThisFit];
	for(int bin=0; bin<numBinsThisFit; bin++){
		double modelPrediction = exp(fitParams[0]*(hNumObserved->GetBinCenter(bin+1)) + fitParams[1]);
		numExpected[bin] = modelPrediction;
	}

	/*
		Now we compute the log-likelihood by hand
	*/

	double logL=0.;
	double numObservedTotal=0.;
	double numExpectedTotal_sum=0;
	for(int bin=0; bin<numBinsThisFit; bin++){
		double thisObserved = hNumObserved->GetBinContent(bin+1);
		double thisExpected = numExpected[bin];
		printf("At bin %d Observed %.2f and Expected %.2f \n", bin, thisObserved, thisExpected );
		numObservedTotal+=thisObserved;
		numExpectedTotal_sum+=thisExpected;
		logL += ReturnLogL_highN( thisObserved,thisExpected );
	}
	printf("The logL is %.3f \n", logL);
	
	// double numExpectedTotal_Integral = (1./(fitParams[0])) * ( exp(fitParams[0]*end_of_fit + fitParams[1]) - exp(fitParams[0]*start_of_fit + fitParams[1]));
	// double numExpectedTotal_Integral = (1./(fitParams[0]*binWidthIntercept)) * ( exp(fitParams[0]*end_of_fit + fitParams[1]) - exp(fitParams[0]*start_of_fit + fitParams[1]));
	double numExpectedTotal_Integral = (1./(fitParams[0]*binWidthIntercept)) * ( exp(fitParams[0]*rightBoundary + fitParams[1]) - exp(fitParams[0]*leftBoundary + fitParams[1]));
	printf("Best fit sum bins: %.2f. Best fit do integral: %.2f. Num observed: %.2f. \n",numExpectedTotal_sum, numExpectedTotal_Integral, numObservedTotal);

	TRandom3 *test_random = new TRandom3();
	cout<<"Poisson number of random events is "<<test_random->Poisson(numExpectedTotal_Integral)<<endl;

	/*
		Now for toy simulations
	*/
	sprintf(this_fit_title,"fCopyFit", slope);
	TF1 *fitCopy = new TF1(this_fit_title, "exp([0]*x+[1])", start_of_fit, end_of_fit);
	fitCopy->SetParameters(fitParams[0], fitParams[1]);
	double less_BestFit_logL = 0.; // values lower than BestFit_logL
	double Total_Toy_logL = 0.; // total logL values from Toy

	int num_Toy = 10000;
	int Toy_logL_bin = 100;
	double min_Toy_logL = 0.;
	double max_Toy_logL = 400.;
	TH1D *hToy_logL = new TH1D("hToy_logL", "", Toy_logL_bin, min_Toy_logL, max_Toy_logL );
	char test_title[400];
	for(int num=0; num<num_Toy; num++){
		// fill this toy pseudo experiment with a poisson fluctuations of the events observed
		sprintf(test_title, "hToy %d ", num);
		TH1D *hToy = new TH1D(test_title,"",numBinsThisFit,leftBoundary,rightBoundary);
		// int Evts_Poisson = test_random->Poisson(numExpectedTotal_sum);
		int Evts_Poisson = test_random->Poisson(numExpectedTotal_Integral);
		hToy->FillRandom("fCopyFit",Evts_Poisson);
		double logL_log_Toy=0.; // get logL for this toy
		for(int bin=0; bin<numBinsThisFit; bin++){
			logL_log_Toy+=ReturnLogL_highN(hToy->GetBinContent(bin+1), numExpected[bin]);
		}
		// TCanvas *cToyHist = new TCanvas ("cToyHist","", 800, 600);
		// cToyHist->cd();
		// 	cToyHist->cd()->SetLogy();
		// 	sprintf( test_title, "Toy hist, evts : %d, -2Ln(L): %.2f", Evts_Poisson, logL_log_Toy );
		// 	hToy->SetTitle(test_title);
		// 	hToy->Draw();
		// 	fit->Draw("same");
		// char this_save_title[400];
		// sprintf(this_save_title,"toy%d.png",num);
		// cToyHist->SaveAs(this_save_title);
		// delete cToyHist;
		delete hToy;
		hToy_logL->Fill(logL_log_Toy);
		if ( logL_log_Toy <= logL )
			less_BestFit_logL += logL_log_Toy;
		Total_Toy_logL += logL_log_Toy;
	}

	// vertical line for log likelihood
	double P_value = less_BestFit_logL / Total_Toy_logL;
	double DataLogL_x[2] = { logL, logL };
	double DataLogL_y[2] = { 0, 5 };
	TGraph *gDataLogL = new TGraph (2, DataLogL_x, DataLogL_y);

	printf(BLUE"About to do background estimation\n"RESET);
	/*
		Now, we must find our background estimate
		Which is where we integrate the exponential above our cut value
		And then use that to get s_up
	*/

	double S_up_array[numSNRbins];
	double S_array[numSNRbins];
	for(int bin=0; bin<numSNRbins; bin++){
		S_up_array[bin]=0.;
		S_array[bin]=0.;
	}
	int startBin = 80;
	for(int bin=startBin; bin<numSNRbins; bin++){
		double cut = intercept[bin];
		double back_estimate = (1./(fitParams[0]*binWidthIntercept)) * (-exp(fitParams[0]*cut + fitParams[1]));
		back_estimate*=10.; // make it 10 times bigger, for switch to 100% sample
		double achieved_alpha;
		// double s_up = GetS_up_TMath(back_estimate,achieved_alpha, 0.9); // compute S_up for this background
		double s_up = GetS_up(back_estimate,achieved_alpha, 0.9); // compute S_up for this background
		S_up_array[bin] = s_up;
		printf("For cut %.2f background estimate is %.3f with sup %.2f \n",cut,back_estimate,S_up_array[bin]);
	}

	// for(double cut=8.; cut<14; cut+=0.1){
	// 	double back_estimate = (1./(fitParams[0]*binWidthIntercept)) * (-exp(fitParams[0]*cut + fitParams[1]));
	// 	back_estimate*=10; // make it 10 times bigger, for switch to 100% sample
	// 	double achieved_alpha;
	// 	double s_up = GetS_up_TMath(back_estimate,achieved_alpha, 0.9); // compute 90% CL UL for this many background
	// 	printf("For cut %.2f background estimate is %.3f with sup %.2f \n",cut,back_estimate,s_up);
	// }

	double select_slope=slope;

	/*
		Now we must loop over data and simulation
	*/

	TChain simVTree("VTree");
	TChain simHTree("HTree");
	TChain simAllTree("AllTree");
	char the_sims[500];
	sprintf(the_sims,"/fs/scratch/PAS0654/ara/sim/ValsForCuts/A%d/c%d/E%d/cutvals_drop_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_*.root",station,config,int(year_or_energy));
	simVTree.Add(the_sims);
	simHTree.Add(the_sims);
	simAllTree.Add(the_sims);
	int numSimEvents = simVTree.GetEntries();
	printf("Num of sim entries is %d \n", numSimEvents);

	TH2D *h2SNRvsCorr_sim[2]; // SNR on Y axis, Corr on X axis, like in the TB
	h2SNRvsCorr_sim[0]=new TH2D("","V Sim",100,0,max,30,0,30);
	h2SNRvsCorr_sim[1]=new TH2D("","H Sim",100,0,max,30,0,30);

	// and now get values out
	{
		double corr_val[2];
		double snr_val[2];
		int WFRMS[2];
		double frac_of_power_notched_V[8];
		double frac_of_power_notched_H[8];
		int Refilt[2];

		simVTree.SetBranchAddress("corr_val_V",&corr_val[0]);
		simVTree.SetBranchAddress("snr_val_V",&snr_val[0]);
		simVTree.SetBranchAddress("wfrms_val_V",&WFRMS[0]);
		simVTree.SetBranchAddress("Refilt_V",&Refilt[0]);
		simHTree.SetBranchAddress("corr_val_H",&corr_val[1]);
		simHTree.SetBranchAddress("snr_val_H",&snr_val[1]);
		simHTree.SetBranchAddress("wfrms_val_H",&WFRMS[1]);
		simHTree.SetBranchAddress("Refilt_H",&Refilt[1]);

		int isCal;
		int isSoft;
		int isShort;
		int isCW;
		int isNewBox;
		int isSurf[2];
		int isBadEvent;
		double weight;
		int isSurfEvent_top[2];
		int unixTime;
		int isFirstFiveEvent;
		int hasBadSpareChanIssue;

		simAllTree.SetBranchAddress("cal",&isCal);
		simAllTree.SetBranchAddress("soft",&isSoft);
		simAllTree.SetBranchAddress("short",&isShort);
		simAllTree.SetBranchAddress("CW",&isCW);
		simAllTree.SetBranchAddress("box",&isNewBox);
		simAllTree.SetBranchAddress("surf_V",&isSurf[0]);
		simAllTree.SetBranchAddress("surf_H",&isSurf[1]);
		simAllTree.SetBranchAddress("bad",&isBadEvent);
		simAllTree.SetBranchAddress("weight",&weight);
		simAllTree.SetBranchAddress("surf_top_V",&isSurfEvent_top[0]);
		simAllTree.SetBranchAddress("surf_top_H",&isSurfEvent_top[1]);
		simAllTree.SetBranchAddress("unixTime",&unixTime);
		simAllTree.SetBranchAddress("isFirstFiveEvent",&isFirstFiveEvent);
		simAllTree.SetBranchAddress("hasBadSpareChanIssue",&hasBadSpareChanIssue);

		stringstream ss;
		for(int i=0; i<8; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			simVTree.SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_V[i]);
		}
		for(int i=8; i<16; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			simHTree.SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_H[i-8]);
		}

		for(int event=0; event<numSimEvents; event++){

			simVTree.GetEvent(event);
			simHTree.GetEvent(event);
			simAllTree.GetEvent(event);
			for(int pol=0; pol<2; pol++){

				h2SNRvsCorr_sim[pol]->Fill(corr_val[pol], snr_val[pol],weight);
				if(pol==1){
					// only vpol for right now
					continue;
				}

				bool failsCWPowerCut=false;
				if(Refilt[pol] && !WFRMS[pol]){
					vector<double> frac;
					for(int i=0; i<8; i++){
						if(pol==0) frac.push_back(frac_of_power_notched_V[i]);
						else if(pol==1) frac.push_back(frac_of_power_notched_H[i]);
					}
					sort(frac.begin(), frac.end(), std::greater<double>());
					if(frac[2]>0.06){
						failsCWPowerCut=true;
					}
				} //refiltered?

				double this_SNR=snr_val[pol];
				double this_corr=corr_val[pol];

				if(!WFRMS[pol] && !failsCWPowerCut){
					if(!isNewBox){
						if(!isSurf[0] && !isSurf[1] && !isSurfEvent_top[pol]){
							// loop over every bin (intercept value), and figure out if this event would have passed or not
							for(int bin=startBin; bin<numSNRbins; bin++){
								double failsRcut=false;
								double thisIntercept = intercept[bin];
								double this_y_val = this_corr*select_slope + thisIntercept;
								if(this_SNR>=this_y_val){
									S_array[bin]+=weight;
								}
							}
						}
					}
				}
			} //loop over polarizations
		}
	}
	double SoverSup[numSNRbins];
	for(int bin=0; bin<numSNRbins; bin++){
		double this_S = S_array[bin];
		double this_Sup = S_up_array[bin];
		double this_SoverSup;
		if(!this_Sup>0)
			this_SoverSup=0.;
		else
			this_SoverSup = this_S/this_Sup;

		SoverSup[bin] = this_SoverSup;
		printf("For bin %d, intercept %.2f, S is %.2f and S_up is %.2f for S/S_up of %.2f  \n", bin,intercept[bin],this_S, this_Sup, this_SoverSup);
	}
	double max_SoverSup=0.;
	double optimal_intercept=0.;
	for(int bin=0; bin<numSNRbins; bin++){
		double this_intercept = intercept[bin];
		if(this_intercept<8. || this_intercept>13.){
			continue;
		}
		else{
			double this_SoverSup=SoverSup[bin];
			if(this_SoverSup>max_SoverSup){
				printf("At bin %d, for intercept %.2f, new S/Sup of %.2f is greater than current %.2f \n",bin, intercept[bin], this_SoverSup,max_SoverSup);
				optimal_intercept = this_intercept;
				max_SoverSup=this_SoverSup;
			}
		}
	}
	TGraph *gSoverSup = new TGraph(numSNRbins,intercept,SoverSup);

	printf("Found optimal intercept at %.2f \n", optimal_intercept);
	double select_inter=optimal_intercept;

	vector <double> x_vals_for_line;
	vector <double> y_vals_for_line;
	for(double x=0; x<0.020; x+=0.00001){
		double y_val = (slope * x ) + optimal_intercept;
		x_vals_for_line.push_back(x);
		y_vals_for_line.push_back(y_val);
	}
	TGraph *cut_line = new TGraph(x_vals_for_line.size(), &x_vals_for_line[0], &y_vals_for_line[0]);


	TCanvas *cRcut = new TCanvas("","",4*850,2*850);
	cRcut->Divide(4,2);
	cRcut->cd(1);
		h2SNRvsCorr[0]->Draw("colz");
		h2SNRvsCorr[0]->GetXaxis()->SetTitle("Correlation Value");
		h2SNRvsCorr[0]->GetYaxis()->SetTitle("SNR");
		gPad->SetLogz();
		cut_line->Draw("same");
		cut_line->SetLineColor(kRed);
	cRcut->cd(2);
		hEventsVsSNR->Draw("");
		hEventsVsSNR->GetXaxis()->SetTitle("SNR Cut (y-intercept value)");
		hEventsVsSNR->GetYaxis()->SetTitle("Number of Events Cut");
		gPad->SetLogy();
		// hEventsVsSNR->GetXaxis()->SetRangeUser(8.6,10.);
		// hEventsVsSNR->GetYaxis()->SetRangeUser(8e1,1e4);
	cRcut->cd(3);
		hNumObserved->Draw("HIST");
		hNumObserved->GetXaxis()->SetTitle("SNR Cut (y-intercept value)");
		hNumObserved->GetYaxis()->SetTitle("Number of Events Cut");
		// hNumObserved->GetXaxis()->SetRangeUser(8.6,10.);
		// hNumObserved->GetYaxis()->SetRangeUser(8e1,1e4);
		fit->Draw("same");
		gPad->SetLogy();
		// hNumObserved->GetYaxis()->SetRangeUser(0.2,2e3);
	cRcut->cd(4);
		sprintf( test_title, "data logL: %.2f, P-value: %f", logL, P_value );
		hToy_logL->SetTitle(test_title);
		hToy_logL->Draw();
		hToy_logL->GetYaxis()->SetTitle("Number of Pseudo Experiments");
		hToy_logL->GetXaxis()->SetTitle("-2log(L)");
		gDataLogL->SetLineColor(kRed);
		gDataLogL->Draw("l");
		gPad->SetLogy();
	cRcut->cd(5);
		h2SNRvsCorr_sim[0]->Draw("colz");
		h2SNRvsCorr_sim[0]->GetXaxis()->SetTitle("Correlation Value");
		h2SNRvsCorr_sim[0]->GetYaxis()->SetTitle("SNR");
		gPad->SetLogz();
		cut_line->Draw("same");
		cut_line->SetLineColor(kRed);
	cRcut->cd(6);
		gSoverSup->Draw("ALP");
		gSoverSup->GetXaxis()->SetRangeUser(8.,13.);
		gSoverSup->GetYaxis()->SetRangeUser(0.,5e3);
		gSoverSup->GetYaxis()->SetTitle("S/S_up");
		gSoverSup->GetYaxis()->SetTitleOffset(1.8);
		gSoverSup->GetXaxis()->SetTitle("SNR Cut (y-intercept value)");
	char save_title[400];
	sprintf(save_title,"%s/optimize/A%d_config%d_DiffDistro_RCutSlope%.4f.png",plotPath,station,config,slope);
	cRcut->SaveAs(save_title);

	printf(BLUE"Now to loop back through the data and compute cut tables\n"RESET);

	/*
		We must now also loop through the data one more time to compute the "as last cut"
		VPol for now, just to get a feeling for what the filter is doing.
	*/
	double num_total_data=0.;
	double remove_bad_runs_and_livetime_data=0.;
	double and_remove_soft_data=0.;
	double and_remove_short_and_glitch_data=0.;
	double and_remove_tagged_cal_data=0.;

	double fails_WFRMS_first_data[2]={0.,0.};
	double fails_box_first_data[2]={0.,0.};
	double fails_surface_first_data[2]={0.,0.};
	double fails_rcut_first_data[2]={0.,0.};

	double fails_WFRMS_last_data[2]={0.,0.};
	double fails_box_last_data[2]={0.,0.};
	double fails_surface_last_data[2]={0.,0.};
	double fails_rcut_last_data[2]={0.,0.};

	double fails_WFRMS_insequence_data[2]={0.,0.};
	double fails_box_insequence_data[2]={0.,0.};
	double fails_surface_insequence_data[2]={0.,0.};
	double fails_rcut_insequence_data[2]={0.,0.};

	// now, we go back through the data a second time to get "as last cut" tables
	for(int file_num=6; file_num<argc; file_num++){

		string chRun = "run";
		string file = string(argv[file_num]);
		size_t foundRun = file.find(chRun);
		string strRunNum = file.substr(foundRun+4,4);
		int runNum = atoi(strRunNum.c_str());

		TFile *inputFile = TFile::Open(argv[file_num]);
		if(!inputFile){
			cout<<"Can't open val for cuts file!"<<endl;
			return -1;
		}
		printf("As last cut loop data file %d: run %d \n", file_num, runNum);

		TTree *trees[3];
		trees[0] = (TTree*) inputFile->Get("VTree");
		trees[1] = (TTree*) inputFile->Get("HTree");
		trees[2] = (TTree*) inputFile->Get("AllTree");

		double corr_val[2];
		double snr_val[2];
		int WFRMS[2];
		double frac_of_power_notched_V[8];
		double frac_of_power_notched_H[8];
		int Refilt[2];

		trees[0]->SetBranchAddress("corr_val_V",&corr_val[0]);
		trees[0]->SetBranchAddress("snr_val_V",&snr_val[0]);
		trees[0]->SetBranchAddress("wfrms_val_V",&WFRMS[0]);
		trees[0]->SetBranchAddress("Refilt_V",&Refilt[0]);
		trees[1]->SetBranchAddress("corr_val_H",&corr_val[1]);
		trees[1]->SetBranchAddress("snr_val_H",&snr_val[1]);
		trees[1]->SetBranchAddress("wfrms_val_H",&WFRMS[1]);
		trees[1]->SetBranchAddress("Refilt_H",&Refilt[1]);

		int isCal;
		int isSoft;
		int isShort;
		int isCW;
		int isNewBox;
		int isSurf[2];
		int isBadEvent;
		double weight;
		int isSurfEvent_top[2];
		int unixTime;
		int isFirstFiveEvent;
		int hasBadSpareChanIssue;

		trees[2]->SetBranchAddress("cal",&isCal);
		trees[2]->SetBranchAddress("soft",&isSoft);
		trees[2]->SetBranchAddress("short",&isShort);
		trees[2]->SetBranchAddress("CW",&isCW);
		trees[2]->SetBranchAddress("box",&isNewBox);
		trees[2]->SetBranchAddress("surf_V",&isSurf[0]);
		trees[2]->SetBranchAddress("surf_H",&isSurf[1]);
		trees[2]->SetBranchAddress("bad",&isBadEvent);
		trees[2]->SetBranchAddress("weight",&weight);
		trees[2]->SetBranchAddress("surf_top_V",&isSurfEvent_top[0]);
		trees[2]->SetBranchAddress("surf_top_H",&isSurfEvent_top[1]);
		trees[2]->SetBranchAddress("unixTime",&unixTime);
		trees[2]->SetBranchAddress("isFirstFiveEvent",&isFirstFiveEvent);
		trees[2]->SetBranchAddress("hasBadSpareChanIssue",&hasBadSpareChanIssue);

		stringstream ss;
		for(int i=0; i<8; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			trees[0]->SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_V[i]);
		}
		for(int i=8; i<16; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			trees[1]->SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_H[i-8]);
		}
		
		int numEntries = trees[0]->GetEntries();
		numTotal+=numEntries;

		if(isBadRun(station,runNum,BadRunList)){
			continue;
		}

		//now to loop over events
		for(int event=0; event<numEntries; event++){

			trees[0]->GetEvent(event);
			trees[1]->GetEvent(event);
			trees[2]->GetEvent(event);

			num_total_data+=weight;


			if(isBadLivetime(station,unixTime) && !isSim){
				continue;
			}
			remove_bad_runs_and_livetime_data+=weight;

			if(isSoft){
				continue;
			}
			and_remove_soft_data+=weight;

			if(
				(isBadEvent
				|| isFirstFiveEvent
				|| hasBadSpareChanIssue
				|| isShort)
				&& !isSim
				){
				continue;
			}
			and_remove_short_and_glitch_data+=weight;

			if(isCal){
				continue;
			}
			and_remove_tagged_cal_data+=weight;
				
			for(int pol=0; pol<2; pol++){
				bool failsCWPowerCut=false;
				if(Refilt[pol] && !WFRMS[pol]){
					vector<double> frac;
					for(int i=0; i<8; i++){
						if(pol==0) frac.push_back(frac_of_power_notched_V[i]);
						else if(pol==1) frac.push_back(frac_of_power_notched_H[i]);
					}
					sort(frac.begin(), frac.end(), std::greater<double>());
					if(frac[2]>0.06){
						failsCWPowerCut=true;
					}
				} //refiltered?

				bool failsRcut=false;
				double this_SNR=snr_val[pol];
				double this_corr=corr_val[pol];
				double this_y_val = this_corr*select_slope + select_inter;
				if(this_SNR < this_y_val){
					failsRcut=true;
				}

				// "as first cut"
					// fail WFRMS first?
					if(WFRMS[pol] || failsCWPowerCut){
						fails_WFRMS_first_data[pol]+=weight;
					}
					// fail box first?
					if(isNewBox){
						fails_box_first_data[pol]+=weight;
					}
					// fail surface first?
					if(isSurf[0] || isSurf[1] || isSurfEvent_top[pol]){
						fails_surface_first_data[pol]+=weight;
					}
					if(failsRcut){
						fails_rcut_first_data[pol]+=weight;
					}


				// "as last cut"
					// fails as last cut with surface?
					// survives WFRMS and box and Rcut, but doesn't survive surface
					if(!WFRMS[pol] && !failsCWPowerCut && !isNewBox && !failsRcut){
						if(isSurf[0] || isSurf[1] || isSurfEvent_top[pol]){
							fails_surface_last_data[pol]+=weight;
						}
					}
					// fails as last cut with WFRMS?
					// survives box and surface, but doesn't survive WFRMS
					if(!isNewBox && !isSurf[0] && !isSurf[1] && !isSurfEvent_top[pol] && !failsRcut){
						if(WFRMS[pol] || failsCWPowerCut){
							fails_WFRMS_last_data[pol]+=weight;
						}
					}
					// fails as last cut with box?
					// survives surface and WFRMS, but not the box
					if(!isSurf[0] && !isSurf[1] && !isSurfEvent_top[pol] && !WFRMS[pol] && !failsCWPowerCut && !failsRcut){
						if(isNewBox){
							fails_box_last_data[pol]+=weight;
						}
					}
					// fails as last cust with Rcut?
					if(!isSurf[0] && !isSurf[1] && !isSurfEvent_top[pol] && !WFRMS[pol] && !failsCWPowerCut && !isNewBox){
						if(failsRcut){
							fails_rcut_last_data[pol]+=weight;
						}
					}

				// "in sequence"
					// fails WFRMS first? (same as "as first" for this cut only)
					if(WFRMS[pol] || failsCWPowerCut){
						fails_WFRMS_insequence_data[pol]+=weight;
					}
					// passes WFRMS, but fails box
					if(!WFRMS[pol] && !failsCWPowerCut && isNewBox){
						fails_box_insequence_data[pol]+=weight;
					}
					// passes WFRMS and box, but fails surface
					if(!WFRMS[pol] && !failsCWPowerCut && !isNewBox && (isSurf[0] || isSurf[1] || isSurfEvent_top[pol])){
						fails_surface_insequence_data[pol]+=weight;
					}
					// passes WFRMS, box, and surface, but fails Rcut (same as "as last" for this cut only)
					if(!WFRMS[pol] && !failsCWPowerCut && !isNewBox && (!isSurf[0] && !isSurf[1] && !isSurfEvent_top[pol]) && failsRcut){
						fails_rcut_insequence_data[pol]+=weight;
					}
			}
		} // loop over events
		inputFile->Close();
		delete inputFile;
	} // loop over files

	printf("Num total          :           %7.1f\n",num_total_data);
	printf("Livetime           :           %7.1f\n",remove_bad_runs_and_livetime_data);
	printf("Soft Trig          :           %7.1f\n",and_remove_soft_data);
	printf("Short and Glitch   :           %7.1f\n",and_remove_short_and_glitch_data);
	printf("Tagged Cal         :           %7.1f\n",and_remove_tagged_cal_data);
	printf("------------------------------------------\n");
	printf("WFRMS              :           %7.1f, %7.1f, %7.1f | %7.1f, %7.1f, %7.1f \n",fails_WFRMS_first_data[0],fails_WFRMS_insequence_data[0],fails_WFRMS_last_data[0],fails_WFRMS_first_data[1],fails_WFRMS_insequence_data[1],fails_WFRMS_last_data[1]);
	printf("Box                :           %7.1f, %7.1f, %7.1f | %7.1f, %7.1f, %7.1f \n",fails_box_first_data[0],fails_box_insequence_data[0],fails_box_last_data[0],fails_box_first_data[1],fails_box_insequence_data[1],fails_box_last_data[1]);
	printf("Surf               :           %7.1f, %7.1f, %7.1f | %7.1f, %7.1f, %7.1f \n",fails_surface_first_data[0],fails_surface_insequence_data[0],fails_surface_last_data[0],fails_surface_first_data[1],fails_surface_insequence_data[1],fails_surface_last_data[1]);
	printf("Rcut               :           %7.1f, %7.1f, %7.1f | %7.1f, %7.1f, %7.1f \n",fails_rcut_first_data[0],fails_rcut_insequence_data[0],fails_rcut_last_data[0],fails_rcut_first_data[1],fails_rcut_insequence_data[1],fails_rcut_last_data[1]);

	printf(RED"Now to try and go over the simulation again\n"RESET);
	/*
		And now we go through the simulation again to get our as last cut table and to compute efficiencies
	*/
	double num_total_sim=0.;

	double fails_WFRMS_first_sim[2]={0.,0.};
	double fails_box_first_sim[2]={0.,0.};
	double fails_surface_first_sim[2]={0.,0.};
	double fails_rcut_first_sim[2]={0.,0.};

	double fails_WFRMS_last_sim[2]={0.,0.};
	double fails_box_last_sim[2]={0.,0.};
	double fails_surface_last_sim[2]={0.,0.};
	double fails_rcut_last_sim[2]={0.,0.};

	double fails_WFRMS_insequence_sim[2]={0.,0.};
	double fails_box_insequence_sim[2]={0.,0.};
	double fails_surface_insequence_sim[2]={0.,0.};
	double fails_rcut_insequence_sim[2]={0.,0.};

	TH1D *all_events[2];
	TH1D *pass_soft_short_cal_wfrms[2];
	TH1D *pass_soft_short_cal_wfrms_box[2];
	TH1D *pass_soft_short_cal_wfrms_box_surf[2];
	TH1D *pass_soft_short_cal_wfrms_box_surf_rcut[2];

	TH1D *eff[2];
	TH1D *eff_soft_short_cal[2];
	TH1D *eff_soft_short_cal_wfrms[2];
	TH1D *eff_soft_short_cal_wfrms_box[2];
	TH1D *eff_soft_short_cal_wfrms_box_surf[2];
	TH1D *eff_soft_short_cal_wfrms_box_surf_rcut[2];

	for(int i=0; i<2; i++){
		all_events[i] = new TH1D("","",30,0,30);
		pass_soft_short_cal_wfrms[i] = new TH1D("","",30,0,30);
		pass_soft_short_cal_wfrms_box[i] = new TH1D("","",30,0,30);
		pass_soft_short_cal_wfrms_box_surf[i] = new TH1D("","",30,0,30);
		pass_soft_short_cal_wfrms_box_surf_rcut[i] = new TH1D("","",30,0,30);

		eff_soft_short_cal_wfrms[i] = new TH1D("","",30,0,30);
		eff_soft_short_cal_wfrms_box[i] = new TH1D("","",30,0,30);
		eff_soft_short_cal_wfrms_box_surf[i] = new TH1D("","",30,0,30);
		eff_soft_short_cal_wfrms_box_surf_rcut[i] = new TH1D("","",30,0,30);
	}


	{
		double corr_val[2];
		double snr_val[2];
		int WFRMS[2];
		double frac_of_power_notched_V[8];
		double frac_of_power_notched_H[8];
		int Refilt[2];

		simVTree.SetBranchAddress("corr_val_V",&corr_val[0]);
		simVTree.SetBranchAddress("snr_val_V",&snr_val[0]);
		simVTree.SetBranchAddress("wfrms_val_V",&WFRMS[0]);
		simVTree.SetBranchAddress("Refilt_V",&Refilt[0]);
		simHTree.SetBranchAddress("corr_val_H",&corr_val[1]);
		simHTree.SetBranchAddress("snr_val_H",&snr_val[1]);
		simHTree.SetBranchAddress("wfrms_val_H",&WFRMS[1]);
		simHTree.SetBranchAddress("Refilt_H",&Refilt[1]);

		int isCal;
		int isSoft;
		int isShort;
		int isCW;
		int isNewBox;
		int isSurf[2];
		int isBadEvent;
		double weight;
		int isSurfEvent_top[2];
		int unixTime;
		int isFirstFiveEvent;
		int hasBadSpareChanIssue;

		simAllTree.SetBranchAddress("cal",&isCal);
		simAllTree.SetBranchAddress("soft",&isSoft);
		simAllTree.SetBranchAddress("short",&isShort);
		simAllTree.SetBranchAddress("CW",&isCW);
		simAllTree.SetBranchAddress("box",&isNewBox);
		simAllTree.SetBranchAddress("surf_V",&isSurf[0]);
		simAllTree.SetBranchAddress("surf_H",&isSurf[1]);
		simAllTree.SetBranchAddress("bad",&isBadEvent);
		simAllTree.SetBranchAddress("weight",&weight);
		simAllTree.SetBranchAddress("surf_top_V",&isSurfEvent_top[0]);
		simAllTree.SetBranchAddress("surf_top_H",&isSurfEvent_top[1]);
		simAllTree.SetBranchAddress("unixTime",&unixTime);
		simAllTree.SetBranchAddress("isFirstFiveEvent",&isFirstFiveEvent);
		simAllTree.SetBranchAddress("hasBadSpareChanIssue",&hasBadSpareChanIssue);

		stringstream ss;
		for(int i=0; i<8; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			simVTree.SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_V[i]);
		}
		for(int i=8; i<16; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			simHTree.SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_H[i-8]);
		}

		for(int event=0; event<numSimEvents; event++){

			simVTree.GetEvent(event);
			simHTree.GetEvent(event);
			simAllTree.GetEvent(event);

			num_total_sim+=weight;

			for(int pol=0; pol<2; pol++){
				bool failsCWPowerCut=false;
				if(Refilt[pol] && !WFRMS[pol]){
					vector<double> frac;
					for(int i=0; i<8; i++){
						if(pol==0) frac.push_back(frac_of_power_notched_V[i]);
						else if(pol==1) frac.push_back(frac_of_power_notched_H[i]);
					}
					sort(frac.begin(), frac.end(), std::greater<double>());
					if(frac[2]>0.06){
						failsCWPowerCut=true;
					}
				} //refiltered?

				bool failsRcut=false;
				double this_SNR=snr_val[pol];
				double this_corr=corr_val[pol];
				double this_y_val = this_corr*select_slope + select_inter;
				if(this_SNR < this_y_val){
					failsRcut=true;
				}

				if (this_SNR>30.) this_SNR=30.;
				all_events[pol]->Fill(this_SNR,weight);
				if(!WFRMS[pol] && !failsCWPowerCut){
					pass_soft_short_cal_wfrms[pol]->Fill(this_SNR,weight);
					if(!isNewBox){
						pass_soft_short_cal_wfrms_box[pol]->Fill(this_SNR,weight);
						if(!isSurf[0] && !isSurf[1] && !isSurfEvent_top[pol]){
							pass_soft_short_cal_wfrms_box_surf[pol]->Fill(this_SNR,weight);
							if(!failsRcut)
								pass_soft_short_cal_wfrms_box_surf_rcut[pol]->Fill(this_SNR,weight);
						}
					}
				}

				// "as first cut"
					// fail WFRMS first?
					if(WFRMS[pol] || failsCWPowerCut){
						fails_WFRMS_first_sim[pol]+=weight;
					}
					// fail box first?
					if(isNewBox){
						fails_box_first_sim[pol]+=weight;
					}
					// fail surface first?
					if(isSurf[0] || isSurf[1] || isSurfEvent_top[pol]){
						fails_surface_first_sim[pol]+=weight;
					}
					if(failsRcut){
						fails_rcut_first_sim[pol]+=weight;
					}


				// "as last cut"
					// fails as last cut with surface?
					// survives WFRMS and box, but doesn't survive surface
					if(!WFRMS[pol] && !failsCWPowerCut && !isNewBox && !failsRcut){
						if(isSurf[0] || isSurf[1] || isSurfEvent_top[pol]){
							fails_surface_last_sim[pol]+=weight;
						}
					}
					// fails as last cut with WFRMS?
					// survives box and surface, but doesn't survive WFRMS
					if(!isNewBox && !isSurf[0] && !isSurf[1] && !isSurfEvent_top[pol] && !failsRcut){
						if(WFRMS[pol] || failsCWPowerCut){
							fails_WFRMS_last_sim[pol]+=weight;
						}
					}
					// fails as last cut with box?
					// survives surface and WFRMS, but not the box
					if(!isSurf[0] && !isSurf[1] && !isSurfEvent_top[pol] && !WFRMS[pol] && !failsCWPowerCut && !failsRcut){
						if(isNewBox){
							fails_box_last_sim[pol]+=weight;
						}
					}
					// fails as last cust with Rcut?
					if(!isSurf[0] && !isSurf[1] && !isSurfEvent_top[pol] && !WFRMS[pol] && !failsCWPowerCut && !isNewBox){
						if(failsRcut){
							fails_rcut_last_sim[pol]+=weight;
						}
					}

				// "in sequence"
					// fails WFRMS first? (same as "as first" for this cut only)
					if(WFRMS[pol] || failsCWPowerCut){
						fails_WFRMS_insequence_sim[pol]+=weight;
					}
					// passes WFRMS, but fails box
					if(!WFRMS[pol] && !failsCWPowerCut && isNewBox){
						fails_box_insequence_sim[pol]+=weight;
					}
					// passes WFRMS and box, but fails surface
					if(!WFRMS[pol] && !failsCWPowerCut && !isNewBox && (isSurf[0] || isSurf[1] || isSurfEvent_top[pol])){
						fails_surface_insequence_sim[pol]+=weight;
					}
					// passes WFRMS, box, and surface, but fails Rcut (same as "as last" for this cut only)
					if(!WFRMS[pol] && !failsCWPowerCut && !isNewBox && (!isSurf[0] && !isSurf[1] && !isSurfEvent_top[pol]) && failsRcut){
						fails_rcut_insequence_sim[pol]+=weight;
					}
			} // loop over pol
		} //loop over sim events

	}

	printf("Sim Num total          :           %7.1f\n",num_total_sim);
	printf("------------------------------------------\n");
	printf("Sim WFRMS              :           %7.1f, %7.1f, %7.1f | %7.1f, %7.1f, %7.1f \n",fails_WFRMS_first_sim[0],fails_WFRMS_insequence_sim[0],fails_WFRMS_last_sim[0],fails_WFRMS_first_sim[1],fails_WFRMS_insequence_sim[1],fails_WFRMS_last_sim[1]);
	printf("Sim Box                :           %7.1f, %7.1f, %7.1f | %7.1f, %7.1f, %7.1f \n",fails_box_first_sim[0],fails_box_insequence_sim[0],fails_box_last_sim[0],fails_box_first_sim[1],fails_box_insequence_sim[1],fails_box_last_sim[1]);
	printf("Sim Surf               :           %7.1f, %7.1f, %7.1f | %7.1f, %7.1f, %7.1f \n",fails_surface_first_sim[0],fails_surface_insequence_sim[0],fails_surface_last_sim[0],fails_surface_first_sim[1],fails_surface_insequence_sim[1],fails_surface_last_sim[1]);
	printf("Sim Rcut               :           %7.1f, %7.1f, %7.1f | %7.1f, %7.1f, %7.1f \n",fails_rcut_first_sim[0],fails_rcut_insequence_sim[0],fails_rcut_last_sim[0],fails_rcut_first_sim[1],fails_rcut_insequence_sim[1],fails_rcut_last_sim[1]);

	printf("V Fit Parameters are %.2f and %.2f \n", fitParams[0], fitParams[1]);
	printf("Found optimal V intercept at %.2f \n", optimal_intercept);

	int colors [28] = { kBlue, kRed, kGreen, kMagenta, kCyan};
	for(int pol=0; pol<2; pol++){
		for(int bin=0; bin<=all_events[pol]->GetNbinsX(); bin++){
			double thrown = all_events[pol]->GetBinContent(bin);
			double pass_soft_short_cal_wfrms_this = pass_soft_short_cal_wfrms[pol]->GetBinContent(bin);
			double pass_soft_short_cal_wfrms_box_this = pass_soft_short_cal_wfrms_box[pol]->GetBinContent(bin);
			double pass_soft_short_cal_wfrms_box_surf_this = pass_soft_short_cal_wfrms_box_surf[pol]->GetBinContent(bin);
			double pass_soft_short_cal_wfrms_box_surf_rcut_this = pass_soft_short_cal_wfrms_box_surf_rcut[pol]->GetBinContent(bin);
			if(thrown>0.){
				eff_soft_short_cal_wfrms[pol]->SetBinContent(bin,pass_soft_short_cal_wfrms_this/thrown);
				eff_soft_short_cal_wfrms_box[pol]->SetBinContent(bin,pass_soft_short_cal_wfrms_box_this/thrown);
				eff_soft_short_cal_wfrms_box_surf[pol]->SetBinContent(bin,pass_soft_short_cal_wfrms_box_surf_this/thrown);
				eff_soft_short_cal_wfrms_box_surf_rcut[pol]->SetBinContent(bin,pass_soft_short_cal_wfrms_box_surf_rcut_this/thrown);
			}
			else{
				eff_soft_short_cal_wfrms[pol]->SetBinContent(bin,0.);
				eff_soft_short_cal_wfrms_box[pol]->SetBinContent(bin,0.);
				eff_soft_short_cal_wfrms_box_surf[pol]->SetBinContent(bin,0.);
				eff_soft_short_cal_wfrms_box_surf_rcut[pol]->SetBinContent(bin,0.);
			}
		}
	}

	TCanvas *c3 = new TCanvas("mycanv","mycanv",2*850,850);
	c3->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		if(pol==0){
			eff_soft_short_cal_wfrms[pol]->SetTitle("Efficiency VPol");
		}
		if(pol==1){
			eff_soft_short_cal_wfrms[pol]->SetTitle("Efficiency HPol");
		}

		c3->cd(pol+1);
			eff_soft_short_cal_wfrms[pol]->Draw("");
			eff_soft_short_cal_wfrms_box[pol]->Draw("same");
			eff_soft_short_cal_wfrms_box_surf[pol]->Draw("same");
			eff_soft_short_cal_wfrms_box_surf_rcut[pol]->Draw("same");


			eff_soft_short_cal_wfrms[pol]->SetLineColor(colors[0]);
			eff_soft_short_cal_wfrms_box[pol]->SetLineColor(colors[1]);
			eff_soft_short_cal_wfrms_box_surf[pol]->SetLineColor(colors[2]);
			eff_soft_short_cal_wfrms_box_surf_rcut[pol]->SetLineColor(colors[3]);

			eff_soft_short_cal_wfrms[pol]->SetLineWidth(2.);
			eff_soft_short_cal_wfrms_box[pol]->SetLineWidth(2.);
			eff_soft_short_cal_wfrms_box_surf[pol]->SetLineWidth(2.);
			eff_soft_short_cal_wfrms_box_surf_rcut[pol]->SetLineWidth(2.);
			
			// if(pol+1==1){
			// 	TLegend *leg = new TLegend(0.48,0.6,0.9,0.9);
			// 	leg->AddEntry(eff_soft_short_cal_wfrms[pol],"Cut WFMRS","l");
			// 	leg->AddEntry(eff_soft_short_cal_wfrms_box[pol],"+Cut Cal Pulser Reco","l");
			// 	leg->AddEntry(eff_soft_short_cal_wfrms_box_surf[pol],"+Cut Surface","l");
			// 	leg->AddEntry(eff_soft_short_cal_wfrms_box_surf_rcut[pol],"+Cut Peak/Corr","l");
			// 	leg->Draw();
			// }

		eff_soft_short_cal_wfrms[pol]->GetXaxis()->SetTitle("3rd Highest Vpeak/RMS");
		eff_soft_short_cal_wfrms[pol]->GetYaxis()->SetTitle("Efficiency (weighted)");

	}
	char efficiency_title[400];
	sprintf(efficiency_title,
			 "%s/optimize/A%d_config%d_E%2.1f_Efficiency.png",plotPath,station,config,224.);
			 // "%s/optimize/%d.%d.%d_A%d_config%d_E%2.1f_Efficiency.png",plotPath,year_now, month_now, day_now,station,config,224.);
	c3->SaveAs(efficiency_title);
	delete c3;
}