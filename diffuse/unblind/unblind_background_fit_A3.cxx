////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	unblind_background_fit.cxx
////	unblind the fits to the background region
////
////	Oct 2019
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
#include "TTimeStamp.h"

//AraRoot includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraQualCuts.h"

// analysis custom
#include "tools_Cuts.h"
#include "tools_Stats.h"
#include "tools_CommandLine.h"
#include "tools_outputObjects.h"
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"

// for ray trace correlator
#include "Settings.h"
#include "Event.h"
#include "Detector.h"
#include "Report.h"
#include "RayTraceCorrelator.h"

AraAntPol::AraAntPol_t Vpol = AraAntPol::kVertical;
AraAntPol::AraAntPol_t Hpol = AraAntPol::kHorizontal;

using namespace std;
int PlotThisEvent(int station, int runNum, int event, int problempol);

int main(int argc, char **argv)
{

	time_t time_now = time(0); //get the time now
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	stringstream ss;
	gStyle->SetOptStat(0);

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;


	if(argc<5){
		cout<< "Usage\n" << argv[0] << " <1-station> <2-config> <100 or 10> <isOrg> <output_location>"<<endl;;
		return -1;
	}
	int station = atoi(argv[1]);
	int config = atoi(argv[2]);
	int full_or_partial = atoi(argv[3]);
	int isOrg = atoi(argv[4]);
	string outputLocation = argv[5];

	if(station!=2 && station!=3){
		printf("No good! You asked for station %d, but this code only works for stations 2 and 3 \n",station);
		return -1;
	}

	double side_cut_corr[2] = {0.005, 0.004}; // in corr, side region is <=0.005 for V and <=0.004 in H
	double side_cut_snr[2] = {4.,4.}; // in snr, side region is <=4 in V and H


	// get the slope and intercept of the slanted line for V and H
	double vslope, hslope;
	double vintercept, hintercept;
	getRCutValues(station, config, 0, vslope, vintercept);
	getRCutValues(station, config, 1, hslope, hintercept);
	double selected_slopes[2] = {vslope, hslope};
	double selected_intercepts[2] = {vintercept, hintercept};

	printf("V slope and intercept: %.2f, %.2f \n", vslope, vintercept);
	printf("H slope and intercept: %.2f , %.2f \n", hslope, hintercept);

	// let's see if we can do it with arrays....
	int numSNRbins=300; //let this go all the way up to 30
	double numEventsPassed[2][numSNRbins]; // number of events passing in a polarization, and for that snr bin
	double numEventsPassed_diff[2][numSNRbins]; // the difference between the number of events cut between this bin and the next

	double dRcut=0.1; // SNR bin incremets (bins in 0.1 SNR units)
	double Rcutmin=0.; //mininum value is 0
	double intercept[2][numSNRbins]; // what numerical value corresponds to that bin in the snr (intercept) scan
	for(int pol=0; pol<2; pol++){
		for(int bin=0; bin<numSNRbins; bin++){
			numEventsPassed[pol][bin]=0.; //initialize...
			numEventsPassed_diff[pol][bin]=0.;
			intercept[pol][bin] = Rcutmin + double(bin)*dRcut;
		}
	}

	vector<int> BadRunList=BuildBadRunList(station);
	vector<int> BadSurfaceRunList=BuildSurfaceRunList(station);

	int numTotal=0;

	// need to be able to make the final 2D distribution
	double max=0.05;
	TH2D *h2SNRvsCorr[2]; // SNR on Y axis, Corr on X axis, like in the TB
	h2SNRvsCorr[0]=new TH2D("2DDistroV","2DDistroV",100,0,max,300,0,30);
	h2SNRvsCorr[1]=new TH2D("2DDistroH","2DDistroH",100,0,max,300,0,30);

	char outputFileName[500];
	sprintf(outputFileName,"%s/2d_cut_values_A%d_c%d_sample%d.root",outputLocation.c_str(),station,config,full_or_partial);
	TFile *outFile = TFile::Open(outputFileName,"RECREATE");
	TTree *outTree = new TTree("outTree","outTree");
	int hist_this_pol[2];
	double corr_val_out[2];
	double snr_val_out[2];
	int runNum_out;
	int eventNum_out;
	outTree->Branch("passes_this_pol_V",&hist_this_pol[0]);
	outTree->Branch("passes_this_pol_H",&hist_this_pol[1]);
	outTree->Branch("corr_val_V",&corr_val_out[0]);
	outTree->Branch("corr_val_H",&corr_val_out[1]);
	outTree->Branch("snr_val_V",&snr_val_out[0]);
	outTree->Branch("snr_val_H",&snr_val_out[1]);
	outTree->Branch("runNum_out",&runNum_out);
	outTree->Branch("eventNum_out",&eventNum_out);

	TChain dataVTree("VTree");
	TChain dataHTree("HTree");
	TChain dataAllTree("AllTree");
	TChain dataRecoTree("OutputTreeReco");
	TChain dataFilterTree("OutputTree");
	char the_data[500];

	if(full_or_partial==100){
		if (config==1) sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_CWThresh2.0_snrbins_0_1_wfrmsvals_-1.2_-1.3_run_*.root",station,config);
  	if (config==2) sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_CWThresh2.0_snrbins_0_1_wfrmsvals_-1.3_-1.4_run_*.root",station,config);
  	if (config==3 || config==4) sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_CWThresh2.0_snrbins_0_1_wfrmsvals_-1.0_-1.1_run_*.root",station,config);
  	if (config==5) sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_CWThresh2.0_snrbins_0_1_wfrmsvals_-0.7_-0.8_run_*.root",station,config);
	}
	if(full_or_partial==10 && isOrg==1){
		// use the *original* 10pct sample from optimization
		if (config==1) sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/ValsForCuts_v2/A%d/c%d/cutvals_drop_FiltSurface_CWThresh2.0_snrbins_0_1_wfrmsvals_-1.2_-1.3_run_*.root",station,config);
		if (config==2) sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/ValsForCuts_v2/A%d/c%d/cutvals_drop_FiltSurface_CWThresh2.0_snrbins_0_1_wfrmsvals_-1.3_-1.4_run_*.root",station,config);
		if (config==3 || config==4) sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/ValsForCuts_v2/A%d/c%d/cutvals_drop_FiltSurface_CWThresh2.0_snrbins_0_1_wfrmsvals_-1.0_-1.1_run_*.root",station,config);
		if (config==5) sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/ValsForCuts_v2/A%d/c%d/cutvals_drop_FiltSurface_CWThresh2.0_snrbins_0_1_wfrmsvals_-0.7_-0.8_run_*.root",station,config);
	}
	if(full_or_partial==10 && isOrg==0){
		// use the *new* 10pct sample from the slightly revised code base
		// sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/10pct_verify_try2/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_*.root",station,config);
	}

	printf("The data: %s\n", the_data);
	dataVTree.Add(the_data);
	dataHTree.Add(the_data);
	dataAllTree.Add(the_data);
	dataFilterTree.Add(the_data);

	int numDataEvents = dataVTree.GetEntries();
	// numDataEvents=100;

	// do this inside brackets for scoping power and re-use of identical variable names when it comes time for simulation to happen
	{

		double corr_val[2];
		double snr_val[2];
		int WFRMS[2];
		int theta_300[2];
		int phi_300[2];
		int theta_41[2];
		int phi_41[2];

		int Refilt[2];
		dataVTree.SetBranchAddress("Refilt_V",&Refilt[0]);
		dataHTree.SetBranchAddress("Refilt_H",&Refilt[1]);

		dataVTree.SetBranchAddress("corr_val_V_new",&corr_val[0]);
		dataVTree.SetBranchAddress("snr_val_V_new",&snr_val[0]);
		dataVTree.SetBranchAddress("wfrms_val_V_new",&WFRMS[0]);
		dataVTree.SetBranchAddress("theta_300_V_new",&theta_300[0]);
		dataVTree.SetBranchAddress("theta_41_V_new",&theta_41[0]);
		dataVTree.SetBranchAddress("phi_300_V_new",&phi_300[0]);
		dataVTree.SetBranchAddress("phi_41_V_new",&phi_41[0]);

		dataHTree.SetBranchAddress("corr_val_H_new",&corr_val[1]);
		dataHTree.SetBranchAddress("snr_val_H_new",&snr_val[1]);
		dataHTree.SetBranchAddress("wfrms_val_H_new",&WFRMS[1]);
		dataHTree.SetBranchAddress("theta_300_H_new",&theta_300[1]);
		dataHTree.SetBranchAddress("theta_41_H_new",&theta_41[1]);
		dataHTree.SetBranchAddress("phi_300_H_new",&phi_300[1]);
		dataHTree.SetBranchAddress("phi_41_H_new",&phi_41[1]);

		int isCal;
		int isSoft;
		int isShort;
		int isCW;
		int isNewBox;

		dataAllTree.SetBranchAddress("cal",&isCal);
		dataAllTree.SetBranchAddress("soft",&isSoft);
		dataAllTree.SetBranchAddress("short",&isShort);
		dataAllTree.SetBranchAddress("CW",&isCW);
		dataAllTree.SetBranchAddress("box",&isNewBox);

		int isSurf[2]; // a surface event after filtering?
		int isSurfEvent_top[2]; // a top event?

		dataAllTree.SetBranchAddress("surf_V_new2",&isSurf[0]);
		dataAllTree.SetBranchAddress("surf_H_new2",&isSurf[1]);

		dataAllTree.SetBranchAddress("surf_top_V",&isSurfEvent_top[0]);
		dataAllTree.SetBranchAddress("surf_top_H",&isSurfEvent_top[1]);

		int isBadEvent;
		double weight;
		int unixTime;
		int isFirstFiveEvent;
		int hasBadSpareChanIssue;
		int hasBadSpareChanIssue2;
		int runNum;
		int eventNumber;
		bool isSpikey;
    bool isCliff;
    bool OutofBandIssue;
    bool bad_v2;
    bool isRFEvent;
    bool isPayloadBlast2;
    int box300;

		dataAllTree.SetBranchAddress("bad",&isBadEvent);
		dataAllTree.SetBranchAddress("weight",&weight);
		dataAllTree.SetBranchAddress("unixTime",&unixTime);
		dataAllTree.SetBranchAddress("isFirstFiveEvent",&isFirstFiveEvent);
		dataAllTree.SetBranchAddress("hasBadSpareChanIssue",&hasBadSpareChanIssue);
		dataAllTree.SetBranchAddress("hasBadSpareChanIssue2",&hasBadSpareChanIssue2);
		dataAllTree.SetBranchAddress("runNum",&runNum);
		dataAllTree.SetBranchAddress("eventNumber",&eventNumber);
		dataAllTree.SetBranchAddress("isSpikey",&isSpikey);
		dataAllTree.SetBranchAddress("isCliff",&isCliff);
		dataAllTree.SetBranchAddress("OutofBandIssue",&OutofBandIssue);
		dataAllTree.SetBranchAddress("bad_v2",&bad_v2);
		dataAllTree.SetBranchAddress("isRFEvent",&isRFEvent);
		dataAllTree.SetBranchAddress("isPayloadBlast2",&isPayloadBlast2);
		dataAllTree.SetBranchAddress("box300",&box300);

		double frac_of_power_notched_V[8];
		double frac_of_power_notched_H[8];

		stringstream ss;
		for(int i=0; i<8; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			dataVTree.SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_V[i]);
		}
		for(int i=8; i<16; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			dataHTree.SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_H[i-8]);
		}

		double VPeakOverRMS[16];
		int waveformLength[16];
		dataFilterTree.SetBranchAddress("VPeakOverRMS",&VPeakOverRMS);
		dataFilterTree.SetBranchAddress("waveformLength",&waveformLength);

		int numEntries = dataVTree.GetEntries();
		numTotal+=numEntries;

		dataAllTree.GetEvent(0);
		int currentRunNum = runNum;
		bool isThisABadRun = isBadRun(station,runNum,BadRunList);
		bool isThisASoftDomRun = isSoftwareDominatedRun("/users/PCON0003/cond0068/ARA/AraRoot/analysis/a23_analysis_tools", station, runNum);

		if(!isThisABadRun){
			isThisABadRun = isBadRun(station, runNum, BadSurfaceRunList);
		}

		for(int event=0; event<numDataEvents; event++){
			dataVTree.GetEvent(event);
			dataHTree.GetEvent(event);
			dataAllTree.GetEvent(event);
			dataFilterTree.GetEvent(event);

			if(runNum!=currentRunNum){
				currentRunNum=runNum;
				isThisABadRun = isBadRun(station,runNum, BadRunList);
				isThisASoftDomRun = isSoftwareDominatedRun("/users/PCON0003/cond0068/ARA/AraRoot/analysis/a23_analysis_tools", station, runNum);

				if(!isThisABadRun){
					isThisABadRun = isBadRun(station,runNum, BadSurfaceRunList);
				}
				if(isThisABadRun){
					printf(RED"*"RESET);
					// printf("     Yup, run %d is bad \n",runNum);
				}
				else{
					printf(GREEN"*"RESET);
				}
			}

			// set or re-set this stuff
			for(int pol=0; pol<2; pol++){
				hist_this_pol[pol]=0;
				corr_val_out[pol]=-10000.;
				snr_val_out[pol]=-10000.;
			}
			runNum_out=runNum;
			eventNum_out=eventNumber;

			// continue;
			if( isSoft || isBadEvent || hasBadSpareChanIssue || hasBadSpareChanIssue2 || isFirstFiveEvent || isShort || isCal || isThisABadRun || isThisASoftDomRun|| isSpikey || isCliff || OutofBandIssue || bad_v2 || isPayloadBlast2 || box300){
				continue;
			}
			if(isBadLivetime(station,unixTime)){
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
						bool this_pass_R_cut = passesRCut(station, config, pol, snr_val[pol], corr_val[pol]);

						// careful now, only check for events *failing* the R-cut
						if(!this_pass_R_cut){

							// if(pol==1 && runNum==1454 && snr_val[pol]>8){
							// 	printf("Run %4d, Event %6d, Pol %d, unixTime %d, theta %3d, phi %3d, coor %.4f, CW status %d\n", runNum, eventNumber, 1, unixTime, theta_300[1], phi_300[1], corr_val[1], Refilt[1]);
							// 	PlotThisEvent(station, runNum, eventNumber, 1);
							// }

							h2SNRvsCorr[pol]->Fill(corr_val[pol],snr_val[pol],weight);

							// set these correctly
							hist_this_pol[pol]=1;
							corr_val_out[pol]=corr_val[pol];
							snr_val_out[pol]=snr_val[pol];
							runNum_out=runNum;
							eventNum_out=eventNumber;
							outTree->Fill(); // save them to the file

							// scan over the choices of intercept, and record if this event would be cut by this choice of intercept
							for(int bin=0; bin<numSNRbins; bin++){ //check all potential intercepts
								double this_y_val = (selected_slopes[pol] * corr_val[pol] ) + intercept[pol][bin]; // compute the SNR to pass at this intercept
								// printf("For bin %d, with intercept %.2f, SNR to pass is %.2f \n", bin, intercept[bin], this_y_val);
								if(snr_val[pol]>=this_y_val){ // does this event pass?
									// printf("     This event has SNR %.2f, so it passes!\n",snr_val[pol]);
									// printf("          Old number of events in this bin %5f \n", numEventsPassed[bin]);
									numEventsPassed[pol][bin]+=1.;
									// printf(".         New number of events in this bin %5f \n", numEventsPassed[bin]);
								} // does event pass Rcut
							} // loop over SNR cuts
						}
					}// not failing CW power cut?
				}// passes rest of analysis (not WFRMS, box, surface)
			}// loop over polarizations
		}// loop over events
	}

	outFile->Write(); // write this down
	outFile->Close(); // close it up

	// now, loop over all the bins, and make the "differetial distribution"
	// that is, how many events are cut by a given movement in the snr cut
	// that's why it's a differential
	for(int pol=0; pol<2; pol++){
		for(int bin=0; bin<numSNRbins-1; bin++){
			numEventsPassed_diff[pol][bin] = numEventsPassed[pol][bin] - numEventsPassed[pol][bin+1];
			// printf("Pol %d Bin %d at cut %.2f has %5.f events passing, and next bin has %5f, so diff is %.5f \n", pol, bin, intercept[pol][bin],numEventsPassed[pol][bin],numEventsPassed[pol][bin+1],numEventsPassed_diff[pol][bin]);
		}
	}

	TH1D *hEventsVsSNR[2];
	hEventsVsSNR[0] = new TH1D("DiffDistroV","DiffDistroV",numSNRbins,0,30);
	hEventsVsSNR[1] = new TH1D("DiffDistroH","DiffDistroH",numSNRbins,0,30);

	// now, we have to do the exponential fit
	// which I chose to code as half-way down the distribution between the maximum bin and the last filled bin

	int max_bin[2];
	int last_filled_bin_above_2[2];
	int fit_start_bin[2];
	double start_of_fit[2];
	int last_filled_bin[2];
	double end_of_fit[2];

	for(int pol=0; pol<2; pol++){
		for(int bin=0; bin<numSNRbins; bin++){
			hEventsVsSNR[pol]->SetBinContent(bin+1,numEventsPassed_diff[pol][bin]);
		}
		max_bin[pol] = hEventsVsSNR[pol]->GetMaximumBin();
		// cout<<"Max bin is "<<max_bin[pol]<<endl;
		last_filled_bin_above_2[pol] = hEventsVsSNR[pol]->FindLastBinAbove(2.,1) + 2;
		// cout<<"Last filled bin above 2 is "<<last_filled_bin_above_2[pol]<<endl;
		fit_start_bin[pol] = int((last_filled_bin_above_2[pol] - max_bin[pol])/2) + max_bin[pol]; //start half-way between the peak bin and the last filled bin
		// cout<<"Fit start bin is "<<fit_start_bin[pol]<<endl;
		start_of_fit[pol] = hEventsVsSNR[pol]->GetBinCenter(fit_start_bin[pol]);
		// cout<<"Start of fit is "<<start_of_fit[pol]<<endl;
		last_filled_bin[pol] = hEventsVsSNR[pol]->FindLastBinAbove(0.,1);
		end_of_fit[pol] = hEventsVsSNR[pol]->GetBinCenter(last_filled_bin[pol]+2.); //go two bins more just to make sure fit is over
		// printf("Pol %d Start of fit is %.2f and end of fit is %.2f \n", pol, start_of_fit[pol], end_of_fit[pol]);

		printf("Pol %d: Last filled bin is bin %d and value %.2f \n", pol, last_filled_bin[pol], hEventsVsSNR[pol]->GetBinCenter(last_filled_bin[pol]));
		printf("Pol %d: Max bin is bin %d and value %.2f \n", pol, max_bin[pol], hEventsVsSNR[pol]->GetBinCenter(max_bin[pol]));
		printf("Pol %d: Proposed start of fit is bin %d and value %.2f \n", pol, fit_start_bin[pol], hEventsVsSNR[pol]->GetBinCenter(fit_start_bin[pol]));
	}

	char fileTitle[300];
	sprintf(fileTitle,"/users/PCON0003/cond0068/ARA/AraRoot/analysis/unblind/background_fit/A%d_c%d_sample%d.root",station,config,full_or_partial);
	TFile *fOut = new TFile(fileTitle, "RECREATE");
	fOut->cd();
	hEventsVsSNR[0]->Write();
	hEventsVsSNR[1]->Write();
	h2SNRvsCorr[0]->Write();
	h2SNRvsCorr[1]->Write();
	fOut->Close();

	// now we actually do the exponential fit
	char equation[150];
	sprintf(equation,"exp([0]*x+[1])");
	char equation_name[2][150];
	TF1 *fit[2];
	int status[2];
	double fitParams[2][2];
	double fitParamErrors[2][2];
	for(int pol=0; pol<2; pol++){
		sprintf(equation_name[pol],"ExpFit%d",pol);
		fit[pol] = new TF1(equation_name[pol],equation,start_of_fit[pol],end_of_fit[pol]);
		status[pol] = hEventsVsSNR[pol]->Fit(equation_name[pol],"LL,R");
		fitParams[pol][0] = fit[pol]->GetParameter(0);
		fitParams[pol][1] = fit[pol]->GetParameter(1);
		fitParamErrors[pol][0] = fit[pol]->GetParError(0);
		fitParamErrors[pol][1] = fit[pol]->GetParError(1);
		printf("Pol %d Fit Parameters are %.2f and %.2f \n", pol, fitParams[pol][0], fitParams[pol][1]);
	}

	// and we record some of the information about the fit
	// like the name, the number of obsered events, etc.
	double binWidthIntercept[2];
	double leftBoundary[2];
	double rightBoundary[2];
	double numBinsThisFit[2];
	char this_fit_title[2][400];
	TH1D *hNumObserved[2];
	for(int pol=0; pol<2; pol++){
		binWidthIntercept[pol] = hEventsVsSNR[pol]->GetBinWidth(1);
		leftBoundary[pol] = start_of_fit[pol] - binWidthIntercept[pol]/2.;
		rightBoundary[pol] = end_of_fit[pol] + binWidthIntercept[pol]/2.;
		numBinsThisFit[pol] = (rightBoundary[pol] - leftBoundary[pol])/binWidthIntercept[pol] + 1;
		sprintf(this_fit_title[pol],"Fit_Pol%d_Slope%.2f",pol,selected_slopes[pol]);
		hNumObserved[pol] = new TH1D(this_fit_title[pol],"",numBinsThisFit[pol],leftBoundary[pol],rightBoundary[pol]);
		for(int bin=0; bin<numBinsThisFit[pol]; bin++){
			double originalContent = hEventsVsSNR[pol]->GetBinContent(bin+fit_start_bin[pol]);
			hNumObserved[pol]->SetBinContent(bin+1,originalContent);
		}
	}

	char fit_title_words[2][400];

	TCanvas *cRcut = new TCanvas("","",6*850,2*850);
	cRcut->Divide(6,2);
	for(int pol=0; pol<2; pol++){
		cRcut->cd(pol+1+(pol==0 ? 0 : 5)); // for corr vs snr, data
			// h2SNRvsCorr[pol]->Draw("colz");
			// h2SNRvsCorr[pol]->GetXaxis()->SetTitle("Correlation Value");
			// h2SNRvsCorr[pol]->GetYaxis()->SetTitle("SNR");
			// gPad->SetLogz();
			// cut_lines[pol]->Draw("same");
			// cut_lines[pol]->SetLineColor(kRed);
		cRcut->cd(pol+2+(pol==0 ? 0 : 5)); // for differential distribution, zoom out
			hEventsVsSNR[pol]->Draw("");
			hEventsVsSNR[pol]->GetXaxis()->SetTitle("SNR Cut (y-intercept value)");
			hEventsVsSNR[pol]->GetYaxis()->SetTitle("Number of Events Cut");
			hEventsVsSNR[pol]->SetTitle("Differential Distribution");
			gPad->SetLogy();
		cRcut->cd(pol+3+(pol==0 ? 0 : 5)); // for differential distribution, zoom in
			hNumObserved[pol]->Draw("HIST");
			hNumObserved[pol]->GetXaxis()->SetTitle("SNR Cut (y-intercept value)");
			hNumObserved[pol]->GetYaxis()->SetTitle("Number of Events Cut");
			sprintf(fit_title_words[pol],"Fit exp((%.2f +- %.2f)x + (%.2f +- %.2f)",fitParams[pol][0], fitParamErrors[pol][0],fitParams[pol][1], fitParamErrors[pol][1]);
			hNumObserved[pol]->SetTitle(fit_title_words[pol]);
			fit[pol]->Draw("same");
			fit[pol]->SetLineColor(kRed);
			gPad->SetLogy();

	char save_title[400];
	// sprintf(save_title,"%s/optimize/A%d_config%d_Final_VSlope_%.2f_HSlope_%.2f_VInt_%.2f_Hint_%.2f.png",
	sprintf(save_title,"%s/unblind/background_fit/%d.%d.%d_A%d_config%d_sample%d_Final_VSlope_%.2f_HSlope_%.2f_VInt_%.2f_Hint_%.2f.png",
						plotPath,
						year_now, month_now, day_now,
						station,
						config,
						full_or_partial,
						selected_slopes[0],
						selected_slopes[1],
						selected_intercepts[0],
						selected_intercepts[1]);
	// cRcut->SaveAs(save_title);
	}
}

int PlotThisEvent(int station, int runNum, int event, int problempol){
	time_t time_now = time(0); //get the time now
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *DataDirPath(getenv("DATA_DIR_100"));
	if (DataDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

	char run_file_name[400];
	sprintf(run_file_name,"%s/RawData/A%d/all_runs/event%d.root",DataDirPath,station,runNum);
	TFile *mapFile = TFile::Open(run_file_name);
	if(!mapFile){
		cout<<"Can't open data file for map!"<<endl;
		return -1;
	}
	TTree *eventTree = (TTree*) mapFile-> Get("eventTree");
	if(!eventTree){
		cout<<"Can't find eventTree for map"<<endl;
		return -1;
	}

	RawAtriStationEvent *rawPtr =0;
	eventTree->SetBranchAddress("event",&rawPtr);
	eventTree->GetEvent(event);

	int stationID = rawPtr->stationId;
	char ped_file_name[400];
	sprintf(ped_file_name,"%s/run_specific_peds/A%d/all_peds/event%d_specificPeds.dat",PedDirPath,station,runNum);
	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	calibrator->setAtriPedFile(ped_file_name,stationID); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist

	AraQualCuts *qualCut = AraQualCuts::Instance();
	UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawPtr,AraCalType::kLatestCalib);
	printf("Run %d, Event %d \n", runNum, realAtriEvPtr->eventNumber);
	printf("	Is Quality Event? %d \n", qualCut->isGoodEvent(realAtriEvPtr));

	int unixTime = (int)rawPtr->unixTime;
	int unixTimeUs =(int)rawPtr->unixTimeUs;
	int timeStamp = (int)rawPtr->timeStamp;
	printf("	Unixtime is %d \n", unixTime);
	printf("	Unixtime microsecond is %d \n", unixTimeUs);
	printf("	timeStamp is %d \n", timeStamp);

	stringstream ss1;
	string xLabel, yLabel;
	xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
	vector<string> titlesForGraphs;
	for (int i = 0; i < 16; i++){
		ss1.str("");
		ss1 << "Channel " << i;
		titlesForGraphs.push_back(ss1.str());
	}

	vector <TGraph*> waveforms = makeGraphsFromRF(realAtriEvPtr,16,xLabel,yLabel,titlesForGraphs);
	vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(waveforms, 0.5, xLabel, yLabel, titlesForGraphs);
	vector<TGraph*> grWaveformsPadded = makePaddedGraphs(grWaveformsInt, 0, xLabel, yLabel, titlesForGraphs);
	xLabel = "Frequency (Hz)"; yLabel = "Power Spectral Density (mV/Hz)";
	vector<TGraph*> grWaveformsPowerSpectrum = makePowerSpectrumGraphs(grWaveformsPadded, xLabel, yLabel, titlesForGraphs);


	bool doMapsAtAll=false;
	if(doMapsAtAll){
		// get the run summary information, if it exists yet
		// and remember, because it's the users job to pass the location of the filter files
		// this should work for simulated events just fine
		char filter_file_name[400];
		// gonna try all of them cuz we lazy af lol
		sprintf(filter_file_name,"%s/processed_station_%d_run_%d_filter.root","/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ProcessedFile/A2/2015",station,runNum);
		bool hasFilterFile = false;
		TFile *filterFile = TFile::Open(filter_file_name);
		if(!filterFile){
			sprintf(filter_file_name,"%s/processed_station_%d_run_%d_filter.root","/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ProcessedFile/A2/2013",station,runNum);
		}
		filterFile = TFile::Open(filter_file_name);
		if(!filterFile){
			sprintf(filter_file_name,"%s/processed_station_%d_run_%d_filter.root","/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ProcessedFile/A2/2014",station,runNum);
		}
		filterFile = TFile::Open(filter_file_name);
		if(!filterFile){
			sprintf(filter_file_name,"%s/processed_station_%d_run_%d_filter.root","/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ProcessedFile/A2/2016",station,runNum);
		}
		filterFile = TFile::Open(filter_file_name);
		if(!filterFile){
			return -1;
		}

		TTree *filterTree;
		double VPeakOverRMS[16];
		if(filterFile){
			printf("Successfully found filter file information \n");
			hasFilterFile=true;
			filterTree = (TTree*) filterFile->Get("OutputTree");
			if(!filterTree) {
				std::cout << "Can't find filterTree\n";
				return -1;
			}
			filterTree->SetBranchAddress("VPeakOverRMS", &VPeakOverRMS);
			filterFile->cd();
		}


		filterTree->GetEvent(event);

		vector<double> chan_SNRs;
		if(hasFilterFile){
			for(int i=0; i<16; i++){
				chan_SNRs.push_back(VPeakOverRMS[i]);
			}
		}

		vector <int> chan_list_V;
		vector <int> chan_list_H;
		for(int chan=0; chan<=7; chan++){
			chan_list_V.push_back(chan);
			chan_list_H.push_back(chan+8);
		}

		if(station==2){
			//for station 2, we need to exclude channel 15 from the analysis
			chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
		}
		else if(station==3){
			// // for station 3 years 2014, 2015, 2016, we need to drop string 4 (channels 3, 7, 11, 15) altogether above some run
			// if(
			// 	(!isSimulation && runNum>getA3BadRunBoundary())
			// 	||
			// 	(isSimulation && yearConfig>2)

			// ){			// drop string four
			// 	chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 3), chan_list_V.end());
			// 	chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 7), chan_list_V.end());
			// 	chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 11), chan_list_H.end());
			// 	chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
			// }
		}

		Settings *settings = new Settings();
		string setupfile = "setup.txt";
		settings->ReadFile(setupfile);
		cout << "Read " << setupfile << " file!" << endl;
		settings->NOFZ=1;
		Detector *detector = 0;

		RayTraceCorrelator *theCorrelator = new RayTraceCorrelator(station, 300, settings, 1, RTTestMode);

		int solNum = 0;
		TH2D *map_V_raytrace = theCorrelator->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V, chan_SNRs, solNum);
		TH2D *map_H_raytrace = theCorrelator->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H, chan_SNRs, solNum);

		int new_theta;
		int new_phi;
		double new_corr;

		getCorrMapPeak(map_V_raytrace,new_theta, new_phi, new_corr);
		char title_for_map[300];
		sprintf(title_for_map,"VMap: peak theta %d, phi %d, corr %.4f",new_theta, new_phi,new_corr);
		map_V_raytrace->SetTitle(title_for_map);

		int new_theta_H;
		int new_phi_H;
		double new_corr_H;

		getCorrMapPeak(map_H_raytrace,new_theta_H, new_phi_H, new_corr_H);
		char title_for_map_H[300];
		sprintf(title_for_map_H,"HMap: peak theta %d, phi %d, corr %.4f",new_theta_H, new_phi_H,new_corr_H);
		map_H_raytrace->SetTitle(title_for_map_H);

		bool print_maps = true;
		if(print_maps){
			gStyle->SetOptStat(0);
			TCanvas *cMaps = new TCanvas("","",2*1100,850);
			// TCanvas *cMaps = new TCanvas("","",1100,1.1*850);
			cMaps->Divide(2,1);
				cMaps->cd(1);
				map_V_raytrace->Draw("colz");
				gPad->SetRightMargin(0.15);
				cMaps->cd(2);
				map_H_raytrace->Draw("colz");
				gPad->SetRightMargin(0.15);
			char save_temp_title[400];
			sprintf(save_temp_title,"/users/PCON0003/cond0068/ARA/AraRoot/analysis/unblind/surface/trouble_events/surface_%d.%d.%d_Run%d_Ev%d_Maps.png",year_now,month_now,day_now,runNum,event);
			cMaps->SaveAs(save_temp_title);
			delete cMaps;
		}

		filterFile->Close();
	} // do maps at all

	// bool doContributingMaps=true;
	// if(doContributingMaps){
	// 	stringstream ss1;
	// 	vector<string> titlesForGraphs;
	// 	vector <TH2D*> individuals;

	// 	double SNR_scaling=0.;

	// 	for(int i=0; i<7; i++){
	// 		for(int j=i+1; j<8; j++){
	// 			ss1.str("");
	// 			ss1<<"Pair "<<i<<" and "<<j;
	// 			titlesForGraphs.push_back(ss1.str());
	// 			TH2D *map = theCorrelator->getInterferometricMap_RT_NewNormalization_PairSelect(settings, detector, realAtriEvPtr, Vpol, isSimulation, i, j, solNum);

	// 			// now do the SNR scaling
	// 			double this_snr_product = chan_SNRs[i] * chan_SNRs[j];
	// 			// printf("Weighting term for %d, %d is %.3f \n", i, j, this_snr_product);
	// 			map->Scale(this_snr_product);

	// 			SNR_scaling+=this_snr_product;

	// 			individuals.push_back(map);
	// 		}
	// 	}

	// 	for(int i=0; i<individuals.size(); i++){
	// 		individuals[i]->Scale(1./SNR_scaling);
	// 	}

	// 	// compute the average myself manually
	// 	TH2D *average = (TH2D*) individuals[0]->Clone();
	// 	for(int i=0; i<individuals.size(); i++){
	// 		average->Add(individuals[i]);
	// 	}
	// 	// average->Scale(1./28);
	// 	average->SetTitle("Summed Maps");

	// 	vector<double> mins;
	// 	vector<double> maxs;
	// 	for(int i=0; i<individuals.size(); i++){
	// 		mins.push_back(individuals[i]->GetMinimum());
	// 		maxs.push_back(individuals[i]->GetMaximum());
	// 	}
	// 	std::sort(mins.begin(), mins.end()); //sort smallest to largest
	// 	std::sort(maxs.begin(), maxs.end()); //sort smallest to largest
	// 	std::reverse(maxs.begin(), maxs.end()); //reverse order to get largest to smallest

	// 	TCanvas *cMaps2 = new TCanvas("","",8*850,4*850);
	// 	cMaps2->Divide(7,5);
	// 	for(int i=0; i<individuals.size(); i++){
	// 		cMaps2->cd(i+1);
	// 		individuals[i]->Draw("colz");
	// 		individuals[i]->GetZaxis()->SetRangeUser(mins[0],maxs[0]);
	// 		individuals[i]->SetTitle(titlesForGraphs[i].c_str());
	// 		gStyle->SetTitleFontSize(0.07);
	// 	}
	// 	cMaps2->cd(35);
	// 		average->Draw("colz");
	// 	char save_temp_title[400];
	// 	sprintf(save_temp_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/unblind/surface/surface_%d.%d.%d_Run%d_Ev%d_AllMaps.png",year_now,month_now,day_now,runNum,event);
	// 	cMaps2->SaveAs(save_temp_title);
	// 	delete cMaps2;
	// }


	// need dummy canvases becaus root can't make things pretty to save its life
	vector<TGraph*> dummy;
	for(int i=0; i<16; i++){
		vector<double> thisX;
		vector<double> thisY;
		thisX.push_back(-200);
		thisX.push_back(800);
		thisY.push_back(-700);
		thisY.push_back(700);
		dummy.push_back(new TGraph(2,&thisX[0], &thisY[0]));
	}

	char save_temp_title[300];
	sprintf(save_temp_title,"%s/unblind/background_fit/peculiar_events/A%d_Run%d_Ev%d_ProblemPol%d_Waveforms.png",plotPath,station,runNum,event,problempol);
	TCanvas *cWave = new TCanvas("","",4*1100,4*850);
	cWave->Divide(4,4);
	for(int i=0; i<16; i++){
		cWave->cd(i+1);
		dummy[i]->Draw("AP");
		dummy[i]->SetLineColor(kWhite);
		dummy[i]->GetXaxis()->SetRangeUser(-200.,700.);
		dummy[i]->GetYaxis()->SetRangeUser(-300.,300.);

		waveforms[i]->Draw("sameL");
		waveforms[i]->SetLineWidth(3);
	}
	cWave->SaveAs(save_temp_title);
	delete cWave;

	sprintf(save_temp_title,"%s/unblind/background_fit/peculiar_events/A%d_Run%d_Ev%d_ProblemPol%d_Spectra.png",plotPath,station,runNum,event,problempol);
	TCanvas *cSpec = new TCanvas("","",4*1100,4*850);
	cSpec->Divide(4,4);
	for(int i=0; i<16; i++){
		cSpec->cd(i+1);
		grWaveformsPowerSpectrum[i]->Draw("AL");
		grWaveformsPowerSpectrum[i]->SetLineWidth(3);
		gPad->SetLogy();
	}
	cSpec->SaveAs(save_temp_title);
	delete cSpec;

	for(int i=0; i<16; i++){
		delete waveforms[i];
		delete grWaveformsInt[i];
		delete grWaveformsPadded[i];
		delete grWaveformsPowerSpectrum[i];
	}
	delete realAtriEvPtr;
	mapFile->Close();
	delete mapFile;
	return 0;
}
