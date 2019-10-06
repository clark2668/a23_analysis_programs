////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	unblind_surface.cxx
////	unblind the surface cut
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
#include "TGraph.h"

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

void configure(TH1D *gr);

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


	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <1-station> <2-config> <100 or 10> <isOrg>"<<endl;;
		return -1;
	}
	int station = atoi(argv[1]);
	int config = atoi(argv[2]);
	int full_or_partial = atoi(argv[3]);
	int isOrg = atoi(argv[4]);

	if(station!=2 && station!=3){
		printf("No good! You asked for station %d, but this code only works for stations 2 and 3 \n",station);
		return -1;
	}

	// a place to print things to file
	char title_txt[200];
	sprintf(title_txt,"%s/unblind/surface/A%d_c%d_%dSample_EventListsigmavsfreq_ch%d.txt", plotPath,station,config,full_or_partial);

	vector<int> BadRunList=BuildBadRunList(station);

	int numTotal=0;

	// we also want to work up a plot for the surface rate
	TTimeStamp start;
	TTimeStamp stop;
	int numBins;

	// run 6635
	start.Set(2015, 12, 29, 18, 00, 0, 0, true, 0);
	stop.Set(2015, 12, 29, 22, 00, 0, 0, true, 0);
	numBins=100;

	// run 6655
	start.Set(2016, 01, 02, 18, 40, 0, 0, true, 0);
	stop.Set(2016, 01, 02, 18, 42, 0, 0, true, 0);

	// run 6554
	start.Set(2015, 12, 13, 00, 00, 0, 0, true, 0);
	stop.Set(215, 12, 14, 00, 00, 0, 0, true, 0);	

	int start_bin = start.GetSec();
	int stop_bin = stop.GetSec();

	TH2D *events_vs_time[2];
	for(int i=0; i<2; i++){
		if(i==0){
			ss.str("");
			ss<<"VPol";
		}
		else if(i==1){
			ss.str("");
			ss<<"HPol";
		}
		events_vs_time[i] = new TH2D(ss.str().c_str(),ss.str().c_str(),numBins, start_bin, stop_bin, 360, -180, 180);
		events_vs_time[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		events_vs_time[i]->GetXaxis()->SetTimeFormat("%H:%M:%S");
		events_vs_time[i]->GetXaxis()->SetTimeOffset(0.,"GMT");
	}

	vector<int> interesting_unixTimes_V;
	vector<int> interesting_Phis_V;

	vector<int> interesting_unixTimes_H;
	vector<int> interesting_Phis_H;

	TChain dataVTree("VTree");
	TChain dataHTree("HTree");
	TChain dataAllTree("AllTree");
	TChain dataRecoTree("OutputTreeReco");
	TChain dataFilterTree("OutputTree");
	char the_data[500];

	if(full_or_partial==100){
		// use the 100pct sample
		// sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_*.root",station,config);

	}
	if(full_or_partial==10 && isOrg==1){
		// use the *original* 10pct sample from optimization
		// sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_*.root",station,config);
	}
	if(full_or_partial==10 && isOrg==0){
		// use the *new* 10pct sample from the slightly revised code base
		// sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/10pct_verify_try2/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_*.root",station,config);
	}

	// printf("The data: %s\n", the_data);
	// dataVTree.Add(the_data);
	// dataHTree.Add(the_data);
	// dataAllTree.Add(the_data);
	// dataFilterTree.Add(the_data);

	// if(full_or_partial==10 && isOrg==1){ // to make plotting just a little easier

	// 	vector<int> runs_to_loop_over;
	// 	if(config==1){
	// 		runs_to_loop_over.push_back(2636);
	// 		runs_to_loop_over.push_back(2662);
	// 		runs_to_loop_over.push_back(2678);
	// 		runs_to_loop_over.push_back(2990);
	// 		runs_to_loop_over.push_back(3206);
	// 	}
	// 	if(config==2){
	// 		runs_to_loop_over.push_back(1830);
	// 		runs_to_loop_over.push_back(1835);
	// 		runs_to_loop_over.push_back(1860);
	// 		runs_to_loop_over.push_back(2090);
	// 		runs_to_loop_over.push_back(2091);
	// 		runs_to_loop_over.push_back(2155);
	// 	}
	// 	if(config==4){
	// 		runs_to_loop_over.push_back(4398);
	// 		runs_to_loop_over.push_back(4497);
	// 		runs_to_loop_over.push_back(4777);
	// 		runs_to_loop_over.push_back(4817);
	// 		runs_to_loop_over.push_back(4837);
	// 		runs_to_loop_over.push_back(4842);
	// 		runs_to_loop_over.push_back(4979);
	// 		runs_to_loop_over.push_back(5004);
	// 		runs_to_loop_over.push_back(5007);
	// 		runs_to_loop_over.push_back(5032);
	// 		runs_to_loop_over.push_back(5174);
	// 		runs_to_loop_over.push_back(5369);
	// 		runs_to_loop_over.push_back(5515);
	// 		runs_to_loop_over.push_back(5516);
	// 		runs_to_loop_over.push_back(5517);
	// 		runs_to_loop_over.push_back(5586);
	// 		runs_to_loop_over.push_back(5615);
	// 		runs_to_loop_over.push_back(5649);
	// 		runs_to_loop_over.push_back(5660);
	// 		runs_to_loop_over.push_back(5664);
	// 		runs_to_loop_over.push_back(5666);
	// 		runs_to_loop_over.push_back(5670);
	// 		runs_to_loop_over.push_back(5680);
	// 		runs_to_loop_over.push_back(5699);
	// 		runs_to_loop_over.push_back(5701);
	// 		runs_to_loop_over.push_back(5702);
	// 		runs_to_loop_over.push_back(5704);
	// 		runs_to_loop_over.push_back(5775);
	// 		runs_to_loop_over.push_back(5849);
	// 		runs_to_loop_over.push_back(6080);
	// 		runs_to_loop_over.push_back(6445);
	// 	}
	// 	if(config==5){
	// 		runs_to_loop_over.push_back(6532);
	// 		runs_to_loop_over.push_back(6536);
	// 		runs_to_loop_over.push_back(6542);
	// 		runs_to_loop_over.push_back(6554);
	// 		runs_to_loop_over.push_back(6577);
	// 		runs_to_loop_over.push_back(6635);
	// 		runs_to_loop_over.push_back(6655);
	// 		runs_to_loop_over.push_back(6669);
	// 		runs_to_loop_over.push_back(6674);
	// 		runs_to_loop_over.push_back(6679);
	// 		runs_to_loop_over.push_back(6694);
	// 		runs_to_loop_over.push_back(6705);
	// 		runs_to_loop_over.push_back(6733);
	// 		runs_to_loop_over.push_back(6818);
	// 		runs_to_loop_over.push_back(7836);
	// 		runs_to_loop_over.push_back(8074);
	// 	}

	// 	for(int i=0; i<runs_to_loop_over.size(); i++){
	// 		sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_%d.root",station,config, runs_to_loop_over[i]);
	// 		dataVTree.Add(the_data);
	// 		dataHTree.Add(the_data);
	// 		dataAllTree.Add(the_data);
	// 		dataFilterTree.Add(the_data);
	// 	}
	// }


	sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_6554.root",station,config);
	// sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_6554.root",station,config);
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

		dataAllTree.SetBranchAddress("surf_V_new",&isSurf[0]);
		dataAllTree.SetBranchAddress("surf_H_new",&isSurf[1]);

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

		dataAllTree.SetBranchAddress("bad",&isBadEvent);
		dataAllTree.SetBranchAddress("weight",&weight);
		dataAllTree.SetBranchAddress("unixTime",&unixTime);
		dataAllTree.SetBranchAddress("isFirstFiveEvent",&isFirstFiveEvent);
		dataAllTree.SetBranchAddress("hasBadSpareChanIssue",&hasBadSpareChanIssue);
		dataAllTree.SetBranchAddress("hasBadSpareChanIssue2",&hasBadSpareChanIssue2);
		dataAllTree.SetBranchAddress("runNum",&runNum);
		dataAllTree.SetBranchAddress("eventNumber",&eventNumber);

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

		for(int event=0; event<numDataEvents; event++){
			dataVTree.GetEvent(event);
			dataHTree.GetEvent(event);
			dataAllTree.GetEvent(event);
			dataFilterTree.GetEvent(event);

			if(runNum!=currentRunNum){
				currentRunNum=runNum;
				isThisABadRun = isBadRun(station,runNum, BadRunList);
				if(isThisABadRun){
					printf(RED"*"RESET);
					// printf("     Yup, run %d is bad \n",runNum);
				}
				else{
					printf(GREEN"*"RESET);
				}
			}

			// continue;
			if( isSoft || isBadEvent || hasBadSpareChanIssue || hasBadSpareChanIssue2 || isFirstFiveEvent || isShort || isCal || isThisABadRun){
				continue;
			}
			if(isBadLivetime(station,unixTime)){
				continue;
			}
			for(int pol=0; pol<2; pol++){
				if(!WFRMS[pol] 
					&& !isNewBox 
					&& (isSurf[0] || isSurf[1] || isSurfEvent_top[pol])
				){
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

						double which_corr_to_use = corr_val[pol];
						bool this_pass_R_cut = passesRCut(station, config, pol, snr_val[pol], which_corr_to_use);

						if(this_pass_R_cut){

							if(pol==0){
								interesting_unixTimes_V.push_back(unixTime);
								interesting_Phis_V.push_back(phi_300[pol]);
							}
							else if(pol==1){
								interesting_unixTimes_H.push_back(unixTime);
								interesting_Phis_H.push_back(phi_300[pol]);
							}	

							// FILE *fout = fopen(title_txt, "a");
							// fprintf(fout,"Run %4d, Event %6d, Pol %1d, unixTime %d, theta %3d, phi %3d, coor %.4f, refilt %d, surfv %d, surfh %d, surftopppol %d \n",runNum, eventNumber, pol, unixTime, theta_300[pol], phi_300[pol], corr_val[pol], Refilt[pol], isSurf[0], isSurf[1], isSurfEvent_top[pol]);
							// fclose(fout);//close sigmavsfreq.txt file

							// printf(BLUE"\n Run %d, Event %d, pol %d, unixTime %d, reco theta %d, phi %d, Refilt status %d, corr is %.4f \n"RESET, runNum, eventNumber, pol, unixTime, theta_300[pol], phi_300[pol], Refilt[pol], corr_val[pol]);
							// PlotThisEvent(station, config, runNum, eventNumber, pol);
						}

					}// not failing CW power cut?
				}// passes rest of analysis (not WFRMS, box, surface)
			}// loop over polarizations
		}// loop over events
	}
	std::cout<<endl;


	TGraph *grPhiTracking[2];
	grPhiTracking[0] = new TGraph(interesting_unixTimes_V.size(), &interesting_unixTimes_V[0], &interesting_Phis_V[0]);
	grPhiTracking[1] = new TGraph(interesting_unixTimes_H.size(), &interesting_unixTimes_H[0], &interesting_Phis_H[0]);

	TCanvas *cSurfaceRate = new TCanvas("","",2.1*850,2.1*850);
	cSurfaceRate->Divide(1,2);
	for(int i=0; i<2; i++){
		cSurfaceRate->cd(i+1);
		// num_samps[i]->GetYaxis()->SetRangeUser(0.1,2e7);
		// num_samps[i]->GetXaxis()->SetRangeUser(0,750);

		// gPad->SetTopMargin(0.1);
		// gPad->SetRightMargin(0.03);
		// gPad->SetLeftMargin(0.05);
		// gPad->SetBottomMargin(0.11);

		// configure(events_vs_time[i]);

		events_vs_time[i]->Draw("");
		events_vs_time[i]->GetXaxis()->SetTitle("UnixTime");
		events_vs_time[i]->GetYaxis()->SetTitle("Phi [deg]");
		grPhiTracking[i]->Draw("samepl");
		grPhiTracking[i]->SetMarkerStyle(kCircle);
	}
	char thistitle[300];
	sprintf(thistitle, "%s/unblind/surface/%d.%d.%d_A%d_c%d_Surface_%dSample_PhiVsTime.png",plotPath,year_now,month_now,day_now,station,config,full_or_partial);
	cSurfaceRate->SaveAs(thistitle);
	delete cSurfaceRate;

}


void configure(TH1D *gr){
	gr->GetXaxis()->SetLabelSize(0.05);
	gr->GetXaxis()->SetTitleSize(0.05);
	
	gr->GetYaxis()->SetLabelSize(0.05);
	gr->GetYaxis()->SetTitleSize(0.05);
	// gr->GetYaxis()->SetTitleOffset(1.1);
	
	gStyle->SetTitleFontSize(0.05);
	gr->GetXaxis()->SetNdivisions(12,0,0,false);
	gr->SetLineWidth(2);
}