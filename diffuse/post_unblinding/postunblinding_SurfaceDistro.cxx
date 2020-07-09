////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	postunblinding_SurfaceDistro.cxx
////	look at the theta distribution
////	this saves the histograms and theta values out to some file
////
////	Nov 2019
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

double getCosZenith(double theta){
	// cout<<"theta is "<<theta<<endl;
	double this_theta = abs(theta-90.);
	// cout<<"after rotation "<<this_theta<<endl;
	this_theta*=TMath::DegToRad();
	// cout<<"after convert to radians "<<this_theta<<endl;
	double retVal = TMath::Cos(this_theta);
	// cout<<"return value is "<<retVal<<endl;
	return retVal;
}

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

	vector<int> BadRunList=BuildBadRunList(station);
	vector<int> BadSurfaceRunList=BuildSurfaceRunList(station);

	int numTotal=0;

	TChain dataVTree("VTree");
	TChain dataHTree("HTree");
	TChain dataAllTree("AllTree");
	TChain dataRecoTree("OutputTreeReco");
	TChain dataFilterTree("OutputTree");
	char the_data[500];

	if(full_or_partial==100 && isOrg==1){
		// use the 100pct sample
		// sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_*.root",station,config);
		// sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_*.root",station,config);
	}
	if(full_or_partial==10 && isOrg==1){
		// use the *original* 10pct sample from optimization
		// sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_*.root",station,config);
		// sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_*.root",station,config);
		sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_*.root",station,config);
	}

	printf("The Data file: %s\n", the_data);
	dataVTree.Add(the_data);
	dataHTree.Add(the_data);
	dataAllTree.Add(the_data);
	dataFilterTree.Add(the_data);

	TH1D *h1_theta[2];
	h1_theta[0] = new TH1D("Vdist_theta","Vdist_theta",180,-90,90);
	h1_theta[1] = new TH1D("Hdist_theta","Hdist_theta",180,-90,90);

	TH1D *h1_costheta[2];
	h1_costheta[0] = new TH1D("Vdist_costheta","Vdist_costheta",180,-1,1);
	h1_costheta[1] = new TH1D("Hdist_costheta","Hdist_costheta",180,-1,1);

	char outputFileName[500];
	sprintf(outputFileName,"%s/reco_and_2d_cut_values_A%d_c%d_10sample.root",outputLocation.c_str(),station,config);
	TFile *outFile = TFile::Open(outputFileName,"RECREATE");
	TTree *outTree = new TTree("outTree","outTree");
	int hist_this_pol[2];
	int theta_val_300[2];
	int phi_val_300[2];
	double corr_val_300[2];
	double snr_val_300[2];
	int runNum_out;
	int eventNum_out;
	int unixTime_out;
	bool isThisABadRun_out;
	bool isThisBadLivetime_out;
	int surf_out[2];
	int surf_top_out[2];
	outTree->Branch("passes_this_pol_V",&hist_this_pol[0]);
	outTree->Branch("passes_this_pol_H",&hist_this_pol[1]);
	outTree->Branch("theta_val_300_V",&theta_val_300[0]);
	outTree->Branch("theta_val_300_H",&theta_val_300[1]);
	outTree->Branch("phi_val_300_V",&phi_val_300[0]);
	outTree->Branch("phi_val_300_H",&phi_val_300[1]);
	outTree->Branch("corr_val_300_V",&corr_val_300[0]);
	outTree->Branch("corr_val_300_H",&corr_val_300[1]);
	outTree->Branch("snr_val_300_V",&snr_val_300[0]);
	outTree->Branch("snr_val_300_H",&snr_val_300[1]);
	outTree->Branch("runNum_out",&runNum_out);
	outTree->Branch("eventNum_out",&eventNum_out);
	outTree->Branch("unixTime_out",&unixTime_out);
	outTree->Branch("isThisABadRun_out",&isThisABadRun_out);
	outTree->Branch("isThisBadLivetime_out",&isThisBadLivetime_out);
	outTree->Branch("surf_V",&surf_out[0]);
	outTree->Branch("surf_h",&surf_out[1]);
	outTree->Branch("surf_top_out_V",&surf_top_out[0]);
	outTree->Branch("surf_top_out_H",&surf_top_out[1]);


	int numDataEvents = dataVTree.GetEntries();

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
		
		// check both
		bool isThisABadRun = isBadRun(station,runNum,BadRunList);
		if(!isThisABadRun){
			// isThisABadRun = isBadRun(station, runNum, BadSurfaceRunList);
		}

		cout<<"Num data entries is "<<numDataEvents<<endl;
		for(int event=0; event<numDataEvents; event++){
			dataVTree.GetEvent(event);
			dataHTree.GetEvent(event);
			dataAllTree.GetEvent(event);
			dataFilterTree.GetEvent(event);

			if(runNum!=currentRunNum){
				currentRunNum=runNum;
				isThisABadRun = isBadRun(station,runNum, BadRunList);
				if(!isThisABadRun){
					// isThisABadRun = isBadRun(station,runNum, BadSurfaceRunList);
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
				hist_this_pol[pol]=-1000;
				theta_val_300[pol]=-1000;
				phi_val_300[pol]=-1000;
				corr_val_300[pol]=-1000.;
				snr_val_300[pol]=-1000.;
				surf_out[pol]=0;
				surf_top_out[pol]=0;
			}
			runNum_out=runNum;
			eventNum_out=eventNumber;
			unixTime_out=unixTime;
			isThisABadRun_out=isThisABadRun;
			// continue;
			if( isSoft || isBadEvent || hasBadSpareChanIssue || hasBadSpareChanIssue2 || isFirstFiveEvent || isShort || isCal){
			// if( isSoft || isBadEvent || hasBadSpareChanIssue || hasBadSpareChanIssue2 || isFirstFiveEvent || isShort || isCal || isThisABadRun){
				continue;
			}
			bool local_badlivetime = isBadLivetime(station,unixTime);
			// if(){
			// 	// continue;
			// }
			isThisBadLivetime_out=local_badlivetime;

			for(int pol=0; pol<2; pol++){
				// if(!WFRMS[pol] && !isNewBox && !isSurf[0] && !isSurf[1] && !isSurfEvent_top[pol]){
				if(!WFRMS[pol] && !isNewBox){
				// if(!isNewBox){
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

						hist_this_pol[pol]=1;

						h1_theta[pol]->Fill(double(theta_300[pol]));
						// cout<<"Theta is "<<theta_300[pol]<<" and function return is "<<getCosZenith(double(theta_300[pol]))<<endl;
						h1_costheta[pol]->Fill(getCosZenith(double(theta_300[pol])));
						
					}// not failing CW power cut?
				}// passes rest of analysis (not WFRMS, box, surface)
			}// loop over polarizations

			if(hist_this_pol[0] || hist_this_pol[1]){
				// set or re-set this stuff
				for(int pol2=0; pol2<2; pol2++){
					theta_val_300[pol2]=theta_300[pol2];
					phi_val_300[pol2]=phi_300[pol2];
					corr_val_300[pol2]=corr_val[pol2];
					snr_val_300[pol2]=snr_val[pol2];
					surf_out[pol2]=isSurf[pol2];
					surf_top_out[pol2]=isSurfEvent_top[pol2];
				}
				runNum_out=runNum;
				eventNum_out=eventNumber;
				unixTime_out=unixTime;
				outTree->Fill();
			}

		}// loop over events
	}

	outTree->Write();
	outFile->Close();


	char fileTitle[300];
	sprintf(fileTitle,"/users/PAS0654/osu0673/A23_analysis_new2/results/post_unblinding/A%d_c%d_SurfaceDistro.root",station,config);
	TFile *fOut = new TFile(fileTitle, "RECREATE");
	fOut->cd();
	h1_theta[0]->Write();
	h1_theta[1]->Write();
	h1_costheta[0]->Write();
	h1_costheta[1]->Write();
	fOut->Close();

	double bot = 25.;
	double top = 55.;
	double bot_coszen = getCosZenith(bot);
	double top_coszen = getCosZenith(top);

	gStyle->SetOptStat(0);
	TCanvas *cDistro = new TCanvas("","",2*1100,2*850);
	cDistro->Divide(2,2);
	for(int pol=0; pol<2; pol++){
		cDistro->cd(pol+1);
			h1_theta[pol]->Draw("");
			h1_theta[pol]->GetXaxis()->SetRangeUser(bot,top);
			h1_theta[pol]->GetXaxis()->SetTitle("#theta [deg]");
			h1_theta[pol]->GetYaxis()->SetTitle("Number of Events");
			gPad->SetLogy();
		cDistro->cd(pol+3);
			h1_costheta[pol]->Draw("");
			h1_costheta[pol]->GetXaxis()->SetRangeUser(bot_coszen, top_coszen);
			h1_costheta[pol]->GetXaxis()->SetTitle("cos(#theta)");
			h1_costheta[pol]->GetYaxis()->SetTitle("Number of Events");
			gPad->SetLogy();
	}
	char save_temp_title[400];
	sprintf(save_temp_title,"%s/post_unblinding/A%d_c%d_ThetaDistributions_%dEvents.png",plotPath,station,config,numTotal);
	cDistro->SaveAs(save_temp_title);
	delete cDistro;
}