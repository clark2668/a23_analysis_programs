////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	unblind_surface_SaveToTree.cxx
////	unblind the surface cut
////	if an event passes surface cut, save it to a reduced root file
////
////	Sep 2019
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

using namespace std;

int main(int argc, char **argv)
{
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *toolsPath(getenv("TOOLS_DIR"));
	if (toolsPath == NULL) {
		std::cout << "Warning! $TOOLS_DIR is not set!" << endl;
		return -1;
	}

	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <1-station> <2-config> <3-save_vals_file_1> <4-save_vals_file_2> .... <n-save_vals_file_n>"<<endl;;
		return -1;
	}
	int station = atoi(argv[1]);
	int config = atoi(argv[2]);

	if(station!=2 && station!=3){
		printf("No good! You asked for station %d, but this code only works for stations 2 and 3 \n",station);
		return -1;
	}
	
	vector<int> BadRunList=BuildBadRunList(station);

	for(int file_num=3; file_num<argc; file_num++){
		
		TFile *inputFile = TFile::Open(argv[file_num],"read");
		if(!inputFile){
			cout<<"Can't open vals for cuts file!"<<endl;			
			return -1;
		}
		TTree *inputVTree = (TTree*) inputFile->Get("VTree");
		TTree *inputHTree = (TTree*) inputFile->Get("HTree");
		TTree *inputAllTree = (TTree*) inputFile->Get("AllTree");
		TTree *inputRecoTree = (TTree*) inputFile->Get("OutputTreeReco");
		TTree *inputFilterTree = (TTree*) inputFile->Get("OutputTree");

		double corr_val[2];
		double snr_val[2];
		int WFRMS[2];
		int theta_300[2];
		int phi_300[2];
		int theta_41[2];
		int phi_41[2];

		int Refilt[2];
		inputVTree->SetBranchAddress("Refilt_V",&Refilt[0]);
		inputHTree->SetBranchAddress("Refilt_H",&Refilt[1]);

		inputVTree->SetBranchAddress("corr_val_V_new",&corr_val[0]);
		inputVTree->SetBranchAddress("snr_val_V_new",&snr_val[0]);
		inputVTree->SetBranchAddress("wfrms_val_V_new",&WFRMS[0]);
		inputVTree->SetBranchAddress("theta_300_V_new",&theta_300[0]);
		inputVTree->SetBranchAddress("theta_41_V_new",&theta_41[0]);
		inputVTree->SetBranchAddress("phi_300_V_new",&phi_300[0]);
		inputVTree->SetBranchAddress("phi_41_V_new",&phi_41[0]);

		inputHTree->SetBranchAddress("corr_val_H_new",&corr_val[1]);
		inputHTree->SetBranchAddress("snr_val_H_new",&snr_val[1]);
		inputHTree->SetBranchAddress("wfrms_val_H_new",&WFRMS[1]);
		inputHTree->SetBranchAddress("theta_300_H_new",&theta_300[1]);
		inputHTree->SetBranchAddress("theta_41_H_new",&theta_41[1]);
		inputHTree->SetBranchAddress("phi_300_H_new",&phi_300[1]);
		inputHTree->SetBranchAddress("phi_41_H_new",&phi_41[1]);

		int isCal;
		int isSoft;
		int isShort;
		int isCW;
		int isNewBox;

		inputAllTree->SetBranchAddress("cal",&isCal);
		inputAllTree->SetBranchAddress("soft",&isSoft);
		inputAllTree->SetBranchAddress("short",&isShort);
		inputAllTree->SetBranchAddress("CW",&isCW);
		inputAllTree->SetBranchAddress("box",&isNewBox);

		int isSurf[2]; // a surface event after filtering?
		int isSurfEvent_top[2]; // a top event?

		inputAllTree->SetBranchAddress("surf_V_new2",&isSurf[0]);
		inputAllTree->SetBranchAddress("surf_H_new2",&isSurf[1]);

		inputAllTree->SetBranchAddress("surf_top_V",&isSurfEvent_top[0]);
		inputAllTree->SetBranchAddress("surf_top_H",&isSurfEvent_top[1]);

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

		inputAllTree->SetBranchAddress("bad",&isBadEvent);
		inputAllTree->SetBranchAddress("weight",&weight);
		inputAllTree->SetBranchAddress("unixTime",&unixTime);
		inputAllTree->SetBranchAddress("isFirstFiveEvent",&isFirstFiveEvent);
		inputAllTree->SetBranchAddress("hasBadSpareChanIssue",&hasBadSpareChanIssue);
		inputAllTree->SetBranchAddress("hasBadSpareChanIssue2",&hasBadSpareChanIssue2);
		inputAllTree->SetBranchAddress("runNum",&runNum);
		inputAllTree->SetBranchAddress("eventNumber",&eventNumber);
		inputAllTree->SetBranchAddress("isSpikey",&isSpikey);
		inputAllTree->SetBranchAddress("isCliff",&isCliff);
		inputAllTree->SetBranchAddress("OutofBandIssue",&OutofBandIssue);
		inputAllTree->SetBranchAddress("bad_v2",&bad_v2);
		inputAllTree->SetBranchAddress("isRFEvent",&isRFEvent);
		inputAllTree->SetBranchAddress("isPayloadBlast2",&isPayloadBlast2);
		inputAllTree->SetBranchAddress("box300",&box300);

		double frac_of_power_notched_V[8];
		double frac_of_power_notched_H[8];

		stringstream ss;
		for(int i=0; i<8; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			inputVTree->SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_V[i]);
		}
		for(int i=8; i<16; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			inputHTree->SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_H[i-8]);
		}

		double VPeakOverRMS[16];
		int waveformLength[16];
		inputFilterTree->SetBranchAddress("VPeakOverRMS",&VPeakOverRMS);
		inputFilterTree->SetBranchAddress("waveformLength",&waveformLength);

		int numEntries = inputVTree->GetEntries();

		inputAllTree->GetEvent(0);

		bool isThisABadRun = isBadRun(station,runNum,BadRunList);
		bool isThisASoftDomRun = isSoftwareDominatedRun(toolsPath, station, runNum);

		if(isThisABadRun || isThisASoftDomRun){
			inputFile->Close();
			delete inputFile;
		}

		// prepare the output file
		// just write it straight to disk because idc anymore, and the "write to disk rate" should be small
		char outputFileName[400];
		sprintf(outputFileName,"/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/unblind_surface_redo/A%d/c%d/surface_events_run%d.root",station, config, runNum);
		TFile *outputFile = TFile::Open(outputFileName, "RECREATE");
		TTree *trees[6];
		trees[0] = inputVTree->CloneTree(0);
		trees[1] = inputHTree->CloneTree(0);
		trees[2] = inputAllTree->CloneTree(0);
		trees[3] = inputRecoTree->CloneTree(0);
		trees[4] = inputFilterTree->CloneTree(0);
		trees[5]= new TTree("SurfaceTree","SurfaceTree");
		bool passesV;
		bool passesH;
		bool passesBoth;
		trees[5]->Branch("passesV",&passesV,"passesV/O");
		trees[5]->Branch("passesH",&passesH,"passesH/O");
		trees[5]->Branch("passesBoth",&passesBoth,"passesBoth/O");

		int numFoundThisRun=0;

		for(int event=0; event<numEntries; event++){
			inputVTree->GetEvent(event);
			inputHTree->GetEvent(event);
			inputAllTree->GetEvent(event);
			inputRecoTree->GetEvent(event);
			inputFilterTree->GetEvent(event);

			passesV=false;
			passesH=false;
			passesBoth=false;

			// continue;
			if( isSoft || isBadEvent || hasBadSpareChanIssue || hasBadSpareChanIssue2 || isFirstFiveEvent || isShort || isCal){
				continue;
			}
			// new, since we learned more about A3
			if(!isRFEvent || isThisASoftDomRun || isSpikey || isCliff || OutofBandIssue || bad_v2 || isPayloadBlast2){
				continue;
			}
			if(isBadLivetime(station,unixTime)){
				continue;
			}
			bool passes_this_pol[2] = {false};
			for(int pol=0; pol<2; pol++){
				if(!WFRMS[pol] 
					&& !isNewBox && !box300
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
							passes_this_pol[pol]=true;
						}
					}// not failing CW power cut?
				}// passes rest of analysis (not WFRMS, box, surface)
			}// loop over polarizations

			if(passes_this_pol[0]){
				passesV=true;
			}
			if(passes_this_pol[1]){
				passesH=true;
			}
			if(passes_this_pol[0] && passes_this_pol[1]){
				passesBoth=true;
			}

			if(passes_this_pol[0] || passes_this_pol[1]){
				trees[0]->Fill();
				trees[1]->Fill();
				trees[2]->Fill();
				trees[3]->Fill();
				trees[4]->Fill();
				trees[5]->Fill();
				numFoundThisRun++;
			}
		}// loop over events

		inputFile->Close();
		delete inputFile;

		outputFile->Write();
		outputFile->Close();
		delete outputFile;

		printf("In run %d, found %d events \n", runNum, numFoundThisRun);

		// if we didn't find anything, remove the output file to decrease clutter (this should be the vast majority of cases)
		if(numFoundThisRun==0){
			remove(outputFileName);
		}

	} // loop over input files
}