////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	add_branch_updategeocut.cxx
////	Add a branch to an existing root file, adjusting geo cut values
////
////	Nov 2019
////    Jorge Torres
////    updated by baclark to be "sim" friendly
////	and then updated again to 
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <memory>
#include <vector>
#include <cmath>
#include <ctime>
#include <TSystem.h>
#include <sys/stat.h>


//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"

#include "RawAraStationEvent.h"
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraEventCalibrator.h"
#include "AraStationInfo.h"
#include "AraQualCuts.h"
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_inputParameters.h"
#include "tools_outputObjects.h"
#include "tools_Cuts.h"
#include "tools_WaveformFns.h"

using namespace std;

int main(int argc, char **argv){

	if(argc<5){
		cout<< "Usage\n" << argv[0] << " <1-isSim?> <2-station> <3-config> <4-ValForCuts filename 1...> <5-ValForCuts filename 2... >"<<endl;;
		return -1;
	}
	bool isSimulation = atoi(argv[1]);
	int station = atoi(argv[2]);
	int config = atoi(argv[3]);
	bool drop_chan=false;
	if(config>2) drop_chan=true;

	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL){
		std::cout << "Warning! $PED_DIR is not set!" << endl;
		return -1;
	}
	AraQualCuts *qualCut = AraQualCuts::Instance();

	for(int file_num=4; file_num<argc; file_num++){

		if(gSystem->AccessPathName(argv[file_num])){
			std::cout << "file does not exist" << std::endl;
			return -1;
		}

		string chRun = "run";
		string file = string(argv[file_num]);
		size_t foundRun = file.find(chRun);
		string strRunNum = file.substr(foundRun+4,4);
		int runNum_in = atoi(strRunNum.c_str());
		printf("On file %d, run number %d \n", file_num-2, runNum_in);

		TFile *f = new TFile(argv[file_num],"update");

		TTree *VTree = (TTree*)f->Get("VTree");
		TTree *HTree = (TTree*)f->Get("HTree");
		TTree *AllTree = (TTree*)f->Get("AllTree");

		int theta_300[2];
		int phi_300[2];
		VTree->SetBranchAddress("theta_300_V_new",&theta_300[0]);
		VTree->SetBranchAddress("phi_300_V_new",&phi_300[0]);
		HTree->SetBranchAddress("theta_300_H_new",&theta_300[1]);
		HTree->SetBranchAddress("phi_300_H_new",&phi_300[1]);
		int eventNumber;
		AllTree->SetBranchAddress("eventNumber",&eventNumber);
		int isCal;
		int box;
		AllTree->SetBranchAddress("cal",&isCal);
		AllTree->SetBranchAddress("box",&box);

		// now, the new things we want to set up as outputs
		// first the surface cut
		int surf_out[2];
		TBranch *newSurfVBranch = AllTree->Branch("surf_V_new2",&surf_out[0]);
		TBranch *newSurfHBranch = AllTree->Branch("surf_H_new2",&surf_out[1]);
		int isNewBox_out;
		TBranch *newBoxBranch = AllTree->Branch("box300",&isNewBox_out);

		int numEntries = AllTree->GetEntries();
		// cout << "Number of events is " << numEntries << endl;
		// numEntries=4000;

		for(int event=0; event<numEntries; event++){

			surf_out[0]=0;
			surf_out[1]=0;
			isNewBox_out=0;

			VTree->GetEvent(event);
			HTree->GetEvent(event);
			AllTree->GetEvent(event);

			// first, update the surface cut
			if(theta_300[0]>17){
				surf_out[0]=1;

			}
			if(theta_300[1]>17){
				surf_out[1]=1;
			}
			// printf("Event %d, thetaV %d, is surf now %d \n", eventNumber, theta_300[0], surf_out[0]);
			// printf("Event %d, thetaH %d, is surf now %d \n", eventNumber, theta_300[1], surf_out[1]);
			// printf("\n");

			// and now update the cal pulser cut
			bool isCP5=false;
			bool isCP6=false;
			for(int pol=0; pol<2; pol++){
				identifyCalPulser(station,config, theta_300[pol], phi_300[pol], isCP5, isCP6);
			}
			if(isCP5 || isCP6){
				isNewBox_out=1;
			}
			// printf("Event %3d, isCal %d, isOldBox %d, isNewBox_out %d\n", eventNumber, isCal, box, isNewBox_out);

			newSurfVBranch->Fill();
			newSurfHBranch->Fill();
			newBoxBranch->Fill();

		}
		f->cd();
		AllTree->Write(0,TObject::kOverwrite);
		delete f;
	}
}
