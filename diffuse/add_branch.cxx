////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	add_branch.cxx
////	Add a branch to an existing root file
////
////	Nov 2019
////    Jorge Torres
////    updated by baclark to be "sim" friendly
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <memory>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fftw3.h>
#include <ctime>
#include "time.h" // for time convert
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
#include "AraGeomTool.h"
#include "FFTtools.h"
#include "AraStationInfo.h"
#include "AraQualCuts.h"
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_inputParameters.h"
#include "tools_outputObjects.h"
#include "tools_Cuts.h"
#include "tools_WaveformFns.h"

TStyle* RootStyle();

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
		TTree *T = (TTree*)f->Get("AllTree");
		bool isSpikey;
		bool isCliff;
		bool OutofBandIssue;
		bool bad_v2;
		bool isRFEvent;
		bool isPayloadBlast;

		int eventNumber;
		T->SetBranchAddress("eventNumber",&eventNumber);

		TBranch *isSpikeyBranch = T->Branch("isSpikey",&isSpikey);
		TBranch *isCliffBranch = T->Branch("isCliff",&isCliff);
		TBranch *hasOutofBandIssueBranch = T->Branch("OutofBandIssue",&OutofBandIssue);
		TBranch *bad_v2Branch = T->Branch("bad_v2",&bad_v2);
		TBranch *isRFEventBranch = T->Branch("isRFEvent",&isRFEvent);
		TBranch *isPayloadBlastBranch = T->Branch("isPayloadBlast2",&isPayloadBlast);


		char dataFName[70];
		if(!isSimulation){
			sprintf(dataFName,"/fs/project/PAS0654/ARA_DATA/A23/10pct/RawData/A3/all_runs/event%i.root",runNum_in);
		}
		else if(isSimulation){
			sprintf(dataFName,"/fs/project/PAS0654/ARA_DATA/A23/sim/RawSim/A%d/c%d/E224/AraOut.A%d_c%d_E224.txt.run%d.root",station,config,station,config,runNum_in);
		}

		TFile *fp = TFile::Open(dataFName);
		if(!fp) {
			std::cout << "Can't open file\n";
			return -1;
		}
		TTree *eventTree;
		eventTree= (TTree*) fp->Get("eventTree");
		if(!eventTree) {
			std::cout << "Can't find eventTree\n";
			return -1;
		}
		RawAtriStationEvent *rawEvPtr=0;
		UsefulAtriStationEvent *realAtriEvPtr_fullcalib=0;
		if(!isSimulation){
			eventTree->SetBranchAddress("event",&rawEvPtr);
		}
		else{
			eventTree->SetBranchAddress("UsefulAtriStationEvent",&realAtriEvPtr_fullcalib);
		}

		if(!isSimulation){
			AraEventCalibrator *calib = AraEventCalibrator::Instance(); //make a calibrator
			char ped_file_name[400];
			sprintf(ped_file_name,"%s/run_specific_peds/A%d/all_peds/event%d_specificPeds.dat",PedDirPath,station,runNum_in);
			calib->setAtriPedFile(ped_file_name,station); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist
		}

		int numEntries = T->GetEntries();
		cout << "Number of events is " << numEntries << endl;

		for(int event=0; event<numEntries; event++){
			// reset some variables
			isSpikey=false;
			isCliff=false;
			OutofBandIssue=false;
			bad_v2=false;
			isRFEvent=false;
			isPayloadBlast=false;

			eventTree->GetEvent(event);
			T->GetEvent(event);

			// only if it's a *data* event do we need to build the real event, and, do qual cuts and isRFTrigger on it
			if(!isSimulation){
				realAtriEvPtr_fullcalib = new UsefulAtriStationEvent(rawEvPtr, AraCalType::kLatestCalib); //make the event
				if(!qualCut->isGoodEvent(realAtriEvPtr_fullcalib)) bad_v2=true;
				isRFEvent=rawEvPtr->isRFTrigger();
				int evt_num = rawEvPtr->eventNumber;//event number
			}
			vector<TGraph*> graphs;
			for (int i = 0; i < 16; i++){
				TGraph* gr = realAtriEvPtr_fullcalib->getGraphFromRFChan(i);
				graphs.push_back(gr);
			}

			if(!isSimulation){
				// only call the deletion on this pointer if it's data, and therefore was created locally
				delete realAtriEvPtr_fullcalib;
			}

			if(isCliffEvent(graphs)) isCliff=true;
			if(isSpikeyStringEvent(station,drop_chan,graphs,config)) isSpikey=true;
			if(hasOutofBandIssue(graphs,drop_chan)) OutofBandIssue=true;
			if(isHighPowerStringEvent(graphs, station, config)) isPayloadBlast=true;

			for (int i=0; i < graphs.size(); i++){
				delete graphs[i];
			}
			graphs.clear();

			isSpikeyBranch->Fill();
			isCliffBranch->Fill();
			hasOutofBandIssueBranch->Fill();
			bad_v2Branch->Fill();
			isRFEventBranch->Fill();
			isPayloadBlastBranch->Fill();
		}
		f->cd();
		T->Write(0,TObject::kOverwrite);
		delete f;
		delete fp;
	}
}
