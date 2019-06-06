////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  2018.12.23-LowFreqPower.cxx 
////  store fraction of spectrum power below 75 MHz
////
////  Dec 2018
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <string>
#include <sstream>

//AraRoot Includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraQualCuts.h"
#include "FFTtools.h"

#include "tools_PlottingFns.h"
#include "tools_WaveformFns.h"
#include "tools_Cuts.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"

using namespace std;

int main(int argc, char **argv)
{

	if(argc<4) {
		std::cout << "Usage\n" << argv[0] << " <station> <input_file> <output_location> "<<endl;
		return -1;
	}
	int station = atoi(argv[1]);

	/*
	arguments
	0: exec
	1: station
	2: input data file
	3: output location
	*/
	
	TFile *fpIn = TFile::Open(argv[2]);
	if(!fpIn) {
		std::cout << "Can't open file\n";
		return -1;
	}
	TTree *eventTree = (TTree*) fpIn->Get("eventTree");
	if(!eventTree) {
		std::cout << "Can't find eventTree\n";
		return -1;
	}
	RawAtriStationEvent *rawAtriEvPtr=0;
	eventTree->SetBranchAddress("event",&rawAtriEvPtr);
	int run;
	eventTree->SetBranchAddress("run",&run);
	eventTree->GetEntry(0);
	printf("Run Number %d \n", run);

	char ped_file_name[400];
	sprintf(ped_file_name,"/fs/scratch/PAS0654/ara/peds/run_specific_peds/A%d/all_peds/event%d_specificPeds.dat",station,run);
	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	calibrator->setAtriPedFile(ped_file_name,station); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist

	AraQualCuts *qualCut = AraQualCuts::Instance(); //we also need a qual cuts tool

	char outfile_name[400];
	sprintf(outfile_name,"%s/spare_channel_stats_run%d.root",argv[3],run);

	TFile *fpOut = TFile::Open(outfile_name, "RECREATE");
	TTree* outTree = new TTree("outTree", "outTree");
	bool isCal;
	bool isSoft;
	bool isShort;
	double spareChannelMaxSamp[4];
	double spareChannelMinSamp[4];
	double spareChannelRMS[4];
	int runNum;
	int unixTime;
	int unixTimeUs;
	int timeStamp;
	int eventNumber;
	bool hasDigitizerError=false;
	bool hasA3S4Issue;
	outTree->Branch("isCal", &isCal, "isCal/O");
	outTree->Branch("isSoft", &isSoft, "isSoft/O");
	outTree->Branch("isShort", &isShort, "isShort/O");
	outTree->Branch("unixTime",&unixTime, "unixTime/I");
	outTree->Branch("unixTimeUs",&unixTimeUs, "unixTimeUs/I");
	outTree->Branch("timeStamp",&timeStamp, "timeStamp/I");
	outTree->Branch("eventNumber",&eventNumber, "eventNumber/I");
	outTree->Branch("hasDigitizerError", &hasDigitizerError, "hasDigitizerError/O");
	outTree->Branch("spareChannelMaxSamp", &spareChannelMaxSamp, "spareChannelMaxSamp[4]/D");
	outTree->Branch("spareChannelMinSamp", &spareChannelMinSamp, "spareChannelMinSamp[4]/D");
	outTree->Branch("spareChannelRMS", &spareChannelRMS, "spareChannelRMS[4]/D");
	runNum=run;
	outTree->Branch("run",&runNum);

	Long64_t numEntries=eventTree->GetEntries();
	int start=0;
	for(Long64_t event=start;event<numEntries;event++) {
		cout<<"On event "<<event<<endl;
		eventTree->GetEntry(event);
		isCal = rawAtriEvPtr->isCalpulserEvent();
		isSoft = rawAtriEvPtr->isSoftwareTrigger();
		unixTime = (int)rawAtriEvPtr->unixTime;
		unixTimeUs= (int)rawAtriEvPtr->unixTimeUs;
		timeStamp = (int)rawAtriEvPtr->timeStamp;
		eventNumber = (int)rawAtriEvPtr->eventNumber;

		UsefulAtriStationEvent *ev = new UsefulAtriStationEvent(rawAtriEvPtr,AraCalType::kLatestCalib);
		// hasDigitizerError = !(qualCut->isGoodEvent(ev));

		vector<TGraph*> electChansGraphs;
		electChansGraphs.push_back(ev->getGraphFromElecChan(6));
		electChansGraphs.push_back(ev->getGraphFromElecChan(14));
		electChansGraphs.push_back(ev->getGraphFromElecChan(22));
		electChansGraphs.push_back(ev->getGraphFromElecChan(30));
		for(int i=0; i<4; i++){
			if(electChansGraphs[i]->GetN()<64) isShort=true;
			spareChannelMaxSamp[i]=electChansGraphs[i]->GetMaximum();
			spareChannelMinSamp[i]=electChansGraphs[i]->GetMinimum();
			spareChannelRMS[i]=electChansGraphs[i]->GetRMS(2);
		}
		outTree->Fill();
		deleteGraphVector(electChansGraphs);
		delete ev;
	} //loop over events
	
	fpOut->Write();
	fpOut->Close();
	
	fpIn->Close();
	delete fpIn;
}