////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  2019.11.05-PowerContent.cxx 
////  work on power ratio cut
////
////  Nov 2019
////  trying to work out the ratio of power in strings relative to eachother
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
	// RawAtriStationEvent *rawAtriEvPtr=0;
	// eventTree->SetBranchAddress("event",&rawAtriEvPtr);
	UsefulAtriStationEvent *ev;
	eventTree->SetBranchAddress("UsefulAtriStationEvent",&ev);
	double weight;
	eventTree->SetBranchAddress("weight",&weight);
	int run;
	// eventTree->SetBranchAddress("run",&run);

	// for simulation only (stupid, stupid)
	run=getrunNum(argv[2]);
	// eventTree->GetEntry(0);
	printf("Starting Run %d \n", run);

	// char ped_file_name[400];
	// sprintf(ped_file_name,"/fs/project/PAS0654/ARA_DATA/A23/peds/run_specific_peds/A%d/all_peds/event%d_specificPeds.dat",station,run);
	// AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	// calibrator->setAtriPedFile(ped_file_name,station); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist

	// AraQualCuts *qualCut = AraQualCuts::Instance(); //we also need a qual cuts tool
	// vector<int> BadRunList=BuildBadRunList(station);

	char outfile_name[400];
	sprintf(outfile_name,"%s/diagnostic_info_run%d.root",argv[3],run);

	TFile *fpOut = TFile::Open(outfile_name, "RECREATE");
	TTree* outTree = new TTree("outTree", "outTree");
	bool isCal;
	bool isSoft;
	bool isShort;
	bool isSuperShort;
	bool hasSpareChannelIssuev1;
	bool hasSpareChannelIssuev2;
	bool hasDigitizerError;
	bool isKnownBadLivetime;
	outTree->Branch("isCal", &isCal, "isCal/O");
	outTree->Branch("isSoft", &isSoft, "isSoft/O");
	outTree->Branch("isShort", &isShort, "isShort/O");
	outTree->Branch("isSuperShort", &isSuperShort, "isSuperShort/O");
	outTree->Branch("hasSpareChannelIssuev1", &hasSpareChannelIssuev1, "hasSpareChannelIssuev1/O");
	outTree->Branch("hasSpareChannelIssuev2", &hasSpareChannelIssuev2, "hasSpareChannelIssuev2/O");
	outTree->Branch("hasDigitizerError", &hasDigitizerError, "hasDigitizerError/O");
	outTree->Branch("isKnownBadLivetime", &isKnownBadLivetime, "isKnownBadLivetime/O");

	int unixTime;
	int unixTimeUs;
	int timeStamp;
	int eventNumber;
	double deepChannelRMS[16];
	double spareChannelRMS[4];
	outTree->Branch("unixTime",&unixTime, "unixTime/I");
	outTree->Branch("unixTimeUs",&unixTimeUs, "unixTimeUs/I");
	outTree->Branch("timeStamp",&timeStamp, "timeStamp/I");
	outTree->Branch("eventNumber",&eventNumber, "eventNumber/I");
	outTree->Branch("deepChannelRMS",&deepChannelRMS,"deepChannelRMS[16]/D");
	outTree->Branch("spareChannelRMS", &spareChannelRMS, "spareChannelRMS[4]/D");

	double powerPerChannel[16];
	double powerPerString[4];
	double powerRatio;
	outTree->Branch("powerPerChannel",&powerPerChannel,"powerPerChannel[16]/D");
	outTree->Branch("powerPerString", &powerPerString, "powerPerString[4]/D");	
	outTree->Branch("powerRatio", &powerRatio, "powerRatio/D");	

	int runNum;	
	bool isKnownBadRun;
	outTree->Branch("runNum",&runNum,"runNum/I");
	outTree->Branch("isKnownBadRun", &isKnownBadRun, "isKnownBadRun/O");

	double weight_out;
	outTree->Branch("weight_out",&weight_out,"weight_out/D");

	runNum=run;
	// isKnownBadRun = isBadRun(station,runNum,BadRunList);
	isKnownBadRun=false;

	Long64_t numEntries=eventTree->GetEntries();
	int start=0;
	for(Long64_t event=start;event<numEntries;event++) {
		// cout<<"Event is "<<event<<endl;
		eventTree->GetEntry(event);
		
		isCal=false;
		isSoft=false;
		unixTime=0;
		unixTimeUs=0;
		timeStamp=0;
		weight_out=weight;

		// isCal = rawAtriEvPtr->isCalpulserEvent();
		// isSoft = rawAtriEvPtr->isSoftwareTrigger();
		// unixTime = (int)rawAtriEvPtr->unixTime;
		// unixTimeUs= (int)rawAtriEvPtr->unixTimeUs;
		// timeStamp = (int)rawAtriEvPtr->timeStamp;
		// eventNumber = (int)rawAtriEvPtr->eventNumber;
		eventNumber=event;
		// isKnownBadLivetime = isBadLivetime(station,unixTime);

		// get info on the deep (rf) channels
		// UsefulAtriStationEvent *ev = new UsefulAtriStationEvent(rawAtriEvPtr,AraCalType::kLatestCalib);
		// hasDigitizerError = !(qualCut->isGoodEvent(ev));
		hasDigitizerError=false;

		vector<TGraph*> rfChanGraphs;
		for(int i=0; i<16; i++){
			rfChanGraphs.push_back(ev->getGraphFromRFChan(i));
			deepChannelRMS[i] = rfChanGraphs[i]->GetRMS(2);
			// cout<<"RMS for chan "<<i<<" is "<<rfChanGraphs[i]->GetRMS(2)<<endl;
			if(rfChanGraphs[i]->GetN()<500){
				isShort=true;
				if(rfChanGraphs[i]->GetN()<64){
					isSuperShort=true;
				}
			}
		}

		// vector<TGraph*> spareElecChanGraphs;
		// spareElecChanGraphs.push_back(ev->getGraphFromElecChan(6));
		// spareElecChanGraphs.push_back(ev->getGraphFromElecChan(14));
		// spareElecChanGraphs.push_back(ev->getGraphFromElecChan(22));
		// spareElecChanGraphs.push_back(ev->getGraphFromElecChan(30));
		// for(int i=0; i<4; i++){
		// 	spareChannelRMS[i]=spareElecChanGraphs[i]->GetRMS(2);
		// }
		// hasSpareChannelIssuev1 = hasSpareChannelIssue(spareElecChanGraphs);
		// hasSpareChannelIssuev2 = hasSpareChannelIssue_v2(spareElecChanGraphs, station);
		hasSpareChannelIssuev1=false;
		hasSpareChannelIssuev2=false;


		for(int chan=0; chan<16; chan++){
			powerPerChannel[chan]=0.; // set this, or god knows what will happen
			Double_t *yVals = rfChanGraphs[chan]->GetY();
			int N = rfChanGraphs[chan]->GetN();
			for(int samp=0; samp<N; samp++) powerPerChannel[chan]+=yVals[samp]*yVals[samp];
		}

		vector<double> powerPerString_local;
		for(int string=0; string<4; string++){
			powerPerString[string] = 0.;
			double thisPowerPerString = powerPerChannel[string]+powerPerChannel[string+4]+powerPerChannel[string+8]+powerPerChannel[string+12];
			powerPerString[string] = thisPowerPerString;
			powerPerString_local.push_back(thisPowerPerString);
		}
		std::sort(powerPerString_local.begin(), powerPerString_local.end()); //sort smallest to largest
		powerRatio=powerPerString_local[3]/powerPerString_local[0];
		cout<<"power ratio is "<<powerRatio<<endl;

		// fill output
		outTree->Fill();

		// clean-up
		// deleteGraphVector(spareElecChanGraphs);
		deleteGraphVector(rfChanGraphs);

		// delete ev;
	} //loop over events
	
	fpOut->Write();
	fpOut->Close();
	
	fpIn->Close();
	delete fpIn;
	printf("Done Run %d \n", run);


}