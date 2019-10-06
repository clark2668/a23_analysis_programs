////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  2019.09.25-TraceLength.cxx 
////  save out the number of samples and trace length in ns
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

double MaxMeanBlock(TGraph *grIn);

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

	char outfile_name[400];
	sprintf(outfile_name,"%s/tracelength_run%d.root",argv[3],run);

	TFile *fpOut = TFile::Open(outfile_name, "RECREATE");
	TTree* outTree = new TTree("outTree", "outTree");
	bool isCal;
	bool isSoft;
	bool isShort;
	int runNum;
	int unixTime;
	int unixTimeUs;
	int timeStamp;
	int numSamples[16];
	double traceLength[16];

	runNum=run;

	outTree->Branch("isCal", &isCal, "isCal/O");
	outTree->Branch("isSoft", &isSoft, "isSoft/O");
	outTree->Branch("isShort", &isShort, "isShort/O");
	outTree->Branch("unixTime",&unixTime, "unixTime/I");
	outTree->Branch("unixTime",&unixTime, "unixTime/I");
	outTree->Branch("unixTimeUs",&unixTimeUs, "unixTimeUs/I");
	outTree->Branch("timeStamp",&timeStamp, "timeStamp/I");
	outTree->Branch("traceLength", &numSamples, "numSamples[16]/I");
	outTree->Branch("numSamples", &traceLength, "traceLength[16]/D");
	outTree->Branch("run",&runNum);

	//now get the waveforms
	stringstream ss1;
	string xLabel, yLabel;
	xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
	vector<string> titlesForGraphs;
	for (int i = 0; i < 16; i++){
		ss1.str("");
		ss1 << "Channel " << i;
		titlesForGraphs.push_back(ss1.str());
	}

	Long64_t numEntries=eventTree->GetEntries();
	int start=0;
	int stop = numEntries;
	stop=10;
	for(Long64_t event=start;event<numEntries;event++) {
		eventTree->GetEntry(event);
		isCal = rawAtriEvPtr->isCalpulserEvent();
		isSoft = rawAtriEvPtr->isSoftwareTrigger();
		unixTime = (int)rawAtriEvPtr->unixTime;
		unixTimeUs= (int)rawAtriEvPtr->unixTimeUs;
		timeStamp = (int)rawAtriEvPtr->timeStamp;

		UsefulAtriStationEvent *ev = new UsefulAtriStationEvent(rawAtriEvPtr,AraCalType::kLatestCalib);
		vector <TGraph*> grWaveformsRaw = makeGraphsFromRF(ev,16,xLabel,yLabel,titlesForGraphs);

		for(int i=0; i<16; i++){
			int this_numSamples = grWaveformsRaw[i]->GetN();
			numSamples[i] = this_numSamples;
			traceLength[i] = grWaveformsRaw[i]->GetX()[this_numSamples-1] - grWaveformsRaw[i]->GetX()[0];
		}
		for(int i=0; i<16; i++){
			printf("Cal %d, Soft %d, Event %d: Chan %d, Samps %d, Length %.2f \n", isCal, isSoft, event, i, numSamples[i], traceLength[i]);
		}
		outTree->Fill();
		deleteGraphVector(grWaveformsRaw);
		delete ev;
	} //loop over events
	
	fpOut->Write();
	fpOut->Close();
	
	fpIn->Close();
	delete fpIn;
}