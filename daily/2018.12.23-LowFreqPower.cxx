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
#include "FFTtools.h"

#include "tools_PlottingFns.h"
#include "tools_WaveformFns.h"
#include "tools_Cuts.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"

using namespace std;

double cumulativePowerBelow(TGraph *grIn, double freq);


int main(int argc, char **argv)
{

	if(argc<3) {
		std::cout << "Usage\n" << argv[0] << " <station> <year> <input_file> <output_location> "<<endl;
		return -1;
	}
	int station = atoi(argv[1]);
	int year = atoi(argv[2]);

	/*
	arguments
	0: exec
	1: station
	2: year
	3: input data file
	4: output location
	*/
	
	TFile *fpIn = TFile::Open(argv[3]);
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
	printf("Filter Run Number %d \n", run);

	char ped_file_name[400];
	if(year==2013){
		sprintf(ped_file_name,"/fs/scratch/PAS0654/ara/peds/run_specific_peds/A%d/%d/event%d_specificPeds.dat",station,year,run);
	}
	else if(year==2014 || year==2015 || year==2016){
		sprintf(ped_file_name,"/fs/scratch/PAS0654/ara/peds/run_specific_peds/A%d/%d/event00%d_specificPeds.dat",station,year,run);
	}
	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	calibrator->setAtriPedFile(ped_file_name,station); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist

	char outfile_name[400];
	sprintf(outfile_name,"%s/low_freq_power_run%d.root",argv[4],run);

	TFile *fpOut = TFile::Open(outfile_name, "RECREATE");
	TTree* outTree = new TTree("outTree", "outTree");
	bool isCal;
	bool isSoft;
	int waveform_length[16];
	double frac_power[16];
	int runNum;
	bool hasDigitizerError;
	outTree->Branch("isCal", &isCal, "isCal/O");
	outTree->Branch("isSoft", &isSoft, "isSoft/O");
	outTree->Branch("hasDigitizerError", &hasDigitizerError, "hasDigitizerError/O");
	outTree->Branch("waveform_length", &waveform_length, "waveform_length[16]/I");
	outTree->Branch("frac_power", &frac_power, "frac_power[16]/D");
	runNum=run;
	outTree->Branch("run",&runNum);

	Long64_t numEntries=eventTree->GetEntries();

	for(Long64_t event=0;event<numEntries;event++) {
		eventTree->GetEntry(event);
		isCal = rawAtriEvPtr->isCalpulserEvent();
		isSoft = rawAtriEvPtr->isSoftwareTrigger();

		UsefulAtriStationEvent *ev = new UsefulAtriStationEvent(rawAtriEvPtr,AraCalType::kLatestCalib);

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
		vector <TGraph*> grWaveformsRaw = makeGraphsFromRF(ev,16,xLabel,yLabel,titlesForGraphs);
		hasDigitizerError = hasDigitizerIssue(grWaveformsRaw);
		if(hasDigitizerError){ //cleanup
			outTree->Fill();
			deleteGraphVector(grWaveformsRaw); //cleanup
			delete ev;
			continue;
		}
		for(int i=0; i<16; i++){
			waveform_length[i]=grWaveformsRaw[i]->GetN();
			frac_power[i]=cumulativePowerBelow(grWaveformsRaw[i],75.);
		}
		outTree->Fill();

		//cleanup
		deleteGraphVector(grWaveformsRaw); //cleanup
		delete ev;
	} //loop over events
	
	fpOut->Write();
	fpOut->Close();
	
	fpIn->Close();
	delete fpIn;
}