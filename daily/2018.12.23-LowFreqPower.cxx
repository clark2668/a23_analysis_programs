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

	if(argc<3) {
		std::cout << "Usage\n" << argv[0] << " <station> <input_file> <output_location> "<<endl;
		return -1;
	}
	int station = atoi(argv[1]);

	/*
	arguments
	0: exec
	1: station
	3: input data file
	4: output location
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
	printf("Filter Run Number %d \n", run);

	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;

	char ped_file_name[400];
	sprintf(ped_file_name,"%s/run_specific_peds/A%d/all_peds/event%d_specificPeds.dat",PedDirPath,station,run);
	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	calibrator->setAtriPedFile(ped_file_name,station); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist

	AraQualCuts *qualCut = AraQualCuts::Instance(); //we also need a qual cuts tool

	char outfile_name[400];
	sprintf(outfile_name,"%s/low_freq_power_run%d.root",argv[3],run);

	TFile *fpOut = TFile::Open(outfile_name, "RECREATE");
	TTree* outTree = new TTree("outTree", "outTree");
	bool isCal;
	bool isSoft;
	int waveform_length[16];
	double frac_power_75[16];
	double frac_power_110[16];
	double frac_power_150[16];
	double frac_above_850[16];
	double frac_between[16];
	double power_75[16];
	double power_110[16];
	double power_150[16];
	double above_850[16];
	double between[16];
	int runNum;
	bool hasDigitizerError;
	outTree->Branch("isCal", &isCal, "isCal/O");
	outTree->Branch("isSoft", &isSoft, "isSoft/O");
	outTree->Branch("hasDigitizerError", &hasDigitizerError, "hasDigitizerError/O");
	outTree->Branch("waveform_length", &waveform_length, "waveform_length[16]/I");
	outTree->Branch("frac_power_75", &frac_power_75, "frac_power_75[16]/D");
	outTree->Branch("frac_power_110", &frac_power_110, "frac_power_110[16]/D");
	outTree->Branch("frac_power_150", &frac_power_150, "frac_power_150[16]/D");
	outTree->Branch("frac_above_850", &frac_above_850, "frac_above_850[16]/D");
	outTree->Branch("frac_between", &frac_between, "frac_between[16]/D");
	outTree->Branch("power_75", &power_75, "power_75[16]/D");
	outTree->Branch("power_110", &power_110, "power_110[16]/D");
	outTree->Branch("power_150", &power_150, "power_150[16]/D");
	outTree->Branch("above_850", &above_850, "above_850[16]/D");
	outTree->Branch("between", &between, "between[16]/D");
	runNum=run;
	outTree->Branch("run",&runNum);

	Long64_t numEntries=eventTree->GetEntries();
	// numEntries=20;

	for(Long64_t event=0;event<numEntries;event++) {
		eventTree->GetEntry(event);
		isCal = rawAtriEvPtr->isCalpulserEvent();
		isSoft = rawAtriEvPtr->isSoftwareTrigger();

		UsefulAtriStationEvent *ev = new UsefulAtriStationEvent(rawAtriEvPtr,AraCalType::kLatestCalib);
		hasDigitizerError = (qualCut->hasBlockGap(ev) || qualCut->hasTooFewBlocks(ev) || qualCut->hasTimingError(ev));
		// hasDigitizerError = !(qualCut->isGoodEvent(ev));
		if(hasDigitizerError){ //cleanup
			outTree->Fill();
			delete ev;
			continue;
		}

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

		vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(grWaveformsRaw, 0.5, xLabel, yLabel, titlesForGraphs);
		vector<TGraph*> grWaveformsPadded = makePaddedGraphs(grWaveformsInt, 0, xLabel, yLabel, titlesForGraphs);
		vector<TGraph*> grSpectra = makePowerSpectrumGraphs(grWaveformsPadded, xLabel, yLabel, titlesForGraphs);

		for(int i=0; i<16; i++){
			waveform_length[i]=grWaveformsRaw[i]->GetN();
			frac_power_75[i]=cumulativePowerBelowfromSpectrum(grSpectra[i],75.,power_75[i]);
			frac_power_110[i]=cumulativePowerBelowfromSpectrum(grSpectra[i],110.,power_110[i]);
			frac_power_150[i]=cumulativePowerBelowfromSpectrum(grSpectra[i],150.,power_150[i]);
			frac_above_850[i]=cumulativePowerAbovefromSpectrum(grSpectra[i],850.,above_850[i]);
			frac_between[i]=cumulativePowerBetweenfromSpectrum(grSpectra[i],75.,850.,between[i]);
		}
		outTree->Fill();

		//cleanup
		deleteGraphVector(grWaveformsRaw);
		deleteGraphVector(grWaveformsInt);
		deleteGraphVector(grWaveformsPadded);
		deleteGraphVector(grSpectra);
		delete ev;
	} //loop over events
	
	fpOut->Write();
	fpOut->Close();
	
	fpIn->Close();
	delete fpIn;
}