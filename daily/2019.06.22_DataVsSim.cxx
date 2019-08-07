////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  v2_analysis_filter.cxx 
////  A23 diffuse, filter events
////
////  Nov 2018
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <iomanip>
#include <sstream>

//AraRoot Includes
#include "RawIcrrStationEvent.h"
#include "RawAtriStationEvent.h"
#include "UsefulAraStationEvent.h"
#include "UsefulIcrrStationEvent.h"
#include "UsefulAtriStationEvent.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"

RawAtriStationEvent *rawAtriEvPtr;
UsefulAtriStationEvent *realAtriEvPtr;

#include "Event.h"
#include "Detector.h"
#include "Report.h"

#include "AraGeomTool.h"
#include "AraQualCuts.h"
#include "AraAntennaInfo.h"

#include "tools_inputParameters.h"
#include "tools_outputObjects.h"
#include "tools_runSummaryObjects.h"
#include "tools_WaveformFns.h"
#include "tools_PlottingFns.h"
#include "tools_Constants.h"
#include "tools_RecoFns.h"
#include "tools_filterEvent.h"
#include "tools_Cuts.h"

int main(int argc, char **argv)
{

	if(argc<4) {
		std::cout << "Usage\n" << argv[0] << " <simulation_flag> <station> <input_file> \n";
		return -1;
	}

	isSimulation=atoi(argv[1]);
	int station_num=atoi(argv[2]);
		
	stringstream ss;
	string xLabel, yLabel;
	vector<string> titlesForGraphs;
	for (int i = 0; i < nGraphs; i++){
		ss.str("");
		ss << "Channel " << i;
		titlesForGraphs.push_back(ss.str());
	}
	
	AraGeomTool * geomTool = new AraGeomTool();
	AraQualCuts *qualCut = AraQualCuts::Instance();
	
	TFile *fp = TFile::Open(argv[3]);
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
	
	double weight;
	int unixTime;
	int unixTimeUs;
	int eventNumber;
	double maxPeakVfromSim;
	double PeakVfromSim[16][2];
	int Trig_Pass[16] = {0};

	if(isSimulation){
		eventTree->SetBranchAddress("UsefulAtriStationEvent", &realAtriEvPtr);
		eventTree->SetBranchAddress("weight", &weight);
		printf("Simulation; load useful event tree straight away \n");
	}
	else{
		eventTree->SetBranchAddress("event",&rawAtriEvPtr);
		printf("Data; load raw event tree \n");
	}
	
	Long64_t numEntries=eventTree->GetEntries();
	Long64_t starEvery=numEntries/80;
	if(starEvery==0) starEvery++;
		
	int runNum = getrunNum(argv[6]);
	printf("Filter Run Number %d \n", runNum);
	printf("Num entries is %d \n", numEntries);
	numEntries=5;

	for(Long64_t event=0;event<numEntries;event++) {

		if(event%starEvery==0) {
			std::cout << "*";       
		}
		
		eventTree->GetEntry(event);
		if (isSimulation == false){
			unixTime=(int)rawAtriEvPtr->unixTime;
			unixTimeUs=(int)rawAtriEvPtr->unixTimeUs;
			eventNumber=(int)rawAtriEvPtr->eventNumber;
		} else{
			eventNumber = event;
		}

		if(!isSimulation){
			realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
		}
			
		if (isSimulation){
			isCalpulser = false;
			isSoftTrigger = false;
		} else{
			isCalpulser = rawAtriEvPtr->isCalpulserEvent();
			isSoftTrigger = rawAtriEvPtr->isSoftwareTrigger();
			weight = 1.;
		}

		bool analyzeEvent = false;
		if(!isSoftTrigger) analyzeEvent=true;

		if (analyzeEvent == true){

			weight_out = weight;
			if(!isSimulation)
				hasDigitizerError = !(qualCut->isGoodEvent(realAtriEvPtr));
			else
				hasDigitizerError=false;
			xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
			vector<TGraph*> grWaveformsRaw = makeGraphsFromRF(realAtriEvPtr, nGraphs, xLabel, yLabel, titlesForGraphs);
			ss.str("");

			for (int i = 0; i < 16; i++){
				waveformLength[i] = grWaveformsRaw[i]->GetN();
				printf("Graph %d length is %d\n", i, waveformLength[i]);
			}
	
			xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
			vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(grWaveformsRaw, interpolationTimeStep, xLabel, yLabel, titlesForGraphs);
	
			for(int i=0; i<16; i++){
				printf("Interpolated graph %d length is %d \n", i, grWaveformsInt[i]->GetN());
			}

			deleteGraphVector(grWaveformsInt);
			deleteGraphVector(grWaveformsRaw);
			if (isSimulation == false) {
				delete realAtriEvPtr;
			}
		} //analyze event?
	} //loop over events
	
	fp->Close();
	delete fp;
	printf("Done! Run Number %d \n", runNum);

}