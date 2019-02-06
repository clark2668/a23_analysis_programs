////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	v2_analysis_run_summary.cxx 
////	A23 diffuse, make a run summary
////
////	Nov 2018
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>

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

#include "tools_inputParameters.h"
#include "tools_outputObjects.h"
#include "tools_runSummaryObjects.h"
#include "tools_WaveformFns.h"
#include "tools_PlottingFns.h"
#include "tools_Cuts.h"

int main(int argc, char **argv)
{

	if(argc<5) {
		std::cout << "Usage\n" << argv[0] << " <simulation_flag> <station> <output directory> <input file> <pedestal file> \n";
		return -1;
	}
	isSimulation=atoi(argv[1]);
	int station_num=atoi(argv[2]);

	numEvents = 0;
	numSoftTriggers = 0;
	numCalpulsers = 0;
	numRFTriggers = 0;

	vector<double> RMS_SoftTrigger_total;
	RMS_SoftTrigger_total.resize(nGraphs);
	vector<double> RMS_Calpulser_total;
	RMS_Calpulser_total.resize(nGraphs);
	vector<double> RMS_RFTrigger_total;
	RMS_RFTrigger_total.resize(nGraphs);
	vector<double> RMS_All_total;
	RMS_All_total.resize(nGraphs);

	vector<TGraph*> PowerSpectrumAverage_RFTrigger;
	vector<TGraph*> PowerSpectrumAverage_Calpulser;
	vector<TGraph*> PowerSpectrumAverage_SoftTrigger;
	PowerSpectrumAverage_RFTrigger.resize(nGraphs);
	PowerSpectrumAverage_Calpulser.resize(nGraphs);
	PowerSpectrumAverage_SoftTrigger.resize(nGraphs);

	for (int i = 0; i < nGraphs; i++){
		PowerSpectrumAverage_RFTrigger[i] = new TGraph(1024);
		PowerSpectrumAverage_Calpulser[i] = new TGraph(1024);
		PowerSpectrumAverage_SoftTrigger[i] = new TGraph(1024);
	}

	stringstream ss;
	string xLabel, yLabel;
	vector<string> titlesForGraphs;
	for (int i = 0; i < nGraphs; i++){
		ss.str("");
		ss << "Channel " << i;
		titlesForGraphs.push_back(ss.str());
	}

	TFile *fp = TFile::Open(argv[4]);
	if(!fp) {
		std::cerr << "Can't open file\n";
		return -1;
	}
	TTree *eventTree; 
		eventTree= (TTree*) fp->Get("eventTree");
		if(!eventTree) {
		std::cerr << "Can't find eventTree\n";
		return -1;
	}
	AraEventCalibrator * calibrator = AraEventCalibrator::Instance();
	if(argc==6){
		calibrator->setAtriPedFile(argv[5], station_num);
	}
	else{
		calibrator->setAtriPedFile("", station_num);
	}

	if(isSimulation){
		eventTree->SetBranchAddress("UsefulAtriStationEvent", &realAtriEvPtr);
		printf("Simulation; load useful event tree straight away \n");
	}
	else{
		eventTree->SetBranchAddress("event",&rawAtriEvPtr);
		printf("Data; load raw event tree \n");
	}

	Long64_t numEntries=eventTree->GetEntries();
	Long64_t starEvery=numEntries/100;
	if(starEvery==0) starEvery++;

	int runNum = getrunNum(argv[4]);
	printf("Run Number %d \n", runNum);
	runNumber = runNum;

	string runSummaryFilename = getRunSummaryFilename(station_num, argv[3], argv[4]);
	TFile *OutputFile = TFile::Open(runSummaryFilename.c_str(), "RECREATE");
	TTree* OutputTree=new TTree("SummaryTree", "SummaryTree");

	OutputTree->Branch("runNumber", &runNumber, "runNumber/I");
	OutputTree->Branch("numEvents", &numEvents, "numEvents/I");
	OutputTree->Branch("numRFTriggers", &numRFTriggers, "numRFTriggers/I");
	OutputTree->Branch("numSoftTriggers", &numSoftTriggers, "numSoftTriggers/I");
	OutputTree->Branch("numCalpulsers", &numCalpulsers, "numCalpulsers/I");

	OutputTree->Branch("RMS_SoftTrigger", &RMS_SoftTrigger, "RMS_SoftTrigger[20]/D");
	OutputTree->Branch("RMS_RFTrigger", &RMS_RFTrigger, "RMS_RFTrigger[20]/D");
	OutputTree->Branch("RMS_Calpulser", &RMS_Calpulser, "RMS_Calpulser[20]/D");
	OutputTree->Branch("RMS_All", &RMS_All, "RMS_All[20]/D");

	for (int i = 0; i < 16; i++){
		ss.str("");	ss << "PowerSpectrumAverage_RFTrigger_" << i;
		OutputTree->Branch(ss.str().c_str(), &PowerSpectrumAverage_RFTrigger[i]);
		ss.str("");	ss << "PowerSpectrumAverage_Calpulser_" << i;
		OutputTree->Branch(ss.str().c_str(), &PowerSpectrumAverage_Calpulser[i]);
		ss.str("");	ss << "PowerSpectrumAverage_SoftTrigger_" << i;
		OutputTree->Branch(ss.str().c_str(), &PowerSpectrumAverage_SoftTrigger[i]);
	}
	for(Long64_t event=0;event<numEntries;event++) {
		if(event%starEvery==0) {
			std::cerr << "*";       
		}

		eventTree->GetEntry(event);

		if (isSimulation == false){
			realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
		}

		if (isSimulation){
			isCalpulser = false;
			isSoftTrigger = false;
		} else{
			isCalpulser = rawAtriEvPtr->isCalpulserEvent();
			isSoftTrigger = rawAtriEvPtr->isSoftwareTrigger();
		}

		xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
		vector<TGraph*> grWaveformsRaw = makeGraphsFromRF(realAtriEvPtr, nGraphs, xLabel, yLabel, titlesForGraphs);

		hasDigitizerError = hasDigitizerIssue(grWaveformsRaw); 
		//if the event has a  digitizer error, skip it
		//note that exiting this far up will prevent any errors with the averaging
		//because we don't count contributions to the average until the numEvents++ later
		if(hasDigitizerError){
			deleteGraphVector(grWaveformsRaw); //cleanup
			continue; //skip this event
		}
	   
		ss.str("");

		xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
		vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(grWaveformsRaw, interpolationTimeStep, xLabel, yLabel, titlesForGraphs);
		
		ss.str("");

		vector<double> vWaveformRMS;
		getRMS(grWaveformsInt, vWaveformRMS, 0);

		if (isSoftTrigger == true){
			transform(RMS_SoftTrigger_total.begin(), RMS_SoftTrigger_total.end(), vWaveformRMS.begin(),
			RMS_SoftTrigger_total.begin(), std::plus<double>());
			numSoftTriggers++;
		}
		if (isCalpulser == true){
			transform(RMS_Calpulser_total.begin(), RMS_Calpulser_total.end(), vWaveformRMS.begin(),
			RMS_Calpulser_total.begin(), std::plus<double>());
			numCalpulsers++;
		}
		if (isCalpulser == false && isSoftTrigger == false){
			transform(RMS_RFTrigger_total.begin(), RMS_RFTrigger_total.end(), vWaveformRMS.begin(),
			RMS_RFTrigger_total.begin(), std::plus<double>());
			numRFTriggers++;
		}
		transform(RMS_All_total.begin(), RMS_All_total.end(), vWaveformRMS.begin(),
		RMS_All_total.begin(), std::plus<double>());
		numEvents++;

		vector<TGraph*> grWaveformsPadded = makePaddedGraphs(grWaveformsInt, isSoftTrigger, xLabel, yLabel, titlesForGraphs);
		vector<TGraph*> grWaveformsPowerSpectrum = makePowerSpectrumGraphs(grWaveformsPadded, xLabel, yLabel, titlesForGraphs);

		for (int i = 0; i < nGraphs; i++){
			for (int point = 0; point < grWaveformsPowerSpectrum[i]->GetN(); point++){
				double x1,y1;
				grWaveformsPowerSpectrum[i]->GetPoint(point, x1,y1);
				double x2,y2;
				if (isCalpulser == true){
					PowerSpectrumAverage_Calpulser[i]->GetPoint(point,x2,y2);
					PowerSpectrumAverage_Calpulser[i]->SetPoint(point,x1,y1+y2);
				}
				if (isSoftTrigger == true){
					PowerSpectrumAverage_SoftTrigger[i]->GetPoint(point,x2,y2);
					PowerSpectrumAverage_SoftTrigger[i]->SetPoint(point,x1,y1+y2);
				}
				if (isSoftTrigger == false && isCalpulser == false){
					PowerSpectrumAverage_RFTrigger[i]->GetPoint(point,x2,y2);
					PowerSpectrumAverage_RFTrigger[i]->SetPoint(point,x1,y1+y2);
				}
			}
		}

		deleteGraphVector(grWaveformsPowerSpectrum);
		deleteGraphVector(grWaveformsPadded);
		deleteGraphVector(grWaveformsInt);
		deleteGraphVector(grWaveformsRaw);
		if (isSimulation == false) {
			delete realAtriEvPtr;
		}	 
	}

	cout << numSoftTriggers << " : ";
	for (int i = 0; i < nGraphs; i++){
		RMS_RFTrigger_total[i] = RMS_RFTrigger_total[i]/(double)numRFTriggers;
		cout << RMS_SoftTrigger_total[i] << " : "; 
		RMS_SoftTrigger_total[i] = RMS_SoftTrigger_total[i]/(double)numSoftTriggers;
		RMS_Calpulser_total[i] = RMS_Calpulser_total[i]/(double)numCalpulsers;
		RMS_All_total[i] = RMS_All_total[i]/(double)numEvents;
	}
	cout << endl;

	copy(RMS_All_total.begin(), RMS_All_total.begin()+16, RMS_All);
	copy(RMS_RFTrigger_total.begin(), RMS_RFTrigger_total.begin()+16, RMS_RFTrigger);
	copy(RMS_SoftTrigger_total.begin(), RMS_SoftTrigger_total.begin()+16, RMS_SoftTrigger);
	copy(RMS_Calpulser_total.begin(), RMS_Calpulser_total.begin()+16, RMS_Calpulser);

	for (int i = 0; i < nGraphs; i++){
		for (int point = 0; point < PowerSpectrumAverage_Calpulser[i]->GetN(); point++){
			double x,y;
			PowerSpectrumAverage_Calpulser[i]->GetPoint(point,x,y);
			if (numCalpulsers != 0){
				PowerSpectrumAverage_Calpulser[i]->SetPoint(point,x,y/(double)numCalpulsers);
			} else {
				PowerSpectrumAverage_Calpulser[i]->SetPoint(point,x,0);
			}


			PowerSpectrumAverage_SoftTrigger[i]->GetPoint(point,x,y);
			if (numSoftTriggers != 0){
				PowerSpectrumAverage_SoftTrigger[i]->SetPoint(point,x,y/(double)numSoftTriggers);
			} else {
				PowerSpectrumAverage_SoftTrigger[i]->SetPoint(point,x,0);
			}

			PowerSpectrumAverage_RFTrigger[i]->GetPoint(point,x,y);
			if (numRFTriggers != 0){
				PowerSpectrumAverage_RFTrigger[i]->SetPoint(point,x,y/(double)numRFTriggers);
			} else{
				PowerSpectrumAverage_RFTrigger[i]->SetPoint(point,x,0);
			}
		}
	}

	OutputTree->Fill();
	OutputFile->Write();
	OutputFile->Close();
	fp->Close();
	delete fp;
	printf("Done! Run Number %d \n", runNum);
}