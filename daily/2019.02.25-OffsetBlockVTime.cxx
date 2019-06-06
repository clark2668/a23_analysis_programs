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

	AraQualCuts *qualCut = AraQualCuts::Instance(); //we also need a qual cuts tool

	char outfile_name[400];
	sprintf(outfile_name,"%s/offset_block_run%d.root",argv[3],run);

	TFile *fpOut = TFile::Open(outfile_name, "RECREATE");
	TTree* outTree = new TTree("outTree", "outTree");
	bool isCal;
	bool isSoft;
	bool isShort;
	double first_block_mean[16];
	double maxMeans[16];
	double minMeans[16];
	double absMaxBlockMeans[16];
	double RMS[16];
	int runNum;
	int unixTime;
	int unixTimeUs;
	int timeStamp;
	bool hasDigitizerError;
	bool hasA3S4Issue;
	outTree->Branch("isCal", &isCal, "isCal/O");
	outTree->Branch("isSoft", &isSoft, "isSoft/O");
	outTree->Branch("isShort", &isShort, "isShort/O");
	outTree->Branch("unixTime",&unixTime, "unixTime/I");
	outTree->Branch("unixTime",&unixTime, "unixTime/I");
	outTree->Branch("unixTimeUs",&unixTimeUs, "unixTimeUs/I");
	outTree->Branch("timeStamp",&timeStamp, "timeStamp/I");
	outTree->Branch("hasDigitizerError", &hasDigitizerError, "hasDigitizerError/O");
	outTree->Branch("first_block_mean", &first_block_mean, "first_block_mean[16]/D");
	outTree->Branch("maxMeans", &maxMeans, "maxMeans[16]/D");
	outTree->Branch("minMeans", &minMeans, "minMeans[16]/D");
	outTree->Branch("absMaxBlockMeans", &absMaxBlockMeans, "absMaxBlockMeans[16]/D");
	outTree->Branch("RMS", &RMS, "RMS[16]/D");
	outTree->Branch("hasA3S4Issue", &hasA3S4Issue, "hasA3S4Issue/O");
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

		UsefulAtriStationEvent *ev = new UsefulAtriStationEvent(rawAtriEvPtr,AraCalType::kLatestCalib);
		hasDigitizerError = !(qualCut->isGoodEvent(ev));
		if(hasDigitizerError){ //cleanup
			outTree->Fill();
			delete ev;
			continue;
		}

		// hasA3S4Issue=qualCut->hasA3S4OffsetBlocks(ev);

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
		// for(int chan=0; chan<16; chan++){
		// 	double first_block_mean_this=0.;
		// 	int num_samps_in_first_block=0;
		// 	double *oldX = grWaveformsRaw[chan]->GetX();
		// 	double *oldY = grWaveformsRaw[chan]->GetY();
		// 	double first_time = oldX[0];
		// 	for(int samp=0; samp<grWaveformsRaw[chan]->GetN(); samp++){
		// 		if(oldX[samp]<first_time+20.){
		// 			first_block_mean_this+=oldY[samp];
		// 			num_samps_in_first_block++;
		// 		}
		// 	}
		// 	first_block_mean[chan]=first_block_mean_this/double(num_samps_in_first_block);
		// }
		vector<TGraph*> rollingMeans;
		rollingMeans.resize(16);
		for(int i=0; i<16; i++){
			if(grWaveformsRaw[i]->GetN()<64){
		 		isShort=true;
		 		continue;
		 	}
		 	TGraph *grInt = FFTtools::getInterpolatedGraph(grWaveformsRaw[i],0.5);
		 	if(grInt->GetN()<101){
		 		isShort=true;
		 		delete grInt;
		 		continue;
		 	}
		 	// rollingMeans.push_back(qualCut->getRollingMean(grInt,160));
		 	rollingMeans[i]=(qualCut->getRollingMean(grInt,100));
		 	maxMeans[i] = abs(TMath::MaxElement(rollingMeans[i]->GetN(), rollingMeans[i]->GetY()));
		 	minMeans[i] = abs(TMath::MinElement(rollingMeans[i]->GetN(), rollingMeans[i]->GetY()));
		}		
		deleteGraphVector(rollingMeans);

		for(int chan=0; chan<16; chan++){
		 	absMaxBlockMeans[chan]=MaxMeanBlock(grWaveformsRaw[chan]);
		}
		//for(int chan=0; chan<16; chan++){
		//	RMS[chan] = grWaveformsRaw[chan]->GetRMS(2);
		//}

		outTree->Fill();

		//cleanup
		deleteGraphVector(grWaveformsRaw);
		delete ev;
	} //loop over events
	
	fpOut->Write();
	fpOut->Close();
	
	fpIn->Close();
	delete fpIn;
}

double MaxMeanBlock(TGraph *grIn){

	int n_input = grIn->GetN();
	double *oldX = grIn->GetX();
	double *oldY = grIn->GetY();

	deque <double> inX;
	deque <double> inY;	

	for(int samp=0; samp<n_input; samp++){
		inX.push_back(oldX[samp]);
		inY.push_back(oldY[samp]);
	}

	vector <double> means;
	while(inX.size()>0){
		vector <double> sub_X;
		vector <double> sub_Y;
		double first_time = inX[0];
		int num_to_pop=0;
		double mean_this_section=0.;
		for(int samp=0; samp<inX.size(); samp++){
			if(inX[samp]<=first_time+50.){
				sub_X.push_back(inX[samp]);
				sub_Y.push_back(inY[samp]);
				num_to_pop++;
				mean_this_section+=inY[samp];
			}
		}
		mean_this_section/=double(sub_X.size());
		means.push_back(abs(mean_this_section));
		for(int samp=0; samp<sub_X.size(); samp++){ //now fix the mean
			sub_Y[samp]=mean_this_section;
		}
		for(int iPop=0; iPop<num_to_pop; iPop++){
			inX.pop_front();
			inY.pop_front();
		}
	}

	std::sort(means.begin(), means.end(), std::greater<double>());
	return means[0];
}
