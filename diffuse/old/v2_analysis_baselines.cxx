////////////////////////////////////////////////////////////////////////////////
////	v2_analysis_baselines.cxx 
////	A23 diffuse, make FFT baselines for all runs
////
////	Nov 2018
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <sstream>
#include <string>

//AraRoot Includes
#include "RawAtriStationEvent.h"
#include "UsefulAraStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraQualCuts.h"

RawAtriStationEvent *rawAtriEvPtr;
UsefulAtriStationEvent *realAtriEvPtr;

//ROOT Includes
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "FFTtools.h"

#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"

#include "tools_inputParameters.h"
#include "tools_outputObjects.h"
#include "tools_runSummaryObjects.h"
#include "tools_WaveformFns.h"
#include "tools_PlottingFns.h"
#include "tools_Cuts.h"

using namespace std;
TGraph *customInterpolation(TGraph *grIn, Double_t deltaT);

int main(int argc, char **argv)
{

	if(argc<3) {
		std::cout << "Usage\n" << argv[0] << " <simulation_flag> <station> <output_location> <input_file> <pedestal_file> \n";
		return -1;
	}

	int isSimulation = atoi(argv[1]);
	int station_num = atoi(argv[2]);
	string output_location = argv[3];

	int numAnts=16;
	double lowFreqLimit=120;
	double highFreqLimit=900.;
	int WaveformLength=2048;
	int newLength=(WaveformLength/2)+1;

	TFile *fp = TFile::Open(argv[4]);
	if(!fp) {
		std::cout << "Can't open file\n";
		return -1;
	}
	TTree *eventTree= (TTree*) fp->Get("eventTree");
	if(!eventTree) {
		std::cout << "Can't find eventTree\n";
		return -1;
	}

	int runNum = getrunNum(argv[4]);
	printf("Run Number %d \n", runNum);
	runNumber = runNum;

	AraQualCuts *qualCut = AraQualCuts::Instance(); //we also need a qual cuts tool


	if(argc==6){
		//only if they gave us a pedestal should we fire up the calibrator
		AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
		calibrator->setAtriPedFile(argv[5],station_num);
	}

	if(isSimulation){
		eventTree->SetBranchAddress("UsefulAtriStationEvent", &realAtriEvPtr);
		printf("Simulation; load useful event tree straight away \n");
	}
	else{
		eventTree->SetBranchAddress("event",&rawAtriEvPtr);
		printf("Data; load raw event tree \n");
	}

	double frequencyArray[newLength];
	double FFTre[16][newLength];
	double FFTim[16][newLength];
	double FFT[16][newLength];
	for(int j=0; j<newLength; j++){
		for(int i=0; i<16; i++){
			FFTre[i][j]=0;
			FFTim[i][j]=0;
			FFT[i][j]=0;
		}
		frequencyArray[j]=0;
	}
	int eventsIncludedCtr=0;

	Long64_t numEntries=eventTree->GetEntries();

	for(int event=0; event<eventTree->GetEntries();event++){
		eventTree->GetEvent(event);

		if (isSimulation == false){
			realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
		}
		bool isCalpulser;
		bool isSoftTrigger;

		if (isSimulation){
			isCalpulser = false;
			isSoftTrigger = false;
		} else{
			isCalpulser = rawAtriEvPtr->isCalpulserEvent();
			isSoftTrigger = rawAtriEvPtr->isSoftwareTrigger();
		}

		if(!isCalpulser && !isSoftTrigger && qualCut->isGoodEvent(realAtriEvPtr)){

			for(int chan=0; chan<numAnts; chan++){
				TGraph *gr = realAtriEvPtr->getGraphFromRFChan(chan); //get waveform				
				TGraph *grInt = customInterpolation(gr,0.6);
				TGraph *grPad = FFTtools::padWaveToLength(grInt,WaveformLength); //pad
				double *getX = grPad->GetX();
				double deltaT = getX[1]-getX[0];
				while(grPad->GetN()<WaveformLength){
					double lastX = 0.;
					double lastY = 0.;
					grPad->GetPoint(grPad->GetN()-1,lastX,lastY);
					grPad->SetPoint(grPad->GetN(),lastX+deltaT,0);
				}
				double *getY = grPad->GetY();
				int length = grPad->GetN();

				FFTWComplex *theFFT = FFTtools::doFFT(length,getY);
				double deltaF = 1/(deltaT * length);
				deltaF*=1e3;

				for(int samp=0; samp<newLength; samp++){
					FFTre[chan][samp]+=pow(theFFT[samp].re,2.); //the division by 1E5 is here to control the size of this number in memory
					FFTim[chan][samp]+=pow(theFFT[samp].im,2.); //the division by 1E5 is here to control the size of this number in memory

					if(chan==0){
						if(samp==0) frequencyArray[samp]=0.;
						if(samp>0) frequencyArray[samp]=frequencyArray[samp-1]+deltaF;
					}
				}
				delete [] theFFT;
				delete grPad;
				delete grInt;
				delete gr;
			}

			eventsIncludedCtr++;
		}
		if(isSimulation==false){
			delete realAtriEvPtr;
		}
	}

	if(eventsIncludedCtr>0){

		for(int chan=0; chan<numAnts; chan++){
			for(int samp=0; samp<newLength;samp++){
				FFTre[chan][samp]=FFTre[chan][samp]/double(eventsIncludedCtr);
				FFTim[chan][samp]=FFTre[chan][samp]/double(eventsIncludedCtr);
				FFT[chan][samp]=sqrt(FFTre[chan][samp] + FFTim[chan][samp]);
				if(FFT[chan][samp]>0.) FFT[chan][samp] = 10*log10(FFT[chan][samp]);
				if(frequencyArray[samp]<lowFreqLimit || frequencyArray[samp]>highFreqLimit){
					FFT[chan][samp]=0;
				}
			}
		}

		char run_file_name[400];
		sprintf(run_file_name,"%s/baseline_station_%d_run_%d.root",output_location.c_str(),station_num, runNum);
		TFile *outFile = TFile::Open(run_file_name,"RECREATE");
		TTree *BaselineTree = new TTree("BaselineTree","BaselineTree");

		vector <TGraph*> average;
		average.resize(16);

		stringstream ss;
		stringstream ss1;

		for(int i=0; i<16; i++){
			ss.str(""); ss<<"baselines_RF_chan_"<<i;
			BaselineTree->Branch(ss.str().c_str(),&average[i]);
		}

		for(int chan=0; chan<16; chan++){
			vector <double> freq;
			vector <double> amps;
			for(int samp=0; samp<newLength; samp++){
				freq.push_back(frequencyArray[samp]);
				amps.push_back(FFT[chan][samp]);
			}
			average[chan]=new TGraph(newLength,&freq[0],&amps[0]);
		}

		BaselineTree->Fill();
		outFile->Write();
		outFile->Close();
	}

	fp->Close();
	delete fp;
	printf("Done! Run Number %d \n", runNum);
}

TGraph *customInterpolation(TGraph *grIn, Double_t deltaT)
{
	std::vector<double> tVec;
	std::vector<double> vVec;

	Int_t numIn=grIn->GetN();
	Double_t tIn,vIn;

	Double_t startTime=0;
	Double_t lastTime=0;
	for (int samp=0;samp<numIn;samp++) {
		grIn->GetPoint(samp,tIn,vIn);
		tVec.push_back(tIn);
		vVec.push_back(vIn);
		if(samp==0) startTime=tIn;
		lastTime=tIn;
	}
	if(tVec.size()<1) {
		std::cout << "Insufficent points for interpolation\n";
		return NULL;
	}

	ROOT::Math::Interpolator chanInterp(tVec,vVec,ROOT::Math::Interpolation::kAKIMA);
	Int_t roughPoints=Int_t((lastTime-startTime)/deltaT);

	Double_t *newTimes = new Double_t[roughPoints+1]; //Will change this at some point, but for now
	Double_t *newVolts = new Double_t[roughPoints+1]; //Will change this at some point, but for now
	Int_t numPoints=0;
	for(Double_t time=startTime;time<=lastTime;time+=deltaT) {
		newTimes[numPoints]=time;
		newVolts[numPoints]=chanInterp.Eval(time);
		numPoints++;
	}

	TGraph *grInt = new TGraph(numPoints,newTimes,newVolts);
	delete [] newTimes;
	delete [] newVolts;
	return grInt;
}