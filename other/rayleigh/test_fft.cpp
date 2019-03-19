////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	save_fft.cxx
////
////	March 2019,  clark.2668@osu.edu
////	Make an FFT of every soft trigger and save the spectral coefficient
////////////////////////////////////////////////////////////////////////////////

// C/C++ Includes
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <algorithm>

//AraRoot Includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "FFTtools.h"
#include "AraGeomTool.h"
#include "AraQualCuts.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"

using namespace std;

TGraph *makeFFTPlot(TGraph *grWave);
TGraph *makeFreqV_MilliVoltsNanoSeconds ( TGraph *grWave );

int main(int argc, char **argv)
{
	if(argc<3) {  // Check to make sure there are enough arguments to do something meaningful
		std::cout << "Usage requires you to provide input parameter of the form " << basename(argv[0]) << " <output_dir> <data_file>" << std::endl;
		return -1;
	}
	
	TFile *fpIn = new TFile(argv[2], "OLD"); //we're going to open the data file
	if(!fpIn){
		std::cerr<< "Can not open the old file: " <<argv[2]<<endl;
		return -1;
	} //throw a warning if you can't open it
	
	fpIn->cd(); //go into that file
	TTree *eventTree = (TTree*) fpIn->Get("eventTree"); //load in the event free for this file
	if(!eventTree){
		std::cerr<<"Can't find eventTree in file" <<argv[2]<<endl;
		return -1;
	} //throw a warning if you can't open it
	RawAtriStationEvent *rawAtriEvPtr=0;
	eventTree->SetBranchAddress("event",&rawAtriEvPtr);
	int runNum;
	eventTree->SetBranchAddress("run",&runNum);
	eventTree->GetEntry(0);
	int station=rawAtriEvPtr->stationId;
	int unixtime = rawAtriEvPtr->unixTime;
	time_t test_time = unixtime;
	tm *time = gmtime( &test_time );
	int year = time->tm_year+1900;

	printf("Filter Run Number %d \n", runNum);

	char ped_file_name[400];
	if(year==2013){
		sprintf(ped_file_name,"/fs/scratch/PAS0654/ara/peds/run_specific_peds/A%d/%d/event%d_specificPeds.dat",station,year,runNum);
	}
	else if(year==2014 || year==2015 || year==2016){
		sprintf(ped_file_name,"/fs/scratch/PAS0654/ara/peds/run_specific_peds/A%d/%d/event00%d_specificPeds.dat",station,year,runNum);
	}
	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	calibrator->setAtriPedFile(ped_file_name,station); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist

	double numEntries = eventTree -> GetEntries(); //get the number of entries in this file

	AraQualCuts *qual = AraQualCuts::Instance();
	int num_found=0;
	for(int event=0; event<numEntries; event++){ //loop over those entries
		if(num_found>2) break;
		eventTree->GetEntry(event); //get the event
		bool isSoftTrigger = rawAtriEvPtr->isSoftwareTrigger();
		if(isSoftTrigger){
			UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
			if(qual->isGoodEvent(realAtriEvPtr)){
				for(int chan=0; chan<16; chan++){
					TGraph *grRaw = realAtriEvPtr->getGraphFromRFChan(chan);
					TGraph *grInt = FFTtools::getInterpolatedGraph(grRaw, 0.5);
					TGraph *grPad = FFTtools::padWaveToLength(grInt,1024);
					// TGraph *spec = makeFFTPlot(grPad);
					TGraph *spec = makeFreqV_MilliVoltsNanoSeconds(grPad);
					delete spec;
					delete grPad;
					delete grInt;
					delete grRaw;
				}
				num_found++;
			}
			delete realAtriEvPtr;
		}
	}

	fpIn->Close();
	delete fpIn;
}//close the main program

TGraph *makeFFTPlot(TGraph *grWave)
{
	double *oldY = grWave->GetY();
	double *oldX = grWave->GetX();
	double deltaT=oldX[1]-oldX[0];
	int length=grWave->GetN();
	FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);
     
	int newLength=(length/2)+1;
     
	double *newY = new double [newLength];
	double *newX = new double [newLength];
     
	double deltaF=1/(deltaT*length); //Hz
	deltaF*=1e3; //MHz
       
	double tempF=0;
	for(int i=0;i<newLength;i++) {
		float power=FFTtools::getAbs(theFFT[i]);//converting from mV to V;
		if(i>0 && i<newLength-1) power*=2; //account for symmetry
		power*=sqrt(deltaT)/(length); //For time-integral squared amplitude
		power/=sqrt(deltaF);//Just to normalise bin-widths
		newX[i]=tempF;
		newY[i]=power; //Units should be mV/MhZ for the y-axis.
		tempF+=deltaF;
	}
   
	TGraph *grPower = new TGraph(newLength,newX,newY);
	delete [] theFFT;
	delete [] newY;
	delete [] newX;
	return grPower;
}

TGraph *makeFreqV_MilliVoltsNanoSeconds ( TGraph *grWave ) {
	double *oldY = grWave->GetY(); // mV
	double *oldX = grWave->GetX(); // ns
	int length=grWave->GetN();
	double deltaT = (oldX[1]-oldX[0]) * 1.e-9; // deltaT in s
	FFTWComplex *theFFT=FFTtools::doFFT(length,oldY); // FFT with mV unit
	int newLength=(length/2)+1;
	double *newY = new double [newLength];
	double *newX = new double [newLength];
	double deltaF=1./(deltaT*(double)length); //Hz
	deltaF*=1e-6; //from Hz to MHz
	double tempF=0;
	for(int i=0;i<newLength;i++) {
		newY[i] = sqrt(2.)*FFTtools::getAbs(theFFT[i]); // from mV to V		
		// newY[i] = sqrt(pow(FFTtools::getAbs(theFFT[i]),2.)*2.); // from mV to V
		newX[i]=tempF;
		tempF+=deltaF;
	}
	TGraph *grPower = new TGraph(newLength,newX,newY);
	delete [] theFFT;
	delete [] newY;
	delete [] newX;

	double time_power=0.;
	for(int i=0; i<grWave->GetN(); i++){
		time_power+=pow(grWave->GetY()[i],2.);
	}
	double freq_power=0.;
	for(int i=0; i<grPower->GetN(); i++){
		freq_power+=pow(grPower->GetY()[i],2.);
	}
	freq_power/=double(grWave->GetN());
	printf("Time power is %.2f \n", time_power);
	printf("Freq power is %.2f \n", freq_power);
	printf("Ratio is %.4f \n", time_power/freq_power);

	return grPower;
}

