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

#include "diode_include.h"

using namespace std;

int main(int argc, char **argv)
{

	if(argc<2) {  // Check to make sure there are enough arguments to do something meaningful
		std::cout << "Run like " << basename(argv[0]) << " <input data file>" << std::endl;
		return -1;
	}
	
	TFile *fpIn = new TFile(argv[1], "OLD"); //we're going to open the data file
	if(!fpIn){
		std::cerr<< "Can not open file: " <<argv[1]<<endl;
		return -1;
	} //throw a warning if you can't open it
	fpIn->cd(); //go into that file
	TTree *eventTree = (TTree*) fpIn->Get("eventTree"); //load in the event free for this file
	if(!eventTree){
		std::cerr<<"Can't find eventTree in file" <<argv[1]<<endl;
		return -1;
	} //throw a warning if you can't open it
	RawAtriStationEvent *rawAtriEvPtr=0;
	eventTree->SetBranchAddress("event",&rawAtriEvPtr);
	int runNum;
	eventTree->SetBranchAddress("run",&runNum);
	
	int numEntries = eventTree -> GetEntries(); //get the number of entries in this file

	AraQualCuts *qual = AraQualCuts::Instance();
	bool found=false;
	for(int event=0; event<numEntries; event++){ //loop over those entries
		if(found) break;
		eventTree->GetEntry(event); //get the event
		bool isCalpulser = rawAtriEvPtr->isCalpulserEvent();
		bool isSoftwareTrigger = rawAtriEvPtr->isSoftwareTrigger();
		if(!isSoftwareTrigger) continue;
		UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
		bool this_qual = qual->isGoodEvent(realAtriEvPtr); //get the quality
		if(!this_qual){
			delete realAtriEvPtr;
			continue;
		}
		found=true;
	
		//now, we'll get the waveform from channel 2
		TGraph *waveform = realAtriEvPtr->getGraphFromRFChan(0);
		TGraph *interpolated_waveform = FFTtools::getInterpolatedGraph(waveform, 0.5); //get an interpolated waveform with 0.5 ns interpolation
		TGraph *padded_waveform = FFTtools::padWaveToLength(interpolated_waveform,1024);

		TGraph *grOut = doConvolve(padded_waveform);
		TCanvas *c = new TCanvas("","",1100,850);
		c->Divide(1,2);
		c->cd(1);
			padded_waveform->Draw("ALP");
			padded_waveform->SetTitle("padded_waveform");
		c->cd(2);
			grOut->Draw("ALP");
			grOut->SetTitle("diode output");
		c->SaveAs("waveform_and_diode.png");
		delete c;

		//now do some cleanup
		delete interpolated_waveform;
		delete padded_waveform;
		delete waveform;
		delete realAtriEvPtr;
	}
}