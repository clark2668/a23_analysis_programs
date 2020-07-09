////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	analysis.cxx
////
////	January 2019,  clark.2668@osu.edu
////	This is an example of how you analyze ARA data
////	We will learn how to get a waveform
////	How to make a spectrum
////	We will run this as *compiled* code, NOT run as a ROOT macro
////	This code executes over raw (L0) data
////
////	We will also learn how to use the ARA geom tool
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

int main(int argc, char **argv)
{
	if(argc<2) {  // Check to make sure there are enough arguments to do something meaningful
		std::cout << "Usage requires you to provide input parameter of the form " << basename(argv[0]) << " <input data file>" << std::endl;
		return -1;
	}
	
	TFile *fpIn = new TFile(argv[1], "OLD"); //we're going to open the data file
	if(!fpIn){
		std::cerr<< "Can not open the old file: " <<argv[1]<<endl;
		return -1;
	} //throw a warning if you can't open it
	
	fpIn->cd(); //go into that file
	TTree *eventTree = (TTree*) fpIn->Get("eventTree"); //load in the event free for this file
	if(!eventTree){
		std::cerr<<"Can't find eventTree in file" <<argv[1]<<" with filename " <<argv[1]<<endl;
		return -1;
	} //throw a warning if you can't open it
	 //set the tree address to access our raw data type
	RawAtriStationEvent *rawAtriEvPtr=0;
	eventTree->SetBranchAddress("event",&rawAtriEvPtr);
	int runNum;
	//we can also get the run number
	eventTree->SetBranchAddress("run",&runNum);
	
	double numEntries = eventTree -> GetEntries(); //get the number of entries in this file

	AraQualCuts *qual = AraQualCuts::Instance();
	// numEntries=3;
	// for(int event=0; event<numEntries; event++){ //loop over those entries
	for(int event=12610-10; event<12624+10; event++){ //loop over those entries
		
		eventTree->GetEntry(event); //get the event

		if(rawAtriEvPtr->eventNumber!=131864) continue;

		// if(rawAtriEvPtr->eventNumber==131864) cout<<"Event is "<<event<<endl;
		// cout<<"Event is "<<event<<endl;

		// //we can see if it has cal pulser or software trigger timing
		// bool isCalpulser = rawAtriEvPtr->isCalpulserEvent();
		// bool isSoftTrigger = rawAtriEvPtr->isSoftwareTrigger();

		// //and also at what unixtime and unixtime microsecond the event was read out by the SBC
		// int unixTime=(int)rawAtriEvPtr->unixTime;
		// int unixTimeUs=(int)rawAtriEvPtr->unixTimeUs;

		// //we can also get the (real) firmware event number the event was stamped with
		// int eventNumber=(int)rawAtriEvPtr->eventNumber;
		
		// //make a *useful* event out of the *raw* event, which functionally just calibrates it
		UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);

		// vector<string> conditioning = realAtriEvPtr->fConditioningList;
		// for(int i=0; i<conditioning.size(); i++){
		// 	cout<<"conditioning value "<<i<<" is "<<conditioning[i]<<endl;
		// }

		// bool this_qual = qual->isGoodEvent(realAtriEvPtr); //get the quality
		// printf("Iterator event %d, is real event number %d and has qual %d \n", event, eventNumber, this_qual);
		// if(!this_qual) continue;
		// cout<<"Qual is "<<this_qual<<endl;

		vector<TGraph*> waveforms;
		vector<TGraph*> interpolated;
		vector<TGraph*> padded;
		vector<TGraph*> spectra;
		for(int i=0; i<16; i++){
			waveforms.push_back(realAtriEvPtr->getGraphFromRFChan(i));
			interpolated.push_back(FFTtools::getInterpolatedGraph(waveforms[i],0.5));
			padded.push_back(FFTtools::padWaveToLength(interpolated[i],2048));
			spectra.push_back(FFTtools::makePowerSpectrumMilliVoltsNanoSeconds(padded[i]));

		}
		for(int i=0; i<16; i++){};


		string xLabel = "Time (ns)"; string yLabel = "Voltage (mV)";
		string xLabel_spec = "Frequency (MHz)"; string yLabel_spec = "Power";
		vector<string> titlesForGraphs;
		stringstream ss;
		for (int i = 0; i < 16; i++){
			ss.str("");
			ss << "Channel " << i;
			titlesForGraphs.push_back(ss.str());
		}
		vector<TGraph*> dummy;
		vector<TGraph*> dummySpectra;
		for(int i=0; i<16; i++){
			vector<double> Xs;
			vector<double> Ys;
			Xs.push_back(-500.);
			Xs.push_back(500.);
			Ys.push_back(-1000.);
			Ys.push_back(1000.);
			dummy.push_back(new TGraph(Xs.size(), &Xs[0], &Ys[0]));
			dummy[i]->GetXaxis()->SetTitle(xLabel.c_str());
			dummy[i]->GetYaxis()->SetTitle(yLabel.c_str());
			dummy[i]->SetTitle(titlesForGraphs[i].c_str());

			vector<double> Xs_spec;
			vector<double> Ys_spec;
			Xs_spec.push_back(0.);
			Xs_spec.push_back(1000.);
			Ys_spec.push_back(1e0);
			Ys_spec.push_back(2e5);
			dummySpectra.push_back(new TGraph(Xs_spec.size(), &Xs_spec[0], &Ys_spec[0]));
			dummySpectra[i]->GetXaxis()->SetTitle(xLabel_spec.c_str());
			dummySpectra[i]->GetYaxis()->SetTitle(yLabel_spec.c_str());
			dummySpectra[i]->SetTitle(titlesForGraphs[i].c_str());
		}

		TCanvas *c = new TCanvas("","",4*1100,4*850);
		c->Divide(4,4);
		for(int i=0; i<16; i++){
			c->cd(i+1);
			dummy[i]->Draw("AP");
			waveforms[i]->Draw("same");
			waveforms[i]->SetLineWidth(2);
			dummy[i]->GetXaxis()->SetRangeUser(-200.,500.);
		}
		char some_titles[400];
		sprintf(some_titles,"run%d_event%d_waveforms.png",runNum,rawAtriEvPtr->eventNumber);
		c->SaveAs(some_titles);

		TCanvas *c2 = new TCanvas("","",4*1100,4*850);
		c2->Divide(4,4);
		for(int i=0; i<16; i++){
			c2->cd(i+1);
			dummySpectra[i]->Draw("AP");
			gPad->SetLogy();
			spectra[i]->Draw("same");
			spectra[i]->SetLineWidth(2);
		}
		sprintf(some_titles,"run%d_event%d_spectra.png",runNum,rawAtriEvPtr->eventNumber);
		c2->SaveAs(some_titles);
	
		// //now, we'll get the waveform from channel 2
		// TGraph *waveform = realAtriEvPtr->getGraphFromRFChan(0);
		// TGraph *interpolated_waveform = FFTtools::getInterpolatedGraph(waveform, 0.5); //get an interpolated waveform with 0.5 ns interpolation
		// TGraph *padded_waveform = FFTtools::padWaveToLength(interpolated_waveform,2048);
		// TGraph *spectrum = FFTtools::makePowerSpectrumMilliVoltsNanoSeconds(padded_waveform); //now make a spectrum
	
		// //now do some cleanup
		// delete spectrum;
		// delete interpolated_waveform;
		// delete padded_waveform;
		// delete waveform;
		// delete realAtriEvPtr;
	}


	// /*
	// Now we can also see how to use the geom tool to get measurement antenna information
	// */

	// int station=2;
	// AraGeomTool *araGeom = AraGeomTool::Instance();

	// for(int i=0; i<16; i++){ //loop over antennas
	// 	int pol = (int) araGeom->getStationInfo(station)->getAntennaInfo(i)->polType; //polarization
	// 	double X = araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[0]; //antenna X location
	// 	double Y = araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[1]; //antenna Y location
	// 	double Z = araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[2]; //antenna Z location
	// 	double delay = araGeom->getStationInfo(station)->getCableDelay(i); //the associated cable delay
	// }

	// /*
	// Now we can also see how to use the geom tool to get cal-pulser antenna information
	// */

	// for(int i=0; i<araGeom->getStationInfo(station)->getNumCalAnts(); i++){ //loop over number of cal antennas
	// 	double X = araGeom->getStationInfo(station)->getCalAntennaInfo(i)->antLocation[0];
	// 	double Y = araGeom->getStationInfo(station)->getCalAntennaInfo(i)->antLocation[1];
	// 	double Z = araGeom->getStationInfo(station)->getCalAntennaInfo(i)->antLocation[2];
	// 	string locName(&araGeom->getStationInfo(station)->getCalAntennaInfo(i)->locationName[0]);
	// 	string antName(araGeom->getStationInfo(station)->getCalAntennaInfo(i)->getCalAntName());
	// }
	
}//close the main program
