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
#include "araSoft.h"

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

	bool found=false;
	for(int event=2; event<3; event++){ //loop over those entries
		// if(found) break;
		eventTree->GetEntry(event); //get the event
		cout<<"Event is "<<event<<endl;

		//we can see if it has cal pulser or software trigger timing
		bool isCalpulser = rawAtriEvPtr->isCalpulserEvent();
		bool isSoftTrigger = rawAtriEvPtr->isSoftwareTrigger();
		cout<<"is software? " <<isSoftTrigger<<endl;

		int length = rawAtriEvPtr->blockVec.size();
		vector <int> single_chan_blocks;

		for(int i=0; i<length; i++){
			// cout<<rawAtriEvPtr->blockVec[i].getBlock()<<endl;
			if(i%4==0) single_chan_blocks.push_back((int) rawAtriEvPtr->blockVec[i].getBlock());
		}
		for(int i=0; i<single_chan_blocks.size(); i++){
			printf("Block %d is %d \n", i, single_chan_blocks[i]);
		}

		for(int i=0; i<MAX_TRIG_BLOCKS; i++){
			cout<<unsigned(rawAtriEvPtr->triggerInfo[i])<<endl;
		}

		if(isCalpulser) found=true;

		//make a *useful* event out of the *raw* event, which functionally just calibrates it
		UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
		vector<TGraph*> graphs;
		for(int i=0; i<16; i++){
			graphs.push_back(realAtriEvPtr->getGraphFromRFChan(i));
		}
		TCanvas *c = new TCanvas("","",1100,850);
		c->Divide(4,4);
		for(int i=0; i<16; i++){
			c->cd(i+1);
			graphs[i]->Draw("alp");
		}
		c->SaveAs("quick.png");
		delete realAtriEvPtr;
	}
	
}//close the main program
