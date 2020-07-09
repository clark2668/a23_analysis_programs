////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	num trig chans check
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
	numEntries=101;
	for(int event=0; event<numEntries; event++){ //loop over those entries
		
		eventTree->GetEntry(event); //get the event
		if(rawAtriEvPtr->isSoftwareTrigger() || rawAtriEvPtr->isCalpulserEvent())
		// if(rawAtriEvPtr->isSoftwareTrigger())	
			continue;
		printf("Event %d has %d chans high and has cal status %d \n", event, rawAtriEvPtr->numTriggerChansHigh(), rawAtriEvPtr->isCalpulserEvent());
		for(int i=0; i<16 ;i++){
			if(i==0 || i==1 || i==4 || i==5 || i==8 || i==9 || i==12 || i==13)
				printf("		VPol Chan %d has high status %d \n", i, rawAtriEvPtr->isTriggerChanHigh(i));
		}
		for(int i=0; i<16 ;i++){
			if(i!=0 && i!=1 && i!=4 && i!=5 &&i!=8 && i!=9 && i!=12 && i!=13)
				printf("		Hpol Chan %d has high status %d \n", i, rawAtriEvPtr->isTriggerChanHigh(i));
		}
	}	
}//close the main program
