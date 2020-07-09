////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	print_waveform_out.cpp
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
	RawAtriStationEvent *rawAtriEvPtr=0;
	eventTree->SetBranchAddress("event",&rawAtriEvPtr);
	int runNum;
	eventTree->SetBranchAddress("run",&runNum);
	
	double numEntries = eventTree -> GetEntries(); //get the number of entries in this file

	for(int event=0; event<numEntries; event++){ //loop over those entries
		
		eventTree->GetEntry(event); //get the event

		if(rawAtriEvPtr->eventNumber!=39072) continue;

		UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
		vector<TGraph*> waveforms;
		vector<TGraph*> interpolated;
		for(int i=0; i<16; i++){
			waveforms.push_back(realAtriEvPtr->getGraphFromRFChan(i));
			interpolated.push_back(FFTtools::getInterpolatedGraph(waveforms[i],0.5));
		}
		for(int i=0; i<16; i++){
				char title_txt[200];
				sprintf(title_txt,"run%d_event%d_waveform_ch%d.txt",runNum,rawAtriEvPtr->eventNumber,i);
				FILE *fout = fopen(title_txt, "a");
				for(int samp=0; samp<interpolated[i]->GetN(); samp++){
					fprintf(fout,"%.3f, %.3f  \n", interpolated[i]->GetX()[samp], interpolated[i]->GetY()[samp]);
				}
				fclose(fout);//close sigmavsfreq.txt file
		}
		break;
	}
}//close the main program
