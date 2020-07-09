////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	save_sim_RPR.cpp
////	save the RPR for a bunch of simulation
////	locatoin: other/sim_number_of_hits/save_sim_RPR.cxx
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

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"

//AraRoot Includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "FFTtools.h"
#include "AraQualCuts.h"
#include "AraRecoHandler.h"

//AraSim Includes
#include "Position.h"
#include "Report.h"
#include "Event.h"
#include "Detector.h"

//custom tools includes
#include "tools_PlottingFns.h"


using namespace std;
UsefulAtriStationEvent *realAtriEvPtr;

int main(int argc, char **argv)
{

	if(argc<4) {  // Check to make sure there are enough arguments to do something meaningful
		std::cout << "Usage requires you to provide input parameter of the form " << basename(argv[0]) << " <station> <output_directory> <input_sim_file>" << std::endl;
		return -1;
	}
	int station = atoi(argv[1]);

	double interpV = 0.4;
	double interpH = 0.625;
	int nIntSamp_V = int(25./interpV);
	int nIntSamp_H = int(25./interpH);

	AraRecoHandler *RecoHandler = new AraRecoHandler();

	// now start "main"
	////////////////////
	////////////////////
	
	for(int file=3; file<argc; file++){

		TFile *fpIn = new TFile(argv[file], "OLD"); //we're going to open the data file
			if(!fpIn){ std::cerr<< "Can not open the old file: " <<argv[file]<<endl; return -1; } //throw a warning if you can't open it
		TTree *eventTree = (TTree*) fpIn->Get("eventTree"); //load in the event free for this file
			if(!eventTree){ std::cerr<<"Can't find eventTree in file" <<argv[file]<<endl; return -1; } //throw a warning if you can't open it
		TTree *AraTree2 = (TTree*) fpIn->Get("AraTree2");
			if(!AraTree2){  std::cerr<<"Can't find AraTree2  in file" <<argv[file]<<endl; return -1; }

		int runNum = getrunNum(argv[file]);
		char outfile_name[400];
		sprintf(outfile_name,"%s/rpr_run%d.root",argv[2],runNum);
		TFile *fpOut = TFile::Open(outfile_name,"RECREATE");
		TTree *outTree = new TTree("outTree", "outTree");
		double weight_out;
		double RPR[16];
		double pnu;
		double times[16][3]; // sixteen waveforms; 3 entries per waveform (start time, end time, and hit time)
		outTree->Branch("weight",&weight_out,"weight/D");
		outTree->Branch("RPR",&RPR,"RPR[16]/D");
		outTree->Branch("pnu",&pnu,"pnu/D");
		outTree->Branch("times",&times,"times[16][3]/D");

		double weight;
		eventTree->SetBranchAddress("UsefulAtriStationEvent",&realAtriEvPtr);
		eventTree->SetBranchAddress("weight",&weight);

		Report *reportPtr;
		Event *eventPtr;
		AraTree2->SetBranchAddress("report",&reportPtr);
		AraTree2->SetBranchAddress("event",&eventPtr);

		cout<<"Num AraTree2 entries is "<<AraTree2->GetEntries()<<endl;
		cout<<"Num eventTree entries is "<<eventTree->GetEntries()<<endl;

		double numEntries = AraTree2 -> GetEntries(); //get the number of entries in this file
		// double numEntries = eventTree -> GetEntries(); //get the number of entries in this file

		int start=0;
		int eventTrig=0;
		for(int event=start; event<numEntries; event++){ //loop over those entries		
			AraTree2->GetEvent(event);
			// only try to access eventTree information if we know the event globally triggered
			int globalpass = reportPtr->stations[0].Global_Pass;
			if(globalpass){
				eventTree->GetEntry(eventTrig); //get the event
				pnu=eventPtr->pnu;
				vector<TGraph*> waveforms_int;
				for(int i=0; i<16; i++){
					TGraph *gr = realAtriEvPtr->getGraphFromRFChan(i);
					times[i][0] = gr->GetX()[0]; // time of first sample
					times[i][1] = gr->GetX()[gr->GetN()-1]; // time of last sample
					TGraph *grInt = FFTtools::getInterpolatedGraph(gr,i<8?interpV:interpH);
					waveforms_int.push_back(grInt);
					delete gr;
				}

				float RPRs[16];
				float hitTimes[16];
				RecoHandler->getChannelSlidingV2SNR_UW(waveforms_int, nIntSamp_V, nIntSamp_H, RPRs, hitTimes);
				
				for(int i=0; i<16; i++){
					RPR[i]=RPRs[i];
					times[i][2]=hitTimes[i];
				}
				weight_out=weight;
				outTree->Fill();

				// int selecChan=1;
				// printf("Event %d, RPR %.2f, start time %.2f, end time %.2f, hit time %.2f \n", event, RPR[selecChan], times[selecChan][0], times[selecChan][1], times[selecChan][2]);

				// clean up and increment triggered event counter
				for(int i=0; i<16; i++) delete waveforms_int[i];
				eventTrig++; //advance the "triggered" counter
			}
		}

		fpOut->Write();
		fpOut->Close();
		delete fpOut;

		fpIn->Close();
		delete fpIn;
	}
}//close the main program














