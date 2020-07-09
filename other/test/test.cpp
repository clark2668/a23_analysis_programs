////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  find_hit_time.cxx
////
////  Nov 2018
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <iomanip>
#include <sstream>

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
#include "FFTtools.h"
#include "TLine.h"
#include "TCanvas.h"

RawAtriStationEvent *rawAtriEvPtr;

#include "AraGeomTool.h"
#include "AraQualCuts.h"
#include "AraAntennaInfo.h"

using namespace std;

int main(int argc, char **argv)
{

	if(argc<3) {
		std::cout << "Usage\n" << argv[0] << " <station> <data_file> \n";
		return -1;
	}

	AraQualCuts *qual = AraQualCuts::Instance();
		
	TFile *fp = TFile::Open(argv[2]);
	if(!fp) {
		std::cout << "Can't open file\n";
		return -1;
	}
	TTree *eventTree; 
	eventTree= (TTree*) fp->Get("eventTree");
	if(!eventTree) {
		std::cout << "Can't find eventTree\n";
		return -1;
	}
	int run;
	eventTree->SetBranchAddress("event",&rawAtriEvPtr);
	eventTree->SetBranchAddress("run",&run);

	int numEntries = eventTree->GetEntries();
	cout<<"numEntries is "<<numEntries<<endl;
	// numEntries=10000;
	int num_bad=0;
	for(int event=0; event<numEntries; event++){
		eventTree->GetEntry(event);
		cout<<"Event number is "<<rawAtriEvPtr->eventNumber<<endl;
		//cout<<"UnixTime is "<<rawAtriEvPtr->unixTime<<endl;
		  // // if(rawAtriEvPtr->eventNumber!=43) continue;
		  // // printf("Event number %d \n", rawAtriEvPtr->eventNumber);
		  // if(1!=1)
		  // 	int a=2;
		  // else{
		  //   // printf("Station ID is %d \n", rawAtriEvPtr->stationId);
		  //   UsefulAtriStationEvent *real = new UsefulAtriStationEvent(rawAtriEvPtr,AraCalType::kLatestCalib);
		  //   // printf("Event %d: hasOffsetBlocksFlag flag is %d \n", rawAtriEvPtr->eventNumber,  qual->hasOffsetBlocks(real));
		  //   bool quality=qual->isGoodEvent(real);
		  //   printf("Event %5d: is quality event %d \n", rawAtriEvPtr->eventNumber, quality);
		  //   if(!quality) num_bad++;
		  //   // break;
		  // }
		// UsefulAtriStationEvent *real = new UsefulAtriStationEvent(rawAtriEvPtr,AraCalType::kLatestCalib);
		// vector<TGraph*> graphs;
		// for(int i=0; i<16; i++) graphs.push_back(real->getGraphFromRFChan(i));
		// TCanvas *c = new TCanvas("","",4*1100,4*850);
		// c->Divide(4,4);
		// for(int i=0; i<16; i++){
		// 	c->cd(i+1);
		// 	graphs[i]->Draw("AL");
		// }
		// char save_title[150];
		// sprintf(save_title,"run%d_event%d.png",run,rawAtriEvPtr->eventNumber);
		// c->SaveAs(save_title);
		// delete c;
		// for(int i=0; i<16; i++) delete graphs[i];
		// delete real;

	}

	// printf("Num bad and num total are %d and %d \n", num_bad, numEntries);
}
