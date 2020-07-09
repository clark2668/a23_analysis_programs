////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  sync_sim.cxx 
////  code to synchronize eventTree and AraTree2
////
////  Nov 2019
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

RawAtriStationEvent *rawAtriEvPtr;
UsefulAtriStationEvent *realAtriEvPtr;

// #include "/home/brianclark/A23_analysis_new2/AraRoot/include/Event.h"
// #include "/home/brianclark/A23_analysis_new2/AraRoot/include/Detector.h"
// #include "/home/brianclark/A23_analysis_new2/AraRoot/include/Report.h"

#include "Event.h"
// #include "Detector.h"
#include "Report.h"

#include "AraGeomTool.h"
#include "AraQualCuts.h"
#include "AraAntennaInfo.h"

#include "tools_inputParameters.h"
#include "tools_outputObjects.h"
#include "tools_runSummaryObjects.h"
#include "tools_WaveformFns.h"
#include "tools_PlottingFns.h"
#include "tools_Constants.h"
#include "tools_RecoFns.h"
#include "tools_filterEvent.h"
#include "tools_Cuts.h"
#include "tools_CommandLine.h"

int main(int argc, char **argv)
{

	if(argc<6) {
		std::cout << "Usage\n" << argv[0] << " <1-station> <2-config> <3-energy> <4-output_directory> <5-input_file> \n";
		return -1;
	}

	/*
	arguments
	0: exec
	1: station
	2: config
	3: energy
	4: output directory
	5: input file
	*/

	int station = atoi(argv[1]);
	int config = atoi(argv[2]);
	int energy = atoi(argv[3]);
	string output_directory = argv[4];

	for(int file=5; file<argc; file++){
			
		TFile *fp = TFile::Open(argv[file],"read");
		if(!fp) {
			std::cout << "Can't open file\n";
			return -1;
		}
		int runNum = getrunNum(argv[file]);
		
		// first, the eventTree

			TTree *org_eventTree; 
			org_eventTree = (TTree*) fp->Get("eventTree");

			// if nothing passes the event tree (no triggered events) then we have nothing to do
			if(!org_eventTree) {
				std::cout << "Can't find eventTree\n";
				return -1;
			}

			int N_eventTree = org_eventTree->GetEntries();

			if(N_eventTree<1){
				cout<<"There are no events in the eventTree (no events triggered) so do nothing!"<<endl;
				return -1;
			}
			double weight_eventTree;
			org_eventTree->SetBranchAddress("weight",&weight_eventTree);

		// second, the AraTree2

			TTree *org_AraTree2;
			org_AraTree2 = (TTree*) fp->Get("AraTree2");
			if(!org_AraTree2) {
				std::cout << "Can't find AraTree2\n";
				return -1;
			}
			int N_AraTree2 = org_AraTree2->GetEntries();

			Event *eventPtr = 0;
			Report *reportPtr = 0;
			org_AraTree2->SetBranchAddress("report",&reportPtr);
			org_AraTree2->SetBranchAddress("event",&eventPtr);

		// third, the AraTree

			// TTree *org_AraTree;
			// org_AraTree = (TTree*) fp->Get("AraTree");
			// if(!org_AraTree) {
			// 	std::cout << "Can't find AraTree\n";
			// 	return -1;
			// }

		char outputFileName[500];
		sprintf(outputFileName,"%s/A%d_c%d_E%d_run%d.root",output_directory.c_str(),station,config,energy,runNum);
		TFile *outFile = TFile::Open(outputFileName,"RECREATE");
		if(!outFile){
			cout<<"I can't open the output file "<<outputFileName<<" so I'm quitting!"<<endl;
			return -1;
		}
		TTree *clone_eventTree = org_eventTree->CloneTree(0);
		// printf(RED"GOT HERE L%d !\n"RESET,__LINE__);
		TTree *clone_AraTree2 = org_AraTree2->CloneTree(0);
		// printf(RED"GOT HERE L%d !\n"RESET,__LINE__);
		// TTree *clone_AraTree = org_AraTree->CloneTree(0);
		// printf(RED"GOT HERE L%d !\n"RESET,__LINE__);

		// first, just copy the AraTree
		// org_AraTree->GetEvent(0);
		// clone_AraTree->Fill();

		// printf("Num events eventTree %d, Num Events AraTree2 %d \n", N_eventTree, N_AraTree2);

		int doVersion2=true;
		if(doVersion2){
			// printf(BLUE"I'm doing version 2 \n"RESET);
			vector<int> eventList_AraTree2;
			for(int j_AraTree2=0; j_AraTree2<N_AraTree2; j_AraTree2++){
				org_AraTree2->GetEvent(j_AraTree2);
				if(reportPtr->stations[0].Global_Pass>0){
					eventList_AraTree2.push_back(j_AraTree2);
				}
			}

			for(int i_eventTree=0; i_eventTree<N_eventTree; i_eventTree++){
				org_eventTree->GetEvent(i_eventTree);
				int j_AraTree2 = eventList_AraTree2[i_eventTree];
				org_AraTree2->GetEvent(j_AraTree2);

				double weight_AraTree2 = eventPtr->Nu_Interaction[0].weight;
				// printf("    eventTree event i_eventTree %d with weight %.6f matches AraTree2 j_AraTree2 %d with weight %.6f \n", i_eventTree, weight_eventTree, j_AraTree2, weight_AraTree2);

				clone_eventTree->Fill();
				clone_AraTree2->Fill();
			}
		}

		
		// int doVersion1=false;
		// if(doVersion1){
		// 	int numNeedToFind = N_eventTree;

		// 	for(int i_eventTree=0; i_eventTree<N_eventTree; i_eventTree++){
		// 		org_eventTree->GetEvent(i_eventTree);
		// 		// cout<<"eventTree entry i_eventTree "<<i_eventTree<<" has weight "<<weight_eventTree<<endl;

		// 		int numFound=0;
		// 		int j_AraTree2=0;

		// 		for(int j_AraTree2=0; j_AraTree2<N_AraTree2; j_AraTree2++){
		// 			org_AraTree2->GetEvent(j_AraTree2);
		// 			if(reportPtr->stations[0].Global_Pass>0){
		// 				double weight_AraTree2 = eventPtr->Nu_Interaction[0].weight;
		// 				numFound++;
		// 				if(numFound==i_eventTree+1){
		// 					// cout<<"    AraTree2 entri j_AraTree2 "<<j_AraTree2<<" has weight "<<weight_AraTree2<<endl;
		// 					printf("    eventTree event i_eventTree %d with weight %.6f matches AraTree2 j_AraTree2 %d with weight %.6f \n", i_eventTree, weight_eventTree, j_AraTree2, weight_AraTree2);

		// 					// once we synchronize them, write them down!
		// 					clone_eventTree->Fill();
		// 					clone_AraTree2->Fill();

		// 				}
		// 			}
		// 		}
		// 	}
		// }

		outFile->Write();
		outFile->Close();
		delete outFile;

		fp->Close();
		delete fp;

		printf("Done run %d \n", runNum);
	}
}
