////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	print first four events
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
#include "TF1.h"

#include "AraQualCuts.h"

using namespace std;

int main(int argc, char **argv)
{
	if(argc<2) {  // Check to make sure there are enough arguments to do something meaningful
		std::cout << "Usage requires you to provide input parameter of the form " << basename(argv[0]) << " <station> <input data file 1 > <input data file 2> ..." << std::endl;
		return -1;
	}
	AraQualCuts *qualCut = AraQualCuts::Instance();
	qualCut->_OffsetBlocksTimeWindowCut=100.;
	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();

	TH1D *MaximumRollingMean[16];
	for(int i=0; i<16; i++){
		MaximumRollingMean[i] = new TH1D("","",100,0,100);
	}

	int station = atoi(argv[1]);
	for(int file=2; file<argc; file++){
		TFile *fpIn = new TFile(argv[file], "OLD"); //we're going to open the data file
		if(!fpIn){
			std::cerr<< "Can not open the old file: " <<argv[file]<<endl;
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
		eventTree->GetEntry(0);
		double numEntries = eventTree -> GetEntries(); //get the number of entries in this file
		char *PedDirPath(getenv("PED_DIR"));
		if (PedDirPath == NULL) std::cout << "Warning! $PED_DIR is not set!" << endl;
		char ped_file_name[400];
		sprintf(ped_file_name,"%s/run_specific_peds/A%d/all_peds/event%d_specificPeds.dat",PedDirPath,station,runNum);
		calibrator->setAtriPedFile(ped_file_name,station); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist
		int start=0;
		// numEntries=13566;
		numEntries=2000;
		for(int event=start; event<numEntries; event++){ //loop over those entries
			eventTree->GetEntry(event); //get the event
			printf("On event %d \n", event);
			int eventNumber=(int)rawAtriEvPtr->eventNumber;
			UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
			vector<TGraph*> graphs;
			vector<TGraph*> rollingMeans;
			bool PrintThisEvent=false;
			for(int i=0; i<16; i++){
				graphs.push_back(realAtriEvPtr->getGraphFromRFChan(i));
				if(graphs[i]->GetN()<65) continue;
				TGraph *grInt = FFTtools::getInterpolatedGraph(graphs[i],0.5);
				rollingMeans.push_back(qualCut->getRollingMean(grInt,64));
				for(int samp=0; samp<rollingMeans[i]->GetN(); samp++){
					rollingMeans[i]->GetX()[samp]+=graphs[i]->GetX()[0];
				}
				delete grInt;
				double this_max = abs(TMath::MaxElement(rollingMeans[i]->GetN(), rollingMeans[i]->GetY()));
				double this_min = abs(TMath::MinElement(rollingMeans[i]->GetN(), rollingMeans[i]->GetY()));
				double abs_max;
				if(this_max>this_min) abs_max=this_max;
				else if (this_min>this_max) abs_max=this_min;
				if(abs_max>=98) abs_max=98.;

				// separations[i]->Fill(this_max-this_min);
				MaximumRollingMean[i]->Fill(abs_max);
				// if(this_max-this_min>800 && (i!=3 && i!=7 && i!=11 && i!=15)) PrintThisEvent=true;
				// printf("Max and min are %.2f and %.2f with separation %.2f \n", this_max, this_min, this_max-this_min);
			}
			// bool hasOffsetBlock = qualCut->hasOffsetBlocks(realAtriEvPtr);
			// if(hasOffsetBlock) PrintThisEvent=false;
			if(PrintThisEvent){
				TCanvas *c = new TCanvas("","",1000,1000);
				c->Divide(4,4);
				for(int i=0; i<16; i++){
					c->cd(i+1);
					graphs[i]->Draw("ALP");
					rollingMeans[i]->Draw("Lsame");
					rollingMeans[i]->SetLineColor(kRed);
				}
				char this_save_name[400];
				sprintf(this_save_name,"./outputs/run%d_event%d.png",runNum,event);
				c->SaveAs(this_save_name);
				delete c;
			}

			// bool isGoodEvent = qualCut->isGoodEvent(realAtriEvPtr);
			// printf("Event %2d is Good Events %d and hasOffsetBlocks is %d \n", event, isGoodEvent, hasOffsetBlock);
			printf("------------------------------------------------------------\n");

			for(int i=0; i<16; i++) delete graphs[i];
			delete realAtriEvPtr;
		}
		fpIn->Close();
	}

	TCanvas *c = new TCanvas("","",2*1000,2*1000);
	c->Divide(4,4);
	for(int i=0; i<16; i++){
		c->cd(i+1);
		MaximumRollingMean[i]->Draw("");
		MaximumRollingMean[i]->GetXaxis()->SetTitle("Absolute Largest Value in 64 Sample Rolling Mean");
		MaximumRollingMean[i]->GetYaxis()->SetTitle("Number of Events");
		MaximumRollingMean[i]->SetLineWidth(2);
		MaximumRollingMean[i]->GetXaxis()->SetRangeUser(20,50);
		gPad->SetLogy();
	}
	char fit_equation[150];
	sprintf(fit_equation,"gaus");
	TF1 *fit = new TF1("ExpoFit",fit_equation,20,30);
	MaximumRollingMean[0]->Fit("ExpoFit","R");
	printf("Chi-Square/NDF %.2f / %.2f \n",fit->GetChisquare(),double(fit->GetNDF()));
	char this_save_name[400];
	sprintf(this_save_name,"./outputs/distribution_of_maximum_rolling_means.png");
	c->SaveAs(this_save_name);
	delete c;
}//close the main program
