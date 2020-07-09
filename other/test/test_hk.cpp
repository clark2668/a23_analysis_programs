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
#include <typeinfo>


//AraRoot Includes
#include "RawIcrrStationEvent.h"
#include "RawAtriStationEvent.h"
#include "UsefulAraStationEvent.h"
#include "UsefulIcrrStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AtriEventHkData.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "FFTtools.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TTimeStamp.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLine.h"

RawAtriStationEvent *rawAtriEvPtr;
AtriEventHkData *eventHkPtr;

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

	vector<double> thresh;
	vector<double> unixTime;

	TTimeStamp start;
	TTimeStamp stop;
	// start.Set(2016, 11, 05, 00, 00,0,0,true,0);
	// stop.Set(2016, 11, 10, 00, 00,0,0,true,0);

	// start.Set(2013, 11, 16, 00, 00,0,0,true,0);
	// stop.Set(2013, 11, 21, 00, 00,0,0,true,0);

	start.Set(2013, 10, 02, 21, 30,0,0,true,0);
	stop.Set(2013, 10, 02, 22, 00,0,0,true,0);

	int start_bin = start.GetSec();
	int stop_bin = stop.GetSec();
	int numBins=stop_bin - start_bin;
	numBins/=10;

	// int start_bin=7750;
	// int stop_bin=7780;

	// int start_bin=1790;
	// int stop_bin=1820;

	// int numBins=stop_bin - start_bin;

	int runEndunixTime;

	// TProfile *rate = new TProfile("","",numBins,start_bin,stop_bin);
	TH2D *rate = new TH2D("","",numBins,start_bin,stop_bin,100,0,50000);

	for(int file=2; file<argc; file++){
		
		TFile *fp = TFile::Open(argv[file]);
		if(!fp) {
			std::cout << "Can't open file\n";
			return -1;
		}
		cout<<"Working on "<<argv[file]<<endl;
		TTree *eventHkTree; 
		eventHkTree= (TTree*) fp->Get("eventHkTree");
		if(!eventHkTree) {
			std::cout << "Can't find eventHkTree\n";
			return -1;
		}
		int runNum;
		eventHkTree->SetBranchAddress("run",&runNum);
		eventHkTree->SetBranchAddress("eventHk",&eventHkPtr);

		int numEntries = eventHkTree->GetEntries();
		// numEntries=100;

		for(int event=0; event<numEntries; event++){
			eventHkTree->GetEntry(event);
			// thresh.push_back(double(eventHkPtr->getSingleChannelThreshold(0,0)));
			// unixTime.push_back(double(eventHkPtr->unixTime));
			// printf("Event %d, thresh %f \n", double(eventHkPtr->getSingleChannelThreshold(0,0)));
			// cout<<"Type is "<<typeid(eventHkPtr->getSingleChannelThreshold(0,0))<<endl;
			int thing = eventHkPtr->getSingleChannelThreshold(0,0);
			int unixTime = double(eventHkPtr->unixTime);
			// thresh.push_back(double(thing));
			// unixTime.push_back(double(eventHkPtr->unixTime));
			rate->Fill(double(unixTime),double(thing));
			cout<<"Threshold is "<<thing<<endl;
			// rate->Fill(runNum,thing);
			// for(int dda=0; dda<4; dda++){
			// 	printf("DDA %d\n",dda );
			// 	for(int chan=0; chan<8; chan++){
			// 		printf("     Chan %d: %d \n",dda, eventHkPtr->getSingleChannelThreshold(dda,chan));
			// 	}
			// }
		}

		if(runNum==3419){
			eventHkTree->GetEntry(numEntries-2);
			runEndunixTime=eventHkPtr->unixTime;
		}

		fp->Close();
	}

	TLine *line = new TLine(1380750395,0,1380750395,31000);
	TLine *line_run_end = new TLine(runEndunixTime,0,runEndunixTime,31000);
	// TLine *cut_start = new TLine(1413550800,0,1413550800,31000);
	// TLine *cut_end = new TLine(1413552600,0,1413552600,31000);


	rate->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
	// rate->GetXaxis()->SetTimeFormat("%m/%d %H");
	// rate->GetXaxis()->SetTimeFormat("%m/%d");
	rate->GetXaxis()->SetTimeOffset(0.,"GMT");
	gStyle->SetOptStat(0);

	// TGraph *gr = new TGraph(unixTime.size(), &unixTime[0], &thresh[0]);
	// gr->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
	// gr->GetXaxis()->SetTimeFormat("%m/%d %H");
	// gr->GetXaxis()->SetTimeOffset(0.,"GMT");
	TCanvas *c = new TCanvas("","",1100,850);
	rate->Draw("colz");
	rate->GetXaxis()->SetTitle("Date (Month/Day)");
	// rate->GetXaxis()->SetTitle("runNum");
	rate->GetYaxis()->SetTitle("L1 Trigger Threshold");
	line->Draw("same");
		line->SetLineWidth(2);
		line->SetLineStyle(9);
	line_run_end->Draw("same");
		line_run_end->SetLineWidth(2);
		line_run_end->SetLineStyle(9);
		line_run_end->SetLineColor(kBlue);
	// cut_start->Draw("same");
	// 	cut_start->SetLineWidth(2);
	// 	cut_start->SetLineStyle(9);
	// 	cut_start->SetLineColor(kGreen);
	// cut_end->Draw("same");
	// 	cut_end->SetLineWidth(2);
	// 	cut_end->SetLineStyle(9);
	// 	cut_end->SetLineColor(kGreen);
	// gr->Draw("AP");
	// gr->SetMarkerStyle(kOpenCircle);
	c->SaveAs("thresh_vs_time_run1569_CWevent.png");
}
