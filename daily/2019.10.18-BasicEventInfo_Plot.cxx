////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  2019.10.17-BasicEventInfo.cxx 
////  store the basic event information
////
////  Oct 2019
////  basic diagnostic plotting bonanza
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <string>
#include <sstream>

//AraRoot Includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraQualCuts.h"
#include "FFTtools.h"

#include "tools_PlottingFns.h"
#include "tools_WaveformFns.h"
#include "tools_Cuts.h"
#include "tools_CommandLine.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TTimeStamp.h"
#include "TH2D.h"
#include "TH1D.h"

using namespace std;

int main(int argc, char **argv)
{
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;


	gStyle->SetOptStat(0);
	gStyle->SetTimeOffset(0);

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

	if(argc<4) {
		std::cout << "Usage\n" << argv[0] << " <station> <config> <basic_info_file> "<<endl;
		return -1;
	}
	int station = atoi(argv[1]);
	int config = atoi(argv[2]);

	/*
		Place to set up all the plots we want
	*/

	TDatime start(2013, 01, 01, 00, 00,0);
	int start_bin = start.Convert();
	TDatime stop(2013, 12, 31, 24, 00,0);
	int stop_bin = stop.Convert();

	TH2D *h2_rms_vs_time[16];
	for(int i=0; i<16; i++){
		stringstream ss1;
		ss1<<"Channel "<<i;
		h2_rms_vs_time[i] = new TH2D(ss1.str().c_str(),ss1.str().c_str(),365, start_bin, stop_bin, 500,0,500);
		h2_rms_vs_time[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		h2_rms_vs_time[i]->GetXaxis()->SetTimeFormat("%y/%m");
		// h2_rms_vs_time[i]->GetXaxis()->SetNdivisions(12,0,0,false);
	}

	int num_total=0;

	for(int file_num=3; file_num<argc; file_num++){

		cout<<argv[file_num]<<endl;

		TFile *fpIn = TFile::Open(argv[file_num]);
		if(!fpIn) { std::cout << "Can't open file\n"; return -1; }
		TTree *inTree = (TTree*) fpIn->Get("outTree");
		if(!inTree){ cout<<"Can't open filter tree"<<endl; return -1; }
		bool isCal;
		bool isSoft;
		bool isShort;
		bool isSuperShort;
		bool hasSpareChannelIssuev1;
		bool hasSpareChannelIssuev2;
		bool hasDigitizerError;
		bool isKnownBadLivetime;
		inTree->SetBranchAddress("isCal", &isCal);
		inTree->SetBranchAddress("isSoft", &isSoft);
		inTree->SetBranchAddress("isShort", &isShort);
		inTree->SetBranchAddress("isSuperShort", &isSuperShort);
		inTree->SetBranchAddress("hasSpareChannelIssuev1", &hasSpareChannelIssuev1);
		inTree->SetBranchAddress("hasSpareChannelIssuev2", &hasSpareChannelIssuev2);
		inTree->SetBranchAddress("hasDigitizerError", &hasDigitizerError);
		inTree->SetBranchAddress("isKnownBadLivetime", &isKnownBadLivetime);

		int unixTime;
		int unixTimeUs;
		int timeStamp;
		int eventNumber;
		double deepChannelRMS[16];
		double spareChannelRMS[4];
		inTree->SetBranchAddress("unixTime",&unixTime);
		inTree->SetBranchAddress("unixTimeUs",&unixTimeUs);
		inTree->SetBranchAddress("timeStamp",&timeStamp);
		inTree->SetBranchAddress("eventNumber",&eventNumber);
		inTree->SetBranchAddress("deepChannelRMS",&deepChannelRMS);
		inTree->SetBranchAddress("spareChannelRMS", &spareChannelRMS);

		int runNum;	
		bool isKnownBadRun;
		inTree->SetBranchAddress("runNum",&runNum);
		inTree->SetBranchAddress("isKnownBadRun", &isKnownBadRun);

		int numEntries = inTree->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;
		inTree->GetEvent(0);

		if(isKnownBadRun){
			fpIn->Close();
			delete fpIn;
			continue;
		}

		//now to loop over events
		for(int event=0; event<numEntries; event++){
			inTree->GetEvent(event);
			num_total++;
			if(isKnownBadLivetime){
				continue;
			}
			for(int chan=0; chan<16; chan++){
				h2_rms_vs_time[chan]->Fill(unixTime,deepChannelRMS[chan]);
			}
		}
		fpIn->Close();
		delete fpIn;
	} //end loop over input files

	printf(GREEN"TIME TO PLOT\n"RESET);

	TCanvas *c = new TCanvas("","",8*850,4*850);
	c->Divide(4,4);
	for(int i=0; i<16; i++){
		c->cd(i+1);
		h2_rms_vs_time[i]->Draw("colz");
		h2_rms_vs_time[i]->GetYaxis()->SetTitle("RMS");
		h2_rms_vs_time[i]->GetXaxis()->SetTitle("unixTime [YY/MM]");
		h2_rms_vs_time[i]->GetXaxis()->SetTitleOffset(1.1);
		h2_rms_vs_time[i]->GetYaxis()->SetTitleOffset(1.1);
		h2_rms_vs_time[i]->GetZaxis()->SetTitleOffset(1.1);
		h2_rms_vs_time[i]->GetXaxis()->SetTitleSize(0.06);
		h2_rms_vs_time[i]->GetYaxis()->SetTitleSize(0.06);
		h2_rms_vs_time[i]->GetZaxis()->SetTitleSize(0.06);
		h2_rms_vs_time[i]->GetXaxis()->SetLabelSize(0.06);
		h2_rms_vs_time[i]->GetYaxis()->SetLabelSize(0.06);
		h2_rms_vs_time[i]->GetZaxis()->SetLabelSize(0.06);
		// hRMS[i]->GetXaxis()->SetRangeUser(10,40);
		// if(i+1==1){
		// 	TLegend *leg = new TLegend(0.48,0.6,0.9,0.9);
		// 	leg->AddEntry(hRMS[i],"All Events","l");
		// 	leg->AddEntry(hRMS_cal[i],"Tagged Cal Pulsers","l");
		// 	// leg->AddEntry(hRMS_stragglers[i],"Straggling Events","l");
		// 	leg->Draw();
		// }
	}
	char save_plot_title[400];
	sprintf(save_plot_title,"%s/basic_info/%d.%d.%d_RMS_vs_Time_A%d_c%d_%dEvents.png",plotPath,year_now,month_now,day_now,station,config,num_total);
	c->SaveAs(save_plot_title);
	delete c;
}