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

void prepTimeAxis(TH1D *h1){
	h1->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
	h1->GetXaxis()->SetTimeFormat("%y/%m");
	h1->GetXaxis()->SetTimeOffset(0.,"GMT");
}

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

	if(argc<5) {
		std::cout << "Usage\n" << argv[0] << " <station> <config> <is_simulation> <basic_info_file> "<<endl;
		return -1;
	}
	int station = atoi(argv[1]);
	int config = atoi(argv[2]);
	int is_simulation = atoi(argv[3]);

	vector<int> BadRunList=BuildBadRunList(station);

	/*
		Place to set up all the plots we want
	*/

	int numBins=500;
	double start_bin=0.;
	double stop_bin=5.;
	TH1D *h1_all_events = new TH1D("all_events","all_events",numBins, start_bin,stop_bin);
	TH1D *h1_tagcal_events = new TH1D ("tagged_cal_events","tagged_cal_events",numBins,start_bin,stop_bin);

	int num_total=0;

	for(int file_num=4; file_num<argc; file_num++){

		cout<<argv[file_num]<<endl;

		TFile *fpIn = TFile::Open(argv[file_num]);
		if(!fpIn) { std::cout << "Can't open file\n"; fpIn->Close(); continue; /*return -1;*/ }
		TTree *inTree = (TTree*) fpIn->Get("outTree");
		if(!inTree){ cout<<"Can't open filter tree"<<endl; fpIn->Close();  continue; /*return -1;*/ }
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

		double powerPerChannel[16];
		double powerPerString[4];
		double powerRatio;
		inTree->SetBranchAddress("powerPerChannel",&powerPerChannel);
		inTree->SetBranchAddress("powerPerString", &powerPerString);	
		inTree->SetBranchAddress("powerRatio", &powerRatio);	

		int runNum;	
		bool isKnownBadRun;
		inTree->SetBranchAddress("runNum",&runNum);
		inTree->SetBranchAddress("isKnownBadRun", &isKnownBadRun);

		double weight_out;
		if(is_simulation)
			inTree->SetBranchAddress("weight_out",&weight_out);

		int numEntries = inTree->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;
		inTree->GetEvent(0);

		// check this again as the bad run list evolves
		bool secondBadRunCheck = isBadRun(station,runNum,BadRunList);

		// if(isKnownBadRun || secondBadRunCheck){
		if(isKnownBadRun){
			fpIn->Close();
			delete fpIn;
			continue;
		}

		//now to loop over events
		for(int event=0; event<numEntries; event++){
			inTree->GetEvent(event);

			if(isKnownBadLivetime){
				continue;
			}
			if(hasSpareChannelIssuev1 || hasSpareChannelIssuev2) continue;

			if(!is_simulation)
				weight_out=1.;

			num_total++;

			vector<double> average_of_antennas_on_string;
			average_of_antennas_on_string.push_back((powerPerChannel[0]+powerPerChannel[4]+powerPerChannel[8]+powerPerChannel[12])/4.);
			average_of_antennas_on_string.push_back((powerPerChannel[1]+powerPerChannel[5]+powerPerChannel[9]+powerPerChannel[13])/4.);
			average_of_antennas_on_string.push_back((powerPerChannel[2]+powerPerChannel[6]+powerPerChannel[10]+powerPerChannel[14])/4.);
			average_of_antennas_on_string.push_back((powerPerChannel[3]+powerPerChannel[7]+powerPerChannel[11])/3.);
			std::sort(average_of_antennas_on_string.begin(), average_of_antennas_on_string.end()); //sort smallest to largest
			double ratio = average_of_antennas_on_string[3]/average_of_antennas_on_string[2];

			if(ratio>10.){
				// printf("Run %d, Event %d \n", runNum, eventNumber);
				// printf("Run %d, Event %d, Glitch Detect 1 %d, Glitch Detect 2 %d \n", runNum, eventNumber, hasSpareChannelIssuev1, hasSpareChannelIssuev2);
			}

			h1_all_events->Fill(ratio,weight_out);

			if(isCal)
				h1_tagcal_events->Fill(ratio,weight_out);
		}
		fpIn->Close();
		delete fpIn;
	} //end loop over input files

	printf(GREEN"TIME TO PLOT\n"RESET);


	double total=0;
	total = h1_all_events->Integral();
	int binMin = h1_all_events->FindBin(5.);
	int binMax = h1_all_events->FindBin(50);
	double above_cut = h1_all_events->Integral(binMin,binMax);
	printf("Fraction above is %.4f/%.4f = %.5f \n", above_cut,total,above_cut/total);

	// gStyle->SetOptStat(101100);

	// print out the "high software rate" runs to a file and make a histogram
	TH1D *h1_dist_soft_frac = new TH1D("dist_soft_frac","dist_soft_frac",300,0,30);

	TCanvas *c = new TCanvas("","",1100,850);
	h1_all_events->Draw("");
		h1_all_events->GetXaxis()->SetTitle("Ratio");
		h1_all_events->GetYaxis()->SetTitle("Number of Events");
		h1_all_events->SetLineWidth(2);
		gPad->SetLogy();
	h1_tagcal_events->Draw("same");
		h1_tagcal_events->SetLineColor(kRed);
		h1_tagcal_events->SetLineWidth(2);
	if(is_simulation)
		h1_all_events->GetYaxis()->SetRangeUser(0.5,1e4);
	char save_plot_title[500];
	sprintf(save_plot_title,"%s/high_power_events/%d.%d.%d_Ratio_AveragePower_A%d_c%d_sim%d.png",plotPath,year_now,month_now,day_now,station,config,is_simulation);
	c->SaveAs(save_plot_title);
	delete c;

}