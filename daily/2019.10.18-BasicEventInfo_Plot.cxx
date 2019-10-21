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
		std::cout << "Usage\n" << argv[0] << " <station> <config> <date_or_run_mode> <basic_info_file> "<<endl;
		return -1;
	}
	int station = atoi(argv[1]);
	int config = atoi(argv[2]);
	int date_or_run_mode = atoi(argv[3]);

	vector<int> BadRunList=BuildBadRunList(station);


	/*
		Place to set up all the plots we want
	*/

	TDatime start(2013, 01, 01, 00, 00,0);
	TDatime stop(2013, 12, 31, 24, 00,0);
	if(config==3){
		TDatime start(2015, 01, 01, 00, 00,0);
		TDatime stop(2015, 12, 31, 24, 00,0);
	}

	if(config==5){
		TDatime start(2014, 01, 01, 00, 00,0);
		TDatime stop(2014, 12, 31, 24, 00,0);
	}
	int start_bin = start.Convert();
	int stop_bin = stop.Convert();
	int numBins = 365;

	if(date_or_run_mode==1){
		printf(BLUE"RUN MODE\n"RESET);
		start_bin=0;
		stop_bin=8000;
		numBins=8000;
	}

	TH2D *h2_rms_vs_time[16];
	for(int i=0; i<16; i++){
		stringstream ss1;
		ss1<<"Channel "<<i;
		h2_rms_vs_time[i] = new TH2D(ss1.str().c_str(),ss1.str().c_str(),numBins, start_bin, stop_bin, 1000,0,1000);
		if(date_or_run_mode==0){
			h2_rms_vs_time[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
			h2_rms_vs_time[i]->GetXaxis()->SetTimeFormat("%y/%m");
			// h2_rms_vs_time[i]->GetXaxis()->SetNdivisions(12,0,0,false);
		}
	}

	TH1D *h1_all_events = new TH1D("all_events","all_events",numBins, start_bin,stop_bin);
	TH1D *h1_soft_events = new TH1D("soft_events","soft_events",numBins,start_bin,stop_bin);
	TH1D *h1_tagcal_events = new TH1D ("tagged_cal_events","tagged_cal_events",numBins,start_bin,stop_bin);
	TH1D *h1_soft_fraction = new TH1D("soft_event_fraction","soft_event_fraction",numBins,start_bin,stop_bin);
	TH1D *h1_tagcal_fraction = new TH1D("tagged_cal_fraction","tagged_cal_fraction",numBins,start_bin,stop_bin);
	if(date_or_run_mode==0){
		prepTimeAxis(h1_all_events);
		prepTimeAxis(h1_soft_events);
		prepTimeAxis(h1_tagcal_events);
		prepTimeAxis(h1_soft_fraction);
		prepTimeAxis(h1_tagcal_fraction);
	}

	int num_total=0;

	for(int file_num=4; file_num<argc; file_num++){

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

		// check this again as the bad run list evolves
		bool secondBadRunCheck = isBadRun(station,runNum,BadRunList);

		if(isKnownBadRun || secondBadRunCheck){
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

			num_total++;

			int var_to_fill;
			if(date_or_run_mode==0) {
				var_to_fill=unixTime;
			}
			else if(date_or_run_mode==1){
				var_to_fill=runNum;
			}

			h1_all_events->Fill(var_to_fill);
			if(isCal){
				h1_tagcal_events->Fill(var_to_fill);
			}
			if(isSoft){
				h1_soft_events->Fill(var_to_fill);
			}

			for(int chan=0; chan<16; chan++){
				h2_rms_vs_time[chan]->Fill(var_to_fill,deepChannelRMS[chan]);
			}
		}
		fpIn->Close();
		delete fpIn;
	} //end loop over input files

	printf(GREEN"TIME TO PLOT\n"RESET);

	double low_edge;
	double high_edge;
	if(config==1){
		low_edge=1400;
		high_edge=2000;
	}
	if(config==2){
		low_edge=500;
		high_edge=1400;
	}
	if(config==3){
		low_edge=3000;
		high_edge=8000;
	}
	if(config==4){
		low_edge=6000;
		high_edge=8000;
	}
	if(config==5){
		low_edge=1800;
		high_edge=3200;
	}

	TCanvas *c = new TCanvas("","",8*850,4*850);
	c->Divide(4,4);
	for(int i=0; i<16; i++){
		c->cd(i+1);
		h2_rms_vs_time[i]->Draw("colz");
		h2_rms_vs_time[i]->GetYaxis()->SetTitle("RMS");
		h2_rms_vs_time[i]->GetXaxis()->SetTitle("unixTime [YY/MM]");
		if(date_or_run_mode==1){
			h2_rms_vs_time[i]->GetXaxis()->SetTitle("Run Number");
			h2_rms_vs_time[i]->GetXaxis()->SetRangeUser(low_edge,high_edge);
		}
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
	sprintf(save_plot_title,"%s/basic_info/%d.%d.%d_RMS_vs_Time_A%d_c%d.png",plotPath,year_now,month_now,day_now,station,config);
	if(date_or_run_mode){
		sprintf(save_plot_title,"%s/basic_info/%d.%d.%d_RMS_vs_Run_A%d_c%d.png",plotPath,year_now,month_now,day_now,station,config);
	}
	c->SaveAs(save_plot_title);
	delete c;

	// make fractions
	for(int bin=1; bin<h1_all_events->GetNbinsX(); bin++){
		int total_events = h1_all_events->GetBinContent(bin);
		if(total_events>0){
			int soft_events  = h1_soft_events->GetBinContent(bin);
			double soft_frac = double(soft_events)/double(total_events);
			int cal_events = h1_tagcal_events->GetBinContent(bin);
			double tagcal_frac = double(cal_events)/double(total_events);
			h1_soft_fraction->SetBinContent(bin,soft_frac);
			h1_tagcal_fraction->SetBinContent(bin,tagcal_frac);
		}
		else{
			h1_soft_fraction->SetBinContent(bin,-0.1);
			h1_tagcal_fraction->SetBinContent(bin,-0.1);
		}
	}
// mode = ksiourmen  (default = 000001111)
// k = 1;  kurtosis printed
// k = 2;  kurtosis and kurtosis error printed

// s = 1;  skewness printed
// s = 2;  skewness and skewness error printed

// i = 1;  integral of bins printed

// o = 1;  number of overflows printed

// u = 1;  number of underflows printed

// r = 1;  rms printed
// r = 2;  rms and rms error printed

// m = 1;  mean value printed
// m = 2;  mean and mean error values printed

// e = 1;  number of entries printed

// n = 1;  name of histogram is printed


	gStyle->SetOptStat(101100);

	// print out the "high software rate" runs to a file and make a histogram
	TH1D *h1_dist_soft_frac = new TH1D("dist_soft_frac","dist_soft_frac",300,0,30);

	char title_txt[200];
	sprintf(title_txt,"%s/basic_info/soft_fraction_A%d_c%d.txt",plotPath,station,config);
	FILE *fout = fopen(title_txt, "a");
	for(int bin=1; bin<h1_soft_fraction->GetNbinsX(); bin++){
		double fraction = h1_soft_fraction->GetBinContent(bin);
		h1_dist_soft_frac->Fill(fraction*100.);
		if(fraction>.19){
			// printf("Run %4d, Rate %.2f \n", bin-1, fraction);
			fprintf(fout,"Run %4d, Rate %.2f \n",bin-1, fraction);
		}
	}
	fclose(fout);//close file

	TCanvas *cDistSoftFrac = new TCanvas("","",1100,850);
	h1_dist_soft_frac->Draw("");
		h1_dist_soft_frac->GetXaxis()->SetTitle("Fraction of Events that are Software Triggers [%]");
		h1_dist_soft_frac->GetYaxis()->SetTitle("Number of Runs");
		h1_dist_soft_frac->SetLineWidth(2);
		gPad->SetLogy();
	sprintf(save_plot_title,"%s/basic_info/%d.%d.%d_DistSoftFraction_A%d_c%d.png",plotPath,year_now,month_now,day_now,station,config);
	cDistSoftFrac->SaveAs(save_plot_title);
	delete cDistSoftFrac;

	gStyle->SetOptStat(0);

	TLine *line;
	line = new TLine(start_bin,0,stop_bin,0);
	line->SetLineColor(kBlack);
	line->SetLineStyle(9);

	TCanvas *c2 = new TCanvas("","",2*850,2*850);
	c2->Divide(1,2);
	c2->cd(1);
		h1_soft_fraction->Draw("");
		h1_soft_fraction->GetXaxis()->SetTitle("unixTime [YY/MM]");
		h1_soft_fraction->GetYaxis()->SetTitle("Fraction of Software Triggers");
		if(date_or_run_mode==1){
			h1_soft_fraction->GetXaxis()->SetTitle("Run Number");
			h1_soft_fraction->GetXaxis()->SetRangeUser(low_edge,high_edge);
		}
		h1_soft_fraction->GetYaxis()->SetRangeUser(-0.2,1.2);
		line->Draw("sameL");
	c2->cd(2);
		h1_tagcal_fraction->Draw("");
		h1_tagcal_fraction->GetXaxis()->SetTitle("unixTime [YY/MM]");
		h1_tagcal_fraction->GetYaxis()->SetTitle("Fraction of Tagged Calpulsers");
		if(date_or_run_mode==1){
			h1_tagcal_fraction->GetXaxis()->SetTitle("Run Number");
			h1_tagcal_fraction->GetXaxis()->SetRangeUser(low_edge,high_edge);
		}
		h1_tagcal_fraction->GetYaxis()->SetRangeUser(-0.2,1.2);
		line->Draw("sameL");
	sprintf(save_plot_title,"%s/basic_info/%d.%d.%d_SoftCalFraction_vs_Time_A%d_c%d.png",plotPath,year_now,month_now,day_now,station,config);
	if(date_or_run_mode){
		sprintf(save_plot_title,"%s/basic_info/%d.%d.%d_SoftCalFraction_vs_Run_A%d_c%d.png",plotPath,year_now,month_now,day_now,station,config);
	}
	c2->SaveAs(save_plot_title);
	delete c2;


	// plot together sometimes
	TCanvas *c3 = new TCanvas("","",2*850,3*850);
	c3->Divide(1,3);
	c3->cd(1);
		h2_rms_vs_time[0]->Draw("colz");
	c3->cd(3);
		h1_tagcal_fraction->Draw("");
		line->Draw("sameL");
	c3->cd(2);
		h1_soft_fraction->Draw("");
		line->Draw("sameL");
	sprintf(save_plot_title,"%s/basic_info/%d.%d.%d_RMS_and_SoftCalFrac_vs_Time_A%d_c%d.png",plotPath,year_now,month_now,day_now,station,config);
	if(date_or_run_mode){
		sprintf(save_plot_title,"%s/basic_info/%d.%d.%d_RMS_and_SoftCalFrac_vs_Run_A%d_c%d.png",plotPath,year_now,month_now,day_now,station,config);
	}
	c3->SaveAs(save_plot_title);
	delete c3;
}