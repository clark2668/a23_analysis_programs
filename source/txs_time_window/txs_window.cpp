////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	txs_window.cxx
//// 	plot ARA events relative to TXS
////////////////////////////////////////////////////////////////////////////////

//C++
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <sys/stat.h>

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTimeStamp.h"
#include "TH1D.h"
#include "TLine.h"
#include "TF1.h"

using namespace std;

int getThisRunNum(char *thearg){
	string chRun = "run";
	string file = string(thearg);
	size_t foundRun = file.find(chRun);
	string strRunNum = file.substr(foundRun+4,4);
	int runNum = atoi(strRunNum.c_str());
	return runNum;
}

int main(int argc, char **argv)
{

	if(argc<4){
		cout<< "Usage\n" << argv[0] << " <station> <output_location> <filename>"<<endl;;
		return -1;
	}

	int station = atoi(argv[1]);
	string output_location = argv[2];


	TTimeStamp start;
	TTimeStamp stop;
	start.Set(2014, 7, 00, 00, 00,0,0,true,0);
	stop.Set(2015, 6, 00, 00, 00,0,0,true,0);

	int start_bin = start.GetSec();
	int stop_bin = stop.GetSec();

	int num_bins=(stop_bin-start_bin)/60/60/24; // convert to days


	TH1D *h1_ev_vs_time = new TH1D("","Ev vs Time",num_bins,start_bin,stop_bin);
	h1_ev_vs_time->GetXaxis()->SetTimeDisplay(1);
	// h1_ev_vs_time->GetXaxis()->SetTimeFormat("%y/%m");
	h1_ev_vs_time->GetXaxis()->SetTimeFormat("%b '%y");
	h1_ev_vs_time->GetXaxis()->SetTimeOffset(0.,"GMT");

	for(int file_num=3; file_num<argc; file_num++){

		cout<<"On file "<<argv[file_num]<<endl;

		// int runNum=getThisRunNum(argv[file_num]);
		// cout<<"runNum is "<<runNum<<endl;
		
		TFile *inputFile = TFile::Open(argv[file_num]);
		if(!inputFile){
			cout<<"Can't open joined file!"<<endl;
			continue;
			// return -1;
		}
		TTree *theTree = (TTree*) inputFile->Get("OutputTree");
		if(!theTree){
			cout<<"Can't get eventTree!"<<endl;
			inputFile->Close();
			continue;
			// return -1;
		}

		int unixTime;
		theTree->SetBranchAddress("unixTime",&unixTime);

		int numEntries = theTree->GetEntries();

		theTree->GetEntry(0);
		int first_event_unixtime = unixTime;
		h1_ev_vs_time->Fill(first_event_unixtime);

		//now to loop over events
		// for(int event=0; event<numEntries; event++){
		// 	// if(event%10000==0) cout<<"on event "<<event<<endl;
		// 	// cout<<"On event"<<event<<endl;
		// 	theTree->GetEntry(event);
		// 	h1_ev_vs_time->Fill(unixTime);
		// }

		inputFile->Close();
		delete inputFile;
	}

	// go through the histogram
	// and just make the bin "1" if we have data there
	for(int bin=0; bin<h1_ev_vs_time->GetNbinsX(); bin++){
		double count = h1_ev_vs_time->GetBinContent(bin);
		if(count>0) count=1.;
		h1_ev_vs_time->SetBinContent(bin,count);
	}

	// normalize the data real quick
	// h1_ev_vs_time->Sumw2();
	// Double_t scale = 1/h1_ev_vs_time->Integral("width");
	// h1_ev_vs_time->Scale(scale);
	
	TTimeStamp txs_center_gaus;
	txs_center_gaus.Set(2014,12,13,00,00,0,0,true,0);
	int txs_center_gaus_time = txs_center_gaus.GetSec();
	TLine gaus_time(txs_center_gaus_time, 0, txs_center_gaus_time, 1);
	double txs_gaus_duration_days = 55;
	double txs_gaus_plot_range_start = txs_center_gaus_time - (200*24*60*60);
	double txs_gaus_plot_range_stop = txs_center_gaus_time + (200*24*60*60);

	TTimeStamp txs_center_box;
	txs_center_box.Set(2014,12,26,00,00,0,0,true,0);
	int txs_center_box_time = txs_center_box.GetSec();
	TLine box_time(txs_center_box_time, 0, txs_center_box_time, 1);
	double txs_box_duration_days = 79;
	double txs_box_plot_range_start = txs_center_box_time - (txs_box_duration_days*24*60*60);
	double txs_box_plot_range_stop = txs_center_box_time + (txs_box_duration_days*24*60*60);
	TLine box_start(txs_box_plot_range_start, 0, txs_box_plot_range_start, 1.);
	TLine box_stop(txs_box_plot_range_stop, 0, txs_box_plot_range_stop, 1.);

	TF1 *f1 = new TF1("f1","TMath::Gaus(x,[0],[1],0)",txs_gaus_plot_range_start,txs_gaus_plot_range_stop);
	f1->SetParameters(double(txs_center_gaus_time),txs_gaus_duration_days*24*60*60);
	f1->GetXaxis()->SetTimeDisplay(1);
	f1->GetXaxis()->SetTimeFormat("%y/%m/%d");
	f1->GetXaxis()->SetTimeOffset(0.,"GMT");

	gStyle->SetOptStat(0);

	TCanvas *c = new TCanvas("","",1100,850);
	h1_ev_vs_time->Draw("");
	// h1_ev_vs_time->GetYaxis()->SetTitle("ARA Has Data");
	// h1_ev_vs_time->SetFillColor(kBlue-10);
	h1_ev_vs_time->SetFillColor(kGray);
	h1_ev_vs_time->SetLineColor(kGray);
	h1_ev_vs_time->SetTitle("A3");

	gaus_time.Draw("same");
	gaus_time.SetLineColor(kRed);
	gaus_time.SetLineWidth(2.);
	box_time.Draw("same");
	box_time.SetLineColor(kBlue);
	box_time.SetLineWidth(2.);


	f1->Draw("same");

	box_start.Draw("same");
	box_start.SetLineColor(kBlue);
	box_start.SetLineWidth(2.);
	box_stop.Draw("same");
	box_stop.SetLineColor(kBlue);	
	box_stop.SetLineWidth(2.);

	char title[150];
	sprintf(title,"%s/events_vs_time_w_txs_A%d.png",output_location.c_str(), station);
	c->SaveAs(title);
	delete c;

}
