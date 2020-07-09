////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	plot_sim_RPR.cxx
////
////	April 2020,  baclark@msu.edu
////	plot distributions of RPR
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
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"

using namespace std;

int main(int argc, char **argv)
{
	if(argc<2) {  // Check to make sure there are enough arguments to do something meaningful
		std::cout << "Usage requires you to provide input parameter of the form " << basename(argv[0]) << " <rpr_file_1> <rpr_file_2> ..." << std::endl;
		return -1;
	}

	TH1D *h1_num_hit_E17 = new TH1D("E17","E17",17,-0.5,16.5);
	TH1D *h1_num_hit_E18 = new TH1D("E18","E18",17,-0.5,16.5);
	TH1D *h1_num_hit_E19 = new TH1D("E19","E19",17,-0.5,16.5);

	TH1D *h1_time_to_hits[16];
	stringstream ss;
	for(int i=0; i<16; i++){
		ss.str("");
		ss<<"Channel "<<i;
		h1_time_to_hits[i] = new TH1D(ss.str().c_str(), ss.str().c_str(), 50, -250, 250);
	}

	for(int file=1; file<argc; file++){
		TFile *fpIn = new TFile(argv[file], "OLD");
		if(!fpIn){ cout<<"Cannnot open file "<<argv[file]<<endl; delete fpIn; continue;}
		fpIn->cd();
		TTree *inTree = (TTree*) fpIn->Get("outTree");
		if(!inTree){ cout<<"Cannot get outTree in file "<<argv[file]<<endl; delete fpIn; continue;}
		double weight;
		double RPR[16];
		double pnu;
		double times[16][3];
		inTree->SetBranchAddress("weight",&weight);
		inTree->SetBranchAddress("RPR",&RPR);
		inTree->SetBranchAddress("pnu",&pnu);
		inTree->SetBranchAddress("times",&times);
		int numEntries = inTree->GetEntries();
		for(int event=0; event<numEntries; event++){
			inTree->GetEntry(event);
			int numHits=0;
			for(int i=0; i<16; i++){
				if(RPR[i]>=8.) numHits++;
			}
			int logE = int(TMath::Log10(pnu));
			if(logE==17){
				h1_num_hit_E17->Fill(numHits,weight);
			}
			if(logE==18){
				h1_num_hit_E18->Fill(numHits,weight);
			}
			if(logE==19){
				h1_num_hit_E19->Fill(numHits,weight);
			}
			if(numHits<6){
				for(int i=0; i<16; i++){
					if(RPR[i]>=8.){
						double time_from_start = times[i][0] - times[i][2];
						double time_to_stop = times[i][1] - times[i][2];
						double which_is_smaller = time_to_stop;
						if(abs(time_from_start)<abs(time_to_stop)) which_is_smaller=time_from_start;
						// cout<<"which is larger "<<which_is_larger<<endl;
						// printf("Chan %d: t-tstart is %.2f, tstop-t %.2f, closest approach %.2f \n", time_from_start, time_to_stop, which_is_smaller);
						h1_time_to_hits[i]->Fill(which_is_smaller, weight);
						// h1_time_to_hits[i]->Fill(times[i][2], weight);
					}
				}
			}
		}
		fpIn->Close();
		delete fpIn;
	}

	gStyle->SetOptStat(0);

	// plots of delays between hits and end of the waveform

	TCanvas *c_hit_times = new TCanvas("c_hit_times","c_hit_times",4*1100,4*850);
	c_hit_times->Divide(4,4);
	for(int i=0; i<16; i++){
		c_hit_times->cd(i+1);
		h1_time_to_hits[i]->Draw("hist");
		h1_time_to_hits[i]->GetXaxis()->SetTitle("Time to Nearest Waveform Edge (ns)");
		h1_time_to_hits[i]->GetYaxis()->SetTitle("Weighted Events");
		h1_time_to_hits[i]->SetLineWidth(3);
	}
	char save_title[150];
	sprintf(save_title,"dist_hit_times.png");
	c_hit_times->SaveAs(save_title);



	// plots of histogram of hits

	h1_num_hit_E19->Sumw2();
	Double_t scale = 1/h1_num_hit_E19->Integral("width");
	h1_num_hit_E19->Scale(scale);
	
	h1_num_hit_E18->Sumw2();
	scale = 1/h1_num_hit_E18->Integral("width");
	h1_num_hit_E18->Scale(scale);

	h1_num_hit_E17->Sumw2();
	scale = 1/h1_num_hit_E17->Integral("width");
	h1_num_hit_E17->Scale(scale);

	TH1 *h1_cumulative_hit_E19 = h1_num_hit_E19->GetCumulative();
	TH1 *h1_cumulative_hit_E18 = h1_num_hit_E18->GetCumulative();
	TH1 *h1_cumulative_hit_E17 = h1_num_hit_E17->GetCumulative();

	sprintf(save_title,"dist_hits.png");
	TCanvas *c = new TCanvas("","",2*1100,850);
	c->Divide(2,1);
	c->cd(1);
		h1_num_hit_E19->Draw("hist");
			h1_num_hit_E19->SetLineWidth(2);
			h1_num_hit_E19->SetTitle("Number of hits w/ RPR>=8; Number of Hits; Normalized Weighted Events");
			h1_num_hit_E19->SetLineColor(kBlue);
		h1_num_hit_E18->Draw("histsame");
			h1_num_hit_E18->SetLineWidth(2);
			h1_num_hit_E18->SetLineColor(kRed);
		h1_num_hit_E17->Draw("histsame");
			h1_num_hit_E17->SetLineWidth(2);
			h1_num_hit_E17->SetLineColor(kGreen);
		TLegend *legend = new TLegend(0.6,0.5,0.8,0.7);
			legend->AddEntry(h1_num_hit_E19,"E19","l");
			legend->AddEntry(h1_num_hit_E18,"E18","l");
			legend->AddEntry(h1_num_hit_E17,"E17","l");
			// legend->SetBorderSize(0);  //no border for legend
			legend->SetFillColor(0);  //fill color is white
			legend->Draw("same");
	c->cd(2);
		h1_cumulative_hit_E19->Draw("hist");
			h1_cumulative_hit_E19->SetLineWidth(2);
			h1_cumulative_hit_E19->SetTitle("  ; Number of Hits; Cumulative Fraction");
			h1_cumulative_hit_E19->SetLineColor(kBlue);		
		h1_cumulative_hit_E18->Draw("histsame");
			h1_cumulative_hit_E18->SetLineWidth(2);
			h1_cumulative_hit_E18->SetLineColor(kRed);
		h1_cumulative_hit_E17->Draw("histsame");
			h1_cumulative_hit_E17->SetLineWidth(2);
			h1_cumulative_hit_E17->SetLineColor(kGreen);
	c->SaveAs(save_title);


}//close the main program