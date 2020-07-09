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

#include "AraGeomTool.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"

using namespace std;

int main(int argc, char **argv)
{
	if(argc<4) {  // Check to make sure there are enough arguments to do something meaningful
		std::cout << "Usage requires you to provide input parameter of the form " << basename(argv[0]) << "<station> <config> <rpr_file_1> <rpr_file_2> ..." << std::endl;
		return -1;
	}

	int station=atoi(argv[1]);
	int config=atoi(argv[2]);
	AraGeomTool *araGeom = AraGeomTool::Instance();
	
	double mostDelay=-100;
	double delays[16];
	for(int i=0; i<16; i++){
		delays[i] = (double) araGeom->getStationInfo(station)->getCableDelay(i);
		printf("Chan %d, delay %.2f \n", i, delays[i]);
		if(delays[i]>mostDelay) mostDelay=delays[i];
	}
	printf("Most delay is %.2f \n", mostDelay);

	int preTriggerBlocks=0;
	int readoutBlocks=0;
	if(station==2){
		if(config==1){
			preTriggerBlocks=12; //8+4
			readoutBlocks=20;
		}

	}

	double window_stop[16];
	double window_start[16];
	for(int i=0; i<16; i++){
		double T0 = 0;
		T0+=40.;
		if(station==2){
			if(config==1 || config==3 || config==4){
				T0+=mostDelay;
			}
		}
		window_start[i] = T0 - (preTriggerBlocks*20.);
		window_stop[i] = window_start[i] + (readoutBlocks*20.);
		printf("Chan %d, start %.2f, stop %.2f \n", i, window_start[i], window_stop[i]);
	}

	TH1D *h1_num_hit = new TH1D("raw","raw",17,-0.5, 16.5);
	TH1D *h1_num_hit_cable = new TH1D("fix cable delays","fix cable delays",17,-0.5, 16.5);
	TH1D *h1_num_hit_redotrig = new TH1D("redo trigger", "redo trigger", 17, -0.5, 16.5);

	int logE=0;

	for(int file=3; file<argc; file++){
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
			logE = int(TMath::Log10(pnu));
			
			// first compute hits straight from AraSim
			int numHits=0;
			for(int i=0; i<16; i++){
				if(RPR[i]>=8.){
					numHits++;
				}
			}
			h1_num_hit->Fill(numHits, weight);
			

			// now compute hits if we adjust position in buffers from cable delays
			int numHits_afterCorrectCableDelay=0;
			for(int i=0; i<16; i++){
				if(RPR[i]>=8.){
					// do hit time + cable delay, which is the location in the buffer
					double corrected_hit_time = times[i][2] + delays[i];
					// figure out if location in the buffer is inside readout window
					// if so, increment up the number of hits 
					if(corrected_hit_time<times[i][1] && corrected_hit_time>times[i][0]){
						numHits_afterCorrectCableDelay++;
					}
				}
			}
			h1_num_hit_cable->Fill(numHits_afterCorrectCableDelay, weight);


			// now compute hits if we had recalculated the trigger window from first principles
			int numHits_afterReformtrigger=0;
			for(int i=0; i<16; i++){
				if(RPR[i]>8.){
					// do hit time + cable delay, which is the location in the buffer
					double corrected_hit_time = times[i][2] + delays[i];
					// figure out if location in the buffer is inside first principles trigger window
					// if so, increment up the number of hits
					if(corrected_hit_time<window_stop[i] && corrected_hit_time>window_start[i]){
						numHits_afterReformtrigger++;
					}
				}
			}
			h1_num_hit_redotrig->Fill(numHits_afterReformtrigger, weight);
		}
		fpIn->Close();
		delete fpIn;
	}

	gStyle->SetOptStat(0);


	// normalize
	h1_num_hit->Sumw2();
	Double_t scale = 1/h1_num_hit->Integral("width");
	h1_num_hit->Scale(scale);

	h1_num_hit_cable->Sumw2();
	scale = 1/h1_num_hit_cable->Integral("width");
	h1_num_hit_cable->Scale(scale);

	h1_num_hit_redotrig ->Sumw2();
	scale = 1/h1_num_hit_redotrig->Integral("width");
	h1_num_hit_redotrig->Scale(scale);

	// get CDF
	TH1 *h1_cumulative_hit = h1_num_hit->GetCumulative();
	TH1 *h1_cumulative_hit_cable = h1_num_hit_cable->GetCumulative();
	TH1 *h1_cumulative_hit_redotrig = h1_num_hit_redotrig->GetCumulative();

	double upperBound;
	if(logE==19){
		upperBound=0.35;
	}
	else if(logE==18){
		upperBound=0.20;
	}
	else if(logE==17){
		upperBound=0.2;
	}
	stringstream ss;
	TH2D *dummy_hists[3];
	for(int i=0; i<3; i++){
		ss.str("");
		ss<<"thing "<<i;
		dummy_hists[i] = new TH2D(ss.str().c_str(), ss.str().c_str(), 17,-0.5, 16.5, 10, 0, upperBound);
	}

	TCanvas *c = new TCanvas("","",2*550,3*425);
	c->Divide(2,3);

	// first, "Raw" from AraSim
	c->cd(1);
		dummy_hists[0]->Draw("");
			dummy_hists[0]->SetTitle("unaltered simulation; Number of Hits; Normalized Weighted Events");
			h1_num_hit->Draw("samehist");
			h1_num_hit->SetLineWidth(2);
			h1_num_hit->SetTitle("unaltered simulation; Number of Hits; Normalized Weighted Events");
	c->cd(2);
		h1_cumulative_hit->Draw("hist");
			h1_cumulative_hit->SetLineWidth(2);
			h1_cumulative_hit->SetTitle("unaltered simulation; Number of Hits; Cumulative Fraction");

	// second, correcting only cable delays
	c->cd(3);
		dummy_hists[1]->Draw("");
			dummy_hists[1]->SetTitle("unaltered simulation; Number of Hits; Normalized Weighted Events");
			h1_num_hit_cable->Draw("samehist");
			h1_num_hit_cable->SetLineWidth(2);
			h1_num_hit_cable->SetTitle("accounting for cable delays; Number of Hits; Normalized Weighted Events");
	c->cd(4);
		h1_cumulative_hit_cable->Draw("hist");
			h1_cumulative_hit_cable->SetLineWidth(2);
			h1_cumulative_hit_cable->SetTitle("accounting for cable delays; Number of Hits; Cumulative Fraction");

	// third, reforming the trigger from "scratch"
	c->cd(5);
		dummy_hists[2]->Draw("");
			dummy_hists[2]->SetTitle("unaltered simulation; Number of Hits; Normalized Weighted Events");
			h1_num_hit_redotrig->Draw("samehist");
			h1_num_hit_redotrig->SetLineWidth(2);
			h1_num_hit_redotrig->SetTitle("redefine trigger window; Number of Hits; Normalized Weighted Events");
	c->cd(6);
		h1_cumulative_hit_redotrig->Draw("hist");
			h1_cumulative_hit_redotrig->SetLineWidth(2);
			h1_cumulative_hit_redotrig->SetTitle("redefine trigger window; Number of Hits; Cumulative Fraction");

	cout<<h1_cumulative_hit->GetBinCenter(5)<<endl;
	printf("Value at 4 hits: default %.2f, cable delay %.2f, redo trig window %.2f \n", 
		h1_cumulative_hit->GetBinContent(5), h1_cumulative_hit_cable->GetBinContent(5),
		h1_cumulative_hit_redotrig->GetBinContent(5));

	char save_title[150];
	// sprintf(save_title,"dist_hits_adjust_delays_and_window_E%d.png",logE);
	// c->SaveAs(save_title);

	TCanvas *c2 = new TCanvas("","",2*550,425);
	c2->Divide(2,1);
	c2->cd(1);
		dummy_hists[0]->Draw("");
			dummy_hists[0]->SetTitle("Number of hits w/ RPR>=8; Number of Hits; Normalized Weighted Events");
		h1_num_hit->Draw("histsame");
		h1_num_hit_cable->Draw("samehist");
			h1_num_hit_cable->SetLineColor(kRed);
		h1_num_hit_redotrig->Draw("histsame");
			h1_num_hit_redotrig->SetLineColor(kGreen);
	c2->cd(2);
		h1_cumulative_hit->Draw("hist");
			h1_cumulative_hit->SetTitle("  ; Number of Hits; Cumulative Fraction");
		h1_cumulative_hit_cable->Draw("samehist");
			h1_cumulative_hit_cable->SetLineColor(kRed);
		h1_cumulative_hit_redotrig->Draw("histsame");
			h1_cumulative_hit_redotrig->SetLineColor(kGreen);		


		TLegend *legend = new TLegend(0.1,0.7,0.4,0.9);
			legend->AddEntry(h1_cumulative_hit,"unaltered","l");
			legend->AddEntry(h1_cumulative_hit_cable,"cable delay","l");
			legend->AddEntry(h1_cumulative_hit_redotrig,"redo trig","l");
			// legend->SetBorderSize(0);  //no border for legend
			legend->SetFillColor(0);  //fill color is white
			legend->Draw("same");
	sprintf(save_title,"dist_hits_adjust_delays_and_window_overlay_E%d.png",logE);
	c2->SaveAs(save_title);

}//close the main program






















