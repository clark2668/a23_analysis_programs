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
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TChain.h"
#include "TTimeStamp.h"
#include "TMath.h"
#include "TF1.h"

// custom includes
#include "tools_Stats.h"
#include "tools_Cuts.h"

using namespace std;

int main(int argc, char **argv){

	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	stringstream ss;
	gStyle->SetOptStat(0);

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;


	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <1-station> <2-config>"<<endl;;
		return -1;
	}
	int station = atoi(argv[1]);
	int config = atoi(argv[2]);

	if(station!=2 && station!=3){
		printf("No good! You asked for station %d, but this code only works for stations 2 and 3 \n",station);
		return -1;
	}

	char fileHundredName[400];
	sprintf(fileHundredName,"/users/PAS0654/osu0673/A23_analysis_new2/results/unblind/background_fit/A%d_c%d_sample100.root",station,config);
	TFile* fileHundred = new TFile(fileHundredName,"READ");
	TH1D *h1Hundred[2];
	h1Hundred[0] = (TH1D*)fileHundred->Get("DiffDistroV");
	h1Hundred[1] = (TH1D*)fileHundred->Get("DiffDistroH");

	TH2D *h2Hundred_SNRvsCorr[2];
	h2Hundred_SNRvsCorr[0] = (TH2D*)fileHundred->Get("2DDistroV");
	h2Hundred_SNRvsCorr[1] = (TH2D*)fileHundred->Get("2DDistroH");
	h2Hundred_SNRvsCorr[0]->SetTitle("2D Distro V - 100%");
	h2Hundred_SNRvsCorr[1]->SetTitle("2D Distro H - 100%");

	char fileHundredRootName[400];
	sprintf(fileHundredRootName,"/users/PAS0654/osu0673/A23_analysis_new2/results/unblind/background_fit/2d_cut_values_A%d_c%d_sample100.root",station,config);
	TFile *outFile = TFile::Open(fileHundredRootName,"READ");
	TTree *outTree = (TTree*) outFile->Get("outTree");
	int hist_this_pol[2];
	double corr_val_out[2];
	double snr_val_out[2];
	int runNum_out;
	int eventNum_out;
	outTree->SetBranchAddress("passes_this_pol_V",&hist_this_pol[0]);
	outTree->SetBranchAddress("passes_this_pol_H",&hist_this_pol[1]);
	outTree->SetBranchAddress("corr_val_V",&corr_val_out[0]);
	outTree->SetBranchAddress("corr_val_H",&corr_val_out[1]);
	outTree->SetBranchAddress("snr_val_V",&snr_val_out[0]);
	outTree->SetBranchAddress("snr_val_H",&snr_val_out[1]);
	outTree->SetBranchAddress("runNum_out",&runNum_out);
	outTree->SetBranchAddress("eventNum_out",&eventNum_out);
	int numEntries = outTree->GetEntries();

	TH2D *redo[2];
	redo[0] = new TH2D("vredo","vredo",100,0,0.05,300,0,30);
	redo[1] = new TH2D("hredo","hredo",100,0,0.05,300,0,30);

	for(int i=0; i<numEntries; i++){
		outTree->GetEntry(i);
		if(hist_this_pol[0]){
			redo[0]->Fill(corr_val_out[0],snr_val_out[0]);
		}
		if(hist_this_pol[1]){
			redo[1]->Fill(corr_val_out[1],snr_val_out[1]);
			if(config==2 && snr_val_out[1]>8){
				printf("Run %d, Event %d \n", runNum_out, eventNum_out);
			}
		}
	}

	vector<TGraph*> cut_lines;
	for(int pol=0; pol<2; pol++){
		vector <double> x_vals_for_line;
		vector <double> y_vals_for_line;
		double rcut_slope;
		double rcut_intercept;
		getRCutValues(station, config, pol, rcut_slope, rcut_intercept);
		for(double x=0; x<0.020; x+=0.00001){
			double y_val = (rcut_slope * x ) + rcut_intercept;
			x_vals_for_line.push_back(x);
			y_vals_for_line.push_back(y_val);
		}
		cut_lines.push_back(new TGraph(x_vals_for_line.size(), &x_vals_for_line[0], &y_vals_for_line[0]));
	}

	TCanvas *c2 = new TCanvas("","",2*850,2*850);
	c2->Divide(2,2);
	for(int pol=0; pol<2; pol++){
		c2->cd(pol+1);
		h2Hundred_SNRvsCorr[pol]->Draw("colz");
		h2Hundred_SNRvsCorr[pol]->GetXaxis()->SetTitle("Correlation Value");
		h2Hundred_SNRvsCorr[pol]->GetYaxis()->SetTitle("SNR");
		cut_lines[pol]->Draw("same");
		cut_lines[pol]->SetLineColor(kRed);
		gPad->SetLogz();
	}
	for(int pol=0; pol<2; pol++){
		c2->cd(pol+3);
		redo[pol]->Draw("colz)");
		gPad->SetLogz();
	}
	char saveTitle[400];
	sprintf(saveTitle,"./A%d_c%d_Hundred2DDistro.png",station,config);
	c2->SaveAs(saveTitle);

	outFile->Close();
	fileHundred->Close();

}