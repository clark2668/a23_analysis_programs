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
#include "TLine.h"

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

	char fileTenName[400];
	sprintf(fileTenName,"/users/PCON0003/cond0068/ARA/AraRoot/analysis/unblind/background_fit/A%d_c%d_sample10.root",station,config);
	TFile* fileTen = new TFile(fileTenName,"READ");
	TH1D *h1Ten[2];
	h1Ten[0] = (TH1D*)fileTen->Get("DiffDistroV");
	h1Ten[1] = (TH1D*)fileTen->Get("DiffDistroH");

	TH1D *h1TenClone[2];
	h1TenClone[0] = (TH1D*) h1Ten[0]->Clone();
	h1TenClone[1] = (TH1D*) h1Ten[1]->Clone();

	char fileHundredName[400];
	sprintf(fileHundredName,"/users/PCON0003/cond0068/ARA/AraRoot/analysis/unblind/background_fit/A%d_c%d_sample100.root",station,config);
	TFile* fileHundred = new TFile(fileHundredName,"READ");
	TH1D *h1Hundred[2];
	h1Hundred[0] = (TH1D*)fileHundred->Get("DiffDistroV");
	h1Hundred[1] = (TH1D*)fileHundred->Get("DiffDistroH");

	TH1D *h1HundredClone[2];
	h1HundredClone[0] = (TH1D*) h1Hundred[0]->Clone();
	h1HundredClone[1] = (TH1D*) h1Hundred[1]->Clone();

	TH2D *h2Hundred_SNRvsCorr[2];
	h2Hundred_SNRvsCorr[0] = (TH2D*)fileHundred->Get("2DDistroV");
	h2Hundred_SNRvsCorr[1] = (TH2D*)fileHundred->Get("2DDistroH");
	h2Hundred_SNRvsCorr[0]->SetTitle("2D Distro V - 100%");
	h2Hundred_SNRvsCorr[1]->SetTitle("2D Distro H - 100%");

	// h1Ten[0]->Scale(10.);
	// h1Ten[1]->Scale(10.);

	/*
		Refit the 10%, scaled by 100
	*/
	int max_bin[2];
	int last_filled_bin_above_2[2];
	int fit_start_bin[2];
	double start_of_fit[2];
	int last_filled_bin[2];
	double end_of_fit[2];

	for(int pol=0; pol<2; pol++){
		max_bin[pol] = h1Ten[pol]->GetMaximumBin();
		// cout<<"Max bin is "<<max_bin[pol]<<endl;
		last_filled_bin_above_2[pol] = h1Ten[pol]->FindLastBinAbove(2.,1) + 2;
		// cout<<"Last filled bin above 2 is "<<last_filled_bin_above_2[pol]<<endl;
		fit_start_bin[pol] = int((last_filled_bin_above_2[pol] - max_bin[pol])/2) + max_bin[pol]; //start half-way between the peak bin and the last filled bin
		// cout<<"Fit start bin is "<<fit_start_bin[pol]<<endl;
		start_of_fit[pol] = h1Ten[pol]->GetBinCenter(fit_start_bin[pol]);
		// cout<<"Start of fit is "<<start_of_fit[pol]<<endl;
		last_filled_bin[pol] = h1Ten[pol]->FindLastBinAbove(0.,1);
		end_of_fit[pol] = h1Ten[pol]->GetBinCenter(last_filled_bin[pol]+2.); //go two bins more just to make sure fit is over
		// printf("Pol %d Start of fit is %.2f and end of fit is %.2f \n", pol, start_of_fit[pol], end_of_fit[pol]);

		// printf("Pol %d: Last filled bin is bin %d and value %.2f \n", pol, last_filled_bin[pol], h1Ten[pol]->GetBinCenter(last_filled_bin[pol]));
		// printf("Pol %d: Max bin is bin %d and value %.2f \n", pol, max_bin[pol], h1Ten[pol]->GetBinCenter(max_bin[pol]));
		// printf("Pol %d: Proposed start of fit is bin %d and value %.2f \n", pol, fit_start_bin[pol], h1Ten[pol]->GetBinCenter(fit_start_bin[pol]));
	}

	// now we actually do the exponential fit
	char equation[150];
	sprintf(equation,"exp([0]*x+[1])");
	char equation_name[2][150];
	TF1 *fit[2];
	int status[2];
	double fitParams[2][2];
	double fitParamErrors[2][2];
	for(int pol=0; pol<2; pol++){
		sprintf(equation_name[pol],"ExpFit%d",pol);
		fit[pol] = new TF1(equation_name[pol],equation,start_of_fit[pol],end_of_fit[pol]);
		status[pol] = h1Ten[pol]->Fit(equation_name[pol],"LL,R");
		fitParams[pol][0] = fit[pol]->GetParameter(0);
		fitParams[pol][1] = fit[pol]->GetParameter(1);
		fitParamErrors[pol][0] = fit[pol]->GetParError(0);
		fitParamErrors[pol][1] = fit[pol]->GetParError(1);
		printf("Pol %d Fit Parameters are %.2f and %.2f \n", pol, fitParams[pol][0], fitParams[pol][1]);
	}

	double fitParams_scale[2][2];
	for(int pol=0; pol<2; pol++){
		fitParams_scale[pol][0] = fitParams[pol][0];
		fitParams_scale[pol][1] = TMath::Log(10.*TMath::Exp(fitParams[pol][1]));
		// fitParams_scale[pol][1] = TMath::Log(1.*TMath::Exp(fitParams[pol][1]));
	}

	// printf("Original amplitude %.2f and scale %.2f \n", fitParams[0][1], fitParams_scale[0][1]);

	// and we record some of the information about the fit
	// like the name, the number of obsered events, etc.
	double binWidthIntercept[2];
	double leftBoundary[2];
	double rightBoundary[2];
	int numBinsThisFit[2];
	char this_fit_title[2][400];
	TH1D *hNumObserved[2];
	for(int pol=0; pol<2; pol++){
		binWidthIntercept[pol] = h1Ten[pol]->GetBinWidth(1);
		leftBoundary[pol] = start_of_fit[pol] - binWidthIntercept[pol]/2.;
		rightBoundary[pol] = end_of_fit[pol] + binWidthIntercept[pol]/2.;
		numBinsThisFit[pol] = (rightBoundary[pol] - leftBoundary[pol])/binWidthIntercept[pol] + 1;
		sprintf(this_fit_title[pol],"Fit_Pol%d",pol);
		hNumObserved[pol] = new TH1D(this_fit_title[pol],"",numBinsThisFit[pol],leftBoundary[pol],rightBoundary[pol]);
		for(int bin=0; bin<numBinsThisFit[pol]; bin++){

			// use the 100%!
			// double originalContent = h1Hundred[pol]->GetBinContent(bin+fit_start_bin[pol]);
			double originalContent = h1Ten[pol]->GetBinContent(bin+fit_start_bin[pol]);
			hNumObserved[pol]->SetBinContent(bin+1,originalContent);
		}
	}

	TF1 *fitCopy[2];

	sprintf(this_fit_title[0],"fCopyFitV");
	fitCopy[0] = new TF1(this_fit_title[0], "exp([0]*x+[1])", start_of_fit[0], end_of_fit[0]);
	fitCopy[0]->SetParameters(fitParams_scale[0][0], fitParams_scale[0][1]);

	sprintf(this_fit_title[1],"fCopyFitH");
	fitCopy[1] = new TF1(this_fit_title[1], "exp([0]*x+[1])", start_of_fit[1], end_of_fit[1]);
	fitCopy[1]->SetParameters(fitParams_scale[1][0], fitParams_scale[1][1]);

	vector<TGraph*> cut_lines;
	double selected_intercepts[2];
	for(int pol=0; pol<2; pol++){
		vector <double> x_vals_for_line;
		vector <double> y_vals_for_line;
		double rcut_slope;
		double rcut_intercept;
		getRCutValues(station, config, pol, rcut_slope, rcut_intercept);
		selected_intercepts[pol]=rcut_intercept;
		for(double x=0; x<0.020; x+=0.00001){
			double y_val = (rcut_slope * x ) + rcut_intercept;
			x_vals_for_line.push_back(x);
			y_vals_for_line.push_back(y_val);
		}
		cut_lines.push_back(new TGraph(x_vals_for_line.size(), &x_vals_for_line[0], &y_vals_for_line[0]));
	}

	TCanvas *c2 = new TCanvas("","",2*850,850);
	c2->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		c2->cd(pol+1);
		h2Hundred_SNRvsCorr[pol]->Draw("colz");
		h2Hundred_SNRvsCorr[pol]->GetXaxis()->SetTitle("Correlation Value");
		h2Hundred_SNRvsCorr[pol]->GetYaxis()->SetTitle("SNR");
		cut_lines[pol]->Draw("same");
		cut_lines[pol]->SetLineColor(kRed);
		gPad->SetLogz();
	}
	char saveTitle[400];
	sprintf(saveTitle,"./A%d_c%d_Hundred2DDistro.png",station,config);
	c2->SaveAs(saveTitle);

	TLine *lines[2];
	lines[0] = new TLine(selected_intercepts[0],0.5,selected_intercepts[0],1e4);
	lines[1] = new TLine(selected_intercepts[1],0.5,selected_intercepts[1],1e4);

	TCanvas *c = new TCanvas("","",2*850,850);
	c->Divide(2,1);
	c->cd(1);
			h1Hundred[0]->Draw("");
			h1Hundred[0]->GetXaxis()->SetRangeUser(11.,18.);
			fitCopy[0]->Draw("same");
			fitCopy[0]->SetLineColor(kBlue);
			fitCopy[0]->SetLineStyle(9);
			h1Hundred[0]->GetXaxis()->SetTitle("y-intercept value");
			h1Hundred[0]->GetYaxis()->SetTitle("differential number of events cut");
			h1Ten[0]->Draw("same");
			h1Ten[0]->SetLineColor(kRed);
			fit[0]->Draw("same");
			fit[0]->SetLineStyle(9);
			lines[0]->Draw("same");
			lines[0]->SetLineColor(kBlack);
			lines[0]->SetLineStyle(9);
			lines[0]->SetLineWidth(2);
			gPad->SetLogy();
		{
			TLegend *leg = new TLegend(0.48,0.6,0.9,0.9);
			leg->AddEntry(h1Ten[0],"10% Data","l");
			leg->AddEntry(fit[0],"10% Fit","l");
			leg->AddEntry(h1Hundred[0],"100% Data","l");
			leg->AddEntry(fitCopy[0],"Scaled 10% Fit","l");
			leg->AddEntry(lines[0],"Cut","l");
			leg->Draw("same");
		}
	c->cd(2);
		h1Hundred[1]->Draw("");
			h1Hundred[1]->GetXaxis()->SetRangeUser(8.,16.);
			fitCopy[1]->Draw("same");
			fitCopy[1]->SetLineColor(kBlue);
			h1Hundred[1]->GetXaxis()->SetTitle("y-intercept value");
			h1Hundred[1]->GetYaxis()->SetTitle("differential number of events cut");
		h1Ten[1]->Draw("same");
			h1Ten[1]->SetLineColor(kRed);
		lines[1]->Draw("same");
			lines[1]->SetLineColor(kBlack);
			lines[1]->SetLineStyle(9);
			lines[1]->SetLineWidth(2);
		gPad->SetLogy();
	sprintf(saveTitle,"./A%d_c%d_compare_TenHundred.png",station,config);
	c->SaveAs(saveTitle);

	fileTen->Close();
	fileHundred->Close();

}
