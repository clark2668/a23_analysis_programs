// this was just a first pass at doing the fits

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

	char fileName[400];
	sprintf(fileName,"/users/PAS0654/osu0673/A23_analysis_new2/results/post_unblinding/A%d_c%d_SurfaceDistro.root",station,config);
	TFile* fpIn = new TFile(fileName,"READ");
	TH1D *h1_theta[2];
	h1_theta[0] = (TH1D*)fpIn->Get("Vdist_theta");
	h1_theta[1] = (TH1D*)fpIn->Get("Hdist_theta");

	int max_bin[2];
	int last_filled_bin_above_2[2];
	double start_of_fit[2];
	double end_of_fit[2];

	start_of_fit[0] = start_of_fit[1] = 40;
	end_of_fit[0] = end_of_fit[1] = 45;

	// now we actually do the exponential fit
	char equation[150];
	sprintf(equation,"gaus");
	char equation_name[2][150];
	TF1 *fit[2];
	int status[2];
	double fitParams[2][3];
	double fitParamErrors[2][3];
	for(int pol=0; pol<2; pol++){
		sprintf(equation_name[pol],"ExpFit%d",pol);
		fit[pol] = new TF1(equation_name[pol],equation,start_of_fit[pol],end_of_fit[pol]);
		status[pol] = h1_theta[pol]->Fit(equation_name[pol],"LL,R");
		fitParams[pol][0] = fit[pol]->GetParameter(0);
		fitParams[pol][1] = fit[pol]->GetParameter(1);
		fitParams[pol][2] = fit[pol]->GetParameter(2);
		// fitParamErrors[pol][0] = fit[pol]->GetParError(0);
		// fitParamErrors[pol][1] = fit[pol]->GetParError(1);
		// printf("Pol %d Fit Parameters are %.2f and %.2f \n", pol, fitParams[pol][0], fitParams[pol][1]);
	}

	TF1 *fitCopy[2];

	char this_fit_title[2][400];
	sprintf(this_fit_title[0],"fCopyFitV");
	fitCopy[0] = new TF1(this_fit_title[0], "gaus", start_of_fit[0]-15., end_of_fit[0]);
	fitCopy[0]->SetParameters(fitParams[0][0], fitParams[0][1], fitParams[0][2]);

	sprintf(this_fit_title[1],"fCopyFitH");
	fitCopy[1] = new TF1(this_fit_title[1], "gaus", start_of_fit[1]-15., end_of_fit[1]);
	fitCopy[1]->SetParameters(fitParams[1][0], fitParams[1][1], fitParams[1][2]);

	double bot = 30.;
	double top = 45.;

	TCanvas *c = new TCanvas("","",2*850,850);
	c->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		c->cd(pol+1);
		h1_theta[pol]->Draw("");
		gPad->SetLogy();
		h1_theta[pol]->GetXaxis()->SetRangeUser(bot,top);
		h1_theta[pol]->GetXaxis()->SetTitle("#theta [deg]");
		h1_theta[pol]->GetYaxis()->SetTitle("Number of Events");
		// h1_theta[pol]->GetYaxis()->SetRangeUser(0.001,4e4);
		fitCopy[pol]->Draw("same");
		fitCopy[pol]->SetLineColor(kRed);
	}
	char saveTitle[500];
	sprintf(saveTitle,"./A%d_c%d_ThetaDistro.png",station,config);
	c->SaveAs(saveTitle);

	fpIn->Close();

}