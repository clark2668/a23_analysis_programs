// here, we can either redo the fits
// or, we can remake the histograms with binning we like ourselves

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

	bool doVersion1=false;
	if(doVersion1){
		char fileName[400];
		sprintf(fileName,"/users/PAS0654/osu0673/A23_analysis_new2/results/post_unblinding/A%d_c%d_SurfaceDistro.root",station,config);
		TFile* fpIn = new TFile(fileName,"READ");
		TH1D *h1_theta[2];
		h1_theta[0] = (TH1D*)fpIn->Get("Vdist_theta");
		h1_theta[1] = (TH1D*)fpIn->Get("Hdist_theta");

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
	}

	bool doVersion2=true;
	if(doVersion2){

		double vslope, hslope;
		double vintercept, hintercept;
		getRCutValues(station, config, 0, vslope, vintercept);
		getRCutValues(station, config, 1, hslope, hintercept);
		double selected_slopes[2] = {vslope, hslope};
		double selected_intercepts[2] = {vintercept, hintercept};
		// selected_intercepts[0]+=50.;
		// selected_intercepts[1]+=50.;
		// selected_intercepts[0]-=0.0001;
		// selected_intercepts[1]-=1; // bring these down

		char fileName[400];
		sprintf(fileName,"/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/other_studies/surface_distro/reco_and_2d_cut_values_A%d_c%d_10sample.root",station,config);
		TFile *fpIn = TFile::Open(fileName,"READ");
		TTree *outTree = (TTree*) fpIn->Get("outTree");
		int hist_this_pol[2];
		int theta_val_300[2];
		int phi_val_300[2];
		double corr_val_300[2];
		double snr_val_300[2];
		int runNum_out;
		int eventNum_out;
		int unixTime_out;
		bool isThisABadRun_out;
		bool isThisBadLivetime_out;
		outTree->SetBranchAddress("passes_this_pol_V",&hist_this_pol[0]);
		outTree->SetBranchAddress("passes_this_pol_H",&hist_this_pol[1]);
		outTree->SetBranchAddress("theta_val_300_V",&theta_val_300[0]);
		outTree->SetBranchAddress("theta_val_300_H",&theta_val_300[1]);
		outTree->SetBranchAddress("phi_val_300_V",&phi_val_300[0]);
		outTree->SetBranchAddress("phi_val_300_H",&phi_val_300[1]);
		outTree->SetBranchAddress("corr_val_300_V",&corr_val_300[0]);
		outTree->SetBranchAddress("corr_val_300_H",&corr_val_300[1]);
		outTree->SetBranchAddress("snr_val_300_V",&snr_val_300[0]);
		outTree->SetBranchAddress("snr_val_300_H",&snr_val_300[1]);
		outTree->SetBranchAddress("runNum_out",&runNum_out);
		outTree->SetBranchAddress("eventNum_out",&eventNum_out);
		outTree->SetBranchAddress("unixTime_out",&unixTime_out);
		outTree->SetBranchAddress("isThisABadRun_out",&isThisABadRun_out);
		outTree->SetBranchAddress("isThisBadLivetime_out",&isThisBadLivetime_out);
		int numEntries = outTree->GetEntries();

		TH1D *h1_costheta[2];
		TH1D *h1_costheta_passingRcut[2];
		h1_costheta[0] = new TH1D("Vdist_costheta","Vdist_costheta",180,-1,1);
		h1_costheta[1] = new TH1D("Hdist_costheta","Hdist_costheta",180,-1,1);
		h1_costheta_passingRcut[0] = new TH1D("Vdist_costheta_reduced","Vdist_costheta_reduced",180,-1,1);
		h1_costheta_passingRcut[1] = new TH1D("Hdist_costheta_reduced","Hdist_costheta_reduced",180,-1,1);

		for(int i=0; i<numEntries; i++){
			outTree->GetEntry(i);
			// if(runNum_out!=2884) continue;
			if(isThisABadRun_out || isThisBadLivetime_out) continue;
			for(int pol=0; pol<2; pol++){
				
				double thisSinTheta=TMath::Sin(double(theta_val_300[pol])*TMath::DegToRad());
				h1_costheta[pol]->Fill(thisSinTheta);
				// if(hist_this_pol[pol]){
					// if(snr_val_300[pol]>5.5 && hist_this_pol[pol]){
					// 	h1_costheta_passingRcut[pol]->Fill(thisSinTheta);
					// }
					// double this_y_val = (selected_slopes[pol] * corr_val_300[pol] ) + selected_intercepts[pol];
					// cout<<"Corr is "<<corr_val_300[pol]<<" and snr is "<<snr_val_300[pol]<<" and this_y_val is "<<this_y_val<<endl;
					// cout<<"Selected slopes is "<<selected_slopes[pol]<<" and selected intercepts is "<<selected_intercepts[pol]<<endl;
					// cout<<"SNR cut is "<<this_y_val<<" and the SNR is "<<snr_val_300[pol]<<endl;
					// if(snr_val_300[pol]>=this_y_val){
					// 	cout<<"Run "<<runNum_out<<" with an SNR of "<<snr_val_300[pol]<<"has bad run flag of "<<isThisABadRun_out<<endl;
					// 	h1_costheta_passingRcut[pol]->Fill(thisSinTheta);
					// }
				// }
			}
		}

		double start_of_fit[2];
		double end_of_fit[2];

		start_of_fit[0] = start_of_fit[1] = 0.66;
		end_of_fit[0] = end_of_fit[1] = 0.73;

		// now we actually do the exponential fit
		char equation[150];
		sprintf(equation,"gaus");
		char equation_name[2][150];
		TF1 *fit[2];
		int status[2];
		double fitParams[2][3];
		double fitParamErrors[2][3];
		for(int pol=0; pol<2; pol++){
			sprintf(equation_name[pol],"GausFit%d",pol);
			fit[pol] = new TF1(equation_name[pol],equation,start_of_fit[pol],end_of_fit[pol]);
			status[pol] = h1_costheta[pol]->Fit(equation_name[pol],"LL,R");
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
		fitCopy[0] = new TF1(this_fit_title[0], "gaus", start_of_fit[0]-0.35, end_of_fit[0]);
		fitCopy[0]->SetParameters(fitParams[0][0], fitParams[0][1], fitParams[0][2]);

		sprintf(this_fit_title[1],"fCopyFitH");
		fitCopy[1] = new TF1(this_fit_title[1], "gaus", start_of_fit[1]-0.35, end_of_fit[1]);
		fitCopy[1]->SetParameters(fitParams[1][0], fitParams[1][1], fitParams[1][2]);	

		double bot = 25.;
		double top = 55.;
		double bot_coszen = TMath::Sin(TMath::DegToRad()*double(bot));
		double top_coszen = TMath::Sin(TMath::DegToRad()*double(bot));

		double surface = 37.;
		double surface_sintheta = TMath::Sin(TMath::DegToRad()*double(surface));
		TLine *line = new TLine(surface_sintheta,1e4, surface_sintheta,1e5);

		char vtitle[400];
		char htitle[400];

		sprintf(vtitle,"const %.2f, mean %.2f, sigma %.2f",fitParams[0][0], fitParams[0][1], fitParams[0][2]);
		sprintf(htitle,"const %.2f, mean %.2f, sigma %.2f",fitParams[1][0], fitParams[1][1], fitParams[1][2]);

		h1_costheta[0]->SetTitle(vtitle);
		h1_costheta[1]->SetTitle(htitle);

		gStyle->SetOptStat(0);
		TCanvas *cDistro = new TCanvas("","",2*1100,850);
		cDistro->Divide(2,1);
		for(int pol=0; pol<2; pol++){
			// cDistro->cd(pol+1);
			// 	h1_theta[pol]->Draw("");
			// 	h1_theta[pol]->GetXaxis()->SetRangeUser(bot,top);
			// 	h1_theta[pol]->GetXaxis()->SetTitle("#theta [deg]");
			// 	h1_theta[pol]->GetYaxis()->SetTitle("Number of Events");
			// 	gPad->SetLogy();
			cDistro->cd(pol+1);
				h1_costheta[pol]->Draw("");
					h1_costheta[pol]->GetXaxis()->SetRangeUser(0.3,0.8);
					h1_costheta[pol]->GetXaxis()->SetTitle("sin(#theta)");
					h1_costheta[pol]->GetYaxis()->SetTitle("Number of Events");
					// h1_costheta[pol]->GetYaxis()->SetRangeUser(0.01,1e5);
				// h1_costheta_passingRcut[pol]->Draw("same");
				// 	h1_costheta_passingRcut[pol]->SetLineColor(kRed);
				fitCopy[pol]->Draw("same");
					fitCopy[pol]->SetLineColor(kRed);
				// line->Draw("same");
				// 	line->SetLineColor(kGreen);
				gPad->SetLogy();
		}
		char save_temp_title[400];
		sprintf(save_temp_title,"./A%d_c%d_SinThetaDistributions.png",station,config);
		cDistro->SaveAs(save_temp_title);
		delete cDistro;

		fpIn->Close();

	}
}