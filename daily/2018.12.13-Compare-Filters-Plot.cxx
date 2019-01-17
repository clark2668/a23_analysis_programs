////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	v2_analysis_save_vals.cxx 
////	A23 diffuse, save values for cuts
////
////	Nov 2018
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
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"

//AraRoot includes
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_inputParameters.h"

using namespace std;

int main(int argc, char **argv)
{
	gStyle->SetOptStat(0);
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	stringstream ss;
	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	
	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <thresh bin> <wfrms cut val> <station> <year> <joined filename 1> <joined filename 2 > ... <joined filename x>"<<endl;
		return 0;
	}
	int selected_bin = atoi(argv[1]);
	double selected_cut = double(atof(argv[2]));
	int station = atoi(argv[3]);
	int year = atoi(argv[4]);

	//just to have the cut parameters up front and easy to find

	// int thresholdBin_pol[]={3,5}; //bin 3 = 2.3, bin 5 = 2.5 //what is the faceRMS inclusion threshold?
	// double wavefrontRMScut[]={-1.5, -1.5}; //event wavefrontRMS < this value

	int thresholdBin_pol[]={selected_bin, selected_bin}; //bin 3 = 2.3, bin 5 = 2.5 //what is the faceRMS inclusion threshold?
	double wavefrontRMScut[]={selected_cut, selected_cut}; //event wavefrontRMS < this value

	TH2D *wfrms_plots[6];
	wfrms_plots[0] = new TH2D("Vpol_12face","Vpol_12face",100,-5,5,40,0,40);
	wfrms_plots[1] = new TH2D("Hpol_12face","Hpol_12face",100,-5,5,40,0,40);
	wfrms_plots[2] = new TH2D("Vpol_3face","Vpol_3face",100,-5,5,40,0,40);
	wfrms_plots[3] = new TH2D("Hpol_3face","Hpol_3face",100,-5,5,40,0,40);

	TH2D *wfrms_plots_cal[6];
	wfrms_plots_cal[0] = new TH2D("Vpol_12face_cal","Vpol_12face_cal",100,-5,5,40,0,40);
	wfrms_plots_cal[1] = new TH2D("Hpol_12face_cal","Hpol_12face_cal",100,-5,5,40,0,40);
	wfrms_plots_cal[2] = new TH2D("Vpol_3face_cal","Vpol_3face_cal",100,-5,5,40,0,40);
	wfrms_plots_cal[3] = new TH2D("Hpol_3face_cal","Hpol_3face_cal",100,-5,5,40,0,40);

	int num_total=0;
	int num_thermal=0;
	int num_passing[] = {0,0};
	int num_passing_alt[] = {0,0};

	for(int file_num=5; file_num<argc; file_num++){

		string chRun = "run";
		string file = string(argv[file_num]);
		size_t foundRun = file.find(chRun);
		string strRunNum = file.substr(foundRun+4,4);
		int runNum = atoi(strRunNum.c_str());
		printf("Run Number %d \n", runNum);

		// cout << "Run " << file_num << " :: " << argv[file_num] << endl;
		
		//first, load in the data file; this shoud be a "joined" file
		//meaning it should contain "filter" trees and "reco" trees
		TFile *inputFile = TFile::Open(argv[file_num]);
		if(!inputFile){
			cout<<"Can't open joined file!"<<endl;
			return -1;
		}
		
		//next, we need to load the filter tree
		ss.str("");
		ss << "OutputTree";
		TTree *inputTree_filter = (TTree*) inputFile->Get(ss.str().c_str());
		if(!inputTree_filter){
			cout<<"Can't open filter tree"<<endl;
			return -1;
		}
		double thirdVPeakOverRMS[3]; //the third highest vpeak over RMS
		double rms_pol_thresh_face[2][15][12];
		double rms_pol_thresh_face_alternate[2][15][3];
		bool isCalPulser;
		bool isSoftTrigger;
		int waveformLength[16];
		inputTree_filter->SetBranchAddress("thirdVPeakOverRMS", &thirdVPeakOverRMS);
		inputTree_filter->SetBranchAddress("rms_pol_thresh_face", &rms_pol_thresh_face);
		inputTree_filter->SetBranchAddress("rms_pol_thresh_face_alternate", &rms_pol_thresh_face_alternate);
		inputTree_filter->SetBranchAddress("isCalpulser",&isCalPulser);
		inputTree_filter->SetBranchAddress("isSoftTrigger",&isSoftTrigger);
		inputTree_filter->SetBranchAddress("waveformLength",&waveformLength);

		int numEntries = inputTree_filter->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;

		//now to loop over events
		for(int event=0; event<numEntries; event++){
			num_total++;

			inputTree_filter->GetEvent(event);

			bool isShort=false;
			bool failWavefrontRMS[2];
			failWavefrontRMS[0]=false;
			failWavefrontRMS[1]=false;

			for(int i=0;i<16;i++){ if(waveformLength[i]<550) isShort=true; }

			//filter associated parameters
			double SNRs[2];
			SNRs[0] = thirdVPeakOverRMS[0];
			SNRs[1] = thirdVPeakOverRMS[1];
			if(SNRs[0]>29.) SNRs[0]=29.;
			if(SNRs[1]>29.) SNRs[1]=29.;

			vector <double>  rms_faces_V;
			rms_faces_V.resize(12);
			vector <double> rms_faces_H;
			rms_faces_H.resize(12);

			vector <double>  rms_faces_V_alt;
			rms_faces_V_alt.resize(3);
			vector <double> rms_faces_H_alt;
			rms_faces_H_alt.resize(3);

			//now, we must loop over the faces
			for(int i=0; i<12; i++){
				rms_faces_V[i] = rms_pol_thresh_face[0][thresholdBin_pol[0]][i];  //this is right RMS for this polarization, threshold requirement, and face
				rms_faces_H[i] = rms_pol_thresh_face[1][thresholdBin_pol[1]][i];
			}
			for(int i=0; i<3; i++){
				rms_faces_V_alt[i] = rms_pol_thresh_face_alternate[0][thresholdBin_pol[0]][i];  //this is right RMS for this polarization, threshold requirement, and face
				rms_faces_H_alt[i] = rms_pol_thresh_face_alternate[1][thresholdBin_pol[1]][i];
			}

			//now to sort them smallest to largest; lowest RMS is best
			sort(rms_faces_V.begin(), rms_faces_V.end());
			sort(rms_faces_H.begin(), rms_faces_H.end());
			sort(rms_faces_V_alt.begin(), rms_faces_V_alt.end());
			sort(rms_faces_H_alt.begin(), rms_faces_H_alt.end());

			double bestFaceRMS[2];
			bestFaceRMS[0]=rms_faces_V[0];
			bestFaceRMS[1]=rms_faces_H[0];

			double bestFaceRMS_alt[2];
			bestFaceRMS_alt[0]=rms_faces_V_alt[0];
			bestFaceRMS_alt[1]=rms_faces_H_alt[0];

			if(!isCalPulser && !isShort){
				num_thermal++;
			}

			for(int pol=0; pol<2; pol++){
				if(!isShort){
					wfrms_plots[pol]->Fill(TMath::Log10(bestFaceRMS[pol]), SNRs[pol]);
					wfrms_plots[pol+2]->Fill(TMath::Log10(bestFaceRMS_alt[pol]), SNRs[pol]);
					if(isCalPulser){
						wfrms_plots_cal[pol]->Fill(TMath::Log10(bestFaceRMS[pol]), SNRs[pol]);
						wfrms_plots_cal[pol+2]->Fill(TMath::Log10(bestFaceRMS_alt[pol]), SNRs[pol]);
					}

					// if(TMath::Log10(bestFaceRMS[pol]) < wavefrontRMScut[pol]){
					// 	num_passing[pol]++;
					// }
					// if(TMath::Log10(bestFaceRMS[pol+2]) < wavefrontRMScut[pol]){
					// 	num_passing_alt[pol+2]++;
					// }
				}
			}//loop over polarization
		
		}//loop over events
		
		inputFile->Close();
		delete inputFile;
	} //end loop over input files

	TH1D *projections[4];
	TH1D *projections_cal[4];

	for(int pol=0; pol<2; pol++){
		projections[pol] = wfrms_plots[pol]->ProjectionX();
		projections[pol+2] = wfrms_plots[pol+2]->ProjectionX();
		projections_cal[pol] = wfrms_plots_cal[pol]->ProjectionX();
		projections_cal[pol+2] = wfrms_plots_cal[pol+2]->ProjectionX();
	}

	for(int pol=0; pol<2; pol++){
		projections[pol] = wfrms_plots[pol]->ProjectionX();
		projections[pol+2] = wfrms_plots[pol+2]->ProjectionX();
	}

	TCanvas *c = new TCanvas("","",2.1*850,3.1*850);
	c->Divide(2,4);
	for(int pol=0; pol<2; pol++){
		c->cd(pol+1);
			wfrms_plots[pol]->Draw("colz");
			wfrms_plots[pol]->GetXaxis()->SetTitle("log10(Wavefront RMS)");
			wfrms_plots[pol]->GetYaxis()->SetTitle("3rd Highest VPeak/RMS");
			wfrms_plots[pol]->GetZaxis()->SetRangeUser(1.,7e6);
			gPad->SetLogz();
		c->cd(pol+3);		
			wfrms_plots[pol+2]->Draw("colz");
			wfrms_plots[pol+2]->GetXaxis()->SetTitle("log10(Wavefront RMS)");
			wfrms_plots[pol+2]->GetYaxis()->SetTitle("3rd Highest VPeak/RMS");
			wfrms_plots[pol+2]->GetZaxis()->SetRangeUser(1.,7e6);
			gPad->SetLogz();
		c->cd(pol+5);		
			projections[pol]->Draw();
			projections[pol]->SetLineWidth(3);
			if(pol==0){
				projections_cal[pol]->Draw("same");
				projections_cal[pol]->SetLineWidth(2);
				projections_cal[pol]->SetLineColor(kRed);
			}
			projections[pol]->GetYaxis()->SetTitle("Counts");
			projections[pol]->GetXaxis()->SetTitle("log10(Wavefront RMS)");
			projections[pol]->GetYaxis()->SetRangeUser(1.,1e7);
			gPad->SetLogy();
		c->cd(pol+7);		
			projections[pol+2]->Draw();
			projections[pol+2]->SetLineWidth(3);
			if(pol==0){
				projections_cal[pol+2]->Draw("same");
				projections_cal[pol+2]->SetLineWidth(2);
				projections_cal[pol+2]->SetLineColor(kRed);
			}
			projections[pol+2]->GetYaxis()->SetTitle("Counts");
			projections[pol+2]->GetXaxis()->SetTitle("log10(Wavefront RMS)");
			projections[pol+2]->GetYaxis()->SetRangeUser(1.,1e7);
			gPad->SetLogy();
	}
	char title[300];
	sprintf(title, "/users/PAS0654/osu0673/A23_analysis_new2/results/%d.%d.%d_A%d_%d_%dEvents_CompareFilters_Vpol%.1f_Hpol%.1f.png",year_now, month_now, day_now,station,year,num_total,0.1*double(thresholdBin_pol[0]) + 2.0,0.1*double(thresholdBin_pol[1])+2.0);
	c->SaveAs(title);
	delete c;
	delete wfrms_plots[0]; delete wfrms_plots[1]; delete wfrms_plots[2]; delete wfrms_plots[3];
	for(int i=0; i<4; i++) delete projections[i];
}
