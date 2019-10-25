////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	2019.10.21_WFRMS_Study.cxx 
////	deep dive into how the wfrms parameter is getting put together, esp in A3
////
////	Oct 2010
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
#include "TLegend.h"

//AraRoot includes
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_inputParameters.h"
#include "tools_outputObjects.h"
#include "tools_Cuts.h"
#include "tools_CommandLine.h"

using namespace std;

int main(int argc, char **argv)
{
	gStyle->SetOptStat(0);
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;
	char *ToolsDirPath(getenv("TOOLS_DIR"));
	if(ToolsDirPath == NULL){
		std::cout << "Warning! $TOOLS_DIR is not set! This is needed for finding where your list of runs without cal pulses is "<< endl;
		return -1;
	}

	stringstream ss;
	
	if(argc<5){
		cout<< "Usage\n" << argv[0] << " <isSim?> <station> <config> <joined filename 1> <joined filename 2 > ... <joined filename x>"<<endl;
		return 0;
	}
	int isSim = atoi(argv[1]);
	int station = atoi(argv[2]);
	int config = atoi(argv[3]);

	//just to have the cut parameters up front and easy to find

	// int thresholdBin_pol[]={3,5}; //bin 3 = 2.3, bin 5 = 2.5 //what is the faceRMS inclusion threshold?
	// double wavefrontRMScut[]={-1.5, -1.5}; //event wavefrontRMS < this value

	int thresholdBin_pol[]={0, 0}; //bin 3 = 2.3, bin 5 = 2.5 //what is the faceRMS inclusion threshold?
	double wavefrontRMScut[]={2, 2}; //event wavefrontRMS < this value

	TH1D *h1_face_RMSs_V[12];
	TH1D *h1_face_RMSs_best_V[12];
	TH1D *h1_face_RMSs_H[12];
	TH1D *h1_face_RMSs_best_H[12];
	for(int i=0; i<12; i++){
		stringstream ss1;
		ss1<<"V Face "<<i;
		h1_face_RMSs_V[i] = new TH1D(ss1.str().c_str(),ss1.str().c_str(),200,-4,4);
		ss1.str("");
		ss1<<"V Face best"<<i;
		h1_face_RMSs_best_V[i] = new TH1D(ss1.str().c_str(),ss1.str().c_str(),200,-4,4);
		ss1.str("");
		ss1<<"H Face "<<i;
		h1_face_RMSs_H[i] = new TH1D(ss1.str().c_str(),ss1.str().c_str(),200,-4,4);
		ss1.str("");
		ss1<<"H Face best"<<i;
		h1_face_RMSs_best_H[i] = new TH1D(ss1.str().c_str(),ss1.str().c_str(),200,-4,4);
	}

	TH1D *h1_best_face_RMS[2];
	h1_best_face_RMS[0] = new TH1D("v_best_face_rms","v_best_face_rms",200,-4,4);
	h1_best_face_RMS[1] = new TH1D("h_best_face_rms","h_best_face_rms",200,-4,4);

	TH1D *h1_best_face[2];
	h1_best_face[0] = new TH1D("v_best_face","v_best_face_rms",12,0,12);
	h1_best_face[1] = new TH1D("h_best_face","h_best_face_rms",12,0,12);

	vector<int> BadRunList=BuildBadRunList(station);
	int num_total=0;

	for(int file_num=4; file_num<argc; file_num++){

		string chRun = "run";
		string file = string(argv[file_num]);
		size_t foundRun = file.find(chRun);
		string strRunNum = file.substr(foundRun+4,4);
		int runNum = atoi(strRunNum.c_str());
		printf("Run Number %d \n", runNum);

		bool doesRunHaveUntagedCalPulses = hasUntaggedCalpul(ToolsDirPath, station, config, runNum);
		if((isBadRun(station,runNum,BadRunList) || doesRunHaveUntagedCalPulses) && !isSim){
			printf("     Skipping run %d \n",runNum);
			continue;
		}

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
		ss << "OutputTree_filter";
		TTree *inputTree_filter = (TTree*) inputFile->Get(ss.str().c_str());
		if(!inputTree_filter){
			cout<<"Can't open filter tree"<<endl;
			return -1;
		}
		double thirdVPeakOverRMS[3]; //the third highest vpeak over RMS
		int numFaces_new_V;
		int numFaces_new_H;
		if(station==2){
			numFaces_new_V = numFaces;
			numFaces_new_H = numFaces_A2_drop;
		}
		else if(station==3){
			numFaces_new_V = numFaces_A3_drop;
			numFaces_new_H = numFaces_A3_drop;
		}
		double rms_pol_thresh_face_alternate_V[thresholdSteps][numFaces_new_V];
		double rms_pol_thresh_face_alternate_H[thresholdSteps][numFaces_new_H];
		bool isCalPulser;
		bool isSoftTrigger;
		int waveformLength[16];
		double weight;
		double SNR_theory;
		int Trig_Pass[16] = {0};
		int unixTime;

		inputTree_filter->SetBranchAddress("thirdVPeakOverRMS", &thirdVPeakOverRMS);
		inputTree_filter->SetBranchAddress("rms_pol_thresh_face_V", &rms_pol_thresh_face_V);
		inputTree_filter->SetBranchAddress("rms_pol_thresh_face_H", &rms_pol_thresh_face_H);
		inputTree_filter->SetBranchAddress("rms_pol_thresh_face_alternate_V", &rms_pol_thresh_face_alternate_V);
		inputTree_filter->SetBranchAddress("rms_pol_thresh_face_alternate_H", &rms_pol_thresh_face_alternate_H);
		inputTree_filter->SetBranchAddress("isCalpulser",&isCalPulser);
		inputTree_filter->SetBranchAddress("isSoftTrigger",&isSoftTrigger);
		inputTree_filter->SetBranchAddress("waveformLength",&waveformLength);
		inputTree_filter->SetBranchAddress("weight",&weight);
		inputTree_filter->SetBranchAddress("maxPeakVfromSim", &SNR_theory);
		inputTree_filter->SetBranchAddress("Trig_Pass", &Trig_Pass);
		inputTree_filter->SetBranchAddress("unixTime",&unixTime);

		int numEntries = inputTree_filter->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;


		//now to loop over events
		// numEntries=15;
		for(int event=0; event<numEntries; event++){
			inputTree_filter->GetEvent(event);

			if(isCalpulser) cout<<"Hi, I'm a tagged cal pulser event"<<endl;

			if(isBadLivetime(station,unixTime) && !isSim){
				continue;
			}

			bool isShort=false;
			bool failWavefrontRMS[2];
			failWavefrontRMS[0]=false;
			failWavefrontRMS[1]=false;

			for(int i=0;i<16;i++){ if(waveformLength[i]<64) isShort=true; }

			if(isSoftTrigger || isShort || isCalPulser) continue;

			num_total++;

			//filter associated parameters
			double SNRs[2];
			SNRs[0] = thirdVPeakOverRMS[0];
			SNRs[1] = thirdVPeakOverRMS[1];
			if(SNRs[0]>29.) SNRs[0]=29.;
			if(SNRs[1]>29.) SNRs[1]=29.;

			vector <double> rms_faces_V;
			vector <double> rms_faces_H;

			int num_faces_for_V_loop;
			int num_faces_for_H_loop;
			if(station==2){
				rms_faces_V.resize(numFaces);
				num_faces_for_V_loop=numFaces;
				rms_faces_H.resize(numFaces_A2_drop);
				num_faces_for_H_loop=numFaces_A2_drop;
			}
			else if(station==3){
				if(runNum>getA3BadRunBoundary()){ //it's 2014+, drop string four
					rms_faces_V.resize(numFaces_A3_drop);
					num_faces_for_V_loop=numFaces_A3_drop;
					rms_faces_H.resize(numFaces_A3_drop);
					num_faces_for_H_loop=numFaces_A3_drop;
					// printf(RED"Run is beyond bad boundary!\n");
				}
				else{ //it's 2013-, keep string four
					rms_faces_V.resize(numFaces);
					num_faces_for_V_loop=numFaces;
					rms_faces_H.resize(numFaces);
					num_faces_for_H_loop=numFaces;
				}
			}

			for(int i=0; i<num_faces_for_V_loop; i++){
				if (num_faces_for_V_loop==12){
					rms_faces_V[i] = rms_pol_thresh_face_V[thresholdBin_pol[0]][i];
				}
 				else{
 					rms_faces_V[i] = rms_pol_thresh_face_alternate_V[thresholdBin_pol[0]][i];
 				}
 				rms_faces_V[i] = TMath::Log10(rms_faces_V[i]);
			}
			for(int i=0; i<num_faces_for_H_loop; i++){
				if (num_faces_for_H_loop==12){
					rms_faces_H[i] = rms_pol_thresh_face_H[thresholdBin_pol[1]][i];
				}
				else{
					rms_faces_H[i] = rms_pol_thresh_face_alternate_H[thresholdBin_pol[1]][i];
				}
				rms_faces_H[i] = TMath::Log10(rms_faces_H[i]);
			}

			int bestFace[2]={40,40};
			double bestWFRMS[2]={3.e3,3.e3};
			// cout<<"bestFace values are "<<bestFace[0]<<" and "<<bestFace[1]<<endl;

			// printf("Event %d \n", event);
			for(int i=0; i<num_faces_for_V_loop; i++){
				// printf("     V face %d RMS is %.4f \n", i, rms_faces_V[i]);
				h1_face_RMSs_V[i]->Fill(rms_faces_V[i]);
				if(rms_faces_V[i]<bestWFRMS[0]){
					bestWFRMS[0]=rms_faces_V[i];
					int faceToUse = i;

					if(station==3){
						if(config>2){
							if(faceToUse==0){
								faceToUse=1;
							}
							else if(faceToUse==1){
								faceToUse=3;
							}
							else if(faceToUse==2){
								faceToUse=7;
							}
						}
					}
					bestFace[0]=faceToUse;
				}
			}
			// printf("     Best V face is %d and best RMS is %.4f \n", bestWFRMS[0], bestFace[0]);
			if(bestWFRMS[0]<2){
				h1_face_RMSs_best_V[bestFace[0]]->Fill(bestWFRMS[0]);
				h1_best_face[0]->Fill(bestFace[0]);
			}

			for(int i=0; i<num_faces_for_H_loop; i++){
				// printf("     H face %d RMS is %.4f \n", i, rms_faces_H[i]);
				h1_face_RMSs_H[i]->Fill(rms_faces_H[i]);
				if(rms_faces_H[i]<bestWFRMS[1]){
					bestWFRMS[1]=rms_faces_H[i];
					int faceToUse=i;
					
					if(station==3){
						if(config>2){
							if(faceToUse==0){
								faceToUse=1;
							}
							else if(faceToUse==1){
								faceToUse=3;
							}
							else if(faceToUse==2){
								faceToUse=7;
							}
						}
					}
					bestFace[1]=faceToUse;
				}
			}
			// printf("     Best H face is %d and best RMS is %.4f \n", bestWFRMS[1], bestFace[1]);
			if(bestWFRMS[1]<2){
				h1_face_RMSs_best_H[bestFace[1]]->Fill(bestWFRMS[1]);
				h1_best_face[1]->Fill(bestFace[1]);
			}

			//now to sort them smallest to largest; lowest RMS is best
			sort(rms_faces_V.begin(), rms_faces_V.end());
			sort(rms_faces_H.begin(), rms_faces_H.end());

			double bestFaceRMS[2];
			bestFaceRMS[0]=rms_faces_V[0];
			bestFaceRMS[1]=rms_faces_H[0];

			h1_best_face_RMS[0]->Fill(bestFaceRMS[0]);
			h1_best_face_RMS[1]->Fill(bestFaceRMS[1]);

		}//loop over events
		
		inputFile->Close();
		delete inputFile;
	} //end loop over input files


	TCanvas *c = new TCanvas("","",4*850,3*850);
	c->Divide(3,4);
	for(int i=0; i<12; i++){
		c->cd(i+1);
		h1_face_RMSs_V[i]->Draw("");
		h1_face_RMSs_V[i]->GetXaxis()->SetTitle("log10(Wavefront RMS)");
		h1_face_RMSs_V[i]->GetYaxis()->SetTitle("Number of Events");
		h1_face_RMSs_V[i]->SetLineColor(kBlue);
		h1_face_RMSs_V[i]->SetLineWidth(3);
		h1_face_RMSs_H[i]->Draw("same");
		h1_face_RMSs_H[i]->SetLineColor(kRed);
		h1_face_RMSs_H[i]->SetLineWidth(3);

		h1_face_RMSs_V[i]->GetZaxis()->SetRangeUser(0.1,5e6);
		
		// h1_face_RMSs_best_V[i]->Draw("same");
		// h1_face_RMSs_best_V[i]->SetLineWidth(3);
		// h1_face_RMSs_best_V[i]->SetLineColor(kBlue-7);
		// h1_face_RMSs_best_H[i]->Draw("same");
		// h1_face_RMSs_best_H[i]->SetLineWidth(3);
		// h1_face_RMSs_best_H[i]->SetLineColor(kRed-7);
		gPad->SetLogy();
		if(i+1==1){
			TLegend *leg = new TLegend(0.2,0.6,0.4,0.9);
			leg->AddEntry(h1_face_RMSs_V[i],"VPol","l");
			// leg->AddEntry(h1_face_RMSs_best_V[i],"VPol Best","l");
			leg->AddEntry(h1_face_RMSs_H[i],"HPol","l");
			// leg->AddEntry(h1_face_RMSs_best_H[i],"HPol Best","l");			
			leg->Draw();
		}
	}
	char title[300];
	sprintf(title, "%s/filter_cut/%d.%d.%d_data_A%d_c%d_%dEvents_AllFaces_WFRMSDistro.png",plotPath,year_now, month_now, day_now,station,config, num_total);
	c->SaveAs(title);
	delete c;

	TCanvas *c2 = new TCanvas("","",2*1100,850);
	c2->Divide(2,1);
	c2->cd(1);
		h1_best_face_RMS[0]->Draw("");
		h1_best_face_RMS[0]->GetXaxis()->SetTitle("best log10(Wavefront RMS)");
		h1_best_face_RMS[0]->GetYaxis()->SetTitle("Number of Events");
		h1_best_face_RMS[0]->SetLineColor(kBlue);
		h1_best_face_RMS[0]->SetLineWidth(3);
		h1_best_face_RMS[0]->SetTitle("");

		h1_best_face_RMS[1]->Draw("same");
		h1_best_face_RMS[1]->SetLineColor(kRed);
		h1_best_face_RMS[1]->SetLineWidth(3);
		gPad->SetLogy();

		{
			TLegend *leg = new TLegend(0.15,0.7,0.4,0.9);
			leg->AddEntry(h1_best_face_RMS[0],"VPol","l");
			leg->AddEntry(h1_best_face_RMS[1],"HPol","l");
			leg->Draw();
		}
	c2->cd(2);
		h1_best_face[0]->Scale(1/h1_best_face[0]->Integral());
		h1_best_face[1]->Scale(1/h1_best_face[1]->Integral());

		h1_best_face[0]->Draw("");
		h1_best_face[0]->GetXaxis()->SetTitle("face with best log10(Wavefront RMS)");
		h1_best_face[0]->GetYaxis()->SetTitle("Fraction of Events");
		h1_best_face[0]->SetLineColor(kBlue);
		h1_best_face[0]->SetLineWidth(3);
		h1_best_face[0]->GetYaxis()->SetRangeUser(0.,.25);
		h1_best_face[0]->GetXaxis()->SetNdivisions(12,0,0,false);
		h1_best_face[0]->SetTitle("");
		h1_best_face[0]->GetYaxis()->SetTitleOffset(1.2);

		h1_best_face[1]->Draw("same");
		h1_best_face[1]->SetLineColor(kRed);
		h1_best_face[1]->SetLineWidth(3);
	sprintf(title, "%s/filter_cut/%d.%d.%d_data_A%d_c%d_%dEvents_Best_WFRMSDistro.png",plotPath,year_now, month_now, day_now,station,config, num_total);
	c2->SaveAs(title);
	delete c2;
}
