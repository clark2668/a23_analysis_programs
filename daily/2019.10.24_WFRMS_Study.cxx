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
	
	if(argc<9){
		cout<< "Usage\n" << argv[0] << " <isSim?> <station> <config> <vbin> <hbin> <vcut> <hcut> <filter filename 1> <filter filename 2 > ... <filter filename x>"<<endl;
		return 0;
	}
	int isSim = atoi(argv[1]);
	int station = atoi(argv[2]);
	int config = atoi(argv[3]);

	//just to have the cut parameters up front and easy to find
	int thresholdBin_pol[] = {atoi(argv[4]), atoi(argv[5])}; //what is the faceRMS inclusion threshold?
	double wavefrontRMScut[]={atof(argv[6]), atof(argv[7])}; //event wavefrontRMS < this value

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

	double num_thermal[2] = {0, 0};
	double num_passing[] = {0.,0.};
	double num_passing_either=0.;

	// double num_passing[2][5][11] = {{{0.}}}; // [pol][snr_bin][wfrms_cut]
	// int numWFRMSlevels=11;
	// double wfrms_level_cuts[11]={-1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5};

	TH1D *h1_SNR_dist[16][2];
	for(int chan=0; chan<16; chan++){
		for(int peak=0; peak<2; peak++){
			stringstream ss;
			ss<<"Chan "<<chan<<" peak "<<peak;
			h1_SNR_dist[chan][peak] = new TH1D(ss.str().c_str(), ss.str().c_str(), 100,0,5);
		}
	}

	vector<int> BadRunList=BuildBadRunList(station);
	int num_total=0;

	for(int file_num=8; file_num<argc; file_num++){

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
		ss << "OutputTree";
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

		double SNR_for_WFRMS[16][2];
		double hitTimes_for_WFRMS[16][2];
		inputTree_filter->SetBranchAddress("SNR_for_WFRMS", &SNR_for_WFRMS);
		inputTree_filter->SetBranchAddress("hitTimes_for_WFRMS",&hitTimes_for_WFRMS);

		int numEntries = inputTree_filter->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;

		//now to loop over events
		// numEntries=15;
		for(int event=0; event<numEntries; event++){
			inputTree_filter->GetEvent(event);
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
			num_thermal[0]+=weight;
			num_thermal[1]+=weight;

			//filter associated parameters
			double SNRs[2];
			SNRs[0] = thirdVPeakOverRMS[0];
			SNRs[1] = thirdVPeakOverRMS[1];
			if(SNRs[0]>29.) SNRs[0]=29.;
			if(SNRs[1]>29.) SNRs[1]=29.;

			vector <double> rms_faces_V;
			vector <double> rms_faces_H;
			vector <double> rms_faces_V_scan;
			vector <double> rms_faces_H_scan;

			int num_faces_for_V_loop;
			int num_faces_for_H_loop;
			if(station==2){
				rms_faces_V.resize(numFaces);
				rms_faces_V_scan.resize(numFaces);
				num_faces_for_V_loop=numFaces;
				rms_faces_H.resize(numFaces_A2_drop);
				rms_faces_H_scan.resize(numFaces);
				num_faces_for_H_loop=numFaces_A2_drop;
			}
			else if(station==3){
				if(runNum>getA3BadRunBoundary()){ //it's 2014+, drop string four
					rms_faces_V.resize(numFaces_A3_drop);
					rms_faces_V_scan.resize(numFaces_A3_drop);
					num_faces_for_V_loop=numFaces_A3_drop;
					rms_faces_H.resize(numFaces_A3_drop);
					rms_faces_H_scan.resize(numFaces_A3_drop);
					num_faces_for_H_loop=numFaces_A3_drop;
					// printf(RED"Run is beyond bad boundary!\n");
				}
				else{ //it's 2013-, keep string four
					rms_faces_V.resize(numFaces);
					rms_faces_V_scan.resize(numFaces);
					num_faces_for_V_loop=numFaces;
					rms_faces_H.resize(numFaces);
					rms_faces_H_scan.resize(numFaces);
					num_faces_for_H_loop=numFaces;
				}
			}

			// // scan over threshold bins for 2D passing rate plot
			// for(int binThresh=0; binThresh<5; binThresh++){
			// 	for(int i=0; i<num_faces_for_V_loop; i++){
			// 		if (num_faces_for_V_loop==12){
			// 			rms_faces_V_scan[i] = rms_pol_thresh_face_V[binThresh][i];
			// 		}
	 	// 			else{
	 	// 				rms_faces_V_scan[i] = rms_pol_thresh_face_alternate_V[binThresh][i];
	 	// 			}
	 	// 			rms_faces_V_scan[i] = TMath::Log10(rms_faces_V_scan[i]);
			// 	}
			// 	for(int i=0; i<num_faces_for_H_loop; i++){
			// 		if (num_faces_for_H_loop==12){
			// 			rms_faces_H_scan[i] = rms_pol_thresh_face_H[binThresh][i];
			// 		}
			// 		else{
			// 			rms_faces_H_scan[i] = rms_pol_thresh_face_alternate_H[binThresh][i];
			// 		}
			// 		rms_faces_H_scan[i] = TMath::Log10(rms_faces_H_scan[i]);
			// 	}

			// 	sort(rms_faces_V_scan.begin(), rms_faces_V_scan.end());
			// 	sort(rms_faces_H_scan.begin(), rms_faces_H_scan.end());

			// 	for(int binWFRMS=0; binWFRMS<numWFRMSlevels; binWFRMS++){
			// 		if(rms_faces_V_scan[0]<wfrms_level_cuts[binWFRMS]){
			// 			num_passing[0][binThresh][binWFRMS]++;
			// 		}
			// 		if(rms_faces_V_scan[1]<wfrms_level_cuts[binWFRMS]){
			// 			num_passing[1][binThresh][binWFRMS]++;
			// 		}
			// 	}
			// }

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


				int faceToUse = i;
				if(station==3){
					if(config>2){
						if(i==0){
							faceToUse=0;
						}
						else if(i==1){
							faceToUse=3;
						}
						else if(i==2){
							faceToUse=7;
						}
					}
				}

				h1_face_RMSs_V[faceToUse]->Fill(rms_faces_V[i]);
				if(rms_faces_V[i]<bestWFRMS[0]){
					bestWFRMS[0]=rms_faces_V[i];
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
				
				int faceToUse=i;
				if(station==3){
					if(config>2){
						if(i==0){
							faceToUse=0;
						}
						else if(i==1){
							faceToUse=3;
						}
						else if(i==2){
							faceToUse=7;
						}
					}
				}
				if(station==2){
					if(i==0){
						faceToUse=0;
					}
					else if(i==1){
						faceToUse=3;
					}
					else if(i==2){
						faceToUse=4;
					}
					else if(i==3){
						faceToUse=7;
					}
					else if(i==4){
						faceToUse=9;
					}
					else if(i==5){
						faceToUse=11;
					}
				}

				h1_face_RMSs_H[faceToUse]->Fill(rms_faces_H[i]);
				if(rms_faces_H[i]<bestWFRMS[1]){
					bestWFRMS[1]=rms_faces_H[i];
					bestFace[1]=faceToUse;
				}
			}
			// printf("     Best H face is %d and best RMS is %.4f \n", bestWFRMS[1], bestFace[1]);
			if(bestWFRMS[1]<2){
				h1_face_RMSs_best_H[bestFace[1]]->Fill(bestWFRMS[1]);
				h1_best_face[1]->Fill(bestFace[1]);
			}

			for(int chan=0; chan<16; chan++){
				for(int peak=0; peak<2; peak++){
					h1_SNR_dist[chan][peak]->Fill(SNR_for_WFRMS[chan][peak]);
				}
			}

			//now to sort them smallest to largest; lowest RMS is best
			sort(rms_faces_V.begin(), rms_faces_V.end());
			sort(rms_faces_H.begin(), rms_faces_H.end());

			double bestFaceRMS[2];
			bestFaceRMS[0]=rms_faces_V[0];
			bestFaceRMS[1]=rms_faces_H[0];

			h1_best_face_RMS[0]->Fill(bestFaceRMS[0]);
			h1_best_face_RMS[1]->Fill(bestFaceRMS[1]);

			bool thisPasses[] = {false,false};
			for(int pol=0; pol<2; pol++){
				if(bestFaceRMS[pol]<wavefrontRMScut[pol]){
					num_passing[pol]+=weight;
					thisPasses[pol]=true;
				}
			}
			if(thisPasses[0] || thisPasses[1]){
				num_passing_either+=weight;
			}
		}//loop over events
		inputFile->Close();
		delete inputFile;
	} //end loop over input files

	// turn all of these into rates
	// for(int pol=0; pol<2; pol++){
	// 	for(int binSNR=0; binSNR<5; binSNR++){
	// 		for(int binWFRMS=0; binWFRMS<numWFRMSlevels; binWFRMS++){
	// 			num_passing[pol][binSNR][binWFRMS]/=num_thermal;
	// 		}
	// 	}
	// }

	for(int pol=0; pol<2; pol++){
		printf("Pol %d \n", pol);
		printf("-----------------------\n");
		printf("	Num thermal passing pol %d: %.3f/%.3f events, %.3f rate \n", pol, num_passing[pol],num_thermal[pol], 100.*num_passing[pol]/num_thermal[pol]);
	}
	printf("-----------------------\n");
	printf("Num passing either: %.3f/%.3f events, %.3f rate \n",num_passing_either, num_thermal[0], 100.*num_passing_either/num_thermal[0]);


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
	sprintf(title, "%s/wavefront_rms_teardown/%d.%d.%d_data_A%d_c%d_VBin%d_HBin%d_%dEvents_AllFaces_WFRMSDistro.png",plotPath,year_now, month_now, day_now,station,config,thresholdBin_pol[0], thresholdBin_pol[1], num_total);
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
	sprintf(title, "%s/wavefront_rms_teardown/%d.%d.%d_data_A%d_c%d_VBin%d_HBin%d_%dEvents_Best_WFRMSDistro.png",plotPath,year_now, month_now, day_now,station,config,thresholdBin_pol[0], thresholdBin_pol[1], num_total);
	c2->SaveAs(title);
	delete c2;

	vector<double> maxs;
	for(int i=0; i<16; i++){
		maxs.push_back(h1_SNR_dist[i][1]->GetMaximum());
	}
	std::sort(maxs.begin(), maxs.end()); //sort smallest to largest
	std::reverse(maxs.begin(), maxs.end()); //reverse order to get largest to smallest


	TCanvas *c3 = new TCanvas("","",4*1100,4*850);
	c3->Divide(4,4);
	for(int chan=0; chan<16; chan++){
		c3->cd(chan+1);
		h1_SNR_dist[chan][0]->Draw("");
		h1_SNR_dist[chan][0]->GetYaxis()->SetRangeUser(0,maxs[0]*1.1);
		h1_SNR_dist[chan][1]->Draw("same");
		h1_SNR_dist[chan][1]->SetLineColor(kRed);
		h1_SNR_dist[chan][1]->SetLineWidth(2);
		h1_SNR_dist[chan][0]->SetLineWidth(2);
		h1_SNR_dist[chan][0]->GetXaxis()->SetTitle("SNR");
		h1_SNR_dist[chan][0]->GetYaxis()->SetTitle("Number of Events");
	}
	sprintf(title, "%s/wavefront_rms_teardown/%d.%d.%d_data_A%d_c%d_%dEvents_SNR_Distributions.png",plotPath,year_now, month_now, day_now,station,config, num_total);
	c3->SaveAs(title);

	// TCanvas *c4 = new TCanvas("","",850,850);
	// h1_SNR_dist[2][0]->Draw("");
	// h1_SNR_dist[14][0]->Draw("same");
	// h1_SNR_dist[14][0]->SetLineColor(kGreen);
	// sprintf(title, "%s/wavefront_rms_teardown/%d.%d.%d_data_A%d_c%d_%dEvents_SNR_Distributions_Compare_1_and_13.png",plotPath,year_now, month_now, day_now,station,config, num_total);
	// c4->SaveAs(title);

	TCanvas *c5 = new TCanvas("","",4*1100,2*850);
	c5->Divide(4,2);
	for(int chan=0; chan<8; chan++){
		c5->cd(chan+1);
		h1_SNR_dist[chan][0]->Draw("");
		h1_SNR_dist[chan][0]->GetYaxis()->SetRangeUser(0,maxs[0]*1.1);
		h1_SNR_dist[chan+8][0]->Draw("same");
		h1_SNR_dist[chan+8][0]->SetLineColor(kRed);
		h1_SNR_dist[chan+8][0]->SetLineWidth(2);
		h1_SNR_dist[chan][0]->SetLineWidth(2);
		h1_SNR_dist[chan][0]->GetXaxis()->SetTitle("SNR");
		h1_SNR_dist[chan][0]->GetYaxis()->SetTitle("Number of Events");
		if(chan==0){
			h1_SNR_dist[chan][0]->SetTitle("String 1 Top Ants");
		}
		if(chan==1){
			h1_SNR_dist[chan][0]->SetTitle("String 2 Top Ants");
		}
		if(chan==2){
			h1_SNR_dist[chan][0]->SetTitle("String 3 Top Ants");
		}
		if(chan==3){
			h1_SNR_dist[chan][0]->SetTitle("String 4 Top Ants");
		}
		if(chan==4){
			h1_SNR_dist[chan][0]->SetTitle("String 1 Bottom Ants");
		}
		if(chan==5){
			h1_SNR_dist[chan][0]->SetTitle("String 2 Bottom Ants");
		}
		if(chan==6){
			h1_SNR_dist[chan][0]->SetTitle("String 3 Bottom Ants");
		}
		if(chan==7){
			h1_SNR_dist[chan][0]->SetTitle("String 4 Bottom Ants");
		}
		if(chan==0){
			TLegend *leg = new TLegend(0.15,0.7,0.4,0.9);
			leg->AddEntry(h1_SNR_dist[chan][0],"VPol","l");
			leg->AddEntry(h1_SNR_dist[chan][1],"HPol","l");
			leg->Draw();
		}
	}
	sprintf(title, "%s/wavefront_rms_teardown/%d.%d.%d_data_A%d_c%d_%dEvents_V_vs_H_SNR_Distributions.png",plotPath,year_now, month_now, day_now,station,config, num_total);
	c5->SaveAs(title);

}
