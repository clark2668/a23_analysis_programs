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

	stringstream ss;
	
	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <isSim?> <thresh bin> <wfrms cut val> <station> <config> <joined filename 1> <joined filename 2 > ... <joined filename x>"<<endl;
		return 0;
	}
	int isSim = atoi(argv[1]);
	int selected_bin = atoi(argv[2]);
	double selected_cut = double(atof(argv[3]));
	int station = atoi(argv[4]);
	int config = atoi(argv[5]);

	//just to have the cut parameters up front and easy to find

	// int thresholdBin_pol[]={3,5}; //bin 3 = 2.3, bin 5 = 2.5 //what is the faceRMS inclusion threshold?
	// double wavefrontRMScut[]={-1.5, -1.5}; //event wavefrontRMS < this value

	int thresholdBin_pol[]={selected_bin, selected_bin}; //bin 3 = 2.3, bin 5 = 2.5 //what is the faceRMS inclusion threshold?
	double wavefrontRMScut[]={selected_cut, selected_cut}; //event wavefrontRMS < this value

	TH2D *wfrms_plots[2];
	wfrms_plots[0] = new TH2D("Vpol","Vpol_org",200,-5,5,30,0,30);
	wfrms_plots[1] = new TH2D("Hpol","Hpol_org",200,-5,5,30,0,30);

	TH2D *wfrms_plots_cal[2];
	wfrms_plots_cal[0] = new TH2D("Vpol_cal","Vpol_cal",200,-5,5,30,0,30);
	wfrms_plots_cal[1] = new TH2D("Hpol_cal","Hpol_cal",200,-5,5,30,0,30);

	TH2D *wfrms_plots_rf[2];
	wfrms_plots_rf[0] = new TH2D("Vpol_rf","Vpol_rf",200,-5,5,30,0,30);
	wfrms_plots_rf[1] = new TH2D("Hpol_rf","Hpol_rf",200,-5,5,30,0,30);

	TH1D *h1_best_face[2];
	h1_best_face[0] = new TH1D("vpol_best_face","vpol_best_face",12,0,12);
	h1_best_face[1] = new TH1D("hpol_best_face","hpol_best_face",12,0,12);

	TH1D *h1_best_face_thresh05[2];
	h1_best_face_thresh05[0] = new TH1D("vpol_best_face","vpol_best_face",12,0,12);
	h1_best_face_thresh05[1] = new TH1D("hpol_best_face","hpol_best_face",12,0,12);

	TH1D *h1_best_face_thresh1[2];
	h1_best_face_thresh1[0] = new TH1D("vpol_best_face","vpol_best_face",12,0,12);
	h1_best_face_thresh1[1] = new TH1D("hpol_best_face","hpol_best_face",12,0,12);

	TH1D *h1_best_face_thresh15[2];
	h1_best_face_thresh15[0] = new TH1D("vpol_best_face","vpol_best_face",12,0,12);
	h1_best_face_thresh15[1] = new TH1D("hpol_best_face","hpol_best_face",12,0,12);

	TH2D *dummy[2];
	dummy[0] = new TH2D("","",100,-5,5,2,0,100.);
	dummy[1] = new TH2D("","",100,-5,5,2,0,100.);

	TH1D *all_events[2];
	TH1D *passed_events[2];
	TH1D *eff[2];
	for(int i=0; i<2; i++){
		all_events[i] = new TH1D("","",30,0,30);
		passed_events[i] = new TH1D("","",30,0,30);
		eff[i] = new TH1D("","",30,0,30);
	}

	double num_total=0.;
	double num_thermal[] = {0., 0.};
	double num_passing[] = {0.,0.};
	double num_passing_either=0.;

	vector<int> BadRunList=BuildBadRunList(station);

	for(int file_num=6; file_num<argc; file_num++){

		string chRun = "run";
		string file = string(argv[file_num]);
		size_t foundRun = file.find(chRun);
		string strRunNum = file.substr(foundRun+4,4);
		int runNum = atoi(strRunNum.c_str());
		printf("Run Number %d \n", runNum);

		if(isBadRun(station,runNum,BadRunList) && !isSim){
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

			if(isBadLivetime(station,unixTime) && !isSim){
				continue;
			}

			num_total+=weight;

			bool isShort=false;
			bool failWavefrontRMS[2];
			failWavefrontRMS[0]=false;
			failWavefrontRMS[1]=false;

			for(int i=0;i<16;i++){ if(waveformLength[i]<64) isShort=true; }

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
			}
			for(int i=0; i<num_faces_for_H_loop; i++){
				if (num_faces_for_H_loop==12){
					rms_faces_H[i] = rms_pol_thresh_face_H[thresholdBin_pol[1]][i];
				}
				else{
					rms_faces_H[i] = rms_pol_thresh_face_alternate_H[thresholdBin_pol[1]][i];
				}
			}

			// printf("Num faces for V loop is %d \n", num_faces_for_V_loop);
			// printf("Num faces for H loop is %d \n", num_faces_for_H_loop);

			// //now we loop over the faces
			// for(int i=0; i<num_faces_for_V_loop; i++){
			// 	rms_faces_V[i] = rms_pol_thresh_face_alternate_V[thresholdBin_pol[0]][i];
			// }
			// for(int i=0; i<num_faces_for_H_loop; i++){
			// 	rms_faces_H[i] = rms_pol_thresh_face_alternate_H[thresholdBin_pol[1]][i];
			// }

			int bestFace[2]={40,40};
			double bestWFRMS[2]={3.e3,3.e3};
			// cout<<"bestFace values are "<<bestFace[0]<<" and "<<bestFace[1]<<endl;

			// printf("Event %d \n", event);
			for(int i=0; i<num_faces_for_V_loop; i++){
				// printf("     V face %d RMS is %.4f \n", i, rms_faces_V[i]);
				if(rms_faces_V[i]<bestWFRMS[0]){
					bestWFRMS[0]=rms_faces_V[i];
					int faceToUse = i;
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
					bestFace[0]=faceToUse;
				}
			}
			// printf("     Best V face is %d and best RMS is %.4f \n", bestWFRMS[0], bestFace[0]);

			for(int i=0; i<num_faces_for_H_loop; i++){
				// printf("     H face %d RMS is %.4f \n", i, rms_faces_H[i]);
				if(rms_faces_H[i]<bestWFRMS[1]){
					bestWFRMS[1]=rms_faces_H[i];
					int faceToUse=i;
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
					bestFace[1]=faceToUse;
				}
			}
			// printf("     Best H face is %d and best RMS is %.4f \n", bestWFRMS[1], bestFace[1]);

			//now to sort them smallest to largest; lowest RMS is best
			sort(rms_faces_V.begin(), rms_faces_V.end());
			sort(rms_faces_H.begin(), rms_faces_H.end());

			double bestFaceRMS[2];
			bestFaceRMS[0]=rms_faces_V[0];
			bestFaceRMS[1]=rms_faces_H[0];


			// bool trig_pol[2];
			// trig_pol[0]=false;
			// trig_pol[1]=false;

			// int num_trig_V=0;
			// int num_trig_H=0;
			// for(int string=0; string<4; string++){
			// 	for(int ant=0; ant<4; ant++){
			// 		int chan = string + 4*ant;
			// 		// printf("Trig pass for string %d ant %d is %d \n", string, ant, Trig_Pass[chan]);
			// 		if(Trig_Pass[chan]>0){
			// 			if(ant==0 || ant==2)
			// 				num_trig_V++;
			// 			else if(ant==1 || ant==3)
			// 				num_trig_H++;
			// 		}
			// 	}
			// }
			// if(num_trig_V>2){
			// 	trig_pol[0]=true;
			// }
			// if(num_trig_H>2){
			// 	trig_pol[1]=true;
			// }

			if(!isCalPulser && !isShort && !isSoftTrigger){
				num_thermal[0]+=weight;
				num_thermal[1]+=weight;
				if(isSim){
					// all_events[0]->Fill(SNR_theory/0.035,weight);
					// all_events[1]->Fill(SNR_theory/0.035,weight);
					all_events[0]->Fill(SNRs[0],weight);
					all_events[1]->Fill(SNRs[1],weight);
				}

				// if(!isSim){
				// 	num_thermal[0]+=weight;
				// 	num_thermal[1]+=weight;
				// }
				// else if(isSim){
				// 	if(trig_pol[0]){
				// 		num_thermal[0]+=weight;
				// 		all_events[0]->Fill(SNR_theory/0.035,weight);
				// 	}
				// 	else if (trig_pol[1]){
				// 		num_thermal[1]+=weight;
				// 		all_events[1]->Fill(SNR_theory/0.035,weight);
				// 	}
				// }
			}

			bool thisPasses[2] = {false};

			for(int pol=0; pol<2; pol++){
				if(!isShort && !isSoftTrigger){
					wfrms_plots[pol]->Fill(TMath::Log10(bestFaceRMS[pol]), SNRs[pol], weight);
					if(isCalPulser){
						wfrms_plots_cal[pol]->Fill(TMath::Log10(bestFaceRMS[pol]), SNRs[pol], weight);
					}
					else if(!isCalPulser){
						wfrms_plots_rf[pol]->Fill(TMath::Log10(bestFaceRMS[pol]), SNRs[pol], weight);
						if(TMath::Log10(bestFaceRMS[pol]) < wavefrontRMScut[pol]){
							num_passing[pol]+=(weight);
							thisPasses[pol]=true;
							if(isSim){
								passed_events[pol]->Fill(SNRs[pol],weight);
							}
						}
						if(TMath::Log10(bestFaceRMS[pol]) < 3){
							// store info about the best face
							h1_best_face[pol]->Fill(bestFace[pol]);
							// cout<<"Event "<<event<<" best face is "<<bestFace[pol]<<endl;
						}

						if(TMath::Log10(bestFaceRMS[pol]) < -0.5){
							// store info about the best face
							h1_best_face_thresh05[pol]->Fill(bestFace[pol]);
						}

						if(TMath::Log10(bestFaceRMS[pol]) < -1){
							// store info about the best face
							h1_best_face_thresh1[pol]->Fill(bestFace[pol]);
						}


						if(TMath::Log10(bestFaceRMS[pol]) < -1.5){
							// store info about the best face
							h1_best_face_thresh15[pol]->Fill(bestFace[pol]);
						}
						// if(TMath::Log10(bestFaceRMS_alt[pol]) < wavefrontRMScut[pol] && !trig_pol[pol]){
						// 	cout<<"Passes WFRMS but didn't trigger in pol "<<pol<<" with weight "<<weight<<endl;
						// }

					}
				}
			}//loop over polarization
			if(thisPasses[0] || thisPasses[1]){
				num_passing_either+=weight;
			}
		}//loop over events
		
		inputFile->Close();
		delete inputFile;
	} //end loop over input files

	// printf("Total Events: %.2f \n", num_total);
	// printf("-----------------------\n");
	// printf("-----------------------\n");
	// for(int pol=0; pol<2; pol++){
	// 	printf("Pol %d \n", pol);
	// 	printf("-----------------------\n");
	// 	printf("	Org : Num passing pol %d: %.3f events, %.3f rate \n", pol, num_passing[pol], 100.*num_passing[pol]/num_total);
	// 	printf("	Alt  face: Num passing pol %d: %.3f events, %.3f rate \n", pol, num_passing_alt[pol], 100.*num_passing_alt[pol]/num_total);
	// }
	// cout<<""<<endl;
	// cout<<""<<endl;
	// cout<<""<<endl;

	// printf("Total Thermal Events: %.2f \n", num_thermal);
	// printf("-----------------------\n");
	// printf("-----------------------\n");
	for(int pol=0; pol<2; pol++){
		printf("Pol %d \n", pol);
		printf("-----------------------\n");
		printf("	Num thermal passing pol %d: %.3f/%.3f events, %.3f rate \n", pol, num_passing[pol],num_thermal[pol], 100.*num_passing[pol]/num_thermal[pol]);
	}
	printf("-----------------------\n");
	printf("Num passing either: %.3f/%.3f events, %.3f rate \n",num_passing_either, num_thermal[0], 100.*num_passing_either/num_thermal[0]);

	TH1D *projections[2];
	TH1D *projections_cal[2];
	TH1D *projections_rf[2];

	for(int pol=0; pol<2; pol++){
		projections[pol] = wfrms_plots[pol]->ProjectionX();
		projections_cal[pol] = wfrms_plots_cal[pol]->ProjectionX();
		projections_rf[pol] = wfrms_plots_rf[pol]->ProjectionX();
	}


	/*
		Now, make plot of passing rate as function of WFRMS
	*/
	char title_txt[200];
	sprintf(title_txt,"%s/filter_cut/filter_rates_A%d_c%d_V%d_H%d.txt",plotPath,station,config,selected_bin,selected_bin);
	// FILE *fout = fopen(title_txt, "a");

	TGraph *passing_rate[2];	
	for(int pol=0; pol<2; pol++){
		if(pol==0){
			// fprintf(fout, "WFRMS, V_rate\n");
		}
		if(pol==1){
			// fprintf(fout, "--------------------------------------\n");
			// fprintf(fout, "WFRMS, H_rate\n");
		}
		// now compute the passing rate
		passing_rate[pol] = new TGraph();
		int numEventsTotal = projections_rf[pol]->Integral();
		int numPointsInGraph=0;
		double binWidth = projections_rf[pol]->GetBinWidth(2);
		for(int binX=0; binX<projections_rf[pol]->GetNbinsX(); binX++){
			double binEdge = projections_rf[pol]->GetBinCenter(binX) + binWidth/2.;
			double thisPassRate = 100.*projections_rf[pol]->Integral(0,binX)/numEventsTotal;
			passing_rate[pol]->SetPoint(numPointsInGraph,binEdge,thisPassRate);
			numPointsInGraph++;
		}
		for(int binX=0; binX<passing_rate[pol]->GetN(); binX++){
			double thisWFRMS = passing_rate[pol]->GetX()[binX];
			double thisRate = passing_rate[pol]->GetY()[binX];
			// fprintf(fout, "%1.2f, %0.5f\n", thisWFRMS, thisRate);
		}
	}
	// fclose(fout);//close sigmavsfreq.txt file


	// for(int pol=0; pol<2; pol++){
	// 	projections[pol] = wfrms_plots[pol]->ProjectionX();
	// 	projections[pol+2] = wfrms_plots[pol+2]->ProjectionX();
	// }

	// first, print out everything

	TCanvas *c = new TCanvas("","",2.1*850,2.1*850);
	c->Divide(2,2);
	for(int pol=0; pol<2; pol++){
		c->cd(pol+1);
			wfrms_plots[pol]->Draw("colz");
			wfrms_plots[pol]->GetXaxis()->SetTitle("log10(Wavefront RMS)");
			wfrms_plots[pol]->GetYaxis()->SetTitle("3rd Highest VPeak/RMS");
			wfrms_plots[pol]->GetZaxis()->SetRangeUser(0.1,7e6);
			gPad->SetLogz();
		c->cd(pol+3);		
			projections[pol]->Draw();
			projections[pol]->SetLineWidth(3);
			projections[pol]->GetYaxis()->SetTitle("Counts");
			projections[pol]->GetXaxis()->SetTitle("log10(Wavefront RMS)");
			projections[pol]->GetYaxis()->SetRangeUser(0.1,1e7);
			gPad->SetLogy();
	}
	char title[300];
	if(isSim)
		sprintf(title, "%s/filter_cut/%d.%d.%d_sim_A%d_c%d_%dEvents_Filter_AllEvents_Vpol%.1f_Hpol%.1f.png",plotPath,year_now, month_now, day_now,station,config,int(num_total),0.1*double(thresholdBin_pol[0]) + 2.0,0.1*double(thresholdBin_pol[1])+2.0);
	else
		sprintf(title, "%s/filter_cut/%d.%d.%d_data_A%d_c%d_%dEvents_Filter_AllEvents_Vpol%.1f_Hpol%.1f.png",plotPath,year_now, month_now, day_now,station,config,int(num_total),0.1*double(thresholdBin_pol[0]) + 2.0,0.1*double(thresholdBin_pol[1])+2.0);
	// c->SaveAs(title);
	delete c;

	TCanvas *c2 = new TCanvas("","",2.1*850,3.1*850);
	c2->Divide(2,3);
	for(int pol=0; pol<2; pol++){
		c2->cd(pol+1);
			wfrms_plots_rf[pol]->Draw("colz");
			wfrms_plots_rf[pol]->GetXaxis()->SetTitle("log10(Wavefront RMS)");
			wfrms_plots_rf[pol]->GetYaxis()->SetTitle("3rd Highest VPeak/RMS");
			wfrms_plots_rf[pol]->GetZaxis()->SetRangeUser(0.1,7e6);
			gPad->SetLogz();
		c2->cd(pol+3);		
			projections_rf[pol]->Draw();
			projections_rf[pol]->SetLineWidth(3);
			projections_rf[pol]->GetYaxis()->SetTitle("Counts");
			projections_rf[pol]->GetXaxis()->SetTitle("log10(Wavefront RMS)");
			projections_rf[pol]->GetYaxis()->SetRangeUser(0.1,1e7);
			// if(pol==0){
			// 	projections_cal[pol]->Draw("same");
			// 	projections_cal[pol]->SetLineColor(kRed);
			// 	projections_cal[pol]->SetLineWidth(3);
			// }
			gPad->SetLogy();
		c2->cd(pol+5);
			dummy[pol]->Draw();
			passing_rate[pol]->GetYaxis()->SetTitle("Passing Rate");
			passing_rate[pol]->GetXaxis()->SetTitle("log10(Wavefront RMS)");
			passing_rate[pol]->Draw("LP");
			// dummy[pol]->GetYaxis()->SetRangeUser(1.e-2,1.);
			// gPad->SetLogy();
	}
	if(isSim)
		sprintf(title, "%s/filter_cut/%d.%d.%d_sim_A%d_c%d_%dEvents_Filter_RFEvents_Vpol%.1f_Hpol%.1f.png",plotPath,year_now, month_now, day_now,station,config,int(num_total),0.1*double(thresholdBin_pol[0]) + 2.0,0.1*double(thresholdBin_pol[1])+2.0);
	else
		sprintf(title, "%s/filter_cut/%d.%d.%d_data_A%d_c%d_%dEvents_Filter_RFEvents_Vpol%.1f_Hpol%.1f.png",plotPath,year_now, month_now, day_now,station,config,int(num_total),0.1*double(thresholdBin_pol[0]) + 2.0,0.1*double(thresholdBin_pol[1])+2.0);
	// c2->SaveAs(title);
	delete c2;

	for(int i=0; i<2; i++){
		h1_best_face[i]->Scale(1/h1_best_face[i]->Integral());
		h1_best_face_thresh05[i]->Scale(1/h1_best_face_thresh05[i]->Integral());
		h1_best_face_thresh1[i]->Scale(1/h1_best_face_thresh1[i]->Integral());
		h1_best_face_thresh15[i]->Scale(1/h1_best_face_thresh15[i]->Integral());
	}

	TCanvas *c3 = new TCanvas("","",2*850,850);
	c3->Divide(2,1);
	for(int i=0; i<2; i++){
		c3->cd(i+1);
		h1_best_face[i]->Draw("");
		h1_best_face[i]->GetXaxis()->SetTitle("Face with best log10(WFRMS)");
		h1_best_face[i]->GetYaxis()->SetTitle("Fraction of Events");
		h1_best_face[i]->GetYaxis()->SetTitleOffset(1.6);
		// h1_best_face[i]->GetYaxis()->SetRangeUser(1e-1,2e3);
		h1_best_face[i]->GetYaxis()->SetRangeUser(0,0.4);
		
		if(config>2){
			h1_best_face[i]->GetYaxis()->SetRangeUser(0,1);
		}

		h1_best_face[i]->GetXaxis()->SetNdivisions(12,0,0,false);
		h1_best_face[i]->SetLineColor(kBlack);
		h1_best_face[i]->SetLineWidth(2);

		h1_best_face_thresh05[i]->Draw("same");
		h1_best_face_thresh05[i]->SetLineColor(kRed);
		h1_best_face_thresh05[i]->SetLineWidth(2);

	
		h1_best_face_thresh1[i]->Draw("same");
		h1_best_face_thresh1[i]->SetLineColor(kBlue);
		h1_best_face_thresh1[i]->SetLineWidth(2);

		h1_best_face_thresh15[i]->Draw("same");
		h1_best_face_thresh15[i]->SetLineColor(kGreen);
		h1_best_face_thresh15[i]->SetLineWidth(2);

		if(i+1==1){
			TLegend *leg = new TLegend(0.48,0.8,0.9,0.9);
			leg->AddEntry(h1_best_face[i],"log10(WFRMS)<3","l");
			leg->AddEntry(h1_best_face_thresh05[i],"log10(WFRMS)<-0.5","l");
			leg->AddEntry(h1_best_face_thresh1[i],"log10(WFRMS)<-1","l");
			leg->AddEntry(h1_best_face_thresh15[i],"log10(WFRMS)<-1.5","l");
			leg->Draw("same");
		}
		// gPad->SetLogy();
	}
	if(isSim)
		sprintf(title, "%s/filter_cut/%d.%d.%d_sim_A%d_c%d_%dEvents_Filter_DistBestFaces_Vpol%.1f_Hpol%.1f.png",plotPath,year_now, month_now, day_now,station,config,int(num_total),0.1*double(thresholdBin_pol[0]) + 2.0,0.1*double(thresholdBin_pol[1])+2.0);
	else
		sprintf(title, "%s/filter_cut/%d.%d.%d_data_A%d_c%d_%dEvents_Filter_DistBestFaces_Vpol%.1f_Hpol%.1f.png",plotPath,year_now, month_now, day_now,station,config,int(num_total),0.1*double(thresholdBin_pol[0]) + 2.0,0.1*double(thresholdBin_pol[1])+2.0);
	c3->SaveAs(title);
	delete c3;


	for(int i=0; i<2; i++){
		delete wfrms_plots[i];
		delete wfrms_plots_cal[i];
		delete wfrms_plots_rf[i];
		delete projections[i];
		delete projections_rf[i];
		delete projections_cal[i];
	}




	if(isSim){
		for(int pol=0; pol<2; pol++){
			for(int i=0; i<passed_events[pol]->GetNbinsX(); i++){
				double thrown = all_events[pol]->GetBinContent(i);
				double passed = passed_events[pol] -> GetBinContent(i);
				if(passed > 1E-6)
					eff[pol]->SetBinContent(i, passed/thrown);
				else
					eff[pol]->SetBinContent(i,0.);
			}
		}
		TCanvas *c2 = new TCanvas("","",1.5*1100,2*850);
		c2->Divide(2,3);
		for(int pol=0; pol<2; pol++){
			//pol 0 -> canvas 1,3,5
			//pol 1 -> canvas 2, 4, 6

			if(pol==0){
				all_events[pol]->SetTitle("Triggered Events VPol");
				passed_events[pol]->SetTitle("Passing WFRMS Filter VPol");
				eff[pol]->SetTitle("Efficiency VPol");
			}

			if(pol==1){
				all_events[pol]->SetTitle("Triggered Events HPol");
				passed_events[pol]->SetTitle("Passing WFRMS Filter HPol");
				eff[pol]->SetTitle("Efficiency HPol");
			}

			c2->cd(pol+1);
				all_events[pol]->Draw("");
				all_events[pol]->SetLineWidth(3.);
				all_events[pol]->GetYaxis()->SetRangeUser(1,1e3);
				SetAxisLabels(all_events[pol],"3rd Highest VPeak/RMS", "Counts");
				gPad->SetLogy();
			c2->cd(pol+3);
				passed_events[pol]->Draw("");
				passed_events[pol]->SetLineWidth(3.);
				passed_events[pol]->GetYaxis()->SetRangeUser(1,1e3);
				SetAxisLabels(passed_events[pol],"3rd Highest VPeak/RMS", "Counts");
				gPad->SetLogy();
			c2->cd(pol+5);
				eff[pol]->Draw("");
				eff[pol]->SetLineWidth(3.);
				SetAxisLabels(eff[pol],"3rd Highest VPeak/RMS", "Efficiency");
		}
		sprintf(title,"%s/filter_cut/%d.%d.%d_sim_A%d_c%d_%dEvents_Filter_Efficiency_Vpol%.1f_Hpol%.1f.png",plotPath,year_now, month_now, day_now,station,config,int(num_total),0.1*double(thresholdBin_pol[0]) + 2.0,0.1*double(thresholdBin_pol[1])+2.0);
		c2->SaveAs(title);
	}

}
