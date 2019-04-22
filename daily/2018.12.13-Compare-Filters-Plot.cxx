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
#include "tools_outputObjects.h"

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
		cout<< "Usage\n" << argv[0] << " <isSim?> <thresh bin> <wfrms cut val> <station> <year> <joined filename 1> <joined filename 2 > ... <joined filename x>"<<endl;
		return 0;
	}
	int isSim = atoi(argv[1]);
	int selected_bin = atoi(argv[2]);
	double selected_cut = double(atof(argv[3]));
	int station = atoi(argv[4]);
	int year = atoi(argv[5]);

	//just to have the cut parameters up front and easy to find

	// int thresholdBin_pol[]={3,5}; //bin 3 = 2.3, bin 5 = 2.5 //what is the faceRMS inclusion threshold?
	// double wavefrontRMScut[]={-1.5, -1.5}; //event wavefrontRMS < this value

	int thresholdBin_pol[]={selected_bin, selected_bin}; //bin 3 = 2.3, bin 5 = 2.5 //what is the faceRMS inclusion threshold?
	double wavefrontRMScut[]={selected_cut, selected_cut}; //event wavefrontRMS < this value

	TH2D *wfrms_plots[4];
	wfrms_plots[0] = new TH2D("Vpol_org","Vpol_org",100,-5,5,40,0,40);
	wfrms_plots[1] = new TH2D("Hpol_org","Hpol_org",100,-5,5,40,0,40);
	wfrms_plots[2] = new TH2D("Vpol_alt","Vpol_alt",100,-5,5,40,0,40);
	wfrms_plots[3] = new TH2D("Hpol_alt","Hpol_alt",100,-5,5,40,0,40);

	TH2D *wfrms_plots_cal[4];
	wfrms_plots_cal[0] = new TH2D("Vpol_org_cal","Vpol_org_cal",100,-5,5,40,0,40);
	wfrms_plots_cal[1] = new TH2D("Hpol_org_cal","Hpol_org_cal",100,-5,5,40,0,40);
	wfrms_plots_cal[2] = new TH2D("Vpol_alt_cal","Vpol_alt_cal",100,-5,5,40,0,40);
	wfrms_plots_cal[3] = new TH2D("Hpol_alt_cal","Hpol_alt_cal",100,-5,5,40,0,40);

	TH2D *wfrms_plots_rf[4];
	wfrms_plots_rf[0] = new TH2D("Vpol_org_rf","Vpol_org_rf",100,-5,5,40,0,40);
	wfrms_plots_rf[1] = new TH2D("Hpol_org_rf","Hpol_org_rf",100,-5,5,40,0,40);
	wfrms_plots_rf[2] = new TH2D("Vpol_alt_rf","Vpol_alt_rf",100,-5,5,40,0,40);
	wfrms_plots_rf[3] = new TH2D("Hpol_alt_rf","Hpol_alt_rf",100,-5,5,40,0,40);

	TH1D *all_events[2];
	TH1D *passed_events[2];
	TH1D *eff[2];
	for(int i=0; i<2; i++){
		all_events[i] = new TH1D("","",50,0,50);
		passed_events[i] = new TH1D("","",50,0,50);
		eff[i] = new TH1D("","",50,0,50);
	}

	double num_total=0.;
	double num_thermal[] = {0., 0.};
	double num_passing[] = {0., 0.};
	double num_passing_alt[] = {0.,0.};

	for(int file_num=6; file_num<argc; file_num++){

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

		int numEntries = inputTree_filter->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;

		//now to loop over events
		// numEntries=2;
		for(int event=0; event<numEntries; event++){
			inputTree_filter->GetEvent(event, weight);

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
			// if(SNRs[0]>29.) SNRs[0]=29.;
			// if(SNRs[1]>29.) SNRs[1]=29.;

			vector <double>  rms_faces_V;
			rms_faces_V.resize(12);
			vector <double> rms_faces_H;
			rms_faces_H.resize(12);

			vector <double> rms_faces_V_alt;
			vector <double> rms_faces_H_alt;

			rms_faces_V_alt.resize(numFaces_new_V);
			rms_faces_H_alt.resize(numFaces_new_H);

			//now, we must loop over the faces
			for(int i=0; i<12; i++){
				rms_faces_V[i] = rms_pol_thresh_face_V[thresholdBin_pol[0]][i];  //this is right RMS for this polarization, threshold requirement, and face
				rms_faces_H[i] = rms_pol_thresh_face_H[thresholdBin_pol[1]][i];
			}
			for(int i=0; i<numFaces_new_V; i++){
				rms_faces_V_alt[i] = rms_pol_thresh_face_alternate_V[thresholdBin_pol[0]][i];  //this is right RMS for this polarization, threshold requirement, and face
			}
			for(int i=0; i<numFaces_new_H; i++){
				rms_faces_H_alt[i] = rms_pol_thresh_face_alternate_H[thresholdBin_pol[1]][i];  //this is right RMS for this polarization, threshold requirement, and face
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

			bool trig_pol[2];
			trig_pol[0]=false;
			trig_pol[1]=false;

			int num_trig_V=0;
			int num_trig_H=0;
			for(int string=0; string<4; string++){
				for(int ant=0; ant<4; ant++){
					int chan = string + 4*ant;
					// printf("Trig pass for string %d ant %d is %d \n", string, ant, Trig_Pass[chan]);
					if(Trig_Pass[chan]>0){
						if(ant==0 || ant==2)
							num_trig_V++;
						else if(ant==1 || ant==3)
							num_trig_H++;
					}
				}
			}
			if(num_trig_V>2){
				trig_pol[0]=true;
			}
			if(num_trig_H>2){
				trig_pol[1]=true;
			}

			if(!isCalPulser && !isShort){
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

			for(int pol=0; pol<2; pol++){
				if(!isShort && !isSoftTrigger){
					if(isCalPulser){
						wfrms_plots_cal[pol]->Fill(TMath::Log10(bestFaceRMS[pol]), SNRs[pol], weight);
						wfrms_plots_cal[pol+2]->Fill(TMath::Log10(bestFaceRMS_alt[pol]), SNRs[pol], weight);
					}
					if(!isCalPulser){
						wfrms_plots[pol]->Fill(TMath::Log10(bestFaceRMS[pol]), SNRs[pol], weight);
						wfrms_plots[pol+2]->Fill(TMath::Log10(bestFaceRMS_alt[pol]), SNRs[pol], weight);
						wfrms_plots_rf[pol]->Fill(TMath::Log10(bestFaceRMS[pol]), SNRs[pol], weight);
						wfrms_plots_rf[pol+2]->Fill(TMath::Log10(bestFaceRMS[pol]), SNRs[pol], weight);						
						if(TMath::Log10(bestFaceRMS[pol]) < wavefrontRMScut[pol]){
							num_passing[pol]+=(weight);
						}
						// if(TMath::Log10(bestFaceRMS_alt[pol]) < wavefrontRMScut[pol] && !trig_pol[pol]){
						// 	cout<<"Passes WFRMS but didn't trigger in pol "<<pol<<" with weight "<<weight<<endl;
						// }
						if(TMath::Log10(bestFaceRMS_alt[pol]) < wavefrontRMScut[pol]){
							num_passing_alt[pol]+=(weight);
							if(isSim)
								passed_events[pol]->Fill(SNRs[pol],weight);
								// passed_events[pol]->Fill(SNR_theory/0.035, weight);
						}
					}
				}
			}//loop over polarization
		
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
		// printf("	Org face: Num thermal passing pol %d: %.3f events, %.3f rate \n", pol, num_passing[pol], 100.*num_passing[pol]/num_thermal);
		printf("	Alt  face: Num thermal passing pol %d: %.3f/%.3f events, %.3f rate \n", pol, num_passing_alt[pol],num_thermal[pol], 100.*num_passing_alt[pol]/num_thermal[pol]);
	}

	TH1D *projections[4];
	TH1D *projections_cal[4];
	TH1D *projections_rf[4];

	for(int pol=0; pol<2; pol++){
		projections[pol] = wfrms_plots[pol]->ProjectionX();
		projections[pol+2] = wfrms_plots[pol+2]->ProjectionX();
		projections_cal[pol] = wfrms_plots_cal[pol]->ProjectionX();
		projections_cal[pol+2] = wfrms_plots_cal[pol+2]->ProjectionX();
		projections_rf[pol] = wfrms_plots_rf[pol]->ProjectionX();
		projections_rf[pol+2] = wfrms_plots_rf[pol+2]->ProjectionX();
	}

	// for(int pol=0; pol<2; pol++){
	// 	projections[pol] = wfrms_plots[pol]->ProjectionX();
	// 	projections[pol+2] = wfrms_plots[pol+2]->ProjectionX();
	// }

	TCanvas *c = new TCanvas("","",2.1*850,3.1*850);
	c->Divide(2,4);
	for(int pol=0; pol<2; pol++){
		c->cd(pol+1);
			wfrms_plots[pol]->Draw("colz");
			wfrms_plots[pol]->GetXaxis()->SetTitle("log10(Wavefront RMS)");
			wfrms_plots[pol]->GetYaxis()->SetTitle("3rd Highest VPeak/RMS");
			wfrms_plots[pol]->GetZaxis()->SetRangeUser(0.1,7e6);
			gPad->SetLogz();
		c->cd(pol+3);		
			wfrms_plots[pol+2]->Draw("colz");
			wfrms_plots[pol+2]->GetXaxis()->SetTitle("log10(Wavefront RMS)");
			wfrms_plots[pol+2]->GetYaxis()->SetTitle("3rd Highest VPeak/RMS");
			wfrms_plots[pol+2]->GetZaxis()->SetRangeUser(0.1,7e6);
			gPad->SetLogz();
		c->cd(pol+5);		
			projections[pol]->Draw();
			projections[pol]->SetLineWidth(3);
			// if(pol==0 || pol==1){
			// 	projections_cal[pol]->Draw("same");
			// 	projections_cal[pol]->SetLineWidth(2);
			// 	projections_cal[pol]->SetLineColor(kRed);
			// }
			// projections_rf[pol]->Draw("same");
			// projections_rf[pol]->SetLineWidth(2);
			// projections_rf[pol]->SetLineColor(kRed);
			projections[pol]->GetYaxis()->SetTitle("Counts");
			projections[pol]->GetXaxis()->SetTitle("log10(Wavefront RMS)");
			projections[pol]->GetYaxis()->SetRangeUser(0.1,1e7);
			gPad->SetLogy();
		c->cd(pol+7);		
			projections[pol+2]->Draw();
			projections[pol+2]->SetLineWidth(3);
			// if(pol==0 || pol==1){
			// 	projections_cal[pol+2]->Draw("same");
			// 	projections_cal[pol+2]->SetLineWidth(2);
			// 	projections_cal[pol+2]->SetLineColor(kRed);
			// }
			// projections_rf[pol+2]->Draw("same");
			// projections_rf[pol+2]->SetLineWidth(2);
			// projections_rf[pol+2]->SetLineColor(kRed);
			projections[pol+2]->GetYaxis()->SetTitle("Counts");
			projections[pol+2]->GetXaxis()->SetTitle("log10(Wavefront RMS)");
			projections[pol+2]->GetYaxis()->SetRangeUser(0.1,1e7);
			gPad->SetLogy();
	}
	char title[300];
	if(isSim) sprintf(title, "/users/PAS0654/osu0673/A23_analysis_new2/results/%d.%d.%d_sim_A%d_%d_%dEvents_CompareFilters_Vpol%.1f_Hpol%.1f.png",year_now, month_now, day_now,station,year,int(num_total),0.1*double(thresholdBin_pol[0]) + 2.0,0.1*double(thresholdBin_pol[1])+2.0);
	sprintf(title, "/users/PAS0654/osu0673/A23_analysis_new2/results/%d.%d.%d_data_A%d_%d_%dEvents_CompareFilters_Vpol%.1f_Hpol%.1f.png",year_now, month_now, day_now,station,year,int(num_total),0.1*double(thresholdBin_pol[0]) + 2.0,0.1*double(thresholdBin_pol[1])+2.0);
	c->SaveAs(title);
	delete c;
	delete wfrms_plots[0]; delete wfrms_plots[1]; delete wfrms_plots[2]; delete wfrms_plots[3];
	for(int i=0; i<4; i++) delete projections[i];

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
		c2->SaveAs("/users/PAS0654/osu0673/A23_analysis_new2/results/all_events.png");
	}

}
