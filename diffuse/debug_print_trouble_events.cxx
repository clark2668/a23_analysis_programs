////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	v2_final_plots.cxx 
////	A23 diffuse, make plots of the final cut parameter space
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
#include "TF1.h"

//AraRoot includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "Settings.h"
#include "Detector.h"
#include "Report.h"
#include "RayTraceCorrelator.h"
AraAntPol::AraAntPol_t Vpol = AraAntPol::kVertical;
AraAntPol::AraAntPol_t Hpol = AraAntPol::kHorizontal;
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_Cuts.h"

using namespace std;

int PlotThisEvent(int station, int year, int runNum, int event, Settings *settings, Detector *detector, RayTraceCorrelator *theCorrelators[2]);

int main(int argc, char **argv)
{

	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <station> <year> <ValForCuts filename>"<<endl;;
		return -1;
	}
	int station = atoi(argv[1]);
	int year = atoi(argv[2]);
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

	stringstream ss;

	gStyle->SetOptStat(11);

	//set up the ray tracer
	Settings *settings = new Settings();
	string setupfile = "setup.txt";
	settings->ReadFile(setupfile);
	cout << "Read " << setupfile << " file!" << endl;
	settings->NOFZ=1;
	Detector *detector=0;
	RayTraceCorrelator *theCorrelators[2];
	theCorrelators[0] =  new RayTraceCorrelator(station, 41., settings, 1, 4); //41 m, cal puser
	theCorrelators[1] =  new RayTraceCorrelator(station, 300., settings, 1, 4);//300 m, far reco
	
	TH2D *PeakCorr_vs_SNR_all[2];
	PeakCorr_vs_SNR_all[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_all[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *PeakCorr_vs_SNR_cutCal[2];
	PeakCorr_vs_SNR_cutCal[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_cutCal[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft[2];
	PeakCorr_vs_SNR_cutCal_cutSoft[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_cutCal_cutSoft[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[2];
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[2];
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[2];
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[2];
	// PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[0]=new TH2D("","V",90,0,30,500,0,1);
	// PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[1]=new TH2D("","H",90,0,30,500,0,1);
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *special[2];
	special[0]=new TH2D("","V",90,0,30,500,0,1);
	special[1]=new TH2D("","H",90,0,30,500,0,1);

	TH1D *fracs_power_cut[2];
	fracs_power_cut[0]=new TH1D("","V",100,0,1);
	fracs_power_cut[1]=new TH1D("","H",100,0,1);
	
	int num_total=0;
	int num_in_final_plot=0;
	int num_refilt=0;

	for(int file_num=3; file_num<argc; file_num++){

		string chRun = "run";
		string file = string(argv[file_num]);
		size_t foundRun = file.find(chRun);
		string strRunNum = file.substr(foundRun+4,4);
		int runNum = atoi(strRunNum.c_str());

		if(isBadRun(station,year,runNum)) continue;

		TFile *inputFile = TFile::Open(argv[file_num]);
		if(!inputFile){
			cout<<"Can't open joined file!"<<endl;
			return -1;
		}
		cout << "Run " << file_num << " :: " << argv[file_num] << endl;

		TTree *trees[3];
		trees[0] = (TTree*) inputFile->Get("VTree");
		trees[1] = (TTree*) inputFile->Get("HTree");
		trees[2] = (TTree*) inputFile->Get("AllTree");

		double corr_val[2];
		double snr_val[2];
		int WFRMS[2];
		double frac_of_power_notched_V[8];
		double frac_of_power_notched_H[8];
		int Refilt[2];

		trees[0]->SetBranchAddress("corr_val_V",&corr_val[0]);
		trees[0]->SetBranchAddress("snr_val_V",&snr_val[0]);
		trees[0]->SetBranchAddress("wfrms_val_V",&WFRMS[0]);
		trees[0]->SetBranchAddress("Refilt_V",&Refilt[0]);
		trees[1]->SetBranchAddress("corr_val_H",&corr_val[1]);
		trees[1]->SetBranchAddress("snr_val_H",&snr_val[1]);
		trees[1]->SetBranchAddress("wfrms_val_H",&WFRMS[1]);
		trees[1]->SetBranchAddress("Refilt_H",&Refilt[1]);

		int isCal;
		int isSoft;
		int isShort;
		int isCW;
		int isNewBox;
		int isSurf;
		int isBadEvent;

		trees[2]->SetBranchAddress("cal",&isCal);
		trees[2]->SetBranchAddress("soft",&isSoft);
		trees[2]->SetBranchAddress("short",&isShort);
		trees[2]->SetBranchAddress("CW",&isCW);
		trees[2]->SetBranchAddress("box",&isNewBox);
		trees[2]->SetBranchAddress("surf",&isSurf);
		trees[2]->SetBranchAddress("bad",&isBadEvent);

		stringstream ss;
		for(int i=0; i<8; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			trees[0]->SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_V[i]);
		}
		for(int i=8; i<16; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			trees[1]->SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_H[i-8]);
		}
		
		int numEntries = trees[0]->GetEntries();

		//now to loop over events
		for(int event=0; event<numEntries; event++){

			trees[0]->GetEvent(event);
			trees[1]->GetEvent(event);
			trees[2]->GetEvent(event);

			num_total++;

			if(isBadEvent){
				continue;
			}

			for(int pol=0; pol<2; pol++){
				PeakCorr_vs_SNR_all[pol]->Fill(snr_val[pol],corr_val[pol]);
				
				if(!isCal){ //cut cal pulsers
					PeakCorr_vs_SNR_cutCal[pol]->Fill(snr_val[pol],corr_val[pol]);
					
					if(!isSoft){ //cut software triggers 
						PeakCorr_vs_SNR_cutCal_cutSoft[pol]->Fill(snr_val[pol],corr_val[pol]);
						
						if(!isShort){ //cut short
							PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->Fill(snr_val[pol],corr_val[pol]);
							
							if(!WFRMS[pol]){ //cut WRMS
								PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->Fill(snr_val[pol],corr_val[pol]);
								
								if(!isNewBox){ //cut cal box
									PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->Fill(snr_val[pol],corr_val[pol]);

									if(!isSurf){

										bool condition = false;
										if(snr_val[pol]>=28.) condition=true;

										if(Refilt[pol]){
											num_refilt++;

											vector<double> frac;
											for(int i=0; i<8; i++){
												if(pol==0) frac.push_back(frac_of_power_notched_V[i]);
												else if(pol==1) frac.push_back(frac_of_power_notched_H[i]);
											}
											sort(frac.begin(), frac.end(), std::greater<double>());
											fracs_power_cut[pol]->Fill(frac[2]);
											if(frac[2]<=0.06){
												PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->Fill(snr_val[pol],corr_val[pol]);
												if(condition){
													// printf("		Event %d is refiltered in pol %d \n", event,pol);
													// cout<<"			Frac of power notched is "<<frac[2]<<endl;
													PlotThisEvent(station,year,runNum,event, settings, detector, theCorrelators);
												}
											}
										} //refiltered?
										else{
											PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->Fill(snr_val[pol],corr_val[pol]);
											if(condition){
												// printf("		Event %d is NOT refiltered in pol %d \n", event,pol);
												PlotThisEvent(station,year,runNum,event, settings, detector, theCorrelators);
											}
										}
										num_in_final_plot++;
									}
								}
							}
						}
					}
				}
			}
		}
		inputFile->Close();
		delete inputFile;
	}
	cout<<"Num total is "<<num_total<<endl;
	cout<<"Num in final plot "<<num_in_final_plot<<endl;
	cout<<"Num re-filtered is "<<num_refilt<<endl;

	gStyle->SetOptStat(11);
	gStyle->SetStatY(0.9);
	gStyle->SetStatX(0.9);
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.2);

	//save out SNR vs WavefrontRMS plot
	char graph_title[2][300];
	char title[300];

	int cal=0;
	int soft=0;
	int Short=0;
	int wrms=0;
	int box = 0;
	int surf=0;
	int cw=0;

	
	//save out the Corr vs SNR plot for all 
	// sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	// sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	// TCanvas *c2 = new TCanvas("","",2.1*850,850);
	// c2->Divide(2,1);
	// for(int pol=0; pol<2; pol++){
	// 	c2->cd(pol+1);
	// 	PeakCorr_vs_SNR_all[pol]->Draw("colz");
	// 	PeakCorr_vs_SNR_all[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
	// 	PeakCorr_vs_SNR_all[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
	// 	PeakCorr_vs_SNR_all[pol]->SetTitle(graph_title[pol]);
	// 	gPad->SetLogz();
	// }
	// sprintf(title, "%s/%d.%d.%d_A%d_%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,year,num_total,cal,soft,Short,wrms,box,surf);
	// c2->SaveAs(title);
	// delete c2;
	// delete PeakCorr_vs_SNR_all[0]; delete PeakCorr_vs_SNR_all[1];

	// //turn on cal
	// cal=1;
	// sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	// sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	// TCanvas *c3 = new TCanvas("","",2.1*850,850);
	// c3->Divide(2,1);
	// for(int pol=0; pol<2; pol++){
	// 	c3->cd(pol+1);
	// 	PeakCorr_vs_SNR_cutCal[pol]->Draw("colz");
	// 	PeakCorr_vs_SNR_cutCal[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
	// 	PeakCorr_vs_SNR_cutCal[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
	// 	PeakCorr_vs_SNR_cutCal[pol]->SetTitle(graph_title[pol]);
	// 	gPad->SetLogz();
	// }
	// sprintf(title, "%s/%d.%d.%d_A%d_%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,year,num_total,cal,soft,Short,wrms,box,surf);
	// c3->SaveAs(title);
	// delete c3;
	// delete PeakCorr_vs_SNR_cutCal[0]; delete PeakCorr_vs_SNR_cutCal[1];

	// //turn on cal, soft
	// soft=1;
	// sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	// sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	// TCanvas *c4 = new TCanvas("","",2.1*850,850);
	// c4->Divide(2,1);
	// for(int pol=0; pol<2; pol++){
	// 	c4->cd(pol+1);
	// 	PeakCorr_vs_SNR_cutCal_cutSoft[pol]->Draw("colz");
	// 	PeakCorr_vs_SNR_cutCal_cutSoft[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
	// 	PeakCorr_vs_SNR_cutCal_cutSoft[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
	// 	PeakCorr_vs_SNR_cutCal_cutSoft[pol]->SetTitle(graph_title[pol]);
	// 	gPad->SetLogz();
	// }
	// sprintf(title, "%s/%d.%d.%d_A%d_%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,year,num_total,cal,soft,Short,wrms,box,surf);
	// c4->SaveAs(title);
	// delete c4;
	// delete PeakCorr_vs_SNR_cutCal_cutSoft[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft[1];

	// //turn on cal, soft, short
	// Short=1;
	// sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	// sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	// TCanvas *c5 = new TCanvas("","",2.1*850,850);
	// c5->Divide(2,1);
	// for(int pol=0; pol<2; pol++){
	// 	c5->cd(pol+1);
	// 	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->Draw("colz");
	// 	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
	// 	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
	// 	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->SetTitle(graph_title[pol]);
	// 	gPad->SetLogz();
	// }
	// sprintf(title, "%s/%d.%d.%d_A%d_%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,year,num_total,cal,soft,Short,wrms,box,surf);
	// c5->SaveAs(title);
	// delete c5;
	// delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[1];

	// //turn on cal, soft, short, wmrs
	// wrms=1;
	// sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	// sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	// TCanvas *c6 = new TCanvas("","",2.1*850,850);
	// c6->Divide(2,1);
	// for(int pol=0; pol<2; pol++){
	// 	c6->cd(pol+1);
	// 	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->Draw("colz");
	// 	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
	// 	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
	// 	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->SetTitle(graph_title[pol]);
	// 	gPad->SetLogz();
	// }
	// sprintf(title, "%s/%d.%d.%d_A%d_%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,year,num_total,cal,soft,Short,wrms,box,surf);
	// c6->SaveAs(title);
	// delete c6;
	// delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[1];

	// //turn on cal, soft, short, wmrs, box
	// box=1;
	// sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	// sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	// TCanvas *c7 = new TCanvas("","",2.1*850,850);
	// c7->Divide(2,1);
	// for(int pol=0; pol<2; pol++){
	// 	c7->cd(pol+1);
	// 	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->Draw("colz");
	// 	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
	// 	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
	// 	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->SetTitle(graph_title[pol]);
	// 	gPad->SetLogz();
	// }
	// sprintf(title, "%s/%d.%d.%d_A%d_%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,year,num_total,cal,soft,Short,wrms,box,surf);
	// c7->SaveAs(title);
	// delete c7;
	// delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[1];

	// surf=1;
	// sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	// sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	// TCanvas *c8 = new TCanvas("","",2.1*850,850);
	// c8->Divide(2,1);
	// for(int pol=0; pol<2; pol++){
	// 	c8->cd(pol+1);
	// 	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->Draw("colz");
	// 	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
	// 	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
	// 	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->SetTitle(graph_title[pol]);
	// 	gPad->SetLogz();
	// 	// PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->GetXaxis()->SetRangeUser(0,10);
	// 	// PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->GetYaxis()->SetRangeUser(0,0.5);
	// }
	// sprintf(title, "%s/%d.%d.%d_A%d_%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,year,num_total,cal,soft,Short,wrms,box,surf);
	// c8->SaveAs(title);
	// delete c8;
	// delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[1];
	

}

int PlotThisEvent(int station, int year, int runNum, int event, Settings *settings, Detector *detector, RayTraceCorrelator *theCorrelators[2]){
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *DataDirPath(getenv("DATA_DIR"));
	if (DataDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

	char run_file_name[400];
	if(year==2013){
		sprintf(run_file_name,"%s/RawData/A%d/%d/run%d/event%d.root",DataDirPath,station,year,runNum,runNum);
	}
	else if(year==2014 || year==2015 || year==2016){
		sprintf(run_file_name,"%s/RawData/A%d/%d/sym_links/event00%d.root",DataDirPath,station,year,runNum,runNum);
	}
	TFile *mapFile = TFile::Open(run_file_name);
	if(!mapFile){
		cout<<"Can't open data file for map!"<<endl;
		return -1;
	}
	TTree *eventTree = (TTree*) mapFile-> Get("eventTree");
	if(!eventTree){
		cout<<"Can't find eventTree for map"<<endl;
		return -1;
	}

	RawAtriStationEvent *rawPtr =0;
	eventTree->SetBranchAddress("event",&rawPtr);
	eventTree->GetEvent(event);

	int stationID = rawPtr->stationId;
	char ped_file_name[400];

	if(year==2013){
		sprintf(ped_file_name,"%s/run_specific_peds/A%d/%d/event%d_specificPeds.dat",PedDirPath,station,year,runNum);
	}
	else if(year==2014 || year==2015 || year==2016){
		sprintf(ped_file_name,"%s/run_specific_peds/A%d/%d/event00%d_specificPeds.dat",PedDirPath,station,year,runNum);
	}
	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	calibrator->setAtriPedFile(ped_file_name,stationID); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist

	UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawPtr,AraCalType::kLatestCalib);

	int unixTime = (int)rawPtr->unixTime;
	int unixTimeUs =(int)rawPtr->unixTimeUs;
	int timeStamp = (int)rawPtr->timeStamp;
	printf("Unixtime is %d \n", unixTime);
	printf("Unixtime microsecond is %d \n", unixTimeUs);
	printf("timeStamp is %d \n", timeStamp);

	stringstream ss1;
	string xLabel, yLabel;
	xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
	vector<string> titlesForGraphs;
	for (int i = 0; i < 16; i++){
		ss1.str("");
		ss1 << "Channel " << i;
		titlesForGraphs.push_back(ss1.str());
	}

	char cw_file_name[400];
	sprintf(cw_file_name,"%s/CWID/A%d/%d/CWID_station_%d_run_%d.root",DataDirPath,station,year,station,runNum);
	TFile *NewCWFile = TFile::Open(cw_file_name);
	if(!NewCWFile) {
		std::cerr << "Can't open new CW file\n";
		return -1;
	}
	TTree* NewCWTree = (TTree*) NewCWFile->Get("NewCWTree");   
	if(!NewCWTree) {
		std::cerr << "Can't find NewCWTree\n";
		return -1;
	}
	vector<vector<double> > *badFreqs_fwd =0;
	vector<vector<double> > *badFreqs_back=0;
	vector<vector<double> > *badSigmas_fwd=0;
	vector<vector<double> > *badSigmas_back=0;
	vector<vector<double> > *badFreqs_baseline=0;

	NewCWTree->SetBranchAddress("badFreqs_fwd",&badFreqs_fwd);
	NewCWTree->SetBranchAddress("badSigmas_fwd",&badSigmas_fwd);
	NewCWTree->SetBranchAddress("badFreqs_back",&badFreqs_back);
	NewCWTree->SetBranchAddress("badSigmas_back",&badSigmas_back);
	NewCWTree->SetBranchAddress("badFreqs_baseline",&badFreqs_baseline);

	//deal w/ CW cut
	//inputTree_CW->GetEntry(event);
	NewCWTree->GetEntry(event);

	bool isCutonCW_fwd[2]; isCutonCW_fwd[0]=false; isCutonCW_fwd[1]=false;
	bool isCutonCW_back[2]; isCutonCW_back[0]=false; isCutonCW_back[1]=false;
	bool isCutonCW_baseline[2]; isCutonCW_baseline[0]=false; isCutonCW_baseline[1]=false;
	
	for(int pol=0; pol<badFreqs_baseline->size(); pol++){
		vector<double> badFreqListLocal_baseline = badFreqs_baseline->at(pol);
		if(badFreqListLocal_baseline.size()>0) isCutonCW_baseline[pol]=true;
	}

	double threshCW = 1.0;
	vector<double> badFreqList_fwd;
	vector<double> badSigmaList_fwd;
	for(int pol=0; pol<badFreqs_fwd->size(); pol++){
		badFreqList_fwd=badFreqs_fwd->at(pol);
		badSigmaList_fwd=badSigmas_fwd->at(pol);
		for(int iCW=0; iCW<badFreqList_fwd.size(); iCW++){
			if(
				badSigmaList_fwd[iCW] > threshCW 
				&& 
				abs(300. - badFreqList_fwd[iCW]) > 2.
				&&
				abs(500. - badFreqList_fwd[iCW]) > 2.
				&&
				abs(125. - badFreqList_fwd[iCW]) > 2.
			){
				isCutonCW_fwd[pol] = true;
			}
		}
	}
	vector<double> badFreqList_back;
	vector<double> badSigmaList_back;
	for(int pol=0; pol<badFreqs_back->size(); pol++){
		badFreqList_back=badFreqs_back->at(pol);
		badSigmaList_back=badSigmas_back->at(pol);
		for(int iCW=0; iCW<badFreqList_back.size(); iCW++){
			if(
				badSigmaList_back[iCW] > threshCW 
				&& 
				abs(300. - badFreqList_back[iCW]) > 2.
				&&
				abs(500. - badFreqList_back[iCW]) > 2.
				&&
				abs(125. - badFreqList_back[iCW]) > 2.
			){
				isCutonCW_back[pol] = true;
			}
		}
	}
	for(int pol=0; pol<2; pol++){
		if(isCutonCW_fwd[pol] || isCutonCW_back[pol] || isCutonCW_baseline[pol]){
			printf("Has CW issue in pol %d \n", pol);
			printf("CW in FWD %d, BWD %d, or baseline %d? \n", isCutonCW_fwd[pol], isCutonCW_back[pol], isCutonCW_baseline[pol]);
			//get the frequencies to notch
			vector<double> badFreqListLocal_fwd;
			vector <double> badFreqListLocal_back;
			vector <double> mergedFreqList;

			//merge the two lists of frequencies
			//if it's cut going both forward and backward
			if(isCutonCW_fwd[pol] && isCutonCW_back[pol]){
				badFreqListLocal_fwd=badFreqs_fwd->at(pol);
				badFreqListLocal_back=badFreqs_back->at(pol);
				for(int iFreq=0; iFreq<badFreqListLocal_fwd.size(); iFreq++){
					mergedFreqList.push_back(badFreqListLocal_fwd[iFreq]);
				}
				for(int iFreq=0; iFreq<badFreqListLocal_back.size(); iFreq++){
					double new_freq=badFreqListLocal_back[iFreq];
					for(int iFreqOld=0; iFreqOld<badFreqListLocal_fwd.size(); iFreqOld++){
						if(abs(new_freq-mergedFreqList[iFreqOld])>0.1){
							mergedFreqList.push_back(new_freq);
						}
					}
				}
			}
			//if it's cut only going forward
			else if(isCutonCW_fwd[pol] && !isCutonCW_back[pol]){
				badFreqListLocal_fwd=badFreqs_fwd->at(pol);
				for(int iFreq=0; iFreq<badFreqListLocal_fwd.size(); iFreq++){
					mergedFreqList.push_back(badFreqListLocal_fwd[iFreq]);
				}
			}
			//if it's cut only going backward
			else if(!isCutonCW_fwd[pol] && isCutonCW_back[pol]){
				badFreqListLocal_back=badFreqs_back->at(pol);
				for(int iFreq=0; iFreq<badFreqListLocal_back.size(); iFreq++){
					mergedFreqList.push_back(badFreqListLocal_back[iFreq]);
				}
			}

			vector<double> more_freqs_to_add;
			vector<double> badFreqListLocal_baseline = badFreqs_baseline->at(pol);
			if(mergedFreqList.size()>0){ //do we already have frequencies to check against?
				//loop over everything identified by the CW baseline cut
				for(int newFreq=0; newFreq<badFreqListLocal_baseline.size(); newFreq++){
					double new_freq = badFreqListLocal_baseline[newFreq];
					//now, loop over everything already in the list
					for(int oldFreq=0; oldFreq<mergedFreqList.size(); oldFreq++){
						//if there's a genuinely new frequency, add it to the list of things to be adde
						if(abs(new_freq-mergedFreqList[oldFreq])>0.1){
							more_freqs_to_add.push_back(new_freq);
						}
					}
				}
			}
			else{ //otherwise we take only those found by the CW ID cut
				for(int newFreq=0; newFreq<badFreqListLocal_baseline.size(); newFreq++){
					double new_freq = badFreqListLocal_baseline[newFreq];
					more_freqs_to_add.push_back(new_freq);
				}
			}

			//now actually add it to the merged freq list
			for(int iFreq=0; iFreq<more_freqs_to_add.size(); iFreq++){
				mergedFreqList.push_back(more_freqs_to_add[iFreq]);
			}

			//they need to be in smaller -> larger order for notching
			sort(mergedFreqList.begin(), mergedFreqList.end());
			vector <double> uniqueNotchFreqs;
			vector <double> uniqueNotchBands;			
			theCorrelators[0]->pickFreqsAndBands(mergedFreqList,uniqueNotchFreqs,uniqueNotchBands);
			for (int i = 0; i < uniqueNotchFreqs.size(); ++i)
			{
				printf("Unique freq to be notched is %.2f with width %.2f \n", uniqueNotchFreqs[i],uniqueNotchBands[i]);
			}
			vector <TGraph*> grWaveformsRaw = makeGraphsFromRF(realAtriEvPtr,16,xLabel,yLabel,titlesForGraphs);
			vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(grWaveformsRaw, 0.6, xLabel, yLabel, titlesForGraphs);
			vector<TGraph*> grWaveformsPadded = makePaddedGraphs(grWaveformsInt, 0, xLabel, yLabel, titlesForGraphs);
			vector <TGraph*> grNotched;
			for(int i=0; i<16; i++){
				TGraph *grNotchAmp = theCorrelators[0]->applyAdaptiveFilter_singleAnt_FiltMany(grWaveformsPadded[i],uniqueNotchFreqs,uniqueNotchBands);
				grNotched.push_back(theCorrelators[0]->GeometricFilter(grNotchAmp,uniqueNotchFreqs,uniqueNotchBands,uniqueNotchFreqs));
				delete grNotchAmp;
			}
			vector<TGraph*> grWaveformsPowerSpectrum = makePowerSpectrumGraphs(grWaveformsPadded, xLabel, yLabel, titlesForGraphs);
			vector<TGraph*> grWaveformsPowerSpectrum_notched = makePowerSpectrumGraphs(grNotched, xLabel, yLabel, titlesForGraphs);
			
			char save_temp_title[300];
			sprintf(save_temp_title,"%s/trouble_events/%d.%d.%d_Run%d_Ev%d_WaveformsNotch.png",plotPath,year_now,month_now,day_now,runNum,event);
			TCanvas *cWave = new TCanvas("","",4*1100,4*850);
			cWave->Divide(4,4);
			for(int i=0; i<16; i++){
				cWave->cd(i+1);
				grWaveformsRaw[i]->Draw("AL");
				grWaveformsRaw[i]->SetLineWidth(3);
				grNotched[i]->Draw("same");
				grNotched[i]->SetLineWidth(3);
				grNotched[i]->SetLineColor(kRed);
			}
			cWave->SaveAs(save_temp_title);
			delete cWave;

			sprintf(save_temp_title,"%s/trouble_events/%d.%d.%d_Run%d_Ev%d_SpectraNotch.png",plotPath,year_now,month_now,day_now,runNum,event);
			TCanvas *cSpec = new TCanvas("","",4*1100,4*850);
			cSpec->Divide(4,4);
			for(int i=0; i<16; i++){
				cSpec->cd(i+1);
				grWaveformsPowerSpectrum[i]->Draw("AL");
				grWaveformsPowerSpectrum[i]->SetLineWidth(3);
				grWaveformsPowerSpectrum_notched[i]->Draw("same");
				grWaveformsPowerSpectrum_notched[i]->SetLineWidth(3);
				grWaveformsPowerSpectrum_notched[i]->SetLineColor(kRed);
				gPad->SetLogy();
			}
			cSpec->SaveAs(save_temp_title);
			delete cSpec;
			deleteGraphVector(grWaveformsRaw);
			deleteGraphVector(grWaveformsInt);
			deleteGraphVector(grWaveformsPadded);
			deleteGraphVector(grNotched);
			deleteGraphVector(grWaveformsPowerSpectrum);
			deleteGraphVector(grWaveformsPowerSpectrum_notched);
		}
	}

	vector <TGraph*> waveforms = makeGraphsFromRF(realAtriEvPtr,16,xLabel,yLabel,titlesForGraphs);
	vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(waveforms, 0.5, xLabel, yLabel, titlesForGraphs);
	vector<TGraph*> grWaveformsPadded = makePaddedGraphs(grWaveformsInt, 0, xLabel, yLabel, titlesForGraphs);
	xLabel = "Frequency (Hz)"; yLabel = "Power Spectral Density (mV/Hz)";
	vector<TGraph*> grWaveformsPowerSpectrum = makePowerSpectrumGraphs(grWaveformsPadded, xLabel, yLabel, titlesForGraphs);

	bool do_reco=true;
	if(do_reco){
		TH2D *map_30m_V;
		TH2D *map_300m_V;
		TH2D *map_30m_H;
		TH2D *map_300m_H;
		TH2D *map_30m_V_select;

		vector <int> chan_list_V;
		vector <int> chan_list_H;
		for(int chan=0; chan<=7; chan++){
			chan_list_V.push_back(chan);
			chan_list_H.push_back(chan+8);
		}

		if(station==2){
			//for station 2, we need to exclude channel 15 from the analysis
			chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
		}
		else if(station==3){
			//for station 3 years 2014 and 2015, we need to drop string 4 (channels 3, 7, 11, 15) altogether
			if(year==2014 || year==2015 || year==2016){
				chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 3), chan_list_V.end());
				chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 7), chan_list_V.end());

				chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 11), chan_list_H.end());
				chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
			}
		}

		// chan_list_V.clear();
		// chan_list_V.push_back(0);
		// chan_list_V.push_back(1);
		// chan_list_V.push_back(2);
		// chan_list_V.push_back(3);

		// chan_list_H.clear();
		// chan_list_H.push_back(8);
		// chan_list_H.push_back(9);
		// chan_list_H.push_back(10);
		// chan_list_H.push_back(11);

		map_30m_V = theCorrelators[0]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Vpol, false,chan_list_V) ;
		map_300m_V = theCorrelators[1]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Vpol, false,chan_list_V);
		map_30m_H = theCorrelators[0]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Hpol, false,chan_list_H);
		map_300m_H = theCorrelators[1]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Hpol, false,chan_list_H);

		int PeakTheta_Recompute_30m_H;
		int PeakTheta_Recompute_300m_H;
		int PeakPhi_Recompute_30m_H;
		int PeakPhi_Recompute_300m_H;
		double PeakCorr_Recompute_30m_H;
		double PeakCorr_Recompute_300m_H;
		int PeakTheta_Recompute_30m_V;
		int PeakTheta_Recompute_300m_V;
		int PeakPhi_Recompute_30m_V;
		int PeakPhi_Recompute_300m_V;
		double PeakCorr_Recompute_30m_V;
		double PeakCorr_Recompute_300m_V;
		getCorrMapPeak(map_30m_H,PeakTheta_Recompute_30m_H,PeakPhi_Recompute_30m_H,PeakCorr_Recompute_30m_H);
		getCorrMapPeak(map_300m_H,PeakTheta_Recompute_300m_H,PeakPhi_Recompute_300m_H,PeakCorr_Recompute_300m_H);
		getCorrMapPeak(map_30m_V,PeakTheta_Recompute_30m_V,PeakPhi_Recompute_30m_V,PeakCorr_Recompute_30m_V);
		getCorrMapPeak(map_300m_V,PeakTheta_Recompute_300m_V,PeakPhi_Recompute_300m_V,PeakCorr_Recompute_300m_V);

		printf("30m H theta and phi %d and %d \n", PeakTheta_Recompute_30m_H, PeakPhi_Recompute_30m_H);
		stringstream ss30H;
		ss30H<<" Peak Theta, Phi is "<<PeakTheta_Recompute_30m_H<<" , "<<PeakPhi_Recompute_30m_H;
		map_30m_H->SetTitle(ss30H.str().c_str());
		printf("300m H theta and phi %d and %d \n", PeakTheta_Recompute_300m_H, PeakPhi_Recompute_300m_H);
		stringstream ss300H;
		ss30H<<" Peak Theta, Phi is "<<PeakTheta_Recompute_300m_H<<" , "<<PeakPhi_Recompute_300m_H;
		map_300m_H->SetTitle(ss300H.str().c_str());
		printf("30m V theta and phi %d and %d \n", PeakTheta_Recompute_30m_V, PeakPhi_Recompute_30m_V);
		stringstream ss30V;
		ss30V<<" Peak Theta, Phi is "<<PeakTheta_Recompute_30m_V<<" , "<<PeakPhi_Recompute_30m_V;
		map_30m_V->SetTitle(ss30V.str().c_str());
		printf("300m V theta and phi %d and %d \n", PeakTheta_Recompute_300m_V, PeakPhi_Recompute_300m_V);
		stringstream ss300V;
		ss300V<<" Peak Theta, Phi is "<<PeakTheta_Recompute_300m_V<<" , "<<PeakPhi_Recompute_300m_V;
		map_300m_V->SetTitle(ss300V.str().c_str());

		TCanvas *cMaps = new TCanvas("","",2*1100,2*850);
		cMaps->Divide(2,2);
			cMaps->cd(3);
			map_30m_V->Draw("colz");
			cMaps->cd(4);
			map_30m_H->Draw("colz");
			cMaps->cd(1);
			map_300m_V->Draw("colz");
			cMaps->cd(2);
			map_300m_H->Draw("colz");
		char save_temp_title[400];		
		sprintf(save_temp_title,"%s/trouble_events/%d.%d.%d_Run%d_Ev%d_Maps.png",plotPath,year_now,month_now,day_now,runNum,event);
		cMaps->SaveAs(save_temp_title);
		delete cMaps;
		delete map_30m_V; delete map_300m_V; delete map_30m_H; delete map_300m_H; 
	}


	char save_temp_title[300];
	sprintf(save_temp_title,"%s/trouble_events/%d.%d.%d_Run%d_Ev%d_Waveforms.png",plotPath,year_now,month_now,day_now,runNum,event);
	TCanvas *cWave = new TCanvas("","",4*1100,4*850);
	cWave->Divide(4,4);
	for(int i=0; i<16; i++){
		cWave->cd(i+1);
		waveforms[i]->Draw("AL");
		waveforms[i]->SetLineWidth(3);
	}
	cWave->SaveAs(save_temp_title);
	delete cWave;

	sprintf(save_temp_title,"%s/trouble_events/%d.%d.%d_Run%d_Ev%d_Spectra.png",plotPath,year_now,month_now,day_now,runNum,event);
	TCanvas *cSpec = new TCanvas("","",4*1100,4*850);
	cSpec->Divide(4,4);
	for(int i=0; i<16; i++){
		cSpec->cd(i+1);
		grWaveformsPowerSpectrum[i]->Draw("AL");
		grWaveformsPowerSpectrum[i]->SetLineWidth(3);
		gPad->SetLogy();
	}
	cSpec->SaveAs(save_temp_title);
	delete cSpec;
	for(int i=0; i<16; i++){
		delete waveforms[i];
		delete grWaveformsInt[i];
		delete grWaveformsPadded[i];
		delete grWaveformsPowerSpectrum[i];
	}
	delete realAtriEvPtr;
	mapFile->Close();
	delete mapFile;
	return 0;
}
