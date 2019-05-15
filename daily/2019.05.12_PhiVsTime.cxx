////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	2019.05.12_PhiVsTime.cxx 
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
#include "TLine.h"
#include "TTimeStamp.h"

//AraRoot includes
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_Cuts.h"
#include "tools_CW.h"

using namespace std;

int main(int argc, char **argv)
{

	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <station> <ValForCuts filename>"<<endl;;
		return -1;
	}
	int station = atoi(argv[1]);
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

	stringstream ss;

	gStyle->SetOptStat(0);

	TTimeStamp start;
	TTimeStamp stop;
	start.Set(2014,01,10,19,15,0,0,true,0);
	stop.Set(2014,01,10,20,15,0,0,true,0);

	// run 4775
	// start.Set(2015, 01, 03, 15, 30,0,0,true,0);
	// stop.Set(2015, 01, 03, 16, 00,0,0,true,0);

	// run 6705
	// start.Set(2016, 01, 12, 00, 00,0,0,true,0);
	// stop.Set(2016, 01, 12, 24, 00,0,0,true,0);

	int start_bin = start.GetSec();
	int stop_bin = stop.GetSec();

	TH2D *phi_vs_event[2];
	TH2D *theta_vs_event[2];
	phi_vs_event[0] = new TH2D("","V Phi Vs Time",stop_bin-start_bin, start_bin, stop_bin, 360,-180,180);
	theta_vs_event[0] = new TH2D("","V Theta Vs Time",stop_bin-start_bin, start_bin, stop_bin, 180,-90,90);
	phi_vs_event[1] = new TH2D("","H Phi Vs Time",stop_bin-start_bin, start_bin, stop_bin, 360,-180,180);
	theta_vs_event[1] = new TH2D("","H Theta Vs Time",stop_bin-start_bin, start_bin, stop_bin, 180,-90,90);
	for(int i=0; i<2; i++){
		phi_vs_event[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		phi_vs_event[i]->GetXaxis()->SetTimeFormat("%H:%M");
		theta_vs_event[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		theta_vs_event[i]->GetXaxis()->SetTimeFormat("%H:%M");
		phi_vs_event[i]->GetXaxis()->SetTimeOffset(0.,"GMT");
		theta_vs_event[i]->GetXaxis()->SetTimeOffset(0.,"GMT");
	}

	vector<vector<double> > unixTimes;
	vector<vector<double> > phis;
	vector<vector<double> > thetas;
	unixTimes.resize(2);
	phis.resize(2);
	thetas.resize(2);
	vector<vector<double> > all_unixTimes;
	vector<vector<double> > all_phis;
	vector<vector<double> > all_thetas;
	all_unixTimes.resize(2);
	all_phis.resize(2);
	all_thetas.resize(2);	

	for(int file_num=2; file_num<argc; file_num++){

		string chRun = "run";
		string file = string(argv[file_num]);
		size_t foundRun = file.find(chRun);
		string strRunNum = file.substr(foundRun+4,4);
		int runNum = atoi(strRunNum.c_str());
		int isThisBadABadRun = isBadRun(station,runNum);

		TFile *inputFile = TFile::Open(argv[file_num]);
		if(!inputFile){
			cout<<"Can't open joined file!"<<endl;
			return -1;
		}
		// cout << "Run " << file_num << " :: " << argv[file_num] << endl;

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
		int theta_300[2];
		int phi_300[2];
		int theta_41[2];
		int phi_41[2];

		trees[0]->SetBranchAddress("corr_val_V",&corr_val[0]);
		trees[0]->SetBranchAddress("snr_val_V",&snr_val[0]);
		trees[0]->SetBranchAddress("wfrms_val_V",&WFRMS[0]);
		trees[0]->SetBranchAddress("Refilt_V",&Refilt[0]);
		trees[0]->SetBranchAddress("theta_300_V",&theta_300[0]);
		trees[0]->SetBranchAddress("theta_41_V",&theta_41[0]);
		trees[0]->SetBranchAddress("phi_300_V",&phi_300[0]);
		trees[0]->SetBranchAddress("phi_41_V",&phi_41[0]);
		
		trees[1]->SetBranchAddress("corr_val_H",&corr_val[1]);
		trees[1]->SetBranchAddress("snr_val_H",&snr_val[1]);
		trees[1]->SetBranchAddress("wfrms_val_H",&WFRMS[1]);
		trees[1]->SetBranchAddress("Refilt_H",&Refilt[1]);
		trees[1]->SetBranchAddress("theta_300_H",&theta_300[1]);
		trees[1]->SetBranchAddress("theta_41_H",&theta_41[1]);
		trees[1]->SetBranchAddress("phi_300_H",&phi_300[1]);
		trees[1]->SetBranchAddress("phi_41_H",&phi_41[1]);

		int isCal;
		int isSoft;
		int isShort;
		int isCW;
		int isNewBox;
		int isSurf[2];
		int isBadEvent;
		int isSurfEvent_top[2];
		int unixTime;

		trees[2]->SetBranchAddress("cal",&isCal);
		trees[2]->SetBranchAddress("soft",&isSoft);
		trees[2]->SetBranchAddress("short",&isShort);
		trees[2]->SetBranchAddress("CW",&isCW);
		trees[2]->SetBranchAddress("box",&isNewBox);
		trees[2]->SetBranchAddress("surf_V",&isSurf[0]);
		trees[2]->SetBranchAddress("surf_H",&isSurf[1]);
		trees[2]->SetBranchAddress("bad",&isBadEvent);
		trees[2]->SetBranchAddress("surf_top_V",&isSurfEvent_top[0]);
		trees[2]->SetBranchAddress("surf_top_H",&isSurfEvent_top[1]);
		trees[2]->SetBranchAddress("unixTime",&unixTime);

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
		int start=0;
		// numEntries=10000;
		for(int event=start; event<numEntries; event++){

			trees[0]->GetEvent(event);
			trees[1]->GetEvent(event);
			trees[2]->GetEvent(event);

			if(isBadEvent){
				continue;
			}

			if(!isCal && !isSoft && !isShort && !isNewBox){
				for(int pol=0; pol<2; pol++){
					phi_vs_event[pol]->Fill(unixTime, phi_300[pol]);
					theta_vs_event[pol]->Fill(unixTime, theta_300[pol]);
					all_unixTimes[pol].push_back(unixTime);
					all_phis[pol].push_back(phi_300[pol]);
					all_thetas[pol].push_back(theta_300[pol]);

					if(runNum==4775 && pol==0 && unixTime==1420317755 && event==9520){
						unixTimes[0].push_back(double(unixTime)); phis[0].push_back(double(phi_300[0])); thetas[0].push_back(double(theta_300[0]));
					}
					if(runNum==4775 && pol==1 && unixTime==1420317755 && event==9520){
						unixTimes[1].push_back(double(unixTime)); phis[1].push_back(double(phi_300[1])); thetas[1].push_back(double(theta_300[1]));
					}
					if(runNum==4775 && pol==0 && unixTime==1420317765 && event==9533){
						unixTimes[0].push_back(double(unixTime)); phis[0].push_back(double(phi_300[0])); thetas[0].push_back(double(theta_300[0]));
					}
					if(runNum==4775 && pol==1 && unixTime==1420317765 && event==9533){
						unixTimes[1].push_back(double(unixTime)); phis[1].push_back(double(phi_300[1])); thetas[1].push_back(double(theta_300[1]));
					}
					if(runNum==4775 && pol==0 && unixTime==1420317765 && event==9536){
						unixTimes[0].push_back(double(unixTime)); phis[0].push_back(double(phi_300[0])); thetas[0].push_back(double(theta_300[0]));
					}
					if(runNum==4775 && pol==1 && unixTime==1420317765 && event==9536){
						unixTimes[1].push_back(double(unixTime)); phis[1].push_back(double(phi_300[1])); thetas[1].push_back(double(theta_300[1]));
					}
					if(runNum==4775 && pol==1 && unixTime==1420317816 && event==9570){
						unixTimes[1].push_back(double(unixTime)); phis[1].push_back(double(phi_300[1])); thetas[1].push_back(double(theta_300[1]));
					}
					if(runNum==4775 && pol==0 && unixTime==1420317831 && event==9596){
						unixTimes[0].push_back(double(unixTime)); phis[0].push_back(double(phi_300[0])); thetas[0].push_back(double(theta_300[0]));
					}
					if(runNum==4775 && pol==0 && unixTime==1420317832 && event==9598){
						unixTimes[0].push_back(double(unixTime)); phis[0].push_back(double(phi_300[0])); thetas[0].push_back(double(theta_300[0]));
					}
					if(runNum==4775 && pol==0 && unixTime==1420317835 && event==9604){
						unixTimes[0].push_back(double(unixTime)); phis[0].push_back(double(phi_300[0])); thetas[0].push_back(double(theta_300[0]));
					}
					if(runNum==4775 && pol==0 && unixTime==1420317859 && event==9615){
						unixTimes[0].push_back(double(unixTime)); phis[0].push_back(double(phi_300[0])); thetas[0].push_back(double(theta_300[0]));
					}

					if(runNum==2868 && pol==0 && unixTime==1389382075 && event==928){
						unixTimes[0].push_back(double(unixTime)); phis[0].push_back(double(phi_300[0])); thetas[0].push_back(double(theta_300[0]));
					}
					if(runNum==2868 && pol==0 && unixTime==1389383363 && event==1774){
						unixTimes[0].push_back(double(unixTime)); phis[0].push_back(double(phi_300[0])); thetas[0].push_back(double(theta_300[0]));
					}
					if(runNum==2868 && pol==0 && unixTime==1389383366 && event==1777){
						unixTimes[0].push_back(double(unixTime)); phis[0].push_back(double(phi_300[0])); thetas[0].push_back(double(theta_300[0]));
					}

					if(WFRMS[pol])
						continue;
				}	
			}
		}
		inputFile->Close();
		delete inputFile;

		TGraph *gr_outlier_phi_V = new TGraph(unixTimes[0].size(),&unixTimes[0][0],&phis[0][0]);
		gr_outlier_phi_V->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		TGraph *gr_outlier_thetas_V = new TGraph(unixTimes[0].size(),&unixTimes[0][0],&thetas[0][0]);
		gr_outlier_thetas_V->GetXaxis()->SetTimeDisplay(1); //turn on a time axis

		TGraph *gr_outlier_phi_H = new TGraph(unixTimes[1].size(),&unixTimes[1][0],&phis[1][0]);
		gr_outlier_phi_H->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		TGraph *gr_outlier_thetas_H = new TGraph(unixTimes[1].size(),&unixTimes[1][0],&thetas[1][0]);
		gr_outlier_thetas_H->GetXaxis()->SetTimeDisplay(1); //turn on a time axis

		vector<TGraph*> g_phi_vs_event;
		vector<TGraph*> g_theta_vs_event;
		for(int pol=0; pol<2; pol++){
			g_phi_vs_event.push_back(new TGraph(all_unixTimes[pol].size(), &all_unixTimes[pol][0], &all_phis[pol][0]));
			g_theta_vs_event.push_back(new TGraph(all_unixTimes[pol].size(), &all_unixTimes[pol][0], &all_thetas[pol][0]));
			
			g_phi_vs_event[pol]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
			g_theta_vs_event[pol]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
			g_phi_vs_event[pol]->GetXaxis()->SetTimeFormat("%H:%M");
			g_theta_vs_event[pol]->GetXaxis()->SetTimeFormat("%H:%M");
			g_phi_vs_event[pol]->GetXaxis()->SetTimeOffset(0.,"GMT");
			g_theta_vs_event[pol]->GetXaxis()->SetTimeOffset(0.,"GMT");
		}

		char title[300];
		gStyle->SetOptStat(111111);
		TCanvas *c_phi_vs_event = new TCanvas("","",2*1100,2*850);
		c_phi_vs_event->Divide(2,2);
		c_phi_vs_event->cd(1);
			g_phi_vs_event[0]->Draw("AP");
			g_phi_vs_event[0]->GetXaxis()->SetTitle("Unixtime");
			g_phi_vs_event[0]->GetYaxis()->SetTitle("Phi");
			g_phi_vs_event[0]->SetMarkerSize(0.5);
			g_phi_vs_event[0]->SetMarkerStyle(31);
			g_phi_vs_event[0]->GetXaxis()->SetRangeUser(start_bin, stop_bin);
			gr_outlier_phi_V->Draw("sameP");
			gr_outlier_phi_V->SetMarkerColor(kRed);
			gr_outlier_phi_V->SetMarkerStyle(4);
			gr_outlier_phi_V->SetMarkerSize(2);
		c_phi_vs_event->cd(2);
			g_phi_vs_event[1]->Draw("AP");
			g_phi_vs_event[1]->GetXaxis()->SetTitle("Unixtime");
			g_phi_vs_event[1]->GetYaxis()->SetTitle("Phi");
			g_phi_vs_event[1]->SetMarkerSize(0.5);
			g_phi_vs_event[1]->SetMarkerStyle(31);
			g_phi_vs_event[1]->GetXaxis()->SetRangeUser(start_bin, stop_bin);
			gr_outlier_phi_H->Draw("sameP");
			gr_outlier_phi_H->SetMarkerColor(kRed);
			gr_outlier_phi_H->SetMarkerStyle(4);
			gr_outlier_phi_H->SetMarkerSize(2);
		c_phi_vs_event->cd(3);
			g_theta_vs_event[0]->Draw("AP");
			g_theta_vs_event[0]->GetXaxis()->SetTitle("Unixtime");
			g_theta_vs_event[0]->GetYaxis()->SetTitle("Theta");
			g_theta_vs_event[0]->SetMarkerSize(0.5);
			g_theta_vs_event[0]->SetMarkerStyle(31);
			g_theta_vs_event[0]->GetXaxis()->SetRangeUser(start_bin, stop_bin);
			gr_outlier_thetas_V->Draw("sameP");
			gr_outlier_thetas_V->SetMarkerColor(kRed);
			gr_outlier_thetas_V->SetMarkerStyle(4);
			gr_outlier_thetas_V->SetMarkerSize(2);			
		c_phi_vs_event->cd(4);
			g_theta_vs_event[1]->Draw("AP");
			g_theta_vs_event[1]->GetXaxis()->SetTitle("Unixtime");
			g_theta_vs_event[1]->GetYaxis()->SetTitle("Theta");
			g_theta_vs_event[1]->SetMarkerSize(0.5);
			g_theta_vs_event[1]->SetMarkerStyle(31);
			g_theta_vs_event[1]->GetXaxis()->SetRangeUser(start_bin, stop_bin);
			gr_outlier_thetas_H->Draw("sameP");
			gr_outlier_thetas_H->SetMarkerColor(kRed);
			gr_outlier_thetas_H->SetMarkerStyle(4);
			gr_outlier_thetas_H->SetMarkerSize(2);	
		sprintf(title, "%s/%d.%d.%d_A%d_Run%d_RecoVsTime.png",plotPath,year_now, month_now, day_now,station,runNum);
		c_phi_vs_event->SaveAs(title);
		delete c_phi_vs_event;
		delete g_phi_vs_event[0]; delete g_theta_vs_event[1];

		// char title[300];
		// gStyle->SetOptStat(111111);
		// TCanvas *c_phi_vs_event = new TCanvas("","",2*1100,2*850);
		// c_phi_vs_event->Divide(2,2);
		// c_phi_vs_event->cd(1);
		// 	phi_vs_event[0]->Draw("");
		// 	phi_vs_event[0]->GetXaxis()->SetTitle("Unixtime");
		// 	phi_vs_event[0]->GetYaxis()->SetTitle("Phi");
		// 	gr_outlier_phi_V->Draw("sameP");
		// 	gr_outlier_phi_V->SetMarkerColor(kRed);
		// 	gr_outlier_phi_V->SetMarkerStyle(27);
		// 	gr_outlier_phi_V->SetMarkerSize(2);
		// c_phi_vs_event->cd(2);
		// 	phi_vs_event[1]->Draw("");
		// 	phi_vs_event[1]->GetXaxis()->SetTitle("Unixtime");
		// 	phi_vs_event[1]->GetYaxis()->SetTitle("Phi");
		// 	gr_outlier_phi_H->Draw("sameP");
		// 	gr_outlier_phi_H->SetMarkerColor(kRed);
		// 	gr_outlier_phi_H->SetMarkerStyle(27);
		// 	gr_outlier_phi_H->SetMarkerSize(2);
		// c_phi_vs_event->cd(3);
		// 	theta_vs_event[0]->Draw("");
		// 	theta_vs_event[0]->GetXaxis()->SetTitle("Unixtime");
		// 	theta_vs_event[0]->GetYaxis()->SetTitle("Theta");
		// 	gr_outlier_thetas_V->Draw("sameP");
		// 	gr_outlier_thetas_V->SetMarkerColor(kRed);
		// 	gr_outlier_thetas_V->SetMarkerStyle(27);
		// 	gr_outlier_thetas_V->SetMarkerSize(2);			
		// c_phi_vs_event->cd(4);
		// 	theta_vs_event[1]->Draw("");
		// 	theta_vs_event[1]->GetXaxis()->SetTitle("Unixtime");
		// 	theta_vs_event[1]->GetYaxis()->SetTitle("Theta");
		// 	gr_outlier_thetas_H->Draw("sameP");
		// 	gr_outlier_thetas_H->SetMarkerColor(kRed);
		// 	gr_outlier_thetas_H->SetMarkerStyle(27);
		// 	gr_outlier_thetas_H->SetMarkerSize(2);	
		// sprintf(title, "%s/%d.%d.%d_A%d_Run%d_RecoVsTime.png",plotPath,year_now, month_now, day_now,station,runNum);
		// c_phi_vs_event->SaveAs(title);
		// delete c_phi_vs_event;
		// delete phi_vs_event[0]; delete phi_vs_event[1];	
	}

}