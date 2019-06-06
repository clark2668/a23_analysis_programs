////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	cuts_refine_cal_box
////	Reads recco_save files and refines the cal pule cut box
////
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
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLine.h"
#include "TTimeStamp.h"

//AraRoot includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"

//custom analysis includes
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_Cuts.h"

char *DataDirPath(getenv("DATA_DIR"));
char *PedDirPath(getenv("PED_DIR"));
char *plotPath(getenv("PLOT_PATH"));

using namespace std;

int main(int argc, char **argv)
{
	gStyle->SetOptStat(1111);
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <station> <config> <reco_val_save_file_1> <reco_val_save_file_2 > ... <reco_val_save_file_x>"<<endl;
		return 0;
	}
	int station = atoi(argv[1]);
	int config = atoi(argv[2]);

	TTimeStamp start;
	TTimeStamp stop;
	start.Set(2013, 01, 00, 00, 00,0,0,true,0);
	stop.Set(2016, 12, 31, 24, 00,0,0,true,0);

	int start_bin = start.GetSec();
	int stop_bin = stop.GetSec();

	TH2D *phi_vs_time[2];
	TH2D *theta_vs_time[2];
	int num_bins=(stop_bin-start_bin)/60/60;
	phi_vs_time[0] = new TH2D("","V Phi Vs Time",num_bins, start_bin, stop_bin, 360,-180,180);
	theta_vs_time[0] = new TH2D("","V Theta Vs Time",num_bins, start_bin, stop_bin, 180,-90,90);
	phi_vs_time[1] = new TH2D("","H Phi Vs Time",num_bins, start_bin, stop_bin, 360,-180,180);
	theta_vs_time[1] = new TH2D("","H Theta Vs Time",num_bins, start_bin, stop_bin, 180,-90,90);
	for(int i=0; i<2; i++){
		phi_vs_time[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		phi_vs_time[i]->GetXaxis()->SetTimeFormat("%y/%m");
		theta_vs_time[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		theta_vs_time[i]->GetXaxis()->SetTimeFormat("%y/%m");
		phi_vs_time[i]->GetXaxis()->SetTimeOffset(0.,"GMT");
		theta_vs_time[i]->GetXaxis()->SetTimeOffset(0.,"GMT");
	}

	stringstream ss;
	int num_total=0;

	for(int file_num=3; file_num<argc; file_num++){

		cout << "Run " << file_num << " :: " << argv[file_num] << endl;

		string chRun = "run";
		string file = string(argv[file_num]);
		size_t foundRun = file.find(chRun);
		string strRunNum = file.substr(foundRun+4,4);
		int runNum = atoi(strRunNum.c_str());

		// if(isBadRun(station,runNum)) continue;

		TFile *inputFile = TFile::Open(argv[file_num]);
		if(!inputFile){
			cout<<"Can't open joined file!"<<endl;
			return -1;
		}
		TTree *trees = (TTree*) inputFile->Get("RecoVals");
		if(!trees){
			cout<<"Can't open reco tree"<<endl;
			return -1;
		}
		double phi_41_V;
		double theta_41_V;
		double phi_41_H;
		double theta_41_H;
		double corr_val_V;
		double corr_val_H;
		int isCal;
		int isBad;
		int unixTime;
		trees->SetBranchAddress("phi_41_V",&phi_41_V);
		trees->SetBranchAddress("theta_41_V",&theta_41_V);
		trees->SetBranchAddress("phi_41_H",&phi_41_H);
		trees->SetBranchAddress("theta_41_H",&theta_41_H);
		trees->SetBranchAddress("corr_val_V",&corr_val_V);
		trees->SetBranchAddress("corr_val_H",&corr_val_H);
		trees->SetBranchAddress("cal",&isCal);
		trees->SetBranchAddress("isBad",&isBad);
		trees->SetBranchAddress("unixTime",&unixTime);

		int numEntries = trees->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;

		//now to loop over events
		for(int event=0; event<numEntries; event++){

			trees->GetEvent(event);

			if(!isCal || isBad) continue;

			num_total++;
			if(corr_val_V>corr_val_H){

				phi_vs_time[0]->Fill(unixTime, phi_41_V);
				theta_vs_time[0]->Fill(unixTime, theta_41_V);
			}
			else if(corr_val_H>corr_val_V){
				phi_vs_time[1]->Fill(unixTime, phi_41_H);
				theta_vs_time[1]->Fill(unixTime, theta_41_H);
			}
		}//loop over events
		inputFile->Close();
		delete inputFile;
	} //end loop over input files

	char title[300];
	gStyle->SetOptStat(0);
	TCanvas *c_phi_vs_time = new TCanvas("","",2*1100,2*850);
	c_phi_vs_time->Divide(2,2);
	c_phi_vs_time->cd(1);
		phi_vs_time[0]->Draw("colz");
		phi_vs_time[0]->GetXaxis()->SetTitle("Unixtime");
		phi_vs_time[0]->GetYaxis()->SetTitle("Phi");
	c_phi_vs_time->cd(2);
		phi_vs_time[1]->Draw("colz");
		phi_vs_time[1]->GetXaxis()->SetTitle("Unixtime");
		phi_vs_time[1]->GetYaxis()->SetTitle("Phi");
	c_phi_vs_time->cd(3);
		theta_vs_time[0]->Draw("colz");
		theta_vs_time[0]->GetXaxis()->SetTitle("Unixtime");
		theta_vs_time[0]->GetYaxis()->SetTitle("Theta");		
	c_phi_vs_time->cd(4);
		theta_vs_time[1]->Draw("colz");
		theta_vs_time[1]->GetXaxis()->SetTitle("Unixtime");
		theta_vs_time[1]->GetYaxis()->SetTitle("Theta");
	sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_RecoVsTime.png",plotPath,year_now, month_now, day_now,station,config,num_total);
	c_phi_vs_time->SaveAs(title); delete c_phi_vs_time;
	delete phi_vs_time[0]; delete phi_vs_time[1];	

}