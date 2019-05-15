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

using namespace std;

int getThisRunNum(char *thearg){
	string chRun = "run";
	string file = string(thearg);
	size_t foundRun = file.find(chRun);
	string strRunNum = file.substr(foundRun+4,4);
	int runNum = atoi(strRunNum.c_str());
	return runNum;
}

int main(int argc, char **argv)
{

	if(argc<4){
		cout<< "Usage\n" << argv[0] << " <station> <output_location> <ValForCuts filename>"<<endl;;
		return -1;
	}
	int station = atoi(argv[1]);
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	// char *plotPath(getenv("PLOT_PATH"));
	// if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

	gStyle->SetOptStat(0);

	int startRun = getThisRunNum(argv[3]);
	int stopRun = getThisRunNum(argv[argc-1]);

	char title[300];
	// TCanvas *c_start = new TCanvas("start","start",2*850,2*850);
	// c_start->Divide(2,2);
	// sprintf(title, "%s/%d.%d.%d_A%d_Run%d_to_Run%d_RecoVsTime.pdf(",argv[2],year_now, month_now, day_now,station,startRun,stopRun);
	// c_start->Print(title,"pdf");

	TH1D *surface_events_per_minute[2];
	surface_events_per_minute[0]= new TH1D("","VPol Surface Events Per Minute",200,0,200);
	surface_events_per_minute[1]= new TH1D("","Hpol Surface Events Per Minute",200,0,200);

	for(int file_num=3; file_num<argc; file_num++){

		int runNum=getThisRunNum(argv[file_num]);
		
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

		int theta_300[2];
		int phi_300[2];

		trees[0]->SetBranchAddress("theta_300_V",&theta_300[0]);
		trees[0]->SetBranchAddress("phi_300_V",&phi_300[0]);
		
		trees[1]->SetBranchAddress("theta_300_H",&theta_300[1]);
		trees[1]->SetBranchAddress("phi_300_H",&phi_300[1]);

		int isCal;
		int isSoft;
		int isShort;
		int isNewBox;
		int isBadEvent;
		int unixTime;

		trees[2]->SetBranchAddress("cal",&isCal);
		trees[2]->SetBranchAddress("soft",&isSoft);
		trees[2]->SetBranchAddress("short",&isShort);
		trees[2]->SetBranchAddress("box",&isNewBox);
		trees[2]->SetBranchAddress("bad",&isBadEvent);
		trees[2]->SetBranchAddress("unixTime",&unixTime);

		int numEntries = trees[0]->GetEntries();

		trees[2]->GetEvent(0);
		int unixStart = unixTime;
		trees[2]->GetEvent(numEntries-1);
		int unixStop = unixTime;

		TTimeStamp start(unixStart);
		TTimeStamp stop(unixStop);
		int start_bin = start.GetSec();
		int stop_bin = stop.GetSec();
		// gStyle->SetTimeOffset(0.,"GMT");

		TH2D *phi_vs_event[2];
		TH2D *theta_vs_event[2];
		int num_bins=int((stop_bin-start_bin)/60.);
		phi_vs_event[0] = new TH2D("","V Phi Vs Time",num_bins, start_bin, stop_bin, 360,-180,180);
		theta_vs_event[0] = new TH2D("","V Theta Vs Time",num_bins, start_bin, stop_bin, 180,-90,90);
		phi_vs_event[1] = new TH2D("","H Phi Vs Time",num_bins, start_bin, stop_bin, 360,-180,180);
		theta_vs_event[1] = new TH2D("","H Theta Vs Time",num_bins, start_bin, stop_bin, 180,-90,90);
		for(int i=0; i<2; i++){
			phi_vs_event[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
			phi_vs_event[i]->GetXaxis()->SetTimeFormat("%H:%M");
			phi_vs_event[i]->GetXaxis()->SetTimeOffset(0.,"GMT");
			theta_vs_event[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
			theta_vs_event[i]->GetXaxis()->SetTimeFormat("%H:%M");
			theta_vs_event[i]->GetXaxis()->SetTimeOffset(0.,"GMT");
		}

		//now to loop over events
		for(int event=0; event<numEntries; event++){
			trees[0]->GetEvent(event);
			trees[1]->GetEvent(event);
			trees[2]->GetEvent(event);

			if(isBadEvent) continue;

			if(!isCal && !isSoft && !isShort && !isNewBox){
				for(int pol=0; pol<2; pol++){
					phi_vs_event[pol]->Fill(unixTime, phi_300[pol]);
					theta_vs_event[pol]->Fill(unixTime, theta_300[pol]);
				}	
			}
		}

		gStyle->SetOptStat(0);
		gStyle->SetPalette(55,0);
		TCanvas *c_temp = new TCanvas("","",2*850,3*850);
		c_temp->Divide(2,3);
		c_temp->cd(1);
			phi_vs_event[0]->Draw("colz");
			phi_vs_event[0]->GetXaxis()->SetTitle("Unixtime");
			phi_vs_event[0]->GetYaxis()->SetTitle("VPol Phi");
			char run_info[150];
			sprintf(run_info,"run %d",runNum);
			phi_vs_event[0]->SetTitle(run_info);
			// phi_vs_event[0]->GetZaxis()->SetRangeUser(1,6);
		c_temp->cd(2);
			phi_vs_event[1]->Draw("colz");
			phi_vs_event[1]->GetXaxis()->SetTitle("Unixtime");
			phi_vs_event[1]->GetYaxis()->SetTitle("Hpol Phi");
			// phi_vs_event[1]->GetZaxis()->SetRangeUser(1,6);
		c_temp->cd(3);
			theta_vs_event[0]->Draw("colz");
			theta_vs_event[0]->GetXaxis()->SetTitle("Unixtime");
			theta_vs_event[0]->GetYaxis()->SetTitle("VPol Theta");
			char start_string[150];
			sprintf(start_string,"start: %s",start.AsString("s"));
			theta_vs_event[0]->SetTitle(start_string);
			// theta_vs_event[0]->GetZaxis()->SetRangeUser(1,6);
		c_temp->cd(4);
			theta_vs_event[1]->Draw("colz");
			theta_vs_event[1]->GetXaxis()->SetTitle("Unixtime");
			theta_vs_event[1]->GetYaxis()->SetTitle("Hpol Theta");
			char stop_string[150];
			sprintf(stop_string,"stop: %s",stop.AsString("s"));
			theta_vs_event[1]->SetTitle(stop_string);

		TH1D *projectTheta[2];
		for(int pol=0; pol<2; pol++){
			projectTheta[pol] = (TH1D*) theta_vs_event[pol]->ProjectionX("",127,180);
			for(int bin=0; bin<projectTheta[pol]->GetNbinsX(); bin++){
				surface_events_per_minute[pol]->Fill(projectTheta[pol]->GetBinContent(bin));
			}
		}
		// c_temp->cd(5);
		// 	projectTheta[0]->Draw();
		// 	projectTheta[0]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		// 	projectTheta[0]->GetXaxis()->SetTimeFormat("%H:%M");
		// 	projectTheta[0]->GetXaxis()->SetTimeOffset(0.,"GMT");
		// c_temp->cd(6);
		// 	projectTheta[1]->Draw();
		// 	projectTheta[1]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		// 	projectTheta[1]->GetXaxis()->SetTimeFormat("%H:%M");
		// 	projectTheta[1]->GetXaxis()->SetTimeOffset(0.,"GMT");
		// sprintf(title, "%s/%d.%d.%d_A%d_Run%d_to_Run%d_RecoVsTime.png",argv[2],year_now, month_now, day_now,station,startRun,stopRun);
		sprintf(title, "%s/%d.%d.%d_A%d_Run%d_RecoVsTime.png",argv[2],year_now, month_now, day_now,station,runNum);			
		// c_temp->Print(title,"png");
		delete c_temp;
		delete phi_vs_event[0]; delete phi_vs_event[1];	delete theta_vs_event[0]; delete theta_vs_event[1];
		delete projectTheta[0]; delete projectTheta[1];

		// close these after for dumb root permissions reasons
		inputFile->Close();
		delete inputFile;
	}

	TCanvas *c_surfevpermin = new TCanvas("","",2*850,850);
	c_surfevpermin->Divide(2,1);
	c_surfevpermin->cd(1);
		surface_events_per_minute[0]->Draw();
		surface_events_per_minute[0]->GetXaxis()->SetTitle("Surface Events per minute");
		surface_events_per_minute[0]->GetYaxis()->SetTitle("Number of Minutes");
		gPad->SetLogy();
		surface_events_per_minute[0]->SetTitleOffset(1.1);
	c_surfevpermin->cd(2);
		surface_events_per_minute[1]->Draw();
		surface_events_per_minute[1]->GetXaxis()->SetTitle("Surface Events per minute");
		surface_events_per_minute[1]->GetYaxis()->SetTitle("Number of Minutes");
		gPad->SetLogy();
	sprintf(title, "%s/%d.%d.%d_A%d_SurfaceEventsPerMin_Run%d_to_Run%d.png",argv[2],year_now, month_now, day_now,station,startRun,stopRun);			
	c_surfevpermin->Print(title,"png");	

	// TCanvas *c_stop= new TCanvas("stop","stop",2*850,2*850);
	// c_stop->Divide(2,2);
	// sprintf(title, "%s/%d.%d.%d_A%d_Run%d_to_Run%d_RecoVsTime.pdf)",argv[2],year_now, month_now, day_now,station,startRun,stopRun);
	// c_stop->Print(title,"pdf");

}