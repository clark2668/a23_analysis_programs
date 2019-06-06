////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	v2_analysis_reco.cxx 
////	A23 diffuse, do reconstruction
////
////	Nov 2018
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <iomanip>
#include <sstream>

//AraRoot Includes
#include "RawAtriStationEvent.h"
#include "UsefulAraStationEvent.h"
#include "UsefulAtriStationEvent.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
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
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;
	
	if(argc<3) {
		std::cout << "Usage\n" << argv[0] << " <station> <ValsForCut_Filename> \n";
		return -1;
	}
	int station=atoi(argv[1]);
	
	TTimeStamp start;
	TTimeStamp stop;
	start.Set(2013,1,01,00,00,00,0,true,0);
	stop.Set(2016,12,31,24,00,00,0,true,0);
	int start_bin = start.GetSec();
	int stop_bin = stop.GetSec();
	int numBins=(stop_bin-start_bin)/60/60;

	TH2D *phi_vs_time[2];
	TH2D *theta_vs_time[2];
	phi_vs_time[0] = new TH2D("","V Phi Vs Time",numBins, start_bin, stop_bin, 360,-180,180);
	theta_vs_time[0] = new TH2D("","V Theta Vs Time",numBins, start_bin, stop_bin, 180,-90,90);
	phi_vs_time[1] = new TH2D("","H Phi Vs Time",numBins, start_bin, stop_bin, 360,-180,180);
	theta_vs_time[1] = new TH2D("","H Theta Vs Time",numBins, start_bin, stop_bin, 180,-90,90);
	for(int i=0; i<2; i++){
		phi_vs_time[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		phi_vs_time[i]->GetXaxis()->SetTimeFormat("%y/%m");
		phi_vs_time[i]->GetXaxis()->SetTimeOffset(0.,"GMT");
		phi_vs_time[i]->GetXaxis()->SetNdivisions(4,12,0,false);
		theta_vs_time[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		theta_vs_time[i]->GetXaxis()->SetTimeFormat("%y/%m");
		theta_vs_time[i]->GetXaxis()->SetTimeOffset(0.,"GMT");
		theta_vs_time[i]->GetXaxis()->SetNdivisions(4,12,0,false);

	}
	TH2D *theta_vs_phi[2];
	theta_vs_phi[0] = new TH2D("","300m V Phi vs Theta ", 360,-180,180,180,-90,90);
	theta_vs_phi[1] = new TH2D("","300m H Phi vs Theta ", 360,-180,180,180,-90,90);
	int numUsed=0;

	for(int file_num=2; file_num<argc; file_num++){

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

		//now to loop over events
		for(int event=0; event<numEntries; event++){
			trees[0]->GetEvent(event);
			trees[1]->GetEvent(event);
			trees[2]->GetEvent(event);

			if(isBadEvent){
				continue;
			}

			if(!isCal && !isSoft && !isShort && !isNewBox){
				numUsed++;
				for(int pol=0; pol<2; pol++){
					phi_vs_time[pol]->Fill(unixTime, phi_300[pol]);
					theta_vs_time[pol]->Fill(unixTime, theta_300[pol]);
					theta_vs_phi[pol]->Fill(phi_300[pol],theta_300[pol]);
				}	
			}
		}

		// close these after for dumb root permissions reasons
		inputFile->Close();
		delete inputFile;
	}



	gStyle->SetOptStat(0);
	gStyle->SetPalette(55,0);
	
	char title[300];
	TCanvas *cRecoVsTime = new TCanvas("","",2*850,2*850);
	cRecoVsTime->Divide(2,2);
	cRecoVsTime->cd(1);
		gPad->SetLogz();
		phi_vs_time[0]->Draw("colz");
		phi_vs_time[0]->GetXaxis()->SetTitle("Unixtime (Year/Month)");
		phi_vs_time[0]->GetYaxis()->SetTitle("VPol Phi");
		phi_vs_time[0]->GetZaxis()->SetTitle("Number of Events");
		gPad->SetLogz();
		gPad->SetRightMargin(0.12);
	cRecoVsTime->cd(2);
		gPad->SetLogz();
		phi_vs_time[1]->Draw("colz");
		phi_vs_time[1]->GetXaxis()->SetTitle("Unixtime (Year/Month)");
		phi_vs_time[1]->GetYaxis()->SetTitle("Hpol Phi");
		phi_vs_time[1]->GetZaxis()->SetTitle("Number of Events");
		gPad->SetLogz();
		gPad->SetRightMargin(0.12);
	cRecoVsTime->cd(3);
		gPad->SetLogz();
		theta_vs_time[0]->Draw("colz");
		theta_vs_time[0]->GetXaxis()->SetTitle("Unixtime (Year/Month)");
		theta_vs_time[0]->GetYaxis()->SetTitle("VPol Theta");
		theta_vs_time[0]->GetZaxis()->SetTitle("Number of Events");
		gPad->SetLogz();
		gPad->SetRightMargin(0.12);
	cRecoVsTime->cd(4);
		gPad->SetLogz();
		theta_vs_time[1]->Draw("colz");
		theta_vs_time[1]->GetXaxis()->SetTitle("Unixtime (Year/Month)");
		theta_vs_time[1]->GetYaxis()->SetTitle("Hpol Theta");
		theta_vs_time[1]->GetZaxis()->SetTitle("Number of Events");
		gPad->SetLogz();
		gPad->SetRightMargin(0.12);
	sprintf(title, "%s/maps_stacking/%d.%d.%d_A%d_%dEvents_RecoVsTime.png",plotPath,year_now, month_now, day_now,station,numUsed);
	cRecoVsTime->SaveAs(title);
	delete cRecoVsTime;
	delete phi_vs_time[0]; delete phi_vs_time[1];	delete theta_vs_time[0]; delete theta_vs_time[1];

	TCanvas *cThetaVsPhi = new TCanvas("","",2*850,850);
	cThetaVsPhi->Divide(2,1);
	cThetaVsPhi->cd(1);
		theta_vs_phi[0]->Draw("colz");
		theta_vs_phi[0]->GetXaxis()->SetTitle("Phi (deg)");
		theta_vs_phi[0]->GetYaxis()->SetTitle("Theta (deg)");
		theta_vs_phi[0]->GetZaxis()->SetTitle("Number of Events");
		gPad->SetLogz();
		gPad->SetRightMargin(0.12);
	cThetaVsPhi->cd(2);
		theta_vs_phi[1]->Draw("colz");
		theta_vs_phi[1]->GetXaxis()->SetTitle("Phi (deg)");
		theta_vs_phi[1]->GetYaxis()->SetTitle("Theta (deg)");
		theta_vs_phi[1]->GetZaxis()->SetTitle("Number of Events");
		gPad->SetLogz();
		gPad->SetRightMargin(0.12);
	sprintf(title, "%s/maps_stacking/%d.%d.%d_A%d_%dEvents_ThetaVsPhi.png",plotPath,year_now, month_now, day_now,station,numUsed);	
	cThetaVsPhi->SaveAs(title);
	delete cThetaVsPhi;

	TH1D *phi_projection[2];
	TH1D *theta_projection[2];
	for(int pol=0; pol<2; pol++){
		phi_projection[pol] = (TH1D*) theta_vs_phi[pol]->ProjectionX()->Clone();
		theta_projection[pol] = (TH1D*) theta_vs_phi[pol]->ProjectionY()->Clone();
	}

	TCanvas *cProjection = new TCanvas("","",2*850,2*850);
	cProjection->Divide(2,2);
	cProjection->cd(1);
		phi_projection[0]->Draw("");
		phi_projection[0]->SetTitle("V Phi Peak Locations");
		phi_projection[0]->GetXaxis()->SetTitle("Phi (deg)");
		phi_projection[0]->GetYaxis()->SetTitle("Number of Events");
		phi_projection[0]->GetYaxis()->SetTitleOffset(1.2);
		phi_projection[0]->SetLineWidth(3);
		gPad->SetLogy();
		gPad->SetRightMargin(0.12);
	cProjection->cd(2);
		phi_projection[1]->Draw("");
		phi_projection[1]->SetTitle("H Phi Peak Locations");
		phi_projection[1]->GetYaxis()->SetTitleOffset(1.2);
		phi_projection[1]->GetXaxis()->SetTitle("Phi (deg)");
		phi_projection[1]->GetYaxis()->SetTitle("Number of Events");
		phi_projection[1]->GetYaxis()->SetTitleOffset(1.2);
		phi_projection[1]->SetLineWidth(3);
		gPad->SetLogy();
		gPad->SetRightMargin(0.12);
	cProjection->cd(3);
		theta_projection[0]->Draw("");
		theta_projection[0]->SetTitle("V Theta Peak Locations");
		theta_projection[0]->GetXaxis()->SetTitle("Theta (deg)");
		theta_projection[0]->GetYaxis()->SetTitle("Number of Events");
		theta_projection[0]->GetYaxis()->SetTitleOffset(1.2);
		theta_projection[0]->SetLineWidth(3);
		gPad->SetLogy();
		gPad->SetRightMargin(0.12);
	cProjection->cd(4);
		theta_projection[1]->Draw("");
		theta_projection[1]->SetTitle("H Theta Peak Locations");
		theta_projection[1]->GetXaxis()->SetTitle("Theta (deg)");
		theta_projection[1]->GetYaxis()->SetTitle("Number of Events");
		theta_projection[1]->GetYaxis()->SetTitleOffset(1.2);
		theta_projection[1]->SetLineWidth(3);
		gPad->SetLogy();
		gPad->SetRightMargin(0.12);
	sprintf(title, "%s/maps_stacking/%d.%d.%d_A%d_%dEvents_ThetaPhiProjections.png",plotPath,year_now, month_now, day_now,station,numUsed);
	cProjection->SaveAs(title);
	delete cProjection;
	// delete phi_projection[0]; delete phi_projection[1];	delete theta_projection[0]; delete theta_projection[1];
	// delete theta_vs_phi[0]; delete theta_vs_phi[1];

}