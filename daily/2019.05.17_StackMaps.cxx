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

RawAtriStationEvent *rawAtriEvPtr;
UsefulAtriStationEvent *realAtriEvPtr;

#include "Settings.h"
#include "Event.h"
#include "Detector.h"
#include "Report.h"

#include "AraAntennaInfo.h"
#include "AraQualCuts.h"
#include "RayTraceCorrelator.h"

#include "tools_inputParameters.h"
#include "tools_outputObjects.h"
#include "tools_runSummaryObjects.h"
#include "tools_WaveformFns.h"
#include "tools_PlottingFns.h"
#include "tools_Constants.h"
#include "tools_RecoFns.h"
#include "tools_Cuts.h"

AraAntPol::AraAntPol_t Vpol = AraAntPol::kVertical;
AraAntPol::AraAntPol_t Hpol = AraAntPol::kHorizontal;

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
		std::cout << "Usage\n" << argv[0] << " <station> <input reco file> \n";
		return -1;
	}
	int station=atoi(argv[1]);

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
		TTree *OutputTree;
		OutputTree = (TTree*) inputFile->Get("OutputTree");
		vector<TH2D*> maps;
		maps.resize(2);
		int unixTime;
		int isShort;
		OutputTree->SetBranchAddress("VMap",&maps[0]);
		OutputTree->SetBranchAddress("HMap",&maps[1]);
		OutputTree->SetBranchAddress("hasDigitizerError",&hasDigitizerError);
		OutputTree->SetBranchAddress("isCalpulser", &isCalpulser);
		OutputTree->SetBranchAddress("isSoftTrigger", &isSoftTrigger);
		OutputTree->SetBranchAddress("unixTime", &unixTime);
		OutputTree->SetBranchAddress("isShort",&isShort);

		Long64_t numEntries=OutputTree->GetEntries();
		Long64_t starEvery=numEntries/80;
		if(starEvery==0) starEvery++;
		bool found_good_start=false;
		int good_start_index;
		TH2D *stacked_maps[2];
		int numStacked=0;
		// numEntries=50;
		for(int event=0; event<numEntries; event++){
			OutputTree->GetEntry(event);
			if(!hasDigitizerError && !isCalpulser && !isSoftTrigger && !isShort){
				found_good_start=true;
				stacked_maps[0]=(TH2D*) maps[0]->Clone();
				stacked_maps[1]=(TH2D*) maps[1]->Clone();
				numStacked++;
				good_start_index=event;
			}
			if(found_good_start)
				break;
		}
		for(int event=good_start_index+1; event<numEntries; event++){
			OutputTree->GetEntry(event);
			if(!hasDigitizerError && !isCalpulser && !isSoftTrigger && !isShort){
				cout<<"On event "<<event<<endl;
				stacked_maps[0]->Add(maps[0]);
				stacked_maps[1]->Add(maps[1]);
				numStacked++;
			}
		}
		stacked_maps[0]->Scale(1./double(numStacked));
		stacked_maps[1]->Scale(1./double(numStacked));

		beautify_TH2D();
		gStyle->SetOptStat(0);
		TCanvas *c = new TCanvas("","",2.1*850,850);
		c->Divide(2,1);
		c->cd(1);
			stacked_maps[0]->Draw("colz");
			stacked_maps[0]->GetXaxis()->SetTitle("Phi (deg)");
			stacked_maps[0]->GetYaxis()->SetTitle("Theta (deg)");
			stacked_maps[0]->GetZaxis()->SetRangeUser(0.04,0.075);
		c->cd(2);
			stacked_maps[1]->Draw("colz");
			stacked_maps[1]->GetXaxis()->SetTitle("Phi (deg)");
			stacked_maps[1]->GetYaxis()->SetTitle("Theta (deg)");
			stacked_maps[1]->GetZaxis()->SetRangeUser(0.04,0.075);
		char title[300];
		sprintf(title, "%s/maps_stacking/%d.%d.%d_A%d_Run%d_%dStackedMap.png",plotPath,year_now, month_now, day_now,station,runNum,numStacked);
		c->SaveAs(title); 
		delete c;

		cout<<"Done making an averaged map..."<<endl;

		TH2D *before_distro[2];
		TH2D *after_distro[2];
		for(int pol=0; pol<2; pol++){
			before_distro[pol] = new TH2D("","",360,-180,180,180,-90,90);
			after_distro[pol] = new TH2D("","",360,-180,180,180,-90,90);
		}

		bool doIndividualSubtraction=true;
		if(doIndividualSubtraction){
			for(int event=good_start_index+1; event<numEntries; event++){
				cout<<"Fixing event "<<event<<endl;
				OutputTree->GetEntry(event);
				if(!hasDigitizerError && !isSoftTrigger && !isShort){
				// if(!hasDigitizerError && !isCalpulser && !isSoftTrigger && !isShort){
					TH2D *subtracted[2];
					int before_theta[2];
					int before_phi[2];
					int after_theta[2];
					int after_phi[2];
					double before_corr[2];
					double after_corr[2];
					for(int pol=0; pol<2; pol++){
						
						//get original map
						subtracted[pol] = (TH2D*) maps[pol]->Clone();
						
						//compute stuff about the original
						getCorrMapPeak(subtracted[pol], before_theta[pol], before_phi[pol], before_corr[pol]);
						// printf("Event %d before location is %d and %d \n", event, before_theta[pol], before_phi[pol]);
						before_distro[pol]->Fill(before_phi[pol], before_theta[pol]);

						//background subtract
						subtracted[pol]->Add(stacked_maps[pol],-1.);
						getCorrMapPeak(subtracted[pol], after_theta[pol], after_phi[pol], after_corr[pol]);
						// printf("Event %d after location is %d and %d \n", event, after_theta[pol], after_phi[pol]);
						after_distro[pol]->Fill(after_phi[pol], after_theta[pol]);
					}

					bool doThisPrint=false;
					if(doThisPrint){
						stringstream V_original;
							V_original<<" VPol Original Peak (#phi, #theta)=("<<before_phi[0]<<" , "<<before_theta[0]<<")";
							maps[0]->SetTitle(V_original.str().c_str());
						stringstream V_subtracted;
							V_subtracted<<" VPol Subtracted Peak (#phi, #theta)=("<<after_phi[0]<<" , "<<after_theta[0]<<")";
							subtracted[0]->SetTitle(V_subtracted.str().c_str());
						stringstream H_original;
							H_original<<" HPol Original Peak (#phi, #theta)=("<<before_phi[1]<<" , "<<before_theta[1]<<")";
							maps[1]->SetTitle(H_original.str().c_str());
						stringstream H_subtracted;
							H_subtracted<<" HPol Subtracted Peak (#phi, #theta)=("<<after_phi[1]<<" , "<<after_theta[1]<<")";
							subtracted[1]->SetTitle(H_subtracted.str().c_str());

						TCanvas *c2del = new TCanvas("","",3*850,2*850);
						c2del->Divide(3,2);
						// vpol
							c2del->cd(1);
								maps[0]->Draw("colz");
								maps[0]->GetXaxis()->SetTitle("Phi (deg)");
								maps[0]->GetYaxis()->SetTitle("Theta (deg)");
							c2del->cd(2);
								stacked_maps[0]->Draw("colz");
								stacked_maps[0]->SetTitle("VPol Averaged Stacked Map");
								stacked_maps[0]->GetXaxis()->SetTitle("Phi (deg)");
								stacked_maps[0]->GetYaxis()->SetTitle("Theta (deg)");
							c2del->cd(3);
								subtracted[0]->Draw("colz");
								subtracted[0]->GetXaxis()->SetTitle("Phi (deg)");
								subtracted[0]->GetYaxis()->SetTitle("Theta (deg)");
						// hpol
							c2del->cd(4);
								maps[1]->Draw("colz");
								maps[1]->GetXaxis()->SetTitle("Phi (deg)");
								maps[1]->GetYaxis()->SetTitle("Theta (deg)");
							c2del->cd(5);
								stacked_maps[1]->Draw("colz");
								stacked_maps[1]->SetTitle("HPol Averaged Stacked Map");
								stacked_maps[1]->GetXaxis()->SetTitle("Phi (deg)");
								stacked_maps[1]->GetYaxis()->SetTitle("Theta (deg)");
							c2del->cd(6);
								subtracted[1]->Draw("colz");
								subtracted[1]->GetXaxis()->SetTitle("Phi (deg)");
								subtracted[1]->GetYaxis()->SetTitle("Theta (deg)");
						char this_title[400];
						sprintf(this_title, "%s/maps_stacking/%d.%d.%d_Run%d_Event%d_OriginalAndSubtracted.png",plotPath,year_now, month_now, day_now,runNum,event);
						c2del->SaveAs(this_title); 
						delete c2del;
					}
					delete subtracted[0]; delete subtracted[1];
				}
			}
		}

		TCanvas *c2 = new TCanvas("","",2.1*850,2.1*850);
		c2->Divide(2,2);
		c2->cd(1);
			before_distro[0]->SetTitle("VPol Event Distribution, before subtraction");
			before_distro[0]->Draw("colz");
			before_distro[0]->GetXaxis()->SetTitle("Phi (deg)");
			before_distro[0]->GetYaxis()->SetTitle("Theta (deg)");
			before_distro[0]->GetZaxis()->SetRangeUser(0,1);
		c->cd(2);
			after_distro[0]->SetTitle("VPol Event Distribution, after subtraction");
			after_distro[0]->Draw("colz");
			after_distro[0]->GetXaxis()->SetTitle("Phi (deg)");
			after_distro[0]->GetYaxis()->SetTitle("Theta (deg)");
			after_distro[0]->GetZaxis()->SetRangeUser(0,1);
		c2->cd(3);
			before_distro[1]->SetTitle("HPol Event Distribution, before subtraction");
			before_distro[1]->Draw("colz");
			before_distro[1]->GetXaxis()->SetTitle("Phi (deg)");
			before_distro[1]->GetYaxis()->SetTitle("Theta (deg)");
			before_distro[1]->GetZaxis()->SetRangeUser(0,1);
		c->cd(4);
			after_distro[1]->SetTitle("HPol Event Distribution, after subtraction");
			after_distro[1]->Draw("colz");
			after_distro[1]->GetXaxis()->SetTitle("Phi (deg)");
			after_distro[1]->GetYaxis()->SetTitle("Theta (deg)");
			after_distro[1]->GetZaxis()->SetRangeUser(0,1);
		sprintf(title, "%s/maps_stacking/%d.%d.%d_A%d_Run%d_%dEvents_BeforeAndAfterSubtraction.png",plotPath,year_now, month_now, day_now,station,runNum,numStacked);
		c2->SaveAs(title); 
		delete c2;

		delete stacked_maps[0]; delete stacked_maps[1];
		inputFile->Close();
	}
}