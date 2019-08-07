////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  2019.02.25-OffsetBlockVTime-Plot.cxx 
////  plot histogram of fractions of spectrum power below 75 MHz
////
////  Dec 2018
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
#include "TDatime.h"
#include "TLegend.h"
#include "TLine.h"

#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraQualCuts.h"

#include "Settings.h"
#include "Detector.h"
#include "Report.h"
#include "RayTraceCorrelator.h"
AraAntPol::AraAntPol_t Vpol = AraAntPol::kVertical;
AraAntPol::AraAntPol_t Hpol = AraAntPol::kHorizontal;

#include "tools_PlottingFns.h"
#include "tools_Cuts.h"
#include "tools_RecoFns.h"
#include "tools_outputObjects.h"


using namespace std;

int main(int argc, char **argv)
{
	gStyle->SetOptStat(110011);
	gStyle->SetTimeOffset(0);
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	stringstream ss;

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;
	char *DataDirPath(getenv("DATA_DIR"));
	if (DataDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	
	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <station> <year> <runNum> <event>"<<endl;
		return 0;
	}
	int station = atoi(argv[1]);
	int year = atoi(argv[2]);
	int runNum = atoi(argv[3]);
	int event = atoi(argv[4]);

	char run_file_name[400];
	sprintf(run_file_name,"%s/RawData/A%d/all_runs/event%d.root",DataDirPath,station,runNum);
	TFile *mapFile = TFile::Open(run_file_name);
	if(!mapFile){
		cout<<"Can't open data file for map!"<<endl;
	}
	TTree *eventTree = (TTree*) mapFile-> Get("eventTree");
	if(!eventTree){
		cout<<"Can't find eventTree for map"<<endl;
	}

	RawAtriStationEvent *rawPtr =0;
	eventTree->SetBranchAddress("event",&rawPtr);
	eventTree->GetEvent(event);

	int unixTime = (int)rawPtr->unixTime;
	int unixTimeUs =(int)rawPtr->unixTimeUs;
	int timeStamp = (int)rawPtr->timeStamp;

	int stationID = rawPtr->stationId;
	char ped_file_name[400];
	sprintf(ped_file_name,"%s/run_specific_peds/A%d/all_peds/event%d_specificPeds.dat",PedDirPath,station,runNum);
	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	calibrator->setAtriPedFile(ped_file_name,stationID); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist
	
	UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawPtr,AraCalType::kLatestCalib);


	Settings *settings = new Settings();
	string setupfile = "setup.txt";
	settings->ReadFile(setupfile);
	cout << "Read " << setupfile << " file!" << endl;
	settings->NOFZ=1;
	Detector *detector=0;
	RayTraceCorrelator *theCorrelators[2];
	theCorrelators[0] =  new RayTraceCorrelator(station, 41., settings, 1, 4); //41 m, cal puser
	theCorrelators[1] =  new RayTraceCorrelator(station, 300., settings, 1, 4);//300 m, far reco

	char filter_file_name[400];
	sprintf(filter_file_name,"%s/ProcessedFile/A%d/%d/processed_station_%d_run_%d_filter.root",DataDirPath,station,year,station,runNum);
	TFile *filterFile = TFile::Open(filter_file_name);
	TTree *filterTree = (TTree*) filterFile->Get("OutputTree");
	filterTree->SetBranchAddress("VPeakOverRMS", &VPeakOverRMS);
	filterTree->GetEvent(event);

	vector <int> chan_list_V;
	vector <int> chan_list_H;
	vector<double> chan_SNRs;
	for(int chan=0; chan<=7; chan++){
		chan_list_V.push_back(chan);
		chan_list_H.push_back(chan+8);
	}
	for(int i=0; i<16; i++){
		// chan_SNRs.push_back(1.);
		chan_SNRs.push_back(VPeakOverRMS[i]);
		// printf("SNR is %.2f \n", VPeakOverRMS[i]);
	}
	filterFile->Close();

	if(station==2){
		//for station 2, we need to exclude channel 15 from the analysis
		chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
	}
	else if(station==3){
		//for station 3 years 2014 and 2015, we need to drop string 4 (channels 3, 7, 11, 15) altogether
		if(runNum>getA3BadRunBoundary()){
			chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 3), chan_list_V.end());
			chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 7), chan_list_V.end());
			chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 11), chan_list_H.end());
			chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
		}
	}

	TH2D *map_41m_V;
	TH2D *map_300m_V;
	TH2D *map_41m_H;
	TH2D *map_300m_H;

	bool AraSim=false;
	int solNum=0;
	map_41m_V = theCorrelators[0]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Vpol, AraSim, chan_list_V, chan_SNRs, solNum);
	map_300m_V = theCorrelators[1]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Vpol, AraSim, chan_list_V, chan_SNRs,  solNum);
	map_41m_H = theCorrelators[0]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Hpol, AraSim, chan_list_H, chan_SNRs, solNum);
	map_300m_H = theCorrelators[1]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Hpol, AraSim, chan_list_H, chan_SNRs, solNum);

	// map_41m_V = theCorrelators[0]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Vpol, AraSim, chan_list_V);
	// map_300m_V = theCorrelators[1]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Vpol, AraSim, chan_list_V);
	// map_41m_H = theCorrelators[0]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Hpol, AraSim, chan_list_H);
	// map_300m_H = theCorrelators[1]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Hpol, AraSim, chan_list_H);

	int PeakTheta_Recompute_41m_H;
	int PeakTheta_Recompute_300m_H;
	int PeakPhi_Recompute_41m_H;
	int PeakPhi_Recompute_300m_H;
	double PeakCorr_Recompute_41m_H;
	double PeakCorr_Recompute_300m_H;
	int PeakTheta_Recompute_41m_V;
	int PeakTheta_Recompute_300m_V;
	int PeakPhi_Recompute_41m_V;
	int PeakPhi_Recompute_300m_V;
	double PeakCorr_Recompute_41m_V;
	double PeakCorr_Recompute_300m_V;
	getCorrMapPeak(map_41m_H,PeakTheta_Recompute_41m_H,PeakPhi_Recompute_41m_H,PeakCorr_Recompute_41m_H);
	getCorrMapPeak(map_300m_H,PeakTheta_Recompute_300m_H,PeakPhi_Recompute_300m_H,PeakCorr_Recompute_300m_H);
	getCorrMapPeak(map_41m_V,PeakTheta_Recompute_41m_V,PeakPhi_Recompute_41m_V,PeakCorr_Recompute_41m_V);
	getCorrMapPeak(map_300m_V,PeakTheta_Recompute_300m_V,PeakPhi_Recompute_300m_V,PeakCorr_Recompute_300m_V);

	printf("	Rconstruction Information\n");
	printf("		41m H theta and phi %d and %d \n", PeakTheta_Recompute_41m_H, PeakPhi_Recompute_41m_H);
	stringstream ss30H;
	ss30H<<" 41m H Peak Theta, Phi is "<<PeakTheta_Recompute_41m_H<<" , "<<PeakPhi_Recompute_41m_H;
	map_41m_H->SetTitle(ss30H.str().c_str());
	printf("		300m H theta and phi %d and %d \n", PeakTheta_Recompute_300m_H, PeakPhi_Recompute_300m_H);
	stringstream ss300H;
	ss300H<<" 300m H Peak Theta, Phi is "<<PeakTheta_Recompute_300m_H<<" , "<<PeakPhi_Recompute_300m_H;
	map_300m_H->SetTitle(ss300H.str().c_str());
	printf("		41m V theta and phi %d and %d \n", PeakTheta_Recompute_41m_V, PeakPhi_Recompute_41m_V);
	stringstream ss30V;
	ss30V<<" 41m V Peak Theta, Phi is "<<PeakTheta_Recompute_41m_V<<" , "<<PeakPhi_Recompute_41m_V;
	map_41m_V->SetTitle(ss30V.str().c_str());
	printf("		300m V theta and phi %d and %d \n", PeakTheta_Recompute_300m_V, PeakPhi_Recompute_300m_V);
	stringstream ss300V;
	ss300V<<" 300m V Peak Theta, Phi is "<<PeakTheta_Recompute_300m_V<<" , "<<PeakPhi_Recompute_300m_V;
	map_300m_V->SetTitle(ss300V.str().c_str());

	beautify_TH2D();
	TCanvas *cMaps = new TCanvas("","",2*850,2*850);
	cMaps->Divide(2,2);
		cMaps->cd(3);
		map_41m_V->Draw("colz");
		// map_41m_V->GetXaxis()->SetRangeUser(50,80);
		// map_41m_V->GetYaxis()->SetRangeUser(-10,20);
		// map_41m_V->GetZaxis()->SetRangeUser(0.005,0.04);
		cMaps->cd(4);
		map_41m_H->Draw("colz");
		cMaps->cd(1);
		map_300m_V->Draw("colz");
		cMaps->cd(2);
		map_300m_H->Draw("colz");
	char save_temp_title[400];		
	sprintf(save_temp_title,"%s/single_events/%d.%d.%d_Run%d_Ev%d_Maps_SNRweighted.png",plotPath,year_now,month_now,day_now,runNum,event);
	cMaps->SaveAs(save_temp_title);
	delete cMaps;
	gStyle->SetOptStat(0);

	TH2D *map_41_V_clone_draw = (TH2D*) map_41m_V->Clone();
	TH2D *map_41_V_clone_adjust = (TH2D*) map_41m_V->Clone();

	int maxX;
	int maxY;
	int maxZ;
	map_41_V_clone_adjust->GetMaximumBin(maxX, maxY, maxZ);
	
	// switch to degree space
	maxY = maxY - 90;
	maxX = maxX - 180;

	int flexBoundary=5;
	bool isCloseToPulser=false;
	int config=6;
	identifyCalPulser(station,config,maxY,maxX,isCloseToPulser,isCloseToPulser,flexBoundary);

	printf("This is close to the pulser box! Check it further...\n");

	// switch back to bin space
	maxY+=90;
	maxX+=180;
	for(int binX=maxX-3; binX<maxX+3; binX++){
		for(int binY=maxY-3; binY<maxY+3; binY++){
			map_41_V_clone_adjust->SetBinContent(binX,binY,0.);
		}
	}

	int PeakTheta_Recompute_41m_V_mask;
	int PeakPhi_Recompute_41m_V_mask;
	double PeakCorr_Recompute_41m_V_mask;
	getCorrMapPeak(map_41_V_clone_adjust,PeakTheta_Recompute_41m_V_mask,PeakPhi_Recompute_41m_V_mask,PeakCorr_Recompute_41m_V_mask);
	stringstream ss30V_mask;
	ss30V_mask<<" Masked 41m V Peak Theta, Phi is "<<PeakTheta_Recompute_41m_V_mask<<" , "<<PeakPhi_Recompute_41m_V_mask;
	map_41_V_clone_adjust->SetTitle(ss30V_mask.str().c_str());

	bool newCP5, newCP6 = false;
	identifyCalPulser(station,config,PeakTheta_Recompute_41m_V_mask,PeakPhi_Recompute_41m_V_mask,newCP5, newCP6);
	printf("New CP5 and CP6 status is %d and %d \n", newCP5, newCP6);

	double real_phi, real_theta;
	getRealLocation(station,5,1,real_theta,real_phi);
	TLine theta_line_2D(-40,double(real_theta),-10,double(real_theta));
	TLine phi_line_2D(real_phi,-40,real_phi,-10);

	TCanvas *cMaps2 = new TCanvas("","",3*850,850);
	cMaps2->Divide(3,1);
		cMaps2->cd(1);
			map_41m_V->Draw("colz");
		cMaps2->cd(2);
			map_41_V_clone_draw->Draw("colz");
			map_41_V_clone_draw->GetXaxis()->SetRangeUser(-40,-10);
			map_41_V_clone_draw->GetYaxis()->SetRangeUser(-40,-10);
			map_41_V_clone_draw->GetZaxis()->SetRangeUser(0.005,0.018);
		cMaps2->cd(3);
			map_41_V_clone_adjust->Draw("colz");
			map_41_V_clone_adjust->GetXaxis()->SetRangeUser(-40,-10);
			map_41_V_clone_adjust->GetYaxis()->SetRangeUser(-40,-10);
			map_41_V_clone_adjust->GetZaxis()->SetRangeUser(0.005,0.018);
			phi_line_2D.Draw("");
			theta_line_2D.Draw("");
			phi_line_2D.SetLineStyle(9);
			theta_line_2D.SetLineStyle(9);
			theta_line_2D.SetLineWidth(2);
			phi_line_2D.SetLineWidth(2);
	sprintf(save_temp_title,"%s/single_events/%d.%d.%d_Run%d_Ev%d_Maps_SNRweighted_Zoom.png",plotPath,year_now,month_now,day_now,runNum,event);
	cMaps2->SaveAs(save_temp_title);
	delete cMaps2;

	delete map_41m_V; delete map_300m_V; delete map_41m_H; delete map_300m_H;

	delete realAtriEvPtr;
	mapFile->Close();
	delete mapFile;
	return 0;
}