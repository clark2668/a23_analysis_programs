////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	2019.08.26_MapForJorge.cxx 
////	make map for Jorge
////
////	Aug 2019
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <iomanip>
#include <sstream>

//AraRoot Includes
#include "RawAtriStationEvent.h"
#include "UsefulAraStationEvent.h"
#include "UsefulAtriStationEvent.h"

//ROOT Includesf
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

	char *DataDirPath(getenv("DATA_DIR"));
	if (DataDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL) std::cout << "Warning! $PED_DIR is not set!" << endl;
	
	if(argc<5) {
		std::cout << "Usage\n" << argv[0] << " <1-simulation_flag> <2-station> <3-filter_file_dir> <4-input file>\n";
		return -1;
	}

	/*
	arguments
	0: exec
	1: simulation (yes/no)
	2: station num (2/3)
	3: filter file dir
	4: input file
	*/

	isSimulation=atoi(argv[1]);
	int station_num=atoi(argv[2]);
	int yearConfig=atoi(argv[3]);
	
	TFile *fp = TFile::Open(argv[4]);
	if(!fp) {
		std::cout << "Can't open file\n";
		return -1;
	}
	TTree *eventTree; 
	eventTree= (TTree*) fp->Get("eventTree");
	if(!eventTree) {
		std::cout << "Can't find eventTree\n";
		return -1;
	}
			
	Settings *settings = new Settings();
	string setupfile = "setup.txt";
	settings->ReadFile(setupfile);
	cout << "Read " << setupfile << " file!" << endl;
	settings->NOFZ=1;
	Detector *detector = 0;

	RayTraceCorrelator *theCorrelator = new RayTraceCorrelator(station_num, 41, settings, 1, RTTestMode);
	
	double weight;
	int unixTime;
	int unixTimeUs;
	int eventNumber;
	double maxPeakVfromSim;
	double PeakVfromSim[16][2];
	int runNum;
	
	if(isSimulation){
		eventTree->SetBranchAddress("UsefulAtriStationEvent", &realAtriEvPtr);
		eventTree->SetBranchAddress("weight", &weight);
		printf("Simulation; load useful event tree straight away \n");
		runNum = getrunNum(argv[7]);
	}
	else{
		eventTree->SetBranchAddress("event",&rawAtriEvPtr);
		eventTree->SetBranchAddress("run",&runNum);
		printf("Data; load raw event tree \n");

	}
	
	Long64_t numEntries=eventTree->GetEntries();
	Long64_t starEvery=numEntries/80;
	if(starEvery==0) starEvery++;

	eventTree->GetEntry(0); //just to get runNum
	printf("Reco Run Number %d \n", runNum);

	// load pedestals, if they exist yet
	// this will perform fine for simulation; it'll evaluate to garbage (or whatever), but it's not needed, so cool
	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	char ped_file_name[400];
	sprintf(ped_file_name,"%s/run_specific_peds/A%d/all_peds/event%d_specificPeds.dat",PedDirPath,station_num,runNum);
	calibrator->setAtriPedFile(ped_file_name,station_num); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist

	// get the run summary information, if it exists yet
	// and remember, because it's the users job to pass the location of the filter files
	// this should work for simulated events just fine
	char filter_file_name[400];
	sprintf(filter_file_name,"%s/processed_station_%d_run_%d_filter.root",argv[3],station_num,runNum);
	bool hasFilterFile = false;
	TFile *filterFile = TFile::Open(filter_file_name);
	TTree *filterTree;
	if(filterFile){
		printf("Successfully found filter file information \n");
		hasFilterFile=true;
		filterTree = (TTree*) filterFile->Get("OutputTree");
		if(!filterTree) {
			std::cout << "Can't find filterTree\n";
			return -1;
		}
		filterTree->SetBranchAddress("VPeakOverRMS", &VPeakOverRMS);
		filterFile->cd();
	}
	
	int eventSim = 0;
	cerr<<"Run "<<runNum<<" has a starEvery of "<<starEvery<<" and "<<numEntries<<" total events"<<endl;
	numEntries=20;
	for(Long64_t event=0;event<numEntries;event++) {
		if(event%starEvery==0) {
			std::cout << "*";     
		}

		eventTree->GetEntry(event);
		if(hasFilterFile){
			filterTree->GetEntry(event);
		}

		if (isSimulation == false){
			unixTime=(int)rawAtriEvPtr->unixTime;
			unixTimeUs=(int)rawAtriEvPtr->unixTimeUs;
			eventNumber=(int)rawAtriEvPtr->eventNumber;
		}else {
			eventNumber = event;
		}
	
		if(!isSimulation){
			realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
		}
		
		if (isSimulation){
			isCalpulser = false;
			isSoftTrigger = false;
		} else{
			isCalpulser = rawAtriEvPtr->isCalpulserEvent();
			isSoftTrigger = rawAtriEvPtr->isSoftwareTrigger();
		}
		if(!isCalpulser){
			continue;
		}

		bool analyzeEvent = true;
		if (analyzeEvent == true){

			vector<double> chan_SNRs;
			if(hasFilterFile){
				for(int i=0; i<16; i++){
					chan_SNRs.push_back(VPeakOverRMS[i]);
				}
			}

			vector <int> chan_list_V;
			vector <int> chan_list_H;
			for(int chan=0; chan<=7; chan++){
				chan_list_V.push_back(chan);
				chan_list_H.push_back(chan+8);
			}

			if(station_num==2){
				//for station 2, we need to exclude channel 15 from the analysis
				chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
			}
			else if(station_num==3){
				// for station 3 years 2014, 2015, 2016, we need to drop string 4 (channels 3, 7, 11, 15) altogether above some run
				if( 
					(!isSimulation && runNum>getA3BadRunBoundary())
					||
					(isSimulation && yearConfig>2)

				){			// drop string four
					chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 3), chan_list_V.end());
					chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 7), chan_list_V.end());
					chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 11), chan_list_H.end());
					chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
				}
			}
			int solNum = 0;
			TH2D *map_V_raytrace = theCorrelator->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V, chan_SNRs, solNum);
			TH2D *map_H_raytrace = theCorrelator->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H, chan_SNRs, solNum);

			// getCorrMapPeak_wStats(map_V_raytrace, peakTheta_single[0], peakPhi_single[0], peakCorr_single[0], minCorr_single[0], meanCorr_single[0], rmsCorr_single[0], peakSigma_single[0]);
			// getCorrMapPeak_wStats(map_H_raytrace, peakTheta_single[1], peakPhi_single[1], peakCorr_single[1], minCorr_single[1], meanCorr_single[1], rmsCorr_single[1], peakSigma_single[1]);

			map_V_raytrace->GetZaxis()->SetRangeUser(0, 0.04);
			map_H_raytrace->GetZaxis()->SetRangeUser(0, 0.006);

			bool print_maps = true;
			if(print_maps){
				gStyle->SetOptStat(0);
				TCanvas *cMaps = new TCanvas("","",2*1100,2*850);
				cMaps->Divide(2,2);
					cMaps->cd(1);
					map_V_raytrace->Draw("colz");
					cMaps->cd(2);
					map_H_raytrace->Draw("colz");
				char save_temp_title[400];		
				sprintf(save_temp_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/%d.%d.%d_Run%d_Ev%d_Maps_FromRecoCode_forJorge.png",year_now,month_now,day_now,runNum,event);
				cMaps->SaveAs(save_temp_title);
				delete cMaps;
			}

			delete map_V_raytrace;
			delete map_H_raytrace;
				
			if (isSimulation == false) {
				delete realAtriEvPtr;
			}
		}
	}
	
	if(hasFilterFile)
		filterFile->Close();
	fp->Close();

	delete theCorrelator;
	delete settings;
}
