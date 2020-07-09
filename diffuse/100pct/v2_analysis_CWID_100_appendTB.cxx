////////////////////////////////////////////////////////////////////////////////
////	v2_CWID.cxx 
////	A23 diffuse, identify CW frequency
////	We only want to do CW on events that pas the filter
////	So we will need the filter file
////
////	Nov 2018
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <iomanip>
#include <sstream>
#include <deque>
#include <vector>

//AraRoot Includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraGeomTool.h"
#include "AraAntennaInfo.h"
#include "AraQualCuts.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"

#include "tools_inputParameters.h"
#include "tools_WaveformFns.h"
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_CW.h"
#include "tools_CommandLine.h"
#include "tools_outputObjects.h"

RawAtriStationEvent *rawAtriEvPtr;
UsefulAtriStationEvent *realAtriEvPtr;

int main(int argc, char **argv)
{

	if(argc<9) {
		std::cout << "Usage\n" << argv[0] << " <1-simulation_flag> <2-station> <3-year/config> <4-drop_bad_chans> <5-filter_file_dir> <6-run summary directory> <7-output directory> <8-input file> <9-pedestal file> \n";
		return -1;
	}

	// check for TOOLS_DIR, will need this to get default baselines
	char *toolsPath(getenv("TOOLS_DIR"));
	if (toolsPath == NULL) {
		std::cout << "Warning! $TOOLS_DIR is not set!" << endl;
		return -1;
	}


	/*
	arguments
	0: exec
	1: simulation (yes/no)
	2: station num (2/3)
	3: year/ config
	4: drop bad chans
	5: filter file directory
	6: run summary directory
	7: output directory
	8: input file
	9: pedestal file
	*/

	int isSimulation=atoi(argv[1]);
	int station_num=atoi(argv[2]);
	int yearConfig=atoi(argv[3]);
	int drop_bad_chans=atoi(argv[4]);
	drop_bad_chans=1; //always drop bad channels...

	//open the file
	TFile *fp = TFile::Open(argv[8]);
	if(!fp) {
		std::cerr << "Can't open file\n";
		return -1;
	}
	TTree *eventTree= (TTree*) fp->Get("eventTree");
	if(!eventTree) {
		std::cerr << "Can't find eventTree\n";
		return -1;
	}
	int runNum = getrunNum(argv[8]);
	int config = getConfig(station_num, runNum);
	printf("Filter Run Number %d \n", runNum);

	/*
		Adjust filter cut settings
		in A2, we were happy with snr in 0,0 and -1.3, -1.4 in V and H
		in A3, we want to do something a little more complicated and just thresholds
		by configuration
	*/
	int thresholdBin_pol[]={0,0};
	double wavefrontRMScut[]={2., 2.};

	int thisConfiguration = yearConfig;
	if(!isSimulation){
		thisConfiguration = getConfig(station_num, runNum);
	}
	getWFRMSCutValues(station_num, thisConfiguration, thresholdBin_pol[0], thresholdBin_pol[1], wavefrontRMScut[0], wavefrontRMScut[1]);
	// printf("Wavefront RMS cut values for config %d are vBin %d and hBin %d, vCut %.2f, hCut %.2f \n", thisConfiguration, thresholdBin_pol[0], thresholdBin_pol[1], wavefrontRMScut[0], wavefrontRMScut[1]);

	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();	
	if(argc==10){
		//only if they gave us a pedestal should we fire up the calibrator
		calibrator->setAtriPedFile(argv[9],station_num);
	}
	else{
		calibrator->setAtriPedFile("", station_num);
	}
	AraQualCuts *qualCut = AraQualCuts::Instance(); //we also need a qual cuts tool

	if(isSimulation){
		eventTree->SetBranchAddress("UsefulAtriStationEvent", &realAtriEvPtr);
		printf("Simulation; load useful event tree straight away \n");
	}
	else{
		eventTree->SetBranchAddress("event",&rawAtriEvPtr);
		printf("Data; load raw event tree \n");
	}
  
  	Long64_t numEntries=eventTree->GetEntries();
	Long64_t starEvery=numEntries/80;
	if(starEvery==0) starEvery++;
	printf("Num events is %d \n", numEntries);
	// numEntries=3000;
	cout<<"This file has a starEvery of "<<starEvery<<endl;

	//first, let's get the baselines loaded in
	string runSummaryFilename;
	if (!isSimulation){
		runSummaryFilename = getRunSummaryFilename(station_num, argv[6], argv[8]);
	}
	else {
		if(station_num==2){
			runSummaryFilename = "/fs/scratch/PAS0654/ara/sim/RunSummary/run_summary_station_2_run_20.root";
		}
		else if(station_num==3){
			runSummaryFilename = "/fs/scratch/PAS0654/ara/sim/RunSummary/run_summary_station_3_run_30.root";
		}
	}
	TFile *SummaryFile = TFile::Open(runSummaryFilename.c_str(),"read");
	if(!SummaryFile) {
		std::cerr << "Can't open summary file\n";
		return -1;
	}
	TTree* SummaryTree = (TTree*) SummaryFile->Get("BaselineTree");   
	if(!SummaryTree) {
		std::cerr << "Can't find SummaryTree\n";
		return -1;
	}
	vector <TGraph*> average;
	average.resize(16);
	stringstream ss1;
	for(int i=0; i<16; i++){
		ss1.str(""); ss1<<"baselines_RF_chan_"<<i;
		SummaryTree->SetBranchAddress(ss1.str().c_str(),&average[i]);
	}
	SummaryTree->GetEntry(0);
	vector<TGraph*> average_copy;
	for(int i=0; i<16; i++){
		average_copy.push_back((TGraph*)average[i]->Clone());
	}


	/*
		Now, we test and see if these averages are "good"/worth using
	*/
	bool baselinesGood = areBaselinesGood(average_copy,station_num,runNum);
	deleteGraphVector(average_copy); // clean this up, we don't need the copies anymore

	if(!baselinesGood){
		printf(RED"Warning! Run %d has bad baeslines, preparing for swap to averages...\n"RESET,runNum);
		char default_filename[500];
		sprintf(default_filename,"%s/data/average_baseline_station%d_config%d.root",toolsPath,station_num,config);
		TFile *DefaultFile = TFile::Open(default_filename);
		if(!DefaultFile) { std::cout << "Can't open default baseline file, and run "<<runNum<<" needs it\n"; return -1;}
		TTree* AverageBaselineTree = (TTree*) DefaultFile->Get("AverageBaselineTree");   
		if(!AverageBaselineTree) { std::cerr << "Can't find AverageBaselineTree for run "<<runNum<<" which needs it.\n"; return -1;}
		vector <TGraph*> config_averages;
		config_averages.resize(16);
		stringstream ss2;
		for(int i=0; i<16; i++){
			ss2.str(""); ss2<<"baselines_average_RF_chan_"<<i;
			AverageBaselineTree->SetBranchAddress(ss2.str().c_str(),&config_averages[i]);
		}
		AverageBaselineTree->GetEntry(0);
		for(int i=0; i<16; i++){
			delete average[i]; //delete
		}
		average.clear(); //clear this
		for(int i=0; i<16; i++){
			average.push_back((TGraph*)config_averages[i]->Clone()); //push back the new thing
		}
		// now that we've *cloned* the result, we can close up shop right away
		// and now "average" will have been replaced with the default tree, which is exactly what we wanted
		DefaultFile->Close();
		delete DefaultFile;
	}

	AraGeomTool * geomTool = new AraGeomTool();
	int nGraphs=16;

	vector<int> chan_exclusion_list;
	if(drop_bad_chans){
		if(station_num==2){
			//drop the last hpol channel
			chan_exclusion_list.push_back(15);
		}
		else if(station_num==3){
			if( 
				(!isSimulation && runNum>getA3BadRunBoundary())
				||
				(isSimulation && yearConfig>2)

			){			// drop string four
				printf("Yes, we need to drop channels!\n");
				chan_exclusion_list.push_back(3);
				chan_exclusion_list.push_back(7);
				chan_exclusion_list.push_back(11);
				chan_exclusion_list.push_back(15);
			}
		}
	}


	// we do *not* need to reload the filter file in this append operation
	// w were smart and included a "does/doesn't" pass WFRMS flag already
	// so we can just check that

	//now set up the outputs
	string output_location = argv[7];
	char run_file_name[400];
	sprintf(run_file_name,"%s/CWID_station_%d_run_%d.root",output_location.c_str(),station_num, runNum);
	TFile *outFile = TFile::Open(run_file_name,"update");
	TTree *NewCWTree = (TTree*)outFile->Get("NewCWTree");
	int numEntries_org_file = NewCWTree->GetEntries();
	// numEntries_org_file=3000;
	if(numEntries_org_file!=numEntries){
		cout<<"Warning! Not the same number of events in both. Halt!"<<endl;
		outFile->Close();
		return -1;
	}
	// the file should *already* have this does/doesn't pass wavefront RMS flag
	// so we can just call it!
	bool passesWavefrontRMS[2];
	NewCWTree->SetBranchAddress("passesWavefrontRMS", &passesWavefrontRMS);

	// and, set the output
	vector<vector<double> > badFreqs_baseline;
	TBranch *newBranch = NewCWTree->Branch("badFreqs_baseline_TB",&badFreqs_baseline);

	printf(BLUE"About to loop over events\n"RESET);

	//now, to loop over events!
	for(Long64_t event=0;event<numEntries;event++){

		badFreqs_baseline.clear();
    
		if(event%starEvery==0) { std::cerr << "*"; }


		NewCWTree->GetEntry(event);
		eventTree->GetEntry(event); //get the event

		bool isCalPulser;
		bool isSoftTrigger;
		if(isSimulation==false){
			isCalPulser = rawAtriEvPtr->isCalpulserEvent();
			isSoftTrigger = rawAtriEvPtr->isSoftwareTrigger();
		}
		else if (isSimulation){
			isCalPulser=false;
			isSoftTrigger=false;;
		}


		if(isCalPulser || isSoftTrigger){
			// printf("No need to reco event %d this event! WFRMS are %.2f, %.2f \n",event, TMath::Log10(bestFaceRMS[0]), TMath::Log10(bestFaceRMS[1]));
			newBranch->Fill(); //fill this anyway with garbage
			continue; //don't do any further processing on this event
		}

		// if it doesn't pass WFRMS filter in one of the two pols, no need to do reconstruction!
		if(!passesWavefrontRMS[0] && !passesWavefrontRMS[1]){
			// printf("No need to reco event %d this event! WFRMS are %.2f, %.2f \n",event, TMath::Log10(bestFaceRMS[0]), TMath::Log10(bestFaceRMS[1]));
			newBranch->Fill(); //fill this anyway with garbage
			continue; //don't do any further processing on this event
		}

		// printf("Event %d good for CWID! CalPul %d, Soft %d, WFRMS are %.2f, %.2f \n", event, isCalPulser, isSoftTrigger, TMath::Log10(bestFaceRMS[0]), TMath::Log10(bestFaceRMS[1]));

		if (isSimulation == false){
			realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
		}

		bool hasError=false;

		if (isSimulation == false){
			hasError = !(qualCut->isGoodEvent(realAtriEvPtr));
			if(realAtriEvPtr->eventNumber<5){
				hasError=true;
			}
			vector<TGraph*> spareElecChanGraphs;
			spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(6));
			spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(14));
			spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(22));
			spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(30));
			bool hasBadSpareChanIssue = hasSpareChannelIssue(spareElecChanGraphs);
			deleteGraphVector(spareElecChanGraphs);
			if(hasBadSpareChanIssue){
				hasError=true;
			}
		}
		else if(isSimulation){
			hasError=false;
		}

		// if(!hasError){
		// for station 3, skip cal pulsers
		// so only proceed if it doen't have an error
		// and if it's not (A3 and a cal pulser)
		// this is done for backwards compatibility with A2, which did not exclude cal pulsers
		// from phase variance search
		if(!hasError && !(isCalPulser && station_num==3)){

			cout<<"Redoing event "<<rawAtriEvPtr->eventNumber<<endl;

			stringstream ss;
			string xLabel, yLabel;
			xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
			vector<string> titlesForGraphs;
			for (int i = 0; i < nGraphs; i++){
				ss.str("");
				ss << "Channel " << i;
				titlesForGraphs.push_back(ss.str());
			}


			vector<TGraph*> grWaveformsRaw = makeGraphsFromRF(realAtriEvPtr, 16, xLabel, yLabel, titlesForGraphs);


			//before we do the phase variance, we should check for baseline violations	
			vector<double> baseline_CW_cut_V = CWCut_TB(grWaveformsRaw, average, 0, 6., 5.5, station_num, 3, chan_exclusion_list, runNum, event,false);
			vector<double> baseline_CW_cut_H = CWCut_TB(grWaveformsRaw, average, 1, 6., 5.5, station_num, 3, chan_exclusion_list, runNum, event,false);
			
			
			// for(int i=0; i<baseline_CW_cut_V.size(); i++){
			// 	printf(CYAN"	V: Event %d Baseline CW Cut %.2f \n"RESET, event, baseline_CW_cut_V[i]);
			// }
			// for(int i=0; i<baseline_CW_cut_H.size(); i++){
			// 	printf(CYAN"	H: Event %d Baseline CW Cut %.2f \n"RESET, event, baseline_CW_cut_H[i]);
			// }
			

			badFreqs_baseline.push_back(baseline_CW_cut_V);
			badFreqs_baseline.push_back(baseline_CW_cut_H);

			deleteGraphVector(grWaveformsRaw);
		}
		newBranch->Fill();
		if(isSimulation==false){
			delete realAtriEvPtr;
		}
	}//loop over events

	outFile->Write(0,TObject::kOverwrite);
	outFile->Close();
	fp->Close();
	SummaryFile->Close();

	printf("Done! Run Number %d \n", runNum);

}//end main
