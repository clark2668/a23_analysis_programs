////////////////////////////////////////////////////////////////////////////////
////	v2_CWID.cxx 
////	A23 diffuse, identify CW freequency
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
#include "TF1.h"

#include "tools_inputParameters.h"
#include "tools_WaveformFns.h"
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_CW.h"
#include "tools_Cuts.h"
#include "tools_CommandLine.h"

RawAtriStationEvent *rawAtriEvPtr;
UsefulAtriStationEvent *realAtriEvPtr;

int main(int argc, char **argv)
{

	if(argc<8) {
		std::cout << "Usage\n" << argv[0] << " <simulation_flag> <station> <year> <drop_bad_chans> <run summary directory> <output directory> <input file> <pedestal file> \n";
		return -1;
	}

	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;

	/*
	arguments
	0: exec
	1: simulation (yes/no)
	2: station num (2/3)
	3: year
	4: drop bad chans
	5: run summary directory
	6: output directory
	7: input file
	8: pedestal file
	*/

	int isSimulation=atoi(argv[1]);
	int station_num=atoi(argv[2]);
	int year=atoi(argv[3]);
	int drop_bad_chans=atoi(argv[4]);

	//open the file
	TFile *fp = TFile::Open(argv[7]);
	if(!fp) {
		std::cerr << "Can't open file\n";
		return -1;
	}
	TTree *eventTree= (TTree*) fp->Get("eventTree");
	if(!eventTree) {
		std::cerr << "Can't find eventTree\n";
		return -1;
	}
	int runNum = getrunNum(argv[7]);
	printf("Run Number %d \n", runNum);

	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();	
	if(argc==9){
		//only if they gave us a pedestal should we fire up the calibrator
		calibrator->setAtriPedFile(argv[8],station_num);
	}
	else{
		if(!isSimulation){
			char ped_file_name[400];
			sprintf(ped_file_name,"%s/run_specific_peds/A%d/all_peds/event%d_specificPeds.dat",PedDirPath,station_num,runNum);
			calibrator->setAtriPedFile(ped_file_name,station_num); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist
		}
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
	cout<<"This file has a starEvery of "<<starEvery<<endl;

	//first, let's get the baselines loaded in
	string runSummaryFilename;
	if (!isSimulation){
		runSummaryFilename = getRunSummaryFilename(station_num, argv[5], argv[7]);
	}
	else {
		if(station_num==2){
			runSummaryFilename = "/fs/scratch/PAS0654/ara/sim/RunSummary/run_summary_station_2_run_20.root";
		}
		else if(station_num==3){
			runSummaryFilename = "/fs/scratch/PAS0654/ara/sim/RunSummary/run_summary_station_3_run_30.root";
		}
	}
	TFile *SummaryFile = TFile::Open(runSummaryFilename.c_str());
	// TFile *SummaryFile = TFile::Open("/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/RunSummary/A3/2013/run_summary_station_3_run_1550.root");
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

	bool baselinesGood = areBaselinesGood(average);
	cout<<"Baselines are good? : "<<baselinesGood<<endl;

	if(!baselinesGood){
		cout<<"The baselines are bad! Pull a default!"<<endl;

		// first, close up shop


		int numTries=0;
		int runToTry=runNum-1;
		while(numTries<10){
			
		}
	}

	char equation[150];
	sprintf(equation,"([3]*x*x*x + [2]*x*x + [0]*x + [1])");
	char equation_name[16][150];
	TF1 *fit[16];
	double start_of_fit[16];
	double end_of_fit[16];
	for(int pol=0; pol<16; pol++){
		start_of_fit[pol]=150.;
		end_of_fit[pol]=850.;
		sprintf(equation_name[pol],"Fit%d",pol);
		fit[pol] = new TF1(equation_name[pol],equation,start_of_fit[pol],end_of_fit[pol]);
		average[pol]->Fit(equation_name[pol],"Q,R");
		// printf("Pol %d Fit Parameters are %.2f and %.2f \n", pol, fitParams[pol][0], fitParams[pol][1]);
	}
	// average[0]->Fit(equation_name[0],"LL,R");
	// fitParams[pol][0] = fit[pol]->GetParameter(0);
	// fitParams[pol][1] = fit[pol]->GetParameter(1);
	// fitParamErrors[pol][0] = fit[pol]->GetParError(0);
	// fitParamErrors[pol][1] = fit[pol]->GetParError(1);
	// TF1 *f1 = new TF1("f1", "[0]+[1]*x", 500, 750);
	// average[0]->Fit(f1);

	TH1D *violations = new TH1D("","",50,-10,10);

	for(int chan=0; chan<16; chan++){
		int numSamps = average[chan]->GetN();
		Double_t *theseY = average[chan]->GetY();
		Double_t *theseX = average[chan]->GetX();
		for(int samp=0; samp<numSamps; samp++){
			double thisX = theseX[samp];
			double thisY = theseY[samp];
			// if(thisX>150. && thisX<850. && (thisX-300.>1.) && (thisX-500.>1.)){
			if(abs(thisX-300.)<2.) continue;
			if(abs(thisX-500.)<2.) continue;
			if(thisX>150. && thisX<800.){
				double expectedY = fit[chan]->Eval(thisX);
				double violation = thisY - expectedY;
				if(chan==0)
					violations->Fill(violation);
				if( violation > 2.){
					printf("Chan %d, Freq %.2f has %.2f violation \n", chan, thisX, violation);
				}
			}
		}
	}

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;
	char save_temp_title[300];

	TCanvas *cViolation = new TCanvas("","",1100,850);
		violations->Draw("");
		violations->GetXaxis()->SetTitle("Fit - Baseline");
		violations->GetYaxis()->SetTitle("Number of Frequency Bins");
	sprintf(save_temp_title,"%s/trouble_events/Run%d_DiffBetweenFitAndBaseline.png",plotPath,runNum);
	gPad->SetLogy();
	cViolation->SaveAs(save_temp_title);
	delete cViolation;

	TCanvas *cAvg = new TCanvas("","",4*1100,4*850);
	cAvg->Divide(4,4);
	for(int i=0; i<16; i++){
		cAvg->cd(i+1);
		average[i]->Draw("AL");
		average[i]->GetYaxis()->SetRangeUser(15,40);
		average[i]->GetYaxis()->SetTitle("Power (dB)");
		average[i]->GetXaxis()->SetTitle("Frequency (MHz)");
		if(i==1){
			TLegend *leg = new TLegend(0.5,0.4,0.9,0.2);
			leg->AddEntry(average[0],"Baseline","l");
			leg->AddEntry(fit[0],"Poly Fit","l");
			leg->Draw();
		}
	}
	sprintf(save_temp_title,"%s/trouble_events/Run%d_CWBaselineOnly.png",plotPath,runNum);
	cAvg->SaveAs(save_temp_title);
	delete cAvg;


	// return -1;

	AraGeomTool * geomTool = new AraGeomTool();
	int nGraphs=16;

	vector<int> chan_exclusion_list;
	if(drop_bad_chans){
		if(station_num==2){
			//drop the last hpol channel
			chan_exclusion_list.push_back(15);
		}
		if(station_num==3 && runNum>2901){
			// drop string four
			chan_exclusion_list.push_back(3);
			chan_exclusion_list.push_back(7);
			chan_exclusion_list.push_back(11);
			chan_exclusion_list.push_back(15);
		}
	}

	//now, to loop over events!
	for(Long64_t event=0;event<numEntries;event++){
		// cout<<"On event "<<event<<endl;
   
		if(event%starEvery==0) { std::cerr << "*"; }
    
		eventTree->GetEntry(event); //get the event
		int eventNumber = rawAtriEvPtr->eventNumber;
		if(eventNumber!=394) continue;
		if(eventNumber==394){
			printf(RED"Found event %d, unixTime %d \n"RESET, rawAtriEvPtr->eventNumber, rawAtriEvPtr->unixTime);
		}

		if (isSimulation == false){
			realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
		}

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
		vector<vector<double> > badFreqs_baseline;

		if(1==1){

			//before we do the phase variance, we should check for baseline violations	
			vector<double> baseline_CW_cut_V = CWCut_TB(grWaveformsRaw, average, 0, 6., 5.5, station_num, 3, chan_exclusion_list, runNum, eventNumber, true);
			vector<double> baseline_CW_cut_H = CWCut_TB(grWaveformsRaw, average, 1, 6., 5.5, station_num, 3, chan_exclusion_list, runNum, eventNumber, true);
						
			for(int i=0; i<baseline_CW_cut_V.size(); i++){
				printf("V: Event %d Baseline CW Cut %.2f \n", event, baseline_CW_cut_V[i]);
			}
			for(int i=0; i<baseline_CW_cut_H.size(); i++){
				printf("H: Event %d Baseline CW Cut %.2f \n", event, baseline_CW_cut_H[i]);
			}
			

			badFreqs_baseline.push_back(baseline_CW_cut_V);
			badFreqs_baseline.push_back(baseline_CW_cut_H);
		}

		deleteGraphVector(grWaveformsRaw);
		if(isSimulation==false){
			delete realAtriEvPtr;
		}
	}//loop over events


	SummaryFile->Close();

	printf("Done! Run Number %d \n", runNum);

}//end main