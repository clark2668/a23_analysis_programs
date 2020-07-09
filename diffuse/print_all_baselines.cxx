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
#include "tools_CommandLine.h"

RawAtriStationEvent *rawAtriEvPtr;
UsefulAtriStationEvent *realAtriEvPtr;

int main(int argc, char **argv)
{

	if(argc<8) {
		std::cout << "Usage\n" << argv[0] << " <config> <station> <year> <drop_bad_chans> <run summary directory> <output directory> <input file> <input file 1> \n";
		return -1;
	}

	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

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

	int config=atoi(argv[1]);
	int station_num=atoi(argv[2]);
	int year=atoi(argv[3]);
	int drop_bad_chans=atoi(argv[4]);

	vector<TGraph*> mega_average;
	bool foundMegaAverage=false;
	int numInAverage=0;
	char save_temp_title[300];

	for(int file=7; file<argc; file++){

		//open the file
		TFile *fp = TFile::Open(argv[file]);
		if(!fp) {
			std::cerr << "Can't open file\n";
			continue;
			// return -1;
		}
		TTree *eventTree= (TTree*) fp->Get("eventTree");
		if(!eventTree) {
			std::cerr << "Can't find eventTree\n";
			fp->Close();
			continue;
			// return -1;
		}
		int runNum;
		eventTree->SetBranchAddress("run",&runNum);
		eventTree->GetEntry(0);

		// int runNum = getrunNum(argv[file]);
		// printf("Run Number %d \n", runNum);

		//first, let's get the baselines loaded in
		string runSummaryFilename;
		// runSummaryFilename = getRunSummaryFilename(station_num, argv[5], argv[file]);
		char summaryFileName[500];
		// sprintf(summaryFileName,"/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/RunSummary/A3/2013/run_summary_station_3_run_%d.root",runNum);
		sprintf(summaryFileName,"/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/RunSummary/A3/by_config/c%d/run_summary_station_%d_run_%d.root",config,station_num,runNum);
		// TFile *SummaryFile = TFile::Open("/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/RunSummary/A3/2013/run_summary_station_3_run_1550.root");
		// TFile *SummaryFile = TFile::Open(runSummaryFilename.c_str());
		TFile *SummaryFile = TFile::Open(summaryFileName);
		if(!SummaryFile) {
			std::cerr << "Can't open summary file\n";
			continue;
			// return -1;
		}
		TTree* SummaryTree = (TTree*) SummaryFile->Get("BaselineTree");   
		if(!SummaryTree) {
			std::cerr << "Can't find SummaryTree\n";
			SummaryFile->Close();
			continue;
			// return -1s;
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

		// TCanvas *cAvg = new TCanvas("","",4*1100,4*850);
		// cAvg->Divide(4,4);
		// for(int i=0; i<16; i++){
		// 	cAvg->cd(i+1);
		// 	average[i]->Draw("AL");
		// 	average[i]->GetYaxis()->SetRangeUser(15,40);
		// 	average[i]->GetYaxis()->SetTitle("Power (dB)");
		// 	average[i]->GetXaxis()->SetTitle("Frequency (MHz)");
		// }
		// // sprintf(save_temp_title,"%s/trouble_events/Run%d_Baselines_Good%d.png",plotPath,runNum,baselinesGood);
		// sprintf(save_temp_title,"/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/other_studies/baselines/A%d/c%d/Run%d_Baselines_Good%d.png",station_num,config,runNum,1);
		// cAvg->SaveAs(save_temp_title);
		// delete cAvg;

		bool baselinesGood = areBaselinesGood(average_copy,station_num,runNum);
		// bool baselinesGood=true;
		for(int i=0; i<16; i++) delete average_copy[i];
		printf("Baselines run %d are good %d? (%s)\n", runNum, baselinesGood, argv[file]);

		if(baselinesGood && foundMegaAverage==false){
			for(int i=0; i<16; i++){
				mega_average.push_back((TGraph*)average[i]->Clone());
			}
			numInAverage++;
			foundMegaAverage=true;
		}
		else if(baselinesGood && foundMegaAverage==true){
			for(int i=0; i<16; i++){
				int numSamps = mega_average[i]->GetN();
				for(int samp=0; samp<numSamps; samp++){
					mega_average[i]->GetY()[samp]+=average[i]->GetY()[samp];
				}
			}
			numInAverage++;
		}

		SummaryFile->Close();
	}
	cout<<"num in average is "<<numInAverage<<endl;

	for(int i=0; i<16; i++){
		int numSamps= mega_average[i]->GetN();
		for(int samp=0; samp<numSamps; samp++){
			mega_average[i]->GetY()[samp]/=double(numInAverage);
		}
	}

	TCanvas *cAvg = new TCanvas("","",4*1100,4*850);
	cAvg->Divide(4,4);
	for(int i=0; i<16; i++){
		cAvg->cd(i+1);
		mega_average[i]->Draw("AL");
		mega_average[i]->GetYaxis()->SetRangeUser(15,40);
		mega_average[i]->GetYaxis()->SetTitle("Power (dB)");
		mega_average[i]->GetXaxis()->SetTitle("Frequency (MHz)");
	}
	sprintf(save_temp_title,"%s/trouble_events/MegaAverage_station%d_config%d.png",plotPath,station_num,config);
	cAvg->SaveAs(save_temp_title);
	delete cAvg;

	char outputFileName[500];
	sprintf(outputFileName,"%s/average_baseline_station%d_config%d.root",plotPath,station_num,config);
	TFile *OutputFile = TFile::Open(outputFileName, "RECREATE");
	TTree* AverageBaselineTree=new TTree("AverageBaselineTree", "AverageBaselineTree");
	vector<TGraph*> average_to_save;
	average_to_save.resize(16);
	stringstream ss2;
	for(int i=0; i<16; i++){
		ss2.str(""); ss2<<"baselines_average_RF_chan_"<<i;
		AverageBaselineTree->Branch(ss2.str().c_str(), &average_to_save[i]);
	}
	for(int i=0; i<16; i++){
		average_to_save[i]=(TGraph*)mega_average[i]->Clone();
	}
	AverageBaselineTree->Fill();
	OutputFile->Write();
	OutputFile->Close();
	delete OutputFile;

}//end main