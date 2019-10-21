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

#include "tools_PlottingFns.h"
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraQualCuts.h"
#include "tools_Cuts.h"

using namespace std;

void PlotThisEvent(int station, int runNum, int event);

int main(int argc, char **argv)
{
	// gStyle->SetOptStat(110011);
	gStyle->SetOptStat(0);
	gStyle->SetTimeOffset(0);
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	stringstream ss;
	
	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <station> <config> <offset block file>"<<endl;
		return 0;
	}
	int station = atoi(argv[1]);
	int config = atoi(argv[2]);

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

	TH1D *all_events = new TH1D("all_events","all_events",8000,0,8000);
	TH1D *soft_events = new TH1D("soft_events","soft_events",8000,0,8000);
	TH1D *fraction = new TH1D("fraction", "fraction", 8000, 0, 8000);

	vector<int> BadRunList=BuildBadRunList(station);
	int numBadEventsTotal=0;
	int num_total=0;

	for(int file_num=3; file_num<argc; file_num++){

		cout<<argv[file_num]<<endl;

		TFile *fpIn = TFile::Open(argv[file_num]);
		if(!fpIn) {
			std::cout << "Can't open file\n";
			return -1;
		}
		stringstream ss;
		ss << "outTree";
		TTree *inTree = (TTree*) fpIn->Get(ss.str().c_str());
		if(!inTree){
			cout<<"Can't open filter tree"<<endl;
			return -1;
		}
		bool isCal;
		bool isSoft;
		bool isShort;
		double spareChannelMaxSamp[4];
		double spareChannelMinSamp[4];
		double spareChannelRMS[4];
		int runNum;
		int unixTime;
		int unixTimeUs;
		int timeStamp;
		bool hasDigitizerError=false;
		int eventNumber;
		inTree->SetBranchAddress("isCal", &isCal);
		inTree->SetBranchAddress("isSoft", &isSoft);
		inTree->SetBranchAddress("isShort", &isShort);
		inTree->SetBranchAddress("run",&runNum);
		inTree->SetBranchAddress("hasDigitizerError",&hasDigitizerError);
		inTree->SetBranchAddress("unixTime",&unixTime);
		inTree->SetBranchAddress("spareChannelMaxSamp",&spareChannelMaxSamp);
		inTree->SetBranchAddress("spareChannelMinSamp",&spareChannelMinSamp);
		inTree->SetBranchAddress("spareChannelRMS",&spareChannelRMS);
		inTree->SetBranchAddress("run",&runNum);
		inTree->SetBranchAddress("eventNumber",&eventNumber);

		int numEntries = inTree->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;
		inTree->GetEvent(0);

		// TH1D *event_type = new TH1D("","",2300,0,2300);
		int num_local=0;
		for(int event=0; event<numEntries; event++){
			num_total++;
			num_local++;
			inTree->GetEvent(event);
			all_events->Fill(runNum);
			if(isSoft){
				soft_events->Fill(runNum);
			}

			// if(isSoft){
			// 	event_type->Fill(event,2);
			// }
			// else{
			// 	event_type->Fill(event,1);
			// }
		}

		// TCanvas *cLocal = new TCanvas("","",2*850,850);
		// 	event_type->Draw("");
		// 	event_type->GetXaxis()->SetTitle("10pct iterator event index");
		// 	event_type->GetYaxis()->SetTitle("Is Software Trigger Yes/No 2/1");
		// 	event_type->SetLineWidth(2);
		// char save_plot_title[500];
		// sprintf(save_plot_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/%d.%d.%d_SoftwareTrigList_A%d_c%d_Run%d_%dEvents.png",year_now,month_now,day_now,station,config,runNum,num_local);
		// cLocal->SaveAs(save_plot_title);
		// delete cLocal;

		fpIn->Close();
		delete fpIn;

	} //end loop over input files
	for(int bin=1; bin<all_events->GetNbinsX(); bin++){
		int total = all_events->GetBinContent(bin);
		if(total>0){
			int soft  = soft_events->GetBinContent(bin);
			double this_frac=double(soft)/double(total);
			fraction->SetBinContent(bin,this_frac);
			if(bin>1750 && bin<1850){
				printf("Bin %d, Run %d, Fraction is %.2f\n", bin, bin-1, this_frac);
			}
		}
		else{
			fraction->SetBinContent(bin,-0.1);
		}
	}
	TCanvas *c = new TCanvas("","",2*850,850);
	fraction->Draw("");
	fraction->GetXaxis()->SetTitle("Run Number");
	fraction->GetYaxis()->SetTitle("Fraction of Software Triggers");
	TLine *line;
	if(config==1){
		fraction->GetXaxis()->SetRangeUser(1400,2000);
		line = new TLine(1400,0,2000,0);
	}
	if(config==2){
		fraction->GetXaxis()->SetRangeUser(450,1500);
		line = new TLine(500,0,1500,0);
	}
	if(config==5){
		fraction->GetXaxis()->SetRangeUser(450,1500);
		line = new TLine(500,0,1500,0);
	}
	fraction->GetYaxis()->SetRangeUser(-0.2,1);
	line->Draw("sameL");
	line->SetLineColor(kBlack);
	line->SetLineStyle(9);
	char save_plot_title[400];
	sprintf(save_plot_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/%d.%d.%d_Distro_FracSoft_A%d_c%d_%dEvents.png",year_now,month_now,day_now,station,config,num_total);
	c->SaveAs(save_plot_title);
	delete c;
}

void PlotThisEvent(int station, int runNum, int event){
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *DataDirPath(getenv("DATA_DIR"));
	if (DataDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;

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
	// printf("Unixtime is %d \n", unixTime);
	// printf("Unixtime microsecond is %d \n", unixTimeUs);
	// printf("TimeStamp is %d \n", timeStamp);

	bool hasA3S4Issue=false;

	bool do_print=true;
	if(do_print){

		int stationID = rawPtr->stationId;
		char ped_file_name[400];
		sprintf(ped_file_name,"%s/run_specific_peds/A%d/all_peds/event%d_specificPeds.dat",PedDirPath,station,runNum);
		AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
		calibrator->setAtriPedFile(ped_file_name,stationID); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist
	
		UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawPtr,AraCalType::kLatestCalib);
		// AraQualCuts *qualCut = AraQualCuts::Instance();

		bool print_graphs=false;
		if(print_graphs){
			stringstream ss1;
			string xLabel, yLabel;
			xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
			vector<string> titlesForGraphs;
			for (int i = 0; i < 16; i++){
				ss1.str("");
				ss1 << "Channel " << i;
				titlesForGraphs.push_back(ss1.str());
			}
			vector <TGraph*> waveforms = makeGraphsFromRF(realAtriEvPtr,16,xLabel,yLabel,titlesForGraphs);

			char save_temp_title[300];
			sprintf(save_temp_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/graphs/%d.%d.%d_Station%d_Run%d_Ev%d_Waveforms.png",year_now,month_now,day_now,station,runNum,event);
			TCanvas *cWave = new TCanvas("","",4*1100,4*850);
			cWave->Divide(4,4);
			for(int i=0; i<16; i++){
				cWave->cd(i+1);
				waveforms[i]->Draw("AL");
				waveforms[i]->SetLineWidth(3);
			}
			cWave->SaveAs(save_temp_title);
			delete cWave;	

			bool print_spectrum=false;
			if(print_spectrum){
				vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(waveforms, 0.6, xLabel, yLabel, titlesForGraphs);
				vector<TGraph*> grWaveformsPadded = makePaddedGraphs(grWaveformsInt, 0, xLabel, yLabel, titlesForGraphs);
				xLabel = "Frequency (Hz)"; yLabel = "Power Spectral Density (mV/Hz)";
				vector<TGraph*> grWaveformsPowerSpectrum = makePowerSpectrumGraphs(grWaveformsPadded, xLabel, yLabel, titlesForGraphs);
				sprintf(save_temp_title,"/users/PAS0654/osu0673/A23_analysis/results/debug/%d.%d.%d_Run%d_Ev%d_Spectra.png",year_now,month_now,day_now,runNum,event);
				TCanvas *cSpec = new TCanvas("","",4*1100,4*850);
				cSpec->Divide(4,4);
				for(int i=0; i<16; i++){
					cSpec->cd(i+1);
					grWaveformsPowerSpectrum[i]->Draw("AL");
					grWaveformsPowerSpectrum[i]->SetLineWidth(3);
					gPad->SetLogy();
				}
				cSpec->SaveAs(save_temp_title);
				delete cSpec;
				for(int i=0; i<16; i++){
					delete grWaveformsInt[i];
					delete grWaveformsPadded[i];
					delete grWaveformsPowerSpectrum[i];
				}
			}

			for(int i=0; i<16; i++){
				delete waveforms[i];
			}
		}

		bool print_new_graphs=true;
		if(print_new_graphs){
			vector<string> titlesForGraphs;
			titlesForGraphs.clear();
			stringstream ss1;
			for (int i = 0; i < 32; i++){
				ss1.str("");
				ss1 << "Channel " << i;
				titlesForGraphs.push_back(ss1.str());
			}

			vector<TGraph*> waveformsAll;
			for(int i=0; i<32; i++){
				waveformsAll.push_back(realAtriEvPtr->getGraphFromElecChan(i));
				// printf("RMS of Elec Chan %d is %3.2f\n",i,waveformsAll[i]->GetRMS(2));
				waveformsAll[i]->SetTitle(titlesForGraphs[i].c_str());
			}
			char save_temp_title[400];
			sprintf(save_temp_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/graphs/%d.%d.%d_Run%d_Ev%d_AllWaveforms.png",year_now,month_now,day_now,runNum,event);
			TCanvas *cAllWave = new TCanvas("","",8*1100,4*850);
			cAllWave->Divide(8,4);
			for(int i=0; i<32; i++){
				cAllWave->cd(i+1);
				waveformsAll[i]->Draw("ALsame");
				waveformsAll[i]->SetLineWidth(2);
				// waveforms[i]->SetMarkerStyle(kFullCircle);
				// waveforms[i]->SetMarkerSize(2);
			}
			cAllWave->SaveAs(save_temp_title);
			delete cAllWave;
		}

		delete realAtriEvPtr;

	}
	mapFile->Close();
	delete mapFile;
}
