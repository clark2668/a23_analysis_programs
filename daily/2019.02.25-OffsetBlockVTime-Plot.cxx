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

#include "tools_PlottingFns.h"
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"

using namespace std;

double PlotThisEvent(int station, int year, int runNum, int event);

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
	
	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <station> <year> <offset block file>"<<endl;
		return 0;
	}
	int station = atoi(argv[1]);
	int year = atoi(argv[2]);

	TDatime start(year, 01, 01, 00, 00,0);
	int start_bin = start.Convert();
	TDatime stop(year, 12, 31, 24, 00,0);
	int stop_bin = stop.Convert();
	TDatime start_long(2012, 12, 31, 24, 00,00);
	int start_bin_long = start_long.Convert();
	TDatime stop_long(2016, 12, 31, 24, 00,00);
	int stop_bin_long = stop_long.Convert();

	TH1D *stamp_timing = new TH1D("timeStamp","timeStamp",110,-100,1e8);

	TH1D *distro[16];
	TH2D *distro_time[16];
	TH2D *distro_time_long[16];
	for(int i=0; i<16; i++){
		distro[i] = new TH1D("","",100,-100,100);
		distro_time[i] = new TH2D("","",365, start_bin, stop_bin, 100,-100,100);
		distro_time[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		distro_time[i]->GetXaxis()->SetTimeFormat("%m");
		distro_time[i]->GetXaxis()->SetNdivisions(12,0,0,false);

		distro_time_long[i] = new TH2D("","",1460, start_bin_long, stop_bin_long, 100,-100,100);
		distro_time_long[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		distro_time_long[i]->GetXaxis()->SetTimeFormat("%b'%y");
		distro_time_long[i]->GetXaxis()->SetNdivisions(8,6,0,false);
	}

	int num_total=0;
	int num_total_glitch=0;

	int glitch_number[16]={0};

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
		int runNum;
		bool hasDigitizerError;
		double first_block_mean[16];
		int unixTime;

		inTree->SetBranchAddress("isCal", &isCal);
		inTree->SetBranchAddress("isSoft", &isSoft);
		inTree->SetBranchAddress("run",&runNum);
		inTree->SetBranchAddress("hasDigitizerError",&hasDigitizerError);
		inTree->SetBranchAddress("unixTime",&unixTime);
		inTree->SetBranchAddress("first_block_mean",&first_block_mean);

		int numEntries = inTree->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;
		inTree->GetEvent(0);

		//now to loop over events
		for(int event=0; event<numEntries; event++){
			inTree->GetEvent(event);
			if(isCal) PlotThisEvent(station, year, runNum, event);

			// if(isCal || isSoft || hasDigitizerError) continue;
			// // if(runNum>=1308 && runNum<=1407) continue;
			// num_total++;
			// for(int i=0; i<16; i++){
			// 	distro[i]->Fill(first_block_mean[i]);
			// 	distro_time[i]->Fill(unixTime,first_block_mean[i]);
			// 	distro_time_long[i]->Fill(unixTime,first_block_mean[i]);
			// }
			// if(first_block_mean[0]<-30){
			// 	cout<<"First block mean is "<<first_block_mean[0]<<endl;
			// 	double stamp_timping_this = PlotThisEvent(station, year, runNum, event);
			// 	// stamp_timing->Fill((double)stamp_timping_this);
			// }

		}
		fpIn->Close();
		delete fpIn;
	} //end loop over input files

	stringstream ss1;
	vector<string> titlesForGraphs;
	for (int i = 0; i < 16; i++){
		ss1.str("");
		ss1 << "Channel " << i;
		titlesForGraphs.push_back(ss1.str());
	}

	gStyle->SetOptStat(0);
	TCanvas *c = new TCanvas("","",4*1100,4*850);
	c->Divide(4,4);
	for(int i=0; i<16; i++){
		c->cd(i+1);
		distro[i]->Draw("");
		distro[i]->SetLineWidth(5);
		distro[i]->SetLineColor(kBlue);
		gPad->SetLogy();
		// distro[i]->GetYaxis()->SetRangeUser(1.,1e8);
		// distro[i]->GetYaxis()->SetRangeUser(0.00001,200.);		
		distro[i]->GetYaxis()->SetTitle("Number of Events");
		distro[i]->GetXaxis()->SetTitle("Mean Voltage of First Block (mV)");
		distro[i]->SetTitle(titlesForGraphs[i].c_str());
		distro[i]->GetXaxis()->SetTitleOffset(1.1);
		distro[i]->GetYaxis()->SetTitleOffset(1.1);
		distro[i]->GetZaxis()->SetTitleOffset(1.1);
		distro[i]->GetXaxis()->SetTitleSize(0.06);
		distro[i]->GetYaxis()->SetTitleSize(0.06);
		distro[i]->GetZaxis()->SetTitleSize(0.06);
		distro[i]->GetXaxis()->SetLabelSize(0.06);
		distro[i]->GetYaxis()->SetLabelSize(0.06);
		distro[i]->GetZaxis()->SetLabelSize(0.06);
	}
	char save_plot_title[400];
	sprintf(save_plot_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/%d.%d.%d_Distro_FirstBlockMean_A%d_%d_%dEvents.png",year_now,month_now,day_now,station,year,num_total);
	c->SaveAs(save_plot_title);
	delete c;

	TCanvas *c2 = new TCanvas("","",4*1100,4*850);
	c2->Divide(4,4);
	for(int i=0; i<16; i++){
		c2->cd(i+1);
		distro_time[i]->Draw("colz");
		gPad->SetLogz();
		distro_time[i]->GetZaxis()->SetRangeUser(1.,1e4);
		// distro[i]->GetYaxis()->SetRangeUser(0.00001,200.);		
		distro_time[i]->GetYaxis()->SetTitle("Mean Voltage of First Block (mV)");
		distro_time[i]->GetXaxis()->SetTitle("Time");
		distro_time[i]->GetZaxis()->SetTitle("Number of Events");
		distro_time[i]->SetTitle(titlesForGraphs[i].c_str());
		distro_time[i]->GetXaxis()->SetTitleOffset(1.1);
		distro_time[i]->GetYaxis()->SetTitleOffset(1.1);
		// distro_time[i]->GetZaxis()->SetTitleOffset(1.1);
		gPad->SetRightMargin(0.15);
		distro_time[i]->GetXaxis()->SetTitleSize(0.06);
		distro_time[i]->GetYaxis()->SetTitleSize(0.06);
		distro_time[i]->GetZaxis()->SetTitleSize(0.06);
		distro_time[i]->GetXaxis()->SetLabelSize(0.06);
		distro_time[i]->GetYaxis()->SetLabelSize(0.06);
		distro_time[i]->GetZaxis()->SetLabelSize(0.06);
	}
	sprintf(save_plot_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/%d.%d.%d_Distro_FirstBlockMean_A%d_%d_%dEvents_2D.png",year_now,month_now,day_now,station,year,num_total);
	c2->SaveAs(save_plot_title);
	delete c2;


	TCanvas *c3 = new TCanvas("","",8*1100,4*850);
	c3->Divide(4,4);
	for(int i=0; i<16; i++){
		c3->cd(i+1);
		distro_time_long[i]->Draw("colz");
		gPad->SetLogz();
		distro_time_long[i]->GetZaxis()->SetRangeUser(1.,1e5);
		// distro[i]->GetYaxis()->SetRangeUser(0.00001,200.);		
		distro_time_long[i]->GetYaxis()->SetTitle("Mean Voltage of First Block (mV)");
		distro_time_long[i]->GetXaxis()->SetTitle("Time");
		distro_time_long[i]->GetZaxis()->SetTitle("Number of Events");
		distro_time_long[i]->SetTitle(titlesForGraphs[i].c_str());
		distro_time_long[i]->GetXaxis()->SetTitleOffset(1.1);
		distro_time_long[i]->GetYaxis()->SetTitleOffset(1.1);
		// distro_time[i]->GetZaxis()->SetTitleOffset(1.1);
		gPad->SetRightMargin(0.15);
		distro_time_long[i]->GetXaxis()->SetTitleSize(0.06);
		distro_time_long[i]->GetYaxis()->SetTitleSize(0.06);
		distro_time_long[i]->GetZaxis()->SetTitleSize(0.06);
		distro_time_long[i]->GetXaxis()->SetLabelSize(0.06);
		distro_time_long[i]->GetYaxis()->SetLabelSize(0.06);
		distro_time_long[i]->GetZaxis()->SetLabelSize(0.06);
	}
	sprintf(save_plot_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/%d.%d.%d_Distro_FirstBlockMean_MultiYear_A%d_%dEvents_MuliYear.png",year_now,month_now,day_now,station,num_total);
	c3->SaveAs(save_plot_title);
	delete c3;

	TCanvas *c4 = new TCanvas("","",1100,850);
	stamp_timing->Draw();
	stamp_timing->SetLineWidth(2);
	stamp_timing->SetLineColor(kBlack);
	stamp_timing->GetYaxis()->SetTitle("Number of Events");
	stamp_timing->GetXaxis()->SetTitle("Clock Cylces");
	// gPad->SetLogy();
	// stamp_timing->GetYaxis()->SetRangeUser(1e0,3e6);
	// stamp_timing->GetXaxis()->SetRangeUser(-1e2,1e5);
	sprintf(save_plot_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/%d.%d.%d_UntaggedCalTiming.png",year_now,month_now,day_now);
	c4->SaveAs(save_plot_title);
	delete c4;

}

double PlotThisEvent(int station, int year, int runNum, int event){
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char run_file_name[400];
	if(year==2013){
		sprintf(run_file_name,"/fs/scratch/PAS0654/ara/10pct/RawData/A%d/%d/run%d/event%d.root",station,year,runNum,runNum);
	}
	else if(year==2014 || year==2015 || year==2016){
		sprintf(run_file_name,"/fs/scratch/PAS0654/ara/10pct/RawData/A%d/%d/sym_links/event00%d.root",station,year,runNum,runNum);
	}
	TFile *mapFile = TFile::Open(run_file_name);
	if(!mapFile){
		cout<<"Can't open data file for map!"<<endl;
		return -1;
	}
	TTree *eventTree = (TTree*) mapFile-> Get("eventTree");
	if(!eventTree){
		cout<<"Can't find eventTree for map"<<endl;
		return -1;
	}

	RawAtriStationEvent *rawPtr =0;
	eventTree->SetBranchAddress("event",&rawPtr);
	eventTree->GetEvent(event);

	int unixTime = (int)rawPtr->unixTime;
	int unixTimeUs =(int)rawPtr->unixTimeUs;
	int timeStamp = (int)rawPtr->timeStamp;
	printf("Unixtime is %d \n", unixTime);
	printf("Unixtime microsecond is %d \n", unixTimeUs);
	printf("TimeStamp is %d \n", timeStamp);

	bool do_print=true;
	if(do_print){

		int stationID = rawPtr->stationId;
		char ped_file_name[400];
		if(year==2013){
			sprintf(ped_file_name,"/fs/scratch/PAS0654/ara/peds/run_specific_peds/A%d/%d/event%d_specificPeds.dat",station,year,runNum);
		}
		else if(year==2014 || year==2015 || year==2016){
			sprintf(ped_file_name,"/fs/scratch/PAS0654/ara/peds/run_specific_peds/A%d/%d/event00%d_specificPeds.dat",station,year,runNum);
		}
		AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
		calibrator->setAtriPedFile(ped_file_name,stationID); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist
	
		UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawPtr,AraCalType::kLatestCalib);

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
		sprintf(save_temp_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/%d.%d.%d_Station%d_Run%d_Ev%d_Waveforms.png",year_now,month_now,day_now,station,runNum,event);
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
		delete realAtriEvPtr;
	}
	mapFile->Close();
	delete mapFile;
	return double(timeStamp);
}