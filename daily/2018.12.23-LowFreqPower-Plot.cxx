////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  2018.12.23-LowFreqPower-Plot.cxx 
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
#include "TCanvas.h"
#include "TStyle.h"

#include "tools_PlottingFns.h"
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"

using namespace std;

int PlotThisEvent(int station, int year, int runNum, int event);

int main(int argc, char **argv)
{
	gStyle->SetOptStat(110011);
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	stringstream ss;
	
	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <station> <year> <low freq power file>"<<endl;
		return 0;
	}
	int station = atoi(argv[1]);
	int year = atoi(argv[2]);

	TH1D *distro[16];
	TH1D *distro_75[16];
	TH1D *glitch[16];
	for(int i=0; i<16; i++){
		distro[i] = new TH1D("","",100,0,1); 
		distro_75[i] = new TH1D("","",100,0,1); 
		glitch[i] = new TH1D("","",100,0,1);
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
		int waveform_length[16];
		double frac_power_75[16];
		double frac_power_110[16];
		double frac_power_150[16];
		int runNum;
		bool hasDigitizerError;
		inTree->SetBranchAddress("isCal", &isCal);
		inTree->SetBranchAddress("isSoft", &isSoft);
		inTree->SetBranchAddress("waveform_length", &waveform_length);
		inTree->SetBranchAddress("frac_power_75", &frac_power_75);
		inTree->SetBranchAddress("frac_power_110", &frac_power_110);
		inTree->SetBranchAddress("frac_power_150", &frac_power_150);
		inTree->SetBranchAddress("run",&runNum);
		inTree->SetBranchAddress("hasDigitizerError",&hasDigitizerError);

		int numEntries = inTree->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;

		//now to loop over events
		for(int event=0; event<numEntries; event++){
			inTree->GetEvent(event);

			if(isCal || isSoft || hasDigitizerError) continue;
			num_total++;
			bool need_to_print=false;
			int this_glitch_counter=0;
			for(int i=0; i<16; i++){
				distro[i]->Fill(frac_power_150[i]);
				distro_75[i]->Fill(frac_power_75[i]);
				if(frac_power_75[i]>0.1 && i!=11){
					// if(station==3 && i!=3 && i!=7 && i!=11 && i!=15) this_glitch_counter++;
					glitch_number[i]++;
					this_glitch_counter++;
				}
			}
			if(this_glitch_counter>2){
				need_to_print=false;
				num_total_glitch++;
			}
			if(need_to_print) PlotThisEvent(station, year, runNum ,event);
		}
		fpIn->Close();
		delete fpIn;
	} //end loop over input files

	for(int i=0; i<16; i++){
		printf("Chan %d fraction: %.6f \n", i, double(glitch_number[i])/double(num_total));
	}

	stringstream ss1;
	vector<string> titlesForGraphs;
	for (int i = 0; i < 16; i++){
		ss1.str("");
		ss1 << "Channel " << i;
		titlesForGraphs.push_back(ss1.str());
	}

	printf("Total number cut is %d/%d = %.6f \n", num_total_glitch, num_total, double(num_total_glitch)/double(num_total));
	TH1 *cumulative[16];
	for(int i=0; i<16; i++){
		// distro[i]->Scale(1./double(distro[i]->GetEntries())*100.);
		cumulative[i] = distro[i]->GetCumulative();
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
		distro[i]->GetYaxis()->SetRangeUser(1.,1e8);
		// distro[i]->GetYaxis()->SetRangeUser(0.00001,200.);		
		distro[i]->GetYaxis()->SetTitle("Number of Events");
		distro[i]->GetXaxis()->SetTitle("Fraction of Power Below X MHz");
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

		distro_75[i]->Draw("same");
		distro_75[i]->SetLineWidth(5);
		distro_75[i]->SetLineColor(kRed);
	}
	char save_plot_title[400];
	sprintf(save_plot_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/%d.%d.%d_Distro_Glitch_A%d_%d_%dEvents.png",year_now,month_now,day_now,station,year,num_total);
	c->SaveAs(save_plot_title);
	delete c;
}

int PlotThisEvent(int station, int year, int runNum, int event){
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

	int unixTime = (int)rawPtr->unixTime;
	int unixTimeUs =(int)rawPtr->unixTimeUs;
	printf("Unixtime is %d \n", unixTime);
	printf("Unixtime microsecond is %d \n", unixTimeUs);

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
	mapFile->Close();
	delete mapFile;
	return 0;
}