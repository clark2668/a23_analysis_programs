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
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TChain.h"
#include "TTimeStamp.h"
#include "TMath.h"
#include "TF1.h"

// AraRoot includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "FFTtools.h"
#include "AraGeomTool.h"
#include "AraQualCuts.h"

// analysis custom
#include "tools_Cuts.h"
#include "tools_Stats.h"
#include "tools_CommandLine.h"
#include "tools_outputObjects.h"
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"

using namespace std;
int PlotThisEvent(int station, int config, int runNum, int event, int problempol);

int main(int argc, char **argv){

	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	stringstream ss;
	gStyle->SetOptStat(0);

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;


	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <1-station> <2-config>"<<endl;;
		return -1;
	}
	int station = atoi(argv[1]);
	int config = atoi(argv[2]);

	if(station!=2 && station!=3){
		printf("No good! You asked for station %d, but this code only works for stations 2 and 3 \n",station);
		return -1;
	}

	char fileHundredName[400];
	sprintf(fileHundredName,"/users/PAS0654/osu0673/A23_analysis_new2/results/unblind/background_fit/A%d_c%d_sample100.root",station,config);
	TFile* fileHundred = new TFile(fileHundredName,"READ");
	TH1D *h1Hundred[2];
	h1Hundred[0] = (TH1D*)fileHundred->Get("DiffDistroV");
	h1Hundred[1] = (TH1D*)fileHundred->Get("DiffDistroH");

	TH2D *h2Hundred_SNRvsCorr[2];
	h2Hundred_SNRvsCorr[0] = (TH2D*)fileHundred->Get("2DDistroV");
	h2Hundred_SNRvsCorr[1] = (TH2D*)fileHundred->Get("2DDistroH");
	h2Hundred_SNRvsCorr[0]->SetTitle("2D Distro V - 100%");
	h2Hundred_SNRvsCorr[1]->SetTitle("2D Distro H - 100%");

	char fileHundredRootName[400];
	sprintf(fileHundredRootName,"/users/PAS0654/osu0673/A23_analysis_new2/results/unblind/background_fit/2d_cut_values_A%d_c%d_sample100.root",station,config);
	TFile *outFile = TFile::Open(fileHundredRootName,"READ");
	TTree *outTree = (TTree*) outFile->Get("outTree");
	int hist_this_pol[2];
	double corr_val_out[2];
	double snr_val_out[2];
	int runNum_out;
	int eventNum_out;
	outTree->SetBranchAddress("passes_this_pol_V",&hist_this_pol[0]);
	outTree->SetBranchAddress("passes_this_pol_H",&hist_this_pol[1]);
	outTree->SetBranchAddress("corr_val_V",&corr_val_out[0]);
	outTree->SetBranchAddress("corr_val_H",&corr_val_out[1]);
	outTree->SetBranchAddress("snr_val_V",&snr_val_out[0]);
	outTree->SetBranchAddress("snr_val_H",&snr_val_out[1]);
	outTree->SetBranchAddress("runNum_out",&runNum_out);
	outTree->SetBranchAddress("eventNum_out",&eventNum_out);
	int numEntries = outTree->GetEntries();

	TH2D *redo[2];
	redo[0] = new TH2D("vredo","vredo",100,0,0.05,300,0,30);
	redo[1] = new TH2D("hredo","hredo",100,0,0.05,300,0,30);

	for(int i=0; i<numEntries; i++){
		outTree->GetEntry(i);
		if(hist_this_pol[0]){
			redo[0]->Fill(corr_val_out[0],snr_val_out[0]);
		}
		if(hist_this_pol[1]){
			redo[1]->Fill(corr_val_out[1],snr_val_out[1]);
			// if(config==2 && snr_val_out[1]>8){
				// printf("Run %d, Event %d \n", runNum_out, eventNum_out);
				// PlotThisEvent(station, config, runNum_out, eventNum_out, 1);
			// }
			if(config==4 && corr_val_out[1]>0.009){
				printf("on i %d, Run %d, Event %d \n", i, runNum_out, eventNum_out);
				// PlotThisEvent(station, config, runNum_out, eventNum_out, 1);
			}
		}
	}

	vector<TGraph*> cut_lines;
	for(int pol=0; pol<2; pol++){
		vector <double> x_vals_for_line;
		vector <double> y_vals_for_line;
		double rcut_slope;
		double rcut_intercept;
		getRCutValues(station, config, pol, rcut_slope, rcut_intercept);
		for(double x=0; x<0.020; x+=0.00001){
			double y_val = (rcut_slope * x ) + rcut_intercept;
			x_vals_for_line.push_back(x);
			y_vals_for_line.push_back(y_val);
		}
		cut_lines.push_back(new TGraph(x_vals_for_line.size(), &x_vals_for_line[0], &y_vals_for_line[0]));
	}

	TCanvas *c2 = new TCanvas("","",2*850,850);
	c2->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		c2->cd(pol+1);
		h2Hundred_SNRvsCorr[pol]->Draw("colz");
		h2Hundred_SNRvsCorr[pol]->GetXaxis()->SetTitle("Correlation Value");
		h2Hundred_SNRvsCorr[pol]->GetYaxis()->SetTitle("SNR");
		cut_lines[pol]->Draw("same");
		cut_lines[pol]->SetLineColor(kRed);
		gPad->SetLogz();
	}
	// for(int pol=0; pol<2; pol++){
	// 	c2->cd(pol+3);
	// 	redo[pol]->Draw("colz)");
	// 	gPad->SetLogz();
	// }
	char saveTitle[400];
	sprintf(saveTitle,"/users/PAS0654/osu0673/A23_analysis_new2/results/unblind/background_fit/A%d_c%d_Hundred2DDistro.png",station,config);
	c2->SaveAs(saveTitle);

	outFile->Close();
	fileHundred->Close();
}

int PlotThisEvent(int station, int config, int runNum, int event, int problempol){
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *DataDirPath(getenv("DATA_DIR_100"));
	if (DataDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

	char run_file_name[400];
	sprintf(run_file_name,"%s/RawData/A%d/all_runs/event%d.root",DataDirPath,station,runNum);
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
	sprintf(ped_file_name,"%s/run_specific_peds/A%d/all_peds/event%d_specificPeds.dat",PedDirPath,station,runNum);
	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	calibrator->setAtriPedFile(ped_file_name,stationID); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist

	AraQualCuts *qualCut = AraQualCuts::Instance();
	UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawPtr,AraCalType::kLatestCalib);
	printf("Run %d, Event %d \n", runNum, realAtriEvPtr->eventNumber);
	printf("	Is Quality Event? %d \n", qualCut->isGoodEvent(realAtriEvPtr));

	int unixTime = (int)rawPtr->unixTime;
	int unixTimeUs =(int)rawPtr->unixTimeUs;
	int timeStamp = (int)rawPtr->timeStamp;
	printf("	Unixtime is %d \n", unixTime);
	printf("	Unixtime microsecond is %d \n", unixTimeUs);
	printf("	timeStamp is %d \n", timeStamp);

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
	vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(waveforms, 0.5, xLabel, yLabel, titlesForGraphs);
	vector<TGraph*> grWaveformsPadded = makePaddedGraphs(grWaveformsInt, 0, xLabel, yLabel, titlesForGraphs);
	xLabel = "Frequency (Hz)"; yLabel = "Power Spectral Density (mV/Hz)";
	vector<TGraph*> grWaveformsPowerSpectrum = makePowerSpectrumGraphs(grWaveformsPadded, xLabel, yLabel, titlesForGraphs);

	char save_temp_title[300];
	sprintf(save_temp_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/unblind/background_fit/peculiar_events/%d.%d.%d_Run%d_Ev%d_ProblemPol%d_Waveforms.png",plotPath,year_now,month_now,day_now,runNum,event,problempol);
	TCanvas *cWave = new TCanvas("","",4*1100,4*850);
	cWave->Divide(4,4);
	for(int i=0; i<16; i++){
		cWave->cd(i+1);
		// dummy[i]->Draw("AL");
		// dummy[i]->SetLineColor(kWhite);
		// dummy[i]->GetXaxis()->SetRangeUser(300.,500.);

		waveforms[i]->Draw("AL");
		waveforms[i]->SetLineWidth(3);
		// waveforms[i]->GetXaxis()->SetRangeUser(300.,500.);
	}
	cWave->SaveAs(save_temp_title);
	delete cWave;

	// vector<TGraph*> electChansGraphs;
	// electChansGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(6));
	// electChansGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(14));
	// electChansGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(22));
	// electChansGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(30));

	// TCanvas *cWave_spare = new TCanvas("","",4*1100,850);
	// cWave_spare->Divide(4,1);
	// for(int i=0; i<4; i++){
	// 	cWave_spare->cd(i+1);
	// 	electChansGraphs[i]->Draw("AL");
	// 	waveforms[i]->SetLineWidth(3);
	// }
	// sprintf(save_temp_title,"%s/trouble_events/sideregion_%d.%d.%d_Run%d_Ev%d_ProblemPol%d_Waveforms_SpareChans.png",plotPath,year_now,month_now,day_now,runNum,event,problempol);
	// cWave_spare->SaveAs(save_temp_title);
	// delete cWave_spare;


	sprintf(save_temp_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/unblind/background_fit/peculiar_events/%d.%d.%d_Run%d_Ev%d_ProblemPol%d_Spectra.png",plotPath,year_now,month_now,day_now,runNum,event,problempol);
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
		delete waveforms[i];
		delete grWaveformsInt[i];
		delete grWaveformsPadded[i];
		delete grWaveformsPowerSpectrum[i];
	}
	delete realAtriEvPtr;
	mapFile->Close();
	delete mapFile;
	return 0;

}