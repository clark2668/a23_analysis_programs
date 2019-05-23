////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	v2_print_event.cxx 
////	Print waveforms, spectra, and maps for a single event
////
////	Nov 20187
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
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"

//AraRoot includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "Settings.h"
#include "Detector.h"
#include "Report.h"
#include "RayTraceCorrelator.h"
AraAntPol::AraAntPol_t Vpol = AraAntPol::kVertical;
AraAntPol::AraAntPol_t Hpol = AraAntPol::kHorizontal;
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_Cuts.h"

using namespace std;

int main(int argc, char **argv)
{

	if(argc<5){
		cout<< "Usage\n" << argv[0] << " <station> <year> <run num> <event>"<<endl;;
		return -1;
	}
	int station = atoi(argv[1]);
	int year = atoi(argv[2]);
	int runNum = atoi(argv[3]);
	int event = atoi(argv[4]);
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;
	char *DataDirPath(getenv("DATA_DIR"));
	if (DataDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;

	gStyle->SetOptStat(11);

	char run_file_name[400];
	sprintf(run_file_name,"%s/RawData/A%d/%d/sym_links/event%d.root",DataDirPath,station,year,runNum,runNum);
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
		sprintf(ped_file_name,"%s/run_specific_peds/A%d/%d/event%d_specificPeds.dat",PedDirPath,station,year,runNum);
	}
	else if(year==2014 || year==2015 || year==2016){
		sprintf(ped_file_name,"%s/run_specific_peds/A%d/%d/event00%d_specificPeds.dat",PedDirPath,station,year,runNum);
	}
	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	calibrator->setAtriPedFile(ped_file_name,stationID); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist
	
	UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawPtr,AraCalType::kLatestCalib);

	int unixTime = (int)rawPtr->unixTime;
	int unixTimeUs =(int)rawPtr->unixTimeUs;
	int timeStamp = (int)rawPtr->timeStamp;
	printf("Unixtime is %d \n", unixTime);
	printf("Unixtime microsecond is %d \n", unixTimeUs);
	printf("timeStamp is %d \n", timeStamp);
	printf("Event Number is %d \n", realAtriEvPtr->eventNumber);

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
	for(int i=0; i<16; i++){
		// printf("Chan %d has %d samples\n", i, waveforms[i]->GetN());
		// for(int j=0; j<realAtriEvPtr->blockVec.size(); j++){
		// 	printf("Block %d is %d \n", j, realAtriEvPtr->blockVec[j].getBlock());
		// }

	}

	vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(waveforms, 0.5, xLabel, yLabel, titlesForGraphs);
	vector<TGraph*> grWaveformsPadded = makePaddedGraphs(grWaveformsInt, 0, xLabel, yLabel, titlesForGraphs);
	xLabel = "Frequency (Hz)"; yLabel = "Power Spectral Density (mV/Hz)";
	vector<TGraph*> grWaveformsPowerSpectrum = makePowerSpectrumGraphs(grWaveformsPadded, xLabel, yLabel, titlesForGraphs);

	gStyle->SetOptStat(0);
	bool do_reco=false;
	if(do_reco){
		//set up the ray tracer
		Settings *settings = new Settings();
		string setupfile = "setup.txt";
		settings->ReadFile(setupfile);
		cout << "Read " << setupfile << " file!" << endl;
		settings->NOFZ=1;
		Detector *detector=0;
		RayTraceCorrelator *theCorrelators[2];
		theCorrelators[0] =  new RayTraceCorrelator(station, 41., settings, 1, 4); //41 m, cal puser
		theCorrelators[1] =  new RayTraceCorrelator(station, 300., settings, 1, 4);//300 m, far reco

		TH2D *map_41m_V;
		TH2D *map_300m_V;
		TH2D *map_41m_H;
		TH2D *map_300m_H;
		TH2D *map_41m_V_select;

		map_41m_V = theCorrelators[0]->getInterferometricMap_RT(settings, detector, realAtriEvPtr, Vpol, 0, 0);
		map_300m_V = theCorrelators[1]->getInterferometricMap_RT(settings, detector, realAtriEvPtr, Vpol, 0, 0);
		map_41m_H = theCorrelators[0]->getInterferometricMap_RT(settings, detector, realAtriEvPtr, Hpol, 0, 0);
		map_300m_H = theCorrelators[1]->getInterferometricMap_RT(settings, detector, realAtriEvPtr, Hpol, 0, 0);

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

		TCanvas *cMaps = new TCanvas("","",2*1100,2*850);
		cMaps->Divide(2,2);
			cMaps->cd(3);
			map_41m_V->Draw("colz");
			cMaps->cd(4);
			map_41m_H->Draw("colz");
			cMaps->cd(1);
			map_300m_V->Draw("colz");
			cMaps->cd(2);
			map_300m_H->Draw("colz");
		char save_temp_title[400];		
		sprintf(save_temp_title,"%s/single_events/%d.%d.%d_Run%d_Ev%d_Maps.png",plotPath,year_now,month_now,day_now,runNum,event);
		cMaps->SaveAs(save_temp_title);
		delete cMaps;
		delete map_41m_V; delete map_300m_V; delete map_41m_H; delete map_300m_H; 
		// delete map_41m_V_select;
	}

	bool do_reco_snrweighted=false;
	if(do_reco_snrweighted){
		//set up the ray tracer
		Settings *settings = new Settings();
		string setupfile = "setup.txt";
		settings->ReadFile(setupfile);
		cout << "Read " << setupfile << " file!" << endl;
		settings->NOFZ=1;
		Detector *detector=0;
		RayTraceCorrelator *theCorrelators[2];
		theCorrelators[0] =  new RayTraceCorrelator(station, 41., settings, 1, 4); //41 m, cal puser
		theCorrelators[1] =  new RayTraceCorrelator(station, 300., settings, 1, 4);//300 m, far reco

		TH2D *map_41m_V;
		TH2D *map_300m_V;
		TH2D *map_41m_H;
		TH2D *map_300m_H;
		TH2D *map_41m_V_select;

		map_41m_V = theCorrelators[0]->getInterferometricMap_RT_SNRweighted(settings, detector, realAtriEvPtr, Vpol, 0, 0);
		map_300m_V = theCorrelators[1]->getInterferometricMap_RT_SNRweighted(settings, detector, realAtriEvPtr, Vpol, 0, 0);
		map_41m_H = theCorrelators[0]->getInterferometricMap_RT_SNRweighted(settings, detector, realAtriEvPtr, Hpol, 0, 0);
		map_300m_H = theCorrelators[1]->getInterferometricMap_RT_SNRweighted(settings, detector, realAtriEvPtr, Hpol, 0, 0);

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

		TCanvas *cMaps = new TCanvas("","",2*1100,2*850);
		cMaps->Divide(2,2);
			cMaps->cd(3);
			map_41m_V->Draw("colz");
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
		delete map_41m_V; delete map_300m_V; delete map_41m_H; delete map_300m_H; 
		// delete map_41m_V_select;
	}

	vector<TGraph*> dummy;
	for(int i=0; i<16; i++){
		vector<double> thisX;
		vector<double> thisY;
		thisX.push_back(-200.);
		thisX.push_back(400.);
		thisY.push_back(-100.);
		thisY.push_back(100.);
		dummy.push_back(new TGraph(thisX.size(), &thisX[0], &thisY[0]));
		dummy[i]->GetXaxis()->SetTitle("Time (ns)");
		dummy[i]->GetYaxis()->SetTitle("Voltage (mV)");
		dummy[i]->GetXaxis()->SetLabelSize(0.07);
		dummy[i]->GetYaxis()->SetLabelSize(0.07);
		dummy[i]->GetXaxis()->SetTitleSize(0.07);
		dummy[i]->GetYaxis()->SetTitleSize(0.07);
	}

	char save_temp_title[300];
	sprintf(save_temp_title,"%s/single_events/%d.%d.%d_Run%d_Ev%d_Waveforms.png",plotPath,year_now,month_now,day_now,runNum,event);
	TCanvas *cWave = new TCanvas("","",4*1100,4*850);
	cWave->Divide(4,4);
	for(int i=0; i<16; i++){
		cWave->cd(i+1);
		dummy[i]->Draw("AP");
		waveforms[i]->Draw("Lsame");
		waveforms[i]->SetLineWidth(2);
	}
	cWave->SaveAs(save_temp_title);
	delete cWave;

	sprintf(save_temp_title,"%s/single_events/%d.%d.%d_Run%d_Ev%d_Spectra.png",plotPath,year_now,month_now,day_now,runNum,event);
	TCanvas *cSpec = new TCanvas("","",4*1100,4*850);
	cSpec->Divide(4,4);
	for(int i=0; i<16; i++){
		cSpec->cd(i+1);
		grWaveformsPowerSpectrum[i]->Draw("AL");
		grWaveformsPowerSpectrum[i]->SetLineWidth(3);
		gPad->SetLogy();
		grWaveformsPowerSpectrum[i]->GetYaxis()->SetRangeUser(10.,1e7);
	}
	cSpec->SaveAs(save_temp_title);
	delete cSpec;

	vector <TGraph*> more_spec;
	for(int i=0; i<16; i++){
		more_spec.push_back(FFTtools::makePowerSpectrumMilliVoltsNanoSecondsdB(grWaveformsPadded[i]));
		more_spec[i]->GetXaxis()->SetTitle("Frequency");
		more_spec[i]->GetYaxis()->SetTitle("Power (dB)");
		more_spec[i]->SetLineWidth(3);
		more_spec[i]->GetYaxis()->SetRangeUser(0,80);
		more_spec[i]->SetTitle("");
	}
	TCanvas *cnow = new TCanvas("","",3*1100,850);
	cnow->Divide(3,1);
	for(int i=1; i<4; i++){
		cnow->cd(i);
		more_spec[i]->Draw("AL");
	}
	sprintf(save_temp_title,"%s/single_events/%d.%d.%d_Run%d_Ev%d_Spectra_Special.png",plotPath,year_now,month_now,day_now,runNum,event);
	cnow->SaveAs(save_temp_title);
	delete cnow;

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