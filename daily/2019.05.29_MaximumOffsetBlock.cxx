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
#include "AraQualCuts.h"

using namespace std;

void PlotThisEvent(int station, int runNum, int event);

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
		cout<< "Usage\n" << argv[0] << " <station> <config> <offset block file>"<<endl;
		return 0;
	}
	int station = atoi(argv[1]);
	int config = atoi(argv[2]);

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

	TDatime start(2013, 01, 01, 00, 00,0);
	int start_bin = start.Convert();
	TDatime stop(2013, 12, 31, 24, 00,0);
	int stop_bin = stop.Convert();
	TH2D *distro_time[16];
	for(int i=0; i<16; i++){
		distro_time[i] = new TH2D("","",365, start_bin, stop_bin, 200,0,200);
		distro_time[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		distro_time[i]->GetXaxis()->SetTimeFormat("%m");
		distro_time[i]->GetXaxis()->SetNdivisions(12,0,0,false);
	}

	TH1D *MaximumRollingMean[16];
	TH1D *MaximumRollingMeanCal[16];
	TH1D *MaximumRollingS4Glitch[16];
	TH1D *stragglers[16];
	for(int i=0; i<16; i++){
		stringstream ss1;
		ss1<<"Ch"<<i;
		MaximumRollingMean[i] = new TH1D(ss1.str().c_str(),"",401,0,401);
		ss1.str("");
		ss1<<"ChCal"<<i;
		MaximumRollingMeanCal[i] = new TH1D(ss1.str().c_str(),"",401,0,401);
		ss1.str("");
		ss1<<"ChGlitch"<<i;
		MaximumRollingS4Glitch[i] = new TH1D(ss1.str().c_str(),"",401,0,401);
		ss1.str("");
		ss1<<"ChStragglers"<<i;
		stragglers[i] = new TH1D(ss1.str().c_str(),"",401,0,401);
	}
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
		int runNum;
		bool hasDigitizerError;
		double first_block_mean[16];
		double maxMeans[16];
		double minMeans[16];
		double absMaxBlockMeans[16];
		double power_75[16];
		double power_110[16];
		double power_150[16];
		double above_850[16];
		double between[16];
		double RMS[16];
		int unixTime;
		bool hasA3S4Issue;
		inTree->SetBranchAddress("isCal", &isCal);
		inTree->SetBranchAddress("isSoft", &isSoft);
		inTree->SetBranchAddress("isShort", &isShort);
		inTree->SetBranchAddress("run",&runNum);
		inTree->SetBranchAddress("hasDigitizerError",&hasDigitizerError);
		inTree->SetBranchAddress("unixTime",&unixTime);
		inTree->SetBranchAddress("first_block_mean",&first_block_mean);
		inTree->SetBranchAddress("maxMeans",&maxMeans);
		inTree->SetBranchAddress("minMeans",&minMeans);
		inTree->SetBranchAddress("hasA3S4Issue", &hasA3S4Issue);
		inTree->SetBranchAddress("absMaxBlockMeans",&absMaxBlockMeans);
		//inTree->Branch("power_75", &power_75);
		//inTree->Branch("power_110", &power_110);
		//inTree->Branch("power_150", &power_150);
		//inTree->Branch("above_850", &above_850);
		//inTree->Branch("between", &between);
		//inTree->SetBranchAddress("RMS",&RMS);

		int numEntries = inTree->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;
		inTree->GetEvent(0);

		//now to loop over events
		for(int event=0; event<numEntries; event++){
			inTree->GetEvent(event);
			if(runNum==1862 && event==1){
				printf("Run %d Event %d Digitizer %d and Software %d \n", runNum,event,hasDigitizerError,isSoft);
			}
			if(hasDigitizerError || isSoft) continue;
			// if(hasDigitizerError || isSoft || isShort || hasA3S4Issue) continue;
			bool TroubleEvent=false;
	
			// config 1
			if(runNum==1795 && event==5617) TroubleEvent=true;
			if(runNum==1795 && event==5617) TroubleEvent=true;
			if(runNum==1805 && event==70) TroubleEvent=true;
			if(runNum==1806 && event==41) TroubleEvent=true;
			if(runNum==1806 && event==41) TroubleEvent=true;
			if(runNum==1806 && event==99) TroubleEvent=true;
			if(runNum==1810 && event==52) TroubleEvent=true;
			if(runNum==1810 && event==144) TroubleEvent=true;
			if(runNum==1810 && event==144) TroubleEvent=true;
			// more subtle glitch
			// if(runNum==1862 && event==1) TroubleEvent=true;
			// if(runNum==1862 && event==1) TroubleEvent=true;

			// config 2
			// if(runNum==1075 && event==2953) TroubleEvent=true;
			// if(runNum==476 && event==223) TroubleEvent=true;
			// if(runNum==476 && event==223) TroubleEvent=true;
			// if(runNum==476 && event==506) TroubleEvent=true;
			// if(runNum==476 && event==726) TroubleEvent=true;
			// if(runNum==553 && event==1) TroubleEvent=true;
			// if(runNum==704 && event==4855) TroubleEvent=true;
			// if(runNum==711 && event==0) TroubleEvent=true;
			// if(runNum==711 && event==0) TroubleEvent=true;
			// if(runNum==808 && event==9238) TroubleEvent=true;

			// config 3
			if(runNum==3419 && event==13502) TroubleEvent=true;
			if(runNum==3419 && event==13502) TroubleEvent=true;
			// if(runNum==3419 && event==13676) TroubleEvent=true;
			if(runNum==3421 && event==81) TroubleEvent=true;
			if(runNum==3421 && event==81) TroubleEvent=true;
			if(runNum==3422 && event==49) TroubleEvent=true;
			if(runNum==3423 && event==19) TroubleEvent=true;
			if(runNum==3423 && event==30) TroubleEvent=true;
			if(runNum==3424 && event==107) TroubleEvent=true;
			if(runNum==3424 && event==107) TroubleEvent=true;
			// if(runNum==3426 && event==6) TroubleEvent=true;
			if(runNum==3437 && event==10823) TroubleEvent=true;
			// if(runNum==3681 && event==0) TroubleEvent=true;
			// if(runNum==3841 && event==5) TroubleEvent=true;
			// if(runNum==3979 && event==1) TroubleEvent=true;
			// if(runNum==3979 && event==1) TroubleEvent=true;
			// if(runNum==3983 && event==6) TroubleEvent=true;
			// if(runNum==3983 && event==6) TroubleEvent=true;
			// if(runNum==3990 && event==7) TroubleEvent=true;
			// if(runNum==3994 && event==4) TroubleEvent=true;
			// if(runNum==3994 && event==4) TroubleEvent=true;
			// if(runNum==4099 && event==1) TroubleEvent=true;
			// if(runNum==4101 && event==0) TroubleEvent=true;
			// if(runNum==4101 && event==0) TroubleEvent=true;
			// if(runNum==4102 && event==1) TroubleEvent=true;
			// if(runNum==4129 && event==1) TroubleEvent=true;
			// if(runNum==4129 && event==1) TroubleEvent=true;
			// if(runNum==4140 && event==6) TroubleEvent=true;
			// if(runNum==4140 && event==6) TroubleEvent=true;
			// if(runNum==4141 && event==0) TroubleEvent=true;
			// if(runNum==4150 && event==4) TroubleEvent=true;
			// if(runNum==4150 && event==4) TroubleEvent=true;
			// if(runNum==4167 && event==1) TroubleEvent=true;
			// if(runNum==4175 && event==1) TroubleEvent=true;
			// if(runNum==4175 && event==1) TroubleEvent=true;
			// if(runNum==4214 && event==1) TroubleEvent=true;
			// if(runNum==4271 && event==0) TroubleEvent=true;
			// if(runNum==4274 && event==1) TroubleEvent=true;
			if(runNum==4921 && event==76) TroubleEvent=true;
			if(runNum==4921 && event==76) TroubleEvent=true;
			if(runNum==4926 && event==78) TroubleEvent=true;
			if(runNum==4926 && event==78) TroubleEvent=true;
			if(runNum==4929 && event==97) TroubleEvent=true;
			if(runNum==4931 && event==104) TroubleEvent=true;
			if(runNum==4931 && event==104) TroubleEvent=true;
			if(runNum==4932 && event==69) TroubleEvent=true;
			if(runNum==4932 && event==69) TroubleEvent=true;
			if(runNum==4933 && event==76) TroubleEvent=true;
			if(runNum==4933 && event==89) TroubleEvent=true;
			if(runNum==4933 && event==89) TroubleEvent=true;
			if(runNum==4938 && event==84) TroubleEvent=true;
			if(runNum==4938 && event==84) TroubleEvent=true;
			if(runNum==4939 && event==0) TroubleEvent=true;
			if(runNum==4939 && event==0) TroubleEvent=true;
			if(runNum==4947 && event==13) TroubleEvent=true;
			if(runNum==4947 && event==13) TroubleEvent=true;
			if(runNum==4948 && event==0) TroubleEvent=true;
			if(runNum==4948 && event==0) TroubleEvent=true;
			if(runNum==4948 && event==10) TroubleEvent=true;
			if(runNum==4948 && event==10) TroubleEvent=true;
			if(runNum==4951 && event==14) TroubleEvent=true;
			if(runNum==4951 && event==14) TroubleEvent=true;
			if(runNum==4956 && event==6) TroubleEvent=true;
			if(runNum==4957 && event==84) TroubleEvent=true;
			if(runNum==4957 && event==84) TroubleEvent=true;
			if(runNum==7760 && event==99) TroubleEvent=true;

			// config 4
			if(runNum==6684 && event==1) TroubleEvent=true;

			// config 5
			if(runNum==2001 && event==0) TroubleEvent=true;
			if(runNum==2466 && event==2232) TroubleEvent=true;
			if(runNum==2466 && event==2232) TroubleEvent=true;
			if(runNum==2472 && event==41) TroubleEvent=true;

			num_total++;
			double thisAbs;
			bool printThisEvent=false;
			for(int i=0; i<16; i++){
				if(abs(maxMeans[i])>abs(minMeans[i]))
					thisAbs=abs(maxMeans[i]);
				else if(abs(minMeans[i])>abs(maxMeans[i]))
					thisAbs=abs(minMeans[i]);
				// thisAbs=abs(maxMeans[i]-minMeans[i]);
				// thisAbs=RMS[i];
				// thisAbs=absMaxBlockMeans[i];
				if(thisAbs>399.) thisAbs=399.;
				MaximumRollingMean[i]->Fill(thisAbs);
				if(thisAbs>180 && i==4) printThisEvent=true;
				// distro_time[i]->Fill(unixTime,thisAbs);
				if(isCal)
					MaximumRollingMeanCal[i]->Fill(thisAbs);
				if(hasA3S4Issue){
					MaximumRollingS4Glitch[i]->Fill(thisAbs);
					// printThisEvent=true;
				}
				if(TroubleEvent){
					stragglers[i]->Fill(thisAbs);
					printf("Trouble event! This abs is %.2f \n", thisAbs);
				}
			}
			// printThisEvent=false;
			if(printThisEvent){
				PlotThisEvent(station,runNum,event);
			}
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
	TCanvas *c = new TCanvas("","",3*1100,3*850);
	c->Divide(4,4);
	for(int i=0; i<16; i++){
		c->cd(i+1);
		MaximumRollingMean[i]->Draw("");
		MaximumRollingMean[i]->SetLineWidth(4);
		MaximumRollingMean[i]->SetLineColor(kBlue);
		MaximumRollingMeanCal[i]->Draw("same");
		MaximumRollingMeanCal[i]->SetLineColor(kGreen);
		MaximumRollingMeanCal[i]->SetLineWidth(3);
		// MaximumRollingS4Glitch[i]->Draw("same");
		// MaximumRollingS4Glitch[i]->SetLineColor(kCyan);
		// MaximumRollingS4Glitch[i]->SetLineWidth(3);
		stragglers[i]->Draw("same");
		stragglers[i]->SetLineColor(kMagenta);
		stragglers[i]->SetLineWidth(3);
		gPad->SetLogy();
		// MaximumRollingMean[i]->GetYaxis()->SetRangeUser(1.,1e8);
		// MaximumRollingMean[i]->GetYaxis()->SetRangeUser(0.00001,200.);		
		MaximumRollingMean[i]->GetYaxis()->SetTitle("Number of Events");
		// MaximumRollingMean[i]->GetXaxis()->SetTitle("RMS");
		MaximumRollingMean[i]->GetXaxis()->SetTitle("Abs Max Rolling Mean");
		MaximumRollingMean[i]->SetTitle(titlesForGraphs[i].c_str());
		MaximumRollingMean[i]->GetXaxis()->SetTitleOffset(1.1);
		MaximumRollingMean[i]->GetYaxis()->SetTitleOffset(1.1);
		MaximumRollingMean[i]->GetZaxis()->SetTitleOffset(1.1);
		MaximumRollingMean[i]->GetXaxis()->SetTitleSize(0.06);
		MaximumRollingMean[i]->GetYaxis()->SetTitleSize(0.06);
		MaximumRollingMean[i]->GetZaxis()->SetTitleSize(0.06);
		MaximumRollingMean[i]->GetXaxis()->SetLabelSize(0.06);
		MaximumRollingMean[i]->GetYaxis()->SetLabelSize(0.06);
		MaximumRollingMean[i]->GetZaxis()->SetLabelSize(0.06);
		// MaximumRollingMean[i]->GetXaxis()->SetRangeUser(10,40);		
	}

	// char fit_equation[150];
	// sprintf(fit_equation,"gaus");
	// TF1 *fit = new TF1("fit",fit_equation,20,30);
	// MaximumRollingMeanCal[0]->Fit("fit","R");
	// printf("Chi-Square/NDF %.2f / %.2f \n",fit->GetChisquare(),double(fit->GetNDF()));

	char save_plot_title[400];
	// sprintf(save_plot_title,"%s/glitch_detect/%d.%d.%d_Distro_MaxBlockMeanAverageOverBlocksOnly_A%d_c%d_%dEvents.png",plotPath,year_now,month_now,day_now,station,config,num_total);
	sprintf(save_plot_title,"%s/glitch_detect/%d.%d.%d_Distro_AbsMaxRollingMean_A%d_c%d_%dEvents.png",plotPath,year_now,month_now,day_now,station,config,num_total);
	// sprintf(save_plot_title,"%s/glitch_detect/%d.%d.%d_Distro_RMS_A%d_c%d_%dEvents.png",plotPath,year_now,month_now,day_now,station,config,num_total);
	c->SaveAs(save_plot_title);
	delete c;

	// TCanvas *c2 = new TCanvas("","",4*1100,4*850);
	// c2->Divide(4,4);
	// for(int i=0; i<16; i++){
	// 	c2->cd(i+1);
	// 	distro_time[i]->Draw("colz");
	// 	gPad->SetLogz();
	// 	// distro_time[i]->GetZaxis()->SetRangeUser(1.,1e4);
	// 	// distro[i]->GetYaxis()->SetRangeUser(0.00001,200.);		
	// 	distro_time[i]->GetYaxis()->SetTitle("Largest Rolling Mean");
	// 	distro_time[i]->GetXaxis()->SetTitle("Time");
	// 	distro_time[i]->GetZaxis()->SetTitle("Number of Events");
	// 	distro_time[i]->SetTitle(titlesForGraphs[i].c_str());
	// 	distro_time[i]->GetXaxis()->SetTitleOffset(1.1);
	// 	distro_time[i]->GetYaxis()->SetTitleOffset(1.1);
	// 	// distro_time[i]->GetZaxis()->SetTitleOffset(1.1);
	// 	gPad->SetRightMargin(0.15);
	// 	distro_time[i]->GetXaxis()->SetTitleSize(0.06);
	// 	distro_time[i]->GetYaxis()->SetTitleSize(0.06);
	// 	distro_time[i]->GetZaxis()->SetTitleSize(0.06);
	// 	distro_time[i]->GetXaxis()->SetLabelSize(0.06);
	// 	distro_time[i]->GetYaxis()->SetLabelSize(0.06);
	// 	distro_time[i]->GetZaxis()->SetLabelSize(0.06);
	// }
	// sprintf(save_plot_title,"%s/glitch_detect/%d.%d.%d_Distro_MaxBlockMeanAverageOverBlocksOnly_A%d_%d_%dEvents_2D.png",plotPath,year_now,month_now,day_now,station,config,num_total);
	// c2->SaveAs(save_plot_title);
	// delete c2;
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

		bool print_graphs=true;
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
		}

		delete realAtriEvPtr;
	}
	mapFile->Close();
	delete mapFile;
}
