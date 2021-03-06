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
	TDatime stop(2016, 12, 31, 24, 00,0);
	int stop_bin = stop.Convert();
	// TDatime start_long(2012, 12, 31, 24, 00,00);
	// int start_bin_long = start_long.Convert();
	// TDatime stop_long(2016, 12, 31, 24, 00,00);
	// int stop_bin_long = stop_long.Convert();

	// TH1D *distro[16];
	TH2D *distro_time[4];
	// TH2D *distro_time_long[16];
	for(int i=0; i<4; i++){
		// distro[i] = new TH1D("","",100,-100,100);
		distro_time[i] = new TH2D("","",365*4, start_bin, stop_bin, 100,0,100);
		distro_time[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		// distro_time[i]->GetXaxis()->SetTimeFormat("%m");
		// distro_time[i]->GetXaxis()->SetNdivisions(12,0,0,false);

		distro_time[i]->GetXaxis()->SetTimeFormat("%b'%y");
		distro_time[i]->GetXaxis()->SetNdivisions(8,6,0,false);

		// distro_time_long[i] = new TH2D("","",1460, start_bin_long, stop_bin_long, 100,-100,100);
		// distro_time_long[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		// distro_time_long[i]->GetXaxis()->SetTimeFormat("%b'%y");
		// distro_time_long[i]->GetXaxis()->SetNdivisions(8,6,0,false);
	}

	TH2D *distro_2dcut = new TH2D("","",50,0,100,5,-0.5,4.5);
	TGraph *TroubleEventOVerlay = new TGraph();
	int num_in_overlay=0;

	TH1D *hMeanOverStdDev[4];
	TH1D *hMeanOverStdDev_cal[4];
	TH1D *hMeanOverStdDev_stragglers[4];
	for(int i=0; i<4; i++){
		stringstream ss1;
		ss1<<"Ch"<<i;
		hMeanOverStdDev[i] = new TH1D(ss1.str().c_str(),"",100,0,100);
		ss1.str("");
		ss1<<"ChCal"<<i;
		hMeanOverStdDev_cal[i] = new TH1D(ss1.str().c_str(),"",100,0,100);
		ss1.str("");
		ss1<<"ChStragglers"<<i;
		hMeanOverStdDev_stragglers[i] = new TH1D(ss1.str().c_str(),"",401,0,401);
	}
	int num_total=0;

	vector<int> BadRunList=BuildBadRunList(station);
	int numBadEventsTotal=0;

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
		double spareChannelMaxMeanOverStdDev[4];
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
		inTree->SetBranchAddress("spareChannelMaxMeanOverStdDev",&spareChannelMaxMeanOverStdDev);
		inTree->SetBranchAddress("run",&runNum);
		inTree->SetBranchAddress("eventNumber",&eventNumber);

		int numEntries = inTree->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;
		inTree->GetEvent(0);

		//now to loop over events
		for(int event=0; event<numEntries; event++){
			// cout<<"On event"<<endl;
			inTree->GetEvent(event);
			// if(isSoft || eventNumber<5 || runNum==480 || runNum==479 || runNum==485 || runNum==486 || runNum==484 || isBadRun(station,runNum,BadRunList)) continue;
			if(isSoft || eventNumber<5) continue;
			// if(isSoft || isShort ) continue;
			// if(hasDigitizerError || isSoft || isShort || hasA3S4Issue) continue;
			bool TroubleEvent=false;
			// config 1
			// if(runNum==1795 && event==5617) TroubleEvent=true;
			// if(runNum==1795 && event==5617) TroubleEvent=true;
			// if(runNum==1805 && event==70) TroubleEvent=true;
			// if(runNum==1806 && event==41) TroubleEvent=true;
			// if(runNum==1806 && event==41) TroubleEvent=true;
			// if(runNum==1806 && event==99) TroubleEvent=true;
			// if(runNum==1810 && event==52) TroubleEvent=true;
			// if(runNum==1810 && event==144) TroubleEvent=true;
			// if(runNum==1810 && event==144) TroubleEvent=true;
			// if(runNum==1862 && event==1) TroubleEvent=true;
			// if(runNum==1862 && event==1) TroubleEvent=true;

			// // config 2
			// // if(runNum==1075 && event==2953) TroubleEvent=true;
			// if(runNum==476 && event==223) TroubleEvent=true;
			// if(runNum==476 && event==506) TroubleEvent=true;
			// if(runNum==476 && event==726) TroubleEvent=true;
			// if(runNum==553 && event==1) TroubleEvent=true;
			// if(runNum==711 && event==0) TroubleEvent=true;

			// // config 3
			// if(runNum==3419 && event==13502) TroubleEvent=true;
			// if(runNum==3419 && event==13502) TroubleEvent=true;
			// if(runNum==3419 && event==13676) TroubleEvent=true;
			// if(runNum==3421 && event==81) TroubleEvent=true;
			// if(runNum==3421 && event==81) TroubleEvent=true;
			// if(runNum==3422 && event==49) TroubleEvent=true;
			// if(runNum==3423 && event==19) TroubleEvent=true;
			// if(runNum==3423 && event==30) TroubleEvent=true;
			// if(runNum==3424 && event==107) TroubleEvent=true;
			// if(runNum==3424 && event==107) TroubleEvent=true;
			// if(runNum==3426 && event==6) TroubleEvent=true;
			// // if(runNum==3437 && event==10823) TroubleEvent=true;
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
			// if(runNum==4921 && event==76) TroubleEvent=true;
			// if(runNum==4921 && event==76) TroubleEvent=true;
			// if(runNum==4926 && event==78) TroubleEvent=true;
			// if(runNum==4926 && event==78) TroubleEvent=true;
			// if(runNum==4929 && event==97) TroubleEvent=true;
			// if(runNum==4931 && event==104) TroubleEvent=true;
			// if(runNum==4931 && event==104) TroubleEvent=true;
			// if(runNum==4932 && event==69) TroubleEvent=true;
			// if(runNum==4932 && event==69) TroubleEvent=true;
			// if(runNum==4933 && event==76) TroubleEvent=true;
			// if(runNum==4933 && event==89) TroubleEvent=true;
			// if(runNum==4933 && event==89) TroubleEvent=true;
			// if(runNum==4938 && event==84) TroubleEvent=true;
			// if(runNum==4938 && event==84) TroubleEvent=true;
			// if(runNum==4939 && event==0) TroubleEvent=true;
			// if(runNum==4939 && event==0) TroubleEvent=true;
			// if(runNum==4947 && event==13) TroubleEvent=true;
			// if(runNum==4947 && event==13) TroubleEvent=true;
			// if(runNum==4948 && event==0) TroubleEvent=true;
			// if(runNum==4948 && event==0) TroubleEvent=true;
			// if(runNum==4948 && event==10) TroubleEvent=true;
			// if(runNum==4948 && event==10) TroubleEvent=true;
			// if(runNum==4951 && event==14) TroubleEvent=true;
			// if(runNum==4951 && event==14) TroubleEvent=true;
			// if(runNum==4956 && event==6) TroubleEvent=true;
			// if(runNum==4957 && event==84) TroubleEvent=true;
			// if(runNum==4957 && event==84) TroubleEvent=true;
			// if(runNum==7760 && event==99) TroubleEvent=true;

			// // config 4
			// if(runNum==6684 && event==1) TroubleEvent=true;

			// // config 5
			// if(runNum==2001 && event==0) TroubleEvent=true;
			// if(runNum==2466 && event==2232) TroubleEvent=true;
			// if(runNum==2466 && event==2232) TroubleEvent=true;
			// if(runNum==2472 && event==41) TroubleEvent=true;
			if(runNum==3663 && event==6){
				TroubleEvent=true;
			}

			// what about A2
			num_total++;
			bool printThisEvent=false;
			int numBad=0;
			int numReallBad=0;

			for(double level=0.; level<100.; level+=2.){
				int numBad=0;
				for(int i=0; i<4; i++){
					if(spareChannelMaxMeanOverStdDev[i]>level)
						numBad++;
				}
				distro_2dcut->Fill(level,numBad);
			}

			if(TroubleEvent){
				// TroubleEventOverlay->SetPoint(numTroubleFound,);
				num_in_overlay++;
			}

			for(int i=0; i<4; i++){
				double thisRMS = spareChannelMaxMeanOverStdDev[i];
				hMeanOverStdDev[i]->Fill(thisRMS);
				distro_time[i]->Fill(unixTime, thisRMS);
				// if(thisRMS>99.) thisRMS=99.;
				if(isCal)
					hMeanOverStdDev_cal[i]->Fill(thisRMS);
				if(TroubleEvent){
					hMeanOverStdDev_stragglers[i]->Fill(thisRMS);
					// printf("Trouble event! This RMS is %.2f \n", thisRMS);
				}
				if(thisRMS>50)
					numBad++;
				// printf("Channel %d RMS is %.2f \n", i, thisRMS);
			}
			if(numBad>1)
				numBadEventsTotal++;
			printThisEvent=false;
			if(printThisEvent){
				PlotThisEvent(station,runNum,event);
			}
		}
		fpIn->Close();
		delete fpIn;
	} //end loop over input files

	printf("Num bad total is %d out of %d = %2.5f\n", numBadEventsTotal, num_total, double(numBadEventsTotal)/double(num_total));

	stringstream ss1;
	vector<string> titlesForGraphs;
	ss1.str("");
	ss1<<"Elec Channel "<<6;
	titlesForGraphs.push_back(ss1.str());
	ss1.str("");
	ss1<<"Elec Channel "<<14;
	titlesForGraphs.push_back(ss1.str());
	ss1.str("");
	ss1<<"Elec Channel "<<22;
	titlesForGraphs.push_back(ss1.str());
	ss1.str("");
	ss1<<"Elec Channel "<<30;
	titlesForGraphs.push_back(ss1.str());

	gStyle->SetOptStat(0);
	TCanvas *c = new TCanvas("","",4*850,850);
	c->Divide(4,1);
	for(int i=0; i<4; i++){
		c->cd(i+1);
		// hRMS[i]->Scale(1./hRMS[i]->GetMaximum());
		// hRMS_cal[i]->Scale(1./hRMS_cal[i]->GetMaximum());
		// hRMS_stragglers[i]->Scale(1./hRMS_stragglers[i]->GetMaximum());
		hMeanOverStdDev[i]->Draw("");
		hMeanOverStdDev[i]->SetLineWidth(4);
		hMeanOverStdDev[i]->SetLineColor(kBlue);
		hMeanOverStdDev_cal[i]->Draw("same");
		hMeanOverStdDev_cal[i]->SetLineColor(kGreen);
		hMeanOverStdDev_cal[i]->SetLineWidth(4);
		hMeanOverStdDev_stragglers[i]->Draw("same");
		hMeanOverStdDev_stragglers[i]->SetLineColor(kRed);
		hMeanOverStdDev_stragglers[i]->SetLineWidth(4);
		gPad->SetLogy();
		hMeanOverStdDev[i]->GetYaxis()->SetRangeUser(0.1,1e8);
		// hRMS[i]->GetYaxis()->SetRangeUser(hRMS[i]->GetMinimum(), hRMS[i]->GetMaximum());
		hMeanOverStdDev[i]->GetYaxis()->SetTitle("Number of Events");
		hMeanOverStdDev[i]->GetXaxis()->SetTitle("RMS");
		hMeanOverStdDev[i]->SetTitle(titlesForGraphs[i].c_str());
		hMeanOverStdDev[i]->GetXaxis()->SetTitleOffset(1.1);
		hMeanOverStdDev[i]->GetYaxis()->SetTitleOffset(1.1);
		hMeanOverStdDev[i]->GetZaxis()->SetTitleOffset(1.1);
		hMeanOverStdDev[i]->GetXaxis()->SetTitleSize(0.06);
		hMeanOverStdDev[i]->GetYaxis()->SetTitleSize(0.06);
		hMeanOverStdDev[i]->GetZaxis()->SetTitleSize(0.06);
		hMeanOverStdDev[i]->GetXaxis()->SetLabelSize(0.06);
		hMeanOverStdDev[i]->GetYaxis()->SetLabelSize(0.06);
		hMeanOverStdDev[i]->GetZaxis()->SetLabelSize(0.06);
		// hRMS[i]->GetXaxis()->SetRangeUser(10,40);
		if(i+1==1){
			TLegend *leg = new TLegend(0.48,0.6,0.9,0.9);
			leg->AddEntry(hMeanOverStdDev[i],"All Events","l");
			leg->AddEntry(hMeanOverStdDev_cal[i],"Tagged Cal Pulsers","l");
			// leg->AddEntry(hRMS_stragglers[i],"Straggling Events","l");
			leg->Draw();
		}
		double total = hMeanOverStdDev[i]->Integral();
		double above = hMeanOverStdDev[i]->Integral(50,100);
		printf("Number above 50 is %d/%d = %2.5f\n",int(above),int(total),double(above/total));
	}
	char save_plot_title[400];
	sprintf(save_plot_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/%d.%d.%d_Distro_SpareChannelMeanOverStdDev_A%d_c%d_%dEvents.png",year_now,month_now,day_now,station,config,num_total);
	c->SaveAs(save_plot_title);
	delete c;

	TCanvas *c2 = new TCanvas("","",4*850,2*850);
	c2->Divide(2,2);
	for(int i=0; i<4; i++){
		c2->cd(i+1);
		distro_time[i]->Draw("colz");
		gPad->SetLogz();
		distro_time[i]->GetZaxis()->SetRangeUser(1.,1e4);
		// distro[i]->GetYaxis()->SetRangeUser(0.00001,200.);		
		distro_time[i]->GetYaxis()->SetTitle("RMS");
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
	sprintf(save_plot_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/%d.%d.%d_Distro_SpareChannelMeanOverStdDevvsTime_A%d_c%d_%dEvents.png",year_now,month_now,day_now,station,config,num_total);
	c2->SaveAs(save_plot_title);
	delete c2;
	
	beautify_TH2D();
	distro_2dcut->Scale(1./distro_2dcut->GetMaximum());
	TCanvas *c3 = new TCanvas("","",1100,850);
		distro_2dcut->Draw("colz");
		gPad->SetLogz();
		distro_2dcut->GetZaxis()->SetRangeUser(1e-6,1);
		distro_2dcut->GetZaxis()->SetRangeUser(distro_2dcut->GetMinimum(),distro_2dcut->GetMaximum());
		distro_2dcut->GetYaxis()->SetTitle("Number of Channels");
		distro_2dcut->GetXaxis()->SetTitle("#mu/(#sigma/#sqrt(N))");
		distro_2dcut->GetZaxis()->SetTitle("Fraction of Events");
		distro_2dcut->GetXaxis()->SetTitleOffset(1.1);
		distro_2dcut->GetYaxis()->SetTitleOffset(1.1);
		// distro_time[i]->GetZaxis()->SetTitleOffset(1.1);
		gPad->SetRightMargin(0.15);
		distro_2dcut->GetXaxis()->SetTitleSize(0.04);
		distro_2dcut->GetYaxis()->SetTitleSize(0.04);
		distro_2dcut->GetZaxis()->SetTitleSize(0.04);
		distro_2dcut->GetXaxis()->SetLabelSize(0.04);
		distro_2dcut->GetYaxis()->SetLabelSize(0.04);
		distro_2dcut->GetZaxis()->SetLabelSize(0.04);
	sprintf(save_plot_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/%d.%d.%d_ThresholdMeanOverStdDev_NumChannel_Scan_A%d_c%d_%dEvents.png",year_now,month_now,day_now,station,config,num_total);
	c3->SaveAs(save_plot_title);
	delete c3;

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
