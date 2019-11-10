
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  print_event.cxx
////
////  Nov 2019
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <iomanip>
#include <sstream>

//AraRoot Includes
#include "RawIcrrStationEvent.h"
#include "RawAtriStationEvent.h"
#include "UsefulAraStationEvent.h"
#include "UsefulIcrrStationEvent.h"
#include "UsefulAtriStationEvent.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "FFTtools.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"

RawAtriStationEvent *rawAtriEvPtr;

using namespace std;

int main(int argc, char **argv)
{


	char *DataDirPath(getenv("DATA_DIR_100"));
	if (DataDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;

	gStyle->SetOptStat(11);

	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	if(argc<4){
		cout<<"Use like "<<argv[0]<<" <station> <runNum> <eventNumber> "<<endl;
		return -1;
	}

	int station = atoi(argv[1]);
	int runNum = atoi(argv[2]);
	int eventNumber = atoi(argv[3]);

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
	eventTree->GetEvent(eventNumber);

	int stationID = rawPtr->stationId;
	char ped_file_name[400];
	sprintf(ped_file_name,"%s/run_specific_peds/A%d/all_peds/event%d_specificPeds.dat",PedDirPath,station,runNum);
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

	vector<TGraph*> waveforms;
	for(int i=0; i<16; i++){
		TGraph *gr = realAtriEvPtr->getGraphFromRFChan(i);
		gr->GetXaxis()->SetTitle(xLabel.c_str());
		gr->GetYaxis()->SetTitle(yLabel.c_str());
		gr->SetTitle(titlesForGraphs[i].c_str());
		waveforms.push_back(gr);
	}

	double powerInside[16];
	for(int chan=0; chan<16; chan++) powerInside[chan]=0.;
	for(int chan=0; chan<16; chan++){
		Double_t *yVals = waveforms[chan]->GetY();
		int N = waveforms[chan]->GetN();
		for(int samp=0; samp<N; samp++) powerInside[chan]+=yVals[samp]*yVals[samp];
		printf("Chan %2d has power %.2f \n", chan, powerInside[chan]);
	}
	
	vector<double> powerStrings;
	for(int string=0; string<4; string++){
	  if(string!=4)
	    powerStrings.push_back(powerInside[string]+powerInside[string+4]+powerInside[string+8]+powerInside[string+12]);
	  else if(string==3)
	    powerStrings.push_back(powerInside[string]+powerInside[string+4]+powerInside[string+8]);
	}
	for(int string=0; string<4; string++){
		cout<<"String 0 has power "<<powerStrings[string]<<endl;
	}
	  powerStrings[0]/=4.;
	  powerStrings[1]/=4.;
	  powerStrings[2]/=4.;
	  powerStrings[3]/=3.;
	std::sort(powerStrings.begin(), powerStrings.end()); //sort smallest to largest
	cout<<"Average Second Lowest "<<powerStrings[1]<<endl;
	cout<<"Average Highest "<<powerStrings[3]<<endl;
	cout<<"Ratio average highest to second average lowest "<<powerStrings[3]/powerStrings[2]<<endl;


	vector<TGraph*> dummy;
	for(int i=0; i<16; i++){
		vector<double> thisX;
		vector<double> thisY;
		thisX.push_back(-200.);
		thisX.push_back(700.);
		// thisX.push_back(45);
		// thisX.push_back(55);
		thisY.push_back(-1500.);
		thisY.push_back(1500.);
		// thisY.push_back(-3000.);
		// thisY.push_back(3000.);
		dummy.push_back(new TGraph(thisX.size(), &thisX[0], &thisY[0]));
		dummy[i]->GetXaxis()->SetTitle("Time (ns)");
		dummy[i]->GetYaxis()->SetTitle("Voltage (mV)");
		dummy[i]->GetXaxis()->SetLabelSize(0.07);
		dummy[i]->GetYaxis()->SetLabelSize(0.07);
		dummy[i]->GetXaxis()->SetTitleSize(0.07);
		dummy[i]->GetYaxis()->SetTitleSize(0.07);
		dummy[i]->SetTitle(titlesForGraphs[i].c_str());
		dummy[i]->GetYaxis()->SetTitleOffset(1.2);
	}

	char save_temp_title[300];
	sprintf(save_temp_title,"./%d.%d.%d_Run%d_Ev%d_Waveforms.png",year_now,month_now,day_now,runNum,eventNumber);
	TCanvas *cWave = new TCanvas("","",4*1100,4*850);
	cWave->Divide(4,4);
	for(int i=0; i<16; i++){
		cWave->cd(i+1);
		dummy[i]->Draw("AP");
		waveforms[i]->Draw("LPsame");
		waveforms[i]->SetLineWidth(2);
		// waveforms[i]->SetMarkerStyle(kFullCircle);
		// waveforms[i]->SetMarkerSize(2);
	}
	cWave->SaveAs(save_temp_title);
	delete cWave;
}
