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
#include <algorithm>
#include <vector>
#include <deque>

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "FFTtools.h"

//AraRoot includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraQualCuts.h"
using namespace std;

double MaxMeanBlock(TGraph *grIn, int evt_num, int chan, vector<double> &means, bool print);

int countExcursions(vector<double> maxMeans){
	if(maxMeans.size()!=16){
		cout<<"Something wrong! More than 16 elements!"<<endl;
		return 0;
	}
	int numExcursions=0;
	for(int i=0; i<16; i++){
		if((i+1)%4==0) continue;
		if(maxMeans[i]>150)
			numExcursions++;
	}
	return numExcursions;
}


int main(int argc, char **argv)
{

	if(argc<4){
		cout<< "Usage\n" << argv[0] << " <station> <run num> <event>"<<endl;;
		return -1;
	}
	int station = atoi(argv[1]);
	int runNum = atoi(argv[2]);
	int event = atoi(argv[3]);
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *DataDirPath(getenv("DATA_DIR"));
	if (DataDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;


	gStyle->SetOptStat(11);

	char run_file_name[400];
	sprintf(run_file_name,"%s/RawData/A%d/all_runs/event%d.root",DataDirPath,station,runNum,runNum);
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
	
	UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawPtr,AraCalType::kLatestCalib);

	int unixTime = (int)rawPtr->unixTime;
	int unixTimeUs =(int)rawPtr->unixTimeUs;
	int timeStamp = (int)rawPtr->timeStamp;
	printf("Unixtime is %d \n", unixTime);
	printf("Unixtime microsecond is %d \n", unixTimeUs);
	printf("timeStamp is %d \n", timeStamp);
	printf("Event Number is %d \n", realAtriEvPtr->eventNumber);

	AraQualCuts *qualCut = AraQualCuts::Instance();
	bool quality = qualCut->isGoodEvent(realAtriEvPtr);
	printf("Quality is %d \n", quality);

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
		waveforms.push_back(realAtriEvPtr->getGraphFromRFChan(i));
		waveforms[i]->GetXaxis()->SetTitle(xLabel.c_str());
		waveforms[i]->GetYaxis()->SetTitle(yLabel.c_str());
		waveforms[i]->SetTitle(titlesForGraphs[i].c_str());
	}

	// rolling mean test
	// vector<TGraph*> rollingMeans;
	// vector<double> maxMeans;
	// rollingMeans.resize(16);
	// for(int i=0; i<16; i++){
	// 	if(waveforms[i]->GetN()<64){
	//  		continue;
	//  	}
	//  	TGraph *grInt = FFTtools::getInterpolatedGraph(waveforms[i],0.5);
	//  	if(grInt->GetN()<101){
	//  		delete grInt;
	//  		continue;
	//  	}
	//  	rollingMeans[i]=(qualCut->getRollingMean(grInt,100));
	//  	double thisMax = abs(TMath::MaxElement(rollingMeans[i]->GetN(), rollingMeans[i]->GetY()));
	//  	double thisMin = abs(TMath::MinElement(rollingMeans[i]->GetN(), rollingMeans[i]->GetY()));
	//  	double thisAbs;
	//  	if(abs(thisMax)>abs(thisMin))
	//  		thisAbs=abs(thisMax);
	//  	else
	//  		thisAbs=abs(thisMin);
	//  	maxMeans.push_back(thisAbs);
	//  	printf("Chan %d max abs is %4.2f\n",i,thisAbs);
	// }
	// for(int i=0; i<16; i++) delete rollingMeans[i];
	// int numExcursions = countExcursions(maxMeans);
	// printf("Found %d excursions\n",numExcursions);

	vector< vector<double> > vvmeans;
	for(int i=0; i<16; i++){
		// double thisRMS = waveforms[i]->GetRMS(2);
		// printf("Graph %d RMS is %.2f \n",i,thisRMS);
		vector<double> vmeans;
		double this_max = MaxMeanBlock(waveforms[i], event, i, vmeans, true);
		vvmeans.push_back(vmeans);
	}
	for(int i=0; i<16; i++){
		
	}

	// for(int i=0; i<vvmeans[0].size(); i++){
	// 	bool allSameSignAfter=true;
	// 	double startVal = vvmeans[0][i];
	// 	for(int j=i+1; j<vvmeans[0].size(); j++){
	// 		bool sameSign = signbit(startVal)==signbit(vvmeans[0][j]);
	// 		if(!sameSign) allSameSignAfter=false;
	// 		// printf("	Block %d value is %3.2f which is %d sign \n", j, vvmeans[0][j], sameSign);
	// 	}
	// 	printf("All Same Sign after block %d is %d \n", i, allSameSignAfter);
	// }

	// for(int i=0; i<16; i++){
	// 	if(i==2 || i==3 || i==6 || i==7 || i==10 || i==11 || i==14 || i==15) continue;
	// 	printf("Chan %2d : ",i);
	// 	for(int j=0; j<vvmeans[i].size();j++)
	// 		printf("%4.f : ", vvmeans[i][j]);
	// 	cout<<endl;
	// }
	vector<TGraph*> dummy;
	for(int i=0; i<16; i++){
		vector<double> thisX;
		vector<double> thisY;
		thisX.push_back(-200.);
		thisX.push_back(400.);
		thisY.push_back(-200.);
		thisY.push_back(200.);
		dummy.push_back(new TGraph(thisX.size(), &thisX[0], &thisY[0]));
		dummy[i]->GetXaxis()->SetTitle("Time (ns)");
		dummy[i]->GetYaxis()->SetTitle("Voltage (mV)");
		dummy[i]->GetXaxis()->SetLabelSize(0.07);
		dummy[i]->GetYaxis()->SetLabelSize(0.07);
		dummy[i]->GetXaxis()->SetTitleSize(0.07);
		dummy[i]->GetYaxis()->SetTitleSize(0.07);
	}

	char save_temp_title[300];
	sprintf(save_temp_title,"outputs/%d.%d.%d_Run%d_Ev%d_Waveforms.png",year_now,month_now,day_now,runNum,event);
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

	for(int i=0; i<16; i++){
		delete waveforms[i];
	}
	delete realAtriEvPtr;
	mapFile->Close();
	delete mapFile;
	return 0;
}

double MaxMeanBlock(TGraph *grIn, int evt_num, int chan, vector<double> &means, bool print){

	int n_input = grIn->GetN();
	double *oldX = grIn->GetX();
	double *oldY = grIn->GetY();

	deque <double> inX;
	deque <double> inY;	

	for(int samp=0; samp<n_input; samp++){
		inX.push_back(oldX[samp]);
		inY.push_back(oldY[samp]);
	}

	// vector <double> means;
	vector <TGraph*> grPieces;

	while(inX.size()>0){
		vector <double> sub_X;
		vector <double> sub_Y;
		double first_time = inX[0];
		int num_to_pop=0;
		double mean_this_section=0.;
		for(int samp=0; samp<inX.size(); samp++){
			if(inX[samp]<=first_time+20.){
				sub_X.push_back(inX[samp]);
				sub_Y.push_back(inY[samp]);
				num_to_pop++;
				mean_this_section+=inY[samp];
			}
		}
		mean_this_section/=double(sub_X.size());
		means.push_back((mean_this_section));
		for(int samp=0; samp<sub_X.size(); samp++){ //now fix the mean
			sub_Y[samp]=mean_this_section;
		}
		grPieces.push_back(new TGraph(sub_X.size(),&sub_X[0],&sub_Y[0]));

		for(int iPop=0; iPop<num_to_pop; iPop++){
			inX.pop_front();
			inY.pop_front();
		}
	}

	int colors [28] = { kBlue, kSpring, kYellow, kTeal, kMagenta, kAzure, kRed, kCyan, kViolet, kGreen, kOrange, kPink, kBlue, kSpring, kYellow,kTeal, kMagenta, kAzure, kRed, kCyan, kViolet, kGreen, kOrange, kPink, kBlue, kSpring, kYellow, kTeal}; 

	if(print){
		TCanvas *c = new TCanvas("","",1100,850);
		grIn->Draw("AL");
		for(int i=0; i<grPieces.size(); i++){
			grPieces[i]->Draw("same");
			grPieces[i]->SetLineColor(colors[i]);
			grPieces[i]->SetLineWidth(3);
		}
		char title[400];
		sprintf(title,"outputs/GlitchWaveform_ev%d_chan%d.png",evt_num,chan);
		c->SaveAs(title);	
		delete c;
	}
	for(int i=0; i<grPieces.size(); i++){
		delete grPieces[i];
	}
	// for(int i=0; i<means.size(); i++){
	// 	printf("Block %d has mean %.2f \n", i, means[i]);
	// }
	std::sort(means.begin(), means.end(), std::greater<double>());
	// printf("Biggest value is %.2f \n", means[0]);
	return means[0];
}