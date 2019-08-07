////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  2018.12.23-LowFreqPower.cxx 
////  store fraction of spectrum power below 75 MHz
////
////  Dec 2018
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <string>
#include <sstream>

//AraRoot Includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraQualCuts.h"
#include "FFTtools.h"

#include "tools_PlottingFns.h"
#include "tools_WaveformFns.h"
#include "tools_Cuts.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"

using namespace std;

double MyFunc(TGraph *grIn, int evt_num, int chan, vector<double> &means, vector<double> &stdDev, vector<double> MeanOverstdDev);
void MyFunc2(TGraph *grIn, vector<double> &means);

int main(int argc, char **argv)
{

	if(argc<4) {
		std::cout << "Usage\n" << argv[0] << " <station> <input_file> <output_location> "<<endl;
		return -1;
	}
	int station = atoi(argv[1]);

	/*
	arguments
	0: exec
	1: station
	2: input data file
	3: output location
	*/
	
	TFile *fpIn = TFile::Open(argv[2]);
	if(!fpIn) {
		std::cout << "Can't open file\n";
		return -1;
	}
	TTree *eventTree = (TTree*) fpIn->Get("eventTree");
	if(!eventTree) {
		std::cout << "Can't find eventTree\n";
		return -1;
	}
	RawAtriStationEvent *rawAtriEvPtr=0;
	eventTree->SetBranchAddress("event",&rawAtriEvPtr);
	int run;
	eventTree->SetBranchAddress("run",&run);
	eventTree->GetEntry(0);
	printf("Run Number %d \n", run);

	char ped_file_name[400];
	sprintf(ped_file_name,"/fs/scratch/PAS0654/ara/peds/run_specific_peds/A%d/all_peds/event%d_specificPeds.dat",station,run);
	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	calibrator->setAtriPedFile(ped_file_name,station); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist

	AraQualCuts *qualCut = AraQualCuts::Instance(); //we also need a qual cuts tool

	char outfile_name[400];
	sprintf(outfile_name,"%s/spare_channel_stats_wMeanOverStdDev_wSpareBlockMeans_run%d.root",argv[3],run);

	TFile *fpOut = TFile::Open(outfile_name, "RECREATE");
	TTree* outTree = new TTree("outTree", "outTree");
	bool isCal;
	bool isSoft;
	bool isShort;
	double spareChannelMaxSamp[4];
	double spareChannelMinSamp[4];
	double spareChannelRMS[4];
	double spareChannelMaxMeanOverStdDev[4];
	// vector< vector< double> > 
	vector<vector< double > > spareChannelBlockMeans;
	int runNum;
	int unixTime;
	int unixTimeUs;
	int timeStamp;
	int eventNumber;
	bool hasDigitizerError=false;
	bool hasA3S4Issue;
	outTree->Branch("isCal", &isCal, "isCal/O");
	outTree->Branch("isSoft", &isSoft, "isSoft/O");
	outTree->Branch("isShort", &isShort, "isShort/O");
	outTree->Branch("unixTime",&unixTime, "unixTime/I");
	outTree->Branch("unixTimeUs",&unixTimeUs, "unixTimeUs/I");
	outTree->Branch("timeStamp",&timeStamp, "timeStamp/I");
	outTree->Branch("eventNumber",&eventNumber, "eventNumber/I");
	outTree->Branch("hasDigitizerError", &hasDigitizerError, "hasDigitizerError/O");
	outTree->Branch("spareChannelMaxSamp", &spareChannelMaxSamp, "spareChannelMaxSamp[4]/D");
	outTree->Branch("spareChannelMinSamp", &spareChannelMinSamp, "spareChannelMinSamp[4]/D");
	outTree->Branch("spareChannelRMS", &spareChannelRMS, "spareChannelRMS[4]/D");
	outTree->Branch("spareChannelMaxMeanOverStdDev", &spareChannelMaxMeanOverStdDev, "spareChannelMaxMeanOverStdDev[4]/D");
	outTree->Branch("spareChannelBlockMeans", &spareChannelBlockMeans);

	runNum=run;
	outTree->Branch("run",&runNum);

	Long64_t numEntries=eventTree->GetEntries();
	int start=0;
	for(Long64_t event=start;event<numEntries;event++) {
		cout<<"On event "<<event<<endl;
		eventTree->GetEntry(event);
		isCal = rawAtriEvPtr->isCalpulserEvent();
		isSoft = rawAtriEvPtr->isSoftwareTrigger();
		unixTime = (int)rawAtriEvPtr->unixTime;
		unixTimeUs= (int)rawAtriEvPtr->unixTimeUs;
		timeStamp = (int)rawAtriEvPtr->timeStamp;
		eventNumber = (int)rawAtriEvPtr->eventNumber;

		UsefulAtriStationEvent *ev = new UsefulAtriStationEvent(rawAtriEvPtr,AraCalType::kLatestCalib);
		// hasDigitizerError = !(qualCut->isGoodEvent(ev));

		vector<TGraph*> electChansGraphs;
		electChansGraphs.push_back(ev->getGraphFromElecChan(6));
		electChansGraphs.push_back(ev->getGraphFromElecChan(14));
		electChansGraphs.push_back(ev->getGraphFromElecChan(22));
		electChansGraphs.push_back(ev->getGraphFromElecChan(30));
		for(int i=0; i<4; i++){
			if(electChansGraphs[i]->GetN()<64) isShort=true;
			spareChannelMaxSamp[i]=electChansGraphs[i]->GetMaximum();
			spareChannelMinSamp[i]=electChansGraphs[i]->GetMinimum();
			spareChannelRMS[i]=electChansGraphs[i]->GetRMS(2);
			vector<double> means;
			vector<double> stdDevs;
			vector<double> MeanOverstdDev;			
			spareChannelMaxMeanOverStdDev[i] = MyFunc(electChansGraphs[i],event,i,means,stdDevs,MeanOverstdDev);

			vector<double> theseMeans;
			MyFunc2(electChansGraphs[i], theseMeans);
			spareChannelBlockMeans.push_back(theseMeans);
		}
		outTree->Fill();
		deleteGraphVector(electChansGraphs);
		spareChannelBlockMeans.clear();
		delete ev;
	} //loop over events
	
	fpOut->Write();
	fpOut->Close();
	
	fpIn->Close();
	delete fpIn;
}

void MyFunc2(TGraph *grIn, vector<double> &means){
	int n_input = grIn->GetN();
	double *oldX = grIn->GetX();
	double *oldY = grIn->GetY();

	deque <double> inX;
	deque <double> inY;

	for(int samp=0; samp<n_input; samp++){
		inX.push_back(oldX[samp]);
		inY.push_back(oldY[samp]);
	}

	while(inX.size()>63){
		vector <double> sub_X;
		vector <double> sub_Y;
		int num_to_pop=0;
		for(int samp=0; samp<64; samp++){
			sub_X.push_back(inX[samp]);
			sub_Y.push_back(inY[samp]);
			num_to_pop++;
		}
		double average = std::accumulate( sub_Y.begin(), sub_Y.end(), 0.0)/sub_Y.size();
		means.push_back(average);
		for(int iPop=0; iPop<num_to_pop; iPop++){
			inX.pop_front();
			inY.pop_front();
		}
	}
}


double MyFunc(TGraph *grIn, int evt_num, int chan, vector<double> &means, vector<double> &stdDev, vector<double> MeanOverstdDev){

	int n_input = grIn->GetN();
	double *oldX = grIn->GetX();
	double *oldY = grIn->GetY();

	deque <double> inX;
	deque <double> inY;	

	for(int samp=0; samp<n_input; samp++){
		inX.push_back(oldX[samp]);
		inY.push_back(oldY[samp]);
	}

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

		double forStdDev=0.;
		for(int samp=0; samp<sub_X.size(); samp++){
			double thisY = sub_Y[samp];
			forStdDev+=pow(thisY - mean_this_section,2.);
		}
		forStdDev/=double(sub_X.size());
		forStdDev=sqrt(forStdDev);

		means.push_back(abs(mean_this_section));
		stdDev.push_back(forStdDev);	
		MeanOverstdDev.push_back(abs(mean_this_section/(forStdDev/sqrt(double(sub_X.size())))));
		// MeanOverstdDev.push_back(abs(mean_this_section/(forStdDev/double(sub_X.size()))));

		for(int iPop=0; iPop<num_to_pop; iPop++){
			inX.pop_front();
			inY.pop_front();
		}
	}
	std::sort(MeanOverstdDev.begin(), MeanOverstdDev.end(), std::greater<double>());
	std::sort(means.begin(), means.end(), std::greater<double>());
	// return MeanOverstdDev[0];
	return means[0];
}