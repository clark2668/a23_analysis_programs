////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	print_distro.cxx
//// 	print reco vs time
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
#include "TCanvas.h"
#include "TStyle.h"
#include "TTimeStamp.h"
#include "TH1D.h"

//AraRoot includes
#include "RawAtriStationEvent.h"

using namespace std;

int getThisRunNum(char *thearg){
	string chRun = "run";
	string file = string(thearg);
	size_t foundRun = file.find(chRun);
	string strRunNum = file.substr(foundRun+4,4);
	int runNum = atoi(strRunNum.c_str());
	return runNum;
}

int main(int argc, char **argv)
{

	if(argc<4){
		cout<< "Usage\n" << argv[0] << " <station> <output_location> <filename>"<<endl;;
		return -1;
	}
	int station = atoi(argv[1]);
	string output_location = argv[2];

	for(int file_num=3; file_num<argc; file_num++){

		cout<<"On file "<<argv[file_num]<<endl;

		// int runNum=getThisRunNum(argv[file_num]);
		
		TFile *inputFile = TFile::Open(argv[file_num]);
		if(!inputFile){
			cout<<"Can't open joined file!"<<endl;
			return -1;
		}
		TTree *eventTree = (TTree*) inputFile->Get("eventTree");
		if(!eventTree){
			cout<<"Can't get eventTree!"<<endl;
			return -1;
		}
		RawAtriStationEvent *rawAtriEvPtr=0;
		eventTree->SetBranchAddress("event",&rawAtriEvPtr);
		int runNum;
		eventTree->SetBranchAddress("run",&runNum);

		int numEntries = eventTree->GetEntries();

		eventTree->GetEntry(0);
		int first_event_unixtime = (int) rawAtriEvPtr->unixTime;

		eventTree->GetEntry(numEntries-1);
		int last_event_unixtime = (int) rawAtriEvPtr->unixTime;

		TTimeStamp start(first_event_unixtime-600);
		TTimeStamp stop(last_event_unixtime+600);
		int start_bin = start.GetSec();
		int stop_bin = stop.GetSec();

		TH1D *h1_ev_vs_time = new TH1D("","",stop_bin-start_bin,start_bin,stop_bin);
			h1_ev_vs_time->GetXaxis()->SetTimeDisplay(1);
			h1_ev_vs_time->GetXaxis()->SetTimeFormat("%H:%M");
			h1_ev_vs_time->GetXaxis()->SetTimeOffset(0.,"GMT");

		//now to loop over events
		for(int event=0; event<numEntries; event++){
			if(event%10000==0) cout<<"on event "<<event<<endl;
			// cout<<"On event"<<event<<endl;
			eventTree->GetEntry(event);
			h1_ev_vs_time->Fill(rawAtriEvPtr->unixTime);
		}

		int maxBin = h1_ev_vs_time->GetMaximumBin();
		double maxValue = h1_ev_vs_time->GetBinContent(maxBin);
		cout<<"maxValue is "<<maxValue<<endl;

		char title_txt[200];
		sprintf(title_txt,"/users/PAS0654/osu0673/A23_analysis_new2/AraRoot/analysis/a23_analysis_programs/other/trigger_rate/rates_A%d.txt",station);
		FILE *fout = fopen(title_txt, "a");
		fprintf(fout,"%d, %d \n",runNum, int(maxValue));
		fclose(fout);//close sigmavsfreq.txt file

		// TCanvas *c = new TCanvas("","",1100,850);
		// h1_ev_vs_time->Draw("");
		// h1_ev_vs_time->GetYaxis()->SetTitle("Number of Events");
		// char title[150];
		// sprintf(title,"%s/rate_run%d.png",output_location.c_str(),runNum);
		// c->SaveAs(title);
		// delete c;

		inputFile->Close();
		delete inputFile;
	}
}
