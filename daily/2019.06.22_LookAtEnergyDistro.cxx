////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	2019.06.22_LookAtCorrDistro.cxx 
////	Just look at distro of corr peaks, really straightforward
////
////	June 2019
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

//AraRoot includes
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_inputParameters.h"
#include "tools_outputObjects.h"
#include "tools_Cuts.h"
#include "tools_CommandLine.h"

int main(int argc, char **argv)
{
	// gStyle->SetOptStat(0);
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

	if(argc<5){
		cout<< "Usage\n" << argv[0] << " <isSim?> <station> <config> <reco filename 1 > <reco filename 2 > ... <reco filename n>"<<endl;
		return 0;
	}
	int isSimulation = atoi(argv[1]);
	int station = atoi(argv[2]);
	int config = atoi(argv[3]);

	TH1D *hEnergy = new TH1D("","",10,16,21);

	for(int file_num=4; file_num<argc; file_num++){

		string file = string(argv[file_num]);
		string wordRun = "run_";
		size_t foundRun = file.find(wordRun);
		string wordFilter = "_joined";
		size_t foundFilter = file.find(wordFilter);
		size_t diff=(foundFilter-wordRun.length())-foundRun;
		string strRunNum = file.substr(foundRun+4,diff);
		int runNum = atoi(strRunNum.c_str());


		double energy;
		double weight;
		TFile *inputFile = TFile::Open(argv[file_num]);
		if(!inputFile){
			cout<<"Can't open joined file!"<<endl;			
			return -1;
		}
		TTree *inTree = (TTree*) inputFile->Get("OutputTree");
		if(!inTree){
			cout<<"Can't find the in tree"<<endl;
			return -1;
		}
		inTree->SetBranchAddress("energy",&energy);
		inTree->SetBranchAddress("weight",&weight);

		int numEntries = inTree->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;
		// cout<<"Star every is "<<starEvery<<endl;

		for(int event=0; event<numEntries; event++){
			if(event%starEvery==0) {
				// std::cout<<"*";
			}
			inTree->GetEvent(event);
			hEnergy->Fill(TMath::Log10(energy),weight);
		}
	}


	TCanvas *c = new TCanvas("","",850,850);
	hEnergy->Draw();
	gPad->SetLogy();
	char save_title[400];
	sprintf(save_title,"%s/corr_study/%d.%d.%d_A%d_c%d_Sim%d_EnergyDistro.png",plotPath,year_now,month_now,day_now,station,config,isSimulation);
	c->SaveAs(save_title);
}