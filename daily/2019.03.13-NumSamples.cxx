////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	num_samples.cxx 
////	make histograms of the number of samples in A23 data set
////
////	Mar 2019
////////////////////////////////////////////////////////////////////////////////

//C++
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace std;
void configure(TH1D *gr);

int main(int argc, char **argv)
{
	gStyle->SetOptStat(0);
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

	stringstream ss;
	
	if(argc<4){
		cout<< "Usage\n" << argv[0] << " <station> <year> <joined filename 1> <joined filename 2 > ... <joined filename x>"<<endl;
		return 0;
	}
	int station = atoi(argv[1]);
	int year = atoi(argv[2]);
	
	TH1D *num_samps[16];
	for(int i=0; i<16; i++){
		num_samps[i] = new TH1D("","",200,0,2000);
	}
	int num_total=0;

	for(int file_num=3; file_num<argc; file_num++){

		cout << "Run " << file_num << " :: " << argv[file_num] << endl;
		
		TFile *inputFile = TFile::Open(argv[file_num]);
		if(!inputFile){
			cout<<"Can't open joined file!"<<endl;
			return -1;
		}
		ss.str("");
		ss << "OutputTree_filter";
		TTree *inputTree_filter = (TTree*) inputFile->Get(ss.str().c_str());
		if(!inputTree_filter){
			cout<<"Can't open filter tree"<<endl;
			return -1;
		}
		bool isCalPulser;
		bool isSoftTrigger;
		int waveformLength[16];
		bool hasDigitizerError;
		inputTree_filter->SetBranchAddress("isCalpulser",&isCalPulser);
		inputTree_filter->SetBranchAddress("isSoftTrigger",&isSoftTrigger);
		inputTree_filter->SetBranchAddress("waveformLength",&waveformLength);
		inputTree_filter->SetBranchAddress("hasDigitizerError",&hasDigitizerError);

		int numEntries = inputTree_filter->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;

		for(int event=0; event<numEntries; event++){
			// if(event%starEvery==0) { std::cout << "	On event "<<event<<endl;}
			inputTree_filter->GetEvent(event);
			num_total++;
			for(int i=0; i<16; i++){
				num_samps[i]->Fill(waveformLength[i]);
			}
		}//loop over events
		inputFile->Close();
		delete inputFile;
	} //end loop over input files

	gStyle->SetOptStat(0);
	gStyle->SetStatY(0.9);
	gStyle->SetStatX(0.9);
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.2);

	TCanvas *c = new TCanvas("","",2.1*850,2.1*850);
	c->Divide(4,4);
	for(int i=0; i<16; i++){
		c->cd(i+1);
		configure(num_samps[i]);

		gPad->SetTopMargin(0.1);
		gPad->SetRightMargin(0.03);
		gPad->SetLeftMargin(0.12);
		gPad->SetBottomMargin(0.11);

		num_samps[i]->Draw("");
		num_samps[i]->GetXaxis()->SetTitle("Numer of Samples");
		num_samps[i]->GetYaxis()->SetTitle("Counts");
		gPad->SetLogy();
	}
	char title[300];
	sprintf(title, "%s/A%d_%d_%dEvents_NumSamps.png",plotPath,station,year,int(num_total));
	c->SaveAs(title);
	delete c;
}

void configure(TH1D *gr){
	gr->GetXaxis()->SetLabelSize(0.07);
	gr->GetXaxis()->SetTitleSize(0.07);
	
	gr->GetYaxis()->SetLabelSize(0.07);
	gr->GetYaxis()->SetTitleSize(0.07);
	gr->GetYaxis()->SetTitleOffset(1.2);
	
	gStyle->SetTitleFontSize(0.07);
	gr->GetXaxis()->SetNdivisions(4,0,0,false);
	gr->SetLineWidth(3);
}