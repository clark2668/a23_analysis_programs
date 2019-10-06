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

#include "TTimeStamp.h"


using namespace std;
void configure(TH1D *gr);

int main(int argc, char **argv)
{
	gStyle->SetOptStat(11);
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

	stringstream ss;
	
	if(argc<4){
		cout<< "Usage\n" << argv[0] << " <station> <config> <save_vals_file_1> <save_vals_file_2> ... <save_vals_file_x>"<<endl;
		return 0;
	}
	int station = atoi(argv[1]);
	int config = atoi(argv[2]);
	
	TH1D *num_samps[16];
	TH1D *num_samps_special[16];

	for(int i=0; i<16; i++){
		num_samps[i] = new TH1D("","",300,0,3000);
		num_samps_special[i] = new TH1D("","",300,0,3000);
	}

	TTimeStamp start;
	TTimeStamp stop;

	start.Set(2015, 01, 00, 00, 00,0,0,true,0);
	 stop.Set(2015, 12, 31, 00, 00,0,0,true,0);

	int start_bin = start.GetSec();
	int stop_bin = stop.GetSec();

	TH2D *samps_vs_time[16];
	for(int i=0; i<16; i++){
		ss.str("");
		ss << "Chan "<<i;
		// samps_vs_time[i] = new TH2D(ss.str().c_str(),ss.str().c_str(),(stop_bin-start_bin)/8760, start_bin, stop_bin, 100,0,4000);
		samps_vs_time[i] = new TH2D("","",(stop_bin-start_bin)/8760, start_bin, stop_bin, 100,0,3000);
		samps_vs_time[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		samps_vs_time[i]->GetXaxis()->SetTimeFormat("%y/%m");
		samps_vs_time[i]->GetXaxis()->SetTimeOffset(0.,"GMT");
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
		ss << "OutputTree";
		TTree *inputTree_filter = (TTree*) inputFile->Get(ss.str().c_str());
		if(!inputTree_filter){
			cout<<"Can't open filter tree"<<endl;
			return -1;
		}
		bool isCalPulser;
		bool isSoftTrigger;
		int waveformLength[16];
		bool hasDigitizerError;
		int eventNumber;
		int unixTime;
		inputTree_filter->SetBranchAddress("isCalpulser",&isCalPulser);
		inputTree_filter->SetBranchAddress("isSoftTrigger",&isSoftTrigger);
		inputTree_filter->SetBranchAddress("waveformLength",&waveformLength);
		inputTree_filter->SetBranchAddress("hasDigitizerError",&hasDigitizerError);
		inputTree_filter->SetBranchAddress("eventNumber",&eventNumber);
		inputTree_filter->SetBranchAddress("unixTime",&unixTime);

		ss.str("");
		ss << "AllTree";
		TTree *inputTree_All = (TTree*) inputFile->Get(ss.str().c_str());
		if(!inputTree_All){
			cout<<"Can't open all tree"<<endl;
			return -1;
		}
		int runNum;
		inputTree_All->SetBranchAddress("runNum",&runNum);

		int numEntries = inputTree_filter->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;

		int start=0;
		int stop = numEntries;
		// stop = 20;

		for(int event=0; event<stop; event++){
			// if(event%starEvery==0) { std::cout << "	On event "<<event<<endl;}
			inputTree_filter->GetEvent(event);
			inputTree_All->GetEvent(event);
			num_total++;
			for(int i=0; i<16; i++){
				num_samps[i]->Fill(waveformLength[i]);
				samps_vs_time[i]->Fill(unixTime, waveformLength[i]);
				if((runNum==5679 && eventNumber==30424)
					|| (runNum==4111 && eventNumber==57506)
					|| (runNum==4837 && eventNumber==72008)
					|| (runNum==5660 && eventNumber==107807)
					|| (runNum==5669 && eventNumber==72868)
					|| (runNum==5676 && eventNumber==90691)
					|| (runNum==6171 && eventNumber==49022)
					){
					num_samps_special[i]->Fill(waveformLength[i]);
					printf("Run %d, Event %d, Chan %d, Samps %d\n", runNum, eventNumber, i, waveformLength[i]);
				}
			}

			// for(int i=0; i<1; i++){
			// 	printf("Cal %d, Soft %d, Event %5d: Chan %4d, Samps %4d \n", isCalPulser, isSoftTrigger, event, i, waveformLength[i]);
			// }

		}//loop over events
		inputFile->Close();
		delete inputFile;
	} //end loop over input files

	gStyle->SetOptStat(111);
	// gStyle->SetStatY(0.9);
	// gStyle->SetStatX(0.9);
	// gStyle->SetStatW(0.2);
	// gStyle->SetStatH(0.2);

	stringstream ss1;
	vector<string> titlesForGraphs;
	for (int i = 0; i < 16; i++){
		ss1.str("");
		ss1 << "Channel " << i;
		titlesForGraphs.push_back(ss1.str());
	}

	TCanvas *c = new TCanvas("","",4.1*850,2.1*850);
	c->Divide(4,2);
	for(int i=0; i<8; i++){
		c->cd(i+1);
		configure(num_samps[i]);
		configure(num_samps_special[i]);
		num_samps[i]->GetYaxis()->SetRangeUser(0.1,2e8);
		// num_samps[i]->GetXaxis()->SetRangeUser(0,750);

		gPad->SetTopMargin(0.1);
		gPad->SetRightMargin(0.03);
		gPad->SetLeftMargin(0.12);
		gPad->SetBottomMargin(0.11);

		num_samps[i]->Draw("");
		num_samps[i]->GetXaxis()->SetTitle("Number of Samples");
		num_samps[i]->GetYaxis()->SetTitle("Number of Events");
		gPad->SetLogy();

		num_samps_special[i]->Draw("same");
		num_samps_special[i]->SetLineColor(kRed);

		double total_num = num_samps[i]->Integral();
		TAxis *xaxis = num_samps[i]->GetXaxis();
		Int_t binx_start = xaxis->FindBin(2100);
		Int_t binx_end = xaxis->FindBin(3500);

		double above = num_samps[i]->Integral(binx_start, binx_end);
		printf("Total %7.1f, Above %7.1f, percentage %1.4f \n", total_num, above, 100.*above/total_num);
		num_samps[i]->SetTitle(titlesForGraphs[i].c_str());
	}
	char title[300];
	sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_NumSamps.png",plotPath,year_now,month_now,day_now,station,config,int(num_total));
	c->SaveAs(title);
	delete c;

	// TCanvas *c2 = new TCanvas("","",2.1*850,2.1*850);
	// c2->Divide(4,4);
	// for(int i=0; i<16; i++){
	// 	c2->cd(i+1);
	// 	configure(num_samps[i]);
	// 	// num_samps[i]->GetYaxis()->SetRangeUser(0.1,2e7);
	// 	// num_samps[i]->GetXaxis()->SetRangeUser(0,750);

	// 	gPad->SetTopMargin(0.1);
	// 	gPad->SetRightMargin(0.03);
	// 	gPad->SetLeftMargin(0.12);
	// 	gPad->SetBottomMargin(0.11);


	// 	samps_vs_time[i]->Draw("colz");
	// 	samps_vs_time[i]->GetXaxis()->SetTitle("unixTime");
	// 	samps_vs_time[i]->GetYaxis()->SetTitle("Number of Samples");
	// 	samps_vs_time[i]->GetZaxis()->SetTitle("Number of Events");
	// 	gPad->SetLogz();

	// 	samps_vs_time[i]->SetTitle(titlesForGraphs[i].c_str());
	// }
	// sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_NumSamps_vs_Time.png",plotPath,year_now,month_now,day_now,station,config,int(num_total));
	// c2->SaveAs(title);
	// delete c2;


}

void configure(TH1D *gr){
	gr->GetXaxis()->SetLabelSize(0.07);
	gr->GetXaxis()->SetTitleSize(0.07);
	
	gr->GetYaxis()->SetLabelSize(0.07);
	gr->GetYaxis()->SetTitleSize(0.07);
	gr->GetYaxis()->SetTitleOffset(1.2);
	
	gStyle->SetTitleFontSize(0.07);
	gr->GetXaxis()->SetNdivisions(4,0,0,false);
	gr->SetLineWidth(2);
}