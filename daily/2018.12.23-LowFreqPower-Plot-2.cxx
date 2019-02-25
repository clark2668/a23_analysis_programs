////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  2018.12.23-LowFreqPower-Plot.cxx 
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

#include "tools_PlottingFns.h"
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"

using namespace std;

int PlotThisEvent(int station, int year, int runNum, int event);

int main(int argc, char **argv)
{
	gStyle->SetOptStat(110011);
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	stringstream ss;
	
	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <station> <year> <low freq power file>"<<endl;
		return 0;
	}
	int station = atoi(argv[1]);
	int year = atoi(argv[2]);

	TH1D *distro[16];
	TH1D *glitch[16];
	TH2D *distro_2d = new TH2D("","",100,0.005,1.005,18,-1.5,16.5);
	for(int i=0; i<16; i++){
		distro[i] = new TH1D("","",100,0,1); 
		glitch[i] = new TH1D("","",100,0,1);
	}

	int num_total=0;
	int num_total_glitch=0;

	int glitch_number[16]={0};

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
		int waveform_length[16];
		double frac_power_75[16];
		double frac_power_110[16];
		double frac_power_150[16];
		int runNum;
		bool hasDigitizerError;
		inTree->SetBranchAddress("isCal", &isCal);
		inTree->SetBranchAddress("isSoft", &isSoft);
		inTree->SetBranchAddress("waveform_length", &waveform_length);
		inTree->SetBranchAddress("frac_power_75", &frac_power_75);
		inTree->SetBranchAddress("frac_power_110", &frac_power_110);
		inTree->SetBranchAddress("frac_power_150", &frac_power_150);
		inTree->SetBranchAddress("run",&runNum);
		inTree->SetBranchAddress("hasDigitizerError",&hasDigitizerError);

		int numEntries = inTree->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;

		//now to loop over events
		for(int event=0; event<numEntries; event++){
			inTree->GetEvent(event);

			if(isCal || isSoft || hasDigitizerError) continue;
			num_total++;
			for(int i=0; i<16; i++) distro[i]->Fill(frac_power_75[i]);

			for(double power=0.; power<1.; power+=0.01){
				int this_glitch_counter=0;
				for(int i=0; i<16; i++){
					// if(frac_power_75[i]>power && i!=3 && i!=7 && i!=11 && i!=15){
					if(frac_power_150[i]>power){
						this_glitch_counter++;
					}
				}
				distro_2d->Fill(power,this_glitch_counter);
			}
		}
		fpIn->Close();
		delete fpIn;
	} //end loop over input files
	distro_2d->Scale(1/double(num_total));

	// for(int i=0; i<distro_2d->GetNbinsX(); i++){
	// 	cout<<"Bin "<<i<<" center is "<<distro_2d->GetXaxis()->GetBinCenter(i)<<" and width is "<<distro_2d->GetXaxis()->GetBinWidth(i)<<endl;
	// }
	
	TCanvas *c = new TCanvas("","",2*1100,2*850);
	c->Divide(4,4);
	for(int i=0; i<16; i++){
		c->cd(i+1);
		distro[i]->Draw("");
		distro[i]->SetLineWidth(2);
		distro[i]->SetLineColor(kBlue);
		gPad->SetLogy();
		distro[i]->GetYaxis()->SetRangeUser(1.,1e8);
		distro[i]->GetYaxis()->SetTitle("Number of Events");
		distro[i]->GetXaxis()->SetTitle("Fraction of Power Below 75 MHz");
	}
	char save_plot_title[400];
	sprintf(save_plot_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/%d.%d.%d_Distro_Glitch_A%d_%d_%dEvents.png",year_now,month_now,day_now,station,year,num_total);
	c->SaveAs(save_plot_title);
	delete c;

	gStyle->SetOptStat(0);
	TCanvas *c2 = new TCanvas("","",2*1100,2*850);
	c2->SetGrid();
	distro_2d->Draw("colz");
	distro_2d->GetYaxis()->SetNdivisions(16,0,0);
	distro_2d->GetYaxis()->SetTitle("Number of Channels");
	distro_2d->GetXaxis()->SetTitle("Fraction of Power Below 150 MHz");
	distro_2d->GetZaxis()->SetTitle("Fraction of Events");
	distro_2d->GetZaxis()->SetTitleOffset(1.2);
	gPad->SetLogz();
	gPad->SetRightMargin(0.15);
	distro_2d->GetZaxis()->SetRangeUser(1e-8,1.);		
	sprintf(save_plot_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/%d.%d.%d_2D_Distro_Glitch_A%d_%d_%dEvents.png",year_now,month_now,day_now,station,year,num_total);
	c2->SaveAs(save_plot_title);
	delete c2;
}