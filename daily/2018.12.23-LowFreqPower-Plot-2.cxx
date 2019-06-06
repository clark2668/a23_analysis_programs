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
#include "TLegend.h"

#include "tools_PlottingFns.h"
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"

using namespace std;

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
		cout<< "Usage\n" << argv[0] << " <station> <config> <low freq power file>"<<endl;
		return 0;
	}
	int station = atoi(argv[1]);
	int config = atoi(argv[2]);

	TH1D *distro[16];
	TH1D *glitch[16];
	TH1D *cal[16];
	TH2D *distro_2d = new TH2D("","",100,0.005,1.005,18,-1.5,16.5);
	TH2D *power_oob_vs_ib[16];
	TH1D *distro_above[16];
	TH1D *glitch_above[16];
	TH1D *cal_above[16];
	for(int i=0; i<16; i++){
		distro[i] = new TH1D("","",100,0,1); 
		glitch[i] = new TH1D("","",100,0,1);
		cal[i] = new TH1D("","",100,0,1);
		power_oob_vs_ib[i] = new TH2D("","",1e3,0,1e6,1e3,0,1e5);
		distro_above[i] = new TH1D("","",100,0,1); 
		glitch_above[i] = new TH1D("","",100,0,1);
		cal_above[i] = new TH1D("","",100,0,1);
	}

	TGraph *TroubleEventOverlay[16];
	for(int i=0; i<16; i++){
		TroubleEventOverlay[i]=new TGraph();
	}
	int numTroubleFound=0;

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
		double frac_above_850[16];
		double frac_between[16];
		double power_75[16];
		double power_110[16];
		double power_150[16];
		double above_850[16];
		double between[16];
		int runNum;
		bool hasDigitizerError;
		inTree->SetBranchAddress("isCal", &isCal);
		inTree->SetBranchAddress("isSoft", &isSoft);
		inTree->SetBranchAddress("waveform_length", &waveform_length);
		inTree->SetBranchAddress("frac_power_75", &frac_power_75);
		// inTree->SetBranchAddress("frac_power_110", &frac_power_110);
		// inTree->SetBranchAddress("frac_power_150", &frac_power_150);
		inTree->SetBranchAddress("frac_above_850", &frac_above_850);
		// inTree->SetBranchAddress("frac_between", &frac_between);
		inTree->SetBranchAddress("run",&runNum);
		inTree->SetBranchAddress("hasDigitizerError",&hasDigitizerError);
		inTree->SetBranchAddress("power_75", &power_75);
		inTree->SetBranchAddress("power_110", &power_110);
		inTree->SetBranchAddress("power_150", &power_150);
		inTree->SetBranchAddress("above_850", &above_850);
		inTree->SetBranchAddress("between", &between);

		int numEntries = inTree->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;

		//now to loop over events
		for(int event=0; event<numEntries; event++){
			inTree->GetEvent(event);

			if(isSoft || hasDigitizerError) continue;
			num_total++;

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

			// config 2
			// if(runNum==1075 && event==2953) TroubleEvent=true;
			if(runNum==476 && event==223) TroubleEvent=true;
			if(runNum==476 && event==506) TroubleEvent=true;
			if(runNum==476 && event==726) TroubleEvent=true;
			// if(runNum==553 && event==1) TroubleEvent=true;
			// if(runNum==704 && event==4855) TroubleEvent=true;
			// if(runNum==711 && event==0) TroubleEvent=true;
			// if(runNum==711 && event==0) TroubleEvent=true;
			// if(runNum==808 && event==9238) TroubleEvent=true;

			// config 3
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

			// config 4
			if(runNum==6684 && event==1) TroubleEvent=true;

			// config 5
			// if(runNum==2001 && event==0) TroubleEvent=true;
			// if(runNum==2466 && event==2232) TroubleEvent=true;
			// if(runNum==2466 && event==2232) TroubleEvent=true;
			// if(runNum==2472 && event==41) TroubleEvent=true;

			// frac_power_75[i];
			// frac_above_850[i];

			for(int i=0; i<16; i++){
				distro[i]->Fill(frac_power_75[i]);
				distro_above[i]->Fill(frac_above_850[i]);
				if(TroubleEvent){
					glitch[i]->Fill(frac_power_75[i]);
					glitch_above[i]->Fill(frac_above_850[i]);
					printf("Trouble event %4d: ch %2d frac below 75 %1.2f \n", event,i,frac_power_75[i]);
				}
				if(isCal){
					cal[i]->Fill(frac_power_75[i]);
					cal_above[i]->Fill(frac_above_850[i]);
				}

				power_oob_vs_ib[i]->Fill(between[i],power_75[i]+above_850[i]);
				if(TroubleEvent){
					TroubleEventOverlay[i]->SetPoint(numTroubleFound,between[i],power_75[i]+above_850[i]);
					numTroubleFound++;
					printf("Trouble event %4d: ch %2d frac below 75 %1.2f \n", event,i,frac_power_75[i]);
				}

				// power_oob_vs_ib[i]->Fill(10.*TMath::Log10(between[i]),10.*TMath::Log10(power_75[i]+above_850[i]));
				// if(TroubleEvent){
				// 	TroubleEventOverlay[i]->SetPoint(numTroubleFound,10.*TMath::Log10(between[i]),10.*TMath::Log10(power_75[i]+above_850[i]));
				// 	numTroubleFound++;
				// }
			}

			// for(double power=0.; power<1.; power+=0.01){
			// 	int this_glitch_counter=0;
			// 	for(int i=0; i<16; i++){
			// 		// if(frac_power_75[i]>power && i!=3 && i!=7 && i!=11 && i!=15){
			// 		if(frac_power_75[i]>power){
			// 			this_glitch_counter++;
			// 		}
			// 	}
			// 	distro_2d->Fill(power,this_glitch_counter);
			// }
		}
		fpIn->Close();
		delete fpIn;
	} //end loop over input files
	distro_2d->Scale(1/double(num_total));
	// for(int i=0; i<16; i++) power_oob_vs_ib[i]->Scale(1./power_oob_vs_ib[i]->GetMaximum());

	// for(int i=0; i<distro_2d->GetNbinsX(); i++){
	// 	cout<<"Bin "<<i<<" center is "<<distro_2d->GetXaxis()->GetBinCenter(i)<<" and width is "<<distro_2d->GetXaxis()->GetBinWidth(i)<<endl;
	// }

	stringstream ss1;
	vector<string> titlesForGraphs;
	for (int i = 0; i < 16; i++){
		ss1.str("");
		ss1 << "Channel " << i;
		titlesForGraphs.push_back(ss1.str());
	}
	
	TCanvas *c = new TCanvas("","",2*1100,2*850);
	c->Divide(4,4);
	for(int i=0; i<16; i++){
		c->cd(i+1);
		distro[i]->Draw("");
		distro[i]->SetLineWidth(3);
		distro[i]->SetLineColor(kBlue);
		cal[i]->Draw("same");
		cal[i]->SetLineWidth(3);
		cal[i]->SetLineColor(kGreen);
		glitch[i]->Draw("same");
		glitch[i]->SetLineWidth(3);
		glitch[i]->SetLineColor(kRed);
		gPad->SetLogy();
		distro[i]->GetYaxis()->SetRangeUser(0.1,1e7);
		distro[i]->GetYaxis()->SetTitle("Number of Events");
		distro[i]->GetXaxis()->SetTitle("Fraction of Power Below 75 MHz");
		distro[i]->SetTitle(titlesForGraphs[i].c_str());
		if(i+1==1){
			TLegend *leg = new TLegend(0.48,0.6,0.9,0.9);
			leg->AddEntry(distro_above[i],"All Events (no soft, no bad qual)","l");
			leg->AddEntry(cal_above[i],"Tagged Cal Pulsers","l");
			leg->AddEntry(glitch_above[i],"Straggling Events","l");
			leg->Draw();
		}
	}
	char save_plot_title[400];
	sprintf(save_plot_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/%d.%d.%d_Distro_PowerBelow75_A%d_c%d_%dEvents.png",year_now,month_now,day_now,station,config,num_total);
	c->SaveAs(save_plot_title);
	delete c;

	TCanvas *c4 = new TCanvas("","",2*1100,2*850);
	c4->Divide(4,4);
	for(int i=0; i<16; i++){
		c4->cd(i+1);
		distro_above[i]->Draw("");
		distro_above[i]->SetLineWidth(3);
		distro_above[i]->SetLineColor(kBlue);
		cal_above[i]->Draw("same");
		cal_above[i]->SetLineWidth(3);
		cal_above[i]->SetLineColor(kGreen);
		glitch_above[i]->Draw("same");
		glitch_above[i]->SetLineWidth(3);
		glitch_above[i]->SetLineColor(kRed);
		gPad->SetLogy();
		distro_above[i]->GetYaxis()->SetRangeUser(0.1,1e7);
		distro_above[i]->GetYaxis()->SetTitle("Number of Events");
		// distro[i]->GetXaxis()->SetTitle("Fraction of Power Below 75 MHz");
		distro_above[i]->GetXaxis()->SetTitle("Fraction of Power Above 850 MHz");
		distro_above[i]->SetTitle(titlesForGraphs[i].c_str());
		if(i+1==1){
			TLegend *leg = new TLegend(0.48,0.6,0.9,0.9);
			leg->AddEntry(distro_above[i],"All Events (no soft, no bad qual)","l");
			leg->AddEntry(cal_above[i],"Tagged Cal Pulsers","l");
			leg->AddEntry(glitch_above[i],"Straggling Events","l");
			leg->Draw();
		}
	}
	sprintf(save_plot_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/%d.%d.%d_Distro_PowerAbove850_A%d_c%d_%dEvents.png",year_now,month_now,day_now,station,config,num_total);
	c4->SaveAs(save_plot_title);
	delete c4;

	printf("Size of overlay is %d \n", TroubleEventOverlay[0]->GetN());

	gStyle->SetOptStat(1111111);
	TCanvas *c3 = new TCanvas("","",2*1100,2*850);
	c3->Divide(4,4);
	for(int i=0; i<16; i++){
		c3->cd(i+1);
		power_oob_vs_ib[i]->Draw("colz");
		gPad->SetLogz();
		TroubleEventOverlay[i]->Draw("Psame");
		TroubleEventOverlay[i]->SetMarkerSize(20);
		// gPad->SetLogy();
		// gPad->SetLogx();
		power_oob_vs_ib[i]->GetXaxis()->SetTitle("Power In Band");
		power_oob_vs_ib[i]->GetYaxis()->SetTitle("Power Out of Band");
		power_oob_vs_ib[i]->GetZaxis()->SetTitle("Fraction of Events");
		gPad->SetRightMargin(0.15);
		// power_oob_vs_ib[i]->GetZaxis()->SetRangeUser(1e-3,1.);
		double min = power_oob_vs_ib[i]->GetMinimum();
		double max= power_oob_vs_ib[i]->GetMaximum();
		// power_oob_vs_ib[i]->GetZaxis()->SetRangeUser(1e-1,1.);
	}
	sprintf(save_plot_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/%d.%d.%d_Distro_Glitch_PowerOOBvsIB_A%d_c%d_%dEvents.png",year_now,month_now,day_now,station,config,num_total);
	c3->SaveAs(save_plot_title);
	delete c3;

	// gStyle->SetOptStat(0);
	// TCanvas *c2 = new TCanvas("","",2*1100,2*850);
	// c2->SetGrid();
	// distro_2d->Draw("colz");
	// distro_2d->GetYaxis()->SetNdivisions(16,0,0);
	// distro_2d->GetYaxis()->SetTitle("Number of Channels");
	// distro_2d->GetXaxis()->SetTitle("Fraction of Power Below 75 MHz");
	// distro_2d->GetZaxis()->SetTitle("Fraction of Events");
	// distro_2d->GetZaxis()->SetTitleOffset(1.2);
	// gPad->SetLogz();
	// gPad->SetRightMargin(0.15);
	// distro_2d->GetZaxis()->SetRangeUser(1e-8,1.);		
	// sprintf(save_plot_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/%d.%d.%d_2D_Distro_Glitch_A%d_c%d_%dEvents.png",year_now,month_now,day_now,station,config,num_total);
	// c2->SaveAs(save_plot_title);
	// delete c2;
}