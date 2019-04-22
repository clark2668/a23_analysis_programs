///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////		2019.02.25_OffsetBlock-Demo.cxx
//////		June 2018, Brian Clark 
//////		Demo the offset block problem
////////////////////////////////////////////////////////////////////////////////

//C++
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

//ROOT Includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TTree.h"
#include "TStyle.h"

#include "FFTtools.h"

//AraRoot includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraQualCuts.h"

#include "tools_PlottingFns.h"
#include "tools_inputParameters.h"
#include "tools_runSummaryObjects.h"
#include "tools_outputObjects.h"
#include "tools_Cuts.h"

using namespace std;

double MaxMeanBlock(TGraph *grIn, int evt_num, int chan, vector<double> &means, bool print);

int main(int argc, char **argv)
{
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	stringstream ss;
	gStyle->SetOptStat(11);
	
	if(argc<1){
		cout<< "Usage\n" << argv[0] << " <data_file> <ped_file>"<<endl;
		return 0;
	}
	TH1D *distro[16];
	TH1D *distro_cal[16];
	TH1D *glitch[16];
	for(int i=0; i<16; i++){
		distro[i] = new TH1D("","",100,-600,600);
		distro_cal[i] = new TH1D("","",100,-600,600);
		glitch[i] = new TH1D("","",100,-600,600);
	}

	//now we need to get the data
	
	TFile *fpIn = TFile::Open(argv[1]);
	if(!fpIn){
		cout<<"Can't open data file for map!"<<endl;
		return -1;
	}
	TTree *eventTree = (TTree*) fpIn-> Get("eventTree");
	if(!eventTree){
		cout<<"Can't find eventTree for map"<<endl;
		return -1;
	}
	RawAtriStationEvent *rawPtr =0;
	int runNum;
	eventTree->SetBranchAddress("event",&rawPtr);
	eventTree->SetBranchAddress("run",&runNum);
	eventTree->GetEvent(0);
	int stationID = rawPtr->stationId;


	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	calibrator->setAtriPedFile(argv[2], stationID);
	AraQualCuts *qual = AraQualCuts::Instance();

	int numEntries = eventTree->GetEntries();
	cout<<"Num entries is "<<numEntries<<endl;
	Long64_t starEvery=numEntries/100;
	if(starEvery==0) starEvery++;
	numEntries=4;
	for(int event=3; event<numEntries; event++){

		if(event%starEvery==0) {
			std::cout << "	On event "<<event<<endl;
		}

		eventTree->GetEvent(event);
		UsefulAtriStationEvent *ev = new UsefulAtriStationEvent(rawPtr,AraCalType::kLatestCalib);
		if(!(qual->isGoodEvent(ev))){
			delete ev;
			continue;
		}
		cout<<"Event is "<<ev->eventNumber<<endl;

		//now get the waveforms
		stringstream ss1;
		string xLabel, yLabel;
		xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
		vector<string> titlesForGraphs;
		for (int i = 0; i < 16; i++){
			ss1.str("");
			ss << "Channel " << i;
			titlesForGraphs.push_back(ss1.str());
		}
		vector <TGraph*> grWaveformsRaw = makeGraphsFromRF(ev,16,xLabel,yLabel,titlesForGraphs);
		int has_short = hasShortWaveformMultiGraph(grWaveformsRaw);
		if(has_short>0){
			deleteGraphVector(grWaveformsRaw);
			delete ev;
			continue;
		}


		//now to do glitch detection
		double this_max_chan=0;
		for(int chan=0; chan<16; chan++){
			vector<double> means;
			double this_max = MaxMeanBlock(grWaveformsRaw[chan], event, chan, means, true);
			// if(chan==3) this_max_chan=this_max;
			if(!ev->isCalpulserEvent() && !ev->isSoftwareTrigger()){
				distro[chan]->Fill(this_max);
				// for(int block=0; block<means.size(); block++){
				// 	distro[chan]->Fill(means[block]);
				// }
			}


			// if(ev->isSoftwareTrigger()) distro_cal[chan]->Fill(TMath::Log10(this_max));
			// if(event==40) glitch[chan]->Fill(TMath::Log10(this_max));
		}
		// if(TMath::Log10(this_max_chan)<1.1) plot=true;
		// bool plot=false;
		// if(plot){
		// 	TCanvas *c = new TCanvas("","",2*1100,2*850);
		// 	c->Divide(4,4);
		// 	for(int i=0; i<16; i++){
		// 		c->cd(i+1);
		// 		grWaveformsRaw[i]->Draw("al");
		// 		grWaveformsRaw[i]->SetLineWidth(2);
		// 	}
		// 	char save_plot_title[400];
		// 	sprintf(save_plot_title,"/users/PAS0654/osu0673/A23_analysis/results/glitch_detect/%d.%d.%d_Glitch_Run%d_Ev%d.png",year_now,month_now,day_now,runNum,event);
		// 	c->SaveAs(save_plot_title);
		// 	delete c;
		// }
		deleteGraphVector(grWaveformsRaw);
		delete ev;

	}
	
	gStyle->SetOptStat(111111);
	TCanvas *c = new TCanvas("","",1100,850);
	c->Divide(4,4);
	for(int i=0; i<16; i++){
		c->cd(i+1);
		// distro[i]->Draw();
		// distro_cal[i]->GetXaxis()->SetTitle("Log10( Highest Abs Maximum of Mean Among Blocks)");
		// distro_cal[i]->GetYaxis()->SetTitle("Number of Events");
		// distro[i]->SetLineWidth(1);
		// distro[i]->SetLineColor(kBlack);
		// distro_cal[i]->Draw();
		distro[i]->Draw();
		distro[i]->SetLineColor(kBlack);
		// distro_cal[i]->SetLineWidth(1);
		// distro_cal[i]->SetLineColor(kGreen);
		// distro_cal[i]->SetLineStyle(9);
		// glitch[i]->Draw("same");
		// glitch[i]->SetLineWidth(1);
		// glitch[i]->SetLineColor(kGreen);
		// glitch[i]->SetLineStyle(10);
		gPad->SetLogy();
		// gPad->SetLogx();
		// distro_cal[i]->GetYaxis()->SetRangeUser(0.1,2000.);
		// distro[i]->GetXaxis()->SetRangeUser(0.01,3000.);
	}
	char save_plot_title[400];
	sprintf(save_plot_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/%d.%d.%d_GlitchRemoval_Distros.pdf",year_now,month_now,day_now);
	c->SaveAs(save_plot_title);
	delete c;

	fpIn->Close();
	
}

double MaxMeanBlock(TGraph *grIn, int evt_num, int chan, vector<double> &means, bool print){

	int n_input = grIn->GetN();
	double *oldX = grIn->GetX();
	double *oldY = grIn->GetY();

	deque <double> inX;
	deque <double> inY;

	vector <double> newX;
	vector <double> newY;

	double totalYValues=0.;
	double originalmean = 0.;

	for(int samp=0; samp<n_input; samp++){
		inX.push_back(oldX[samp]);
		inY.push_back(oldY[samp]);
		newX.push_back(oldX[samp]);
		newY.push_back(oldY[samp]);
		totalYValues+=oldY[samp];
	}
	printf("Total y values in channel %d is %.2f/%d = %.2f \n", chan, totalYValues, n_input, totalYValues/double(n_input));
	originalmean = totalYValues/double(n_input);
	//do mean subtraction
	double totalYValues_aftersubtract=0;
	for(int samp=0; samp<inX.size(); samp++){
		newY[samp]-=originalmean;
		totalYValues_aftersubtract+=newY[samp];
	}
	TGraph *grNew = new TGraph(newY.size(), &newX[0], &newY[0]);
	printf("Total y values after subtract in channel %d is %.2f/%d = %.2f \n", chan, totalYValues_aftersubtract, n_input, totalYValues_aftersubtract/double(n_input));


	double newtotalYValues=0.;
	int numfornewtotal=0.;
	double first_time=inX[0];
	for(int samp=0; samp<n_input; samp++){
		if(inX[samp]>first_time+20.){
			newtotalYValues+=inY[samp];
			numfornewtotal++;
		}
	}
	printf("Total y values in channel %d is %.2f/%d = %.2f  \n", chan, newtotalYValues, numfornewtotal, newtotalYValues/double(numfornewtotal));

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
		means.push_back(abs(mean_this_section));
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
		grIn->GetYaxis()->SetRangeUser(-120.,120.);
		grIn->SetLineWidth(1);
		for(int i=0; i<grPieces.size(); i++){
			grPieces[i]->Draw("same");
			grPieces[i]->SetLineColor(colors[i]);
			grPieces[i]->SetLineWidth(3);
		}
		char title[400];
		sprintf(title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/GlitchWaveform_ev%d_chan%d.png",evt_num,chan);
		c->SaveAs(title);	
		delete c;

		TGraph *interpolated_waveform = FFTtools::getInterpolatedGraph(grNew, 0.6); //get an interpolated waveform with 0.5 ns interpolation
		TGraph *pad = FFTtools::padWaveToLength(interpolated_waveform,2048);
		TGraph *spectrum = FFTtools::makePowerSpectrumMilliVoltsNanoSecondsdB(pad); //now make a spectrum
		TCanvas *c2 = new TCanvas("","",1100,850);
				spectrum->Draw("AL");
				spectrum->SetLineWidth(1);
				spectrum->GetYaxis()->SetTitle("Power (dB)");
				spectrum->GetXaxis()->SetTitle("Frequency (MHz)");
				spectrum->GetYaxis()->SetRangeUser(0,40);
				// spectrum->GetXaxis()->SetRangeUser(0,200);
				// gPad->SetLogx();
		sprintf(title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/GlitchWaveformSpec_ev%d_chan%d.png",evt_num,chan);
		c2->SaveAs(title);
		delete c2;
		delete spectrum;
		delete pad;
		delete interpolated_waveform;
	}
	for(int i=0; i<grPieces.size(); i++){
		delete grPieces[i];
	}

	for(int i=0; i<means.size(); i++){
		printf("Chan %d Block %d has mean %.2f \n", chan, i, means[i]);
	}
	std::sort(means.begin(), means.end(), std::greater<double>());
	// printf("Biggest value is %.2f \n", means[0]);
	return means[0];
}