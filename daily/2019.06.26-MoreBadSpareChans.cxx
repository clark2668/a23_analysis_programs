///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////		2018.12.16_Glitch_Detect_Distro.cxx
//////		June 2018, Brian Clark 
//////		Distributions to do glitch detection
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
#include "tools_CommandLine.h"

using namespace std;

double MyFunc(TGraph *grIn, int evt_num, int chan, vector<double> &means, vector<double> &stdDev, vector<double> MeanOverstdDev);
void MyFunc2(TGraph *grIn, vector<double> &means);

int main(int argc, char **argv)
{
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	stringstream ss;
	gStyle->SetOptStat(11);
	
	if(argc!=3){
		cout<< "Usage\n" << argv[0] << " <data_file> <ped_file>"<<endl;
		return 0;
	}
	TH1D *distro[4];
	TH1D *distro_cal[4];
	TH1D *glitch[4];
	for(int i=0; i<4; i++){
		distro[i] = new TH1D("","",100,0,100);
		distro_cal[i] = new TH1D("","",100,0,100);
		glitch[i] = new TH1D("","",100,0,100);
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
	numEntries=7;
	for(int event=5; event<numEntries; event++){

		if(event%starEvery==0) {
			std::cout << "	On event "<<event<<endl;
		}

		eventTree->GetEvent(event);
		UsefulAtriStationEvent *ev = new UsefulAtriStationEvent(rawPtr,AraCalType::kLatestCalib);
		// if(!(qual->isGoodEvent(ev))){
		// 	delete ev;
		// 	continue;
		// }
		// cout<<"On event"<<event<<endl;

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
		// int has_short = hasShortWaveformMultiGraph(grWaveformsRaw);
		// if(has_short>0){
		// 	// cout<<"Has short!"<<endl;
		// 	deleteGraphVector(grWaveformsRaw);
		// 	delete ev;
		// 	continue;
		// }

		vector<TGraph*> electChansGraphs;
		electChansGraphs.push_back(ev->getGraphFromElecChan(6));
		electChansGraphs.push_back(ev->getGraphFromElecChan(14));
		electChansGraphs.push_back(ev->getGraphFromElecChan(22));
		electChansGraphs.push_back(ev->getGraphFromElecChan(30));

		vector<vector< double > > means;
		// means.resize(4);
		int shortest=300;
		for(int i=0; i<4; i++){
			vector<double> theseMeans;
			MyFunc2(electChansGraphs[i], theseMeans);
			// printf("For elec chan %d, found %d means \n",i,theseMeans.size());
			if(theseMeans.size()<shortest){
				shortest=theseMeans.size();
			}
			for(int j=0; j<theseMeans.size(); j++){
				// printf("		For elec chan %d, means are %.2f \n", i, theseMeans[j]);
			}
			means.push_back(theseMeans);
		}
		// printf("Shortest now is %d \n", shortest);
		for(int i=0; i<4; i++){
			// printf("For chan %d, length is %d \n",i, means[i].size());
			while(means[i].size()>shortest){
				means[i].pop_back();
			}
		}
		for(int i=0; i<4; i++){
			for(int block=0; block<means[i].size(); block++){
				distro[i]->Fill(abs(means[i][block]));
				if(runNum==3663 && event==6 && block<4){
					// cout<<"Yup, filling this!"<<endl;
					glitch[i]->Fill(abs(means[i][block]));
				}
			}
			// printf("Now, for channel %d, I have %d means \n",i, means[i].size());
		}

		int numViolatingBlocks=0;
		for(int block=0; block<shortest; block++){
			int numViolating=0;
			for(int chan=0; chan<4; chan++){
				if(abs(means[chan][block])>20){
					numViolating++;
				}
			}
			if(numViolating>1){
				numViolatingBlocks++;
			}
			if(numViolatingBlocks>1){
				// that's enough for the cut, be done
				break;
			}
		}

		bool hasViolation=false;
		if(numViolatingBlocks>1){
			hasViolation=true;
		}

		if(hasViolation){
			printf(RED"Event %d has a coincident offset block violation with %d violating blocks!\n"RESET,event,numViolatingBlocks);
		}



		// int nChanBelowThresh=0;
		// vector<double> maxTimeVec;
		// for(int i=0; i<4; i++){
		// 	printf("Num samples in chan %d is %d \n", i, electChansGraphs[i]->GetN());
		// 	TGraph *grInt = FFTtools::getInterpolatedGraph(electChansGraphs[i],0.4);
		// 	TGraph *grMean = qual->getRollingMean(grInt,64);
		// 	double maxTime;
		// 	double meanMax = qual->getMax(grMean,&maxTime);
		// 	printf("Channel %d meanMax is %.2f at time %.2f \n", i, meanMax, maxTime);
		// 	if(-1.*fabs(meanMax)<-20.){
		// 		nChanBelowThresh++;
		// 		maxTimeVec.push_back(maxTime);
		// 	}
		// 	delete grInt;
		// 	delete grMean;
		// }
		// if(nChanBelowThresh>1){
		// 	double timeRange = *max_element(maxTimeVec.begin(), maxTimeVec.end()) - *min_element(maxTimeVec.begin(), maxTimeVec.end());
		// 	printf("Timee range is timeRange %.2f \n", timeRange);
		// }

		for(int i=0; i<4; i++){
			vector<double> means;
			vector<double> stdDevs;
			vector<double> MeanOverstdDev;
		

			// double thisMax = MyFunc(electChansGraphs[i],event,i,means,stdDevs,MeanOverstdDev);
			// // if(!ev->isCalpulserEvent() && !ev->isSoftwareTrigger()){
			// if(!ev->isSoftwareTrigger()){
			// 	distro[i]->Fill(thisMax);
			// }
			// if(runNum==3663 && event==6){
			// 	glitch[i]->Fill(thisMax);
			// }
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
		deleteGraphVector(electChansGraphs);
		delete ev;

	}
	
	gStyle->SetOptStat(111111);
	TCanvas *c = new TCanvas("","",4*850,850);
	c->Divide(4,1);
	for(int i=0; i<4; i++){
		c->cd(i+1);
		// distro[i]->Draw();
		// distro_cal[i]->GetXaxis()->SetTitle("Log10( Highest Abs Maximum of Mean Among Blocks)");
		// distro_cal[i]->GetYaxis()->SetTitle("Number of Events");
		// distro[i]->SetLineWidth(1);
		// distro[i]->SetLineColor(kBlack);
		// distro_cal[i]->Draw();
		distro[i]->Draw();
		distro[i]->SetLineColor(kBlack);
		// distro_cal[i]->SetLineWidth(3);
		// distro_cal[i]->SetLineColor(kGreen);
		// distro_cal[i]->SetLineStyle(9);
		glitch[i]->Draw("same");
		glitch[i]->SetLineWidth(3);
		glitch[i]->SetLineColor(kRed);
		glitch[i]->SetLineStyle(10);
		gPad->SetLogy();
		// gPad->SetLogx();
		// distro_cal[i]->GetYaxis()->SetRangeUser(0.1,2000.);
		// distro[i]->GetXaxis()->SetRangeUser(0.01,3000.);
	}
	char save_plot_title[400];
	sprintf(save_plot_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/%d.%d.%d_GlitchRemoval_Distros_AllSpareChanBlockMeans.png",year_now,month_now,day_now);
	c->SaveAs(save_plot_title);
	delete c;

	fpIn->Close();
	
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

		// double forStdDev=0.;
		// for(int samp=0; samp<sub_X.size(); samp++){
		// 	double thisY = sub_Y[samp];
		// 	forStdDev+=pow(thisY - mean_this_section,2.);
		// }
		// forStdDev/=double(sub_X.size());
		// forStdDev=sqrt(forStdDev);
		double forStdDev=6.;

		means.push_back(abs(mean_this_section));
		stdDev.push_back(forStdDev);	
		MeanOverstdDev.push_back(abs(mean_this_section/(forStdDev/sqrt(double(sub_X.size())))));
		// MeanOverstdDev.push_back(abs(mean_this_section/(forStdDev/double(sub_X.size()))));

		for(int iPop=0; iPop<num_to_pop; iPop++){
			inX.pop_front();
			inY.pop_front();
		}
	}

	for(int i=0; i<means.size(); i++){
		// printf("Block %d has mean %.2f and sigma %.2f and mean over sigma thing of %.3f \n", i, means[i], stdDev[i], MeanOverstdDev[i]);
	}
	std::sort(MeanOverstdDev.begin(), MeanOverstdDev.end(), std::greater<double>());
	std::sort(means.begin(), means.end(), std::greater<double>());
	// return MeanOverstdDev[0];
	return means[0];
}