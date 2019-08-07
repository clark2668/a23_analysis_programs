////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	cuts_refine_cal_box
////	Reads recco_save files and refines the cal puler cut box
////
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
#include "TF1.h"
#include "TLine.h"

//AraRoot includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "Settings.h"
#include "Detector.h"
#include "Report.h"
#include "RayTraceCorrelator.h"
#include "AraGeomTool.h"
#include "TTimeStamp.h"
AraAntPol::AraAntPol_t Vpol = AraAntPol::kVertical;
AraAntPol::AraAntPol_t Hpol = AraAntPol::kHorizontal;

//custom analysis includes
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_Cuts.h"
#include "tools_CW.h"
#include "tools_outputObjects.h"

bool doRezero=0;

int PlotThisEvent(int station, int config, int runNum, int event, int problempol, Settings *settings, Detector *detector, RayTraceCorrelator *theCorrelators[2]);

char *DataDirPath(getenv("DATA_DIR"));
char *PedDirPath(getenv("PED_DIR"));
char *plotPath(getenv("PLOT_PATH"));

using namespace std;

int main(int argc, char **argv)
{
	gStyle->SetOptStat(1111);
	gStyle->SetOptStat(0);
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <station> <config> <pulser> <pol> <reco_val_save_file_1> <reco_val_save_file_2 > ... <reco_val_save_file_x>"<<endl;
		return 0;
	}
	int station = atoi(argv[1]);
	int config = atoi(argv[2]);
	int pulser = atoi(argv[3]);
	int pol = atoi(argv[4]);

	//set up the ray tracer
	Settings *settings = new Settings();
	string setupfile = "setup.txt";
	settings->ReadFile(setupfile);
	cout << "Read " << setupfile << " file!" << endl;
	settings->NOFZ=1;
	Detector *detector=0;
	RayTraceCorrelator *theCorrelators[2];
	theCorrelators[0] =  new RayTraceCorrelator(station, 41., settings, 1, 4); //41 m, cal puser
	theCorrelators[1] =  new RayTraceCorrelator(station, 300., settings, 1, 4);//300 m, far reco

	stringstream ss;

	TH2D *locations[2];
	locations[0]=new TH2D("","V",360,-180,180,180,-90,90);
	locations[1]=new TH2D("","H",360,-180,180,180,-90,90);

	TH2D *locations_zoom[2];
	locations_zoom[0]=new TH2D("","V",360,-180,180,180,-90,90);
	locations_zoom[1]=new TH2D("","H",360,-180,180,180,-90,90);

	TH1D *locations_phi_slice[2];
	locations_phi_slice[0]=new TH1D("","V",360,-180,180);
	locations_phi_slice[1]=new TH1D("","H",360,-180,180);

	TH1D *locations_theta_slice[2];
	locations_theta_slice[0]=new TH1D("","V",360,-180,180);
	locations_theta_slice[1]=new TH1D("","H",360,-180,180);

	TTimeStamp start;
	TTimeStamp stop;
	start.Set(2013, 01, 00, 00, 00,0,0,true,0);
	stop.Set(2016, 12, 31, 24, 00,0,0,true,0);

	int start_bin = start.GetSec();
	int stop_bin = stop.GetSec();

	TH2D *phi_vs_time[2];
	TH2D *theta_vs_time[2];
	int num_bins=(stop_bin-start_bin)/60/60;
	phi_vs_time[0] = new TH2D("","V Phi Vs Time",num_bins, start_bin, stop_bin, 360,-180,180);
	theta_vs_time[0] = new TH2D("","V Theta Vs Time",num_bins, start_bin, stop_bin, 180,-90,90);
	phi_vs_time[1] = new TH2D("","H Phi Vs Time",num_bins, start_bin, stop_bin, 360,-180,180);
	theta_vs_time[1] = new TH2D("","H Theta Vs Time",num_bins, start_bin, stop_bin, 180,-90,90);
	for(int i=0; i<2; i++){
		phi_vs_time[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		phi_vs_time[i]->GetXaxis()->SetTimeFormat("%y/%m");
		theta_vs_time[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		theta_vs_time[i]->GetXaxis()->SetTimeFormat("%y/%m");
		phi_vs_time[i]->GetXaxis()->SetTimeOffset(0.,"GMT");
		theta_vs_time[i]->GetXaxis()->SetTimeOffset(0.,"GMT");
	}

	int num_total=0;

	vector<int> BadRunList=BuildBadRunList(station);

	for(int file_num=5; file_num<argc; file_num++){

		cout << "Run " << file_num << " :: " << argv[file_num] << endl;

		string chRun = "run";
		string file = string(argv[file_num]);
		size_t foundRun = file.find(chRun);
		string strRunNum = file.substr(foundRun+4,4);
		int runNum = atoi(strRunNum.c_str());

		int isThisBadABadRun = isBadRun(station,runNum,BadRunList);

		if(isThisBadABadRun) continue;
		if(runNum==2428 || runNum==3442 || runNum==7100) continue;
		if(runNum>=8100) continue;

		TFile *inputFile = TFile::Open(argv[file_num]);
		if(!inputFile){
			cout<<"Can't open joined file!"<<endl;
			return -1;
		}
		TTree *trees = (TTree*) inputFile->Get("RecoVals");
		if(!trees){
			cout<<"Can't open reco tree"<<endl;
			return -1;
		}
		double phi_41_V;
		double theta_41_V;
		double phi_41_H;
		double theta_41_H;
		double corr_val_V;
		double corr_val_H;
		int isCal;
		int isBad;
		int unixTime;
		trees->SetBranchAddress("phi_41_V",&phi_41_V);
		trees->SetBranchAddress("theta_41_V",&theta_41_V);
		trees->SetBranchAddress("phi_41_H",&phi_41_H);
		trees->SetBranchAddress("theta_41_H",&theta_41_H);
		trees->SetBranchAddress("corr_val_V",&corr_val_V);
		trees->SetBranchAddress("corr_val_H",&corr_val_H);
		trees->SetBranchAddress("cal",&isCal);
		trees->SetBranchAddress("isBad",&isBad);
		trees->SetBranchAddress("unixTime",&unixTime);


		int numEntries = trees->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;

		//now to loop over events
		for(int event=0; event<numEntries; event++){

			trees->GetEvent(event);

			if(!isCal || isBad) continue;
			if(isBadLivetime(station,unixTime)){
				continue;
			}

			num_total++;
			if(corr_val_V>corr_val_H){
				locations[0]->Fill(phi_41_V,theta_41_V);
				locations_zoom[0]->Fill(phi_41_V,theta_41_V);
				locations_phi_slice[0]->Fill(phi_41_V);
				locations_theta_slice[0]->Fill(theta_41_V);
				phi_vs_time[0]->Fill(unixTime, phi_41_V);
				theta_vs_time[0]->Fill(unixTime, theta_41_V);

				double this_phi = phi_41_V;
				double this_theta = theta_41_V;
				// if(this_theta<-30 && this_phi<50){
				// if(this_phi>50 && this_theta<0){
				if(this_phi>100){
				// if(this_theta<20 && this_theta>-10 && this_phi<0){
				// if(this_theta>20 && this_phi<0 && this_phi>-50){
				// if(this_phi<-40){
				// if(this_phi>-50 && this_phi<0 && this_theta<-30){
				  //PlotThisEvent(station, config, runNum, event, 0, settings, detector, theCorrelators);
				}

			}
			else if(corr_val_H>corr_val_V){
				locations[1]->Fill(phi_41_H,theta_41_H);
				locations_zoom[1]->Fill(phi_41_H,theta_41_H);
				locations_phi_slice[1]->Fill(phi_41_H);
				locations_theta_slice[1]->Fill(theta_41_H);
				phi_vs_time[1]->Fill(unixTime, phi_41_H);
				theta_vs_time[1]->Fill(unixTime, theta_41_H);
			}
		}//loop over events
		inputFile->Close();
		delete inputFile;
	} //end loop over input files
	
	char graph_title[2][300];
	char title[300];

	sprintf(graph_title[0],"VPol");
	sprintf(graph_title[1],"HPol");
	TCanvas *c = new TCanvas("","",2.1*850,850);
	c->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		c->cd(pol+1);
		locations[pol]->Draw("colz");
		locations[pol]->GetYaxis()->SetTitle("Theta (deg)");
		locations[pol]->GetXaxis()->SetTitle("Phi (deg)");
		locations[pol]->SetTitle(graph_title[pol]);
		gPad->SetLogz();
	}
	sprintf(title, "%s/cal_cut/%d.%d.%d_A%d_c%d_%dEvents_DistroCalPulses.png",plotPath,year_now, month_now, day_now,station,config,num_total);
	c->SaveAs(title);
	delete c;
	delete locations[0]; delete locations[1];



	int projectX_low, projectX_high, projectY_low, projectY_high;
	getCalCutPlotBoundary(station,config,pulser,pol,false,projectX_low,projectX_high,projectY_low,projectY_high);


		projectX_low=-35;
		projectX_high=-15;
		projectY_low=20;
		projectY_high=50;


	double real_phi, real_theta;
	cout<<"Pulser is "<<pulser<<endl;
	getRealLocation(station,pulser,pol,real_theta,real_phi);
	printf("Theta, Phi %.2f %.2f\n", real_theta, real_phi);
	TLine theta_line_2D(projectX_low,double(real_theta),projectX_high,double(real_theta));
	TLine phi_line_2D(real_phi,projectY_low,real_phi,projectY_high);

	TCanvas *c2 = new TCanvas("","",3*850,850);
	c2->Divide(3,1);
	c2->cd(1);
		locations_zoom[pol]->Draw("colz");
		locations_zoom[pol]->GetYaxis()->SetTitle("Theta (deg)");
		locations_zoom[pol]->GetXaxis()->SetTitle("Phi (deg)");
		locations_zoom[pol]->SetTitle(graph_title[0]);
		gPad->SetLogz();
		locations_zoom[pol]->GetXaxis()->SetRangeUser(projectX_low,projectX_high);
		locations_zoom[pol]->GetYaxis()->SetRangeUser(projectY_low,projectY_high);

	TH2D *locations_zoom_clone[2];
	locations_zoom_clone[0] = (TH2D*) locations_zoom[0]->Clone();
	locations_zoom_clone[1] = (TH2D*) locations_zoom[1]->Clone();

	getCalCutPlotBoundary(station,config,pulser,pol,true,projectX_low,projectX_high,projectY_low,projectY_high);


		projectX_low=-35;
	projectX_high=-15;
	projectY_low=20;
		projectY_high=50;

		projectX_low+=180;
		projectX_high+=180;
		projectY_low+=90;
		projectY_high+=90;

	TH1D *project_V[2];
	project_V[0] = (TH1D*) locations_zoom_clone[pol]->ProjectionX("",projectY_low,projectY_high)->Clone();
	project_V[1] = (TH1D*) locations_zoom_clone[pol]->ProjectionY("",projectX_low,projectX_high)->Clone();
	project_V[0]->SetTitle("Phi Projection near peak");
	project_V[1]->SetTitle("Theta Projection near peak");
	Double_t scale = 1/project_V[0]->Integral();
	project_V[0]->Scale(scale);
	scale = 1/project_V[1]->Integral();
	project_V[1]->Scale(scale);

	double startPhiFit, stopPhiFit, startThetaFit, stopThetaFit;
	getCalFitRange(station, config, pulser,pol, startPhiFit, stopPhiFit, startThetaFit, stopThetaFit);

	char equation_phi[150];
	sprintf(equation_phi,"gaus");
	TF1 *fit_phi = new TF1("GausFit_phi",equation_phi,startPhiFit,stopPhiFit);
	project_V[0]->Fit("GausFit_phi");
	printf("Chi-Square/NDF %.2f / %.2f \n",fit_phi->GetChisquare(),double(fit_phi->GetNDF()));
	char fit_graph_title_phi[150];
	// sprintf(fit_graph_title_phi,"Vpol %.2f*exp(-0.5((x-%.2f)/%.2f)^2)     #chi^{2}/NDF = %.2f / %d \n",fit_phi->GetParameter(0), fit_phi->GetParameter(1),fit_phi->GetParameter(2),fit_phi->GetChisquare(),fit_phi->GetNDF());
	// sprintf(fit_graph_title_phi,"Vpol %.2f*exp(-0.5((x-%.2f)/%.2f)^2)\n",fit_phi->GetParameter(0), fit_phi->GetParameter(1),fit_phi->GetParameter(2));
	sprintf(fit_graph_title_phi,"Phi Fit #mu=%.2f, #sigma=%.2f. Expected: %.2f. ",fit_phi->GetParameter(1),fit_phi->GetParameter(2),real_phi);


	char equation_theta[150];
	sprintf(equation_theta,"gaus");
	TF1 *fit_theta= new TF1("GausFit_theta",equation_theta,startThetaFit,stopThetaFit);
	project_V[1]->Fit("GausFit_theta");
	printf("Chi-Square/NDF %.2f / %.2f \n",fit_theta->GetChisquare(),double(fit_theta->GetNDF()));
	char fit_graph_title_theta[150];
	// sprintf(fit_graph_title_theta,"Vpol %.2f*exp(-0.5((x-%.2f)/%.2f)^2)     #chi^{2}/NDF = %.2f / %d \n",fit_theta->GetParameter(0), fit_theta->GetParameter(1),fit_theta->GetParameter(2),fit_theta->GetChisquare(),fit_theta->GetNDF());
	// sprintf(fit_graph_title_theta,"Vpol %.2f*exp(-0.5((x-%.2f)/%.2f)^2)\n",fit_theta->GetParameter(0), fit_theta->GetParameter(1),fit_theta->GetParameter(2));  
	sprintf(fit_graph_title_theta,"Theta Fit #mu=%.2f, #sigma=%.2f. Expected %.2f ",fit_theta->GetParameter(1),fit_theta->GetParameter(2),real_theta);

	printf("Six sigma point phi %.2f \n",6.*fit_phi->GetParameter(2));
	printf("Six sigma points phi %.2f to %.2f \n",(fit_phi->GetParameter(1) - 6.*fit_phi->GetParameter(2)),(fit_phi->GetParameter(1) + 6.*fit_phi->GetParameter(2)));	
	printf("Six sigma point theta %.2f to %.2f \n",(fit_theta->GetParameter(1) - 6.*fit_theta->GetParameter(2)),(fit_theta->GetParameter(1) + 6.*fit_theta->GetParameter(2)));

	/*
		We live in ROOT's world where we now need to redraw things...
	*/

	locations_zoom[0]->Draw("colz");
	phi_line_2D.Draw("");
	theta_line_2D.Draw("");
	phi_line_2D.SetLineStyle(9);
	theta_line_2D.SetLineStyle(9);
	theta_line_2D.SetLineWidth(2);
	phi_line_2D.SetLineWidth(2);

	// double max_amp;
	// max_amp=project_V[0]->GetMaximum()*7;
	// if(project_V[1]->GetMaximum()*7 > max_amp)
	// 	max_amp=project_V[1]->GetMaximum()*7;

	c2->cd(2);
		project_V[0]->Draw();
		project_V[0]->GetYaxis()->SetTitle("Normalized Number of Events");
		project_V[0]->GetXaxis()->SetTitle("Phi (deg)");
		project_V[0]->SetTitle(fit_graph_title_phi);
		gPad->SetLogy();
		project_V[0]->GetYaxis()->SetTitleOffset(1.5);
		// project_V[0]->GetYaxis()->SetRangeUser(0.1,max_amp);
		// TLine phi_proj_line(real_phi,0.1,real_phi,max_amp);
		project_V[0]->GetYaxis()->SetRangeUser(1e-7,3);
		TLine phi_proj_line(real_phi,1e-7,real_phi,3);
		phi_proj_line.Draw("");
		phi_proj_line.SetLineStyle(9);
		phi_proj_line.SetLineWidth(2);
	c2->cd(3);
		project_V[1]->Draw();
		project_V[1]->GetYaxis()->SetTitle("Normalized Number of Events");
		project_V[1]->GetXaxis()->SetTitle("Theta (deg)");
		project_V[1]->SetTitle(graph_title[0]);
		project_V[1]->SetTitle(fit_graph_title_theta);
		gPad->SetLogy();
		project_V[1]->GetYaxis()->SetTitleOffset(1.5);
		// project_V[1]->GetYaxis()->SetRangeUser(0.1,max_amp);
		// TLine theta_proj_line(real_theta,0.1,real_theta,max_amp);		
		project_V[1]->GetYaxis()->SetRangeUser(1e-7,3);
		TLine theta_proj_line(real_theta,1e-7,real_theta,3);
		theta_proj_line.Draw("");
		theta_proj_line.SetLineStyle(9);
		theta_proj_line.SetLineWidth(2);
	sprintf(title, "%s/cal_cut/%d.%d.%d_A%d_c%d_%dEvents_DistroCalPulses_Pol%d_CP%d_Zoom.png",plotPath, year_now, month_now, day_now,station,config,num_total,pol,pulser);
	c2->SaveAs(title);
	delete c2;

	// for plotting reco tracking
	gStyle->SetOptStat(0);
	TCanvas *c_phi_vs_time = new TCanvas("","",2*1100,2*850);
	c_phi_vs_time->Divide(2,2);
	c_phi_vs_time->cd(1);
		phi_vs_time[0]->Draw("colz");
		phi_vs_time[0]->GetXaxis()->SetTitle("Unixtime");
		phi_vs_time[0]->GetYaxis()->SetTitle("Phi");
	c_phi_vs_time->cd(2);
		// phi_vs_time[1]->Draw("colz");
		phi_vs_time[1]->GetXaxis()->SetTitle("Unixtime");
		phi_vs_time[1]->GetYaxis()->SetTitle("Phi");
	c_phi_vs_time->cd(3);
		theta_vs_time[0]->Draw("colz");
		theta_vs_time[0]->GetXaxis()->SetTitle("Unixtime");
		theta_vs_time[0]->GetYaxis()->SetTitle("Theta");		
	c_phi_vs_time->cd(4);
		// theta_vs_time[1]->Draw("colz");
		theta_vs_time[1]->GetXaxis()->SetTitle("Unixtime");
		theta_vs_time[1]->GetYaxis()->SetTitle("Theta");
	sprintf(title, "%s/cal_cut/%d.%d.%d_A%d_c%d_%dEvents_RecoVsTime.png",plotPath,year_now, month_now, day_now,station,config,num_total);
	c_phi_vs_time->SaveAs(title); delete c_phi_vs_time;
	delete phi_vs_time[0]; delete phi_vs_time[1];	


}

int PlotThisEvent(int station, int config, int runNum, int event, int problempol, Settings *settings, Detector *detector, RayTraceCorrelator *theCorrelators[2]){
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	// beautify_TH2D();

	char run_file_name[400];
	sprintf(run_file_name,"%s/RawData/A%d/by_config/c%d/event%d.root",DataDirPath,station,config,runNum);
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
	printf("Unixtime is %d \n", unixTime);
	printf("Unixtime microsecond is %d \n", unixTimeUs);

	stringstream ss1;
	string xLabel, yLabel;
	xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
	vector<string> titlesForGraphs;
	for (int i = 0; i < 16; i++){
		ss1.str("");
		ss1 << "Channel " << i;
		titlesForGraphs.push_back(ss1.str());
	}
	vector <TGraph*> waveforms = makeGraphsFromRF(realAtriEvPtr,16,xLabel,yLabel,titlesForGraphs);
	vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(waveforms, 0.5, xLabel, yLabel, titlesForGraphs);
	vector<TGraph*> grWaveformsPadded = makePaddedGraphs(grWaveformsInt, 0, xLabel, yLabel, titlesForGraphs);
	xLabel = "Frequency (Hz)"; yLabel = "Power Spectral Density (mV/Hz)";
	vector<TGraph*> grWaveformsPowerSpectrum = makePowerSpectrumGraphs(grWaveformsPadded, xLabel, yLabel, titlesForGraphs);

	vector<TGraph*> spareElecChanGraphs;
	spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(6));
	spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(14));
	spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(22));
	spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(30));
	int hasBadSpareChanIssue = hasSpareChannelIssue(spareElecChanGraphs);
	if(hasSpareChannelIssue){
		cout<<"Has bad spare chan issue! Like, yes, actually"<<endl;
	}
	deleteGraphVector(spareElecChanGraphs);


	char cw_file_name[400];
	sprintf(cw_file_name,"%s/CWID/A%d/all_runs/CWID_station_%d_run_%d.root",DataDirPath,station,station,runNum);
	TFile *NewCWFile = TFile::Open(cw_file_name);
	if(!NewCWFile) {
		std::cerr << "Can't open new CW file\n";
		return -1;
	}
	TTree* NewCWTree = (TTree*) NewCWFile->Get("NewCWTree");   
	if(!NewCWTree) {
		std::cerr << "Can't find NewCWTree\n";
		return -1;
	}
	vector<vector<double> > *badFreqs_fwd =0;
	vector<vector<double> > *badFreqs_back=0;
	vector<vector<double> > *badSigmas_fwd=0;
	vector<vector<double> > *badSigmas_back=0;
	vector<vector<double> > *badFreqs_baseline=0;

	NewCWTree->SetBranchAddress("badFreqs_fwd",&badFreqs_fwd);
	NewCWTree->SetBranchAddress("badSigmas_fwd",&badSigmas_fwd);
	NewCWTree->SetBranchAddress("badFreqs_back",&badFreqs_back);
	NewCWTree->SetBranchAddress("badSigmas_back",&badSigmas_back);
	NewCWTree->SetBranchAddress("badFreqs_baseline",&badFreqs_baseline);

	//deal w/ CW cut
	//inputTree_CW->GetEntry(event);
	NewCWTree->GetEntry(event);

	bool isCutonCW_fwd[2]; isCutonCW_fwd[0]=false; isCutonCW_fwd[1]=false;
	bool isCutonCW_back[2]; isCutonCW_back[0]=false; isCutonCW_back[1]=false;
	bool isCutonCW_baseline[2]; isCutonCW_baseline[0]=false; isCutonCW_baseline[1]=false;
	
	for(int pol=0; pol<badFreqs_baseline->size(); pol++){
		vector<double> badFreqListLocal_baseline = badFreqs_baseline->at(pol);
		if(badFreqListLocal_baseline.size()>0) isCutonCW_baseline[pol]=true;
	}
	for(int pol=0; pol<2; pol++){
		char run_summary_filename[400];
		sprintf(run_summary_filename,"%s/RunSummary/A%d/all_runs/run_summary_station_%d_run_%d.root",DataDirPath,station,station,runNum);
		TFile *SummaryFile = TFile::Open(run_summary_filename);
		if(!SummaryFile) {
			std::cerr << "Can't open summary file\n";
			return -1;
		}
		TTree* SummaryTree = (TTree*) SummaryFile->Get("BaselineTree");   
		if(!SummaryTree) {
			std::cerr << "Can't find SummaryTree\n";
			return -1;
		}
		vector <TGraph*> average;
		average.resize(16);
		stringstream ss1;
		for(int i=0; i<16; i++){
			ss1.str(""); ss1<<"baselines_RF_chan_"<<i;
			SummaryTree->SetBranchAddress(ss1.str().c_str(),&average[i]);
		}
		SummaryTree->GetEntry(0);
		// char *plotPath(getenv("PLOT_PATH"));
		// if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;
		// TCanvas *c = new TCanvas("","",1100,850);
		// c->Divide(4,4);
		// for(int i=0; i<16; i++){
		// 	c->cd(i+1);
		// 	average[i]->Draw("ALP");
		// }
		// char save_temp_title[300];
		// sprintf(save_temp_title,"%s/trouble_events/Run%d_JustBaseline.png",plotPath,runNum);
		// c->SaveAs(save_temp_title);
		// delete c;
		vector<int> chan_exclusion_list;
		vector<double> baseline_CW_freqs = CWCut_TB(waveforms, average, pol, 6., 5.5, station, 3, chan_exclusion_list, runNum, event, false);
		cout<<"Recomputed baseline CW freqs is "<<baseline_CW_freqs.size()<<endl;
	}

	double threshCW = 1.5;
	vector<double> badFreqList_fwd;
	vector<double> badSigmaList_fwd;
	for(int pol=0; pol<badFreqs_fwd->size(); pol++){
		badFreqList_fwd=badFreqs_fwd->at(pol);
		badSigmaList_fwd=badSigmas_fwd->at(pol);
		for(int iCW=0; iCW<badFreqList_fwd.size(); iCW++){
			if(
				badSigmaList_fwd[iCW] > threshCW 
				&& 
				abs(300. - badFreqList_fwd[iCW]) > 2.
				&&
				abs(500. - badFreqList_fwd[iCW]) > 2.
				&&
				abs(125. - badFreqList_fwd[iCW]) > 2.
				&&
				badFreqList_fwd[iCW] < 850.
			){
				isCutonCW_fwd[pol] = true;
			}
		}
	}
	vector<double> badFreqList_back;
	vector<double> badSigmaList_back;
	for(int pol=0; pol<badFreqs_back->size(); pol++){
		badFreqList_back=badFreqs_back->at(pol);
		badSigmaList_back=badSigmas_back->at(pol);
		for(int iCW=0; iCW<badFreqList_back.size(); iCW++){
			if(
				badSigmaList_back[iCW] > threshCW 
				&& 
				abs(300. - badFreqList_back[iCW]) > 2.
				&&
				abs(500. - badFreqList_back[iCW]) > 2.
				&&
				abs(125. - badFreqList_back[iCW]) > 2.
				&&
				badFreqList_fwd[iCW] < 850.
			){
				isCutonCW_back[pol] = true;
			}
		}
	}

	bool skipCW=true;
	for(int pol=0; pol<2; pol++){
		if(skipCW || pol!=problempol) continue;
		if(isCutonCW_fwd[pol] || isCutonCW_back[pol] || isCutonCW_baseline[pol]){
			printf("	Has CW issue in pol %d \n", pol);
			printf("		CW in FWD %d, BWD %d, or baseline %d? \n", isCutonCW_fwd[pol], isCutonCW_back[pol], isCutonCW_baseline[pol]);
			//get the frequencies to notch
			vector<double> badFreqListLocal_fwd;
			vector <double> badFreqListLocal_back;
			vector <double> mergedFreqList;

			//merge the two lists of frequencies
			//if it's cut going both forward and backward
			if(isCutonCW_fwd[pol] && isCutonCW_back[pol]){
				badFreqListLocal_fwd=badFreqs_fwd->at(pol);
				badFreqListLocal_back=badFreqs_back->at(pol);
				for(int iFreq=0; iFreq<badFreqListLocal_fwd.size(); iFreq++){
					mergedFreqList.push_back(badFreqListLocal_fwd[iFreq]);
				}
				for(int iFreq=0; iFreq<badFreqListLocal_back.size(); iFreq++){
					double new_freq=badFreqListLocal_back[iFreq];
					for(int iFreqOld=0; iFreqOld<badFreqListLocal_fwd.size(); iFreqOld++){
						if(abs(new_freq-mergedFreqList[iFreqOld])>0.1){
							mergedFreqList.push_back(new_freq);
						}
					}
				}
			}
			//if it's cut only going forward
			else if(isCutonCW_fwd[pol] && !isCutonCW_back[pol]){
				badFreqListLocal_fwd=badFreqs_fwd->at(pol);
				for(int iFreq=0; iFreq<badFreqListLocal_fwd.size(); iFreq++){
					mergedFreqList.push_back(badFreqListLocal_fwd[iFreq]);
				}
			}
			//if it's cut only going backward
			else if(!isCutonCW_fwd[pol] && isCutonCW_back[pol]){
				badFreqListLocal_back=badFreqs_back->at(pol);
				for(int iFreq=0; iFreq<badFreqListLocal_back.size(); iFreq++){
					mergedFreqList.push_back(badFreqListLocal_back[iFreq]);
				}
			}

			vector<double> more_freqs_to_add;
			vector<double> badFreqListLocal_baseline = badFreqs_baseline->at(pol);
			if(mergedFreqList.size()>0){ //do we already have frequencies to check against?
				//loop over everything identified by the CW baseline cut
				for(int newFreq=0; newFreq<badFreqListLocal_baseline.size(); newFreq++){
					double new_freq = badFreqListLocal_baseline[newFreq];
					//now, loop over everything already in the list
					for(int oldFreq=0; oldFreq<mergedFreqList.size(); oldFreq++){
						//if there's a genuinely new frequency, add it to the list of things to be adde
						if(abs(new_freq-mergedFreqList[oldFreq])>0.1){
							more_freqs_to_add.push_back(new_freq);
						}
					}
				}
			}
			else{ //otherwise we take only those found by the CW ID cut
				for(int newFreq=0; newFreq<badFreqListLocal_baseline.size(); newFreq++){
					double new_freq = badFreqListLocal_baseline[newFreq];
					more_freqs_to_add.push_back(new_freq);
				}
			}

			//now actually add it to the merged freq list
			for(int iFreq=0; iFreq<more_freqs_to_add.size(); iFreq++){
				mergedFreqList.push_back(more_freqs_to_add[iFreq]);
			}

			//they need to be in smaller -> larger order for notching
			sort(mergedFreqList.begin(), mergedFreqList.end());
			vector <double> uniqueNotchFreqs;
			vector <double> uniqueNotchBands;			
			theCorrelators[0]->pickFreqsAndBands(mergedFreqList,uniqueNotchFreqs,uniqueNotchBands);
			for (int i = 0; i < uniqueNotchFreqs.size(); ++i)
			{
				printf("			Unique freq to be notched is %.2f with width %.2f \n", uniqueNotchFreqs[i],uniqueNotchBands[i]);
			}
			vector <TGraph*> grWaveformsRaw = makeGraphsFromRF(realAtriEvPtr,16,xLabel,yLabel,titlesForGraphs);
			vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(grWaveformsRaw, 0.5, xLabel, yLabel, titlesForGraphs);
			vector<TGraph*> grWaveformsPadded = makePaddedGraphs(grWaveformsInt, 0, xLabel, yLabel, titlesForGraphs);
			vector <TGraph*> grNotched;
			for(int i=0; i<16; i++){
				TGraph *grNotchAmp = theCorrelators[0]->applyAdaptiveFilter_singleAnt_FiltMany(grWaveformsPadded[i],uniqueNotchFreqs,uniqueNotchBands);
				grNotched.push_back(theCorrelators[0]->GeometricFilter(grNotchAmp,uniqueNotchFreqs,uniqueNotchBands,uniqueNotchFreqs));
				delete grNotchAmp;
			}
			vector<TGraph*> grWaveformsPowerSpectrum = makePowerSpectrumGraphs(grWaveformsPadded, xLabel, yLabel, titlesForGraphs);
			vector<TGraph*> grWaveformsPowerSpectrum_notched = makePowerSpectrumGraphs(grNotched, xLabel, yLabel, titlesForGraphs);
			
			char save_temp_title[300];
			sprintf(save_temp_title,"%s/cal_reco_events/%d.%d.%d_Run%d_Ev%d_ProblemPol%d_WaveformsNotch.png",plotPath,year_now,month_now,day_now,runNum,event,problempol);
			TCanvas *cWave = new TCanvas("","",4*1100,4*850);
			cWave->Divide(4,4);
			for(int i=0; i<16; i++){
				cWave->cd(i+1);
				grWaveformsRaw[i]->Draw("AL");
				grWaveformsRaw[i]->SetLineWidth(3);
				grNotched[i]->Draw("same");
				grNotched[i]->SetLineWidth(3);
				grNotched[i]->SetLineColor(kRed);
			}
			cWave->SaveAs(save_temp_title);
			delete cWave;

			sprintf(save_temp_title,"%s/cal_reco_events/%d.%d.%d_Run%d_Ev%d_ProblemPol%d_SpectraNotch.png",plotPath,year_now,month_now,day_now,runNum,event,problempol);
			TCanvas *cSpec = new TCanvas("","",4*1100,4*850);
			cSpec->Divide(4,4);
			for(int i=0; i<16; i++){
				cSpec->cd(i+1);
				grWaveformsPowerSpectrum[i]->Draw("AL");
				grWaveformsPowerSpectrum[i]->SetLineWidth(3);
				grWaveformsPowerSpectrum_notched[i]->Draw("same");
				grWaveformsPowerSpectrum_notched[i]->SetLineWidth(3);
				grWaveformsPowerSpectrum_notched[i]->SetLineColor(kRed);
				gPad->SetLogy();
			}
			cSpec->SaveAs(save_temp_title);
			delete cSpec;

			TH2D *map_30m_V;
			TH2D *map_300m_V;
			TH2D *map_30m_H;
			TH2D *map_300m_H;
			vector <int> chan_list_V;
			vector <int> chan_list_H;
			for(int chan=0; chan<=7; chan++){
				chan_list_V.push_back(chan);
				chan_list_H.push_back(chan+8);
			}

			vector<double> chan_SNRs;

			if(station==2){
				//for station 2, we need to exclude channel 15 from the analysis
				chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
			}
			else if(station==3){
				//for station 3 years 2014 and 2015, we need to drop string 4 (channels 3, 7, 11, 15) altogether
				if(runNum>getA3BadRunBoundary()){
					chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 3), chan_list_V.end());
					chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 7), chan_list_V.end());

					chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 11), chan_list_H.end());
					chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
				}
			}

			map_30m_V = theCorrelators[0]->getInterferometricMap_RT_FiltMany_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Vpol, false,chan_list_V,chan_SNRs,0,-1,uniqueNotchFreqs,uniqueNotchBands) ;
			map_300m_V = theCorrelators[1]->getInterferometricMap_RT_FiltMany_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Vpol, false,chan_list_V,chan_SNRs,0,-1,uniqueNotchFreqs,uniqueNotchBands);
			map_30m_H = theCorrelators[0]->getInterferometricMap_RT_FiltMany_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Hpol, false,chan_list_H,chan_SNRs,0,-1,uniqueNotchFreqs,uniqueNotchBands);
			map_300m_H = theCorrelators[1]->getInterferometricMap_RT_FiltMany_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Hpol, false,chan_list_H,chan_SNRs,0,-1,uniqueNotchFreqs,uniqueNotchBands);

			int PeakTheta_Recompute_30m_H;
			int PeakTheta_Recompute_300m_H;
			int PeakPhi_Recompute_30m_H;
			int PeakPhi_Recompute_300m_H;
			double PeakCorr_Recompute_30m_H;
			double PeakCorr_Recompute_300m_H;
			int PeakTheta_Recompute_30m_V;
			int PeakTheta_Recompute_300m_V;
			int PeakPhi_Recompute_30m_V;
			int PeakPhi_Recompute_300m_V;
			double PeakCorr_Recompute_30m_V;
			double PeakCorr_Recompute_300m_V;
			getCorrMapPeak(map_30m_H,PeakTheta_Recompute_30m_H,PeakPhi_Recompute_30m_H,PeakCorr_Recompute_30m_H);
			getCorrMapPeak(map_300m_H,PeakTheta_Recompute_300m_H,PeakPhi_Recompute_300m_H,PeakCorr_Recompute_300m_H);
			getCorrMapPeak(map_30m_V,PeakTheta_Recompute_30m_V,PeakPhi_Recompute_30m_V,PeakCorr_Recompute_30m_V);
			getCorrMapPeak(map_300m_V,PeakTheta_Recompute_300m_V,PeakPhi_Recompute_300m_V,PeakCorr_Recompute_300m_V);

			printf("	Rconstruction Information\n");
			printf("		30m H theta and phi %d and %d \n", PeakTheta_Recompute_30m_H, PeakPhi_Recompute_30m_H);
			stringstream ss30H;
			ss30H<<" Peak Theta, Phi is "<<PeakTheta_Recompute_30m_H<<" , "<<PeakPhi_Recompute_30m_H;
			map_30m_H->SetTitle(ss30H.str().c_str());
			printf("		300m H theta and phi %d and %d \n", PeakTheta_Recompute_300m_H, PeakPhi_Recompute_300m_H);
			stringstream ss300H;
			ss300H<<" Peak Theta, Phi is "<<PeakTheta_Recompute_300m_H<<" , "<<PeakPhi_Recompute_300m_H;
			map_300m_H->SetTitle(ss300H.str().c_str());
			printf("		30m V theta and phi %d and %d \n", PeakTheta_Recompute_30m_V, PeakPhi_Recompute_30m_V);
			stringstream ss30V;
			ss30V<<" Peak Theta, Phi is "<<PeakTheta_Recompute_30m_V<<" , "<<PeakPhi_Recompute_30m_V;
			map_30m_V->SetTitle(ss30V.str().c_str());
			printf("		300m V theta and phi %d and %d \n", PeakTheta_Recompute_300m_V, PeakPhi_Recompute_300m_V);
			stringstream ss300V;
			ss300V<<" Peak Theta, Phi is "<<PeakTheta_Recompute_300m_V<<" , "<<PeakPhi_Recompute_300m_V;
			map_300m_V->SetTitle(ss300V.str().c_str());

			TCanvas *cMaps = new TCanvas("","",2*1100,2*850);
			cMaps->Divide(2,2);
				cMaps->cd(3);
				map_30m_V->Draw("colz");
				cMaps->cd(4);
				map_30m_H->Draw("colz");
				cMaps->cd(1);
				map_300m_V->Draw("colz");
				cMaps->cd(2);
				map_300m_H->Draw("colz");
			sprintf(save_temp_title,"%s/cal_reco_events/%d.%d.%d_Run%d_Ev%d_ProblemPol%d_FilteredMaps.png",plotPath,year_now,month_now,day_now,runNum,event,problempol);
			cMaps->SaveAs(save_temp_title);
			delete cMaps;
			delete map_30m_V; delete map_300m_V; delete map_30m_H; delete map_300m_H;

			deleteGraphVector(grWaveformsRaw);
			deleteGraphVector(grWaveformsInt);
			deleteGraphVector(grWaveformsPadded);
			deleteGraphVector(grNotched);
			deleteGraphVector(grWaveformsPowerSpectrum);
			deleteGraphVector(grWaveformsPowerSpectrum_notched);
		}
	}


	vector <int> chan_list_V;
	vector <int> chan_list_H;
	for(int chan=0; chan<=7; chan++){
		chan_list_V.push_back(chan);
		chan_list_H.push_back(chan+8);
	}

	bool do_reco=false;
	if(do_reco){
		TH2D *map_30m_V;
		TH2D *map_300m_V;
		TH2D *map_30m_H;
		TH2D *map_300m_H;
		// TH2D *map_30m_V_select;

		map_30m_V = theCorrelators[0]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Vpol, 0, chan_list_V);
		map_300m_V = theCorrelators[1]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Vpol, 0, chan_list_V);
		map_30m_H = theCorrelators[0]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Hpol, 0, chan_list_H);
		map_300m_H = theCorrelators[1]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Hpol, 0, chan_list_H);

		int PeakTheta_Recompute_30m;
		int PeakTheta_Recompute_300m;
		int PeakPhi_Recompute_30m;
		int PeakPhi_Recompute_300m;
		double PeakCorr_Recompute_30m;
		double PeakCorr_Recompute_300m;
		double MinCorr_Recompute_30m;
		double MinCorr_Recompute_300m;
		double MeanCorr_Recompute_30m;
		double MeanCorr_Recompute_300m;
		double RMSCorr_Recompute_30m;
		double RMSCorr_Recompute_300m;
		double PeakSigma_Recompute_30m;
		double PeakSigma_Recompute_300m;
		getCorrMapPeak_wStats(map_30m_V,PeakTheta_Recompute_30m,PeakPhi_Recompute_30m,PeakCorr_Recompute_30m,MinCorr_Recompute_30m,MeanCorr_Recompute_30m,RMSCorr_Recompute_30m,PeakSigma_Recompute_30m);
		getCorrMapPeak_wStats(map_300m_V,PeakTheta_Recompute_300m,PeakPhi_Recompute_300m,PeakCorr_Recompute_300m,MinCorr_Recompute_300m,MeanCorr_Recompute_300m,RMSCorr_Recompute_300m,PeakSigma_Recompute_300m);

		printf("30m theta and phi %d and %d \n", PeakTheta_Recompute_30m, PeakPhi_Recompute_30m);
		printf("300m theta and phi %d and %d \n", PeakTheta_Recompute_300m, PeakPhi_Recompute_300m);

		// vector <int> chan_list;
		// chan_list.push_back(4);
		// chan_list.push_back(5);
		// chan_list.push_back(6);
		// chan_list.push_back(7);
		// map_30m_V_select = theCorrelators[0]->getInterferometricMap_RT_select(settings,detector,realAtriEvPtr,Vpol,0,chan_list,0);


		TCanvas *cMaps = new TCanvas("","",2*1100,3*850);
		cMaps->Divide(2,2);
			cMaps->cd(3);
			map_30m_V->Draw("colz");
			cMaps->cd(4);
			map_30m_H->Draw("colz");
			cMaps->cd(1);
			map_300m_V->Draw("colz");
			cMaps->cd(2);
			map_300m_H->Draw("colz");
			// cMaps->cd(5);
			// map_30m_V_select->Draw("colz");
		char save_temp_title[400];		
		sprintf(save_temp_title,"%s/cal_reco_events/%d.%d.%d_Run%d_Ev%d_Maps.png",plotPath,year_now,month_now,day_now,runNum,event);
		cMaps->SaveAs(save_temp_title);
		delete cMaps;
		delete map_30m_V; delete map_300m_V; delete map_30m_H; delete map_300m_H; 
		// delete map_30m_V_select;
	}

	bool do_reco_snrweighted_newnormalization_select=true;
	if(do_reco_snrweighted_newnormalization_select){

		char filter_file_name[400];
		sprintf(filter_file_name,"%s/ProcessedFile/A%d/all_runs/processed_station_%d_run_%d_filter.root",DataDirPath,station,station,runNum);
		TFile *filterFile = TFile::Open(filter_file_name);
		TTree *filterTree = (TTree*) filterFile->Get("OutputTree");
		filterTree->SetBranchAddress("VPeakOverRMS", &VPeakOverRMS);
		filterTree->GetEvent(event);

		vector <int> chan_list_V;
		vector <int> chan_list_H;
		vector<double> chan_SNRs;
		for(int chan=0; chan<=7; chan++){
			chan_list_V.push_back(chan);
			chan_list_H.push_back(chan+8);
		}
		for(int i=0; i<16; i++){
			// chan_SNRs.push_back(1.);
			chan_SNRs.push_back(VPeakOverRMS[i]);
			// printf("SNR is %.2f \n", VPeakOverRMS[i]);
		}
		filterFile->Close();

		if(station==2){
			//for station 2, we need to exclude channel 15 from the analysis
			chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
		}
		else if(station==3){
			//for station 3 years 2014 and 2015, we need to drop string 4 (channels 3, 7, 11, 15) altogether
			if(runNum>getA3BadRunBoundary()){
				chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 3), chan_list_V.end());
				chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 7), chan_list_V.end());
				chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 11), chan_list_H.end());
				chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
			}
		}

		TH2D *map_41m_V;
		TH2D *map_300m_V;
		TH2D *map_41m_H;
		TH2D *map_300m_H;

		bool AraSim=false;
		int solNum=0;
		map_41m_V = theCorrelators[0]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Vpol, AraSim, chan_list_V, chan_SNRs, solNum);
		map_300m_V = theCorrelators[1]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Vpol, AraSim, chan_list_V, chan_SNRs,  solNum);
		map_41m_H = theCorrelators[0]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Hpol, AraSim, chan_list_H, chan_SNRs, solNum);
		map_300m_H = theCorrelators[1]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Hpol, AraSim, chan_list_H, chan_SNRs, solNum);

		// map_41m_V = theCorrelators[0]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Vpol, AraSim, chan_list_V);
		// map_300m_V = theCorrelators[1]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Vpol, AraSim, chan_list_V);
		// map_41m_H = theCorrelators[0]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Hpol, AraSim, chan_list_H);
		// map_300m_H = theCorrelators[1]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Hpol, AraSim, chan_list_H);

		int PeakTheta_Recompute_41m_H;
		int PeakTheta_Recompute_300m_H;
		int PeakPhi_Recompute_41m_H;
		int PeakPhi_Recompute_300m_H;
		double PeakCorr_Recompute_41m_H;
		double PeakCorr_Recompute_300m_H;
		int PeakTheta_Recompute_41m_V;
		int PeakTheta_Recompute_300m_V;
		int PeakPhi_Recompute_41m_V;
		int PeakPhi_Recompute_300m_V;
		double PeakCorr_Recompute_41m_V;
		double PeakCorr_Recompute_300m_V;
		getCorrMapPeak(map_41m_H,PeakTheta_Recompute_41m_H,PeakPhi_Recompute_41m_H,PeakCorr_Recompute_41m_H);
		getCorrMapPeak(map_300m_H,PeakTheta_Recompute_300m_H,PeakPhi_Recompute_300m_H,PeakCorr_Recompute_300m_H);
		getCorrMapPeak(map_41m_V,PeakTheta_Recompute_41m_V,PeakPhi_Recompute_41m_V,PeakCorr_Recompute_41m_V);
		getCorrMapPeak(map_300m_V,PeakTheta_Recompute_300m_V,PeakPhi_Recompute_300m_V,PeakCorr_Recompute_300m_V);

		printf("	Rconstruction Information\n");
		printf("		41m H theta and phi %d and %d \n", PeakTheta_Recompute_41m_H, PeakPhi_Recompute_41m_H);
		stringstream ss30H;
		ss30H<<" 41m H Peak Theta, Phi is "<<PeakTheta_Recompute_41m_H<<" , "<<PeakPhi_Recompute_41m_H;
		map_41m_H->SetTitle(ss30H.str().c_str());
		printf("		300m H theta and phi %d and %d \n", PeakTheta_Recompute_300m_H, PeakPhi_Recompute_300m_H);
		stringstream ss300H;
		ss300H<<" 300m H Peak Theta, Phi is "<<PeakTheta_Recompute_300m_H<<" , "<<PeakPhi_Recompute_300m_H;
		map_300m_H->SetTitle(ss300H.str().c_str());
		printf("		41m V theta and phi %d and %d \n", PeakTheta_Recompute_41m_V, PeakPhi_Recompute_41m_V);
		stringstream ss30V;
		ss30V<<" 41m V Peak Theta, Phi is "<<PeakTheta_Recompute_41m_V<<" , "<<PeakPhi_Recompute_41m_V;
		map_41m_V->SetTitle(ss30V.str().c_str());
		printf("		300m V theta and phi %d and %d \n", PeakTheta_Recompute_300m_V, PeakPhi_Recompute_300m_V);
		stringstream ss300V;
		ss300V<<" 300m V Peak Theta, Phi is "<<PeakTheta_Recompute_300m_V<<" , "<<PeakPhi_Recompute_300m_V;
		map_300m_V->SetTitle(ss300V.str().c_str());

		// beautify_TH2D();
		TCanvas *cMaps = new TCanvas("","",2*1100,2*850);
		cMaps->Divide(2,2);
			cMaps->cd(3);
			map_41m_V->Draw("colz");
			// map_41m_V->GetXaxis()->SetRangeUser(50,80);
			// map_41m_V->GetYaxis()->SetRangeUser(-10,20);
			// map_41m_V->GetZaxis()->SetRangeUser(0.005,0.04);
			cMaps->cd(4);
			map_41m_H->Draw("colz");
			cMaps->cd(1);
			map_300m_V->Draw("colz");
			cMaps->cd(2);
			map_300m_H->Draw("colz");
		char save_temp_title[400];		
		sprintf(save_temp_title,"%s/single_events/%d.%d.%d_Run%d_Ev%d_Maps_SNRweighted.png",plotPath,year_now,month_now,day_now,runNum,event);
		cMaps->SaveAs(save_temp_title);
		delete cMaps;
		delete map_41m_V; delete map_300m_V; delete map_41m_H; delete map_300m_H; 
		// delete map_41m_V_select;
	}



	char save_temp_title[300];
	sprintf(save_temp_title,"%s/cal_reco_events/%d.%d.%d_Run%d_Ev%d_Waveforms.png",plotPath,year_now,month_now,day_now,runNum,event);
	TCanvas *cWave = new TCanvas("","",4*1100,4*850);
	cWave->Divide(4,4);
	for(int i=0; i<16; i++){
		cWave->cd(i+1);
		waveforms[i]->Draw("AL");
		waveforms[i]->SetLineWidth(3);
	}
	cWave->SaveAs(save_temp_title);
	delete cWave;

	sprintf(save_temp_title,"%s/cal_reco_events/%d.%d.%d_Run%d_Ev%d_Spectra.png",plotPath,year_now,month_now,day_now,runNum,event);
	TCanvas *cSpec = new TCanvas("","",4*1100,4*850);
	cSpec->Divide(4,4);
	for(int i=0; i<16; i++){
		cSpec->cd(i+1);
		grWaveformsPowerSpectrum[i]->Draw("AL");
		grWaveformsPowerSpectrum[i]->SetLineWidth(3);
		gPad->SetLogy();
	}
	cSpec->SaveAs(save_temp_title);
	delete cSpec;
	for(int i=0; i<16; i++){
		delete waveforms[i];
		delete grWaveformsInt[i];
		delete grWaveformsPadded[i];
		delete grWaveformsPowerSpectrum[i];
	}
	delete realAtriEvPtr;
	mapFile->Close();
	delete mapFile;
	return 0;
}
