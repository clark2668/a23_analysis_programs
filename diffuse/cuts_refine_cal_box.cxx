////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	cuts_refine_cal_box
////	Reads recco_save files and refines the cal pule cut box
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
AraAntPol::AraAntPol_t Vpol = AraAntPol::kVertical;
AraAntPol::AraAntPol_t Hpol = AraAntPol::kHorizontal;

//custom analysis includes
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_Cuts.h"

bool doRezero=0;

int PlotThisEvent(int station, int config, int runNum, int event, Settings *settings, Detector *detector, RayTraceCorrelator *theCorrelators[2]);

char *DataDirPath(getenv("DATA_DIR"));
char *PedDirPath(getenv("PED_DIR"));
char *plotPath(getenv("PLOT_PATH"));

using namespace std;

int main(int argc, char **argv)
{
	gStyle->SetOptStat(1111);
	
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
		if(runNum==2428 || runNum==3442) continue;

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
		trees->SetBranchAddress("phi_41_V",&phi_41_V);
		trees->SetBranchAddress("theta_41_V",&theta_41_V);
		trees->SetBranchAddress("phi_41_H",&phi_41_H);
		trees->SetBranchAddress("theta_41_H",&theta_41_H);
		trees->SetBranchAddress("corr_val_V",&corr_val_V);
		trees->SetBranchAddress("corr_val_H",&corr_val_H);
		trees->SetBranchAddress("cal",&isCal);
		trees->SetBranchAddress("isBad",&isBad);


		int numEntries = trees->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;

		//now to loop over events
		for(int event=0; event<numEntries; event++){

			trees->GetEvent(event);

			if(!isCal || isBad) continue;

			num_total++;
			if(corr_val_V>corr_val_H){
				locations[0]->Fill(phi_41_V,theta_41_V);
				locations_zoom[0]->Fill(phi_41_V,theta_41_V);
				locations_phi_slice[0]->Fill(phi_41_V);
				locations_theta_slice[0]->Fill(theta_41_V);

				// double this_phi = phi_41_V;
				// double this_theta = theta_41_V;
				// if(this_phi>-50 && this_phi<0 && this_theta<-30){
				// 	PlotThisEvent(station, config, runNum, event, settings, detector, theCorrelators);
				// }

			}
			else if(corr_val_H>corr_val_V){
				locations[1]->Fill(phi_41_H,theta_41_H);
				locations_zoom[1]->Fill(phi_41_H,theta_41_H);
				locations_phi_slice[1]->Fill(phi_41_H);
				locations_theta_slice[1]->Fill(theta_41_H);
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

}

int PlotThisEvent(int station, int config, int runNum, int event, Settings *settings, Detector *detector, RayTraceCorrelator *theCorrelators[2]){
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

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

	vector <int> chan_list_V;
	vector <int> chan_list_H;
	for(int chan=0; chan<=7; chan++){
		chan_list_V.push_back(chan);
		chan_list_H.push_back(chan+8);
	}

	bool do_reco=true;
	if(do_reco){
		TH2D *map_30m_V;
		TH2D *map_300m_V;
		TH2D *map_30m_H;
		TH2D *map_300m_H;
		TH2D *map_30m_V_select;

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

		vector <int> chan_list;
		chan_list.push_back(4);
		chan_list.push_back(5);
		chan_list.push_back(6);
		chan_list.push_back(7);
		map_30m_V_select = theCorrelators[0]->getInterferometricMap_RT_select(settings,detector,realAtriEvPtr,Vpol,0,chan_list,0);


		TCanvas *cMaps = new TCanvas("","",2*1100,3*850);
		cMaps->Divide(2,3);
			cMaps->cd(3);
			map_30m_V->Draw("colz");
			cMaps->cd(4);
			map_30m_H->Draw("colz");
			cMaps->cd(1);
			map_300m_V->Draw("colz");
			cMaps->cd(2);
			map_300m_H->Draw("colz");
			cMaps->cd(5);
			map_30m_V_select->Draw("colz");
		char save_temp_title[400];		
		sprintf(save_temp_title,"%s/cal_reco_events/%d.%d.%d_Run%d_Ev%d_Maps.png",plotPath,year_now,month_now,day_now,runNum,event);
		cMaps->SaveAs(save_temp_title);
		delete cMaps;
		delete map_30m_V; delete map_300m_V; delete map_30m_H; delete map_300m_H; 
		delete map_30m_V_select;
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