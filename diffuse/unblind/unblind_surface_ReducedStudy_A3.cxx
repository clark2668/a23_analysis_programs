////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	unblind_surface_SaveToTree.cxx
////	unblind the surface cut
////	if an event passes surface cut, save it to a reduced root file
////
////	Sep 2019
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
#include "TF1.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TChain.h"
#include "TTimeStamp.h"
#include "TLegend.h"

//AraRoot includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraQualCuts.h"

// analysis custom
#include "tools_Cuts.h"
#include "tools_Stats.h"
#include "tools_CommandLine.h"
#include "tools_outputObjects.h"
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"

// for ray trace correlator
#include "Settings.h"
#include "Event.h"
#include "Detector.h"
#include "Report.h"
#include "RayTraceCorrelator.h"
#include "RaySolver.h"
#include "Position.h"
#include "Vector.h"
#include "Settings.h"

using namespace std;

AraAntPol::AraAntPol_t Vpol = AraAntPol::kVertical;
AraAntPol::AraAntPol_t Hpol = AraAntPol::kHorizontal;

// this function is simple enough, just declare it up front
// otherwise it's at the end and clogs things up
void configure(TH1D *gr){
	gr->GetXaxis()->SetLabelSize(0.05);
	gr->GetXaxis()->SetTitleSize(0.05);
	gr->GetYaxis()->SetLabelSize(0.05);
	gr->GetYaxis()->SetTitleSize(0.05);
	// gr->GetYaxis()->SetTitleOffset(1.1);
	gStyle->SetTitleFontSize(0.05);
	gr->SetLineWidth(2);
}

int PlotThisEvent(int station, int runNum, int event, int problempol);

double getAntennaArrivalTheta(int station, Settings *settings, Position target, RaySolver *raysolver, double stationCenter[3], int thetaIn, int phiIn){
	double phiWave = double(phiIn)*TMath::DegToRad();
	double thetaWave = double(thetaIn)*TMath::DegToRad();
	double R = 300;
	
	Double_t xs = R*TMath::Cos(thetaWave)*TMath::Cos(phiWave);
	Double_t ys = R*TMath::Cos(thetaWave)*TMath::Sin(phiWave);
	Double_t zs = R*TMath::Sin(thetaWave);

	xs = xs + stationCenter[0];
	ys = ys + stationCenter[1];
	zs = zs + stationCenter[2];
	Position source;
	source.SetXYZ(xs, ys, zs);

	vector < vector <double> > outputs; //place for the answer
	raysolver->Solve_Ray_org(source, target, outputs, settings);
	double thetaOut=-1000;
	if(outputs.size()>0){
		thetaOut=outputs[2][0]*TMath::RadToDeg() -90.;
	}
	thetaOut = abs(thetaOut-90.);
	return thetaOut;
}

int main(int argc, char **argv)
{
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;


	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <1-station> <2-config> <3-save_vals_file_1> <4-save_vals_file_2> .... <n-save_vals_file_n>"<<endl;;
		return -1;
	}
	int station = atoi(argv[1]);
	int config = atoi(argv[2]);

	if(station!=2 && station!=3){
		printf("No good! You asked for station %d, but this code only works for stations 2 and 3 \n",station);
		return -1;
	}

	/*
		We're gonna need some ray tracing tools to get arrival angle information
	*/

	Settings *settings = new Settings();
	string setupfile = "setup.txt";
	settings->ReadFile(setupfile);
	settings->NOFZ=1; //yes, variable depth index of refraction
	AraGeomTool *araGeom = AraGeomTool::Instance(); // need a geom tool
	double stationCenter[3]={0.};
	for(int i=0; i<16; i++){
		for(int ii=0; ii<3; ii++){
			stationCenter[ii]+=(araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[ii]);
		}
	}
	for(int ii=0; ii<3; ii++){
		stationCenter[ii]/=16.;
	}
	Position target;
	double x1 = araGeom->getStationInfo(station)->getAntennaInfo(0)->antLocation[0];
	double y1 = araGeom->getStationInfo(station)->getAntennaInfo(0)->antLocation[1];
	double z1 = araGeom->getStationInfo(station)->getAntennaInfo(0)->antLocation[2];
	target.SetXYZ(x1,y1,z1);
	RaySolver *raysolver = new RaySolver; // invoke the ray solver


	/*
		Plot declaration
	*/
	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

	TTimeStamp start;
	TTimeStamp stop;
	int numBins;

	start.Set(2013, 01, 01, 00, 00, 0, 0, true, 0);
	stop.Set(2016, 12, 31, 24, 00, 0, 0, true, 0);
	numBins=4*365;

	int start_bin = start.GetSec();
	int stop_bin = stop.GetSec();

	TH1D *h1_events_vs_time[3];
	TH1D *h1_events_vs_time_finer[3];
	TH1D *h1_events_vs_run[3];

	stringstream ss;
	for(int i=0; i<3; i++){
		if(i==0){
			ss.str("");
			ss<<"VPol Only Events vs Time";
		}
		else if(i==1){
			ss.str("");
			ss<<"HPol Only Events vs Time";
		}
		else if(i==2){
			ss.str("");
			ss<<"Combo Events vs Time";			
		}
		h1_events_vs_time[i] = new TH1D(ss.str().c_str(),ss.str().c_str(),numBins, start_bin, stop_bin);
		h1_events_vs_time[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		h1_events_vs_time[i]->GetXaxis()->SetTimeOffset(0.,"GMT");

		if(i==0){
			ss.str("");
			ss<<"VPol Only Events vs Time Zoom";
		}
		else if(i==1){
			ss.str("");
			ss<<"HPol Only Events vs Time Zoom";
		}
		else if(i==2){
			ss.str("");
			ss<<"Combo Events vs Time Zoom";			
		}

		h1_events_vs_time_finer[i] = new TH1D(ss.str().c_str(),ss.str().c_str(),numBins*24, start_bin, stop_bin);
		h1_events_vs_time_finer[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		h1_events_vs_time_finer[i]->GetXaxis()->SetTimeOffset(0.,"GMT");
		
		if(i==0){
			ss.str("");
			ss<<"VPol Only Events vs Run";
		}
		else if(i==1){
			ss.str("");
			ss<<"HPol Only Events vs Run";
		}
		else if(i==2){
			ss.str("");
			ss<<"Combo Events vs Run";
		}
		h1_events_vs_run[i] = new TH1D(ss.str().c_str(),ss.str().c_str(),9000, 0, 9000);
	}

	TH1D *h1_PhiDist[4]; // 0 is vpol singlet, 1 is hpol singlet, 2 is vhcombo v, 3 is vhcombo h
	h1_PhiDist[0] = new TH1D("'Isolated' VPol Only Distribution Phi","'Isolated' VPol Only Distribution Phi",360, -180, 180);
	h1_PhiDist[1] = new TH1D("Isolated HPol Only Distribution Phi","Isolated HPol Only Distribution Phi",360, -180, 180);
	h1_PhiDist[2] = new TH1D("Isolated Combo Distribution, V Phi","Isolated Combo Distribution, V Phi",360, -180, 180);
	h1_PhiDist[3] = new TH1D("Isolated Combo Distribution, H Phi","Isolated Combo Distribution, H Phi",360, -180, 180);

	TH1D *h1_ThetaDist[4]; // 0 is vpol singlet, 1 is hpol singlet, 2 is vhcombo v, 3 is vhcombo h
	h1_ThetaDist[0] = new TH1D("'Isolated' VPol Only Distribution Theta","'Isolated' VPol Only Distribution Theta",180, -90, 90);
	h1_ThetaDist[1] = new TH1D("Isolated HPol Only Distribution Theta","Isolated HPol Only Distribution Theta",180, -90, 90);
	h1_ThetaDist[2] = new TH1D("Isolated Combo Distribution, V Theta","Isolated Combo Distribution, V Theta",180, -90, 90);
	h1_ThetaDist[3] = new TH1D("Isolated Combo Distribution, H Theta","Isolated Combo Distribution, H Theta",180, -90, 90);

	vector<int> v_hpol_isol_unixtime;
	vector<int> v_hpol_isol_runs;
	vector<int> v_hpol_isol_entries;

	vector<int> v_combo_isol_unixtime;
	vector<int> v_combo_isol_runs;
	vector<int> v_combo_isol_entries;

	vector<int> BadRunList=BuildBadRunList(station);
	vector<int> BadSurfaceRunList=BuildSurfaceRunList(station);

	// for distribution of number of events
	TH1D *distro_of_events = new TH1D("","",25,0,50);



	for(int file_num=3; file_num<argc; file_num++){
		
		TFile *inputFile = TFile::Open(argv[file_num]);
		if(!inputFile){
			cout<<"Can't open vals for cuts file!"<<endl;			
			return -1;
		}
		TTree *inputVTree = (TTree*) inputFile->Get("VTree");
		TTree *inputHTree = (TTree*) inputFile->Get("HTree");
		TTree *inputAllTree = (TTree*) inputFile->Get("AllTree");
		TTree *inputRecoTree = (TTree*) inputFile->Get("OutputTreeReco");
		TTree *inputFilterTree = (TTree*) inputFile->Get("OutputTree");
		TTree *inputSurfaceTree = (TTree*) inputFile->Get("SurfaceTree");

		double corr_val[2];
		double snr_val[2];
		int WFRMS[2];
		int theta_300[2];
		int phi_300[2];
		int theta_41[2];
		int phi_41[2];

		int Refilt[2];
		inputVTree->SetBranchAddress("Refilt_V",&Refilt[0]);
		inputHTree->SetBranchAddress("Refilt_H",&Refilt[1]);

		inputVTree->SetBranchAddress("corr_val_V_new",&corr_val[0]);
		inputVTree->SetBranchAddress("snr_val_V_new",&snr_val[0]);
		inputVTree->SetBranchAddress("wfrms_val_V_new",&WFRMS[0]);
		inputVTree->SetBranchAddress("theta_300_V_new",&theta_300[0]);
		inputVTree->SetBranchAddress("theta_41_V_new",&theta_41[0]);
		inputVTree->SetBranchAddress("phi_300_V_new",&phi_300[0]);
		inputVTree->SetBranchAddress("phi_41_V_new",&phi_41[0]);

		inputHTree->SetBranchAddress("corr_val_H_new",&corr_val[1]);
		inputHTree->SetBranchAddress("snr_val_H_new",&snr_val[1]);
		inputHTree->SetBranchAddress("wfrms_val_H_new",&WFRMS[1]);
		inputHTree->SetBranchAddress("theta_300_H_new",&theta_300[1]);
		inputHTree->SetBranchAddress("theta_41_H_new",&theta_41[1]);
		inputHTree->SetBranchAddress("phi_300_H_new",&phi_300[1]);
		inputHTree->SetBranchAddress("phi_41_H_new",&phi_41[1]);

		int isCal;
		int isSoft;
		int isShort;
		int isCW;
		int isNewBox;

		inputAllTree->SetBranchAddress("cal",&isCal);
		inputAllTree->SetBranchAddress("soft",&isSoft);
		inputAllTree->SetBranchAddress("short",&isShort);
		inputAllTree->SetBranchAddress("CW",&isCW);
		inputAllTree->SetBranchAddress("box",&isNewBox);

		int isSurf[2]; // a surface event after filtering?
		int isSurfEvent_top[2]; // a top event?

		inputAllTree->SetBranchAddress("surf_V_new2",&isSurf[0]);
		inputAllTree->SetBranchAddress("surf_H_new2",&isSurf[1]);

		inputAllTree->SetBranchAddress("surf_top_V",&isSurfEvent_top[0]);
		inputAllTree->SetBranchAddress("surf_top_H",&isSurfEvent_top[1]);

		int isBadEvent;
		double weight;
		int unixTime;
		int isFirstFiveEvent;
		int hasBadSpareChanIssue;
		int hasBadSpareChanIssue2;
		int runNum;
		int eventNumber;

		inputAllTree->SetBranchAddress("bad",&isBadEvent);
		inputAllTree->SetBranchAddress("weight",&weight);
		inputAllTree->SetBranchAddress("unixTime",&unixTime);
		inputAllTree->SetBranchAddress("isFirstFiveEvent",&isFirstFiveEvent);
		inputAllTree->SetBranchAddress("hasBadSpareChanIssue",&hasBadSpareChanIssue);
		inputAllTree->SetBranchAddress("hasBadSpareChanIssue2",&hasBadSpareChanIssue2);
		inputAllTree->SetBranchAddress("runNum",&runNum);
		inputAllTree->SetBranchAddress("eventNumber",&eventNumber);

		double frac_of_power_notched_V[8];
		double frac_of_power_notched_H[8];

		stringstream ss;
		for(int i=0; i<8; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			inputVTree->SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_V[i]);
		}
		for(int i=8; i<16; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			inputHTree->SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_H[i-8]);
		}

		double VPeakOverRMS[16];
		int waveformLength[16];
		inputFilterTree->SetBranchAddress("VPeakOverRMS",&VPeakOverRMS);
		inputFilterTree->SetBranchAddress("waveformLength",&waveformLength);

		bool passesV;
		bool passesH;
		bool passesBoth;
		inputSurfaceTree->SetBranchAddress("passesV",&passesV);
		inputSurfaceTree->SetBranchAddress("passesH",&passesH);
		inputSurfaceTree->SetBranchAddress("passesBoth",&passesBoth);

		int numEntries = inputVTree->GetEntries();

		inputAllTree->GetEvent(0);

		bool isThisABadRun = isBadRun(station,runNum,BadRunList);
		if(isThisABadRun){
			inputFile->Close();
			delete inputFile;
			continue;
		}
		bool isThisABadSurfaceRun = isBadRun(station,runNum,BadSurfaceRunList);
		if(isThisABadSurfaceRun){
			inputFile->Close();
			delete inputFile;
			continue;
		}		

		int numHPolIsol=0;
		int numVPolIsol=0;
		int numCombo=0;
		int numEither=0;

		for(int event=0; event<numEntries; event++){
			inputVTree->GetEvent(event);
			inputHTree->GetEvent(event);
			inputAllTree->GetEvent(event);
			inputRecoTree->GetEvent(event);
			inputFilterTree->GetEvent(event);
			inputSurfaceTree->GetEvent(event);

			if(passesV && !passesH){
				numVPolIsol++;
				h1_events_vs_time[0]->Fill(unixTime);
				h1_events_vs_time_finer[0]->Fill(unixTime);
				h1_events_vs_run[0]->Fill(runNum);
				// PlotThisEvent(station, runNum, eventNumber, 0);
			}
			if(passesH && !passesV){
				numHPolIsol++;
				h1_events_vs_time[1]->Fill(unixTime);
				h1_events_vs_time_finer[1]->Fill(unixTime);
				h1_events_vs_run[1]->Fill(runNum);
				// PlotThisEvent(station, runNum, eventNumber, 1);
			}

			if(passesBoth){
				h1_events_vs_time[2]->Fill(unixTime);
				h1_events_vs_time_finer[2]->Fill(unixTime);
				h1_events_vs_run[2]->Fill(runNum);
				numCombo++;
				// PlotThisEvent(station, runNum, eventNumber, 2);
			}

			if(passesV || passesH){
				numEither++;
			}

		}// loop over events

		// if(numEither>=14){
			// printf("Run %4d, Num VPol Isol, Num HPol Isol, Num Combo, Num Either: %3d, %3d, %3d, %3d \n", runNum, numVPolIsol, numHPolIsol, numCombo, numEither);
			// printf("%4d, %3d \n", runNum, numEither);
		// }
		// printf("Run %4d, Num VPol Isol, Num HPol Isol, Num Combo, Num Either: %3d, %3d, %3d, %3d \n", runNum, numVPolIsol, numHPolIsol, numCombo, numEither);
		// distro_of_events->Fill(numEither);


		if(numHPolIsol==0 && numVPolIsol==1 && numCombo==0){
			
			// now, loop over again and plot
			for(int event=0; event<numEntries; event++){
				inputAllTree->GetEvent(event);
				inputVTree->GetEvent(event);
				inputHTree->GetEvent(event);
				// v_hpol_isol_unixtime.push_back(unixTime);
				// v_hpol_isol_entries.push_back(1);
				// v_hpol_isol_runs.push_back(runNum);
				// PlotThisEvent(station, runNum, eventNumber, 1);

				h1_PhiDist[0]->Fill(phi_300[0]);
				double thisZenith = getAntennaArrivalTheta(station, settings, target, raysolver, stationCenter, theta_300[0], phi_300[0]);
				// h1_ThetaDist[0]->Fill(theta_300[0]);
				h1_ThetaDist[0]->Fill(thisZenith);
				
				// v_hpol_isol_phi.push_back(phi_300[1]);
				// v_hpol_isol_theta.push_back(theta_300[1]);
			}
		}

		if(numHPolIsol==1 && numVPolIsol==0 && numCombo==0){
			// printf("HPol Isol: Run %d, V, H, Combo: %d, %d, %d \n", runNum, numVPolIsol, numHPolIsol, numCombo);
			
			// now, loop over again and plot
			for(int event=0; event<numEntries; event++){
				// printf("%d, %d \n", runNum, eventNumber);
				// printf("Run %4d, Event %6d, Pol %d, unixTime %d, theta %3d, phi %3d, coor %.4f, CW status %d\n", runNum, eventNumber, 1, unixTime, theta_300[1], phi_300[1], corr_val[1], Refilt[1]);
				// fprintf(fout,"Run %4d, Event %6d, Pol %1d, unixTime %d, theta %3d, phi %3d, coor %.4f, refilt %d, surfv %d, surfh %d, surftopppol %d \n",runNum, eventNumber, pol, unixTime, theta_300[pol], phi_300[pol], corr_val[pol], Refilt[pol], isSurf[0], isSurf[1], isSurfEvent_top[pol]);
				inputAllTree->GetEvent(event);
				inputVTree->GetEvent(event);
				inputHTree->GetEvent(event);
				v_hpol_isol_unixtime.push_back(unixTime);
				v_hpol_isol_entries.push_back(1);
				v_hpol_isol_runs.push_back(runNum);
				// PlotThisEvent(station, runNum, eventNumber, 1);

				h1_PhiDist[1]->Fill(phi_300[1]);
				double thisZenith = getAntennaArrivalTheta(station, settings, target, raysolver, stationCenter, theta_300[1], phi_300[1]);
				// h1_ThetaDist[1]->Fill(theta_300[1]);
				h1_ThetaDist[1]->Fill(thisZenith);

			}
		}

		if(numHPolIsol==0 && numVPolIsol==0 && numCombo==1){
			// printf("Combo Isol: Run %d, V, H, Combo: %d, %d, %d \n", runNum, numVPolIsol, numHPolIsol, numCombo);
			
			// now, loop over again and plot
			for(int event=0; event<numEntries; event++){
				// printf("%d, %d \n", runNum, eventNumber);
				// printf("Run %4d, Event %6d, unixTime %d, V theta %3d, V phi %3d, V coor %.4f, H theta %3d, H phi %3d, H corr %.4f, V CW status %d, H CW Status %d, Surf V %d, Surf H %d, Top Surf V %d, Top Surf H %d \n", runNum, eventNumber, unixTime, theta_300[0], phi_300[0], corr_val[0], theta_300[1], phi_300[1], corr_val[1], Refilt[0] , Refilt[1], isSurf[0], isSurf[1], isSurfEvent_top[0], isSurfEvent_top[1]);
				inputAllTree->GetEvent(event);
				v_combo_isol_unixtime.push_back(unixTime);
				v_combo_isol_entries.push_back(1);
				v_combo_isol_runs.push_back(runNum);
				// PlotThisEvent(station, runNum, eventNumber, 2);

				h1_PhiDist[2]->Fill(phi_300[0]);
				h1_PhiDist[3]->Fill(phi_300[1]);

				double thisZenithV = getAntennaArrivalTheta(station, settings, target, raysolver, stationCenter, theta_300[0], phi_300[0]);
				double thisZenithH = getAntennaArrivalTheta(station, settings, target, raysolver, stationCenter, theta_300[1], phi_300[1]);

				h1_ThetaDist[2]->Fill(thisZenithV);
				h1_ThetaDist[3]->Fill(thisZenithH);
			}
		}

		inputFile->Close();
		delete inputFile;

	} // loop over input files

	gStyle->SetOptStat(101110);


	char thistitle[300];
	
	bool plotDistroOfNumSurfaceEvents=false;
	if(plotDistroOfNumSurfaceEvents){
		Int_t nq = 20;
		Double_t xq[nq];  // position where to compute the quantiles in [0,1]
		Double_t yq[nq];  // array to contain the quantiles
		for (Int_t i=0;i<nq;i++) xq[i] = Float_t(i+1)/nq;
		distro_of_events->GetQuantiles(nq,yq,xq);
		for(int i=0; i<nq; i++){
			printf("%.2f quantile is %.2f \n", xq[i], yq[i]);
		}
		distro_of_events->Fit("expo");

		TCanvas *c = new TCanvas("","",850,850);
			distro_of_events->Draw("");
			distro_of_events->GetXaxis()->SetTitle("Number of Surface Events");
			distro_of_events->GetYaxis()->SetTitle("Number of Runs");
			distro_of_events->GetXaxis()->SetRangeUser(0,50);
			distro_of_events->GetYaxis()->SetTitleOffset(1.2);
			// distro_of_events->SetLineWidth(3);
			gPad->SetLogy();
		char title[150];
		sprintf(title,"/users/PAS0654/osu0673/A23_analysis_new2/results/unblind/surface/A%d_NumSurfEventsDistro.png",station);
		// c->SaveAs(title);
	}

	bool plotSurfaceRate=false;
	if(plotSurfaceRate){

		TGraph *g_HpolIsol_vsRun = new TGraph(v_hpol_isol_runs.size(), &v_hpol_isol_runs[0], &v_hpol_isol_entries[0]);
		TGraph *g_ComboIsol_vsRun = new TGraph(v_combo_isol_runs.size(), &v_combo_isol_runs[0], &v_combo_isol_entries[0]);

		TGraph *g_HpolIsol_vsTime = new TGraph(v_hpol_isol_unixtime.size(), &v_hpol_isol_unixtime[0], &v_hpol_isol_entries[0]);
		TGraph *g_ComboIsol_vsTime = new TGraph(v_combo_isol_unixtime.size(), &v_combo_isol_unixtime[0], &v_combo_isol_entries[0]);

		g_HpolIsol_vsTime->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		g_HpolIsol_vsTime->GetXaxis()->SetTimeOffset(0.,"GMT");
		g_HpolIsol_vsTime->GetXaxis()->SetTimeFormat("%y/%m");

		g_ComboIsol_vsTime->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		g_ComboIsol_vsTime->GetXaxis()->SetTimeOffset(0.,"GMT");
		g_ComboIsol_vsTime->GetXaxis()->SetTimeFormat("%y/%m");

		TCanvas *cSurfaceRate_vsTime = new TCanvas("","",4.1*850,3.1*850);
		cSurfaceRate_vsTime->Divide(1,3);
		for(int i=0; i<3; i++){
			cSurfaceRate_vsTime->cd(i+1);
			// num_samps[i]->GetYaxis()->SetRangeUser(0.1,2e7);
			// num_samps[i]->GetXaxis()->SetRangeUser(0,750);

			// gPad->SetTopMargin(0.1);
			// gPad->SetRightMargin(0.03);
			// gPad->SetLeftMargin(0.05);
			// gPad->SetBottomMargin(0.11);

			configure(h1_events_vs_time[i]);
			configure(h1_events_vs_time_finer[i]);

			h1_events_vs_time[i]->Draw("");
			h1_events_vs_time[i]->GetXaxis()->SetTitle("UnixTime");
			h1_events_vs_time[i]->GetYaxis()->SetTitle("Number of Events");
			h1_events_vs_time[i]->GetXaxis()->SetTimeFormat("%y/%m");
			h1_events_vs_time[i]->GetXaxis()->SetNdivisions(24,0,0,false);
			gPad->SetLogy();
			
			h1_events_vs_time_finer[i]->GetXaxis()->SetTitle("UnixTime");
			h1_events_vs_time_finer[i]->GetYaxis()->SetTitle("Number of Events");

			if(i==1){
				g_HpolIsol_vsTime->Draw("samep");
				g_HpolIsol_vsTime->SetMarkerStyle(kFullCircle);
				g_HpolIsol_vsTime->SetMarkerColor(kRed);
				g_HpolIsol_vsTime->SetMarkerSize(3);
			}
			if(i==2){
				g_ComboIsol_vsTime->Draw("samep");
				g_ComboIsol_vsTime->SetMarkerStyle(kFullCircle);
				g_ComboIsol_vsTime->SetMarkerColor(kRed);
				g_ComboIsol_vsTime->SetMarkerSize(3);
			}

		}
		sprintf(thistitle, "%s/unblind/surface/%d.%d.%d_A%d_c%d_Surface_EventsVsTime.png",plotPath,year_now,month_now,day_now,station,config);
		// cSurfaceRate_vsTime->SaveAs(thistitle);
		delete cSurfaceRate_vsTime;


		TCanvas *cSurfaceRate_vsRun = new TCanvas("","",4.1*850,3.1*850);
		cSurfaceRate_vsRun->Divide(1,3);
		for(int i=0; i<3; i++){
			cSurfaceRate_vsRun->cd(i+1);

			configure(h1_events_vs_run[i]);

			h1_events_vs_run[i]->Draw("");
			h1_events_vs_run[i]->GetXaxis()->SetTitle("runNumber");
			h1_events_vs_run[i]->GetYaxis()->SetTitle("Number of Events");
			gPad->SetLogy();

			if(i==1){
				g_HpolIsol_vsRun->Draw("samep");
				g_HpolIsol_vsRun->SetMarkerStyle(kFullCircle);
				g_HpolIsol_vsRun->SetMarkerColor(kRed);
				g_HpolIsol_vsRun->SetMarkerSize(3);
			}
			if(i==2){
				g_ComboIsol_vsRun->Draw("samep");
				g_ComboIsol_vsRun->SetMarkerStyle(kFullCircle);
				g_ComboIsol_vsRun->SetMarkerColor(kRed);
				g_ComboIsol_vsRun->SetMarkerSize(3);
			}
		}
		sprintf(thistitle, "%s/unblind/surface/%d.%d.%d_A%d_c%d_Surface_EventsVsRun.png",plotPath,year_now,month_now,day_now,station,config);
		// cSurfaceRate_vsRun->SaveAs(thistitle);
		delete cSurfaceRate_vsRun;

		bool doVsTime=false;
		if(doVsTime){
			TCanvas *cSurfaceRate_vsRun_vsTime_Zoom = new TCanvas("","",8.1*850,3.1*850);
			cSurfaceRate_vsRun_vsTime_Zoom->Divide(2,3);
			
			for(int i=0; i<3; i++){
				h1_events_vs_time_finer[i]->GetXaxis()->SetTimeFormat("%m/%d");
				h1_events_vs_time_finer[i]->GetXaxis()->SetNdivisions(14,0,0,false);
				h1_events_vs_time_finer[i]->GetYaxis()->SetTitle("Number of Events per Hour");
			}

			for(int i=0; i<3; i++){
				h1_events_vs_time[i]->SetLineWidth(3);
				h1_events_vs_run[i]->SetLineWidth(3);
				h1_events_vs_time_finer[i]->SetLineWidth(3);
			}

			// the left column, which is "vs run"
			cSurfaceRate_vsRun_vsTime_Zoom->cd(1);
				h1_events_vs_run[0]->Draw("");
				gPad->SetLogy();
			cSurfaceRate_vsRun_vsTime_Zoom->cd(3);
				h1_events_vs_run[1]->Draw("");
				gPad->SetLogy();
				g_HpolIsol_vsRun->Draw("samep");
			cSurfaceRate_vsRun_vsTime_Zoom->cd(5);
				h1_events_vs_run[2]->Draw("");
				gPad->SetLogy();
				g_ComboIsol_vsRun->Draw("samep");	
			
			// the right column, which is "vs time"
			cSurfaceRate_vsRun_vsTime_Zoom->cd(2);
				h1_events_vs_time_finer[0]->Draw("");
				gPad->SetLogy();
			cSurfaceRate_vsRun_vsTime_Zoom->cd(4);
				h1_events_vs_time_finer[1]->Draw("");
				gPad->SetLogy();
				g_HpolIsol_vsTime->Draw("samep");
			cSurfaceRate_vsRun_vsTime_Zoom->cd(6);
				h1_events_vs_time_finer[2]->Draw("");
				gPad->SetLogy();
				g_ComboIsol_vsTime->Draw("samep");	

			int numSecDay = 86400;
			int intervalWidth=7*numSecDay;
			
			/*
				Now we loop and maked "zoomed in" versions of all these plots around the events of interest
			*/
			for(int vrun=0; vrun<v_hpol_isol_runs.size(); vrun++){
				int theRun = v_hpol_isol_runs[vrun];
				int theUnixtime = v_hpol_isol_unixtime[vrun];
				for(int i=0; i<3; i++){
					h1_events_vs_run[i]->GetXaxis()->SetRangeUser(theRun-50, theRun+50);
					h1_events_vs_time_finer[i]->GetXaxis()->SetRangeUser(theUnixtime-intervalWidth, theUnixtime+intervalWidth);
					h1_events_vs_time_finer[i]->GetXaxis()->SetTimeFormat("%b %d '%y");
					h1_events_vs_time_finer[i]->GetXaxis()->SetNdivisions(14,0,0,false);
				}
				sprintf(thistitle, "%s/unblind/surface/A%d_HPolOnly_Run%d_NearbySurfaceEvents.png",plotPath,station,theRun);
				// cSurfaceRate_vsRun_vsTime_Zoom->SaveAs(thistitle);
			}

			/*
				Now we loop and maked "zoomed in" versions of all these plots around the events of interest
			*/
			for(int vrun=0; vrun<v_combo_isol_runs.size(); vrun++){
				int theRun = v_combo_isol_runs[vrun];
				int theUnixtime = v_combo_isol_unixtime[vrun];
				for(int i=0; i<3; i++){
					h1_events_vs_run[i]->GetXaxis()->SetRangeUser(theRun-50, theRun+50);
					h1_events_vs_time_finer[i]->GetXaxis()->SetRangeUser(theUnixtime-intervalWidth, theUnixtime+intervalWidth);
					h1_events_vs_time_finer[i]->GetXaxis()->SetTimeFormat("%b %d '%y");
					h1_events_vs_time_finer[i]->GetXaxis()->SetNdivisions(14,0,0,false);
				}
				sprintf(thistitle, "%s/unblind/surface/A%d_HVCombo_Run%d_NearbySurfaceEvents.png",plotPath,station,theRun);
				cSurfaceRate_vsRun_vsTime_Zoom->SaveAs(thistitle);
			}

			delete cSurfaceRate_vsRun_vsTime_Zoom;
		}
	}

	bool doSpatial=true;
	if(doSpatial){

		printf("Num HPol Isolated Events %d \n", (int)h1_ThetaDist[0]->Integral());
		printf("Num VPOl Isolated Events %d \n", (int)h1_ThetaDist[1]->Integral());

		TCanvas *cSpatial_Isol = new TCanvas("","",2.1*850, 2.1*850);
		cSpatial_Isol->Divide(2,3);
		cSpatial_Isol->cd(1);
			h1_ThetaDist[0]->Draw("");
			h1_ThetaDist[0]->SetLineWidth(3);
			h1_ThetaDist[0]->GetXaxis()->SetTitle("Zenith [deg]");
			h1_ThetaDist[0]->GetYaxis()->SetTitle("Number of Events");
			h1_ThetaDist[0]->GetXaxis()->SetRangeUser(0,90);
			gPad->SetLogy();
		cSpatial_Isol->cd(3);
			h1_ThetaDist[1]->Draw("");
			h1_ThetaDist[1]->SetLineWidth(3);
			h1_ThetaDist[1]->GetXaxis()->SetTitle("Zenith [deg]");
			h1_ThetaDist[1]->GetYaxis()->SetTitle("Number of Events");
			h1_ThetaDist[1]->GetXaxis()->SetRangeUser(0,90);
		cSpatial_Isol->cd(5);
			h1_ThetaDist[2]->Draw("");
				h1_ThetaDist[2]->SetLineWidth(3);
				h1_ThetaDist[2]->GetXaxis()->SetTitle("Zenith [deg]");
				h1_ThetaDist[2]->GetYaxis()->SetTitle("Number of Events");
				h1_ThetaDist[2]->GetXaxis()->SetRangeUser(0,90);
				h1_ThetaDist[2]->SetTitle("Isolated H&V Combination Events Theta");
			h1_ThetaDist[3]->Draw("same");
				h1_ThetaDist[3]->SetLineWidth(3);
				h1_ThetaDist[3]->SetLineColor(kRed);
			{
				TLegend *leg = new TLegend(0.48,0.6,0.9,0.9);
				leg->AddEntry(h1_ThetaDist[2],"VPol","l");
				leg->AddEntry(h1_ThetaDist[3],"HPol","l");
				// leg->AddEntry(hRMS_stragglers[i],"Straggling Events","l");
				leg->Draw();
			}

		cSpatial_Isol->cd(2);
			h1_PhiDist[0]->Draw("");
			h1_PhiDist[0]->SetLineWidth(3);
			h1_PhiDist[0]->GetXaxis()->SetTitle("Phi [deg]");
			h1_PhiDist[0]->GetYaxis()->SetTitle("Number of Events");
			gPad->SetLogy();
		cSpatial_Isol->cd(4);
			h1_PhiDist[1]->Draw("");
			h1_PhiDist[1]->SetLineWidth(3);
			h1_PhiDist[1]->GetXaxis()->SetTitle("Phi [deg]");
			h1_PhiDist[1]->GetYaxis()->SetTitle("Number of Events");
		cSpatial_Isol->cd(6);
			h1_PhiDist[2]->Draw("");
				h1_PhiDist[2]->SetLineWidth(3);
				h1_PhiDist[2]->GetXaxis()->SetTitle("Phi [deg]");
				h1_PhiDist[2]->GetYaxis()->SetTitle("Number of Events");
				h1_PhiDist[2]->SetTitle("Isolated H&V Combination Events Phi");
			h1_PhiDist[3]->Draw("same");
				h1_PhiDist[3]->SetLineWidth(3);
				h1_PhiDist[3]->SetLineColor(kRed);
		sprintf(thistitle, "%s/unblind/surface/A%d_SurfaceSpatialDistribution_IsolatedEvents_ExcludingBursts.png",plotPath,station);
		// cSpatial_Isol->SaveAs(thistitle);
		delete cSpatial_Isol;
	}
}

int PlotThisEvent(int station, int runNum, int event, int problempol){
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *DataDirPath(getenv("DATA_DIR_100"));
	if (DataDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

	char run_file_name[400];
	sprintf(run_file_name,"%s/RawData/A%d/all_runs/event%d.root",DataDirPath,station,runNum);
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

	AraQualCuts *qualCut = AraQualCuts::Instance();
	UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawPtr,AraCalType::kLatestCalib);
	printf("Run %d, Event %d \n", runNum, realAtriEvPtr->eventNumber);
	printf("	Is Quality Event? %d \n", qualCut->isGoodEvent(realAtriEvPtr));

	int unixTime = (int)rawPtr->unixTime;
	int unixTimeUs =(int)rawPtr->unixTimeUs;
	int timeStamp = (int)rawPtr->timeStamp;
	printf("	Unixtime is %d \n", unixTime);
	printf("	Unixtime microsecond is %d \n", unixTimeUs);
	printf("	timeStamp is %d \n", timeStamp);

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

	
	bool doMapsAtAll=false;
	if(doMapsAtAll){
		// get the run summary information, if it exists yet
		// and remember, because it's the users job to pass the location of the filter files
		// this should work for simulated events just fine
		char filter_file_name[400];
		// gonna try all of them cuz we lazy af lol
		sprintf(filter_file_name,"%s/processed_station_%d_run_%d_filter.root","/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ProcessedFile/A3/2015",station,runNum);
		bool hasFilterFile = false;
		TFile *filterFile = TFile::Open(filter_file_name);
		if(!filterFile){
			sprintf(filter_file_name,"%s/processed_station_%d_run_%d_filter.root","/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ProcessedFile/A3/2013",station,runNum);
		}
		filterFile = TFile::Open(filter_file_name);
		if(!filterFile){
			sprintf(filter_file_name,"%s/processed_station_%d_run_%d_filter.root","/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ProcessedFile/A3/2014",station,runNum);
		}
		filterFile = TFile::Open(filter_file_name);
		if(!filterFile){
			sprintf(filter_file_name,"%s/processed_station_%d_run_%d_filter.root","/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ProcessedFile/A3/2016",station,runNum);
		}
		filterFile = TFile::Open(filter_file_name);
		if(!filterFile){
			return -1;
		}

		TTree *filterTree;
		double VPeakOverRMS[16];
		if(filterFile){
			printf("Successfully found filter file information \n");
			hasFilterFile=true;
			filterTree = (TTree*) filterFile->Get("OutputTree");
			if(!filterTree) {
				std::cout << "Can't find filterTree\n";
				return -1;
			}
			filterTree->SetBranchAddress("VPeakOverRMS", &VPeakOverRMS);
			filterFile->cd();
		}


		filterTree->GetEvent(event);

		vector<double> chan_SNRs;
		if(hasFilterFile){
			for(int i=0; i<16; i++){
				chan_SNRs.push_back(VPeakOverRMS[i]);
			}
		}

		vector <int> chan_list_V;
		vector <int> chan_list_H;
		for(int chan=0; chan<=7; chan++){
			chan_list_V.push_back(chan);
			chan_list_H.push_back(chan+8);
		}
		
		isSimulation=false;
		if(station==2){
			//for station 2, we need to exclude channel 15 from the analysis
			chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
		}
		else if(station==3){
			// for station 3 years 2014, 2015, 2016, we need to drop string 4 (channels 3, 7, 11, 15) altogether above some run
			if( 
				(!isSimulation && runNum>getA3BadRunBoundary())

			){			// drop string four
				chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 3), chan_list_V.end());
				chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 7), chan_list_V.end());
				chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 11), chan_list_H.end());
				chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
			}
		}

		Settings *settings = new Settings();
		string setupfile = "setup.txt";
		settings->ReadFile(setupfile);
		cout << "Read " << setupfile << " file!" << endl;
		settings->NOFZ=1;
		Detector *detector = 0;

		RayTraceCorrelator *theCorrelator = new RayTraceCorrelator(station, 300, settings, 1, RTTestMode);

		int solNum = 0;
		TH2D *map_V_raytrace = theCorrelator->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V, chan_SNRs, solNum);
		TH2D *map_H_raytrace = theCorrelator->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H, chan_SNRs, solNum);

		int new_theta;
		int new_phi;
		double new_corr;

		getCorrMapPeak(map_V_raytrace,new_theta, new_phi, new_corr);
		char title_for_map[300];
		sprintf(title_for_map,"VMap: peak theta %d, phi %d, corr %.4f",new_theta, new_phi,new_corr);
		map_V_raytrace->SetTitle(title_for_map);

		int new_theta_H;
		int new_phi_H;
		double new_corr_H;

		getCorrMapPeak(map_H_raytrace,new_theta_H, new_phi_H, new_corr_H);
		char title_for_map_H[300];
		sprintf(title_for_map_H,"HMap: peak theta %d, phi %d, corr %.4f",new_theta_H, new_phi_H,new_corr_H);
		map_H_raytrace->SetTitle(title_for_map_H);

		bool print_maps = true;
		if(print_maps){
			gStyle->SetOptStat(0);
			TCanvas *cMaps = new TCanvas("","",2*1100,850);
			// TCanvas *cMaps = new TCanvas("","",1100,1.1*850);
			cMaps->Divide(2,1);
				cMaps->cd(1);
				map_V_raytrace->Draw("colz");
				gPad->SetRightMargin(0.15);
				cMaps->cd(2);
				map_H_raytrace->Draw("colz");
				gPad->SetRightMargin(0.15);
			char save_temp_title[400];
			sprintf(save_temp_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/unblind/surface/all_events/pol%d/A%d_Run%d_Ev%d_ProblemPol%d_Maps.png",problempol,station,runNum,event, problempol);
			cMaps->SaveAs(save_temp_title);
			delete cMaps;
		}

		filterFile->Close();
	} // do maps at all

	// bool doContributingMaps=true;
	// if(doContributingMaps){
	// 	stringstream ss1;
	// 	vector<string> titlesForGraphs;
	// 	vector <TH2D*> individuals;

	// 	double SNR_scaling=0.;

	// 	for(int i=0; i<7; i++){
	// 		for(int j=i+1; j<8; j++){
	// 			ss1.str("");
	// 			ss1<<"Pair "<<i<<" and "<<j;
	// 			titlesForGraphs.push_back(ss1.str());
	// 			TH2D *map = theCorrelator->getInterferometricMap_RT_NewNormalization_PairSelect(settings, detector, realAtriEvPtr, Vpol, isSimulation, i, j, solNum);
				
	// 			// now do the SNR scaling
	// 			double this_snr_product = chan_SNRs[i] * chan_SNRs[j];
	// 			// printf("Weighting term for %d, %d is %.3f \n", i, j, this_snr_product);
	// 			map->Scale(this_snr_product);

	// 			SNR_scaling+=this_snr_product;

	// 			individuals.push_back(map);
	// 		}
	// 	}

	// 	for(int i=0; i<individuals.size(); i++){
	// 		individuals[i]->Scale(1./SNR_scaling);
	// 	}

	// 	// compute the average myself manually
	// 	TH2D *average = (TH2D*) individuals[0]->Clone();
	// 	for(int i=0; i<individuals.size(); i++){
	// 		average->Add(individuals[i]);
	// 	}
	// 	// average->Scale(1./28);
	// 	average->SetTitle("Summed Maps");

	// 	vector<double> mins;
	// 	vector<double> maxs;
	// 	for(int i=0; i<individuals.size(); i++){
	// 		mins.push_back(individuals[i]->GetMinimum());
	// 		maxs.push_back(individuals[i]->GetMaximum());
	// 	}
	// 	std::sort(mins.begin(), mins.end()); //sort smallest to largest
	// 	std::sort(maxs.begin(), maxs.end()); //sort smallest to largest
	// 	std::reverse(maxs.begin(), maxs.end()); //reverse order to get largest to smallest

	// 	TCanvas *cMaps2 = new TCanvas("","",8*850,4*850);
	// 	cMaps2->Divide(7,5);
	// 	for(int i=0; i<individuals.size(); i++){
	// 		cMaps2->cd(i+1);
	// 		individuals[i]->Draw("colz");
	// 		individuals[i]->GetZaxis()->SetRangeUser(mins[0],maxs[0]);
	// 		individuals[i]->SetTitle(titlesForGraphs[i].c_str());
	// 		gStyle->SetTitleFontSize(0.07);
	// 	}
	// 	cMaps2->cd(35);
	// 		average->Draw("colz");
	// 	char save_temp_title[400];
	// 	sprintf(save_temp_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/unblind/surface/surface_%d.%d.%d_Run%d_Ev%d_AllMaps.png",year_now,month_now,day_now,runNum,event);
	// 	cMaps2->SaveAs(save_temp_title);
	// 	delete cMaps2;
	// }

	
	// need dummy canvases becaus root can't make things pretty to save its life
	vector<TGraph*> dummy;
	for(int i=0; i<16; i++){
		vector<double> thisX;
		vector<double> thisY;
		thisX.push_back(-200);
		thisX.push_back(800);
		thisY.push_back(-800);
		thisY.push_back(800);
		dummy.push_back(new TGraph(2,&thisX[0], &thisY[0]));
	}

	char save_temp_title[300];
	sprintf(save_temp_title,"%s/unblind/surface/all_events/pol%d/A%d_Run%d_Ev%d_ProblemPol%d_Waveforms.png",plotPath,problempol,station,runNum,event,problempol);
	TCanvas *cWave = new TCanvas("","",4*1100,4*850);
	cWave->Divide(4,4);
	for(int i=0; i<16; i++){
		cWave->cd(i+1);
		dummy[i]->Draw("AP");
		dummy[i]->SetLineColor(kWhite);
		dummy[i]->GetXaxis()->SetRangeUser(-200.,700.);
		dummy[i]->GetXaxis()->SetRangeUser(-700.,700.);

		waveforms[i]->Draw("sameL");
		waveforms[i]->SetLineWidth(3);
	}
	cWave->SaveAs(save_temp_title);
	delete cWave;

	sprintf(save_temp_title,"%s/unblind/surface/all_events/pol%d/A%d_Run%d_Ev%d_ProblemPol%d_Spectra.png",plotPath,problempol,station,runNum,event,problempol);
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