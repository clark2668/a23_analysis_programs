////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	proposal_passing_distro.cxx
////	make theta and phi distro for v vs h isolate for proposal
////
////	Nov 2019
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
#include "TArrow.h"
#include "TLatex.h"

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

double getAntennaArrivalTheta(int station, Settings *settings, Position target, RaySolver *raysolver, double stationCenter[3], int thetaIn, int phiIn){
	double phiWave = (double(phiIn))*TMath::DegToRad();
	double thetaWave = (double(thetaIn))*TMath::DegToRad();
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

	cout<<"Station depth is "<<stationCenter[2]<<endl;

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

	stringstream ss;

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

	vector<int> BadRunList=BuildBadRunList(station);
	vector<int> BadSurfaceRunList=BuildSurfaceRunList(station);

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

		inputAllTree->SetBranchAddress("surf_V_new",&isSurf[0]);
		inputAllTree->SetBranchAddress("surf_H_new",&isSurf[1]);

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

		for(int event=0; event<numEntries; event++){
			inputVTree->GetEvent(event);
			inputHTree->GetEvent(event);
			inputAllTree->GetEvent(event);
			inputRecoTree->GetEvent(event);
			inputFilterTree->GetEvent(event);
			inputSurfaceTree->GetEvent(event);

			if(passesV && !passesH){
				numVPolIsol++;
				// h1_events_vs_time[0]->Fill(unixTime);
				// h1_events_vs_time_finer[0]->Fill(unixTime);
				// h1_events_vs_run[0]->Fill(runNum);
				// PlotThisEvent(station, runNum, eventNumber, 0);
			}
			if(passesH && !passesV){
				numHPolIsol++;
				// h1_events_vs_time[1]->Fill(unixTime);
				// h1_events_vs_time_finer[1]->Fill(unixTime);
				// h1_events_vs_run[1]->Fill(runNum);
				// PlotThisEvent(station, runNum, eventNumber, 1);
			}

			if(passesBoth){
				// h1_events_vs_time[2]->Fill(unixTime);
				// h1_events_vs_time_finer[2]->Fill(unixTime);
				// h1_events_vs_run[2]->Fill(runNum);
				numCombo++;
				// PlotThisEvent(station, runNum, eventNumber, 2);
			}
		}// loop over events

		// printf("Run %4d, Num VPol Isol, Num HPol Isol, Num Combo: %3d, %3d, %3d \n", runNum, numVPolIsol, numHPolIsol, numCombo);


		if(numHPolIsol==0 && numVPolIsol==1 && numCombo==0){
			
			// now, loop over again and plot
			for(int event=0; event<numEntries; event++){
				inputAllTree->GetEvent(event);
				inputVTree->GetEvent(event);
				inputHTree->GetEvent(event);

				h1_PhiDist[0]->Fill(phi_300[0]);
				double thisZenith = getAntennaArrivalTheta(station, settings, target, raysolver, stationCenter, theta_300[0], phi_300[0]);
				// h1_ThetaDist[0]->Fill(theta_300[0]);
				h1_ThetaDist[0]->Fill(thisZenith);
				
				// v_hpol_isol_phi.push_back(phi_300[1]);
				// v_hpol_isol_theta.push_back(theta_300[1]);
			}
		}

		if(numHPolIsol==1 && numVPolIsol==0 && numCombo==0){
			cout<<"I found an hpol isolated event"<<endl;
			// printf("HPol Isol: Run %d, V, H, Combo: %d, %d, %d \n", runNum, numVPolIsol, numHPolIsol, numCombo);
			
			// now, loop over again and plot
			for(int event=0; event<numEntries; event++){
				// printf("%d, %d \n", runNum, eventNumber);
				// printf("Run %4d, Event %6d, Pol %d, unixTime %d, theta %3d, phi %3d, coor %.4f, CW status %d\n", runNum, eventNumber, 1, unixTime, theta_300[1], phi_300[1], corr_val[1], Refilt[1]);
				// fprintf(fout,"Run %4d, Event %6d, Pol %1d, unixTime %d, theta %3d, phi %3d, coor %.4f, refilt %d, surfv %d, surfh %d, surftopppol %d \n",runNum, eventNumber, pol, unixTime, theta_300[pol], phi_300[pol], corr_val[pol], Refilt[pol], isSurf[0], isSurf[1], isSurfEvent_top[pol]);
				inputAllTree->GetEvent(event);
				inputVTree->GetEvent(event);
				inputHTree->GetEvent(event);

				h1_PhiDist[1]->Fill(phi_300[1]);
				double thisZenith = getAntennaArrivalTheta(station, settings, target, raysolver, stationCenter, theta_300[1], phi_300[1]);
				// h1_ThetaDist[1]->Fill(theta_300[1]);
				h1_ThetaDist[1]->Fill(thisZenith);

			}
		}

		// if(numHPolIsol==0 && numVPolIsol==0 && numCombo==1){
		// 	// printf("Combo Isol: Run %d, V, H, Combo: %d, %d, %d \n", runNum, numVPolIsol, numHPolIsol, numCombo);
			
		// 	// now, loop over again and plot
		// 	for(int event=0; event<numEntries; event++){
		// 		// printf("%d, %d \n", runNum, eventNumber);
		// 		// printf("Run %4d, Event %6d, unixTime %d, V theta %3d, V phi %3d, V coor %.4f, H theta %3d, H phi %3d, H corr %.4f, V CW status %d, H CW Status %d, Surf V %d, Surf H %d, Top Surf V %d, Top Surf H %d \n", runNum, eventNumber, unixTime, theta_300[0], phi_300[0], corr_val[0], theta_300[1], phi_300[1], corr_val[1], Refilt[0] , Refilt[1], isSurf[0], isSurf[1], isSurfEvent_top[0], isSurfEvent_top[1]);
		// 		inputAllTree->GetEvent(event);
		// 		inputVTree->GetEvent(event);
		// 		inputHTree->GetEvent(event);

		// 		h1_PhiDist[2]->Fill(phi_300[0]);
		// 		h1_PhiDist[3]->Fill(phi_300[1]);

		// 		double thisZenithV = getAntennaArrivalTheta(station, settings, target, raysolver, stationCenter, theta_300[0], phi_300[0]);
		// 		double thisZenithH = getAntennaArrivalTheta(station, settings, target, raysolver, stationCenter, theta_300[1], phi_300[1]);

		// 		h1_ThetaDist[2]->Fill(thisZenithV);
		// 		h1_ThetaDist[3]->Fill(thisZenithH);
		// 	}
		// }

		inputFile->Close();
		delete inputFile;

	} // loop over input files

	char thistitle[300];

	// bool plotSurfaceRate=false;
	// if(plotSurfaceRate){

	// 	TGraph *g_HpolIsol_vsRun = new TGraph(v_hpol_isol_runs.size(), &v_hpol_isol_runs[0], &v_hpol_isol_entries[0]);
	// 	TGraph *g_ComboIsol_vsRun = new TGraph(v_combo_isol_runs.size(), &v_combo_isol_runs[0], &v_combo_isol_entries[0]);

	// 	TGraph *g_HpolIsol_vsTime = new TGraph(v_hpol_isol_unixtime.size(), &v_hpol_isol_unixtime[0], &v_hpol_isol_entries[0]);
	// 	TGraph *g_ComboIsol_vsTime = new TGraph(v_combo_isol_unixtime.size(), &v_combo_isol_unixtime[0], &v_combo_isol_entries[0]);

	// 	g_HpolIsol_vsTime->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
	// 	g_HpolIsol_vsTime->GetXaxis()->SetTimeOffset(0.,"GMT");
	// 	g_HpolIsol_vsTime->GetXaxis()->SetTimeFormat("%y/%m");

	// 	g_ComboIsol_vsTime->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
	// 	g_ComboIsol_vsTime->GetXaxis()->SetTimeOffset(0.,"GMT");
	// 	g_ComboIsol_vsTime->GetXaxis()->SetTimeFormat("%y/%m");

	// 	TCanvas *cSurfaceRate_vsTime = new TCanvas("","",4.1*850,3.1*850);
	// 	cSurfaceRate_vsTime->Divide(1,3);
	// 	for(int i=0; i<3; i++){
	// 		cSurfaceRate_vsTime->cd(i+1);
	// 		// num_samps[i]->GetYaxis()->SetRangeUser(0.1,2e7);
	// 		// num_samps[i]->GetXaxis()->SetRangeUser(0,750);

	// 		// gPad->SetTopMargin(0.1);
	// 		// gPad->SetRightMargin(0.03);
	// 		// gPad->SetLeftMargin(0.05);
	// 		// gPad->SetBottomMargin(0.11);

	// 		configure(h1_events_vs_time[i]);
	// 		configure(h1_events_vs_time_finer[i]);

	// 		h1_events_vs_time[i]->Draw("");
	// 		h1_events_vs_time[i]->GetXaxis()->SetTitle("UnixTime");
	// 		h1_events_vs_time[i]->GetYaxis()->SetTitle("Number of Events");
	// 		h1_events_vs_time[i]->GetXaxis()->SetTimeFormat("%y/%m");
	// 		h1_events_vs_time[i]->GetXaxis()->SetNdivisions(24,0,0,false);
	// 		gPad->SetLogy();
			
	// 		h1_events_vs_time_finer[i]->GetXaxis()->SetTitle("UnixTime");
	// 		h1_events_vs_time_finer[i]->GetYaxis()->SetTitle("Number of Events");

	// 		if(i==1){
	// 			g_HpolIsol_vsTime->Draw("samep");
	// 			g_HpolIsol_vsTime->SetMarkerStyle(kFullCircle);
	// 			g_HpolIsol_vsTime->SetMarkerColor(kRed);
	// 			g_HpolIsol_vsTime->SetMarkerSize(3);
	// 		}
	// 		if(i==2){
	// 			g_ComboIsol_vsTime->Draw("samep");
	// 			g_ComboIsol_vsTime->SetMarkerStyle(kFullCircle);
	// 			g_ComboIsol_vsTime->SetMarkerColor(kRed);
	// 			g_ComboIsol_vsTime->SetMarkerSize(3);
	// 		}

	// 	}
	// 	sprintf(thistitle, "%s/unblind/surface/%d.%d.%d_A%d_c%d_Surface_EventsVsTime.png",plotPath,year_now,month_now,day_now,station,config);
	// 	// cSurfaceRate_vsTime->SaveAs(thistitle);
	// 	delete cSurfaceRate_vsTime;


	// 	TCanvas *cSurfaceRate_vsRun = new TCanvas("","",4.1*850,3.1*850);
	// 	cSurfaceRate_vsRun->Divide(1,3);
	// 	for(int i=0; i<3; i++){
	// 		cSurfaceRate_vsRun->cd(i+1);

	// 		configure(h1_events_vs_run[i]);

	// 		h1_events_vs_run[i]->Draw("");
	// 		h1_events_vs_run[i]->GetXaxis()->SetTitle("runNumber");
	// 		h1_events_vs_run[i]->GetYaxis()->SetTitle("Number of Events");
	// 		gPad->SetLogy();

	// 		if(i==1){
	// 			g_HpolIsol_vsRun->Draw("samep");
	// 			g_HpolIsol_vsRun->SetMarkerStyle(kFullCircle);
	// 			g_HpolIsol_vsRun->SetMarkerColor(kRed);
	// 			g_HpolIsol_vsRun->SetMarkerSize(3);
	// 		}
	// 		if(i==2){
	// 			g_ComboIsol_vsRun->Draw("samep");
	// 			g_ComboIsol_vsRun->SetMarkerStyle(kFullCircle);
	// 			g_ComboIsol_vsRun->SetMarkerColor(kRed);
	// 			g_ComboIsol_vsRun->SetMarkerSize(3);
	// 		}
	// 	}
	// 	sprintf(thistitle, "%s/unblind/surface/%d.%d.%d_A%d_c%d_Surface_EventsVsRun.png",plotPath,year_now,month_now,day_now,station,config);
	// 	// cSurfaceRate_vsRun->SaveAs(thistitle);
	// 	delete cSurfaceRate_vsRun;

	// 	bool doVsTime=false;
	// 	if(doVsTime){
	// 		TCanvas *cSurfaceRate_vsRun_vsTime_Zoom = new TCanvas("","",8.1*850,3.1*850);
	// 		cSurfaceRate_vsRun_vsTime_Zoom->Divide(2,3);
			
	// 		for(int i=0; i<3; i++){
	// 			h1_events_vs_time_finer[i]->GetXaxis()->SetTimeFormat("%m/%d");
	// 			h1_events_vs_time_finer[i]->GetXaxis()->SetNdivisions(14,0,0,false);
	// 			h1_events_vs_time_finer[i]->GetYaxis()->SetTitle("Number of Events per Hour");
	// 		}

	// 		for(int i=0; i<3; i++){
	// 			h1_events_vs_time[i]->SetLineWidth(3);
	// 			h1_events_vs_run[i]->SetLineWidth(3);
	// 			h1_events_vs_time_finer[i]->SetLineWidth(3);
	// 		}

	// 		// the left column, which is "vs run"
	// 		cSurfaceRate_vsRun_vsTime_Zoom->cd(1);
	// 			h1_events_vs_run[0]->Draw("");
	// 			gPad->SetLogy();
	// 		cSurfaceRate_vsRun_vsTime_Zoom->cd(3);
	// 			h1_events_vs_run[1]->Draw("");
	// 			gPad->SetLogy();
	// 			g_HpolIsol_vsRun->Draw("samep");
	// 		cSurfaceRate_vsRun_vsTime_Zoom->cd(5);
	// 			h1_events_vs_run[2]->Draw("");
	// 			gPad->SetLogy();
	// 			g_ComboIsol_vsRun->Draw("samep");	
			
	// 		// the right column, which is "vs time"
	// 		cSurfaceRate_vsRun_vsTime_Zoom->cd(2);
	// 			h1_events_vs_time_finer[0]->Draw("");
	// 			gPad->SetLogy();
	// 		cSurfaceRate_vsRun_vsTime_Zoom->cd(4);
	// 			h1_events_vs_time_finer[1]->Draw("");
	// 			gPad->SetLogy();
	// 			g_HpolIsol_vsTime->Draw("samep");
	// 		cSurfaceRate_vsRun_vsTime_Zoom->cd(6);
	// 			h1_events_vs_time_finer[2]->Draw("");
	// 			gPad->SetLogy();
	// 			g_ComboIsol_vsTime->Draw("samep");	

	// 		int numSecDay = 86400;
	// 		int intervalWidth=7*numSecDay;
			
	// 		/*
	// 			Now we loop and maked "zoomed in" versions of all these plots around the events of interest
	// 		*/
	// 		for(int vrun=0; vrun<v_hpol_isol_runs.size(); vrun++){
	// 			int theRun = v_hpol_isol_runs[vrun];
	// 			int theUnixtime = v_hpol_isol_unixtime[vrun];
	// 			for(int i=0; i<3; i++){
	// 				h1_events_vs_run[i]->GetXaxis()->SetRangeUser(theRun-50, theRun+50);
	// 				h1_events_vs_time_finer[i]->GetXaxis()->SetRangeUser(theUnixtime-intervalWidth, theUnixtime+intervalWidth);
	// 				h1_events_vs_time_finer[i]->GetXaxis()->SetTimeFormat("%b %d '%y");
	// 				h1_events_vs_time_finer[i]->GetXaxis()->SetNdivisions(14,0,0,false);
	// 			}
	// 			sprintf(thistitle, "%s/unblind/surface/A%d_HPolOnly_Run%d_NearbySurfaceEvents.png",plotPath,station,theRun);
	// 			// cSurfaceRate_vsRun_vsTime_Zoom->SaveAs(thistitle);
	// 		}

	// 		/*
	// 			Now we loop and maked "zoomed in" versions of all these plots around the events of interest
	// 		*/
	// 		for(int vrun=0; vrun<v_combo_isol_runs.size(); vrun++){
	// 			int theRun = v_combo_isol_runs[vrun];
	// 			int theUnixtime = v_combo_isol_unixtime[vrun];
	// 			for(int i=0; i<3; i++){
	// 				h1_events_vs_run[i]->GetXaxis()->SetRangeUser(theRun-50, theRun+50);
	// 				h1_events_vs_time_finer[i]->GetXaxis()->SetRangeUser(theUnixtime-intervalWidth, theUnixtime+intervalWidth);
	// 				h1_events_vs_time_finer[i]->GetXaxis()->SetTimeFormat("%b %d '%y");
	// 				h1_events_vs_time_finer[i]->GetXaxis()->SetNdivisions(14,0,0,false);
	// 			}
	// 			sprintf(thistitle, "%s/unblind/surface/A%d_HVCombo_Run%d_NearbySurfaceEvents.png",plotPath,station,theRun);
	// 			cSurfaceRate_vsRun_vsTime_Zoom->SaveAs(thistitle);
	// 		}

	// 		delete cSurfaceRate_vsRun_vsTime_Zoom;
	// 	}
	// }

	bool doSpatial=true;
	if(doSpatial){

		printf("Num VPol Isolated is %d \n", int(h1_ThetaDist[0]->Integral()));
		printf("Num HPOl Isolated is %d \n", int(h1_ThetaDist[1]->Integral()));

		TCanvas *cSpatial_Isol = new TCanvas("","",2.1*850, 2.1*850);
		cSpatial_Isol->Divide(2,2);
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
		// cSpatial_Isol->cd(5);
		// 	h1_ThetaDist[2]->Draw("");
		// 		h1_ThetaDist[2]->SetLineWidth(3);
		// 		h1_ThetaDist[2]->GetXaxis()->SetTitle("Zenith [deg]");
		// 		h1_ThetaDist[2]->GetYaxis()->SetTitle("Number of Events");
		// 		h1_ThetaDist[2]->GetXaxis()->SetRangeUser(0,90);
		// 		h1_ThetaDist[2]->SetTitle("Isolated H&V Combination Events Theta");
		// 	h1_ThetaDist[3]->Draw("same");
		// 		h1_ThetaDist[3]->SetLineWidth(3);
		// 		h1_ThetaDist[3]->SetLineColor(kRed);
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
		// cSpatial_Isol->cd(6);
		// 	h1_PhiDist[2]->Draw("");
		// 		h1_PhiDist[2]->SetLineWidth(3);
		// 		h1_PhiDist[2]->GetXaxis()->SetTitle("Phi [deg]");
		// 		h1_PhiDist[2]->GetYaxis()->SetTitle("Number of Events");
		// 		h1_PhiDist[2]->SetTitle("Isolated H&V Combination Events Phi");
		// 	h1_PhiDist[3]->Draw("same");
		// 		h1_PhiDist[3]->SetLineWidth(3);
		// 		h1_PhiDist[3]->SetLineColor(kRed);
		// sprintf(thistitle, "/users/PAS0654/osu0673/A23_analysis_new2/AraRoot/analysis/a23_analysis_programs/proposal2019/A%d_SurfaceSpatialDistribution_IsolatedEvents_ExcludingBursts.png",station);
		// cSpatial_Isol->SaveAs(thistitle);
		delete cSpatial_Isol;

		gStyle->SetOptStat(0);
		h1_ThetaDist[0]->SetStats(false);


		int firebrick = TColor::GetColor("#B22222");
		int royalblue = TColor::GetColor("#4169E1");

		h1_ThetaDist[0]->SetLineColor(royalblue);
		h1_ThetaDist[1]->SetLineColor(firebrick);
		h1_PhiDist[0]->SetLineColor(royalblue);
		h1_PhiDist[1]->SetLineColor(firebrick);


		TCanvas *cPlot = new TCanvas("","",2*500,2*500);
		cPlot->Divide(1,2);
		cPlot->cd(1);
			h1_ThetaDist[0]->Draw("");
			h1_ThetaDist[0]->GetXaxis()->SetRangeUser(20,40);
			h1_ThetaDist[0]->GetYaxis()->SetRangeUser(0.5,2e2);
			h1_ThetaDist[0]->SetTitle("");
			gStyle->SetOptStat(0);
			h1_ThetaDist[0]->SetLineWidth(4);
			h1_ThetaDist[0]->GetXaxis()->SetTitle("Zenith [deg]");

			h1_ThetaDist[1]->Draw("same");
			h1_ThetaDist[1]->SetLineColor(kRed);
			h1_ThetaDist[1]->SetLineWidth(4);

			// TLine *TIR_line = new TLine(36,0.5,36,2e2);
			// TIR_line->Draw("same");
			// TIR_line->SetLineColor(kBlack);
			// TIR_line->SetLineStyle(9);
			// TIR_line->SetLineWidth(3);

			TArrow *TIR_arrow = new TArrow(36.,2e2,36.,8e1 ,0.03,">");

			TIR_arrow->Draw("");
			TIR_arrow->SetLineColor(kBlack);
			TIR_arrow->SetLineWidth(3);

			TLatex TIR_label(36.5,8e1,"TIR");
			TIR_label.Draw();
   			TIR_label.SetTextSize(0.06);
			
			gPad->SetLogy();
			gPad->SetRightMargin(0.05);
			gPad->SetTopMargin(0.05);
			gPad->SetBottomMargin(0.15);
			gPad->SetLeftMargin(0.1);

			h1_ThetaDist[0]->GetXaxis()->SetTitleSize(0.08);
			h1_ThetaDist[0]->GetXaxis()->SetLabelSize(0.07);
			h1_ThetaDist[0]->GetYaxis()->SetTitleSize(0.08);
			h1_ThetaDist[0]->GetYaxis()->SetLabelSize(0.07);
			h1_ThetaDist[0]->GetYaxis()->SetTitleOffset(0.65);


			{
				TLegend *leg = new TLegend(0.15,0.6,0.4,0.85);
				leg->AddEntry(h1_ThetaDist[0],"VPol","l");
				leg->AddEntry(h1_ThetaDist[1],"HPol","l");
				leg->Draw();
				leg->SetBorderSize(0);  //no border for legend                                                                                                                                                                                         
				leg->SetFillColor(0);  //fill color is white      
			}

			gPad->Modified();
			gPad->Update();
		cPlot->cd(2);
			h1_PhiDist[0]->Draw("");
			h1_PhiDist[0]->SetTitle("");
			h1_PhiDist[0]->GetYaxis()->SetRangeUser(0.5,3e1);
			gStyle->SetOptStat(0);
			h1_PhiDist[0]->SetLineWidth(4);
			h1_PhiDist[0]->GetXaxis()->SetTitle("Azimuth [deg]");

			h1_PhiDist[1]->Draw("same");
			h1_PhiDist[1]->SetLineColor(kRed);
			h1_PhiDist[1]->SetLineWidth(4);

			// TLine *SP_line = new TLine(-106.,0.5,-106.,3e1);
			// SP_line->Draw("same");
			// SP_line->SetLineColor(kBlack);
			// SP_line->SetLineStyle(9);
			// SP_line->SetLineWidth(3);

			TArrow *SP_arrow = new TArrow(-106.,30,-106,15 ,0.03,">");

			SP_arrow->Draw("");
			SP_arrow->SetLineColor(kBlack);
			SP_arrow->SetLineWidth(3);

			TLatex SP_label(-95,15,"South Pole");
			SP_label.Draw();
   			SP_label.SetTextSize(0.06);

			gPad->SetLogy();
			gPad->SetRightMargin(0.05);
			gPad->SetTopMargin(0.05);
			gPad->SetBottomMargin(0.15);
			gPad->SetLeftMargin(0.1);

			h1_PhiDist[0]->GetXaxis()->SetTitleSize(0.08);
			h1_PhiDist[0]->GetXaxis()->SetLabelSize(0.07);
			h1_PhiDist[0]->GetYaxis()->SetTitleSize(0.08);
			h1_PhiDist[0]->GetYaxis()->SetLabelSize(0.07);
			h1_PhiDist[0]->GetYaxis()->SetTitleOffset(0.65);

			gPad->Modified();
			gPad->Update();

		// sprintf(thistitle, "/users/PAS0654/osu0673/A23_analysis_new2/AraRoot/analysis/a23_analysis_programs/proposal2019/a%d_isol_sample_dist.png",station);
		cPlot->SaveAs(thistitle);
		delete cPlot;	
	}
}