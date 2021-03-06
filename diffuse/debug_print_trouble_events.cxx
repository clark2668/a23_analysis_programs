////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	v2_final_plots.cxx 
////	A23 diffuse, make plots of the final cut parameter space
////
////	Nov 2018
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
#include "TLine.h"

//AraRoot includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "Settings.h"
#include "Detector.h"
#include "Report.h"
#include "RayTraceCorrelator.h"
#include "AraQualCuts.h"
AraAntPol::AraAntPol_t Vpol = AraAntPol::kVertical;
AraAntPol::AraAntPol_t Hpol = AraAntPol::kHorizontal;
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_Cuts.h"
#include "tools_CW.h"

using namespace std;

int PlotThisEvent(int station, int config, int runNum, int event, int problempol, Settings *settings, Detector *detector, RayTraceCorrelator *theCorrelators[2]);
bool ReconsiderThisEventForGlitch(int station,int runNum,int event,Settings *settings, Detector *detector, RayTraceCorrelator *theCorrelators[2]);

int main(int argc, char **argv)
{

	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <station> <config> <ValForCuts filename>"<<endl;;
		return -1;
	}
	int station = atoi(argv[1]);
	int config = atoi(argv[2]);

	if(station!=2 && station!=3){
		printf("No good! You asked for station %d, but this code only works for stations 2 and 3 \n",station);
		return -1;
	}
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

	stringstream ss;

	gStyle->SetOptStat(11);

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
	
	// double max=1.;
	double max=0.05;

	TH2D *PeakCorr_vs_SNR_all[2];
	PeakCorr_vs_SNR_all[0]=new TH2D("","V",30,0,30,100,0,max);
	PeakCorr_vs_SNR_all[1]=new TH2D("","H",30,0,30,100,0,max);

	TH2D *PeakCorr_vs_SNR_cutCal[2];
	PeakCorr_vs_SNR_cutCal[0]=new TH2D("","V",30,0,30,100,0,max);
	PeakCorr_vs_SNR_cutCal[1]=new TH2D("","H",30,0,30,100,0,max);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft[2];
	PeakCorr_vs_SNR_cutCal_cutSoft[0]=new TH2D("","V",30,0,30,100,0,max);
	PeakCorr_vs_SNR_cutCal_cutSoft[1]=new TH2D("","H",30,0,30,100,0,max);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[2];
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[0]=new TH2D("","V",30,0,30,100,0,max);
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[1]=new TH2D("","H",30,0,30,100,0,max);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[2];
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[0]=new TH2D("","V",30,0,30,100,0,max);
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[1]=new TH2D("","H",30,0,30,100,0,max);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[2];
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[0]=new TH2D("","V",30,0,30,100,0,max);
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[1]=new TH2D("","H",30,0,30,100,0,max);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[2];
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[0]=new TH2D("","V",30,0,30,100,0,max);
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[1]=new TH2D("","H",30,0,30,100,0,max);

	TH2D *special[2];
	special[0]=new TH2D("","V",90,0,30,500,0,1);
	special[1]=new TH2D("","H",90,0,30,500,0,1);

	TH1D *fracs_power_cut[2];
	fracs_power_cut[0]=new TH1D("","V",100,0,1);
	fracs_power_cut[1]=new TH1D("","H",100,0,1);

	TH2D *spatial_distro_remaining[4];
	spatial_distro_remaining[0]=new TH2D("","V41",360,-180,180,180,-90,90);
	spatial_distro_remaining[1]=new TH2D("","H41",360,-180,180,180,-90,90);
	spatial_distro_remaining[2]=new TH2D("","V300",360,-180,180,180,-90,90);
	spatial_distro_remaining[3]=new TH2D("","H300",360,-180,180,180,-90,90);

	TH1D *surface_distro[2];
	surface_distro[0]=new TH1D("","SurfaceV",91,0,91);
	surface_distro[1]=new TH1D("","SurfaceH",91,0,91);
	TH1D *surface_distro_good[2];
	surface_distro_good[0]=new TH1D("","SurfaceVGood",91,0,91);
	surface_distro_good[1]=new TH1D("","SurfaceHGood",91,0,91);

	TH1D *h_num_surface_events = new TH1D("","",500,0,500);
	
	int num_total=0;
	int num_in_final_plot=0;
	int num_refilt=0;

	vector<int> BadRunList=BuildBadRunList(station);

	for(int file_num=3; file_num<argc; file_num++){

		string chRun = "run";
		string file = string(argv[file_num]);
		size_t foundRun = file.find(chRun);
		string strRunNum = file.substr(foundRun+4,4);
		int runNum = atoi(strRunNum.c_str());
		int isThisBadABadRun = isBadRun(station,runNum,BadRunList);

		if(isThisBadABadRun)
			continue;

		TFile *inputFile = TFile::Open(argv[file_num]);
		if(!inputFile){
			cout<<"Can't open joined file "<<argv[file_num]<<endl;
			return -1;
		}
		printf("File %d: run %d \n", file_num, runNum);

		TTree *trees[3];
		trees[0] = (TTree*) inputFile->Get("VTree");
		trees[1] = (TTree*) inputFile->Get("HTree");
		trees[2] = (TTree*) inputFile->Get("AllTree");

		double corr_val[2];
		double snr_val[2];
		int WFRMS[2];
		double frac_of_power_notched_V[8];
		double frac_of_power_notched_H[8];
		int Refilt[2];
		int theta_300[2];
		int phi_300[2];
		int theta_41[2];
		int phi_41[2];

		trees[0]->SetBranchAddress("corr_val_V",&corr_val[0]);
		trees[0]->SetBranchAddress("snr_val_V",&snr_val[0]);
		trees[0]->SetBranchAddress("wfrms_val_V",&WFRMS[0]);
		trees[0]->SetBranchAddress("Refilt_V",&Refilt[0]);
		trees[0]->SetBranchAddress("theta_300_V",&theta_300[0]);
		trees[0]->SetBranchAddress("theta_41_V",&theta_41[0]);
		trees[0]->SetBranchAddress("phi_300_V",&phi_300[0]);
		trees[0]->SetBranchAddress("phi_41_V",&phi_41[0]);
		
		trees[1]->SetBranchAddress("corr_val_H",&corr_val[1]);
		trees[1]->SetBranchAddress("snr_val_H",&snr_val[1]);
		trees[1]->SetBranchAddress("wfrms_val_H",&WFRMS[1]);
		trees[1]->SetBranchAddress("Refilt_H",&Refilt[1]);
		trees[1]->SetBranchAddress("theta_300_H",&theta_300[1]);
		trees[1]->SetBranchAddress("theta_41_H",&theta_41[1]);
		trees[1]->SetBranchAddress("phi_300_H",&phi_300[1]);
		trees[1]->SetBranchAddress("phi_41_H",&phi_41[1]);

		int isCal;
		int isSoft;
		int isShort;
		int isCW;
		int isNewBox;
		int isSurf[2];
		int isBadEvent;
		int isSurfEvent_top[2];
		int unixTime;
		int isFirstFiveEvent;
		int hasBadSpareChanIssue;

		trees[2]->SetBranchAddress("cal",&isCal);
		trees[2]->SetBranchAddress("soft",&isSoft);
		trees[2]->SetBranchAddress("short",&isShort);
		trees[2]->SetBranchAddress("CW",&isCW);
		trees[2]->SetBranchAddress("box",&isNewBox);
		trees[2]->SetBranchAddress("surf_V",&isSurf[0]);
		trees[2]->SetBranchAddress("surf_H",&isSurf[1]);
		trees[2]->SetBranchAddress("bad",&isBadEvent);
		trees[2]->SetBranchAddress("surf_top_V",&isSurfEvent_top[0]);
		trees[2]->SetBranchAddress("surf_top_H",&isSurfEvent_top[1]);
		trees[2]->SetBranchAddress("unixTime",&unixTime);
		trees[2]->SetBranchAddress("isFirstFiveEvent",&isFirstFiveEvent);
		trees[2]->SetBranchAddress("hasBadSpareChanIssue",&hasBadSpareChanIssue);

		stringstream ss;
		for(int i=0; i<8; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			trees[0]->SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_V[i]);
		}
		for(int i=8; i<16; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			trees[1]->SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_H[i-8]);
		}
		
		int numEntries = trees[0]->GetEntries();
		int num_surf_this_run=0;

		//now to loop over events
		for(int event=0; event<numEntries; event++){

			trees[0]->GetEvent(event);
			trees[1]->GetEvent(event);
			trees[2]->GetEvent(event);

			num_total++;

			if(isBadEvent|| isFirstFiveEvent || hasBadSpareChanIssue){
				continue;
			}
			if(isBadLivetime(station,unixTime)){
				continue;
			}

			// if(!isCal && !isSoft && !isShort && !isNewBox){
				
			// 	//count number of surface events
			// 	if(isSurf[0] || isSurf[1])
			// 		num_surf_this_run++;

			// 	for(int pol=0; pol<2; pol++){
			// 		if(WFRMS[pol])
			// 			continue;
			// 		surface_distro[pol]->Fill(theta_300[pol]);
			// 		if(!isThisBadABadRun)
			// 			surface_distro_good[pol]->Fill(theta_300[pol]);
			// 	}	
			// }

			for(int pol=0; pol<2; pol++){

				PeakCorr_vs_SNR_all[pol]->Fill(snr_val[pol],corr_val[pol]);
				
				if(!isCal){ //cut cal pulsers
					PeakCorr_vs_SNR_cutCal[pol]->Fill(snr_val[pol],corr_val[pol]);
					
					if(!isSoft){ //cut software triggers 
						PeakCorr_vs_SNR_cutCal_cutSoft[pol]->Fill(snr_val[pol],corr_val[pol]);
						
						if(!isShort){ //cut short
							PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->Fill(snr_val[pol],corr_val[pol]);
							
							if(!WFRMS[pol]){ //cut WRMS
								PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->Fill(snr_val[pol],corr_val[pol]);
								
								if(!isNewBox){ //cut cal box
									PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->Fill(snr_val[pol],corr_val[pol]);


									if((!isSurf[0] && !isSurf[1])  && !isSurfEvent_top[pol]){
									// if((!isSurf[0] && !isSurf[1])  && !isSurfEvent_top[0] && !isSurfEvent_top[1]){

										bool condition = false;
										// if(snr_val[pol]>=8.) condition=true;
										// if(corr_val[pol]>=0.15) condition=true;
										// if(snr_val[pol]>=7 && corr_val[pol]>=0.12) condition=true;
										// if(corr_val[pol]<0.003 || snr_val[pol]>=7.) condition=true;
										// if(corr_val[pol]>0.01) condition=true;
										// if(corr_val[pol]>0.01) condition=true;
										// if(snr_val[pol]>=7.) condition=true;
										// if(snr_val[pol]>=8.) condition=true;

										if(Refilt[pol]){
											num_refilt++;

											vector<double> frac;
											for(int i=0; i<8; i++){
												if(pol==0) frac.push_back(frac_of_power_notched_V[i]);
												else if(pol==1) frac.push_back(frac_of_power_notched_H[i]);
											}
											sort(frac.begin(), frac.end(), std::greater<double>());
											fracs_power_cut[pol]->Fill(frac[2]);
											if(frac[2]<=0.06){
												// if(!condition){
												// 	// printf("Run %d, Event %d: Corr %.2f and SNR %.2f \n", runNum, event, corr_val[pol], snr_val[pol]);
												// 	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->Fill(snr_val[pol],corr_val[pol]);
												// }
												// else{											
												// 	printf("Run %d, Event %d: Corr %.2f and SNR %.2f and V theta, phi = %d, %d  and H theta, phi = %d, %d \n", runNum, event, corr_val[pol], snr_val[pol],theta_300[0], phi_300[0],theta_300[1], phi_300[1]);
												// 	printf("	Run %d, Event %d Surface status in pol %d is %d and in pol %d is %d \n",runNum, event, 0,isSurf[0], 1, isSurf[1]);
												// 	// printf("Reconsiering for Glitch Run %d Event %d in Pol %d\n",runNum,event,pol);
												// 	bool failsReconsideration=false;
												// 	failsReconsideration=ReconsiderThisEventForGlitch(station,runNum,event,settings, detector, theCorrelators);
												// 	if(!failsReconsideration){
												// 		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->Fill(snr_val[pol],corr_val[pol]);
												// 		PlotThisEvent(station,config,runNum,event, pol, settings, detector, theCorrelators);
												// 	}
												// }
												PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->Fill(snr_val[pol],corr_val[pol]);
												if(condition){

													// printf("		Event %d is refiltered in pol %d \n", event,pol);
													// cout<<"			Frac of power notched is "<<frac[2]<<endl;
													// printf("Event has condition in pol %d \n", pol);
													// printf("Event %d Unixtime is %d \n", event, unixTime);
													
													// printf("if(runNum==%d && pol==%d && unixTime==%d && event==%d){\n\tunixTimes[%d].push_back(double(unixTime)); phis[%d].push_back(double(phi_300[%d])); thetas[%d].push_back(double(theta_300[%d]));\n}\n",runNum,pol,unixTime,event,pol,pol,pol,pol,pol);

													spatial_distro_remaining[pol]->Fill(phi_41[pol],theta_41[pol]);
													spatial_distro_remaining[pol+2]->Fill(phi_300[pol], theta_300[pol]);
													// printf("if(runNum==%d && event==%d) TroubleEvent=true;\n",runNum,event);
													// printf("Run %d, Event %d: Corr %.2f and SNR %.2f \n", runNum, event, corr_val[pol], snr_val[pol]);
													PlotThisEvent(station,config,runNum,event, pol, settings, detector, theCorrelators);
												}
											}
										} //refiltered?
										else{
											// if(!condition){
											// 	// printf("Run %d, Event %d: Corr %.2f and SNR %.2f \n", runNum, event, corr_val[pol], snr_val[pol]);
											// 	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->Fill(snr_val[pol],corr_val[pol]);
											// }
											// else{
											// 	printf("Run %d, Event %d: Corr 4%.2f and SNR %.2f \n", runNum, event, corr_val[pol], snr_val[pol]);
											// 	// printf("Reconsiering for Glitch Run %d Event %d in Pol %d\n",runNum,event,pol);
											// 	bool failsReconsideration=false;
											// 	failsReconsideration = ReconsiderThisEventForGlitch(station,runNum,event,settings, detector, theCorrelators);
											// 	if(!failsReconsideration){
											// 		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->Fill(snr_val[pol],corr_val[pol]);
											// 		PlotThisEvent(station,config,runNum,event, pol, settings, detector, theCorrelators);
											// 	}
											// }
											PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->Fill(snr_val[pol],corr_val[pol]);
											if(condition){
												// printf("		Event %d is NOT refiltered in pol %d \n", event,pol);
												spatial_distro_remaining[pol]->Fill(phi_41[pol],theta_41[pol]);
												spatial_distro_remaining[pol+2]->Fill(phi_300[pol], theta_300[pol]);
												// printf("Event %d Unixtime is %d \n", event, unixTime);
												// printf("if(runNum==%d && pol==%d && unixTime==%d && event==%d){\n\tunixTimes[%d].push_back(double(unixTime)); phis[%d].push_back(double(phi_300[%d])); thetas[%d].push_back(double(theta_300[%d]));\n}\n",runNum,pol,unixTime,event,pol,pol,pol,pol,pol);
												// printf("if(runNum==%d && event==%d) TroubleEvent=true;\n",runNum,event);
												// printf("Run %d, Event %d: Corr %.2f and SNR %.2f \n", runNum, event, corr_val[pol], snr_val[pol]);
												PlotThisEvent(station,config,runNum,event, pol, settings, detector, theCorrelators);
											}
										}
										num_in_final_plot++;
									}
								}
							}
						}
					}
				}
			}
		}
		h_num_surface_events->Fill(num_surf_this_run);
		// if(num_surf_this_run>500)
		
		// char title_txt[200];
		// sprintf(title_txt,"%s/surf_channels_per_run_config%d.txt",plotPath,config);
		// FILE *fout = fopen(title_txt, "a");
		// fprintf(fout,"%d, %d \n",runNum,num_surf_this_run);
		// fclose(fout);

		inputFile->Close();
		delete inputFile;
	}

	char title[300];


	bool print_surface_stuff=false;
	if(print_surface_stuff){

		TLegend *leg = new TLegend(0.52,0.7,0.75,0.9);
		leg->AddEntry(surface_distro[0],"All Events","l");
		leg->AddEntry(surface_distro_good[0],"Events in 'Good' Runs","l");
		leg->SetTextSize(0.02);

		TLine l(37,0.1,37,surface_distro[0]->GetMaximum()*1.2);
		l.SetLineStyle(9);

		TCanvas *c_spatial_distro = new TCanvas("","",2.1*850,2.1*850);
		c_spatial_distro->Divide(2,2);
		for(int pol=0; pol<2; pol++){
			c_spatial_distro->cd(pol+1);
				spatial_distro_remaining[pol]->Draw("colz");
				spatial_distro_remaining[pol]->GetZaxis()->SetRangeUser(1,3);
			c_spatial_distro->cd(pol+3);
				spatial_distro_remaining[pol+2]->Draw("colz");
				spatial_distro_remaining[pol+2]->GetZaxis()->SetRangeUser(1,3);
		}
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_SpatialDistroRemainingEvents.png",plotPath,year_now, month_now, day_now,station,config,num_total);
		c_spatial_distro->SaveAs(title);
		delete c_spatial_distro;

		TH1D *project_theta[2];
		project_theta[0]=(TH1D*) spatial_distro_remaining[2]->ProjectionY("")->Clone();
		project_theta[1]=(TH1D*) spatial_distro_remaining[3]->ProjectionY("")->Clone();
		TH1D *project_phi[2];
		project_phi[0]=(TH1D*) spatial_distro_remaining[2]->ProjectionX("")->Clone();
		project_phi[1]=(TH1D*) spatial_distro_remaining[3]->ProjectionX("")->Clone();

		TCanvas *c_spatial_distro_project = new TCanvas("","",2.1*850,850);
		c_spatial_distro_project->Divide(2,2);
		for(int pol=0; pol<2; pol++){
			c_spatial_distro_project->cd(pol+1);
				project_theta[pol]->Draw("");
				project_theta[pol]->GetXaxis()->SetTitle("Theta");
				project_theta[pol]->GetYaxis()->SetTitle("Number of Events");
			c_spatial_distro_project->cd(pol+3);
				project_phi[pol]->Draw("");
				project_phi[pol]->GetXaxis()->SetTitle("Phi");
				project_phi[pol]->GetYaxis()->SetTitle("Number of Events");
			if(pol==0){
				project_theta[pol]->SetTitle("VPol Theta Distribution of Stragglers");
				project_phi[pol]->SetTitle("VPol Phi Distribution of Stragglers");
			}
			else{
				project_theta[pol]->SetTitle("HPol Theta Distribution of Stragglers");
				project_phi[pol]->SetTitle("Hpol Phi Distribution of Stragglers");
			}

		}
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_ThetaPhiDistroRemainingEvents.png",plotPath,year_now, month_now, day_now,station,config,num_total);
		c_spatial_distro_project->SaveAs(title);
		delete c_spatial_distro_project;
		delete project_theta[0]; delete project_theta[1];
		delete spatial_distro_remaining[0]; delete spatial_distro_remaining[1]; delete spatial_distro_remaining[2]; delete spatial_distro_remaining[3];

		gStyle->SetOptStat(111111);
		TCanvas *c_num_surface_per_run = new TCanvas("","",850,850);
			h_num_surface_events->Draw();
			h_num_surface_events->GetYaxis()->SetTitle("Number of Runs");
			h_num_surface_events->GetXaxis()->SetTitle("Number of Surface Events per Run");
			h_num_surface_events->GetYaxis()->SetTitleOffset(1.3);
			// gPad->SetLogx();
			gPad->SetLogy();
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_NumSurfEventsPerRun.png",plotPath,year_now, month_now, day_now,station,config,num_total);
		c_num_surface_per_run->SaveAs(title);
		delete c_num_surface_per_run;
		delete h_num_surface_events;

		// char equation_phi[150];
		// sprintf(equation_phi,"expo");
		// TF1 *fit_phi = new TF1("GausFit_surface",equation_phi,25.,37.);
		// surface_distro[0]->Fit("GausFit_surface","R");
		// printf("Chi-Square/NDF %.2f / %.2f \n",fit_phi->GetChisquare(),double(fit_phi->GetNDF()));

		TCanvas *c_surface_event_distro = new TCanvas("","",2*850,850);
		c_surface_event_distro->Divide(2,1);
		c_surface_event_distro->cd(1);
			surface_distro[0]->Draw();
			surface_distro[0]->GetXaxis()->SetTitle("Theta (deg)");
			surface_distro[0]->GetYaxis()->SetTitle("Number of Events");
			surface_distro[0]->GetYaxis()->SetTitleOffset(1.3);
			surface_distro_good[0]->Draw("same");
			surface_distro_good[0]->SetLineColor(kRed);
			gPad->SetLogy();
			leg->Draw();
			l.Draw("same");
		c_surface_event_distro->cd(2);
			surface_distro[1]->Draw();
			surface_distro[1]->GetXaxis()->SetTitle("Theta (deg)");
			surface_distro[1]->GetYaxis()->SetTitle("Number of Events");
			surface_distro[1]->GetYaxis()->SetTitleOffset(1.3);
			surface_distro_good[1]->Draw("same");
			surface_distro_good[1]->SetLineColor(kRed);
			gPad->SetLogy();
			l.Draw("same");
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_SurfaceEventDistro.png",plotPath,year_now, month_now, day_now,station,config,num_total);
		c_surface_event_distro->SaveAs(title);
		delete c_surface_event_distro;
		delete surface_distro[0]; delete surface_distro[1];
	}

	bool print_summary=true;
	if(print_summary){

		gStyle->SetOptStat(11);


		cout<<"Num total is "<<num_total<<endl;
		cout<<"Num in final plot "<<num_in_final_plot<<endl;
		cout<<"Num re-filtered is "<<num_refilt<<endl;

		gStyle->SetOptStat(11);
		gStyle->SetStatY(0.9);
		gStyle->SetStatX(0.9);
		gStyle->SetStatW(0.2);
		gStyle->SetStatH(0.2);

		//save out SNR vs WavefrontRMS plot
		char graph_title[2][300];
		char title[300];

		int cal=0;
		int soft=0;
		int Short=0;
		int wrms=0;
		int box = 0;
		int surf=0;
		int cw=0;

		//save out the Corr vs SNR plot for all 
		sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
		sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
		TCanvas *c2 = new TCanvas("","",2.1*850,850);
		c2->Divide(2,1);
		for(int pol=0; pol<2; pol++){
			c2->cd(pol+1);
			PeakCorr_vs_SNR_all[pol]->Draw("colz");
			PeakCorr_vs_SNR_all[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
			PeakCorr_vs_SNR_all[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
			PeakCorr_vs_SNR_all[pol]->SetTitle(graph_title[pol]);
			gPad->SetLogz();
		}
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,num_total,cal,soft,Short,wrms,box,surf);
		c2->SaveAs(title);
		delete c2;
		delete PeakCorr_vs_SNR_all[0]; delete PeakCorr_vs_SNR_all[1];

		//turn on cal
		cal=1;
		sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
		sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
		TCanvas *c3 = new TCanvas("","",2.1*850,850);
		c3->Divide(2,1);
		for(int pol=0; pol<2; pol++){
			c3->cd(pol+1);
			PeakCorr_vs_SNR_cutCal[pol]->Draw("colz");
			PeakCorr_vs_SNR_cutCal[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
			PeakCorr_vs_SNR_cutCal[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
			PeakCorr_vs_SNR_cutCal[pol]->SetTitle(graph_title[pol]);
			gPad->SetLogz();
		}
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,num_total,cal,soft,Short,wrms,box,surf);
		c3->SaveAs(title);
		delete c3;
		delete PeakCorr_vs_SNR_cutCal[0]; delete PeakCorr_vs_SNR_cutCal[1];

		//turn on cal, soft
		soft=1;
		sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
		sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
		TCanvas *c4 = new TCanvas("","",2.1*850,850);
		c4->Divide(2,1);
		for(int pol=0; pol<2; pol++){
			c4->cd(pol+1);
			PeakCorr_vs_SNR_cutCal_cutSoft[pol]->Draw("colz");
			PeakCorr_vs_SNR_cutCal_cutSoft[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
			PeakCorr_vs_SNR_cutCal_cutSoft[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
			PeakCorr_vs_SNR_cutCal_cutSoft[pol]->SetTitle(graph_title[pol]);
			gPad->SetLogz();
		}
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,num_total,cal,soft,Short,wrms,box,surf);
		c4->SaveAs(title);
		delete c4;
		delete PeakCorr_vs_SNR_cutCal_cutSoft[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft[1];

		//turn on cal, soft, short
		Short=1;
		sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
		sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
		TCanvas *c5 = new TCanvas("","",2.1*850,850);
		c5->Divide(2,1);
		for(int pol=0; pol<2; pol++){
			c5->cd(pol+1);
			PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->Draw("colz");
			PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
			PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
			PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->SetTitle(graph_title[pol]);
			gPad->SetLogz();
		}
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,num_total,cal,soft,Short,wrms,box,surf);
		c5->SaveAs(title);
		delete c5;
		delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[1];

		//turn on cal, soft, short, wmrs
		wrms=1;
		sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
		sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
		TCanvas *c6 = new TCanvas("","",2.1*850,850);
		c6->Divide(2,1);
		for(int pol=0; pol<2; pol++){
			c6->cd(pol+1);
			PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->Draw("colz");
			PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
			PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
			PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->SetTitle(graph_title[pol]);
			gPad->SetLogz();
		}
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,num_total,cal,soft,Short,wrms,box,surf);
		c6->SaveAs(title);
		delete c6;
		delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[1];

		//turn on cal, soft, short, wmrs, box
		box=1;
		sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
		sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
		TCanvas *c7 = new TCanvas("","",2.1*850,850);
		c7->Divide(2,1);
		for(int pol=0; pol<2; pol++){
			c7->cd(pol+1);
			PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->Draw("colz");
			PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
			PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
			PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->SetTitle(graph_title[pol]);
			gPad->SetLogz();
		}
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,num_total,cal,soft,Short,wrms,box,surf);
		c7->SaveAs(title);
		delete c7;
		delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[1];

		gStyle->SetOptStat(0);
		surf=1;
		sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
		sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
		TCanvas *c8 = new TCanvas("","",2.1*850,850);
		c8->Divide(2,1);
		for(int pol=0; pol<2; pol++){
			c8->cd(pol+1);
			PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->Draw("colz");
			PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
			PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
			PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->SetTitle(graph_title[pol]);
			gPad->SetLogz();
			// PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->GetXaxis()->SetRangeUser(0,10);
			// PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->GetYaxis()->SetRangeUser(0,0.5);
		}
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,num_total,cal,soft,Short,wrms,box,surf);
		c8->SaveAs(title);
		delete c8;
		delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[1];
	}	

}

int PlotThisEvent(int station, int config, int runNum, int event, int problempol, Settings *settings, Detector *detector, RayTraceCorrelator *theCorrelators[2]){
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *DataDirPath(getenv("DATA_DIR"));
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
	}

	double threshCW = 1.0;
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
			){
				isCutonCW_back[pol] = true;
			}
		}
	}
	bool skipCW=false;
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
			sprintf(save_temp_title,"%s/trouble_events/%d.%d.%d_Run%d_Ev%d_ProblemPol%d_WaveformsNotch.png",plotPath,year_now,month_now,day_now,runNum,event,problempol);
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

			sprintf(save_temp_title,"%s/trouble_events/%d.%d.%d_Run%d_Ev%d_ProblemPol%d_SpectraNotch.png",plotPath,year_now,month_now,day_now,runNum,event,problempol);
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
			sprintf(save_temp_title,"%s/trouble_events/%d.%d.%d_Run%d_Ev%d_ProblemPol%d_FilteredMaps.png",plotPath,year_now,month_now,day_now,runNum,event,problempol);
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

	// printf("Trying to reconstruct theta\n");
	// AraGeomTool *araGeom = AraGeomTool::Instance();
	// double dist=abs(araGeom->getStationInfo(station)->getAntennaInfo(3)->antLocation[2] - araGeom->getStationInfo(station)->getAntennaInfo(7)->antLocation[2]);
	// TGraph *corr = FFTtools::getInterpolatedCorrelationGraph(waveforms[7],waveforms[3],0.5);
	// int peak_bin=FFTtools::getPeakBin(corr);
	// double delay=corr->GetX()[peak_bin];
	// cout<<"Delay is "<<delay<<endl;
	// delay*=1e-9;
	// double theta = TMath::ASin(3e8*delay/1.78/dist)*TMath::RadToDeg();
	// printf("Delay is %.2f, and reco theta is %.2f \n", delay,theta);
	// delete corr;

	bool do_reco=true;
	if(do_reco){
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

		vector<double> chan_SNRs;

		map_30m_V = theCorrelators[0]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Vpol, false,chan_list_V,chan_SNRs,0,-1) ;
		map_300m_V = theCorrelators[1]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Vpol, false,chan_list_V,chan_SNRs,0,-1);
		map_30m_H = theCorrelators[0]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Hpol, false,chan_list_H,chan_SNRs,0,-1);
		map_300m_H = theCorrelators[1]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Hpol, false,chan_list_H,chan_SNRs,0,-1);

		// map_30m_V = theCorrelators[0]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Vpol, false,chan_list_V) ;
		// map_300m_V = theCorrelators[1]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Vpol, false,chan_list_V);
		// map_30m_H = theCorrelators[0]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Hpol, false,chan_list_H);
		// map_300m_H = theCorrelators[1]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Hpol, false,chan_list_H);

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
		char save_temp_title[400];		
		sprintf(save_temp_title,"%s/trouble_events/%d.%d.%d_Run%d_Ev%d_ProblemPol%d_Maps.png",plotPath,year_now,month_now,day_now,runNum,event,problempol);
		cMaps->SaveAs(save_temp_title);
		delete cMaps;
		delete map_30m_V; delete map_300m_V; delete map_30m_H; delete map_300m_H;

		chan_list_V.clear();
		chan_list_V.push_back(0);
		chan_list_V.push_back(1);
		chan_list_V.push_back(2);
		if(!(station==3 && runNum>getA3BadRunBoundary())){ //if dropping bad chans and station 3, don't keep fourth string
			chan_list_V.push_back(3);
		}

		chan_list_H.clear();
		chan_list_H.push_back(8);
		chan_list_H.push_back(9);
		chan_list_H.push_back(10);
		if(!(station==3 && runNum>getA3BadRunBoundary())){ //if dropping bad chans and station 3, don't keep fourth string
			chan_list_H.push_back(11);
		}

		TH2D *map_300m_top_V = theCorrelators[1]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Vpol, false, chan_list_V);
		TH2D *map_300m_top_H = theCorrelators[1]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Hpol, false, chan_list_H);

		int PeakTheta_Recompute_300m_top_V;
		int PeakPhi_Recompute_300m_top_V;
		double PeakCorr_Recompute_300m_top_V;
		int PeakTheta_Recompute_300m_top_H;
		int PeakPhi_Recompute_300m_top_H;
		double PeakCorr_Recompute_300m_top_H;
		getCorrMapPeak(map_300m_top_V,PeakTheta_Recompute_300m_top_V,PeakPhi_Recompute_300m_top_V,PeakCorr_Recompute_300m_top_V);
		getCorrMapPeak(map_300m_top_H,PeakTheta_Recompute_300m_top_H,PeakPhi_Recompute_300m_top_H,PeakCorr_Recompute_300m_top_H);

		stringstream ss300H_top;
		ss300H_top<<" Peak Theta, Phi is "<<PeakTheta_Recompute_300m_top_H<<" , "<<PeakPhi_Recompute_300m_top_H;
		map_300m_top_H->SetTitle(ss300H_top.str().c_str());

		stringstream ss300V_top;
		ss300V_top<<" Peak Theta, Phi is "<<PeakTheta_Recompute_300m_top_V<<" , "<<PeakPhi_Recompute_300m_top_V;
		map_300m_top_V->SetTitle(ss300V_top.str().c_str());

		TCanvas *cMaps_top = new TCanvas("","",2*850,850);
		cMaps_top->Divide(2,1);
			cMaps_top->cd(1);
			map_300m_top_V->Draw("colz");
			cMaps_top->cd(2);
			map_300m_top_H->Draw("colz");
		sprintf(save_temp_title,"%s/trouble_events/%d.%d.%d_Run%d_Ev%d_ProblemPol%d_MapsTop.png",plotPath,year_now,month_now,day_now,runNum,event,problempol);
		cMaps_top->SaveAs(save_temp_title);
		delete cMaps_top;
		delete map_300m_top_V; delete map_300m_top_H;

		chan_list_V.clear();
		chan_list_H.clear();
		chan_list_V.push_back(4);
		chan_list_V.push_back(5);
		chan_list_V.push_back(6);
		chan_list_V.push_back(7);

		chan_list_H.push_back(12);
		chan_list_H.push_back(13);
		chan_list_H.push_back(14);
		chan_list_H.push_back(15);

		TH2D *map_300m_bottom_V = theCorrelators[1]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Vpol, false, chan_list_V);
		TH2D *map_300m_bottom_H = theCorrelators[1]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Hpol, false, chan_list_H);

		int PeakTheta_Recompute_300m_bottom_V;
		int PeakPhi_Recompute_300m_bottom_V;
		double PeakCorr_Recompute_300m_bottom_V;
		int PeakTheta_Recompute_300m_bottom_H;
		int PeakPhi_Recompute_300m_bottom_H;
		double PeakCorr_Recompute_300m_bottom_H;
		getCorrMapPeak(map_300m_bottom_V,PeakTheta_Recompute_300m_bottom_V,PeakPhi_Recompute_300m_bottom_V,PeakCorr_Recompute_300m_bottom_V);
		getCorrMapPeak(map_300m_bottom_H,PeakTheta_Recompute_300m_bottom_H,PeakPhi_Recompute_300m_bottom_H,PeakCorr_Recompute_300m_bottom_H);

		stringstream ss300H_bottom;
		ss300H_bottom<<" Peak Theta, Phi is "<<PeakTheta_Recompute_300m_bottom_H<<" , "<<PeakPhi_Recompute_300m_bottom_H;
		map_300m_bottom_H->SetTitle(ss300H_bottom.str().c_str());

		stringstream ss300V_bottom;
		ss300V_bottom<<" Peak Theta, Phi is "<<PeakTheta_Recompute_300m_bottom_V<<" , "<<PeakPhi_Recompute_300m_bottom_V;
		map_300m_bottom_V->SetTitle(ss300V_bottom.str().c_str());

		TCanvas *cMaps_bottom = new TCanvas("","",2*850,850);
		cMaps_bottom->Divide(2,1);
			cMaps_bottom->cd(1);
			map_300m_bottom_V->Draw("colz");
			cMaps_bottom->cd(2);
			map_300m_bottom_H->Draw("colz");
		sprintf(save_temp_title,"%s/trouble_events/%d.%d.%d_Run%d_Ev%d_ProblemPol%d_MapsBottom.png",plotPath,year_now,month_now,day_now,runNum,event,problempol);
		cMaps_bottom->SaveAs(save_temp_title);
		delete cMaps_bottom;
		delete map_300m_bottom_V; delete map_300m_bottom_H;
	}

	// vector<TGraph*> dummy;
	// for(int i=0; i<16; i++){
	// 	vector<double> thisX;
	// 	vector<double> thisY;
	// 	thisX.push_back(300);
	// 	thisX.push_back(600);
	// 	thisY.push_back(-700);
	// 	thisY.push_back(700);
	// 	dummy.push_back(new TGraph(2,&thisX[0], &thisY[0]));
	// }


	char save_temp_title[300];
	sprintf(save_temp_title,"%s/trouble_events/%d.%d.%d_Run%d_Ev%d_ProblemPol%d_Waveforms.png",plotPath,year_now,month_now,day_now,runNum,event,problempol);
	TCanvas *cWave = new TCanvas("","",4*1100,4*850);
	cWave->Divide(4,4);
	for(int i=0; i<16; i++){
		cWave->cd(i+1);
		// dummy[i]->Draw("AL");
		// dummy[i]->SetLineColor(kWhite);
		// dummy[i]->GetXaxis()->SetRangeUser(300.,500.);

		waveforms[i]->Draw("AL");
		waveforms[i]->SetLineWidth(3);
		// waveforms[i]->GetXaxis()->SetRangeUser(300.,500.);
	}
	cWave->SaveAs(save_temp_title);
	delete cWave;

	sprintf(save_temp_title,"%s/trouble_events/%d.%d.%d_Run%d_Ev%d_ProblemPol%d_Spectra.png",plotPath,year_now,month_now,day_now,runNum,event,problempol);
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


double MaxMeanBlock(TGraph *grIn, int evt_num, int chan, vector<double> &means, bool print){

	int n_input = grIn->GetN();
	double *oldX = grIn->GetX();
	double *oldY = grIn->GetY();

	deque <double> inX;
	deque <double> inY;	

	for(int samp=0; samp<n_input; samp++){
		inX.push_back(oldX[samp]);
		inY.push_back(oldY[samp]);
	}

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
		means.push_back((mean_this_section));
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
		for(int i=0; i<grPieces.size(); i++){
			grPieces[i]->Draw("same");
			grPieces[i]->SetLineColor(colors[i]);
			grPieces[i]->SetLineWidth(3);
		}
		char title[400];
		sprintf(title,"/users/PAS0654/osu0673/A23_analysis_new2/results/glitch_detect/GlitchWaveform_ev%d_chan%d.png",evt_num,chan);
		c->SaveAs(title);	
		delete c;
	}
	for(int i=0; i<grPieces.size(); i++){
		delete grPieces[i];
	}
	// for(int i=0; i<means.size(); i++){
	// 	printf("Block %d has mean %.2f \n", i, means[i]);
	// }
	std::sort(means.begin(), means.end(), std::greater<double>());
	// printf("Biggest value is %.2f \n", means[0]);
	return means[0];
}

int countExcursions(vector<double> maxMeans){
	if(maxMeans.size()!=16){
		cout<<"Something wrong! More than 16 elements!"<<endl;
		return 0;
	}
	int numExcursions=0;
	for(int i=0; i<16; i++){
		if((i+1)%4==0) continue;
		if(maxMeans[i]>150)
			numExcursions++;
	}
	return numExcursions;
}

bool ReconsiderThisEventForGlitch(int station,int runNum, int event,Settings *settings, Detector *detector, RayTraceCorrelator *theCorrelators[2]){
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *DataDirPath(getenv("DATA_DIR"));
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
	}
	TTree *eventTree = (TTree*) mapFile-> Get("eventTree");
	if(!eventTree){
		cout<<"Can't find eventTree for map"<<endl;
	}

	RawAtriStationEvent *rawPtr =0;
	eventTree->SetBranchAddress("event",&rawPtr);
	eventTree->GetEvent(event);
	int realEventNumber = rawPtr->eventNumber;
	bool hasRealEventNumberIssue = false;
	if(realEventNumber<=4)
		hasRealEventNumberIssue=true;

	int stationID = rawPtr->stationId;
	char ped_file_name[400];
	sprintf(ped_file_name,"%s/run_specific_peds/A%d/all_peds/event%d_specificPeds.dat",PedDirPath,station,runNum);
	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	calibrator->setAtriPedFile(ped_file_name,stationID); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist

	AraQualCuts *qualCut = AraQualCuts::Instance();
	UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawPtr,AraCalType::kLatestCalib);
	bool isGoodEventNew = qualCut->isGoodEvent(realAtriEvPtr);
	printf("Run %d, Event %d \n", runNum, realEventNumber);

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


	/*
		Rolling means excursion check
	*/

	// vector<TGraph*> rollingMeans;
	// vector<double> maxMeans;
	// rollingMeans.resize(16);
	// for(int i=0; i<16; i++){
	// 	if(waveforms[i]->GetN()<64){
	//  		continue;
	//  	}
	//  	TGraph *grInt = FFTtools::getInterpolatedGraph(waveforms[i],0.5);
	//  	if(grInt->GetN()<101){
	//  		delete grInt;
	//  		continue;
	//  	}
	//  	rollingMeans[i]=(qualCut->getRollingMean(grInt,100));
	//  	double thisMax = abs(TMath::MaxElement(rollingMeans[i]->GetN(), rollingMeans[i]->GetY()));
	//  	double thisMin = abs(TMath::MinElement(rollingMeans[i]->GetN(), rollingMeans[i]->GetY()));
	//  	double thisAbs;
	//  	if(abs(thisMax)>abs(thisMin))
	//  		thisAbs=abs(thisMax);
	//  	else
	//  		thisAbs=abs(thisMin);
	//  	maxMeans.push_back(thisAbs);
	//  	// printf("Chan %d max abs is %4.2f\n",i,thisAbs);
	// }
	// for(int i=0; i<16; i++) delete rollingMeans[i];
	// int numExcursions = countExcursions(maxMeans);
	// printf("Num excursions is %d \n", numExcursions);

	bool cutOnExcursions=false;
	// if(numExcursions>2)
	// 	cutOnExcursions=true;

	/*
		Number of non-zero blocks cut
	*/

	// vector< vector<double> > vvmeans;
	// vector<int> vNumViolatingBlocksPerChannel;
	// for(int i=0; i<16; i++){
	// 	vector<double> vmeans;
	// 	double this_max = MaxMeanBlock(waveforms[i], event, i, vmeans, false);
	// 	vvmeans.push_back(vmeans);
	// }
	// for(int i=0; i<16; i++){
	// 	int numViolatingBlocks=0;
	// 	for(int j=0; j<vvmeans[i].size();j++){
	// 		if(abs(vvmeans[i][j])>20)
	// 			numViolatingBlocks++;
	// 	}
	// 	vNumViolatingBlocksPerChannel.push_back(numViolatingBlocks);
	// 	// printf("Chan %2d has %2d num violating blocks \n", i, numViolating);
	// }
	// int numChansWithManyViolations=0;
	// for(int i=0; i<16; i++){
	// 	if(vNumViolatingBlocksPerChannel[i]>10)
	// 		numChansWithManyViolations++;
	// }
	bool cutOnManyOffsetBlocks=false;
	// if(numChansWithManyViolations>=4)
	// 	cutOnManyOffsetBlocks=true;
	// printf("Num chans with many block violations is %d \n", numChansWithManyViolations);


	/*
		Spare channels cut
	*/

	vector<TGraph*> electChansGraphs;
	electChansGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(6));
	electChansGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(14));
	electChansGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(22));
	electChansGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(30));
	vector<double> spareRMS;
	for(int i=0; i<4; i++){
		spareRMS.push_back(electChansGraphs[i]->GetRMS(2));
	}
	for(int i=0; i<4; i++)
		delete electChansGraphs[i];
	int numBadSpareChans=0;
	int numReallyBadSpareChans=0;
	for(int i=0; i<4; i++){
		printf("Eventd %d Spare RMS string %d is %3.2f\n",event,i,spareRMS[i]);
		if(spareRMS[i]>20 && i!=3){
			numBadSpareChans++;
		}
		if(spareRMS[i]>60 && i!=3){
			numReallyBadSpareChans++;
		}
	}
	bool hasBadSpareChansIssue=false;
	if(numBadSpareChans>1 || numReallyBadSpareChans>0)
		hasBadSpareChansIssue=true;


	for(int i=0; i<16; i++){
		delete waveforms[i];
	}


	/*
		SNR Weighted Reconstruction check for rejecting any surface events we'll get rid of later
	*/
	bool isSNRWeightedSurfaceEvent=false;
	// only do these reconstructions if removing the event in every other way has failed...
	// if(!cutOnManyOffsetBlocks 
	// 	&& !hasBadSpareChansIssue 
	// 	&& !hasRealEventNumberIssue){

	// 	vector <int> chan_list_V;
	// 	vector <int> chan_list_H;
	// 	vector<double> chan_SNRs;
	// 	for(int chan=0; chan<=7; chan++){
	// 		chan_list_V.push_back(chan);
	// 		chan_list_H.push_back(chan+8);
	// 	}
	// 	for(int i=0; i<16; i++){
	// 		// chan_SNRs.push_back(1.);
	// 		// chan_SNRs.push_back(VPeakOverRMS[i]);
	// 		// printf("SNR is %.2f \n", VPeakOverRMS[i]);
	// 	}
	// 	// filterFile->Close();

	// 	if(station==2){
	// 		//for station 2, we need to exclude channel 15 from the analysis
	// 		chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
	// 	}
	// 	else if(station==3){
	// 		//for station 3 years 2014 and 2015, we need to drop string 4 (channels 3, 7, 11, 15) altogether
	// 		if(runNum>getA3BadRunBoundary()){
	// 			chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 3), chan_list_V.end());
	// 			chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 7), chan_list_V.end());
	// 			chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 11), chan_list_H.end());
	// 			chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
	// 		}
	// 	}
	// 	TH2D *map_300m_V;
	// 	TH2D *map_300m_H;
	// 	bool AraSim=false;
	// 	int solNum=0;
	// 	map_300m_V = theCorrelators[1]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Vpol, AraSim, chan_list_V, chan_SNRs,  solNum);
	// 	map_300m_H = theCorrelators[1]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Hpol, AraSim, chan_list_H, chan_SNRs, solNum);
	// 	int PeakTheta_Recompute_300m_H;
	// 	int PeakPhi_Recompute_300m_H;
	// 	double PeakCorr_Recompute_300m_H;
	// 	int PeakTheta_Recompute_300m_V;
	// 	int PeakPhi_Recompute_300m_V;
	// 	double PeakCorr_Recompute_300m_V;
	// 	getCorrMapPeak(map_300m_H,PeakTheta_Recompute_300m_H,PeakPhi_Recompute_300m_H,PeakCorr_Recompute_300m_H);
	// 	getCorrMapPeak(map_300m_V,PeakTheta_Recompute_300m_V,PeakPhi_Recompute_300m_V,PeakCorr_Recompute_300m_V);

	// 	if(PeakTheta_Recompute_300m_V >=37 || PeakTheta_Recompute_300m_H>=37){
	// 		isSNRWeightedSurfaceEvent=true;
	// 		cout<<"Hey! This one reconstructs to the surface!"<<endl;
	// 	}

	// 	delete map_300m_V;
	// 	delete map_300m_H;
	// }




	// delete realAtriEvPtr;
	mapFile->Close();
	delete mapFile;
	// if( cutOnExcursions || cutOnManyOffsetBlocks || !isGoodEventNew)
	// if(cutOnManyOffsetBlocks 
	// 	|| hasBadSpareChansIssue 
	// 	|| hasRealEventNumberIssue
	// 	|| isSNRWeightedSurfaceEvent)
	if(hasRealEventNumberIssue || hasBadSpareChansIssue){
		cout<<"Rejected by reconsideration!"<<endl;
		return true;
	}
	else
		return false;
}
