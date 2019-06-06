////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	v2_analysis_reco.cxx 
////	A23 diffuse, do reconstruction
////
////	Nov 2018
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <iomanip>
#include <sstream>

//AraRoot Includes
#include "RawAtriStationEvent.h"
#include "UsefulAraStationEvent.h"
#include "UsefulAtriStationEvent.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TH1D.h"

RawAtriStationEvent *rawAtriEvPtr;
UsefulAtriStationEvent *realAtriEvPtr;

#include "Settings.h"
#include "Event.h"
#include "Detector.h"
#include "Report.h"

#include "AraAntennaInfo.h"
#include "AraQualCuts.h"
#include "RayTraceCorrelator.h"

#include "tools_inputParameters.h"
#include "tools_outputObjects.h"
#include "tools_runSummaryObjects.h"
#include "tools_WaveformFns.h"
#include "tools_PlottingFns.h"
#include "tools_Constants.h"
#include "tools_RecoFns.h"
#include "tools_Cuts.h"

AraAntPol::AraAntPol_t Vpol = AraAntPol::kVertical;
AraAntPol::AraAntPol_t Hpol = AraAntPol::kHorizontal;

int main(int argc, char **argv)
{
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;
	char *DataDir(getenv("DATA_DIR"));
	if (DataDir == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	
	if(argc<3) {
		std::cout << "Usage\n" << argv[0] << " <station> <year> <SNRlow> <SNRhigh> <input reco file> \n";
		return -1;
	}
	int station=atoi(argv[1]);
	int year=atoi(argv[2]);
	double SNRlow = double(atof(argv[3]));
	double SNRhigh = double(atof(argv[4]));

	TH2D *h2RecoDirections[2];
	TH1D *h1SNR[2];
	for(int pol=0; pol<2; pol++){
		h2RecoDirections[pol] = new TH2D("","",360,-180,180,180,-90,90);
			h2RecoDirections[pol]->GetYaxis()->SetTitle("Theta (deg)");
			h2RecoDirections[pol]->GetXaxis()->SetTitle("Phi (deg)");
		h1SNR[pol] = new TH1D("","",100,0,10);
			h1SNR[pol]->GetXaxis()->SetTitle("SNR");
			h1SNR[pol]->GetYaxis()->SetTitle("Number of Events");
			h1SNR[pol]->GetYaxis()->SetTitleOffset(1.2);
	}

	vector<int> BadRunList=BuildBadRunList(station);

	for(int file_num=5; file_num<argc; file_num++){
		string chRun = "run";
		string file = string(argv[file_num]);
		size_t foundRun = file.find(chRun);
		string strRunNum = file.substr(foundRun+4,4);
		int runNum = atoi(strRunNum.c_str());
		int isThisBadABadRun = isBadRun(station,runNum,BadRunList);

		TFile *inputFile = TFile::Open(argv[file_num]);
		if(!inputFile){
			cout<<"Can't open joined file!"<<endl;
			return -1;
		}
		TTree *OutputTree;
		OutputTree = (TTree*) inputFile->Get("OutputTree");
		vector<TH2D*> maps;
		maps.resize(2);
		int unixTime;
		int isShort;
		OutputTree->SetBranchAddress("VMap",&maps[0]);
		OutputTree->SetBranchAddress("HMap",&maps[1]);
		OutputTree->SetBranchAddress("hasDigitizerError",&hasDigitizerError);
		OutputTree->SetBranchAddress("isCalpulser", &isCalpulser);
		OutputTree->SetBranchAddress("isSoftTrigger", &isSoftTrigger);
		OutputTree->SetBranchAddress("unixTime", &unixTime);
		OutputTree->SetBranchAddress("isShort",&isShort);

		Long64_t numEntries=OutputTree->GetEntries();
		Long64_t starEvery=numEntries/80;
		if(starEvery==0) starEvery++;

		char filter_file_name[400];
		sprintf(filter_file_name,"%s/Joined/A%d/%d/processed_station_%d_run_%d_joined_bins_6_19.root",DataDir,station,year,station,runNum);
		// TFile *filterFile = TFile::Open(filter_file_name);
		TTree *filterTree;
		// TTree *filterTree = (TTree*) filterFile->Get("OutputTree_filter");
		// filterTree->SetBranchAddress("thirdVPeakOverRMS", &thirdVPeakOverRMS);

		bool doStackingAnalysis=true;
		if(doStackingAnalysis){
			bool found_good_start=false;
			int good_start_index;
			TH2D *stacked_maps[2];
			int numStacked=0;
			// numEntries=1000;
			for(int event=0; event<numEntries; event++){
				OutputTree->GetEntry(event);
				if(!hasDigitizerError && !isCalpulser && !isSoftTrigger && !isShort){
					found_good_start=true;
					stacked_maps[0]=(TH2D*) maps[0]->Clone();
					stacked_maps[1]=(TH2D*) maps[1]->Clone();
					numStacked++;
					good_start_index=event;
				}
				if(found_good_start)
					break;
			}
			for(int event=good_start_index+1; event<numEntries; event++){
				OutputTree->GetEntry(event);
				if(!hasDigitizerError && !isCalpulser && !isSoftTrigger && !isShort){
					if(event%starEvery==0)
						std:cout<<"*"<<"Stacking Event "<<event<<endl;
					stacked_maps[0]->Add(maps[0]);
					stacked_maps[1]->Add(maps[1]);
					numStacked++;
				}
			}
			stacked_maps[0]->Scale(1./double(numStacked));
			stacked_maps[1]->Scale(1./double(numStacked));

			for(int pol=0; pol<2; pol++){
				int maxTheta, maxPhi, maxCorrBin;
				stacked_maps[pol]->GetMaximumBin(maxPhi, maxTheta, maxCorrBin);
				double maxCorr = stacked_maps[pol]->GetMaximum();
				maxTheta-=90;
				maxPhi-=180;
				printf("Pol %d Max Corr is %.5f at Theta, Phi = %d, %d \n", pol, maxCorr, maxTheta, maxPhi);

				int minTheta, minPhi, minCorrBin;
				stacked_maps[pol]->GetMinimumBin(minPhi, minTheta, minCorrBin);
				double minCorr = stacked_maps[pol]->GetMinimum();
				minTheta-=90;
				minPhi-=180;	
				printf("Pol %d 	Min Corr is %.5f at Theta, Phi = %d, %d \n", pol, minCorr, minTheta, minPhi);

				printf("Pol %d Ratio of Max/Min is %.5f/%.5f = %.5f \n", pol, maxCorr, minCorr, maxCorr/minCorr);
				printf("-------------------------------------------------\n");
			}

			double plot_max=-100.;
			double plot_min = 100.;
			for(int pol=0; pol<2; pol++){
				double this_max = stacked_maps[pol]->GetMaximum();
				double this_min = stacked_maps[pol]->GetMinimum();
				if(this_max>plot_max) plot_max=this_max;
				if(this_min<plot_min) plot_min=this_min;
			}

			beautify_TH2D();
			gStyle->SetOptStat(0);
			TCanvas *c = new TCanvas("","",2.1*850,850);
			c->Divide(2,1);
			c->cd(1);
				stacked_maps[0]->Draw("colz");
				stacked_maps[0]->GetXaxis()->SetTitle("Phi (deg)");
				stacked_maps[0]->GetYaxis()->SetTitle("Theta (deg)");
				stacked_maps[0]->GetZaxis()->SetRangeUser(plot_min, plot_max);
				gPad->SetRightMargin(0.15);
			c->cd(2);
				stacked_maps[1]->Draw("colz");
				stacked_maps[1]->GetXaxis()->SetTitle("Phi (deg)");
				stacked_maps[1]->GetYaxis()->SetTitle("Theta (deg)");
				stacked_maps[1]->GetZaxis()->SetRangeUser(plot_min, plot_max);
				gPad->SetRightMargin(0.15);
			char title[300];
			sprintf(title, "%s/maps_stacking/%d.%d.%d_A%d_Run%d_%dStackedMap.png",plotPath,year_now, month_now, day_now,station,runNum,numStacked);
			c->SaveAs(title); 
			delete c;
			cout<<"Done making an averaged map..."<<endl;

			bool doIndividualSubtraction=true;
			if(doIndividualSubtraction){

				TH2D *before_distro[2];
				TH2D *after_distro[2];
				for(int pol=0; pol<2; pol++){
					before_distro[pol] = new TH2D("","",360,-180,180,180,-90,90);
					after_distro[pol] = new TH2D("","",360,-180,180,180,-90,90);
				}

				for(int event=good_start_index+1; event<numEntries; event++){
					OutputTree->GetEntry(event);
					// if(!hasDigitizerError && !isSoftTrigger && !isShort){
					// if(!hasDigitizerError && !isCalpulser && !isSoftTrigger && !isShort){
					cout<<"Doing indiv event "<<event<<endl;
					if(isCalpulser){
						TH2D *subtracted[2];
						int before_theta[2];
						int before_phi[2];
						int after_theta[2];
						int after_phi[2];
						double before_corr[2];
						double after_corr[2];
						for(int pol=0; pol<2; pol++){
							
							//get original map
							subtracted[pol] = (TH2D*) maps[pol]->Clone();
							
							//compute stuff about the original
							getCorrMapPeak(subtracted[pol], before_theta[pol], before_phi[pol], before_corr[pol]);
							// printf("Event %d before location is %d and %d \n", event, before_theta[pol], before_phi[pol]);
							before_distro[pol]->Fill(before_phi[pol], before_theta[pol]);

							//background subtract
							subtracted[pol]->Add(stacked_maps[pol],-1.);
							getCorrMapPeak(subtracted[pol], after_theta[pol], after_phi[pol], after_corr[pol]);
							// printf("Event %d after location is %d and %d \n", event, after_theta[pol], after_phi[pol]);
							after_distro[pol]->Fill(after_phi[pol], after_theta[pol]);
						}

						bool doThisPrint=true;
						if(doThisPrint){
							stringstream V_original;
								V_original<<" VPol Original Peak (#phi, #theta)=("<<before_phi[0]<<" , "<<before_theta[0]<<")";
								maps[0]->SetTitle(V_original.str().c_str());
							stringstream V_subtracted;
								V_subtracted<<" VPol Subtracted Peak (#phi, #theta)=("<<after_phi[0]<<" , "<<after_theta[0]<<")";
								subtracted[0]->SetTitle(V_subtracted.str().c_str());
							stringstream H_original;
								H_original<<" HPol Original Peak (#phi, #theta)=("<<before_phi[1]<<" , "<<before_theta[1]<<")";
								maps[1]->SetTitle(H_original.str().c_str());
							stringstream H_subtracted;
								H_subtracted<<" HPol Subtracted Peak (#phi, #theta)=("<<after_phi[1]<<" , "<<after_theta[1]<<")";
								subtracted[1]->SetTitle(H_subtracted.str().c_str());

							TCanvas *c2del = new TCanvas("","",3*850,2*850);
							c2del->Divide(3,2);
							// vpol
								c2del->cd(1);
									maps[0]->Draw("colz");
									maps[0]->GetXaxis()->SetTitle("Phi (deg)");
									maps[0]->GetYaxis()->SetTitle("Theta (deg)");
								c2del->cd(2);
									stacked_maps[0]->Draw("colz");
									stacked_maps[0]->SetTitle("VPol Averaged Stacked Map");
									stacked_maps[0]->GetXaxis()->SetTitle("Phi (deg)");
									stacked_maps[0]->GetYaxis()->SetTitle("Theta (deg)");
								c2del->cd(3);
									subtracted[0]->Draw("colz");
									subtracted[0]->GetXaxis()->SetTitle("Phi (deg)");
									subtracted[0]->GetYaxis()->SetTitle("Theta (deg)");
							// hpol
								c2del->cd(4);
									maps[1]->Draw("colz");
									maps[1]->GetXaxis()->SetTitle("Phi (deg)");
									maps[1]->GetYaxis()->SetTitle("Theta (deg)");
								c2del->cd(5);
									stacked_maps[1]->Draw("colz");
									stacked_maps[1]->SetTitle("HPol Averaged Stacked Map");
									stacked_maps[1]->GetXaxis()->SetTitle("Phi (deg)");
									stacked_maps[1]->GetYaxis()->SetTitle("Theta (deg)");
								c2del->cd(6);
									subtracted[1]->Draw("colz");
									subtracted[1]->GetXaxis()->SetTitle("Phi (deg)");
									subtracted[1]->GetYaxis()->SetTitle("Theta (deg)");
							char this_title[400];
							sprintf(this_title, "%s/maps_stacking/%d.%d.%d_Run%d_Event%d_OriginalAndSubtracted.png",plotPath,year_now, month_now, day_now,runNum,event);
							c2del->SaveAs(this_title); 
							delete c2del;
						}
						delete subtracted[0]; delete subtracted[1];
					}
				}


				TCanvas *c2 = new TCanvas("","",2.1*850,2.1*850);
				c2->Divide(2,2);
				c2->cd(1);
					before_distro[0]->SetTitle("VPol Event Distribution, before subtraction");
					before_distro[0]->Draw("colz");
					before_distro[0]->GetXaxis()->SetTitle("Phi (deg)");
					before_distro[0]->GetYaxis()->SetTitle("Theta (deg)");
					before_distro[0]->GetZaxis()->SetRangeUser(0,1);
				c->cd(2);
					after_distro[0]->SetTitle("VPol Event Distribution, after subtraction");
					// after_distro[0]->Draw("colz");
					after_distro[0]->GetXaxis()->SetTitle("Phi (deg)");
					after_distro[0]->GetYaxis()->SetTitle("Theta (deg)");
					after_distro[0]->GetZaxis()->SetRangeUser(0,1);
				c2->cd(3);
					before_distro[1]->SetTitle("HPol Event Distribution, before subtraction");
					before_distro[1]->Draw("colz");
					before_distro[1]->GetXaxis()->SetTitle("Phi (deg)");
					before_distro[1]->GetYaxis()->SetTitle("Theta (deg)");
					before_distro[1]->GetZaxis()->SetRangeUser(0,1);
				c->cd(4);
					after_distro[1]->SetTitle("HPol Event Distribution, after subtraction");
					// after_distro[1]->Draw("colz");
					after_distro[1]->GetXaxis()->SetTitle("Phi (deg)");
					after_distro[1]->GetYaxis()->SetTitle("Theta (deg)");
					after_distro[1]->GetZaxis()->SetRangeUser(0,1);
				sprintf(title, "%s/maps_stacking/%d.%d.%d_A%d_Run%d_%dEvents_BeforeAndAfterSubtraction.png",plotPath,year_now, month_now, day_now,station,runNum,numStacked);
				c2->SaveAs(title); 
				delete c2;


				TH1D *projections[2];
				double max=0;
				for(int pol=0; pol<2; pol++){
					projections[pol]=(TH1D*)before_distro[pol]->ProjectionY()->Clone();
					projections[pol]->GetXaxis()->SetTitle("Theta (deg)");
					projections[pol]->GetYaxis()->SetTitle("Number of Events");
					if(projections[pol]->GetMaximum()>max) max=projections[pol]->GetMaximum();
				}
				for(int pol=0; pol<2; pol++){
					projections[pol]->GetYaxis()->SetRangeUser(0,max*1.1);
					double total = projections[pol]->Integral();
					double inrange = projections[pol]->Integral(-45+90, -39+90);
					printf("Pol %d In range over total: %.2f/%.2f = %.2f \n", pol, inrange,total, inrange/total);
				}
				TCanvas *c3 = new TCanvas("","",2.1*850,850);
				c3->Divide(2,1);
				c3->cd(1);
					projections[0]->Draw("");
					projections[0]->SetTitle("VPol Distribution, Theta Projection");
				c3->cd(2);
					projections[1]->Draw("");
					projections[1]->SetTitle("HPol Distribution, Theta Projection");
				sprintf(title,"%s/maps_stacking/%d.%d.%d_A%d_Run%d_%dEvents_BeforeProjectionX.png",plotPath,year_now, month_now, day_now,station,runNum,numStacked);
				c3->SaveAs(title);
				delete c3;
			}
			delete stacked_maps[0]; delete stacked_maps[1];
		}

		bool doAnisotropyAnalysis=false;
		if(doAnisotropyAnalysis){
			int numUsed=0;
			// numEntries=10000;
			for(int event=0; event<numEntries; event++){
				if(event%starEvery==0)
					std::cout<<"*"<<"Event "<<event<<endl;
				OutputTree->GetEntry(event);
				// filterTree->GetEntry(event);
				// if(!hasDigitizerError && !isSoftTrigger && !isShort){
				if(!hasDigitizerError && !isCalpulser && !isSoftTrigger && !isShort){
					numUsed++;
					TH2D *maps_copy[2];
					int thetas[2];
					int phis[2];
					double corrs[2];
					for(int pol=0; pol<2; pol++){
						if(thirdVPeakOverRMS[pol]>=SNRlow && thirdVPeakOverRMS[pol]<SNRhigh){
							maps_copy[pol] = (TH2D*) maps[pol]->Clone();
							getCorrMapPeak(maps_copy[pol], thetas[pol], phis[pol], corrs[pol]);
							h2RecoDirections[pol]->Fill(phis[pol], thetas[pol]);
							h1SNR[pol]->Fill(thirdVPeakOverRMS[pol]);
							delete maps_copy[pol];
						}
					}
				}
			}

			gStyle->SetOptStat(0);
			TCanvas *c2 = new TCanvas("","",2.1*850,850);
			c2->Divide(2,1);
			c2->cd(1);
				h2RecoDirections[0]->SetTitle("VPol Event Distribution");
				h2RecoDirections[0]->Draw("colz");
				h2RecoDirections[0]->GetXaxis()->SetTitle("Phi (deg)");
				h2RecoDirections[0]->GetYaxis()->SetTitle("Theta (deg)");
				h2RecoDirections[0]->GetZaxis()->SetRangeUser(0,1);
			c2->cd(2);
				h2RecoDirections[1]->SetTitle("HPol Event Distribution");
				h2RecoDirections[1]->Draw("colz");
				h2RecoDirections[1]->GetXaxis()->SetTitle("Phi (deg)");
				h2RecoDirections[1]->GetYaxis()->SetTitle("Theta (deg)");
				h2RecoDirections[1]->GetZaxis()->SetRangeUser(0,1);
			char title[300];
			sprintf(title, "%s/maps_stacking/%d.%d.%d_A%d_Run%d_New_%dEvents_DistroRecoDirections_SNR%d_to_%d.png",plotPath,year_now, month_now, day_now,station,runNum,numUsed,int(SNRlow),int(SNRhigh));
			c2->SaveAs(title); 
			delete c2;

			TCanvas *cSNR = new TCanvas("","",2.1*850,850);
			cSNR->Divide(2,1);
			cSNR->cd(1);
				h1SNR[0]->Draw("");
				h1SNR[0]->SetTitle("VPol SNR Distribution");
			cSNR->cd(2);
				h1SNR[1]->Draw("");
				h1SNR[1]->SetTitle("HPol SNR Distribution");
			sprintf(title,"%s/maps_stacking/%d.%d.%d_A%d_Run%d_New_%dEvents_DistroSNR_SNR%d_to_%d.png",plotPath,year_now, month_now, day_now,station,runNum,numUsed,int(SNRlow),int(SNRhigh));
			cSNR->SaveAs(title);


			TH1D *projections[2];
			double max=0;
			for(int pol=0; pol<2; pol++){
				projections[pol]=(TH1D*)h2RecoDirections[pol]->ProjectionY()->Clone();
				projections[pol]->GetXaxis()->SetTitle("Theta (deg)");
				projections[pol]->GetYaxis()->SetTitle("Number of Events");
				if(projections[pol]->GetMaximum()>max) max=projections[pol]->GetMaximum();
			}
			for(int pol=0; pol<2; pol++){
				projections[pol]->GetYaxis()->SetRangeUser(0,max*1.1);
				double total = projections[pol]->Integral();
				double inrange = projections[pol]->Integral(-45+90, -39+90);
				printf("Pol %d In range over total: %.2f/%.2f = %.2f \n", pol, inrange,total, inrange/total);
			}
			TCanvas *c3 = new TCanvas("","",2.1*850,850);
			c3->Divide(2,1);
			c3->cd(1);
				projections[0]->Draw("");
				projections[0]->SetTitle("VPol Distribution, Theta Projection");
			c3->cd(2);
				projections[1]->Draw("");
				projections[1]->SetTitle("HPol Distribution, Theta Projection");
			sprintf(title,"%s/maps_stacking/%d.%d.%d_A%d_Run%d_New_%dEvents_Projection_SNR%d_to_%d.png",plotPath,year_now, month_now, day_now,station,runNum,numUsed,int(SNRlow),int(SNRhigh));
			c3->SaveAs(title);
			delete c3;
		}
		inputFile->Close();
	}
}