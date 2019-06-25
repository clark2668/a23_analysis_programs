////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	2019.06.23_FilterAngularAcceptance.cxx 
////	see what angular distribution of events passing filter is
////////////////////////////////////////////////////////////////////////////////

//C++
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"

//AraRoot includes
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_inputParameters.h"
#include "tools_outputObjects.h"
#include "tools_Cuts.h"
#include "tools_CommandLine.h"

using namespace std;

int main(int argc, char **argv)
{
	gStyle->SetOptStat(0);
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *DataDirPath(getenv("DATA_DIR"));
	if (DataDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *SimDirPath(getenv("SIM_DIR"));
	if (SimDirPath == NULL) std::cout << "Warning! $SIM_DIR is not set!" << endl;
	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

	stringstream ss;
	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	
	if(argc<12){
		cout<< "Usage\n" << argv[0] << " <isSim?> <station> <config> <year_or_energy (as float, eg 17.0 or 18.5)> <drop_bad_chan> <output_location> <V SNR bin> <H SNR bin> <V WFRMS val> <H WFRMS val> <joined filename 1> <joined filename 2 > ... <joined filename x>"<<endl;
		return 0;
	}
	int isSimulation = atoi(argv[1]);
	int station = atoi(argv[2]);
	int config = atoi(argv[3]);
	double year_or_energy = double(atof(argv[4]));
	int dropBadChans = atoi(argv[5]);
	string output_location = argv[6];

	//just to have the cut parameters up front and easy to find
	int thresholdBin_pol[]={atoi(argv[7]), atoi(argv[8])}; //bin 3 = 2.3, bin 5 = 2.5 //what is the faceRMS inclusion threshold?
	double wavefrontRMScut[]={atof(argv[9]),atof(argv[10])}; //event wavefrontRMS < this value

	TH1D *filter_acceptance[2];
	filter_acceptance[0] = new TH1D("V","V",90,-180,180);
	filter_acceptance[1] = new TH1D("H","H",90,-180,180);

	for(int file_num=11; file_num<argc; file_num++){

		string file = string(argv[file_num]);
		string wordRun = "run_";
		size_t foundRun = file.find(wordRun);
		string wordFilter = "_joined";
		size_t foundFilter = file.find(wordFilter);
		size_t diff=(foundFilter-wordRun.length())-foundRun;
		string strRunNum = file.substr(foundRun+4,diff);
		int runNum = atoi(strRunNum.c_str());		

		cout << "Run " << file_num << " :: " << argv[file_num] << endl;
		
		//first, load in the data file; this shoud be a "joined" file
		//meaning it should contain "filter" trees and "reco" trees
		TFile *inputFile = TFile::Open(argv[file_num]);
		if(!inputFile){
			cout<<"Can't open joined file!"<<endl;			
			return -1;
		}
		
		//next, we need to load the filter tree
		ss.str("");
		ss << "OutputTree_filter";
		TTree *inputTree_filter = (TTree*) inputFile->Get(ss.str().c_str());
		if(!inputTree_filter){
			cout<<"Can't open filter tree"<<endl;
			return -1;
		}
		double inweight;
		inputTree_filter->SetBranchAddress("weight", &inweight);

		int numFaces_new_V;
		int numFaces_new_H;
		if(station==2){
			numFaces_new_V = numFaces;
			numFaces_new_H = numFaces_A2_drop;
		}
		else if(station==3){
			numFaces_new_V = numFaces_A3_drop;
			numFaces_new_H = numFaces_A3_drop;
		}
		double rms_pol_thresh_face_alternate_V[thresholdSteps][numFaces_new_V];
		double rms_pol_thresh_face_alternate_H[thresholdSteps][numFaces_new_H];
		inputTree_filter->SetBranchAddress("rms_pol_thresh_face_alternate_V", &rms_pol_thresh_face_alternate_V);
		inputTree_filter->SetBranchAddress("rms_pol_thresh_face_alternate_H", &rms_pol_thresh_face_alternate_H);

		//next, load the reco tree
		TTree *inputTree_reco[35];
		double peakCorr[35][2];
		int peakTheta[35][2];
		int peakPhi[35][2];
		int recoBinSelect = 19; //300 m map
		int recoBinCalpulser = 6; //41 m map
		for(int i=0; i<35; i++){
			if(i==recoBinSelect||i==recoBinCalpulser){
				ss.str("");
				ss << "OutputTree_recoRadius_" << i;
				inputTree_reco[i] = (TTree*) inputFile->Get(ss.str().c_str());
				if(!inputTree_reco[i]) {
					std::cout << "Can't find OutputTree: " << i << "\n";
					return -1;
				}
				inputTree_reco[i]->SetBranchAddress("peakCorr_single", &peakCorr[i]);
				inputTree_reco[i]->SetBranchAddress("peakTheta_single", &peakTheta[i]);
				inputTree_reco[i]->SetBranchAddress("peakPhi_single", &peakPhi[i]);
			}
		}

		int numEntries = inputTree_filter->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;

		int start=0;
		for(int event=start; event<numEntries; event++){
			if(event%starEvery==0) {
				// std::cout<<"*";
			}
			inputTree_filter->GetEvent(event);

			bool failWavefrontRMS[2];
			failWavefrontRMS[0]=false;
			failWavefrontRMS[1]=false;

			for (int i = 0; i < 35; i++){
				if (i == recoBinSelect || i == recoBinCalpulser){
					inputTree_reco[i]->GetEntry(event);
				}
			}

			//figure out which reconstruction map (vpol or hpol) is best
			//in the present analysis, this only matters for the 300m bin
			double bestCorr[] = {0., 0., 0.};
			int bestCorrRadiusBin[3];
			int bestPol = 2;
			int bestTheta[3];
			int bestPhi[3];

			for(int pol=0; pol<2; pol++){
				for(int i=0; i<35; i++){
					if(i==recoBinSelect){
						if(peakCorr[i][pol] > bestCorr[pol]){
							bestCorr[pol]=peakCorr[i][pol];
							bestCorrRadiusBin[pol]=i;
							bestTheta[pol]=peakTheta[i][pol];
							bestPhi[pol]=peakPhi[i][pol];
						}
						if(peakCorr[i][pol] > bestCorr[2]){
							bestCorr[2]=peakCorr[i][pol];
							bestCorrRadiusBin[2]=i;
							bestTheta[2]=peakTheta[i][pol];
							bestPhi[2]=peakPhi[i][pol];
							bestPol=pol;
						}
					}//300m bin check
				}//loop over reco bins
			}//loop over polarizations

			bool isSurf[2] = {false};

			for(int pol=0; pol<2; pol++){
				// printf("Pol %d has theta %d \n",pol,bestTheta[pol]);
				if(bestTheta[pol] >= 37){
					isSurf[pol]=true;
				}
			}

			//figure out which reconstruction map (vpol or hpol) is best
			//for the 41m bin
			double bestCorr_pulser[] = {0., 0., 0.};
			int bestCorrRadiusBin_pulser[3];
			int bestPol_pulser = 2;
			int bestTheta_pulser[3];
			int bestPhi_pulser[3];

			for(int pol=0; pol<2; pol++){
				for(int i=0; i<35; i++){
					if (i == recoBinCalpulser){
						if (peakCorr[i][pol] > bestCorr_pulser[pol]){
							bestCorr_pulser[pol] = peakCorr[i][pol];
							bestCorrRadiusBin_pulser[pol] = i;
							bestTheta_pulser[pol] = peakTheta[i][pol];
							bestPhi_pulser[pol] = peakPhi[i][pol];
						}
						if (peakCorr[i][pol] > bestCorr_pulser[2]){
							bestCorr_pulser[2] = peakCorr[i][pol];
							bestCorrRadiusBin_pulser[2] = i;
							bestTheta_pulser[2] = peakTheta[i][pol];
							bestPhi_pulser[2] = peakPhi[i][pol];
							bestPol_pulser = pol;
						}
					}//cal pulser (41m) bin check
				}//loop over reco bins
			}//loop over polarizations
			
			int theta_300[2];
			int phi_300[2];
			int theta_41[2];
			int phi_41[2];

			for(int pol=0; pol<2; pol++){
				theta_300[pol]=bestTheta[pol];
				phi_300[pol]=bestPhi[pol];
				theta_41[pol]=bestTheta_pulser[pol];
				phi_41[pol]=bestPhi_pulser[pol];
			}

			// so, now we have the reco directions
			vector <double> rms_faces_V;
			vector <double> rms_faces_H;

			if(dropBadChans){
				int num_faces_for_V_loop;
				int num_faces_for_H_loop;
				if(station==2){
					rms_faces_V.resize(numFaces);
					num_faces_for_V_loop=numFaces;
					rms_faces_H.resize(numFaces_A2_drop);
					num_faces_for_H_loop=numFaces_A2_drop;
				}
				else if(station==3){
					if(
						(!isSimulation && runNum>getA3BadRunBoundary())
						||
						(isSimulation && config>2)
					){ //it's 2014+, drop string four
						rms_faces_V.resize(numFaces_A3_drop);
						num_faces_for_V_loop=numFaces_A3_drop;
						rms_faces_H.resize(numFaces_A3_drop);
						num_faces_for_H_loop=numFaces_A3_drop;
					}
					else{ //it's 2013-, keep string four
						rms_faces_V.resize(numFaces);
						num_faces_for_V_loop=numFaces;
						rms_faces_H.resize(numFaces);
						num_faces_for_H_loop=numFaces;
					}
				}
				//now we loop over the faces
				for(int i=0; i<num_faces_for_V_loop; i++){
					rms_faces_V[i] = rms_pol_thresh_face_alternate_V[thresholdBin_pol[0]][i];
				}
				for(int i=0; i<num_faces_for_H_loop; i++){
					rms_faces_H[i] = rms_pol_thresh_face_alternate_H[thresholdBin_pol[1]][i];
				}
			}
			else{
				rms_faces_V.resize(numFaces);
				rms_faces_H.resize(numFaces);
				//now, we must loop over the faces
				for(int i=0; i<numFaces; i++){
					rms_faces_V[i] = rms_pol_thresh_face_V[thresholdBin_pol[0]][i];  //this is right RMS for this polarization, threshold requirement, and face
					rms_faces_H[i] = rms_pol_thresh_face_H[thresholdBin_pol[1]][i];
				}
			}

			//now to sort them smallest to largest; lowest RMS is best
			sort(rms_faces_V.begin(), rms_faces_V.end());
			sort(rms_faces_H.begin(), rms_faces_H.end());

			double bestFaceRMS[2];
			bestFaceRMS[0]=rms_faces_V[0];
			bestFaceRMS[1]=rms_faces_H[0];

			if(log(bestFaceRMS[0])/log(10) >= wavefrontRMScut[0]){
				failWavefrontRMS[0]=true;
			}
			if(log(bestFaceRMS[1])/log(10) >= wavefrontRMScut[1]){
				failWavefrontRMS[1]=true;
			}

			if(!failWavefrontRMS[0]){
				filter_acceptance[0]->Fill(phi_300[0],inweight);
			}
			if(!failWavefrontRMS[1]){
				filter_acceptance[1]->Fill(phi_300[1],inweight);
			}

		}//loop over events
		inputFile->Close();
		delete inputFile;
	} //end loop over input files

	// double scale[2]={filter_acceptance[0]->Integral("width"), filter_acceptance[1]->Integral("width")};
	// filter_acceptance[0]->Scale(1./scale[0]);
	// filter_acceptance[1]->Scale(1./scale[1]);

	TCanvas *cAcceptance = new TCanvas("","",2*850,850);
	cAcceptance->Divide(2,1);
	cAcceptance->cd(1);
		filter_acceptance[0]->Draw("");
		filter_acceptance[0]->GetYaxis()->SetRangeUser(0.,140.);
		filter_acceptance[0]->GetYaxis()->SetTitle("Number of Events (weighted)");
		filter_acceptance[0]->GetXaxis()->SetTitle("Azimuth (deg");
		filter_acceptance[0]->GetYaxis()->SetTitleOffset(1.4);
	cAcceptance->cd(2);
		filter_acceptance[1]->Draw("");
		filter_acceptance[1]->GetYaxis()->SetRangeUser(0.,140.);
		filter_acceptance[1]->GetYaxis()->SetTitle("Number of Events (weighted)");
		filter_acceptance[1]->GetXaxis()->SetTitle("Azimuth (deg");
		filter_acceptance[1]->GetYaxis()->SetTitleOffset(1.4);
	char save_title[400];
	sprintf(save_title,"%s/filter_cut/A%d_C%d_FilterAngularAcceptance.png",plotPath,station,config);
	cAcceptance->SaveAs(save_title);

}
