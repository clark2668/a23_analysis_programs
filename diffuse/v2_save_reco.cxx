////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	v2_analysis_save_vals.cxx 
////	A23 diffuse, save values for cuts
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

using namespace std;

int main(int argc, char **argv)
{

	stringstream ss;
	
	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <station> <year> <output_location> <joined filename 1> <joined filename 2 > ... <joined filename x>";
		return 0;
	}
	int station = atoi(argv[1]);
	int year = atoi(argv[2]);
	string output_location = argv[3];

	for(int file_num=4; file_num<argc; file_num++){

		string chRun = "run";
		string file = string(argv[file_num]);
		size_t foundRun = file.find(chRun);
		string strRunNum = file.substr(foundRun+4,4);
		int runNum = atoi(strRunNum.c_str());
		printf("Run Number %d \n", runNum);

		char outfile_name[300];
		sprintf(outfile_name,"%s/vals_for_cut_run_%d.root",output_location.c_str(),runNum);
		TFile *fpOut = new TFile(outfile_name,"recreate");
		TTree *trees = new TTree("RecoVals","RecoVals");

		double phi_41_V;
		double theta_41_V;
		double phi_41_H;
		double theta_41_H;
		double corr_val_V;
		double corr_val_H;
		int isCal;
		trees->Branch("phi_41_V",&phi_41_V);
		trees->Branch("theta_41_V",&theta_41_V);
		trees->Branch("phi_41_H",&phi_41_H);
		trees->Branch("theta_41_H",&theta_41_H);
		trees->Branch("corr_val_V",&corr_val_V);
		trees->Branch("corr_val_H",&corr_val_H);
		trees->Branch("cal",&isCal);

		cout << "Run " << file_num << " :: " << argv[file_num] << endl;
		
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
		bool isCalPulser;
		inputTree_filter->SetBranchAddress("isCalpulser",&isCalPulser);

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

		//now to loop over events
		for(int event=0; event<numEntries; event++){

			// if(event%starEvery==0) {
			// 	std::cout << "	On event "<<event<<endl;
			// }

			isCal=0;
			phi_41_V=-10000;
			theta_41_V-10000;
			phi_41_H-10000;
			theta_41_H-10000;
			corr_val_V-10000;;
			corr_val_H-10000;;


			inputTree_filter->GetEvent(event);
			if(isCalPulser) isCal=1;

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

			for(int pol=0; pol<2; pol++){
				if(pol==0){
					corr_val_V=bestCorr_pulser[0];
					phi_41_V = bestPhi_pulser[0];
					theta_41_V = bestTheta_pulser[0];
				}
				if(pol==1){
					corr_val_H=bestCorr_pulser[1];
					phi_41_H = bestPhi_pulser[1];
					theta_41_H = bestTheta_pulser[1];
				}
			}
			trees->Fill();
		}//loop over events
		inputFile->Close();
		delete inputFile;

		fpOut->Write();
		fpOut->Close();
		delete fpOut;
		printf("Done! Run Number %d", runNum);
	} //end loop over input files
}
