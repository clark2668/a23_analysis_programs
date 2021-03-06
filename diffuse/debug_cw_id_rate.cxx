////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	debug_cw_id_rate.cxx 
////	strip down save vals code to just look at wtf is going on with CWID...
////
////	Oct 2019
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
#include "TTimeStamp.h"

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
#include "tools_inputParameters.h"
#include "tools_outputObjects.h"
#include "tools_Cuts.h"
#include "tools_CommandLine.h"
#include "tools_CW.h"

using namespace std;

int PlotThisEvent(int station, int config, int runNum, int event, int problempol);

int main(int argc, char **argv)
{
	gStyle->SetOptStat(0);
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *AuxDirPath(getenv("AUX_DIR"));
	if (AuxDirPath == NULL){
		std::cout << "Warning! $AUX_DIR is not set! You need this for CWID and RunSummary files" << endl;
		return -1;
	}
	char *DataDirPath_Project(getenv("DATA_DIR_PROJECT"));
	if (DataDirPath_Project == NULL){
		std::cout << "Warning! $DATA_DIR_PROJECT is not set!" << endl;
		return -1;
	}
	char *SimDirPath(getenv("SIM_DIR"));
	if (SimDirPath == NULL){
		std::cout << "Warning! $SIM_DIR is not set!" << endl;
		return -1;
	}
	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL){
		std::cout << "Warning! $PED_DIR is not set!" << endl;
		return -1;
	}
	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL){
		std::cout << "Warning! $PLOT_PATH is not set!" << endl;
		return -1;
	}
	char *ToolsDirPath(getenv("TOOLS_DIR"));
	if(ToolsDirPath == NULL){
		std::cout << "Warning! $TOOLS_DIR is not set! This is needed for finding where your list of runs without cal pulses is "<< endl;
		return -1;
	}


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

	TH1D *filter_freqs = new TH1D("","",1024,0,1000);

	TTimeStamp start;
	TTimeStamp stop;
	int numBins=365*4;

	start.Set(2013, 01, 01, 00, 00,0,0,true,0);
	stop.Set(2013, 12, 31, 00, 00,0,0,true,0);

	int start_bin = start.GetSec();
	int stop_bin = stop.GetSec();

	TH2D *h2_filt_freq_vs_time = new TH2D("","",365,start_bin,stop_bin,1024,0,1000);
	h2_filt_freq_vs_time->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
	h2_filt_freq_vs_time->GetXaxis()->SetTimeFormat("%y/%m/%d");
	h2_filt_freq_vs_time->GetXaxis()->SetTimeOffset(0.,"GMT");

	TH2D *h2_filt_freq_vs_run = new TH2D("","",8000,0,8000,1024,0,1000);

	for(int file_num=11; file_num<argc; file_num++){

		int num_total[2]={0,0};
		int num_baseline_filt[2]={0,0};

		string file = string(argv[file_num]);
		string wordRun = "run_";
		size_t foundRun = file.find(wordRun);
		string wordFilter = "_joined";
		size_t foundFilter = file.find(wordFilter);
		size_t diff=(foundFilter-wordRun.length())-foundRun;
		string strRunNum = file.substr(foundRun+4,diff);
		int runNum = atoi(strRunNum.c_str());

		if(!isSimulation){
			//we're almost certainly going to need the calibrator, so let's just load it now
			char ped_file_name[400];
			sprintf(ped_file_name,"%s/run_specific_peds/A%d/all_peds/event%d_specificPeds.dat",PedDirPath,station,runNum);
			calibrator->setAtriPedFile(ped_file_name,station); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist
		}

		// need "org" and "new" variables so we know what's happening before and after filtering
		// the way I'm doing this is so, so dumb. Brian, don't ever do it this way again

		double corr_val_org[2];
		double snr_val_org[2];
		int WFRMS_org[2];
		int theta_300_org[2];
		int phi_300_org[2];
		int theta_41_org[2];
		int phi_41_org[2];

		double corr_val_new[2];
		double snr_val_new[2];
		int WFRMS_new[2];
		int theta_300_new[2];
		int phi_300_new[2];
		int theta_41_new[2];
		int phi_41_new[2];

		int Refilt[2];

		double frac_of_power_notched_V[8];
		double frac_of_power_notched_H[8];

		// mark these as "out" so we can see them clearly

		int isCal_out;
		int isSoft_out;
		int isShortWave_out;
		int isCW_out;
		int isNewBox_out;
		int isSurfEvent_org_out[2]; // originally a surf event?
		int isSurfEvent_new_out[2]; // a surface event after filtering?
		int isSurfEvent_top[2]; // a top event?

		int isBadEvent_out;
		double outweight;
		int Trig_Pass_out[16];
		int unixTime_out;
		int hasBadSpareChanIssue_out;
		int hasBadSpareChanIssue2_out;
		int isFirstFiveEvent_out;
		int eventNumber_out;		

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
		bool isCalPulser_in;
		bool isSoftTrigger_in;
		int waveformLength_in[16];
		bool hasDigitizerError_in;
		double inweight;
		int Trig_Pass_in[16];
		int unixTime_in;
		int eventNumber_in;

		inputTree_filter->SetBranchAddress("VPeakOverRMS",&VPeakOverRMS); //get SNR's of all channels
		inputTree_filter->SetBranchAddress("thirdVPeakOverRMS", &thirdVPeakOverRMS);
		inputTree_filter->SetBranchAddress("rms_pol_thresh_face_V", &rms_pol_thresh_face_V);
		inputTree_filter->SetBranchAddress("rms_pol_thresh_face_H", &rms_pol_thresh_face_H);
		inputTree_filter->SetBranchAddress("weight", &inweight);
		inputTree_filter->SetBranchAddress("unixTime",&unixTime_in);
		inputTree_filter->SetBranchAddress("eventNumber",&eventNumber_in);
		if(isSimulation)
			inputTree_filter->SetBranchAddress("Trig_Pass", &Trig_Pass_in);

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
		
		inputTree_filter->SetBranchAddress("isCalpulser",&isCalPulser_in);
		inputTree_filter->SetBranchAddress("isSoftTrigger",&isSoftTrigger_in);
		inputTree_filter->SetBranchAddress("waveformLength",&waveformLength_in);
		inputTree_filter->SetBranchAddress("hasDigitizerError",&hasDigitizerError_in);

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

		char summary_file_name[400];
		if(isSimulation){
			if(year_or_energy<25)
				sprintf(summary_file_name,"%s/CWID/A%d/c%d/E%2.2f/CWID_station_%d_run_%d.root",SimDirPath,station,config,year_or_energy,station,runNum);
			else
				sprintf(summary_file_name,"%s/CWID/A%d/c%d/E%d/CWID_station_%d_run_%d.root",SimDirPath,station,config,int(year_or_energy),station,runNum);
		}
		else{
			sprintf(summary_file_name,"%s/CWID/A%d/all_runs/CWID_station_%d_run_%d.root",AuxDirPath,station,station,runNum);
		}
		TFile *NewCWFile = TFile::Open(summary_file_name);
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

		// start out with a problem threshold of 1.5 for A2 and most of A3
		// if the run has untagged cal pulses though, lift the threshold to 2.0
		double threshCW=1.5;
		bool doesRunHaveUntagedCalPulses = hasUntaggedCalpul(ToolsDirPath, station, config, runNum);
		if(station==3 && doesRunHaveUntagedCalPulses){
			threshCW=2.0;
		}

		int numEntries = inputTree_filter->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;
		cout<<"Num entries is "<<numEntries<<endl;
		cout<<"Star every is "<<starEvery<<endl;

		int start=0;
		//now to loop over events
		for(int event=start; event<numEntries; event++){
			if(event%starEvery==0) {
				// std::cout<<"*";
			}

			// reset the variables that we are passing *out*

			// things that don't require loops
			isCal_out=0;
			isSoft_out=0;
			isShortWave_out=0;
			isCW_out=0;
			isNewBox_out=0;
			isFirstFiveEvent_out=false;
			hasBadSpareChanIssue_out=false;
			hasBadSpareChanIssue2_out=false;

			for(int pol=0; pol<2; pol++){
				isSurfEvent_org_out[pol]=0;
				isSurfEvent_new_out[pol]=0;
				Refilt[pol]=0;
				corr_val_org[pol]=-10000;
				snr_val_org[pol]=-10000;
				corr_val_new[pol]=-10000;
				snr_val_new[pol]=-10000;
				WFRMS_org[pol]=0;
				WFRMS_new[pol]=0;
				for(int i=0; i<8; i++){
					frac_of_power_notched_V[i]=0.;
					frac_of_power_notched_H[i]=0.;
				}
				isSurfEvent_top[pol]=0;
			}

			inputTree_filter->GetEvent(event);

			unixTime_out=unixTime_in; //copy over the unixtime
			eventNumber_out=eventNumber_in;
			if(eventNumber_out<5 && !isSimulation){ //eep, never check this for simulation, will be huge efficiency hit!
				isFirstFiveEvent_out=true;
			}

			bool isShort=false;
			bool isSurf[2];
			isSurf[0]=false;
			isSurf[1]=false;
			// bool isSurf[2]={false};
			bool isCP5=false;
			bool isCP6=false;
			bool failWavefrontRMS[2];
			failWavefrontRMS[0]=false;
			failWavefrontRMS[1]=false;

			outweight=inweight;
			isBadEvent_out=hasDigitizerError_in;

			for(int i=0;i<16;i++){
				if(waveformLength_in[i]<500){
					isShort=true;
				}
				Trig_Pass_out[i]=Trig_Pass_in[i];
			}

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


			for(int pol=0; pol<2; pol++){
				if(bestTheta[pol] >= 37){
					// printf("Pol %d has theta %d \n",pol,bestTheta[pol]);
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
			
			//draw a box around the cal pulser
			for (int pol = 0; pol < 2; pol++){
				identifyCalPulser(station,config, bestTheta_pulser[pol], bestPhi_pulser[pol], isCP5, isCP6);
			}
			for(int pol=0; pol<2; pol++){
				theta_300_org[pol]=bestTheta[pol];
				phi_300_org[pol]=bestPhi[pol];
				theta_41_org[pol]=bestTheta_pulser[pol];
				phi_41_org[pol]=bestPhi_pulser[pol];
				
				// we also need to populate the "new" arrays
				theta_300_new[pol]=theta_300_org[pol];
				phi_300_new[pol]=phi_300_org[pol];
				theta_41_new[pol]=theta_41_org[pol];
				phi_41_new[pol]=phi_41_org[pol];
			}

			//deal w/ CW cut
			//inputTree_CW->GetEntry(event);
			NewCWTree->GetEntry(event);

			bool isCutonCW_fwd[2]; isCutonCW_fwd[0]=false; isCutonCW_fwd[1]=false;
			bool isCutonCW_back[2]; isCutonCW_back[0]=false; isCutonCW_back[1]=false;
			bool isCutonCW_baseline[2]; isCutonCW_baseline[0]=false; isCutonCW_baseline[1]=false;
			
			for(int pol=0; pol<badFreqs_baseline->size(); pol++){
				vector<double> badFreqListLocal_baseline = badFreqs_baseline->at(pol);
				if(badFreqListLocal_baseline.size()>0){
					for(int freq=0; freq<badFreqListLocal_baseline.size(); freq++){
						// if(badFreqListLocal_baseline[freq]>140.)
						if(badFreqListLocal_baseline[freq]>0.)
							isCutonCW_baseline[pol]=true;
					}
				}
			}

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
						badFreqList_back[iCW] < 850.
					){
						isCutonCW_back[pol] = true;
						// cout<<"The bad frequency mode is "<<badFreqList_back[iCW]<<endl;
					}
				}
			}

			//filter associated parameters
			double SNRs[2];
			SNRs[0] = thirdVPeakOverRMS[0];
			SNRs[1] = thirdVPeakOverRMS[1];
			vector <double> rms_faces_V;
			vector <double> rms_faces_H;

			if(dropBadChans){
				int num_faces_for_V_loop;
				int num_faces_for_H_loop;
				
				bool needToSwitchToAltWFRMSArrays=false;
				// we need an actual flag for if we need to switch to the _alternate_ arrays

				if(station==2){
					rms_faces_V.resize(numFaces);
					num_faces_for_V_loop=numFaces;
					rms_faces_H.resize(numFaces_A2_drop);
					num_faces_for_H_loop=numFaces_A2_drop;
					needToSwitchToAltWFRMSArrays=true;
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
						needToSwitchToAltWFRMSArrays=true;
					}
					else{ //it's 2013-, keep string four
						rms_faces_V.resize(numFaces);
						num_faces_for_V_loop=numFaces;
						rms_faces_H.resize(numFaces);
						num_faces_for_H_loop=numFaces;
					}
				}
				if(needToSwitchToAltWFRMSArrays){
					//now we loop over the faces in the *alternate* arrays
					for(int i=0; i<num_faces_for_V_loop; i++){
						rms_faces_V[i] = rms_pol_thresh_face_alternate_V[thresholdBin_pol[0]][i];
					}
					for(int i=0; i<num_faces_for_H_loop; i++){
						rms_faces_H[i] = rms_pol_thresh_face_alternate_H[thresholdBin_pol[1]][i];
					}
				}
				else{
					//now we loop over the faces in the not alternate arrays
					for(int i=0; i<num_faces_for_V_loop; i++){
						rms_faces_V[i] = rms_pol_thresh_face_V[thresholdBin_pol[0]][i];
					}
					for(int i=0; i<num_faces_for_H_loop; i++){
						rms_faces_H[i] = rms_pol_thresh_face_H[thresholdBin_pol[1]][i];
					}
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

			if(isCalPulser_in){
				isCal_out=1;
			}
			if(isSoftTrigger_in){
				isSoft_out=1;
			}
			if(isShort){
				isShortWave_out=1;
			}
			if(isCP5 || isCP6 ){
				isNewBox_out=1;
			}

			for(int pol=0; pol<2; pol++){
				if(failWavefrontRMS[pol]){
					WFRMS_org[pol]=1;
				}
				if(isSurf[pol]){
					isSurfEvent_org_out[pol]=1;
				}
				// and, we need to set the "new" parameters to the "old" parameters
				// so that if they don't get altered again after filtering, they are set correctly
				isSurfEvent_new_out[pol] = isSurfEvent_org_out[pol];
				WFRMS_new[pol] = WFRMS_org[pol];
			}

			for(int pol=0; pol<1; pol++){
			// for(int pol=0; pol<2; pol++){

					num_total[pol]++;
					if(isCutonCW_baseline[pol] && !isSoftTrigger_in && !isShort ){
					// if(isCutonCW_baseline[pol]){
						num_baseline_filt[pol]++;
					}
					// printf("CW Status for Run %4d Event %5d, Pol %d, Cal %d, CP5 %d, CP6 %d, Soft %d, Short %d, badSpareChan %d: CW fwd %d, CW back %d, CW Baseline %d \n", runNum, event, pol, isCalPulser_in, isCP5, isCP6, isSoftTrigger_in, isShort, hasBadSpareChanIssue_out, isCutonCW_fwd[pol], isCutonCW_back[pol], isCutonCW_baseline[pol]);

					if(isCutonCW_baseline[pol] && pol==0 && !isSoftTrigger_in && !isShort){
					// if(isCutonCW_baseline[pol] && pol==0){
						vector<double> badFreqListLocal_baseline = badFreqs_baseline->at(pol);
						for(int newFreq=0; newFreq<badFreqListLocal_baseline.size(); newFreq++){
							double new_freq = badFreqListLocal_baseline[newFreq];
							// cout<<"     unixTime is "<< unixTime_in<<" and Problem freq is "<<new_freq<<endl;
							filter_freqs->Fill(new_freq);
							h2_filt_freq_vs_time->Fill(unixTime_in,new_freq);
							h2_filt_freq_vs_run->Fill(runNum,new_freq);
						}
						// PlotThisEvent(station, config, runNum, event, pol);
					}


				corr_val_org[pol]=bestCorr[pol];
				snr_val_org[pol]=SNRs[pol];

				// copy these over
				corr_val_new[pol] = corr_val_org[pol];
				snr_val_new[pol] = snr_val_org[pol];

				// printf(BLUE"Run %4d, Event %5d/%5d, Pol %d : isCal %d, isSoft %d, isShort %d, does Fail WFRMS %d, isCP5 %d, isCP6 %d, isSurf %d and %d, isBad %d, isFirstFive %d \n"RESET
				// 	,runNum, event, eventNumber_out, pol, isCalPulser_in, isSoftTrigger, isShort, failWavefrontRMS[pol], isCP5, isCP6, isSurf[0], isSurf[1], isBadEvent_out, isFirstFiveEvent_out);

				if(1==2
				// if(!isCalPulser_in
				// 	&& !isSoftTrigger_in
				// 	&& !isShort
				// 	&& !failWavefrontRMS[pol]
				// 	&& !isCP5 && !isCP6
				// 	&& !isBadEvent_out
				// 	&& !isFirstFiveEvent_out
					// && 2==3
					// now, we will let all events which are surface be filtered too. this is a  major change.
					// && !isSurf[0] && !isSurf[1] //check both pols for surface
				){ //cut cal pulsers

					// printf(RED"Time to do math on event %d \n"RESET,event);

					// load in the data for the event
					char run_file_name[400];
					if(isSimulation)
						if(year_or_energy<25)
							sprintf(run_file_name,"%s/RawSim/A%d/c%d/E%2.1f/AraOut.A%d_c%d_E%2.1f.txt.run%d.root",SimDirPath,station,config,year_or_energy,station,config,year_or_energy,runNum);
						else{
							// should just be Kotera
							sprintf(run_file_name,"%s/RawSim/A%d/c%d/E%d/AraOut.A%d_c%d_E%d.txt.run%d.root",SimDirPath,station,config,int(year_or_energy),station,config,int(year_or_energy),runNum);
						}
					else{
						// sprintf(run_file_name,"%s/RawData/A%d/by_config/c%d/event%d.root",DataDirPath,station,config,runNum);
						// FIX ME: you should never hard code this, but I'm hacking it in...
						sprintf(run_file_name,"/fs/scratch/PAS0654/ara/10pct/RawData/A%d/by_config/c%d/event%d.root",station,config,runNum);
						// sprintf(run_file_name,"%s/RawData/A%d/by_config/c%d/event%d.root",DataDirPath_Project,station,config,runNum);
					}
					TFile *mapFile = TFile::Open(run_file_name,"READ");
					if(!mapFile){
						cout<<"Can't open data file for map!"<<endl;
						return -1;
					}
					TTree *eventTree = (TTree*) mapFile-> Get("eventTree");
					if(!eventTree){
						cout<<"Can't find eventTree for map"<<endl;
						return -1;
					}

					UsefulAtriStationEvent *realAtriEvPtr=0;
					RawAtriStationEvent *rawPtr =0;

					if(isSimulation){
						eventTree->SetBranchAddress("UsefulAtriStationEvent", &realAtriEvPtr);
						eventTree->GetEvent(event);
					}
					else if(!isSimulation){
						eventTree->SetBranchAddress("event",&rawPtr);
						eventTree->GetEvent(event);
						realAtriEvPtr = new UsefulAtriStationEvent(rawPtr,AraCalType::kLatestCalib);
					}

					// check to see if this event is experiencing a digitizer glitch
					// can only do this with data
					if(!isSimulation){
						vector<TGraph*> spareElecChanGraphs;
						spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(6));
						spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(14));
						spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(22));
						spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(30));
						hasBadSpareChanIssue_out = hasSpareChannelIssue(spareElecChanGraphs);
						hasBadSpareChanIssue2_out = hasSpareChannelIssue_v2(spareElecChanGraphs, station);

						if(pol==0 && 1==2){
							vector<TGraph*> dummy;
							for(int i=0; i<4; i++){
								vector<double> thisX;
								vector<double> thisY;
								thisY.push_back(-200);
								thisY.push_back(300);
								thisX.push_back(0);
								thisX.push_back(400);
								dummy.push_back(new TGraph(thisX.size(), &thisX[0], &thisY[0]));
								dummy[i]->GetXaxis()->SetTitle("Time (ns)");
								dummy[i]->GetYaxis()->SetTitle("Voltage (mV)");
								if(i==0)
									dummy[i]->SetTitle("Chan 6");
								if(i==1)
									dummy[i]->SetTitle("Chan 14");
								if(i==2)
									dummy[i]->SetTitle("Chan 22");
								if(i==3)
									dummy[i]->SetTitle("Chan 30");

							}

							TCanvas *cSpare = new TCanvas("","",4*850,850);
							cSpare->Divide(4,1);
							for(int i=0; i<4; i++){
								cSpare->cd(i+1);
								dummy[i]->Draw("AP");
								dummy[i]->SetLineColor(kWhite);
								spareElecChanGraphs[i]->Draw("sameL");
							}
							char save_temp_title[300];
							sprintf(save_temp_title,"%s/trouble_events/CWIDissue_%d.%d.%d_Run%d_Ev%d_ProblemPol%d_SpareElecChans.png",plotPath,year_now,month_now,day_now,runNum,event,pol);
							cSpare->SaveAs(save_temp_title);
							delete cSpare;
						}

						deleteGraphVector(spareElecChanGraphs);
					}

					stringstream ss1;
					string xLabel, yLabel;
					xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
					vector<string> titlesForGraphs;
					for (int i = 0; i < 16; i++){
						ss1.str("");
						ss << "Channel " << i;
						titlesForGraphs.push_back(ss1.str());
					}


					// vector <TGraph*> grWaveformsRaw = makeGraphsFromRF(realAtriEvPtr,16,xLabel,yLabel,titlesForGraphs);
					// for(int i=0; i<16; i++){
					// 	printf("Number of points is %d \n", grWaveformsRaw[i]->GetN());
					// }

					// printf(GREEN"     Event %d has bad channel %d and %d  \n"RESET,event,hasBadSpareChanIssue_out, hasBadSpareChanIssue2_out);

					/*
					if it's not in need of re-filtering, check the "top" reco again
					*/

					// printf("CW Status for Run %4d Event %5d, Pol %d: CW fwd %d, CW back %d, CW Baseline %d \n", runNum, event, pol, isCutonCW_fwd[pol], isCutonCW_back[pol], isCutonCW_baseline[pol]);
					// if((!isCutonCW_fwd[pol] && !isCutonCW_back[pol] && !isCutonCW_baseline[pol]) && !hasBadSpareChanIssue_out && !hasBadSpareChanIssue2_out){
					// // if((!isCutonCW_fwd[pol] && !isCutonCW_back[pol] && !isCutonCW_baseline[pol]) && !hasBadSpareChanIssue && 2==3){

					// 	// printf(BLUE"          Okay, so, no CW, yay\n"RESET);
					// 	vector <int> chan_list_V;
					// 	vector <int> chan_list_H;
						
					// 	chan_list_V.clear();
					// 	chan_list_V.push_back(0);
					// 	chan_list_V.push_back(1);
					// 	chan_list_V.push_back(2);

					// 	chan_list_H.clear();
					// 	chan_list_H.push_back(8);
					// 	chan_list_H.push_back(9);
					// 	chan_list_H.push_back(10);

					// 	// this looks weird here because we're aiming for what channels to *include*

					// 	if(
					// 		!(
					// 			dropBadChans
					// 			&& station==3
					// 			&& (
					// 				(!isSimulation && runNum>getA3BadRunBoundary())
					// 				||
					// 				(isSimulation && config>2)
					// 			)
					// 		)
					// 	){
					// 		chan_list_V.push_back(3);
					// 		chan_list_H.push_back(11);
					// 	}


					// 	vector<double> chan_SNRs;
					// 	for(int i=0; i<16; i++){
					// 		chan_SNRs.push_back(VPeakOverRMS[i]);
					// 	}

					// 	int solNum=0;

					// 	TH2D *map_300m_top;
					// 	if(pol==0){
					// 		// map_300m_top = theCorrelators[1]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V);
					// 		map_300m_top = theCorrelators[1]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V, chan_SNRs, solNum);
					// 	}
					// 	if(pol==1){
					// 		// map_300m_top = theCorrelators[1]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H);
					// 		map_300m_top = theCorrelators[1]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H, chan_SNRs, solNum);
					// 	}

					// 	bool printReco=false;
					// 	if(printReco){
					// 		TCanvas *cMaps = new TCanvas("","",2*1100,2*850);
					// 		map_300m_top->Draw("colz");
					// 		char save_temp_title[400];		
					// 		sprintf(save_temp_title,"%s/trouble_events/%d.%d.%d_Run%d_Ev%d_Maps_InsideSaveVals.png",plotPath,year_now,month_now,day_now,runNum,event);
					// 		cMaps->SaveAs(save_temp_title);
					// 		delete cMaps;
					// 	}
					// 	int PeakTheta_Recompute_300m_top;
					// 	int PeakPhi_Recompute_300m_top;
					// 	double PeakCorr_Recompute_300m_top;
					// 	getCorrMapPeak(map_300m_top,PeakTheta_Recompute_300m_top,PeakPhi_Recompute_300m_top,PeakCorr_Recompute_300m_top);
					// 	// cout<<"                    New theta peak in map is "<<PeakTheta_Recompute_300m_top<<endl;

					// 	if(PeakTheta_Recompute_300m_top>=37){
					// 		// cout<<"                         Ah-ha! This is a top surface event!"<<endl;
					// 		isSurfEvent_top[pol]=1;
					// 	}

					// 	delete map_300m_top;
					// }

					// and now to do *filtering*
					// if((isCutonCW_fwd[pol] || isCutonCW_back[pol] || isCutonCW_baseline[pol]) && !hasBadSpareChanIssue_out && !hasBadSpareChanIssue2_out){
					// // if((isCutonCW_fwd[pol] || isCutonCW_back[pol] || isCutonCW_baseline[pol]) && !hasBadSpareChanIssue && 2==3){
					// 	isCW_out=1;
					// 	Refilt[pol]=1;

					// 	printf(RED"	Need to filter event %d in pol %d \n"RESET,event,pol);

					// 	//get the frequencies to notch
					// 	vector<double> badFreqListLocal_fwd;
					// 	vector <double> badFreqListLocal_back;
					// 	vector <double> mergedFreqList;

					// 	//merge the two lists of frequencies
					// 	//if it's cut going both forward and backward
					// 	if(isCutonCW_fwd[pol] && isCutonCW_back[pol]){
					// 		badFreqListLocal_fwd=badFreqs_fwd->at(pol);
					// 		badFreqListLocal_back=badFreqs_back->at(pol);
					// 		for(int iFreq=0; iFreq<badFreqListLocal_fwd.size(); iFreq++){
					// 			mergedFreqList.push_back(badFreqListLocal_fwd[iFreq]);
					// 		}
					// 		for(int iFreq=0; iFreq<badFreqListLocal_back.size(); iFreq++){
					// 			double new_freq=badFreqListLocal_back[iFreq];
					// 			for(int iFreqOld=0; iFreqOld<badFreqListLocal_fwd.size(); iFreqOld++){
					// 				if(abs(new_freq-mergedFreqList[iFreqOld])>0.1){
					// 					mergedFreqList.push_back(new_freq);
					// 				}
					// 			}
					// 		}
					// 	}
					// 	//if it's cut only going forward
					// 	else if(isCutonCW_fwd[pol] && !isCutonCW_back[pol]){
					// 		badFreqListLocal_fwd=badFreqs_fwd->at(pol);
					// 		for(int iFreq=0; iFreq<badFreqListLocal_fwd.size(); iFreq++){
					// 			mergedFreqList.push_back(badFreqListLocal_fwd[iFreq]);
					// 		}
					// 	}
					// 	//if it's cut only going backward
					// 	else if(!isCutonCW_fwd[pol] && isCutonCW_back[pol]){
					// 		badFreqListLocal_back=badFreqs_back->at(pol);
					// 		for(int iFreq=0; iFreq<badFreqListLocal_back.size(); iFreq++){
					// 			mergedFreqList.push_back(badFreqListLocal_back[iFreq]);
					// 		}
					// 	}

					// 	vector<double> more_freqs_to_add;
					// 	vector<double> badFreqListLocal_baseline = badFreqs_baseline->at(pol);
					// 	if(mergedFreqList.size()>0){ //do we already have frequencies to check against?
					// 		//loop over everything identified by the CW baseline cut
					// 		for(int newFreq=0; newFreq<badFreqListLocal_baseline.size(); newFreq++){
					// 			double new_freq = badFreqListLocal_baseline[newFreq];
					// 			//now, loop over everything already in the list
					// 			for(int oldFreq=0; oldFreq<mergedFreqList.size(); oldFreq++){
					// 				//if there's a genuinely new frequency, add it to the list of things to be adde
					// 				if(abs(new_freq-mergedFreqList[oldFreq])>0.1){
					// 					more_freqs_to_add.push_back(new_freq);
					// 				}
					// 			}
					// 		}
					// 	}
					// 	else{ //otherwise we take only those found by the CW ID cut
					// 		for(int newFreq=0; newFreq<badFreqListLocal_baseline.size(); newFreq++){
					// 			double new_freq = badFreqListLocal_baseline[newFreq];
					// 			more_freqs_to_add.push_back(new_freq);
					// 		}
					// 	}

					// 	//now actually add it to the merged freq list
					// 	for(int iFreq=0; iFreq<more_freqs_to_add.size(); iFreq++){
					// 		mergedFreqList.push_back(more_freqs_to_add[iFreq]);
					// 	}

					// 	//they need to be in smaller -> larger order for notching
					// 	sort(mergedFreqList.begin(), mergedFreqList.end());

					// 	//identify the unique center frequencies and the bandwidths around them
					// 	vector <double> uniqueNotchFreqs;
					// 	vector <double> uniqueNotchBands;
					// 	for(int iFreq=0; iFreq<mergedFreqList.size(); iFreq++){
					// 		// cout<<"Frequency "<<iFreq<<" to be notched is "<<mergedFreqList[iFreq]<<endl;
					// 	}

					// 	theCorrelators[0]->pickFreqsAndBands(mergedFreqList,uniqueNotchFreqs,uniqueNotchBands);
					// 	for (int i = 0; i < uniqueNotchFreqs.size(); ++i)
					// 	{
					// 		// printf("Unique freq to be notched is %.2f with width %.2f \n", uniqueNotchFreqs[i],uniqueNotchBands[i]);
					// 	}

					// 	// for(int iFreq=0; iFreq<uniqueNotchFreqs.size(); iFreq++)
					// 	// 	printf("				Unique freq %d is %.2f with band %.2f\n",iFreq,uniqueNotchFreqs[iFreq],uniqueNotchBands[iFreq]);

					// 	/*
					// 		First we must re-do the SNR calculation (so much code...)
					// 	*/

					// 	stringstream ss1;
					// 	string xLabel, yLabel;
					// 	xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
					// 	vector<string> titlesForGraphs;
					// 	for (int i = 0; i < 16; i++){
					// 		ss1.str("");
					// 		ss << "Channel " << i;
					// 		titlesForGraphs.push_back(ss1.str());
					// 	}

					// 	vector <TGraph*> grWaveformsRaw = makeGraphsFromRF(realAtriEvPtr,16,xLabel,yLabel,titlesForGraphs);
					// 	vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(grWaveformsRaw, interpolationTimeStep, xLabel, yLabel, titlesForGraphs);
					// 	vector<TGraph*> grWaveformsPadded = makePaddedGraphs(grWaveformsInt, 0, xLabel, yLabel, titlesForGraphs);
					// 	vector <TGraph*> grNotched;
					// 	for(int i=0; i<16; i++){
					// 		TGraph *grNotchAmp = theCorrelators[0]->applyAdaptiveFilter_singleAnt_FiltMany(grWaveformsPadded[i],uniqueNotchFreqs,uniqueNotchBands);
					// 		grNotched.push_back(theCorrelators[0]->GeometricFilter(grNotchAmp,uniqueNotchFreqs,uniqueNotchBands,uniqueNotchFreqs));
					// 		delete grNotchAmp;
					// 	}
					// 	vector<TGraph*> grWaveformsPowerSpectrum = makePowerSpectrumGraphs(grWaveformsPadded, xLabel, yLabel, titlesForGraphs);
					// 	vector<TGraph*> grWaveformsPowerSpectrum_notched = makePowerSpectrumGraphs(grNotched, xLabel, yLabel, titlesForGraphs);

					// 	for(int i=0; i<16; i++){
					// 		double original_power=0.;
					// 		double final_power=0.;
					// 		for(int samp=0; samp<grWaveformsPowerSpectrum[i]->GetN(); samp++){
					// 			original_power+=grWaveformsPowerSpectrum[i]->GetY()[samp];
					// 		}
					// 		for(int samp=0; samp<grWaveformsPowerSpectrum_notched[i]->GetN(); samp++){
					// 			final_power+=grWaveformsPowerSpectrum_notched[i]->GetY()[samp];
					// 		}
					// 		if(i<8)
					// 			frac_of_power_notched_V[i]=(original_power-final_power)/original_power;
					// 		else
					// 			frac_of_power_notched_H[i-8]=(original_power-final_power)/original_power;
					// 	}

					// 	//okay, now to do the filtering
					// 	vector<double> vVPeak_new;
					// 	vector<double> vVPeakTimes_new;
					// 	double VPeak_new[16];
					// 	getAbsMaximum(grNotched, vVPeakTimes_new, vVPeak_new);
					// 	copy(vVPeak_new.begin(), vVPeak_new.begin()+16, VPeak_new); //copy these results into the arrays, because ROOT prefers arrays
					// 	double VPeakTimes_new[16];
					// 	copy(vVPeakTimes_new.begin(), vVPeakTimes_new.begin()+16, VPeakTimes_new); //copy these results into the arrays, because ROOT prefers arrays
					// 	vector<double> vWaveformRMS_new;
					// 	vWaveformRMS_new.resize(16);
						
					// 	char run_summary_name[400];
					// 	if (isSimulation == false){
					// 		sprintf(run_summary_name,"%s/RunSummary/A%d/by_config/c%d/run_summary_station_%d_run_%d.root",AuxDirPath,station,config,station,runNum);
					// 	}
					// 	else {
					// 		if(station==2){
					// 			sprintf(run_summary_name,"/fs/scratch/PAS0654/ara/sim/RunSummary/run_summary_station_2_run_20.root");
					// 		}
					// 		else if(station==3){
					// 			sprintf(run_summary_name,"/fs/scratch/PAS0654/ara/sim/RunSummary/run_summary_station_3_run_30.root");
					// 		}
					// 	}

					// 	TFile *summaryFile = TFile::Open(run_summary_name);
					// 	if(!summaryFile){
					// 		cout<<"Can't open summary file!"<<endl;
					// 		return -1;
					// 	}
					// 	TTree *SummaryTree = (TTree*) summaryFile-> Get("SummaryTree");
					// 	if(!SummaryTree){
					// 		cout<<"Can't find summaryTree for map"<<endl;
					// 		return -1;
					// 	}
						
					// 	double RMS_SoftTrigger[16];
					// 	double RMS_RFTrigger[16];
					// 	SummaryTree->SetBranchAddress("RMS_RFTrigger", &RMS_RFTrigger);
					// 	SummaryTree->SetBranchAddress("RMS_SoftTrigger", &RMS_SoftTrigger);
					// 	SummaryTree->GetEntry(0);

					// 	int nGraphs=16;
					// 	for (int i = 0; i < nGraphs; i++){ //loop over graphs
					// 		//the RMS_SoftTrigger comes out of the run summary
					// 		//so, what we want to do is see if the RMS of the software triggers was computed successfully
					// 		if (RMS_SoftTrigger[i] == RMS_SoftTrigger[i]){ //check to make sure it's not a nan
					// 			vWaveformRMS_new[i] = RMS_SoftTrigger[i];
					// 		} else { //if it was a nan, then instead we'll look at the RF trigger version
					// 			if (RMS_RFTrigger[i] == RMS_RFTrigger[i]){ //make sure it's not a nan
					// 				vWaveformRMS_new[i] = RMS_RFTrigger[i];
					// 			}
					// 		}
					// 	}
					// 	double waveformRMS_new[16];
					// 	copy(vWaveformRMS_new.begin(), vWaveformRMS_new.begin()+16, waveformRMS_new); //copy into the array
					// 	vector<double> vWaveformRMS_50ns_new;
					// 	int first50ns = (int)(50./interpolationTimeStep);
					// 	getRMS(grNotched, vWaveformRMS_50ns_new, first50ns);
					// 	double waveformRMS_50ns_new[16];
					// 	copy(vWaveformRMS_50ns_new.begin(), vWaveformRMS_50ns_new.begin()+16, waveformRMS_50ns_new); //copy those results into an array
					// 	vector<double> vVPeakOverRMS_new;
					// 	vVPeakOverRMS_new.resize(16);
					// 	for (int i = 0; i < 16; i++){
					// 		vVPeakOverRMS_new[i] = VPeak_new[i]/waveformRMS_new[i];
					// 		vVPeakOverRMS_new[i] = VPeak_new[i]/waveformRMS_new[i];
					// 	}
					// 	AraGeomTool * geomTool = new AraGeomTool();
					// 	vector<int> polarizations;
					// 	vector<int> antenna_numbers;
					// 	polarizations.resize(16);
					// 	antenna_numbers.resize(16);
					// 	vector< vector <double> > ant_loc; //will be 16x3 vector of the x,y,z's the 16 antennas
					// 	ant_loc.resize(16);
					// 	for (int i = 0; i < 16; i++){
					// 		ant_loc[i].resize(3);
					// 		ant_loc[i][0] = geomTool->getStationInfo(station)->getAntennaInfo(i)->antLocation[0];
					// 		ant_loc[i][1] = geomTool->getStationInfo(station)->getAntennaInfo(i)->antLocation[1];
					// 		ant_loc[i][2] = geomTool->getStationInfo(station)->getAntennaInfo(i)->antLocation[2];
					// 		polarizations[i] = (int)geomTool->getStationInfo(station)->getAntennaInfo(i)->polType;
					// 		antenna_numbers[i]=i;
					// 	}

					// 	vector<int> chan_exclusion_list;
					// 	if(dropBadChans){
					// 		if(station==2){
					// 			// hpol channel 15
					// 			chan_exclusion_list.push_back(15);
					// 		}
					// 		else if(station==3){
					// 			if( 
					// 				(!isSimulation && runNum>getA3BadRunBoundary())
					// 				||
					// 				(isSimulation && config>2)

					// 			){								// vpol sring 4
					// 				chan_exclusion_list.push_back(3);
					// 				chan_exclusion_list.push_back(7);
					// 				chan_exclusion_list.push_back(11);
					// 				chan_exclusion_list.push_back(15);
					// 			}
					// 		}
					// 	}

					// 	vector<double> vThirdVPeakOverRMS_new;
					// 	double thirdVPeakOverRMS_new[3];
					// 	getThirdVPeakOverRMS(vVPeakOverRMS_new, polarizations, antenna_numbers, chan_exclusion_list, vThirdVPeakOverRMS_new);
					// 	for (int i = 0 ; i< 3; i++){ //pull out the first three entries
					// 		thirdVPeakOverRMS_new[i] = vThirdVPeakOverRMS_new[i];
					// 	}

					// 	xLabel = "Time (ns)"; yLabel = "Integrated Power (arb units)";
					// 	int numBinsToIntegrate = (int)(5./interpolationTimeStep);
					// 	vector<TGraph*> grIntPower = makeIntegratedBinPowerGraphs(grNotched, numBinsToIntegrate, xLabel, yLabel, titlesForGraphs);

					// 	vector<double> hitTimes_new; //what are the hit times
					// 	vector<double> peakIntPowers_new; //what are the powers at those hit times?
					// 	getAbsMaximum(grIntPower, hitTimes_new, peakIntPowers_new);
					// 	vector<vector<double> > vvHitTimes_new; //vector of vector of hit times
					// 	vector<vector<double> > vvPeakIntPowers_new; //vector of vector of power at the those hit times
					// 	int numSearchPeaks = 2;
					// 	const int numFaces = 12;
					// 	getAbsMaximum_N(grIntPower, numSearchPeaks, 5.0, vvHitTimes_new, vvPeakIntPowers_new);
					// 	vector<double> peakIntRMS_new;       
					// 	for (int i = 0; i < peakIntPowers_new.size(); i++){
					// 		peakIntRMS_new.push_back(sqrt(peakIntPowers_new[i]/numBinsToIntegrate));
					// 	}
					// 	double avgPeakPower_5ns_new[16];
					// 	double peakPowerTimes_new[16];
					// 	for (int i = 0; i < 16; i++){
					// 		avgPeakPower_5ns_new[i] = peakIntPowers_new[i]/numBinsToIntegrate;
					// 		peakPowerTimes_new[i] = hitTimes_new[i];
					// 	}
					// 	vector<double> RMS_10overRMS_new;
					// 	for (int i = 0; i < 16; i++){
					// 		RMS_10overRMS_new.push_back(sqrt(avgPeakPower_5ns_new[i])/waveformRMS_new[i]);
					// 	}
					// 	vector<vector<double> > vvRMS_10overRMS_new;
					// 	vvRMS_10overRMS_new.resize(16);
					// 	for (int i = 0; i < 16; i++){
					// 		vvRMS_10overRMS_new[i].resize(vvPeakIntPowers_new[i].size());
					// 		for (int j = 0; j < vvPeakIntPowers_new[i].size(); j++){
					// 			vvRMS_10overRMS_new[i][j] = sqrt(vvPeakIntPowers_new[i][j]/numBinsToIntegrate)/waveformRMS_new[i];	
					// 		}
					// 	}

					// 	vector< vector< int > > pairs_V_new;
					// 	vector< vector< int > > pairs_H_new;
					// 	setupCorrelationPairs(station, pairs_V_new, pairs_H_new); //just sets up the pairs (like, 0,1, 0,2 etc) that go into the correlation

					// 	vector<double> bestTimes_V_new;
					// 	vector<double> bestCorrs_V_new;
					// 	vector<double> bestTimes_H_new;
					// 	vector<double> bestCorrs_H_new;

					// 	//now, to set up all the pairs that contribute to the faces
					// 	vector<vector<vector<vector<int> > > > faces = setupFaces(station);

					// 	//loop over the thresholds that decide if a face is allowed to contribute
					// 	const int thresholdSteps = 15;
					// 	double thresholdMin = 2.0;
					// 	double thresholdStep = 0.1;
					// 	double rms_pol_thresh_face_new_V[15][12];
					// 	double rms_pol_thresh_face_new_H[15][12];
					// 	for (int thresholdBin = 0; thresholdBin < thresholdSteps; thresholdBin++){
					// 		double threshold = thresholdMin + thresholdStep*(double)thresholdBin;
							
					// 		//get the RMS of all the faces at this threshold bin
					// 		vector<double> rms_faces_V_new = getRms_Faces_Thresh_N(vvHitTimes_new, vvRMS_10overRMS_new, threshold, 0, faces, ant_loc);
					// 		vector<double> rms_faces_H_new = getRms_Faces_Thresh_N(vvHitTimes_new, vvRMS_10overRMS_new, threshold, 1, faces, ant_loc);

					// 		for (int i = 0; i < numFaces; i++){ //loop over the faces, and store the RMS for that polarization, threshold bin, and face
					// 			rms_pol_thresh_face_new_V[thresholdBin][i] = rms_faces_V_new[i];
					// 			rms_pol_thresh_face_new_H[thresholdBin][i] = rms_faces_H_new[i];
					// 		}	
					// 	} // end threshold scan

					// 	//now the dropped channel case
					// 	vector<vector<vector<vector<int> > > > faces_drop = setupFaces(station, dropBadChans);
					// 	double rms_pol_thresh_face_new_V_drop[15][numFaces_new_V];
					// 	double rms_pol_thresh_face_new_H_drop[15][numFaces_new_H];
					// 	for (int thresholdBin = 0; thresholdBin < thresholdSteps; thresholdBin++){
					// 		double threshold = thresholdMin + thresholdStep*(double)thresholdBin;
							
					// 		//get the RMS of all the faces at this threshold bin
					// 		vector<double> rms_faces_V_new_drop = getRms_Faces_Thresh_N(vvHitTimes_new, vvRMS_10overRMS_new, threshold, 0, faces_drop, ant_loc);
					// 		vector<double> rms_faces_H_new_drop = getRms_Faces_Thresh_N(vvHitTimes_new, vvRMS_10overRMS_new, threshold, 1, faces_drop, ant_loc);

					// 		for (int i = 0; i < numFaces_new_V; i++){
					// 			rms_pol_thresh_face_new_V_drop[thresholdBin][i] = rms_faces_V_new_drop[i];
					// 		}
					// 		for (int i = 0; i < numFaces_new_H; i++){
					// 			rms_pol_thresh_face_new_H_drop[thresholdBin][i] = rms_faces_H_new_drop[i];
					// 		}
					// 	} // end threshold scan

					// 	int thresholdBin_pol_new[]={thresholdBin_pol[0],thresholdBin_pol[1]};

					// 	vector <double> rms_faces_V_new;
					// 	vector <double> rms_faces_H_new;

					// 	if(dropBadChans){
					// 		int num_faces_for_V_loop;
					// 		int num_faces_for_H_loop;

					// 		bool needToSwitchToAltWFRMSArrays_local=false;
					// 		// we need an actual flag for if we need to switch to the _alternate_ arrays

					// 		if(station==2){
					// 			rms_faces_V_new.resize(numFaces);
					// 			num_faces_for_V_loop=numFaces;
					// 			rms_faces_H_new.resize(numFaces_A2_drop);
					// 			num_faces_for_H_loop=numFaces_A2_drop;
					// 			needToSwitchToAltWFRMSArrays_local=true;
					// 		}
					// 		else if(station==3){
					// 			if(
					// 				(!isSimulation && runNum>getA3BadRunBoundary())
					// 				||
					// 				(isSimulation && config>2)
					// 			){ //it's 2014+, drop string four
					// 				rms_faces_V_new.resize(numFaces_A3_drop);
					// 				num_faces_for_V_loop=numFaces_A3_drop;
					// 				rms_faces_H_new.resize(numFaces_A3_drop);
					// 				num_faces_for_H_loop=numFaces_A3_drop;
					// 				needToSwitchToAltWFRMSArrays_local=true;
					// 			}
					// 			else{ //it's 2013-, keep string four
					// 				rms_faces_V_new.resize(numFaces);
					// 				num_faces_for_V_loop=numFaces;
					// 				rms_faces_H_new.resize(numFaces);
					// 				num_faces_for_H_loop=numFaces;
					// 			}
					// 		}

					// 		if(needToSwitchToAltWFRMSArrays_local){
					// 			//now we loop over the faces
					// 			for(int i=0; i<num_faces_for_V_loop; i++){
					// 				rms_faces_V_new[i] = rms_pol_thresh_face_new_V_drop[thresholdBin_pol_new[0]][i];
					// 			}
					// 			for(int i=0; i<num_faces_for_H_loop; i++){
					// 				rms_faces_H_new[i] = rms_pol_thresh_face_new_H_drop[thresholdBin_pol_new[1]][i];
					// 			}
					// 		}
					// 		else{
					// 			//now we loop over the faces
					// 			for(int i=0; i<num_faces_for_V_loop; i++){
					// 				rms_faces_V_new[i] = rms_pol_thresh_face_new_V[thresholdBin_pol_new[0]][i];
					// 			}
					// 			for(int i=0; i<num_faces_for_H_loop; i++){
					// 				rms_faces_H_new[i] = rms_pol_thresh_face_new_H[thresholdBin_pol_new[1]][i];
					// 			}
					// 		}
					// 	}
					// 	else{
					// 		rms_faces_V_new.resize(12);
					// 		rms_faces_H_new.resize(12);
					// 		//now, we must loop over the faces
					// 		for(int i=0; i<12; i++){
					// 			rms_faces_V_new[i] = rms_pol_thresh_face_new_V[thresholdBin_pol_new[0]][i];  //this is right RMS for this polarization, threshold requirement, and face
					// 			rms_faces_H_new[i] = rms_pol_thresh_face_new_H[thresholdBin_pol_new[1]][i];
					// 		}
					// 	}

					// 	sort(rms_faces_V_new.begin(), rms_faces_V_new.end());
					// 	sort(rms_faces_H_new.begin(), rms_faces_H_new.end());
						
					// 	double bestFaceRMS_new[2];
					// 	bestFaceRMS_new[0]=rms_faces_V_new[0];
					// 	bestFaceRMS_new[1]=rms_faces_H_new[0];

					// 	double SNRs_new[2];
					// 	SNRs_new[0] = vThirdVPeakOverRMS_new[0];
					// 	SNRs_new[1] = vThirdVPeakOverRMS_new[1];

					// 	/*
					// 		Now we must re-do the reconstructions
					// 	*/
						
					// 	vector <int> chan_list_V;
					// 	vector <int> chan_list_H;
					// 	for(int chan=0; chan<=7; chan++){
					// 		chan_list_V.push_back(chan);
					// 		chan_list_H.push_back(chan+8);
					// 	}
					// 	if(dropBadChans){
					// 		if(station==2){
					// 			//for station 2, we need to exclude channel 15 from the analysis
					// 			chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
					// 		}
					// 		else if(station==3){
					// 			//for station 3, remove for data after getA3BadRunBoundary, or for sim after config 2
					// 			if( 
					// 				(!isSimulation && runNum>getA3BadRunBoundary())
					// 				||
					// 				(isSimulation && config>2)

					// 			){
					// 				chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 3), chan_list_V.end());
					// 				chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 7), chan_list_V.end());

					// 				chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 11), chan_list_H.end());
					// 				chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
					// 			}
					// 		}
					// 	}

					// 	// get those filtered SNR's we just computed
					// 	vector<double> chan_SNRs;
					// 	for(int i=0; i<16; i++){
					// 		chan_SNRs.push_back(vVPeakOverRMS_new[i]);
					// 	}

					// 	TH2D *map_30m;
					// 	TH2D *map_300m;
					// 	int solNum = 0;
					// 	if(pol==0){
					// 		// map_30m = theCorrelators[0]->getInterferometricMap_RT_FiltMany_select(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V, solNum,-1,uniqueNotchFreqs,uniqueNotchBands);
					// 		// map_300m = theCorrelators[1]->getInterferometricMap_RT_FiltMany_select(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V, solNum,-1,uniqueNotchFreqs,uniqueNotchBands);

					// 		map_30m = theCorrelators[0]->getInterferometricMap_RT_FiltMany_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V, chan_SNRs, solNum,-1,uniqueNotchFreqs,uniqueNotchBands);
					// 		map_300m = theCorrelators[1]->getInterferometricMap_RT_FiltMany_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V, chan_SNRs, solNum,-1,uniqueNotchFreqs,uniqueNotchBands);

					// 	}
					// 	if(pol==1){
					// 		// map_30m = theCorrelators[0]->getInterferometricMap_RT_FiltMany_select(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H, solNum,-1,uniqueNotchFreqs,uniqueNotchBands);
					// 		// map_300m = theCorrelators[1]->getInterferometricMap_RT_FiltMany_select(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H, solNum,-1,uniqueNotchFreqs,uniqueNotchBands);

					// 		map_30m = theCorrelators[0]->getInterferometricMap_RT_FiltMany_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H, chan_SNRs, solNum,-1,uniqueNotchFreqs,uniqueNotchBands);
					// 		map_300m = theCorrelators[1]->getInterferometricMap_RT_FiltMany_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H, chan_SNRs, solNum,-1,uniqueNotchFreqs,uniqueNotchBands);

					// 	}
					// 	int PeakTheta_Recompute_30m;
					// 	int PeakTheta_Recompute_300m;
					// 	int PeakPhi_Recompute_30m;
					// 	int PeakPhi_Recompute_300m;
					// 	double PeakCorr_Recompute_30m;
					// 	double PeakCorr_Recompute_300m;
					// 	getCorrMapPeak(map_30m,PeakTheta_Recompute_30m,PeakPhi_Recompute_30m,PeakCorr_Recompute_30m);
					// 	getCorrMapPeak(map_300m,PeakTheta_Recompute_300m,PeakPhi_Recompute_300m,PeakCorr_Recompute_300m);

					// 	// theta_300[pol]=PeakTheta_Recompute_300m;
					// 	// phi_300[pol]=PeakPhi_Recompute_300m;
					// 	// theta_41[pol]=PeakTheta_Recompute_30m;
					// 	// phi_41[pol]=PeakPhi_Recompute_30m;

					// 	// update--only change the _new verions
					// 	theta_300_new[pol]=PeakTheta_Recompute_300m;
					// 	phi_300_new[pol]=PeakPhi_Recompute_300m;
					// 	theta_41_new[pol]=PeakTheta_Recompute_30m;
					// 	phi_41_new[pol]=PeakPhi_Recompute_30m;						
						
					// 	chan_list_V.clear();
					// 	chan_list_V.push_back(0);
					// 	chan_list_V.push_back(1);
					// 	chan_list_V.push_back(2);

					// 	chan_list_H.clear();
					// 	chan_list_H.push_back(8);
					// 	chan_list_H.push_back(9);
					// 	chan_list_H.push_back(10);

					// 	// this looks weird here because we're aiming for what channels to *include*

					// 	if(
					// 		!(
					// 			dropBadChans
					// 			&& station==3
					// 			&& (
					// 				(!isSimulation && runNum>getA3BadRunBoundary())
					// 				||
					// 				(isSimulation && config>2)
					// 			)
					// 		)
					// 	){
					// 		chan_list_V.push_back(3);
					// 		chan_list_H.push_back(11);
					// 	}

					// 	TH2D *map_300m_top;
					// 	if(pol==0){
					// 		map_300m_top = theCorrelators[1]->getInterferometricMap_RT_FiltMany_select(settings, detector, realAtriEvPtr, Vpol, 0, chan_list_V, 0,-1,uniqueNotchFreqs,uniqueNotchBands);
					// 	}
					// 	if(pol==1){
					// 		map_300m_top = theCorrelators[1]->getInterferometricMap_RT_FiltMany_select(settings, detector, realAtriEvPtr, Hpol, 0, chan_list_H, 0,-1,uniqueNotchFreqs,uniqueNotchBands);
					// 	}

					// 	int PeakTheta_Recompute_300m_top;
					// 	int PeakPhi_Recompute_300m_top;
					// 	double PeakCorr_Recompute_300m_top;
					// 	getCorrMapPeak(map_300m_top,PeakTheta_Recompute_300m_top,PeakPhi_Recompute_300m_top,PeakCorr_Recompute_300m_top);

					// 	//cleanup
					// 	deleteGraphVector(grWaveformsRaw);
					// 	deleteGraphVector(grWaveformsInt);
					// 	deleteGraphVector(grWaveformsPadded);
					// 	deleteGraphVector(grIntPower);
					// 	deleteGraphVector(grNotched);
					// 	deleteGraphVector(grWaveformsPowerSpectrum);
					// 	deleteGraphVector(grWaveformsPowerSpectrum_notched);

					// 	printf("		old vs new logrms calc in pol %d: %.2f vs %.2f \n",pol,log(bestFaceRMS[pol])/log(10),log(bestFaceRMS_new[pol])/log(10));
					// 	printf("		old vs new snr in pol %d: %.2f vs %.2f \n",pol,SNRs[pol],SNRs_new[pol] );
					// 	printf("		old vs new corr in pol %d: %.4f vs %.4f \n",pol,corr_val_org[pol],PeakCorr_Recompute_300m);


					// 	// re-check top face reco
					// 	if(PeakTheta_Recompute_300m_top>=37)
					// 		isSurfEvent_top[pol]=1;

					// 	// re-check surface cut
					// 	if(PeakTheta_Recompute_300m>=37){
					// 		// isSurfEvent[pol]=1;
					// 		isSurfEvent_new_out[pol]=1;
					// 	}
					// 	else{
					// 		// isSurfEvent[pol]=0;
					// 		isSurfEvent_new_out[pol]=0;
					// 	}

					// 	// re-check wavefront RMS cut
					// 	if(log(bestFaceRMS_new[pol])/log(10) < wavefrontRMScut[pol]){
					// 		// WFRMS[pol]=0;
					// 		WFRMS_new[pol]=0;
					// 	}
					// 	else{
					// 		// WFRMS[pol]=1;
					// 		WFRMS_new[pol]=1;
					// 	}

					// 	// assign the newly computed corr and snr values
					// 	// corr_val[pol]=PeakCorr_Recompute_300m;
					// 	// snr_val[pol]=SNRs_new[pol];

					// 	corr_val_new[pol]=PeakCorr_Recompute_300m;
					// 	snr_val_new[pol]=SNRs_new[pol];

					// 	delete map_300m;
					// 	delete map_300m_top;
					// 	delete map_30m;

					// 	// yes, the summary file close and delete needs to be down here
					// 	// to make ROOT's silly ownership thingy work out and not cause segfault
					// 	summaryFile->Close();
					// 	delete summaryFile;
					// } //if any frequencies are flagged for filtering
					if(!isSimulation) delete realAtriEvPtr;
					mapFile->Close();
					delete mapFile;
				} //cal pulser
			}//loop over polarization
		}//loop over events
		inputFile->Close();
		NewCWFile->Close();
		delete inputFile;

		for(int pol=0; pol<2; pol++){
			printf("Num baseline filt over num total: %6d/%6d = %.2f \n", num_baseline_filt[pol], num_total[pol], double(num_baseline_filt[pol])/double(num_total[pol]));
		}
		printf("Done! Run Number %d \n", runNum);
	} //end loop over input files

	TCanvas *cFilterFreqs = new TCanvas("","",850,850);
	filter_freqs->Draw("");
	filter_freqs->GetXaxis()->SetTitle("Frequency (MHz)");
	filter_freqs->GetYaxis()->SetTitle("Number of Events");
	gPad->SetLogy();
	char this_title[400];
	sprintf(this_title,"%s/trouble_events/CWIDissue_%d.%d.%d_A%d_c%d_DistOfFilterFreqs.png",plotPath,year_now,month_now,day_now,station,config);
	cFilterFreqs->SaveAs(this_title);

	TCanvas *cFilterFreqsVsTime = new TCanvas("","",2*850,850);
	h2_filt_freq_vs_time->Draw("colz");
		h2_filt_freq_vs_time->GetXaxis()->SetTitle("UnixTime [YY/MM]");
		h2_filt_freq_vs_time->GetYaxis()->SetTitle("Contaminated Frequency [MHz]");
		h2_filt_freq_vs_time->GetYaxis()->SetRangeUser(100,900);
	sprintf(this_title,"%s/trouble_events/CWIDissue_%d.%d.%d_A%d_c%d_DistOfFilterFreqs_vs_Time.png",plotPath,year_now,month_now,day_now,station,config);
	cFilterFreqsVsTime->SaveAs(this_title);

	double low_edge;
	double high_edge;
	if(config==1){
		low_edge=1400;
		high_edge=2000;
	}
	if(config==2){
		low_edge=500;
		high_edge=1400;
	}
	if(config==3){
		low_edge=3000;
		high_edge=8000;
	}
	if(config==4){
		low_edge=6000;
		high_edge=8000;
	}
	if(config==5){
		low_edge=1800;
		high_edge=3200;
	}

	TCanvas *cFilterFreqsVsRun = new TCanvas("","",2*850,850);
	h2_filt_freq_vs_run->Draw("colz");
		h2_filt_freq_vs_run->GetXaxis()->SetTitle("Run Number");
		h2_filt_freq_vs_run->GetYaxis()->SetTitle("Contaminated Frequency [MHz]");
		h2_filt_freq_vs_run->GetYaxis()->SetRangeUser(100,900);
		h2_filt_freq_vs_run->GetXaxis()->SetRangeUser(low_edge,high_edge);
	sprintf(this_title,"%s/trouble_events/CWIDissue_%d.%d.%d_A%d_c%d_DistOfFilterFreqs_vs_Run.png",plotPath,year_now,month_now,day_now,station,config);
	cFilterFreqsVsRun->SaveAs(this_title);	
}

int PlotThisEvent(int station, int config, int runNum, int event, int problempol){
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char *DataDirPath(getenv("DATA_DIR_PROJECT"));
	if (DataDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;
	char *AuxDirPath(getenv("AUX_DIR"));
	if (AuxDirPath == NULL){
		std::cout << "Warning! $AUX_DIR is not set! You need this for CWID and RunSummary files" << endl;
		return -1;
	}

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

	for(int pol=0; pol<2; pol++){
		char run_summary_filename[400];
		sprintf(run_summary_filename,"%s/RunSummary/A%d/all_runs/run_summary_station_%d_run_%d.root",AuxDirPath,station,station,runNum);
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
		vector<double> baseline_CW_freqs = CWCut_TB(waveforms, average, pol, 6., 5.5, station, 3, chan_exclusion_list, runNum, event, true);
	}




	vector<TGraph*> dummy;
	for(int i=0; i<16; i++){
		vector<double> thisX;
		vector<double> thisY;
		thisY.push_back(-400);
		thisY.push_back(400);
		thisX.push_back(-200);
		thisX.push_back(700);
		dummy.push_back(new TGraph(thisX.size(), &thisX[0], &thisY[0]));
	}

	char save_temp_title[300];
	sprintf(save_temp_title,"%s/trouble_events/CWIDissue_%d.%d.%d_Run%d_Ev%d_ProblemPol%d_Waveforms.png",plotPath,year_now,month_now,day_now,runNum,event,problempol);
	TCanvas *cWave = new TCanvas("","",4*1100,4*850);
	cWave->Divide(4,4);
	for(int i=0; i<16; i++){
		cWave->cd(i+1);
		dummy[i]->Draw("AL");
		dummy[i]->SetLineColor(kWhite);
		// dummy[i]->GetXaxis()->SetRangeUser(300.,500.);

		waveforms[i]->Draw("sameL");
		waveforms[i]->SetLineWidth(3);
		// waveforms[i]->GetXaxis()->SetRangeUser(300.,500.);
	}
	cWave->SaveAs(save_temp_title);
	delete cWave;

	sprintf(save_temp_title,"%s/trouble_events/CWIDissue_%d.%d.%d_Run%d_Ev%d_ProblemPol%d_Spectra.png",plotPath,year_now,month_now,day_now,runNum,event,problempol);
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