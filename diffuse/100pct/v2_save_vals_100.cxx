////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	v2_analysis_save_vals_100.cxx 
////	A23 diffuse, save values for cuts for 100pct sample
////	Should rely on "reco" files and "CWID" files only
////
////	Aug 2019
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

//AraRoot includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "Settings.h"
#include "Detector.h"
#include "Report.h"
#include "RayTraceCorrelator.h"
AraAntPol::AraAntPol_t Vpol = AraAntPol::kVertical;
AraAntPol::AraAntPol_t Hpol = AraAntPol::kHorizontal;
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

	char *SimDirPath(getenv("SIM_DIR"));
	if (SimDirPath == NULL){
		std::cout << "Warning! $SIM_DIR is not set!" << endl;
		return -1;
	}
	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL){
		std::cout << "Warning! $PED_DIR is not set! This is needed to find pedestal files." << endl;
		return -1;
	}
	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL){
		std::cout << "Warning! $PLOT_PATH is not set!" << endl;
		return -1;
	}

	stringstream ss;
	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	
	if(argc<15){
	  cout<< "Usage is: " << argv[0] << " <isSim?> <station> <config> <year_or_energy (as float, eg 17.0 or 18.5)> <drop_bad_chan> <output_location> <data_directory> <cwid_directory> <run_summary_directory> <V SNR bin> <H SNR bin> <V WFRMS val> <H WFRMS val> <100pct reco filename 1> <100pct reco filename 2 > ... <100pct reco filename x>"<<endl;
	  return 0;	
	}
	int isSimulation = atoi(argv[1]);
	int station = atoi(argv[2]);
	int config = atoi(argv[3]);
	double year_or_energy = double(atof(argv[4]));
	int dropBadChans = atoi(argv[5]);	
	string output_location = argv[6];
	string data_directory = argv[7];
	string cw_directory = argv[8];
	string runsummary_directory = argv[9];

	//just to have the cut parameters up front and easy to find
	int thresholdBin_pol[]={atoi(argv[10]), atoi(argv[11])}; //bin 0 = 2.0, bin 0 = 2.0 //what is the faceRMS SNR inclusion threshold?
	double wavefrontRMScut[]={atof(argv[12]),atof(argv[13])}; //event wavefrontRMS < this value

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

	for(int file_num=14; file_num<argc; file_num++){

		string file = string(argv[file_num]);
		string wordRun = "run_";
		size_t foundRun = file.find(wordRun);
		string wordFilter = "_reco";
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

		char outfile_name[300];
		sprintf(outfile_name,"%s/cutvals_snrbins_%d_%d_wfrmsvals_%.1f_%.1f_run_%d.root",output_location.c_str(),thresholdBin_pol[0], thresholdBin_pol[1], (wavefrontRMScut[0]),(wavefrontRMScut[1]),runNum);
		if(dropBadChans){
			sprintf(outfile_name,"%s/cutvals_drop_FiltSurface_snrbins_%d_%d_wfrmsvals_%.1f_%.1f_run_%d.root",output_location.c_str(),thresholdBin_pol[0], thresholdBin_pol[1], (wavefrontRMScut[0]), (wavefrontRMScut[1]),runNum);
			// sprintf(outfile_name,"%s/cutvals_drop_snrbins_%s%d_%s%d_wfrmsvals_%.1f_%.1f_run_%d.root",output_location.c_str(),thresholdBin_pol[0], thresholdBin_pol[1], wavefrontRMScut[0]<0?'-':'',(wavefrontRMScut[0]), wavefrontRMScut[0]<0?'-':'',(wavefrontRMScut[1]),runNum);
		}
		TFile *fpOut = new TFile(outfile_name,"recreate");
		TTree *trees[5];
		trees[0]= new TTree("VTree","VTree");
		trees[1]= new TTree("HTree","HTree");
		trees[2]= new TTree("AllTree","AllTree");

		// need "org" and "new" variables so we know what's happening before and after filtering
		// the way I'm doing this is so, so dumb. Brian, don't ever do it this way again
		// nothing here really needs to change (I think)

		double corr_val_org[2];
		double snr_val_org[2];
		int WFRMS_org[2];
		int theta_300_org[2];
		int phi_300_org[2];
		int theta_41_org[2];
		int phi_41_org[2];

		trees[0]->Branch("corr_val_V_org",&corr_val_org[0]);
		trees[0]->Branch("snr_val_V_org",&snr_val_org[0]);
		trees[0]->Branch("wfrms_val_V_org",&WFRMS_org[0]);
		trees[0]->Branch("theta_300_V_org",&theta_300_org[0]);
		trees[0]->Branch("theta_41_V_org",&theta_41_org[0]);
		trees[0]->Branch("phi_300_V_org",&phi_300_org[0]);
		trees[0]->Branch("phi_41_V_org",&phi_41_org[0]);

		trees[1]->Branch("corr_val_H_org",&corr_val_org[1]);
		trees[1]->Branch("snr_val_H_org",&snr_val_org[1]);
		trees[1]->Branch("wfrms_val_H_org",&WFRMS_org[1]);
		trees[1]->Branch("theta_300_H_org",&theta_300_org[1]);
		trees[1]->Branch("theta_41_H_org",&theta_41_org[1]);
		trees[1]->Branch("phi_300_H_org",&phi_300_org[1]);
		trees[1]->Branch("phi_41_H_org",&phi_41_org[1]);

		double corr_val_new[2];
		double snr_val_new[2];
		int WFRMS_new[2];
		int theta_300_new[2];
		int phi_300_new[2];
		int theta_41_new[2];
		int phi_41_new[2];

		trees[0]->Branch("corr_val_V_new",&corr_val_new[0]);
		trees[0]->Branch("snr_val_V_new",&snr_val_new[0]);
		trees[0]->Branch("wfrms_val_V_new",&WFRMS_new[0]);
		trees[0]->Branch("theta_300_V_new",&theta_300_new[0]);
		trees[0]->Branch("theta_41_V_new",&theta_41_new[0]);
		trees[0]->Branch("phi_300_V_new",&phi_300_new[0]);
		trees[0]->Branch("phi_41_V_new",&phi_41_new[0]);

		trees[1]->Branch("corr_val_H_new",&corr_val_new[1]);
		trees[1]->Branch("snr_val_H_new",&snr_val_new[1]);
		trees[1]->Branch("wfrms_val_H_new",&WFRMS_new[1]);
		trees[1]->Branch("theta_300_H_new",&theta_300_new[1]);
		trees[1]->Branch("theta_41_H_new",&theta_41_new[1]);
		trees[1]->Branch("phi_300_H_new",&phi_300_new[1]);
		trees[1]->Branch("phi_41_H_new",&phi_41_new[1]);

		int Refilt[2];	
		trees[0]->Branch("Refilt_V",&Refilt[0]);
		trees[1]->Branch("Refilt_H",&Refilt[1]);

		double frac_of_power_notched_V[8];
		double frac_of_power_notched_H[8];

		for(int i=0; i<8; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			trees[0]->Branch(ss.str().c_str(),&frac_of_power_notched_V[i]);
		}
		for(int i=8; i<16; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			trees[1]->Branch(ss.str().c_str(),&frac_of_power_notched_H[i-8]);
		}

		// mark these as "out" so we can see them clearly

		int isCal_out;
		int isSoft_out;
		int isShortWave_out;
		int isCW_out;
		int isNewBox_out;
		int isSurfEvent_org_out[2]; // originally a surf event?
		int isSurfEvent_new_out[2]; // a surface event after filtering?
		int isSurfEvent_top[2]; // a top event?
		trees[2]->Branch("cal",&isCal_out);
		trees[2]->Branch("soft",&isSoft_out);
		trees[2]->Branch("short",&isShortWave_out);
		trees[2]->Branch("CW",&isCW_out);
		trees[2]->Branch("box",&isNewBox_out);
		trees[2]->Branch("surf_V_org",&isSurfEvent_org_out[0]);
		trees[2]->Branch("surf_H_org",&isSurfEvent_org_out[1]);
		trees[2]->Branch("surf_V_new",&isSurfEvent_new_out[0]);
		trees[2]->Branch("surf_H_new",&isSurfEvent_new_out[1]);
		trees[2]->Branch("surf_top_V",&isSurfEvent_top[0]);
		trees[2]->Branch("surf_top_H",&isSurfEvent_top[1]);

		int isBadEvent_out;
		double outweight;
		int Trig_Pass_out[16];
		int unixTime_out;
		int hasBadSpareChanIssue_out;
		int hasBadSpareChanIssue2_out;
		int isFirstFiveEvent_out;
		int eventNumber_out;
		int runNum_out;

		trees[2]->Branch("bad",&isBadEvent_out);
		trees[2]->Branch("weight",&outweight);
		trees[2]->Branch("unixTime",&unixTime_out);
		trees[2]->Branch("hasBadSpareChanIssue",&hasBadSpareChanIssue_out);
		trees[2]->Branch("hasBadSpareChanIssue2",&hasBadSpareChanIssue2_out);
		trees[2]->Branch("isFirstFiveEvent",&isFirstFiveEvent_out);
		trees[2]->Branch("eventNumber",&eventNumber_out);
		trees[2]->Branch("runNum",&runNum_out);
		if(isSimulation)
			trees[2]->Branch("Trig_Pass", &Trig_Pass_out, "Trig_Pass_out[16]/I");

		cout << "Run " << file_num << " :: " << argv[file_num] << endl;

		// now, for the 100pct sample, we are just drawing information in from the
		// reco file, which will contain info about the filter and info about the reco all in one go
		// this should actually look *simpler* than the 10pct case 
		// because the output is more sensibly organized (theoretically)
		
		TFile *inputFile = TFile::Open(argv[file_num]);
		if(!inputFile){
			cout<<"Can't open joined file!"<<endl;			
			return -1;
		}
		
		// in the 100pct files, these are just called "OutputTree"
		// these are really the filter trees of course
		ss.str("");
		ss << "OutputTree";
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

		// copy over the whole filter tree
		fpOut->cd();
		trees[3]= inputTree_filter->CloneTree(0);

		// next, load the reco tree
		// this actually looks *simpler* than it did before because all the reco happened at one
		// like it should have been in the 10pct if Brian had refactored it sooner, but oh well...
		ss.str("");
		ss << "OutputTreeReco";
		TTree *inputTree_reco = (TTree*) inputFile->Get(ss.str().c_str());
		if(!inputTree_reco){
			cout<<"Can't open reco tree"<<endl;
			return -1;
		}
		// copy over the whole reco tree also, because why lose information right?
		trees[4] = inputTree_reco->CloneTree(0);
		// 2x2 array
		// first element is for radius (0=41m and 1=300)
		// second element is polarization (0=V and 1=H)

		double peakCorr[2][2];
		int peakTheta[2][2];
		int peakPhi[2][2];
		// firt correlation
		inputTree_reco->SetBranchAddress("peakCorr_41m",peakCorr[0]);
		inputTree_reco->SetBranchAddress("peakCorr_300m",peakCorr[1]);
		// then peak theta
		inputTree_reco->SetBranchAddress("peakTheta_41m",peakTheta[0]);
		inputTree_reco->SetBranchAddress("peakTheta_300m",peakTheta[1]);
		// then peak phi
		inputTree_reco->SetBranchAddress("peakPhi_41m",peakPhi[0]);
		inputTree_reco->SetBranchAddress("peakPhi_300m",peakPhi[1]);

		int recoBinSelect = 19; //300 m map
		int recoBinCalpulser = 6; //41 m map

		runNum_out=runNum; //get the run number in here formally finally

		// now to open the CW file
		// which should work just the same as it did before (yay)
		char summary_file_name[400];
		if(isSimulation){
			if(year_or_energy<25)
				sprintf(summary_file_name,"%s/CWID/A%d/c%d/E%2.2f/CWID_station_%d_run_%d.root",SimDirPath,station,config,year_or_energy,station,runNum);
			else
				sprintf(summary_file_name,"%s/CWID/A%d/c%d/E%d/CWID_station_%d_run_%d.root",SimDirPath,station,config,int(year_or_energy),station,runNum);
		}
		else{
			sprintf(summary_file_name,"%s/CWID/A%d/all_runs/CWID_station_%d_run_%d.root",cw_directory.c_str(),station,station,runNum);
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

		int numEntries = inputTree_filter->GetEntries();
		Long64_t starEvery=numEntries/200;
		if(starEvery==0) starEvery++;
		cout<<"Num entries is "<<numEntries<<endl;
		cout<<"Star every is "<<starEvery<<endl;
		//numEntries=100;

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
			bool isSurf[2]={false};
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

			// probably should check the WFRMS filter here so that we don't try to do things with garbage events
			// FIX ME!!!!
			// bool didYouFixEarlyWFRMSCheck=false;
			// if(didYouFixEarlyWFRMSCheck==false){
			// 	cout<<"Brian! You need to fix the check for the WFRMS earlier in the code!"<<endl;
			// 	return -1;
			// }

			// our peak finding in the 100pct code is "simpler" than in 10pct, but still, be careful...

			// get the event
			inputTree_reco->GetEntry(event);

			//figure out which reconstruction map (vpol or hpol) is best
			//in the present analysis, this only matters for the 300m bin
			double bestCorr[] = {0., 0., 0.};
			int bestCorrRadiusBin[3];
			int bestPol = 2;
			int bestTheta[3];
			int bestPhi[3];

			// "bestTheta" is the 300m radius, so peakTheta[1][x], peakPhi[1][x], etc.
			// where x=0 for VPol and x=1 for Hpol
			// we never use the third element of bestTheta or bestPhi or bestCorr, so let's just leave those unfilled

			bestTheta[0] = peakTheta[1][0];
			bestTheta[1] = peakTheta[1][1];
			bestPhi[0] = peakPhi[1][0];
			bestPhi[1] = peakPhi[1][1];
			bestCorr[0] = peakCorr[1][0];
			bestCorr[1] = peakCorr[1][1];

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

			// "bestTheta_pulser" is the 300m radius, so peakTheta[0][x], peakPhi[0][x], etc.
			// where x=0 for VPol and x=1 for Hpol
			// we never use the third element of bestTheta_pulser or bestPhi_pulser or bestCorr_pulser, so let's just leave those unfilled
			bestTheta_pulser[0] = peakTheta[0][0];
			bestTheta_pulser[1] = peakTheta[0][1];
			bestPhi_pulser[0] = peakPhi[0][0];
			bestPhi_pulser[1] = peakPhi[0][1];
			bestCorr_pulser[0] = peakCorr[0][0];
			bestCorr_pulser[1] = peakCorr[0][1];
			
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
				if(badFreqListLocal_baseline.size()>0) isCutonCW_baseline[pol]=true;
			}

			double threshCW=10;
			if(station==2){
				threshCW = 1.5;
			}
			else if(station==3){
				threshCW = 2.0;
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

			for(int pol=0; pol<2; pol++){
				corr_val_org[pol]=bestCorr[pol];
				snr_val_org[pol]=SNRs[pol];

				// copy these over
				corr_val_new[pol] = corr_val_org[pol];
				snr_val_new[pol] = snr_val_org[pol];

				// printf(BLUE"Run %4d, Event %5d/%5d, Pol %d : isCal %d, isSoft %d, isShort %d, does Fail WFRMS %d, isCP5 %d, isCP6 %d, isSurf %d and %d, isBad %d, isFirstFive %d \n"RESET
				// 	,runNum, event, eventNumber_out, pol, isCalPulser, isSoftTrigger, isShort, failWavefrontRMS[pol], isCP5, isCP6, isSurf[0], isSurf[1], isBadEvent, isFirstFiveEvent);
				

				// this is a sad and desperate mix of "input" and "output variables"
				// that is objectively a bad idea; I should really clean this up.

				// printf("Status: %d, %d, %d, %d, %d, %d, %d, %d \n"
				// 				,isCalPulser_in
				// 				,isSoftTrigger_in
				// 				,isShort
				// 				,failWavefrontRMS[pol]
				// 				,isCP5
				// 				,isCP6
				// 				,isBadEvent_out
				// 				,isFirstFiveEvent_out
				// 				 );

				// if(!failWavefrontRMS[pol])
				// 	cout<<"Event "<<event<<" doesn't fail the WFRMS filter in pol "<<pol<<" with a value of "<<TMath::Log10(bestFaceRMS[pol])<<endl;

				if(!isCalPulser_in
					&& !isSoftTrigger_in
					&& !isShort
					&& !failWavefrontRMS[pol]
					&& !isCP5 && !isCP6
					&& !isBadEvent_out
					&& !isFirstFiveEvent_out
					// now, we will let all events which are surface be filtered too. this is a  major change.
					// && !isSurf[0] && !isSurf[1] //check both pols for surface
				){ //cut cal pulsers

					// printf(RED"Time to do math on event %d \n"RESET,event);

					// load in the data for the event
					char run_file_name[400];
					if(isSimulation){
						if(year_or_energy<25)
							sprintf(run_file_name,"%s/RawSim/A%d/c%d/E%2.1f/AraOut.A%d_c%d_E%2.1f.txt.run%d.root",SimDirPath,station,config,year_or_energy,station,config,year_or_energy,runNum);
						else{
							// should just be Kotera
							sprintf(run_file_name,"%s/RawSim/A%d/c%d/E%d/AraOut.A%d_c%d_E%d.txt.run%d.root",SimDirPath,station,config,int(year_or_energy),station,config,int(year_or_energy),runNum);
						}
					}
					else{
						sprintf(run_file_name,"%s/RawData/A%d/by_config/c%d/event%d.root",data_directory.c_str(),station,config,runNum);
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

					// printf(GREEN"     Event %d has bad channel %d and %d  \n"RESET,event,hasBadSpareChanIssue_out, hasBadSpareChanIssue2_out);

					/*
					if it's not in need of re-filtering, check the "top" reco again
					*/

					// printf("CW Status for Run %4d Event %5d, Pol %d: CW fwd %d, CW back %d, CW Baseline %d \n", runNum, event, pol, isCutonCW_fwd[pol], isCutonCW_back[pol], isCutonCW_baseline[pol]);
					if(
						(!isCutonCW_fwd[pol] && !isCutonCW_back[pol] && !isCutonCW_baseline[pol]) 
						&& !hasBadSpareChanIssue_out 
						&& !hasBadSpareChanIssue2_out)
					{

						// printf(BLUE"          Okay, so, no CW, yay\n"RESET);
						vector <int> chan_list_V;
						vector <int> chan_list_H;
						
						chan_list_V.clear();
						chan_list_V.push_back(0);
						chan_list_V.push_back(1);
						chan_list_V.push_back(2);

						chan_list_H.clear();
						chan_list_H.push_back(8);
						chan_list_H.push_back(9);
						chan_list_H.push_back(10);

						// this looks weird here because we're aiming for what channels to *include*

						if(
							!(
								dropBadChans
								&& station==3
								&& (
						 			(!isSimulation && runNum>getA3BadRunBoundary())
									||
									(isSimulation && config>2)
								)
							)
						){
							chan_list_V.push_back(3);
							chan_list_H.push_back(11);
						}

						vector<double> chan_SNRs;
						for(int i=0; i<16; i++){
							chan_SNRs.push_back(VPeakOverRMS[i]);
						}

						int solNum=0;

						TH2D *map_300m_top;
						if(pol==0){
							// map_300m_top = theCorrelators[1]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V);
							map_300m_top = theCorrelators[1]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V, chan_SNRs, solNum);
						}
						if(pol==1){
							// map_300m_top = theCorrelators[1]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H);
							map_300m_top = theCorrelators[1]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H, chan_SNRs, solNum);
						}

						bool printReco=false;
						if(printReco){
							TCanvas *cMaps = new TCanvas("","",2*1100,2*850);
							map_300m_top->Draw("colz");
							char save_temp_title[400];		
							sprintf(save_temp_title,"%s/trouble_events/%d.%d.%d_Run%d_Ev%d_Maps_InsideSaveVals.png",plotPath,year_now,month_now,day_now,runNum,event);
							cMaps->SaveAs(save_temp_title);
							delete cMaps;
						}
						int PeakTheta_Recompute_300m_top;
						int PeakPhi_Recompute_300m_top;
						double PeakCorr_Recompute_300m_top;
						getCorrMapPeak(map_300m_top,PeakTheta_Recompute_300m_top,PeakPhi_Recompute_300m_top,PeakCorr_Recompute_300m_top);
						// cout<<"                    New theta peak in map is "<<PeakTheta_Recompute_300m_top<<endl;

						if(PeakTheta_Recompute_300m_top>=37){
							// cout<<"                         Ah-ha! This is a top surface event!"<<endl;
							isSurfEvent_top[pol]=1;
						}

						delete map_300m_top;
					}

					// and now to do *filtering*
					if(
						(isCutonCW_fwd[pol] || isCutonCW_back[pol] || isCutonCW_baseline[pol]) 
						&& !hasBadSpareChanIssue_out 
						&& !hasBadSpareChanIssue2_out
					)
					{
						isCW_out=1;
						Refilt[pol]=1;

						// printf(RED"	Need to filter event %d in pol %d \n"RESET,event,pol);

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

						//identify the unique center frequencies and the bandwidths around them
						vector <double> uniqueNotchFreqs;
						vector <double> uniqueNotchBands;
						for(int iFreq=0; iFreq<mergedFreqList.size(); iFreq++){
							// cout<<"Frequency "<<iFreq<<" to be notched is "<<mergedFreqList[iFreq]<<endl;
						}

						theCorrelators[0]->pickFreqsAndBands(mergedFreqList,uniqueNotchFreqs,uniqueNotchBands);
						for (int i = 0; i < uniqueNotchFreqs.size(); ++i)
						{
							// printf("Unique freq to be notched is %.2f with width %.2f \n", uniqueNotchFreqs[i],uniqueNotchBands[i]);
						}

						// for(int iFreq=0; iFreq<uniqueNotchFreqs.size(); iFreq++)
						// 	printf("				Unique freq %d is %.2f with band %.2f\n",iFreq,uniqueNotchFreqs[iFreq],uniqueNotchBands[iFreq]);

						/*
							First we must re-do the SNR calculation (so much code...)
						*/

						stringstream ss1;
						string xLabel, yLabel;
						xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
						vector<string> titlesForGraphs;
						for (int i = 0; i < 16; i++){
							ss1.str("");
							ss << "Channel " << i;
							titlesForGraphs.push_back(ss1.str());
						}

						vector <TGraph*> grWaveformsRaw = makeGraphsFromRF(realAtriEvPtr,16,xLabel,yLabel,titlesForGraphs);
						vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(grWaveformsRaw, interpolationTimeStep, xLabel, yLabel, titlesForGraphs);
						vector<TGraph*> grWaveformsPadded = makePaddedGraphs(grWaveformsInt, 0, xLabel, yLabel, titlesForGraphs);
						vector <TGraph*> grNotched;
						for(int i=0; i<16; i++){
							TGraph *grNotchAmp = theCorrelators[0]->applyAdaptiveFilter_singleAnt_FiltMany(grWaveformsPadded[i],uniqueNotchFreqs,uniqueNotchBands);
							grNotched.push_back(theCorrelators[0]->GeometricFilter(grNotchAmp,uniqueNotchFreqs,uniqueNotchBands,uniqueNotchFreqs));
							delete grNotchAmp;
						}
						vector<TGraph*> grWaveformsPowerSpectrum = makePowerSpectrumGraphs(grWaveformsPadded, xLabel, yLabel, titlesForGraphs);
						vector<TGraph*> grWaveformsPowerSpectrum_notched = makePowerSpectrumGraphs(grNotched, xLabel, yLabel, titlesForGraphs);

						for(int i=0; i<16; i++){
							double original_power=0.;
							double final_power=0.;
							for(int samp=0; samp<grWaveformsPowerSpectrum[i]->GetN(); samp++){
								original_power+=grWaveformsPowerSpectrum[i]->GetY()[samp];
							}
							for(int samp=0; samp<grWaveformsPowerSpectrum_notched[i]->GetN(); samp++){
								final_power+=grWaveformsPowerSpectrum_notched[i]->GetY()[samp];
							}
							if(i<8)
								frac_of_power_notched_V[i]=(original_power-final_power)/original_power;
							else
								frac_of_power_notched_H[i-8]=(original_power-final_power)/original_power;
						}

						//okay, now to do the filtering
						vector<double> vVPeak_new;
						vector<double> vVPeakTimes_new;
						double VPeak_new[16];
						getAbsMaximum(grNotched, vVPeakTimes_new, vVPeak_new);
						copy(vVPeak_new.begin(), vVPeak_new.begin()+16, VPeak_new); //copy these results into the arrays, because ROOT prefers arrays
						double VPeakTimes_new[16];
						copy(vVPeakTimes_new.begin(), vVPeakTimes_new.begin()+16, VPeakTimes_new); //copy these results into the arrays, because ROOT prefers arrays
						vector<double> vWaveformRMS_new;
						vWaveformRMS_new.resize(16);
						
						char run_summary_name[400];
						if (isSimulation == false){
							sprintf(run_summary_name,"%s/RunSummary/A%d/by_config/c%d/run_summary_station_%d_run_%d.root",runsummary_directory.c_str(),station,config,station,runNum);
						}
						else {
							if(station==2){
								sprintf(run_summary_name,"/fs/scratch/PAS0654/ara/sim/RunSummary/run_summary_station_2_run_20.root");
							}
							else if(station==3){
								sprintf(run_summary_name,"/fs/scratch/PAS0654/ara/sim/RunSummary/run_summary_station_3_run_30.root");
							}
						}

						TFile *summaryFile = TFile::Open(run_summary_name);
						if(!summaryFile){
							cout<<"Can't open summary file!"<<endl;
							return -1;
						}
						TTree *SummaryTree = (TTree*) summaryFile-> Get("SummaryTree");
						if(!SummaryTree){
							cout<<"Can't find summaryTree for map"<<endl;
							return -1;
						}
						
						double RMS_SoftTrigger[16];
						double RMS_RFTrigger[16];
						SummaryTree->SetBranchAddress("RMS_RFTrigger", &RMS_RFTrigger);
						SummaryTree->SetBranchAddress("RMS_SoftTrigger", &RMS_SoftTrigger);
						SummaryTree->GetEntry(0);

						int nGraphs=16;
						for (int i = 0; i < nGraphs; i++){ //loop over graphs
							//the RMS_SoftTrigger comes out of the run summary
							//so, what we want to do is see if the RMS of the software triggers was computed successfully
							if (RMS_SoftTrigger[i] == RMS_SoftTrigger[i]){ //check to make sure it's not a nan
								vWaveformRMS_new[i] = RMS_SoftTrigger[i];
							} else { //if it was a nan, then instead we'll look at the RF trigger version
								if (RMS_RFTrigger[i] == RMS_RFTrigger[i]){ //make sure it's not a nan
									vWaveformRMS_new[i] = RMS_RFTrigger[i];
								}
							}
						}
						double waveformRMS_new[16];
						copy(vWaveformRMS_new.begin(), vWaveformRMS_new.begin()+16, waveformRMS_new); //copy into the array
						vector<double> vWaveformRMS_50ns_new;
						int first50ns = (int)(50./interpolationTimeStep);
						getRMS(grNotched, vWaveformRMS_50ns_new, first50ns);
						double waveformRMS_50ns_new[16];
						copy(vWaveformRMS_50ns_new.begin(), vWaveformRMS_50ns_new.begin()+16, waveformRMS_50ns_new); //copy those results into an array
						vector<double> vVPeakOverRMS_new;
						vVPeakOverRMS_new.resize(16);
						for (int i = 0; i < 16; i++){
							vVPeakOverRMS_new[i] = VPeak_new[i]/waveformRMS_new[i];
							vVPeakOverRMS_new[i] = VPeak_new[i]/waveformRMS_new[i];
						}
						AraGeomTool * geomTool = new AraGeomTool();
						vector<int> polarizations;
						vector<int> antenna_numbers;
						polarizations.resize(16);
						antenna_numbers.resize(16);
						vector< vector <double> > ant_loc; //will be 16x3 vector of the x,y,z's the 16 antennas
						ant_loc.resize(16);
						for (int i = 0; i < 16; i++){
							ant_loc[i].resize(3);
							ant_loc[i][0] = geomTool->getStationInfo(station)->getAntennaInfo(i)->antLocation[0];
							ant_loc[i][1] = geomTool->getStationInfo(station)->getAntennaInfo(i)->antLocation[1];
							ant_loc[i][2] = geomTool->getStationInfo(station)->getAntennaInfo(i)->antLocation[2];
							polarizations[i] = (int)geomTool->getStationInfo(station)->getAntennaInfo(i)->polType;
							antenna_numbers[i]=i;
						}

						vector<int> chan_exclusion_list;
						if(dropBadChans){
							if(station==2){
								// hpol channel 15
								chan_exclusion_list.push_back(15);
							}
							else if(station==3){
								if( 
									(!isSimulation && runNum>getA3BadRunBoundary())
									||
									(isSimulation && config>2)

								){								// vpol sring 4
									chan_exclusion_list.push_back(3);
									chan_exclusion_list.push_back(7);
									chan_exclusion_list.push_back(11);
									chan_exclusion_list.push_back(15);
								}
							}
						}

						vector<double> vThirdVPeakOverRMS_new;
						double thirdVPeakOverRMS_new[3];
						getThirdVPeakOverRMS(vVPeakOverRMS_new, polarizations, antenna_numbers, chan_exclusion_list, vThirdVPeakOverRMS_new);
						for (int i = 0 ; i< 3; i++){ //pull out the first three entries
							thirdVPeakOverRMS_new[i] = vThirdVPeakOverRMS_new[i];
						}

						xLabel = "Time (ns)"; yLabel = "Integrated Power (arb units)";
						int numBinsToIntegrate = (int)(5./interpolationTimeStep);
						vector<TGraph*> grIntPower = makeIntegratedBinPowerGraphs(grNotched, numBinsToIntegrate, xLabel, yLabel, titlesForGraphs);

						vector<double> hitTimes_new; //what are the hit times
						vector<double> peakIntPowers_new; //what are the powers at those hit times?
						getAbsMaximum(grIntPower, hitTimes_new, peakIntPowers_new);
						vector<vector<double> > vvHitTimes_new; //vector of vector of hit times
						vector<vector<double> > vvPeakIntPowers_new; //vector of vector of power at the those hit times
						int numSearchPeaks = 2;
						const int numFaces = 12;
						getAbsMaximum_N(grIntPower, numSearchPeaks, 5.0, vvHitTimes_new, vvPeakIntPowers_new);
						vector<double> peakIntRMS_new;       
						for (int i = 0; i < peakIntPowers_new.size(); i++){
							peakIntRMS_new.push_back(sqrt(peakIntPowers_new[i]/numBinsToIntegrate));
						}
						double avgPeakPower_5ns_new[16];
						double peakPowerTimes_new[16];
						for (int i = 0; i < 16; i++){
							avgPeakPower_5ns_new[i] = peakIntPowers_new[i]/numBinsToIntegrate;
							peakPowerTimes_new[i] = hitTimes_new[i];
						}
						vector<double> RMS_10overRMS_new;
						for (int i = 0; i < 16; i++){
							RMS_10overRMS_new.push_back(sqrt(avgPeakPower_5ns_new[i])/waveformRMS_new[i]);
						}
						vector<vector<double> > vvRMS_10overRMS_new;
						vvRMS_10overRMS_new.resize(16);
						for (int i = 0; i < 16; i++){
							vvRMS_10overRMS_new[i].resize(vvPeakIntPowers_new[i].size());
							for (int j = 0; j < vvPeakIntPowers_new[i].size(); j++){
								vvRMS_10overRMS_new[i][j] = sqrt(vvPeakIntPowers_new[i][j]/numBinsToIntegrate)/waveformRMS_new[i];	
							}
						}

						vector< vector< int > > pairs_V_new;
						vector< vector< int > > pairs_H_new;
						setupCorrelationPairs(station, pairs_V_new, pairs_H_new); //just sets up the pairs (like, 0,1, 0,2 etc) that go into the correlation

						vector<double> bestTimes_V_new;
						vector<double> bestCorrs_V_new;
						vector<double> bestTimes_H_new;
						vector<double> bestCorrs_H_new;

						//now, to set up all the pairs that contribute to the faces
						vector<vector<vector<vector<int> > > > faces = setupFaces(station);

						//loop over the thresholds that decide if a face is allowed to contribute
						const int thresholdSteps = 15;
						double thresholdMin = 2.0;
						double thresholdStep = 0.1;
						double rms_pol_thresh_face_new_V[15][12];
						double rms_pol_thresh_face_new_H[15][12];
						for (int thresholdBin = 0; thresholdBin < thresholdSteps; thresholdBin++){
							double threshold = thresholdMin + thresholdStep*(double)thresholdBin;
							
							//get the RMS of all the faces at this threshold bin
							vector<double> rms_faces_V_new = getRms_Faces_Thresh_N(vvHitTimes_new, vvRMS_10overRMS_new, threshold, 0, faces, ant_loc);
							vector<double> rms_faces_H_new = getRms_Faces_Thresh_N(vvHitTimes_new, vvRMS_10overRMS_new, threshold, 1, faces, ant_loc);

							for (int i = 0; i < numFaces; i++){ //loop over the faces, and store the RMS for that polarization, threshold bin, and face
								rms_pol_thresh_face_new_V[thresholdBin][i] = rms_faces_V_new[i];
								rms_pol_thresh_face_new_H[thresholdBin][i] = rms_faces_H_new[i];
							}	
						} // end threshold scan

						//now the dropped channel case
						vector<vector<vector<vector<int> > > > faces_drop = setupFaces(station, dropBadChans);
						double rms_pol_thresh_face_new_V_drop[15][numFaces_new_V];
						double rms_pol_thresh_face_new_H_drop[15][numFaces_new_H];
						for (int thresholdBin = 0; thresholdBin < thresholdSteps; thresholdBin++){
							double threshold = thresholdMin + thresholdStep*(double)thresholdBin;
							
							//get the RMS of all the faces at this threshold bin
							vector<double> rms_faces_V_new_drop = getRms_Faces_Thresh_N(vvHitTimes_new, vvRMS_10overRMS_new, threshold, 0, faces_drop, ant_loc);
							vector<double> rms_faces_H_new_drop = getRms_Faces_Thresh_N(vvHitTimes_new, vvRMS_10overRMS_new, threshold, 1, faces_drop, ant_loc);

							for (int i = 0; i < numFaces_new_V; i++){
								rms_pol_thresh_face_new_V_drop[thresholdBin][i] = rms_faces_V_new_drop[i];
							}
							for (int i = 0; i < numFaces_new_H; i++){
								rms_pol_thresh_face_new_H_drop[thresholdBin][i] = rms_faces_H_new_drop[i];
							}
						} // end threshold scan

						int thresholdBin_pol_new[]={thresholdBin_pol[0],thresholdBin_pol[1]};

						vector <double> rms_faces_V_new;
						vector <double> rms_faces_H_new;

						if(dropBadChans){
							int num_faces_for_V_loop;
							int num_faces_for_H_loop;
							if(station==2){
								rms_faces_V_new.resize(numFaces);
								num_faces_for_V_loop=numFaces;
								rms_faces_H_new.resize(numFaces_A2_drop);
								num_faces_for_H_loop=numFaces_A2_drop;
							}
							else if(station==3){
								if(
									(!isSimulation && runNum>getA3BadRunBoundary())
									||
									(isSimulation && config>2)
								){ //it's 2014+, drop string four
									rms_faces_V_new.resize(numFaces_A3_drop);
									num_faces_for_V_loop=numFaces_A3_drop;
									rms_faces_H_new.resize(numFaces_A3_drop);
									num_faces_for_H_loop=numFaces_A3_drop;
								}
								else{ //it's 2013-, keep string four
									rms_faces_V_new.resize(numFaces);
									num_faces_for_V_loop=numFaces;
									rms_faces_H_new.resize(numFaces);
									num_faces_for_H_loop=numFaces;
								}
							}

							//now we loop over the faces
							for(int i=0; i<num_faces_for_V_loop; i++){
								rms_faces_V_new[i] = rms_pol_thresh_face_new_V_drop[thresholdBin_pol_new[0]][i];
							}
							for(int i=0; i<num_faces_for_H_loop; i++){
								rms_faces_H_new[i] = rms_pol_thresh_face_new_H_drop[thresholdBin_pol_new[1]][i];
							}
						}
						else{
							rms_faces_V_new.resize(12);
							rms_faces_H_new.resize(12);
							//now, we must loop over the faces
							for(int i=0; i<12; i++){
								rms_faces_V_new[i] = rms_pol_thresh_face_new_V[thresholdBin_pol_new[0]][i];  //this is right RMS for this polarization, threshold requirement, and face
								rms_faces_H_new[i] = rms_pol_thresh_face_new_H[thresholdBin_pol_new[1]][i];
							}
						}

						sort(rms_faces_V_new.begin(), rms_faces_V_new.end());
						sort(rms_faces_H_new.begin(), rms_faces_H_new.end());
						
						double bestFaceRMS_new[2];
						bestFaceRMS_new[0]=rms_faces_V_new[0];
						bestFaceRMS_new[1]=rms_faces_H_new[0];

						double SNRs_new[2];
						SNRs_new[0] = vThirdVPeakOverRMS_new[0];
						SNRs_new[1] = vThirdVPeakOverRMS_new[1];

						/*
							Now we must re-do the reconstructions
						*/
						
						vector <int> chan_list_V;
						vector <int> chan_list_H;
						for(int chan=0; chan<=7; chan++){
							chan_list_V.push_back(chan);
							chan_list_H.push_back(chan+8);
						}
						if(dropBadChans){
							if(station==2){
								//for station 2, we need to exclude channel 15 from the analysis
								chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
							}
							else if(station==3){
								//for station 3, remove for data after getA3BadRunBoundary, or for sim after config 2
								if( 
									(!isSimulation && runNum>getA3BadRunBoundary())
									||
									(isSimulation && config>2)

								){
									chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 3), chan_list_V.end());
									chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 7), chan_list_V.end());

									chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 11), chan_list_H.end());
									chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
								}
							}
						}

						// get those filtered SNR's we just computed
						vector<double> chan_SNRs;
						for(int i=0; i<16; i++){
							chan_SNRs.push_back(vVPeakOverRMS_new[i]);
						}

						TH2D *map_30m;
						TH2D *map_300m;
						int solNum = 0;
						if(pol==0){
							// map_30m = theCorrelators[0]->getInterferometricMap_RT_FiltMany_select(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V, solNum,-1,uniqueNotchFreqs,uniqueNotchBands);
							// map_300m = theCorrelators[1]->getInterferometricMap_RT_FiltMany_select(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V, solNum,-1,uniqueNotchFreqs,uniqueNotchBands);

							map_30m = theCorrelators[0]->getInterferometricMap_RT_FiltMany_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V, chan_SNRs, solNum,-1,uniqueNotchFreqs,uniqueNotchBands);
							map_300m = theCorrelators[1]->getInterferometricMap_RT_FiltMany_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V, chan_SNRs, solNum,-1,uniqueNotchFreqs,uniqueNotchBands);

						}
						if(pol==1){
							// map_30m = theCorrelators[0]->getInterferometricMap_RT_FiltMany_select(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H, solNum,-1,uniqueNotchFreqs,uniqueNotchBands);
							// map_300m = theCorrelators[1]->getInterferometricMap_RT_FiltMany_select(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H, solNum,-1,uniqueNotchFreqs,uniqueNotchBands);

							map_30m = theCorrelators[0]->getInterferometricMap_RT_FiltMany_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H, chan_SNRs, solNum,-1,uniqueNotchFreqs,uniqueNotchBands);
							map_300m = theCorrelators[1]->getInterferometricMap_RT_FiltMany_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H, chan_SNRs, solNum,-1,uniqueNotchFreqs,uniqueNotchBands);

						}
						int PeakTheta_Recompute_30m;
						int PeakTheta_Recompute_300m;
						int PeakPhi_Recompute_30m;
						int PeakPhi_Recompute_300m;
						double PeakCorr_Recompute_30m;
						double PeakCorr_Recompute_300m;
						getCorrMapPeak(map_30m,PeakTheta_Recompute_30m,PeakPhi_Recompute_30m,PeakCorr_Recompute_30m);
						getCorrMapPeak(map_300m,PeakTheta_Recompute_300m,PeakPhi_Recompute_300m,PeakCorr_Recompute_300m);

						// theta_300[pol]=PeakTheta_Recompute_300m;
						// phi_300[pol]=PeakPhi_Recompute_300m;
						// theta_41[pol]=PeakTheta_Recompute_30m;
						// phi_41[pol]=PeakPhi_Recompute_30m;

						// update--only change the _new verions
						theta_300_new[pol]=PeakTheta_Recompute_300m;
						phi_300_new[pol]=PeakPhi_Recompute_300m;
						theta_41_new[pol]=PeakTheta_Recompute_30m;
						phi_41_new[pol]=PeakPhi_Recompute_30m;						
						
						chan_list_V.clear();
						chan_list_V.push_back(0);
						chan_list_V.push_back(1);
						chan_list_V.push_back(2);

						chan_list_H.clear();
						chan_list_H.push_back(8);
						chan_list_H.push_back(9);
						chan_list_H.push_back(10);

						// this looks weird here because we're aiming for what channels to *include*

						if(
							!(
								dropBadChans
								&& station==3
								&& (
									(!isSimulation && runNum>getA3BadRunBoundary())
									||
									(isSimulation && config>2)
								)
							)
						){
							chan_list_V.push_back(3);
							chan_list_H.push_back(11);
						}

						TH2D *map_300m_top;
						if(pol==0){
							map_300m_top = theCorrelators[1]->getInterferometricMap_RT_FiltMany_select(settings, detector, realAtriEvPtr, Vpol, 0, chan_list_V, 0,-1,uniqueNotchFreqs,uniqueNotchBands);
						}
						if(pol==1){
							map_300m_top = theCorrelators[1]->getInterferometricMap_RT_FiltMany_select(settings, detector, realAtriEvPtr, Hpol, 0, chan_list_H, 0,-1,uniqueNotchFreqs,uniqueNotchBands);
						}

						int PeakTheta_Recompute_300m_top;
						int PeakPhi_Recompute_300m_top;
						double PeakCorr_Recompute_300m_top;
						getCorrMapPeak(map_300m_top,PeakTheta_Recompute_300m_top,PeakPhi_Recompute_300m_top,PeakCorr_Recompute_300m_top);

						//cleanup
						deleteGraphVector(grWaveformsRaw);
						deleteGraphVector(grWaveformsInt);
						deleteGraphVector(grWaveformsPadded);
						deleteGraphVector(grIntPower);
						deleteGraphVector(grNotched);
						deleteGraphVector(grWaveformsPowerSpectrum);
						deleteGraphVector(grWaveformsPowerSpectrum_notched);

						printf("		old vs new logrms calc in pol %d: %.2f vs %.2f \n",pol,log(bestFaceRMS[pol])/log(10),log(bestFaceRMS_new[pol])/log(10));
						printf("		old vs new snr in pol %d: %.2f vs %.2f \n",pol,SNRs[pol],SNRs_new[pol] );
						printf("		old vs new corr in pol %d: %.4f vs %.4f \n",pol,corr_val_org[pol],PeakCorr_Recompute_300m);


						// re-check top face reco
						if(PeakTheta_Recompute_300m_top>=37)
							isSurfEvent_top[pol]=1;

						// re-check surface cut
						if(PeakTheta_Recompute_300m>=37){
							// isSurfEvent[pol]=1;
							isSurfEvent_new_out[pol]=1;
						}
						else{
							// isSurfEvent[pol]=0;
							isSurfEvent_new_out[pol]=0;
						}

						// re-check wavefront RMS cut
						if(log(bestFaceRMS_new[pol])/log(10) < wavefrontRMScut[pol]){
							// WFRMS[pol]=0;
							WFRMS_new[pol]=0;
						}
						else{
							// WFRMS[pol]=1;
							WFRMS_new[pol]=1;
						}

						// assign the newly computed corr and snr values
						// corr_val[pol]=PeakCorr_Recompute_300m;
						// snr_val[pol]=SNRs_new[pol];

						corr_val_new[pol]=PeakCorr_Recompute_300m;
						snr_val_new[pol]=SNRs_new[pol];

						delete map_300m;
						delete map_300m_top;
						delete map_30m;

						// yes, the summary file close and delete needs to be down here
						// to make ROOT's silly ownership thingy work out and not cause segfault
						summaryFile->Close();
						delete summaryFile;
					} //if any frequencies are flagged for filtering
					if(!isSimulation) delete realAtriEvPtr;
					mapFile->Close();
					delete mapFile;
				} //cal pulser
				trees[pol]->Fill();
			}//loop over polarization
			trees[2]->Fill();
			trees[3]->Fill();
			trees[4]->Fill();
		}//loop over events
		inputFile->Close();
		NewCWFile->Close();
		delete inputFile;
		delete NewCWFile;

		fpOut->Write();
		fpOut->Close();
		delete fpOut;
		printf("Done! Run Number %d \n", runNum);
	} //end loop over input files
}
