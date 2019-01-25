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

using namespace std;

int PlotThisEvent(int station, int year, int runNum, int event, Settings *settings, Detector *detector, RayTraceCorrelator *theCorrelators[2]);

int main(int argc, char **argv)
{
	gStyle->SetOptStat(0);
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	stringstream ss;
	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	
	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <station> <year> <output_location> <joined filename 1> <joined filename 2 > ... <joined filename x>"<<endl;
		return 0;
	}
	int station = atoi(argv[1]);
	int year = atoi(argv[2]);
	string output_location = argv[3];

	//just to have the cut parameters up front and easy to find
	int thresholdBin_pol[]={3,5}; //bin 3 = 2.3, bin 5 = 2.5 //what is the faceRMS inclusion threshold?
	double wavefrontRMScut[]={-1.5, -1.5}; //event wavefrontRMS < this value

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
		TTree *trees[3];
		trees[0]= new TTree("VTree","VTree");
		trees[1]= new TTree("HTree","HTree");
		trees[2]= new TTree("AllTree","AllTree");

		double corr_val[2];
		double snr_val[2];
		int WFRMS[2];
		double frac_of_power_notched_V[8];
		double frac_of_power_notched_H[8];

		int Refilt[2];
		trees[0]->Branch("corr_val_V",&corr_val[0]);
		trees[0]->Branch("snr_val_V",&snr_val[0]);
		trees[0]->Branch("wfrms_val_V",&WFRMS[0]);
		trees[0]->Branch("Refilt_V",&Refilt[0]);
		trees[1]->Branch("corr_val_H",&corr_val[1]);
		trees[1]->Branch("snr_val_H",&snr_val[1]);
		trees[1]->Branch("wfrms_val_H",&WFRMS[1]);
		trees[0]->Branch("Refilt_H",&Refilt[1]);
		for(int i=0; i<8; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			trees[0]->Branch(ss.str().c_str(),&frac_of_power_notched_V[i]);
		}
		for(int i=8; i<16; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			trees[1]->Branch(ss.str().c_str(),&frac_of_power_notched_H[i]);
		}

		int isCal;
		int isSoft;
		int isShortWave;
		int isCW;
		int isNewBox;
		int isSurfEvent;

		trees[2]->Branch("cal",&isCal);
		trees[2]->Branch("soft",&isSoft);
		trees[2]->Branch("short",&isShortWave);
		trees[2]->Branch("CW",&isCW);
		trees[2]->Branch("box",&isNewBox);
		trees[2]->Branch("surf",&isSurfEvent);

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
		double thirdVPeakOverRMS[3]; //the third highest vpeak over RMS
		double rms_pol_thresh_face[2][15][12];
		bool isCalPulser;
		bool isSoftTrigger;
		int waveformLength[16];
		inputTree_filter->SetBranchAddress("thirdVPeakOverRMS", &thirdVPeakOverRMS);
		inputTree_filter->SetBranchAddress("rms_pol_thresh_face", &rms_pol_thresh_face);
		inputTree_filter->SetBranchAddress("isCalpulser",&isCalPulser);
		inputTree_filter->SetBranchAddress("isSoftTrigger",&isSoftTrigger);
		inputTree_filter->SetBranchAddress("waveformLength",&waveformLength);

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
		sprintf(summary_file_name,"/fs/scratch/PAS0654/ara/10pct/CWID/A%d/%d/CWID_station_%d_run_%d.root",station,year,station,runNum);
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

		//now to loop over events
		for(int event=0; event<numEntries; event++){

			// if(event%starEvery==0) {
			// 	std::cout << "	On event "<<event<<endl;
			// }

			isCal=0;
			isSoft=0;
			isShortWave=0;
			isCW=0;
			isNewBox=0;
			isSurfEvent=0;
			Refilt[0]=0;
			Refilt[1]=0;
			corr_val[0]=0.;
			corr_val[1]=0.;
			snr_val[0]=0.;
			snr_val[1]=0.;
			WFRMS[0]=0;
			WFRMS[1]=0;
			for(int i=0; i<8; i++){
				frac_of_power_notched_V[i]=0.;
				frac_of_power_notched_H[i]=0.;
			}

			inputTree_filter->GetEvent(event);

			bool isShort=false;
			bool isSurf=false;
			bool isCP5=false;
			bool isCP6=false;
			bool failWavefrontRMS[2];
			failWavefrontRMS[0]=false;
			failWavefrontRMS[1]=false;

			for(int i=0;i<16;i++){ if(waveformLength[i]<550) isShort=true; }

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
				if(bestTheta[pol] > 37) isSurf=true;
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
				if (bestPhi_pulser[pol] > -30 && bestPhi_pulser[pol] < -20 && bestTheta_pulser[pol] > -25 && bestTheta_pulser[pol] < -10){
					isCP5=true;
				}
				//if (bestPhi_pulser[pol] > 60 && bestPhi_pulser[pol] < 70 && bestTheta_pulser[pol] > 10 && bestTheta_pulser[pol] < 25){
				if (bestPhi_pulser[pol] > 60 && bestPhi_pulser[pol] < 70 && bestTheta_pulser[pol] > 0 && bestTheta_pulser[pol] < 15){
					isCP6=true;
				}
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

			//filter associated parameters
			double SNRs[2];
			SNRs[0] = thirdVPeakOverRMS[0];
			SNRs[1] = thirdVPeakOverRMS[1];
			if(SNRs[0]>29.) SNRs[0]=29.;
			if(SNRs[1]>29.) SNRs[1]=29.;

			vector <double>  rms_faces_V;
			rms_faces_V.resize(12);
			vector <double> rms_faces_H;
			rms_faces_H.resize(12);

			//now, we must loop over the faces
			for(int i=0; i<12; i++){
				rms_faces_V[i] = rms_pol_thresh_face[0][thresholdBin_pol[0]][i];  //this is right RMS for this polarization, threshold requirement, and face
				rms_faces_H[i] = rms_pol_thresh_face[1][thresholdBin_pol[1]][i];
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

			if(isCalPulser) isCal=1;
			if(isSoftTrigger) isSoft=1;
			if(isShort) isShortWave=1;
			if(failWavefrontRMS[0]) WFRMS[0]=1;
			if(failWavefrontRMS[1]) WFRMS[1]=1;
			if(isCP5 || isCP6 ) isNewBox=1;
			if(isSurf) isSurfEvent=1;

			for(int pol=0; pol<2; pol++){
				corr_val[pol]=bestCorr[pol];
				snr_val[pol]=SNRs[pol];
				
				if(!isCalPulser
					&& !isSoftTrigger
					&& !isShort
					&& !failWavefrontRMS[pol]
					&& !isCP5 && !isCP6
					&& !isSurf
				){ //cut cal pulsers

					//what happens if we act like it's a *cut*
					//we want all of the clean events in both histograms (as cut as after filter)
					if(!isCutonCW_fwd[pol] && !isCutonCW_back[pol] && !isCutonCW_baseline[pol]){
						PlotThisEvent(station,year,runNum,event, settings, detector, theCorrelators);
					} //not cut by any CW

					//and now to do *filtering*
					if(isCutonCW_fwd[pol] || isCutonCW_back[pol] || isCutonCW_baseline[pol]){
						isCW=1;
						Refilt[pol]=1;

						cout<<"			Need to filter event "<<event<<endl;

						char run_file_name[400];
						if(year==2013){
							sprintf(run_file_name,"/fs/scratch/PAS0654/ara/10pct/RawData/A%d/%d/run%d/event%d.root",station,year,runNum,runNum);
						}
						else if(year==2014 || year==2015 || year==2016){
							sprintf(run_file_name,"/fs/scratch/PAS0654/ara/10pct/RawData/A%d/%d/sym_links/event00%d.root",station,year,runNum,runNum);
						}
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

						if(year==2013){
							sprintf(ped_file_name,"/fs/scratch/PAS0654/ara/peds/run_specific_peds/A%d/%d/event%d_specificPeds.dat",station,year,runNum);
						}
						else if(year==2014 || year==2015 || year==2016){
							sprintf(ped_file_name,"/fs/scratch/PAS0654/ara/peds/run_specific_peds/A%d/%d/event00%d_specificPeds.dat",station,year,runNum);
						}
						calibrator->setAtriPedFile(ped_file_name,stationID); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist
						
						UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawPtr,AraCalType::kLatestCalib);

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

						// for(int iFreq=0; iFreq<uniqueNotchFreqs.size(); iFreq++) printf("				Unique freq %d is %.2f with band %.2f\n",iFreq,uniqueNotchFreqs[iFreq],uniqueNotchBands[iFreq]);

						//now, we must re-do the interferometry
						TH2D *map_30m;
						TH2D *map_300m;
						if(pol==0){
							map_30m = theCorrelators[0]->getInterferometricMap_RT_FiltMany(settings, detector, realAtriEvPtr, Vpol, 0, 0,-1,uniqueNotchFreqs,uniqueNotchBands);
							map_300m = theCorrelators[1]->getInterferometricMap_RT_FiltMany(settings, detector, realAtriEvPtr, Vpol, 0, 0,-1,uniqueNotchFreqs,uniqueNotchBands);
						}
						if(pol==1){
							map_30m = theCorrelators[0]->getInterferometricMap_RT_FiltMany(settings, detector, realAtriEvPtr, Hpol, 0, 0,-1,uniqueNotchFreqs,uniqueNotchBands);
							map_300m = theCorrelators[1]->getInterferometricMap_RT_FiltMany(settings, detector, realAtriEvPtr, Hpol, 0, 0,-1,uniqueNotchFreqs,uniqueNotchBands);
						}
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
						getCorrMapPeak_wStats(map_30m,PeakTheta_Recompute_30m,PeakPhi_Recompute_30m,PeakCorr_Recompute_30m,MinCorr_Recompute_30m,MeanCorr_Recompute_30m,RMSCorr_Recompute_30m,PeakSigma_Recompute_30m);
						getCorrMapPeak_wStats(map_300m,PeakTheta_Recompute_300m,PeakPhi_Recompute_300m,PeakCorr_Recompute_300m,MinCorr_Recompute_300m,MeanCorr_Recompute_300m,RMSCorr_Recompute_300m,PeakSigma_Recompute_300m);

						//and we must also redo the WRMS calculation
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

						int istart, istop;
						if(pol==0){ istart=0; istop=8; }
						else if(pol==1){ istart=8; istop=16; }

						for(int i=istart; i<istop; i++){
							double original_power=0.;
							double final_power=0.;
							for(int samp=0; samp<grWaveformsPowerSpectrum[i]->GetN(); samp++){
								original_power+=grWaveformsPowerSpectrum[i]->GetY()[samp];
							}
							for(int samp=0; samp<grWaveformsPowerSpectrum_notched[i]->GetN(); samp++){
								final_power+=grWaveformsPowerSpectrum_notched[i]->GetY()[samp];
							}
							if(pol==0)
								frac_of_power_notched_V[i]=(original_power-final_power)/original_power;
							else if(pol==1)
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

						//get run summary information
						
						char run_summary_name[400];
						sprintf(run_summary_name,"/fs/scratch/PAS0654/ara/10pct/RunSummary/A%d/%d/run_summary_station_%d_run_%d.root",station,year,station,runNum);
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
						polarizations.resize(16);
						vector< vector <double> > ant_loc; //will be 16x3 vector of the x,y,z's the 16 antennas
						ant_loc.resize(16);
						for (int i = 0; i < 16; i++){
							ant_loc[i].resize(3);
							ant_loc[i][0] = geomTool->getStationInfo(stationID)->getAntennaInfo(i)->antLocation[0];
							ant_loc[i][1] = geomTool->getStationInfo(stationID)->getAntennaInfo(i)->antLocation[1];
							ant_loc[i][2] = geomTool->getStationInfo(stationID)->getAntennaInfo(i)->antLocation[2];
							polarizations[i] = (int)geomTool->getStationInfo(stationID)->getAntennaInfo(i)->polType;
						}

						vector<double> vThirdVPeakOverRMS_new;
						double thirdVPeakOverRMS_new[3];
						getThirdVPeakOverRMS(vVPeakOverRMS_new, polarizations, vThirdVPeakOverRMS_new);
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
						setupCorrelationPairs(2, pairs_V_new, pairs_H_new); //just sets up the pairs (like, 0,1, 0,2 etc) that go into the correlation

						vector<double> bestTimes_V_new;
						vector<double> bestCorrs_V_new;
						vector<double> bestTimes_H_new;
						vector<double> bestCorrs_H_new;

						//now, to set up all the pairs that contribute to the faces
						vector<vector<vector<vector<int> > > > faces = setupFaces(stationID);

						//loop over the thresholds that decide if a face is allowed to contribute
						const int thresholdSteps = 15;
						double thresholdMin = 2.0;
						double thresholdStep = 0.1;
						double rms_pol_thresh_face_new[2][15][12];
						for (int thresholdBin = 0; thresholdBin < thresholdSteps; thresholdBin++){
							double threshold = thresholdMin + thresholdStep*(double)thresholdBin;
							
							//get the RMS of all the faces at this threshold bin
							vector<double> rms_faces_V_new = getRms_Faces_Thresh_N(vvHitTimes_new, vvRMS_10overRMS_new, threshold, 0, faces, ant_loc);
							vector<double> rms_faces_H_new = getRms_Faces_Thresh_N(vvHitTimes_new, vvRMS_10overRMS_new, threshold, 1, faces, ant_loc);


							for (int i = 0; i < numFaces; i++){ //loop over the faces, and store the RMS for that polarization, threshold bin, and face
								rms_pol_thresh_face_new[0][thresholdBin][i] = rms_faces_V_new[i];
								rms_pol_thresh_face_new[1][thresholdBin][i] = rms_faces_H_new[i];
							}
									
								//sort, in increasing order, the RMSs for this threshold
							sort(rms_faces_V_new.begin(), rms_faces_V_new.end());
							sort(rms_faces_H_new.begin(), rms_faces_H_new.end());
						} // end threshold scan

						int thresholdBin_pol_new[]={3,5};

						vector <double>  rms_faces_V_new;
						rms_faces_V_new.resize(12);
						vector <double> rms_faces_H_new;
						rms_faces_H_new.resize(12);
						//now, we must loop over the faces
						for(int i=0; i<12; i++){
							rms_faces_V_new[i] = rms_pol_thresh_face_new[0][thresholdBin_pol_new[0]][i];  //this is right RMS for this polarization, threshold requirement, and face
							rms_faces_H_new[i] = rms_pol_thresh_face_new[1][thresholdBin_pol_new[1]][i];
						}
						sort(rms_faces_V_new.begin(), rms_faces_V_new.end());
						sort(rms_faces_H_new.begin(), rms_faces_H_new.end());
						double bestFaceRMS_new[2];
						bestFaceRMS_new[0]=rms_faces_V_new[0];
						bestFaceRMS_new[1]=rms_faces_H_new[0];

						//cout<<"New vpol log rms value "<<log(bestFaceRMS_new[0])/log(10)<<endl;
						//cout<<"New third highest vpeak over rms "<<vThirdVPeakOverRMS_new[0]<<endl;

						double SNRs_new[2];
						SNRs_new[0] = vThirdVPeakOverRMS_new[0];
						SNRs_new[1] = vThirdVPeakOverRMS_new[1];
						if(SNRs_new[0]>29.) SNRs_new[0]=29.;
						if(SNRs_new[1]>29.) SNRs_new[1]=29.;

						//cleanup
						for(int i=0; i<16; i++){
							delete grWaveformsRaw[i];
							delete grWaveformsInt[i];
							delete grWaveformsPadded[i];
							delete grIntPower[i];
							delete grNotched[i];
							grWaveformsPowerSpectrum[i];
							grWaveformsPowerSpectrum_notched[i];
						}
						summaryFile->Close();
						delete summaryFile;
						// printf("				old vs new logrms calc in pol %d: %.2f vs %.2f \n",pol,log(bestFaceRMS[pol])/log(10),log(bestFaceRMS_new[pol])/log(10));
						// printf("				old vs new snr in pol %d: %.2f vs %.2f \n",pol,SNRs[pol],SNRs_new[pol] );


						isSurfEvent=1; //assume again it's surface
						if(PeakTheta_Recompute_300m<=37){ //recheck for surface
							isSurfEvent=0;  //mark it as not a surface event

							//recheck wrms and use the recomputed SNR
							WFRMS[pol]=1; //assume it will fail
							if(log(bestFaceRMS_new[pol])/log(10) < wavefrontRMScut[pol]){ //recheck if it *passes* the WRMS cut
								// cout<<"new wavefront RMS is "<<log(bestFaceRMS_new[pol])/log(10)<<endl;
								WFRMS[pol]=0; //actually, it passes!

								//save this out for use in optimization later
								corr_val[pol]=PeakCorr_Recompute_300m;
								snr_val[pol]=SNRs_new[pol];

								//printf("Still passes log rms cut of %.2f \n", wavefrontRMScut[pol]);
								// if(PeakCorr_Recompute_300m>=0.14){
								// 	stringstream ss1;
								// 	string xLabel, yLabel;
								// 	xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
								// 	vector<string> titlesForGraphs;
								// 	for (int i = 0; i < 16; i++){
								// 		ss1.str("");
								// 		ss << "Channel " << i;
								// 		titlesForGraphs.push_back(ss1.str());
								// 	}
								// 	vector <TGraph*> waveforms = makeGraphsFromRF(realAtriEvPtr,16,xLabel,yLabel,titlesForGraphs);
								// 	vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(waveforms, interpolationTimeStep, xLabel, yLabel, titlesForGraphs);
								// 	vector<TGraph*> grWaveformsPadded = makePaddedGraphs(grWaveformsInt, 0, xLabel, yLabel, titlesForGraphs);
								// 	xLabel = "Frequency (Hz)"; yLabel = "Power Spectral Density (mV/Hz)";
								// 	vector<TGraph*> grWaveformsPowerSpectrum = makePowerSpectrumGraphs(grWaveformsPadded, xLabel, yLabel, titlesForGraphs);

								// 	char save_temp_title[300];
								// 	sprintf(save_temp_title,"/home/brianclark/results/A23/%d.%d.%d_YesCWIssues_Waveforms_Run%d_Ev%d_pol%d.png",year_now,month_now,day_now,runNum,event,pol);
								// 	TCanvas *cWave = new TCanvas("","",4*1100,4*850);
								// 	cWave->Divide(4,4);
								// 	for(int i=0; i<16; i++){
								// 		cWave->cd(i+1);
								// 		waveforms[i]->Draw("AL");
								// 		waveforms[i]->SetLineWidth(3);
								// 	}
								// 	cWave->SaveAs(save_temp_title);
								// 	delete cWave;

								// 	TCanvas *cMaps = new TCanvas("","",2*1100,850);
								// 	cMaps->Divide(2,1);
								// 		cMaps->cd(1);
								// 		map_30m->Draw("colz");
								// 		cMaps->cd(2);
								// 		map_300m->Draw("colz");
								// 	sprintf(save_temp_title,"/home/brianclark/results/A23/%d.%d.%d_YesCWIssues_Maps_Run%d_Ev%d_pol%d.png",year_now,month_now,day_now,runNum,event,pol);
								// 	cMaps->SaveAs(save_temp_title);
								// 	delete cMaps;

								// 	sprintf(save_temp_title,"/home/brianclark/results/A23/%d.%d.%d_YesCWIssues_Spectra_Run%d_Ev%d_pol%d.png",year_now,month_now,day_now,runNum,event,pol);
								// 	TCanvas *cSpec = new TCanvas("","",4*1100,4*850);
								// 	cSpec->Divide(4,4);
								// 	for(int i=0; i<16; i++){
								// 		cSpec->cd(i+1);
								// 		grWaveformsPowerSpectrum[i]->Draw("AL");
								// 		grWaveformsPowerSpectrum[i]->SetLineWidth(3);
								// 		gPad->SetLogy();
								// 	}
								// 	cSpec->SaveAs(save_temp_title);
								// 	delete cSpec;
								// 	for(int i=0; i<16; i++){
								// 		delete waveforms[i];
								// 		delete grWaveformsInt[i];
								// 		delete grWaveformsPadded[i];
								// 		delete grWaveformsPowerSpectrum[i];
								// 	}
								// } //check if it's still in the signal box for some reason
							} //WFRMS cut on new event
						} //recheck the surface cut
						delete map_300m;
						delete map_30m;
						delete realAtriEvPtr;
						mapFile->Close();
						delete mapFile;
					} //if any frequencies are flagged for filtering
				} //cal pulser
				trees[pol]->Fill();
			}//loop over polarization
			trees[2]->Fill();
		}//loop over events
		inputFile->Close();
		NewCWFile->Close();
		delete inputFile;

		fpOut->Write();
		fpOut->Close();
		delete fpOut;
		printf("Done! Run Number %d", runNum);
	} //end loop over input files
}


int PlotThisEvent(int station, int year, int runNum, int event, Settings *settings, Detector *detector, RayTraceCorrelator *theCorrelators[2]){
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	char run_file_name[400];
	if(year==2013){
		sprintf(run_file_name,"/fs/scratch/PAS0654/ara/10pct/RawData/A%d/%d/run%d/event%d.root",station,year,runNum,runNum);
	}
	else if(year==2014 || year==2015 || year==2016){
		sprintf(run_file_name,"/fs/scratch/PAS0654/ara/10pct/RawData/A%d/%d/sym_links/event00%d.root",station,year,runNum,runNum);
	}
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

	if(year==2013){
		sprintf(ped_file_name,"/fs/scratch/PAS0654/ara/peds/run_specific_peds/A%d/%d/event%d_specificPeds.dat",station,year,runNum);
	}
	else if(year==2014 || year==2015 || year==2016){
		sprintf(ped_file_name,"/fs/scratch/PAS0654/ara/peds/run_specific_peds/A%d/%d/event00%d_specificPeds.dat",station,year,runNum);
	}
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

	bool do_reco=true;
	if(do_reco){
		TH2D *map_30m_V;
		TH2D *map_300m_V;
		TH2D *map_30m_H;
		TH2D *map_300m_H;
		TH2D *map_30m_V_select;

		map_30m_V = theCorrelators[0]->getInterferometricMap_RT_Rezero(settings, detector, realAtriEvPtr, Vpol, 0, 0,-1);
		map_300m_V = theCorrelators[1]->getInterferometricMap_RT_Rezero(settings, detector, realAtriEvPtr, Vpol, 0, 0,-1);
		map_30m_H = theCorrelators[0]->getInterferometricMap_RT_Rezero(settings, detector, realAtriEvPtr, Hpol, 0, 0,-1);
		map_300m_H = theCorrelators[1]->getInterferometricMap_RT_Rezero(settings, detector, realAtriEvPtr, Hpol, 0, 0,-1);

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
		// chan_list.push_back(5);
		// chan_list.push_back(6);
		// chan_list.push_back(7);
		// map_30m_V_select = theCorrelators[0]->getInterferometricMap_RT_select(settings,detector,realAtriEvPtr,Vpol,0,chan_list,0);


		TCanvas *cMaps = new TCanvas("","",2*1100,2*850);
		cMaps->Divide(2,2);
			cMaps->cd(1);
			map_30m_V->Draw("colz");
			cMaps->cd(2);
			map_30m_H->Draw("colz");
			cMaps->cd(3);
			map_300m_V->Draw("colz");
			cMaps->cd(4);
			map_300m_H->Draw("colz");
			// cMaps->cd(5);
			// map_30m_V_select->Draw("colz");
		char save_temp_title[400];		
		sprintf(save_temp_title,"/users/PAS0654/osu0673/A23_analysis/results/trouble_events/%d.%d.%d_Run%d_Ev%d_Maps_.png",year_now,month_now,day_now,runNum,event);
		cMaps->SaveAs(save_temp_title);
		delete cMaps;
		delete map_30m_V; delete map_300m_V; delete map_30m_H; delete map_300m_H; 
		// delete map_30m_V_select;
	}


	char save_temp_title[300];
	sprintf(save_temp_title,"/users/PAS0654/osu0673/A23_analysis/results/trouble_events/%d.%d.%d_Run%d_Ev%d_Waveforms.png",year_now,month_now,day_now,runNum,event);
	TCanvas *cWave = new TCanvas("","",4*1100,4*850);
	cWave->Divide(4,4);
	for(int i=0; i<16; i++){
		cWave->cd(i+1);
		waveforms[i]->Draw("AL");
		waveforms[i]->SetLineWidth(3);
	}
	cWave->SaveAs(save_temp_title);
	delete cWave;

	sprintf(save_temp_title,"/users/PAS0654/osu0673/A23_analysis/results/trouble_events/%d.%d.%d_Run%d_Ev%d_Spectra.png",year_now,month_now,day_now,runNum,event);
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