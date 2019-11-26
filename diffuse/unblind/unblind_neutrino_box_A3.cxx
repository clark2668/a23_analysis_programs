////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	unblind_neutrino_box.cxx
////	unblinding the neutrino signal region!
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

//AraRoot includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraQualCuts.h"

// analysis custom
// #include "tools_Cuts.h"
#include "tools_Stats.h"
#include "tools_CommandLine.h"
#include "tools_outputObjects.h"
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_CW.h"

// for ray trace correlator
#include "Settings.h"
#include "Event.h"
#include "Detector.h"
#include "Report.h"
#include "RayTraceCorrelator.h"

AraAntPol::AraAntPol_t Vpol = AraAntPol::kVertical;
AraAntPol::AraAntPol_t Hpol = AraAntPol::kHorizontal;

using namespace std;

int PlotThisEvent(int station, int runNum, int event, int problempol);

int main(int argc, char **argv)
{

	time_t time_now = time(0); //get the time now
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	stringstream ss;
	gStyle->SetOptStat(0);

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;


	if(argc<5){
		cout<< "Usage\n" << argv[0] << " <1-station> <2-config> <100 or 10> <isOrg> <output_location>"<<endl;;
		return -1;
	}
	int station = atoi(argv[1]);
	int config = atoi(argv[2]);
	int full_or_partial = atoi(argv[3]);
	int isOrg = atoi(argv[4]);
	string outputLocation = argv[5];

	if(station!=2 && station!=3){
		printf("No good! You asked for station %d, but this code only works for stations 2 and 3 \n",station);
		return -1;
	}

	vector<int> BadRunList=BuildBadRunList(station);
	vector<int> BadSurfaceRunList=BuildSurfaceRunList(station);

	int numTotal=0;

	// char outputFileName[500];
	// sprintf(outputFileName,"%s/neutrino_box_2d_cut_values_A%d_c%d_sample%d.root",outputLocation.c_str(),station,config,full_or_partial);
	// TFile *outFile = TFile::Open(outputFileName,"RECREATE");
	// TTree *outTree = new TTree("outTree","outTree");
	// int hist_this_pol[2];
	// double corr_val_out[2];
	// double snr_val_out[2];
	// int runNum_out;
	// int eventNum_out;
	// outTree->Branch("passes_this_pol_V",&hist_this_pol[0]);
	// outTree->Branch("passes_this_pol_H",&hist_this_pol[1]);
	// outTree->Branch("corr_val_V",&corr_val_out[0]);
	// outTree->Branch("corr_val_H",&corr_val_out[1]);
	// outTree->Branch("snr_val_V",&snr_val_out[0]);
	// outTree->Branch("snr_val_H",&snr_val_out[1]);
	// outTree->Branch("runNum_out",&runNum_out);
	// outTree->Branch("eventNum_out",&eventNum_out);

	TChain dataVTree("VTree");
	TChain dataHTree("HTree");
	TChain dataAllTree("AllTree");
	TChain dataRecoTree("OutputTreeReco");
	TChain dataFilterTree("OutputTree");
	char the_data[500];

	if(full_or_partial==100){
		// use the 100pct sample
    if (config==1) sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_CWThresh2.0_snrbins_0_1_wfrmsvals_-1.2_-1.3_run_*.root",station,config);
  	if (config==2) sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_CWThresh2.0_snrbins_0_1_wfrmsvals_-1.3_-1.4_run_*.root",station,config);
  	if (config==3 || config==4) sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_CWThresh2.0_snrbins_0_1_wfrmsvals_-1.0_-1.1_run_*.root",station,config);
  	if (config==5) sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_CWThresh2.0_snrbins_0_1_wfrmsvals_-0.7_-0.8_run_*.root",station,config);
	}
	if(full_or_partial==10 && isOrg==1){
		// use the *original* 10pct sample from optimization
		// sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_*.root",station,config);
		sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_*.root",station,config);
	}
	if(full_or_partial==10 && isOrg==0){
		// use the *new* 10pct sample from the slightly revised code base
		// sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/10pct_verify_try2/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_*.root",station,config);
	}

	printf("The data: %s\n", the_data);
	dataVTree.Add(the_data);
	dataHTree.Add(the_data);
	dataAllTree.Add(the_data);
	dataFilterTree.Add(the_data);

	int numDataEvents = dataVTree.GetEntries();
	// numDataEvents=100;

	// do this inside brackets for scoping power and re-use of identical variable names when it comes time for simulation to happen
	{

		double corr_val[2];
		double snr_val[2];
		int WFRMS[2];
		int theta_300[2];
		int phi_300[2];
		int theta_41[2];
		int phi_41[2];

		int Refilt[2];
		dataVTree.SetBranchAddress("Refilt_V",&Refilt[0]);
		dataHTree.SetBranchAddress("Refilt_H",&Refilt[1]);

		dataVTree.SetBranchAddress("corr_val_V_new",&corr_val[0]);
		dataVTree.SetBranchAddress("snr_val_V_new",&snr_val[0]);
		dataVTree.SetBranchAddress("wfrms_val_V_new",&WFRMS[0]);
		dataVTree.SetBranchAddress("theta_300_V_new",&theta_300[0]);
		dataVTree.SetBranchAddress("theta_41_V_new",&theta_41[0]);
		dataVTree.SetBranchAddress("phi_300_V_new",&phi_300[0]);
		dataVTree.SetBranchAddress("phi_41_V_new",&phi_41[0]);

		dataHTree.SetBranchAddress("corr_val_H_new",&corr_val[1]);
		dataHTree.SetBranchAddress("snr_val_H_new",&snr_val[1]);
		dataHTree.SetBranchAddress("wfrms_val_H_new",&WFRMS[1]);
		dataHTree.SetBranchAddress("theta_300_H_new",&theta_300[1]);
		dataHTree.SetBranchAddress("theta_41_H_new",&theta_41[1]);
		dataHTree.SetBranchAddress("phi_300_H_new",&phi_300[1]);
		dataHTree.SetBranchAddress("phi_41_H_new",&phi_41[1]);

		int isCal;
		int isSoft;
		int isShort;
		int isCW;
		int isNewBox;

		dataAllTree.SetBranchAddress("cal",&isCal);
		dataAllTree.SetBranchAddress("soft",&isSoft);
		dataAllTree.SetBranchAddress("short",&isShort);
		dataAllTree.SetBranchAddress("CW",&isCW);
		dataAllTree.SetBranchAddress("box",&isNewBox);

		int isSurf[2]; // a surface event after filtering?
		int isSurfEvent_top[2]; // a top event?

		dataAllTree.SetBranchAddress("surf_V_new2",&isSurf[0]);
		dataAllTree.SetBranchAddress("surf_H_new2",&isSurf[1]);

		dataAllTree.SetBranchAddress("surf_top_V",&isSurfEvent_top[0]);
		dataAllTree.SetBranchAddress("surf_top_H",&isSurfEvent_top[1]);

		int isBadEvent;
		double weight;
		int unixTime;
		int isFirstFiveEvent;
		int hasBadSpareChanIssue;
		int hasBadSpareChanIssue2;
		int runNum;
		int eventNumber;
		bool isSpikey;
    bool isCliff;
    bool OutofBandIssue;
    bool bad_v2;
    bool isRFEvent;
    bool isPayloadBlast2;
    int box300;

		dataAllTree.SetBranchAddress("bad",&isBadEvent);
		dataAllTree.SetBranchAddress("weight",&weight);
		dataAllTree.SetBranchAddress("unixTime",&unixTime);
		dataAllTree.SetBranchAddress("isFirstFiveEvent",&isFirstFiveEvent);
		dataAllTree.SetBranchAddress("hasBadSpareChanIssue",&hasBadSpareChanIssue);
		dataAllTree.SetBranchAddress("hasBadSpareChanIssue2",&hasBadSpareChanIssue2);
		dataAllTree.SetBranchAddress("runNum",&runNum);
		dataAllTree.SetBranchAddress("eventNumber",&eventNumber);
		dataAllTree.SetBranchAddress("isSpikey",&isSpikey);
    dataAllTree.SetBranchAddress("isCliff",&isCliff);
    dataAllTree.SetBranchAddress("OutofBandIssue",&OutofBandIssue);
    dataAllTree.SetBranchAddress("bad_v2",&bad_v2);
    dataAllTree.SetBranchAddress("isRFEvent",&isRFEvent);
    dataAllTree.SetBranchAddress("isPayloadBlast2",&isPayloadBlast2);
    dataAllTree.SetBranchAddress("box300",&box300);

		double frac_of_power_notched_V[8];
		double frac_of_power_notched_H[8];

		stringstream ss;
		for(int i=0; i<8; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			dataVTree.SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_V[i]);
		}
		for(int i=8; i<16; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			dataHTree.SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_H[i-8]);
		}

		double VPeakOverRMS[16];
		int waveformLength[16];
		dataFilterTree.SetBranchAddress("VPeakOverRMS",&VPeakOverRMS);
		dataFilterTree.SetBranchAddress("waveformLength",&waveformLength);

		int numEntries = dataVTree.GetEntries();
		numTotal+=numEntries;

		dataAllTree.GetEvent(0);
		int currentRunNum = runNum;
		bool isThisABadRun = isBadRun(station,runNum,BadRunList);
		bool isThisASoftDomRun = isSoftwareDominatedRun("/users/PCON0003/cond0068/ARA/AraRoot/analysis/a23_analysis_tools", station, runNum);

		if(!isThisABadRun){
			isThisABadRun = isBadRun(station, runNum, BadSurfaceRunList);
		}

		for(int event=0; event<numDataEvents; event++){
			dataVTree.GetEvent(event);
			dataHTree.GetEvent(event);
			dataAllTree.GetEvent(event);
			dataFilterTree.GetEvent(event);

			if(runNum!=currentRunNum){
				currentRunNum=runNum;
				isThisABadRun = isBadRun(station,runNum, BadRunList);
				isThisASoftDomRun = isSoftwareDominatedRun("/users/PCON0003/cond0068/ARA/AraRoot/analysis/a23_analysis_tools", station, runNum);

				if(!isThisABadRun){
					isThisABadRun = isBadRun(station,runNum, BadSurfaceRunList);
				}
				if(isThisABadRun){
					printf(RED"*"RESET);
					// printf("     Yup, run %d is bad \n",runNum);
				}
				else{
					printf(GREEN"*"RESET);
				}
			}

			// set or re-set this stuff
			// for(int pol=0; pol<2; pol++){
			// 	hist_this_pol[pol]=0;
			// 	corr_val_out[pol]=-10000.;
			// 	snr_val_out[pol]=-10000.;
			// }
			// runNum_out=runNum;
			// eventNum_out=eventNumber;

			// continue;
			if( isSoft || isBadEvent || hasBadSpareChanIssue || hasBadSpareChanIssue2 || isFirstFiveEvent || isShort || isCal || isThisABadRun || isThisASoftDomRun|| isSpikey || isCliff || OutofBandIssue || bad_v2 || isPayloadBlast2 || box300){
				continue;
			}
			if(isBadLivetime(station,unixTime)){
				continue;
			}
			for(int pol=0; pol<2; pol++){
				if(!WFRMS[pol] && !isNewBox && !isSurf[0] && !isSurf[1] && !isSurfEvent_top[pol]){
					bool failsCWPowerCut=false;
					if(Refilt[pol] && !WFRMS[pol]){
						vector<double> frac;
						for(int i=0; i<8; i++){
							if(pol==0) frac.push_back(frac_of_power_notched_V[i]);
							else if(pol==1) frac.push_back(frac_of_power_notched_H[i]);
						}
						sort(frac.begin(), frac.end(), std::greater<double>());
						if(frac[2]>0.06){
							failsCWPowerCut=true;
						}
					} //refiltered?
					if(!failsCWPowerCut){
						bool this_pass_R_cut = passesRCut(station, config, pol, snr_val[pol], corr_val[pol]);

						// it's time to unblind
						// check events *passing* the cut
						if(this_pass_R_cut){


							printf(RED"Run %d, Event Number %d, pol %d, snr %.2f, corr-val %.2f, surfV %d, surfH %d, surfTopV %d, surfTopH %d \n"RESET, runNum, eventNumber, pol, snr_val[pol], corr_val[pol], isSurf[0], isSurf[1], isSurfEvent_top[0], isSurfEvent_top[1]);
							printf(RED"Run %d, Event Number %d, snrV %.2f, corrV %.4f, snrH %.2f, corrH %.4f\n"RESET, runNum, eventNumber, pol, snr_val[0], corr_val[0], snr_val[1], corr_val[1]);
							printf(RED"Run %d, Event Number %d, wfrmsV %d, wfrmsH %d \n"RESET, runNum, eventNumber, WFRMS[0], WFRMS[1]);

							// PlotThisEvent(station, runNum, eventNumber, pol);

							// set these correctly
							// hist_this_pol[pol]=1;
							// corr_val_out[pol]=corr_val[pol];
							// snr_val_out[pol]=snr_val[pol];
							// runNum_out=runNum;
							// eventNum_out=eventNumber;
							// outTree->Fill(); // save them to the file

						}
					}// not failing CW power cut?
				}// passes rest of analysis (not WFRMS, box, surface)
			}// loop over polarizations
		}// loop over events
	}

	// outFile->Write(); // write this down
	// outFile->Close(); // close it up
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

	char cw_file_name[400];
	sprintf(cw_file_name,"/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/CWID/A%d/all_runs/CWID_station_%d_run_%d.root",station,station,runNum);
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
		if(badFreqListLocal_baseline.size()>0){
			isCutonCW_baseline[pol]=true;
			// printf(RED"Cut on CW baseline in pol %d ?\n"RESET, pol);
		}
		for(int freq=0; freq<badFreqListLocal_baseline.size(); freq++){
			// printf(RED"    Freq %d is %.2f in pol %d \n"RESET, freq, badFreqListLocal_baseline[freq], pol);
		}
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
		// vector<int> chan_exclusion_list;
		// vector<double> baseline_CW_freqs = CWCut_TB(waveforms, average, pol, 6., 5.5, station, 3, chan_exclusion_list, runNum, event, true);
		// for(int i=0; i<baseline_CW_freqs.size(); i++){
		// 	printf("Pol %d, Freq %d, is %.2f \n", pol, i, baseline_CW_freqs[i]);
		// }
	}

	double threshCW = 1.5;
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
			}
		}
	}

	Settings *settings = new Settings();
	string setupfile = "setup.txt";
	settings->ReadFile(setupfile);
	cout << "Read " << setupfile << " file!" << endl;
	settings->NOFZ=1;
	Detector *detector=0;
	RayTraceCorrelator *theCorrelators[2];
	theCorrelators[0] =  new RayTraceCorrelator(station, 41., settings, 1, 4); //41 m, cal puser
	theCorrelators[1] =  new RayTraceCorrelator(station, 300., settings, 1, 4);//300 m, far reco
	int solNum = 0;
	bool isSimulation=false;

	bool skipCW=false;
	for(int pol=0; pol<2; pol++){
		if(skipCW || pol!=problempol) continue;
		if(isCutonCW_fwd[pol] || isCutonCW_back[pol] || isCutonCW_baseline[pol]){
			// printf(RED"	Has CW issue in pol %d \n"RESET, pol);
			// printf("		CW in FWD %d, BWD %d, or baseline %d? \n", isCutonCW_fwd[pol], isCutonCW_back[pol], isCutonCW_baseline[pol]);
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
			sprintf(save_temp_title,"%s/unblind/neutrino_box/passing_events/A%d_Run%d_Ev%d_ProblemPol%d_WaveformsNotch.png",plotPath,station,runNum,event,problempol);
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

			sprintf(save_temp_title,"%s/unblind/neutrino_box/passing_events/A%d_Run%d_Ev%d_ProblemPol%d_SpectraNotch.png",plotPath,station,runNum,event,problempol);
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
			sprintf(save_temp_title,"%s/unblind/neutrino_box/passing_events/A%d_Run%d_Ev%d_ProblemPol%d_FilteredMaps.png",plotPath,station,runNum,event,problempol);
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



	bool doMapsAtAll=false;
	if(doMapsAtAll){
		// get the run summary information, if it exists yet
		// and remember, because it's the users job to pass the location of the filter files
		// this should work for simulated events just fine
		char filter_file_name[400];
		// gonna try all of them cuz we lazy af lol
		sprintf(filter_file_name,"%s/processed_station_%d_run_%d_filter.root","/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ProcessedFile/A2/2015",station,runNum);
		bool hasFilterFile = false;
		TFile *filterFile = TFile::Open(filter_file_name);
		if(!filterFile){
			sprintf(filter_file_name,"%s/processed_station_%d_run_%d_filter.root","/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ProcessedFile/A2/2013",station,runNum);
		}
		filterFile = TFile::Open(filter_file_name);
		if(!filterFile){
			sprintf(filter_file_name,"%s/processed_station_%d_run_%d_filter.root","/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ProcessedFile/A2/2014",station,runNum);
		}
		filterFile = TFile::Open(filter_file_name);
		if(!filterFile){
			sprintf(filter_file_name,"%s/processed_station_%d_run_%d_filter.root","/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ProcessedFile/A2/2016",station,runNum);
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

		if(station==2){
			//for station 2, we need to exclude channel 15 from the analysis
			chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
			chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 13), chan_list_H.end());

			chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 5), chan_list_V.end());

		}
		else if(station==3){
			// // for station 3 years 2014, 2015, 2016, we need to drop string 4 (channels 3, 7, 11, 15) altogether above some run
			// if(
			// 	(!isSimulation && runNum>getA3BadRunBoundary())
			// 	||
			// 	(isSimulation && yearConfig>2)

			// ){			// drop string four
			// 	chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 3), chan_list_V.end());
			// 	chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 7), chan_list_V.end());
			// 	chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 11), chan_list_H.end());
			// 	chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
			// }
		}

		// Settings *settings = new Settings();
		// string setupfile = "setup.txt";
		// settings->ReadFile(setupfile);
		// cout << "Read " << setupfile << " file!" << endl;
		// settings->NOFZ=1;
		// Detector *detector = 0;

		// RayTraceCorrelator *theCorrelator = new RayTraceCorrelator(station, 41, settings, 1, RTTestMode);

		// int solNum = 0;
		TH2D *map_V_raytrace = theCorrelators[1]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V, chan_SNRs, solNum);
		TH2D *map_H_raytrace = theCorrelators[1]->getInterferometricMap_RT_select_NewNormalization_SNRweighted(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H, chan_SNRs, solNum);

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
			sprintf(save_temp_title,"/users/PCON0003/cond0068/ARA/AraRoot/analysis/unblind/neutrino_box/A%d_Run%d_Ev%d_Maps_300m_Exclude5and13.png",station,runNum,event);
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
		thisY.push_back(-700);
		thisY.push_back(700);
		dummy.push_back(new TGraph(2,&thisX[0], &thisY[0]));
	}

	char save_temp_title[300];
	sprintf(save_temp_title,"%s/unblind/neutrino_box/passing_events/A%d_Run%d_Ev%d_ProblemPol%d_Waveforms.png",plotPath,station,runNum,event,problempol);
	TCanvas *cWave = new TCanvas("","",4*1100,4*850);
	cWave->Divide(4,4);
	for(int i=0; i<16; i++){
		cWave->cd(i+1);
		dummy[i]->Draw("AP");
		dummy[i]->SetLineColor(kWhite);
		dummy[i]->GetXaxis()->SetRangeUser(-200.,700.);
		dummy[i]->GetYaxis()->SetRangeUser(-700.,700.);

		waveforms[i]->Draw("sameL");
		waveforms[i]->SetLineWidth(3);
	}
	cWave->SaveAs(save_temp_title);
	delete cWave;

	sprintf(save_temp_title,"%s/unblind/neutrino_box/passing_events/A%d_Run%d_Ev%d_ProblemPol%d_Spectra.png",plotPath,station,runNum,event,problempol);
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

	TCanvas *cSpec2 = new TCanvas("","",2*1100,850);
	cSpec2->Divide(2,1);
	cSpec2->cd(1);
		grWaveformsPowerSpectrum[0]->Draw("AL");
		grWaveformsPowerSpectrum[4]->Draw("sameL");
			grWaveformsPowerSpectrum[4]->SetLineColor(kRed);
		grWaveformsPowerSpectrum[8]->Draw("sameL");
			grWaveformsPowerSpectrum[8]->SetLineColor(kBlue);
		gPad->SetLogy();
		grWaveformsPowerSpectrum[0]->GetXaxis()->SetRangeUser(0.,450.);
		grWaveformsPowerSpectrum[0]->SetTitle("String 1 TV, BV, TH antenna spectra");

	cSpec2->cd(2);
		grWaveformsPowerSpectrum[3]->Draw("AL");
		grWaveformsPowerSpectrum[7]->Draw("sameL");
			grWaveformsPowerSpectrum[7]->SetLineColor(kRed);
		grWaveformsPowerSpectrum[11]->Draw("sameL");
			grWaveformsPowerSpectrum[11]->SetLineColor(kBlue);
		gPad->SetLogy();
		grWaveformsPowerSpectrum[3]->GetXaxis()->SetRangeUser(0.,450.);
		grWaveformsPowerSpectrum[3]->SetTitle("String 4 TV, BV, TH antenna spectra");
	sprintf(save_temp_title,"%s/unblind/neutrino_box/passing_events/A%d_Run%d_Ev%d_ProblemPol%d_SpectraString4.png",plotPath,station,runNum,event,problempol);
	cSpec2->SaveAs(save_temp_title);
	delete cSpec2;

	vector<TGraph*> spareElecChanGraphs;
	spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(6));
	spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(14));
	spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(22));
	spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(30));
	for(int i=0; i<4; i++){
		printf("Spar Elec chan index %d has RMS %.3f \n", i, spareElecChanGraphs[i]->GetRMS(2));
	}
	vector<TGraph*> dummy2;
	for(int i=0; i<4; i++){
		vector<double> thisX;
		vector<double> thisY;
		thisY.push_back(-200);
		thisY.push_back(200);
		thisX.push_back(0);
		thisX.push_back(600);
		dummy2.push_back(new TGraph(thisX.size(), &thisX[0], &thisY[0]));
		dummy2[i]->GetXaxis()->SetTitle("Time (ns)");
		dummy2[i]->GetYaxis()->SetTitle("Voltage (mV)");
		if(i==0)
			dummy2[i]->SetTitle("Chan 6");
		if(i==1)
			dummy2[i]->SetTitle("Chan 14");
		if(i==2)
			dummy2[i]->SetTitle("Chan 22");
		if(i==3)
			dummy2[i]->SetTitle("Chan 30");
		}
		TCanvas *cSpare = new TCanvas("","",4*850,850);
		cSpare->Divide(4,1);
		for(int i=0; i<4; i++){
			cSpare->cd(i+1);
			dummy2[i]->Draw("AP");
			dummy2[i]->SetLineColor(kWhite);
			spareElecChanGraphs[i]->Draw("sameL");
		}
		sprintf(save_temp_title,"%s/unblind/neutrino_box/passing_events/A%d_Run%d_Ev%d_ProblemPol%d_SpareChans.png",plotPath,station,runNum,event,problempol);
		cSpare->SaveAs(save_temp_title);
		delete cSpare;

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
