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
	
	if(argc<7) {
		std::cout << "Usage\n" << argv[0] << " <simulation_flag> <station> <year> <radius_bin> <output directory> <input file> <pedestal file> \n";
		return -1;
	}

	/*
	arguments
	0: exec
	1: simulation (yes/no)
	2: station num (2/3)
	3: year (2013-2016)
	4: radius bin
	5: output directory
	6: input file
	7: pedestal file
	*/

	isSimulation=atoi(argv[1]);
	int station_num=atoi(argv[2]);
	int year = atoi(argv[3]);
	calpulserRunMode=0;
	int radiusBin = atoi(argv[4]);
	
	int numRadiiScanned = 35;
	int startingRadiusBin = radiusBin;

	AraQualCuts *qualCut = AraQualCuts::Instance(); //we also need a qual cuts tool
	
	stringstream ss;
	string xLabel, yLabel;
	vector<string> titlesForGraphs;
	for (int i = 0; i < nGraphs; i++){
	ss.str("");
	ss << "Channel " << i;
		titlesForGraphs.push_back(ss.str());
	}
		
	TFile *fp = TFile::Open(argv[6]);
	if(!fp) {
		std::cout << "Can't open file\n";
		return -1;
	}
	TTree *eventTree; 
	eventTree= (TTree*) fp->Get("eventTree");
	if(!eventTree) {
		std::cout << "Can't find eventTree\n";
		return -1;
	}
	
	TTree *simTree;
	TTree *simSettingsTree;
	Event *eventPtr = 0;
	Detector *detector = 0;
	Report *reportPtr = 0;
	
	if (isSimulation == true){
		simSettingsTree=(TTree*) fp->Get("AraTree");
		if (!simSettingsTree) {
			std::cout << "Can't find AraTree\n";
			return -1;
		}
	
		simSettingsTree->SetBranchAddress("detector", &detector);
		simSettingsTree->GetEntry(0);
	
		for (int i = 0; i < detector->stations.size(); i++){
			int n_antennas = 0;
			for (int ii = 0; ii < detector->stations[i].strings.size(); ii++){
				for (int iii = 0; iii < detector->stations[i].strings[ii].antennas.size(); iii++){
					detectorCenter[0] += detector->stations[i].strings[ii].antennas[iii].GetX();
					detectorCenter[1] += detector->stations[i].strings[ii].antennas[iii].GetY();
					detectorCenter[2] += detector->stations[i].strings[ii].antennas[iii].GetZ();
					n_antennas++;
				}
			}
			cout << "Detector Center: ";
			for (int ii = 0; ii < 3; ii++){
				detectorCenter[ii] = detectorCenter[ii]/(double)n_antennas;
				cout << detectorCenter[ii] << " : ";
			}
			cout << endl;
		}
	
	
		simTree=(TTree*) fp->Get("AraTree2");
		if (!simTree) {
			std::cout << "Can't find AraTree2\n";
			return -1;
		}
		simTree->SetBranchAddress("event", &eventPtr);
		simTree->SetBranchAddress("report", &reportPtr);
		simTree->GetEvent(0);
	}
		
	Settings *settings = new Settings();
	string setupfile = "setup.txt";
	settings->ReadFile(setupfile);
	cout << "Read " << setupfile << " file!" << endl;
	settings->NOFZ=1;
	
	numRadiiScanned=1;
	RayTraceCorrelator *theCorrelators[numRadiiScanned];
	for (int i = 0; i < numRadiiScanned; i++){
		theCorrelators[i] = 0;
	}

	for (int i = 0; i < numRadiiScanned; i++){
		double radius_temp = radii[i+startingRadiusBin];
		cout << "Setup RTCorr : " << radius_temp << endl;
		theCorrelators[i] = new RayTraceCorrelator(station_num, radius_temp, settings, angularBinSize, RTTestMode);
	}	

	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	
	if (argc == 8){
		cout << "Trying to load named pedestal" << endl;
		calibrator->setAtriPedFile(argv[7], station_num);
		cout << "Loaded named pedestal" << endl;
	} else {
		cout << "Trying to load blank pedestal" << endl;
		calibrator->setAtriPedFile("", station_num);
		cout << "Loaded blank pedestal" << endl;
	}
	
	double weight;
	int unixTime;
	int unixTimeUs;
	int eventNumber;
	double maxPeakVfromSim;
	double PeakVfromSim[16][2];

	eventTree->ResetBranchAddresses();
	
	if(isSimulation){
		eventTree->SetBranchAddress("UsefulAtriStationEvent", &realAtriEvPtr);
		eventTree->SetBranchAddress("weight", &weight);
		printf("Simulation; load useful event tree straight away \n");
	}
	else{
		eventTree->SetBranchAddress("event",&rawAtriEvPtr);
		printf("Data; load raw event tree \n");
	}
	
	Long64_t numEntries=eventTree->GetEntries();
	Long64_t starEvery=numEntries/80;
	if(starEvery==0) starEvery++;
	
	int runNum = getrunNum(argv[6]);
	printf("Reco Run Number %d \n", runNum);

	vector<TH2D*> maps;
	maps.resize(2);

	string processedFilename = getProcessedFilename_recoRadius(station_num, argv[5], argv[6], radii[radiusBin]);
	TFile *OutputFile = TFile::Open(processedFilename.c_str(), "RECREATE");
	TTree* OutputTree=new TTree("OutputTree", "OutputTree");
	// reconstruction information
	OutputTree->Branch("runReconstruction", &runReconstruction, "runReconstruction/O");

	OutputTree->Branch("peakCorr_single", &peakCorr_single, "peakCorr_single[2]/D");
	OutputTree->Branch("peakTheta_single", &peakTheta_single, "peakTheta_single[2]/I");
	OutputTree->Branch("peakPhi_single", &peakPhi_single, "peakPhi_single[2]/I");
	OutputTree->Branch("minCorr_single", &minCorr_single, "minCorr_single[2]/D");
	OutputTree->Branch("meanCorr_single", &meanCorr_single, "meanCorr_single[2]/D");
	OutputTree->Branch("rmsCorr_single", &rmsCorr_single, "rmsCorr_single[2]/D");
	OutputTree->Branch("peakSigma_single", &peakSigma_single, "peakSigma_single[2]/D");

	OutputTree->Branch("VMap",&maps[0]);
	OutputTree->Branch("HMap",&maps[1]);

	OutputTree->Branch("hasDigitizerError",&hasDigitizerError);
	OutputTree->Branch("isCalpulser", &isCalpulser, "isCalpulser/O");
	OutputTree->Branch("isSoftTrigger", &isSoftTrigger, "isSoftTrigger/O");
	OutputTree->Branch("unixTime", &unixTime);
	OutputTree->Branch("unixTimeUs", &unixTimeUs);
	OutputTree->Branch("eventNumber", &eventNumber);
	OutputTree->Branch("maxPeakVfromSim", &maxPeakVfromSim);
	OutputTree->Branch("PeakVfromSim", &PeakVfromSim, "peakVfromSim[16][2]/D");
	int isShort;
	OutputTree->Branch("isShort",&isShort);

	// simulation parameters
	OutputTree->Branch("weight", &weight_out, "weight/D");
	OutputTree->Branch("flavor", &flavor, "flavor/I");
	OutputTree->Branch("nu_nubar", &nu_nubar, "nu_nubar/I");
	OutputTree->Branch("energy", &energy, "energy/D");
	OutputTree->Branch("posnu", &posnu, "posnu[3]/D");
	OutputTree->Branch("viewAngle", &viewAngle, "viewAngle[16][2]/D");
	OutputTree->Branch("viewAngleAvg", &viewAngleAvg, "viewAngleAvg[2]/D");

	int eventSim = 0;
	numEntries=150;
	cerr<<"Run "<<runNum<<" has a starEvery of "<<starEvery<<endl;
	for(Long64_t event=100;event<numEntries;event++) {
		cout<<"On event "<<event<<endl;
		if(event%starEvery==0) {
			std::cerr << "*";     
		}
	
		eventTree->GetEntry(event);
		if (isSimulation == false){
			unixTime=(int)rawAtriEvPtr->unixTime;
			unixTimeUs=(int)rawAtriEvPtr->unixTimeUs;
			eventNumber=(int)rawAtriEvPtr->eventNumber;
		}else {
			eventNumber = event;
		}
	
		if (isSimulation == true){
			bool foundNextSimEvent = false;
		
			while (foundNextSimEvent == false){
				simTree->GetEntry(eventSim);
				if (reportPtr->stations[0].Global_Pass != 0 ){
					flavor = eventPtr->nuflavorint;
					nu_nubar = eventPtr->nu_nubar;
					energy = eventPtr->pnu;
					posnu[0] = eventPtr->Nu_Interaction[0].posnu.GetX();
					posnu[1] = eventPtr->Nu_Interaction[0].posnu.GetY();
					posnu[2] = eventPtr->Nu_Interaction[0].posnu.GetZ();
					weight = eventPtr->Nu_Interaction[0].weight;       
					maxPeakVfromSim = reportPtr->stations[0].max_PeakV;
					for (int i = 0; i < 4; i++){
						for (int ii = 0; ii < 4; ii++){
							int chan = ii +4*i;
							for (int j = 0; j < reportPtr->stations[0].strings[ii].antennas[i].PeakV.size(); j++){
								PeakVfromSim[chan][j] = reportPtr->stations[0].strings[ii].antennas[i].PeakV[j];
							}
						}
					}

					int avgCounter[2];
					avgCounter[0] = 0;       avgCounter[1] = 0;
					viewAngleAvg[0] = 0.;        viewAngleAvg[1] = 0.;
					for (int i = 0; i < 16; i++){
						for (int ii = 0; ii < 2; ii++){
							viewAngle[i][ii] = 0.; 
						}
					}
					for (int i = 0; i < reportPtr->stations[0].strings.size(); i++){
						for (int ii = 0; ii < reportPtr->stations[0].strings[i].antennas.size(); ii++){
							int channel = 4*i+ii;
							for (int iii = 0; iii < reportPtr->stations[0].strings[i].antennas[ii].view_ang.size(); iii++){
								viewAngleAvg[iii] += reportPtr->stations[0].strings[i].antennas[ii].view_ang[iii];
								avgCounter[iii]++;
								viewAngle[channel][iii] = reportPtr->stations[0].strings[i].antennas[ii].view_ang[iii];
							}
						}
					}
					for (int i = 0; i < 2; i++){
						if (avgCounter[i] == 0) {
							viewAngleAvg[i] = 0.;
						} else {
							viewAngleAvg[i] = viewAngleAvg[i]/(double)avgCounter[i];
						}
					}
					foundNextSimEvent=true;
				}
				eventSim++;
			}
		} else {
			posnu[0] = -10000000;
			posnu[1] = -10000000;
			posnu[2] = -10000000;
			flavor = -1;
			nu_nubar = -1;
			energy = -1.;
		}
	
		if(!isSimulation){
			realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
		}
		
		if (isSimulation){
			isCalpulser = false;
			isSoftTrigger = false;
		} else{
			isCalpulser = rawAtriEvPtr->isCalpulserEvent();
			isSoftTrigger = rawAtriEvPtr->isSoftwareTrigger();
			weight = 1.;
			isShort=0;
		}

		bool analyzeEvent = false;
		if (calpulserRunMode == 0) { analyzeEvent = true; } // analyze all events
		if (calpulserRunMode == 1 && isCalpulser == false && isSoftTrigger == false) { analyzeEvent = true; } // analyze only RF-triggered, non-calpulser events
		if (calpulserRunMode == 2 && isCalpulser == true) { analyzeEvent = true; } // analyze only calpulser events
		if (calpulserRunMode == 3 && isSoftTrigger == true) { analyzeEvent = true; } // analyze only software triggered  events
		if (analyzeEvent == true){

			weight_out = weight;
			if(!isSimulation)
				hasDigitizerError = !(qualCut->isGoodEvent(realAtriEvPtr));
			else
				hasDigitizerError=false;

			xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
			vector<TGraph*> grWaveformsRaw = makeGraphsFromRF(realAtriEvPtr, nGraphs, xLabel, yLabel, titlesForGraphs);
			ss.str("");
			for (int i = 0; i < 16; i++){
				if(grWaveformsRaw[i]->GetN()<500)
					isShort=1;
			}
			deleteGraphVector(grWaveformsRaw);

			//if the event has a  digitizer error, skip it
			if(hasDigitizerError){
				OutputTree->Fill(); //fill this anyway with garbage
				if (isSimulation == false) {
					delete realAtriEvPtr;
				}
				continue; //don't do any further processing on this event
			}
			// deleteGraphVector(grWaveformsRaw);

			int radiusBin_adjusted = radiusBin-startingRadiusBin;

			vector <int> chan_list_V;
			vector <int> chan_list_H;
			for(int chan=0; chan<=7; chan++){
				chan_list_V.push_back(chan);
				chan_list_H.push_back(chan+8);
			}

			if(station_num==2){
				//for station 2, we need to exclude channel 15 from the analysis
				chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
			}
			else if(station_num==3){
				//for station 3 years 2014 and 2015, we need to drop string 4 (channels 3, 7, 11, 15) altogether
				if(year==2014 || year==2015 || year==2016){

					chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 3), chan_list_V.end());
					chan_list_V.erase(remove(chan_list_V.begin(), chan_list_V.end(), 7), chan_list_V.end());

					chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 11), chan_list_H.end());
					chan_list_H.erase(remove(chan_list_H.begin(), chan_list_H.end(), 15), chan_list_H.end());
				}
			}
			// maps[0] = theCorrelators[radiusBin_adjusted]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V);
			// maps[1] = theCorrelators[radiusBin_adjusted]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H);
			// maps[0] = theCorrelators[radiusBin_adjusted]->getInterferometricMap_RT_select_SNRweighted(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V);
			// maps[1] = theCorrelators[radiusBin_adjusted]->getInterferometricMap_RT_select_SNRweighted(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H);
			maps[0] = theCorrelators[radiusBin_adjusted]->getInterferometricMap_RT_select_NewNormalization(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V);
			maps[1] = theCorrelators[radiusBin_adjusted]->getInterferometricMap_RT_select_NewNormalization(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H);
			// TH2D *map_V_raytrace = theCorrelators[radiusBin_adjusted]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Vpol, isSimulation, chan_list_V);
			// TH2D *map_H_raytrace = theCorrelators[radiusBin_adjusted]->getInterferometricMap_RT_select(settings, detector, realAtriEvPtr, Hpol, isSimulation, chan_list_H);

			getCorrMapPeak_wStats(maps[0], peakTheta_single[0], peakPhi_single[0], peakCorr_single[0], minCorr_single[0], meanCorr_single[0], rmsCorr_single[0], peakSigma_single[0]);
			getCorrMapPeak_wStats(maps[1], peakTheta_single[1], peakPhi_single[1], peakCorr_single[1], minCorr_single[1], meanCorr_single[1], rmsCorr_single[1], peakSigma_single[1]);
				
			OutputTree->Fill();

			delete maps[0];
			delete maps[1];

			if (isSimulation == false) {
				delete realAtriEvPtr;
			}
		}
	}
	
	OutputFile->Write();
	OutputFile->Close();
	fp->Close();
	
	for (int i = startingRadiusBin; i < startingRadiusBin-numRadiiScanned; i++){
		delete theCorrelators[i];
	}

	delete settings;

	cout<<endl;
	printf("Done! Run Number %d \n", runNum);
}
