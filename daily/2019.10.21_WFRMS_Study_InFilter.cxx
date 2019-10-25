////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  v2_analysis_filter.cxx 
////  A23 diffuse, filter events
////
////  Nov 2018
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <iomanip>
#include <sstream>

//AraRoot Includes
#include "RawIcrrStationEvent.h"
#include "RawAtriStationEvent.h"
#include "UsefulAraStationEvent.h"
#include "UsefulIcrrStationEvent.h"
#include "UsefulAtriStationEvent.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"

RawAtriStationEvent *rawAtriEvPtr;
UsefulAtriStationEvent *realAtriEvPtr;

#include "Event.h"
#include "Detector.h"
#include "Report.h"

#include "AraGeomTool.h"
#include "AraQualCuts.h"
#include "AraAntennaInfo.h"

#include "tools_inputParameters.h"
#include "tools_outputObjects.h"
#include "tools_runSummaryObjects.h"
#include "tools_WaveformFns.h"
#include "tools_PlottingFns.h"
#include "tools_Constants.h"
#include "tools_RecoFns.h"
#include "tools_filterEvent.h"
#include "tools_Cuts.h"
#include "tools_CommandLine.h"

int main(int argc, char **argv)
{

	if(argc<7) {
		std::cout << "Usage\n" << argv[0] << " <simulation_flag> <station> <year/config> <run summary directory> <output directory> <input file> <pedestal file> \n";
		return -1;
	}

	/*
	arguments
	0: exec
	1: simulation (yes/no)
	2: station num (2/3)
	3: year/config
	4: run summary directory
	5: output directory
	6: input file
	7: pedestal file
	*/

	isSimulation=atoi(argv[1]);
	int station_num=atoi(argv[2]);
	calpulserRunMode=0; //analyze all events
	int yearConfig=atoi(argv[3]);

	if(isSimulation && yearConfig>5){
		std:cout<<"Warning! You have called for simulation with a config larger than 5, which cannot be true. Are you sure you didnt' accidentally use year instead of config?"<<endl;
		return -1;
	}

	TH1D *h1_SNR_dist[16][2];
	for(int chan=0; chan<16; chan++){
		for(int peak=0; peak<2; peak++){
			stringstream ss;
			ss<<"Chan "<<chan<<" peak "<<peak;
			h1_SNR_dist[chan][peak] = new TH1D(ss.str().c_str(), ss.str().c_str(), 100,0,5);
		}
	}

		
	stringstream ss;
	string xLabel, yLabel;
	vector<string> titlesForGraphs;
	for (int i = 0; i < nGraphs; i++){
		ss.str("");
		ss << "Channel " << i;
		titlesForGraphs.push_back(ss.str());
	}
	
	AraGeomTool * geomTool = new AraGeomTool();
	AraQualCuts *qualCut = AraQualCuts::Instance();
	
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
	
	vector<int> polarizations;
	vector<int> antenna_numbers;
	polarizations.resize(16);
	antenna_numbers.resize(16);
	vector< vector <double> > ant_loc;
	ant_loc.resize(16);
	for (int i = 0; i < 16; i++){
		ant_loc[i].resize(3);
		ant_loc[i][0] = geomTool->getStationInfo(station_num)->getAntennaInfo(i)->antLocation[0];
		ant_loc[i][1] = geomTool->getStationInfo(station_num)->getAntennaInfo(i)->antLocation[1];
		ant_loc[i][2] = geomTool->getStationInfo(station_num)->getAntennaInfo(i)->antLocation[2];
		polarizations[i] = (int)geomTool->getStationInfo(station_num)->getAntennaInfo(i)->polType;
		antenna_numbers[i]=i;
	}

	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	
	if (argc==8){
		cout << "Trying to load named pedestal" << endl;
		calibrator->setAtriPedFile(argv[7], station_num);
		cout << "Loaded named pedestal" << endl;
	} else {
		cout << "Trying to load blank pedestal" << endl;
		calibrator->setAtriPedFile("", station_num);
		cout << "Loaded blank pedestal" << endl;
	}
	
	vector<vector<vector<vector<int> > > > Pairs = setupPairs(station_num); // face type, polarization, pair #, 2 antennas
	double weight;
	int unixTime;
	int unixTimeUs;
	int eventNumber;
	double maxPeakVfromSim;
	double PeakVfromSim[16][2];
	int Trig_Pass[16] = {0};

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
	printf("Filter Run Number %d \n", runNum);
	printf("Num entries is %d \n", numEntries);
	numEntries=2;

	string runSummaryFilename;
	if (isSimulation == false){
		runSummaryFilename = getRunSummaryFilename(station_num, argv[4], argv[6]);
	}
	else {
	}
	TFile *SummaryFile = TFile::Open(runSummaryFilename.c_str());
	if(!SummaryFile) {
		std::cout << "Can't open summary file\n";
		return -1;
	}
	TTree* SummaryTree;
	SummaryTree = (TTree*) SummaryFile->Get("SummaryTree");   
	if(!SummaryTree) {
		std::cout << "Can't find SummaryTree\n";
		return -1;
	}
	
	SummaryTree->SetBranchAddress("RMS_RFTrigger", &RMS_RFTrigger);
	SummaryTree->SetBranchAddress("RMS_SoftTrigger", &RMS_SoftTrigger);
	SummaryTree->GetEntry(0);
	
	string processedFilename = getProcessedFilename_filter(station_num, argv[5], argv[6]);
	cout<<"processed filename is "<<processedFilename<<endl;
	TFile *OutputFile = TFile::Open(processedFilename.c_str(), "RECREATE");
	TTree* OutputSettingsTree = new TTree("OutputSettingsTree", "OutputSettingsTree");
	OutputSettingsTree->Branch("detectorCenter", &detectorCenter, "detectorCenter[3]/D");
	OutputSettingsTree->Branch("calpulserRunMode", &calpulserRunMode, "calpulserRunMode/I");
	OutputSettingsTree->Branch("numFaces", &numFaces_save, "numFaces/I");
	OutputSettingsTree->Branch("numFaces_A2_drop", &numFaces_A2_drop_save, "numFaces_A2_drop/I");
	OutputSettingsTree->Branch("numFaces_A3_drop", &numFaces_A3_drop_save, "numFaces_A3_drop/I");
	OutputSettingsTree->Branch("numSearchPeaks", &numSearchPeaks, "numSearchPeaks/I");
	OutputSettingsTree->Branch("thresholdMin", &thresholdMin, "thresholdMin/I");
	OutputSettingsTree->Branch("thresholdStep", &thresholdStep, "thresholdStep/D");
	OutputSettingsTree->Branch("thresholdSteps", &thresholdSteps_v, "thresholdSteps/D");
	OutputSettingsTree->Branch("interpolationTimeStep", &interpolationTimeStep, "interpolationTimeSteps/D");
	OutputSettingsTree->Branch("numBinsToIntegrate", &numBinsToIntegrate, "numBinsToIntegrate/I");
	OutputSettingsTree->Fill();
	
	TTree* OutputTree=new TTree("OutputTree", "OutputTree");
	
	// event summary information
	OutputTree->Branch("isCalpulser", &isCalpulser, "isCalpulser/O");
	OutputTree->Branch("isSoftTrigger", &isSoftTrigger, "isSoftTrigger/O");
	OutputTree->Branch("thirdVPeakOverRMS", &thirdVPeakOverRMS, "thirdVPeakOverRMS[3]/D");
	OutputTree->Branch("hasDigitizerError", &hasDigitizerError, "hasDigitizerError/O");
	OutputTree->Branch("unixTime", &unixTime);
	OutputTree->Branch("unixTimeUs", &unixTimeUs);
	OutputTree->Branch("eventNumber", &eventNumber);
	OutputTree->Branch("maxPeakVfromSim", &maxPeakVfromSim);
	OutputTree->Branch("PeakVfromSim", &PeakVfromSim, "peakVfromSim[16][2]/D");
	OutputTree->Branch("Trig_Pass", &Trig_Pass, "Trig_Pass[16]/I");


	// simulation parameters
	OutputTree->Branch("weight", &weight_out, "weight/D");
	OutputTree->Branch("flavor", &flavor, "flavor/I");
	OutputTree->Branch("nu_nubar", &nu_nubar, "nu_nubar/I");
	OutputTree->Branch("energy", &energy, "energy/D");
	OutputTree->Branch("posnu", &posnu, "posnu[3]/D");
	OutputTree->Branch("viewAngle", &viewAngle, "viewAngle[16][2]/D");
	OutputTree->Branch("viewAngleAvg", &viewAngleAvg, "viewAngleAvg[2]/D");
	
	// event channel summary information
	OutputTree->Branch("VPeak", &VPeak, "VPeak[16]/D");
	OutputTree->Branch("waveformRMS", &waveformRMS, "waveformRMS[16]/D");
	OutputTree->Branch("waveformRMS_50ns", &waveformRMS_50ns, "waveformRMS_50ns[16]/D");
	OutputTree->Branch("VPeakOverRMS", &VPeakOverRMS, "VPeakOverRMS[16]/D");
	OutputTree->Branch("avgPeakPower_5ns", &avgPeakPower_5ns, "avgPeakPower_5ns[16]/D");
	OutputTree->Branch("waveformLength", &waveformLength, "waveformLength[16]/I");
	
	
	// event filter information   
	OutputTree->Branch("TSQualParam", &TSQualParam, "TSQualParam/D");   
	OutputTree->Branch("rms_pol_thresh_face_V", &rms_pol_thresh_face_V, "rms_pol_thresh_face_V[15][12]/D");
	OutputTree->Branch("rms_pol_thresh_face_H", &rms_pol_thresh_face_H, "rms_pol_thresh_face_H[15][12]/D");


	int dropBadChans=1;
	int numFaces_new_V;
	int numFaces_new_H;
	if(station_num==2){
		numFaces_new_V = numFaces;
		numFaces_new_H = numFaces_A2_drop;
	}
	else if(station_num==3){
		numFaces_new_V = numFaces_A3_drop;
		numFaces_new_H = numFaces_A3_drop;
	}

	// list of what channels to exclude in finding the SNR of the event
	vector<int> chan_exclusion_list;
	if(station_num==2){
		// hpol channel 15
		chan_exclusion_list.push_back(15);
		printf("Station 2: Dropping ch 15\n");
	}
	else if(station_num==3){
		if( 
			(!isSimulation && runNum>getA3BadRunBoundary())
			||
			(isSimulation && yearConfig>2)

		){
			printf("Station 3: Dropping ch 3, 7, 11, 15\n");
			chan_exclusion_list.push_back(3);
			chan_exclusion_list.push_back(7);
			chan_exclusion_list.push_back(11);
			chan_exclusion_list.push_back(15);
		}
	}

	double rms_pol_thresh_face_alternate_V[thresholdSteps][numFaces_new_V];
	double rms_pol_thresh_face_alternate_H[thresholdSteps][numFaces_new_H];
	char rms_title_V[300];
	char rms_title_H[300];
	sprintf(rms_title_V,"rms_pol_thresh_face_alternate_V[15][%d]/D",numFaces_new_V); 
	sprintf(rms_title_H,"rms_pol_thresh_face_alternate_H[15][%d]/D",numFaces_new_H); 
	OutputTree->Branch("rms_pol_thresh_face_alternate_V", &rms_pol_thresh_face_alternate_V,rms_title_V);
	OutputTree->Branch("rms_pol_thresh_face_alternate_H", &rms_pol_thresh_face_alternate_H,rms_title_H);
	
	// polarization parameters
	OutputTree->Branch("polarizationRatio", &polarizationRatio, "polarizationRatio/D");   

	int eventSim = 0;
	int start=0;

	for(Long64_t event=start;event<numEntries;event++) {

		printf("On event %d \n", event);

		if(event%starEvery==0) {
			std::cout << "*";       
		}
		
		eventTree->GetEntry(event);
		if (isSimulation == false){
			unixTime=(int)rawAtriEvPtr->unixTime;
			unixTimeUs=(int)rawAtriEvPtr->unixTimeUs;
			eventNumber=(int)rawAtriEvPtr->eventNumber;
		} else{
			eventNumber = event;
		}

	
		if (isSimulation == true){
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
		}


		bool analyzeEvent = false;
		calpulserRunMode=1;
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
				waveformLength[i] = grWaveformsRaw[i]->GetN();
			}

			if(hasDigitizerError){
				OutputTree->Fill(); //fill this anyway with garbage
				deleteGraphVector(grWaveformsRaw);
				if (isSimulation == false) {
					delete realAtriEvPtr;
				}
				continue; //don't do any further processing on this event
			}
	
			double qualArray[4];
			filterEvent * filterEventPtr = new filterEvent();
			TSQualParam = -1.;
			TSQualParam = filterEventPtr->getQualityParameter(grWaveformsRaw, ant_loc, station_num, qualArray);
			delete filterEventPtr;
	
			xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
			vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(grWaveformsRaw, interpolationTimeStep, xLabel, yLabel, titlesForGraphs);


			bool printStuff=true;
			if(printStuff){
				TCanvas *c = new TCanvas("","",4*1100,4*850);
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
				sprintf(save_temp_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/wavefront_rms_teardown/Run%d_Ev%d_InterpolatedGraphs.png",runNum,event);
				TCanvas *cWave = new TCanvas("","",4*1100,4*850);
				cWave->Divide(4,4);
				for(int i=0; i<16; i++){
					cWave->cd(i+1);
					dummy[i]->Draw("AL");
					dummy[i]->SetLineColor(kWhite);
					// dummy[i]->GetXaxis()->SetRangeUser(300.,500.);

					grWaveformsInt[i]->Draw("sameL");
					grWaveformsInt[i]->SetLineWidth(3);
					// waveforms[i]->GetXaxis()->SetRangeUser(300.,500.);
				}
				cWave->SaveAs(save_temp_title);
				delete cWave;
			}

			vector<double> vVPeak;
			vector<double> vVPeakTimes;
	
			getAbsMaximum(grWaveformsInt, vVPeakTimes, vVPeak);
			copy(vVPeak.begin(), vVPeak.begin()+16, VPeak);
			copy(vVPeakTimes.begin(), vVPeakTimes.begin()+16, VPeakTimes);

			for(int i=0; i<vVPeak.size(); i++){
				// printf("	Chan %d: vVPeak %.2f and vVPeakTimes %.2f \n", i, vVPeak[i], vVPeakTimes[i]);
			}

			vector<double> vWaveformRMS;
			vWaveformRMS.resize(nGraphs);
			for (int i = 0; i < nGraphs; i++){
				if (RMS_SoftTrigger[i] == RMS_SoftTrigger[i] && 1==2){
					vWaveformRMS[i] = RMS_SoftTrigger[i];
				} else {
					if (RMS_RFTrigger[i] == RMS_RFTrigger[i]){
						vWaveformRMS[i] = RMS_RFTrigger[i];
					}
				}
				// printf("Chan %d: soft RMS %.3f and triggered RMS %.3f \n", i, RMS_SoftTrigger[i], RMS_RFTrigger[i]);
			}       
			copy(vWaveformRMS.begin(), vWaveformRMS.begin()+16, waveformRMS);
			for(int i=0; i<vVPeak.size(); i++){
				// printf("	Chan %d RMS %.2f \n", i, vWaveformRMS[i]);
			}

	
			vector<double> vWaveformRMS_50ns;
			getRMS(grWaveformsInt, vWaveformRMS_50ns, first50ns);
			copy(vWaveformRMS_50ns.begin(), vWaveformRMS_50ns.begin()+16, waveformRMS_50ns);
	
			vector<double> vVPeakOverRMS;
			vVPeakOverRMS.resize(16);
			for (int i = 0; i < 16; i++){
				VPeakOverRMS[i] = VPeak[i]/waveformRMS[i];
				vVPeakOverRMS[i] = VPeak[i]/waveformRMS[i];
			}
	
			vector<double> vThirdVPeakOverRMS;
			getThirdVPeakOverRMS(vVPeakOverRMS, polarizations, antenna_numbers, chan_exclusion_list, vThirdVPeakOverRMS);
			for (int i = 0 ; i< 3; i++){
				thirdVPeakOverRMS[i] = vThirdVPeakOverRMS[i];
			}
	
			xLabel = "Time (ns)"; yLabel = "Integrated Power (arb units)";
			vector<TGraph*> grIntPower = makeIntegratedBinPowerGraphs(grWaveformsInt, numBinsToIntegrate, xLabel, yLabel, titlesForGraphs);
			ss.str("");
	
			vector<double> hitTimes;
			vector<double> peakIntPowers;
			getAbsMaximum(grIntPower, hitTimes, peakIntPowers);

			for(int i=0; i<vVPeak.size(); i++){
				// printf("	Chan %d: peakIntPowers %.2f and hitTimes %.2f \n", i, peakIntPowers[i], hitTimes[i]);
			}
			printStuff=true;
			if(printStuff){
				TCanvas *c = new TCanvas("","",4*1100,4*850);
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
				sprintf(save_temp_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/wavefront_rms_teardown/Run%d_Ev%d_IntegratedPowerGraphs.png",runNum,event);
				TCanvas *cWave = new TCanvas("","",4*1100,4*850);
				cWave->Divide(4,4);
				for(int i=0; i<16; i++){
					cWave->cd(i+1);
					// dummy[i]->Draw("AL");
					// dummy[i]->SetLineColor(kWhite);
					// dummy[i]->GetXaxis()->SetRangeUser(300.,500.);

					grIntPower[i]->Draw("AL");
					grIntPower[i]->SetLineWidth(3);
					// waveforms[i]->GetXaxis()->SetRangeUser(300.,500.);
				}
				cWave->SaveAs(save_temp_title);
				delete cWave;
			}

			vector<vector<double> > vvHitTimes;
			vector<vector<double> > vvPeakIntPowers;
			getAbsMaximum_N(grIntPower, numSearchPeaks, 5.0, vvHitTimes, vvPeakIntPowers);
			for(int i=0; i<vvHitTimes.size(); i++){
				// printf("	Chan %d \n", i);
				for(int j=0; j<vvHitTimes[i].size(); j++){
					// printf("		Peak %d, power %.2f and time %.2f \n", j, vvPeakIntPowers[i][j], vvHitTimes[i][j]);
				}
			}

			vector<double> peakIntRMS;       
			for (int i = 0; i < peakIntPowers.size(); i++){
				peakIntRMS.push_back(sqrt(peakIntPowers[i]/numBinsToIntegrate));
			}
	
			for (int i = 0; i < 16; i++){
				avgPeakPower_5ns[i] = peakIntPowers[i]/numBinsToIntegrate;
				peakPowerTimes[i] = hitTimes[i];
			}
	
			vector<double> RMS_10overRMS;
			for (int i = 0; i < 16; i++){
				RMS_10overRMS.push_back(sqrt(avgPeakPower_5ns[i])/waveformRMS[i]);
			}
	
			vector<vector<double> > vvRMS_10overRMS;
			vvRMS_10overRMS.resize(16);
			for (int i = 0; i < 16; i++){
				vvRMS_10overRMS[i].resize(vvPeakIntPowers[i].size());
				for (int j = 0; j < vvPeakIntPowers[i].size(); j++){
					vvRMS_10overRMS[i][j] = sqrt(vvPeakIntPowers[i][j]/numBinsToIntegrate)/waveformRMS[i];
				}
			}
			for(int i=0; i<vvHitTimes.size(); i++){
				// printf("	Chan %d \n", i);
				for(int j=0; j<vvHitTimes[i].size(); j++){
					// printf("		Peak %d, SNR %.2f and time %.2f \n", j, vvRMS_10overRMS[i][j], vvHitTimes[i][j]);
					h1_SNR_dist[i][j]->Fill(vvRMS_10overRMS[i][j]);
				}
			}
		
			// polarizationRatio = getPolarizationRatio(grWaveformsRaw, polarizations);

			// //first apply the filter with *no* dropped channels
			// vector<vector<vector<vector<int> > > > faces = setupFaces(station_num);
			// for (int thresholdBin = 0; thresholdBin < thresholdSteps; thresholdBin++){
			// 	double threshold = thresholdMin + thresholdStep*(double)thresholdBin;
			// 	// if(thresholdBin==0) printf("     4-String Threshold is %.2f \n", threshold);
		
			// 	vector<double> rms_faces_V = getRms_Faces_Thresh_N(vvHitTimes, vvRMS_10overRMS, threshold, 0, faces, ant_loc);
			// 	vector<double> rms_faces_H = getRms_Faces_Thresh_N(vvHitTimes, vvRMS_10overRMS, threshold, 1, faces, ant_loc);
		
			// 	for (int i = 0; i < numFaces; i++){
			// 		// if(thresholdBin==0) printf("          4-string face %d value is  is %.4f \n", i, rms_faces_V[i]);
			// 		rms_pol_thresh_face_V[thresholdBin][i] = rms_faces_V[i];
			// 		rms_pol_thresh_face_H[thresholdBin][i] = rms_faces_H[i];
			// 	}
		
			// } // end threshold scan

			// //then apply the filter over again *with* dropped channels
			// vector<vector<vector<vector<int> > > > faces_drop = setupFaces(station_num, dropBadChans);
			// for (int thresholdBin = 0; thresholdBin < thresholdSteps; thresholdBin++){
			// 	double threshold = thresholdMin + thresholdStep*(double)thresholdBin;
			// 	// if(thresholdBin==0) printf("     3-String Threshold is %.2f \n", threshold);
		
			// 	vector<double> rms_faces_V_alternate = getRms_Faces_Thresh_N(vvHitTimes, vvRMS_10overRMS, threshold, 0, faces_drop, ant_loc);
			// 	vector<double> rms_faces_H_alternate = getRms_Faces_Thresh_N(vvHitTimes, vvRMS_10overRMS, threshold, 1, faces_drop, ant_loc);
		
			// 	for (int i = 0; i < numFaces_new_V; i++){
			// 		rms_pol_thresh_face_alternate_V[thresholdBin][i] = rms_faces_V_alternate[i];
			// 		// if(thresholdBin==0) printf("          3-string face %d value is %.4f \n", i, rms_faces_V_alternate[i]);
			// 	}
			// 	for (int i = 0; i < numFaces_new_H; i++){
			// 		rms_pol_thresh_face_alternate_H[thresholdBin][i] = rms_faces_H_alternate[i];
			// 	}
			// } // end threshold scan
				
			OutputTree->Fill();

			deleteGraphVector(grIntPower);
			deleteGraphVector(grWaveformsInt);
			deleteGraphVector(grWaveformsRaw);
			if (isSimulation == false) {
				delete realAtriEvPtr;
			}
		} //analyze event?
	} //loop over events

	gStyle->SetOptStat(1111111);

	TCanvas *c = new TCanvas("","",4*1100,4*850);
	c->Divide(4,4);
	for(int chan=0; chan<16; chan++){
		c->cd(chan+1);
		h1_SNR_dist[chan][0]->Draw("");
		h1_SNR_dist[chan][1]->Draw("same");
		h1_SNR_dist[chan][1]->SetLineColor(kRed);
		h1_SNR_dist[chan][1]->SetLineWidth(2);
		h1_SNR_dist[chan][0]->SetLineWidth(2);
		h1_SNR_dist[chan][0]->GetXaxis()->SetTitle("SNR");
		h1_SNR_dist[chan][0]->GetYaxis()->SetTitle("Number of Events");
		h1_SNR_dist[chan][0]->GetYaxis()->SetRangeUser(0,200);
	}
	char save_temp_title[300];
	sprintf(save_temp_title,"/users/PAS0654/osu0673/A23_analysis_new2/results/wavefront_rms_teardown/Run%d_SNRDist.png",runNum);
	c->SaveAs(save_temp_title);
	
	OutputFile->Write();
	OutputFile->Close();
	fp->Close();
	delete fp;
	cout<<endl;
	printf("Done! Run Number %d \n", runNum);

}