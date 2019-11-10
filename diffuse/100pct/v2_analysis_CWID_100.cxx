////////////////////////////////////////////////////////////////////////////////
////	v2_CWID.cxx 
////	A23 diffuse, identify CW frequency
////	We only want to do CW on events that pas the filter
////	So we will need the filter file
////
////	Nov 2018
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <iomanip>
#include <sstream>
#include <deque>
#include <vector>

//AraRoot Includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraGeomTool.h"
#include "AraAntennaInfo.h"
#include "AraQualCuts.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"

#include "tools_inputParameters.h"
#include "tools_WaveformFns.h"
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_CW.h"
#include "tools_Cuts.h"
#include "tools_CommandLine.h"
#include "tools_outputObjects.h"

RawAtriStationEvent *rawAtriEvPtr;
UsefulAtriStationEvent *realAtriEvPtr;

int main(int argc, char **argv)
{

	if(argc<9) {
		std::cout << "Usage\n" << argv[0] << " <1-simulation_flag> <2-station> <3-year/config> <4-drop_bad_chans> <5-filter_file_dir> <6-run summary directory> <7-output directory> <8-input file> <9-pedestal file> \n";
		return -1;
	}

	int thresholdBin_pol[]={0,0};
	double wavefrontRMScut[]={2., 2.};

	/*
	arguments
	0: exec
	1: simulation (yes/no)
	2: station num (2/3)
	3: year/ config
	4: drop bad chans
	5: filter file directory
	6: run summary directory
	7: output directory
	8: input file
	9: pedestal file
	*/

	int isSimulation=atoi(argv[1]);
	int station_num=atoi(argv[2]);
	int yearConfig=atoi(argv[3]);
	int drop_bad_chans=atoi(argv[4]);
	drop_bad_chans=1; //always drop bad channels...

	//open the file
	TFile *fp = TFile::Open(argv[8]);
	if(!fp) {
		std::cerr << "Can't open file\n";
		return -1;
	}
	TTree *eventTree= (TTree*) fp->Get("eventTree");
	if(!eventTree) {
		std::cerr << "Can't find eventTree\n";
		return -1;
	}
	int runNum = getrunNum(argv[8]);
	printf("Filter Run Number %d \n", runNum);

	/*
		Adjust filter cut settings
		in A2, we were happy with snr in 0,0 and -1.3, -1.4 in V and H
		in A3, we want to do something a little more complicated and just thresholds
		by configuration
	*/

	int thisConfiguration = yearConfig;
	if(!isSimulation){
		thisConfiguration = getConfig(station_num, runNum);
	}
	getWFRMSCutValues(station_num, thisConfiguration, thresholdBin_pol[0], thresholdBin_pol[1], wavefrontRMScut[0], wavefrontRMScut[1]);
	// printf("Wavefront RMS cut values for config %d are vBin %d and hBin %d, vCut %.2f, hCut %.2f \n", thisConfiguration, thresholdBin_pol[0], thresholdBin_pol[1], wavefrontRMScut[0], wavefrontRMScut[1]);

	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();	
	if(argc==10){
		//only if they gave us a pedestal should we fire up the calibrator
		calibrator->setAtriPedFile(argv[9],station_num);
	}
	else{
		calibrator->setAtriPedFile("", station_num);
	}
	AraQualCuts *qualCut = AraQualCuts::Instance(); //we also need a qual cuts tool

	if(isSimulation){
		eventTree->SetBranchAddress("UsefulAtriStationEvent", &realAtriEvPtr);
		printf("Simulation; load useful event tree straight away \n");
	}
	else{
		eventTree->SetBranchAddress("event",&rawAtriEvPtr);
		printf("Data; load raw event tree \n");
	}
  
  	Long64_t numEntries=eventTree->GetEntries();
	Long64_t starEvery=numEntries/80;
	if(starEvery==0) starEvery++;
	printf("Num events is %d \n", numEntries);
	//numEntries=500;
	cout<<"This file has a starEvery of "<<starEvery<<endl;

	//first, let's get the baselines loaded in
	string runSummaryFilename;
	if (!isSimulation){
		runSummaryFilename = getRunSummaryFilename(station_num, argv[6], argv[8]);
	}
	else {
		if(station_num==2){
			runSummaryFilename = "/fs/scratch/PAS0654/ara/sim/RunSummary/run_summary_station_2_run_20.root";
		}
		else if(station_num==3){
			runSummaryFilename = "/fs/scratch/PAS0654/ara/sim/RunSummary/run_summary_station_3_run_30.root";
		}
	}
	TFile *SummaryFile = TFile::Open(runSummaryFilename.c_str(),"read");
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

	AraGeomTool * geomTool = new AraGeomTool();
	int nGraphs=16;

	vector<int> chan_exclusion_list;
	if(drop_bad_chans){
		if(station_num==2){
			//drop the last hpol channel
			chan_exclusion_list.push_back(15);
		}
		else if(station_num==3){
			if( 
				(!isSimulation && runNum>getA3BadRunBoundary())
				||
				(isSimulation && yearConfig>2)

			){			// drop string four
				printf("Yes, we need to drop channels!\n");
				chan_exclusion_list.push_back(3);
				chan_exclusion_list.push_back(7);
				chan_exclusion_list.push_back(11);
				chan_exclusion_list.push_back(15);
			}
		}
	}

	// now, we need the filter file
	// get the run summary information, if it exists yet
	// and remember, because it's the users job to pass the location of the filter files
	// this should work for simulated events just fine
	// note that also, we need to have information about the wavefront RMS filter
	// so that we can only do reconstructions on events which pass the WFRMS filter
	char filter_file_name[400];
	sprintf(filter_file_name,"%s/processed_station_%d_run_%d_filter.root",argv[5],station_num,runNum);
	TFile *filterFile = TFile::Open(filter_file_name,"READ"); //read only please
	if(!filterFile){
		std::cout << "Can't open filter file for run "<<runNum<<endl;
		return -1;
	}
	TTree *filterTree = (TTree*) filterFile->Get("OutputTree");
	if(!filterTree){
		std::cout << "Can't find filter tree\n";
		return -1;
	}
	// check for number of events mis-match (shouldn't happen, but you never now)
	int numEntriesFilter = filterTree->GetEntries();
	//numEntriesFilter=500;
	if(numEntries!=numEntriesFilter){
		printf("There is a mismatch between the number of entries in the data (%d) and the number in the fitler file (%d). Abort!\n", numEntries, numEntriesFilter);
		return -1;
	}
	filterTree->SetBranchAddress("VPeakOverRMS", &VPeakOverRMS);
	filterTree->SetBranchAddress("rms_pol_thresh_face_V", &rms_pol_thresh_face_V);
	filterTree->SetBranchAddress("rms_pol_thresh_face_H", &rms_pol_thresh_face_H);
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
	double rms_pol_thresh_face_alternate_V[thresholdSteps][numFaces_new_V];
	double rms_pol_thresh_face_alternate_H[thresholdSteps][numFaces_new_H];
	filterTree->SetBranchAddress("rms_pol_thresh_face_alternate_V", &rms_pol_thresh_face_alternate_V);
	filterTree->SetBranchAddress("rms_pol_thresh_face_alternate_H", &rms_pol_thresh_face_alternate_H);
	filterFile->cd();	

	//now set up the outputs
	string output_location = argv[7];
	char run_file_name[400];
	sprintf(run_file_name,"%s/CWID_station_%d_run_%d.root",output_location.c_str(),station_num, runNum);
	TFile *outFile = TFile::Open(run_file_name,"RECREATE");
	TTree *NewCWTree = new TTree("NewCWTree","NewCWTree");
	vector<vector<double> > badFreqs_fwd;
	vector<vector<double> > badFreqs_back;
	vector<vector<double> > badSigmas_fwd;
	vector<vector<double> > badSigmas_back;
	vector<vector<double> > badFreqs_baseline;
	bool passesWavefrontRMS[2];
	double bestFaceRMS[2];
	NewCWTree->Branch("badFreqs_fwd",&badFreqs_fwd);
	NewCWTree->Branch("badSigmas_fwd",&badSigmas_fwd);
	NewCWTree->Branch("badFreqs_back",&badFreqs_back);
	NewCWTree->Branch("badSigmas_back",&badSigmas_back);
	NewCWTree->Branch("badFreqs_baseline",&badFreqs_baseline);
	NewCWTree->Branch("passesWavefrontRMS", &passesWavefrontRMS, "passesWavefrontRMS[2]/O");
	NewCWTree->Branch("bestFaceRMS", &bestFaceRMS, "bestFaceRMS[2]/D");

	/*
	In this re-factorization, we want to do all the event FFT's up front
	and write them to a local ROOT file which we can delete later
	This way we dont' have to constantly re-do the FFT's
	*/

	printf(BLUE"About to make FFTs\n"RESET);

	char del_me_file_name[400];
	sprintf(del_me_file_name,"%s/delme_run%d.root",output_location.c_str(),runNum);
	TFile *tempFile = TFile::Open(del_me_file_name,"RECREATE");
	tempFile->cd();
	TTree *tempTree = new TTree("tempTree","tempTree");
	bool hasError=0;
	bool isCalPulser=0;
	tempTree->Branch("hasError",&hasError);
	tempTree->Branch("isCalPulser", &isCalPulser);
	vector<TGraph*> temp_phs;
	temp_phs.resize(16);
	stringstream ss;
	for(int i=0; i<16; i++){
		ss.str(""); ss<<"temp_phs_"<<i;
		tempTree->Branch(ss.str().c_str(),&temp_phs[i]);
	}
	for(int event=0; event<numEntries; event++){
		eventTree->GetEntry(event);
		if (isSimulation == false){
			realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
			hasError = !(qualCut->isGoodEvent(realAtriEvPtr));
			isCalPulser=rawAtriEvPtr->isCalpulserEvent();
			if(realAtriEvPtr->eventNumber<5){
				hasError=true;
			}
			if(station_num==3){ // for some reason, only activate this for A3 (oops)
				vector<TGraph*> spareElecChanGraphs;
				spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(6));
				spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(14));
				spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(22));
				spareElecChanGraphs.push_back(realAtriEvPtr->getGraphFromElecChan(30));
				bool hasBadSpareChanIssue = hasSpareChannelIssue(spareElecChanGraphs);
				deleteGraphVector(spareElecChanGraphs);
				if(hasBadSpareChanIssue){
					hasError=true;
				}
			}
		}
		else if(isSimulation){
			hasError=false;
			isCalPulser=false;
		}

		stringstream ss;
		string xLabel, yLabel;
		xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
		vector<string> titlesForGraphs;
		for (int i = 0; i < nGraphs; i++){
			ss.str("");
			ss << "Channel " << i;
			titlesForGraphs.push_back(ss.str());
		}
		vector<TGraph*> grWaveformsRaw = makeGraphsFromRF(realAtriEvPtr, 16, xLabel, yLabel, titlesForGraphs);
		if(hasError){
			//if it has a digitizer error, just push back junk
			for(int i=0; i<16; i++){
				temp_phs[i] = new TGraph();
			}
		}
		else{
			// otherwise, interpolate, pad, and get the phase
			vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(grWaveformsRaw, interpolationTimeStep, xLabel, yLabel, titlesForGraphs);
			vector<TGraph*> grWaveformsPadded = makePaddedGraphs(grWaveformsInt, 0, xLabel, yLabel, titlesForGraphs);
			for(int chan=0; chan<16; chan++){
				temp_phs[chan] = getFFTPhase(grWaveformsPadded[chan],120.,1000.);
			}
			deleteGraphVector(grWaveformsInt);
			deleteGraphVector(grWaveformsPadded);
		}
		tempTree->Fill(); //fill the tree
		deleteGraphVector(grWaveformsRaw);
		for(int i=0; i<16; i++) delete temp_phs[i];
		if(!isSimulation) delete realAtriEvPtr;
	}
	tempFile->Write();
	printf(GREEN"Done making FFTs.\n"RESET);

	temp_phs.clear();
	temp_phs.resize(16);

	printf(BLUE"About to loop over events\n"RESET);

	//now, to loop over events!
	for(Long64_t event=0;event<numEntries;event++){
		// cout<<"On event "<<event<<endl;

		badFreqs_fwd.clear();
		badSigmas_fwd.clear();
		badFreqs_back.clear();
		badSigmas_back.clear();
		badFreqs_baseline.clear();

		// always start by assuming they pass
		// NB: passing the WFRMS filter therefore means that the event might have other issues
		// so it shouldn't be taken as an exclusive way to reject the event
		passesWavefrontRMS[0] = true;
		passesWavefrontRMS[1] = true;
		// assigning these to 100 should make them stand out
		// because the value for actually fails the filter is 3
		bestFaceRMS[0]=100; 
		bestFaceRMS[1]=100;
    
		if(event%starEvery==0) { std::cerr << "*"; }

		// here we need to check if it passes the WFRMS filter
		filterTree->GetEntry(event);
		tempTree->GetEntry(event);

		vector <double> rms_faces_V;
		vector <double> rms_faces_H;

		// upgrade in Oct 2019 to deal string dropping in A3 correctly
		bool needToSwitchToAltWFRMSArrays=false;
		// we need an actual flag for if we need to switch to the _alternate_ arrays


		int num_faces_for_V_loop;
		int num_faces_for_H_loop;
		if(station_num==2){
			rms_faces_V.resize(numFaces);
			num_faces_for_V_loop=numFaces;
			rms_faces_H.resize(numFaces_A2_drop);
			num_faces_for_H_loop=numFaces_A2_drop;
			needToSwitchToAltWFRMSArrays=true;
		}
		else if(station_num==3){
			if(
				(!isSimulation && runNum>getA3BadRunBoundary())
				||
				(isSimulation && yearConfig>2)
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

		//now to sort them smallest to largest; lowest RMS is best
		sort(rms_faces_V.begin(), rms_faces_V.end());
		sort(rms_faces_H.begin(), rms_faces_H.end());

		bestFaceRMS[0]=rms_faces_V[0];
		bestFaceRMS[1]=rms_faces_H[0];
		
		// if it's log10(wavefront RMS is <-1.3 in V and <-1.4 in H, we're good to go!)
		if(log(bestFaceRMS[0])/log(10) >= wavefrontRMScut[0]){
			passesWavefrontRMS[0]=false;
		}
		if(log(bestFaceRMS[1])/log(10) >= wavefrontRMScut[1]){
			passesWavefrontRMS[1]=false;
		}

		// printf("Event %d: CalPul %d, Soft %d, WFRMS are %.2f, %.2f \n", event, isCalPulser, isSoftTrigger, TMath::Log10(bestFaceRMS[0]), TMath::Log10(bestFaceRMS[1]));

		// if it doesn't pass WFRMS filter in one of the two pols, no need to do reconstruction!
		if(isCalPulser || isSoftTrigger){
			// printf("No need to reco event %d this event! WFRMS are %.2f, %.2f \n",event, TMath::Log10(bestFaceRMS[0]), TMath::Log10(bestFaceRMS[1]));
			NewCWTree->Fill(); //fill this anyway with garbage
			continue; //don't do any further processing on this event
		}

		// if it doesn't pass WFRMS filter in one of the two pols, no need to do reconstruction!
		if(!passesWavefrontRMS[0] && !passesWavefrontRMS[1]){
			// printf("No need to reco event %d this event! WFRMS are %.2f, %.2f \n",event, TMath::Log10(bestFaceRMS[0]), TMath::Log10(bestFaceRMS[1]));
			NewCWTree->Fill(); //fill this anyway with garbage
			continue; //don't do any further processing on this event
		}

		// printf("Event %d good for CWID! CalPul %d, Soft %d, WFRMS are %.2f, %.2f \n", event, isCalPulser, isSoftTrigger, TMath::Log10(bestFaceRMS[0]), TMath::Log10(bestFaceRMS[1]));

		eventTree->GetEntry(event); //get the event

		if (isSimulation == false){
			realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
		}
		bool isSoftTrigger = rawAtriEvPtr->isSoftwareTrigger();

		stringstream ss;
		string xLabel, yLabel;
		xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
		vector<string> titlesForGraphs;
		for (int i = 0; i < nGraphs; i++){
			ss.str("");
			ss << "Channel " << i;
			titlesForGraphs.push_back(ss.str());
		}

		vector<TGraph*> grWaveformsRaw = makeGraphsFromRF(realAtriEvPtr, 16, xLabel, yLabel, titlesForGraphs);

		// if(!hasError){
		// for station 3, skip cal pulsers
		// so only proceed if it doen't have an error
		// and if it's not (A3 and a cal pulser)
		// this is done for backwards compatibility with A2, which did not exclude cal pulsers
		// from phase variance search
		if(!hasError && !(isCalPulser && station_num==3)){

			//before we do the phase variance, we should check for baseline violations	
			vector<double> baseline_CW_cut_V = CWCut_TB(grWaveformsRaw, average, 0, 6., 5.5, station_num, 3, chan_exclusion_list, runNum, event,false);
			vector<double> baseline_CW_cut_H = CWCut_TB(grWaveformsRaw, average, 1, 6., 5.5, station_num, 3, chan_exclusion_list, runNum, event,false);
			
			/*
			for(int i=0; i<baseline_CW_cut_V.size(); i++){
				printf(CYAN"	V: Event %d Baseline CW Cut %.2f \n"RESET, event, baseline_CW_cut_V[i]);
			}
			for(int i=0; i<baseline_CW_cut_H.size(); i++){
				printf(CYAN"	H: Event %d Baseline CW Cut %.2f \n"RESET, event, baseline_CW_cut_H[i]);
			}
			*/

			badFreqs_baseline.push_back(baseline_CW_cut_V);
			badFreqs_baseline.push_back(baseline_CW_cut_H);

			const int numPols = 2; //how many polarization do we want to think about
			const int numEventsForPhaseVariance = 15; //how many events do we need for the phase variance technique?
			const int numPairs = 28; // N(N-1)/2 unique pairs

			int numPairs_pol[2];

			if(drop_bad_chans){
				if(station_num==2){
					numPairs_pol[0]=numPairs;
					numPairs_pol[1]=21; // drop channel 15, only 21 pairs left
				}
				else if(station_num==3){
					if( 
						(!isSimulation && runNum>getA3BadRunBoundary())
						||
						(isSimulation && yearConfig>2)

					){  // drop string four (two channels per polarization)
						// only 15 pairs left
						numPairs_pol[0]=15;
						numPairs_pol[1]=15;
					}
					else{
						// drop no channels, keep them all
						numPairs_pol[0]=numPairs;
						numPairs_pol[1]=numPairs;
					}
				}
			}
			else{
				numPairs_pol[0]=numPairs;
				numPairs_pol[1]=numPairs;
			}

			/*
				So, when we want to drop channels, we need to change the number of pairs we will scan over with getPhaseVariance
				E.g., reduce to 21 pairs 
				But the only iterator we have is the global pair iterator, which runs from 0->28
				So, we make the size of the container 21 (or whatever)
				And just skip the global iterator pair when we find a pair we don't want
				That's the pair_in_use variable below, instead of using pairIndex directly
			*/


			vector<vector<deque<TGraph*> > > vvdGrPhaseDiff_fwd;
			vector<vector<deque<TGraph*> > > vvdGrPhaseDiff_back;
			vvdGrPhaseDiff_fwd.resize(numPols); //need an entry per polarization
			vvdGrPhaseDiff_back.resize(numPols); //need an entry per polarization
			for (int i = 0 ; i < numPols; i++){
				vvdGrPhaseDiff_fwd[i].resize(numPairs_pol[i]); //and for the number of pairs in that polarization
				vvdGrPhaseDiff_back[i].resize(numPairs_pol[i]); //and for the number of pairs in that polarization
			}

			//forward case
			//assume the initial event is the "0th" entry, and try to go forward
			//if there aren't enough good events ahead of us (no errors) to do the forward case,
			//then we declare it a failure and leave the vector of bad freqs and such empty!

			vector<vector<TGraph*> > phases_forward;
			vector <TGraph*> first_event_in_sequence_phases_forward;
			for(int chan=0; chan<16; chan++){
				first_event_in_sequence_phases_forward.push_back((TGraph*) temp_phs[chan]->Clone());
			}
			phases_forward.push_back(first_event_in_sequence_phases_forward);

			//okay, now we need to try and move forward
			int found_events_forward=0;
			for(int event_next=event+1; event_next<numEntries;event_next++){
				// printf("			Trying to move forwards to event %d \n",event_next);
				// printf("			I've found %d good events \n",found_events_forward);
				if(found_events_forward==14) break; //after you've collected 15 events (0->14), we're good to go.
				tempTree->GetEntry(event_next);
				if(!hasError && !(isCalPulser && station_num==3)){
					found_events_forward++;
					vector<TGraph*> this_event_phases;
					for(int chan=0; chan<16; chan++){
						this_event_phases.push_back((TGraph*) temp_phs[chan]->Clone());
					}
					phases_forward.push_back(this_event_phases);
				}
			}
			//if we have enough events to conduct the CW check
			if(found_events_forward==14){
				// printf("	We have sufficient number of events to do phase variance calculation in forward direction by event %d \n",event);
				int chan1, chan2;
				for(int use_event=0; use_event<15; use_event++){ //loop over the events that we stored
					for(int pol=0; pol<numPols; pol++){ //loop over polarizations
						int pair_in_use=0;
						for(int pairIndex = 0; pairIndex < numPairs; pairIndex++){ //loop over pairs for that event and polarization
							// cout<<"Working on global pair index "<<pairIndex<<endl;
							getChansfromPair(geomTool,station_num,pol,pairIndex,chan1,chan2); //get chan numbers for this pair and pol
							//now, make sure the fetch didn't fail, and that neither pair is in the "channel exclusion" list
							// printf("pairIndex %2d: Chan %2d, %2d \n", pairIndex, chan1,chan2);
							if(
								(std::find(chan_exclusion_list.begin(), chan_exclusion_list.end(), chan1) != chan_exclusion_list.end())
								||
								(std::find(chan_exclusion_list.begin(), chan_exclusion_list.end(), chan2) != chan_exclusion_list.end())
							){
								// if this global pair number contains a bad channel, we should skip it
								// cout<<"		Skipping global pair "<<pairIndex<<" which has channels "<<chan1<<" and "<<chan2<<endl;
								continue;
							}

							if (chan1 != -1
								&&
								chan2 != -1
								){
									vvdGrPhaseDiff_fwd[pol][pair_in_use].push_back(getPhaseDifference(phases_forward[use_event][chan1], phases_forward[use_event][chan2]));
									pair_in_use++;
							}
						}
					}
				}

				//printf("	Got phase difference; on to phase variance\n");
				vector<TGraph*> vGrSigmaVarianceAverage_fwd;
				vGrSigmaVarianceAverage_fwd.resize(numPols);
				for(int pol=0; pol<numPols; pol++){
					vGrSigmaVarianceAverage_fwd[pol] = getPhaseVariance(vvdGrPhaseDiff_fwd[pol], runNum, event, pol, false);
					vector<double> badFreqs_temp;
					vector<double> badSigmas_temp;
					double threshold = 1.0;
					getPeaksAboveThreshold(vGrSigmaVarianceAverage_fwd[pol], threshold, badFreqs_temp, badSigmas_temp);
					badFreqs_fwd.push_back(badFreqs_temp);
					badSigmas_fwd.push_back(badSigmas_temp);
					for(int i=0; i<badFreqs_temp.size(); i++){
						// cout<<"Forward event "<<event<<" :: pol " <<pol<<" :: freq "<<badFreqs_temp[i]<<", sigma "<<badSigmas_temp[i]<<endl;
					}
					delete vGrSigmaVarianceAverage_fwd[pol];
				}
				for(int use_event=0; use_event<15; use_event++){
					for(int pol=0; pol<numPols; pol++){
						for(int pairIndex=0; pairIndex<numPairs_pol[pol]; pairIndex++){
							delete vvdGrPhaseDiff_fwd[pol][pairIndex][use_event];
						}
					}
				}
			}
			for(int use_event=0; use_event<phases_forward.size();use_event++){
				for(int chan=0; chan<16; chan++){
					delete phases_forward[use_event][chan];
				}
			}
			// continue;

			//reverse case
			//assume the initial event is the "15th" entry, and try to go backwards
			//if there aren't enough good events ahead of us (no errors) to do the forward case,
			//then we declare it a failure and leave the vector of bad freqs and such empty!

			tempTree->GetEntry(event);

			vector<vector<TGraph*> > phases_backward;
			vector <TGraph*> first_event_in_sequence_phases_backwards;
			for(int chan=0; chan<16; chan++){
				first_event_in_sequence_phases_backwards.push_back((TGraph*) temp_phs[chan]->Clone());
			}
			phases_backward.push_back(first_event_in_sequence_phases_backwards);
			//okay, now we need to try and move backwards
			int found_events_backwards=14;
			for(int event_next=event-1; event_next>=0;event_next--){
				// printf("			Trying to move backwards to event %d \n",event_next);
				if(found_events_backwards==0) break;
				tempTree->GetEntry(event_next);
				if(!hasError && !(isCalPulser && station_num==3)){
					found_events_backwards--;
					vector<TGraph*> this_event_phases;
					for(int chan=0; chan<16; chan++){
						this_event_phases.push_back((TGraph*) temp_phs[chan]->Clone());
					}
					phases_backward.push_back(this_event_phases);
				}
			}
			//if we have enough events to conduct the CW check
			if(found_events_backwards==0){
				// printf("	We have sufficient number of events to do phase variance calculation in backward direction by event %d \n",event);
				int chan1, chan2;
				for(int use_event=0; use_event<15; use_event++){ //loop over the events that we stored
					for(int pol=0; pol<numPols; pol++){ //loop over polarizations
						int pair_in_use=0;
						for(int pairIndex = 0; pairIndex < numPairs; pairIndex++){ //loop over pairs for that event and polarization
							// cout<<"Working on global pair index "<<pairIndex<<endl;
							getChansfromPair(geomTool,station_num,pol,pairIndex,chan1,chan2); //get chan numbers for this pair and pol
							// printf("pairIndex %2d: Chan %2d, %2d \n", pairIndex, chan1,chan2);
							//now, make sure the fetch didn't fail, and that neither pair is in the "channel exclusion" list
							if(
								(std::find(chan_exclusion_list.begin(), chan_exclusion_list.end(), chan1) != chan_exclusion_list.end())
								||
								(std::find(chan_exclusion_list.begin(), chan_exclusion_list.end(), chan2) != chan_exclusion_list.end())
							){
								// if this global pair number contains a bad channel, we should skip it
								// cout<<"		Skipping global pair "<<pairIndex<<" which has channels "<<chan1<<" and "<<chan2<<endl;
								continue;
							}
							if (chan1 != -1
								&&
								chan2 != -1
								){
									vvdGrPhaseDiff_back[pol][pair_in_use].push_back(getPhaseDifference(phases_backward[use_event][chan1], phases_backward[use_event][chan2]));
									pair_in_use++;
							}
						}
					}
				}

				vector<TGraph*> vGrSigmaVarianceAverage_back;
				vGrSigmaVarianceAverage_back.resize(numPols);
				for(int pol=0; pol<numPols; pol++){
					vGrSigmaVarianceAverage_back[pol] = getPhaseVariance(vvdGrPhaseDiff_back[pol], runNum, event, pol, false);
					vector<double> badFreqs_temp;
					vector<double> badSigmas_temp;
					double threshold = 1.0;
					getPeaksAboveThreshold(vGrSigmaVarianceAverage_back[pol], threshold, badFreqs_temp, badSigmas_temp);
					badFreqs_back.push_back(badFreqs_temp);
					badSigmas_back.push_back(badSigmas_temp);
					for(int i=0; i<badFreqs_temp.size(); i++){
						// cout<<"Backward event "<<event<<" :: " <<pol<<" :: freq "<<badFreqs_temp[i]<<", sigma "<<badSigmas_temp[i]<<endl;
					}					
					delete vGrSigmaVarianceAverage_back[pol];
				}
				for(int use_event=0; use_event<15; use_event++){
					for(int pol=0; pol<numPols; pol++){
						for(int pairIndex=0; pairIndex<numPairs_pol[pol]; pairIndex++){
							delete vvdGrPhaseDiff_back[pol][pairIndex][use_event];
						}
					}
				}
			}
			for(int use_event=0; use_event<phases_backward.size();use_event++){
				for(int chan=0; chan<16; chan++){
					delete phases_backward[use_event][chan];
				}
			}
		}
		NewCWTree->Fill();
		deleteGraphVector(grWaveformsRaw);
		if(isSimulation==false){
			delete realAtriEvPtr;
		}
	}//loop over events

	outFile->Write();
	outFile->Close();
	fp->Close();
	filterFile->Close();
	SummaryFile->Close();

	tempFile->Close();
	remove(del_me_file_name); // remove CW file

	printf("Done! Run Number %d \n", runNum);

}//end main
