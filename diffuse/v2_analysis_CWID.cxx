////////////////////////////////////////////////////////////////////////////////
////	v2_CWID.cxx 
////	A23 diffuse, identify CW freequency
////
////	Nov 2018
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <iomanip>
#include <sstream>
#include <deque>

//AraRoot Includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "AraGeomTool.h"
#include "AraAntennaInfo.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"

#include "tools_WaveformFns.h"
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_CW.h"
#include "tools_Cuts.h"

RawAtriStationEvent *rawAtriEvPtr;
UsefulAtriStationEvent *realAtriEvPtr;

int main(int argc, char **argv)
{

	if(argc<6) {
		std::cout << "Usage\n" << argv[0] << " <simulation_flag> <station> <year> <output directory> <input file> <pedestal file> \n";
		return -1;
	}
	int isSimulation=atoi(argv[1]);
	int station_num=atoi(argv[2]);
	int year=atoi(argv[3]);

	//now, to be as general as possible, we need to check for what kind of event we're up against
	TFile *fp = TFile::Open(argv[5]);
	if(!fp) {
		std::cerr << "Can't open file\n";
		return -1;
	}
	TTree *eventTree= (TTree*) fp->Get("eventTree");
	if(!eventTree) {
		std::cerr << "Can't find eventTree\n";
		return -1;
	}
	int runNum = getrunNum(argv[5]);
	printf("Run Number %d \n", runNum);

	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();	
	if(argc==7){
		//only if they gave us a pedestal should we fire up the calibrator
		calibrator->setAtriPedFile(argv[6],station_num);
	}

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
  
	//first, let's get the baselines loaded in
	char summary_file_name[400];
	sprintf(summary_file_name,"/fs/scratch/PAS0654/ara/10pct/Baselines/A%d/%d/baseline_station_%d_run_%d.root",station_num,year,station_num, runNum);
	TFile *SummaryFile = TFile::Open(summary_file_name);
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

	//now set up the outputs
	string output_location = argv[4];
	char run_file_name[400];
	sprintf(run_file_name,"%s/CWID_station_%d_run_%d.root",output_location.c_str(),station_num, runNum);
	TFile *outFile = TFile::Open(run_file_name,"RECREATE");
	TTree *NewCWTree = new TTree("NewCWTree","NewCWTree");
	vector<vector<double> > badFreqs_fwd;
	vector<vector<double> > badFreqs_back;
	vector<vector<double> > badSigmas_fwd;
	vector<vector<double> > badSigmas_back;
	vector<vector<double> > badFreqs_baseline;
	NewCWTree->Branch("badFreqs_fwd",&badFreqs_fwd);
	NewCWTree->Branch("badSigmas_fwd",&badSigmas_fwd);
	NewCWTree->Branch("badFreqs_back",&badFreqs_back);
	NewCWTree->Branch("badSigmas_back",&badSigmas_back);
	NewCWTree->Branch("badFreqs_baseline",&badFreqs_baseline);

	/*
	In this re-factorization, we want to do all the event FFT's up front
	and write them to a local ROOT file which we can delete later
	This way we dont' have to constantly re-do the FFT's
	*/

	numEntries=300;

	char del_me_file_name[400];
	sprintf(del_me_file_name,"%s/delme_run%d.root",output_location.c_str(),runNum);
	TFile *tempFile = TFile::Open(del_me_file_name,"RECREATE");
	TTree *tempTree = new TTree("tempTree","tempTree");
	bool hasError=0;
	tempTree->Branch("hasError",&hasError);

	TGraph *temp_phs[16];
	tempTree->Branch("temp_phs",&temp_phs, "temp_phs[16]");
	
	// vector<TGraph*> temp_phs(16);
	// tempTree->Branch("temp_phs",&temp_phs);
	for(int event=0; event<numEntries; event++){
		eventTree->GetEntry(event);
		if (isSimulation == false){
			realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
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
		hasError = hasDigitizerIssue(grWaveformsRaw);
		if(hasError){
			//if it has a digitizer error, just push back junk
			for(int i=0; i<16; i++){
				temp_phs[i] = new TGraph();
				// temp_phs.push_back(new TGraph());
			}
		}
		else{
			//otherwise, interpolate, pad, and get the phase
			vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(grWaveformsRaw, 0.6, xLabel, yLabel, titlesForGraphs);
			vector<TGraph*> grWaveformsPadded = makePaddedGraphs(grWaveformsInt, 0, xLabel, yLabel, titlesForGraphs);
			for(int chan=0; chan<16; chan++){
				temp_phs[chan] = getFFTPhase(grWaveformsPadded[chan],120.,1000.);
				// temp_phs.push_back(getFFTPhase(grWaveformsPadded[chan],120.,1000.));
			}
			deleteGraphVector(grWaveformsInt);
			deleteGraphVector(grWaveformsPadded);
		}
		tempTree->Fill(); //fill the tree
		tempFile->Write();
		deleteGraphVector(grWaveformsRaw);
		// for(int i=0; i<16; i++) delete temp_phs[i];
		// deleteGraphVector(temp_phs);
		if(!isSimulation) delete realAtriEvPtr;
	}

	//now, to loop over events!
	for(Long64_t event=0;event<numEntries;event++){

		badFreqs_fwd.clear();
		badSigmas_fwd.clear();
		badFreqs_back.clear();
		badSigmas_back.clear();
		badFreqs_baseline.clear();
    
		if(event%starEvery==0) { std::cerr << "*"; }
    
		eventTree->GetEntry(event); //get the event

		if (isSimulation == false){
			realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
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

		tempTree->GetEntry(event);
	
		if(!hasError){

			//before we do the phase variance, we should check for baseline violations	
			vector<double> baseline_CW_cut_V = CWCut_TB(grWaveformsRaw, average, 0, 6., 5.5, station_num, 3);
			vector<double> baseline_CW_cut_H = CWCut_TB(grWaveformsRaw, average, 1, 6., 5.5, station_num, 3);
			for(int i=0; i<baseline_CW_cut_V.size(); i++){
				// printf("Event %d Baseline CW Cut %.2f \n", event, baseline_CW_cut_V[i]);
			}		
			badFreqs_baseline.push_back(baseline_CW_cut_V);
			badFreqs_baseline.push_back(baseline_CW_cut_H);

			const int numPols = 2; //how many polarization do we want to think about
			const int numEventsForPhaseVariance = 15; //how many events do we need for the phase variance technique?
			const int numPairs = 28; //7+6+5+4+3+2+1 pairs
  	
			vector<vector<deque<TGraph*> > > vvdGrPhaseDiff_fwd;
			vector<vector<deque<TGraph*> > > vvdGrPhaseDiff_back;
			vvdGrPhaseDiff_fwd.resize(numPols); //need an entry per polarization
			vvdGrPhaseDiff_back.resize(numPols); //need an entry per polarization
			for (int i = 0 ; i < numPols; i++){
				vvdGrPhaseDiff_fwd[i].resize(numPairs); //and for the number of pairs in that polarization
				vvdGrPhaseDiff_back[i].resize(numPairs); //and for the number of pairs in that polarization
			}

			// vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(grWaveformsRaw, 0.6, xLabel, yLabel, titlesForGraphs);
			// vector<TGraph*> grWaveformsPadded = makePaddedGraphs(grWaveformsInt, 0, xLabel, yLabel, titlesForGraphs);

			//forward case
			//assume the initial event is the "0th" entry, and try to go forward
			//if there aren't enough good events ahead of us (no errors) to do the forward case,
			//then we declare it a failure and leave the vector of bad freqs and such empty!

			vector<vector<TGraph*> > phases_forward;
			vector <TGraph*> first_event_in_sequence_phases_forward;
			for(int chan=0; chan<16; chan++){
				first_event_in_sequence_phases_forward.push_back((TGraph*) temp_phs[chan]->Clone());
				// first_event_in_sequence_phases_forward.push_back(getFFTPhase(grWaveformsPadded[chan],120.,1000.));
			}
			phases_forward.push_back(first_event_in_sequence_phases_forward);

			//okay, now we need to try and move forward
			int found_events_forward=0;
			for(int event_next=event+1; event_next<numEntries;event_next++){
				// printf("			Trying to move forwards to event %d \n",event_next);
				// printf("			I've found %d good events \n",found_events_forward);
				if(found_events_forward==14) break; //after you've collected 15 events (0->14), we're good to go.
				tempTree->GetEntry(event_next);
				// eventTree->GetEntry(event_next);
				// UsefulAtriStationEvent *realAtriEvPtr_next = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
				// vector<TGraph*> grWaveformsRaw_next = makeGraphsFromRF(realAtriEvPtr_next, 16, xLabel, yLabel, titlesForGraphs);
				if(!hasError){
					found_events_forward++;
					// vector<TGraph*> grWaveformsInt_next = makeInterpolatedGraphs(grWaveformsRaw_next, 0.6, xLabel, yLabel, titlesForGraphs);
					// vector<TGraph*> grWaveformsPadded_next = makePaddedGraphs(grWaveformsInt_next, 0, xLabel, yLabel, titlesForGraphs);
					vector<TGraph*> this_event_phases;
					for(int chan=0; chan<16; chan++){
						this_event_phases.push_back((TGraph*) temp_phs[chan]->Clone());
						// this_event_phases.push_back(getFFTPhase(grWaveformsPadded_next[chan],120.,1000.));
					}
					phases_forward.push_back(this_event_phases);
					// for(int chan=0; chan<16; chan++) {delete grWaveformsInt_next[chan]; delete grWaveformsPadded_next[chan];}  //cleanup
				}
				// for(int chan=0; chan<16; chan++) delete grWaveformsRaw_next[chan]; //cleanup
				// delete realAtriEvPtr_next;
			}
			//if we have enough events to conduct the CW check
			if(found_events_forward==14){
				// printf("	We have sufficient number of events to do phase variance calculation in forward direction\n");
				int chan1, chan2;
				for(int use_event=0; use_event<15; use_event++){ //loop over the events that we stored
					for(int pol=0; pol<numPols; pol++){ //loop over polarizations
						for(int pairIndex = 0; pairIndex < numPairs; pairIndex++){ //loop over pairs for that event and polarization
							getChansfromPair(geomTool,station_num,pol,pairIndex,chan1,chan2); //get chan numbers for this pair and pol
							if (chan1 != -1 && chan2 != -1){
								vvdGrPhaseDiff_fwd[pol][pairIndex].push_back(getPhaseDifference(phases_forward[use_event][chan1], phases_forward[use_event][chan2]));
							}
						}
					}
				}
				//printf("	Got phase difference; on to phase variance\n");
				vector<TGraph*> vGrSigmaVarianceAverage_fwd;
				vGrSigmaVarianceAverage_fwd.resize(numPols);
				for(int pol=0; pol<numPols; pol++){
					vGrSigmaVarianceAverage_fwd[pol] = getPhaseVariance(vvdGrPhaseDiff_fwd[pol]);
					vector<double> badFreqs_temp;
					vector<double> badSigmas_temp;
					double threshold = 1.0;
					getPeaksAboveThreshold(vGrSigmaVarianceAverage_fwd[pol], threshold, badFreqs_temp, badSigmas_temp);
					badFreqs_fwd.push_back(badFreqs_temp);
					badSigmas_fwd.push_back(badSigmas_temp);
					for(int i=0; i<badFreqs_temp.size(); i++){
						// cout<<"event "<<event<<" :: " <<pol<<" :: freq "<<badFreqs_temp[i]<<", sigma "<<badSigmas_temp[i]<<endl;
					}
					delete vGrSigmaVarianceAverage_fwd[pol];
				}
				for(int use_event=0; use_event<15; use_event++){
					for(int pol=0; pol<numPols; pol++){
						for(int pairIndex=0; pairIndex<numPairs; pairIndex++){
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
			
			//reverse case
			//assume the initial event is the "15th" entry, and try to go backwards
			//if there aren't enough good events ahead of us (no errors) to do the forward case,
			//then we declare it a failure and leave the vector of bad freqs and such empty!


			tempTree->GetEntry(event);

			vector<vector<TGraph*> > phases_backward;
			vector <TGraph*> first_event_in_sequence_phases_backwards;
			for(int chan=0; chan<16; chan++){
				first_event_in_sequence_phases_backwards.push_back((TGraph*) temp_phs[chan]->Clone());
				// first_event_in_sequence_phases_backwards.push_back(getFFTPhase(grWaveformsPadded[chan],120.,1000.));
			}
			phases_backward.push_back(first_event_in_sequence_phases_backwards);
			//okay, now we need to try and move backwards
			int found_events_backwards=14;
			for(int event_next=event-1; event_next>=0;event_next--){
				//printf("			Trying to move backwards to event %d \n");
				if(found_events_backwards==0) break;
				tempTree->GetEntry(event_next);
				// eventTree->GetEntry(event_next);
				// UsefulAtriStationEvent *realAtriEvPtr_prev = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
				// vector<TGraph*> grWaveformsRaw_prev = makeGraphsFromRF(realAtriEvPtr_prev, 16, xLabel, yLabel, titlesForGraphs);
				if(!hasError){
					found_events_backwards--;
					// vector<TGraph*> grWaveformsInt_prev = makeInterpolatedGraphs(grWaveformsRaw_prev, 0.6, xLabel, yLabel, titlesForGraphs);
					// vector<TGraph*> grWaveformsPadded_prev = makePaddedGraphs(grWaveformsInt_prev, 0, xLabel, yLabel, titlesForGraphs);
					vector<TGraph*> this_event_phases;
					for(int chan=0; chan<16; chan++){
						this_event_phases.push_back((TGraph*) temp_phs[chan]->Clone());
						// this_event_phases.push_back(getFFTPhase(grWaveformsPadded_prev[chan],120.,1000.));
					}
					phases_backward.push_back(this_event_phases);
					// for(int chan=0; chan<16; chan++) {delete grWaveformsInt_prev[chan]; delete grWaveformsPadded_prev[chan];}  //cleanup
				}
				// for(int chan=0; chan<16; chan++) delete grWaveformsRaw_prev[chan]; //cleanup
				// delete realAtriEvPtr_prev;
			}
			//if we have enough events to conduct the CW check
			if(found_events_backwards==0){
				//printf("	We have sufficient number of events to do phase variance calculation in backward direction\n");
				int chan1, chan2;
				for(int use_event=0; use_event<15; use_event++){ //loop over the events that we stored
					for(int pol=0; pol<numPols; pol++){ //loop over polarizations
						for(int pairIndex = 0; pairIndex < numPairs; pairIndex++){ //loop over pairs for that event and polarization
							getChansfromPair(geomTool,station_num,pol,pairIndex,chan1,chan2); //get chan numbers for this pair and pol
							if (chan1 != -1 && chan2 != -1){
								vvdGrPhaseDiff_back[pol][pairIndex].push_back(getPhaseDifference(phases_backward[use_event][chan1], phases_backward[use_event][chan2]));
							}
						}
					}
				}

				vector<TGraph*> vGrSigmaVarianceAverage_back;
				vGrSigmaVarianceAverage_back.resize(numPols);
				for(int pol=0; pol<numPols; pol++){
					vGrSigmaVarianceAverage_back[pol] = getPhaseVariance(vvdGrPhaseDiff_back[pol]);
					vector<double> badFreqs_temp;
					vector<double> badSigmas_temp;
					double threshold = 1.0;
					getPeaksAboveThreshold(vGrSigmaVarianceAverage_back[pol], threshold, badFreqs_temp, badSigmas_temp);
					badFreqs_back.push_back(badFreqs_temp);
					badSigmas_back.push_back(badSigmas_temp);
					for(int i=0; i<badFreqs_temp.size(); i++){
						// cout<<"event "<<event<<" :: " <<pol<<" :: freq "<<badFreqs_temp[i]<<", sigma "<<badSigmas_temp[i]<<endl;
					}					
					delete vGrSigmaVarianceAverage_back[pol];
				}
				for(int use_event=0; use_event<15; use_event++){
					for(int pol=0; pol<numPols; pol++){
						for(int pairIndex=0; pairIndex<numPairs; pairIndex++){
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
			// //cleanup these intermediate waveforms
			// for(int i=0; i<16; i++){
			// 	delete grWaveformsPadded[i];
			// 	delete grWaveformsInt[i];
			// }

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
	SummaryFile->Close();

	tempFile->Close();

	printf("Done! Run Number %d \n", runNum);

}//end main