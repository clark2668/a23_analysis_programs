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

RawAtriStationEvent *rawAtriEvPtr;
UsefulAtriStationEvent *realAtriEvPtr;

int main(int argc, char **argv)
{

	if(argc<7) {
		std::cout << "Usage\n" << argv[0] << " <simulation_flag> <station> <year> <run summary directory> <output directory> <input file> <pedestal file> \n";
		return -1;
	}

	/*
	arguments
	0: exec
	1: simulation (yes/no)
	2: station num (2/3)
	3: year
	4: run summary directory
	5: output directory
	6: input file
	7: pedestal file
	*/

	int isSimulation=atoi(argv[1]);
	int station_num=atoi(argv[2]);
	int year=atoi(argv[3]);

	//open the file
	TFile *fp = TFile::Open(argv[6]);
	if(!fp) {
		std::cerr << "Can't open file\n";
		return -1;
	}
	TTree *eventTree= (TTree*) fp->Get("eventTree");
	if(!eventTree) {
		std::cerr << "Can't find eventTree\n";
		return -1;
	}
	int runNum = getrunNum(argv[6]);
	printf("Run Number %d \n", runNum);

	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();	
	if(argc==8){
		//only if they gave us a pedestal should we fire up the calibrator
		calibrator->setAtriPedFile(argv[7],station_num);
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

	//first, let's get the baselines loaded in
	string runSummaryFilename;
	if (!isSimulation){
		runSummaryFilename = getRunSummaryFilename(station_num, argv[4], argv[6]);
	}
	else {
		if(station_num==2){
			runSummaryFilename = "/fs/scratch/PAS0654/ara/sim/RunSummary/run_summary_station_2_run_20.root";
		}
		else if(station_num==3){
			runSummaryFilename = "/fs/scratch/PAS0654/ara/sim/RunSummary/run_summary_station_3_run_0.root";
		}
	}
	TFile *SummaryFile = TFile::Open(runSummaryFilename.c_str());
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
	char del_me_file_name[400];
	sprintf(del_me_file_name,"%s/delme_run%d.root",output_location.c_str(),runNum);
	TFile *tempFile = TFile::Open(del_me_file_name,"RECREATE");
	tempFile->cd();
	TTree *tempTree = new TTree("tempTree","tempTree");
	bool hasError=0;
	tempTree->Branch("hasError",&hasError);
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
		}
		else if(isSimulation){
			hasError=false;
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
				// temp_phs.push_back(new TGraph());
				temp_phs[i] = new TGraph();
			}
		}
		else{
			//otherwise, interpolate, pad, and get the phase
			vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(grWaveformsRaw, interpolationTimeStep, xLabel, yLabel, titlesForGraphs);
			vector<TGraph*> grWaveformsPadded = makePaddedGraphs(grWaveformsInt, 0, xLabel, yLabel, titlesForGraphs);
			for(int chan=0; chan<16; chan++){
				// temp_phs.push_back(getFFTPhase(grWaveformsPadded[chan],120.,1000.));
				temp_phs[chan] = getFFTPhase(grWaveformsPadded[chan],120.,1000.);
			}
			deleteGraphVector(grWaveformsInt);
			deleteGraphVector(grWaveformsPadded);
		}
		tempTree->Fill(); //fill the tree
		deleteGraphVector(grWaveformsRaw);
		for(int i=0; i<16; i++) delete temp_phs[i];
		// deleteGraphVector(temp_phs);
		if(!isSimulation) delete realAtriEvPtr;
	}
	tempFile->Write();

	temp_phs.clear();
	temp_phs.resize(16);

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

		/*
		vector <TGraph*> grInt;
		vector <TGraph*> grPad;
		vector <TGraph*> grSpectrum;

		for(int i=0; i<16; i++){
			grInt.push_back(FFTtools::getInterpolatedGraph(grWaveformsRaw[i],interpolationTimeStep));
			grPad.push_back(FFTtools::padWaveToLength(grInt[i],2048));
			grSpectrum.push_back(FFTtools::makePowerSpectrumMilliVoltsNanoSecondsdB(grPad[i]));
		}

		TCanvas *c2 = new TCanvas("","",1100,850);
		c2->Divide(4,4);
		for(int i=0; i<16; i++){
			c2->cd(i+1);
			grWaveformsRaw[i]->Draw("ALP");
		}
		char save_title2[300];
		sprintf(save_title2,"testgraphs%d.png",event);
		c2->SaveAs(save_title2);


		TCanvas *c = new TCanvas("","",1100,850);
		c->Divide(4,4);
		for(int i=0; i<16; i++){
			c->cd(i+1);
			grSpectrum[i]->Draw("ALP");
			average[i]->Draw("Lsame");
			average[i]->SetLineColor(kRed);
		}
		char save_title[300];
		sprintf(save_title,"test%d.png",event);
		c->SaveAs(save_title);
		*/
	
		if(!hasError){

			//before we do the phase variance, we should check for baseline violations	
			vector<double> baseline_CW_cut_V = CWCut_TB(grWaveformsRaw, average, 0, 6., 5.5, station_num, 3);
			vector<double> baseline_CW_cut_H = CWCut_TB(grWaveformsRaw, average, 1, 6., 5.5, station_num, 3);
			/*			
			for(int i=0; i<baseline_CW_cut_V.size(); i++){
				printf("V: Event %d Baseline CW Cut %.2f \n", event, baseline_CW_cut_V[i]);
			}
			for(int i=0; i<baseline_CW_cut_H.size(); i++){
				printf("H: Event %d Baseline CW Cut %.2f \n", event, baseline_CW_cut_H[i]);
			}
			*/

			//incorrectly for now (!), if is simulation, clear these out
			if(isSimulation){
				baseline_CW_cut_V.clear();
				baseline_CW_cut_H.clear();
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
				if(!hasError){
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
						// cout<<"Forward event "<<event<<" :: " <<pol<<" :: freq "<<badFreqs_temp[i]<<", sigma "<<badSigmas_temp[i]<<endl;
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
			}
			phases_backward.push_back(first_event_in_sequence_phases_backwards);
			//okay, now we need to try and move backwards
			int found_events_backwards=14;
			for(int event_next=event-1; event_next>=0;event_next--){
				//printf("			Trying to move backwards to event %d \n");
				if(found_events_backwards==0) break;
				tempTree->GetEntry(event_next);
				if(!hasError){
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
						// cout<<"Backward event "<<event<<" :: " <<pol<<" :: freq "<<badFreqs_temp[i]<<", sigma "<<badSigmas_temp[i]<<endl;
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


			//incorrectly for now (!), if is simulation, clear these out
			if(isSimulation){
				badFreqs_fwd.clear();
				badFreqs_back.clear();
				badSigmas_fwd.clear();
				badSigmas_back.clear();
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
	SummaryFile->Close();

	tempFile->Close();
	remove(del_me_file_name);

	printf("Done! Run Number %d \n", runNum);

}//end main