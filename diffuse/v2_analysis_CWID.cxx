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

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"

RawAtriStationEvent *rawAtriEvPtr;
UsefulAtriStationEvent *realAtriEvPtr;

#include "AraGeomTool.h"
#include "AraAntennaInfo.h"
#include "tools_WaveformFns.h"
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_CW.h"

vector<double> CWCut_TB(vector <TGraph*> waveforms, vector <TGraph*> baselines, int pol, double dBCut, double dBCutBroad, int station, int num_coinc);

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
	sprintf(summary_file_name,"/data/user/brianclark/A23Diffuse/Baselines/A%d/%d/baseline_station_%d_run_%d.root",station_num,year,station_num, runNum);
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

		//before we do the phase variance, we should check for baseline violations	
		vector<double> baseline_CW_cut_V = CWCut_TB(grWaveformsRaw, average, 0, 6., 5.5, station_num, 3);
		vector<double> baseline_CW_cut_H = CWCut_TB(grWaveformsRaw, average, 1, 6., 5.5, station_num, 3);
		for(int i=0; i<baseline_CW_cut_V.size(); i++){
			// printf("Event %d Baseline CW Cut %.2f \n", event, baseline_CW_cut_V[i]);
		}		
		badFreqs_baseline.push_back(baseline_CW_cut_V);
		badFreqs_baseline.push_back(baseline_CW_cut_H);
		
		if(1==1){

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

			vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(grWaveformsRaw, 0.6, xLabel, yLabel, titlesForGraphs);
			vector<TGraph*> grWaveformsPadded = makePaddedGraphs(grWaveformsInt, 0, xLabel, yLabel, titlesForGraphs);

			//forward case
			//assume the initial event is the "0th" entry, and try to go forward
			//if there aren't enough good events ahead of us (no cal, no soft, no short) to do the forward case,
			//then we declare it a failure and leave the vector of bad freqs and such empty!

			vector<vector<TGraph*> > phases_forward;
			vector <TGraph*> first_event_in_sequence_phases_forward;
			for(int chan=0; chan<16; chan++){
				first_event_in_sequence_phases_forward.push_back(getFFTPhase(grWaveformsPadded[chan],120.,1000.));
			}
			phases_forward.push_back(first_event_in_sequence_phases_forward);

			//okay, now we need to try and move forward
			int found_events_forward=0;
			for(int event_next=event+1; event_next<numEntries;event_next++){
				//printf("			Trying to move forwards to event %d \n",event_next);
				//printf("			I've found %d good events \n",found_events_forward);
				if(found_events_forward==14) break; //after you've collected 15 events (0->14), we're good to go.
				eventTree->GetEntry(event_next);
				UsefulAtriStationEvent *realAtriEvPtr_next = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
				vector<TGraph*> grWaveformsRaw_next = makeGraphsFromRF(realAtriEvPtr_next, 16, xLabel, yLabel, titlesForGraphs);
				if(1==1){
					found_events_forward++;
					vector<TGraph*> grWaveformsInt_next = makeInterpolatedGraphs(grWaveformsRaw_next, 0.6, xLabel, yLabel, titlesForGraphs);
					vector<TGraph*> grWaveformsPadded_next = makePaddedGraphs(grWaveformsInt_next, 0, xLabel, yLabel, titlesForGraphs);
					vector<TGraph*> this_event_phases;
					for(int chan=0; chan<16; chan++){
						this_event_phases.push_back(getFFTPhase(grWaveformsPadded_next[chan],120.,1000.));
						//phases_forward[found_events_forward][chan] = getFFTPhase(grWaveformsPadded_next[chan], 120, 1000);
					}
					phases_forward.push_back(this_event_phases);
					for(int chan=0; chan<16; chan++) {delete grWaveformsInt_next[chan]; delete grWaveformsPadded_next[chan];}  //cleanup
				}
				for(int chan=0; chan<16; chan++) delete grWaveformsRaw_next[chan]; //cleanup
				delete realAtriEvPtr_next;
			}

			//if we have enough events to conduct the CW check
			if(found_events_forward==14){
				//printf("	We have sufficient number of events to do phase variance calculation in forward direction\n");
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
			//if there aren't enough good events behind us (no cal, no soft, no short) to do the backward case
			//then we declare it a failure and leave the vector of bad freqs and such empty!

			vector<vector<TGraph*> > phases_backward;
			vector <TGraph*> first_event_in_sequence_phases_backwards;
			for(int chan=0; chan<16; chan++){
				first_event_in_sequence_phases_backwards.push_back(getFFTPhase(grWaveformsPadded[chan],120.,1000.));
			}
			phases_backward.push_back(first_event_in_sequence_phases_backwards);

			//okay, now we need to try and move backwards
			int found_events_backwards=14;
			for(int event_next=event-1; event_next>=0;event_next--){
				//printf("			Trying to move backwards to event %d \n");
				if(found_events_backwards==0) break;
				eventTree->GetEntry(event_next);
				UsefulAtriStationEvent *realAtriEvPtr_prev = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
				vector<TGraph*> grWaveformsRaw_prev = makeGraphsFromRF(realAtriEvPtr_prev, 16, xLabel, yLabel, titlesForGraphs);
				if(1==1){
					found_events_backwards--;
					vector<TGraph*> grWaveformsInt_prev = makeInterpolatedGraphs(grWaveformsRaw_prev, 0.6, xLabel, yLabel, titlesForGraphs);
					vector<TGraph*> grWaveformsPadded_prev = makePaddedGraphs(grWaveformsInt_prev, 0, xLabel, yLabel, titlesForGraphs);
					vector<TGraph*> this_event_phases;
					for(int chan=0; chan<16; chan++){
						this_event_phases.push_back(getFFTPhase(grWaveformsPadded_prev[chan],120.,1000.));
					}
					phases_backward.push_back(this_event_phases);
					for(int chan=0; chan<16; chan++) {delete grWaveformsInt_prev[chan]; delete grWaveformsPadded_prev[chan];}  //cleanup
				}
				for(int chan=0; chan<16; chan++) delete grWaveformsRaw_prev[chan]; //cleanup
				delete realAtriEvPtr_prev;
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
			//cleanup these intermediate waveforms
			for(int i=0; i<16; i++){
				delete grWaveformsPadded[i];
				delete grWaveformsInt[i];
			}

		}
		NewCWTree->Fill();
		for(int i=0; i<16; i++){
			delete grWaveformsRaw[i];
		}
		if(isSimulation==false){
			delete realAtriEvPtr;
		}
	}//loop over events

	outFile->Write();
	outFile->Close();
	fp->Close();
	SummaryFile->Close();

	printf("Done! Run Number %d \n", runNum);

}//end main

vector<double> CWCut_TB(vector <TGraph*> waveforms, vector <TGraph*> baselines, int pol, double dBCut, double dBCutBroad, int station, int num_coinc){
	double lowFreqLimit=120.;
	double highFreqLimit=900.;
	double halfrange = (highFreqLimit - lowFreqLimit)/2.;
	double halfway = double(lowFreqLimit + halfrange);
	const int numAnts = 16;
	double deltaTInt = 0.6;

	vector < vector < double> > badFreqs;
	badFreqs.resize(numAnts);
	vector < vector < double> > badFreqsBroad;
	badFreqsBroad.resize(numAnts);

	TGraph *baseline_clone[numAnts];
	vector <TGraph*> newFFTs;
	vector <TGraph*> newBaselines;

	double deltaF_save;

	AraGeomTool *geomTool = AraGeomTool::Instance();

	for(int ant=0; ant<numAnts; ant++){

		double magFFT[2000];
		double magFFTBegin[2000];
		double frequencyArray[2000];
		for(int i=0; i<2000; i++){
			magFFT[i]=0;
			frequencyArray[i]=-1;
		}

		int WaveformLength = 2048; //big, unfortunately...
		//what comes next is a not-so-obvious (imo) way of padding the waveform
		TGraph *chan1Int = FFTtools::getInterpolatedGraph(waveforms[ant],deltaTInt);
		double *getX = chan1Int->GetX();
		double deltaT = getX[1]-getX[0];
		while(chan1Int->GetN() < WaveformLength){
			double lastX = 0.;
			double lastY = 0.;
			chan1Int->GetPoint(chan1Int->GetN()-1,lastX,lastY);
			chan1Int->SetPoint(chan1Int->GetN(), lastX + deltaT, 0 );
		}
		double *getY = chan1Int->GetY();
		int length=chan1Int->GetN();
		FFTWComplex *theFFT=FFTtools::doFFT(length,getY);

		int newLength=(length/2)+1;
		double deltaF=1/(deltaT*length); //Hz
		deltaF*=1e3; //MHz
		deltaF_save=deltaF;

		for(int i=1;i<newLength;i++) {
			if (i==0) frequencyArray[i]=0.;
			if (i>0) frequencyArray[i]=frequencyArray[i-1]+deltaF;
			if (frequencyArray[i]>=lowFreqLimit && frequencyArray[i]<=highFreqLimit){
				magFFT[i]+=theFFT[i].re*theFFT[i].re+theFFT[i].im*theFFT[i].im;
			}
			else{
				magFFT[i]=-1000;
			}
		}
		delete chan1Int;
		delete [] theFFT;

		for(int i=0;i<newLength;i++) {
			if (frequencyArray[i]>=lowFreqLimit && frequencyArray[i]<=highFreqLimit){
				magFFT[i]=10*log10(sqrt(magFFT[i]));
			}
			else magFFT[i]=-1000;
		}
		
		//need to copy the baselines into new graphs
		double xx, yy;
		double bXreal[newLength];
		double bYreal[newLength];

		baseline_clone[ant] = new TGraph();

		for (int i3 = 0; i3 < newLength; i3++){
			baselines[ant]->GetPoint(i3, xx, yy);
			baseline_clone[ant]->SetPoint(i3, xx, yy);
			bXreal[i3] = xx;
			bYreal[i3] = yy;
		}
		double *bY = baselines[ant]->GetY();
		double *bX = baselines[ant]->GetX();
		int n = baselines[ant]->GetN();

		double *bY_unmodified=baseline_clone[ant]->GetY();
		double *bX_unmodified=baseline_clone[ant]->GetX();
		int n_unmodified=baseline_clone[ant]->GetN();

		
		// Calculate mean baseline
		//get baseline average so we can bump FFT around
		double mean=0;		
		double meanBaseline=0;
		int navg=0;
		int nfirsthalfavg=0;
		int nsecondhalfavg=0;
		double firstHalfMeanBaseline=0;
		double secondHalfMeanBaseline=0;
		double firstHalfMean=0;
		double secondHalfMean=0;

		for (int i=0;i<n;i++){
			if (bX[i]>=lowFreqLimit && bX[i]<highFreqLimit){
				meanBaseline+=bY[i];
				navg++;
			}
			if (bX[i]>=lowFreqLimit && bX[i]<halfway){
				firstHalfMeanBaseline+=bY[i];
				nfirsthalfavg++;
			}
			if (bX[i]>=halfway && bX[i]<highFreqLimit){
				secondHalfMeanBaseline+=bY[i];
				nsecondhalfavg++;
			}
		}
		meanBaseline=meanBaseline/double(navg);
		firstHalfMeanBaseline=firstHalfMeanBaseline/double(nfirsthalfavg);
		secondHalfMeanBaseline=secondHalfMeanBaseline/double(nsecondhalfavg);

		navg=0;
		nfirsthalfavg=0;
		nsecondhalfavg=0;

		//get average of graph in question
		for (int i=0;i<newLength;i++){
			if (frequencyArray[i]>=lowFreqLimit && frequencyArray[i]<highFreqLimit){
				mean+=magFFT[i];
				navg++;
			}
			if (frequencyArray[i]>=lowFreqLimit && frequencyArray[i]<halfway){
				firstHalfMean+=magFFT[i];
				nfirsthalfavg++;
			}
			if (frequencyArray[i]>=halfway && frequencyArray[i]<highFreqLimit){
				secondHalfMean+=magFFT[i];
				nsecondhalfavg++;
			}
		}
		mean=mean/double(navg);
		firstHalfMean=firstHalfMean/double(nfirsthalfavg);
		secondHalfMean=secondHalfMean/double(nsecondhalfavg);	

		//now bump the average to the baseline average and apply a tilt correction to baseline
		double deltaMean=mean-meanBaseline;
		double deltaMeanFirst=firstHalfMean-firstHalfMeanBaseline-deltaMean;
		double deltaMeanSecond=secondHalfMean-secondHalfMeanBaseline-deltaMean;
		double slope=(deltaMeanFirst-deltaMeanSecond)/halfrange;

		for (int i=0;i<newLength;i++){
			magFFTBegin[i]=magFFT[i];
			magFFT[i]=magFFT[i]-deltaMean;
		}
		for (int ctr=0;ctr<n;ctr++){
			if (bX[ctr]>=lowFreqLimit && bX[ctr]<highFreqLimit){
				bYreal[ctr]=bYreal[ctr]-slope*(bXreal[ctr]-halfway);
			}
		}
		
		//now see if any peaks are ndB above the baseline.
		double deltaMag[newLength];
		
		int j;
		for (int i=0;i<newLength;i++){
			if (frequencyArray[i]>lowFreqLimit+2 && frequencyArray[i]<highFreqLimit-2){
				for (j=0;j<n;j++){
					if (bX[j]>frequencyArray[i]) break; // finds bin where frequencies match
				}
				deltaMag[i]=magFFT[i]-bYreal[j]; // changed
			}
			else deltaMag[i]=-1000;
		}

		TGraph *magFFTgraph = new TGraph(2000,frequencyArray,magFFT);
		TGraph *newBaseline = new TGraph(2000,frequencyArray,bYreal);
		newFFTs.push_back(magFFTgraph);
		newBaselines.push_back(newBaseline);

		int index;
		for (int bin = 0; bin < newLength; bin++){
			if (deltaMag[bin] > dBCut){
				badFreqs[ant].push_back(frequencyArray[bin]);
			}
			if (deltaMag[bin] > dBCutBroad){
				badFreqsBroad[ant].push_back(frequencyArray[bin]);
			}
		}
	} //loop over antennas

	vector< vector <double> > freqMatching;
	double freqRangeBroad=40.;
	double freqRangeBroad_bins = freqRangeBroad/deltaF_save; //number of bins spanned by 40 MHz
	vector<double> FreqToNotch;

	for(int i=0; i<numAnts-1; i++){
		int i_pol = geomTool->getStationInfo(station)->getAntennaInfo(i)->polType;
		if(i_pol!=pol) continue; //check polarization
		for(int ii=0; ii<int(badFreqs[i].size()); ii++){ //loop over the bad frequencies for this antenna
			int matchedFreqs=1;
			int matchedFreqsFull=0;
			bool broad_band1=false;
			bool broad_band2=false;
			int broad_freqs1=0;
			int broad_freqs2=0;

			//loop over the other frequencies which were tagged as being above the "broad" threshold in this channel
			for(int ii2=0; ii2<int(badFreqsBroad[i].size()); ii2++){
				//if that frequency is w/in freqRangeBroad window
				if((abs(badFreqs[i][ii] - badFreqsBroad[i][ii2]) < freqRangeBroad)){
					broad_freqs1++;
				}
			}

			//if there are frequencies over the broad threshold w/in 40 MHz
			//then we want to know if the number of contaminated bins/number of bins in 40 MHz < 50%
			//if it's larger, then it's broadband, and we shouldn't touch it!
			if( (freqRangeBroad_bins) != 0){
				if(double(broad_freqs1-1)/double(int(freqRangeBroad_bins)-1) > 0.5){
					broad_band1=true;
				}
				else{
					matchedFreqsFull++;
				}
			}
			else{
				matchedFreqsFull++;
			}

			//now loop over all the other antennas in the array
			for(int j=i+1; j<numAnts; j++){
				int j_pol = geomTool->getStationInfo(station)->getAntennaInfo(j)->polType;
				if(j_pol!=i_pol) continue; //check polarization agreement
				bool matched_ant = false;
				//check all of their bad frequencies
				for(int jj=0; (jj<int(badFreqs[j].size()) && matched_ant==false); jj++){
					//if their bad frequencies happens within 5 MHz of the first bad frequency
					//we know something is up
					if((abs(badFreqs[i][ii] - badFreqs[j][jj]) < 5.)){

						broad_freqs2=0;
						matchedFreqs++;
						matched_ant=true;
						//now check all of the bad frequencies in this secondary antenna
                        for (int jj2 = 0; jj2 < int(badFreqsBroad[j].size()); jj2++){
                        	//if it's trouble make frequencies are w/in 40 MHz, then we want to record that
                        	if ((abs(badFreqs[j][jj] - badFreqsBroad[j][jj2]) < freqRangeBroad)){
                        		broad_freqs2++;
                        	}
                        }
                        //if there are frequencies over the broad threshold w/in 40 MHz
						//then we want to know if the number of contaminated bins/number of bins in 40 MHz < 50%
						//if it's larger, then it's broadband, and we shouldn't touch it!
                        if(freqRangeBroad_bins!=0){
                        	if(double(broad_freqs2-1)/double(int(freqRangeBroad_bins)-1) > 0.5){
                        		broad_band2=true;
                        	}
                        }
                        //if the first ant was *not* broadband, and *neither* was this one
                        //then we should say "yeah, we found something narrow, please notch me"
                        if((broad_band1==false) && (broad_band2==false)){
                        	matchedFreqsFull++;
                        } //were the trouble frequencies for both ants independently not broadband
					} //was the second ant's bad frequency within 5 MHz of the first antennas bad frequency?
				} //loop over second antennas bad freqs
			}//loop over second antenna
			if(matchedFreqsFull>=num_coinc){
				double new_freq=badFreqs[i][ii];
				for(int k=0; k<FreqToNotch.size(); k++){
					if(abs(new_freq)-FreqToNotch[k] < 0.01)
						new_freq=false;
				}
				if(new_freq)
					FreqToNotch.push_back(badFreqs[i][ii]);
			}
		}//loop over trouble frequencies for antenna 1
	} //loop over antenna 1

	for(int i=0; i<16; i++){
		delete baseline_clone[i];
		delete newFFTs[i];
		delete newBaselines[i];
	}
	return FreqToNotch;
}