////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  find_hit_time.cxx
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
#include "FFTtools.h"
#include "TLine.h"
#include "TCanvas.h"

RawAtriStationEvent *rawAtriEvPtr;

#include "AraGeomTool.h"
#include "AraQualCuts.h"
#include "AraAntennaInfo.h"

#include "tools.h"

int main(int argc, char **argv)
{

	if(argc<3) {
		std::cout << "Usage\n" << argv[0] << " <station> <data_file> \n";
		return -1;
	}
		
	TFile *fp = TFile::Open(argv[2]);
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
	eventTree->SetBranchAddress("event",&rawAtriEvPtr);

	stringstream ss;
	string xLabel, yLabel;
	vector<string> titlesForGraphs;
	for (int i = 0; i < 16; i++){
		ss.str("");
		ss << "Channel " << i;
		titlesForGraphs.push_back(ss.str());
	}

	// get an event
	eventTree->GetEntry(39072);
	cout<<"Event number is "<<rawAtriEvPtr->eventNumber<<endl;
	UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
	vector<TGraph*> grWaveformsRaw;
	vector<TGraph*> grWaveformsInt;

	for(int i=0; i<16; i++){
		grWaveformsRaw.push_back(realAtriEvPtr->getGraphFromRFChan(i));
		grWaveformsInt.push_back(FFTtools::getInterpolatedGraph(grWaveformsRaw[i], 0.5));
	};

	int numBinsToIntegrate=10.;
	xLabel = "Time (ns)"; yLabel = "Integrated Power (arb units)";
	vector<TGraph*> grIntPower = makeIntegratedBinPowerGraphs(grWaveformsInt, numBinsToIntegrate, xLabel, yLabel, titlesForGraphs);

	vector<vector<double> > vvHitTimes; // vector of vector of hit times (first index is antenna, second index is hit time)
	vector<vector<double> > vvPeakIntPowers; // vector of vector of peak values
	int numSearchPeaks=2; //only look for two peaks
	double peakSeparation=10.0;
	getAbsMaximum_N(grIntPower, numSearchPeaks, peakSeparation, vvHitTimes, vvPeakIntPowers);

	for(int chan=0; chan<vvHitTimes.size(); chan++){
		for(int hit=0; hit<vvHitTimes[chan].size(); hit++){
			printf("%d, %d, value %.2f \n", chan, hit, vvHitTimes[chan][hit]);
		}
	}

	TLine *line1 = new TLine(vvHitTimes[9][0],0,vvHitTimes[9][0],8E6);
	TLine *line2 = new TLine(vvHitTimes[9][1],0,vvHitTimes[9][1],8E6);


	TCanvas *c = new TCanvas("","",2*850,850);
	c->Divide(2,1);
	c->cd(1);
		grWaveformsRaw[9]->Draw("ALP");
		grWaveformsRaw[9]->GetXaxis()->SetTitle("Time (ns)");
		grWaveformsRaw[9]->GetYaxis()->SetTitle("Voltage (mV)");
		grWaveformsRaw[9]->SetTitle("Channel 9");
	c->cd(2);
		grIntPower[9]->Draw("ALP");
		line1->Draw("same");
		line2->Draw("same");
		line1->SetLineColor(kRed);
		line2->SetLineColor(kRed);
		line1->SetLineWidth(2);
		line2->SetLineWidth(2);
		grIntPower[9]->GetYaxis()->SetTitleOffset(2.1);
	c->SaveAs("hit_finding_demo.png");

	fp->Close();
	delete fp;
}