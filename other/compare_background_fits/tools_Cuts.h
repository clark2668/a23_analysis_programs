
//Includes
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <deque>
#include <complex>
#include <algorithm>
#include <numeric>
#include <fstream>

#include "TGraph.h"
#include "TMath.h"

#include "AraGeomTool.h"
#include "tools_PlottingFns.h"
#include "tools_WaveformFns.h"
#include "FFTtools.h"


using namespace std;

void getRCutValues(int station, int config, int pol, double &slope, double &intercept){
	if(station==2){
		// these were updated 2019-09-10 with the results of the optimization after the re-processing with corrected pedestals
		if(pol==0){
			slope=-2080.;
			if(config==1){
				intercept=20.6;
			}
			else if(config==2){
				intercept=22.6;
			}
			else if(config==3){
				intercept=20.8;
			}
			else if(config==4){
				intercept=21.5;
			}
			else if(config==5){
				intercept=21.6;
			}
		}
		else if(pol==1){
			slope=-640.;
			if(config==1){
				intercept=12.7;
			}
			else if(config==2){
				intercept=13.6;
			}
			else if(config==3){
				intercept=12.3;
			}
			else if(config==4){
				intercept=13.5;
			}
			else if(config==5){
				intercept=13.4;
			}
		}
	}
	else if(station==3){
		// these were updated 2019-11-25 with the results of the A3 optimization
		if(pol==0){
			// slope=-1240.;
			slope=-1420.;
			if(config==1){
				// intercept=14.6;
				intercept=15.8;
			}
			else if(config==2){
				// intercept=15.4;
				intercept=16.5;
			}
			else if(config==3){
				// intercept=17.9;
				intercept=19.5;
			}
			else if(config==4){
				// intercept=17.9;
				intercept=19.4;
			}
			else if(config==5){
				// intercept=16.1;
				intercept=17.5;
			}
		}
		else if(pol==1){
			// slope=-600.;
			slope=-680.;
			if(config==1){
				// intercept=13.4;
				intercept=14.1;
			}
			else if(config==2){
				// intercept=14.1;
				intercept=14.6;
			}
			else if(config==3){
				// intercept=14.0;
				intercept=14.8;
			}
			else if(config==4){
				// intercept=13.7;
				intercept=14.3;
			}
			else if(config==5){
				// intercept=12.3;
				intercept=13.0;
			}
		}
	}
}

bool passesRCut(int station, int config, int pol, double SNR, double corr){
	bool passes=false;
	if(station==2 || station==3){
		double slope_val;
		double intercept_val;
		getRCutValues(station, config, pol, slope_val, intercept_val);
		// passes if SNR > (corr * slope) + intercept
		if (SNR >= ((corr * slope_val ) + intercept_val)){
			passes=true;
		}
	}
	return passes;
}

int getA3BadRunBoundary(){
	return 1901;
}

int getConfig(int station, int runNum){
	int config=0;
	if(station==2){
		if(runNum>=0 && runNum<=4){
			config=1;
		}
		else if(runNum>=11 && runNum<=60){
			config=4;
		}
		else if(runNum>=120 && runNum<=2274){
			config=2;
		}
		else if(runNum>=2275 && runNum<=3463){
			config=1;
		}
		else if(runNum>=3465 && runNum<=4027){
			config=3;
		}
		else if(runNum>=4029 && runNum<=6481){
			config=4;
		}
		else if(runNum>=6500 && runNum<=8097){
			config=5;
		}
		else if(runNum>=8100 && runNum<=8246){
			config=4;
		}
	}
	else if(station==3){
		if(runNum>=0 && runNum<=4){
			config=1;
		}
		if(runNum>=470 && runNum<=1448){
			config=2;
		}
		if(runNum>=1449 && runNum<=1901){
			config=1;
		}
		if(runNum>=1902 && runNum<=3103){
			config=5;
		}
		if(runNum>=3104 && runNum<=6004){
			config=3;
		}
		if(runNum>=6005 && runNum<=7653){
			config=4;
		}
		if(runNum>=7658 && runNum<=7808){
			config=3;
		}
	}
	return config;
}

void getWFRMSCutValues(int station, int config, int &vBin, int &hBin, double &vThresh, double &hThresh){

	if(station==2){
		vBin=0;
		hBin=0;
		vThresh=-1.3;
		hThresh=-1.4;
	}
	else if(station==3){
		vBin=0;
		hBin=1;
		if(config==1){
			vThresh=-1.2;
			hThresh=-1.3;
		}
		else if(config==2){
			vThresh=-1.3;
			hThresh=-1.4;
		}
		else if(config==3 || config==4){
			vThresh=-1.0;
			hThresh=-1.1;
		}
		else if(config==5){
			vThresh=-0.7;
			hThresh=-0.8;
		}
	}
}

/*
	input: station, config, pulser, pol, startX, stopX, startY, stopY (the last four are the bounds of the projection region for projecting the TH2D into TH1D)
	output: changes the value of thetas and phis

	function: Returns the bins for plotting the cal pulser cut
				if "inBins" is true, the answer will be reported as a bin number from 0->360/0->180
			if "inBins" is false, the answer will be reported as a number in degrees from -180->180/-90->90
*/

void getCalCutPlotBoundary(int station, int config, int pulser, int pol, bool inBins, int &startX, int &stopX, int &startY, int &stopY){
	if(station==2){
		if(pulser==5){
			if(pol==0){
				if(config==2 || config==4 || config==5){
					startX=-35;
					stopX=-15;
					startY=-35;
					stopY=-15;
				}
			}
		}
		if(pulser==6){
			if(pol==0){
				if(config==1 || config==2 || config==3 || config==4){
					startX=58;
					stopX=72;
					startY=-4;
					stopY=14;
				}
			}
		}

	}
	if(inBins){
		startX+=180;
		stopX+=180;
		startY+=90;
		stopY+=90;
	}
}

/*
	input: station, config, pulser, pol, startPhi, stopPhi, startTheta, stopTheta (the last four are the bounds of the fit region for computing the cut)
	output: changes the value of thetas and phis

	function:  looks up the boundaries for the gaussian fits used to set the cal pulser cut boxes
*/


void getCalFitRange(int station, int config, int pulser, int pol, double &startPhi, double &stopPhi, double &startTheta, double &stopTheta){
	if(station==2){
		if(pulser==5){
			if(pol==0){
				if(config==2 || config==4){
					startPhi=-30.;
					stopPhi=-18.;
					startTheta=-30.;
					stopTheta=-16.;
				}
			}
		}
		if(pulser==6){
			if(pol==0){
				if(config==1 || config==2 || config==3 || config==4){
					startPhi=62.;
					stopPhi=70.;
					startTheta=2.;
					stopTheta=10.;
				}
			}
		}
	}
}

/*
	input: station, pulser (5 or 6), polarization, theta and phi
	output: changes the value of theta and phi to where the cal pulser is

	function: computes the "truth" theta and phi of the cal pulser
*/

void getRealLocation(int station, int pulser, int pol, double &theta, double &phi){

	int cal_ant_index;
	if(pulser==5){
		if(pol==0) cal_ant_index=1;
		else if(pol==1) cal_ant_index=0;
	}
	else if(pulser==6){
		if(pol==0) cal_ant_index=3;
		else if(pol==1) cal_ant_index=4;
	}

	AraGeomTool *araGeom = AraGeomTool::Instance();

	//compute the average depth of the station
	double antenna_average[3]={0.};
	for(int i=0; i<16; i++){
		for(int ii=0; ii<3; ii++){
			antenna_average[ii]+=(araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[ii]);
		}
	}
	for(int ii=0; ii<3; ii++){
		antenna_average[ii]/=16.;
	}
	double X = araGeom->getStationInfo(station)->getCalAntennaInfo(cal_ant_index)->antLocation[0];
	double Y = araGeom->getStationInfo(station)->getCalAntennaInfo(cal_ant_index)->antLocation[1];
	double Z = araGeom->getStationInfo(station)->getCalAntennaInfo(cal_ant_index)->antLocation[2];
	phi = TMath::ATan2(Y-antenna_average[1],X-antenna_average[0])*TMath::RadToDeg();
	double depth_diff = Z-antenna_average[2];
	double horz_dist = sqrt(pow((X-antenna_average[0]),2.0)+pow((Y-antenna_average[1]),2.0));
	theta = TMath::ATan2(depth_diff,horz_dist)*TMath::RadToDeg();
}

void identifyCalPulser(int station, int config, int peakTheta, int peakPhi, bool &isCP5, bool &isCP6, int flex=0){
	if(station==2){
		if (peakPhi >= -30-flex && peakPhi <= -20+flex
		 	&& peakTheta >= -29-flex && peakTheta <= -14+flex)
		{
			isCP5=true;
		}
		if (peakPhi >= 60-flex && peakPhi <= 70+flex
			&& peakTheta >= 0-flex && peakTheta <= 15+flex)
		{
			isCP6=true;
		}
		if(config==5){
			if (peakPhi >= -28-flex && peakPhi <= -18+flex
		 		&& peakTheta >= 30-flex && peakTheta <= 41+flex)
			{
				isCP5=true;
			}
		}
	}
	else if(station==3){
		if (peakPhi >= -28-flex && peakPhi <= -18+flex
			&& peakTheta >= -20-flex && peakTheta <= -5+flex)
		{
			isCP5=true;
		}
		if (peakPhi >= 60-flex && peakPhi <= 70+flex
			&& peakTheta >= -20-flex && peakTheta <= -5+flex)
		{
			isCP6=true;
		}
	}
}



/*
	input: vector of graphs for an event
	output: 0 (has no error), 1 (has less than 550 points error)

	function: checks to see if any waveform is < 550 samples
*/

int hasShortWaveformMultiGraph(vector <TGraph*> grs){
	int event_has_error=0;
	for(int i=0; i<grs.size(); i++){
		if(grs[i]->GetN()<500.) event_has_error++;
	}
	return event_has_error;
}


void GetVectorOfBlockMeans(TGraph *grIn, vector<double> &means){
	int n_input = grIn->GetN();
	double *oldX = grIn->GetX();
	double *oldY = grIn->GetY();

	deque <double> inX;
	deque <double> inY;

	for(int samp=0; samp<n_input; samp++){
		inX.push_back(oldX[samp]);
		inY.push_back(oldY[samp]);
	}

	while(inX.size()>63){
		vector <double> sub_X;
		vector <double> sub_Y;
		int num_to_pop=0;
		for(int samp=0; samp<64; samp++){
			sub_X.push_back(inX[samp]);
			sub_Y.push_back(inY[samp]);
			num_to_pop++;
		}
		double average = std::accumulate( sub_Y.begin(), sub_Y.end(), 0.0)/sub_Y.size();
		means.push_back(average);
		for(int iPop=0; iPop<num_to_pop; iPop++){
			inX.pop_front();
			inY.pop_front();
		}
	}
}

bool hasSpareChannelIssue_v2(vector<TGraph*> electChansGraphs, int station){

	// this is only tuned for A2 currrently (ugh)
	bool hasIssue=false;
	if(station!=2){
		return hasIssue;
	}

	int shortest=300;
	vector<vector< double > > means;
	for(int i=0; i<4; i++){
		vector<double> theseMeans;
		GetVectorOfBlockMeans(electChansGraphs[i], theseMeans);
		if(theseMeans.size()<shortest){
			shortest=theseMeans.size();
		}
		means.push_back(theseMeans);
	}
	for(int i=0; i<4; i++){
		while(means[i].size()>shortest){
			means[i].pop_back();
		}
	}
	int numViolatingBlocks=0;
	for(int block=0; block<shortest; block++){
		int numViolating=0;
		for(int chan=0; chan<4; chan++){
			if(abs(means[chan][block])>20){
				numViolating++;
			}
		}
		if(numViolating>1){
			numViolatingBlocks++;
		}
		if(numViolatingBlocks>0){
			break; //that's enough to call it bad under our current criteria
		}
	}

	if(numViolatingBlocks>0){
		hasIssue=true;
	}
	return hasIssue;
}


/*
	input: vector of graphs of the electronics channels
	output: 0 (has no error), 1 (has spare channel offset error)

	function: checks to see if any two spare channels have RMS> 20 or any one channel has RMS>60
*/

bool hasSpareChannelIssue(vector<TGraph*> electChansGraphs){
	vector<double> spareRMS;
	for(int i=0; i<4; i++){
		spareRMS.push_back(electChansGraphs[i]->GetRMS(2));
		// cout<<"Spare RMS is "<<spareRMS[i]<<endl;
	}
	int numBadSpareChans=0;
	int numReallyBadSpareChans=0;
	for(int i=0; i<4; i++){
		if(spareRMS[i]>20 && i!=3){
			numBadSpareChans++;
			// cout<<"Bad spare chan in chan "<<i<<endl;
		}
		if(spareRMS[i]>60 && i!=3){
			numReallyBadSpareChans++;
			// cout<<"Actually, it's a really bad spare chan issue in chan "<<i<<endl;
		}
	}
	bool hasBadSpareChansIssue=false;
	if(numBadSpareChans>1 || numReallyBadSpareChans>0){
		hasBadSpareChansIssue=true;
	}
	// printf("Num bad %d and very bad %d channels \n", numBadSpareChans, numReallyBadSpareChans);
	return hasBadSpareChansIssue;
}

/*
	input: station, unixTime (UTC)
	output: 0 (is good time), 1 (is bad time)

	function: checks if the livetime is bad
*/

int isBadLivetime(int station, int unixTime){

	bool isBadLivetime=false;
	if(station==2){

		if(
			// Livetime flagged as bad by my me
			(unixTime>=1389381600 && unixTime<=1389384000) || // from run 2868
			(unixTime>=1420317600 && unixTime<=1420318200) || // from run 4775
			// (unixTime>=1449189600 && unixTime<=1449190200) || // from run 6507
			(unixTime>=1449187200 && unixTime<=1449196200) || // from run 6507

			//Livetime flagged as bad by my undergrads
			//config 1
			// (unixTime>=1380234000 && unixTime<=1380236400) || // from run 2428 22 hour balloon launch
			// (unixTime>=1382046000 && unixTime<=1382047500) || // from run 2536 22 hour balloon launch
			(unixTime>=1382712900 && unixTime<=1382713500) || // from run 2575
			(unixTime>=1382972700 && unixTime<=1382973300) || // from run 2589
			// (unixTime>=1383689100 && unixTime<=1383690900) || // from run 2631 22 hour balloon launch
			(unixTime>=1383884400 && unixTime<=1383886200) || // from run 2642
			(unixTime>=1384060200 && unixTime<=1384061100) || // from run 2652
			(unixTime>=1384487400 && unixTime<=1384489800) || // from run 2677
			(unixTime>=1384489980 && unixTime<=1384491060) || // from run 2678 at start may be glitch or continued from 2677
			(unixTime>=1384856520 && unixTime<=1384856640) || // from run 2698 super zoomed in two minute window
			// (unixTime>=1385674200 && unixTime<=1385675100) || // from run 2744 22 hour balloon launch
			(unixTime>=1389381600 && unixTime<=1389383700) || // from run 2868 first of two from run 2868
			(unixTime>=1389398700 && unixTime<=1389400200) || // from run 2868 second of two from run 2868
			(unixTime>=1389665100 && unixTime<=1389666300) || // from run 2884
			(unixTime>=1393288800 && unixTime<=1393289400) || // from run 3099
			// (unixTime>=1397856600 && unixTime<=1397858400) || // from run 3442 22 hour balloon launch

			//config 2
			(unixTime>=1376731800 && unixTime<=1376733000) || // from run 2235

			//conifg 3
			(unixTime>=1400276700 && unixTime<=1400277300) || // from run 3605 mainly looks like glitch at end

			//config 4
			(unixTime>=1409986500 && unixTime<=1409988000) || // from run 4184
			// (unixTime>=1412026200 && unixTime<=1412027100) || // from run 4301 22 hr balloon
			// (unixTime>=1412285400 && unixTime<=1412288100) || // from run 4316 weird 22hr balloon
			// (unixTime>=1412544600 && unixTime<=1412545500) || // from run 4331 22hr balloon
			// (unixTime>=1412803800 && unixTime<=1412804700) || // from run 4346 22hr balloon
			(unixTime>=1413898200 && unixTime<=1413899100) || // from run 4408
			(unixTime>=1414083900 && unixTime<=1414086000) || // from run 4418
			(unixTime>=1414350300 && unixTime<=1414351200) || // from run 4434 pt 1
			// (unixTime>=1414358700 && unixTime<=1414359900) || // from run 4434 pt 2 22hr balloon
			(unixTime>=1414674300 && unixTime<=1414674780) || // from run 4452
			(unixTime>=1414986600 && unixTime<=1414987200) || // from run 4471
			(unixTime>=1415223000 && unixTime<=1415223900) || // from run 4483
			(unixTime>=1415380500 && unixTime<=1415381400) || // from run 4493
			(unixTime>=1415558100 && unixTime<=1415559000) || // from run 4503
			(unixTime>=1415742300 && unixTime<=1415743800) || // from run 4513
			(unixTime>=1416207000 && unixTime<=1416212100) || // from run 4541
			(unixTime>=1420978200 && unixTime<=1420978800) || // from run 4814
			(unixTime>=1416905100 && unixTime<=1416910500) || // from run 4579 two spikes about an hour apart
			// (unixTime>=1416950700 && unixTime<=1416951600) || // from run 4582 22 hour balloon launch
			(unixTime>=1417677000 && unixTime<=1417678200) || // from run 4621  weird and cool
			(unixTime>=1417836000 && unixTime<=1417837500) || // from run 4631
			(unixTime>=1420097100 && unixTime<=1420098300) || // from run 4763
			(unixTime>=1420293300 && unixTime<=1420294200) || // from run 4774
			(unixTime>=1420317600 && unixTime<=1420318200) || // from run 4775
			(unixTime>=1420978200 && unixTime<=1420978800) || // from run 4814
			(unixTime>=1421024400 && unixTime<=1421025300) || // from run 4817
			(unixTime>=1421713200 && unixTime<=1421718600) || // from run 4872 looks full of errors and not spiky but could have a spiky
			(unixTime>=1421718000 && unixTime<=1421725800) || // from run 4873 definitely an error but also has spiky boy, part 1 of 2
			(unixTime>=1421733300 && unixTime<=1421733900) || // from run 4873 spiky boy alone but in a run with errors, part 2 of 2
			(unixTime>=1421783400 && unixTime<=1421794200) || // from run 4876 definitely an error but not a spikey boy
			// (unixTime>=1428529800 && unixTime<=1428530700) || // from run 5389 22 hour balloon launch
			(unixTime>=1435623000 && unixTime<=1435623600) || // from run 5801
			// (unixTime>=1436394000 && unixTime<=1436395200) || // from run 5845 22 hour balloon launch
			(unixTime>=1437601200 && unixTime<=1437602700) || // from run 5915 looks like error at the start
			// (unixTime>=1439933700 && unixTime<=1439934960) || // from run 6048 22 hour balloon launch
			(unixTime>=1440581700 && unixTime<=1440582480) || // from run 6086
			// (unixTime>=1441489200 && unixTime<=1441490280) || // from run 6137 22 hour balloon launch
			// (unixTime>=1444685400 && unixTime<=1444687080) || // from run 6322 22 hour balloon launch
			// (unixTime>=1445722020 && unixTime<=1445723220) || // from run 6383 22 hour balloon launch
			(unixTime>=1445934900 && unixTime<=1445935500) || // from run 6396
			(unixTime>=1445960400 && unixTime<=1445961000) || // from run 6397
			// (unixTime>=1445982120 && unixTime<=1445982900) || // from run 6398 22 hour balloon launch
			(unixTime>=1446165600 && unixTime<=1446166200) || // from run 6408
			// (unixTime>=1446327300 && unixTime<=1446328200) || // from run 6418 22 hour balloon launch
			(unixTime>=1446607800 && unixTime<=1446608640) || // from run 6433 looks like an error at end
			(unixTime>=1446784200 && unixTime<=1446784800) || // from run 6445
			// (unixTime>=1476739800 && unixTime<=1476741000) || // from run 8100 22 hour balloon launch
			// (unixTime>=1476999000 && unixTime<=1476999900) || // from run 8114 22 hour balloon launch but barely noticeable
			// (unixTime>=1477258200 && unixTime<=1477259100) || // from run 8129 22 hour balloon launch
			(unixTime>=1477511700 && unixTime<=1477512600) || // from run 8143 weird possible balloon launch
			(unixTime>=1477950300 && unixTime<=1477951500) || // from run 8168 22 hour balloon launch
			// (unixTime>=1478033400 && unixTime<=1478034000) || // from run 8173 22 hour balloon launch
			// (unixTime>=1478295300 && unixTime<=1478296200) || // from run 8188 22 hour balloon launch
			// (unixTime>=1478728500 && unixTime<=1478729400) || // from run 8213 22 hour balloon launch
			(unixTime>=1479231900 && unixTime<=1479232500) || // from run 8241

			// config 5
			(unixTime>=1449280500 && unixTime<=1449281100) || // from run 6513
			(unixTime>=1449610200 && unixTime<=1449612000) || // from run 6531
			(unixTime>=1450536000 && unixTime<=1450537200) || // from run 6584
			// (unixTime>=1450906200 && unixTime<=1450907100) || // from run 6606    22hr
			// (unixTime>=1451423700 && unixTime<=1451424600) || // from run 6635   22hr
			(unixTime>=1452008100 && unixTime<=1452009000) || // from run 6669
			// (unixTime>=1452115800 && unixTime<=1452116700) || // from run 6675    22hr
			(unixTime>=1452197700 && unixTime<=1452198600) || // from run 6679
			(unixTime>=1452213600 && unixTime<=1452214200) || // from run 6680
			(unixTime>=1452282000 && unixTime<=1452282600) || // from run 6684
			(unixTime>=1452298200 && unixTime<=1452298800) || // from run 6685    possible error
			(unixTime>=1452385500 && unixTime<=1452386400) || // from run 6690
			// (unixTime>=1452457800 && unixTime<=1452458700) || // from run 6694   22 hr
			(unixTime>=1452494100 && unixTime<=1452495000) || // from run 6696   possible error
			// (unixTime>=1452545100 && unixTime<=1452545880) || // from run 6700    could be error or 22hr
			// (unixTime>=1452636900 && unixTime<=1452637500) || // from run 6705   could be error or 22hr
			(unixTime>=1452715200 && unixTime<=1452716100) || // from run 6709   possible error
			(unixTime>=1452972300 && unixTime<=1452973440) || // from run 6724   possible error
			// (unixTime>=1453325400 && unixTime<=1453326600) || // from run 6743   22 hr
			(unixTime>=1453408500 && unixTime<=1453409400) || // from run 6747
			(unixTime>=1453930200 && unixTime<=1453931400) || // from run 6776
			// (unixTime>=1454535000 && unixTime<=1454536500) || // from run 6818   22 hr
			// (unixTime>=1455746400 && unixTime<=1455747900) || // from run 6889   22 hr
			(unixTime>=1456200900 && unixTime<=1456201800) || // from run 6916
			(unixTime>=1456392600 && unixTime<=1456393800) || // from run 6927
			(unixTime>=1456997400 && unixTime<=1456999200) || // from run 6962
			// (unixTime>=1457559000 && unixTime<=1457560800) || // from run 6994   22 hr
			(unixTime>=1460842800 && unixTime<=1460844600) || // from run 7119   22 hr // has CW contam cal pulsers
			// (unixTime>=1461620100 && unixTime<=1461621900) || // from run 7161   22 hr
			(unixTime>=1463002200 && unixTime<=1463004000) || // from run 7243  22 hr // has CW contam cal pulsers
			(unixTime>=1466501400 && unixTime<=1466503200) || // from run 7474
			(unixTime>=1466721900 && unixTime<=1466724600) || // from run 7486 22 hr // has CW contam cal pulsers
			(unixTime>=1466805600 && unixTime<=1466808300) || // from run 7489 22 hr // has CW contam cal pulsers
			(unixTime>=1466890200 && unixTime<=1466892000) || // from run 7494   22 hr // has CW contam cal pulsers
			(unixTime>=1467927600 && unixTime<=1467929700) || // from run 7552   22 hr
			// (unixTime>=1472333400 && unixTime<=1472335200) || // from run 7831   22 hr
			(unixTime>=1473111300 && unixTime<=1473112800) || // from run 7879    22 hr // has CW contam cal
			// (unixTime>=1473370500 && unixTime<=1473372900) || // from run 7899   22 hr
			// (unixTime>=1475011500 && unixTime<=1475013600) || // from run 7993   22 hr
			(unixTime>=1475185200 && unixTime<=1475187900) || // from run 8003 balloon 22hr // has CW contam cal pulsers
			// (unixTime>=1475358000 && unixTime<=1475359800) || // from run 8013 balloon 22hr
			(unixTime>=1475529900 && unixTime<=1475531400) || // from run 8023 balloon 22hr // has CW contam cal pulsers
			// (unixTime>=1475702700 && unixTime<=1475704200) || // from run 8033 balloon 22hr
			(unixTime>=1476221400 && unixTime<=1476222300) // from run 8069 balloon 22hr // has CW contam cal pulsers
			// (unixTime>=1476479700 && unixTime<=1476481800) || // from run 8084 balloon 22hr

		)
			{
				isBadLivetime=true;
		}
	}
	else if(station==3){

		if(

			// config 1 from undergrads
			(unixTime>=1380234300 && unixTime<=1380235500) || // from run 1538, 22 hour balloon launch
			(unixTime>=1381008600 && unixTime<=1381010400) || // from run 1584, 22 hour balloon launch
			(unixTime>=1382476200 && unixTime<=1382477400) || // from run 1670, 22 hour balloon launch-ish
			(unixTime>=1382687400 && unixTime<=1382688600) || // from run 1682
			(unixTime>=1382712600 && unixTime<=1382713800) || // from run 1684, 15 hour spike
			(unixTime>=1382972700 && unixTime<=1382973300) || // from run 1698, 15 hour spike
			(unixTime>=1383688800 && unixTime<=1383691500) || // from run 1739, 22 hour balloon launch
			(unixTime>=1384060200 && unixTime<=1384060800) || // from run 1761
			(unixTime>=1384208700 && unixTime<=1384209900) || // from run 1770, 22 hour balloon launch
			(unixTime>=1384486200 && unixTime<=1384492800) || // from run 1786, repeated bursts over ~2 hrs
			(unixTime>=1389399600 && unixTime<=1389400800) || // from run 1980
			(unixTime>=1389744000 && unixTime<=1389747600) || // from run 2001, lots of activity, sweeps in phi
			(unixTime>=1390176600 && unixTime<=1390182000) || // from run 2025
			(unixTime>=1391027700 && unixTime<=1391028900) || // from run 2079, 22 hour balloon launch, but early?
			(unixTime>=1393652400 && unixTime<=1393660800) || // from run 2235, repeated bursts over ~2 hrs
			(unixTime>=1394846400 && unixTime<=1394856000) || // from run 2328, repeated bursts over ~2.5 hours
			(unixTime>=1395437400 && unixTime<=1395438600) || // from run 2363, 22 hour balloon launch
			(unixTime>=1397856300 && unixTime<=1397857800) || // from run 2526, 22 hour balloon launch

			// config 2
			(unixTime>=1390176600 && unixTime<=1390182000) || // from run 3533

			// config 3
			(unixTime>=1409954100 && unixTime<=1409956200) || // from run 3216, 22 hour balloon launch
			(unixTime>=1409986800 && unixTime<=1409988600) || // from run 3217
			(unixTime>=1412026200 && unixTime<=1412028000) || // from run 3332
			(unixTime>=1412284920 && unixTime<=1412287020) || // from run 3347, 22 hour balloon launch
			(unixTime>=1412544120 && unixTime<=1412546400) || // from run 3362, 22 hour balloon launch
			(unixTime>=1412803620 && unixTime<=1412805780) || // from run 3377, 22 hour balloon launch
			(unixTime>=1413897900 && unixTime<=1413899100) || // from run 3439
			(unixTime>=1413914400 && unixTime<=1413922200) || // from run 3440 big wide weird above ground
			(unixTime>=1414083600 && unixTime<=1414086300) || // from run 3449 , 2 spikes
			(unixTime>=1413550800 && unixTime<=1413552600) || // from run 3419, end of the run, before a software dominated run starts
			(unixTime>=1414674000 && unixTime<=1414675500) || // from run 3478
			(unixTime>=1415380500 && unixTime<=1415381400) || // from run 3520
			(unixTime>=1415460600 && unixTime<=1415461500) || // from run 3524
			(unixTime>=1415742000 && unixTime<=1415744100) || // from run 3540 22hr balloon
			(unixTime>=1416207300 && unixTime<=1416209700) || // from run 3568 2 small spikes
			(unixTime>=1416457800 && unixTime<=1416459000) || // from run 3579
			(unixTime>=1416909600 && unixTime<=1416910680) || // from run 3605
			(unixTime>=1416951000 && unixTime<=1416952500) || // from run 3608 22hr balloon
			(unixTime>=1417676400 && unixTime<=1417679400) || // from run 3647
			(unixTime>=1417742400 && unixTime<=1417743600) || // from run 3651
			(unixTime>=1417836600 && unixTime<=1417839300) || // from run 3656
			(unixTime>=1420317000 && unixTime<=1420318200) || // from run 3800
			(unixTime>=1420493700 && unixTime<=1420494600) || // from run 3810 22hr balloon
			(unixTime>=1420513200 && unixTime<=1420515000) || // from run 3811
			(unixTime>=1420598700 && unixTime<=1420600500) || // from run 3816
			(unixTime>=1420857900 && unixTime<=1420859700) || // from run 3830
			(unixTime>=1421019000 && unixTime<=1421020200) || // from run 3840 22hr balloon maybe?
			(unixTime>=1421101800 && unixTime<=1421103600) || // from run 3863 22hr balloon
			(unixTime>=1421723400 && unixTime<=1421723940) || // from run 3910
			(unixTime>=1421750700 && unixTime<=1421751720) || // from run 3912
			(unixTime>=1421868600 && unixTime<=1421881200) || // from run 3977 looks intentional
			(unixTime>=1421881200 && unixTime<=1421884680) || // from run 3978 continuation of thing above
			(unixTime>=1422048900 && unixTime<=1422049800) || // from run 3987 , 22 hour balloon launch
			(unixTime>=1422307200 && unixTime<=1422308100) || // from run 3995 22hr balloon
			(unixTime>=1423660800 && unixTime<=1423661700) || // from run 4132
			(unixTime>=1424819880 && unixTime<=1424820720) || // from run 4200
			(unixTime>=1428529500 && unixTime<=1428531000) || // from run 4412, 22 hour balloon launch
			(unixTime>=1429094400 && unixTime<=1429095600) || // from run 4445
			(unixTime>=1429615800 && unixTime<=1429617600) || // from run 4473
			(unixTime>=1429616700 && unixTime<=1429627500) || // from run 4474
			(unixTime>=1429733400 && unixTime<=1429734600) || // from run 4482
			(unixTime>=1431034500 && unixTime<=1431036900) || // from run 4557 , 22 hour balloon launch
			(unixTime>=1433365500 && unixTime<=1433367900) || // from run 4693
			(unixTime>=1435755600 && unixTime<=1435756500) || // from run 4829
			(unixTime>=1435791000 && unixTime<=1435791600) || // from run 4832
			(unixTime>=1436393700 && unixTime<=1436395500) || // from run 4867
			(unixTime>=1476740100 && unixTime<=1476741300) || // from run 7658
			(unixTime>=1477511400 && unixTime<=1477518300) || // from run 7704, big spike followed by nothing at all
			(unixTime>=1477604700 && unixTime<=1477605900) || // from run 7709,  22 hour balloon launch
			(unixTime>=1477950300 && unixTime<=1477951500) || // from run 7729
			(unixTime>=1479231600 && unixTime<=1479235800) || // from run 7802  , big spike followed by nothing at all


			//config 4
			(unixTime>=1448959200 && unixTime<=1448960100) || // from run 6009
			(unixTime>=1449610500 && unixTime<=1449611400) || // from run 6046 22 hour balloon launch
			(unixTime>=1450119900 && unixTime<=1450120500) || // from run 6077 possible 22 hour balloon launch
			(unixTime>=1450536360 && unixTime<=1450536720) || // from run 6098 spike is at end of time
			(unixTime>=1452116100 && unixTime<=1452116700) || // from run 6188 end of time and possible balloon launch
			(unixTime>=1452196800 && unixTime<=1452198600) || // from run 6193 could be balloon
			(unixTime>=1452213600 && unixTime<=1452214200) || // from run 6194
			(unixTime>=1452282300 && unixTime<=1452282900) || // from run 6198 could be balloon
			(unixTime>=1452298500 && unixTime<=1452299100) || // from run 6199 spike is at end of measured time
			(unixTime>=1452385800 && unixTime<=1452386400) || // from run 6203 spike is at end of measured time
			(unixTime>=1452457800 && unixTime<=1452458700) || // from run 6206 spike is at end of measured time, could be balloon
			(unixTime>=1452494100 && unixTime<=1452494700) || // from run 6208 spike is at end of measured time
			(unixTime>=1452544980 && unixTime<=1452545580) || // from run 6212 could be balloon
			(unixTime>=1452561120 && unixTime<=1452561480) || // from run 6213 spike is at end of measured time
			(unixTime>=1452637020 && unixTime<=1452637260) || // from run 6219 spike is at end of measured time, could be balloon
			(unixTime>=1452715320 && unixTime<=1452715680) || // from run 6223 spike is at end of measured time
			(unixTime>=1452972660 && unixTime<=1452973020) || // from run 6239 spike is at end of measured time
			(unixTime>=1453325400 && unixTime<=1453326300) || // from run 6259 could be balloon
			(unixTime>=1453930500 && unixTime<=1453931100) || // from run 6295 could be balloon
			(unixTime>=1454535000 && unixTime<=1454536200) || // from run 6328 could be balloon
			(unixTime>=1454911200 && unixTime<=1454911800) || // from run 6349 spike is at end of measured time could match below
			(unixTime>=1454911200 && unixTime<=1454912100) || // from run 6350 spike is at start of measured time could match above
			(unixTime>=1455746400 && unixTime<=1455747300) || // from run 6397 could be balloon
			(unixTime>=1456374300 && unixTime<=1456374900) || // from run 6433
			(unixTime>=1457559300 && unixTime<=1457560500) || // from run 6501 could be balloon
			(unixTime>=1460843100 && unixTime<=1460844600) || // from run 6618 spike is at start of measured time, could be balloon
			(unixTime>=1467927840 && unixTime<=1467929640) || // from run 7052 could be balloon
			(unixTime>=1473371280 && unixTime<=1473372180) || // from run 7458 could be balloon
			(unixTime>=1475186100 && unixTime<=1475187000) || // from run 7562 could be balloon
			(unixTime>=1475530500 && unixTime<=1475531700) || // from run 7584 could be balloon
			(unixTime>=1476221400 && unixTime<=1476222600) // from run 7625 could be balloon

		)
			{
				isBadLivetime=true;

		}


	}
	return isBadLivetime;

}
/*
	input: station
	output:  vector of bad runs

	function: runs vector of ints of runs
				the we decided had surface activity
				during unblinding process
*/
vector<int> BuildSurfaceRunList(int station){


	vector<int> exclude;
	if(station==2){

		/*
		Runs shared with Ming-Yuan
			http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1889
		*/
		exclude.push_back(2090);
		exclude.push_back(2678);
		exclude.push_back(4777);
		exclude.push_back(5516);
		exclude.push_back(5619);
		exclude.push_back(5649);
		exclude.push_back(5664);
		exclude.push_back(5666);
		exclude.push_back(5670);
		exclude.push_back(5680);
		exclude.push_back(6445);
		exclude.push_back(6536);
		exclude.push_back(6542);
		exclude.push_back(6635);
		exclude.push_back(6655);
		exclude.push_back(6669);
		exclude.push_back(6733);

		/*
		Runs identified independently
		*/
		exclude.push_back(2091);
		exclude.push_back(2155);
		exclude.push_back(2636);
		exclude.push_back(2662);
		exclude.push_back(2784);
		exclude.push_back(4837);
		exclude.push_back(4842);
		exclude.push_back(5675);
		exclude.push_back(5702);
		exclude.push_back(6554);
		exclude.push_back(6818);
		exclude.push_back(6705);
		exclude.push_back(8074);

	}
	else if(station==3){

		/*
		Runs shared with Ming-Yuan
			http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=2041
		*/
		exclude.push_back(977);
		exclude.push_back(1240);
		exclude.push_back(3158);
		exclude.push_back(3431);
		exclude.push_back(3432);
		exclude.push_back(3435);
		exclude.push_back(3437);
		exclude.push_back(3438);
		exclude.push_back(3439);
		exclude.push_back(3440);
		exclude.push_back(3651);
		exclude.push_back(3841);
		exclude.push_back(4472);
		exclude.push_back(4963);
		exclude.push_back(4988);
		exclude.push_back(4989);


		/*
		Runs identified independently
		*/
		exclude.push_back(1745);
		exclude.push_back(3157);
		exclude.push_back(3652);
		exclude.push_back(3800);
		exclude.push_back(6193);
		exclude.push_back(6319);
		exclude.push_back(6426);

		/*
		Runs I am sure we will exclude...
		*/
		exclude.push_back(2000);
		exclude.push_back(2001);


	}

	return exclude;
}

/*
	input: station
	output:  vector of bad runs

	function: returns vector of ints of bad run numbers
*/
vector<int> BuildBadRunList(int station){

	vector<int> exclude;
	if(station==2){


		/*2013*/

		/*2014*/

			/*
			2014 rooftop pulsing
				http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
			*/
			exclude.push_back(3120);
			exclude.push_back(3242);


			/*
			2014 surface pulsing
				originally flagged by 2884, 2895, 2903, 2912, 2916
				going to throw all runs jan 14-20 (<1 week of livetime, oh well...)
			*/
			exclude.push_back(2884); //jan 14 2014 surface pulser runs //actual problem causer
				exclude.push_back(2885); //jan 16 2014 surface pulser runs //exclusion by proximity
				exclude.push_back(2889); //jan 16 2014 surface pulser runs //exclusion by proximity
				exclude.push_back(2890); //jan 16 2014 surface pulser runs //exclusion by proximity
				exclude.push_back(2891); //jan 16 2014 surface pulser runs //exclusion by proximity
				exclude.push_back(2893); //jan 16 2014 surface pulser runs //exclusion by proximity
			exclude.push_back(2895); //jan 16 2014 surface pulser runs //actual problem causer
				exclude.push_back(2898); //jan 16 2014 surface pulser runs //exclusion by proximity
				exclude.push_back(2900); //jan 17 2014 surface pulser runs //exclusion by proximity
				exclude.push_back(2901); //jan 17 2014 surface pulser runs //exclusion by proximity
				exclude.push_back(2902); //jan 17 2014 surface pulser runs //exclusion by proximity
			exclude.push_back(2903); //jan 18 2014 surface pulser runs //actual problem causer
				exclude.push_back(2905); //jan 18 2014 surface pulser runs //exclusion by proximity
				exclude.push_back(2906); //jan 18 2014 surface pulser runs //exclusion by proximity
				exclude.push_back(2907); //jan 18 2014 surface pulser runs //exclusion by proximity
			exclude.push_back(2912); //jan 19 2014 surface pulser runs //actual problem causer
				exclude.push_back(2915); //jan 18 2014 surface pulser runs //exclusion by proximity
			exclude.push_back(2916); //jan 20 2014 surface pulser runs //actual problem causer
				exclude.push_back(2918); //jan 20 2014 surface pulser runs

			exclude.push_back(2938); //surface pulsing from m richman (identified by MYL http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1889 slide 14)
			exclude.push_back(2939); //surface pulsing from m richman (identified by MYL http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1889 slide 14)

			/*
			2014 Cal pulser sweep
				http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
			*/
			for(int i=3139; i<=3162; i++){ exclude.push_back(i); }
			for(int i=3164; i<=3187; i++){ exclude.push_back(i); }
			for(int i=3289; i<=3312; i++){ exclude.push_back(i); }

			/*
			2014 L2 Scaler Masking Issue
				Cal pulsers sysemtatically do not reconstruct correctly, rate is only 1 Hz
				Excluded because configuration was not "science good"
			*/
			for(int i=3464; i<=3504; i++){ exclude.push_back(i); }

			/*
			2014 Trigger Length Window Sweep
				http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
			*/
			for(int i=3578; i<=3598; i++){ exclude.push_back(i); }




		/*2015*/
			exclude.push_back(4004);
			/*
			2015 icecube deep pulsing
				4787 is the "planned" run
				4795,4797-4800 were accidental
			*/
			exclude.push_back(4785); //accidental deep pulser run (http://ara.physics.wisc.edu/docs/0017/001719/003/181001_ARA02AnalysisUpdate.pdf, slide 38)
			exclude.push_back(4787); //deep pulser run (http://ara.physics.wisc.edu/docs/0017/001724/004/181015_ARA02AnalysisUpdate.pdf, slide 29)
			for(int i=4795; i<=4800; i++){ exclude.push_back(i); }

			/*
			2015 noise source tests
				January 2015
				http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2015
			*/
			for(int i=4820; i<=4825; i++){ exclude.push_back(i); }
			for(int i=4850; i<=4854; i++){ exclude.push_back(i); }
			for(int i=4879; i<=4936; i++){ exclude.push_back(i); }
			for(int i=5210; i<=5277; i++){ exclude.push_back(i); }

			/*
			2015 surface pulsing
				January 2015
				http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1339 (slide 5)
			*/
			exclude.push_back(4872);
			exclude.push_back(4873);
			exclude.push_back(4876); // Identified by MYL http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1889 slide 14

			/*
			2015 Pulser Lift
				December 2015
				http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1269 (page 2)
				Run number from private communication with John Kelley
			*/
			exclude.push_back(6513);

			/*
			2015 ICL pulsing
				December 2015
				http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1269 (page 7)
			*/
			exclude.push_back(6527);


		/*2016*/

			/*
			2016 cal pulser sweep
				January 2015
				http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2016
			*/
			for(int i=7625; i<=7686; i++){ exclude.push_back(i); }


		/*Other*/

			/*
			D1 Glitches
				Identified by MYL as having glitches after long periods of downtime
			*/
			exclude.push_back(3);
			exclude.push_back(11);
			exclude.push_back(59);
			exclude.push_back(60);
			exclude.push_back(71);

			/*
			Badly misreconstructing runs
			*/
			for(int i=8100; i<=8246; i++){
				exclude.push_back(i);
			}

	}
	else if(station==3){

		/*2013*/

			// /*
			// 	Misc tests: http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2013
			// */
			// for(int i=22; i<=62; i++){ exclude.push_back(i); }



			// /*
			// 	ICL rooftop
			// 	http://ara.icecube.wisc.edu/wiki/index.php/A23_Diffuse_UW
			// */

			// for(int i=63; i<=70; i++){ exclude.push_back(i); }
			// for(int i=333; i<=341; i++){ exclude.push_back(i); }


			// /*
			// 	Cal sweep
			// 	http://ara.icecube.wisc.edu/wiki/index.php/A23_Diffuse_UW
			// */

			// for(int i=72; i<=297; i++){ exclude.push_back(i); }
			// for(int i=346; i<=473; i++){ exclude.push_back(i); }

			/*
				Eliminate all early data taking (all runs before 508)
			*/
			for(int i=0; i<=508; i++){ exclude.push_back(i); }


			/*
				Cal sweep
				http://ara.icecube.wisc.edu/wiki/index.php/A23_Diffuse_UW
			*/


		/*2014*/

			/*
			2014 Rooftop Pulser
				http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
			*/
			exclude.push_back(2235);
			exclude.push_back(2328);

			/*
			2014 Cal Pulser Sweep
				http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
			*/
			for(int i=2251; i<=2274; i++){ exclude.push_back(i); }
			for(int i=2376; i<=2399; i++){ exclude.push_back(i); }

		/*2015*/

			/*
			2015 surface or deep pulsing
				got through cuts
				happened jan 5-6, some jan 8
				waveforms clearly show double pulses or things consistent with surface pulsing
			*/
			exclude.push_back(3811); //deep pulser run
				exclude.push_back(3810); //elminated by proximity to deep pulser run
				exclude.push_back(3820); //elminated by proximity to deep pulser run
				exclude.push_back(3821); //elminated by proximity to deep pulser run
				exclude.push_back(3822); //elminated by proximity to deep pulser run
			exclude.push_back(3823); //deep pulser, observation of 10% iterator event numbers 496, 518, 674, 985, 1729, 2411

			/*
			2015 noise source tests
				January 2015
				http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2015
			*/
			for(int i=3844; i<=3860; i++){ exclude.push_back(i); }
			for(int i=3881; i<=3891; i++){ exclude.push_back(i); }
			for(int i=3916; i<=3918; i++){ exclude.push_back(i); }
			for(int i=3920; i<=3975; i++){ exclude.push_back(i); }
			for(int i=4009; i<=4073; i++){ exclude.push_back(i); }

			/*
			2015 surface pulsing
				January 2015
				http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1339 (slide 5)
			*/
			exclude.push_back(3977);
			exclude.push_back(3978);

			/*
			2015 ICL pulsing
				December 2015
				http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1269 (page 7)
			*/
			exclude.push_back(6041);

			/*
			2015 station anomaly
				see moni report: http://ara.physics.wisc.edu/cgi-bin/DocDB/ShowDocument?docid=1213
				identified by MYL: http://ara.icecube.wisc.edu/wiki/index.php/A23_Diffuse_UW
			*/
			for(int i=4914; i<=4960; i++){ exclude.push_back(i); }


		/*2016*/

			/*
			2016 Cal Pulser Sweep
				http://ara.icecube.wisc.edu/wiki/index.php/Run_Log_2014
			*/
			for(int i=7126; i<=7253; i++){ exclude.push_back(i); }

			/*
				More events with no RF/deep triggers
				seems to precede coming test
			*/
			exclude.push_back(7125);
	}
	return exclude;
}

/*
	input: station, run number
	output: 0 (is good run), 1 (is bad run)

	function: looks through list of known "bad runs"
			reports if the run you gave it is bad
*/
int isBadRun(int station, int run_num, vector<int>BadRunList){
	if(BadRunList.size()<1){
		BadRunList=BuildBadRunList(station);
	}
	int found = (std::find(BadRunList.begin(), BadRunList.end(), run_num) != BadRunList.end());
	return found;
}

/*
	input: path to "the a23_analysis_tools" directory, station number, config number, and run number
	output: 0/false (does not contained untagged cal pulse), 1/true (does contain untagged cal pulse)

	function: looks through list of runs that are known to have <3% cal pulser content
				which we interpret to mean "has no tagged cal pulsers"
				and return if this run for this station and config is such
*/
bool isSoftwareDominatedRun(std::string pathToToolsDir, int station, int run_num){
	char filename[200];
	sprintf(filename,"%s/data/A%i_software_dominated_list.txt",pathToToolsDir.c_str(), station);
	ifstream infile(filename);
	string line;
	string str;

	bool isSoftwareDominated=false;
	//  Read the file
	while (getline(infile, line))
	{   istringstream ss (line);
		vector <string> record;

		while (getline(ss, str, ','))
		record.push_back(str);
		int runNum;
		std::stringstream(record[0]) >> runNum;
		if (runNum==run_num){
			isSoftwareDominated=true;
		}
	}

	return isSoftwareDominated;
}


/*
	input: path to "the a23_analysis_tools" directory, station number, config number, and run number
	output: 0/false (does not contained untagged cal pulse), 1/true (does contain untagged cal pulse)

	function: looks through list of runs that are known to have <3% cal pulser content
				which we interpret to mean "has no tagged cal pulsers"
				and return if this run for this station and config is such
*/
bool hasUntaggedCalpul(std::string pathToToolsDir, int station, int config, int run_num){
	char filename[200];
	sprintf(filename,"%s/data/A%d_c%i_untagged_calpul.csv",pathToToolsDir.c_str(), station,config);
	ifstream infile(filename);
	string line;
	string str;

	bool hasUntagged=false;

	//  Read the file
	while (getline(infile, line))
	{   istringstream ss (line);
		vector <string> record;

		while (getline(ss, str, ','))
		record.push_back(str);
		// cout << record[0] << endl;
		int runNum;
		std::stringstream(record[0]) >> runNum;
		if (runNum==run_num){
			// cout << "Untagged" << endl;
			hasUntagged=true;
		}
	}

	return hasUntagged;
}

Double_t getPeakSqVal(TGraph *gr, int *index){
  Double_t x,y;
  gr->GetPoint(0,x,y);
  Double_t peakVal=y*y;
  Int_t peakBin=0;
  Int_t numPoints=gr->GetN();
  for(int i=1;i<numPoints;i++) {
     gr->GetPoint(i,x,y);
     if(peakVal<y*y) {
        peakVal=y*y;
        peakBin=i;
     }
  }
  if(index) *index=peakBin;
  return peakVal;
 }

//isSpikeyStringEvent
 bool isSpikeyStringEvent(int stationId, bool dropARA03D4, vector <TGraph*> wf, int config){
	 double spikeyRatio = 0.;
	 if (stationId != 3){ //Only A3 has such problem
		 return false;
	 }
	 double maxV[16];
	 int maxVIdx;
	 std::fill(&maxV[0], &maxV[16], 0.);
	 double avgSNR[4];
	 for(int ch=0; ch<16; ch++){
		 maxV[ch] = sqrt(getPeakSqVal(wf[ch], &maxVIdx));
	 }
	 for (int string=0; string<4; string++){
		 avgSNR[string] = (maxV[string] + maxV[string+4] + maxV[string+12])/3.; //TH don't seem to see the spike, so exclude
	 }
	 if(dropARA03D4) spikeyRatio = 2. * avgSNR[0] / (avgSNR[1] + avgSNR[2]);
	 else            spikeyRatio = 3. * avgSNR[0] / (avgSNR[1] + avgSNR[2] + avgSNR[3]);

	 bool isSpikey=false;
	 // We now need to see if the spikey ratio is greater than the cuts imposed by Ming-Yuan Lu.
	 if(config==1 && spikeyRatio>=2.4055) isSpikey=true;
	 if(config==2 && spikeyRatio>=2.5267) isSpikey=true;
	 if(config==3 && spikeyRatio>=3.5007) isSpikey=true;
	 if(config==4 && spikeyRatio>=3.9877) isSpikey=true;
	 if(config==5 && spikeyRatio>=3.6174) isSpikey=true;

	 return isSpikey;
 }

 /* Check for cliff event. Cut created by Ming-Yuan Lu and adapted by Jorge Torres
 for the OSU A23 analysis.

 /* Check if a cliff string exists. If so, the event is considered a cliff event and discarded.
 /* A cliff string is defined as either string 1,2,3 where all 4 channels on it shows |median difference| > predefined thresholds
 */
 bool isCliffEvent(vector <TGraph*> grInt){
	 int cliff_threshold_A3_string1=100;
	 int cliff_threshold_A3_string2=45;
	 int cliff_threshold_A3_string3=100;
	 int cliffCount_string1, cliffCount_string2, cliffCount_string3;
	 int const IRS2SamplePerBlock=64;
	 bool isCliff;
	 double firstBlockMedian, lastBlockMedian;
	 double firstBlockSamples[IRS2SamplePerBlock], lastBlockSamples[IRS2SamplePerBlock];
	 cliffCount_string1 = cliffCount_string2 = cliffCount_string3 = 0;
	 int len;
	 for (int ch=0; ch<16; ch++){
		 int len = grInt[ch]->GetN();
		 double * voltValues;
		 voltValues = grInt[ch]->GetY();
		 for (int s=0; s<IRS2SamplePerBlock; s++){
			 firstBlockSamples[s] = voltValues[s];
			 lastBlockSamples[s]  = voltValues[len-IRS2SamplePerBlock+s];
		 }
		 firstBlockMedian = TMath::Median(IRS2SamplePerBlock, firstBlockSamples);
		 lastBlockMedian  = TMath::Median(IRS2SamplePerBlock, lastBlockSamples);
		 double medianDiff = fabs(firstBlockMedian - lastBlockMedian);

		 if (ch%4 == 0){ //string 1
			 if (medianDiff > cliff_threshold_A3_string1){
				 isCliff = true;
			 }
		 } else if (ch%4 == 1) { //string 2
			 if (medianDiff > cliff_threshold_A3_string2){
				 isCliff = true;
			 }
		 } else if (ch%4 == 2){ //string 3
			 if (medianDiff > cliff_threshold_A3_string3){
				 isCliff = true;
			 }
		 }
	 }//end of ch

	 return isCliff;
 }

 //hasOutofBandIssue: Cut implemented by Jorge Torres.
 //Inputs: vector of dim=16 containing waveforms for each channels
 //Outputs: returns a false or true depending on whether the event has more than
 //10% of power below 120 MHz.

 //Other version implemented by Ming-Yuan Lu. It tags the event as bad if the peak power bin
 //is out of band.
 bool hasOutofBandIssue(vector <TGraph*> wform, bool dropDDA4){
	 double interpolation_step = 0.5;
	 int counts = 0;
	 bool isGlitch=false;

	 // for (int i = 0; i < 16; i++){
		//  TGraph *Waveform_Interpolated = FFTtools::getInterpolatedGraph(wform[i],interpolation_step);
		//  //	delete gr;
		//  TGraph *Waveform_Padded = FFTtools::padWaveToLength(Waveform_Interpolated, Waveform_Interpolated->GetN()+6000);
		//  delete Waveform_Interpolated;
		//  TGraph *Waveform_Cropped=FFTtools::cropWave(Waveform_Padded,-300.,300.);
		//  delete Waveform_Padded;
		//  TGraph* spectra = FFTtools::makePowerSpectrumMilliVoltsNanoSeconds(Waveform_Cropped);
		//  double outOfBandPower = 0;
		//  double integral = 0;
		//  double fracOutPower = 0;
		//  int num_bins = spectra->GetN();
		//  int upBin = (int) 50*(spectra->GetX()[num_bins-1]-spectra->GetX()[0])/num_bins;//120 MHz is the freq below which it's out of band
		//  // printf("number of bins is:%i, xmin is:%f, xmax is:%f\n",num_bins,spectra->GetX()[0],spectra->GetX()[num_bins-1]);
		//  for(int samp=0; samp<num_bins; samp++){
	 //     integral+=spectra->GetY()[samp];
	 //   }
		//  for(int j = 0; j < upBin+1; j++){
		// 	 outOfBandPower+= spectra->GetY()[j];
		//  }
		//  delete spectra;
		//  fracOutPower = outOfBandPower/integral;
		//  // cout << fracOutPower << endl;
		//  if(fracOutPower>=0.5){
		// 	 isGlitch=true;
		// 	 break;
		//  }
	 // }
	 for (int i = 0; i < 16; i++){
		 TGraph *Waveform_Interpolated = FFTtools::getInterpolatedGraph(wform[i],interpolation_step);
		 //	delete gr;
		 TGraph *Waveform_Padded = FFTtools::padWaveToLength(Waveform_Interpolated, Waveform_Interpolated->GetN()+6000);
		 delete Waveform_Interpolated;
		 TGraph *Waveform_Cropped=FFTtools::cropWave(Waveform_Padded,-300.,300.);
		 delete Waveform_Padded;
		 TGraph* spectra = FFTtools::makePowerSpectrumMilliVoltsNanoSeconds(Waveform_Cropped);
		 int num_bins = spectra->GetN();
		 int upBin = (int) 100/((spectra->GetX()[num_bins-1]-spectra->GetX()[0])/num_bins);//120 MHz is the freq below which it's out of band
		 // printf("number of bins is:%i, xmin is:%f, xmax is:%f\n",num_bins,spectra->GetX()[0],spectra->GetX()[num_bins-1]);
		 int peakBin = FFTtools::getPeakBin(spectra);
		 // cout << peakBin << endl;
		 delete spectra;
		 delete Waveform_Cropped;

		 if(dropDDA4==true && (i==3 || i==7 || i==11 || i==15)) continue;
		 if(peakBin<upBin){
			 // printf("PeakBin is:%i, UpBin is:%i\n",peakBin,upBin);
			 // cout << i << endl;
			 // if (peakBin>0) counts+=1;
			 counts+=1;
		}
	 }
	 if(counts>3) isGlitch=true;
	 // printf("counts is %i\n",counts);
	 return isGlitch;
 }

/*
	input: antenna_powers (take a vector of the power contained in each channel), and dropBadChans (for if to drop bad channels)
	output: the ratio

	function: return the ratio of the average antenna power on the string with the most power
				over the average antenna power on the string with the next most power
*/
double returnStringPowerRatio(vector<double> antenna_powers, int station, bool dropBadChans){

 	// protections against Brian's stupidity
 	vector<double> average_string_powers;
 	double ratio;
 	if(antenna_powers.size()<16){
 		ratio=0.;
 		return ratio;
 	}

 	if(station==2){
 		average_string_powers.push_back((antenna_powers[0]+antenna_powers[4]+antenna_powers[8]+antenna_powers[12])/4.);
 		average_string_powers.push_back((antenna_powers[1]+antenna_powers[5]+antenna_powers[9]+antenna_powers[13])/4.);
 		average_string_powers.push_back((antenna_powers[2]+antenna_powers[6]+antenna_powers[10]+antenna_powers[14])/4.);
 		average_string_powers.push_back((antenna_powers[3]+antenna_powers[7]+antenna_powers[11])/3.);
 	}

 	if(station==3){
 		if(!dropBadChans){
			average_string_powers.push_back((antenna_powers[0]+antenna_powers[4]+antenna_powers[8]+antenna_powers[12])/4.);
 			average_string_powers.push_back((antenna_powers[1]+antenna_powers[5]+antenna_powers[9]+antenna_powers[13])/4.);
 			average_string_powers.push_back((antenna_powers[2]+antenna_powers[6]+antenna_powers[10]+antenna_powers[14])/4.);
 			average_string_powers.push_back((antenna_powers[3]+antenna_powers[7]+antenna_powers[11]+antenna_powers[15])/4.);
 		}
 		else if(dropBadChans){
 			average_string_powers.push_back((antenna_powers[0]+antenna_powers[4]+antenna_powers[8]+antenna_powers[12])/4.);
 			average_string_powers.push_back((antenna_powers[1]+antenna_powers[5]+antenna_powers[9]+antenna_powers[13])/4.);
 			average_string_powers.push_back((antenna_powers[2]+antenna_powers[6]+antenna_powers[10]+antenna_powers[14])/4.);
 			average_string_powers.push_back(0.); // push back a zero for string 4, just to keep the vector the right length
 		}
 	}

 	//sort smallest to largest
 	if(average_string_powers.size()!=4){
 		cout<<"WARNING! This array is too short"<<endl;
 	}
 	std::sort(average_string_powers.begin(), average_string_powers.end());
	ratio = average_string_powers[3]/average_string_powers[2];
	return ratio;
}
/*
	input: waveforms (vector of waveforms; raw, interpolated, don't matter), and station, runNum
	output: does the event have too much power concentrated in one run or not?

	function: check if the event has too much power concentrated in a single run,
				like we saw in the A2 unblindinging
*/
bool isHighPowerStringEvent(vector<TGraph*> waveforms, int station, int config){

	// protections against Brian's stupidity
	bool this_isHighPowerEvent=false;
	if(waveforms.size()<16 || (station!=2 && station!=3)){
		this_isHighPowerEvent=false;
		return this_isHighPowerEvent;
	}

	// compute power in the antennas
	vector<double> powers_on_antennas;
	for(int chan=0; chan<16; chan++){
		Double_t *yVals = waveforms[chan]->GetY();
		int N = waveforms[chan]->GetN();
		double thisPower =0.;
		for(int samp=0; samp<N; samp++){
			thisPower+=yVals[samp]*yVals[samp];
		}
		powers_on_antennas.push_back(thisPower);
	}

	// figure out if we need to drop bad channels based on the runNum
	bool dropBadChans=false;
	if(station==2){
		dropBadChans=true;
	}
	else if(station==3 && config>2){
		dropBadChans=true;
	}

	// get the ratio
	double ratio = returnStringPowerRatio(powers_on_antennas, station, dropBadChans);

	// if the average power in string with highest power is >5x average power in string with next highest power
	// mark it as a high power event
	if(ratio>5.)
		this_isHighPowerEvent=true;

	// return
	return this_isHighPowerEvent;
}


/*
	input: waveforms (vector of waveforms; raw, interpolated, don't matter), freqLimit (below what we consider out of band),
	dropDDA4 (drop string 4?), absPower (vector to be filled)
	output: vector conataining absolute power per channel

	function: calculates the absolute power concentrated below "freqLimit" MHz.
*/

int outOfBandAbsPower(vector <TGraph*> wform, bool dropDDA4, double freqLimit, vector <double> & absPower){
	double interpolation_step = 0.5;
	for (int i = 0; i < 16; i++){
		if(dropDDA4==true && (i==3 || i==7 || i==11 || i==15)){
			absPower.push_back(-1);
			continue;
		}
		TGraph *Waveform_Interpolated = FFTtools::getInterpolatedGraph(wform[i],interpolation_step);
		TGraph *Waveform_Padded = FFTtools::padWaveToLength(Waveform_Interpolated, Waveform_Interpolated->GetN()+6000);
		delete Waveform_Interpolated;
		TGraph *Waveform_Cropped=FFTtools::cropWave(Waveform_Padded,-500.,700.);
		delete Waveform_Padded;
		TGraph* spectra = FFTtools::makePowerSpectrumMilliVoltsNanoSeconds(Waveform_Cropped);
		delete Waveform_Cropped;
		double outOfBandPower = 0.;
		int num_bins = spectra->GetN();
		int upBin = (int) freqLimit/((spectra->GetX()[num_bins-1]-spectra->GetX()[0])/num_bins);//freqLimit is the freq below which it's out of band
		for(int j = 0; j < upBin+1; j++){
			outOfBandPower+= spectra->GetY()[j];
		}
		delete spectra;
		absPower.push_back(outOfBandPower);
	}
	return 0;
}
