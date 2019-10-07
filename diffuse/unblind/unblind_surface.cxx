////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	unblind_surface.cxx
////	unblind the surface cut
////
////	Sep 2019
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
#include "tools_Cuts.h"
#include "tools_Stats.h"
#include "tools_CommandLine.h"
#include "tools_outputObjects.h"
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"

// for ray trace correlator
#include "Settings.h"
#include "Event.h"
#include "Detector.h"
#include "Report.h"
#include "RayTraceCorrelator.h"

AraAntPol::AraAntPol_t Vpol = AraAntPol::kVertical;
AraAntPol::AraAntPol_t Hpol = AraAntPol::kHorizontal;

using namespace std;

void configure(TH1D *gr);

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


	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <1-station> <2-config> <100 or 10> <isOrg>"<<endl;;
		return -1;
	}
	int station = atoi(argv[1]);
	int config = atoi(argv[2]);
	int full_or_partial = atoi(argv[3]);
	int isOrg = atoi(argv[4]);

	if(station!=2 && station!=3){
		printf("No good! You asked for station %d, but this code only works for stations 2 and 3 \n",station);
		return -1;
	}

	double side_cut_corr[2] = {0.005, 0.004}; // in corr, side region is <=0.005 for V and <=0.004 in H
	double side_cut_snr[2] = {4.,4.}; // in snr, side region is <=4 in V and H

	// a place to print things to file
	char title_txt[200];
	if(full_or_partial==100 && isOrg==0){
		sprintf(title_txt,"%s/unblind/surface/reduced_surface/A%d_c%d_%dSample_EventList.txt", plotPath,station,config,full_or_partial);
	}
	else{
		sprintf(title_txt,"%s/unblind/surface/A%d_c%d_%dSample_EventList.txt", plotPath,station,config,full_or_partial);
	}
	
	vector<int> BadRunList=BuildBadRunList(station);
	vector<int> BadSurfaceRunList=BuildSurfaceRunList(station);

	// for(int i=0; i<BadSurfaceRunList.size(); i++){
	// 	printf("Bad run %d \n", BadSurfaceRunList[i]);
	// }
	// return 0;

	int numTotal=0;

	// need to be able to make the final 2D distribution
	double max=0.05;
	TH2D *h2SNRvsCorr[2]; // SNR on Y axis, Corr on X axis, like in the TB
	h2SNRvsCorr[0]=new TH2D("","V Data",100,0,max,300,0,30);
	h2SNRvsCorr[1]=new TH2D("","H Data",100,0,max,300,0,30);


	// we also want to work up a plot for the surface rate
	TTimeStamp start;
	TTimeStamp stop;
	int numBins=365*4;

	start.Set(2013, 01, 00, 00, 00,0,0,true,0);
	stop.Set(2016, 12, 31, 00, 00,0,0,true,0);

	if(config==1 || config==2){
		start.Set(2013, 01, 00, 00, 00,0,0,true,0);
		stop.Set(2014, 12, 31, 00, 00,0,0,true,0);
		numBins=365*2;
	}
	if(config==4){
		start.Set(2014, 01, 00, 00, 00,0,0,true,0);
		stop.Set(2015, 12, 31, 00, 00,0,0,true,0);
		numBins=365*2;
	}
	if(config==5){
		start.Set(2015, 01, 00, 00, 00,0,0,true,0);
		stop.Set(2016, 12, 31, 00, 00,0,0,true,0);
		numBins=365*2;
	}

	TH1D *phi_dist = new TH1D ("","",360, -180, 180);

	TH2D *zoom_burst = new TH2D("","",360,-180,180,180,-90,90);

	int start_bin = start.GetSec();
	int stop_bin = stop.GetSec();

	TH1D *events_vs_time[2];
	for(int i=0; i<2; i++){
		if(i==0){
			ss.str("");
			ss<<"VPol";
		}
		else if(i==1){
			ss.str("");
			ss<<"HPol";
		}
		events_vs_time[i] = new TH1D(ss.str().c_str(),ss.str().c_str(),numBins, start_bin, stop_bin);
		events_vs_time[i]->GetXaxis()->SetTimeDisplay(1); //turn on a time axis
		events_vs_time[i]->GetXaxis()->SetTimeFormat("%y/%m");
		events_vs_time[i]->GetXaxis()->SetTimeOffset(0.,"GMT");
	}

	TChain dataVTree("VTree");
	TChain dataHTree("HTree");
	TChain dataAllTree("AllTree");
	TChain dataRecoTree("OutputTreeReco");
	TChain dataFilterTree("OutputTree");
	char the_data[500];

	if(full_or_partial==100 && isOrg==1){
		// use the 100pct sample
		sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_*.root",station,config);

	}
	if(full_or_partial==10 && isOrg==1){
		// use the *original* 10pct sample from optimization
		// sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_*.root",station,config);
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

	if(full_or_partial==10 && isOrg==1){ // to make plotting just a little easier

		vector<int> runs_to_loop_over;
		if(config==1){
			runs_to_loop_over.push_back(2636);
			runs_to_loop_over.push_back(2662);
			runs_to_loop_over.push_back(2678);
			runs_to_loop_over.push_back(2990);
			runs_to_loop_over.push_back(3206);
		}
		if(config==2){
			runs_to_loop_over.push_back(1830);
			runs_to_loop_over.push_back(1835);
			runs_to_loop_over.push_back(1860);
			runs_to_loop_over.push_back(2090);
			runs_to_loop_over.push_back(2091);
			runs_to_loop_over.push_back(2155);
		}
		if(config==4){
			runs_to_loop_over.push_back(4398);
			runs_to_loop_over.push_back(4497);
			runs_to_loop_over.push_back(4777);
			runs_to_loop_over.push_back(4817);
			runs_to_loop_over.push_back(4837);
			runs_to_loop_over.push_back(4842);
			runs_to_loop_over.push_back(4979);
			runs_to_loop_over.push_back(5004);
			runs_to_loop_over.push_back(5007);
			runs_to_loop_over.push_back(5032);
			runs_to_loop_over.push_back(5174);
			runs_to_loop_over.push_back(5369);
			runs_to_loop_over.push_back(5515);
			runs_to_loop_over.push_back(5516);
			runs_to_loop_over.push_back(5517);
			runs_to_loop_over.push_back(5586);
			runs_to_loop_over.push_back(5615);
			runs_to_loop_over.push_back(5649);
			runs_to_loop_over.push_back(5660);
			runs_to_loop_over.push_back(5664);
			runs_to_loop_over.push_back(5666);
			runs_to_loop_over.push_back(5670);
			runs_to_loop_over.push_back(5680);
			runs_to_loop_over.push_back(5699);
			runs_to_loop_over.push_back(5701);
			runs_to_loop_over.push_back(5702);
			runs_to_loop_over.push_back(5704);
			runs_to_loop_over.push_back(5775);
			runs_to_loop_over.push_back(5849);
			runs_to_loop_over.push_back(6080);
			runs_to_loop_over.push_back(6445);
		}
		if(config==5){
			runs_to_loop_over.push_back(6532);
			runs_to_loop_over.push_back(6536);
			runs_to_loop_over.push_back(6542);
			runs_to_loop_over.push_back(6554);
			runs_to_loop_over.push_back(6577);
			runs_to_loop_over.push_back(6635);
			runs_to_loop_over.push_back(6655);
			runs_to_loop_over.push_back(6669);
			runs_to_loop_over.push_back(6674);
			runs_to_loop_over.push_back(6679);
			runs_to_loop_over.push_back(6694);
			runs_to_loop_over.push_back(6705);
			runs_to_loop_over.push_back(6733);
			runs_to_loop_over.push_back(6818);
			runs_to_loop_over.push_back(7836);
			runs_to_loop_over.push_back(8074);
		}

		for(int i=0; i<runs_to_loop_over.size(); i++){
			sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_%d.root",station,config, runs_to_loop_over[i]);
			dataVTree.Add(the_data);
			dataHTree.Add(the_data);
			dataAllTree.Add(the_data);
			dataFilterTree.Add(the_data);
		}
	}

	if(full_or_partial==100 && isOrg==0){ // to make plotting just a little easier

		vector<int> runs_to_loop_over;
		if(config==1){
			runs_to_loop_over.push_back(2466);
			runs_to_loop_over.push_back(2589);
			runs_to_loop_over.push_back(2626);
			runs_to_loop_over.push_back(2636);
			runs_to_loop_over.push_back(2662);
			runs_to_loop_over.push_back(2678);
			runs_to_loop_over.push_back(2779);
			runs_to_loop_over.push_back(2784);
			runs_to_loop_over.push_back(2849);
			runs_to_loop_over.push_back(2864);
			runs_to_loop_over.push_back(2869);
			runs_to_loop_over.push_back(2871);
			runs_to_loop_over.push_back(2934);
			runs_to_loop_over.push_back(2937);
			runs_to_loop_over.push_back(2990);
			runs_to_loop_over.push_back(3043);
			runs_to_loop_over.push_back(3202);
			runs_to_loop_over.push_back(3206);
			runs_to_loop_over.push_back(3232);
			runs_to_loop_over.push_back(3233);
			runs_to_loop_over.push_back(3325);
			runs_to_loop_over.push_back(3364);
			runs_to_loop_over.push_back(3392);
			runs_to_loop_over.push_back(3412);
		}
		if(config==2){
			runs_to_loop_over.push_back(1455);
			runs_to_loop_over.push_back(1514);
			runs_to_loop_over.push_back(1518);
			runs_to_loop_over.push_back(1571);
			runs_to_loop_over.push_back(1629);
			runs_to_loop_over.push_back(1647);
			runs_to_loop_over.push_back(1689);
			runs_to_loop_over.push_back(1647);
			runs_to_loop_over.push_back(1689);
			runs_to_loop_over.push_back(1776);
			runs_to_loop_over.push_back(1778);
			runs_to_loop_over.push_back(1779);
			runs_to_loop_over.push_back(1830);
			runs_to_loop_over.push_back(1831);
			runs_to_loop_over.push_back(1834);
			runs_to_loop_over.push_back(1835);
			runs_to_loop_over.push_back(1849);
			runs_to_loop_over.push_back(1850);
			runs_to_loop_over.push_back(1855);
			runs_to_loop_over.push_back(1858);
			runs_to_loop_over.push_back(1859);
			runs_to_loop_over.push_back(1860);
			runs_to_loop_over.push_back(1950);
			runs_to_loop_over.push_back(1952);
			runs_to_loop_over.push_back(1953);
			runs_to_loop_over.push_back(2035);
			runs_to_loop_over.push_back(2049);
			runs_to_loop_over.push_back(2087);
			runs_to_loop_over.push_back(2090);
			runs_to_loop_over.push_back(2091);
			runs_to_loop_over.push_back(2155);
			runs_to_loop_over.push_back(2156);
			runs_to_loop_over.push_back(2165);
			runs_to_loop_over.push_back(2173);
			runs_to_loop_over.push_back(2174);
		}
		if(config==3){
			runs_to_loop_over.push_back(3529);
			runs_to_loop_over.push_back(3534);
			runs_to_loop_over.push_back(3543);
			runs_to_loop_over.push_back(3663);
			runs_to_loop_over.push_back(3726);
			runs_to_loop_over.push_back(3766);
			runs_to_loop_over.push_back(3772);
			runs_to_loop_over.push_back(3905);
			runs_to_loop_over.push_back(3906);
			runs_to_loop_over.push_back(3907);
			runs_to_loop_over.push_back(3909);
			runs_to_loop_over.push_back(3917);
			runs_to_loop_over.push_back(4006);
			runs_to_loop_over.push_back(4008);
		}
		if(config==4){
			runs_to_loop_over.push_back(4106);
			runs_to_loop_over.push_back(4114);
			runs_to_loop_over.push_back(4140);
			runs_to_loop_over.push_back(4317);
			runs_to_loop_over.push_back(4336);
			runs_to_loop_over.push_back(4386);
			runs_to_loop_over.push_back(4398);
			runs_to_loop_over.push_back(4399);
			runs_to_loop_over.push_back(4406);
			runs_to_loop_over.push_back(4407);
			runs_to_loop_over.push_back(4408);
			runs_to_loop_over.push_back(4458);
			runs_to_loop_over.push_back(4474);
			runs_to_loop_over.push_back(4486);
			runs_to_loop_over.push_back(4497);
			runs_to_loop_over.push_back(4536);
			runs_to_loop_over.push_back(4572);
			runs_to_loop_over.push_back(4622);
			runs_to_loop_over.push_back(4624);
			runs_to_loop_over.push_back(4651);
			runs_to_loop_over.push_back(4772);
			runs_to_loop_over.push_back(4777);
			runs_to_loop_over.push_back(4817);
			runs_to_loop_over.push_back(4837);
			runs_to_loop_over.push_back(4842);
			runs_to_loop_over.push_back(4855);
			runs_to_loop_over.push_back(4943);
			runs_to_loop_over.push_back(4947);
			runs_to_loop_over.push_back(4979);
			runs_to_loop_over.push_back(5002);
			runs_to_loop_over.push_back(5004);
			runs_to_loop_over.push_back(5007);
			runs_to_loop_over.push_back(5009);
			runs_to_loop_over.push_back(5032);
			runs_to_loop_over.push_back(5042);
			runs_to_loop_over.push_back(5043);
			runs_to_loop_over.push_back(5047);
			runs_to_loop_over.push_back(5049);
			runs_to_loop_over.push_back(5052);
			runs_to_loop_over.push_back(5062);
			runs_to_loop_over.push_back(5094);
			runs_to_loop_over.push_back(5162);
			runs_to_loop_over.push_back(5168);
			runs_to_loop_over.push_back(5169);
			runs_to_loop_over.push_back(5170);
			runs_to_loop_over.push_back(5174);
			runs_to_loop_over.push_back(5185);
			runs_to_loop_over.push_back(5310);
			runs_to_loop_over.push_back(5324);
			runs_to_loop_over.push_back(5356);
			runs_to_loop_over.push_back(5357);
			runs_to_loop_over.push_back(5369);
			runs_to_loop_over.push_back(5405);
			runs_to_loop_over.push_back(5406);
			runs_to_loop_over.push_back(5446);
			runs_to_loop_over.push_back(5505);
			runs_to_loop_over.push_back(5512);
			runs_to_loop_over.push_back(5514);
			runs_to_loop_over.push_back(5515);
			runs_to_loop_over.push_back(5516);
			runs_to_loop_over.push_back(5517);
			runs_to_loop_over.push_back(5539);
			runs_to_loop_over.push_back(5586);
			runs_to_loop_over.push_back(5600);
			runs_to_loop_over.push_back(5610);
			runs_to_loop_over.push_back(5614);
			runs_to_loop_over.push_back(5615);
			runs_to_loop_over.push_back(5616);
			runs_to_loop_over.push_back(5617);
			runs_to_loop_over.push_back(5619);
			runs_to_loop_over.push_back(5620);
			runs_to_loop_over.push_back(5621);
			runs_to_loop_over.push_back(5625);
			runs_to_loop_over.push_back(5645);
			runs_to_loop_over.push_back(5649);
			runs_to_loop_over.push_back(5650);
			runs_to_loop_over.push_back(5652);
			runs_to_loop_over.push_back(5654);
			runs_to_loop_over.push_back(5660);
			runs_to_loop_over.push_back(5661);
			runs_to_loop_over.push_back(5662);
			runs_to_loop_over.push_back(5664);
			runs_to_loop_over.push_back(5666);
			runs_to_loop_over.push_back(5667);
			runs_to_loop_over.push_back(5669);
			runs_to_loop_over.push_back(5670);
			runs_to_loop_over.push_back(5671);
			runs_to_loop_over.push_back(5675);
			runs_to_loop_over.push_back(5676);
			runs_to_loop_over.push_back(5680);
			runs_to_loop_over.push_back(5681);
			runs_to_loop_over.push_back(5683);
			runs_to_loop_over.push_back(5684);
			runs_to_loop_over.push_back(5699);
			runs_to_loop_over.push_back(5700);
			runs_to_loop_over.push_back(5701);
			runs_to_loop_over.push_back(5702);
			runs_to_loop_over.push_back(5704);
			runs_to_loop_over.push_back(5706);
			runs_to_loop_over.push_back(5760);
			runs_to_loop_over.push_back(5764);
			runs_to_loop_over.push_back(5774);
			runs_to_loop_over.push_back(5775);
			runs_to_loop_over.push_back(5780);
			runs_to_loop_over.push_back(5799);
			runs_to_loop_over.push_back(5810);
			runs_to_loop_over.push_back(5811);
			runs_to_loop_over.push_back(5820);
			runs_to_loop_over.push_back(5847);
			runs_to_loop_over.push_back(5849);
			runs_to_loop_over.push_back(5854);
			runs_to_loop_over.push_back(5861);
			runs_to_loop_over.push_back(5881);
			runs_to_loop_over.push_back(5894);
			runs_to_loop_over.push_back(5895);
			runs_to_loop_over.push_back(5896);
			runs_to_loop_over.push_back(5897);
			runs_to_loop_over.push_back(6080); //not in 100pct?
			runs_to_loop_over.push_back(6445);
		}
		if(config==5){
			runs_to_loop_over.push_back(6506);
			runs_to_loop_over.push_back(6511);
			runs_to_loop_over.push_back(6532);
			runs_to_loop_over.push_back(6536);
			runs_to_loop_over.push_back(6542);
			runs_to_loop_over.push_back(6544);
			runs_to_loop_over.push_back(6554);
			runs_to_loop_over.push_back(6561);
			runs_to_loop_over.push_back(6577);
			runs_to_loop_over.push_back(6588);
			runs_to_loop_over.push_back(6593);
			runs_to_loop_over.push_back(6603);
			runs_to_loop_over.push_back(6610);
			runs_to_loop_over.push_back(6628);
			runs_to_loop_over.push_back(6635);
			runs_to_loop_over.push_back(6655);
			runs_to_loop_over.push_back(6657);
			runs_to_loop_over.push_back(6669);
			runs_to_loop_over.push_back(6673);
			runs_to_loop_over.push_back(6674);
			runs_to_loop_over.push_back(6675);
			runs_to_loop_over.push_back(6679);
			runs_to_loop_over.push_back(6687);
			runs_to_loop_over.push_back(6694);
			runs_to_loop_over.push_back(6705);
			runs_to_loop_over.push_back(6733);
			runs_to_loop_over.push_back(6759);
			runs_to_loop_over.push_back(6818);
			runs_to_loop_over.push_back(6820);
			runs_to_loop_over.push_back(6851);
			runs_to_loop_over.push_back(6857);
			runs_to_loop_over.push_back(6860);
			runs_to_loop_over.push_back(6861);
			runs_to_loop_over.push_back(6876);
			runs_to_loop_over.push_back(6943);
			runs_to_loop_over.push_back(7170);
			runs_to_loop_over.push_back(7264);
			runs_to_loop_over.push_back(7418);
			runs_to_loop_over.push_back(7434);
			runs_to_loop_over.push_back(7494);
			runs_to_loop_over.push_back(7501);
			runs_to_loop_over.push_back(7584);
			runs_to_loop_over.push_back(7814);
			runs_to_loop_over.push_back(7824);
			runs_to_loop_over.push_back(7836);
			runs_to_loop_over.push_back(7923);
			runs_to_loop_over.push_back(8052);
			runs_to_loop_over.push_back(8064);
			runs_to_loop_over.push_back(8074);
			runs_to_loop_over.push_back(8085);
		}

		for(int i=0; i<runs_to_loop_over.size(); i++){
			// sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/10pct_redo/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_%d.root",station,config, runs_to_loop_over[i]);
			sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_%d.root",station,config, runs_to_loop_over[i]);
			// printf("/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_%d.root \n",station,config, runs_to_loop_over[i]);
			dataVTree.Add(the_data);
			dataHTree.Add(the_data);
			dataAllTree.Add(the_data);
			dataFilterTree.Add(the_data);
		}
	}


	// sprintf(the_data,"/fs/project/PAS0654/ARA_DATA/A23/100pct_try2/ValsForCuts/A%d/c%d/cutvals_drop_FiltSurface_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_6171.root",station,config);
	// dataVTree.Add(the_data);
	// dataHTree.Add(the_data);
	// dataAllTree.Add(the_data);
	// dataFilterTree.Add(the_data);

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

		dataAllTree.SetBranchAddress("surf_V_new",&isSurf[0]);
		dataAllTree.SetBranchAddress("surf_H_new",&isSurf[1]);

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

		dataAllTree.SetBranchAddress("bad",&isBadEvent);
		dataAllTree.SetBranchAddress("weight",&weight);
		dataAllTree.SetBranchAddress("unixTime",&unixTime);
		dataAllTree.SetBranchAddress("isFirstFiveEvent",&isFirstFiveEvent);
		dataAllTree.SetBranchAddress("hasBadSpareChanIssue",&hasBadSpareChanIssue);
		dataAllTree.SetBranchAddress("hasBadSpareChanIssue2",&hasBadSpareChanIssue2);
		dataAllTree.SetBranchAddress("runNum",&runNum);
		dataAllTree.SetBranchAddress("eventNumber",&eventNumber);

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
		
		// check both
		bool isThisABadRun = isBadRun(station,runNum,BadRunList);
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

			// continue;
			if( isSoft || isBadEvent || hasBadSpareChanIssue || hasBadSpareChanIssue2 || isFirstFiveEvent || isShort || isCal || isThisABadRun){
				continue;
			}
			if(isBadLivetime(station,unixTime)){
				continue;
			}
			for(int pol=0; pol<2; pol++){
				if(!WFRMS[pol] 
					&& !isNewBox 
					&& (isSurf[0] || isSurf[1] || isSurfEvent_top[pol])
				){
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

						double which_corr_to_use = corr_val[pol];
						bool this_pass_R_cut = passesRCut(station, config, pol, snr_val[pol], which_corr_to_use);

						if(this_pass_R_cut){

							events_vs_time[pol]->Fill(unixTime);
							h2SNRvsCorr[pol]->Fill(which_corr_to_use,snr_val[pol],weight);

							if(pol==0){
								phi_dist->Fill(phi_300[pol]);
							}

							// if(runNum==2090 || runNum==2091){
							// 	zoom_burst->Fill(phi_300[pol], theta_300[pol]);
							// }

							// FILE *fout = fopen(title_txt, "a");
							// fprintf(fout,"Run %4d, Event %6d, Pol %1d, unixTime %d, theta %3d, phi %3d, coor %.4f, refilt %d, surfv %d, surfh %d, surftopppol %d \n",runNum, eventNumber, pol, unixTime, theta_300[pol], phi_300[pol], corr_val[pol], Refilt[pol], isSurf[0], isSurf[1], isSurfEvent_top[pol]);
							// fclose(fout);//close sigmavsfreq.txt file

							// printf(BLUE"\n Run %d, Event %d, pol %d, unixTime %d, reco theta %d, phi %d, Refilt status %d, corr is %.4f \n"RESET, runNum, eventNumber, pol, unixTime, theta_300[pol], phi_300[pol], Refilt[pol], corr_val[pol]);
							// PlotThisEvent(station, config, runNum, eventNumber, pol);
						}

					}// not failing CW power cut?
				}// passes rest of analysis (not WFRMS, box, surface)
			}// loop over polarizations
		}// loop over events
	}
	std::cout<<endl;

	// TCanvas *cZoomBurst = new TCanvas("","",1.1*850,850);
	// 	zoom_burst->Draw("colz");
	// 	zoom_burst->GetXaxis()->SetTitle("Phi [deg]");
	// 	zoom_burst->GetYaxis()->SetTitle("Theta [deg]");
	// 	zoom_burst->GetZaxis()->SetTitle("Number of Events");
	// 	zoom_burst->GetXaxis()->SetRangeUser(-136,-124);
	// 	zoom_burst->GetYaxis()->SetRangeUser(36,48);
	// 	gPad->SetRightMargin(0.15);
	// char thistitle_zoomburst[300];
	// sprintf(thistitle_zoomburst, "%s/unblind/surface/%d.%d.%d_A%d_c%d_%dEvents_SurfaceRate_%dSample_ZoomBurst.png",plotPath,year_now,month_now,day_now,station,config,int(numTotal),full_or_partial);
	// cZoomBurst->SaveAs(thistitle_zoomburst);
	// delete cZoomBurst;

	
	if(full_or_partial==100 && isOrg==0){

		TCanvas *cPhiDist = new TCanvas("","",1.1*850,850);
		phi_dist->Draw("");
			phi_dist->GetXaxis()->SetTitle("Phi [deg]");
			phi_dist->GetYaxis()->SetTitle("Number of Events");
			gPad->SetLogy();
		char thistitle[300];
		sprintf(thistitle, "%s/unblind/surface/reduced_surface/%d.%d.%d_A%d_c%d_%dEvents_PhiDist_%dSample.png",plotPath,year_now,month_now,day_now,station,config,int(numTotal),full_or_partial);
		cPhiDist->SaveAs(thistitle);
		delete cPhiDist;
	}

	TCanvas *cSurfaceRate = new TCanvas("","",2.1*850,2.1*850);
	cSurfaceRate->Divide(1,2);
	for(int i=0; i<2; i++){
		cSurfaceRate->cd(i+1);
		// num_samps[i]->GetYaxis()->SetRangeUser(0.1,2e7);
		// num_samps[i]->GetXaxis()->SetRangeUser(0,750);

		// gPad->SetTopMargin(0.1);
		// gPad->SetRightMargin(0.03);
		// gPad->SetLeftMargin(0.05);
		// gPad->SetBottomMargin(0.11);

		configure(events_vs_time[i]);

		events_vs_time[i]->Draw("");
		events_vs_time[i]->GetXaxis()->SetTitle("UnixTime [YY/MM]");
		events_vs_time[i]->GetYaxis()->SetTitle("Number of Events");
		gPad->SetLogy();
	}
	char thistitle[300];
	sprintf(thistitle, "%s/unblind/surface/%d.%d.%d_A%d_c%d_%dEvents_SurfaceRate_%dSample.png",plotPath,year_now,month_now,day_now,station,config,int(numTotal),full_or_partial);
	if(full_or_partial==100 && isOrg==0){
		sprintf(thistitle, "%s/unblind/surface/reduced_surface/%d.%d.%d_A%d_c%d_%dEvents_SurfaceRate_%dSample.png",plotPath,year_now,month_now,day_now,station,config,int(numTotal),full_or_partial);
	}
	cSurfaceRate->SaveAs(thistitle);
	delete cSurfaceRate;


	vector<TGraph*> cut_lines;
	for(int pol=0; pol<2; pol++){
		vector <double> x_vals_for_line;
		vector <double> y_vals_for_line;
		double rcut_slope;
		double rcut_intercept;
		getRCutValues(station, config, pol, rcut_slope, rcut_intercept);
		for(double x=0; x<0.020; x+=0.00001){
			double y_val = (rcut_slope * x ) + rcut_intercept;
			x_vals_for_line.push_back(x);
			y_vals_for_line.push_back(y_val);
		}
		cut_lines.push_back(new TGraph(x_vals_for_line.size(), &x_vals_for_line[0], &y_vals_for_line[0]));
	}

	// now to actually draw the safety cut lines
	double vert_cuts[2] = {side_cut_corr[0], side_cut_corr[1]};
	vector<TGraph*> vert_lines;
	for(int pol=0; pol<2; pol++){
		vector <double> x_vals_for_line;
		vector <double> y_vals_for_line;
		for(double y=0; y<30; y+=1.){
			y_vals_for_line.push_back(y);
			x_vals_for_line.push_back(vert_cuts[pol]);
		}
		vert_lines.push_back(new TGraph(x_vals_for_line.size(), &x_vals_for_line[0], &y_vals_for_line[0]));
	}

	double horz_cuts[2] = {side_cut_snr[0], side_cut_snr[1]};
	vector<TGraph*> horz_lines;
	for(int pol=0; pol<2; pol++){
		vector <double> x_vals_for_line;
		vector <double> y_vals_for_line;
		for(double x=0; x<0.05; x+=0.01){
			y_vals_for_line.push_back(horz_cuts[pol]);
			x_vals_for_line.push_back(x);
		}
		horz_lines.push_back(new TGraph(x_vals_for_line.size(), &x_vals_for_line[0], &y_vals_for_line[0]));
	}

	// okay, now save out the 2D histogram
	TCanvas *cSNRvsCorr = new TCanvas("","",2.1*850, 850);
	cSNRvsCorr->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		cSNRvsCorr->cd(pol+1);
		h2SNRvsCorr[pol]->Draw("colz");
		h2SNRvsCorr[pol]->GetYaxis()->SetTitle("3rd Highest VPeak/RMS");
		h2SNRvsCorr[pol]->GetXaxis()->SetTitle("Peak Corr");
		gPad->SetLogz();
		cut_lines[pol]->Draw("same");
		cut_lines[pol]->SetLineColor(kRed);
		vert_lines[pol]->Draw("same");
		vert_lines[pol]->SetLineColor(kBlue);
		horz_lines[pol]->Draw("same");
		horz_lines[pol]->SetLineColor(kBlue);
	}
	char title[300];
	if(full_or_partial==100 && isOrg==1)
	  sprintf(title, "%s/unblind/surface/%d.%d.%d_A%d_c%d_%dEvents_UnblindSurface_SNRvsCorr_100Sample.png",plotPath,year_now,month_now,day_now,station,config,numTotal);
	else if(full_or_partial==100 && isOrg==0)
	  sprintf(title, "%s/unblind/surface/reduced_surface/%d.%d.%d_A%d_c%d_%dEvents_UnblindSurface_SNRvsCorr_100Sample.png",plotPath,year_now,month_now,day_now,station,config,numTotal);
	else if(full_or_partial==10 && isOrg==0)
	  sprintf(title, "%s/unblind/surface/%d.%d.%d_A%d_c%d_%dEvents_UnblindSurface_SNRvsCorr_10Sample_New.png",plotPath,year_now,month_now,day_now,station,config,numTotal);
	else if(full_or_partial==10 && isOrg==1)
	  sprintf(title, "%s/unblind/surface/%d.%d.%d_A%d_c%d_%dEvents_UnblindSurface_SNRvsCorr_10Sample_Org.png",plotPath,year_now,month_now,day_now,station,config,numTotal);
	cSNRvsCorr->SaveAs(title);
	delete cSNRvsCorr;
}


void configure(TH1D *gr){
	gr->GetXaxis()->SetLabelSize(0.05);
	gr->GetXaxis()->SetTitleSize(0.05);
	
	gr->GetYaxis()->SetLabelSize(0.05);
	gr->GetYaxis()->SetTitleSize(0.05);
	// gr->GetYaxis()->SetTitleOffset(1.1);
	
	gStyle->SetTitleFontSize(0.05);
	gr->GetXaxis()->SetNdivisions(12,0,0,false);
	gr->SetLineWidth(2);
}