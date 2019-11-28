////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	proposal_csw.cxx
////	make csw's and do something (idk what yet) with them
////	this particular case is just for the vpol singlets which do *not* point at SP
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
#include "TLegend.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"

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
#include "RaySolver.h"
#include "Position.h"
#include "Vector.h"
#include "Settings.h"

using namespace std;

AraAntPol::AraAntPol_t Vpol = AraAntPol::kVertical;
AraAntPol::AraAntPol_t Hpol = AraAntPol::kHorizontal;

// this function is simple enough, just declare it up front
// otherwise it's at the end and clogs things up
void configure(TH1D *gr){
	gr->GetXaxis()->SetLabelSize(0.05);
	gr->GetXaxis()->SetTitleSize(0.05);
	gr->GetYaxis()->SetLabelSize(0.05);
	gr->GetYaxis()->SetTitleSize(0.05);
	// gr->GetYaxis()->SetTitleOffset(1.1);
	gStyle->SetTitleFontSize(0.05);
	gr->SetLineWidth(2);
}

double getTravelTime(Settings *settings, RaySolver *raysolver, double stationCenter[3], int thetaIn, int phiIn, Position target){
	// double phiWave = (double(phiIn)+0.5)*TMath::DegToRad();
	// double thetaWave = (double(thetaIn)+0.5)*TMath::DegToRad();
	double phiWave = (double(phiIn))*TMath::DegToRad();
	double thetaWave = (double(thetaIn))*TMath::DegToRad();
	double R = 300;

	Double_t xs = R*TMath::Cos(thetaWave)*TMath::Cos(phiWave);
	Double_t ys = R*TMath::Cos(thetaWave)*TMath::Sin(phiWave);
	Double_t zs = R*TMath::Sin(thetaWave);

	xs = xs + stationCenter[0];
	ys = ys + stationCenter[1];
	zs = zs + stationCenter[2];

	Position source;
	source.SetXYZ(xs, ys, zs);

	vector < vector <double> > outputs; //place for the answer
	raysolver->Solve_Ray_org(source, target, outputs, settings);
	double TOF=-1000;
	if(outputs.size()>0){
		TOF=outputs[3][0];
		// cout<<"Time of flight is "<<TOF<<endl;
	}
	return TOF;
}

double getAntennaArrivalTheta(int station, Settings *settings, Position target, RaySolver *raysolver, double stationCenter[3], int thetaIn, int phiIn){
	double phiWave = (double(phiIn))*TMath::DegToRad();
	double thetaWave = (double(thetaIn))*TMath::DegToRad();
	double R = 300;
	
	Double_t xs = R*TMath::Cos(thetaWave)*TMath::Cos(phiWave);
	Double_t ys = R*TMath::Cos(thetaWave)*TMath::Sin(phiWave);
	Double_t zs = R*TMath::Sin(thetaWave);

	xs = xs + stationCenter[0];
	ys = ys + stationCenter[1];
	zs = zs + stationCenter[2];

	Position source;
	source.SetXYZ(xs, ys, zs);

	vector < vector <double> > outputs; //place for the answer
	raysolver->Solve_Ray_org(source, target, outputs, settings);
	double thetaOut=-1000;
	if(outputs.size()>0){
		thetaOut=outputs[2][0]*TMath::RadToDeg() -90.;
	}
	thetaOut = abs(thetaOut-90.);
	return thetaOut;
}

void makeCSW(int station, int runNum, int eventNumber, int config, int thetaIn, int phiIn);

int main(int argc, char **argv)
{
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;


	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <1-station> <2-config> <3-save_vals_file_1> <4-save_vals_file_2> .... <n-save_vals_file_n>"<<endl;;
		return -1;
	}
	int station = atoi(argv[1]);
	int config = atoi(argv[2]);

	if(station!=2 && station!=3){
		printf("No good! You asked for station %d, but this code only works for stations 2 and 3 \n",station);
		return -1;
	}

	/*
		Plot declaration
	*/
	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

	stringstream ss;

	vector<int> BadRunList=BuildBadRunList(station);
	vector<int> BadSurfaceRunList=BuildSurfaceRunList(station);

	for(int file_num=3; file_num<argc; file_num++){
		
		TFile *inputFile = TFile::Open(argv[file_num]);
		if(!inputFile){
			cout<<"Can't open vals for cuts file!"<<endl;			
			return -1;
		}
		TTree *inputVTree = (TTree*) inputFile->Get("VTree");
		TTree *inputHTree = (TTree*) inputFile->Get("HTree");
		TTree *inputAllTree = (TTree*) inputFile->Get("AllTree");
		TTree *inputRecoTree = (TTree*) inputFile->Get("OutputTreeReco");
		TTree *inputFilterTree = (TTree*) inputFile->Get("OutputTree");
		TTree *inputSurfaceTree = (TTree*) inputFile->Get("SurfaceTree");

		double corr_val[2];
		double snr_val[2];
		int WFRMS[2];
		int theta_300[2];
		int phi_300[2];
		int theta_41[2];
		int phi_41[2];

		int Refilt[2];
		inputVTree->SetBranchAddress("Refilt_V",&Refilt[0]);
		inputHTree->SetBranchAddress("Refilt_H",&Refilt[1]);

		inputVTree->SetBranchAddress("corr_val_V_new",&corr_val[0]);
		inputVTree->SetBranchAddress("snr_val_V_new",&snr_val[0]);
		inputVTree->SetBranchAddress("wfrms_val_V_new",&WFRMS[0]);
		inputVTree->SetBranchAddress("theta_300_V_new",&theta_300[0]);
		inputVTree->SetBranchAddress("theta_41_V_new",&theta_41[0]);
		inputVTree->SetBranchAddress("phi_300_V_new",&phi_300[0]);
		inputVTree->SetBranchAddress("phi_41_V_new",&phi_41[0]);

		inputHTree->SetBranchAddress("corr_val_H_new",&corr_val[1]);
		inputHTree->SetBranchAddress("snr_val_H_new",&snr_val[1]);
		inputHTree->SetBranchAddress("wfrms_val_H_new",&WFRMS[1]);
		inputHTree->SetBranchAddress("theta_300_H_new",&theta_300[1]);
		inputHTree->SetBranchAddress("theta_41_H_new",&theta_41[1]);
		inputHTree->SetBranchAddress("phi_300_H_new",&phi_300[1]);
		inputHTree->SetBranchAddress("phi_41_H_new",&phi_41[1]);

		int isCal;
		int isSoft;
		int isShort;
		int isCW;
		int isNewBox;

		inputAllTree->SetBranchAddress("cal",&isCal);
		inputAllTree->SetBranchAddress("soft",&isSoft);
		inputAllTree->SetBranchAddress("short",&isShort);
		inputAllTree->SetBranchAddress("CW",&isCW);
		inputAllTree->SetBranchAddress("box",&isNewBox);

		int isSurf[2]; // a surface event after filtering?
		int isSurfEvent_top[2]; // a top event?

		inputAllTree->SetBranchAddress("surf_V_new",&isSurf[0]);
		inputAllTree->SetBranchAddress("surf_H_new",&isSurf[1]);

		inputAllTree->SetBranchAddress("surf_top_V",&isSurfEvent_top[0]);
		inputAllTree->SetBranchAddress("surf_top_H",&isSurfEvent_top[1]);

		int isBadEvent;
		double weight;
		int unixTime;
		int isFirstFiveEvent;
		int hasBadSpareChanIssue;
		int hasBadSpareChanIssue2;
		int runNum;
		int eventNumber;

		inputAllTree->SetBranchAddress("bad",&isBadEvent);
		inputAllTree->SetBranchAddress("weight",&weight);
		inputAllTree->SetBranchAddress("unixTime",&unixTime);
		inputAllTree->SetBranchAddress("isFirstFiveEvent",&isFirstFiveEvent);
		inputAllTree->SetBranchAddress("hasBadSpareChanIssue",&hasBadSpareChanIssue);
		inputAllTree->SetBranchAddress("hasBadSpareChanIssue2",&hasBadSpareChanIssue2);
		inputAllTree->SetBranchAddress("runNum",&runNum);
		inputAllTree->SetBranchAddress("eventNumber",&eventNumber);

		double frac_of_power_notched_V[8];
		double frac_of_power_notched_H[8];

		stringstream ss;
		for(int i=0; i<8; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			inputVTree->SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_V[i]);
		}
		for(int i=8; i<16; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			inputHTree->SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_H[i-8]);
		}

		double VPeakOverRMS[16];
		int waveformLength[16];
		inputFilterTree->SetBranchAddress("VPeakOverRMS",&VPeakOverRMS);
		inputFilterTree->SetBranchAddress("waveformLength",&waveformLength);

		bool passesV;
		bool passesH;
		bool passesBoth;
		inputSurfaceTree->SetBranchAddress("passesV",&passesV);
		inputSurfaceTree->SetBranchAddress("passesH",&passesH);
		inputSurfaceTree->SetBranchAddress("passesBoth",&passesBoth);

		int numEntries = inputVTree->GetEntries();

		inputAllTree->GetEvent(0);

		bool isThisABadRun = isBadRun(station,runNum,BadRunList);
		if(isThisABadRun){
			inputFile->Close();
			delete inputFile;
			continue;
		}
		bool isThisABadSurfaceRun = isBadRun(station,runNum,BadSurfaceRunList);
		if(isThisABadSurfaceRun){
			inputFile->Close();
			delete inputFile;
			continue;
		}		

		int numHPolIsol=0;
		int numVPolIsol=0;
		int numCombo=0;

		for(int event=0; event<numEntries; event++){
			inputVTree->GetEvent(event);
			inputHTree->GetEvent(event);
			inputAllTree->GetEvent(event);
			inputRecoTree->GetEvent(event);
			inputFilterTree->GetEvent(event);
			inputSurfaceTree->GetEvent(event);

			if(passesV && !passesH){
				numVPolIsol++;
				// h1_events_vs_time[0]->Fill(unixTime);
				// h1_events_vs_time_finer[0]->Fill(unixTime);
				// h1_events_vs_run[0]->Fill(runNum);
				// PlotThisEvent(station, runNum, eventNumber, 0);
			}
			if(passesH && !passesV){
				numHPolIsol++;
				// h1_events_vs_time[1]->Fill(unixTime);
				// h1_events_vs_time_finer[1]->Fill(unixTime);
				// h1_events_vs_run[1]->Fill(runNum);
				// PlotThisEvent(station, runNum, eventNumber, 1);
			}

			if(passesBoth){
				// h1_events_vs_time[2]->Fill(unixTime);
				// h1_events_vs_time_finer[2]->Fill(unixTime);
				// h1_events_vs_run[2]->Fill(runNum);
				numCombo++;
				// PlotThisEvent(station, runNum, eventNumber, 2);
			}
		}// loop over events

		// printf("Run %4d, Num VPol Isol, Num HPol Isol, Num Combo: %3d, %3d, %3d \n", runNum, numVPolIsol, numHPolIsol, numCombo);

		// if(numHPolIsol==1 && numVPolIsol==0 && numCombo==0){
		if(numHPolIsol==0 && numVPolIsol==1 && numCombo==0){
			// printf("HPol Isol: Run %d, V, H, Combo: %d, %d, %d \n", runNum, numVPolIsol, numHPolIsol, numCombo);
			
			// now, loop over again and plot
			for(int event=0; event<numEntries; event++){
				// printf("%d, %d \n", runNum, eventNumber);
				// printf("Run %4d, Event %6d, Pol %d, unixTime %d, theta %3d, phi %3d, coor %.4f, CW status %d\n", runNum, eventNumber, 1, unixTime, theta_300[1], phi_300[1], corr_val[1], Refilt[1]);
				// fprintf(fout,"Run %4d, Event %6d, Pol %1d, unixTime %d, theta %3d, phi %3d, coor %.4f, refilt %d, surfv %d, surfh %d, surftopppol %d \n",runNum, eventNumber, pol, unixTime, theta_300[pol], phi_300[pol], corr_val[pol], Refilt[pol], isSurf[0], isSurf[1], isSurfEvent_top[pol]);
				inputAllTree->GetEvent(event);
				inputVTree->GetEvent(event);
				inputHTree->GetEvent(event);


				if(phi_300[0]<-140 || phi_300[0]>-90){
					printf("Run %d, Event %d, Theta %d, Phi %d \n", runNum, eventNumber, theta_300[0], phi_300[0]);
					// cout<<"About to call the CSW maker for this phi-selected event"<<endl;
					makeCSW(station, runNum, eventNumber, config, theta_300[0], phi_300[0]);
				}


				// getTravelTime(settings, raysolver, stationCenter, theta_300[1], phi_300[1], target);
				// double getTravelTime(Settings *settings, RaySolver *raysolver, double stationCenter[3], int thetaIn, int phiIn, Position target){
				// double thisZenith = getAntennaArrivalTheta(station, settings, target, raysolver, stationCenter, theta_300[0], phi_300[0]);
			}
		}

		inputFile->Close();
		delete inputFile;

	} // loop over input files

	char thistitle[300];
}


// do the freaking interpolation myself...
// ffttools is so dumb for trying so hard to be helpful...
TGraph *customInterpolation(TGraph *grIn)
{
	std::vector<double> tVec;
	std::vector<double> vVec;

	Int_t numIn=grIn->GetN();
	Double_t tIn,vIn;

	Double_t startTime=0;
	Double_t lastTime=0;
	for (int samp=0;samp<numIn;samp++) {
		grIn->GetPoint(samp,tIn,vIn);
		tVec.push_back(tIn);
		vVec.push_back(vIn);
		if(samp==0) startTime=tIn;
		lastTime=tIn;
	}
	if(tVec.size()<1) {
		std::cout << "Insufficent points for interpolation\n";
		return NULL;
	}

	ROOT::Math::Interpolator chanInterp(tVec,vVec,ROOT::Math::Interpolation::kAKIMA);

	vector<double> newTimes;
	vector<double> newVolts;
	for(double time=-500.; time<1000.; time+=0.2){
		newTimes.push_back(time);
		newVolts.push_back(chanInterp.Eval(time));
	}
	TGraph *grInt = new TGraph(newTimes.size(), &newTimes[0], &newVolts[0]);
	return grInt;

	// Int_t roughPoints=Int_t((lastTime-startTime)/deltaT);
	// Double_t *newTimes = new Double_t[roughPoints+1]; //Will change this at some point, but for now
	// Double_t *newVolts = new Double_t[roughPoints+1]; //Will change this at some point, but for now
	// Int_t numPoints=0;
	// for(Double_t time=startTime;time<=lastTime;time+=deltaT) {
	// 	newTimes[numPoints]=time;
	// 	newVolts[numPoints]=chanInterp.Eval(time);
	// 	numPoints++;
	// }

	// TGraph *grInt = new TGraph(numPoints,newTimes,newVolts);
	// delete [] newTimes;
	// delete [] newVolts;
	// return grInt;
}

void makeCSW(int station, int runNum, int eventNumber, int config, int thetaIn, int phiIn){

	/*
		In this modified version, we want to just do the cross-correlation to get delays
	*/

	/*
		get waveforms
	*/
	char *DataDirPath(getenv("DATA_DIR_100"));
	if (DataDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char run_file_name[400];
	sprintf(run_file_name,"%s/RawData/A%d/all_runs/event%d.root",DataDirPath,station,runNum);
	TFile *mapFile = TFile::Open(run_file_name);
	if(!mapFile){
		cout<<"Can't open data file for map!"<<endl;
		// return -1;
	}
	TTree *eventTree = (TTree*) mapFile-> Get("eventTree");
	if(!eventTree){
		cout<<"Can't find eventTree for map"<<endl;
		// return -1;
	}

	RawAtriStationEvent *rawPtr =0;
	eventTree->SetBranchAddress("event",&rawPtr);
	eventTree->GetEvent(eventNumber);

	int stationID = rawPtr->stationId;
	char ped_file_name[400];
	sprintf(ped_file_name,"%s/run_specific_peds/A%d/all_peds/event%d_specificPeds.dat",PedDirPath,station,runNum);
	AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	calibrator->setAtriPedFile(ped_file_name,stationID); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist

	UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawPtr,AraCalType::kLatestCalib);

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

	// now, we get delays

	double start_x;
	double stop_x;
	int center_chan;
	double extra_shift=0.;
	if(runNum==1514){
		start_x=100;
		stop_x=300.;
		center_chan=2;
	}
	if(runNum==1518){
		start_x=100;
		stop_x=300.;
		center_chan=2;
	}
	if(runNum==1629){
		start_x=0;
		stop_x=100.;
		center_chan=2;
	}
	else if(runNum==1647){
		start_x=0;
		stop_x=100.;
		center_chan=2;
	}
	else if(runNum==3232){
		start_x=-50;
		stop_x=50.;
		center_chan=2;
	}
	else if(runNum==4317 || runNum==4386){
		start_x=100.;
		stop_x=300.;
		center_chan=2;
	}
	else if(runNum==4855){
		start_x=150.;
		stop_x=250.;
		center_chan=2;
	}
	else if(runNum==4947){
		start_x=-100.;
		stop_x=200.;
		center_chan=6;
	}
	else if(runNum==4994 || runNum==4997 || runNum==5002 || runNum==5047 || runNum==5062 || runNum==5094){
		start_x=100.;
		stop_x=300.;
		center_chan=2;
	}
	else if(runNum==5052){
		start_x=100.;
		stop_x=300.;
		center_chan=2;
	}
	else if(runNum==5446 || runNum==5600){
		start_x=-100.;
		stop_x=100.;
		center_chan=2;
	}
	else if(runNum==5799 || runNum==5861){
		start_x=200.;
		stop_x=400.;
		center_chan=2;
	}
	else if(runNum==5946){
		start_x=0.;
		stop_x=200.;
		center_chan=5;
	}
	else if(runNum==6544){
		start_x=0.;
		stop_x=200.;
		center_chan=6;
	}
	else if(runNum==6588){
		start_x=100.;
		stop_x=300.;
		center_chan=5;
	}
	else if(runNum==6628){
		start_x=100.;
		stop_x=300.;
		center_chan=2;
	}
	else if(runNum==6657){
		start_x=200.;
		stop_x=400.;
		center_chan=2;
	}
	else if(runNum==6687){
		start_x=0.;
		stop_x=200.;
		center_chan=2;
	}
	else if(runNum==6759 || runNum==7836){
		start_x=100.;
		stop_x=300.;
		center_chan=2;
	}
	else if(runNum==6820 || runNum==6860 || runNum==7264){
		start_x=200.;
		stop_x=400.;
		center_chan=2;
	}
	else if(runNum==7418){
		start_x=100.;
		stop_x=300.;
		center_chan=4;
	}
	else if(runNum==7494){
		start_x=400.;
		stop_x=600.;
		center_chan=5;
	}
	else{
		start_x=-50;
		stop_x=50.;
		center_chan=2;
	}

	vector<double> delays;
	for(int ant=0; ant<16; ant++){
		TGraph *corr = FFTtools::getInterpolatedCorrelationGraph(waveforms[center_chan],waveforms[ant],0.1);
		int peak_bin = FFTtools::getPeakBin(corr);
		double this_delay = corr->GetX()[peak_bin];
		delays.push_back(this_delay);
		delete corr;
	}

	// now, we translate

	vector <TGraph*> waveforms_translated;
	for(int chan=0; chan<16; chan++){
		waveforms_translated.push_back((TGraph*)waveforms[chan]->Clone());
		int N = waveforms_translated[chan]->GetN();
		for(int samp=0; samp<N; samp++){
			waveforms_translated[chan]->GetX()[samp]+=delays[chan];
			// if(chan==10) waveforms_translated[chan]->GetX()[samp]-=6.797;
		}
	}

	/*
		Try this differently....
	*/
	vector<TGraph*> interpolatedWaveforms;
	for(int chan=0; chan<16; chan++){
		interpolatedWaveforms.push_back(FFTtools::getInterpolatedGraph(waveforms_translated[chan],0.2));
	}
	vector<TGraph*> paddedWaves;
	vector<TGraph*> croppedWaves;
	for(int chan=0; chan<16; chan++){
		paddedWaves.push_back(FFTtools::padWave(interpolatedWaveforms[chan],7));
		croppedWaves.push_back(FFTtools::cropWave(paddedWaves[chan],-1000,1500));
	}
	vector<TGraph*> interpAgain;
	for(int chan=0; chan<16; chan++){
		interpAgain.push_back(customInterpolation(croppedWaves[chan]));
	}

	cout<<"First samp chan 10 "<<interpAgain[10]->GetX()[10]<<endl;
	cout<<"First samp chan 12 "<<interpAgain[12]->GetX()[10]<<endl;

	vector<double> CSW_x;
	vector<double> CSW_y;
	for(double start=-500.; start<1000.; start+=0.2){
		CSW_x.push_back(start);
		CSW_y.push_back(0.);
	}
	for(int samp=0; samp<CSW_x.size(); samp++){
		for(int chan=0; chan<8; chan++){
			CSW_y[samp]+=interpAgain[chan]->GetY()[samp];
		}
		CSW_y[samp]/=8.;
	}
	TGraph *CSW = new TGraph(CSW_x.size(), &CSW_x[0], &CSW_y[0]);

	// and one more time translation
	TGraph *CSW_tosave = (TGraph*)CSW->Clone();
	int peak_bin = FFTtools::getPeakBin(CSW_tosave);
	double to_shift = CSW_tosave->GetX()[peak_bin];
	cout<<"Time to shift is "<<to_shift<<endl;
	int N = CSW_tosave->GetN();
	// shift the CSV back to zero, with an allowed 20ns of padding
	for(int samp=0; samp<N; samp++){
		CSW_tosave->GetX()[samp]-=to_shift-20.-extra_shift;
	}
	char title_txt[200];
	sprintf(title_txt,"/users/PAS0654/osu0673/A23_analysis_new2/AraRoot/analysis/a23_analysis_programs/proposal2019/c%d/V_A%d_run%d_ev%d.csv",config,station,runNum,eventNumber);
	FILE *fout = fopen(title_txt,"w");
	for(int samp=0; samp<N; samp++){
		fprintf(fout,"%.2f,%.2f\n",CSW_tosave->GetX()[samp],CSW_tosave->GetY()[samp]);
	}
	fclose(fout);//close sigmavsfreq.txt file

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
	// sprintf(save_temp_title,"%s/unblind/surface/mixed_events/surface_%d.%d.%d_A%d_Run%d_Ev%d_ProblemPol%d_Waveforms.png",plotPath,year_now,month_now,day_now,station,runNum,event,problempol);
	sprintf(save_temp_title, "/users/PAS0654/osu0673/A23_analysis_new2/AraRoot/analysis/a23_analysis_programs/proposal2019/c%d/V_A%d_Run%d_Even%d_Waveforms.png",config,station, runNum, eventNumber);
	TCanvas *cWave = new TCanvas("","",4*1100,4*850);
	cWave->Divide(4,4);
	for(int i=0; i<16; i++){
		cWave->cd(i+1);
		dummy[i]->Draw("AP");
		dummy[i]->SetLineColor(kWhite);
		dummy[i]->GetXaxis()->SetRangeUser(-200.,700.);
		dummy[i]->GetXaxis()->SetRangeUser(-700.,700.);
		// dummy[i]->GetXaxis()->SetRangeUser(140.,200.);

		waveforms[i]->Draw("sameL");
		waveforms[i]->SetLineWidth(3);

		// paddedWaves[i]->Draw("sameL");
		// paddedWaves[i]->SetLineColor(kRed);

		waveforms_translated[i]->Draw("sameL");
		waveforms_translated[i]->SetLineColor(kRed);
		waveforms_translated[i]->SetLineWidth(3);
	}
	cWave->SaveAs(save_temp_title);
	delete cWave;


	int colors [8] = {kRed, kOrange, kGreen, kBlue, kViolet, kCyan, kMagenta, kGray};

	char this_plot_title[150];
	sprintf(this_plot_title,"Run %d, Event %d",runNum,eventNumber);
	dummy[0]->SetTitle(this_plot_title);

	char this_plot_title_CSW[150];
	sprintf(this_plot_title_CSW,"Run %d, Event %d CSW",runNum,eventNumber);
	dummy[1]->SetTitle(this_plot_title_CSW);

	TCanvas *cstaggered = new TCanvas("","",1100,2*850);
	cstaggered->Divide(1,2);
	cstaggered->cd(1);
		dummy[0]->Draw("AP");
		dummy[0]->GetXaxis()->SetRangeUser(start_x,stop_x);
		dummy[0]->GetYaxis()->SetRangeUser(-500.,500.);
		dummy[0]->GetXaxis()->SetTitle("Time (ns)");
		dummy[0]->GetYaxis()->SetTitle("Voltage (mV)");
		for(int i=0; i<8; i++){
			waveforms_translated[i]->Draw("sameL");
			waveforms_translated[i]->SetLineColor(colors[i]);
			if(i==center_chan){
				waveforms_translated[i]->SetLineColor(kBlack);
			}
		}
		{
			TLegend *leg = new TLegend(0.58,0.75,0.9,0.9);
			leg->AddEntry(waveforms_translated[0],"Chan 0 (TV1)","l");
			leg->AddEntry(waveforms_translated[1],"Chan 1 (TV2)","l");
			leg->AddEntry(waveforms_translated[2],"Chan 2 (TV3)","l");
			leg->AddEntry(waveforms_translated[3],"Chan 3 (TV4)","l");
			leg->AddEntry(waveforms_translated[4],"Chan 4 (BV1)","l");
			leg->AddEntry(waveforms_translated[5],"Chan 5 (BV2)","l");
			leg->AddEntry(waveforms_translated[6],"Chan 6 (BV3)","l");
			leg->AddEntry(waveforms_translated[7],"Chan 7 (BV4)","l");
			leg->Draw();
		}
	cstaggered->cd(2);
		dummy[1]->Draw("AP");
		dummy[1]->GetXaxis()->SetRangeUser(start_x,stop_x);
		dummy[1]->GetYaxis()->SetRangeUser(-500.,500.);
		dummy[1]->GetXaxis()->SetTitle("Time (ns)");
		dummy[1]->GetYaxis()->SetTitle("Voltage (mV)");
		CSW->Draw("sameL");
		CSW->SetLineWidth(2);
	sprintf(save_temp_title, "/users/PAS0654/osu0673/A23_analysis_new2/AraRoot/analysis/a23_analysis_programs/proposal2019/c%d/V_A%d_Run%d_Even%d_WaveformsTimeTranslated.png",config,station, runNum, eventNumber);
	cstaggered->SaveAs(save_temp_title);
	delete cstaggered;


	for(int i=0; i<16; i++) delete waveforms[i];
	delete realAtriEvPtr;

}
