////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	noise_csw_spectra.cxx
////	make spectra of csw
////
////	Apr 2020
////	WARNING BY THE WAY, THIS HAS ALLLL THE MEMORY LEAKS
////	THIS WAS QUICK AND DIRTY AND IS NOT HOW THINGS SHOULD BE DONE!
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
#include "FFTtools.h"

using namespace std;

TGraph* makeCSW(int station, int runNum, int eventNumber, int thetaIn, int phiIn, int &wanToUse);

int main(int argc, char **argv)
{
	if(argc<4){
		cout<< "Usage\n" << argv[0] << " <1-station> <2-runNum> <3-eventNumber>"<<endl;;
		return -1;
	}
	int station = atoi(argv[1]);
	int runNum = atoi(argv[2]);
	int eventNumber = atoi(argv[3]);

	vector<TGraph*> toAverage;
	for(int i=eventNumber-2; i>eventNumber-100; i--){
		int wantToUse=0;
		cout<<"Targetting event "<<eventNumber<<endl;
		TGraph *result = makeCSW(station,runNum, i, 0, 0, wantToUse);
		if(wantToUse){
			toAverage.push_back(result);
		}
	}


	for(int samp=0; samp<toAverage[0]->GetN(); samp++){
		double t1, v1;
		double toAdd=0.;
		toAverage[0]->GetPoint(samp, t1, v1);
		for(int event=1; event<toAverage.size(); event++){
			toAdd+=toAverage[event]->GetY()[samp];
		}
		v1+=toAdd;
		v1/=double(toAverage.size());
		toAverage[0]->SetPoint(samp,t1,v1);
	}
	char title_txt[200];
	sprintf(title_txt,"/users/PAS0654/osu0673/A23_analysis_new2/AraRoot/analysis/a23_analysis_programs/proposal2019/cr_stuff/plot_spectra/noise_spec.csv");
	FILE *fout = fopen(title_txt,"w");
	for(int samp=0; samp<toAverage[0]->GetN(); samp++){
		fprintf(fout,"%.2f,%.2f\n",toAverage[0]->GetX()[samp],toAverage[0]->GetY()[samp]);
	}
	fclose(fout);//close sigmavsfreq.txt file

}	


// do the freaking interpolation myself...
// ffttools is so dumb for trying so hard to be helpful...
TGraph *customInterpolation(TGraph *grIn, double time_step=0.2)
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
	for(double time=-500.; time<1000.; time+=time_step){
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

TGraph *makeCSW(int station, int runNum, int eventNumber, int thetaIn, int phiIn, int &wantToUse){

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
	if(rawPtr->isCalpulserEvent() || rawPtr->isSoftwareTrigger()) wantToUse=0;
	else wantToUse=1;

	stringstream ss1;
	string xLabel, yLabel;
	xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
	vector<string> titlesForGraphs;
	for (int i = 0; i < 16; i++){
		ss1.str("");
		ss1 << "Channel " << i;
		titlesForGraphs.push_back(ss1.str());
	}

	vector<TGraph*> waveforms;
	for(int i=0; i<16; i++){
		waveforms.push_back(realAtriEvPtr->getGraphFromRFChan(i));
	}

	// vector <TGraph*> waveforms = makeGraphsFromRF(realAtriEvPtr,16,xLabel,yLabel,titlesForGraphs);

	// now, we get delays

	double start_x=150.;
	double stop_x=250.;
	int center_chan=12;
	double extra_shift=0.;
	if(runNum==1571){
		start_x=150.;
		stop_x=250.;
		center_chan=12;
	}
	if(runNum==1455){
		start_x=150.;
		stop_x=250.;
		center_chan=11;
	}
	if(runNum==2779){
		start_x=0.;
		stop_x=100.;
		center_chan=10;
	}
	if(runNum==2871){
		start_x=0.;
		stop_x=100.;
		center_chan=10;
	}
	if(runNum==2937){
		start_x=0.;
		stop_x=100.;
		center_chan=10;
	}
	if(runNum==3202){
		start_x=0.;
		stop_x=100.;
		center_chan=10;
	}
	if(runNum==3206){
		start_x=0.;
		stop_x=100.;
		center_chan=10;
	}
	if(runNum==3325){
		start_x=-50.;
		stop_x=50.;
		center_chan=10;
	}
	if(runNum==3392){
		start_x=-50.;
		stop_x=50.;
		center_chan=10;
	}
	if(runNum==3529){
		start_x=0.;
		stop_x=100.;
		center_chan=10;
	}
	if(runNum==6674){
		start_x=200;
		stop_x=400.;
		center_chan=8;
	}
	if(runNum==6861){
		start_x=300;
		stop_x=400.;
		center_chan=9;
	}
	if(runNum==7170){
		start_x=250;
		stop_x=350.;
		center_chan=9;
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
		for(int chan=8; chan<15; chan++){
			// if(chan==8) continue;
			CSW_y[samp]+=interpAgain[chan]->GetY()[samp];
		}
		CSW_y[samp]/=7.;
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
	// char title_txt[200];
	// sprintf(title_txt,"/users/PAS0654/osu0673/A23_analysis_new2/AraRoot/analysis/a23_analysis_programs/proposal2019/cr_stuff/plot_spectra/A%d_run%d_ev%d_wave.csv",station,runNum,eventNumber);
	// FILE *fout = fopen(title_txt,"w");
	// for(int samp=0; samp<N; samp++){
	// 	fprintf(fout,"%.2f,%.2f\n",CSW_tosave->GetX()[samp],CSW_tosave->GetY()[samp]);
	// }
	// fclose(fout);//close sigmavsfreq.txt file

	vector<double> thisX_spec;
	vector<double> thisY_spec;
	thisX_spec.push_back(0.);
	thisX_spec.push_back(1200.);
	thisY_spec.push_back(-10);
	thisY_spec.push_back(40);
	TGraph *dummy_spec = new TGraph(2, &thisX_spec[0], &thisY_spec[0]);

	// plot spectrum
	TGraph *CSW_downsample = customInterpolation(CSW_tosave,0.5);
	TGraph *CSW_downsample_trimmed = FFTtools::cropWave(CSW_downsample, -50,100);
	int length_to_pad_to = CSW_downsample_trimmed->GetN();
	length_to_pad_to = pow(2, ceil(log(length_to_pad_to)/log(2))); // round to power of two
	TGraph *csw_padded = FFTtools::padWaveToLength(CSW_downsample_trimmed,length_to_pad_to);
	TGraph *csw_spec = FFTtools::makePowerSpectrumMilliVoltsNanoSecondsdB(csw_padded);
	TCanvas *cSpec = new TCanvas("spec","spec",1100,850);
		dummy_spec->Draw("AP");
		dummy_spec->SetLineColor(kWhite);
		csw_spec->Draw("sameL");
	char save_temp_title[300];
	sprintf(save_temp_title, "/users/PAS0654/osu0673/A23_analysis_new2/AraRoot/analysis/a23_analysis_programs/proposal2019/cr_stuff/plot_spectra/A%d_Run%d_Even%d_Spec.png",station, runNum, eventNumber);
	// cSpec->SaveAs(save_temp_title);


	// char title_txt[200];
	// sprintf(title_txt,"/users/PAS0654/osu0673/A23_analysis_new2/AraRoot/analysis/a23_analysis_programs/proposal2019/cr_stuff/plot_spectra/A%d_run%d_ev%d_spec.csv",station,runNum,eventNumber);
	// FILE *fout = fopen(title_txt,"w");
	// for(int samp=0; samp<csw_spec->GetN(); samp++){
	// 	fprintf(fout,"%.2f,%.2f\n",csw_spec->GetX()[samp],csw_spec->GetY()[samp]);
	// }
	// fclose(fout);//close sigmavsfreq.txt file

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

	// sprintf(save_temp_title,"%s/unblind/surface/mixed_events/surface_%d.%d.%d_A%d_Run%d_Ev%d_ProblemPol%d_Waveforms.png",plotPath,year_now,month_now,day_now,station,runNum,event,problempol);
	sprintf(save_temp_title, "/users/PAS0654/osu0673/A23_analysis_new2/AraRoot/analysis/a23_analysis_programs/proposal2019/cr_stuff/plot_spectra/A%d_Run%d_Even%d_Waveforms.png",station, runNum, eventNumber);
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
	// cWave->SaveAs(save_temp_title);
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
		for(int i=8; i<15; i++){
			// if(i==8) continue;
			waveforms_translated[i]->Draw("sameL");
			waveforms_translated[i]->SetLineColor(colors[i-8]);
			if(i==center_chan){
				waveforms_translated[i]->SetLineColor(kBlack);
			}
		}
		{
			TLegend *leg = new TLegend(0.58,0.75,0.9,0.9);
			leg->AddEntry(waveforms_translated[8],"Chan 8 (TH1)","l");
			leg->AddEntry(waveforms_translated[9],"Chan 9 (TH2)","l");
			leg->AddEntry(waveforms_translated[10],"Chan 10 (TH3)","l");
			leg->AddEntry(waveforms_translated[11],"Chan 11 (TH4)","l");
			leg->AddEntry(waveforms_translated[12],"Chan 12 (BH1)","l");
			leg->AddEntry(waveforms_translated[13],"Chan 13 (BH2)","l");
			leg->AddEntry(waveforms_translated[14],"Chan 14 (BH3)","l");
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
	sprintf(save_temp_title, "/users/PAS0654/osu0673/A23_analysis_new2/AraRoot/analysis/a23_analysis_programs/proposal2019/cr_stuff/plot_spectra/A%d_Run%d_Even%d_WaveformsTimeTranslated.png",station, runNum, eventNumber);
	// cstaggered->SaveAs(save_temp_title);
	delete cstaggered;


	for(int i=0; i<16; i++) delete waveforms[i];
	delete realAtriEvPtr;

	return csw_spec;

}
