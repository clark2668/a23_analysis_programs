////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	Try to work out the bias in reconstruction from correlation function
////////////////////////////////////////////////////////////////////////////////

//C++
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <sys/stat.h>
#include <time.h>

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TRandom3.h"

//AraRoot includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "Settings.h"
#include "Detector.h"
#include "Report.h"
#include "RayTraceCorrelator.h"
AraAntPol::AraAntPol_t Vpol = AraAntPol::kVertical;
AraAntPol::AraAntPol_t Hpol = AraAntPol::kHorizontal;
#include "tools_PlottingFns.h"
#include "tools_RecoFns.h"
#include "tools_Cuts.h"

using namespace std;

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"


TGraph *getCorrelationGraph_timeDomain_N(TGraph *gr1, TGraph *gr2){
	int length1= gr1->GetN(); //get the number of points in graph 1
	int length2= gr2->GetN(); //get the number of points in graph 2

	if(length1!=length2){cout<<"The Waveforms are not of the same length! This many not work..."<<endl;}

	double x1,y1=0.; //starting values for the first data point
	double x2,y2=0.; //starting values for the second data point
	
	int padLength = 3*length1; // x3 to allow for one whole waveform to go before, middle, and after during the correlation
	
	TGraph *gr1Padded = FFTtools::padWaveToLength(gr1,padLength); //pad them to the same length for ease
	TGraph *gr2Padded = FFTtools::padWaveToLength(gr2,padLength); //pad them to the same length for ease

	gr1Padded->GetPoint(1,x2,y2); //get this first point
	gr1Padded->GetPoint(0,x1,y1); //get this second point
	double deltaT = x2-x1; //find the time step for the waveform
	gr2Padded->GetPoint(0,x2,y2);
	double waveOffset = x1-x2; //get the time offset of wave 2 from wave 1
	//if(abs(waveOffset)>0.0001){cout<<"The two graphs don't have the same starting point, this could be a problem..."<<endl;}
	//double maxDelay = deltaT * length1;
	int maxDelay = length1;
	cout<<"MaxDelay is "<<maxDelay<<endl;
	
	Double_t *oldX1 = gr1Padded->GetX();
	Double_t *oldY1 = gr1Padded->GetY();
	
	Double_t *oldX2 = gr2Padded->GetX();
	Double_t *oldY2 = gr2Padded->GetY();
	
	//Double_t *time_lag[padLength]={};
	//Double_t *correlationValues[padLength]={};

	vector <double> time_lag;
	vector <double> correlationValues;
	vector <double> numContributing;
	
	//cout<<"     MaxDelay is "<<maxDelay<<endl;
	//cout<<"     Length1 is "<<length1<<endl;

	for(int delay = -maxDelay; delay<maxDelay; delay++){
	//for(int delay = -maxDelay; delay<-5998; delay++){
		//cout<<"I'm working on delay "<<delay<<endl;
		double sum=0.;
		int num_contributing_points = 0;
		double total1=0.;
		double total2=0.;
		for(int i=length1; i<2*length1; i++){//integrate over just the center region of waveform 1
		//for(int i=length1; i<6002; i++){//integrate over just the center region of waveform 1
			//cout<<"I'm working on i= "<<i<<endl;
			int j=i+delay; //but, integrate over all possible delayed regions of waveform 2
			//cout<<"     I'm working on j= "<<j<<endl;
			double yVal1 = oldY1[i];
			total1+=oldY1[i]*oldY1[i];
			//cout<<"          First Y value is "<<yVal1<<endl;
			double yVal2 = oldY2[j];
			total2+=oldY2[j]*oldY2[j];
			//cout<<"          Second Y value is "<<yVal2<<endl;
			double product = yVal1*yVal2; //deltaT added 8/5
			//cout<<"          The product of the two is "<<product<<endl;
			//if(abs(yVal1)>0. || abs(yVal2)>0.) num_contributing_points++;
			// if(abs(yVal1)>0.0001) num_contributing_points++;
			// if(abs(yVal2)>0.0001) num_contributing_points++;
			if((abs(yVal1)>0.0001)&&(abs(yVal2)>0.0001)){ num_contributing_points++;} //increase the number of counted overlap
			//cout<<"The number of contributing points now is "<<num_contributing_points<<endl;
			sum+=product;
		}
		//cout<<"For time delay "<<delay<<" I'm working with time lag "<<(-maxDelay*deltaT)<<" and sum value "<<sum/((double) num_contributing_points)<<endl;
		time_lag.push_back(-delay*deltaT);
		// if(num_contributing_points>0){correlationValues.push_back(sum/((double) num_contributing_points)); numContributing.push_back((double) num_contributing_points);}
		if(num_contributing_points>0){correlationValues.push_back(sum/((double) sqrt(total1)*sqrt(total2))); numContributing.push_back((double) num_contributing_points);}
		else {correlationValues.push_back(0.) /*num_contributing_points=1*/; numContributing.push_back(0.);}
		//else{correlationValues.push_back(sum/((double) num_contributing_points));}
	}
	
	if(time_lag.size()!=correlationValues.size()){cout<<"Mismatch in size of time lag and correlation values, something is wrong!"<<endl;}
	
	TGraph *output = new TGraph(time_lag.size(),&time_lag[0],&correlationValues[0]);

	delete gr1Padded;
	delete gr2Padded;

	//delete output;
	return output;
}

double *custom_getCorrelation(int length,double *oldY1, double *oldY2){
	FFTWComplex *theFFT1=FFTtools::doFFT(length,oldY1);
	FFTWComplex *theFFT2=FFTtools::doFFT(length,oldY2);
	int newLength=(length/2)+1;
	FFTWComplex *tempStep = new FFTWComplex [newLength];
	int no2=length>>1;
	for(int i=0;i<newLength;i++) {
		double reFFT1=theFFT1[i].re;
		double imFFT1=theFFT1[i].im;
		double reFFT2=theFFT2[i].re;
		double imFFT2=theFFT2[i].im;
		//Real part of output 
		tempStep[i].re=(reFFT1*reFFT2+imFFT1*imFFT2);
		//Imaginary part of output 
		tempStep[i].im=(imFFT1*reFFT2-reFFT1*imFFT2);
	}
	double *theOutput=FFTtools::doInvFFT(length,tempStep);
	delete [] theFFT1;
	delete [] theFFT2;
	delete [] tempStep;
	return theOutput;
}

TGraph *custom_getCorrelationGraph(TGraph *gr1, TGraph *gr2){
	//Now we'll extend this up to a power of 2
	int length=gr1->GetN();
	int length2=gr2->GetN();

	int N=int(TMath::Power(2,int(TMath::Log2(length))+2));
	if(N<length2)
		N=int(TMath::Power(2,int(TMath::Log2(length2))+2));

	//Will really assume that N's are equal for now
	int firstRealSamp=(N-length)/2;

	double *oldY1 = new double [N];
	double *oldY2 = new double [N];

	double x,y;
	Double_t x2,y2;
	gr1->GetPoint(1,x2,y2);
	gr1->GetPoint(0,x,y);
	double deltaT=x2-x;
	double firstX=x;

	gr2->GetPoint(0,x2,y2);
	double waveOffset=firstX-x2;

	for(int i=0;i<N;i++) {
		if(i<firstRealSamp || i>=firstRealSamp+length)
			y=0;
		else {
			gr1->GetPoint(i-firstRealSamp,x,y);
		}
		oldY1[i]=y;

		if(i<firstRealSamp || i>=firstRealSamp+length2)
			y=0;
		else {
			gr2->GetPoint(i-firstRealSamp,x,y);
		}
		oldY2[i]=y;          
	}

	double *xVals = new double [N];
	double *yVals = new double [N];
	double *corVals=custom_getCorrelation(N,oldY1,oldY2);
	for(int i=0;i<N;i++) {
		if(i<N/2) {
			//Positive
			xVals[i+(N/2)]=(i*deltaT)+waveOffset;
			yVals[i+(N/2)]=corVals[i];
		}
		else {
			//Negative
			xVals[i-(N/2)]=((i-N)*deltaT)+waveOffset;
			yVals[i-(N/2)]=corVals[i];	  
		}
	}

	TGraph *grCor = new TGraph(N,xVals,yVals);
	delete [] oldY1;
	delete [] oldY2;
	delete [] xVals;
	delete [] yVals;
	delete [] corVals;
	return grCor;
}

TGraph *custom_getNormalisedCorrelationGraph(TGraph *gr1, TGraph *gr2, Int_t *zeroOffset) {
  //Will also assume these graphs are zero meaned... may fix this assumption
   //Now we'll extend this up to a power of 2
  int length=gr1->GetN();
  Double_t *y1=gr1->GetY();
  int length2=gr2->GetN();
  Double_t *y2=gr2->GetY();
  Double_t denom=gr1->GetRMS(2)*gr2->GetRMS(2);
  
  int N=int(TMath::Power(2,int(TMath::Log2(length))+2));
  if(N<length2)
    N=int(TMath::Power(2,int(TMath::Log2(length2))+2));
  
  //Will really assume that N's are equal for now
  int firstRealSamp=1+(N-2*length)/2;
  int lastRealSamp=firstRealSamp+2*(length-1);
  TGraph *grCor = FFTtools::getCorrelationGraph(gr1,gr2,zeroOffset);
  // TGraph *grCor = custom_getCorrelationGraph(gr1,gr2);
  Double_t *corVal=grCor->GetY();
  Double_t norm1=0;
  Double_t norm2=0;
  
  for(int i=0;i<N;i++) {
    if(i>=firstRealSamp && i<=lastRealSamp) { //if you are where there's trace
      if(i<=N/2) {
		norm1+=(y1[i-firstRealSamp]*y1[i-firstRealSamp]);
		norm2+=(y2[length-1-(i-firstRealSamp)]*y2[length-1-(i-firstRealSamp)]);
		int effN=1+(i-firstRealSamp);
		corVal[i]/=(sqrt(effN)*denom);
		printf("I is %d and effN is %d \n", i, effN);
		// corVal[i]/=(sqrt(effN)*sqrt(norm1)*sqrt(norm2));
      }
      else if(i<N-1) {
		norm1-=(y1[i-1-(N/2)]*y1[i-1-(N/2)]);
		norm2-=(y2[length-(i-(N/2))]*y2[length-(i-(N/2))]);
		int effN=(1+lastRealSamp-i);
		corVal[i]/=(sqrt(effN)*denom);
		// corVal[i]/=(sqrt(effN)*sqrt(norm1)*sqrt(norm2));
      }
    }
  }

  return grCor;
}


int main(int argc, char **argv)
{

	if(argc<5){
		cout<< "Usage\n" << argv[0] << " <station> <year> <run num> <event>"<<endl;;
		return -1;
	}
	int station = atoi(argv[1]);
	int year = atoi(argv[2]);
	int runNum = atoi(argv[3]);
	int event = atoi(argv[4]);
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time_now_time = localtime(&time_now);
	int year_now = time_now_time -> tm_year + 1900;
	int month_now = time_now_time -> tm_mon + 1;
	int day_now = time_now_time -> tm_mday;

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;
	char *DataDirPath(getenv("DATA_DIR"));
	if (DataDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;
	char *PedDirPath(getenv("PED_DIR"));
	if (PedDirPath == NULL) std::cout << "Warning! $DATA_DIR is not set!" << endl;

	gStyle->SetOptStat(11);

	// char run_file_name[400];
	// sprintf(run_file_name,"%s/RawData/A%d/%d/sym_links/event%d.root",DataDirPath,station,year,runNum,runNum);
	// TFile *mapFile = TFile::Open(run_file_name);
	// if(!mapFile){
	// 	cout<<"Can't open data file for map!"<<endl;
	// 	return -1;
	// }
	// TTree *eventTree = (TTree*) mapFile-> Get("eventTree");
	// if(!eventTree){
	// 	cout<<"Can't find eventTree for map"<<endl;
	// 	return -1;
	// }

	// RawAtriStationEvent *rawPtr =0;
	// eventTree->SetBranchAddress("event",&rawPtr);
	// eventTree->GetEvent(event);

	// int stationID = rawPtr->stationId;
	// char ped_file_name[400];

	// if(year==2013){
	// 	sprintf(ped_file_name,"%s/run_specific_peds/A%d/%d/event%d_specificPeds.dat",PedDirPath,station,year,runNum);
	// }
	// else if(year==2014 || year==2015 || year==2016){
	// 	sprintf(ped_file_name,"%s/run_specific_peds/A%d/%d/event00%d_specificPeds.dat",PedDirPath,station,year,runNum);
	// }
	// AraEventCalibrator *calibrator = AraEventCalibrator::Instance();
	// calibrator->setAtriPedFile(ped_file_name,stationID); //because someone had a brain (!!), this will error handle itself if the pedestal doesn't exist
	
	// UsefulAtriStationEvent *realAtriEvPtr = new UsefulAtriStationEvent(rawPtr,AraCalType::kLatestCalib);

	// int unixTime = (int)rawPtr->unixTime;
	// int unixTimeUs =(int)rawPtr->unixTimeUs;
	// int timeStamp = (int)rawPtr->timeStamp;
	// printf("Unixtime is %d \n", unixTime);
	// printf("Unixtime microsecond is %d \n", unixTimeUs);
	// printf("timeStamp is %d \n", timeStamp);
	// printf("Event Number is %d \n", realAtriEvPtr->eventNumber);

	// stringstream ss1;
	// string xLabel, yLabel;
	// xLabel = "Time (ns)"; yLabel = "Voltage (mV)";
	// vector<string> titlesForGraphs;
	// for (int i = 0; i < 16; i++){
	// 	ss1.str("");
	// 	ss1 << "Channel " << i;
	// 	titlesForGraphs.push_back(ss1.str());
	// }
	// vector <TGraph*> waveforms = makeGraphsFromRF(realAtriEvPtr,16,xLabel,yLabel,titlesForGraphs);

	// vector<TGraph*> grWaveformsInt = makeInterpolatedGraphs(waveforms, 0.5, xLabel, yLabel, titlesForGraphs);
	// vector<TGraph*> grWaveformsPadded = makePaddedGraphs(grWaveformsInt, 0, xLabel, yLabel, titlesForGraphs);
	// xLabel = "Frequency (Hz)"; yLabel = "Power Spectral Density (mV/Hz)";
	// vector<TGraph*> grWaveformsPowerSpectrum = makePowerSpectrumGraphs(grWaveformsPadded, xLabel, yLabel, titlesForGraphs);

	gStyle->SetOptStat(0);
	bool doRecoAtAll=true;
	if(doRecoAtAll){

		//set up the ray tracer
		Settings *settings = new Settings();
		string setupfile = "setup.txt";
		settings->ReadFile(setupfile);
		cout << "Read " << setupfile << " file!" << endl;
		settings->NOFZ=1;
		Detector *detector=0;
		RayTraceCorrelator *theCorrelators[2];
		theCorrelators[0] =  new RayTraceCorrelator(station, 41., settings, 1, 4); //41 m, cal puser
		theCorrelators[1] =  new RayTraceCorrelator(station, 300., settings, 1, 4);//300 m, far reco

		bool doDumbTest=true;
		if(doDumbTest){
			TRandom3 *R = new TRandom3(time(0));
			TRandom3 *R2 = new TRandom3(time(0));
			vector<double> X1, X2;
			vector<double> Y1, Y2;
			for(int i=0; i<100; i++){
				if(i>=10 && i<50){
					Y1.push_back(R->Gaus(0,50));
					X1.push_back(i*0.5);
				}
				// else
				// 	Y1.push_back(0.);
			}
			for(int i=0; i<100; i++){
				if(i>=30 && i<70){
					Y2.push_back(R->Gaus(0,50));
					X2.push_back(i*0.5);
				}
				// else
				// 	Y1.push_back(0.);
			}
			// for(int i=0; i<X1.size(); i++){
			// 	printf("Wave 1 samp %d is %.2f\n", i,X1[i]);
			// }
			// for(int i=0; i<X2.size(); i++){
			// 	printf("Wave 2 samp %d is %.2f\n", i,X2[i]);
			// }
			TGraph *gr1 = new TGraph(X2.size(), &X2[0], &Y1[0]);
			TGraph *gr2 = new TGraph(X1.size(), &X1[0], &Y2[0]);
			// TGraph *gr2 = (TGraph*) gr1->Clone();
			// TGraph *corr = getCorrelationGraph_timeDomain_N(gr1,gr2);
			TGraph *corr2 = theCorrelators[0]->getCorrelationGraph_WFweight(gr1,gr2);
			// TGraph *corr = theCorrelators[0]->getNormalisedCorrelationGraph(gr1,gr2,0);	
			// TGraph *corr2 = custom_getNormalisedCorrelationGraph(gr1,gr2,0);
			TGraph *corr = FFTtools::getCorrelationGraph(gr1,gr2);
			// TGraph *corr = custom_getCorrelationGraph(gr1,gr2);
			// TGraph *corr2 = custom_getCorrelationGraph(gr1,gr2);
			
			// TGraph *corr = custom_getNormalisedCorrelationGraph(gr1,gr2,0);

			bool doBrianDaileyCorr=false;
			if(doBrianDaileyCorr){
				double *times1 = gr1->GetX();
				double *volts1 = gr1->GetY();
				int N = gr1->GetN();
				double *times2 = gr2->GetX();
				double *volts2 = gr2->GetY();
				double delta_t = times1[1] - times1[0];
				double bin_delay = 20./delta_t;
				double offset = (times2[0]-times1[0])/delta_t;
				double corr_Val=0.;
				double v1_sq=0.;
				double v2_sq=0.;
				int bin1;
				double bin2;
				int num_bins=0;
				double total_delay = bin_delay+offset;
				for(int i=-0; i<N; i++){
					bin1=i;
					bin2=round(i-total_delay);
					if(bin2 >0 && (int)bin2<N){
						corr_Val+=volts1[bin1]*volts2[(int)bin2];
					}
				}
			}

			bool testCorr=false;
			if(testCorr){
				//Now we'll extend this up to a power of 2
				int length=gr1->GetN();
				int length2=gr2->GetN();

				int N=int(TMath::Power(2,int(TMath::Log2(length))+2));
				if(N<length2)
					N=int(TMath::Power(2,int(TMath::Log2(length2))+2));

				int firstRealSamp=(N-length)/2;

				printf("N is %d \n", N);
				printf("Length 1 and 2 are %d and %d \n", length, length2);
				printf("First real samp is %d \n", firstRealSamp);

				double x,y;
				Double_t x2,y2;
				gr1->GetPoint(1,x2,y2);
				gr1->GetPoint(0,x,y);
				double deltaT=x2-x; //the interpolation size (assumed same for both)
				double firstX=x; //the first non-zero sample in gr1

				gr2->GetPoint(0,x2,y2);
				double waveOffset=firstX-x2; //time separation between start of gr1 and gr2; if gr2 starts *later*, this value is *negative*
				int OffsetBin = (int)(waveOffset/deltaT); //time separation (waveOffset) in units of bins

				printf("deltaT is %.2f \n", deltaT);
				printf("firstX is %.2f \n", firstX);
				printf("waveOffset is %.2f \n", waveOffset);
				printf("OffsetBin is %d \n", OffsetBin);

				double *oldY1 = new double [N];
				double *oldY2 = new double [N];
				double *oldX = new double [N];

				// this is a super funky way of simulatneously padding the graphs to a factor of 2
				// and aligning the graphs in time before the getCorrelation_NoNorm call
				for(int i=0;i<N;i++) {
					oldX[i] = i*deltaT;
	   
					if(i<firstRealSamp || i>=firstRealSamp+length){
						y=0;
					}
					else {
						gr1->GetPoint(i-firstRealSamp,x,y);
					}
					oldY1[i]=y;
	  
					if(i<firstRealSamp || i>=firstRealSamp+length2)
						y=0;
					else {
						gr2->GetPoint(i-firstRealSamp,x,y);
					}
	  				 oldY2[i]=y;		  
				}
				// for(int i=0; i<N; i++){
				// 	printf("Samp %d: %.2f, %.2f, %.2f \n", i, oldX[i], oldY1[i], oldY2[i]);
				// }
				double *corVals=theCorrelators[0]->getCorrelation_NoNorm(N,oldY1,oldY2);

				double *xVals = new double [N];
				double *yVals = new double [N];

				// now to normalize
				for(int i=0; i<N; i++){
					double Norm1=0.;
					double Norm2=0.;
					int dBin;
					if(i<N/2) { //Positive


						// if we're left of the middle of the new trace
						
						xVals[i+(N/2)]=(i*deltaT)+waveOffset;
						dBin = i+OffsetBin; 
						
						// dBin is now the place in the final correlation graph

						if (dBin<0) {
							// graph 2 *lags* graph 1, i.e., graph 1 *leads* graph 2
							for (int i=-(dBin); i<N; i++) {
								Norm1 += oldY1[i]*oldY1[i];
							}
							for (int i=0; i<N+(dBin); i++) {
								Norm2 += oldY2[i]*oldY2[i];
							}
						}
						else { // dBin >= 0
							// graph 2 *leads* graph 1, i.e., graph 1 *lags* graph 2
							// printf("Samp i %d has dBin %d \n", i, dBin);
							for (int i=0; i<N-(dBin); i++) {
								// printf("		Sub loop: samp %d\n",i);
								Norm1 += oldY1[i]*oldY1[i];
							}
							for (int i=(dBin); i<N; i++) {
								Norm2 += oldY2[i]*oldY2[i];
							}
						}
						printf("For sample i %d, the xVal is %.2f and the norms are %.2f and %.2f \n", i, xVals[i+(N/2)], Norm1, Norm2);

						// cout<<"		Positive! Norm1 : "<<Norm1<<", Norm2 : "<<Norm2<<endl;

						if ( Norm1>0. && Norm2>0. ) 
							yVals[i+(N/2)]=corVals[i] / (sqrt(Norm1)*sqrt(Norm2)); // / (sqrt(Noverlap1)*sqrt(Noverlap2)) ;
						else{
							yVals[i+(N/2)]=corVals[i];
						}
					}
					else {
						//Negative
						xVals[i-(N/2)]=((i-N)*deltaT)+waveOffset;
						dBin = i-N+OffsetBin;

						if (dBin<0) {
							for (int i=-(dBin); i<N; i++) {
								Norm1 += oldY1[i]*oldY1[i];
							}
							for (int i=0; i<N+(dBin); i++) {
								Norm2 += oldY2[i]*oldY2[i];
							}
						}
						else { // dBin >= 0
							for (int i=0; i<N-(dBin); i++) {
								Norm1 += oldY1[i]*oldY1[i];
							}
							for (int i=(dBin); i<N; i++) {
								Norm2 += oldY2[i]*oldY2[i];
							}
						}
						// printf("For sample i %d, the xVal is %.2f and the norms are %.2f and %.2f \n", i, xVals[i-(N/2)], Norm1, Norm2);

						if ( Norm1>0. && Norm2>0. ) 
							yVals[i-(N/2)]=corVals[i] / (sqrt(Norm1)*sqrt(Norm2)); // / (sqrt(Noverlap1)*sqrt(Noverlap2)) ;
						else 
							yVals[i-(N/2)]=corVals[i];	  
					}
				}

				double sum1=0.;
				double sum2=0;
				for(int samp=0; samp<gr1->GetN(); samp++) sum1+=pow(gr1->GetY()[samp],2.);
				for(int samp=0; samp<gr2->GetN(); samp++) sum2+=pow(gr2->GetY()[samp],2.);
				printf("Sum over are %.2f and %.2f \n", sum1, sum2);
			}

			// this website is where I got the idea
			// https://currents.soest.hawaii.edu/ocn_data_analysis/_static/SEM_EDOF.html
			bool otherCorrTest=true;
			if(otherCorrTest){

				TGraph *num_overlap = (TGraph*) corr->Clone();
				TGraph *gr1_integral = (TGraph*) corr->Clone();
				TGraph *gr2_integral = (TGraph*) corr->Clone();
				TGraph *product = (TGraph*) corr->Clone();
				TGraph *corr_copy = (TGraph*) corr->Clone();
				double RMS1 = gr1->GetRMS(2);
				double RMS2 = gr2->GetRMS(2);

				double t1i = gr1->GetX()[0];
				double t1f = gr1->GetX()[gr1->GetN()-1];
				double t2i = gr2->GetX()[0];
				double t2f = gr2->GetX()[gr2->GetN()-1];
				printf("Graph 1 start and stop: %.2f and %.2f \n", t1i, t1f);
				printf("Graph 2 start and stop: %.2f and %.2f \n", t2i, t2f);
				for(int corrsamp=0; corrsamp<corr->GetN(); corrsamp++){
					double lag, corrval;
					corr->GetPoint(corrsamp, lag, corrval);
					double t2i_new = t2i+lag;
					double t2f_new = t2f+lag;
					double integral_start=-1000000;
					double integral_stop=-500000;
					bool do_integral;
					printf(ANSI_COLOR_BLUE" Corr Sample %3d has lag %.2f and unnormalized corr value of %.2f \n "ANSI_COLOR_RESET, corrsamp, lag, corrval);
					printf("    Graph 2 now starts and stop at %.2f and %.2f \n", t2i_new, t2f_new);
					if(
						t2i_new < t1i 
						&& 
						t2f_new < t1f
					)
						{
						// printf("     Case A\n");
						integral_start = t1i;
						integral_stop = t2f_new;
						do_integral=true;
					}
					else if(
						t2i_new > t1i 
						&& 
						t2f_new > t1f
					)
						{
						// printf("     Case B\n");
						integral_start = t2i_new;
						integral_stop = t1f;
						do_integral=true;
					}
					else if(
						t2i_new > t1i
						&&
						t2f_new < t1f
					)
						{
						// printf("     Case C\n");
						integral_start = t2i_new;
						integral_stop = t2f_new;
						do_integral=true;

					}
					else if(
						t2i_new < t1i
						&&
						t2f_new > t1f
					){
						// printf("     Case D\n");
						integral_start = t1i;
						integral_stop = t1f;
						do_integral=true;
					}
					else if(
						(t2i_new-t1i)<0.0001
					){
						// printf("     Case E\n");
						integral_start = t1i;
						integral_stop = t1f;
					}
					double integral_gr1=0.;
					double integral_gr2=0.;
					int n_overlap_1=0;
					int n_overlap_2=0;
					if(do_integral){
						for(int samp1=0; samp1<gr1->GetN(); samp1++){
							double thisX, thisY;
							gr1->GetPoint(samp1,thisX,thisY);
							if(thisX>integral_start && thisX<integral_stop){
								integral_gr1+=thisY*thisY;
								n_overlap_1++;
							}
						}
						for(int samp2=0; samp2<gr2->GetN(); samp2++){
							double thisX, thisY;
							gr2->GetPoint(samp2, thisX, thisY);
							thisX+=lag;
							if(thisX>integral_start && thisX<integral_stop){
								integral_gr2+=thisY*thisY;
								n_overlap_2++;
							}
						}
						if(integral_gr1>0. && integral_gr2>0. && n_overlap_1>1){						
						// if(integral_gr1>0. && integral_gr2>0. && n_overlap_1>2){
							// corrval*=1./(sqrt(integral_gr1)*sqrt(integral_gr2));
							// corrval*=1./n_overlap_1;
							corrval*=1./(sqrt(n_overlap_1)*RMS1*RMS2);
							// corrval*=1./(n_overlap_1);
							printf(ANSI_COLOR_GREEN"     Normalizations are %.2f and %.2f with %d and %d overlapping \n"ANSI_COLOR_RESET, sqrt(integral_gr1), sqrt(integral_gr2), n_overlap_1, n_overlap_2);
						}
						else{
							corrval=0.;
							printf(ANSI_COLOR_RED"     Else case! Sample %d \n"ANSI_COLOR_RESET, corrsamp);
						}
						corr->SetPoint(corrsamp,lag,corrval);
						// printf("     Corrval now is %.3f \n", corrval);

						num_overlap->SetPoint(corrsamp,lag,n_overlap_1);
						gr1_integral->SetPoint(corrsamp,lag,sqrt(integral_gr1));
						gr2_integral->SetPoint(corrsamp,lag,sqrt(integral_gr2));
						product->SetPoint(corrsamp,lag,(sqrt(integral_gr1)*sqrt(integral_gr2)));
					}
					else{
						corr->SetPoint(corrsamp,lag,0.);						
						printf("Not doing integral!\n");
						num_overlap->SetPoint(corrsamp,lag,0.);
						gr1_integral->SetPoint(corrsamp,lag,0.);
						gr2_integral->SetPoint(corrsamp,lag,0.);
						product->SetPoint(corrsamp,lag,0.);
					}
				}

				TCanvas *cOverlap = new TCanvas("","",2*850,4*850);
				cOverlap->Divide(1,6);
				cOverlap->cd(1);
					num_overlap->Draw("ALP");
					num_overlap->GetYaxis()->SetTitle("Number of Overlapping Bins");
					num_overlap->GetXaxis()->SetTitle("Lag (ns)");
				cOverlap->cd(2);
					gr1_integral->Draw("ALP");
					gr1_integral->GetYaxis()->SetTitle("Graph 1 Normalization");
					gr1_integral->GetXaxis()->SetTitle("Lag (ns)");
				cOverlap->cd(3);
					gr2_integral->Draw("ALP");
					gr2_integral->GetYaxis()->SetTitle("Graph 2 Normalization");
					gr2_integral->GetXaxis()->SetTitle("Lag (ns)");
				cOverlap->cd(4);
					product->Draw("ALP");
					product->GetYaxis()->SetTitle("Product of Normalizations");
					product->GetXaxis()->SetTitle("Lag (ns)");
				cOverlap->cd(5);
					corr_copy->Draw("ALP");
					corr_copy->GetYaxis()->SetTitle("Raw Unnormalized Correlation Function");
					corr_copy->GetXaxis()->SetTitle("Lag (ns)");
				cOverlap->cd(6);
					corr->Draw("ALP");
					corr->GetYaxis()->SetTitle("Normalized Correlation Function");
					corr->GetXaxis()->SetTitle("Lag (ns)");
				char save_temp_title[400];
				sprintf(save_temp_title,"%s/corr_study/%d.%d.%d_VerifyOverlap_Dumb.png",plotPath,year_now,month_now,day_now);
				cOverlap->SaveAs(save_temp_title);
				delete cOverlap;
				delete num_overlap;
				delete gr1_integral;
				delete gr2_integral;
			}


			TGraph *rect_corr = FFTtools::rectifyWave(corr);
			int peak_bin = FFTtools::getPeakBin(rect_corr);
			printf("Peak Bin and Value is: %d and %.2f \n", peak_bin,corr->GetX()[peak_bin]);
			TCanvas *cTest = new TCanvas("","",850,2*850);

			vector<TGraph*> dummy;
			for(int i=0; i<3; i++){
				vector<double> thisX;
				vector<double> thisY;

				if(i<2){
					thisY.push_back(-150);
					thisY.push_back(150);
					thisX.push_back(-10.);
					thisX.push_back(70.);
				}
				else{
					thisY.push_back(-1.);
					thisY.push_back(1.);
					thisX.push_back(-50.);
					thisX.push_back(50.);
				}
				dummy.push_back(new TGraph(thisX.size(), &thisX[0], &thisY[0]));
				dummy[i]->GetXaxis()->SetTitle("Time (au)");
				if(i==2)
					dummy[i]->GetYaxis()->SetTitle("Cross Corr (dimless)");
				else
					dummy[i]->GetYaxis()->SetTitle("Amplitude (au)");
			}


			cTest->Divide(1,4);
			cTest->cd(1);
				dummy[0]->Draw("AP");
				gr1->Draw("Lsame");
				gr1->GetYaxis()->SetRangeUser(-150,150);
				dummy[0]->SetTitle("Graph 1");
			cTest->cd(2);
				dummy[1]->Draw("AP");
				gr2->Draw("Lsame");
				gr2->GetYaxis()->SetRangeUser(-150,150);
				dummy[1]->SetTitle("Graph 2");
			cTest->cd(3);
				corr2->Draw("ALP");
				corr2->SetTitle("Original Ray Tracer Normalized Correlation Graph");
				dummy[0]->SetTitle("Graph 1");
			cTest->cd(4);
				// dummy[2]->Draw("AP");
				corr->SetTitle("New Proposed Normalization");
				corr->Draw("ALsame");
				// rect_corr->Draw("Lsame");
			char save_temp_title[400];		
			sprintf(save_temp_title,"%s/corr_study/%d.%d.%d_Test_Correlation_Comparison.png",plotPath,year_now,month_now,day_now);
			cTest->SaveAs(save_temp_title);
			delete gr1;
			delete gr2;
			delete corr;
			delete rect_corr;

		}

		// bool do_reco=false;
		// if(do_reco){
		// 	TH2D *map_41m_V;
		// 	TH2D *map_300m_V;
		// 	TH2D *map_41m_H;
		// 	TH2D *map_300m_H;
		// 	TH2D *map_41m_V_select;

		// 	map_41m_V = theCorrelators[0]->getInterferometricMap_RT(settings, detector, realAtriEvPtr, Vpol, 0, 0);
		// 	map_300m_V = theCorrelators[1]->getInterferometricMap_RT(settings, detector, realAtriEvPtr, Vpol, 0, 0);
		// 	map_41m_H = theCorrelators[0]->getInterferometricMap_RT(settings, detector, realAtriEvPtr, Hpol, 0, 0);
		// 	map_300m_H = theCorrelators[1]->getInterferometricMap_RT(settings, detector, realAtriEvPtr, Hpol, 0, 0);

		// 	int PeakTheta_Recompute_41m_H;
		// 	int PeakTheta_Recompute_300m_H;
		// 	int PeakPhi_Recompute_41m_H;
		// 	int PeakPhi_Recompute_300m_H;
		// 	double PeakCorr_Recompute_41m_H;
		// 	double PeakCorr_Recompute_300m_H;
		// 	int PeakTheta_Recompute_41m_V;
		// 	int PeakTheta_Recompute_300m_V;
		// 	int PeakPhi_Recompute_41m_V;
		// 	int PeakPhi_Recompute_300m_V;
		// 	double PeakCorr_Recompute_41m_V;
		// 	double PeakCorr_Recompute_300m_V;
		// 	getCorrMapPeak(map_41m_H,PeakTheta_Recompute_41m_H,PeakPhi_Recompute_41m_H,PeakCorr_Recompute_41m_H);
		// 	getCorrMapPeak(map_300m_H,PeakTheta_Recompute_300m_H,PeakPhi_Recompute_300m_H,PeakCorr_Recompute_300m_H);
		// 	getCorrMapPeak(map_41m_V,PeakTheta_Recompute_41m_V,PeakPhi_Recompute_41m_V,PeakCorr_Recompute_41m_V);
		// 	getCorrMapPeak(map_300m_V,PeakTheta_Recompute_300m_V,PeakPhi_Recompute_300m_V,PeakCorr_Recompute_300m_V);

		// 	printf("	Rconstruction Information\n");
		// 	printf("		41m H theta and phi %d and %d \n", PeakTheta_Recompute_41m_H, PeakPhi_Recompute_41m_H);
		// 	stringstream ss30H;
		// 	ss30H<<" 41m H Peak Theta, Phi is "<<PeakTheta_Recompute_41m_H<<" , "<<PeakPhi_Recompute_41m_H;
		// 	map_41m_H->SetTitle(ss30H.str().c_str());
		// 	printf("		300m H theta and phi %d and %d \n", PeakTheta_Recompute_300m_H, PeakPhi_Recompute_300m_H);
		// 	stringstream ss300H;
		// 	ss300H<<" 300m H Peak Theta, Phi is "<<PeakTheta_Recompute_300m_H<<" , "<<PeakPhi_Recompute_300m_H;
		// 	map_300m_H->SetTitle(ss300H.str().c_str());
		// 	printf("		41m V theta and phi %d and %d \n", PeakTheta_Recompute_41m_V, PeakPhi_Recompute_41m_V);
		// 	stringstream ss30V;
		// 	ss30V<<" 41m V Peak Theta, Phi is "<<PeakTheta_Recompute_41m_V<<" , "<<PeakPhi_Recompute_41m_V;
		// 	map_41m_V->SetTitle(ss30V.str().c_str());
		// 	printf("		300m V theta and phi %d and %d \n", PeakTheta_Recompute_300m_V, PeakPhi_Recompute_300m_V);
		// 	stringstream ss300V;
		// 	ss300V<<" 300m V Peak Theta, Phi is "<<PeakTheta_Recompute_300m_V<<" , "<<PeakPhi_Recompute_300m_V;
		// 	map_300m_V->SetTitle(ss300V.str().c_str());

		// 	TCanvas *cMaps = new TCanvas("","",2*1100,2*850);
		// 	cMaps->Divide(2,2);
		// 		cMaps->cd(3);
		// 		map_41m_V->Draw("colz");
		// 		cMaps->cd(4);
		// 		map_41m_H->Draw("colz");
		// 		cMaps->cd(1);
		// 		map_300m_V->Draw("colz");
		// 		cMaps->cd(2);
		// 		map_300m_H->Draw("colz");
		// 	char save_temp_title[400];		
		// 	sprintf(save_temp_title,"%s/single_events/%d.%d.%d_Run%d_Ev%d_Maps.png",plotPath,year_now,month_now,day_now,runNum,event);
		// 	cMaps->SaveAs(save_temp_title);
		// 	delete cMaps;
		// 	delete map_41m_V; delete map_300m_V; delete map_41m_H; delete map_300m_H; 
		// 	// delete map_41m_V_select;
		// }


		// bool do_reco_snrweighted=false;
		// if(do_reco_snrweighted){
		// 	TH2D *map_41m_V;
		// 	TH2D *map_300m_V;
		// 	TH2D *map_41m_H;
		// 	TH2D *map_300m_H;
		// 	TH2D *map_41m_V_select;

		// 	map_41m_V = theCorrelators[0]->getInterferometricMap_RT_SNRweighted(settings, detector, realAtriEvPtr, Vpol, 0, 0);
		// 	map_300m_V = theCorrelators[1]->getInterferometricMap_RT_SNRweighted(settings, detector, realAtriEvPtr, Vpol, 0, 0);
		// 	map_41m_H = theCorrelators[0]->getInterferometricMap_RT_SNRweighted(settings, detector, realAtriEvPtr, Hpol, 0, 0);
		// 	map_300m_H = theCorrelators[1]->getInterferometricMap_RT_SNRweighted(settings, detector, realAtriEvPtr, Hpol, 0, 0);

		// 	int PeakTheta_Recompute_41m_H;
		// 	int PeakTheta_Recompute_300m_H;
		// 	int PeakPhi_Recompute_41m_H;
		// 	int PeakPhi_Recompute_300m_H;
		// 	double PeakCorr_Recompute_41m_H;
		// 	double PeakCorr_Recompute_300m_H;
		// 	int PeakTheta_Recompute_41m_V;
		// 	int PeakTheta_Recompute_300m_V;
		// 	int PeakPhi_Recompute_41m_V;
		// 	int PeakPhi_Recompute_300m_V;
		// 	double PeakCorr_Recompute_41m_V;
		// 	double PeakCorr_Recompute_300m_V;
		// 	getCorrMapPeak(map_41m_H,PeakTheta_Recompute_41m_H,PeakPhi_Recompute_41m_H,PeakCorr_Recompute_41m_H);
		// 	getCorrMapPeak(map_300m_H,PeakTheta_Recompute_300m_H,PeakPhi_Recompute_300m_H,PeakCorr_Recompute_300m_H);
		// 	getCorrMapPeak(map_41m_V,PeakTheta_Recompute_41m_V,PeakPhi_Recompute_41m_V,PeakCorr_Recompute_41m_V);
		// 	getCorrMapPeak(map_300m_V,PeakTheta_Recompute_300m_V,PeakPhi_Recompute_300m_V,PeakCorr_Recompute_300m_V);

		// 	printf("	Rconstruction Information\n");
		// 	printf("		41m H theta and phi %d and %d \n", PeakTheta_Recompute_41m_H, PeakPhi_Recompute_41m_H);
		// 	stringstream ss30H;
		// 	ss30H<<" 41m H Peak Theta, Phi is "<<PeakTheta_Recompute_41m_H<<" , "<<PeakPhi_Recompute_41m_H;
		// 	map_41m_H->SetTitle(ss30H.str().c_str());
		// 	printf("		300m H theta and phi %d and %d \n", PeakTheta_Recompute_300m_H, PeakPhi_Recompute_300m_H);
		// 	stringstream ss300H;
		// 	ss300H<<" 300m H Peak Theta, Phi is "<<PeakTheta_Recompute_300m_H<<" , "<<PeakPhi_Recompute_300m_H;
		// 	map_300m_H->SetTitle(ss300H.str().c_str());
		// 	printf("		41m V theta and phi %d and %d \n", PeakTheta_Recompute_41m_V, PeakPhi_Recompute_41m_V);
		// 	stringstream ss30V;
		// 	ss30V<<" 41m V Peak Theta, Phi is "<<PeakTheta_Recompute_41m_V<<" , "<<PeakPhi_Recompute_41m_V;
		// 	map_41m_V->SetTitle(ss30V.str().c_str());
		// 	printf("		300m V theta and phi %d and %d \n", PeakTheta_Recompute_300m_V, PeakPhi_Recompute_300m_V);
		// 	stringstream ss300V;
		// 	ss300V<<" 300m V Peak Theta, Phi is "<<PeakTheta_Recompute_300m_V<<" , "<<PeakPhi_Recompute_300m_V;
		// 	map_300m_V->SetTitle(ss300V.str().c_str());

		// 	TCanvas *cMaps = new TCanvas("","",2*1100,2*850);
		// 	cMaps->Divide(2,2);
		// 		cMaps->cd(3);
		// 		map_41m_V->Draw("colz");
		// 		cMaps->cd(4);
		// 		map_41m_H->Draw("colz");
		// 		cMaps->cd(1);
		// 		map_300m_V->Draw("colz");
		// 		cMaps->cd(2);
		// 		map_300m_H->Draw("colz");
		// 	char save_temp_title[400];		
		// 	sprintf(save_temp_title,"%s/single_events/%d.%d.%d_Run%d_Ev%d_Maps_SNRweighted.png",plotPath,year_now,month_now,day_now,runNum,event);
		// 	cMaps->SaveAs(save_temp_title);
		// 	delete cMaps;
		// 	delete map_41m_V; delete map_300m_V; delete map_41m_H; delete map_300m_H; 
		// 	// delete map_41m_V_select;
		// }

		// bool doContributingMaps=false;
		// if(doContributingMaps){
		// 	stringstream ss1;
		// 	vector<string> titlesForGraphs;
		// 	vector<TH2D*> individuals;
		// 	vector<TH2D*> h2NumOverlappingBins;
		// 	vector<TGraph*> indiv_corrs;
		// 	for(int i=0; i<7; i++){
		// 		for(int j=i+1; j<8; j++){
		// 			ss1.str("");
		// 			ss1<<"Pair "<<i<<" and "<<j;
		// 			titlesForGraphs.push_back(ss1.str());
		// 			h2NumOverlappingBins.push_back(theCorrelators[0]->getNumBinsMap_RT_PairSelect(settings, detector, realAtriEvPtr, Vpol, 0, i, j, 0));
		// 			// individuals.push_back(theCorrelators[0]->getInterferometricMap_RT_PairSelect(settings, detector, realAtriEvPtr, Vpol, 0, i, j, 0));
		// 			// indiv_corrs.push_back(theCorrelators[0]->getCorrelationGraph_WFweight(grWaveformsInt[i],grWaveformsInt[j]));
		// 			// indiv_corrs.push_back(FFTtools::getCorrelationGraph(grWaveformsInt[i],grWaveformsInt[j]));
		// 			// int *blah=0;
		// 			// indiv_corrs.push_back(theCorrelators[0]->getNormalisedCorrelationGraph(grWaveformsInt[i],grWaveformsInt[j],blah));
		// 			// TGraph *waveform1_padded = FFTtools::padWaveToLength(grWaveformsInt[i], grWaveformsInt[i]->GetN()+2000);
		// 			// TGraph *waveform1_cropped =FFTtools::cropWave(waveform1_padded,-300.,500.);
		// 			// TGraph *waveform2_padded = FFTtools::padWaveToLength(grWaveformsInt[j], grWaveformsInt[j]->GetN()+20000);
		// 			// TGraph *waveform2_cropped =FFTtools::cropWave(waveform2_padded,-300.,500.);
		// 			// indiv_corrs.push_back(getCorrelationGraph_timeDomain_N(waveform1_cropped,waveform2_cropped));
		// 			// indiv_corrs.push_back(getCorrelationGraph_timeDomain_N(grWaveformsInt[i],grWaveformsInt[j]));
					
		// 		}
		// 	}

		// 	vector<double> mins;
		// 	vector<double> maxs;
		// 	for(int i=0; i<h2NumOverlappingBins.size(); i++){
		// 		mins.push_back(h2NumOverlappingBins[i]->GetMinimum());
		// 		maxs.push_back(h2NumOverlappingBins[i]->GetMaximum());
		// 	}
		// 	std::sort(mins.begin(), mins.end()); //sort smallest to largest
		// 	std::sort(maxs.begin(), maxs.end()); //sort smallest to largest
		// 	std::reverse(maxs.begin(), maxs.end()); //reverse order to get largest to smallest

		// 	TCanvas *cOverlapBins = new TCanvas("","",8*850,4*850);
		// 	cOverlapBins->Divide(7,4);
		// 	for(int i=0; i<h2NumOverlappingBins.size(); i++){
		// 		cOverlapBins->cd(i+1);
		// 		h2NumOverlappingBins[i]->Draw("colz");
		// 		h2NumOverlappingBins[i]->GetZaxis()->SetRangeUser(mins[0],maxs[0]);
		// 		h2NumOverlappingBins[i]->SetTitle(titlesForGraphs[i].c_str());
		// 		h2NumOverlappingBins[i]->GetXaxis()->SetTitle("Phi (deg)");
		// 		h2NumOverlappingBins[i]->GetYaxis()->SetTitle("Theta (deg)");
		// 		h2NumOverlappingBins[i]->GetZaxis()->SetTitle("Num Overlapping Bins");
		// 	}
		// 	char save_temp_title[400];
		// 	sprintf(save_temp_title,"%s/corr_study/%d.%d.%d_Run%d_Ev%d_NumOverlappingBins_AllMaps.png",plotPath,year_now,month_now,day_now,runNum,event);
		// 	cOverlapBins->SaveAs(save_temp_title);
		// 	delete cOverlapBins;

		// 	mins.clear();
		// 	maxs.clear();
		// 	for(int i=0; i<individuals.size(); i++){
		// 		mins.push_back(individuals[i]->GetMinimum());
		// 		maxs.push_back(individuals[i]->GetMaximum());
		// 	}
		// 	std::sort(mins.begin(), mins.end()); //sort smallest to largest
		// 	std::sort(maxs.begin(), maxs.end()); //sort smallest to largest
		// 	std::reverse(maxs.begin(), maxs.end()); //reverse order to get largest to smallest			

		// 	TCanvas *cMaps = new TCanvas("","",8*850,4*850);
		// 	cMaps->Divide(7,4);
		// 	for(int i=0; i<individuals.size(); i++){
		// 		cMaps->cd(i+1);
		// 		individuals[i]->Draw("colz");
		// 		individuals[i]->GetZaxis()->SetRangeUser(mins[0],maxs[0]);
		// 		individuals[i]->SetTitle(titlesForGraphs[i].c_str());
		// 	}
		// 	char save_temp_title[400];
		// 	sprintf(save_temp_title,"%s/single_events/%d.%d.%d_Run%d_Ev%d_AllMaps.png",plotPath,year_now,month_now,day_now,runNum,event);
		// 	cMaps->SaveAs(save_temp_title);
		// 	delete cMaps;

		// 	vector<TGraph*> dummy;
		// 	for(int i=0; i<individuals.size(); i++){
		// 		vector<double> thisX;
		// 		vector<double> thisY;
		// 		thisX.push_back(-600.);
		// 		thisX.push_back(600.);
		// 		thisY.push_back(-0.20);
		// 		thisY.push_back(0.20);
		// 		dummy.push_back(new TGraph(thisX.size(), &thisX[0], &thisY[0]));
		// 		dummy[i]->GetXaxis()->SetTitle("Time (ns)");
		// 		dummy[i]->GetYaxis()->SetTitle("Cross Corr (dimless)");
		// 		// dummy[i]->GetXaxis()->SetLabelSize(0.07);
		// 		// dummy[i]->GetYaxis()->SetLabelSize(0.07);
		// 		// dummy[i]->GetXaxis()->SetTitleSize(0.07);
		// 		// dummy[i]->GetYaxis()->SetTitleSize(0.07);
		// 	}

		// 	TCanvas *cCorr = new TCanvas("","",8*850,4*850);
		// 	cCorr->Divide(7,4);
		// 	for(int i=0; i<indiv_corrs.size(); i++){
		// 		cMaps->cd(i+1);
		// 		// dummy[i]->Draw("AP");
		// 		// dummy[i]->SetTitle(titlesForGraphs[i].c_str());
		// 		indiv_corrs[i]->Draw("ALP");
		// 	}
		// 	sprintf(save_temp_title,"%s/single_events/%d.%d.%d_Run%d_Ev%d_CorrFunctions.png",plotPath,year_now,month_now,day_now,runNum,event);
		// 	cCorr->SaveAs(save_temp_title);
		// 	delete cCorr;
		// }
	}

	// vector<TGraph*> dummy;
	// for(int i=0; i<16; i++){
	// 	vector<double> thisX;
	// 	vector<double> thisY;
	// 	thisX.push_back(-200.);
	// 	thisX.push_back(400.);
	// 	thisY.push_back(-100.);
	// 	thisY.push_back(100.);
	// 	dummy.push_back(new TGraph(thisX.size(), &thisX[0], &thisY[0]));
	// 	dummy[i]->GetXaxis()->SetTitle("Time (ns)");
	// 	dummy[i]->GetYaxis()->SetTitle("Voltage (mV)");
	// 	dummy[i]->GetXaxis()->SetLabelSize(0.07);
	// 	dummy[i]->GetYaxis()->SetLabelSize(0.07);
	// 	dummy[i]->GetXaxis()->SetTitleSize(0.07);
	// 	dummy[i]->GetYaxis()->SetTitleSize(0.07);
	// }

	// char save_temp_title[300];
	// sprintf(save_temp_title,"%s/single_events/%d.%d.%d_Run%d_Ev%d_Waveforms.png",plotPath,year_now,month_now,day_now,runNum,event);
	// TCanvas *cWave = new TCanvas("","",4*1100,4*850);
	// cWave->Divide(4,4);
	// for(int i=0; i<16; i++){
	// 	cWave->cd(i+1);
	// 	dummy[i]->Draw("AP");
	// 	waveforms[i]->Draw("Lsame");
	// 	waveforms[i]->SetLineWidth(2);
	// }
	// cWave->SaveAs(save_temp_title);
	// delete cWave;

	// sprintf(save_temp_title,"%s/single_events/%d.%d.%d_Run%d_Ev%d_Spectra.png",plotPath,year_now,month_now,day_now,runNum,event);
	// TCanvas *cSpec = new TCanvas("","",4*1100,4*850);
	// cSpec->Divide(4,4);
	// for(int i=0; i<16; i++){
	// 	cSpec->cd(i+1);
	// 	grWaveformsPowerSpectrum[i]->Draw("AL");
	// 	grWaveformsPowerSpectrum[i]->SetLineWidth(3);
	// 	gPad->SetLogy();
	// 	grWaveformsPowerSpectrum[i]->GetYaxis()->SetRangeUser(10.,1e7);
	// }
	// cSpec->SaveAs(save_temp_title);
	// delete cSpec;

	// for(int i=0; i<16; i++){
	// 	delete waveforms[i];
	// 	delete grWaveformsInt[i];
	// 	delete grWaveformsPadded[i];
	// 	delete grWaveformsPowerSpectrum[i];
	// }
	// delete realAtriEvPtr;
	// mapFile->Close();
	// delete mapFile;
	return 0;
}
