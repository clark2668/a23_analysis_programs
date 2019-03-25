#include "TF1.h"
#include "TGraph.h"
using namespace std;

void Exchange( double &a, double &b ) {
	double tmp = a;
	a = b;
	b = tmp;
}

void four1(double *data, const int isign, int nsize) {
	int nn,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;

	int n=nsize/2;

	nn=n << 1;
	j=1;
	for (i=1;i<nn;i+=2) {
		if (j > i) {
			Exchange(data[j-1],data[i-1]);
			Exchange(data[j],data[i]);
		}
		m=n;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (nn > mmax) {
		istep=mmax << 1;
		theta=isign*(TMath::Pi()*2./mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=nn;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j-1]-wi*data[j];
				tempi=wr*data[j]+wi*data[j-1];
				data[j-1]=data[i-1]-tempr;
				data[j]=data[i]-tempi;
				data[i-1] += tempr;
				data[i] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

void realft(double *data, const int isign, int nsize){
	int i, i1, i2, i3, i4;
	double c1=0.5,c2,h1r,h1i,h2r,h2i,wr,wi,wpr,wpi,wtemp,theta;
	theta=TMath::Pi()/(double)(nsize>>1);
	if (isign == 1) {
		c2 = -0.5;
		four1(data,1,nsize);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	for (i=1;i<(nsize>>2);i++) {
		i2=1+(i1=i+i);
		i4=1+(i3=nsize-i1);
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r= -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4]= -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[0] = (h1r=data[0])+data[1];
		data[1] = h1r-data[1];
	} else {
		data[0]=c1*((h1r=data[0])+data[1]);
		data[1]=c1*(h1r-data[1]);
		four1(data,-1,nsize);
	}
}

void getDiodeModel(int NFOUR, double TIMESTEP, vector<double>&fdiode_real, vector<double>&diode_real) {

	double maxt_diode = 70.E-9;
	int maxt_diode_bin = (int)( maxt_diode / TIMESTEP );
	double idelaybeforepeak = (int)(13.E-9 / TIMESTEP);
	int iwindow = (int)(4.E-9 / TIMESTEP);
	int ibinshift = NFOUR/4 - (int)( maxt_diode / TIMESTEP );
	
	//  this is our homegrown diode response function which is a downgoing gaussian followed by an upward step function
	TF1 *fdown1=new TF1("fl_down1","[3]+[0]*exp(-1.*(x-[1])*(x-[1])/(2*[2]*[2]))",-300.E-9,300.E-9);
	fdown1->SetParameter(0,-0.8);
	fdown1->SetParameter(1,15.E-9);
	fdown1->SetParameter(2,2.3E-9);
	fdown1->SetParameter(3,0.);
	
	TF1 *fdown2=new TF1("fl_down2","[3]+[0]*exp(-1.*(x-[1])*(x-[1])/(2*[2]*[2]))",-300.E-9,300.E-9);
	fdown2->SetParameter(0,-0.2);
	fdown2->SetParameter(1,15.E-9);
	fdown2->SetParameter(2,4.0E-9);
	fdown2->SetParameter(3,0.);

	TF1 *f_up=new TF1("f_up","[0]*([3]*(x-[1]))^2*exp(-(x-[1])/[2])",-200.E-9,100.E-9);
	f_up->SetParameter(0,1.);
	f_up->SetParameter(1,18.E-9);
	f_up->SetParameter(2,7.0E-9);
	f_up->SetParameter(3,1.E9);
	f_up->SetParameter(0,-1.*sqrt(2.*TMath::Pi())*(fdown1->GetParameter(0)*fdown1->GetParameter(2)+fdown2->GetParameter(0)*fdown2->GetParameter(2))/(2.*pow(f_up->GetParameter(2),3.)*1.E18));

	// vector<double> diode_real;

	for (int i=0;i<NFOUR/2;i++) {
		diode_real.push_back(0.);   // first puchback 0. value  (this is actually not standard way though works fine)
		if (i<(int)(maxt_diode/TIMESTEP)) { // think this is same with above commented if
			diode_real[i]=fdown1->Eval((double)i*TIMESTEP)+fdown2->Eval((double)i*TIMESTEP);
			if (i>(int)(f_up->GetParameter(1)/TIMESTEP))
				diode_real[i]+=f_up->Eval((double)i*TIMESTEP);
		}
	}
	
	// diode_real is the time domain response of the diode
	// now get f domain response with realft
	double diode_real_fft[NFOUR];  // double sized array for myconvlv	
	for (int i=0; i<NFOUR; i++) {  // 512 bin added for zero padding
		if ( i<(int)(maxt_diode/TIMESTEP) ) {
			diode_real_fft[i] = diode_real[i];
		}
		else {
			diode_real_fft[i] = 0.;
		}
	}	

	// forward FFT
	realft(diode_real_fft,1,NFOUR);

	// save f domain diode response in fdiode_real
	for (int i=0; i<NFOUR; i++) {
		fdiode_real.push_back( diode_real_fft[i] );
	}
	delete fdown1;
	delete fdown2;
	delete f_up;
}

TGraph* doConvolve(TGraph *grIn){

	TGraph *grClone = (TGraph*)grIn->Clone();
	for (int i=0;i<grClone->GetN();i++) grClone->GetX()[i] *= 1.e-9; 

	int length=grClone->GetN();
	int NFOUR=length*2; //NFOUR is 2 x n bins for readout
	int BINSIZE=NFOUR/2;
	double TIMESTEP = grClone->GetX()[1]- grClone->GetX()[0];
	double maxt_diode = 70.E-9;
	int maxt_diode_bin = (int)( maxt_diode / TIMESTEP );

	// setup diode
	vector<double> fdiode_real; //ft of diode time domain response (in realft packaging)
	vector<double> diode_real; //time domain of diode response
	getDiodeModel(NFOUR, TIMESTEP, fdiode_real, diode_real);
	
	//this would print out the time-domain tunnel diode response
	/*
	vector<double> xes;
	for(int i=0; i<diode_real.size(); i++){
		xes.push_back(double(i));
	}
	TGraph *grOut = new TGraph(fdiode_real.size(), &xes[0], &fdiode_real[0]);
	TCanvas *c = new TCanvas("","",1100,850);
	grOut->Draw("ALP");
	c->SaveAs("fdiode_real.png");
	delete c;
	delete grOut;
	*/

	double Zr=50.;
	vector<double> data;
	vector<double> times;
	for(int samp=0; samp<length; samp++){
		data.push_back(grClone->GetY()[samp]);
		times.push_back(grClone->GetX()[samp]);
	}

	double power_noise_copy[length*2];
	for(int i=0; i<length; i++){
		power_noise_copy[i]=(data[i]*data[i])/Zr*TIMESTEP;
	}
	for(int i=length; i<length*2; i++){
		power_noise_copy[i]=0.;
	}
	realft(power_noise_copy, 1, length*2);

	double ans_copy[length*2];
	for(int j=0; j<length; j++){
		ans_copy[2*j]=(power_noise_copy[2*j]*fdiode_real[2*j]-power_noise_copy[2*j+1]*fdiode_real[2*j+1])/((double)length);
		ans_copy[2*j+1]=(power_noise_copy[2*j+1]*fdiode_real[2*j]+power_noise_copy[2*j]*fdiode_real[2*j+1])/((double)length);
	}
	ans_copy[0]=power_noise_copy[0]*fdiode_real[0]/((double)length);
	ans_copy[1]=power_noise_copy[1]*fdiode_real[1]/((double)length);

	realft(ans_copy,-1,length*2);

	vector<double> diodeconv;
	// for(int i=0; i<length+maxt_diode_bin; i++){
	for(int i=0; i<length; i++){
		diodeconv.push_back(ans_copy[i]);
	}

	vector<double> diodeX;
	for(int i=0; i<diodeconv.size(); i++){
		diodeX.push_back(grIn->GetX()[i]);
	}
	TGraph *grDiode = new TGraph(diodeX.size(), &diodeX[0], &diodeconv[0]);
	delete grClone;
	return grDiode;
}

