////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	fit_rayleigh.cxx
////
////	March 2019,  clark.2668@osu.edu
////	Do rayleigh fit
////////////////////////////////////////////////////////////////////////////////

// C/C++ Includes
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <algorithm>

//AraRoot Includes
#include "FFTtools.h"
#include "AraQualCuts.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"

using namespace std;

int main(int argc, char **argv)
{
	if(argc<4) {  // Check to make sure there are enough arguments to do something meaningful
		std::cout << "Usage requires you to provide input parameter of the form " << basename(argv[0]) << " <merged_file> <chan> <freq_bin>" << std::endl;
		return -1;
	}
	int chan=atoi(argv[2]);
	int freq_bin=atoi(argv[3]);
	
	TFile *fpIn = new TFile(argv[1], "OLD"); //we're going to open the data file
	if(!fpIn){
		std::cerr<< "Can not open the old file: " <<argv[2]<<endl;
		return -1;
	} //throw a warning if you can't open it
	
	fpIn->cd(); //go into that file
	TTree *inTree = (TTree*) fpIn->Get("outTree");
	if(!inTree){
		std::cerr<<"Can't find outTree (silly naming...) in file" <<argv[2]<<endl;
		return -1;
	} //throw a warning if you can't open it
	double chan_spec[16][512]={{0}};
	double freqs[16][512]={0};
	inTree->SetBranchAddress("chan_spec", &chan_spec);
	inTree->SetBranchAddress("freqs", &freqs);

	TH1D *h = new TH1D("","",1000,0.,2.);

	double numEntries = inTree -> GetEntries(); //get the number of entries in this file
	for(int event=0; event<numEntries; event++){ //loop over those entries
		inTree->GetEntry(event); //get the event
		h->Fill(chan_spec[chan][freq_bin]);
	}

	h->Sumw2();
	Double_t scale = 1/h->Integral("width");
	h->Scale(scale);

	int binmax = h->GetMaximumBin();//Get bin of max value
	double y_max = h->GetBinContent(binmax); //get value of that bin
	double x_max = h->GetXaxis()->GetBinCenter(binmax);
	TF1 *f = new TF1("f","([0]*x)/([1]*[1])*exp(-[0]*[0]*x*x/(2.*[1]*[1]))",0,h->GetXaxis()->GetBinCenter(binmax));
	f->SetParameters(h->GetMaximum(),h->GetRMS());
	gStyle->SetOptFit(1111);
	h->Fit("f");
	double p0 = f->GetParameter(0);
	double p1 = f->GetParameter(1);
	double chi2 = f->GetChisquare()/f->GetNDF(); //Reduced
	delete f;

	char title_txt[200];
	sprintf(title_txt,"/users/PAS0654/osu0673/A23_analysis_new2/AraRoot/analysis/a23_analysis_programs/other/rayleigh/sigmavsfreq_ch%d.txt",chan);
	FILE *fout = fopen(title_txt, "a");
	fprintf(fout,"%2.4f,%d,%2.4f,%2.4f  \n",freqs[chan][freq_bin],chan,p1,chi2);
	fclose(fout);//close sigmavsfreq.txt file

	TCanvas *c = new TCanvas("","",1100,850);
	h->Draw("E");
	char hist_title[300];
	sprintf(hist_title,"Freq Bin %d: %.2f ",freq_bin,freqs[chan][freq_bin]);
	h->SetTitle(hist_title);
	char save_title[300];
	sprintf(save_title,"fits_ch%d_f%d.png",chan,freq_bin);
	c->SaveAs(save_title);

	// fpIn->Close();
	// delete fpIn;

}//close the main program