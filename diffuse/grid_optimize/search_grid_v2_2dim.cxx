////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////    This time exhaustive search grid w/ S/Sup, where Sup is Sup for sum of backgrounds
////
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

// analysis custom
#include "tools_Cuts.h"
#include "tools_Stats.h"
#include "tools_CommandLine.h"
#include "tools_outputObjects.h"

double func(double *x, double *p){
	double num = x[0];
	double signal = p[0];
	double background = p[1];
	double term1 = TMath::Exp(-(signal+background));
	double term2 = TMath::Power(signal + background,num);
	double term3 = TMath::Gamma(num+1.);
	return term1*term2/term3;
}

// double total_livetime=1142.46; //total livetime of the experiment in days
// double live_frac[5]={
// 	0.15674,
// 	0.12477,
// 	0.08275,
// 	0.38428,
// 	0.25146
// };

double live_frac[5]={
	0.56,
	0.44,
	0.08275,
	0.38428,
	0.25146
};

double getSignal(double Sbin1, double Sbin2){
	Sbin1*=live_frac[0];
	Sbin2*=live_frac[1]*2.; //because Brian was dumb and threw half as many neutrinos in config 2 and config 3
	return Sbin1 + Sbin2;
}

using namespace std;

int main(int argc, char **argv)
{
	
	time_t time_now = time(0); //get the time now                                                                                                                                                                  
	tm *time = localtime(&time_now);
	int year_now = time -> tm_year + 1900;
	int month_now = time -> tm_mon + 1;
	int day_now = time -> tm_mday;

	stringstream ss;
	gStyle->SetOptStat(11);

	char *plotPath(getenv("PLOT_PATH"));
	if (plotPath == NULL) std::cout << "Warning! $PLOT_PATH is not set!" << endl;

	if(argc<4){
		cout<< "Usage\n" << argv[0] << " <1-station> <2-pol> <3-slope>"<<endl;;
		return -1;
	}

	// for(double back=10.; back>0.1; back-=0.1){
	// 	double achieved_alpha;
	// 	cout<<"back is "<<back<<" and sup is "<<GetS_up(back,achieved_alpha,0.9)<<endl;
	// }
	// return 0;

	int station = atoi(argv[1]);
	int pol_select=atoi(argv[2]);
	double slope=double(atof(argv[3]));


	if(station!=2 && station!=3){
		printf("No good! You asked for station %d, but this code only works for stations 2 and 3 \n",station);
		return -1;
	}

	char inputfilename[500];
	sprintf(inputfilename,"A%d_optimize_pol%d_slope%d.root",station,pol_select,int(slope));
	TFile *InputFIle = TFile::Open(inputfilename,"READ");
	TTree *InputTree = (TTree*) InputFIle->Get("OutputTree");

	double intercept[200];
	double backgrounds[200];
	double signal[200];
	double pVal;
	int fitConverge;

	InputTree->SetBranchAddress("intercept",&intercept);
	InputTree->SetBranchAddress("backgrounds",&backgrounds);
	InputTree->SetBranchAddress("signal",&signal);
	InputTree->SetBranchAddress("fitConverge",&fitConverge);
	InputTree->SetBranchAddress("pVal",&pVal);

	// loop over configs, and get the info stored into these larger "mega" arrays, so we don't have to keep asking root to get stuff when we do the grid search
	// c++/the computer has more than enough memory to hold all these numbers at once

	double backgrounds_all[5][200];
	double signal_all[5][200];
	double pVal_all[5];
	int fitConverge_all[5];

	bool allConverge=true;
	bool allGood=true;

	for(int config=0; config<5; config++){

		InputTree->GetEntry(config);
		for(int bin=0; bin<200; bin++){
			backgrounds_all[config][bin] = backgrounds[bin];
			signal_all[config][bin] = signal[bin];
		}
		pVal_all[config] = pVal;
		fitConverge_all[config] = fitConverge;
		if(fitConverge!=1){
			allConverge=false;
		}
		if(pVal<0.05 || pVal>0.99){
			allGood=false;
		}

	}
	// now, we need to check for fit convergence and good p-values in all configurations
	// if not, we stop
	allConverge=true;
	if(!allConverge){
		printf("Not all fits converged! abandon ship! Convergence was: %d, %d, %d, %d, %d \n",fitConverge_all[0], fitConverge_all[1], fitConverge_all[2], fitConverge_all[3], fitConverge_all[4]);
		return -1;
	}
	allGood=true;
	if(!allGood){
		printf("Not all p-vals are reasonable! abandon ship! Values were: %.4f, %.4f, %.4f, %.4f, %.4f \n",pVal_all[0], pVal_all[1], pVal_all[2], pVal_all[3], pVal_all[4]);
		return -1;
	}

	stringstream ss1;

	int best_bins[5] = {0};
	double best_S_over_Sup = -10000.;

	int numBins = 50;

	TH2D *S = new TH2D("","",50,10,15,50,10,15);
	TH2D *Sup = new TH2D("","",50,10,15,50,10,15);
	TH2D *S_over_Sup = new TH2D("","",50,10,15,50,10,15);

	// if all converge and all are good, start grid search (ugh);
	for(int bin1=0; bin1<numBins; bin1++){ // loop over config 1 bins
		
		double back1 = backgrounds_all[0][bin1];
		if(back1>10.) continue; // no background this large in single bin will ever be viable as part of final answer
		// printf("For bin1 %d, at intercept %.2f we find %.6f num backgrounds, and total integral of %.6f \n",bin1,intercept[bin1],back1, total_integrals[0]);
		// printf("For bin1 %d, at intercept %.2f, we find %.10E num backgrounds, and total integral of %.6f with signal %.2f \n", bin1,intercept[bin1],back1,total_integrals[0], signal_all[0][bin1]);
		
		// for(int bin2=0; bin2<numBins+20	; bin2++){ // loop over config 2 bins
		// for(int bin2=0; bin2<numBins+12; bin2++){ // loop over config 2 bins
		for(int bin2=0; bin2<numBins; bin2++){ // loop over config 2 bins


			double back2 = backgrounds_all[1][bin2];
			if(back2>10.) continue; // no background this large in single bin will ever be viable as part of final answer
			// printf("     For bin2 %d, at intercept %.2f, we find %.10E num backgrounds, and total integral of %.6f with signal %.2f \n",bin2,intercept[bin2],back2, total_integrals[1], signal_all[1][bin2]);

			// adaptive step size algorithm for locating the ideal sup
			// this should be much faster than a brute fore scan

			double achieved_alpha;
			double this_Sup = GetS_up(back1+back2,achieved_alpha,0.9);
			// printf("For back %.2f and %.2f, Sup is %.2f \n",back1,back2,this_Sup );

			Sup->SetBinContent(bin1,bin2,this_Sup);

			// cout<<"               I think I found the optimum with an Sup of "<<possibleSup<<" in only "<<counter<<" counts"<<endl;
			double totalSignal = getSignal(signal_all[0][bin1],signal_all[1][bin2]); // this function livetime adjusts the flux and fixes a numerical mistake brian made
			
			// S->SetBinContent(bin1,bin2,(signal_all[0][bin1]*live_frac[0])+(signal_all[1][bin2]*live_frac[1]*2.));
			S->SetBinContent(bin1,bin2,totalSignal);
			double this_S_over_Sup = totalSignal/this_Sup;
			
			S_over_Sup->SetBinContent(bin1,bin2,this_S_over_Sup);

			// printf("               Now, S/Sup = %.2f/%.2f = %.2f\n", totalSignal, possibleSup, this_S_over_Sup);
			if(this_S_over_Sup>best_S_over_Sup){
				// printf(GREEN"     Found new optimum! Old is %.4f and new is S/Sup = %.4f/%.4f = %.4f \n"RESET,best_S_over_Sup, totalSignal, possibleSup, this_S_over_Sup);
				best_S_over_Sup = this_S_over_Sup;
				best_bins[0]=bin1;
				best_bins[1]=bin2;
				// printf(GREEN"          New optimum has intercepts %.2f, %.2f\n"RESET, intercept[bin1], intercept[bin2]);
				// printf("     New optimum has intercepts %.2f, %.2f, %.2f, %.2f, %.2f \n", intercept[bin1], intercept[bin2], intercept[bin3], intercept[bin4], intercept[bin5]);
			}


			// for(int bin3=0; bin3<numBins; bin3++){ // loop over config 3 bins
			// 	double back3 = backgrounds_all[2][bin3];
			// 	if(back3>1.) continue; // no background this large in single bin will ever be viable as part of final answer

			// 	for(int bin4=0; bin4<numBins; bin4++){ // loop over config 4 bins
			// 		double back4 = backgrounds_all[3][bin3];
			// 		if(back4>1.) continue; // no background this large in single bin will ever be viable as part of final answer
					
			// 		for(int bin5=0; bin5<numBins; bin5++){ // loop over config 5 bins
			// 			double back5 = backgrounds_all[4][bin3];
			// 			if(back5>1.) continue; // no background this large in single bin will ever be viable as part of final answer
					
			// 		} // loop over intercepts in config 5
			// 	} // loop over intercepts in config 4
			// } // loop over intercepts in config 3
		} // loop over intercepts in config 2
	} // loop over intercepts in config 1

	printf("I think I have the global S/Sup optimum of %.4f, at intercepts %.2f and %.2f \n",best_S_over_Sup,intercept[best_bins[0]],intercept[best_bins[1]]);

	gStyle->SetOptStat(0);

	TGraph *back_vs_intercept1 = new TGraph(numBins);
	TGraph *back_vs_intercept2 = new TGraph(numBins);
	TGraph *signal_vs_intercept1 = new TGraph(numBins);
	TGraph *signal_vs_intercept2 = new TGraph(numBins);
	for(int bin1=0; bin1<numBins; bin1++){
		back_vs_intercept1->SetPoint(bin1,intercept[bin1],backgrounds_all[0][bin1]);
		back_vs_intercept2->SetPoint(bin1,intercept[bin1],backgrounds_all[1][bin1]);
		signal_vs_intercept1->SetPoint(bin1,intercept[bin1],signal_all[0][bin1]*live_frac[0]);
		signal_vs_intercept2->SetPoint(bin1,intercept[bin1],signal_all[1][bin1]*live_frac[1]*2.);
	}
	TCanvas *c = new TCanvas("","",4*850,2*850);
	c->Divide(4,2);
	c->cd(1);
		back_vs_intercept1->Draw("ALP");
		back_vs_intercept1->SetTitle("Background vs Intercept for Config 1");
		back_vs_intercept1->GetXaxis()->SetTitle("Y-intercept/ cut value");
		back_vs_intercept1->GetYaxis()->SetTitle("Expected Number of Background Events");
	c->cd(2);
		back_vs_intercept2->Draw("ALP");
		back_vs_intercept2->SetTitle("Background vs Intercept for Config 2");
		back_vs_intercept2->GetXaxis()->SetTitle("Y-intercept/ cut value");
		back_vs_intercept2->GetYaxis()->SetTitle("Expected Number of Background Events");
		back_vs_intercept2->GetYaxis()->SetTitleOffset(1.4);		
	c->cd(3);
		signal_vs_intercept1->Draw("ALP");
		signal_vs_intercept1->SetTitle("Signal vs Intercept for Config 1");
		signal_vs_intercept1->GetXaxis()->SetTitle("Y-intercept/ cut value");
		signal_vs_intercept1->GetYaxis()->SetTitle("Number of Weighted Neutrinos Passing All Cuts");
		signal_vs_intercept1->GetYaxis()->SetTitleOffset(1.5);		
	c->cd(4);
		signal_vs_intercept2->Draw("ALP");
		signal_vs_intercept2->SetTitle("Signal vs Intercept for Config 2");
		signal_vs_intercept2->GetXaxis()->SetTitle("Y-intercept/ cut value");
		signal_vs_intercept2->GetYaxis()->SetTitle("Number of Weighted Neutrinos Passing All Cuts");
		signal_vs_intercept2->GetYaxis()->SetTitleOffset(1.5);
	c->cd(5);
		S->Draw("colz");
		S->SetTitle("Summed Signal S from Config 1 and Config 2");
		S->GetXaxis()->SetTitle("Y-intercept Cut in Config 1");
		S->GetYaxis()->SetTitle("Y-intercept Cut in Config 2");
		S->GetYaxis()->SetTitleOffset(1.4);
		S->GetZaxis()->SetTitleOffset(1.6);
		S->GetZaxis()->SetTitle("Summed Weights");
	c->cd(6);
		Sup->Draw("colz");
		Sup->SetTitle("Sup");
		Sup->GetXaxis()->SetTitle("Y-intercept Cut in Config 1");
		Sup->GetYaxis()->SetTitle("Y-intercept Cut in Config 2");
		Sup->GetYaxis()->SetTitleOffset(1.4);
		Sup->GetZaxis()->SetTitleOffset(1.4);
		Sup->GetZaxis()->SetTitle("Sup");
		// gPad->SetLogz();
	c->cd(7);
		S_over_Sup->Draw("colz");
		S_over_Sup->SetTitle("S/Sup");
		S_over_Sup->GetXaxis()->SetTitle("Y-intercept Cut in Config 1");
		S_over_Sup->GetYaxis()->SetTitle("Y-intercept Cut in Config 2");
		S_over_Sup->GetYaxis()->SetTitleOffset(1.4);
		S_over_Sup->GetZaxis()->SetTitleOffset(1.4);
		S_over_Sup->GetZaxis()->SetTitle("S/Sup");
		S_over_Sup->GetZaxis()->SetRangeUser(800,1800);
		gPad->SetLogz();
	char save_title[400];
	sprintf(save_title,"/users/PAS0654/osu0673/A23_analysis_new2/AraRoot/analysis/a23_analysis_programs/diffuse/grid_optimize/plots_v2.png");
	c->SaveAs(save_title);	

}
