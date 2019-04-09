////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	v2_final_plots.cxx 
////	A23 diffuse, make plots of the final cut parameter space
////
////	Nov 2018
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
#include "tools_Cuts.h"

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

	if(argc<6){
		cout<< "Usage\n" << argv[0] << " <isSim> <station> <config> <year_or_energy> <ValForCuts filename>"<<endl;;
		return -1;
	}
	int isSim = atoi(argv[1]);
	int station = atoi(argv[2]);
	int config = atoi(argv[3]);
	double year_or_energy = double(atof(argv[4]));
	
	TH2D *PeakCorr_vs_SNR_all[2];
	PeakCorr_vs_SNR_all[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_all[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *PeakCorr_vs_SNR_cutCal[2];
	PeakCorr_vs_SNR_cutCal[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_cutCal[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft[2];
	PeakCorr_vs_SNR_cutCal_cutSoft[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_cutCal_cutSoft[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[2];
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[2];
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[2];
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[1]=new TH2D("","H",30,0,30,100,0,1);

	TH2D *PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[2];
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[0]=new TH2D("","V",30,0,30,100,0,1);
	PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[1]=new TH2D("","H",30,0,30,100,0,1);

	TH1D *all_events[2];
	TH1D *pass_soft_short_cal[2];
	TH1D *pass_soft_short_cal_wfrms[2];
	TH1D *pass_soft_short_cal_wfrms_box[2];
	TH1D *pass_soft_short_cal_wfrms_box_surf[2];

	TH1D *eff[2];
	TH1D *eff_soft_short_cal[2];
	TH1D *eff_soft_short_cal_wfrms[2];
	TH1D *eff_soft_short_cal_wfrms_box[2];
	TH1D *eff_soft_short_cal_wfrms_box_surf[2];


	for(int i=0; i<2; i++){
		all_events[i] = new TH1D("","",30,0,30);
		pass_soft_short_cal[i] = new TH1D("","",30,0,30);
		pass_soft_short_cal_wfrms[i] = new TH1D("","",30,0,30);
		pass_soft_short_cal_wfrms_box[i] = new TH1D("","",30,0,30);
		pass_soft_short_cal_wfrms_box_surf[i] = new TH1D("","",30,0,30);

		eff_soft_short_cal[i] = new TH1D("","",30,0,30);
		eff_soft_short_cal_wfrms[i] = new TH1D("","",30,0,30);
		eff_soft_short_cal_wfrms_box[i] = new TH1D("","",30,0,30);
		eff_soft_short_cal_wfrms_box_surf[i] = new TH1D("","",30,0,30);
	}

	TH2D *special[2];
	special[0]=new TH2D("","V",90,0,30,500,0,1);
	special[1]=new TH2D("","H",90,0,30,500,0,1);

	TH1D *fracs_power_cut[2];
	fracs_power_cut[0]=new TH1D("","V",100,0,1);
	fracs_power_cut[1]=new TH1D("","H",100,0,1);
	
	int num_total=0;
	int num_in_final_plot=0;
	int num_refilt=0;
	int num_cal=0;
	int num_soft=0;
	int num_short=0;
	int num_box=0;
	int num_surf=0;

	double num_total_trig=0.;
	double num_pass_pol[2]={0.};
	double num_pass_either=0.;

	for(int file_num=5; file_num<argc; file_num++){

		string chRun = "run";
		string file = string(argv[file_num]);
		size_t foundRun = file.find(chRun);
		string strRunNum = file.substr(foundRun+4,4);
		int runNum = atoi(strRunNum.c_str());

		if(isBadRun(station,runNum)) continue;

		TFile *inputFile = TFile::Open(argv[file_num]);
		if(!inputFile){
			cout<<"Can't open val for cuts file!"<<endl;
			return -1;
		}
		cout << "Run " << file_num << " :: " << argv[file_num] << endl;

		TTree *trees[3];
		trees[0] = (TTree*) inputFile->Get("VTree");
		trees[1] = (TTree*) inputFile->Get("HTree");
		trees[2] = (TTree*) inputFile->Get("AllTree");

		double corr_val[2];
		double snr_val[2];
		int WFRMS[2];
		double frac_of_power_notched_V[8];
		double frac_of_power_notched_H[8];
		int Refilt[2];

		trees[0]->SetBranchAddress("corr_val_V",&corr_val[0]);
		trees[0]->SetBranchAddress("snr_val_V",&snr_val[0]);
		trees[0]->SetBranchAddress("wfrms_val_V",&WFRMS[0]);
		trees[0]->SetBranchAddress("Refilt_V",&Refilt[0]);
		trees[1]->SetBranchAddress("corr_val_H",&corr_val[1]);
		trees[1]->SetBranchAddress("snr_val_H",&snr_val[1]);
		trees[1]->SetBranchAddress("wfrms_val_H",&WFRMS[1]);
		trees[1]->SetBranchAddress("Refilt_H",&Refilt[1]);

		int isCal;
		int isSoft;
		int isShort;
		int isCW;
		int isNewBox;
		int isSurf;
		int isBadEvent;
		double weight;

		trees[2]->SetBranchAddress("cal",&isCal);
		trees[2]->SetBranchAddress("soft",&isSoft);
		trees[2]->SetBranchAddress("short",&isShort);
		trees[2]->SetBranchAddress("CW",&isCW);
		trees[2]->SetBranchAddress("box",&isNewBox);
		trees[2]->SetBranchAddress("surf",&isSurf);
		trees[2]->SetBranchAddress("bad",&isBadEvent);
		trees[2]->SetBranchAddress("weight",&weight);

		stringstream ss;
		for(int i=0; i<8; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			trees[0]->SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_V[i]);
		}
		for(int i=8; i<16; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			trees[1]->SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_H[i-8]);
		}
		
		int numEntries = trees[0]->GetEntries();

		//now to loop over events
		for(int event=0; event<trees[0]->GetEntries(); event++){
			trees[0]->GetEvent(event);
			trees[1]->GetEvent(event);
			trees[2]->GetEvent(event);
			if(isBadEvent){
				continue;
			}

			num_total++;
			num_total_trig+=weight;

			if(isCal) num_cal++;
			if(isSoft) num_soft++;
			if(isShort) num_short++;
			if(isNewBox) num_box++;
			if(isSurf) num_surf++;

			if(snr_val[0]>=30.) snr_val[0]=30.;
			if(snr_val[1]>=30.) snr_val[1]=30.;

			if(isSim){
				bool this_pass[2]={false};
				for(int pol=0; pol<2; pol++){
					all_events[pol]->Fill(snr_val[pol],weight);
					if(!isCal && !isSoft && !isShort){
						pass_soft_short_cal[pol]->Fill(snr_val[pol],weight);
						if(!WFRMS[pol]){
							pass_soft_short_cal_wfrms[pol]->Fill(snr_val[pol],weight);
							if(!isNewBox){
								pass_soft_short_cal_wfrms_box[pol]->Fill(snr_val[pol],weight);
								if(!isSurf){
									pass_soft_short_cal_wfrms_box_surf[pol]->Fill(snr_val[pol],weight);
									num_pass_pol[pol]+=weight;
									this_pass[pol]=true;
								}
							}
						}
					}
				}
				if(this_pass[0] || this_pass[1])
					num_pass_either+=weight;
			}


			for(int pol=0; pol<2; pol++){
				PeakCorr_vs_SNR_all[pol]->Fill(snr_val[pol],corr_val[pol],weight);


				
				if(!isCal){ //cut cal pulsers
					PeakCorr_vs_SNR_cutCal[pol]->Fill(snr_val[pol],corr_val[pol],weight);
					
					if(!isSoft){ //cut software triggers 
						PeakCorr_vs_SNR_cutCal_cutSoft[pol]->Fill(snr_val[pol],corr_val[pol],weight);
						
						if(!isShort){ //cut short
							PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->Fill(snr_val[pol],corr_val[pol],weight);
							
							if(!WFRMS[pol]){ //cut WRMS
								PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->Fill(snr_val[pol],corr_val[pol],weight);
								
								if(!isNewBox){ //cut cal box
									PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->Fill(snr_val[pol],corr_val[pol],weight);

									if(!isSurf){

										if(Refilt[pol]){
											num_refilt++;

											vector<double> frac;
											for(int i=0; i<8; i++){
												if(pol==0) frac.push_back(frac_of_power_notched_V[i]);
												else if(pol==1) frac.push_back(frac_of_power_notched_H[i]);
											}
											sort(frac.begin(), frac.end(), std::greater<double>());
											fracs_power_cut[pol]->Fill(frac[2]);
											if(frac[2]<=0.06){ //&& event!=1 && event!=2 && event!=3)
												PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->Fill(snr_val[pol],corr_val[pol],weight);
											}
										} //refiltered?
										else{
											PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->Fill(snr_val[pol],corr_val[pol],weight);
										}
										num_in_final_plot++;
									}
								}
							}
						}
					}
				}
			}
		}
		inputFile->Close();
		delete inputFile;
	}
	cout<<"Num total is "<<num_total<<endl;
	cout<<"Num in final plot "<<num_in_final_plot<<endl;
	cout<<"Num re-filtered is "<<num_refilt<<endl;

	printf("Num cal is %d -- %.2f %\n", num_cal, double(num_cal)/double(num_total)*100.);
	printf("Num box is %d -- %.2f %\n", num_box, double(num_box)/double(num_total)*100.);
	printf("Num soft is %d -- %.2f %\n", num_soft, double(num_soft)/double(num_total)*100.);
	printf("Num short is %d -- %.2f %\n", num_short, double(num_short)/double(num_total)*100.);
	printf("Num surf is %d -- %.2f %\n", num_surf, double(num_surf)/double(num_total)*100.);

	gStyle->SetOptStat(0);
	gStyle->SetStatY(0.9);
	gStyle->SetStatX(0.9);
	gStyle->SetStatW(0.2);
	gStyle->SetStatH(0.2);

	int colors [28] = { kBlue, kRed, kGreen, kMagenta, kCyan};

	if(isSim){
		for(int pol=0; pol<2; pol++){
			for(int bin=0; bin<=all_events[pol]->GetNbinsX(); bin++){
				double thrown = all_events[pol]->GetBinContent(bin);
				double pass_soft_short_cal_this = pass_soft_short_cal[pol]->GetBinContent(bin);
				double pass_soft_short_cal_wfrms_this = pass_soft_short_cal_wfrms[pol]->GetBinContent(bin);
				double pass_soft_short_cal_wfrms_box_this = pass_soft_short_cal_wfrms_box[pol]->GetBinContent(bin);
				double pass_soft_short_cal_wfrms_box_surf_this = pass_soft_short_cal_wfrms_box_surf[pol]->GetBinContent(bin);
				if(thrown>0.){
					eff_soft_short_cal[pol]->SetBinContent(bin,pass_soft_short_cal_this/thrown);
					eff_soft_short_cal_wfrms[pol]->SetBinContent(bin,pass_soft_short_cal_wfrms_this/thrown);
					eff_soft_short_cal_wfrms_box[pol]->SetBinContent(bin,pass_soft_short_cal_wfrms_box_this/thrown);
					eff_soft_short_cal_wfrms_box_surf[pol]->SetBinContent(bin,pass_soft_short_cal_wfrms_box_surf_this/thrown);
				}
				else{
					eff_soft_short_cal[pol]->SetBinContent(bin,0.);
					eff_soft_short_cal_wfrms[pol]->SetBinContent(bin,0.);
					eff_soft_short_cal_wfrms_box[pol]->SetBinContent(bin,0.);
					eff_soft_short_cal_wfrms_box_surf[pol]->SetBinContent(bin,0.);
				}
			}
		}

		TCanvas *c2 = new TCanvas("","",1100,850);
		c2->Divide(2,2);
		for(int pol=0; pol<2; pol++){

			if(pol==0){
				all_events[pol]->SetTitle("Events Passing Cuts VPol");
				eff_soft_short_cal[pol]->SetTitle("Efficiency VPol");
			}

			if(pol==1){
				all_events[pol]->SetTitle("Events Passing Cuts HPol");
				eff_soft_short_cal[pol]->SetTitle("Efficiency HPol");
			}
			all_events[pol]->GetXaxis()->SetTitle("3rd Highest Vpeak/RMS");
			all_events[pol]->GetYaxis()->SetTitle("Events (weighted)");
			eff_soft_short_cal[pol]->GetXaxis()->SetTitle("3rd Highest Vpeak/RMS");
			eff_soft_short_cal[pol]->GetYaxis()->SetTitle("Events (weighted)");
			
			c2->cd(pol+1);
				all_events[pol]->Draw("");
				all_events[pol]->GetYaxis()->SetRangeUser(1,3e3);
				gPad->SetLogy();
				pass_soft_short_cal[pol]->Draw("same");
				pass_soft_short_cal_wfrms[pol]->Draw("same");
				pass_soft_short_cal_wfrms_box[pol]->Draw("same");
				pass_soft_short_cal_wfrms_box_surf[pol]->Draw("same");
				pass_soft_short_cal[pol]->SetLineColor(colors[0]);
				pass_soft_short_cal_wfrms[pol]->SetLineColor(colors[1]);
				pass_soft_short_cal_wfrms_box[pol]->SetLineColor(colors[2]);
				pass_soft_short_cal_wfrms_box_surf[pol]->SetLineColor(colors[3]);
				pass_soft_short_cal[pol]->SetLineWidth(2.);
				pass_soft_short_cal_wfrms[pol]->SetLineWidth(2.);
				pass_soft_short_cal_wfrms_box[pol]->SetLineWidth(2.);
				pass_soft_short_cal_wfrms_box_surf[pol]->SetLineWidth(2.);
				
				if(pol+1==1){
					TLegend *leg = new TLegend(0.48,0.6,0.9,0.9);
					leg->AddEntry(all_events[pol],"All Events","l");
					leg->AddEntry(pass_soft_short_cal[pol],"Cut Soft, Short, and Tagged Cal","l");
					leg->AddEntry(pass_soft_short_cal_wfrms[pol],"+Cut WFMRS","l");
					leg->AddEntry(pass_soft_short_cal_wfrms_box[pol],"+Cut Cal Pulser Reco","l");
					leg->AddEntry(pass_soft_short_cal_wfrms_box_surf[pol],"+Cut Surface","l");
					leg->Draw();
				}

			c2->cd(pol+3);
				eff_soft_short_cal[pol]->Draw("");
				eff_soft_short_cal_wfrms[pol]->Draw("same");
				eff_soft_short_cal_wfrms_box[pol]->Draw("same");
				eff_soft_short_cal_wfrms_box_surf[pol]->Draw("same");
				eff_soft_short_cal[pol]->SetLineColor(colors[0]);
				eff_soft_short_cal_wfrms[pol]->SetLineColor(colors[1]);
				eff_soft_short_cal_wfrms_box[pol]->SetLineColor(colors[2]);
				eff_soft_short_cal_wfrms_box_surf[pol]->SetLineColor(colors[3]);
				eff_soft_short_cal[pol]->SetLineWidth(2.);
				eff_soft_short_cal_wfrms[pol]->SetLineWidth(2.);
				eff_soft_short_cal_wfrms_box[pol]->SetLineWidth(2.);
				eff_soft_short_cal_wfrms_box_surf[pol]->SetLineWidth(2.);
		}

		char efficiency_title[400];
		sprintf(efficiency_title,
				 "%s/%d.%d.%d_A%d_c%d_E%2.1f_%dEvents_Efficiency.png",plotPath,year_now, month_now, day_now,station,config,year_or_energy,num_total);
		c2->SaveAs(efficiency_title);
		delete c2;
	}

	if(isSim){
		for(int pol=0; pol<2; pol++){
			printf("Eff on pol %d is %.3f/%.3f = %.3f \n", pol, num_pass_pol[pol], num_total_trig,num_pass_pol[pol]/num_total_trig);
		}
	}
	printf("Eff on both pols combined is %.3f / %.3f = %.3f \n", num_pass_either, num_total_trig, num_pass_either/num_total_trig);

	//save out SNR vs WavefrontRMS plot
	char graph_title[2][300];
	char title[300];

	int cal=0;
	int soft=0;
	int Short=0;
	int wrms=0;
	int box = 0;
	int surf=0;
	int cw=0;

	//save out the Corr vs SNR plot for all 
	sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	TCanvas *c2 = new TCanvas("","",2.1*850,850);
	c2->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		c2->cd(pol+1);
		PeakCorr_vs_SNR_all[pol]->Draw("colz");

		PeakCorr_vs_SNR_all[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
		PeakCorr_vs_SNR_all[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
		PeakCorr_vs_SNR_all[pol]->SetTitle(graph_title[pol]);
		gPad->SetLogz();
	}
	if(isSim)
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_E%2.1f_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,year_or_energy, num_total,cal,soft,Short,wrms,box,surf);
	else
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,num_total,cal,soft,Short,wrms,box,surf);
	c2->SaveAs(title);
	delete c2;
	delete PeakCorr_vs_SNR_all[0]; delete PeakCorr_vs_SNR_all[1];

	//turn on cal
	cal=1;
	sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	TCanvas *c3 = new TCanvas("","",2.1*850,850);
	c3->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		c3->cd(pol+1);
		PeakCorr_vs_SNR_cutCal[pol]->Draw("colz");
		PeakCorr_vs_SNR_cutCal[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
		PeakCorr_vs_SNR_cutCal[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
		PeakCorr_vs_SNR_cutCal[pol]->SetTitle(graph_title[pol]);
		gPad->SetLogz();
	}
	if(isSim)
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_E%2.1f_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,year_or_energy, num_total,cal,soft,Short,wrms,box,surf);
	else
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,num_total,cal,soft,Short,wrms,box,surf);
	c3->SaveAs(title);
	delete c3;
	delete PeakCorr_vs_SNR_cutCal[0]; delete PeakCorr_vs_SNR_cutCal[1];

	//turn on cal, soft
	soft=1;
	sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	TCanvas *c4 = new TCanvas("","",2.1*850,850);
	c4->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		c4->cd(pol+1);
		PeakCorr_vs_SNR_cutCal_cutSoft[pol]->Draw("colz");
		PeakCorr_vs_SNR_cutCal_cutSoft[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
		PeakCorr_vs_SNR_cutCal_cutSoft[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
		PeakCorr_vs_SNR_cutCal_cutSoft[pol]->SetTitle(graph_title[pol]);
		gPad->SetLogz();
	}
	if(isSim)
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_E%2.1f_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,year_or_energy, num_total,cal,soft,Short,wrms,box,surf);
	else
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,num_total,cal,soft,Short,wrms,box,surf);
	c4->SaveAs(title);
	delete c4;
	delete PeakCorr_vs_SNR_cutCal_cutSoft[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft[1];

	//turn on cal, soft, short
	Short=1;
	sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	TCanvas *c5 = new TCanvas("","",2.1*850,850);
	c5->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		c5->cd(pol+1);
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->Draw("colz");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[pol]->SetTitle(graph_title[pol]);
		gPad->SetLogz();
	}
	if(isSim)
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_E%2.1f_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,year_or_energy, num_total,cal,soft,Short,wrms,box,surf);
	else
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,num_total,cal,soft,Short,wrms,box,surf);
	c5->SaveAs(title);
	delete c5;
	delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort[1];

	//turn on cal, soft, short, wmrs
	wrms=1;
	sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	TCanvas *c6 = new TCanvas("","",2.1*850,850);
	c6->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		c6->cd(pol+1);
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->Draw("colz");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[pol]->SetTitle(graph_title[pol]);
		gPad->SetLogz();
	}
	if(isSim)
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_E%2.1f_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,year_or_energy, num_total,cal,soft,Short,wrms,box,surf);
	else
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,num_total,cal,soft,Short,wrms,box,surf);
	c6->SaveAs(title);
	delete c6;
	delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS[1];

	//turn on cal, soft, short, wmrs, box
	box=1;
	sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	TCanvas *c7 = new TCanvas("","",2.1*850,850);
	c7->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		c7->cd(pol+1);
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->Draw("colz");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[pol]->SetTitle(graph_title[pol]);
		gPad->SetLogz();
	}
	if(isSim)
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_E%2.1f_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,year_or_energy, num_total,cal,soft,Short,wrms,box,surf);
	else
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,num_total,cal,soft,Short,wrms,box,surf);
	c7->SaveAs(title);
	delete c7;
	delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox[1];

	surf=1;
	box=1;
	sprintf(graph_title[0],"VPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	sprintf(graph_title[1],"HPol: cal %d, soft %d ,short %d , wrms  %d , box %d , surf %d, cw %d",cal,soft,Short,wrms,box,surf,cw);
	TCanvas *c8 = new TCanvas("","",2.1*850,850);
	c8->Divide(2,1);
	for(int pol=0; pol<2; pol++){
		c8->cd(pol+1);
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->Draw("colz");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->GetYaxis()->SetTitle("Peak Correlation Value");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->GetXaxis()->SetTitle("3rd Highest VPeak/RMS");
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->SetTitle(graph_title[pol]);
		gPad->SetLogz();
		PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->GetZaxis()->SetRangeUser(1.,5e4);		
		// PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->GetXaxis()->SetRangeUser(0,10);
		// PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[pol]->GetYaxis()->SetRangeUser(0,0.5);
	}
	if(isSim)
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_E%2.1f_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,year_or_energy, num_total,cal,soft,Short,wrms,box,surf);
	else
		sprintf(title, "%s/%d.%d.%d_A%d_c%d_%dEvents_Correlation_vs_SNR_cal%dF_soft%d_short%d_wrms%d_newbox%d_surf%d.png",plotPath,year_now, month_now, day_now,station,config,num_total,cal,soft,Short,wrms,box,surf);
	c8->SaveAs(title);
	delete c8;
	delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[0]; delete PeakCorr_vs_SNR_cutCal_cutSoft_cutShort_cutWRMS_cutBox_cutSurf[1];

}