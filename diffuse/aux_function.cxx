void Optimize(int station, int config, int pol_select, double slope, double &background_out, double &signal_out, double &s_over_sup_out, bool doPrint){

	double select_slope=slope;


	// let's see if we can do it with arrays....
	int numSNRbins=200;
	double numEventsPassed[numSNRbins];
	double numEventsPassed_diff[numSNRbins];

	double dRcut=0.1; // SNR bin incremets (bins in 0.1 SNR units)
	double Rcutmin=0.;
	double intercept[numSNRbins];
	for(int bin=0; bin<numSNRbins; bin++){
		numEventsPassed[bin]=0.;
		numEventsPassed_diff[bin]=0.;
		intercept[bin] = Rcutmin + double(bin)*dRcut;
	}
	// for(int bin=0; bin<numSNRbins; bin++){
	// 	printf("SNR Bin %d has intercept value %.2f \n", bin,intercept[bin]);
	// }
	
	// need to be able to make the final 2D distribution

	/////////////////////////////////////////////////
	/////////////////////////////////////////////////
	/// Load all of the data in
	/////////////////////////////////////////////////
	/////////////////////////////////////////////////

	TChain dataVTree("VTree");
	TChain dataHTree("HTree");
	TChain dataAllTree("AllTree");
	char the_data[500];
	sprintf(the_data,"/fs/scratch/PAS0654/ara/10pct/ValsForCuts/A%d/c%d/cutvals_drop_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_*.root",station,config);
	dataVTree.Add(the_data);
	dataHTree.Add(the_data);
	dataAllTree.Add(the_data);
	int numDataEvents = dataVTree.GetEntries();
	printf("Num of data entries is %d \n", numDataEvents);

	double max=0.05;
	TH2D *h2SNRvsCorr_data[2]; // SNR on Y axis, Corr on X axis, like in the TB
	h2SNRvsCorr_data[0]=new TH2D("","V Data",100,0,max,300,0,30);
	h2SNRvsCorr_data[1]=new TH2D("","H Data",100,0,max,300,0,30);

	// do this inside brackets for scoping power and re-use of identical variable names when it comes time for simulation to happen
	{

		double corr_val[2];
		double snr_val[2];
		int WFRMS[2];
		double frac_of_power_notched_V[8];
		double frac_of_power_notched_H[8];
		int Refilt[2];

		dataVTree.SetBranchAddress("corr_val_V",&corr_val[0]);
		dataVTree.SetBranchAddress("snr_val_V",&snr_val[0]);
		dataVTree.SetBranchAddress("wfrms_val_V",&WFRMS[0]);
		dataVTree.SetBranchAddress("Refilt_V",&Refilt[0]);
		dataHTree.SetBranchAddress("corr_val_H",&corr_val[1]);
		dataHTree.SetBranchAddress("snr_val_H",&snr_val[1]);
		dataHTree.SetBranchAddress("wfrms_val_H",&WFRMS[1]);
		dataHTree.SetBranchAddress("Refilt_H",&Refilt[1]);

		int isCal;
		int isSoft;
		int isShort;
		int isCW;
		int isNewBox;
		int isSurf[2];
		int isBadEvent;
		double weight;
		int isSurfEvent_top[2];
		int unixTime;
		int isFirstFiveEvent;
		int hasBadSpareChanIssue;
		int runNum;
		int badRun;

		dataAllTree.SetBranchAddress("cal",&isCal);
		dataAllTree.SetBranchAddress("soft",&isSoft);
		dataAllTree.SetBranchAddress("short",&isShort);
		dataAllTree.SetBranchAddress("CW",&isCW);
		dataAllTree.SetBranchAddress("box",&isNewBox);
		dataAllTree.SetBranchAddress("surf_V",&isSurf[0]);
		dataAllTree.SetBranchAddress("surf_H",&isSurf[1]);
		dataAllTree.SetBranchAddress("bad",&isBadEvent);
		dataAllTree.SetBranchAddress("weight",&weight);
		dataAllTree.SetBranchAddress("surf_top_V",&isSurfEvent_top[0]);
		dataAllTree.SetBranchAddress("surf_top_H",&isSurfEvent_top[1]);
		dataAllTree.SetBranchAddress("unixTime",&unixTime);
		dataAllTree.SetBranchAddress("isFirstFiveEvent",&isFirstFiveEvent);
		dataAllTree.SetBranchAddress("hasBadSpareChanIssue",&hasBadSpareChanIssue);
		dataAllTree.SetBranchAddress("runNum",&runNum);
		dataAllTree.SetBranchAddress("badRun",&badRun);

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

		for(int event=0; event<numDataEvents; event++){
			dataVTree.GetEvent(event);
			dataHTree.GetEvent(event);
			dataAllTree.GetEvent(event);
			if( isSoft || isBadEvent || hasBadSpareChanIssue || isFirstFiveEvent || isShort || isCal || badRun){
				continue;
			}
			if(isBadLivetime(station,unixTime)){
				continue;
			}
			for(int pol=0; pol<2; pol++){
				if(pol!=pol_select){
					continue;
				}
				if(!WFRMS[pol] && !isNewBox && !isSurf[0] && !isSurf[1] && !isSurfEvent_top[pol]){
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
						h2SNRvsCorr[pol]->Fill(corr_val[pol],snr_val[pol],weight);
						for(int bin=0; bin<numSNRbins; bin++){ //check all potential intercepts
							double this_y_val = (slope * corr_val[pol] ) + intercept[bin]; // compute the SNR to pass at this intercept
							// printf("For bin %d, with intercept %.2f, SNR to pass is %.2f \n", bin, intercept[bin], this_y_val);
							if(snr_val[pol]>=this_y_val){ // does this event pass?
								// printf("     This event has SNR %.2f, so it passes!\n",snr_val[pol]);
								// printf("          Old number of events in this bin %5f \n", numEventsPassed[bin]);
								numEventsPassed[bin]+=1.;
								// printf(".         New number of events in this bin %5f \n", numEventsPassed[bin]);
							} // does event pass Rcut
						} // loop over SNR cuts
					}// not failing CW power cut?
				}// passes rest of analysis (not WFRMS, box, surface)
			}// loop over polarizations
		}// loop over events
	}


	/////////////////////////////////////////////////
	/////////////////////////////////////////////////
	/// Now for differential distribution etc.
	/////////////////////////////////////////////////
	/////////////////////////////////////////////////


	// now to get the differential distribution up and running
	for(int bin=0; bin<numSNRbins-1; bin++){
		numEventsPassed_diff[bin] = numEventsPassed[bin] - numEventsPassed[bin+1];
		// printf("Bin %d at cut %.2f has %5.f events passing, and next bin has %5f, so diff is %.5f \n", bin, intercept[bin],numEventsPassed[bin],numEventsPassed[bin+1],numEventsPassed_diff[bin]);
	}

	TH1D *hEventsVsSNR = new TH1D("","",numSNRbins,0,20);
	for(int bin=0; bin<numSNRbins; bin++){
		// printf("Bin %d I'm in SNR bin %.2f in the histo and I'm going to fill with %.2f \n", bin, hEventsVsSNR->GetBinCenter(bin+1),numEventsPassed_diff[bin]);
		hEventsVsSNR->SetBinContent(bin+1,numEventsPassed_diff[bin]);
	}
	int max_bin = hEventsVsSNR->GetMaximumBin();
	int fit_start_bin = max_bin+14;
	double start_of_fit = hEventsVsSNR->GetBinCenter(fit_start_bin);
	int last_filled_bin = hEventsVsSNR->FindLastBinAbove(0.,1);
	last_filled_bin+=5; // go up two more bins just to make sure the fit is over
	double end_of_fit = hEventsVsSNR->GetBinCenter(last_filled_bin);
	// printf("Start of fit is %.2f and end of fit is %.2f \n", start_of_fit, end_of_fit);

	// now we exponential fit
	char equation[150];
	sprintf(equation,"exp([0]*x+[1])");
	// sprintf(equation,"gaus");
	TF1 *fit = new TF1("ExpFit",equation,start_of_fit,end_of_fit);
	int status = hEventsVsSNR->Fit("ExpFit","LL,R");
	if(status!=0){
		printf("Something went wrong with the fit! Quitting...\n");
		return 0;
	}
	double fitParams[2];
	fitParams[0] = fit->GetParameter(0);
	fitParams[1] = fit->GetParameter(1);
	printf("Fit Parameters are %.2f and %.2f \n", fitParams[0], fitParams[1]);

	double binWidthIntercept = hEventsVsSNR->GetBinWidth(1);
	double leftBoundary = start_of_fit - binWidthIntercept/2.;
	double rightBoundary = end_of_fit + binWidthIntercept/2.;
	int numBinsThisFit = (rightBoundary - leftBoundary)/binWidthIntercept + 1; // how many bins do we need in our histogram to actually do the fitting
	printf("Number of bins in this fit is %d \n", numBinsThisFit);
	char this_fit_title[400];
	sprintf(this_fit_title,"Fit_slope_%.2f ", slope);
	TH1D *hNumObserved = new TH1D(this_fit_title,"",numBinsThisFit,leftBoundary,rightBoundary);
	for(int bin=0; bin<numBinsThisFit; bin++){
		double originalContent = hEventsVsSNR->GetBinContent(bin+fit_start_bin);
		hNumObserved->SetBinContent(bin+1,originalContent);
	}

	/*
		Now we must prepare our *expectation* for the number of events
		in a bin given that we now have our exponential model
	*/

	double numExpected[numBinsThisFit];
	for(int bin=0; bin<numBinsThisFit; bin++){
		double modelPrediction = exp(fitParams[0]*(hNumObserved->GetBinCenter(bin+1)) + fitParams[1]);
		numExpected[bin] = modelPrediction;
	}

	/*
		Now we compute the log-likelihood by hand
	*/

	double logL=0.;
	double numObservedTotal=0.;
	double numExpectedTotal_sum=0;
	for(int bin=0; bin<numBinsThisFit; bin++){
		double thisObserved = hNumObserved->GetBinContent(bin+1);
		double thisExpected = numExpected[bin];
		printf("At bin %d Observed %.2f and Expected %.2f \n", bin, thisObserved, thisExpected );
		numObservedTotal+=thisObserved;
		numExpectedTotal_sum+=thisExpected;
		logL += ReturnLogL_highN( thisObserved,thisExpected );
	}
	printf("The logL is %.3f \n", logL);
	
	// double numExpectedTotal_Integral = (1./(fitParams[0])) * ( exp(fitParams[0]*end_of_fit + fitParams[1]) - exp(fitParams[0]*start_of_fit + fitParams[1]));
	// double numExpectedTotal_Integral = (1./(fitParams[0]*binWidthIntercept)) * ( exp(fitParams[0]*end_of_fit + fitParams[1]) - exp(fitParams[0]*start_of_fit + fitParams[1]));
	double numExpectedTotal_Integral = (1./(fitParams[0]*binWidthIntercept)) * ( exp(fitParams[0]*rightBoundary + fitParams[1]) - exp(fitParams[0]*leftBoundary + fitParams[1]));
	printf("Best fit sum bins: %.2f. Best fit do integral: %.2f. Num observed: %.2f. \n",numExpectedTotal_sum, numExpectedTotal_Integral, numObservedTotal);

	TRandom3 *test_random = new TRandom3();

	/*
		Now for toy simulations
	*/
	sprintf(this_fit_title,"fCopyFit", slope);
	TF1 *fitCopy = new TF1(this_fit_title, "exp([0]*x+[1])", start_of_fit, end_of_fit);
	fitCopy->SetParameters(fitParams[0], fitParams[1]);
	double less_BestFit_logL = 0.; // values lower than BestFit_logL
	double Total_Toy_logL = 0.; // total logL values from Toy

	int num_Toy = 10000;
	int Toy_logL_bin = 150;
	double min_Toy_logL = 0.;
	double max_Toy_logL = 150.;
	TH1D *hToy_logL = new TH1D("hToy_logL", "", Toy_logL_bin, min_Toy_logL, max_Toy_logL );
	char test_title[400];
	for(int num=0; num<num_Toy; num++){
		// fill this toy pseudo experiment with a poisson fluctuations of the events observed
		sprintf(test_title, "hToy %d ", num);
		TH1D *hToy = new TH1D(test_title,"",numBinsThisFit,leftBoundary,rightBoundary);
		int Evts_Poisson = test_random->Poisson(numExpectedTotal_Integral);
		hToy->FillRandom("fCopyFit",Evts_Poisson);
		double logL_log_Toy=0.; // get logL for this toy
		for(int bin=0; bin<numBinsThisFit; bin++){
			logL_log_Toy+=ReturnLogL_highN(hToy->GetBinContent(bin+1), numExpected[bin]);
		}
		// TCanvas *cToyHist = new TCanvas ("cToyHist","", 800, 600);
		// cToyHist->cd();
		// 	cToyHist->cd()->SetLogy();
		// 	sprintf( test_title, "Toy hist, evts : %d, -2Ln(L): %.2f", Evts_Poisson, logL_log_Toy );
		// 	hToy->SetTitle(test_title);
		// 	hToy->Draw();
		// 	fit->Draw("same");
		// char this_save_title[400];
		// sprintf(this_save_title,"toy%d.png",num);
		// cToyHist->SaveAs(this_save_title);
		// delete cToyHist;
		delete hToy;
		hToy_logL->Fill(logL_log_Toy);
		if ( logL_log_Toy <= logL ){
			less_BestFit_logL += logL_log_Toy;
		}
		Total_Toy_logL += logL_log_Toy;
	}

	// vertical line for log likelihood
	double P_value = less_BestFit_logL / Total_Toy_logL;
	double DataLogL_x[2] = { logL, logL };
	double DataLogL_y[2] = { 0, 5 };
	TGraph *gDataLogL = new TGraph (2, DataLogL_x, DataLogL_y);

	printf(BLUE"About to do background estimation\n"RESET);
	/*
		Now, we must find our background estimate
		Which is where we integrate the exponential above our cut value
		And then use that to get s_up
	*/

	double S_up_array[numSNRbins];
	double S_array[numSNRbins];
	double back_estimate[numSNRbins];
	for(int bin=0; bin<numSNRbins; bin++){
		S_up_array[bin]=0.;
		S_array[bin]=0.;
		back_estimate[bin]=0.;
	}
	int startBin = 80;
	for(int bin=startBin; bin<numSNRbins; bin++){
		double cut = intercept[bin];
		double this_back_estimate = (1./(fitParams[0]*binWidthIntercept)) * (-exp(fitParams[0]*cut + fitParams[1]));
		this_back_estimate*=10.; // make it 10 times bigger, for switch to 100% sample
		back_estimate[bin]=this_back_estimate;
		double achieved_alpha;
		// double s_up = GetS_up_TMath(back_estimate,achieved_alpha, 0.9); // compute S_up for this background
		double s_up = GetS_up(back_estimate[bin],achieved_alpha, 0.9); // compute S_up for this background
		S_up_array[bin] = s_up;
		printf("For cut %.2f background estimate is %.3f with sup %.2f \n",cut,back_estimate[bin],S_up_array[bin]);
	}

	/////////////////////////////////////////////////
	/////////////////////////////////////////////////
	/// Now, S over S_up
	/////////////////////////////////////////////////
	/////////////////////////////////////////////////

	/*
		Now loop over simulation to get the S in the S_up
	*/

	TChain simVTree("VTree");
	TChain simHTree("HTree");
	TChain simAllTree("AllTree");
	char the_sims[500];
	sprintf(the_sims,"/fs/scratch/PAS0654/ara/sim/ValsForCuts/A%d/c%d/E%d/cutvals_drop_snrbins_0_0_wfrmsvals_-1.3_-1.4_run_*.root",station,config,int(year_or_energy));
	simVTree.Add(the_sims);
	simHTree.Add(the_sims);
	simAllTree.Add(the_sims);
	int numSimEvents = simVTree.GetEntries();
	printf("Num of sim entries is %d \n", numSimEvents);

	TH2D *h2SNRvsCorr_sim[2]; // SNR on Y axis, Corr on X axis, like in the TB
	h2SNRvsCorr_sim[0]=new TH2D("","V Sim",100,0,max,30,0,30);
	h2SNRvsCorr_sim[1]=new TH2D("","H Sim",100,0,max,30,0,30);
	
	// and now get values out
	{
		double corr_val[2];
		double snr_val[2];
		int WFRMS[2];
		double frac_of_power_notched_V[8];
		double frac_of_power_notched_H[8];
		int Refilt[2];

		simVTree.SetBranchAddress("corr_val_V",&corr_val[0]);
		simVTree.SetBranchAddress("snr_val_V",&snr_val[0]);
		simVTree.SetBranchAddress("wfrms_val_V",&WFRMS[0]);
		simVTree.SetBranchAddress("Refilt_V",&Refilt[0]);
		simHTree.SetBranchAddress("corr_val_H",&corr_val[1]);
		simHTree.SetBranchAddress("snr_val_H",&snr_val[1]);
		simHTree.SetBranchAddress("wfrms_val_H",&WFRMS[1]);
		simHTree.SetBranchAddress("Refilt_H",&Refilt[1]);

		int isCal;
		int isSoft;
		int isShort;
		int isCW;
		int isNewBox;
		int isSurf[2];
		int isBadEvent;
		double weight;
		int isSurfEvent_top[2];
		int unixTime;
		int isFirstFiveEvent;
		int hasBadSpareChanIssue;

		simAllTree.SetBranchAddress("cal",&isCal);
		simAllTree.SetBranchAddress("soft",&isSoft);
		simAllTree.SetBranchAddress("short",&isShort);
		simAllTree.SetBranchAddress("CW",&isCW);
		simAllTree.SetBranchAddress("box",&isNewBox);
		simAllTree.SetBranchAddress("surf_V",&isSurf[0]);
		simAllTree.SetBranchAddress("surf_H",&isSurf[1]);
		simAllTree.SetBranchAddress("bad",&isBadEvent);
		simAllTree.SetBranchAddress("weight",&weight);
		simAllTree.SetBranchAddress("surf_top_V",&isSurfEvent_top[0]);
		simAllTree.SetBranchAddress("surf_top_H",&isSurfEvent_top[1]);
		simAllTree.SetBranchAddress("unixTime",&unixTime);
		simAllTree.SetBranchAddress("isFirstFiveEvent",&isFirstFiveEvent);
		simAllTree.SetBranchAddress("hasBadSpareChanIssue",&hasBadSpareChanIssue);

		stringstream ss;
		for(int i=0; i<8; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			simVTree.SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_V[i]);
		}
		for(int i=8; i<16; i++){
			ss.str("");
			ss<<"PowerNotch_Chan"<<i;
			simHTree.SetBranchAddress(ss.str().c_str(),&frac_of_power_notched_H[i-8]);
		}

		for(int event=0; event<numSimEvents; event++){

			simVTree.GetEvent(event);
			simHTree.GetEvent(event);
			simAllTree.GetEvent(event);
			for(int pol=0; pol<2; pol++){

				h2SNRvsCorr_sim[pol]->Fill(corr_val[pol], snr_val[pol],weight);
				if(pol!=pol_select){
					continue;
				}

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

				double this_SNR=snr_val[pol];
				double this_corr=corr_val[pol];

				if(!WFRMS[pol] && !failsCWPowerCut){
					if(!isNewBox){
						if(!isSurf[0] && !isSurf[1] && !isSurfEvent_top[pol]){
							// loop over every bin (intercept value), and figure out if this event would have passed or not
							for(int bin=startBin; bin<numSNRbins; bin++){
								double failsRcut=false;
								double thisIntercept = intercept[bin];
								double this_y_val = this_corr*select_slope + thisIntercept;
								if(this_SNR>=this_y_val){
									S_array[bin]+=weight;
								} // passes RCut
							} // loop over SNR cut bins
						} // passes surface
					} // passes new box
				} // passes WFRMS and CW power cut
			} //loop over polarizations
		} // loop over events
	} // scoping
		
	double SoverSup[numSNRbins];
	for(int bin=0; bin<numSNRbins; bin++){
		double this_S = S_array[bin];
		double this_Sup = S_up_array[bin];
		double this_SoverSup;
		if(!this_Sup>0)
			this_SoverSup=0.;
		else
			this_SoverSup = this_S/this_Sup;

		SoverSup[bin] = this_SoverSup;
		printf("For bin %d, intercept %.2f, S is %.2f and S_up is %.2f for S/S_up of %.2f  \n", bin,intercept[bin],this_S, this_Sup, this_SoverSup);
	}
	double max_SoverSup=0.;
	double optimal_intercept=0.;
	double optimal_background_estimate=0.;
	double max_signal=0.;
	for(int bin=0; bin<numSNRbins; bin++){
		double this_intercept = intercept[bin];
		if(this_intercept<8. || this_intercept>14.){
			continue;
		}
		else{
			double this_SoverSup=SoverSup[bin];
			if(this_SoverSup>max_SoverSup){
				// printf("At bin %d, for intercept %.2f, new S/Sup of %.2f is greater than current %.2f \n",bin, intercept[bin], this_SoverSup,max_SoverSup);
				optimal_intercept = this_intercept;
				max_SoverSup=this_SoverSup;
				optimal_background_estimate = back_estimate[bin];
				max_signal = S_array[bin];
			}
		}
	}
	background_out = optimal_background_estimate;
	sinal_out = max_signal;
	s_over_sup_out = max_SoverSup;
	TGraph *gSoverSup = new TGraph(numSNRbins,intercept,SoverSup);

	printf("Found optimal intercept at %.2f with %.4f background evnets \n", optimal_intercept, optimal_background_estimate);
	double select_inter=optimal_intercept;

	if(doPrint){
		vector <double> x_vals_for_line;
		vector <double> y_vals_for_line;
		for(double x=0; x<0.020; x+=0.00001){
			double y_val = (slope * x ) + optimal_intercept;
			x_vals_for_line.push_back(x);
			y_vals_for_line.push_back(y_val);
		}
		TGraph *cut_line = new TGraph(x_vals_for_line.size(), &x_vals_for_line[0], &y_vals_for_line[0]);
		
		double optimal_intercept_line_x[2] = { optimal_intercept, optimal_intercept };
		double optimal_intercept_line_y[2] = { 0, 3e3 };
		TGraph *optimal_intercept_line = new TGraph (2, optimal_intercept_line_x, optimal_intercept_line_y);

		TCanvas *cRcut = new TCanvas("","",4*850,2*850);
		cRcut->Divide(4,2);
		cRcut->cd(1);
			h2SNRvsCorr[pol_select]->Draw("colz");
			h2SNRvsCorr[pol_select]->GetXaxis()->SetTitle("Correlation Value");
			h2SNRvsCorr[pol_select]->GetYaxis()->SetTitle("SNR");
			gPad->SetLogz();
			cut_line->Draw("same");
			cut_line->SetLineColor(kRed);
		cRcut->cd(2);
			hEventsVsSNR->Draw("");
			hEventsVsSNR->GetXaxis()->SetTitle("SNR Cut (y-intercept value)");
			hEventsVsSNR->GetYaxis()->SetTitle("Number of Events Cut");
			hEventsVsSNR->SetTitle("Differential Distribution");
			gPad->SetLogy();
			// hEventsVsSNR->GetXaxis()->SetRangeUser(8.6,10.);
			// hEventsVsSNR->GetYaxis()->SetRangeUser(8e1,1e4);
		cRcut->cd(3);
			hNumObserved->Draw("HIST");
			hNumObserved->GetXaxis()->SetTitle("SNR Cut (y-intercept value)");
			hNumObserved->GetYaxis()->SetTitle("Number of Events Cut");
			char fit_title_words[400];
			sprintf(fit_title_words,"Fit exp(%.2fx + %.2f)",fitParams[0], fitParams[1]);
			hNumObserved->SetTitle(fit_title_words);
			// hNumObserved->GetXaxis()->SetRangeUser(8.6,10.);
			// hNumObserved->GetYaxis()->SetRangeUser(8e1,1e4);
			fit->Draw("same");
			gPad->SetLogy();
		cRcut->cd(4);
			sprintf( test_title, "data logL: %.2f, P-value: %f", logL, P_value );
			hToy_logL->SetTitle(test_title);
			hToy_logL->Draw();
			hToy_logL->GetYaxis()->SetTitle("Number of Pseudo Experiments");
			hToy_logL->GetXaxis()->SetTitle("-2log(L)");
			gDataLogL->SetLineColor(kRed);
			gDataLogL->Draw("l");
			gPad->SetLogy();
		cRcut->cd(5);
			h2SNRvsCorr_sim[pol_select]->Draw("colz");
			h2SNRvsCorr_sim[pol_select]->GetXaxis()->SetTitle("Correlation Value");
			h2SNRvsCorr_sim[pol_select]->GetYaxis()->SetTitle("SNR");
			gPad->SetLogz();
			cut_line->Draw("same");
			cut_line->SetLineColor(kRed);
		cRcut->cd(6);
			gSoverSup->Draw("ALP");
			gSoverSup->GetXaxis()->SetRangeUser(8.,13.);
			gSoverSup->GetYaxis()->SetRangeUser(0.,3e3);
			gSoverSup->GetYaxis()->SetTitle("S/S_up");
			gSoverSup->GetYaxis()->SetTitleOffset(1.8);
			gSoverSup->GetXaxis()->SetTitle("SNR Cut (y-intercept value)");
			optimal_intercept_line->Draw("same");
			optimal_intercept_line->SetLineColor(kRed);
			char sup_plot_words[400];
			sprintf(sup_plot_words,"Peak S/S_up %.2f ",optimal_intercept);
			gSoverSup->SetTitle(sup_plot_words);
		// cRcut->cd(7);
		// 	eff_soft_short_cal_wfrms[pol_select]->Draw("");
		// 	eff_soft_short_cal_wfrms_box[pol_select]->Draw("same");
		// 	eff_soft_short_cal_wfrms_box_surf[pol_select]->Draw("same");
		// 	eff_soft_short_cal_wfrms_box_surf_rcut[pol_select]->Draw("same");
		// 		TLegend *leg = new TLegend(0.5,0.4,0.9,0.2);
		// 		leg->AddEntry(eff_soft_short_cal_wfrms[pol_select],"Cut WFMRS","l");
		// 		leg->AddEntry(eff_soft_short_cal_wfrms_box[pol_select],"+Cut Cal Pulser Reco","l");
		// 		leg->AddEntry(eff_soft_short_cal_wfrms_box_surf[pol_select],"+Cut Surface","l");
		// 		leg->AddEntry(eff_soft_short_cal_wfrms_box_surf_rcut[pol_select],"+Cut Peak/Corr","l");
		// 		leg->Draw();
		char save_title[400];
		sprintf(save_title,"%s/optimize/A%d_config%d_Pol%d_Optimization_RCutSlope%.4f.png",plotPath,station,config,pol_select,slope);
		cRcut->SaveAs(save_title);
	}
}