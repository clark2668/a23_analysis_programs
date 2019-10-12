////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////
////		compute geometric info about the SP and other landmarks
////////////////////////////////////////////////////////////////////////////////

// C/C++ Includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <typeinfo>

//ROOT includes
#include "AraGeomTool.h"
#include "TMath.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace std;

int main(int argc, char **argv)
{

	if(argc<2) {  // Check to make sure there are enough arguments to do something meaningful
		std::cout << "Usage requires you to provide input parameter of the form " << basename(argv[0]) << " station" << std::endl;
		return -1;
	}

	int station = atoi(argv[1]);

	TH1D *distro_of_events = new TH1D("","",25,0,50);

	// load the data, which is csv format (why does ROOT make this so difficult)

	for(int config=1; config<6; config++){

		char filename[100];
		sprintf(filename,"A%d_c%i_100Sample_NumSurfaceEventsPerRun.txt",station,config);
		ifstream infile(filename);
		string line;
		string str;

		//  Read the file
		while (getline(infile, line))
		{   istringstream ss (line);
			vector <string> record;
			
			// import order is:
			// 0: run
			// 1: num v
			// 2: num h
			// 3: num either
			// 4: sum of v and h

			while (getline(ss, str, ','))
				record.push_back(str);
			
			int runNum_this;
			std::stringstream(record[0]) >> runNum_this;

			int numEvents_this;
			std::stringstream(record[3]) >> numEvents_this;

			if(runNum_this!=2678
				&& runNum_this!=2090
				&& runNum_this!=4777
				&& runNum_this!=5649
				&& runNum_this!=5660
				&& runNum_this!=5664
				&& runNum_this!=5666
				&& runNum_this!=5670
				&& runNum_this!=5680
				&& runNum_this!=6445
				&& runNum_this!=6536
				&& runNum_this!=6542
				&& runNum_this!=6635
				&& runNum_this!=6655
				&& runNum_this!=6669
				&& runNum_this!=6733
				// && runNum_this!=2636
				// && runNum_this!=6554
				// && runNum_this!=6705
			){
				distro_of_events->Fill(numEvents_this);
				if(numEvents_this>10)
					printf("Run %d, Num Events %d \n", runNum_this, numEvents_this);

			}
		}
	}
	Int_t nq = 20;
	Double_t xq[nq];  // position where to compute the quantiles in [0,1]
	Double_t yq[nq];  // array to contain the quantiles
	for (Int_t i=0;i<nq;i++) xq[i] = Float_t(i+1)/nq;
	distro_of_events->GetQuantiles(nq,yq,xq);
	for(int i=0; i<nq; i++){
		printf("%.2f quantile is %.2f \n", xq[i], yq[i]);
	}
	distro_of_events->Fit("expo");

	gStyle->SetOptStat(1100);

	TCanvas *c = new TCanvas("","",850,850);
		distro_of_events->Draw("");
		distro_of_events->GetXaxis()->SetTitle("Number of Surface Events");
		distro_of_events->GetYaxis()->SetTitle("Number of Runs");
		distro_of_events->GetXaxis()->SetRangeUser(0,50);
		distro_of_events->GetYaxis()->SetTitleOffset(1.2);
		// distro_of_events->SetLineWidth(3);
		gPad->SetLogy();
	char title[150];
	sprintf(title,"A%d_NumSurfEventsDistro.png",station);
	c->SaveAs(title);

}//close the main program
