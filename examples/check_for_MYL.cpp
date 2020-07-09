////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  check_for_MYL.cpp 
////  code to synchronize eventTree and AraTree2
////
////  Nov 2019
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <iomanip>
#include <sstream>


//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"

using namespace std;

int main(int argc, char **argv)
{

	if(argc<2) {
		std::cout << "Usage\n" << argv[0] << " <input_file> \n";
		return -1;
	}

	/*
	arguments
	0: exec
	1: input_file
	*/
			
	TFile *fp = TFile::Open(argv[1],"READ");
	if(!fp) {
		std::cout << "Can't open file\n";
		return -1;
	}
	TTree *dataVTree = (TTree*) fp->Get("VTree");
	TTree *dataHTree = (TTree*) fp->Get("HTree");
	TTree *dataAllTree = (TTree*) fp->Get("AllTree");

	double corr_val[2];
	double snr_val[2];
	int WFRMS[2];
	int theta_300[2];
	int phi_300[2];
	int theta_41[2];
	int phi_41[2];

	int Refilt[2];
	dataVTree->SetBranchAddress("Refilt_V",&Refilt[0]);
	dataHTree->SetBranchAddress("Refilt_H",&Refilt[1]);

	dataVTree->SetBranchAddress("corr_val_V_new",&corr_val[0]);
	dataVTree->SetBranchAddress("snr_val_V_new",&snr_val[0]);
	dataVTree->SetBranchAddress("wfrms_val_V_new",&WFRMS[0]);
	dataVTree->SetBranchAddress("theta_300_V_new",&theta_300[0]);
	dataVTree->SetBranchAddress("theta_41_V_new",&theta_41[0]);
	dataVTree->SetBranchAddress("phi_300_V_new",&phi_300[0]);
	dataVTree->SetBranchAddress("phi_41_V_new",&phi_41[0]);

	dataHTree->SetBranchAddress("corr_val_H_new",&corr_val[1]);
	dataHTree->SetBranchAddress("snr_val_H_new",&snr_val[1]);
	dataHTree->SetBranchAddress("wfrms_val_H_new",&WFRMS[1]);
	dataHTree->SetBranchAddress("theta_300_H_new",&theta_300[1]);
	dataHTree->SetBranchAddress("theta_41_H_new",&theta_41[1]);
	dataHTree->SetBranchAddress("phi_300_H_new",&phi_300[1]);
	dataHTree->SetBranchAddress("phi_41_H_new",&phi_41[1]);

	int isCal;
	int isSoft;
	int isShort;
	int isCW;
	int isNewBox;

	dataAllTree->SetBranchAddress("cal",&isCal);
	dataAllTree->SetBranchAddress("soft",&isSoft);
	dataAllTree->SetBranchAddress("short",&isShort);
	dataAllTree->SetBranchAddress("CW",&isCW);
	dataAllTree->SetBranchAddress("box",&isNewBox);

	int isSurf[2]; // a surface event after filtering?
	int isSurfEvent_top[2]; // a top event?

	dataAllTree->SetBranchAddress("surf_V_new2",&isSurf[0]);
	dataAllTree->SetBranchAddress("surf_H_new2",&isSurf[1]);

	dataAllTree->SetBranchAddress("surf_top_V",&isSurfEvent_top[0]);
	dataAllTree->SetBranchAddress("surf_top_H",&isSurfEvent_top[1]);

	int isBadEvent;
	double weight;
	int unixTime;
	int isFirstFiveEvent;
	int hasBadSpareChanIssue;
	int hasBadSpareChanIssue2;
	int runNum;
	int eventNumber;
	bool isSpikey;
	bool isCliff;
	bool OutofBandIssue;
	bool bad_v2;
	bool isRFEvent;
	bool isPayloadBlast2;
	int box300;

	dataAllTree->SetBranchAddress("bad",&isBadEvent);
	dataAllTree->SetBranchAddress("weight",&weight);
	dataAllTree->SetBranchAddress("unixTime",&unixTime);
	dataAllTree->SetBranchAddress("isFirstFiveEvent",&isFirstFiveEvent);
	dataAllTree->SetBranchAddress("hasBadSpareChanIssue",&hasBadSpareChanIssue);
	dataAllTree->SetBranchAddress("hasBadSpareChanIssue2",&hasBadSpareChanIssue2);
	dataAllTree->SetBranchAddress("runNum",&runNum);
	dataAllTree->SetBranchAddress("eventNumber",&eventNumber);
	dataAllTree->SetBranchAddress("isSpikey",&isSpikey);
    dataAllTree->SetBranchAddress("isCliff",&isCliff);
    dataAllTree->SetBranchAddress("OutofBandIssue",&OutofBandIssue);
    dataAllTree->SetBranchAddress("bad_v2",&bad_v2);
    dataAllTree->SetBranchAddress("isRFEvent",&isRFEvent);
    dataAllTree->SetBranchAddress("isPayloadBlast2",&isPayloadBlast2);
    dataAllTree->SetBranchAddress("box300",&box300);

	dataVTree->GetEvent(144975);
	dataVTree->GetEvent(144975);
	dataAllTree->GetEvent(144975);

	cout <<"Event number is "<<eventNumber<<endl;
	printf("WFRMS values V %d and H %d \n", WFRMS[0], WFRMS[1]);

}