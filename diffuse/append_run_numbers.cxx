////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	add_run_numbers.cxx 
////	Dumb dumb, need to add run numbers to save vals files for later...
////
////	June 2019
////////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <sys/stat.h>

#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include "tools_Cuts.h"

int main(int argc, char **argv){
	
	if(argc<3){
		cout<< "Usage\n" << argv[0] << " <1-station> <2-ValForCuts filename 1...> <3-ValForCuts filename 2... >"<<endl;;
		return -1;
	}
	int station = atoi(argv[1]);

	vector<int> BadRunList=BuildBadRunList(station);

	for(int file_num=2; file_num<argc; file_num++){
		string chRun = "run";
		string file = string(argv[file_num]);
		size_t foundRun = file.find(chRun);
		string strRunNum = file.substr(foundRun+4,4);
		int runNum_in = atoi(strRunNum.c_str());
		printf("On file %d, run number %d \n", file_num-1, runNum_in);

		TFile *f = new TFile(argv[file_num],"update");
		TTree *T = (TTree*)f->Get("AllTree");
		int runNum;
		bool badRun;
		TBranch *runBranch = T->Branch("runNum",&runNum);
		TBranch *badRunBranch = T->Branch("badRun",&badRun);

		badRun = isBadRun(station,runNum,BadRunList);

		int numEntries = T->GetEntries();
		for(int i=0; i<numEntries; i++){
			runNum=runNum_in;
			runBranch->Fill();
			badRunBranch->Fill();
		}
		T->Write();
		delete f;
	}
}

// from Rene Brun
// void upd() { 
// 	TFile *f = new TFile("hs.root","update"); 
// 	TTree *T = (TTree*)f->Get("ntuple"); 
// 	float px,py; 
// 	float pt; 
// 	TBranch *bpt = T->Branch("pt",&pt,"pt/F"); 
// 	T->SetBranchAddress("px",&px); 
// 	T->SetBranchAddress("py",&py); 
// 	Long64_t nentries = T->GetEntries(); 
// 	for (Long64_t i=0;i<nentries;i++) { 
// 		T->GetEntry(i); 
// 		pt = TMath::Sqrt(px*px+py*py); 
// 		bpt->Fill(); 
// 	} 
// 	T->Print(); 
// 	T->Write(); 
// 	delete f; 
// }