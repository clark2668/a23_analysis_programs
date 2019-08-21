////////////////////////////////////////////////////////////////////////////////
////	v2_analysis_audit_100.cxx 
////	Check to see that for this data file, there is a correctly
////	formed filter file, reco file, and CWID file
////	all of this must be true to confidently 
////	proceed with the final stage of the analysis
////
////	Aug 2019
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <string>
#include <sstream>

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "tools_CommandLine.h"

using namespace std;

int main(int argc, char **argv)
{

	stringstream ss;

	if(argc<7) {
		std::cout << "Usage\n" << argv[0] << " <station> <year> <filter_file_loc> <reco_file_location> <cw_file_location> <data_file>"<<endl;
		return -1;
	}

	/*
	arguments
	0: exec
	1: station num (2/3)
	2: year
	3: expected location of the filter file
	4: expected location of the reco file
	5: expected location of the CW file
	6: the data file we want to check in the first place
	*/

	int station = atoi(argv[1]);
	int year = atoi(argv[2]);
	string file_location_filter = argv[3];
	string file_location_reco = argv[4];
	string file_location_cw = argv[5];
	string data_file = argv[6];
	string ped_file;
	if(argc==8){
		ped_file= argv[7];
	}
	else{
		ped_file="";
	}

	int nEvents_data=-5;
	int nEvents_filter=-10;
	int nEvents_reco=-20;
	int nEvents_cw=-30;
	int runNum;

	// to start, we need to know the number of real events
	// and we also need to fetch the run number

	TFile *file_data = TFile::Open(data_file.c_str(),"READ");
	if(!file_data){
		std::cout << "Can't open file : data" << "\n";
	}
	if(file_data){
		TTree *inputTree_data = (TTree*) file_data->Get("eventTree");
		if(!inputTree_data) {
			std::cout << "Can't find data eventTree"  << "\n";
		}
		else{
			nEvents_data = inputTree_data->GetEntries();
			inputTree_data->SetBranchAddress("run",&runNum);
			inputTree_data->GetEntry(0);
			file_data->Close();
		}
	}


	// first, the filter file

	char filename_filter[400];
	sprintf(filename_filter,"%s/processed_station_%d_run_%d_filter.root",file_location_filter.c_str(),station,runNum);
	TFile *file_filter = TFile::Open(filename_filter,"READ");
	if(!file_filter) {
		std::cout << "Can't open file : filter" << "\n";
	}
	if(file_filter){
		TTree *inputTree_filter = (TTree*) file_filter->Get("OutputTree");
		if(!inputTree_filter) {
			std::cout << "Can't find OutputTree: filter"  << "\n";
		}
		else{
			nEvents_filter = inputTree_filter->GetEntries();
		}
		file_filter->Close();
	}

	// then, the reco file

	char filename_reco[400];
	sprintf(filename_reco,"%s/processed_station_%d_run_%d_reco.root",file_location_reco.c_str(),station,runNum);
	TFile *file_reco = TFile::Open(filename_reco,"READ");
	if(!file_reco) {
		std::cout << "Can't open file : reco" << "\n";
	}
	if(file_reco){
		TTree *inputTree_reco = (TTree*) file_reco->Get("OutputTreeReco");
		if(!inputTree_reco) {
			std::cout << "Can't find OutputTreeReco: reco"  << "\n";
		}
		else{
			nEvents_reco = inputTree_reco->GetEntries();
		}
		file_reco->Close();
	}

	// then, the cw file

	char filename_cw[400];
	sprintf(filename_cw,"%s/CWID_station_%d_run_%d.root",file_location_cw.c_str(),station,runNum);
	TFile *file_cw = TFile::Open(filename_cw,"READ");
	if(!file_cw) {
		std::cout << "Can't open file : cw" << "\n";
	}
	if(file_cw){
		TTree *inputTree_cw = (TTree*) file_cw->Get("NewCWTree");
		if(!inputTree_cw) {
			std::cout << "Can't find OutputTree: cw"  << "\n";
		}
		else{
			nEvents_cw = inputTree_cw->GetEntries();
		}
		file_cw->Close();
	}

	int goodFilter=0;
	int goodReco=0;
	int goodCW=0;

	// printf(BLUE"Num events data %d, filter %d, reco %d, cw %d\n"RESET, nEvents_data, nEvents_filter, nEvents_reco, nEvents_cw);

	if(nEvents_filter==nEvents_data){
		goodFilter=1;
	}
	else{
		printf("	Filter Error run %d\n", runNum);
		char logfile_bad[200];
		sprintf(logfile_bad,"/home/brianclark/A23_analysis_new2/a23_analysis_scripts/wipac_scripts/100pct/step5-audit/problems_filter_%d.txt",year);
		FILE *fout = fopen(logfile_bad, "a");
		fprintf(fout,"%s %s\n",data_file.c_str(), ped_file.c_str());
		fclose(fout);//close log file
	}

	if(nEvents_reco==nEvents_data){
		goodReco=1;
	}
	else{
		printf("	Reco Error run %d\n", runNum);
		char logfile_bad[200];
		sprintf(logfile_bad,"/home/brianclark/A23_analysis_new2/a23_analysis_scripts/wipac_scripts/100pct/step5-audit/problems_reco_%d.txt",year);
		FILE *fout = fopen(logfile_bad, "a");
		fprintf(fout,"%s %s\n",data_file.c_str(), ped_file.c_str());
		fclose(fout);//close log file
	}

	if(nEvents_cw==nEvents_data){
		goodCW=1;
	}
	else{
		printf("	CW Error run %d\n", runNum);
		char logfile_bad[200];
		sprintf(logfile_bad,"/home/brianclark/A23_analysis_new2/a23_analysis_scripts/wipac_scripts/100pct/step5-audit/problems_cw_%d.txt",year);
		FILE *fout = fopen(logfile_bad, "a");
		fprintf(fout,"%s %s\n",data_file.c_str(), ped_file.c_str());
		fclose(fout);//close log file
	}

	int retVal=1;

	if(
		goodFilter==1
		&&
		goodReco==1
		&&
		goodCW==1
	){
		printf(GREEN"Good run %d\n"RESET, runNum);
		retVal=0;
	}
	else{
		printf(RED"Bad run %d\n"RESET, runNum);
		retVal=-1;
	}
	return retVal;
}
