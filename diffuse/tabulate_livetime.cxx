////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////
////		Tabulate livetime
////////////////////////////////////////////////////////////////////////////////

// C/C++ Includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <ctime>
#include "time.h" // for time convert

//AraRoot Includes
#include "RawAtriStationEvent.h"

//ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "tools_Cuts.h"

// need this generically
RawAtriStationEvent *rawAtriEvPtr;

int main(int argc, char **argv)
{
	if(argc<2) {  // Check to make sure there are enough arguments to do something meaningful
		std::cout << "Usage requires you to provide input parameter of the form " << basename(argv[0]) << " <station> <config> <input file 1>" << std::endl;
		return -1;
	}

	int station=atoi(argv[1]);
	int config=atoi(argv[2]);
	vector<int> BadRunList=BuildBadRunList(station);

	for(int file=3; file<argc; file++){
		TFile *fp = TFile::Open(argv[file]);
		if(!fp) {
			std::cout << "Can't open file\n";
			return -1;
		}
		TTree *eventTree;
		int runNum;
		eventTree= (TTree*) fp->Get("eventTree");
		if(!eventTree) {
			std::cout << "Can't find eventTree\n";
			return -1;
		}
		eventTree->SetBranchAddress("run",&runNum);
		eventTree->SetBranchAddress("event",&rawAtriEvPtr);

		eventTree->GetEvent(0); // get first event
		cout<<"On run "<<runNum<<endl;

		bool thisIsBadRun = isBadRun(station,runNum,BadRunList);

		eventTree->GetEvent(0); // get first event
		int numEntries=eventTree->GetEntries(); 	
		int unixtime_start = rawAtriEvPtr->unixTime;

		eventTree->GetEvent(numEntries-1); // get last event
		int unixtime_stop = rawAtriEvPtr->unixTime;

		int numBadSecs=0;
		for(int loopSec = unixtime_start; loopSec<unixtime_stop; loopSec++){
			if(isBadLivetime(station,loopSec)){
				numBadSecs++;
			} // do the livetime check
		} // loop seconds inside this thing, looking for bad seconds

		char outfile_name[200];
		sprintf(outfile_name,"/users/PAS0654/osu0673/A23_analysis_new2/livetime/A%d_c%d_livetime.txt",station,config);
		FILE *fout = fopen(outfile_name, "a");
		fprintf(fout,"%d, %d, %d, %d, %d\n",
					runNum,
					thisIsBadRun,
					unixtime_start,
					unixtime_stop,
					numBadSecs
		);	
		fclose(fout);
	} // loop input files
}//close the main program