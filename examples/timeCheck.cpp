////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////
////		Time check; print UTC for the first and last events happen
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
#include "TChain.h"

int main(int argc, char **argv)
{
	if(argc<2) {  // Check to make sure there are enough arguments to do something meaningful
		std::cout << "Usage requires you to provide input parameter of the form " << basename(argv[0]) << " <input file 1>" << std::endl;
		return -1;
	}
		
	TChain chain("eventTree"); //this for the events for the exterior loop
	for(int file=1; file<argc; file++){
		TString fileKey(argv[file]); //a file key
		chain.Add(fileKey); //add files to the chain
	}
	
	RawAtriStationEvent *rawAtriEvPtr=0; //make the raw event pointer
	chain.SetBranchAddress("event",&rawAtriEvPtr); //set the branch address

	chain.GetEvent(0);

	// get  event unixTime
	int unixtime;
	unixtime = rawAtriEvPtr->unixTime;
	time_t test_time = unixtime;
	tm *time = gmtime( &test_time );
	int year = time->tm_year+1900;
	int month = time->tm_mon+1;
	int day = time->tm_mday;
	int hour = time->tm_hour;
	int min = time->tm_min;
	int sec = time->tm_sec;

	printf("First event: %d \n",0);
		printf("	Unixtime: %d \n",unixtime);
		printf("	Year/Month/Day: %d/%d/%d \n",year,month,day);
		printf("	Hour:Min:Sec: %d:%d:%d \n",hour,min,sec);

	chain.GetEvent(chain.GetEntries()-1);
	// get  event unixTime
	unixtime = rawAtriEvPtr->unixTime;
	test_time = unixtime;
	time = gmtime( &test_time );
	year = time->tm_year+1900;
	month = time->tm_mon+1;
	day = time->tm_mday;
	hour = time->tm_hour;
	min = time->tm_min;
	sec = time->tm_sec;

	printf("Last event: %d \n",chain.GetEntries()-1);
		printf("	Unixtime: %d \n",unixtime);
		printf("	Year/Month/Day: %d/%d/%d \n",year,month,day);
		printf("	Hour:Min:Sec: %d:%d:%d \n",hour,min,sec);
	
}//close the main program
