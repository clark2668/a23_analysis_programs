////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  araDataReadout.cxx 
////      Just a very simple example that loops over RawAraEvent objects 
////      calibrating them to make a UsefulAraEvent
////
////    Feb 2016,  meures@icecube.wisc.edu 
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <dirent.h>

using namespace std;

//AraRoot Includes
#include "RawIcrrStationEvent.h"

#include "RawAtriStationEvent.h"
#include "UsefulAraStationEvent.h"
//#include "UsefulIcrrStationEvent.h"
#include "UsefulAtriStationEvent.h"
//#include "AtriEventHkData.h"

//Include FFTtools.h if you want to ask the correlation, etc. tools
#include "FFTtools.h"
#include "AraGeomTool.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
//#include "TGraphErrors.h"
#include "TCanvas.h"
//#include "TH1D.h"
#include "TSystem.h"
#include "TMath.h"
//#include "TF1.h"
#include "TChain.h"
//#include "TH2D.h"
//#include "TLine.h"
#include "time.h"

//selfmade includes
///#include "filter_alg.h"
///#include "waveprocessing.h"
//#include "calibrationToolsVs3.h"
//#include "interferometrytools.h"
//#include "evProcessTools.h"
//#include "filterTools.h"
//#include "filterEvent.h"
//#include "reconstructionTools.h"
//#include "timesequencefilter.h"
//for pdestal file search
#include <sys/types.h>
#include <vector>


static inline void loadBar(int x, int n, int r, int w);
int findPedFile(int stationId, string station_string, int runNumber, char *dir_char);
int analyseRun(char *runChar, char *pedChar, TTree *qualityTree);

//****************************************************************************************************//
//*** This is just some uneccessary visual output, which gives you a loadbar, while processing events.//
//****************************************************************************************************//
static inline void loadBar(int x, int n, int r, int w)
{
       // Only update r times.
       if ( x % (n/r) != 0 ) return;
          
       // Calculuate the ratio of complete-to-incomplete.
       float ratio = x/(float)n;
       int   c     = ratio * w;
                   
       // Show the percentage complete.
       printf("%3d%% [", (int)(ratio*100) );
                          
       // Show the load bar.
       for (int x=0; x<c; x++)
       printf("=");
                       
       for (int x=c; x<w; x++)
       printf(" ");
                                  
       // ANSI Control codes to go back to the
       // previous line and clear it.
       printf("]\n\033[F\033[J");
}



//********************************************************************************************************//
//*** A small module to vertically invert waveforms. *****************************************************//
//*** Some channels in ARA03 seem to have inverted amplifier outputs and need to be corrected by this. ***//
//********************************************************************************************************//
void invertGraph(TGraph *gr)
{
	double t1, v1;
	for(int i=0;i<gr->GetN();i++)
	{
		gr->GetPoint(i,t1,v1);
		gr->SetPoint(i,t1,-v1);
	}
}





int findPedFile(int stationId, string station_string, int runNumber, char *dir_char, char *base_char, string yearStr){
//   char dir_char[200];
   string baseString = string(base_char);   
   int ped_run_no =0;
   DIR *dp;
   for(int pe=0;pe<runNumber;pe++)
   {
      ped_run_no = runNumber - pe-1;
      sprintf(dir_char,(baseString+"/data_files/"+station_string+"/pedestals/"+station_string+"/"+yearStr+"/raw_data/run_%06d").c_str(), ped_run_no );
      if((dp  = opendir(dir_char)) == NULL) {}
      else{
	sprintf(dir_char,(baseString+"/data_files/"+station_string+"/pedestals/"+station_string+"/"+yearStr+"/raw_data/run_%06d/pedestalValues.run%06d.dat").c_str(), ped_run_no, ped_run_no );
      string dir = string(dir_char);
	cout << "Opening ped file " << dir_char << endl;
	return ped_run_no;
     }
     if(ped_run_no<1){
	cout << "Warning: No suitable ped-file found! for station " << stationId  << endl;		
	if(stationId==2){
		sprintf(dir_char,(baseString+"/data_files/"+station_string+"/pedestals/"+station_string+"/"+yearStr+"/raw_data/run_001681/pedestalValues.run001681.dat").c_str());
		return 1681;
	}
	if(stationId==3){
		sprintf(dir_char,(baseString+"/data_files/"+station_string+"/pedestals/"+station_string+"/"+yearStr+"/raw_data/run_000520/pedestalValues.run000520.dat").c_str());
		return 520;
	}
	}
   }
}












//int analyseRun(char *currentDir, char *pedChar, int run)
int analyseRun(char *currentDir, char *pedChar, int run)
{
 

	char txtadd[300];
	//	sprintf(txtadd, "/home/ulatif/txtfiles/ARA03_%i_all.txt", run);
        sprintf(txtadd, "/home/ulatif/src/codeForMaggie/SNRs_%i.txt", run);
       	cout<<txtadd<<endl;
		ofstream aout(txtadd);

   	RawAtriStationEvent *rawAtriEvPtr=0;
   	gSystem->Load("libTree");
   	int pass_evts=0;
  	time_t timer1;
  	time_t timer2;
  	double seconds;
   	char runChar[300];
   	sprintf(runChar, "%s/event%i.root", currentDir, run);
   	TFile *signalFile=new TFile(runChar);
   	TTree *eventTree=(TTree*)signalFile->Get("eventTree");
   	eventTree->SetBranchAddress("event",&rawAtriEvPtr);
   	eventTree->GetEntry(0);

	cerr << "THE CURRENT STATION ID IS: " << rawAtriEvPtr->stationId << endl;
   
   	AraEventCalibrator *calib = AraEventCalibrator::Instance();
   	calib->setAtriPedFile(pedChar, rawAtriEvPtr->stationId);

   	std::cerr << "Opened the evntTree\n";
	UsefulAtriStationEvent *calEvPtr=0;

	int numEntries = eventTree->GetEntries();  

  	Int_t eventType = 0;
  	Int_t eventNumber=0;
  	Double_t unixTime=0;
  	Int_t timeStamp=0;

  	double average[16] = {0};
  	int sampleCount[16] = {0};

  	TGraph *gr_s1[16] = {0};
  	TGraph *gr_s2[16] = {0};
        Int_t nievt=0;

//**************************************************************************************************//
//*****here the station correction file is read. the rows are the different strings,****************//
//*****the 4 coloms are the three coordinate-corrections and then the delay-correction.*************//
//**************************************************************************************************//
//*** The structure of the file is:
//***   4 X-corr
//***   4 Y-corr
//***   4 Z-corr
//***   4 delay-corr
//***	5 calpulser corrections:
//***	   NON-reference(So far this is always D5):
//***         Z-corr
//***         X-corr
//***         Y-corr
//***      REFERENCE-pulser:
//***         Z-corr
//***	      distance-corr ( to be applied as: (1 + distance-corr)*x-coord. or y-coord.  )
//***	1 light speed correction
//***	4 slack corrections.
//*** The last five are not valid!
//
	ifstream ind;
	char posDelayFile[200];
	sprintf(posDelayFile, "/home/meures/analysis/analysisCode/geometryResultsARA%dE.txt", rawAtriEvPtr->stationId );
	ind.open(posDelayFile);
	double posDelayArray[4][4]={{0}};
	double slackArray[4]={0};
	double pulserCorr[6] = {0};
	double corrections;
	if(ind.good()){
	  for(int i=0;i<4;i++){
	    for(int j=0; j<4;j++){
	       	ind >> corrections;
		posDelayArray[j][i] = corrections;   //j and i have to be switched to use the classical correction file!
	    }
	  }
	  for(int i=0;i<6;i++){
		ind >> corrections;
	}
	  ind >> corrections;
	  slackArray[0] = corrections;
	  ind >> corrections;
	  slackArray[1] = corrections;
	  ind >> corrections;
	  slackArray[2] = corrections;
	  ind >> corrections;
	  slackArray[3] = corrections;
	  cerr << "Positions and delays read successfully!" << endl;
	}
	else cout << "couldn't read position and delay correction file!!" << endl;
	ind.close();

	int stationId = rawAtriEvPtr->stationId;


 //***Here the station coordinates are read. the geometry correction is included as well as ***//
 //*** two specific reading errors, which were not corrected in araroot. **********************//
 //*** also is the station center set to 180m under the ice. **********************************//


  	AraGeomTool *geom = AraGeomTool::Instance();
  	std::vector<double> antl;
  	std::vector<std::vector<double> > ant_loc;
	Double_t *antloc=0;


	for(int a=0;a<16;a++){
		antloc = geom->getStationInfo(rawAtriEvPtr->stationId)->getAntennaInfo(a)->getLocationXYZ();		
		if(stationId==2 && a==0) antloc[2] = antloc[2] + 1.68;
		if(stationId==3 && a==10) antloc[2] = antloc[2] + 2.01;
			
		antl.push_back(antloc[0] + posDelayArray[a%4][0]  );
		antl.push_back(antloc[1] + posDelayArray[a%4][1]);
		if((a/4)%2==1)antl.push_back(antloc[2]+180.0 + posDelayArray[a%4][2] + slackArray[a%4] );
		else antl.push_back(antloc[2]+180.0 + posDelayArray[a%4][2]);

		cerr << "antDepth[1]["<<a<<"] = " << antl[2] << endl;
	
		ant_loc.push_back(antl);
		antl.clear();
	}


	int cutWaveAlert = 0;
	double delay = 0;
	double addDelay = 0;
	double times, volts;

	for(int eve=0; eve<numEntries;eve++)
	{
		//This will produce a loading bar at the bottom of your screen. In case you want to monitor the progress.
		if(numEntries>100)loadBar(eve, numEntries, 100, 50);
	
		//Data from the tree is written to the rawAtriEvptr location:
     		eventTree->GetEntry(eve);
		
		//We create a pointer to the calibrated event:
	 	calEvPtr = new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);

		//Initialized to check if an event is corrupted:
		cutWaveAlert = 0;
		
   		for(int a=0;a<16;a++)
   		{//Loop through all 16 channels and write the data to TGraphs
		   
		   //Set up the delays properly(some are not taken care of in the simulation):
		   addDelay = 0.0;
      		   delay= geom->getStationInfo(calEvPtr->stationId)->getCableDelay(a);
		   if(a/4==0){addDelay+=(4.0 + posDelayArray[a%4][3]);}
		   if(a/4==1){addDelay+=(12.0 + posDelayArray[a%4][3]);}
		   if(a/4==2){addDelay+=(0.0 + posDelayArray[a%4][3]);}
		   if(a/4==3){addDelay+=(8.0 + posDelayArray[a%4][3]);}

		   //This gets the graph for the specified RF channel:
		   gr_s1[a] = calEvPtr->getGraphFromRFChan(a);
		   //We need to transcribe the Graphs to another array gr_s2, to adjust a few things:
		   gr_s2[a] = new TGraph();
		   
		   if(calEvPtr->stationId==3 && (a==0||a==4||a==8)){invertGraph(gr_s1[a]);}
		   average[a]=0;
		   //*** The following is to avoid reading corrupted waveforms. ***//
		   //*** I encountered only a few of them, so maybe this is not ***//
		   //*** really neccessary anymore. *******************************//
		   if(gr_s1[a]->GetN()<100 ){ cerr<< "event: " << eve << " wow wow wow " << a << ", points: " << gr_s1[a]->GetN() << endl;cutWaveAlert=1;continue;}
		   //Cut of the first 20 ns in new TGraph, because they are often corrupted.
		   int pc = 0;
		   for(int p=0;p<gr_s1[a]->GetN();p++){
			gr_s1[a]->GetPoint(p, times, volts);
			if(times>20.0 && times - delay - addDelay<9999.0){
				gr_s2[a]->SetPoint(pc, times - delay - addDelay, volts);
				average[a]+=volts;
				pc++;
			}
		   }
		
 		}
		//Part of the corruption check, we need to make sure that all new objects are deleted, before we skip to the next event: 
		if(cutWaveAlert==1){
			cerr << "got alert!" << endl;
			delete calEvPtr;
			for(int i=0;i<16;i++){
				delete gr_s2[i];
				delete gr_s1[i];
			}
			cerr << "skipped event!" << endl;
			continue;
		}



		//********************************************************************************************************************//
		//**************************************ADD YOUR OWN CODE HERE!!!!!***************************************************//
		//*********You can write your Graphs to a file now with the TGraph::GetPoint() command(see usage above).**************//
		//********************************************************************************************************************//
		//*********If you are interested in what kind of an event you are looking at (RF, cal-pulser, software ***************//
		//*********trigger) you can use the following checks:*****************************************************************//

		if(rawAtriEvPtr->isCalpulserEvent()){
		  
		  for(Int_t ich=0; ich<16; ich++){ 
		    
		    gr_s2[ich]= FFTtools::getInterpolatedGraph(gr_s2[ich],0.6);   
		    Double_t *yarr=gr_s2[ich]->GetY();      
		    Double_t *xarr=gr_s2[ich]->GetX(); 
		    
		  }
		}
 
		nievt++;   
	}
	
	
	//The pulser events are tagged, so that is easy:
	if(rawAtriEvPtr->isCalpulserEvent() )
	  {
	    eventType=2;//I used 2 for calpulser, can be anything
	  }
	else{
	  //RF- and software triggers can only be distinguished via the length of the waveform. Software triggers are normally only 160ns longs, while RF events have 500ns.
	  if(calEvPtr->blockVec.size()>50){
	    eventType=0;//For RF events
	  }
	  else{
	    eventType=1;//for Software-triggers
	  }
	}
	
	
	//So you don't get any memory leaks, we have to delete the objects which we created inside the loop:
	delete calEvPtr;
	for(int i=0;i<16;i++){
	  delete gr_s2[i];
	  delete gr_s1[i];
	}
}//loop events
cerr << "Loop finished" << endl;
delete signalFile;

cerr << "Events were looped succesfully!" << endl;
return numEntries;
}




int main(int argc, char **argv)
{
  gSystem->Load("libTree");

  if(argc<5){
    printf("Usage %s <root-file directory: typically something like: /data/exp/ARA/2015/filtered/L0/ARA02/0304/run5144 > <stationId: Integer 2 or 3 (ARA01 doesn't work yet)> <Run number: i.e.: 5144> <pedBaseDir: use '/data/user/meures' for now>\n", argv[0]);
    return 0;
  }

  char currentDir[260];
  char pedBaseDir[260];
//  sprintf(currentDir, "/home/meures/AraCalibration/uhen/ara/data/calibration/ARA03/rootCal/run%d", atoi(argv[3]));
  sprintf(currentDir, "%s", argv[1]);
  string getYear = string(argv[1]);
  Int_t stationId = atoi(argv[2]);
  Int_t runNumber = atoi(argv[3]);
  cout<<"getYear "<<getYear<<" , stationId "<<stationId<<" , runNumber "<<runNumber<<endl; 
 sprintf(pedBaseDir, "%s", argv[4]);
  cerr << "Input stationId: " << stationId << endl;

  cout<<"pedBaseDir "<<pedBaseDir<<", currentDir "<<currentDir<<endl;

	int yearStart = getYear.find("201");
	string yearStr = getYear.substr(yearStart,4);
	int yearNumber = 0;
	std::istringstream(yearStr) >> yearNumber;
	cout << "The year is: " << yearNumber << endl;
	cout<<"yearStr "<<yearStr<<endl;
	
	char station_char[30];
	sprintf(station_char, "ARA%02d", stationId);
	cout<<"station_char "<<station_char<<endl;
     
	string station_string = string(station_char);

	cout<<"station_string "<<station_string<<endl;

   int pedRunNo = 0;
   char pedChar[200];
   char runChar[200];
   int analyzedEvents = 0;
   printf("Opening run %s/event%i.root\n", currentDir, runNumber);
   pedRunNo = findPedFile(stationId, station_string, runNumber, &pedChar[0], &pedBaseDir[0], yearStr);

 // TCanvas *plot = new TCanvas("plot","plot");
 // plot->cd(); 
 // plot->Divide(4,2);


 analyzedEvents = analyseRun(&currentDir[0], &pedChar[0], runNumber);

    return 0; 
}








