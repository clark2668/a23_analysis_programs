////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////	analysis.cxx
////
////	January 2019,  clark.2668@osu.edu
////	This is an example of how you analyze ARA data
////	We will learn how to get a waveform
////	How to make a spectrum
////	We will run this as *compiled* code, NOT run as a ROOT macro
////	This code executes over raw (L0) data
////
////	We will also learn how to use the ARA geom tool
////////////////////////////////////////////////////////////////////////////////

// C/C++ Includes
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <algorithm>

//AraRoot Includes
#include "RawAtriStationEvent.h"
#include "UsefulAtriStationEvent.h"
#include "FFTtools.h"
#include "AraGeomTool.h"
#include "AraQualCuts.h"

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"


using namespace std;

int main(int argc, char **argv)
{
	if(argc<2) {  // Check to make sure there are enough arguments to do something meaningful
		std::cout << "Usage requires you to provide input parameter of the form " << basename(argv[0]) << " <input data file>" << std::endl;
		return -1;
	}

	int station=atoi(argv[1]);
	AraGeomTool *araGeom = AraGeomTool::Instance();

	for(int i=0; i<16; i++){ //loop over antennas
		int pol = (int) araGeom->getStationInfo(station)->getAntennaInfo(i)->polType; //polarization
		double X = araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[0]; //antenna X location
		double Y = araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[1]; //antenna Y location
		double Z = araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[2]; //antenna Z location
		double delay = araGeom->getStationInfo(station)->getCableDelay(i); //the associated cable delay
		string holeName(araGeom->getStationInfo(station)->getAntennaInfo(i)->holeName);
		printf("Hole name for ant %d is %s \n", i, holeName.c_str());
		if(holeName=="BH4")
			printf("	A-ha! BH4 \n");
	}

	/*
	Now we can also see how to use the geom tool to get cal-pulser antenna information
	*/

	for(int i=0; i<araGeom->getStationInfo(station)->getNumCalAnts(); i++){ //loop over number of cal antennas
		double X = araGeom->getStationInfo(station)->getCalAntennaInfo(i)->antLocation[0];
		double Y = araGeom->getStationInfo(station)->getCalAntennaInfo(i)->antLocation[1];
		double Z = araGeom->getStationInfo(station)->getCalAntennaInfo(i)->antLocation[2];
		string locName(&araGeom->getStationInfo(station)->getCalAntennaInfo(i)->locationName[0]);
		string antName(araGeom->getStationInfo(station)->getCalAntennaInfo(i)->getCalAntName());
	}
	
}//close the main program
