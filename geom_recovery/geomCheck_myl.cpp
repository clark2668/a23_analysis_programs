////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////
////		geom check
////		this piece of software walks through the corrections from calibrationTools.cxx
////		which Thomas gave to Ming-Yuan
////		It applies to all antennas (including cal pulsers!)
////		Need to be used with araROOT 2908!
////////////////////////////////////////////////////////////////////////////////

// C/C++ Includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>

//ROOT includes
#include "AraGeomTool.h"

using namespace std;

int main(int argc, char **argv)
{
	if(argc<2) {  // Check to make sure there are enough arguments to do something meaningful
		std::cout << "Usage requires you to provide input parameter of the form " << basename(argv[0]) << " station" << std::endl;
		return -1;
	}

	int station = atoi(argv[1]);

	if(!(station==2 || station==3)){
		std::cout<<"Station "<<station<<" is not supported by this script!"<<endl;
		return -1;
	}


	AraGeomTool *araGeom = AraGeomTool::Instance();


	//antennas
	for(int i=0; i<16; i++){
		int pol = (int) araGeom->getStationInfo(station)->getAntennaInfo(i)->polType;
		double X = araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[0];
		double Y = araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[1];
		double Z = araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[2];
		double delay = araGeom->getStationInfo(station)->getCableDelay(i);

		//the original locations and positions as in the SQL file
		// cout<<i<<","<<X<<","<<Y<<","<<Z<<","<<delay<<endl; //for print out to csv
	}

	for(int i=0; i<araGeom->getStationInfo(station)->getNumCalAnts(); i++){
		double X = araGeom->getStationInfo(station)->getCalAntennaInfo(i)->antLocation[0];
		double Y = araGeom->getStationInfo(station)->getCalAntennaInfo(i)->antLocation[1];
		double Z = araGeom->getStationInfo(station)->getCalAntennaInfo(i)->antLocation[2];
		//the original locations and positions as in the SQL file
		// cout<<i<<","<<X<<","<<Y<<","<<Z<<endl; //for print out to csv
	}

	//first we read in the geometry and delay corrections from Thomas' file

	ifstream ind;
	char posDelayFile[200];
	sprintf(posDelayFile,"from_Thomas/geometryResultsARA%dE.txt",station);
	ind.open(posDelayFile);
	double posDelayArray[4][4]={{0}};
	double pulserCorr[5]={0};	
	double corrections;
	if(ind.good()){
		for(int i=0;i<4;i++){
			for(int j=0; j<4;j++){
				ind >> corrections;
				posDelayArray[j][i] = corrections;   //j and i have to be switched to use the classical correction file!                                              
			}
		}
		for(int i=0;i<5;i++){
			ind >> corrections;
			pulserCorr[i]=corrections;
		}
	}
	else{
		cout << "couldn't read position and delay correction file!!" << endl;
	}
	ind.close();

	//then we correct the geometry

	Double_t *antloc=0;
	for(int a=0;a<16;a++){ //Loop through all 16 channels  

		double myCorrections[3]={0};
		myCorrections[0]=posDelayArray[a%4][0];
		myCorrections[1]=posDelayArray[a%4][1];
		myCorrections[2]=posDelayArray[a%4][2];
		if(station==2 && a==0) myCorrections[2]+=1.68;
		if(station==3 && a==10) myCorrections[2]+=2.01;

		//the geometry corrections
		// cout<<a<<","<<myCorrections[0]<<","<<myCorrections[1]<<","<<myCorrections[2]<<endl; //for print out to csv

		double myFinal[4]={0};
		antloc = araGeom->getStationInfo(station)->getAntennaInfo(a)->getLocationXYZ();
		myFinal[0] = myCorrections[0]+antloc[0];
		myFinal[1] = myCorrections[1]+antloc[1];
		myFinal[2] = myCorrections[2]+antloc[2];

		//the final corrected position values
		// cout<<a<<","<<myFinal[0]<<","<<myFinal[1]<<","<<myFinal[2]<<endl; //for print out to csv
	}
	
	//then we correct the cal pulsers
	
	int num_pulsers = araGeom->getStationInfo(station)->getNumCalAnts();
	// cout<<"Number of cal pulsers for station "<<num_pulsers<<endl;

	for(int pulser=0; pulser<num_pulsers; pulser++){

		//get the name of the cal pulser
		string locName(&araGeom->getStationInfo(station)->getCalAntennaInfo(pulser)->locationName[0]);
		// cout<<"locName: "<<locName<<endl;

		double myOriginal[4]={0};
		antloc = araGeom->getStationInfo(station)->getCalAntennaInfo(pulser)->getLocationXYZ();
		myOriginal[0] = antloc[0];
		myOriginal[1] = antloc[1];
		myOriginal[2] = antloc[2];

		double myCorrections[3]={0};
		double myFinal[3]={0};
		double realCorrections[3]={0}; //the corrections organized in a readable way

		if( locName[locName.length()-1]  == '5' ){

			//the "packing" orrder of pulserCorr is weird!
			//referring to further up in this code in Thomas' comments
			//the D5 corrections are stored (Z-corr, X-corr, Y-corr)

			myCorrections[0] = pulserCorr[1];
			myCorrections[1] = pulserCorr[2];
			myCorrections[2] = pulserCorr[0];

			myFinal[0] = myOriginal[0] + myCorrections[0];
			myFinal[1] = myOriginal[1] + myCorrections[1];
			myFinal[2] = myOriginal[2] + myCorrections[2];

		} else if ( locName[locName.length()-1] == '6'){

			//how this correction works is weird
			//referring to futher up in this code in Thomas' comments
			//the D6 corrections are (1+distance_correction)*x/y_coordinate
			
			myCorrections[0] = 1. + pulserCorr[4];
			myCorrections[1] = 1. + pulserCorr[4];
			myCorrections[2] = pulserCorr[3];

			myFinal[0] = myOriginal[0] * myCorrections[0]; //yes, multiplication! not a mistake!
			myFinal[1] = myOriginal[1] * myCorrections[1]; //yes, multiplication! not a mistake!
			myFinal[2] = myOriginal[2] + myCorrections[2];

		} else {
			cerr<<"Pulser name undefined\n";
			return -1;
		}

		realCorrections[0] = myFinal[0] - myOriginal[0];
		realCorrections[1] = myFinal[1] - myOriginal[1];
		realCorrections[2] = myFinal[2] - myOriginal[2];

		//the corrections, so that they work "normally" like "new = old + correction"
		// cout<<pulser<<","<<realCorrections[0]<<","<<realCorrections[1]<<","<<realCorrections[2]<<endl;

		//the final values
		// cout<<pulser<<","<<myFinal[0]<<","<<myFinal[1]<<","<<myFinal[2]<<endl;
	}
	
}//close the main program
