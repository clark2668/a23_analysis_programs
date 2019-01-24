////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////
////		geom check
////		this piece of software walks through the corrections from mainAnalysis.cxx
////		which Thomas gave to Uzair
////		It only applies to the measurement antennas (no cal pulser corrections here)
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
	AraGeomTool *araGeom = AraGeomTool::Instance();

	for(int i=0; i<16; i++){
		int pol = (int) araGeom->getStationInfo(station)->getAntennaInfo(i)->polType;
		double X = araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[0];
		double Y = araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[1];
		double Z = araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[2];
		double delay = araGeom->getStationInfo(station)->getCableDelay(i);

		//the original locations and positions as in the SQL file
		// cout<<i<<", "<<X<<", "<<Y<<", "<<Z<<", "<<delay<<endl; //for print out to csv
	}

	//first we read in the geometry and delay corrections from Thomas' file

	ifstream ind;
	char posDelayFile[200];
	sprintf(posDelayFile,"geometryResultsARA%dE.txt",station);
	ind.open(posDelayFile);
	double posDelayArray[4][4]={{0}};
	double slackArray[4]={0};
	double pulserCorr[4]={0};
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
		// cerr << "Positions and delays read successfully!" << endl;
	}
	else cout << "couldn't read position and delay correction file!!" << endl;
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
	
	//then we correct the delays

	double delay = 0;
	double addDelay = 0;

	for(int a=0;a<16;a++) { //Loop through all 16 channels                                                                                  

		//Set up the delays properly(some are not taken care of in the simulation):                                                                        
		addDelay = 0.0;
		delay= araGeom->getStationInfo(station)->getCableDelay(a);
		if(a/4==0){
			addDelay+=(4.0 + posDelayArray[a%4][3]);
		}
		if(a/4==1){
			addDelay+=(12.0 + posDelayArray[a%4][3]);
		}
		if(a/4==2){
			addDelay+=(0.0 + posDelayArray[a%4][3]);
		}
		if(a/4==3){
			addDelay+=(8.0 + posDelayArray[a%4][3]);
		}
		
		//the corrections
		// cout<<a<<","<<addDelay<<endl; //for print out to csv
		
		//the final corrected delay values
		// cout<<a<<","<<delay+addDelay<<endl; //for print out to csv
	}
	
}//close the main program
