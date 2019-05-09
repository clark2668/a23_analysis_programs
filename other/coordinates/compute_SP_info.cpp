////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////
////		compute geometric info about the SP
////////////////////////////////////////////////////////////////////////////////

// C/C++ Includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>

//ROOT includes
#include "AraGeomTool.h"
#include "TMath.h"

using namespace std;

int main(int argc, char **argv)
{
	if(argc<2) {  // Check to make sure there are enough arguments to do something meaningful
		std::cout << "Usage requires you to provide input parameter of the form " << basename(argv[0]) << " station" << std::endl;
		return -1;
	}

	int station = atoi(argv[1]);

	/*
		The South Pole
	*/

	AraGeomTool *araGeom = AraGeomTool::Instance();
	Double_t *SP_ARACoord = araGeom->getSouthPole2010_11();
	// printf("SP_ARACoord is %.2f, %.2f, %.2f \n", SP_ARACoord[0],SP_ARACoord[1],SP_ARACoord[2]);

	Double_t SP_A2Coord[3];
	araGeom->convertArrayToStationCoords(station, SP_ARACoord, SP_A2Coord );
	// printf("SP_A2Coord is %.2f, %.2f, %.2f \n", SP_A2Coord[0],SP_A2Coord[1],SP_A2Coord[2]);
	double phi_SP_A2Coord = TMath::ATan2(SP_A2Coord[1],SP_A2Coord[0])*TMath::RadToDeg();
	printf("Phi of SP is %.2f \n", phi_SP_A2Coord);

	/*
		The ICL
	*/

	Double_t *ICL_ARACoord = araGeom->getICLCorner(0);
	for(int i=1; i<=3; i++){
		Double_t *other_corner = araGeom->getICLCorner(i);
		ICL_ARACoord[0]+=other_corner[0];
		ICL_ARACoord[1]+=other_corner[1];
		ICL_ARACoord[2]+=other_corner[2];
	}
	ICL_ARACoord[0]/=4.;
	ICL_ARACoord[1]/=4.;
	ICL_ARACoord[2]/=4.;
	// printf("ICL_ARACoord is %.2f, %.2f, %.2f \n", ICL_ARACoord[0],ICL_ARACoord[1],ICL_ARACoord[2]);

	Double_t ICL_A2Coord[3];
	araGeom->convertArrayToStationCoords(station, ICL_ARACoord, ICL_A2Coord );
	// printf("ICL_A2Coord is %.2f, %.2f, %.2f \n", ICL_A2Coord[0],ICL_A2Coord[1],ICL_A2Coord[2]);
	double phi_ICL_A2Coord = TMath::ATan2(ICL_A2Coord[1],ICL_A2Coord[0])*TMath::RadToDeg();
	printf("Phi of ICL is %.2f \n", phi_ICL_A2Coord);

	
	/*
		WT3
	*/

	Double_t *WT3_ARACoord = araGeom->getWindTurbine(3);
	// printf("WT3_ARACoord is %.2f, %.2f, %.2f \n", WT3_ARACoord[0],WT3_ARACoord[1],WT3_ARACoord[2]);

	Double_t WT3_A2Coord[3];
	araGeom->convertArrayToStationCoords(station, WT3_ARACoord, WT3_A2Coord );
	// printf("WT3_A2Coord is %.2f, %.2f, %.2f \n", WT3_A2Coord[0],WT3_A2Coord[1],WT3_A2Coord[2]);
	double phi_WT3_A2Coord = TMath::ATan2(WT3_A2Coord[1],WT3_A2Coord[0])*TMath::RadToDeg();
	printf("Phi of WT3 is %.2f \n", phi_WT3_A2Coord);


	/*
		The Cal Pulsers
	*/

	//compute the average depth of the station
	double antenna_average[3]={0.};
	for(int i=0; i<16; i++){
		for(int ii=0; ii<3; ii++){
			antenna_average[ii]+=(araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[ii]);
		}
	}
	for(int ii=0; ii<3; ii++){
		antenna_average[ii]/=16.;
	}
	printf("Station center: %.2f, %.2f, %.2f \n", antenna_average[0], antenna_average[1], antenna_average[2]);

	/*
	Now we can also see how to use the geom tool to get cal-pulser antenna information
	*/

	for(int i=0; i<araGeom->getStationInfo(station)->getNumCalAnts(); i++){ //loop over number of cal antennas
		double X = araGeom->getStationInfo(station)->getCalAntennaInfo(i)->antLocation[0];
		double Y = araGeom->getStationInfo(station)->getCalAntennaInfo(i)->antLocation[1];
		double Z = araGeom->getStationInfo(station)->getCalAntennaInfo(i)->antLocation[2];
		string locName(&araGeom->getStationInfo(station)->getCalAntennaInfo(i)->locationName[0]);
		string antName(araGeom->getStationInfo(station)->getCalAntennaInfo(i)->getCalAntName());
		// printf("%s %s: %.2f, %.2f, %.2f  \n", locName.c_str(), antName.c_str(), X, Y , Z);
		double phi = TMath::ATan2(Y-antenna_average[1],X-antenna_average[0])*TMath::RadToDeg();
		double depth_diff = Z-antenna_average[2];
		double horz_dist = sqrt(pow((X-antenna_average[0]),2.0)+pow((Y-antenna_average[1]),2.0));
		double theta = TMath::ATan2(depth_diff,horz_dist)*TMath::RadToDeg();
		// printf("%s %s: %.2f, %.2f, %.2f is %.2f from station center  \n", locName.c_str(), antName.c_str(), X, Y, Z, Z-antenna_average[2]);
		printf("Theta, Phi of %s : %.2f, %.2f \n", antName.c_str(), theta, phi);
	}

}//close the main program
