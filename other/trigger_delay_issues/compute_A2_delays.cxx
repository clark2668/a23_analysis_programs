////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////  v2_analysis_filter.cxx 
////  A23 diffuse, filter events
////
////  Nov 2018
////////////////////////////////////////////////////////////////////////////////

//Includes
#include <iostream>
#include <iomanip>
#include <sstream>

//ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TStyle.h"

// AraRoot
#include "AraGeomTool.h"


// AraSim
#include "Event.h"
#include "Detector.h"
#include "Report.h"
#include "Position.h"
#include "Vector.h"
#include "Settings.h"
#include "RaySolver.h"
#include "IceModel.h"

using namespace std;

int main(int argc, char **argv)
{

	//first, the ray tracer
	Settings *settings = new Settings();
	string setupfile = "setup.txt";
	settings->ReadFile(setupfile);
	settings->NOFZ=1; //yes, variable depth index of refraction
	// settings->RAY_TRACE_ICE_MODEL_PARAMS=0; //use AraSim original default
	int icemodeln = settings->ICE_MODEL + settings->NOFZ*10;
	int crustmodeln = settings->CONSTANTICETHICKNESS * 1000 + settings->CONSTANTCRUST * 100 + settings->FIXEDELEVATION * 10 + 0;
	int moorebayn = settings->MOOREBAY;
	IceModel *icemodel = new IceModel(icemodeln,crustmodeln,moorebayn);// creates
	RaySolver *raysolver = new RaySolver;
	
	//now, the station
	AraGeomTool *araGeom = AraGeomTool::Instance();
	int station=2;

	//setup the cal pulser (D6V)
	Position D6V;
	D6V.SetXYZ(araGeom->getStationInfo(station)->getCalAntennaInfo(1)->antLocation[0],
				araGeom->getStationInfo(station)->getCalAntennaInfo(1)->antLocation[1],
				araGeom->getStationInfo(station)->getCalAntennaInfo(1)->antLocation[2]);
	// printf("D6V=np.array([%.2f, %.2f, %.2f]) \n", D6V.GetX(), D6V.GetY(), D6V.GetZ());
	printf("%.2f, %.2f, %.2f\n", D6V.GetX(), D6V.GetY(), D6V.GetZ());

	for(int i=0; i<8; i++){ //loop over antennas
		double X = araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[0]; //antenna X location
		double Y = araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[1]; //antenna Y location
		double Z = araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[2]; //antenna Z location
		printf("%.2f, %.2f, %.2f\n",i,X,Y,Z);
		// printf("ch%d=np.array([%.2f, %.2f, %.2f])\n",i,X,Y,Z);
		Position target;
		target.SetXYZ(X,Y,Z);
		double T;
		double D;
		int sol;
		raysolver->Solve_Ray_org(D6V.GetX(), D6V.GetY(), D6V.GetZ(),
									target.GetX(), target.GetY(), target.GetZ(),
									T,
									D,
									sol,
									settings);
		// T/=1.e-9; //convert to ns
		// printf("%d %.3f \n", i, T);
	}

}