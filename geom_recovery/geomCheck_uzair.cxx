////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////
////		geom check
////		this piece of software walks through the corrections from mainAnalysis.cxx
////		which Thomas gave to Uzair
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
		// cout<<i<<", "<<X<<", "<<Y<<", "<<Z<<", "<<delay<<endl; //for print out to csv
		// printf("Chan: %d \t X: %.5f \t Y: %.5f \t Z: %.5f \t Delay: %.5f \n", i, X, Y, Z, delay);
		// printf("%d , %.5f , %.5f , %.5f, %.5f \n", i, X, Y, Z, delay);
	}

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

	std::vector<double> antl;
	std::vector<std::vector<double> > ant_loc;
	Double_t *antloc=0;


	for(int a=0;a<16;a++){

		double myDelays[3]={0};
		myDelays[0]=posDelayArray[a%4][0];
		myDelays[1]=posDelayArray[a%4][1];
		myDelays[2]=posDelayArray[a%4][2];
		if(station==2 && a==0) myDelays[2]+=1.68;
		if(station==3 && a==10) myDelays[2]+=2.01;

		// cout<<a<<","<<myDelays[0]<<","<<myDelays[1]<<","<<myDelays[2]<<endl;


		double myFinal[4]={0};
		antloc = araGeom->getStationInfo(station)->getAntennaInfo(a)->getLocationXYZ();
		myFinal[0] = myDelays[0]+antloc[0];
		myFinal[1] = myDelays[1]+antloc[1];
		myFinal[2] = myDelays[2]+antloc[2];

		// cout<<a<<","<<myFinal[0]<<","<<myFinal[1]<<","<<myFinal[2]<<endl;

		// if(station==2 && a==0) antloc[2] = antloc[2] + 1.68;
		// if(station==3 && a==10) antloc[2] = antloc[2] + 2.01;

		// printf(
		// 	"Ant: %d \t  antlocX: %.5f \t antlocY: %.5f \t antlocZ: %.5f posDelayArray [0]: %.5f \t  posDelayArray [1]: %.5f \t  posDelayArray [2]: %.5f \t slackArray %.5f \n",
		// 	a,
		// 	antloc[0],
		// 	antloc[1],
		// 	antloc[2],
		// 	posDelayArray[a%4][0],
		// 	posDelayArray[a%4][1],
		// 	posDelayArray[a%4][2],
		// 	slackArray[a%4]
		// );

		// antl.push_back(antloc[0] + posDelayArray[a%4][0]  );
		// antl.push_back(antloc[1] + posDelayArray[a%4][1]);
		// if((a/4)%2==1){
		// 	antl.push_back(antloc[2]+180.0 + posDelayArray[a%4][2] + slackArray[a%4] );
		// }
		// else antl.push_back(antloc[2]+180.0 + posDelayArray[a%4][2]);

		// cerr << "antDepth[1]["<<a<<"] = " << antl[2] << endl;

		// ant_loc.push_back(antl);
		// antl.clear();
    }
	
	double delay = 0;
	double addDelay = 0;

    //now delay stuff
	for(int a=0;a<16;a++) {//Loop through all 16 channels and write the data to TGraphs                                                                                         

		//Set up the delays properly(some are not taken care of in the simulation):                                                                        
		addDelay = 0.0;
		delay= araGeom->getStationInfo(station)->getCableDelay(a);
		if(a/4==0){
			// cout<<"First condition for ant "<<a<<endl;
			// cout<<"Cable delay is "<<delay<<endl;
			// cout<<"posDelayArray is "<<posDelayArray[a%4][3]<<endl;
			addDelay+=(4.0 + posDelayArray[a%4][3]);
		}
		if(a/4==1){
			// cout<<"Second condition for ant "<<a<<endl;
			// cout<<"Cable delay is "<<delay<<endl;
			// cout<<"posDelayArray is "<<posDelayArray[a%4][3]<<endl;
			addDelay+=(12.0 + posDelayArray[a%4][3]);
		}
		if(a/4==2){
			// cout<<"Third condition for ant "<<a<<endl;
			// cout<<"Cable delay is "<<delay<<endl;
			// cout<<"posDelayArray is "<<posDelayArray[a%4][3]<<endl;
			addDelay+=(0.0 + posDelayArray[a%4][3]);
		}
		if(a/4==3){
			// cout<<"Fourth condition for ant "<<a<<endl;
			// cout<<"Cable delay is "<<delay<<endl;
			// cout<<"posDelayArray is "<<posDelayArray[a%4][3]<<endl;
			addDelay+=(8.0 + posDelayArray[a%4][3]);
		}
		// cout<<a<<","<<addDelay<<endl;
		cout<<a<<","<<delay+addDelay<<endl;
	}



	
}//close the main program
