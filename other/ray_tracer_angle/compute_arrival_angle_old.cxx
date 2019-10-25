#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cstdio>

#include "TMath.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TLine.h"

#include "AraGeomTool.h"

#include "RaySolver.h"
#include "Position.h"
#include "IceModel.h"
#include "Detector.h"
#include "Settings.h"
#include "Vector.h"
#include "Settings.h"


using namespace std;
class EarthModel;

int main(int argc, char **argv)
{

	if(argc<3){
		cout<<"Not enough arguments! Use like "<<argv[0]<<" <station> <antenna> "<<endl;
		return -1;
	}

	int station = atoi(argv[1]);
	int antenna = atoi(argv[2]);

	Settings *settings = new Settings();
	string setupfile = "setup.txt";
	settings->ReadFile(setupfile);
	settings->NOFZ=1; //yes, variable depth index of refraction
	// settings->RAY_TRACE_ICE_MODEL_PARAMS=0; //use AraSim original default

	int icemodeln = settings->ICE_MODEL + settings->NOFZ*10;
	int crustmodeln = settings->CONSTANTICETHICKNESS * 1000 + settings->CONSTANTCRUST * 100 + settings->FIXEDELEVATION * 10 + 0;
	int moorebayn = settings->MOOREBAY;
	IceModel *icemodel = new IceModel(icemodeln,crustmodeln,moorebayn);// creates

	// first, need to compute the detector "core"
	AraGeomTool *araGeom = AraGeomTool::Instance(); // need a geom tool
	double stationRecoCenter_RT[3]={0.};
	for(int i=0; i<16; i++){
		for(int ii=0; ii<3; ii++){
			stationRecoCenter_RT[ii]+=(araGeom->getStationInfo(station)->getAntennaInfo(i)->antLocation[ii]);
		}
	}
	for(int ii=0; ii<3; ii++){
		stationRecoCenter_RT[ii]/=16.;
	}

	// from ray tracer
	// so it wants phi angles from -180 -> 180
	// and theta angles from -90 -> 90
	// in units of radians
	// for(int i=0; i<180; i++){
	// 	double thetaWaveDeg = -90. + 0.5*1. + 1.*i;
	// 	double thetaWaveRad = thetaWaveDeg*TMath::DegToRad();
	// }

	// first, set up the source relative to the station center

	Position source;

	double phiWave = 0.*TMath::DegToRad();
	double thetaWave = -30.*TMath::DegToRad();
	double R = 300;
	printf("Source theta and phi %.2f, %.2f \n", thetaWave, phiWave);

	Double_t xs = R*TMath::Cos(thetaWave)*TMath::Cos(phiWave);
	Double_t ys = R*TMath::Cos(thetaWave)*TMath::Sin(phiWave);
	Double_t zs = R*TMath::Sin(thetaWave);
	
	xs = xs + stationRecoCenter_RT[0];
	ys = ys + stationRecoCenter_RT[1];
	zs = zs + stationRecoCenter_RT[2];
	
	source.SetXYZ(xs, ys, zs);
	printf("Source is at %.2f, %.2f, %.2f \n", xs, ys, zs);

	// now, set up the target

	Position target;
	
	double x1 = araGeom->getStationInfo(station)->getAntennaInfo(antenna)->antLocation[0];
	double y1 = araGeom->getStationInfo(station)->getAntennaInfo(antenna)->antLocation[1];
	double z1 = araGeom->getStationInfo(station)->getAntennaInfo(antenna)->antLocation[2];

	target.SetXYZ(x1,y1,z1);
	printf("Target is at %.2f, %.2f, %.2f \n", x1, y1, z1);


	// invoke the ray solver
	
	vector < vector <double> > outputs; //place for the answer

	/*
		[0][X] = path length
		[1][X] = launch angle
		[2][X] = receipt angle
		[3][X] = reflection angle
		[4][X] = path Time

		X=0 should be direct
		X=1 should be direct/reflected           
	*/

	
	RaySolver *raysolver = new RaySolver;
	
	// raysolver->Solve_Ray(source, target, icemodel, outputs, settings);
	raysolver->Solve_Ray_org(source, target, outputs, settings);

	if(outputs.size()>0){
		cout<<"Direct path length is "<<outputs[0][0]<<endl;
		cout<<"Direct path launch angle is "<<outputs[1][0]*TMath::RadToDeg()<<endl;
		cout<<"Direct path receive angle is "<<outputs[2][0]*TMath::RadToDeg()<<endl;
	}


	vector<double> inTheta;
	vector<double> receiveTheta;

	for(double thetaScan=0.; thetaScan<90.; thetaScan+=1.){

		double this_thetaWave = thetaScan*TMath::DegToRad();
		
		xs = R*TMath::Cos(this_thetaWave)*TMath::Cos(phiWave);
		ys = R*TMath::Cos(this_thetaWave)*TMath::Sin(phiWave);
		zs = R*TMath::Sin(this_thetaWave);
	
		xs = xs + stationRecoCenter_RT[0];
		ys = ys + stationRecoCenter_RT[1];
		zs = zs + stationRecoCenter_RT[2];
		source.SetXYZ(xs, ys, zs);
		
		vector < vector <double> > this_outputs; //place for the answer
		raysolver->Solve_Ray_org(source, target, this_outputs, settings);
		if(this_outputs.size()>0){
			inTheta.push_back(thetaScan);
			receiveTheta.push_back(this_outputs[2][0]*TMath::RadToDeg() -90.);
		}
	}

	TGraph *plot = new TGraph(inTheta.size(),&inTheta[0],&receiveTheta[0]);
	TGraph *plot_straight = new TGraph(inTheta.size(),&inTheta[0],&inTheta[0]);
	TLine *line = new TLine(37.,0,37.,90.);

	TCanvas *c = new TCanvas("","",850,850);
	plot->Draw("AP");
		plot->SetMarkerStyle(kCircle);
		plot->GetXaxis()->SetTitle("Theta Angle in 300m Interferometric Map [deg]");
		plot->GetYaxis()->SetTitle("Theta Angle of Direct Ray Received at the Antenna [deg]");
		plot->SetTitle("Channel 0");
		plot->GetYaxis()->SetRangeUser(0.,90.);
		plot->GetXaxis()->SetRangeUser(0.,90.);
		plot->GetYaxis()->SetTitleOffset(1.2);
	plot_straight->Draw("same");
		plot_straight->SetLineColor(kGray);
		plot_straight->SetLineWidth(2);
	line->Draw("same");
		line->SetLineStyle(9);
		line->SetLineWidth(2);
	TLegend *legend = new TLegend(0.6,0.2,0.9,0.3);
		legend->AddEntry(plot,"Ray Trace Solution","p");
		legend->AddEntry(plot_straight,"1-1 Line","l");
		legend->AddEntry(line,"Surface","l");
		// legend->SetBorderSize(0);  //no border for legend                                                                                                                                                                                         
		// legend->SetFillColor(0);  //fill color is white                                                                                                                                                                                           
		legend->Draw("same");
	char save_title [150];
	sprintf(save_title,"/users/PAS0654/osu0673/A23_analysis_new2/AraRoot/analysis/a23_analysis_programs/other/ray_tracer_angle/map_vs_receive_angle.png");
	c->SaveAs(save_title);

	return 0;
	
}   //end main