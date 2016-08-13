

#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include "TObject.h"
#include "TF1.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1.h"
#include "TH2F.h"
#include "TH2.h"
#include <stdlib.h>
#include <vector>
#include "TString.h"
#include "TCanvas.h"
#include <math.h>
#include "TGraph.h"
#include "TMath.h"
#include <map>

using namespace std;

TF1 * fCBRD;
TF1 * fcross_ax;
TF1 * ftrigger_eff;
TF1 * fcut_eff;

double integrand(double *x, double *par)
{
	double Ea = x[0];			//axion energy
	double E  = par[0];    		// reconstructed energy in detector
	double reso_0 = par[1];		// reconstruced energy resolution
	double reso_a = par[2];		// reconstruced energy resolution
	double pi = TMath::Pi();	//value of pi

	double sigma = sqrt(pow(reso_0,2) + pow(reso_a*E, 2));
	double bolo_reso_func = 1/(sqrt(2*pi)*sigma)*TMath::Exp(-0.5*pow((Ea-E)/sigma,2));

	// Integrand expression in c/keV.kg.d
	return fcross_ax->Eval(Ea) * fCBRD->Eval(Ea) * ftrigger_eff->Eval(E) * fcut_eff->Eval(E) * bolo_reso_func;
}


double convolved_flux(double * x, double * par)
{
	double mass=par[0];
	TF1* f1 = new TF1("integrand",integrand,mass,12,3);	
	f1->SetParameter(0,x[0]);
	f1->SetParameter(1,par[1]);
	f1->SetParameter(2,par[2]);

	int np = 500;
	double *t1=new double[np];
	double *t2=new double[np];
	f1->CalcGaussLegendreSamplingPoints(np,t1,t2,1e-5);

	double Rate=f1->IntegralFast(np,t1,t2,mass,12);
	// Rate is the expected rate in c/kg.d.keV in the detector
	return Rate;
}


int get_convolved_flux(TString bolo_name, vector<float>& mass_vec, map<TString,float>& map_reso_0, map<TString,float>& map_reso_a)
{	


	TFile * file_CBRD = new TFile("../ROOT_files/Axion/CBRD/CBRD_flux_per_cm2daykeV_gAeto1.root", "read");
	fCBRD = (TF1*)file_CBRD->Get("CBRD");	

	TFile * file_trigger_eff = new TFile( TString("../ROOT_files/Axion/Trigger_eff/") + bolo_name + TString("_trigger_eff.root"), "read");
	ftrigger_eff = (TF1*)file_trigger_eff->Get("trigger_eff");	

	// Remove file before starting loop because of "update option"
	std::remove(TString("../ROOT_files/Axion/CBRD_convolved/") + bolo_name + TString("_flux.root"));

   // Loop over masses
   for(std::vector<float>::iterator mass_it = mass_vec.begin(); mass_it != mass_vec.end(); mass_it++) 
   {

      TString smass = Form("%g", *mass_it);
      TFile * file_cut_eff = new TFile( TString("../ROOT_files/Axion/Cut_eff_signal/") + bolo_name + TString("_cut_eff_signal.root"), "read");
      fcut_eff = (TF1*)file_cut_eff->Get(TString("cut_eff_signal_mass_") + smass);  

		TFile * file_cross_ax = new TFile("../ROOT_files/Axion/Axio_electric/axio_electric_cross_section_cm2perkg.root", "read");
		fcross_ax = (TF1*)file_cross_ax->Get("cross_axio_mA_" + smass);	

		TFile * fout = new TFile(TString("../ROOT_files/Axion/CBRD_convolved/") + bolo_name + TString("_flux_perkgdkeV.root"), "update");
		TF1 * fconvolved = new TF1(bolo_name + TString("_flux_mass_") + smass, convolved_flux, 1, 20, 3);
		fconvolved->SetParameter(0,*mass_it);
		fconvolved->SetParameter(1,map_reso_0[bolo_name]);
		fconvolved->SetParameter(2,map_reso_a[bolo_name]);
		fconvolved->SetNpx(500);
		fconvolved->Write();
		delete fconvolved;
		fout->Close();

   }


	return 0;
}

int main() 
{

   vector<TString> bolo_vec;
   map<TString, float> map_reso_0;
   map<TString, float> map_reso_a;
   vector<float> mass_vec ;

   bolo_vec.push_back(TString("FID824"));
   map_reso_0[TString("FID824")] = 0.111282785961;
   map_reso_a[TString("FID824")] = 0.0147253546521;

   bolo_vec.push_back(TString("FID825"));
   map_reso_0[TString("FID825")] = 0.135668311717;
   map_reso_a[TString("FID825")] = 0.0152163196818;

   bolo_vec.push_back(TString("FID827"));
   map_reso_0[TString("FID827")] = 0.133835788942;
   map_reso_a[TString("FID827")] = 0.0151733039577;

   bolo_vec.push_back(TString("FID837"));
   map_reso_0[TString("FID837")] = 0.142294892862;
   map_reso_a[TString("FID837")] = 0.0133482507448;

   bolo_vec.push_back(TString("FID838"));
   map_reso_0[TString("FID838")] = 0.139757549396;
   map_reso_a[TString("FID838")] = 0.0155155288886;

   bolo_vec.push_back(TString("FID839"));
   map_reso_0[TString("FID839")] = 0.167525826292;
   map_reso_a[TString("FID839")] = 0.0175133482679;

   bolo_vec.push_back(TString("FID841"));
   map_reso_0[TString("FID841")] = 0.148667968873;
   map_reso_a[TString("FID841")] = 0.0187514036549;

   bolo_vec.push_back(TString("FID842"));
   map_reso_0[TString("FID842")] = 0.17897114637;
   map_reso_a[TString("FID842")] = 0.013552382694;

   mass_vec.push_back(0.);
   mass_vec.push_back(0.01);
   mass_vec.push_back(0.05);
   mass_vec.push_back(0.1);
   mass_vec.push_back(0.15);
   mass_vec.push_back(0.2);

   // Loop over bolos 
   for(std::vector<TString>::iterator bolo_it = bolo_vec.begin(); bolo_it != bolo_vec.end(); bolo_it++) 
   {
   		cout << "Processing bolo " << *bolo_it << endl;
   		get_convolved_flux(*bolo_it, mass_vec, map_reso_0, map_reso_a);
   	}

	return 0;
}



