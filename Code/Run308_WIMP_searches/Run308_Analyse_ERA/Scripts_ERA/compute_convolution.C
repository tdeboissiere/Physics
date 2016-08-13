
#include <iostream>
#include <fstream>
#include "TF1.h"
#include "TMath.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TFile.h"
using namespace std;


TF1* fbeta;
TF1* fgamma;
TF1* fgamma_target;
TF1* fPb;

double beta_times_gaussian(double * x, double * par) 
{

	double E0=x[0];			// Erec en keV
	double E=par[0];
	double sigma_heat = 0.56/2.3548;
	double pi = TMath::Pi();
	double prob_heat=(1/(sqrt(2*pi)*sigma_heat))*exp(-0.5*pow((E-E0)/sigma_heat,2));
	double rate=prob_heat*fbeta->Eval(E0);

	return rate;
}

double beta_convolution(double * x, double * par)
{

	double E=x[0];

	TF1* f1 = new TF1("beta",beta_times_gaussian,0,70,1);	
	f1->SetParameter(0,E); 

	int np = 1000;
	double *t1=new double[np];
	double *t2=new double[np];
	f1->CalcGaussLegendreSamplingPoints(np,t1,t2,1e-5);
	double Rate=(f1->IntegralFast(np,t1,t2,0.0, 69));

	return Rate;	
}

double Pb_times_gaussian(double * x, double * par) 
{

	double E0=x[0];			// Erec en keV
	double E=par[0];
	double sigma_heat = 0.56/2.3548;
	double pi = TMath::Pi();
	double prob_heat=(1/(sqrt(2*pi)*sigma_heat))*exp(-0.5*pow((E-E0)/sigma_heat,2));
	double rate=prob_heat*fPb->Eval(E0);

	return rate;
}


double Pb_convolution(double * x, double * par)
{

	double E=x[0];

	TF1* f1 = new TF1("Pb",Pb_times_gaussian,0,40,1);	
	f1->SetParameter(0,E); 

	int np = 1000;
	double *t1=new double[np];
	double *t2=new double[np];
	f1->CalcGaussLegendreSamplingPoints(np,t1,t2,1e-5);
	double Rate=(f1->IntegralFast(np,t1,t2,0,39));

	return Rate;	
}

double gamma_times_gaussian(double * x, double * par) 
{

	double E0=x[0];			// Erec en keV
	double E=par[0];
	double sigma_heat = 0.1388;
	double pi = TMath::Pi();
	double prob_heat=(1/(sqrt(2*pi)*sigma_heat))*exp(-0.5*pow((E-E0)/sigma_heat,2));
	double rate=prob_heat*fgamma->Eval(E0);

	return rate;
}


double gamma_convolution(double * x, double * par)
{

	double E=x[0];

	TF1* f1 = new TF1("gamma",gamma_times_gaussian,0,15,1);	
	f1->SetParameter(0,E); 

	int np = 1000;
	double *t1=new double[np];
	double *t2=new double[np];
	f1->CalcGaussLegendreSamplingPoints(np,t1,t2,1e-5);
	double Rate=(f1->IntegralFast(np,t1,t2,0.0, 15));

	return Rate;	
}

int main() 
{
	/*
	Goal: 	-verify the convolution is computed correclty by looking at the gamma
			spectrum.
			-show whether or not the beta/Pb surf spectrum changes a lot after convolution.
			-If it does not, use it as the "original spectrum" (i.e. deconvoluted)
			This will simplify the event simulation for the low mass analyses
	*/


	// //Load functions to convolve: Beta
	// TString surf_path = TString("../Analyse_FID837/ROOT_files/FID837_Surf_spectrum_extrapol.root");
	// TFile filesurf(surf_path);
	// fbeta = (TF1*)filesurf.Get("S1Beta_extra");
	// fPb = (TF1*)filesurf.Get("S1Pb_extra");

	// //Load functions to convolve: Gamma
	// TString gamma_path = TString("../Analyse_FID837/ROOT_files/FID837_Gamma_spectrum_extrapol_deconv.root");
	// TFile filegamma(gamma_path);
	// fgamma = (TF1*)filegamma.Get("FidGamma_deconv");
	// gamma_path = TString("../Analyse_FID837/ROOT_files/FID837_Gamma_spectrum_extrapol.root");
	// TFile filegamma_target(gamma_path);
	// fgamma_target = (TF1*)filegamma_target.Get("FidGamma_extra");

	//Load functions to convolve: Gamma
	TString gamma_path = TString("./gamma.root");
	TFile filegamma(gamma_path);
	fgamma = (TF1*)filegamma.Get("FidGamma_deconv");
	fgamma_target = (TF1*)filegamma.Get("FidGamma");


	// TH1F* h = new TH1F("h", "h", 100, -1, 2);
	// h->SetMaximum(1.5*fbeta->GetMaximum());
	// TF1 * fbeta_conv = new TF1("beta_conv", beta_convolution, -1,60);
	// fbeta_conv->SetLineColor(kBlue);
	// fbeta_conv->SetLineWidth(4);
	// fbeta_conv->SetNpx(500);

	// TH1F* h = new TH1F("h", "h", 100, 0,5);
	// h->SetMaximum(1.5*fPb->GetMaximum());
	// TF1 * fPb_conv = new TF1("Pb_conv", Pb_convolution, 0,39);
	// fPb_conv->SetLineColor(kBlue);
	// fPb_conv->SetLineWidth(4);
	// fPb_conv->SetNpx(500);

	TH1F* h = new TH1F("h", "h", 100, 0,15);
	h->SetMaximum(1.5*fgamma_target->GetMaximum());
	TF1 * fgamma_conv = new TF1("gamma_conv", gamma_convolution, 0,15);
	fgamma_conv->SetLineColor(kBlue);
	fgamma_conv->SetLineWidth(4);
	fgamma_conv->SetNpx(500);

	TCanvas* c = new TCanvas("c", "c");
	h->Draw();
	gPad->SetLogy();
	// fbeta_conv->Draw("same");
	// fbeta->Draw("same");
	// c->Print("beta_and_conv.png");

	// fPb_conv->Draw("same");
	// fPb->Draw("same");
	// c->Print("Pb_and_conv.png");

	fgamma_conv->Draw("same");
	fgamma_target->Draw("same");
	c->Print("gamma_and_conv.png");
	
return 0;

}
