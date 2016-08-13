
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


TF1* f_wave_eff;


double computerate_sfg(double * Erec, double * par) {

double Er0=Erec[0];			// Erec en keV
double Mchi=par[0];		// Mchi en GeV
double Er=par[1];		//ion energy keVee
double Ei=par[2];		//heat energy keVNR	
double sigma_nuc=1E-5;
double A=72;			//A germanium
double delta=0;			// delta en keV		
double sigma_ion= 0.536/2.3548;  //ion res in keV
double sigma_rec= 0.33;   // heat res in keVNR


const double v0=220;
const double vesc=544;
const double rho=0.3;
const double vearth=230;

const double pi=TMath::Pi();
double Mn=0.932*A ;         // mass of target nucleus, 1 amu = 931.46 MeV/c^2 = 0.93146 GeV/c^2
double Mp=0.938 ;           // proton mass in GeV/c^2
Er0=Er0+0.00001;	           // Eviter valeur nulle..

// Units: Er keV; Mn,Mchi GeV; v0,vesc km/s; rho GeV/cm^3; sigma pb
// rate evts/kg.d.keV
double eV=1.6E-19;              // [J]
double keV=eV*1E3;
double GeV=eV*1E9;
double fm=1E-15;
double pb=1E-40 ;              // 1 picobarn=10^-40 m^2
double hbarc=1.054E-34*3.e8;    // [J*m]
double kmpersec=1000;          // [m/s]
double cm3=1E6;                // [m^3]
double c2=9E16;                // [m^2/s^2]

// SPIN-IDEPENDENT COHERENT DIFFUSION, f_n=f_p
double mu=Mchi*Mn/(Mchi+Mn);                  // WIMP-nucleus reduced mass, in real units it would be mu=Mchi*Mn*GeV/((Mchi+Mn)*c2)
double mup=Mchi*Mp/(Mchi+Mp);                 // WIMP-proton  reduced mass, in real units it would be mup=Mchi*Mp*GeV/((Mchi+Mp)*c2)
double sigma=sigma_nuc*pow(mu/mup,2)*pow(A,2);            // [m^2]

// NUCLEAR FORM FACTOR (from Lewin & Smith, 1996)
double q=sqrt(2*Mn*GeV*Er0*keV) ;                   							// [MeV], array
double ss=0.9*fm;											//Lewin p. 98
double aa=0.52*fm	;										//Lewin p. 98
double cc=(1.23*pow(A,1./3.)-0.6)*fm;
double rn=sqrt(pow(cc,2)+(7/3.)*pow(pi,2)*pow(aa,2)-5.*pow(ss,2));                            		// Lewin (4.11)
double qrn=q*rn/hbarc;
double formfactor=3*exp(-pow(q*ss/hbarc,2)/2.)*(sin(qrn)-qrn*cos(qrn))/(pow(qrn,3)) ;  			// Lewin (4.7)


// Vmin CALCULATION
double vmin=sqrt(Mn*Er0*keV*c2/(2*pow(mu,2)*GeV))/kmpersec;
// RATE CALCULATION
double eta;
double rate;
double x=vmin/v0;                          // Adimensional velocities, SFG (12)
double y=vearth/v0;
double z=vesc/v0;

double w1=TMath::Abs(y-z);
double w2=y+z;


double Nesc=TMath::Erf(z)-2.*z*exp(-pow(z,2))/sqrt(pi);							//  Normalisation constant, SFG (8)
if (z< y && x<w1) {eta=1./(v0*y);}
else if (z> y && x<w1) {eta=1./(2*Nesc*v0*y)*(TMath::Erf(x+y)-TMath::Erf(x-y)-4.*y*exp(-pow(z,2))/sqrt(pi));}
else if (w1 < x && x<w2) {eta=1./(2*Nesc*v0*y)*(TMath::Erf(z)-TMath::Erf(x-y)-2.*(y+z-x)*exp(-pow(z,2))/sqrt(pi));}
else if (w2 < x ) {eta = 0; }

double proba_ion=(1/(sqrt(2*pi)*sigma_ion))*exp(-0.5*pow((Ei-0.16*pow(Er0,1.18))/sigma_ion,2));
double proba_rec=(1/(sqrt(2*pi)*sigma_rec))*exp(-0.5*pow((Er-Er0)/sigma_rec,2));
double prob_tot=proba_ion*proba_rec;

rate=prob_tot*sigma*pb*pow(formfactor,2)*rho*cm3*pow(c2,2)*eta*3600*24/(2*Mchi*pow(mu,2)*GeV*1.E6*kmpersec);  // [counts/kg/day/keV]

return rate;

}





double PDF2D(double * x, double * par)
{

TF1* f1 = new TF1("S0",computerate_sfg,0,30.7,3);	


double Er=x[0];
double Ei=x[1];	
	
f1->SetParameter(0,par[0]); 
f1->SetParameter(1,Er);	
f1->SetParameter(2,Ei);

int np = 100;
double *t1=new double[np];
double *t2=new double[np];
f1->CalcGaussLegendreSamplingPoints(np,t1,t2,1e-5);

double sigma_ion= 0.536/2.3548;  //ion res in keV
double sigma_rec= 0.33;   // heat res in ke

double Rate=(f1->IntegralFast(np,t1,t2,0.0, 29));
// We take EE_to_NR(0.5) = 1.26
// We take EE_to_NR(0.65) = 1.625
double eff_ER = 0.5*(1+TMath::Erf((Er-1.625)/(TMath::Sqrt(2)*sigma_rec)));
// double eff_wave_ER = f_wave_eff->Eval(Er);

// //Implement fiducial cut efficiency
//1 sigma on heat only band = EFid>0.256

if (Ei<0*0.256) {return 0;}
else {return Rate*eff_ER;	}



					
}


double simulate_events(vector<int>& Wmass_vec, TString analysis) 
{

	double Er,Ei;

	TFile file_WIMP(TString("./ROOT_files/WIMP_PDF2D_") + analysis + TString(".root"), "read");
	ofstream file_cut(TString("./Text_files/WIMP_PDF_90_and_99_cut_value_") + analysis + TString(".txt"));

	for(std::vector<int>::iterator it_mass = Wmass_vec.begin(); it_mass != Wmass_vec.end(); it_mass++) 
	{
		TString  s_Wmass= Form("%d", *it_mass); 
		TF2* fWIMP2D = (TF2*)file_WIMP.Get(TString("WIMP_") + s_Wmass + TString("_GeV"));

		//Fill histogram to get maximum (root can't find max for 2D functions, only min
		// Limit to 10 in ionisation to avoid crash. It doesn't influence the contour
		// 1.22 keVNR == 0.5 keVee
		// 2.44 keVNR == 1 keVee
		// 3.58 keVNR = 1.5 keVee
		double xmin=1.22 ;double xmax=30.7; double ymin=0 ; double ymax=15.0;
		int binX=2000; int binY=1000;
		float binXf=2000.; float binYf=1000.;

		TH2F* hWIMP2D= new TH2F("h","h",binX, xmin, xmax, binY, ymin, ymax);
		for (int i = 1; i <= binX; i += 1)
		{

			for (int k = 1; k <=binY; k += 1) 
			{
			// if (i>= 1760 ) cout << i  << "," << k<< "," <<xmin+(xmax-xmin)*i/binXf << "," << ymin+(ymax-ymin)*k/binYf  << "," <<fWIMP2D->Eval(27.4456,14) << endl;
			  hWIMP2D->SetBinContent(i,k,fWIMP2D->Eval(xmin+(xmax-xmin)*i/binXf, ymin+(ymax-ymin)*k/binYf));	

			}
		}


		int max_bin = hWIMP2D->GetMaximumBin();
		int max_bin_x, max_bin_y, max_bin_z;
		hWIMP2D->GetBinXYZ(max_bin, max_bin_x, max_bin_y, max_bin_z);
		double x_max, y_max;
		x_max = hWIMP2D->GetXaxis()->GetBinCenter(max_bin_x);
		y_max = hWIMP2D->GetYaxis()->GetBinCenter(max_bin_y);

		TH1F * h1 = new TH1F("h1", "h1", 10000, 0, 2*fWIMP2D->Eval(x_max, y_max));
		TH2F* hfull= new TH2F("hfull","hfull",binX, xmin, xmax, binY, ymin, ymax);
		TH2F* h90= new TH2F("h90","h90",binX, xmin, xmax, binY, ymin, ymax);


		// First get the histogram which we will integrate to find the 90 % contour
		for (int k = 0; k< 200000; k++) 
		{
			// fWIMP2D->GetRandom2(Er,Ei);
			hWIMP2D->GetRandom2(Er,Ei);
			double PDF_val = fWIMP2D->Eval(Er,Ei);
			h1->Fill(PDF_val);
		}

		// Find the value of the PDF (contour level) such that events with a PDF value between it and the max of the PDF
		// make up 90 % of the signal
		double PDF_cut_val_90 = 0;
		for (int i = 1; i<10000; i++)
		{
			if (h1->Integral(i,10000)/h1->Integral()>0.90 && h1->Integral(i+1,10000)/h1->Integral()<0.90)
			{
				PDF_cut_val_90 = h1->GetBinCenter(i) + 0.5*(h1->GetBinCenter(i+1)-h1->GetBinCenter(i));
			}
		}

		// Find the value of the PDF (contour level) such that events with a PDF value between it and the max of the PDF
		// make up 90 % of the signal
		double PDF_cut_val_99 = 0;
		for (int i = 1; i<10000; i++)
		{
			if (h1->Integral(i,10000)/h1->Integral()>0.99 && h1->Integral(i+1,10000)/h1->Integral()<0.99)
			{
				PDF_cut_val_99 = h1->GetBinCenter(i) + 0.5*(h1->GetBinCenter(i+1)-h1->GetBinCenter(i));
			}
		}

		// cout << PDF_cut_val_90 << "   " << PDF_cut_val_99 << endl;

		//Run a check to see we get 90% of events indeed
		double counter = 0;
		for (int k = 0; k< 200000; k++) 
		{
			// fWIMP2D->GetRandom2(Er,Ei);
			hWIMP2D->GetRandom2(Er,Ei);
			double PDF_val = fWIMP2D->Eval(Er,Ei);
			// if (PDF_val>PDF_cut_val)
			//  Add a ion cut for the new heat trick
			if (PDF_val>PDF_cut_val_99 && PDF_val< PDF_cut_val_90 && Ei>0.7)
			{
				h90->Fill(Er,Ei);
				counter+=1;
			}
			hfull->Fill(Er,Ei);
		}

		cout << s_Wmass <<"," << counter/200000. << endl;
		file_cut << s_Wmass <<","<< PDF_cut_val_90 << "," << PDF_cut_val_99 <<"," << counter/200000. << endl;

		TFile fout("./ROOT_files/hist2D.root", "recreate");
		hfull->Write();
		h90->Write();
		fout.Close();

		// // TFile fout("hist.root", "recreate");
		// // h1->Write();
		// // fout.Close();
		
		delete hfull;
		delete h90;
		delete h1;
		delete hWIMP2D;
		delete fWIMP2D;
	}
	return 0;

}



int main() {


	// TString swave_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Wavelet/ROOT_files/FID837/");
	// TFile filou(swave_path + TString("PSA_cut_eff_NR.root"), "read");
	// // f_wave_eff = (TF1*)filou.Get("PSA_eff_NR");


	//Mass vector 
	vector<int> vmass;
	for (int k = 5; k<= 7; k++) {vmass.push_back(k);}
	vmass.push_back(10); 
	vmass.push_back(25); 

	TString analysis = TString("ana_0.5_0_5");

	// ofstream file_num(TString("./Text_files/WIMP_counts_for_1kgday_") + analysis + TString("2D.txt"));

	// TFile fout(TString("./ROOT_files/WIMP_PDF2D_") + analysis + TString(".root"), "recreate");

	// for(std::vector<int>::iterator it = vmass.begin(); it != vmass.end(); it++){

	//  	//0.5 keVee = 1.22 keVNR
	// 	//1 keVee = 2.44 keVNR
	// 	//1.5 keVee = 3.58 keVNR
	// 	//15 keVee = 30.7 keVNR
	// 	double xmin=1.22 ;double xmax=30.7; double ymin=0 ; double ymax=15.0;
	// 	int binX=200; int binY=100;
	// 	float binXf=200.; float binYf=100.;
	// 	TString  s_Wmass= Form("%d", *it);
	// 	//Increase the bounds of TF2 to avoid numerical problems
	// 	TF2* PDF = new TF2(TString("WIMP_") + s_Wmass + TString("_GeV"), PDF2D,  xmin, 32, ymin, 16, 1);
	// 	// cout << PDF << endl;
	// 	// cout << PDF->Eval(5,5) << endl;
	// 	PDF->SetParameter(0,*it);
	// 	PDF->SetNpx(binX);
	// 	PDF->SetNpy(binY);  

	// 	//Also fill a 2D hist
	// 	TH2F* h2DPDF= new TH2F("h","h",binX, xmin, xmax, binY, ymin, ymax);
	// 	for (int i = 1; i <= binX; i += 1) {
	// 	    for (int k = 1; k <=binY; k += 1) {
	// 	    	h2DPDF->SetBinContent(i,k,PDF->Eval(xmin+(xmax-xmin)*i/binXf, ymin+(ymax-ymin)*k/binYf));	
	// 	    }
	// 	}

	// 	//Compute the integral
	// 	double sum=0;
	// 	for (int i = 1; i <= binX; i += 1) {
	// 	    for (int k = 1; k <=binY; k += 1) {
	// 	    	sum+=h2DPDF->GetBinContent(i,k)*((xmax-xmin)/binXf)*((ymax-ymin)/binYf);
	// 	    }
	// 	}

	// 	cout << "Computing for mass: " << *it << endl;
	// 	file_num << *it << ","<<sum << endl;

		
	// 	PDF->Write();
		
	// 	// TCanvas *c = new TCanvas("c", "c");
	// 	// PDF->Draw("cont4z");
	// 	// c->Print("PDF2D.png");

	// 	delete h2DPDF;
	// 	delete PDF;
	// }
	// fout.Close();
	// file_num.close();

	//Compute the 90% contour for WIMPs
	simulate_events(vmass, analysis);

return 0;




}
