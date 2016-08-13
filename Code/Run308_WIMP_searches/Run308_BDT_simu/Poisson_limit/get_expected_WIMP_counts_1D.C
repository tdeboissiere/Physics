
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


// TF1* f_wave_eff;


double computerate_sfg(double * Erec, double * par) {

double Er=Erec[0];			// Erec en keV
double sigma_nuc=1E-5;		// sigma en picobarn
double Mchi=par[0];		// Mchi en GeV
double A=72;
double delta=0;			// delta en keV		if not keyword_set(delta) then delta=0; INELASTIC DARK MATTER with delta!=0

const double v0=220;
const double vesc=544;//par[3];
const double rho=0.3;
const double vearth=230;

const double pi=TMath::Pi();

// Velocity distribution integral from Phys Rev D 74, 043531 (2006), Savage, Freese, Gondolo

//if n_elements(A) ne 1 then print,"pbl with A"
double Mn=0.932*A ;         // mass of target nucleus, 1 amu = 931.46 MeV/c^2 = 0.93146 GeV/c^2
double Mp=0.938 ;           // proton mass in GeV/c^2
Er=Er+0.00001;	           // Eviter valeur nulle..

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
double q=sqrt(2*Mn*GeV*Er*keV) ;                   							// [MeV], array
double ss=0.9*fm;											//Lewin p. 98
double aa=0.52*fm	;										//Lewin p. 98
double cc=(1.23*pow(A,1./3.)-0.6)*fm;
double rn=sqrt(pow(cc,2)+(7/3.)*pow(pi,2)*pow(aa,2)-5.*pow(ss,2));                            		// Lewin (4.11)
double qrn=q*rn/hbarc;
double formfactor=3*exp(-pow(q*ss/hbarc,2)/2.)*(sin(qrn)-qrn*cos(qrn))/(pow(qrn,3)) ;  			// Lewin (4.7)

// Vmin CALCULATION
double vmin=sqrt(Mn*Er*keV*c2/(2*pow(mu,2)*GeV))/kmpersec;

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

rate=sigma*pb*pow(formfactor,2)*rho*cm3*pow(c2,2)*eta*3600*24/(2*Mchi*pow(mu,2)*GeV*1.E6*kmpersec);  // [counts/kg/day/keV]
return rate;

}



int main() {


	// TString swave_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/Wavelet/ROOT_files/FID837/");
	// TFile filou(swave_path + TString("PSA_cut_eff_NR.root"), "read");
	// f_wave_eff = (TF1*)filou.Get("PSA_eff_NR");


	//Mass vector 
	vector<int> vmass;
	// for (int k = 5; k<= 10; k++) {vmass.push_back(k);}
	vmass.push_back(5);
	vmass.push_back(6);
	vmass.push_back(7);
	vmass.push_back(10);
	vmass.push_back(25);

	TString analysis = TString("ana_1.5_0_5");

	ofstream file_num(TString("./Text_files/WIMP_counts_for_1kgday_") + analysis + TString(".txt"));

	for(std::vector<int>::iterator it = vmass.begin(); it != vmass.end(); it++){

		// Integrate 1D PDF from 0 to 30
		TF1 * fWIMP1D = new TF1("WIMP1D", computerate_sfg, 0, 30, 1);
		fWIMP1D->SetParameter(0, *it);

		cout << "Computing for mass: " << *it << endl;
		file_num << *it << ","<<fWIMP1D->Integral(0,30) << endl;

		delete fWIMP1D;
	}
	file_num.close();


return 0;




}
