

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <iostream>
#include <fstream>
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooMCStudy.h"
#include "RooProdPdf.h"
#include "RooNLLVar.h"
#include "RooProfileLL.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include <vector>
#include "TF1.h"
#include "TH1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFile.h"



using namespace RooFit ;
using namespace std;


int get_WIMP_events(float A_target, int Wmass, float vel_0, float vel_esc, float sigma_cross, int num_events, TString bolo_name, TString analysis_type)
{

    ///-----------------------------------///
    ///
    /// Build a list of events 
    /// Arguments:
    ///	A_target 	=> float 	=> target nucleus 
    ///	Mchi		=> int 		=> Wimp mass
    ///	v0 		    => float 	=> earth speed
    ///	vesc		=> float	=> escape velocity
    ///	sigma_cross	=> float	=> cross section
    ///
    /// Ouput:
    ///     a .txt file corresponding to the Wmass which contains the events
    ///-----------------------------------///

    // S e t u p      p d f 
    // ---------------------------------------      
    // Declare variables for the Wimp rate 
    RooRealVar A("A","A",A_target) ;				//Germanium
    RooRealVar E("E","E",0,30);  					//in keV
    RooRealVar Mchi("Mchi","Mchi",Wmass) ;				//in GeV
    RooRealVar v0("v0","v0",vel_0) ;					//in km.s⁻¹
    RooRealVar vesc("vesc","vesc",vel_esc) ;			//in km.s⁻¹
    RooRealVar sigma_nuc("sigma_nuc","sigma_nuc",sigma_cross) ;	//in pb

    // Build WimpRate p.d.f 
    RooWimp WRate("WimpRate","WimpRate",E,Mchi,sigma_nuc,A,v0,vesc) ;
      
    // Generate recoil events 
    // --------------------------------------
    TString  s_Wmass= Form("%d", Wmass); 
    TString filename;
    TString gen_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/");
    RooRandom::randomGenerator()->SetSeed(0);

    //2 simulations: one for training, one for application
    for (int i =0 ; i<2 ; i++)
    {

        //Simulation index
        TString s_index = Form("%d", i);
        filename=gen_path + TString("BDT_") + bolo_name + TString("/") + analysis_type + TString("/WIMP/ROOT_files/recoils_mass_")+s_Wmass+TString("_GeV_") +s_index + TString(".txt");
        RooMCStudy mgr(WRate,RooArgSet(E));
        mgr.generate(1,num_events,kFALSE,filename);
    }

    return 0;
 
}


int intermediate(float A_target, vector<int>& Wmass_vec, float vel_0, float vel_esc, float sigma_cross, vector<int>& num_events, TString bolo_name, TString analysis_type) 
{
    ///-----------------------------------///
    ///
    /// Intermediate function which will generate events 
    /// Arguments:
    ///	A_target 	=> float 	=> target nucleus 
    ///	Mchi		=> int 		=> Wimp mass
    ///	v0 		    => float 	=> earth speed
    ///	vesc		=> float	=> escape velocity
    ///	sigma_cross	=> float	=> cross section
    ///    
    ///
    /// Ouput:
    ///     a .txt file corresponding to the Wmass which contains the events
    ///-----------------------------------///
    for(std::vector<int>::iterator it_mass = Wmass_vec.begin(); it_mass != Wmass_vec.end(); it_mass++) {
        cout << "**********************************************************************" << endl;
        cout << "Now generating events for WIMPS with " << *it_mass << " GeV mass" << endl;
        cout << "**********************************************************************" << endl;
        get_WIMP_events(A_target, *it_mass, vel_0, vel_esc, sigma_cross, num_events[it_mass-Wmass_vec.begin()], bolo_name, analysis_type);
    }

}

void main() 
{

	TString bolo_name = TString("FID837");
	TString analysis_type = TString("ana_0.5_0_5");

    ///-----------------------------------///
    ///
    /// Build a list of events 
    /// Arguments:
    ///    Bolo name
    /// Ouputs:
    ///     .txt files corresponding to given Wmass 
    ///-----------------------------------///


    #ifdef __CINT__
     gROOT->ProcessLineSync(".x /home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/Event_generation/WIMP_generation/RooWimp.cxx+") ;
    #endif


    //------------------------------///
    // Other parameters
    //------------------------------///
    float A_target=72.;
    float vel_0=220.;
    float vel_esc=544.;
    float sigma_cross=1E-6;

    vector<int> Wmass_vec;
    vector<int> num_events;

	Wmass_vec.push_back(5);num_events.push_back(200000);
	Wmass_vec.push_back(6);num_events.push_back(100000);
	Wmass_vec.push_back(7);num_events.push_back(50000);
	Wmass_vec.push_back(10);num_events.push_back(30000);
	Wmass_vec.push_back(25);num_events.push_back(20000);

    intermediate(A_target, Wmass_vec, vel_0, vel_esc, sigma_cross, num_events, bolo_name, analysis_type);

    // exit(1);

    // return 0;

}

