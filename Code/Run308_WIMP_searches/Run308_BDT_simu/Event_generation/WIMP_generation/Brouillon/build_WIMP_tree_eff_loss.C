
#include <string>
#include <iostream>
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TPad.h"
#include "TFile.h"
#include <vector>

using namespace std;


int fill_root_tree (int Wmass, int num_events, float V_fid, float V_surf, float s_C1, float s_C2, float s_IA, float s_IB, float s_IC, float s_ID, TString flag, TString bolo_name)
{

///--------------------------------------------------------------------------------------------------///
///
/// Build a tree 
/// Arguments:
/// Wmass       => int      => Wimp mass
/// num_events  => int      => Number of events in the .txt file
/// V_fid       => float    => Fiduvial bias
/// s_**        => float    => resolution (FWHM/2.3548) for the corresponding channel
///
/// Ouput:
///     a .root file containing EIA/B/C/D and EC1/EC2 simulated from a list of Er
///--------------------------------------------------------------------------------------------------///


// R e a d   f i l e    a n d   c r e a t e     v e c t o r s
// --------------------------------------------------------------------
TString  s_Wmass= Form("%d", Wmass); 
TString filename;



if (flag ==TString("Fid")) {filename=TString("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_") + bolo_name + TString("/standard_resolution/WIMP/Events/recoils_mass_")+s_Wmass+TString("_GeV.txt");}

ifstream file(filename);    
float tempflt;

vector<float> recoil_vec;
vector<float> Q_vec;

for (int k = 0; k < num_events; k += 1)
{
    file >> tempflt;
    recoil_vec.push_back(tempflt); 
    Q_vec.push_back(0.16*pow(tempflt,0.18));   
}   



// C r e a t e     t r e e
// ---------------------------------------

TTree t1("t1","t1");
float  EC1, EC2, EIA, EIB, EIC, EID;
t1.Branch("EC1", &EC1, "EC1/F");
t1.Branch("EC2", &EC2, "EC2/F");
t1.Branch("EIA", &EIA, "EIA/F");
t1.Branch("EIB", &EIB, "EIB/F");
t1.Branch("EIC", &EIC, "EIC/F");
t1.Branch("EID", &EID, "EID/F");



// F i l l     t r e e
// ---------------------------------------

//Prepare random generator
gRandom->SetSeed(0);
TRandom3 u(0);

//Loop over the events

if (flag==TString("Fid")) {
    for (int k=0; k<num_events; k++) {

        float rC1=u.Gaus(0,s_C1);   float rIA=u.Gaus(0,s_IA);   float rIC=u.Gaus(0,s_IC);
        float rC2=u.Gaus(0,s_C2);   float rIB=u.Gaus(0,s_IB);   float rID=u.Gaus(0,s_ID);       

        EIA=rIA;
        EIC=rIC;

        EC1=(1+Q_vec[k]*V_fid/3.)/(1+V_fid/3.)*recoil_vec[k]+rC1;
        EC2=(1+Q_vec[k]*V_fid/3.)/(1+V_fid/3.)*recoil_vec[k]+rC2;       

        EIB=rIB+recoil_vec[k]*Q_vec[k];
        EID=rID+recoil_vec[k]*Q_vec[k];

        if (0.5*(EIB+EID)>1) {t1.Fill();  }
        
    }
}

TString outfilename;

if (flag ==TString("Fid")) {outfilename=TString("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/BDT_") + bolo_name + TString("/standard_resolution/WIMP/Events/") + bolo_name + TString("_WIMP_mass_")+s_Wmass+TString("_ioncut_tree.root");}

ifstream file(filename);

TFile ftree(outfilename, "recreate");
t1.Write();
ftree.Close();


return 0;

}


int intermediate(vector<int>& Wmass_vec, int num_events, float V_fid, float V_surf, float s_C1, float s_C2, float s_IA, float s_IB, float s_IC, float s_ID, TString flag, TString bolo_name) {

///--------------------------------------------------------------------------------------------------///
///
/// Intermediate function 
/// Arguments:
/// Wmass_vec   => vector   => Wimp mass vector
/// num_events  => int      => Number of events in the .txt file
/// V_fid       => float    => Fiduvial bias
/// s_**        => float    => resolution (FWHM/2.3548) for the corresponding channel
///
/// Ouput:
///     a .root file containing EIA/B/C/D and EC1/EC2 simulated from a list of Er
///--------------------------------------------------------------------------------------------------///

for(std::vector<int>::iterator it = Wmass_vec.begin(); it != Wmass_vec.end(); it++) {
    cout << "**********************************************************************" << endl;
    cout << "Now creating tree for WIMPS with " << *it << " GeV mass" << endl;
    cout << "**********************************************************************" << endl;
    fill_root_tree (*it, num_events, V_fid, V_surf, s_C1, s_C2, s_IA, s_IB, s_IC, s_ID, flag, bolo_name);
    }


}






int main(int num_events, TString bolo_name) {

///--------------------------------------------------------------------------------------------------///
///
/// Build a tree for each mass stored in Wmass_vec
/// Arguments:
///     num_events  => int      => the number of events to be simulated
/// Ouput:
///     a .root file containing EIA/B/C/D and EC1/EC2 simulated from a list of Er
///--------------------------------------------------------------------------------------------------///

//------------------------------///
// Other parameters
//------------------------------///

float V_fid=8;
float V_surf=5.5;


float s_C1=0.227195515543;float s_IA=0.405554611857;float s_IC=0.361814166808;
float s_C2=0.204263631731;float s_IB=0.358416850688;float s_ID=0.32656701206;

vector<TString> flag_vec;
flag_vec.push_back(TString("Fid")); 


vector<int> Wmass_vec;


Wmass_vec.push_back(30);

for(std::vector<TString>::iterator it = flag_vec.begin(); it != flag_vec.end(); it++) {

    intermediate(Wmass_vec, num_events, V_fid, V_surf, s_C1, s_C2, s_IA, s_IB, s_IC, s_ID, *it, bolo_name);
    cout << "***********************************" << endl;
    cout << TString("Generating tree for  ") + *it + TString("  as the event type") << endl;
    cout << "***********************************" << endl;
}





return 0;

}













