/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication                                      *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2F.h"

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace std;
using namespace TMVA;

void TMVAClassificationApplication(int Wmass, TString fname, TString fout, TString fout_tree, TString tree_name , TString bolo_name, TString analysis_type) 
{   
#ifdef __CINT__
   gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif


   //General path
   TString gen_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/");


   //---------------------------------------------------------------
   TString  s_Wmass= Form("%d", Wmass);

   // This loads the library
   TMVA::Tools::Instance();
   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;

   // --- Create the Reader object
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   const int num_variables=7;
   
   TString sEC1="EC1";	TString sEIA="EIA";	TString sEIC="EIC";
   TString sEC2="EC2";	TString sEIB="EIB";	TString sEID="EID";
   TString sENR="ENR";

   vector<TString> s_vec;
   s_vec.push_back(sEC1); s_vec.push_back(sEC2); s_vec.push_back(sEIA);	s_vec.push_back(sEIB); 
   s_vec.push_back(sEIC); s_vec.push_back(sEID); s_vec.push_back(sENR);

   Float_t var_val[num_variables];
   
   for(std::vector<TString>::iterator itTS = s_vec.begin(); itTS != s_vec.end(); itTS++) {
   
       if (*itTS != TString("ENR")) {reader->AddVariable(*itTS, &var_val[itTS - s_vec.begin()]);}
       else {reader->AddSpectator(*itTS, &var_val[itTS - s_vec.begin()]);}
   }
   

   // --- Book the MVA methods

   TString weight_path = TString("/weights/mass_")+s_Wmass+TString("_GeV/");
   TString dir    = gen_path + TString("BDT_") + bolo_name + TString("/") + analysis_type + weight_path;
   TString prefix = "TMVAClassification";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = TString(it->first) + TString(" method");
         TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
         reader->BookMVA( methodName, weightfile ); 
      }
   }
   
   // Book output histograms
   UInt_t nbin = 100;
   TH1F   *histLk(0), *histLkD(0), *histLkPCA(0), *histLkKDE(0), *histLkMIX(0), *histPD(0), *histPDD(0);
   TH1F   *histPDPCA(0), *histPDEFoam(0), *histPDEFoamErr(0), *histPDEFoamSig(0), *histKNN(0), *histHm(0);
   TH1F   *histFi(0), *histFiG(0), *histFiB(0), *histLD(0), *histNn(0),*histNnbfgs(0),*histNnbnn(0);
   TH1F   *histNnC(0), *histNnT(0), *histBdt(0), *histBdtG(0), *histBdtD(0), *histRf(0), *histSVMG(0);
   TH1F   *histSVMP(0), *histSVML(0), *histFDAMT(0), *histFDAGA(0), *histCat(0), *histPBdt(0);

   if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -2, 2 );

   //Prepare input tree
   
   TFile *input = new TFile(fname,"READ"); 
   TTree* theTree = (TTree*)input->Get(tree_name);
   std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;

   // --- Write tree to file
   TFile *treetarget=new TFile( fout_tree,"RECREATE" );
   //Prepare additional output tree
   TTree * tout= new TTree("tout","tout");
   float NN_resp;
   tout->Branch("EC1", &var_val[0], "EC1/F");
   tout->Branch("EC2", &var_val[1], "EC2/F");
   tout->Branch("EIA", &var_val[2], "EIA/F");
   tout->Branch("EIB", &var_val[3], "EIB/F");
   tout->Branch("EIC", &var_val[4], "EIC/F");
   tout->Branch("EID", &var_val[5], "EID/F");
   tout->Branch("NN",  &NN_resp   , "NN/F");
   
   std::cout << "--- Select signal sample" << std::endl;
   
   for(std::vector<TString>::iterator itTS2 = s_vec.begin(); itTS2 != s_vec.end(); itTS2++) {
       theTree->SetBranchAddress(*itTS2, &var_val[itTS2 - s_vec.begin()]);
   }
   
   // Efficiency calculator for cut method
   Int_t    nSelCutsGA = 0;
   Double_t effS       = 0.7;

   std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) 
   {
      if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
      //Fill the histogram only with events that pass the BDT cut_mass
      theTree->GetEntry(ievt);
      NN_resp=reader->EvaluateMVA("BDT method");
      tout->Fill();
      // --- Return the MVA outputs and fill into histograms
      if (Use["BDT"          ])   histBdt    ->Fill( reader->EvaluateMVA( "BDT method"           ) );
   }

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   tout->Write();
   treetarget->Close();

   TFile *target  = new TFile( fout,"RECREATE" );
   if (Use["BDT"          ])   histBdt    ->Write();

   target->Close();

   std::cout << "--- Created root file: \"TMVApp.root\" containing the MVA output histograms" << std::endl;
  
   input->Close();
   target->Close();
   // delete reader;
   // delete input;
   // delete target;
   // delete treetarget;
   // delete histBdt;
   // delete tout;

   delete reader;
   delete treetarget;
   delete histBdt; 
   
    
   std::cout << "==> TMVAClassificationApplication is done!" << endl << std::endl;
} 



int main() {

   //General path
   TString gen_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/");

	TString bolo_name = TString("FID837");
	TString analysis_type = TString("ana_0.5_0.3_5");

   vector<int> Wmass_vec;
   Wmass_vec.push_back(3); Wmass_vec.push_back(4); Wmass_vec.push_back(5); Wmass_vec.push_back(6); Wmass_vec.push_back(7); Wmass_vec.push_back(10); Wmass_vec.push_back(25);
   // Wmass_vec.push_back(5); Wmass_vec.push_back(6); Wmass_vec.push_back(7); Wmass_vec.push_back(10); Wmass_vec.push_back(25);

   for(std::vector<int>::iterator it = Wmass_vec.begin(); it != Wmass_vec.end(); it++) 
   {

      TString  s_Wmass= Form("%d", *it);
      TString out_path = gen_path + TString("BDT_") + bolo_name + TString("/") + analysis_type;
            
      TString sS1Beta   = out_path +TString("/Beta_and_Pb/ROOT_files/") + bolo_name +TString("_S1Beta_tree.root");
      TString sS2Beta   = out_path +TString("/Beta_and_Pb/ROOT_files/") + bolo_name +TString("_S2Beta_tree.root");
      TString sS1Pb     = out_path +TString("/Beta_and_Pb/ROOT_files/") + bolo_name +TString("_S1Pb_tree.root");
      TString sS2Pb     = out_path +TString("/Beta_and_Pb/ROOT_files/") + bolo_name +TString("_S2Pb_tree.root");
      TString sFidGamma = out_path +TString("/Gamma/ROOT_files/") + bolo_name +TString("_FidGamma_tree.root");
      TString sS1Gamma  = out_path +TString("/Gamma/ROOT_files/") + bolo_name +TString("_S1Gamma_tree.root");
      TString sS2Gamma  = out_path +TString("/Gamma/ROOT_files/") + bolo_name +TString("_S2Gamma_tree.root");   
      TString sheat     = out_path + TString("/Heatonly/ROOT_files/") + bolo_name +TString("_heatonly_tree.root");
      TString sWIMP     = out_path +TString("/WIMP/ROOT_files/") + bolo_name +TString("_WIMP_mass_") +s_Wmass+ TString("_tree.root");


      TString app_path = out_path + TString("/Application");

      //For the hist output

      TString sBeta_S1out   = app_path + TString("/S1Beta_mass_") + s_Wmass + TString("_hist.root");
      TString sBeta_S2out   = app_path + TString("/S2Beta_mass_") + s_Wmass + TString("_hist.root");
      TString sPb_S1out     = app_path + TString("/S1Pb_mass_") + s_Wmass + TString("_hist.root");
      TString sPb_S2out     = app_path + TString("/S2Pb_mass_") + s_Wmass + TString("_hist.root");
      TString sFidGamma_out = app_path + TString("/FidGamma_mass_") + s_Wmass + TString("_hist.root");
      TString sS1Gamma_out  = app_path + TString("/S1Gamma_mass_") + s_Wmass + TString("_hist.root");
      TString sS2Gamma_out  = app_path + TString("/S2Gamma_mass_") + s_Wmass + TString("_hist.root");
      TString sheat_out     = app_path + TString("/heatonly_mass_") + s_Wmass + TString("_hist.root");
      TString sWIMP_out     = app_path + TString("/WIMP_mass_") + s_Wmass + TString("_hist.root");

      //For the tree output
      TString sBeta_S1out_tree   = app_path + TString("/S1Beta_mass_") + s_Wmass + TString("_tree.root");
      TString sBeta_S2out_tree   = app_path + TString("/S2Beta_mass_") + s_Wmass + TString("_tree.root");
      TString sPb_S1out_tree     = app_path + TString("/S1Pb_mass_") + s_Wmass + TString("_tree.root");
      TString sPb_S2out_tree     = app_path + TString("/S2Pb_mass_") + s_Wmass + TString("_tree.root");
      TString sFidGamma_out_tree = app_path + TString("/FidGamma_mass_") + s_Wmass + TString("_tree.root");
      TString sS1Gamma_out_tree  = app_path + TString("/S1Gamma_mass_") + s_Wmass + TString("_tree.root");
      TString sS2Gamma_out_tree  = app_path + TString("/S2Gamma_mass_") + s_Wmass + TString("_tree.root");
      TString sheat_out_tree     = app_path + TString("/heatonly_mass_") + s_Wmass + TString("_tree.root");
      TString sWIMP_out_tree     = app_path + TString("/WIMP_mass_") + s_Wmass + TString("_tree.root");

      TMVAClassificationApplication(*it, sFidGamma,sFidGamma_out,sFidGamma_out_tree, TString("t_new1"), bolo_name, analysis_type);
      TMVAClassificationApplication(*it, sS1Gamma,sS1Gamma_out,sS1Gamma_out_tree, TString("t_new1"), bolo_name, analysis_type);
      TMVAClassificationApplication(*it, sS2Gamma,sS2Gamma_out, sS2Gamma_out_tree,TString("t_new1"), bolo_name, analysis_type);
      TMVAClassificationApplication(*it, sheat,sheat_out,sheat_out_tree, TString("t_new1"), bolo_name, analysis_type);
      TMVAClassificationApplication(*it, sS1Beta,sBeta_S1out, sBeta_S1out_tree, TString("t_new1"), bolo_name, analysis_type);
      TMVAClassificationApplication(*it, sS2Beta,sBeta_S2out, sBeta_S2out_tree,TString("t_new1"), bolo_name, analysis_type);
      TMVAClassificationApplication(*it, sS1Pb,sPb_S1out, sPb_S1out_tree, TString("t_new1"), bolo_name, analysis_type);
      TMVAClassificationApplication(*it, sS2Pb,sPb_S2out, sPb_S2out_tree, TString("t_new1"), bolo_name, analysis_type);
      TMVAClassificationApplication(*it, sWIMP,sWIMP_out,sWIMP_out_tree, TString("t_new1"), bolo_name, analysis_type);

   }


return 0;

}





