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
 #include <stdio.h>

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

void TMVAClassificationApplication(int Wmass, TString fname, TString fout, TString fout_tree, TString tree_name, TString bolo_name, TString simu_index ) 
{   
#ifdef __CINT__
   gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

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
   const int num_variables=6;
   TString sEC1="EC1";	TString sEIA="EIA";	TString sEIC="EIC";
   TString sEC2="EC2";	TString sEIB="EIB";	TString sEID="EID";
   vector<TString> s_vec;
   s_vec.push_back(sEC1); s_vec.push_back(sEC2); s_vec.push_back(sEIA);	s_vec.push_back(sEIB); 
   s_vec.push_back(sEIC); s_vec.push_back(sEID);
   Float_t var_val[num_variables];
   for(std::vector<TString>::iterator itTS = s_vec.begin(); itTS != s_vec.end(); itTS++) {
   
       reader->AddVariable(*itTS, &var_val[itTS - s_vec.begin()]);
   }
   

   // --- Book the MVA methods
   TString dir    = TString("../BDT_") + bolo_name +TString("/") + analysis_type +  TString("/weights/mass_") + s_Wmass + TString("_GeV/");
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
   TFile *input = new TFile(fname,"read"); 
   TTree* theTree = (TTree*)input->Get(tree_name);
   std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;
   //Prepare additional output tree
   TFile * treetarget= new TFile( fout_tree,"update" );
   TTree * tout= new TTree(TString("tout") + simu_index,TString("tout") + simu_index);
   float NN_resp;
   tout->Branch("EC1", &var_val[0], "EC1/F");
   tout->Branch("EC2", &var_val[1], "EC2/F");
   tout->Branch("EIA", &var_val[2], "EIA/F");
   tout->Branch("EIB", &var_val[3], "EIB/F");
   tout->Branch("EIC", &var_val[4], "EIC/F");
   tout->Branch("EID", &var_val[5], "EID/F");
   tout->Branch("NN",  &NN_resp   , "NN/F");
   
   // --- Event loop
   std::cout << "--- Select signal sample" << std::endl;
   for(std::vector<TString>::iterator itTS2 = s_vec.begin(); itTS2 != s_vec.end(); itTS2++) {
       theTree->SetBranchAddress(*itTS2, &var_val[itTS2 - s_vec.begin()]);
   }

   // Efficiency calculator for cut method
   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

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
   // --- Write histograms
   tout->Write();
   treetarget->Close();


   TFile *target  = new TFile( fout,"update" );
   if (Use["BDT"])  {histBdt->SetName(TString("MVA_BDT") + simu_index); histBdt->Write(); }


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


   TString bolo_name = TString("FID837");
   TString analysis_type = TString("");
   vector<int> Wmass_vec;
   Wmass_vec.push_back(6); Wmass_vec.push_back(7); Wmass_vec.push_back(10); Wmass_vec.push_back(25);

   TString gen_path = TString("../BDT_") + bolo_name + TString("/") + analysis_type;

   //Recreate all files at start of code since we use "UPDATE"
   for(std::vector<int>::iterator it = Wmass_vec.begin(); it != Wmass_vec.end(); it++) {
      
         TString  s_Wmass   = Form("%d", *it);
         TString sfile_hist = gen_path + TString("/Application/data_mass_") + s_Wmass +TString("_hist.root");
         TString sfile_tree = gen_path + TString("/Application/data_mass_") + s_Wmass +TString("_tree.root");
         remove(sfile_hist);
         remove(sfile_tree);
   }

   for(std::vector<int>::iterator it = Wmass_vec.begin(); it != Wmass_vec.end(); it++) {
      for (int k = 0; k < 100 ; ++k)
      {
         TString  s_Wmass= Form("%d", *it);
         TString s_index = Form("%d", k);
         //Input files
         TString sData          = gen_path + TString("/True_events/ROOT_files/") + bolo_name + TString("_true_events_tree.root");
         //For the hist output
         TString sData_out      = gen_path + TString("/Application/data_mass_") + s_Wmass +TString("_hist.root");
         //For the tree output
         TString sData_out_tree = gen_path + TString("/Application/data_mass_") + s_Wmass +TString("_tree.root");

         TMVAClassificationApplication(*it, sData,sData_out,sData_out_tree, TString("t_new") + s_index, bolo_name, s_index);
      }

   
   }

return 0;

}





