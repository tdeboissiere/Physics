//Use only a restricted number of discriminating variables for a first test

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
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Config.h"
#endif

using namespace std;
using namespace TMVA;

void TMVAClassification(int Wmass, TString bolo_name, TString analysis_type TString myMethodList = "" )
{

   //Mass conversion to string
   TString  s_Wmass= Form("%d", Wmass);

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   TString outfileName( TString("../Run308_BDT_simu_corr/BDT_") + bolo_name + TString("/") + analysis_type + TString("/Gui_files/WIMP_mass_")+s_Wmass+TString("_GeV.root"));
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );


   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" );
   //;D;P;G,D

   (TMVA::gConfig().GetIONames()).fWeightFileDir = TString("../BDT_") + bolo_name + TString("/") + analysis_type + TString("/weights/mass_")+s_Wmass+TString("_GeV");


   TString sEC1="EC1";	TString sEIA="EIA";	TString sEIC="EIC";
   TString sEC2="EC2";	TString sEIB="EIB";	TString sEID="EID";
   
   factory->AddVariable(sEC1, sEC1, 'F');
   factory->AddVariable(sEC2, sEC2, 'F');
   
   factory->AddVariable(sEIA, sEIA, 'F');
   factory->AddVariable(sEIB, sEIB, 'F');
   factory->AddVariable(sEIC, sEIC, 'F');
   factory->AddVariable(sEID, sEID, 'F');

   TString gen_path = TString("../Run308_BDT_simu_corr/BDT_") + bolo_name + TString("/") + analysis_type;

   TString sS1Beta   = gen_path + TString("/Beta_and_Pb/ROOT_files/") + bolo_name +TString("_S1Beta_tree.root");
   TString sS2Beta   = gen_path + TString("/Beta_and_Pb/ROOT_files/") + bolo_name +TString("_S2Beta_tree.root");
   TString sS1Pb     = gen_path + TString("/Beta_and_Pb/ROOT_files/") + bolo_name +TString("_S1Pb_tree.root");
   TString sS2Pb     = gen_path + TString("/Beta_and_Pb/ROOT_files/") + bolo_name +TString("_S2Pb_tree.root");
   TString sFidGamma = gen_path + TString("/Gamma/ROOT_files/") + bolo_name +TString("_FidGamma_tree.root");
   TString sS1Gamma  = gen_path + TString("/Gamma/ROOT_files/") + bolo_name +TString("_S1Gamma_tree.root");
   TString sS2Gamma  = gen_path + TString("/Gamma/ROOT_files/") + bolo_name +TString("_S2Gamma_tree.root");
   TString sheat     = gen_path + TString("/Heatonly/ROOT_files/") + bolo_name +TString("_heatonly_tree.root");
   TString sWIMP     = gen_path + TString("/WIMP/Events/") + bolo_name +TString("_WIMP_mass_") +s_Wmass+ TString("_tree.root");
   
   TFile *fS1Beta   = TFile::Open( sS1Beta );
   TFile *fS2Beta   = TFile::Open( sS2Beta );
   TFile *fS1Pb     = TFile::Open( sS1Pb );
   TFile *fS2Pb     = TFile::Open( sS2Pb );    
   TFile *fFidGamma = TFile::Open( sFidGamma );
   TFile *fS1Gamma  = TFile::Open( sS1Gamma );
   TFile *fS2Gamma  = TFile::Open( sS2Gamma );
   TFile *fheat     = TFile::Open( sheat );
   TFile *fWIMP     = TFile::Open( sWIMP );
         
   // --- Register the training and test trees

   TTree *signal        = (TTree*)fWIMP->Get("t_new");
   TTree *bckg_S1Beta   = (TTree*)fS1Beta->Get("t_new");
   TTree *bckg_S2Beta   = (TTree*)fS2Beta->Get("t_new");
   TTree *bckg_S1Pb     = (TTree*)fS1Pb->Get("t_new");
   TTree *bckg_S2Pb     = (TTree*)fS2Pb->Get("t_new");
   TTree *bckg_FidGamma = (TTree*)fFidGamma->Get("t_new");
   TTree *bckg_S1Gamma  = (TTree*)fS1Gamma->Get("t_new");
   TTree *bckg_S2Gamma  = (TTree*)fS2Gamma->Get("t_new");
   TTree* bckg_heat     = (TTree*)fheat->Get("t_new");
   
   factory->AddSignalTree(signal,1);
	factory->AddBackgroundTree(bckg_FidGamma,0.0212424381079);
	factory->AddBackgroundTree(bckg_S1Gamma,0.000170579212848);
	factory->AddBackgroundTree(bckg_S2Gamma,9.44102424779e-05);
	factory->AddBackgroundTree(bckg_S1Beta,0.00260575355008);
	factory->AddBackgroundTree(bckg_S2Beta,0.00190419777459);
	factory->AddBackgroundTree(bckg_S1Pb,0.000648770480249);
	factory->AddBackgroundTree(bckg_S2Pb,0.000817798765804);
	factory->AddBackgroundTree(bckg_heat,0.972516051866);

   TCut mycuts = "";
   TCut mycutb = "";

   factory->PrepareTrainingAndTestTree( mycuts,mycutb, "SplitMode=random:!V:NormMode=EqualNumEvents" );

   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:minNodeSize=0.04%:MaxDepth=5:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );


   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;

}


int main() 
{
   TString bolo_name = TString("FID837");
   TString analysis_type = TString("");
   vector<int> Wmass_vec;
   Wmass_vec.push_back(6); Wmass_vec.push_back(7); Wmass_vec.push_back(10); Wmass_vec.push_back(25);

   for(std::vector<int>::iterator it = Wmass_vec.begin(); it != Wmass_vec.end(); it++) 
   {
      TMVAClassification(*it, bolo_name, analysis_type);

   }

   return 0;
}
