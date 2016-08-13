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

void TMVAClassification(int Wmass, TString bolo_name, TString analysis_type, TString myMethodList = "" )
{

   //General path
   TString gen_path = TString("/home/irfulx204/mnt/tmain/Desktop/Run308_BDT_simu_corr/");

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

   TString gui_path = TString("/Gui_files/WIMP_mass_")+s_Wmass+TString("_GeV.root");
   TString outfileName(gen_path + TString("BDT_") + bolo_name + TString("/") + analysis_type + gui_path);
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );


   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" );
   //;D;P;G,D

   TString weight_path = TString("/weights/mass_")+s_Wmass+TString("_GeV");
   (TMVA::gConfig().GetIONames()).fWeightFileDir = gen_path + TString("BDT_") + bolo_name + TString("/") + analysis_type + weight_path;


   TString sEC1="EC1";	TString sEIA="EIA";	TString sEIC="EIC";
   TString sEC2="EC2";	TString sEIB="EIB";	TString sEID="EID";
   TString sENR="ENR";
   
   factory->AddVariable(sEC1, sEC1, 'F');
   factory->AddVariable(sEC2, sEC2, 'F');
   
   factory->AddVariable(sEIA, sEIA, 'F');
   factory->AddVariable(sEIB, sEIB, 'F');
   factory->AddVariable(sEIC, sEIC, 'F');
   factory->AddVariable(sEID, sEID, 'F');

   factory->AddSpectator(sENR);

   TString out_path = gen_path + TString("BDT_") + bolo_name + TString("/") + analysis_type;

   TString sS1Beta   = out_path + TString("/Beta_and_Pb/ROOT_files/") + bolo_name +TString("_S1Beta_tree.root");
   TString sS2Beta   = out_path + TString("/Beta_and_Pb/ROOT_files/") + bolo_name +TString("_S2Beta_tree.root");
   TString sS1Pb     = out_path + TString("/Beta_and_Pb/ROOT_files/") + bolo_name +TString("_S1Pb_tree.root");
   TString sS2Pb     = out_path + TString("/Beta_and_Pb/ROOT_files/") + bolo_name +TString("_S2Pb_tree.root");
   TString sFidGamma = out_path + TString("/Gamma/ROOT_files/") + bolo_name +TString("_FidGamma_tree.root");
   TString sS1Gamma  = out_path + TString("/Gamma/ROOT_files/") + bolo_name +TString("_S1Gamma_tree.root");
   TString sS2Gamma  = out_path + TString("/Gamma/ROOT_files/") + bolo_name +TString("_S2Gamma_tree.root");
   TString sheat     = out_path + TString("/Heatonly/ROOT_files/") + bolo_name +TString("_heatonly_tree.root");
   TString sWIMP     = out_path + TString("/WIMP/ROOT_files/") + bolo_name +TString("_WIMP_mass_") +s_Wmass+ TString("_tree.root");
   
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

   TTree *signal        = (TTree*)fWIMP->Get("t_new0");
   TTree *bckg_S1Beta   = (TTree*)fS1Beta->Get("t_new0");
   TTree *bckg_S2Beta   = (TTree*)fS2Beta->Get("t_new0");
   TTree *bckg_S1Pb     = (TTree*)fS1Pb->Get("t_new0");
   TTree *bckg_S2Pb     = (TTree*)fS2Pb->Get("t_new0");
   TTree *bckg_FidGamma = (TTree*)fFidGamma->Get("t_new0");
   TTree *bckg_S1Gamma  = (TTree*)fS1Gamma->Get("t_new0");
   TTree *bckg_S2Gamma  = (TTree*)fS2Gamma->Get("t_new0");
   TTree* bckg_heat     = (TTree*)fheat->Get("t_new0");
   
   factory->AddSignalTree(signal,1);
	factory->AddBackgroundTree(bckg_FidGamma,0.07581);
	factory->AddBackgroundTree(bckg_S1Gamma,0.000942006);
	factory->AddBackgroundTree(bckg_S2Gamma,0.000843515);
	factory->AddBackgroundTree(bckg_S1Beta,0.00697338);
	factory->AddBackgroundTree(bckg_S2Beta,0.00633367);
	factory->AddBackgroundTree(bckg_S1Pb,0.000697241);
	factory->AddBackgroundTree(bckg_S2Pb,0.00159634);
	factory->AddBackgroundTree(bckg_heat,0.906804);

   TCut mycuts = "";
   TCut mycutb = "";

   factory->PrepareTrainingAndTestTree( mycuts,mycutb, "SplitMode=random:!V:NormMode=EqualNumEvents" );

   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=200:minNodeSize=0.4%:MaxDepth=6:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );


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
	TString analysis_type = TString("ana_0.5_0.3_5");
   vector<int> Wmass_vec;
   Wmass_vec.push_back(3); Wmass_vec.push_back(4); Wmass_vec.push_back(5); Wmass_vec.push_back(6); Wmass_vec.push_back(7); Wmass_vec.push_back(10); Wmass_vec.push_back(25);
   // Wmass_vec.push_back(5); Wmass_vec.push_back(6); Wmass_vec.push_back(7); Wmass_vec.push_back(10); Wmass_vec.push_back(25);

   for(std::vector<int>::iterator it = Wmass_vec.begin(); it != Wmass_vec.end(); it++) 
   {
      TMVAClassification(*it, bolo_name, analysis_type);

   }

   return 0;
}
