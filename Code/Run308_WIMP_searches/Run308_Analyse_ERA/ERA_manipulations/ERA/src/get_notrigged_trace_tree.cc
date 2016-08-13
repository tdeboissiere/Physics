#include "TF1.h"
#include "TGraph.h"
#include "TF2.h"
#include "TFile.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"
#include "TMath.h"
#include <iostream>
#include "TSystem.h"
#include "TApplication.h"
#include <cstdlib>
#include <vector>
#include <iostream>
#include <numeric>
#include <map>
#include <string>
#include <vector>
#include "TString.h"
#include "TStopwatch.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TCanvas.h"
#include "TH2F.h"
#include <sstream>
#include <fstream>
#include <typeinfo>
#include "TChain.h"
#include "TSpline.h"
// #include "TSpline5.h"
#include <algorithm> 
// #include "/sps/edelweis/tmain/ERA_computations/ERA/src/EdwEvent.h"
// #include "/sps/edelweis/tmain/ERA_computations/ERA/src/EdwPulse.h"
// #include "/sps/edelweis/tmain/ERA_computations/ERA/src/EdwTemplate.h"
// #include "/sps/edelweis/tmain/ERA_computations/ERA/src/EdwUtils.h"
// #include "/sps/edelweis/tmain/ERA_computations/ERA/src/NoiseSpectrum.h"
// #include "/sps/edelweis/tmain/ERA_computations/ERA/src/FitPulse.h"
#include "EdwEvent.h"
#include "EdwPulse.h"
#include "EdwTemplate.h"
#include "EdwUtils.h"
#include "NoiseSpectrum.h"
#include "FitPulse.h"


using namespace std;




int get_notrigged_trace_tree(TString bolo_name, string run_id)
{
	
  // gSystem->Load("/home/irfulx204/mnt/tmain/Desktop/New_ERA/ERA/lib/EraLib.so");
  TString gen_dir("/sps/edelweis/tmain/ERA_computations/");

  TString Run_ID = TString(run_id);

  //Load the list of files
  TString sERA = gen_dir + bolo_name + TString("/Text_files/") + bolo_name +  TString("_ERA_files_") + Run_ID + TString(".txt");
  cout << sERA << endl;
  ifstream fERA(sERA); 
  vector<TString> ERA_vec;
  string temp_str;
  while (fERA >> temp_str) {ERA_vec.push_back(TString(temp_str));} 
  fERA.close();

	//Load tree
  TChain * Tree = new TChain("EdwTree", "EdwTree");
  for(std::vector<TString>::iterator it = ERA_vec.begin(); it != ERA_vec.end(); it++) 
  {
    // cout << "Adding file: " << *it << endl;
    Tree->AddFile(*it);
  }

  EdwEvent * Evt = new EdwEvent();
  Tree->SetBranchAddress("Event", &Evt);
  int nEntries = Tree->GetEntries();

  cout << "Found " << nEntries << " entries"  << endl;


  //Output tree and file
  TString outdir =gen_dir + bolo_name + TString("/") + bolo_name +  TString("_") + Run_ID + TString("_notrigged_trace_tree.root");
  TFile *outf = new TFile(outdir, "recreate");
  TTree* out = new TTree("out", "out");

  int SambaNum;
  char  Run[9];
  int arr_A[1024];
  int arr_B[1024];

  out->Branch("SambaNum", &SambaNum, "SambaNum/I");  out->Branch("Run", &Run, "Run/C"); 
  out->Branch("arr_A",arr_A,"arr_A[1024]/I");
  out->Branch("arr_B",arr_B,"arr_B[1024]/I");



  for (int i = 0; i < nEntries; ++i)
  {
    Tree->GetEntry(i);
    if (Evt->NbPulses() == 6) 
    {

      EdwPulse * pulseA = Evt->Pulse(0);
      EdwPulse * pulseB = Evt->Pulse(1);

      if (pulseA!=NULL && pulseB!=NULL) 
      {
        if (pulseA->IsTrig()==0 && pulseB->IsTrig()==0) 
        {

          strcpy(Run,(Evt->Run()).c_str());
          SambaNum = Evt->SambaNum();

          vector<Short_t> tA = Evt->Pulse(0)->Trace();
          vector<Short_t> tB = Evt->Pulse(1)->Trace();
          for (int k = 0; k<1024; k++) {arr_A[k]=tA[k]; arr_B[k]=tB[k];}

          out->Fill();
        }         
      }
    }

    cout << i << endl;
  }
  out->Write();
  outf->Close();
	return 0;

}

int main(int argc, char *argv[])
{

    TString bolo_name = "FID837";
    if (argc >=2) get_notrigged_trace_tree(bolo_name, argv[1]);
    return 0;
}

