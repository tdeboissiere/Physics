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

#include "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/ERA_manipulations/ERA/src/EdwEvent.h"
#include "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/ERA_manipulations/ERA/src/EdwPulse.h"
#include "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/ERA_manipulations/ERA/src/EdwTemplate.h"
#include "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/ERA_manipulations/ERA/src/EdwUtils.h"
#include "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/ERA_manipulations/ERA/src/NoiseSpectrum.h"
#include "/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/ERA_manipulations/ERA/src/FitPulse.h"

using namespace std;




int merge_notrigged_event_tree(TString bolo_name)
{

	gSystem->Load("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/ERA_manipulations/ERA/lib/EraLib.so");

	TChain *chain=new TChain("out","eionbis");
	TString event_dir = TString("../Event_files/") + bolo_name + TString("/");
	TString file_og = event_dir + bolo_name + TString("_og_notrigged_event_tree.root");
	TString file_oh = event_dir + bolo_name + TString("_oh_notrigged_event_tree.root");
	TString file_oi = event_dir + bolo_name + TString("_oi_notrigged_event_tree.root");
	TString file_oj = event_dir + bolo_name + TString("_oj_notrigged_event_tree.root");
	TString file_ok = event_dir + bolo_name + TString("_ok_notrigged_event_tree.root");

	chain->AddFile(file_og);
	chain->AddFile(file_oh);
	chain->AddFile(file_oi);
	chain->AddFile(file_oj);
	chain->AddFile(file_ok);

	TString output_file =  event_dir + bolo_name + TString("_all_notrigged_event_tree.root");

	chain->Merge(output_file);

}

// int merge_simu_event_tree(TString bolo_name)
// {

// 	gSystem->Load("/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA/ERA_manipulations/ERA/lib/EraLib.so");

// 	TChain *chain=new TChain("out","eionbis");
// 	TString event_dir = TString("../Event_files/") + bolo_name + TString("/");
// 	TString file_og = event_dir + bolo_name + TString("_og_simu_event_tree.root");
// 	TString file_oh = event_dir + bolo_name + TString("_oh_simu_event_tree.root");
// 	TString file_oi = event_dir + bolo_name + TString("_oi_simu_event_tree.root");
// 	TString file_oj = event_dir + bolo_name + TString("_oj_simu_event_tree.root");
// 	TString file_ok = event_dir + bolo_name + TString("_ok_simu_event_tree.root");

// 	chain->AddFile(file_og);
// 	chain->AddFile(file_oh);
// 	chain->AddFile(file_oi);
// 	chain->AddFile(file_oj);
// 	chain->AddFile(file_ok);

// 	TString output_file =  event_dir + bolo_name + TString("_all_simu_event_tree.root");

// 	chain->Merge(output_file);

// }


int main(int argc, char *argv[]) 

{

	TString bolo_name = TString(argv[1]);
	// cout << argv[1] << endl;
	merge_notrigged_event_tree(bolo_name);

}