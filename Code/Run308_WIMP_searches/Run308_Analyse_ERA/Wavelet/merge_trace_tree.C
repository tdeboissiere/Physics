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

using namespace std;




int get_trigged_trace_tree(TString bolo_name)
{

	TChain *chain=new TChain("out","eionbis");
	TString trace_dir = TString("../Trace_files/") + bolo_name + TString("/");
	// TString file_og = trace_dir + bolo_name + TString("_og_trigged_trace_tree.root");
	// TString file_oh = trace_dir + bolo_name + TString("_oh_trigged_trace_tree.root");
	// TString file_oi = trace_dir + bolo_name + TString("_oi_trigged_trace_tree.root");
	// TString file_oj = trace_dir + bolo_name + TString("_oj_trigged_trace_tree.root");
	// TString file_ok = trace_dir + bolo_name + TString("_ok_trigged_trace_tree.root");
	TString file_ol = trace_dir + bolo_name + TString("_ol_trigged_trace_tree.root");
	TString file_pa = trace_dir + bolo_name + TString("_pa_trigged_trace_tree.root");

	// chain->AddFile(file_og);
	// chain->AddFile(file_oh);
	// chain->AddFile(file_oi);
	// chain->AddFile(file_oj);
	// chain->AddFile(file_ok);
	chain->AddFile(file_ol);
	chain->AddFile(file_pa);

	TString output_file =  trace_dir + bolo_name + TString("_ol_pa_trigged_trace_tree.root");

	chain->Merge(output_file);

}

int main() 

{

	TString bolo_name = TString("FID837");
	get_trigged_trace_tree(bolo_name);

}