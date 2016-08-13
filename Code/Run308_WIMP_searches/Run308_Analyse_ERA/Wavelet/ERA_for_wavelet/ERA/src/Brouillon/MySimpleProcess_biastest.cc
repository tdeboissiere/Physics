
#include "EdwUtils.h"
#include "EdwEvent.h"
#include "FitPulse.h"

int FillSimpleAmplTree(TChain* chain, BoloStr bolo, string outfilename) {

  // RESTE A COMPLETER LE BASICPROCESS!!

  // 1) Prepare le ntuple
  TFile* file=new TFile(outfilename.c_str(),"RECREATE");
  TTree* ntp = new TTree(("basicntp_"+bolo.Name).c_str(),"Edelweiss Basic Data Tree");
  ntp->SetDirectory(file);
  int SambaNum, Run, triggerbit1, triggerbit2;
  bool isbolotrigger, saturation;
  ULong64_t DateSec, timestamp, gigastamp;
  float temperature;

  ntp->Branch("Run",&Run,"Run/C");
  ntp->Branch("SambaNum",&SambaNum,"SambaNum/i");
  ntp->Branch("TriggerBit1",&triggerbit1,"TriggerBit1/i");
  ntp->Branch("TriggerBit2",&triggerbit2,"TriggerBit2/i");
  ntp->Branch("IsBoloTrigger",&isbolotrigger,"IsBoloTrigger/O");
  ntp->Branch("Saturation",&saturation,"Saturation/O");
  ntp->Branch("DateSec",&DateSec,"DateSec/l");
  ntp->Branch("TimeStamp",&timestamp,"TimeStamp/l");
  ntp->Branch("GigaStamp",&gigastamp,"GigaStamp/I");
  ntp->Branch("Temperature",&temperature,"Temperature/F");
  // Variables definies par voie
  vector<Float_t> simpleampl(bolo.ChTypes.size(),0);
  vector<Int_t> pileup(bolo.ChTypes.size(),0);
  for (UInt_t p=0; p<simpleampl.size(); p++) {
    ntp->Branch(("SimpleAmpl"+bolo.ChTypes[p]).c_str(),&(simpleampl[p]),("SimpleAmpl"+bolo.ChTypes[p]+"/F").c_str());
    ntp->Branch(("PileUp"+bolo.ChTypes[p]).c_str(),&(pileup[p]),("PileUp"+bolo.ChTypes[p]+"/I").c_str());
  }
  // Variables combinees
  Float_t amplchaltot=0, amplcoltot=0;
  if (bolo.Chal1!="NONE" && bolo.Chal2!="NONE")
    ntp->Branch("SimpleAmplChalTot",&amplchaltot,"SimpleAmplChalTot/F");

  // 2) Boucle sur les evts
  
  chain->SetBranchAddress("Run", &Run);
  chain->SetBranchAddress("SambaNum", &SambaNum);
  chain->SetBranchAddress("TriggerBit1",&triggerbit1);
  chain->SetBranchAddress("TriggerBit2",&triggerbit2);
  chain->SetBranchAddress("IsBoloTrigger",&isbolotrigger);
  chain->SetBranchAddress("TimeStamp",&timestamp);
  chain->SetBranchAddress("GigaStamp",&gigastamp);
  chain->SetBranchAddress("Temperature",&temperature);

  for (UInt_t ievt=0; ievt<chain->GetEntries(); ievt++) {
    chain->GetEntry(ievt);
    
    // 3bis) Boucle sur les pulses
    for (UInt_t p=0; p<simpleampl.size(); p++) {
      edwpulse=evt->Pulse(bolo.ChNames[p]);
      if (edwpulse!=NULL) {
	if (edwpulse->IsTrig()) ntpvar_isbolotrigger=1;
	FitPulse fp(edwpulse);
	fp.BasicPreprocess(1);
	if (fp.IsSaturated()) ntpvar_saturation=1;
	pileup[p]=(Int_t)fp.PeakOutsideWindow();
	simpleampl[p]=fp.SimpleAmpl();
      } else simpleampl[p]=0;
    }
    if (bolo.Chal1!="NONE" && bolo.Chal2!="NONE") {
      edwpulse=evt->Pulse(bolo.Chal1);
      edwpulse2=evt->Pulse(bolo.Chal2);
      if (edwpulse!=NULL && edwpulse2!=NULL) {
	FitPulse fp(edwpulse,edwpulse2,edwpulse->Sign(),edwpulse2->Sign());
	fp.BasicPreprocess();
	amplchaltot=fp.SimpleAmpl();
      } else amplchaltot=0;
    }

    file->cd();
    ntp->Fill();
    evt->Clear(); // A FAIRE: Sinon fuite de memoire!!
    if (ievt % 10000 == 0) cout << "Evt " << ievt << endl;
  }
  file->cd();
  ntp->Write("",TObject::kOverwrite);
  file->Close();

  return 0;
}

int main() {

  // Load parameters

  string AnaDir="/home/irfulx204/mnt/tmain/Desktop/New_traces/OF_bias_test";
  string BoloName="FID837";

  string TraceDir=AnaDir+"/ROOT_files/TraceDir/";
  string AmplDir=AnaDir+ "/ROOT_files/AmplDir/";

  string BoloFile=AnaDir+"/Text_files/liste_bolos.txt";
  string RunFile=AnaDir+"/Text_files/liste_runs_debug.txt";

  BoloStr bolo=Read_liste_bolo(BoloFile,BoloName);
  vector<RunStr> listruns;
  listruns = Read_liste_runs(RunFile,bolo);

  // Boucle sur les runs
  for (UInt_t i_run=0; i_run<listruns.size(); i_run++) {
    RunStr run=listruns[i_run];


    string outputfilename=AmplDir+"/basicntp_"+run.Name+"_"+BoloName+".root";
    
    TChain* ch = new TChain("EdwTree","Edelweiss Raw Data Tree");
    ch->AddFile((TraceDir + BoloName + "_" + run.Name + "_simulated_pulse_tree.root").c_str());

    cout << "Computing ampls for run "<<run.Name<<" ("<<ch->GetEntries()<<" evts)"<<endl;
    FillSimpleAmplTree(ch,bolo,outputfilename);
    ch->Delete("");
    
  }

  return 0;
}
