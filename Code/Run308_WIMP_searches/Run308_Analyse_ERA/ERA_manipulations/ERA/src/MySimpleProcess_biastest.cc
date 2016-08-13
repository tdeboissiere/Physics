
#include "EdwUtils.h"
#include "EdwEvent.h"
#include "FitPulse.h"

int FillSimpleAmplTree(TChain* chain, BoloStr bolo, string outfilename) {

  // RESTE A COMPLETER LE BASICPROCESS!!

  // 1) Prepare le ntuple
  TFile* file=new TFile(outfilename.c_str(),"RECREATE");
  TTree* ntp = new TTree(("basicntp_"+bolo.Name).c_str(),"Edelweiss Basic Data Tree");
  ntp->SetDirectory(file);
  char ntpvar_run[9], ntpvar_partition[4]; // ! size = real size+1 (terminate the char string by zero..)
  UInt_t ntpvar_sambanum, ntpvar_triggerbit1, ntpvar_triggerbit2=0;
  ULong64_t ntpvar_datesec, ntpvar_timestamp=0; Int_t ntpvar_gigastamp=0;
  Bool_t ntpvar_isbolotrigger=0, ntpvar_saturation=0;
  Float_t ntpvar_temperature=0;
  strcpy(ntpvar_run,"aa00a000"); strcpy(ntpvar_partition,"000");
  ntp->Branch("Run",&ntpvar_run,"Run/C");
  ntp->Branch("Partition",&ntpvar_partition,"Partition/C");
  ntp->Branch("SambaNum",&ntpvar_sambanum,"SambaNum/i");
  ntp->Branch("TriggerBit1",&ntpvar_triggerbit1,"TriggerBit1/i");
  ntp->Branch("TriggerBit2",&ntpvar_triggerbit2,"TriggerBit2/i");
  ntp->Branch("IsBoloTrigger",&ntpvar_isbolotrigger,"IsBoloTrigger/O");
  ntp->Branch("Saturation",&ntpvar_saturation,"Saturation/O");
  ntp->Branch("DateSec",&ntpvar_datesec,"DateSec/l");
  ntp->Branch("TimeStamp",&ntpvar_timestamp,"TimeStamp/l");
  ntp->Branch("GigaStamp",&ntpvar_gigastamp,"GigaStamp/I");
  ntp->Branch("Temperature",&ntpvar_temperature,"Temperature/F");
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
  if (bolo.Col1!="NONE" && bolo.Col2!="NONE")
    ntp->Branch("SimpleAmplColTot",&amplcoltot,"SimpleAmplColTot/F");

  // 2) Boucle sur les evts
  EdwEvent* evt = new EdwEvent();
  chain->SetBranchAddress("Event",&evt);
  EdwPulse* edwpulse = new EdwPulse();
  EdwPulse* edwpulse2 = new EdwPulse();
  for (UInt_t ievt=0; ievt<chain->GetEntries(); ievt++) {
    chain->GetEntry(ievt);
    strcpy(ntpvar_run,(evt->Run()).c_str());
    strcpy(ntpvar_partition,(evt->Partition()).c_str());
    ntpvar_sambanum=evt->SambaNum();
    ntpvar_triggerbit1=evt->TriggerBit(1); ntpvar_triggerbit2=evt->TriggerBit(2);
    ntpvar_isbolotrigger=0;
    ntpvar_saturation=0;
    ntpvar_datesec=evt->DateSec();
    ntpvar_timestamp=evt->TimeStamp(); ntpvar_gigastamp=evt->GigaStamp();
    ntpvar_temperature=evt->Temperature();
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
    if (bolo.Col1!="NONE" && bolo.Col2!="NONE") {
      edwpulse=evt->Pulse(bolo.Col1);
      edwpulse2=evt->Pulse(bolo.Col2);
      if (edwpulse!=NULL && edwpulse2!=NULL) {
	FitPulse fp(edwpulse,edwpulse2,edwpulse->Sign(),edwpulse2->Sign());
	fp.BasicPreprocess();
	amplcoltot=fp.SimpleAmpl();
      } else amplcoltot=0;
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

int main(int argc, char *argv[]) {

  // Load parametres...
  string AnaDir="/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA";
  string BoloName=argv[1];

  string TraceDir=AnaDir+"/Event_files/" + BoloName + "/Event_simu/";
  string AmplDir=AnaDir+ "/Amp_files/" + BoloName + "/Amp_simu/";

  string BoloFile=AnaDir+"/ERA_manipulations/Text_files/liste_bolos.txt";
  string RunFile=AnaDir+"/ERA_manipulations/Text_files/" + BoloName + "/liste_runs_bad_runs_removed.txt";

  BoloStr bolo=Read_liste_bolo(BoloFile,BoloName);
  vector<RunStr> listruns;
  listruns = Read_liste_runs(RunFile,bolo);

  // Boucle sur les runs
  for (UInt_t i_run=0; i_run<listruns.size(); i_run++) {
    RunStr run=listruns[i_run];

    string outputfilename=AmplDir+"/basicntp_"+run.Name+"_"+BoloName+".root";
    
    TChain* ch = new TChain("EdwTree","Edelweiss Raw Data Tree");
    ch->AddFile((TraceDir + BoloName + "_" + run.Name + "_simu_event_tree.root").c_str());

    cout << "Computing ampls for run "<<run.Name<<" ("<<ch->GetEntries()<<" evts)"<<endl;
    FillSimpleAmplTree(ch,bolo,outputfilename);
    ch->Delete("");
    
  }

  return 0;
}
