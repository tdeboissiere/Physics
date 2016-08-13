
#include "EdwUtils.h"
#include "NoiseSpectrum.h"
#include "EdwEvent.h"
#include "FitPulse.h"

ULong64_t GetFirstEventAfter(TChain* aEvtChain, ULong64_t aT0) {
  EdwEvent* lEvt=new EdwEvent();
  aEvtChain->SetBranchAddress("Event",&lEvt);

  bool ok=0; ULong64_t b=0, lTime=0;
  while (!ok && b<(ULong64_t)aEvtChain->GetEntries()) {
    aEvtChain->GetEntry(b);
    lTime = lEvt->DateSec();
    if (lTime >= aT0) ok=1;
    lEvt->Clear();
    b+=1;
  }
  if (b==(ULong64_t)aEvtChain->GetEntries()) cerr << "GetFirstEventAfter: Pbl, evt not found.." <<endl;

  delete lEvt;
  return b-1;
}

ULong64_t GetLastEventBefore(TChain* aEvtChain, ULong64_t aT0) {
  EdwEvent* lEvt=new EdwEvent();
  aEvtChain->SetBranchAddress("Event",&lEvt);

  bool ok=0; ULong64_t a=aEvtChain->GetEntries(), lTime=0;
  while (!ok && a>0) {
    a-=1;
    aEvtChain->GetEntry(a);
    lTime = lEvt->DateSec();
    if (lTime <= aT0) ok=1;
    lEvt->Clear();
  }
  if (a==0) cerr << "GetLastEventBefore: Pbl, evt not found.." <<endl;

  delete lEvt;
  return a;
}

vector<Float_t> StackNoiseSpectra(TChain* aChain, string aChannel, ULong64_t aFirstEvt, ULong64_t aLastEvt, UInt_t& NbTracesForNoiseSpectrum, Float_t& Sampling_ms, Int_t aCriterium=FP_STRICT) {
  // aFirstevt et aLastevt sont calcules ailleurs, en fonction de la subdivision en temps
  EdwEvent* lEvt=new EdwEvent();
  // todo= tester que aChain a bien cette branche..
  aChain->SetBranchAddress("Event",&lEvt);
  NbTracesForNoiseSpectrum=0;
  Sampling_ms=0;
  vector<Float_t> lNoiseSpectrum;

  for (ULong64_t ievt=aFirstEvt; ievt<=aLastEvt; ievt++) {
    aChain->GetEntry(ievt);
    EdwPulse* lEdwPulse = lEvt->Pulse(aChannel);
    if (lEdwPulse != NULL) {
      FitPulse lPulse(lEdwPulse);
      UInt_t nbglitches=0;
      // La il faut faire tres attention a la baseline...
      if (aCriterium!=FP_NONE) {
	nbglitches=lPulse.FindGlitches();
	lPulse.FindPeaks(aCriterium,1);
      }
      // sinon normalement peakbins n'est pas rempli! (a bien verifier)
      if (lPulse.NbPeaks() == 0 && nbglitches<5) {
	// !! nouvelle recette sur recherche de gliches: param libre..
	if (lPulse.IsHeat()) lPulse.RemoveLinearBase();
	lPulse.ComputeTraceFFT();
	vector<Float_t> lSpec = lPulse.ProcessedTraceFFT();
	VectMultiply(lSpec,lPulse.ProcessedTraceFFT());
	NbTracesForNoiseSpectrum+=1;
	Sampling_ms=lPulse.Sampling_ms();
	if (lNoiseSpectrum.size() == 0) lNoiseSpectrum = lSpec;
	else VectAdd(lNoiseSpectrum,lSpec);
      }
    }
    lEvt->Clear();
  }
  VectScale(lNoiseSpectrum,((Float_t)1)/(Float_t)NbTracesForNoiseSpectrum);
  return lNoiseSpectrum;
}

vector<Float_t> ComputeNoiseSpectrum(TChain* aChain, string aChannel, ULong64_t aFirstEvt, ULong64_t aLastEvt, Float_t& Sampling_ms) {
  UInt_t NbTracesForNoiseSpectrum=0;
  Sampling_ms=0;
  vector<Float_t> spectrum = StackNoiseSpectra(aChain, aChannel, aFirstEvt,aLastEvt,NbTracesForNoiseSpectrum,Sampling_ms,FP_STRICT);

  if (NbTracesForNoiseSpectrum < 10) {
    cerr << "ComputeNoiseSpectrum : New loop with loose FindPeaks cut" << endl;
    spectrum = StackNoiseSpectra(aChain, aChannel, aFirstEvt,aLastEvt,NbTracesForNoiseSpectrum,Sampling_ms,FP_WEAK);
  }
  if (NbTracesForNoiseSpectrum < 10) {
    cerr << "ComputeNoiseSpectrum : New loop with eXtraWEAK FindPeaks cut."<<endl;
    spectrum = StackNoiseSpectra(aChain, aChannel, aFirstEvt,aLastEvt,NbTracesForNoiseSpectrum,Sampling_ms,FP_XWEAK);
  }
  if (NbTracesForNoiseSpectrum < 5) { // Ultimate round: take all traces!
    cerr << "ComputeNoiseSpectrum :  New loop, desperate case, taking all traces.." << endl;
    spectrum = StackNoiseSpectra(aChain, aChannel, aFirstEvt,aLastEvt,NbTracesForNoiseSpectrum,Sampling_ms,FP_NONE);
  }

  cerr << "ComputeNoiseSpectrum: "<<NbTracesForNoiseSpectrum<<" traces used for "<<aChannel<<endl;
  return spectrum;
}

int main(int argc, char *argv[]) {

  // Load parametres...
  string AnaDir="/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA";
  string BoloName=argv[1];

  string TraceDir=AnaDir+"/Event_files/" + BoloName + "/Event_simu/";
  string NoiseDir=AnaDir+"/Spectra_files/"+ BoloName + "/Spectra_simu/";

  string BoloFile=AnaDir+"/ERA_manipulations/Text_files/liste_bolos.txt";
  string RunFile=AnaDir+"/ERA_manipulations/Text_files/" + BoloName + "/liste_runs_bad_runs_removed.txt";

  BoloStr bolo=Read_liste_bolo(BoloFile,BoloName);
  vector<RunStr> listruns;
  listruns = Read_liste_runs(RunFile,bolo);

  string PeriodFile=AnaDir+"/ERA_manipulations/Text_files/" + BoloName + "/liste_periods.txt";

  // Boucle sur les runs
  for (UInt_t i_run=0; i_run<listruns.size(); i_run++) {
    RunStr run=listruns[i_run];
    string outfilename=NoiseDir+"spectra_"+run.Name+"_"+BoloName+".root";
    vector<PeriodStr> periods=Read_list_periods(PeriodFile,run.Name);

    TChain* ch = new TChain("EdwTree","Edelweiss Raw Data Tree");
    ch->AddFile((TraceDir + BoloName + "_" + run.Name + "_simu_event_tree.root").c_str());

    TFile* RootFile = new TFile(outfilename.c_str(),"RECREATE");
    if (RootFile->IsZombie()) cerr << "Error creating rootfile..." << endl;
    RootFile->cd();
    NoiseSpectrum* noisespec = new NoiseSpectrum();
    TTree SpectrumTree("SpectrumTree","Edelweiss Noise Spectrum Tree");
    SpectrumTree.SetDirectory(RootFile);
    SpectrumTree.Branch("Spectrum","NoiseSpectrum",&noisespec);
    cout << "Computing noise spectra for run "<<run.Name<<" ("<<ch->GetEntries()<<" evts), "<<periods.size()<<" periods"<<endl;
    // Boucle sur les periodes
    for (UInt_t p=0; p<periods.size(); p++) {
      ULong64_t ievt_inf=GetFirstEventAfter(ch,periods[p].Tinf);
      ULong64_t ievt_sup=GetLastEventBefore(ch,periods[p].Tsup);
      // Boucle sur les voies
      for (UInt_t k=0; k<bolo.ChNames.size(); k++) {
	  noisespec->SetChannel(bolo.ChNames[k]);
	  noisespec->SetStartValidity(periods[p].Tinf);
	  noisespec->SetEndValidity(periods[p].Tsup);
	  Float_t sampling_ms=0;
	  vector<Float_t> bruit = ComputeNoiseSpectrum(ch,bolo.ChNames[k],ievt_inf,ievt_sup,sampling_ms);
	  noisespec->SetSpectrum(bruit);
	  noisespec->SetSampling_ms(sampling_ms);
	  RootFile->cd();
	  SpectrumTree.Fill();
      }
    }
    RootFile->cd(); SpectrumTree.Write("",TObject::kOverwrite);
    RootFile->Close();
    ch->Delete("");
  }

  return 0;
}
