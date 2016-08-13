
#include "EdwUtils.h"
#include "EdwEvent.h"
#include "FitPulse.h"
#include "EdwTemplate.h"
#include "NoiseSpectrum.h"

// Voie par voie

// - Pour wienerprocess, contrairement a simpleprocess, on est a un run FIXE!
// - voie par voie: fait un tmpwiener_run_voietype_bolo.root

int FillWienerAmplTree(TChain* chain, BoloStr bolo, UInt_t pchannel, TemplateStr template_voie, TTree* spectrumtree, string outfilename, TTree* synctree=NULL, Int_t syncoffsetrange=0) {

  // 1) Prepare le ntuple
  TFile* file=new TFile(outfilename.c_str(),"RECREATE");
  string ntpname="wienerntp_";
  ntpname+=(bolo.Name+"_"+bolo.ChTypes[pchannel]);
  TTree* ntp = new TTree(ntpname.c_str(),"Edelweiss Wiener Data Tree");
  ntp->SetDirectory(file);
  Float_t wienerampl=0, wienerzeroampl=0, wienerchi2=0, wienertime=0, wienersigmaampl=0;
  Float_t OriginalwienerA=0, OriginalwienerB=0, OriginalHeat=0;
  Float_t wienerfreeampl=0, wienerfreechi2=0, wienersynctime=0; // variables optionnelles

  ntp->Branch("OriWienerAmplA",&OriginalwienerA,"OriWienerAmplA/F");
  ntp->Branch("OriWienerAmplB",&OriginalwienerB,"OriWienerAmplBl/F");
  ntp->Branch("OriHeat",&OriginalHeat,"OriHeat/F");

  ntp->Branch("WienerAmpl",&wienerampl,"WienerAmpl/F");
  ntp->Branch("WienerChi2",&wienerchi2,"WienerChi2/F");
  ntp->Branch("WienerZeroAmpl",&wienerzeroampl,"WienerZeroAmpl/F");
  ntp->Branch("WienerTime",&wienertime,"WienerTime/F");
  ntp->Branch("WienerSigmaAmpl",&wienersigmaampl,"WienerSigmaAmpl/F");
  // Les temps: depend de l'option choisie
  if (synctree!=NULL) {
    synctree->SetBranchAddress("WienerTime",&wienersynctime);
    ntp->Branch("WienerFreeAmpl",&wienerfreeampl,"WienerFreeAmpl/F");
    ntp->Branch("WienerFreeChi2",&wienerfreechi2,"WienerFreeChi2/F");
  }
  // 1bis) Initialise spectre et template
  NoiseSpectrum* spectrum = new NoiseSpectrum();
  spectrumtree->SetBranchAddress("Spectrum",&spectrum);
  EdwTemplate* Template = new EdwTemplate();

  // 2) Boucle sur les evts
  EdwEvent* evt = new EdwEvent();
  chain->SetBranchAddress("WienerAmplA",&OriginalwienerA);
  chain->SetBranchAddress("WienerAmplB",&OriginalwienerB);
  chain->SetBranchAddress("True_heat",&OriginalHeat);

  chain->SetBranchAddress("Event",&evt);
  EdwPulse* edwpulse = new EdwPulse();
  for (UInt_t ievt=0; ievt<chain->GetEntries(); ievt++) {
    chain->GetEntry(ievt);
    edwpulse=evt->Pulse(bolo.ChNames[pchannel]);
    if (edwpulse!=NULL) {
      FitPulse fp(edwpulse);
      // a) Mise a jour template et spectre
      if (spectrum->Channel() != bolo.ChNames[pchannel] || spectrum->StartValidity() > evt->DateSec() ||
	  spectrum->EndValidity() < evt->DateSec()) {
	bool ok=0;
	for (UInt_t i=0; i<spectrumtree->GetEntries(); i++) {
	  spectrumtree->GetEntry(i);
	  if (spectrum->Channel() == bolo.ChNames[pchannel] && spectrum->StartValidity() <= evt->DateSec() && spectrum->EndValidity() >= evt->DateSec()) { 
	    ok=1; break;
	  }
	}
	if (!ok) {cerr << ievt<<" :Not found a good spectrum!!" << endl; exit(-1);}
      }
      if (Template->Channel() != bolo.ChNames[pchannel] || Template->TraceLength() != edwpulse->TraceLength() ||
	  Template->Pretrigger() != edwpulse->Pretrigger()) {
	Int_t offsetmin=fp.FitTimeMin();
	Int_t offsetmax=fp.FitTimeMax();
	Template->Initialize(template_voie.Channel,template_voie.Type,template_voie.Tinf,template_voie.Tsup,template_voie.Coefficients,edwpulse->TraceLength(),edwpulse->Pretrigger(),offsetmin,offsetmax,spectrum->Spectrum());
      }
      if (spectrum->Channel() != bolo.ChNames[pchannel] || spectrum->StartValidity() > evt->DateSec() ||
	  spectrum->EndValidity() < evt->DateSec()) {
	Template->ComputeKernels(spectrum->Spectrum());
      }
      // b) Et maintenant on calcule wiener (offsets en fonction des options)
      fp.BasicPreprocess();
      wienerzeroampl=fp.Sign()*Template->Sign()*(fp.ComputeWienerAmpl(Template->GetOffsetFFT(0),spectrum->Spectrum()))[0];
      vector<Float_t> freetime_result=fp.WienerLoop(Template,spectrum->Spectrum());
      if (synctree==NULL) {
	// i) Aucune synchro : tout est deja fait
	wienerampl=freetime_result[0];
	wienertime=freetime_result[1];
	wienerchi2=freetime_result[2];
	wienersigmaampl=freetime_result[3];
      } else if (synctree!=NULL && syncoffsetrange!=0) {
	// ii) On donne un range de synchronisation
	// calcul range offset (prise en compte modulation de la voie geree)
	wienerfreeampl=freetime_result[0];
	//	wienertime=freetime_result[1]; // A mieux gerer.. plus tard
	wienerfreechi2=freetime_result[2];
	synctree->GetEntry(ievt);
	Float_t scale_fact = (Float_t)fp.ModulationLength();
	if ( (bolo.ChNames[pchannel]).find("slow")!=string::npos ) scale_fact/=100; // gestion voie slow...
	Int_t offsetmin=(Int_t)(wienersynctime/scale_fact-syncoffsetrange/2.0);
	Int_t offsetmax=(Int_t)(wienersynctime/scale_fact+syncoffsetrange/2);
	vector<Float_t> synctime_result=fp.WienerLoop(Template,spectrum->Spectrum(),offsetmin,offsetmax);
	wienerampl=synctime_result[0];
	wienertime=synctime_result[1];
	wienerchi2=synctime_result[2];
	wienersigmaampl=synctime_result[3];
      } else if (synctree!=NULL && syncoffsetrange==0) {
	// iii) On fixe le temps de fit d'apres la synchro
	wienerfreeampl=freetime_result[0];
	wienertime=freetime_result[1];
	wienerfreechi2=freetime_result[2];
	wienersigmaampl=freetime_result[3];
	synctree->GetEntry(ievt);
	Float_t scale_fact = (Float_t)fp.ModulationLength();
	if ( (bolo.ChNames[pchannel]).find("slow")!=string::npos ) scale_fact/=100; // gestion voie slow...
	Int_t offset=(Int_t)(wienersynctime/scale_fact);
	if (offset<Template->MinOffset() || offset>Template->MaxOffset()) {
	  wienerampl=0;
	  wienerchi2=0;
	} else {
	  vector<Float_t> synctime_result=fp.ComputeWienerAmpl(Template->GetOffsetFFT(offset),spectrum->Spectrum(),1);
	  wienerampl=fp.Sign()*Template->Sign()*synctime_result[0];
	  wienerchi2=synctime_result[2];
	}
      }
    } else {
      wienerampl=0; wienertime=0; wienerchi2=0; wienerzeroampl=0;
      wienerfreeampl=0; wienerfreechi2=0; wienersigmaampl=0;
    }
    file->cd();
    ntp->Fill();
    evt->Clear(); // Indispensable contre fuite de memoire...
    if (ievt % 5000 == 0) cout << "Evt " << ievt << endl;
  }
  file->cd();
  ntp->Write("",TObject::kOverwrite);
  file->Close();
  delete spectrum;
  delete Template;
  return 0;
}


int main(int argc, char *argv[]) {

  // Load parametres...
  string AnaDir="/home/irfulx204/mnt/tmain/Desktop/Run308_Analyse_ERA";
  string BoloName=argv[1];

  string TraceDir=AnaDir+"/Event_files/" + BoloName + "/Event_simu/";
  string AmplDir=AnaDir+ "/Amp_files/" + BoloName + "/Amp_simu/";
  string SpectraDir=AnaDir+"/Spectra_files/"+ BoloName + "/Spectra_simu/";

  string BoloFile=AnaDir+"/ERA_manipulations/Text_files/liste_bolos.txt";
  string RunFile=AnaDir+"/ERA_manipulations/Text_files/" + BoloName + "/liste_runs_bad_runs_removed.txt";

  BoloStr bolo=Read_liste_bolo(BoloFile,BoloName);
  vector<RunStr> listruns;
  listruns = Read_liste_runs(RunFile,bolo);

  string PeriodFile=AnaDir+"/ERA_manipulations/Text_files/" + BoloName + "/liste_periods.txt";
  string TemplateFile=AnaDir+"/ERA_manipulations/Text_files/" + BoloName + "/liste_templates.txt";

  string ChannelSelection="Chaleur";
  string SynchroMode="NONE";
  string tmplt_type="NONE";

  // "FromHeat" ou "FromIonFid" [ou "NONE"]
  Bool_t overwrite=1;

  // Boucle sur les runs
  for (UInt_t i_run=0; i_run<listruns.size(); i_run++) {

        RunStr run = listruns[i_run];
        cout << "--> Starting Run " << run.Name << endl;
        TFile* f_spectra = new TFile((SpectraDir+"spectra_"+run.Name+"_"+BoloName+".root").c_str(),"READ");
        TTree* spectrumtree = (TTree*)f_spectra->Get("SpectrumTree");
        if (!f_spectra || !spectrumtree) {cerr << "No spectrum!" <<endl; exit(-1);}
        TChain* ch = new TChain("EdwTree","Edelweiss Raw Data Tree");
        ch->AddFile((TraceDir + BoloName + "_" + run.Name + "_simu_event_tree.root").c_str());
        cout << ch->GetEntries() << " evts to process.." <<endl;

        // Boucle sur les voies
        for (UInt_t p=0; p<bolo.ChTypes.size(); p++) {
            if (ChannelSelection=="Chaleur" && (bolo.ChTypes[p]).find("Chal")==string::npos) continue;
            if (ChannelSelection=="Ionisation" && (bolo.ChTypes[p]).find("Chal")!=string::npos) continue;
            if ( (ChannelSelection=="Col1" || ChannelSelection=="Col2" || ChannelSelection=="Vet1" ||
              ChannelSelection=="Vet2" || ChannelSelection=="Gar1" || ChannelSelection=="Gar2" ||
              ChannelSelection=="Chal1" || ChannelSelection=="Chal2") 
             && bolo.ChTypes[p]!=ChannelSelection) continue;
            TemplateStr tmplt=Read_list_template(TemplateFile,bolo.ChNames[p],tmplt_type);
            // Choix actuel: meme nom de fichier quel que soit le process..
            string fileprefix="tmpwiener";
            if (tmplt_type=="ChalNTD") fileprefix="amplntd";
            string outfile=AmplDir+"/"+fileprefix+"_"+run.Name+"_"+bolo.ChTypes[p]+"_"+BoloName+".root";
            if (!overwrite && file_exists(outfile)) continue; 
            cout << "Calcul amplitudes pour la voie "<<bolo.ChNames[p]<<endl;
            if (SynchroMode=="NONE") {
              FillWienerAmplTree(ch,bolo,p,tmplt,spectrumtree,outfile) ;
            } else if (SynchroMode=="FromHeat") {
              // params en DUR pour l'instant : synchro from Chal1
              string voiesync="Chal1";
              Int_t syncoffsetrange=400; // a verifier que c'est bon...
              // Cas d'une voie slow : on divise a la main par 100...
              if ((bolo.ChNames[p]).find("slow")!=string::npos) syncoffsetrange/=100;
              string synchrofile=AmplDir+"tmpwiener_"+run.Name+"_"+voiesync+"_"+BoloName+".root";
              if (bolo.ChTypes[p]==voiesync || bolo.ChTypes[p]=="Chal2") continue; // bricolo...
              TFile f(synchrofile.c_str(),"READ");
              TTree* tsync = (TTree*)f.Get(("wienerntp_"+BoloName+"_"+voiesync).c_str());
              FillWienerAmplTree(ch,bolo,p,tmplt,spectrumtree,outfile,tsync,syncoffsetrange);
              f.Close();
            } else if (SynchroMode=="FromIonFid") {
              //synchro from Col1 for heats or Vet1, Col2 for Vet2
              string voiesync = (bolo.ChTypes[p]=="Vet2") ? "Col2" : "Col1";
              string synchrofile=AmplDir+"tmpwiener_"+run.Name+"_"+voiesync+"_"+BoloName+".root";
              if (bolo.ChTypes[p]==voiesync || bolo.ChTypes[p]=="Col2") continue; // bricolo a ce stade
              TFile f(synchrofile.c_str(),"READ");
              TTree* tsync = (TTree*)f.Get(("wienerntp_"+BoloName+"_"+voiesync).c_str());
              FillWienerAmplTree(ch,bolo,p,tmplt,spectrumtree,outfile,tsync);
              f.Close();
            } else {cerr<< "Wrong param: SynchroMode"<<endl;exit(-1);}
        }

        f_spectra->Close();
  }
        
  return 0;
}
