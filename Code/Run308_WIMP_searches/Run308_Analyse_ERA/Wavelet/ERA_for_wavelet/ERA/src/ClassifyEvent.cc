
#include "EdwUtils.h"

int EventClass(string infilename, string outfilename, BoloStr bolo) {
  
  TFile infile(infilename.c_str(),"READ");
  TTree* tin=(TTree*)infile.Get(("eionntp_"+bolo.Name).c_str());
  TFile outfile(outfilename.c_str(),"RECREATE");
  TTree* tout = tin->CloneTree(0);
  Float_t ECol1, Chi2Col1, LdbCol1, ToffsetCol1, ECol2, Chi2Col2, LdbCol2, ToffsetCol2;
  Float_t EVet1, Chi2Vet1, LdbVet1, ToffsetVet1, EVet2, Chi2Vet2, LdbVet2, ToffsetVet2;
  Float_t EGar1, Chi2Gar1, LdbGar1, ToffsetGar1, EGar2, Chi2Gar2, LdbGar2, ToffsetGar2;
  if (bolo.Col1!="NONE") {
    tin->SetBranchAddress("ECol1",&ECol1); tin->SetBranchAddress("Chi2Col1",&Chi2Col1);
    tin->SetBranchAddress("LdbCol1",&LdbCol1); tin->SetBranchAddress("ToffsetCol1",&ToffsetCol1);
  }
  if (bolo.Col2!="NONE") {
    tin->SetBranchAddress("ECol2",&ECol2); tin->SetBranchAddress("Chi2Col2",&Chi2Col2);
    tin->SetBranchAddress("LdbCol2",&LdbCol2); tin->SetBranchAddress("ToffsetCol2",&ToffsetCol2);
  }
  if (bolo.Vet1!="NONE") {
    tin->SetBranchAddress("EVet1",&EVet1); tin->SetBranchAddress("Chi2Vet1",&Chi2Vet1);
    tin->SetBranchAddress("LdbVet1",&LdbVet1); tin->SetBranchAddress("ToffsetVet1",&ToffsetVet1);
  }
  if (bolo.Vet2!="NONE") {
    tin->SetBranchAddress("EVet2",&EVet2); tin->SetBranchAddress("Chi2Vet2",&Chi2Vet2);
    tin->SetBranchAddress("LdbVet2",&LdbVet2); tin->SetBranchAddress("ToffsetVet2",&ToffsetVet2);
  }
  if (bolo.Gar1!="NONE") {
    tin->SetBranchAddress("EGar1",&EGar1); tin->SetBranchAddress("Chi2Gar1",&Chi2Gar1);
    tin->SetBranchAddress("LdbGar1",&LdbGar1); tin->SetBranchAddress("ToffsetGar1",&ToffsetGar1);
  }
  if (bolo.Gar2!="NONE") {
    tin->SetBranchAddress("EGar2",&EGar2); tin->SetBranchAddress("Chi2Gar2",&Chi2Gar2);
    tin->SetBranchAddress("LdbGar2",&LdbGar2); tin->SetBranchAddress("ToffsetGar2",&ToffsetGar2);
  }
  Float_t Eion=0, Tion=0; Int_t IonFlags[6]={0,0,0,0,0,0}; Int_t EvtFlag=0;
  tout->Branch("Eion",&Eion,"Eion/F");
  tout->Branch("Tion",&Tion,"Tion/F");
  tout->Branch("IonFlags",IonFlags,"IonFlags[6]/I");
  tout->Branch("EventFlag",&EvtFlag,"EventFlag/I");

  // variables de coupure
  int deltabin_col=3;
  int deltabin_vet=5;
  float sigmacut_efid=2;
  float vetopercentcut=5;
  float sigmacut_delta=2;
  float sigmacut_vet=2;
  float atten=2;
  if (bolo.Col1=="NONE" || bolo.Col2=="NONE") cerr << "Il faut que les cols existent toutes les 2.."<<endl;

  for (UInt_t ievt=0; ievt<tin->GetEntries(); ievt++) {
    tin->GetEntry(ievt);

    //******************************
    // Event Classification algo ICI
    //******************************
    Eion=0; Tion=0; EvtFlag=1;
    for (int i=0; i<6; i++) IonFlags[i]=0;

    // Test bulk
    //***************
    if (fabs(ToffsetCol1-ToffsetCol2) <= deltabin_col) {
      float invdenom=1./(LdbCol1*LdbCol1+LdbCol2*LdbCol2);
      float ldbfid=LdbCol1*LdbCol2*sqrt(invdenom);
      float w1=LdbCol2*LdbCol2*invdenom;
      float w2=LdbCol1*LdbCol1*invdenom;
      float Efid = w1*ECol1+w2*ECol2;
      if (Efid > sigmacut_efid*ldbfid) {
	float BinFid=(ToffsetCol1+ToffsetCol2)/2;
	float truc=0.01*vetopercentcut*Efid;
	float DeltaE = fabs(ECol1-ECol2);
	float Ldb_DeltaE = sqrt(LdbCol1*LdbCol1+LdbCol2*LdbCol2);
	if (DeltaE < sigmacut_delta*Ldb_DeltaE || DeltaE < truc) {
	  // collectrodes are ok, now check other channels
	  bool ok=1;
	  if (bolo.Vet1!="NONE" && LdbVet1!=0) {
            float ecut = LdbVet1*sigmacut_vet;
            if (fabs(BinFid-ToffsetVet1) <= deltabin_vet) ecut /= atten;
            if (EVet1 > ecut && EVet1 > truc) ok=0;
          }
	  if (bolo.Vet2!="NONE" && LdbVet2!=0) {
            float ecut = LdbVet2*sigmacut_vet;
            if (fabs(BinFid-ToffsetVet2) <= deltabin_vet) ecut /= atten;
            if (EVet2 > ecut && EVet2 > truc) ok=0;
          }
	  if (bolo.Gar1!="NONE" && LdbGar1!=0) {
            float ecut = LdbGar1*sigmacut_vet;
            if (fabs(BinFid-ToffsetGar1) <= deltabin_vet) ecut /= atten;
            if (EGar1 > ecut && EGar1 > truc) ok=0;
          }
	  if (bolo.Gar2!="NONE" && LdbGar2!=0) {
            float ecut = LdbGar2*sigmacut_vet;
            if (fabs(BinFid-ToffsetGar2) <= deltabin_vet) ecut /= atten;
            if (EGar2 > ecut && EGar2 > truc) ok=0;
          }
	  if (ok) {
	    // Evt fiduciel. (tfid,Efid) a deja ete calculee..
	    EvtFlag=2;
	    Tion=BinFid;
	    Eion=Efid;
	    IonFlags[0]=1; IonFlags[1]=1;
	  }
	}
      }
    }
    // Fin du test fiduciel
    
    // Test surface haut
    //************
    if (EvtFlag <= 1 && bolo.Vet1!="NONE" && LdbVet1!=0) {
      if (fabs(ToffsetCol1-ToffsetVet1) <= deltabin_col) { // surface haut
	float BinFid = (ToffsetCol1+ToffsetVet1)/2;
	if (ECol1 >= sigmacut_efid*LdbCol1 && EVet1 >= sigmacut_efid*LdbVet1) {
	  bool ok=1;
	  float truc=0.01*vetopercentcut*max(ECol1,EVet1);
	  if (bolo.Vet2!="NONE" && LdbVet2!=0) {
            float ecut = LdbVet2*sigmacut_vet;
            if (fabs(BinFid-ToffsetVet2) <= deltabin_vet) ecut /= atten;
            if (EVet2 > ecut && EVet2 > truc) ok=0;
          }
	  if (bolo.Gar1!="NONE" && LdbGar1!=0) {
            float ecut = LdbGar1*sigmacut_vet;
            if (fabs(BinFid-ToffsetGar1) <= deltabin_vet) ecut /= atten;
            if (EGar1 > ecut && EGar1 > truc) ok=0;
          }
	  if (bolo.Gar2!="NONE" && LdbGar2!=0) {
            float ecut = LdbGar2*sigmacut_vet;
            if (fabs(BinFid-ToffsetGar2) <= deltabin_vet) ecut /= atten;
            if (EGar2 > ecut && EGar2 > truc) ok=0;
          }
	  if (ok) {
	    EvtFlag=3;
	    Tion = BinFid;
	    if (fabs(BinFid-ToffsetCol2) <= deltabin_vet && ECol2 > truc && 
		ECol2 > LdbCol2*sigmacut_vet/atten) { // Evt surface 3 electrodes
	      // La on fait encore un truc simple: on teste pas Ec1=Ec2+Ev1...
	      // et pour l'instant on fait une eion pas bien ponderee
	      Eion=(ECol1+ECol2+EVet1)/2.;
	      IonFlags[0]=1; IonFlags[1]=1; IonFlags[2]=1;
	    } else { // Evt surface 2 electrodes
	      float invdenom=1./(LdbCol1*LdbCol1+LdbVet1*LdbVet1);
	      float w1=LdbVet1*LdbVet1*invdenom;
	      float w2=LdbCol1*LdbCol1*invdenom;
	      Eion=w1*ECol1+w2*EVet1;
	      IonFlags[0]=1; IonFlags[2]=1;
	    }
	  }
	}
      }
    }
    // Fin test surface haut

    // Test surface bas
    //************
    if (EvtFlag <= 1 && bolo.Vet2!="NONE" && LdbVet2!=0) {
      if (fabs(ToffsetCol2-ToffsetVet2) <= deltabin_col) { // surface haut
	float BinFid = (ToffsetCol2+ToffsetVet2)/2;
	if (ECol2 >= sigmacut_efid*LdbCol2 && EVet2 >= sigmacut_efid*LdbVet2) {
	  bool ok=1;
	  float truc=0.01*vetopercentcut*max(ECol2,EVet2);
	  if (bolo.Vet1!="NONE" && LdbVet1!=0) {
            float ecut = LdbVet1*sigmacut_vet;
            if (fabs(BinFid-ToffsetVet1) <= deltabin_vet) ecut /= atten;
            if (EVet1 > ecut && EVet1 > truc) ok=0;
          }
	  if (bolo.Gar1!="NONE" && LdbGar1!=0) {
            float ecut = LdbGar1*sigmacut_vet;
            if (fabs(BinFid-ToffsetGar1) <= deltabin_vet) ecut /= atten;
            if (EGar1 > ecut && EGar1 > truc) ok=0;
          }
	  if (bolo.Gar2!="NONE" && LdbGar2!=0) {
            float ecut = LdbGar2*sigmacut_vet;
            if (fabs(BinFid-ToffsetGar2) <= deltabin_vet) ecut /= atten;
            if (EGar2 > ecut && EGar2 > truc) ok=0;
          }
	  if (ok) {
	    EvtFlag=3;
	    Tion = BinFid;
	    if (fabs(BinFid-ToffsetCol1) <= deltabin_vet && ECol1 > truc && 
		ECol1 > LdbCol1*sigmacut_vet/atten) { // Evt surface 3 electrodes
	      Eion=(ECol2+ECol1+EVet2)/2.;
	      IonFlags[0]=1; IonFlags[1]=1; IonFlags[3]=1;
	    } else { // Evt surface 2 electrodes
	      float invdenom=1./(LdbCol2*LdbCol2+LdbVet2*LdbVet2);
	      float w1=LdbVet2*LdbVet2*invdenom;
	      float w2=LdbCol2*LdbCol2*invdenom;
	      Eion=w1*ECol2+w2*EVet2;
	      IonFlags[1]=1; IonFlags[3]=1;
	    }
	  }
	}
      }
    }
    // Fin test surface bas

    // TODO : test pure garde..

    // Test generique sur toutes les voies ion
    //****************
    if (EvtFlag <= 1) {
      if (bolo.Col1!="NONE" && ECol1>1.5*sigmacut_efid*LdbCol1 && LdbCol1!=0) { // malus 1.5...
        EvtFlag=5;
	Tion=ToffsetCol1;
        Eion+=(ECol1/2.);
        IonFlags[0]=1;
        // ALGO A AMELIORER... il faut tester les bins...
      }
      if (bolo.Col2!="NONE" && ECol2>1.5*sigmacut_efid*LdbCol2 && LdbCol2!=0) {
        EvtFlag=5;
	Tion=ToffsetCol2;
        Eion+=(ECol2/2.);
        IonFlags[1]=1;
      }
      if (bolo.Vet1!="NONE" && EVet1>1.5*sigmacut_efid*LdbVet1 && LdbVet1!=0) {
        EvtFlag=5;
	Tion=ToffsetVet1;
        Eion+=(EVet1/2.);
        IonFlags[2]=1;
      }
      if (bolo.Vet2!="NONE" && EVet2>1.5*sigmacut_efid*LdbVet2 && LdbVet2!=0) {
        EvtFlag=5;
	Tion=ToffsetVet2;
        Eion+=(EVet2/2.);
        IonFlags[3]=1;
      }
      if (bolo.Gar1!="NONE" && EGar1>1.5*sigmacut_efid*LdbGar1 && LdbGar1!=0) {
        EvtFlag=5;
	Tion=ToffsetGar1;
        Eion+=(EGar1/2.);
        IonFlags[4]=1;
      }
      if (bolo.Gar2!="NONE" && EGar2>1.5*sigmacut_efid*LdbGar2 && LdbGar2!=0) {
        EvtFlag=5;
	Tion=ToffsetGar2;
        Eion+=(EGar2/2.);
        IonFlags[5]=1;
      }
    }
    // Fin du test generique

    tout->Fill();
  }
  outfile.cd(); tout->Write("",TObject::kOverwrite);
  outfile.Close();
  infile.Close();

  return 0;
}


int main(int argc, char* argv[]) {

  // Load parametres.
  if (argc!=2) {cerr <<"Wrong nb of arguments" << endl; return 1;}
  string paramfile=string(argv[1]);
  string AnaDir=GetParam(paramfile,"AnaDir");
  string BoloName=GetParam(paramfile,"BoloName");
  if (AnaDir=="NONE" || BoloName=="NONE") {
    cerr << "wrong param file"<<endl; return 1;
  }
  string BoloFile=GetParam(paramfile,"BoloFile");
  if (BoloFile=="NONE") BoloFile=AnaDir+"/liste_bolos.txt";
  string RunFile=GetParam(paramfile,"RunFile");
  if (RunFile=="NONE") RunFile=AnaDir+"/"+BoloName+"/liste_runs.txt";
  Bool_t overwrite=1;
  string dum=GetParam(paramfile,"OverWrite");
  if (dum!="NONE") overwrite=(Bool_t)(atoi(dum.c_str()));

  BoloStr bolo=Read_liste_bolo(BoloFile,BoloName); 
  vector<RunStr> listruns = Read_liste_runs(RunFile,bolo);
  string AmplDir=AnaDir+"/"+BoloName+"/Amplitudes/";

  // Boucle sur les runs
  for (UInt_t i_run=0; i_run<listruns.size(); i_run++) {
    RunStr run=listruns[i_run];
    string infile=AmplDir+"/eion_"+run.Name+"_"+BoloName+".root";
    bool process_run=1;
    if (!overwrite) { // on teste si le fichier eion est deja complet.
      TFile f_test(infile.c_str(),"READ");
      TTree* t_test=(TTree*)f_test.Get(("eionntp_"+BoloName).c_str());
      if (t_test->GetBranch("EventFlag")) process_run=0;
      f_test.Close();
    }
    if (!process_run) continue;
    cout << "Classification for run "<<run.Name<<endl;
    string outfile=AmplDir+"/eionfull_"+run.Name+"_"+BoloName+".root";
    EventClass(infile,outfile,bolo);
    system(("mv -f "+outfile+" "+infile).c_str());
  }

}
