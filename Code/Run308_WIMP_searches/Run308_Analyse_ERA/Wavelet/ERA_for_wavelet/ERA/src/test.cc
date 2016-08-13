
#include "EdwUtils.h"
#include "EdwEvent.h"
#include "FitPulse.h"

int main(int argc, char* argv[]) {
   
  cout << "Hello world" << endl;
  // cout<<SHRT_MAX<<" "<<SHRT_MIN<<endl;

/*
// CPU test...
  for (int i=0; i<1000; i++) {
    if ((i/100)*100==i) cout <<i<<endl;

    EdwEvent* evt = new EdwEvent();
    
    for (int m=0; m<1000; m++) {
      int size=3000;
      double* tableau = new double[size];
      for (int n=0; n<size; n++) {
	tableau[n]=sin(n/(float)size);
      }
      delete[] tableau;
    }
    
    delete evt;
  }
  */
 
 // Test d'un fichier samba... (crosschekc kdata)
  TFile f("/Users/armengau/Edelweiss/Data_Edw/TestERANew/RunTest/FID825/Traces/nj20b002_000_FID825.root","READ");
  TTree* t=(TTree*)f.Get("EdwTree");
    cout << "entries:"<<t->GetEntries() << endl;
   EdwEvent* lEvt = new EdwEvent();
  t->SetBranchAddress("Event",&lEvt);
  int nbpulses_ok=0;
  for (int j=0; j<t->GetEntries(); j++) {
        t->GetEntry(j);
    if (lEvt->Pulse("ionisA FID825")!=NULL) nbpulses_ok++;
    if (j % 10000 == 0) cout << j << endl;
    }
    cout << "pulsesok:"<<nbpulses_ok<<endl;
 
  // Tests structures simples
  /*
  string bolo="ID404";
  string file="/Users/armengau/Edelweiss/Data_Edw/TestERANew/RunA/ID404/liste_runs.txt";
  Bolo b; b.Name="ID404";
  vector<Run> v = Read_liste_runs(file,b);
  cout << v.size()<<endl;
  for (int i=0; i<v.size(); i++) {
    Run toto = v[i];
    cout << toto.Name << endl;
    cout << toto.Vgar2 <<" "<<toto.VrelGar2<< endl;
    }*/

  // Evts - TEST BIT TRIGGER
    /*
  TFile f("/Users/armengau/Edelweiss/Data_Edw/TestERANew/RunA/ID404/Traces/la17e002_000_ID404.root","READ");
  TTree* t=(TTree*)f.Get("EdwTree");
  EdwEvent* lEvt = new EdwEvent();
  t->SetBranchAddress("Event",&lEvt);
  for (int j=0; j<50; j++) {
    t->GetEntry(j);
    int bit1=lEvt->TriggerBit(1); 
    int bit2=lEvt->TriggerBit(2);
    for (int numdetect=0; numdetect<=50; numdetect++)
      if (DetectTrig(bit1,bit2,numdetect)==1) cout << numdetect<<" ";
    cout <<endl;
    }*/

  return 0;
}

