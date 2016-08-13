
#include "EdwUtils.h"

#include "Mixfft.c"

vector<Float_t> EdwRealFFT(vector<Float_t> aVect) {

  Int_t lSize = aVect.size();
  Double_t *in_re = new Double_t[lSize];
  Double_t *in_im = new Double_t[lSize];
  for (Int_t i=0; i<lSize; i++) {
    in_re[i]=(Double_t)(aVect[i]);
    in_im[i]=0;
  }
  Double_t *out_re = new Double_t[lSize];
  Double_t *out_im = new Double_t[lSize];
  fft(lSize,in_re,in_im,out_re,out_im);
  
  // Output in "Root" semi-complex format
  // r0 ... rN/2 i(N/2-1) ... i1
  //  vector<Float_t> lOutput;
  //  for (Int_t i=0; i<= lSize/2; i++) lOutput.push_back((Float_t)(out_re[i]));
  //  for (Int_t i=((lSize/2)-1); i>0; i--) lOutput.push_back((Float_t)(out_im[i]));
  vector<Float_t> lOutput(lSize,0); // eviter les push_back
  for (Int_t i=0; i<= lSize/2; i++) lOutput[i]=(Float_t)(out_re[i]);
  for (Int_t i=((lSize/2)-1); i>0; i--) lOutput[lSize-i]=(Float_t)(out_im[i]);

  delete[] in_re;
  delete[] in_im;
  delete[] out_re;
  delete[] out_im;

  return lOutput;
}

//****************************************************
// Fonctions de manipulation de vecteurs
//****************************************************

Float_t VectMean(const vector<Float_t> & aVect) {
  UInt_t n = aVect.size();
  if (!n) cerr << "vectmean: vector contains nothing." << endl;
  Float_t sum = 0;
  for (UInt_t i=0; i<n; i++) sum+=(aVect[i]);
  return sum/(Float_t)n;
}

Float_t VectMin(const vector<Float_t> & aVect) {
  if (!aVect.size()) cerr << "vectmin: vector contains nothing!" << endl;
  Float_t lMin=aVect[0];
  for (UInt_t i=1; i<aVect.size(); i++) {
    if (aVect[i]<lMin) lMin = aVect[i];
  }
  return lMin;
}

Float_t VectMax(const vector<Float_t> & aVect) {
  if (!aVect.size()) cerr << "vectmax: vector contains nothing!" << endl;
  Float_t lMax=aVect[0];
  for (UInt_t i=1; i<aVect.size(); i++) {
    if (aVect[i]>lMax) lMax = aVect[i];
  }
  return lMax;
}

Float_t VectMedian(const vector<Float_t> & aVect) {
  // For this routine we prefer to use directly the implementation by ROOT.
  UInt_t n = aVect.size();
  if (!n) cerr << "vectmedian: vector contains nothing." << endl;
  Float_t* tab = new Float_t[n];
  for (UInt_t i=0; i<n; i++) tab[i]=aVect[i];
  Float_t lMed = (Float_t)TMath::Median(n,tab);
  delete tab;
  return lMed;
}

Float_t VectRMS(const vector<Float_t> & aVect) {
  UInt_t n = aVect.size();
  if (!n) cerr << "vectrms: vector contains nothing." << endl;
  Float_t sum = 0, sum2 = 0, lVal = 0;
  for (UInt_t i=0; i<n; i++) {
    lVal = aVect[i];
    sum += lVal; sum2 += (lVal*lVal);
  }
  Float_t n1 = 1/(Float_t)n;
  sum *= n1;
  return TMath::Sqrt(TMath::Abs(sum2*n1 - sum*sum));
}

Float_t VectIntegral(const vector<Float_t> & aVect, Int_t xmin, Int_t xmax) {
  if (xmin < 0 || xmax < 0 || xmin >= (Int_t)aVect.size() ||
      xmax >= (Int_t)aVect.size() || xmin > xmax ) {
    cerr << "integral: out of bounds! - exiting" << endl; 
    cerr << xmin<<" "<<xmax<<" size: "<<aVect.size()<<endl; exit(-1);
  }
  Float_t lInt = 0;
  for (Int_t i=xmin; i<=xmax; i++) lInt+=aVect[i];
  return lInt;
}

Float_t VectIntegral(const vector<Float_t> & aVect) {
  Float_t lInt = 0;
  for (UInt_t i=0; i<aVect.size(); i++) lInt+=aVect[i];
  return lInt;
}

UInt_t VectMinBin(const vector<Float_t> & aVect) {
  UInt_t lMinBin = 0;
  Float_t lMin = aVect[0];
  for (UInt_t i=0; i<aVect.size(); i++) {
    if (lMin > aVect[i]) {
      lMinBin = i;
      lMin = aVect[i];
    }
  }
  return lMinBin;
}

UInt_t VectMaxBin(const vector<Float_t> & aVect) {
  UInt_t lMaxBin = 0 ;
  Float_t lMax = aVect[0];
  for (UInt_t i=0; i<aVect.size(); i++) {
    if (lMax < aVect[i]) {
      lMaxBin = i;
      lMax = aVect[i];
    }
  }
  return lMaxBin;
}

void VectMultiply(vector<Float_t> & fVect, const vector<Float_t> & aVect) {
  if (fVect.size() != aVect.size()) {
    cerr << "Multiply with vectors of different sizes. "<< fVect.size()<<" "<<aVect.size()<<" Exit"<< endl; exit(-1);
  }
  for (UInt_t i=0; i<fVect.size(); i++) (fVect[i])*=(aVect[i]);
}

void VectSquare(vector<Float_t> & fVect) {
  for (UInt_t i=0; i<fVect.size(); i++) (fVect[i])*=(fVect[i]);
}

void VectDivide(vector<Float_t> & fVect, const vector<Float_t> & aVect) {
  if (fVect.size() != aVect.size()) {
    cerr << "Divide with vectors of different sizes. "<< fVect.size()<<" "<<aVect.size()<<" Exit"<< endl; exit(-1);
  }
  for (UInt_t i=0; i<fVect.size(); i++) (fVect[i])/=(aVect[i]); // no test on aVect != 0 pour l'instant... ==> si on veut faire il faut tester avec finite() apres..
}

void VectAdd(vector<Float_t> & fVect, const vector<Float_t>& aVect, const Float_t c) {
  if (fVect.size() != aVect.size()) {
    cerr << "Add with vectors of different sizes!!. "<< fVect.size()<<" "<<aVect.size()<<" Exit"<< endl; exit(-1);
  }
  for (UInt_t i=0; i<fVect.size(); i++) fVect[i]+=(c*aVect[i]);
}

void VectScale(vector<Float_t> & fVect, const Float_t& aFactor) {
  for (UInt_t i=0; i<fVect.size(); i++) fVect[i]*=aFactor;
}

vector<Float_t> VectWindowFunction(UInt_t aNbPts, UInt_t aWidth) {
  vector<Float_t> lWindow(aNbPts,1);
  Float_t lOmega = 2*TMath::Pi()/(Float_t)(aNbPts-aWidth);
  for (UInt_t k=0; k<(aNbPts-aWidth)/2; k++)
    lWindow[k] = 0.5-0.5*cos(lOmega*k);
  for (UInt_t k=(aNbPts+aWidth)/2; k<aNbPts; k++)
    lWindow[k] = 0.5-0.5*cos(lOmega*(k-aWidth));  
  return lWindow;
}

//****************************************************
// Fonctions "database"
//****************************************************

// Fonction obsolete... plus besoin de fichier "num voies"
/*
int GetNumVoie(const string aFile, const string voie) {
  int numero=-1;
  ifstream f(aFile.c_str(),ifstream::in);
  if (!f.is_open()) {
    cerr << "GetNumVoie: file not open"<<endl;
    return numero;
  }
  string line;
  while (!f.eof()) {
    getline(f,line);
    if (!line.size()) continue;
    stringstream ss(line);
    string thevoie,dum; int thenum;
    getline(ss,dum,'"'); getline(ss,thevoie,'"'); ss >> thenum;
    if (thevoie==voie) {
      if (numero!=-1) cerr << "GetNumVoie: plusieurs fois la meme voie dans le fichier.." << endl;
      numero=thenum;
    }
  }
  f.close();
  return numero;
}
*/

vector<BoloStr> Read_liste_bolos(const string aFile) {
  vector<BoloStr> list_bolos;
  ifstream f(aFile.c_str(),ifstream::in);
  if (!f.is_open()) cerr << "Read_list_bolos: file not open"<<endl;
  string line,dum;
  getline(f,line); // entete
  bool old_fileversion=0; // Cas d'un fichier avec juste une BB indiquee : par convention BB2="NONE"
  if ( line.find("BB2")==string::npos ) old_fileversion=1;
  while (!f.eof()) {
    getline(f,line);
    if (!line.size()) continue;
    stringstream ss(line);
    BoloStr lbol;
    ss >> lbol.Name;
    ss >> lbol.Mac;
//    ss >> lbol.BB;
    ss >> lbol.BB1; 
    if (old_fileversion==0) ss >> lbol.BB2;
    else lbol.BB2="NONE";
    getline(ss,dum,'"'); getline(ss,lbol.Chal1,'"');
    getline(ss,dum,'"'); getline(ss,lbol.Chal2,'"');
    getline(ss,dum,'"'); getline(ss,lbol.Col1,'"');
    getline(ss,dum,'"'); getline(ss,lbol.Col2,'"');
    getline(ss,dum,'"'); getline(ss,lbol.Vet1,'"');
    getline(ss,dum,'"'); getline(ss,lbol.Vet2,'"');
    getline(ss,dum,'"'); getline(ss,lbol.Gar1,'"');
    getline(ss,dum,'"'); getline(ss,lbol.Gar2,'"');
    if (lbol.Chal1!="NONE") {
      (lbol.ChTypes).push_back("Chal1");
      (lbol.ChNames).push_back(lbol.Chal1);
    }
    if (lbol.Chal2!="NONE") {
      (lbol.ChTypes).push_back("Chal2");
      (lbol.ChNames).push_back(lbol.Chal2);
    }
    if (lbol.Col1!="NONE") {
      (lbol.ChTypes).push_back("Col1");
      (lbol.ChNames).push_back(lbol.Col1);
    }
    if (lbol.Col2!="NONE") {
      (lbol.ChTypes).push_back("Col2");
      (lbol.ChNames).push_back(lbol.Col2);
    }
    if (lbol.Vet1!="NONE") {
      (lbol.ChTypes).push_back("Vet1");
      (lbol.ChNames).push_back(lbol.Vet1);
    }
    if (lbol.Vet2!="NONE") {
      (lbol.ChTypes).push_back("Vet2");
      (lbol.ChNames).push_back(lbol.Vet2);
    }
    if (lbol.Gar1!="NONE") {
      (lbol.ChTypes).push_back("Gar1");
      (lbol.ChNames).push_back(lbol.Gar1);
    }
    if (lbol.Gar2!="NONE") {
      (lbol.ChTypes).push_back("Gar2");
      (lbol.ChNames).push_back(lbol.Gar2);
    }
    list_bolos.push_back(lbol);
  }
  f.close();
  return list_bolos;
}

vector<RunStr> Read_liste_runs(const string aFile, const BoloStr aBolo) {
  vector<RunStr> list_runs;
  ifstream f(aFile.c_str(),ifstream::in);
  if (!f.is_open()) cerr << "Read_list_runs: file not open"<<endl;
  string line,dum;
  getline(f,line); // entete

  while (!f.eof()) {
    getline(f,line);
    if (!line.size()) continue;
    stringstream ss(line);
    RunStr lrun;
    ss >> lrun.Name; ss >> lrun.Type;
    ss >> lrun.ModChal1; ss >> lrun.ModChal2;
    ss >> lrun.Vcol1; ss >> lrun.Vcol2;
    ss >> lrun.Vvet1; ss >> lrun.Vvet2;
    ss >> lrun.Vgar1; ss >> lrun.Vgar2;
    // calcul des tensions relatives utiles pour un (F)ID
    float vmean=0, nbelectrodes=0;
    if (aBolo.Col1!="NONE") {vmean+=lrun.Vcol1;nbelectrodes++;}
    if (aBolo.Col2!="NONE") {vmean+=lrun.Vcol2;nbelectrodes++;}
    if (aBolo.Vet1!="NONE") {vmean+=lrun.Vvet1;nbelectrodes++;}
    if (aBolo.Vet2!="NONE") {vmean+=lrun.Vvet2;nbelectrodes++;}
    if (aBolo.Gar1!="NONE") {vmean+=lrun.Vgar1;nbelectrodes++;}
    if (aBolo.Gar2!="NONE") {vmean+=lrun.Vgar2;nbelectrodes++;}
    if (nbelectrodes) vmean/=nbelectrodes;
    lrun.VrelCol1=lrun.Vcol1-vmean;
    lrun.VrelCol2=lrun.Vcol2-vmean;
    lrun.VrelVet1=lrun.Vvet1-vmean;
    lrun.VrelVet2=lrun.Vvet2-vmean;
    lrun.VrelGar1=lrun.Vgar1-vmean;
    lrun.VrelGar2=lrun.Vgar2-vmean;
    list_runs.push_back(lrun);
  }
  f.close();
  return list_runs;
}

RunStr Read_liste_run(const string aFile, const BoloStr aBolo, const string aRunName) {
  vector<RunStr> list_runs;
  ifstream f(aFile.c_str(),ifstream::in);
  if (!f.is_open()) cerr << "Read_list_runs: file not open"<<endl;
  string line,dum;
  getline(f,line); // entete

  while (!f.eof()) {
    getline(f,line);
    if (!line.size()) continue;
    stringstream ss(line);
    RunStr lrun;
    ss >> lrun.Name; ss >> lrun.Type;
    ss >> lrun.ModChal1; ss >> lrun.ModChal2;
    ss >> lrun.Vcol1; ss >> lrun.Vcol2;
    ss >> lrun.Vvet1; ss >> lrun.Vvet2;
    ss >> lrun.Vgar1; ss >> lrun.Vgar2;
    // calcul des tensions relatives utiles pour un (F)ID
    float vmean=0, nbelectrodes=0;
    if (aBolo.Col1!="NONE") {vmean+=lrun.Vcol1;nbelectrodes++;}
    if (aBolo.Col2!="NONE") {vmean+=lrun.Vcol2;nbelectrodes++;}
    if (aBolo.Vet1!="NONE") {vmean+=lrun.Vvet1;nbelectrodes++;}
    if (aBolo.Vet2!="NONE") {vmean+=lrun.Vvet2;nbelectrodes++;}
    if (aBolo.Gar1!="NONE") {vmean+=lrun.Vgar1;nbelectrodes++;}
    if (aBolo.Gar2!="NONE") {vmean+=lrun.Vgar2;nbelectrodes++;}
    if (nbelectrodes) vmean/=nbelectrodes;
    lrun.VrelCol1=lrun.Vcol1-vmean;
    lrun.VrelCol2=lrun.Vcol2-vmean;
    lrun.VrelVet1=lrun.Vvet1-vmean;
    lrun.VrelVet2=lrun.Vvet2-vmean;
    lrun.VrelGar1=lrun.Vgar1-vmean;
    lrun.VrelGar2=lrun.Vgar2-vmean;
    if (aRunName==lrun.Name) list_runs.push_back(lrun);
  }
  f.close();
  if (list_runs.size()!=1) {cerr << "Not found correct run "<<aRunName<<endl; exit(-1);}
  return list_runs[0];
}

string GetParam(string aFile, string aParam) {
  ifstream lFile(aFile.c_str(),ios::in);
  string s, theparam;
  theparam="NONE";
  while (lFile >> s) {
    if (s == aParam) {
      lFile >> s; 
      if (s != "=") cout << "Syntax error in parameter file."<<endl;
      lFile >> theparam;
    }
  }
  lFile.close();
  return theparam;
}

bool file_exists(string file) {
  ifstream f(file.c_str());
  return (bool)f;
}

BoloStr Read_liste_bolo(const string aFile, const string aBoloName) {
  vector<BoloStr> listebolos = Read_liste_bolos(aFile);
  Int_t i_bolo=-1;
  for (UInt_t p=0; p<listebolos.size(); p++) {
    if (listebolos[p].Name == aBoloName) i_bolo=(Int_t)p;
  }
  if (i_bolo==-1) {
    cerr << "Pas trouve le bolo!" << endl;
    exit(-1);
  }
  BoloStr bolo=listebolos[i_bolo]; 
  return bolo;
}

vector<PeriodStr> Read_list_periods(const string aFile, const string aRunName) {
  vector<PeriodStr> list_periods;
  ifstream f(aFile.c_str(),ifstream::in);
  if (!f.is_open()) cerr << "Read_list_periods: file not open"<<endl;
  string line;
  getline(f,line); // entete
  while (!f.eof()) {
    getline(f,line);
    if (!line.size()) continue;
    stringstream ss(line);
    PeriodStr lper;
    ss >> lper.Tinf;
    ss >> lper.Tsup;
    ss >> lper.Run;
    if (aRunName=="" || lper.Run==aRunName) list_periods.push_back(lper);
  }
  f.close();
  return list_periods;
}

vector<TemplateStr> Read_list_templates(const string aFile){
  vector<TemplateStr> list_templates;
  ifstream f(aFile.c_str(),ifstream::in);
  if (!f.is_open()) cerr << "Read_list_templates: file not open"<<endl;
  string line, dum; float toto;
  getline(f,line); // entete
  while (!f.eof()) {
    getline(f,line);
    if (!line.size()) continue;
    stringstream ss(line);
    TemplateStr tmplt;
    getline(ss,dum,'"'); getline(ss,tmplt.Channel,'"');
    ss >> tmplt.Type;
    ss >> tmplt.Tinf;
    ss >> tmplt.Tsup;
    while (!ss.eof()) {ss >> toto; (tmplt.Coefficients).push_back(toto);}
    list_templates.push_back(tmplt);
  }
  f.close();

  return list_templates;
}

TemplateStr Read_list_template(string aFile, string aChannel, string aType, ULong64_t aTime) {
  ifstream f(aFile.c_str(),ifstream::in);
  if (!f.is_open()) cerr << "Read_list_templates: file not open"<<endl;
  string line, dum; float toto;
  getline(f,line); // entete
  vector<TemplateStr> list_templates;
  while (!f.eof()) {
    getline(f,line);
    if (!line.size()) continue;
    TemplateStr tmplt;
    stringstream ss(line);
    getline(ss,dum,'"'); getline(ss,tmplt.Channel,'"');
    ss >> tmplt.Type;
    ss >> tmplt.Tinf;
    ss >> tmplt.Tsup;
    while (!ss.eof()) {ss >> toto; (tmplt.Coefficients).push_back(toto);}
    // !! Petit bricolage ici.. : si aucun type specifie, on enleve le type "ChalNTD"...
    if (aChannel==tmplt.Channel && (aType==tmplt.Type || (aType=="NONE" && tmplt.Type!="ChalNTD")) && (aTime==0 || (aTime>tmplt.Tinf && aTime<tmplt.Tsup)))
      list_templates.push_back(tmplt);
  }
  f.close();
  if (list_templates.size()!=1) {
    cerr << "Read_list_template: found "<<list_templates.size()<<" templates for "<<aChannel<<endl;
    exit(-1);
  }
  return list_templates.at(0);
}

TChain* ChainFromPartitions(string TraceDir, string Run, string Bolo) {
  DIR* dp = opendir(TraceDir.c_str());
  vector<string> partitions;
  struct dirent *dirp;
  while ((dirp=readdir(dp))!=NULL) {
    string toto=string(dirp->d_name);
    if (toto.find(Run+"_0")!=string::npos && toto.find(Bolo)!=string::npos)
      partitions.push_back(toto);
  }
  if (!partitions.size()) cerr << "No trace file for "<<Run<<endl;
  TChain* ch = new TChain("EdwTree","Edelweiss Raw Data Tree");
  for (UInt_t p=0; p<partitions.size(); p++) 
    ch->AddFile((TraceDir+"/"+partitions[p]).c_str());
  return ch;
}

bool DetectTrig(const vector<UInt_t> bittrigger, const int voienum) {
    // les bit triggers sont des unsigned int sur 32 bits.
  if (voienum<0) return 0 ; // pas de numero...

  //check size vecteur vs num bolo
  if (voienum > 32*bittrigger.size()) {
    cerr << "Detecttrig: pbl taille bittrig "<<voienum<<" "<<bittrigger.size()<<endl;
    return 0;
  }

  Bool_t b=0;
  UInt_t bit=0;
  UInt_t num_detecteur=(UInt_t)voienum; // cast...
  if (num_detecteur<32) {
    bit = ((bittrigger[0] >> num_detecteur) & 1);
  } else if (num_detecteur < 64) {
    UInt_t i=num_detecteur-32;
    bit = ((bittrigger[1] >> i) & 1);
  } else if (num_detecteur < 96) {
    UInt_t i=num_detecteur-64;
    bit = ((bittrigger[2] >> i) & 1);
  } else if (num_detecteur < 128) {
    UInt_t i=num_detecteur-96;
    bit = ((bittrigger[3] >> i) & 1);
  } else if (num_detecteur < 160) {
    UInt_t i=num_detecteur-128;
    bit = ((bittrigger[4] >> i) & 1);
  } else if (num_detecteur < 192) {
    UInt_t i=num_detecteur-160;
    bit = ((bittrigger[5] >> i) & 1);
  } else if (num_detecteur < 224) {
    UInt_t i=num_detecteur-192;
    bit = ((bittrigger[6] >> i) & 1);
  } else if (num_detecteur < 256) {
    UInt_t i=num_detecteur-224;
    bit = ((bittrigger[7] >> i) & 1);
  } else if (num_detecteur < 288) {
    UInt_t i=num_detecteur-256;
    bit = ((bittrigger[8] >> i) & 1);
  } else if (num_detecteur < 320) {
    UInt_t i=num_detecteur-288;
    bit = ((bittrigger[9] >> i) & 1);
  } else if (num_detecteur < 352) {
    UInt_t i=num_detecteur-320;
    bit = ((bittrigger[10] >> i) & 1);
  } else cerr << "Detect trig: Pbl num detecteur? "<<num_detecteur << endl;
  if (bit == 1) b=1;

  return b;
}
