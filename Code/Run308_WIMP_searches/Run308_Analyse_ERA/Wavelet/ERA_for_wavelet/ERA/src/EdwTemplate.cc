
#include "EdwTemplate.h"

ClassImp(EdwTemplate);

EdwTemplate::EdwTemplate() {
  fChannel="NONE";
  fType="NONE";
  fStartValidity=0; fEndValidity=0; fSign=0;
  fPretrigger=0; fTraceLength=0; fNormalisation=0;
}

EdwTemplate::~EdwTemplate() {
  this->Clear();
}

void EdwTemplate::Clear() {
  fChannel="NONE";
  fType="NONE";
  fStartValidity=0; fEndValidity=0; fSign=0;
  fPretrigger=0; fTraceLength=0; fNormalisation=0;
  fCoefficients.clear();
  fTrace.clear();
  fBinsOffsetFFT.clear();
  fDenoms.clear();
  for (UInt_t p=0; p<fTracesOffsetFFT.size(); p++) {
    (fTracesOffsetFFT[p]).clear();
    (fKernels[p]).clear();
  }
  fTracesOffsetFFT.clear();
  fKernels.clear();
}

void EdwTemplate::Initialize(string aChannel, string aType, ULong64_t aTinf, ULong64_t aTsup, vector<Float_t> aCoefs, UInt_t aLength, UInt_t aPretrig, Int_t aOffsetMin, Int_t aOffsetMax, vector<Float_t> aNoise) {
  // remplit les offsetfft a partir des donnees de liste_templates
  fChannel=aChannel;
  fType=aType;
  fStartValidity=aTinf; fEndValidity=aTsup;
  fPretrigger=aPretrig; fTraceLength=aLength;
  fCoefficients=aCoefs; 
  fNormalisation=1;
  vector<Float_t> trace_non_normalisee=this->ComputeTrace(0,0);
  Float_t a=VectMax(trace_non_normalisee), b=VectMin(trace_non_normalisee);
  fNormalisation = (a>-b) ? 1./a : -1./b;
  fSign = (a>-b) ? 1 : -1;
  fTrace = this->ComputeTrace(0,0); // avec la bonne normalisation
  
  for (Int_t k=aOffsetMin; k<=aOffsetMax; k++) {
    fBinsOffsetFFT.push_back(k);
    vector<Float_t> offset_trace=this->ComputeTrace((Float_t)k,1);
    fTracesOffsetFFT.push_back(EdwRealFFT(offset_trace));
  }
  this->ComputeKernels(aNoise); // remplit fDenoms et fKernels
}

Bool_t EdwTemplate::IsTimeValid(ULong64_t aTime) const {
  Bool_t ok=0;
  if (aTime >= fStartValidity && aTime <= fEndValidity) ok=1;
  return ok;
}

vector<Float_t> EdwTemplate::GetOffsetFFT(Int_t aOffset) const {
  vector<Float_t> lOffsetFFT;
  for (UInt_t i=0; i<fBinsOffsetFFT.size(); i++) {
    if (fBinsOffsetFFT[i] == aOffset) lOffsetFFT = fTracesOffsetFFT[i];
  }
  if (lOffsetFFT.size() == 0) {
    cerr<<"No good offset fft found for offset="<<aOffset<< endl;
    exit(-1); // sans pitie!!!
  }

  return lOffsetFFT;
}

vector<Float_t> EdwTemplate::GetKernel(Int_t aOffset) const {
  vector<Float_t> lKernel;
  for (UInt_t i=0; i<fBinsOffsetFFT.size(); i++) {
    if (fBinsOffsetFFT[i] == aOffset) lKernel=fKernels[i]; // ose.. a voir! (la fct est const)
  }
  if (lKernel.size()==0) {
    cerr<<"No good kernel found for offset="<<aOffset<< endl;
    exit(-1); // sans pitie!!!
  }
  return lKernel;
}

Float_t EdwTemplate::GetDenom(Int_t aOffset) const {
  Float_t lDenoms = 0;
  for (UInt_t i=0; i<fBinsOffsetFFT.size(); i++) {
    if (fBinsOffsetFFT[i] == aOffset) lDenoms = fDenoms[i];
  }
  if (lDenoms == 0) {
    cerr << "No good denom found for offset "<<aOffset<<endl;
    exit(-1);
  }
  return lDenoms;
}

void EdwTemplate::ComputeKernels(vector<Float_t> aNoise) {
  if (aNoise.size() == 0) {
    cerr << "ComputeKernels: empty noise. Using a flat noise spectrum as a rustine..."<<endl;
    aNoise.resize(fTraceLength,1);
  }
  if (aNoise.size() != fTraceLength) { cerr << "ComputeKernels: wrong noise size.." << endl; exit(-1);}

  if (fKernels.size() != 0 && fKernels.size() != fBinsOffsetFFT.size()) {
    cerr << "ComputeKernel: pbl with kernel size: "<<fKernels.size()<<" exiting"<<endl; exit(-1);
  }
  if (fDenoms.size() != 0 && fDenoms.size() != fBinsOffsetFFT.size()) {
    cerr << "ComputeKernel: pbl with denom size: "<<fDenoms.size()<<" exiting"<<endl; exit(-1);
  }
  if (fKernels.size() == 0) {
    fKernels.resize(fBinsOffsetFFT.size());
  }
  fDenoms.clear();
  for (UInt_t i=0; i<fKernels.size(); i++) {
    fKernels[i] = fTracesOffsetFFT[i] ;
    vector<Float_t> lK(fKernels[i]);
    vector<Float_t> lDenom(fKernels[i]);
    VectSquare(lDenom);
    VectDivide(lK,aNoise);
    VectDivide(lDenom,aNoise);
    Float_t dede = VectIntegral(lDenom);
    if (dede == 0) dede=1; // cas particulier trace tmplt offset = 00 ==> mettre 1 et ca doit passer
    fDenoms.push_back(1.0/dede); // !!! maj 01/01/2013... pour avoir le bon sigma(A)
    VectScale(lK,fDenoms[i]);
    fKernels[i] = lK;
  }
}

vector<Float_t> EdwTemplate::GetNonIntegerOffsetFFT(Float_t aOffset) const {
  // Methode retenue : interpolation lineaire des offsetfft entiers

  Int_t IntOffset = (Int_t)TMath::Floor(aOffset);
  Float_t alpha = aOffset - (Float_t)IntOffset;
  vector<Float_t> NonIntOffsetFFT(fTraceLength,0);
  // Securite:
  if (IntOffset < fBinsOffsetFFT[0] || (IntOffset+1)> fBinsOffsetFFT[fBinsOffsetFFT.size()-1]) {
    cout << "GetNonIntegerOffsetFFT: out of range!!" << endl;
    exit(-1); // sans pitie
  }
  vector<Float_t> IntOffsetFFT_inf = this->GetOffsetFFT(IntOffset); 
  vector<Float_t> IntOffsetFFT_sup = this->GetOffsetFFT(IntOffset+1); 
  for (UInt_t i=0; i<fTraceLength; i++)
    NonIntOffsetFFT[i]=(1-alpha)*IntOffsetFFT_inf[i]+alpha*IntOffsetFFT_sup[i];
  
  return NonIntOffsetFFT;
}

//****************************************
// Formules analytique des pulses modeles
//****************************************

vector<Float_t> EdwTemplate::ComputeTrace(Float_t offset, Bool_t windowflag) const {
  // "Driver routine"
  // Offset = offset par rapport au pretrigger
  vector<Float_t> trace(fTraceLength,0);
  if (fType=="StandardIon") {
    for (UInt_t p=0; p<fTraceLength; p++)
      trace[p]=this->PulseStandardIon((Float_t)p,offset+(Float_t)fPretrigger);
  } else if (fType=="StandardChal") {
    for (UInt_t p=0; p<fTraceLength; p++)
      trace[p]=this->PulseStandardChal((Float_t)p,offset+(Float_t)fPretrigger);
  } else if (fType=="IonBB2") {
    for (UInt_t p=0; p<fTraceLength; p++)
      trace[p]=this->PulseIonBB2((Float_t)p,offset+(Float_t)fPretrigger);
  } else if (fType=="ChalNTD") {
    for (UInt_t p=0; p<fTraceLength; p++)
      trace[p]=this->PulseChalNTD((Float_t)p,offset+(Float_t)fPretrigger);
  } else {cerr << "ComputeTrace for template:  wrong type"<<endl; exit(-1);}
  VectScale(trace,fNormalisation);

  if (windowflag) {
    Int_t width = 0;
    if ((fType).find("Chal")!=string::npos) width=PLATEAU_WINDOW_HEAT;
    else if ((fType).find("Ion")!=string::npos) width = (Int_t)fTraceLength - 2*RISE_WINDOW_ION;
    else {cerr << "EdwTemplate: tmplt type does not tell if heat or ion..." << endl; exit(-1);}
    vector<Float_t> lWindow=VectWindowFunction(fTraceLength,width);
    VectMultiply(trace,lWindow);
  }

  return trace;
}

Float_t EdwTemplate::PulseStandardIon(Float_t x, Float_t offset) const {
  if (fCoefficients.size()!=3) {
    cerr <<"PulseIonStandard : wrong nb of coefs" << endl; exit(-1);
  }
  Float_t f = 0;
  x-=offset;
  if (x>=0) 
    f = ((1+TMath::Erf(x/fCoefficients[2]))/2.0)*(exp(-x/fCoefficients[0])-fCoefficients[1]);
  return f;
}

Float_t EdwTemplate::PulseIonBB2(Float_t x, Float_t offset) const {
  if (fCoefficients.size()!=1) {
    cerr <<"PulseIonBB2 : wrong nb of coefs" << endl; exit(-1);
  }
  Float_t f = 0;
  x-=offset;
  if (x>=0) 
    f = (1+TMath::Erf(x/fCoefficients[0]))/2.0;
  return f;
}

Float_t EdwTemplate::PulseStandardChal(Float_t x, Float_t offset) const {
  if (fCoefficients.size()!=6) {
    cerr <<"PulseChalStandard : wrong nb of coefs" << endl; exit(-1);
  }
  Double_t f = 0;
  x-=offset;
  if (x >=0)
    f = (1-exp(-x/fCoefficients[0]))
      * (exp(-x/fCoefficients[1]) + fCoefficients[2]*exp(-x/fCoefficients[3])
         + fCoefficients[4]*exp(-x/fCoefficients[5])) ;
  return f;
}

Float_t EdwTemplate::PulseChalNTD(Float_t x, Float_t offset) const {
  if (fCoefficients.size()!=2) {
    cerr <<"PulseChalNTD : wrong nb of coefs" << endl; exit(-1);
  }
  Double_t f = 0;
  x-=offset;
  if (x >=0)
    f = (1-exp(-x/fCoefficients[0]))*exp(-x/fCoefficients[1]) ;
  return f;
}
