
#include "FitPulse.h"

ClassImp(FitPulse);

FitPulse::FitPulse() : EdwPulse() {
  fBasicPreprocessed = 0;
  fWindowed = 0;
}

FitPulse::~FitPulse() {
  fProcessedTrace.clear();
  fProcessedTraceFFT.clear();
  fPattern.clear();
  fPeakBins.clear();
  fTrace.clear();
}

FitPulse::FitPulse(const EdwPulse* aPulse) : EdwPulse() {
  // Copie directe de tout ce que contient aPulse
  fChannel=aPulse->Channel() ;
  fIsHeat = aPulse->IsHeat();
  fIsBBv2 = aPulse->IsBBv2();
  fSampling_ms=aPulse->Sampling_ms() ;
  fPretrigger=aPulse->Pretrigger() ;
  fTrace = aPulse->Trace() ;
  fTraceLength = aPulse->TraceLength();

  fProcessedTrace.resize(fTrace.size(),0);
  for (UInt_t i=0; i<fTrace.size(); i++)
    fProcessedTrace[i]=(Float_t)(fTrace[i]);
  
  fSign=aPulse->Sign();
  fPatternLength=aPulse->PatternLength();
  fModulationLength=aPulse->ModulationLength();

  fBaseStart = 0;
  //  fBaseStop = (Int_t)(PRETRIGGER_SECURITY*(Float_t)fPretrigger);
  // Changement dans cette definition...
  fBaseStop = (Int_t)(PRETRIGGER_SECURITY*(Float_t)(fPretrigger+this->FitTimeMin()));
  fBasicPreprocessed = 0;
  fWindowed = 0;
}

FitPulse::FitPulse(const EdwPulse* aPulse1, const EdwPulse* aPulse2, Float_t lambda1, Float_t lambda2, string comb_name) : EdwPulse() {

  fChannel=comb_name;
  fIsHeat=aPulse1->IsHeat();
  fIsBBv2=aPulse1->IsBBv2();// bon, on gere pas col1 en bbv1 et col2 en bbv2 :)
  fSampling_ms=aPulse1->Sampling_ms();
  fPretrigger=aPulse1->Pretrigger();
  fTraceLength=aPulse1->TraceLength();
  Short_t s1=(lambda1>0) ? aPulse1->Sign() : -aPulse1->Sign() ;
  Short_t s2=(lambda2>0) ? aPulse2->Sign() : -aPulse2->Sign() ;
  fSign=s1;
  if (aPulse2->TraceLength()!=fTraceLength || aPulse2->Pretrigger()!=fPretrigger ||
      aPulse2->Sampling_ms()!=fSampling_ms || aPulse2->IsHeat()!=fIsHeat || s2!=s1)
    cerr << "Fitpulse from 2 edwpulses: coherent?" << endl;
  fPatternLength=aPulse1->PatternLength();
  if (aPulse2->PatternLength() > fPatternLength) fPatternLength=aPulse2->PatternLength();

  // trace = lambda1*trace1 + lambda2*trace2  
  fProcessedTrace.resize(fTraceLength,0);
  fTrace.resize(fTraceLength,0);
  for (UInt_t i=0; i<fTraceLength; i++) {
    fProcessedTrace[i]=(Float_t)(aPulse1->Trace(i))*lambda1+(Float_t)(aPulse2->Trace(i))*lambda2;
    fTrace[i]=(Short_t)fProcessedTrace[i];
  }
  fBaseStart = 0;
  //  fBaseStop = (Int_t)(PRETRIGGER_SECURITY*fPretrigger);
  // Changement dans cette definition...
  fBaseStop = (Int_t)(PRETRIGGER_SECURITY*(Float_t)(fPretrigger+this->FitTimeMin()));
  fBasicPreprocessed = 0;
  fWindowed = 0;
}

//*******************************************************
// A) Processing "de base"
//*******************************************************

// A VOIR: IL Y A UN PBL A LA FIN DE LA TRACE PREPROCESSEE
// DANS LE CAS ON OU FAIT aCheckPeak=1..!!
// hyp: removepattern pas ok si ..??

Int_t FitPulse::FitTimeMin() const {
  if (fIsHeat==1) return FIT_TIMEMIN_HEAT;
  else if (fIsBBv2) return FIT_TIMEMIN_IONBB2;
  else return FIT_TIMEMIN_ION;
}

Int_t FitPulse::FitTimeMax() const {
  if (fIsHeat==1) return FIT_TIMEMAX_HEAT;
  else if (fIsBBv2) return FIT_TIMEMAX_IONBB2;
  else return FIT_TIMEMAX_ION;
}

Bool_t FitPulse::PeakOutsideWindow() {
  // Verifie si parmi les pics detectes, il y en a un hors de la fenetre de recherche naturelle..
  Bool_t peakoutside=0;
  Int_t marge=5; // TODO = A TUNER... (peut-etre en fct de la taille trace etc)
  this->FindPeaks(FP_WEAK,1); // "anysign=1"
  for (UInt_t i=0; i<fPeakBins.size(); i++) {
    if (fPeakBins[i] < fPretrigger+this->FitTimeMin()-marge || 
	fPeakBins[i] > fPretrigger+this->FitTimeMax()+marge)
      peakoutside=1;
  }

  return peakoutside;
}

void FitPulse::BasicPreprocess(Bool_t aCheckPeak) {
  
  if (!this->TraceFullySaturated())
    this->RemoveSaturatedSpikes(); // vraiment a faire? (cas de la somme de trace..)
  this->RemoveBaseline();
  this->RemovePattern(fPatternLength);
  
  if (fIsBBv2 && !fIsHeat) {
    vector<Float_t> lDerivTrace(fProcessedTrace.size()-1,0);
    for (UInt_t i=0; i<fProcessedTrace.size()-1; i++)
      lDerivTrace[i]=fProcessedTrace[i+1]-fProcessedTrace[i];
    if (fSign == -1) {
      fSimpleAmpl = VectMin(lDerivTrace);
      fSimpleAmplBin = VectMinBin(lDerivTrace);
    } else if (fSign == 1) {
      fSimpleAmpl = VectMax(lDerivTrace);
      fSimpleAmplBin = VectMaxBin(lDerivTrace);
    } else cerr <<"BasicPreprocess: Unknown sign.."<<endl;
  } else {
    if (fSign == -1) {
      fSimpleAmpl = VectMin(fProcessedTrace);
      fSimpleAmplBin = VectMinBin(fProcessedTrace);
    } else if (fSign == 1) {
      fSimpleAmpl = VectMax(fProcessedTrace);
      fSimpleAmplBin = VectMaxBin(fProcessedTrace);
    } else {
      Float_t lMin = VectMin(fProcessedTrace);
      fSimpleAmpl = VectMax(fProcessedTrace);
      if (-lMin >= fSimpleAmpl) {
	fSimpleAmpl = lMin;
	fSimpleAmplBin = VectMinBin(fProcessedTrace);
      } else fSimpleAmplBin = VectMaxBin(fProcessedTrace);
    }
  }
  fBasicPreprocessed = 1;

  // if a pulse takes place in the beginning of the trace it might bias the pattern
  if (aCheckPeak) {
    this->FindPeaks(FP_WEAK); // discutable encore..
    if (fPeakBins.size()) {
      Int_t lMinPos=fPeakBins[0];
      for (UInt_t i=1; i<fPeakBins.size(); i++) {
	if (fPeakBins[i]<lMinPos) lMinPos = fPeakBins[i];
      }
      //      Int_t lMinPos = VectMin(fPeakBins);
      if (lMinPos <= fBaseStop && lMinPos > fPatternLength) {
	// An early pulse has biased the pattern!!
	// (note that if the pulse is really really very early I can't do anything...)
        fBaseStop = lMinPos - 1 ;
        if (fBaseStop <= fBaseStart+9) fBaseStop = fBaseStart+10;
        this->BasicPreprocess(0); // no peak check this time!!! (no infinite loop)
      }
    }
  }
}

bool FitPulse::TraceFullySaturated() {
  bool IsFullySaturated=0;
  Float_t nbsatpoints=0;
  for (UInt_t i=0; i<fTrace.size(); i++) {
    if (fTrace[i] >= SHRT_MAX || fTrace[i] <= SHRT_MIN) nbsatpoints++;
  }
  Float_t frac_full_sat=0.8; // param en dur... a tuner... 
  if (nbsatpoints > frac_full_sat*(Float_t)fTraceLength) IsFullySaturated=1;
  return IsFullySaturated;
}

void FitPulse::RemoveSaturatedSpikes(Int_t maxlength) {
  // This step should take place before any other preprocessing (eg. pattern, baseline etc..)
  // Replace saturated points by a linear interpolation of their nearest neighboors.
  // Do this only if the saturation takes place on less than maxlength points.
  // This is a "spike removal" procedure, NOT a processing of saturation due to a real signal.

  // First get the list of points where saturation takes place
  vector<UInt_t> lSatPoints;
  for (UInt_t i=0; i<fTrace.size(); i++) {
    if (fTrace[i] >= SHRT_MAX || fTrace[i] <= SHRT_MIN) lSatPoints.push_back(i);
  }

  // Correction if needed:
  if (lSatPoints.size() != 0) {
    for (UInt_t i=0; i<lSatPoints.size(); i++) {
      Int_t i_pt=(Int_t)(lSatPoints[i]);
      // Find the nearest unsaturated neighboors
      Int_t i_inf=i_pt;
      while (fTrace[i_inf] >= SHRT_MAX || fTrace[i_inf] <= SHRT_MIN) {
        i_inf--;
        if (i_inf<0) break;
      }
      Int_t i_sup=i_pt;
      while (fTrace[i_sup] >= SHRT_MAX || fTrace[i_sup] <= SHRT_MIN) {
        i_sup++;
        if (i_sup >= (Int_t)fTrace.size()) break;
      }
      // Apply the linear interpolation if the "spike length" is not too long
      if ( (i_sup-i_inf-1) <= maxlength ) {
        if (i_inf < 0) {
          fProcessedTrace[i]=fProcessedTrace[i_sup];
        } else if (i_sup >= (Int_t)fTrace.size()) { 
          fProcessedTrace[i]=fProcessedTrace[i_inf];
        } else {
          fProcessedTrace[i]=fProcessedTrace[i_inf]+
            (i-i_inf)*(fProcessedTrace[i_sup]-fProcessedTrace[i_inf])/(Float_t)(i_sup-i_inf);
        }
      }
    }
  }
}

UInt_t FitPulse::FindGlitches(Int_t glitch_ampl) {
  UInt_t nbglitches=0;
  for (UInt_t p=1; p<=fTraceLength; p++) {
    if (abs(fTrace[p]-fTrace[p-1]) > glitch_ampl) nbglitches++;
  }
  return nbglitches;
}

void FitPulse::RemoveBaseline() {

  if (fIsHeat || !fIsBBv2) {
    fMeanBase = VectIntegral(fProcessedTrace,fBaseStart,fBaseStop)/(Float_t)(fBaseStop-fBaseStart+1) ;
    for (UInt_t i=0; i<fProcessedTrace.size(); i++) fProcessedTrace[i]-=fMeanBase;
  } else {
    // Cas d'une ionisation bbv2 : pente constante
    // donc on soustrait une regression lineaire
    // formule y=ax+b avec a=Nsum(xy)-sum(x)sum(y) / Nsum(xx)-(sumx)^2 et b=(sum(y)-asum(x))/N
    // et on applique a x_i=i [??]
    // FORMULE OK SI BASESTART=0
    Float_t sumy=VectIntegral(fProcessedTrace,0,fBaseStop);
    Float_t sumiy=0;
    for (Int_t i=0; i<=fBaseStop; i++) sumiy+=(i*fProcessedTrace[i]);
    Float_t nbpts=(Float_t)fBaseStop+1;
    Float_t a=(12*sumiy-6*(nbpts-1)*sumy)/(nbpts*(nbpts-1)*(nbpts+1));
    Float_t b=sumy/nbpts-0.5*a*(nbpts-1);
    for (UInt_t i=0; i<fProcessedTrace.size(); i++) fProcessedTrace[i]-=(a*i+b);
  }

}

void FitPulse::RemoveLinearBase() {
  // subtract y=ax+b
  // Routine necessaire slt pour le calcul spectres de bruit
  UInt_t n = fProcessedTrace.size();
  Float_t b = fProcessedTrace[0];
  Float_t a = (fProcessedTrace[n-1]-fProcessedTrace[0])/(Float_t)(n-1);
  for (UInt_t k=0; k<n; k++) fProcessedTrace[k] -= (a*k+b);
  this->RemoveBaseline(); // necessary to remove the baseline once again after
}

void FitPulse::RemovePattern(Int_t aNbPatternPts) {
  // Allows to remove a noise with perfectly known period and which remains with a constant
  // phase all over the event (better than single filtering)
  // Especially useful to correct for the heat modulation noise in the ionisation channels

  if (aNbPatternPts != 0) {
    // Construct the pattern
    fPattern.clear();
    fPattern.resize(aNbPatternPts,0);
    Int_t lNbLoops = (Int_t)TMath::Floor((fBaseStop-fBaseStart+1)/(Float_t)aNbPatternPts) ;
    for (Int_t i=0;i<lNbLoops;i++){
      Int_t lIter = aNbPatternPts*i;
      for (Int_t j=0;j<aNbPatternPts;j++) fPattern[j] += (fProcessedTrace[j+lIter]) ;
      // NO MORE factor (+1) = ROOT binning convention
    }
    for (Int_t j=0;j<aNbPatternPts;j++) fPattern[j]/=lNbLoops;
    
    // Substract the pattern
    Int_t lShift, j;
    for (UInt_t i=0; i<fProcessedTrace.size(); i++) {
      lShift = (Int_t)TMath::Floor((i-fBaseStart)/(Float_t)aNbPatternPts); 
      // Need a good drawing!
      // lShift = index in the trace of the nearest beginning of the pattern from point i.
      j = i-(fBaseStart+lShift*aNbPatternPts);
      fProcessedTrace[i] -= fPattern[j]; // NO MORE Root binning
    }
  }
}

void FitPulse::SetBaseFromPeakBin(Int_t aBin) {
  if (aBin != 0) {
    fBaseStart = 0;
    // ERAnew : pour l'instant on a enleve "+FIT_TIMEMIN_HEAT/ION" ... a voir plus tard..
    //    if (this->IsHeat()) fBaseStop = (Int_t)(PRETRIGGER_SECURITY*(aBin)) ; // sign
    //    else 
    fBaseStop = (Int_t)(PRETRIGGER_SECURITY*(aBin)) ;
  }
}

void FitPulse::ApplyWindow(Int_t aWidth) {
  if (aWidth == 0) {
    if (this->IsHeat()) aWidth = PLATEAU_WINDOW_HEAT;
    else aWidth = (Int_t)fProcessedTrace.size() - 2*RISE_WINDOW_ION;
  }
  vector<Float_t> lWindowFunction =VectWindowFunction(fProcessedTrace.size(),aWidth);
  VectMultiply(fProcessedTrace,lWindowFunction);
  fWindowed = 1;
}

void FitPulse::FindPeaks(Int_t aCriterium, Bool_t AnySign, Int_t aLength) {
  if (aCriterium == FP_WEAK && this->IsHeat())
    this->FindPeaksWithParams(FP_SENS_HEAT_WEAK,aLength,FP_ORDER_HEAT_WEAK,AnySign);
  else if (aCriterium == FP_STRICT && this->IsHeat())
    this->FindPeaksWithParams(FP_SENS_HEAT_STRICT,aLength,FP_ORDER_HEAT_STRICT,AnySign);
  else if (aCriterium == FP_WEAK && (!this->IsHeat()))
    this->FindPeaksWithParams(FP_SENS_ION_WEAK,aLength,FP_ORDER_ION_WEAK,AnySign);
  else if (aCriterium == FP_STRICT && (!this->IsHeat()))
    this->FindPeaksWithParams(FP_SENS_ION_STRICT,aLength,FP_ORDER_ION_STRICT,AnySign);
  else if (aCriterium == FP_XWEAK && this->IsHeat())
    this->FindPeaksWithParams(FP_SENS_HEAT_XWEAK,aLength,FP_ORDER_HEAT_XWEAK,AnySign);
  else if (aCriterium == FP_XWEAK && (!this->IsHeat()))
    this->FindPeaksWithParams(FP_SENS_ION_XWEAK,aLength,FP_ORDER_ION_XWEAK,AnySign);
  else cerr << "FindPeaks : incorrect criterium." << endl; 
}

void FitPulse::FindPeaksWithParams(Float_t aSensitivity, Int_t aLength, Int_t aOrder, Bool_t AnySign) {
  fPeakBins.clear();
  if (aSensitivity <= 0 || aLength <= 0 || aOrder <= 0) 
    cerr << "FindPeaks: negative input psrameter!" << endl;
  if (!fBasicPreprocessed) this->BasicPreprocess();

  vector<Float_t> lDerivTrace(fProcessedTrace.size()-aOrder);
  for (UInt_t i=aOrder; i<fProcessedTrace.size(); i++) 
    lDerivTrace[i-aOrder]=fProcessedTrace[i]-fProcessedTrace[i-aOrder];
  Int_t lNbPts=lDerivTrace.size();
  Double_t lThreshold = aSensitivity*VectRMS(lDerivTrace);
  for (Int_t i=1; i<lNbPts; i++) {
    if (( (AnySign || fSign >= 0) && lDerivTrace[i] >= lThreshold) ||
        ( (AnySign || fSign <= 0) &&  lDerivTrace[i] <= -lThreshold)) {
      Int_t j_prev = 0 ;
      if ( fPeakBins.size()!= 0  ) j_prev=fPeakBins[fPeakBins.size()-1];
      if (j_prev == 0 || i > j_prev+aLength) fPeakBins.push_back(i);
    }
  }

}

void FitPulse::ComputeTraceFFT(Bool_t aWindowFlag) {
  if (!fBasicPreprocessed) this->BasicPreprocess();
  if (aWindowFlag && !fWindowed) this->ApplyWindow();
  fProcessedTraceFFT = EdwRealFFT(fProcessedTrace);
}

vector<Float_t> FitPulse::SmoothedTrace(UInt_t aSmoothing) const {
  // top-hat window
  vector<Float_t> lOutput = fProcessedTrace;
  UInt_t N = lOutput.size();
  for (UInt_t bin=0; bin<N; bin++) {
    if (bin >= aSmoothing && bin+aSmoothing < N) { // (car on va tomber hors trace)
      lOutput[bin]=0;
      for (Int_t k=(Int_t)bin-(Int_t)aSmoothing; k<=(Int_t)bin+(Int_t)aSmoothing; k++) lOutput[bin]+=fProcessedTrace[k];
      lOutput[bin] /= ((Float_t)(1+2*aSmoothing));
    }
  }
  return lOutput;
}

vector<Float_t> FitPulse::SmoothedAmpl(const vector<Float_t>& aSmoothTrace) const {
  // Cette routine designee pour le PSA pourrait etre re-designee... cosmetique
  Float_t ampl, amplbin;
  if (fSign == -1) {
    ampl = VectMin(aSmoothTrace);
    amplbin = (Float_t)VectMinBin(aSmoothTrace);
  } else if (fSign == 1) {
    ampl = VectMax(aSmoothTrace); // here maybe need to iterate if fPulseBin < fBaseStop...
    amplbin = (Float_t)VectMaxBin(aSmoothTrace);
  } else {
    Float_t lMin = VectMin(aSmoothTrace);
    ampl = VectMax(aSmoothTrace);
    if (-lMin >= ampl) {
      ampl = lMin;
      amplbin = (Float_t)VectMinBin(aSmoothTrace);
    } else amplbin = (Float_t)VectMaxBin(aSmoothTrace);
  }
  vector<Float_t> out(2,0); out[0]=ampl; out[1]=amplbin;
  return out;
}

vector<Float_t> FitPulse::GetRiseTimeParams(const Float_t aFracLow, const Float_t aFracHigh, UInt_t aSmooth) {
  // Amelioration? : chercher le rise time dans une fenetre reduite... (en fait idem pour fSimpleAmpl eventuellement..)
  if (!fBasicPreprocessed) this->BasicPreprocess();
  if (fWindowed) cerr << "RiseTime: windowed pulse" << endl;
  vector<Float_t> lSmoothTrace = this->SmoothedTrace(aSmooth);
  vector<Float_t> lSmoothAmpl = this->SmoothedAmpl(lSmoothTrace);

  Float_t lTStart = 1;
  Float_t lTEnd = lSmoothAmpl[1];
  if (lTEnd == 1) {
    vector<Float_t> params(4,0);
    return params;
  }
  Float_t lAmplStart = aFracLow*lSmoothAmpl[0];
  Float_t lAmplEnd = aFracHigh*lSmoothAmpl[0];
  // Iteration to find lTStart and lTEnd
  Int_t i = (Int_t)lTEnd;
  Float_t lPtInf = lSmoothTrace[i];
  Float_t lPtSup;
  do {
    i-=1;
    lPtSup = lPtInf;
    lPtInf = lSmoothTrace[i];
    if (fSign < 0) { // if fsign=0, by convention will look rising pulse...
      if (lPtSup <= lAmplEnd && lPtInf >= lAmplEnd)
        lTEnd = i+(lAmplEnd-lPtInf)/(lPtSup-lPtInf);
      if (lPtSup <= lAmplStart && lPtInf >= lAmplStart)
        lTStart = i+(lAmplStart-lPtInf)/(lPtSup-lPtInf);
      
    } else {
      if (lPtSup >= lAmplEnd && lPtInf <= lAmplEnd)
        lTEnd = i+(lAmplEnd-lPtInf)/(lPtSup-lPtInf);
      if (lPtSup >= lAmplStart && lPtInf <= lAmplStart)
        lTStart = i+(lAmplStart-lPtInf)/(lPtSup-lPtInf);
    }
  } while (lTStart == 1 && i!=0);

  vector<Float_t> params;
  params.push_back(lTEnd-lTStart);
  params.push_back(lTStart);
  params.push_back(lAmplStart);
  params.push_back(lAmplEnd);
  return params;
}

Float_t FitPulse::GetRiseTime(const Float_t aFracLow, const Float_t aFracHigh, UInt_t aSmooth) {
  vector<Float_t> params=this->GetRiseTimeParams(aFracLow,aFracHigh,aSmooth);
  return params[0];
}

vector<Float_t> FitPulse::GetFallTimeParams(const Float_t aFracLow, const Float_t aFracHigh, UInt_t aSmooth) {
  if (!fBasicPreprocessed) this->BasicPreprocess();
  if (fWindowed) cerr << "FallTime warning: windowed pulse" << endl;
  vector<Float_t> lSmoothTrace = this->SmoothedTrace(aSmooth);
  vector<Float_t> lSmoothAmpl = this->SmoothedAmpl(lSmoothTrace);
  Float_t lTStart = lSmoothAmpl[1];
  Float_t lLast = fProcessedTrace.size()-1;
  Float_t lTEnd = lLast;
  if (lTStart == lLast) {
    vector<Float_t> params(4,0);
    return params;
  }
  Float_t lAmplStart = aFracHigh*lSmoothAmpl[0];
  Float_t lAmplEnd = aFracLow*lSmoothAmpl[0];
  // Iteration to find lTStart and lTEnd
  Int_t i = (Int_t)lTStart;
  Float_t lPtInf = lSmoothTrace[i];
  Float_t lPtSup;

  do {
    lPtSup = lPtInf;
    lPtInf = lSmoothTrace[i+1];
    if (fSign < 0) {
      if (lPtSup <= lAmplEnd && lPtInf >= lAmplEnd)
        lTEnd = i+(lPtSup-lAmplEnd)/(lPtSup-lPtInf);
      if (lPtSup <= lAmplStart && lPtInf >= lAmplStart)
        lTStart = i+(lPtSup-lAmplStart)/(lPtSup-lPtInf);
    } else {
      if (lPtSup >= lAmplEnd && lPtInf <= lAmplEnd)
        lTEnd = i+(lPtSup-lAmplEnd)/(lPtSup-lPtInf);
      if (lPtSup >= lAmplStart && lPtInf <= lAmplStart)
        lTStart = i+(lPtSup-lAmplStart)/(lPtSup-lPtInf);
    }
    i+=1;
  } while (lTEnd == lLast && i!= lLast);

  vector<Float_t> params;
  params.push_back(lTEnd-lTStart);
  params.push_back(lTStart);
  params.push_back(lAmplStart);
  params.push_back(lAmplEnd);
  return params;
}

Float_t FitPulse::GetFallTime(const Float_t aFracLow, const Float_t aFracHigh, UInt_t aSmooth) {
  vector<Float_t> params=this->GetFallTimeParams(aFracLow,aFracHigh,aSmooth);
  return params[0];
}

vector<Float_t> FitPulse::GetFWHMParams(const Float_t aFracWidth, UInt_t aSmooth) {
  if (!fBasicPreprocessed) this->BasicPreprocess();
  if (fWindowed) cerr << "getfwhm warning: windowed pulse" << endl;
  vector<Float_t> lSmoothTrace = this->SmoothedTrace(aSmooth);
  vector<Float_t> lSmoothAmpl = this->SmoothedAmpl(lSmoothTrace);
  Int_t N = lSmoothTrace.size();
  Float_t t1 = 0;
  Float_t t2 = (Float_t)N-1;
  Float_t lAmplCut = aFracWidth*lSmoothAmpl[0];
  // Iteration to find t1
  Int_t i = (Int_t)lSmoothAmpl[1];
  Float_t lPtInf = lSmoothTrace[i];
  Float_t lPtSup;
  do { 
    i-=1;
    lPtSup = lPtInf;
    lPtInf = lSmoothTrace[i];
    if (fSign < 0) { // if fsign=0, by convention will look rising pulse...
      if (lPtSup <= lAmplCut && lPtInf >= lAmplCut)
        t1 = i+(lAmplCut-lPtInf)/(lPtSup-lPtInf);
    } else {
      if (lPtSup >= lAmplCut && lPtInf <= lAmplCut)
        t1 = i+(lAmplCut-lPtInf)/(lPtSup-lPtInf);
    }
  } while (t1 == 0 && i!=0);
  // Iteration to find t2
  i = (Int_t)lSmoothAmpl[1];
  lPtInf = lSmoothTrace[i];
  do {
    lPtSup = lPtInf;
    lPtInf = lSmoothTrace[i+1];
    if (fSign < 0) {
      if (lPtSup <= lAmplCut && lPtInf >= lAmplCut)
        t2 = i+(lPtSup-lAmplCut)/(lPtSup-lPtInf);
    } else {
      if (lPtSup >= lAmplCut && lPtInf <= lAmplCut)
        t2 = i+(lPtSup-lAmplCut)/(lPtSup-lPtInf);
    }
    i+=1;
  } while ((Int_t)t2 == N-1 && i!= N-1);

  vector<Float_t> params;
  params.push_back(t2-t1);
  params.push_back(t1);
  params.push_back(lAmplCut);
  return params;
}

Float_t FitPulse::GetFWHM(const Float_t aFracWidth, UInt_t aSmooth) {
  vector<Float_t> vect=this->GetFWHMParams(aFracWidth,aSmooth);
  return vect[0];
}



//****************************************************
// B) Filtrage type optimal
// Pour l'instant on reprend pas les fonctions de fit basic utilisant minuit
//****************************************************

vector<Float_t> FitPulse::ComputeWienerAmpl(const vector<Float_t>& aTemplateFFT, const vector<Float_t>& aNoise, Bool_t aChi2Switch) {
  // Optimal filtering in Fourier space
  // aTemplateFFT = fft of the template in "semi-complex" form
  // aNoise = power spectrum of the noise

  // 0) Preliminaires
  vector<Float_t> lOutput(4,0);
  if (aNoise.size() != fTraceLength || aTemplateFFT.size()!=fTraceLength) {
    cerr << "GetWienerAmpl: noise/template has wrong size" << endl;
    return lOutput;
  }
  if (fProcessedTraceFFT.size() == 0) this->ComputeTraceFFT();
  bool zeronoise=0;
  //  if (aNoise.Integral() == 0) zeronoise=1;
  for (UInt_t p=0; p<fTraceLength; p++) {
     if (aNoise[p] == 0) zeronoise=1;
  }
  
  // 1) Calcul amplitude
  vector<Float_t> lDenum = aTemplateFFT;
  VectMultiply(lDenum,aTemplateFFT);
  if (!zeronoise) VectDivide(lDenum,aNoise);
  vector<Float_t> lNumer = aTemplateFFT;
  VectMultiply(lNumer,fProcessedTraceFFT);
  if (!zeronoise) VectDivide(lNumer,aNoise);
  Float_t toto = VectIntegral(lDenum);
  if (toto == 0) toto=1;
  Float_t DenumInt = 1/toto;
  Float_t lAmpl = DenumInt*VectIntegral(lNumer);
  lOutput[0]=lAmpl;
  // 1bis) Sigma(amplitude) .. a etudier!!!
  lOutput[3]=sqrt(DenumInt);

  // 2) En option le Chi2
  if (aChi2Switch) {
    vector<Float_t> lChi2Histo=aTemplateFFT;
    VectScale(lChi2Histo,lAmpl);
    VectAdd(lChi2Histo,fProcessedTraceFFT,-1);
    VectSquare(lChi2Histo);
    if (!zeronoise) VectDivide(lChi2Histo,aNoise);
    lOutput[2]=VectIntegral(lChi2Histo)/fTraceLength;
  }

  return lOutput;
}

vector<Float_t> FitPulse::ComputeWienerFast(const vector<Float_t>& aKernel, const Float_t aDenom) {
  if (fProcessedTraceFFT.size() == 0) this->ComputeTraceFFT();
  vector<Float_t> lConvol = fProcessedTraceFFT;
  VectMultiply(lConvol,aKernel);

  vector<Float_t> lOutput(4,0);
  lOutput[0]=VectIntegral(lConvol);
  lOutput[3]=sqrt(aDenom);
  return lOutput;
}

vector<Float_t> FitPulse::WienerLoop(const EdwTemplate* aTemplate, const vector<Float_t>& aNoise, Int_t OffsetMin, Int_t OffsetMax, Float_t OffsetStep, Bool_t aFast) {

  // A) Teste et corrige les bornes en offset
  if ( OffsetMin == -10000 ) OffsetMin = this->FitTimeMin();
  if ( OffsetMax == -10000 ) OffsetMax = this->FitTimeMax();
  if ( OffsetStep == 0 ) OffsetStep = this->FitTimeStep();
  if (aTemplate->Pretrigger()!=fPretrigger || aTemplate->TraceLength()!=fTraceLength) {
    cerr <<"WienerLoop : template and pulse have different sizes.." <<endl; exit(-1);
  }
  // Ne pas aller au-dela de la position du pulse!
  if (OffsetMin <= -fPretrigger) OffsetMin=-fPretrigger+1;
  Int_t decalage=fTraceLength-fPretrigger;
  if (!fIsHeat) decalage -= RISE_WINDOW_ION;
  if (OffsetMax >= decalage) OffsetMax = decalage-1;
  // Rester dans la zone ou les offsetfft ont ete calcules!
  Int_t lMin = aTemplate->MinOffset();
  Int_t lMax = aTemplate->MaxOffset();
  if (OffsetMin < lMin) OffsetMin = lMin;
  if (OffsetMax < lMin) OffsetMax = lMin;
  if (OffsetMin > lMax) OffsetMin = lMax;
  if (OffsetMax > lMax) OffsetMax = lMax;

  // B) Boucle sur les offset
  vector<Float_t> Output(4,0);
  vector<Float_t> lFFT, kernel;
  Short_t SignChange=fSign*aTemplate->Sign(); // -1, 1 ou 0 (si signe du pulse pas defini..)
  for (Float_t theOffset=(Float_t)OffsetMin; theOffset<=(Float_t)OffsetMax; theOffset+=OffsetStep) {
    vector<Float_t> params;
    if (fabs(floor(OffsetStep)-OffsetStep)<=0.001) { // No subbin accuracy
      if (!aFast) {
	lFFT=aTemplate->GetOffsetFFT((Int_t)theOffset);
	params = this->ComputeWienerAmpl(lFFT,aNoise);
      } else {
	kernel=aTemplate->GetKernel((Int_t)theOffset);
	Float_t denom=aTemplate->GetDenom((Int_t)theOffset);
	params = this->ComputeWienerFast(kernel,denom);
      }
    } else { // Subbin accuracy routine
      lFFT=aTemplate->GetNonIntegerOffsetFFT(theOffset);
      params = this->ComputeWienerAmpl(lFFT,aNoise);
    }
    Float_t theAmpl = params[0];
    if ( (SignChange==1 && theAmpl>=Output[0])||(SignChange==-1 && theAmpl<=Output[0])||Output[0]==0
	 || (SignChange==0 && (fabs(theAmpl)>= fabs(Output[0]))) ) {
      Output[0] = theAmpl;
      Output[1] = theOffset;
      Output[2] = params[2];
      Output[3] = params[3];
    }
  }
  // On termine en calculant le chi2..
  if (fabs(floor(Output[1])-Output[1])<=0.001) lFFT=aTemplate->GetOffsetFFT((Int_t)Output[1]);
  else lFFT=aTemplate->GetNonIntegerOffsetFFT(Output[1]);
  vector<Float_t> toto = this->ComputeWienerAmpl(lFFT,aNoise,1); // chi2 switch ON
  Output[2]=toto[2];
  // .. Et convention pour avoir une amplitude positive!
  if (SignChange==-1) Output[0] *= -1;
  // .. Et convention : prise en compte de la modulation
  Output[1] *= fModulationLength;
  // ... Convention suite : prise en compte voie "slow" a la mains...
  // (principe : tous les temps en base 100kHz...)
  if (fChannel.find("slow")!=string::npos) Output[1] *= 100;

  return Output;
}

