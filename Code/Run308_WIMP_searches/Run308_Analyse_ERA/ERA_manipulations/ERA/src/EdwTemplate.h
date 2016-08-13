
#ifndef _EDWTEMPLATE_H_
#define _EDWTEMPLATE_H_

#include "EdwUtils.h"

class EdwTemplate : public TObject {

 public:
  EdwTemplate() ;
  ~EdwTemplate() ;
  void Clear();
  void Initialize(string aChannel, string aType, ULong64_t aTinf, ULong64_t aTsup, vector<Float_t> aCoefs, UInt_t aLength, UInt_t aPretrig, Int_t aOffsetMin, Int_t aOffsetMax, vector<Float_t> aNoise);

  string Channel() const { return fChannel; } /**< Channel */
  string Type() const { return fType; }
  Short_t Sign() const { return fSign; }
  Int_t Pretrigger() const { return fPretrigger; }
  UInt_t TraceLength() const { return fTraceLength; }
  Bool_t IsTimeValid(ULong64_t aTime) const ;


  Int_t MinOffset() const { return fBinsOffsetFFT[0]; }
  Int_t MaxOffset() const { return fBinsOffsetFFT[fBinsOffsetFFT.size()-1]; }
  vector<Float_t> GetOffsetFFT(Int_t aOffset) const ;
  vector<Float_t> GetNonIntegerOffsetFFT(Float_t aOffset) const ;
  vector<Float_t> GetKernel(Int_t aOffset) const ;
  Float_t GetDenom(Int_t aOffset) const ;

  void ComputeKernels(vector<Float_t> aNoise) ;

  vector<Float_t> ComputeTrace(Float_t offset, Bool_t windowflag=0) const ;
  Float_t PulseStandardIon(Float_t x, Float_t offset) const;
  Float_t PulseIonBB2(Float_t x, Float_t offset) const;
  Float_t PulseStandardChal(Float_t x, Float_t offset) const;
  Float_t PulseChalNTD(Float_t x, Float_t offset) const;

 private:
  string fChannel; /**< Template channel */
  string fType;
  ULong64_t fStartValidity; /**< Period of validity */
  ULong64_t fEndValidity; /**< Period of validity */
  Int_t fPretrigger;
  UInt_t fTraceLength;
  Short_t fSign;
  Float_t fNormalisation;
  vector<Float_t> fCoefficients;
  vector<Float_t> fTrace; /**< Template trace*/
  vector<Int_t> fBinsOffsetFFT; /**< List of bins for which the ffts are saved */
  vector< vector<Float_t> > fTracesOffsetFFT;
  vector<Float_t> fDenoms; /**< kernel denominators */
  vector< vector<Float_t> > fKernels; /**< kernels */


  ClassDef(EdwTemplate,1)
};


#endif
