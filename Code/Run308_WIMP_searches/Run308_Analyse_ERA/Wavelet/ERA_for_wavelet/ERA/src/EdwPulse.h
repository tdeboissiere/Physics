
#ifndef _EDWPULSE_H_
#define _EDWPULSE_H_
#include "EdwUtils.h"

class EdwPulse : public TObject {

 public:
  EdwPulse();
  ~EdwPulse();
  void Clear() { fTrace.clear(); }
  string Channel() const { return fChannel; } /**< Channel */
  vector<Short_t> Trace() const { return fTrace; } /**< Raw trace (short type) */
  UInt_t TraceLength() const { return fTraceLength; } /**< Trace size */
  Short_t Trace(Short_t pt) const { return fTrace.at((Int_t)pt); } /**< Raw trace */
  
  Int_t Pretrigger() const { return fPretrigger; }
  Float_t Sampling_ms() const { return fSampling_ms; } /**< Sampling frequency */
  void SetChannel(string aChannel); /**< Channel */
  void SetSampling_ms(Float_t aSampling) { fSampling_ms = aSampling; } /**< Sampling frequency */
  void SetPretrigger(Int_t pretrig) { fPretrigger = pretrig; }
  void SetTrace(const vector<Short_t> & aTrace) { fTrace = aTrace; fTraceLength=aTrace.size(); } /**< Trace */
  void SetTrace(UShort_t n, const Short_t* aData) ; /**< Trace */

  Bool_t IsSaturated() ; /**< Returns 1 if the trace (of short type) reaches saturation */
  Bool_t IsHeat() const { return fIsHeat; } 
  Bool_t IsBBv2() const { return fIsBBv2; } 
  Bool_t IsTrig() const { return fIsTrig; }

  void SetSign(Short_t aSign) ;
  Short_t Sign() const { return fSign; }
  void SetPatternLength(UShort_t aPattLen) {fPatternLength=aPattLen;}
  UShort_t PatternLength() const { return fPatternLength; }
  void SetModulationLength(UShort_t aMod) { fModulationLength=aMod; }
  UShort_t ModulationLength() const { return fModulationLength; }
  void SetBBv2Flag(Bool_t isbbv2) { fIsBBv2=isbbv2; } 
  void SetTriggerFlag(Bool_t istrig) { fIsTrig=istrig; }

 protected:
  string fChannel; /**< Channel */
  Bool_t fIsHeat;
  Bool_t fIsBBv2; /**< Specific to ion: RC vs Heavyside shape */
  Bool_t fIsTrig;
  Float_t fSampling_ms ; /**< Sampling frequency */
  Int_t fPretrigger;
  vector<Short_t> fTrace; /**< Raw trace */
  UInt_t fTraceLength;
  UShort_t fPatternLength; /**< Pattern */
  Short_t fSign; /**< Pulse sign (defined a priori) */
  UShort_t fModulationLength;

  ClassDef(EdwPulse,1)
};

#endif
