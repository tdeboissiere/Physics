#ifndef _FITPULSE_H_
#define _FITPULSE_H_

#include "EdwPulse.h"
#include "EdwTemplate.h"

// Default parameters for pattern length and pulse signs
#define HEAT_PATTERN_LENGTH 0
#define ION_PATTERN_LENGTH 200
#define HEAT_SIGN -1
#define ION_SIGN 1

#define PRETRIGGER_SECURITY 0.9

#define FP_LENGTH 20 /**< Peak finder parameter */
#define FP_NONE -2
#define FP_XWEAK -1 /**< Extra-weak : for very noisy periods */
#define FP_WEAK 0 /**< "garantit" des vrais pulses */
#define FP_STRICT 1 /**< "garantit" des traces sans pulses */
#define FP_SENS_ION_XWEAK 9 /**< Peak finder parameter */
#define FP_ORDER_ION_XWEAK 2 /**< Peak finder parameter */
#define FP_SENS_ION_WEAK 7 /**< Peak finder parameter */
#define FP_ORDER_ION_WEAK 2 /**< Peak finder parameter */
#define FP_SENS_ION_STRICT 5 /**< Peak finder parameter */
#define FP_ORDER_ION_STRICT 3 /**< Peak finder parameter */
#define FP_SENS_HEAT_XWEAK 5.5 /**< Peak finder parameter */
#define FP_ORDER_HEAT_XWEAK 3 /**< Peak finder parameter */
#define FP_SENS_HEAT_WEAK 4.5 /**< Peak finder parameter */
#define FP_ORDER_HEAT_WEAK 3 /**< Peak finder parameter */
#define FP_SENS_HEAT_STRICT 3.5 /**< Peak finder parameter */
#define FP_ORDER_HEAT_STRICT 6 /**< Peak finder parameter */

#define FIT_TIMEMIN_HEAT -5 /**< Default time window for the heat pulse fit */
#define FIT_TIMEMAX_HEAT 5 /**< Default time window for the heat pulse fit */
#define FIT_TIMESTEP_HEAT 0.2
#define FIT_TIMEMIN_ION -1000 /**< Default time window for the ionization pulse fit */
#define FIT_TIMEMAX_ION 600 /**< Default time window for the ionization pulse fit */
#define FIT_TIMEMIN_IONBB2 -50
#define FIT_TIMEMAX_IONBB2 50
#define FIT_TIMESTEP_ION 1

#define PULSESHAPE_SMOOTHING 2

class FitPulse : public EdwPulse {
 public:
  FitPulse() ;
  FitPulse(const EdwPulse*) ;
  FitPulse(const EdwPulse* aPulse1, const EdwPulse* aPulse2, Float_t lambda1=1, Float_t lambda2=1, string comb_name="Combinaison");

  ~FitPulse() ;
  Int_t FitTimeMin() const;
  Int_t FitTimeMax() const;
  Float_t FitTimeStep() const { return (fIsHeat==1)? FIT_TIMESTEP_HEAT : FIT_TIMESTEP_ION ;}
  Bool_t PeakOutsideWindow() ;

  vector<Float_t> ProcessedTrace() const { return fProcessedTrace; }
  Float_t ProcessedTrace(int i) const { return fProcessedTrace.at(i); }
  vector<Float_t> ProcessedTraceFFT() const { return fProcessedTraceFFT; }
  Float_t MeanBase() const { return fMeanBase; } /**< Baseline level */
  vector<Float_t> Pattern() const { return fPattern ; }
  UInt_t NbPeaks() const { return fPeakBins.size() ; }
  vector<Int_t> PeakBins() const { return fPeakBins; }
  Float_t SimpleAmpl() const { return fSimpleAmpl;}
  Int_t SimpleAmplBin() const { return fSimpleAmplBin; }
	
  void BasicPreprocess(Bool_t aCheckPeak=0);
  void RemoveSaturatedSpikes(Int_t maxlength=3) ;
  UInt_t FindGlitches(Int_t glitch_ampl=1000); // param par defaut en dur..

  bool TraceFullySaturated();
  void RemoveBaseline();
  void RemoveLinearBase();
  void RemovePattern(Int_t aNbPatternPts);
  void SetBaseFromPeakBin(Int_t aBin);
  void ApplyWindow(Int_t aWidth=0) ;
  void FindPeaks(Int_t aCriterium, Bool_t AnySign=0, Int_t aLength=FP_LENGTH) ;
  void FindPeaksWithParams(Float_t aSensitivity, Int_t aLength, Int_t aOrder, Bool_t AnySign=0) ;
  void ComputeTraceFFT(Bool_t aWindowFlag=1) ;
  vector<Float_t> SmoothedTrace(UInt_t aSmoothing=PULSESHAPE_SMOOTHING) const ;
  vector<Float_t> SmoothedAmpl(const vector<Float_t>& aSmoothTrace) const ;
  vector<Float_t> GetRiseTimeParams(const Float_t aFracLow=0.2, const Float_t aFracHigh=0.9, UInt_t aSmooth=PULSESHAPE_SMOOTHING) ; /** returns [risetime,t0,a0,a1] */
  Float_t GetRiseTime(const Float_t aFracLow=0.2, const Float_t aFracHigh=0.9, UInt_t aSmooth=PULSESHAPE_SMOOTHING) ;
  vector<Float_t> GetFallTimeParams(const Float_t aFracLow=0.3, const Float_t aFracHigh=0.9, UInt_t aSmooth=PULSESHAPE_SMOOTHING) ; /** returns [falltime,t0,a0,a1] */
  Float_t GetFallTime(const Float_t aFracLow=0.3, const Float_t aFracHigh=0.9, UInt_t aSmooth=PULSESHAPE_SMOOTHING) ;
  vector<Float_t> GetFWHMParams(const Float_t aFracWidth=0.5, UInt_t aSmooth=PULSESHAPE_SMOOTHING) ; /** returns [fwhm,t_start,ampl_at_fwhm] */
  Float_t GetFWHM(const Float_t aFracWidth=0.5, UInt_t aSmooth=PULSESHAPE_SMOOTHING) ;

  vector<Float_t> ComputeWienerAmpl(const vector<Float_t>& aTemplateFFT, const vector<Float_t>& aNoise, Bool_t aChi2Switch=0) ;
  vector<Float_t> ComputeWienerFast(const vector<Float_t>& aKernel, const Float_t aDenom);
  vector<Float_t> WienerLoop(const EdwTemplate* aTemplate, const vector<Float_t>& aNoise, Int_t aOffsetMin=-10000, Int_t aOffsetMax=-10000, Float_t aOffsetStep=0, Bool_t aFast=1);

 private:
  vector<Float_t> fProcessedTrace;
  vector<Float_t> fProcessedTraceFFT;
  vector<Float_t> fPattern; /**< Pattern */

  Bool_t fBasicPreprocessed; /**< Set to 1 if baseline substraction, depatterning.. were already done for that pulse */
  Int_t fBaseStart; /**< Pretrigger */
  Int_t fBaseStop; /**< Pretrigger */
  Bool_t fWindowed;
  
  Float_t fMeanBase; /**< Average value of baseline from the pretrigger region */
  Int_t fSimpleAmplBin; /**< Raw amplitude */
  Float_t fSimpleAmpl; /**< Raw amplitude */
  vector<Int_t> fPeakBins; /**< Position of peaks detected by the raw peakfinder algorithm */

  ClassDef(FitPulse,1)
};

#endif
