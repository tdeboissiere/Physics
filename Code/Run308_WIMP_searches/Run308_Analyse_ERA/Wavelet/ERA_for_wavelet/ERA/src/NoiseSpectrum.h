
#ifndef _NOISESPECTRUM_H_
#define _NOISESPECTRUM_H_

#include "EdwUtils.h"

class NoiseSpectrum : public TObject {

 public:
  NoiseSpectrum();
  ~NoiseSpectrum();
  NoiseSpectrum(const string channel, const ULong64_t tinf, const ULong64_t tsup) ;

  ULong64_t StartValidity() const { return fStartValidity; } /**< Period of validity of these data */
  ULong64_t EndValidity() const { return fEndValidity; } /**< Period of validity of these data */
  void SetStartValidity(ULong64_t aTime) { fStartValidity = aTime; } /**< Period of validity of these data */
  void SetEndValidity(ULong64_t aTime) { fEndValidity = aTime; } /**< Period of validity of these data */
  string Channel() const { return fChannel; } /**< Channel concerned by the data */
  void SetChannel(string aChannel) { fChannel = aChannel; } /**< Channel */
  void SetSpectrum(vector<Float_t> aNoise) { fSpectrum = aNoise; }
  vector<Float_t> Spectrum() const { return fSpectrum; } 
  Float_t Sampling_ms() const { return fSampling_ms; }
  void SetSampling_ms(Float_t aSampling) { fSampling_ms = aSampling; }

 private:
  string fChannel; /**< Template channel */
  ULong64_t fStartValidity; /**< Period of validity */
  ULong64_t fEndValidity; /**< Period of validity */
  Float_t fSampling_ms; /**< f_max/kHz = 1/2*sampling_ms (Shannon) */
  vector<Float_t> fSpectrum;

  ClassDef(NoiseSpectrum,1)
};

#endif
