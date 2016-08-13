
#include "NoiseSpectrum.h"

ClassImp(NoiseSpectrum);

NoiseSpectrum::NoiseSpectrum() {
  fSpectrum.clear();
  fChannel="NONE";
  fStartValidity=0;
  fEndValidity=0;
  fSampling_ms=0;
}

NoiseSpectrum::~NoiseSpectrum() {
  fSpectrum.clear();
}

NoiseSpectrum::NoiseSpectrum(const string channel, const ULong64_t tinf, const ULong64_t tsup) {
  fChannel=channel;
  fStartValidity=tinf;
  fEndValidity=tsup;
  fSampling_ms=0;
}

// Pour l'instant pas de fonction compliquee a instancier...
// on prefere une classe "nue" comme ca, qu'on va peut-etre completer, a une simple structure c.

