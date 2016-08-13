
#include "EdwPulse.h"

ClassImp(EdwPulse); /**< Root dictionnary */

EdwPulse::EdwPulse() : TObject() {
  fChannel="NONE";
  fIsHeat=0;
  fIsBBv2=0;
  fIsTrig=0;
  fSampling_ms=0;
  fTrace.clear();
  fPretrigger=0;
  fSign=0;
  fTraceLength=0;
  fPatternLength=0;
  fModulationLength=1; // valeur par defaut = pas de modulation

}

EdwPulse::~EdwPulse() {
  fTrace.clear();
}

void EdwPulse::SetTrace(UShort_t n, const Short_t* aData) {
  fTrace.resize(n,0);
  for (unsigned int i=0; i<n; i++) fTrace[i]=aData[i];
  fTraceLength=fTrace.size();
}

Bool_t EdwPulse::IsSaturated() {
  Bool_t is_sat = 0;
  for (UInt_t i=0; i<fTrace.size(); i++) {
    if (fTrace[i] >= SHRT_MAX || fTrace[i] <= SHRT_MIN) is_sat = 1;
  }
  return is_sat;
}

void EdwPulse::SetChannel(string aChannel) {
  fChannel=aChannel;
  fIsHeat=0;
  if (aChannel.substr(0,7) == "chaleur" || aChannel.substr(0,5)=="chalA" ||
      aChannel.substr(0,5)=="chalB") {
    fIsHeat = 1 ;
//    fSign=-1; // PLUS MAINTENANT !
    fPatternLength=0;
  }
}

void EdwPulse::SetSign(Short_t aSign) {
  if (aSign != -1 && aSign != 1 && aSign != 0) cerr << "EdwPulse::SetSign : wrong sign parameter" << endl;
  fSign = aSign;
}
