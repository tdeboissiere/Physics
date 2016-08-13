
#include "EdwEvent.h"

ClassImp(EdwEvent); /**< Root dictionnary */

EdwEvent::EdwEvent() : TObject() {
  fPulses.clear();
  fRun="NONE";
  fSambaNum = 0 ;
  this->SetDateSec(0);
  this->SetDateMuSec(0);
  fTemperature=0;
  
}

EdwEvent::~EdwEvent() {

}

void EdwEvent::Clear() {
  fPulses.clear();
  fRun="NONE";
  fPartition="NONE";
  fSambaNum = 0 ;
  fDateSec=0; fDateMuSec=0;
  fTriggerBit.clear();
  fTimeStamp=0; fGigaStamp=0;
}

UInt_t EdwEvent::TriggerBit(UInt_t num) const {
  UInt_t bit = 0;
  if (num >= fTriggerBit.size()) cerr << "TriggerBit: wrong bit trigger asked"<<endl;
  else bit = fTriggerBit[num];
  return bit;
}

void EdwEvent::SetTriggerBit(UInt_t aNum, UInt_t aBit) {
  if (aNum >= fTriggerBit.size()) fTriggerBit.resize(aNum+1,0);
  fTriggerBit[aNum]=aBit;
}

EdwPulse* EdwEvent::Pulse(string aChannel) {

  EdwPulse* lPulse = NULL;
  for (unsigned int i=0; i<this->NbPulses(); i++) {
    if ( this->Pulse(i)->Channel() == aChannel) {
      if (this->Pulse(i)->TraceLength() != 0) // pulse de longueur nulle donne rien
        lPulse = this->Pulse(i);
    }
  }
  // (we don't check that there could be several pulses with same channel..)
  // When using this function must check if it gives NULL!
  return lPulse;
}

void EdwEvent::Dump(string aFile, string aChannel) {
  if (aFile == "None") {
    cout << "----------------------------------------" << endl;
    /*    cout << "-- EDW Raw Event " << fHeader->Num() << " from run " << fHeader->Run() << endl;
    cout << "Samba delay: " << fHeader->SambaDelay() << endl;
    cout << "Date(s): " << fHeader->DateSec() << " - Date(mus): " << fHeader->DateMuSec() << endl;
    for (UInt_t i=0; i<this->NbPulses(); i++) {
      EdwPulse* lPulse = this->Pulse((Int_t)i);
      if (aChannel == "all" || aChannel == lPulse->Channel()) {
        cout << "------ Channel " << lPulse->Channel()<<" ("<<lPulse->TraceSize()<<" points) ------"<<endl;
        for (UInt_t j=0; j<(UInt_t)(lPulse->TraceSize()); j++) cout << lPulse->Trace(j) << " ";
        cout << endl;
      }
      }*/
    cout << "----------------------------------------" << endl;
  } else {
    ofstream lFile(aFile.c_str(),ios::out);
    lFile << "----------------------------------------" << endl;
    /*   lFile << "-- EDW Raw Event " << fHeader->Num() << " from run " << fHeader->Run() << endl;
    lFile << "Date(s): " << fHeader->DateSec() << " - Date(mus): " << fHeader->DateMuSec() << endl;
    for (UInt_t i=0; i<this->NbPulses(); i++) {
      EdwPulse* lPulse = this->Pulse((Int_t)i);
      if (aChannel == "all" || aChannel == lPulse->Channel()) {
        lFile << "------ Channel " << lPulse->Channel()<<endl;
        for (UInt_t j=0; j<(UInt_t)(lPulse->TraceSize()); j++) 
          lFile << lPulse->Trace(j) << endl;
      }
      }*/
    lFile << "----------------------------------------" << endl;
    lFile.close();
  }

}
