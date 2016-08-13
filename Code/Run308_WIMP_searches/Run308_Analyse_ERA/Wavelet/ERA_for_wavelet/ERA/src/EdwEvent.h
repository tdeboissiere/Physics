
#ifndef _EDWEVENT_H_
#define _EDWEVENT_H_

#include "EdwUtils.h"
#include "EdwPulse.h"

class EdwEvent : public TObject {
 public:
  EdwEvent() ;
  ~EdwEvent() ;

  UInt_t NbPulses() const { return fPulses.size(); } 
  EdwPulse* Pulse(Int_t num) { return &(fPulses.at(num)); }
  /**< Returns a given pulse from its entry number */
  EdwPulse* Pulse(string aChannel);
  /**< Returns a pulse from its channel name. If no pulse with such a name is found in the event, the NULL pointer is returned. */  
  void AddPulse(EdwPulse aPulse) { fPulses.push_back(aPulse); }
  /**< Adds an EdwPulse in the TObjArray structure */

  void Dump(string aFile="None", string aChannel="all") ; /**< Prints header and all traces to standard output */
  void Clear() ; /**< Resets the event */

  string Run() const { return fRun; }
  string Partition() const { return fPartition; }
  UInt_t SambaNum() const { return fSambaNum; }
  ULong64_t DateSec() const { return fDateSec; } 
  Int_t DateMuSec() const { return fDateMuSec; }
  ULong64_t TimeStamp() const { return fTimeStamp; }
  Int_t GigaStamp() const { return fGigaStamp; }
  Float_t Temperature() const { return fTemperature; }
  UInt_t TriggerBit(UInt_t num) const ;
  vector<UInt_t> TriggerBit() const { return fTriggerBit; }
  void SetSambaNum(UInt_t num) { fSambaNum=num; }
  void SetDateSec(ULong64_t sec) { fDateSec = sec; }
  void SetDateMuSec(Int_t musec) { fDateMuSec = musec; }
  void SetRun(string run) { fRun=run; }
  void SetPartition(string aPartition) {fPartition=aPartition;}
  void SetGigaStamp(Int_t gigast) { fGigaStamp = gigast; }
  void SetTimeStamp(ULong64_t stamp) { fTimeStamp = stamp; }
  void SetTriggerBit(UInt_t aNum, UInt_t aBit) ;
  void SetTemperature(Float_t tempe) { fTemperature=tempe; }

 private:
  vector<EdwPulse> fPulses ; /**< Array of pulses (each event may have an arbitrary number of pulses) */
  string fRun ; /**< Run name */
  string fPartition;
  UInt_t fSambaNum ; /**< Samba event num */
  ULong64_t fDateSec ; /**< Event unix time */
  Int_t fDateMuSec; /**< Microsec time */
  Int_t fGigaStamp;
  ULong64_t fTimeStamp;
  vector<UInt_t> fTriggerBit; /**< Bit trigger encoding */
  Float_t fTemperature;

  ClassDef(EdwEvent,1)
};

#endif

