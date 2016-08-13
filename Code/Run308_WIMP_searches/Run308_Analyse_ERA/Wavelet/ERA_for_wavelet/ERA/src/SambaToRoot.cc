
#include "EdwUtils.h"
#include "EdwEvent.h"

enum ByteOrderType { kLittleEndian=0, kBigEndian} ;

ByteOrderType testbyteorder() {
  short word = 0x0001;
  char *byte = (char *) &word;
  return(byte[0] ? kLittleEndian : kBigEndian);
}

Int_t FillRootFile(string aSambaFileName, string aRootFileName, BoloStr aBolo, RunStr aRun) {

  // A) Lecture entete fichier Samba
  string kEnteteEvt = "* Evenement" ; // empirical flags in the raw Edelweiss files
  string kEnteteVoie = "* Voie" ;
  string kSeparator = "* ----------" ;
  short kEnteteEvtSize = kEnteteEvt.size() ;
  short kEnteteVoieSize = kEnteteVoie.size() ;
  short kSeparatorSize = kSeparator.size() ;
  string line, mot1, mot2;
  string version;
  ByteOrderType byte_ordering = kBigEndian; // vieux runs : donnees de samba sans def de l'ordering.. => machines powerpc, big endian
  ByteOrderType machine_ordering = testbyteorder();

  // Get run name, partition, byte ordering
  Int_t length = aSambaFileName.length();
  if (length < 12 || aSambaFileName.compare(length-4,1,"_") )
    cerr << "SambaToRoot: wrong file name" << endl;
  string run_name = aSambaFileName.substr(length-12,8) ;
  string partition = aSambaFileName.substr(length-3,3) ;
  Bool_t IsBBv2 = (aBolo.BB1=="BBv2" || aBolo.BB2!="NONE") ? 1 : 0;
  // NB: si il y a un aBolo.BB2 declare, alors setup edw3 => isbbv2 automatiquement... (que ce soit bb2/3...)
    
  // Teste version du fichier et format (byte ordering)
  ifstream EdwFile(aSambaFileName.c_str(),ios::in) ;
  if (!EdwFile) {
    cerr << "FillRootFile: error opening stream from " << aSambaFileName << endl;
    return 1;
  }
  getline(EdwFile,line);
  if (line!= "* Archive SAMBA") cerr << "FillRootFile: est-ce vraiment une archive samba?" << endl;
  EdwFile >> mot1 >> mot2 >> version;
  if (version!="9.11" && version!="9.29" && version!="8.15" &&
      version!="9.32" && version!="9.33" && version!="9.34" && version!="9.35") 
    cerr << "Version de samba pas verifiee.." << endl;
  while ( line.compare(0,kEnteteEvtSize,kEnteteEvt) ) {
    getline(EdwFile,line);
    if (line.find("Byte-order")!=string::npos) {
      string::size_type eq_pos = line.find("=",0);
      string value = line.substr(eq_pos+1);
      string::size_type diese_pos = line.find("#",0);
      if (diese_pos != string::npos) value = line.substr(eq_pos+1,diese_pos-1-eq_pos) ;
      if (value.find(" little ")!=string::npos) byte_ordering=kLittleEndian;
      else if (value.find(" big ")!=string::npos) byte_ordering=kBigEndian;
      else cerr << "Error reading byte ordering in header " <<value<<"&"<<line<< endl;
    }
    if (!EdwFile.good()) {
      cerr << "FillRootFile: "<< kEnteteEvt<< " not found" << endl;
      return 1;
    }
  }
  
  // B) Creation fichier root
  TFile* RootFile = new TFile(aRootFileName.c_str(),"RECREATE") ;
  if (RootFile->IsZombie()) {
    cerr << "FillRootFile: error opening Root file " << aRootFileName << endl;
    return 1;
  }
  RootFile->cd();
  EdwEvent* Event = new EdwEvent();
  TTree EdwTree("EdwTree","Edelweiss Raw Data Tree");
  EdwTree.SetDirectory(RootFile);
  EdwTree.Branch("Event","EdwEvent",&Event); // default buffer size 32000

  // modulation..
  // Cas de deux chaleurs: on prend la modulation la plus lente
  UShort_t modulation=(aRun.ModChal1>=aRun.ModChal2) ? aRun.ModChal1 : aRun.ModChal2;
  // Float_t reference_sampling=0.01008; // Sampling de reference pour ionisation (100kHz)
  // Il faut faire ca car l'ionsation peut etre moyenneee...

  // C) Boucle sur les evts
  while ( EdwFile.good()) {
    Event->Clear();
    Event->SetRun(run_name);
    Event->SetPartition(partition);
    // C-a) Entete d'evt
    while ( line.compare(0,kEnteteVoieSize,kEnteteVoie) ) {
      getline(EdwFile,line);
      string::size_type diese_pos = line.find("#",0);
      string::size_type eq_pos = line.find("=",0);
      if (eq_pos != string::npos) {
        string field = line.substr(0,eq_pos-1);
        string value = line.substr(eq_pos+1);
        if (diese_pos != string::npos) value = line.substr(eq_pos+1,diese_pos-1-eq_pos) ;
        if (field == "Numero") Event->SetSambaNum( atol(value.c_str()) );
        if (field == "Date.secondes") Event->SetDateSec( atol(value.c_str()) );
        if (field == "Date.microsecs") Event->SetDateMuSec( atoi(value.c_str()) );
        if (field == "Liste:31-0" || field == "Liste:31-00") Event->SetTriggerBit(0,(UInt_t)atoi(value.c_str()));
        if (field == "Liste:63-32") Event->SetTriggerBit(1,(UInt_t)atoi(value.c_str()));
        if (field == "Liste:95-64") Event->SetTriggerBit(2,(UInt_t)atoi(value.c_str()));
        if (field == "Liste:127-96") Event->SetTriggerBit(3,(UInt_t)atoi(value.c_str()));
        if (field == "Liste:159-128") Event->SetTriggerBit(4,(UInt_t)atoi(value.c_str()));
        if (field == "Liste:191-160") Event->SetTriggerBit(5,(UInt_t)atoi(value.c_str()));
        if (field == "Liste:223-192") Event->SetTriggerBit(6,(UInt_t)atoi(value.c_str()));
        if (field == "Liste:255-224") Event->SetTriggerBit(7,(UInt_t)atoi(value.c_str()));
        if (field == "Liste:287-256") Event->SetTriggerBit(8,(UInt_t)atoi(value.c_str()));
        if (field == "Liste:319-288") Event->SetTriggerBit(9,(UInt_t)atoi(value.c_str()));
        if (field == "Liste:351-320") Event->SetTriggerBit(10,(UInt_t)atoi(value.c_str()));
        if (field == "Liste:383-352") Event->SetTriggerBit(11,(UInt_t)atoi(value.c_str()));
        if (field == "GigaStamp") Event->SetGigaStamp( atoi(value.c_str()) );
	    if (field == "TimeStamp") Event->SetTimeStamp( atol(value.c_str()) );
    	if (field == "Temperature" || field=="T_Bolo") {
	    // A partir du run300..: 
        // On teste que c'est pas nul (cas ou T_bolo=0 mais pas Temperature..)
	    Float_t ltempe=(Float_t)atof(value.c_str());
        if (ltempe>0.00001) Event->SetTemperature(ltempe);
	    }
      }
    }
    // C-b) Loop on pulses within a given event
    int nopulse=0;
    
    while (nopulse == 0) {
        EdwPulse pulse;
        pulse.SetChannel(line.substr(8,line.size()-9));
        pulse.SetBBv2Flag(IsBBv2);
        bool pulse_yes=0;
        //   if (!pulse->IsHeat()) cout << pulse->ModulationLength()<<endl;
        //      cout << pulse->Channel()<<" "<<aBolo.Col1<<endl;
        if (pulse.Channel()==aBolo.Col1) {
            pulse.SetSign((aRun.VrelCol1>0 ? -1 : 1));
            // EN DUR: VOIE SLOW ==> PATTERN/=100 (a mettre en param libre...)
            if ( (pulse.Channel()).find("slow")!=string::npos ) pulse.SetPatternLength(modulation/100);
            else pulse.SetPatternLength(modulation);
            pulse_yes=1;
        } else if (pulse.Channel()==aBolo.Col2) {
            pulse.SetSign((aRun.VrelCol2>0 ? -1 : 1));
            if ( (pulse.Channel()).find("slow")!=string::npos ) pulse.SetPatternLength(modulation/100);
            else pulse.SetPatternLength(modulation);
            pulse_yes=1;
        } else if (pulse.Channel()==aBolo.Vet1) {
            pulse.SetSign((aRun.VrelVet1>0 ? -1 : 1));
            if ( (pulse.Channel()).find("slow")!=string::npos ) pulse.SetPatternLength(modulation/100);
            else pulse.SetPatternLength(modulation);
            pulse_yes=1;
        } else if (pulse.Channel()==aBolo.Vet2) {
            pulse.SetSign((aRun.VrelVet2>0 ? -1 : 1));
            if ( (pulse.Channel()).find("slow")!=string::npos ) pulse.SetPatternLength(modulation/100);
            else pulse.SetPatternLength(modulation);
            pulse_yes=1;
        } else if (pulse.Channel()==aBolo.Gar1) {
            pulse.SetSign((aRun.VrelGar1>0 ? -1 : 1));
            if ( (pulse.Channel()).find("slow")!=string::npos ) pulse.SetPatternLength(modulation/100);
            else pulse.SetPatternLength(modulation);
            pulse_yes=1;
        } else if (pulse.Channel()==aBolo.Gar2) {
            pulse.SetSign((aRun.VrelGar2>0 ? -1 : 1));
            if ( (pulse.Channel()).find("slow")!=string::npos ) pulse.SetPatternLength(modulation/100);
            else pulse.SetPatternLength(modulation);
            pulse_yes=1;
            // A VERIFIER!! : signe negatif si BB3 (positif seulement pour BB23)??
        } else if (pulse.Channel()==aBolo.Chal1) {
            pulse.SetSign((aBolo.BB1=="BBv23" ? 1 : -1));
            pulse.SetModulationLength(aRun.ModChal1);
            pulse_yes=1;
        } else if (pulse.Channel()==aBolo.Chal2) {
            pulse.SetSign((aBolo.BB2=="BBv23" ? 1 : -1));
            pulse.SetModulationLength(aRun.ModChal2);
            pulse_yes=1;
        }
        int voienum=-1;
        Short_t pulse_size=0;
        int lPtsFiltre=0;
        while ( line.compare(0,kSeparatorSize,kSeparator) ) { // Pulse header
            getline(EdwFile,line);
            string::size_type diese_pos = line.find("#",0);
            string::size_type eq_pos = line.find("=",0);
            if (eq_pos != string::npos) {
                string field = line.substr(0,eq_pos-1);
                field = field.substr(field.find_first_not_of("\t"));
                string value = line.substr(eq_pos+1);
                if (diese_pos != string::npos) value = line.substr(eq_pos+1,diese_pos-1-eq_pos) ;
                if (field == "Dimension") pulse_size=atoi(value.c_str());
                if (field == "Horloge") pulse.SetSampling_ms( atof(value.c_str()) ) ; // no ms to ns anymore
                if (field == "Filtre.nb") lPtsFiltre=atoi((value).c_str());
                if (field == "Pretrigger") pulse.SetPretrigger(atoi(value.c_str()));
                if (field == "Numero") voienum=atoi((value).c_str());
            } // else cerr << "Strange line for pulse header: " << lLine << endl;
            if (!EdwFile.good()) cerr << "Error reading pulse header" << endl;
        }
        // Rajout : cas de l'ionisation pas a 100kHz!!! On change la modulation - MERDE OBSOLETE.. ?
        //      if (!pulse.IsHeat()) {
        //Float_t rapport=pulse.Sampling_ms()/reference_sampling;
        //if (rapport > 1.2) pulse.SetPatternLength((UShort_t)((Float_t)pulse.PatternLength()/rapport));
        //}
        pulse.SetTriggerFlag(DetectTrig(Event->TriggerBit(),voienum));
        if ( lPtsFiltre ) EdwFile.ignore(8*lPtsFiltre); // filtered data = 8 bytes
        // Finally, check if the pulse is here indeed!
        long lPos = EdwFile.tellg();
        getline(EdwFile,line);
        EdwFile.seekg(lPos);
        if ( (!line.compare(0,kEnteteEvtSize,kEnteteEvt)) || (!line.compare(0,kEnteteVoieSize,kEnteteVoie))) {
            // No pulse, just a header.. We record the pulse in the tree but with no trace.
        } else { // Read the pulse
            Short_t* array = new Short_t[pulse_size];
            EdwFile.read((char*)array,pulse_size*2);
            //if (EdwFile.fail()) cerr << "SambaToRoot::FillRootFile: Error reading a pulse for event " <<Event->SambaNum()<< endl;
            if ((byte_ordering == kBigEndian && machine_ordering == kLittleEndian) || (byte_ordering == kLittleEndian && machine_ordering == kBigEndian)) {
#ifdef DARWIN
                for (Short_t i=0;i<pulse_size;i++) array[i]=ntohs(array[i]);
#else
                for (Short_t i=0;i<pulse_size;i++) array[i]=bswap_16(array[i]);
#endif
            }
            pulse.SetTrace(pulse_size,array);
            delete[] array;
            if (pulse_yes==1) Event->AddPulse(pulse);
            // only if it belongs to the selected bolo
        }
        char toto; // Je comprends pas pourquoi il faut faire ca. Admettons...
        long pos=EdwFile.tellg();
        EdwFile.get(toto);
        if ( EdwFile.good() ) {
            EdwFile.seekg(pos);
            getline(EdwFile,line);
            if ( ! line.compare(0,kEnteteEvtSize,kEnteteEvt) ) nopulse=1;
        } else nopulse=1;
        if (EdwFile.bad()) cerr << "Error reading the edw file.." <<endl;
    }
    RootFile->cd(); EdwTree.Fill();
  }
  RootFile->cd(); EdwTree.Write("",TObject::kOverwrite);
  RootFile->Close() ;
  EdwFile.close();
  return 0;
}

int main(int argc, char* argv[]) {

  // Lecture du fichier de params
  if (argc!=2) {cerr <<"Wrong nb of arguments" << endl; return 1;}
  string paramfile=string(argv[1]);
  string AnaDir=GetParam(paramfile,"AnaDir");
  string SambaDir=GetParam(paramfile,"SambaDir");
  string BoloName=GetParam(paramfile,"BoloName");
  if (AnaDir=="NONE" || SambaDir=="NONE" || BoloName=="NONE") {
    cerr << "wrong param file"<<endl; return 1;
  }
  string BoloFile=GetParam(paramfile,"BoloFile");
  if (BoloFile=="NONE") BoloFile=AnaDir+"/liste_bolos.txt";
  string RunFile=GetParam(paramfile,"RunFile");
  if (RunFile=="NONE") RunFile=AnaDir+"/"+BoloName+"/liste_runs.txt";
  string dum=GetParam(paramfile,"OverWrite");
  Bool_t overwrite=1;
  if (dum!="NONE") overwrite=(Bool_t)(atoi(dum.c_str()));
  string RunName=GetParam(paramfile,"RunName");

  // Recupere le BoloStr et les Runs qu'il faut
  BoloStr bolo=Read_liste_bolo(BoloFile,BoloName); 
  vector<RunStr> listruns;
  if (RunName=="NONE") listruns = Read_liste_runs(RunFile,bolo);
  else listruns.push_back(Read_liste_run(RunFile,bolo,RunName));

  // Boucle sur les runs
  for (UInt_t i_run=0; i_run<listruns.size(); i_run++) {
    RunStr run=listruns[i_run];
    // Recherche des partitions
    DIR* dp = opendir((SambaDir+"/"+run.Name).c_str());
    vector<string> partitions;
    struct dirent *dirp;
    while ((dirp=readdir(dp))!=NULL) {
      string toto=string(dirp->d_name);
      if (toto.find(run.Name+"_0")!=string::npos)
	partitions.push_back(toto);
    }
    if (!partitions.size()) cerr << "No partition for "<<run.Name<<endl;
    // Boucle sur les partitions
    for (UInt_t j_part=0; j_part<partitions.size(); j_part++) {
      string rootfile=AnaDir+"/"+BoloName+"/Traces/"+partitions[j_part]+"_"+BoloName+".root";
      string sambafile=SambaDir+"/"+run.Name+"/"+partitions[j_part];
      if (overwrite || !file_exists(rootfile)) {
	cerr <<"Reading partition "<<partitions[j_part]<<endl;
	FillRootFile(sambafile,rootfile,bolo,run);
      }
    }
  }

  return 0;
}
