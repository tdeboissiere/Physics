
#ifndef _EDWUTILS_H_
#define _EDWUTILS_H_

// Standard library includes
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>
#ifndef __CINT__ // bug apparent de rootcint...
#include <dirent.h>
#endif
#include <byteswap.h>
using namespace std;

// Root includes
#include "TStyle.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TTimeStamp.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TFile.h"
#include "TArrayS.h"
#include "TBranch.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TLegend.h"

// Savoir si on met tous les define ici..??
// Pour ceux-la on en a besoin a la fois chez fitpulse et template (qui doivent etre indptds..)
#define PLATEAU_WINDOW_HEAT 50 /**< default plateau for heat window function */
#define RISE_WINDOW_ION 20 /**< default 'smoothing length' for ion window function */


// Mixfft functions:
void factorize(int n, int *nFact, int fact[]);
void transTableSetup(int sofar[], int actual[], int remain[],int *nFact,int *nPoints);
void permute(int nPoint, int nFact,int fact[], int remain[],double xRe[], double xIm[],double yRe[], double yIm[]);
void initTrig(int radix);
void fft_4(double aRe[], double aIm[]);
void fft_5(double aRe[], double aIm[]);
void fft_8();
void fft_10();
void fft_odd(int radix);
void twiddleTransf(int sofarRadix, int radix, int remainRadix,double yRe[], double yIm[]);
void fft(int n, double xRe[], double xIm[],double yRe[], double yIm[]);

// Interface mixfft
vector<Float_t> EdwRealFFT(vector<Float_t> aVect) ;

// Fonctions pratiques pour la manipulation de vecteurs
Float_t VectMean(const vector<Float_t> & aVect) ; 
Float_t VectMin(const vector<Float_t> & aVect) ;
Float_t VectMax(const vector<Float_t> & aVect) ;
Float_t VectMedian(const vector<Float_t> & aVect) ;
Float_t VectRMS(const vector<Float_t> & aVect) ;
Float_t VectIntegral(const vector<Float_t> & aVect, Int_t xmin, Int_t xmax);
Float_t VectIntegral(const vector<Float_t> & aVect) ;
UInt_t VectMinBin(const vector<Float_t> & aVect) ;
UInt_t VectMaxBin(const vector<Float_t> & aVect) ;
void VectMultiply(vector<Float_t> & fVect, const vector<Float_t> & aVect) ;
void VectSquare(vector<Float_t> & fVect);
void VectDivide(vector<Float_t> & fVect, const vector<Float_t> & aVect) ;
void VectAdd(vector<Float_t> & fVect, const vector<Float_t>& aVect, const Float_t c=1);
void VectScale(vector<Float_t> & fVect, const Float_t& aFactor);

// fonction fenetre
vector<Float_t> VectWindowFunction(UInt_t aNbPts, UInt_t aWidth)  ;

// petites structures utiles pour acceder aux fichiers "liste_trucs.txt"...
// Essayer aussi de pouvoir les utiliser en python, tant qu'a faire...!!!

struct BoloStr {
  string Name;
  string Mac;
//  string BB;
  string BB1, BB2;
  string Chal1, Chal2;
  string Col1, Col2, Vet1, Vet2, Gar1, Gar2;
  // Obsolete : numeros de voies recuperees directement pour chaque evt...
//  int NumCol1, NumCol2, NumVet1, NumVet2, NumGar1, NumGar2, NumChal1, NumChal2;
  vector<string> ChTypes, ChNames; 
};

struct RunStr {
  string Name, Type;
  int ModChal1, ModChal2;
  float Vcol1, Vcol2, Vvet1, Vvet2, Vgar1, Vgar2;
  float VrelCol1, VrelCol2,VrelVet1,VrelVet2,VrelGar1,VrelGar2;
};

struct PeriodStr {
  string Run;
  ULong64_t Tinf, Tsup;
};

struct TemplateStr {
  string Channel;
  string Type;
  ULong64_t Tinf, Tsup;
  vector<Float_t> Coefficients;
};

//int GetNumVoie(const string aFile, const string voie);
vector<BoloStr> Read_liste_bolos(const string aFile);
BoloStr Read_liste_bolo(const string aFile, const string aBoloName);
vector<RunStr> Read_liste_runs(const string aFile, const BoloStr aBolo);
RunStr Read_liste_run(const string aFile, const BoloStr aBolo, const string aRunName);
vector<PeriodStr> Read_list_periods(const string aFile, const string aRunName="");
vector<TemplateStr> Read_list_templates(const string aFile);
TemplateStr Read_list_template(string aFile, string aChannel, string aType="NONE", ULong64_t aTime=0);

string GetParam(string aFile, string aParam) ;
bool file_exists(string file);
TChain* ChainFromPartitions(string TraceDir, string Run, string Bolo);

bool DetectTrig(const vector<UInt_t>, const int voienum);

#endif
