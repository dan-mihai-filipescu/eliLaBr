//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb  6 11:04:50 2017 by ROOT version 5.34/18
// from TTree eliLaBr/eliLaBr
// found on file: nai_9_MeV.root
//////////////////////////////////////////////////////////

#ifndef PileItUp_h
#define PileItUp_h

// Fixed size dimensions of array or collections stored in the TTree if any.
#define NParam 2
#define NDetTypes 4
#define NDetTypes_indv 1
#define NMAXPOINTS 1000
#define NMAXDETECTORS 40
#define BUFFER 20000
//#define PACKCHAN 5. //microseconds to pile-up events
//#define AVGTIMECONTRACT false  //Average the times or not when contracting energy deposition over time spectrum - PileUp)

// include ALL header files needed
//#ifndef __CINT__
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTree.h>
#include <TClassTable.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>

#include "TApplication.h"
#include "TSystem.h"
#include "TClass.h"

#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <math.h>

// Header file for the classes stored in the TTree if any.
#include "vectors.h"

class TSTLvector;
class Tvector;

typedef struct _GaspRecHeader {
  short head[16];
} GaspRecHeader;

typedef struct _GaspEvent {
  //short elem[(NSignals)*(NParam+1)+NDetTypes+2];
  short *elem;
  short evSize;
} GaspEvent;


class PileItUp : public TSelector {

private:

int iGasp, j, z, ii, iChan;
int z_indv, mult_indv;

Int_t NHistMax;
Double_t EHistMax, THistMax;
Double_t EneDivFact, TimeMultFact;
Int_t PartType;
Int_t Option_Time;
Int_t Option_Record;
Int_t Option_OutFiles;
float PackChan;
int AvgTimeContract;
int HistEneTime;

FILE *fp, *fp_indv; 		//  OUTPUT to GASP RunFile (total and individual)
char printfile[50];

int iRecord, iRecord_indv, runnum;

short recSize, recSize_indv;

GaspRecHeader *grh, *grh_indv;
GaspEvent *G, *G_indv;
short *zero;

int Ndet;
int NdetA, NdetB, NdetC;

TH1F* vESingle[NMAXDETECTORS];
TH1F* vEPileUp[NMAXDETECTORS];
TH1F* vESingleInc[NMAXDETECTORS];
TH1F* vTime[NMAXDETECTORS];
TH1F* HistS;
TH1F* HistSp;
TH1F* HistP;
TH1F* HistT;
TH1F* CHistS;
TH1F* CHistSp;
TH1F* CHistP;
TH1F* CHistT;
TH1F* BHistS;
TH1F* BHistSp;
TH1F* BHistP;
TH1F* BHistT;
TH1F* AHistS;
TH1F* AHistSp;
TH1F* AHistP;
TH1F* AHistT;
TH1F* TargetP;
TH2F* ET_TotalProj;

std::stringstream ss, ps;

Tvectors* VectorOfPeaks_Tot;
Tvectors* VectorOfPeaks_R1;
Tvectors* VectorOfPeaks_R2;
Tvectors* VectorOfPeaks_R3;
Tvectors* VectorOfPeaks;
TSTLvector* PileUpVector;
Tvectors* EmptyVector;
Tvectors* PileUpTarget;

//Double_t sum;
double prev_ev, this_ev; 
int ev_index, counter, ev, mult0, Icnt, Jcnt;
float PileUpInp, PoissonAvg, Energy;

public:
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           evNo;
   Int_t           nPoints;
   Double_t        Eini;
   Double_t        Eini_target;
   Double_t        E_target;
   Double_t        time_ini;
   Double_t        tg_time;
   Int_t           mult;
   Int_t           mult_hit;
   Double_t        ev_weight;
   Int_t           det_nb[NMAXDETECTORS];   //[mult]
   Int_t           det_hit[NMAXDETECTORS];   //[mult_hit]
   Double_t        E_br[NMAXDETECTORS];   //[mult]
   Double_t        E_hit[NMAXDETECTORS];   //[mult_hit]
   Double_t        G_time[NMAXDETECTORS];   //[mult]
   Int_t           part_hit[NMAXDETECTORS];   //[mult_hit]
   Double_t        Ene[NMAXDETECTORS],sum[NMAXDETECTORS];
   Bool_t          det_status[NMAXDETECTORS];
   TSTLvector      *PeaksInCountersBranch;
   Tvectors        *VectorOfTargetBranch;

   // List of branches
   TBranch        *b_evNo;   //!
   TBranch        *b_nPoints;   //!
   TBranch        *b_Eini;   //!
   TBranch        *b_Eini_target;   //!
   TBranch        *b_E_target;   //!
   TBranch        *b_time_ini;   //!
   TBranch        *b_tg_time;   //!
   TBranch        *b_mult;   //!
   TBranch        *b_mult_hit;   //!
   TBranch        *b_ev_weight;   //!
   TBranch        *b_det_nb;   //!
   TBranch        *b_det_hit;   //!
   TBranch        *b_E_br;   //!
   TBranch        *b_E_hit;   //!
   TBranch        *b_G_time;   //!
   TBranch        *b_part_hit;   //!
   TBranch        *b_PeaksInCountersBranch;   //!
   TBranch        *b_VectorOfTargetBranch;    //!

   PileItUp(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~PileItUp() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(PileItUp,0);
};

#endif

#ifdef PileItUp_cxx
void PileItUp::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(0);
   PeaksInCountersBranch = new TSTLvector();
   VectorOfTargetBranch = new Tvectors();

   fChain->SetBranchAddress("evNo", &evNo, &b_evNo);
   fChain->SetBranchAddress("nPoints", &nPoints, &b_nPoints);
   fChain->SetBranchAddress("Eini", &Eini, &b_Eini);
   fChain->SetBranchAddress("Eini_target", &Eini_target, &b_Eini_target);
   fChain->SetBranchAddress("E_target", &E_target, &b_E_target);
   fChain->SetBranchAddress("time_ini", &time_ini, &b_time_ini);
   fChain->SetBranchAddress("tg_time", &tg_time, &b_tg_time);
   fChain->SetBranchAddress("mult", &mult, &b_mult);
   fChain->SetBranchAddress("mult_hit", &mult_hit, &b_mult_hit);
   fChain->SetBranchAddress("event_weight", &ev_weight, &b_ev_weight);
   fChain->SetBranchAddress("det_nb", det_nb, &b_det_nb);
   fChain->SetBranchAddress("det_hit", det_hit, &b_det_hit);
   fChain->SetBranchAddress("E_br", E_br, &b_E_br);
   fChain->SetBranchAddress("E_hit", E_hit, &b_E_hit);
   fChain->SetBranchAddress("G_time", G_time, &b_G_time);
   fChain->SetBranchAddress("part_hit", part_hit, &b_part_hit);
   fChain->SetBranchAddress("PeaksInCountersBranch",&PeaksInCountersBranch,&b_PeaksInCountersBranch);
   fChain->SetBranchAddress("VectorOfTargetBranch",&VectorOfTargetBranch,&b_VectorOfTargetBranch);

}

Bool_t PileItUp::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef PileItUp_cxx
