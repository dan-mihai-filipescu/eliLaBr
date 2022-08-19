#ifndef eliLaBr_EventAction_h
#define eliLaBr_EventAction_h 1

#define NMAXPOINTS 100000 //ORIGINAL 1000 changed 2021_09_30
#define PACKCHAN 5. //microseconds to pile-up events
#define AVGTIMECONTRACT false  //Average the times or not when contracting energy deposition over time spectrum - PileUp)
#define NRBIN 5

#include "globals.hh"
#include "G4UserEventAction.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4TrajectoryContainer.hh"
#include "G4SteppingManager.hh"
#include "G4SteppingVerbose.hh"
#include "G4SDManager.hh"
#include "G4THitsCollection.hh"
#include "G4HCofThisEvent.hh"
#include "G4Timer.hh"
#include "G4UnitsTable.hh"

#include "TTimeStamp.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TKey.h"
#include "TArrayD.h"
#include "TArrayI.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom2.h"
#include "TSpectrum.h"

#include <fstream> 
#include <sstream>
#include <vector>

#include "eliLaBr_TrackerHit.hh"
#include "eliLaBr_TrackingAction.hh"
#include "eliLaBr_RunAction.hh"
#include "eliLaBr_HistoManager.hh"
#include "eliLaBr_Analysis.hh"
//#include "eliLaBr_TrackerSD.hh"

#include "vectors.h"

using namespace std;

class G4Event;
class Tvectors;
class TSTLvector;

class eliLaBr_EventAction : public G4UserEventAction
{
public:
	eliLaBr_EventAction();
	~eliLaBr_EventAction();

public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
	void SetSaveThreshold(G4int save);
	void SetEventVerbose(G4int v){fVerbose=v;}
#ifdef OLD_ROOT
	void QuickSort(Float_t *xpeaks, int left, int right);
#else
	void QuickSort(Double_t *xpeaks, int left, int right);
#endif

// private:
	G4int				fSaveThreshold;
	G4int				fVerbose;
    G4double            weight, norm, polAngle, e_x;
    G4int               n_type;
    G4bool              Det_touch;
    G4String            Vname;
    G4String            partName;
    G4int N_Det;
	
// scoring...
    enum { Ndet = 40 }; 		// total nb. of detectors
	double Eini; 
//	double angle [Ndet]; 	
        double br_q  [Ndet]; 	// deposited energy inside bromide...
        double br_time_min [Ndet], br_time_max [Ndet]; //interaction time inside bromide...
        bool br_hit_flag [Ndet];   //flag for detector hit
        bool tg_hit_flag;          //flag for target hit
//	bool hist_writed;

//Histograms of energy deposition over time
Int_t NMaxHit;
Int_t NHistMax;
Double_t THistMax;
std::stringstream ss;
TH1D* vEDepinTime [Ndet];
TH1D* vEDepinTimeTarget;
TSpectrum* s;
Int_t nfound, imin, imax, nrbin;
#ifdef OLD_ROOT
Float_t* xpeaks;
#else
Double_t* xpeaks;
#endif
Tvectors* VectorOfPeaks;
Tvectors* VectorOfTarget;
TSTLvector* PeaksInCounters;

// TTree...
    double E_target, Eini_target;
    double time_ini, tg_time_max, tg_time_min, tg_time;  //initial target hit time and interaction times in target
    int mult,mult_hit; 			// event multiplicity
	int mult_tgt; 
    //int mult_br, mult_nai;
    double ev_weight;
	double init_ang; 
    //int clsN, detN, tgtN, prev_tgt;
    int detN;
    int det_nb[Ndet], det_hit[Ndet];
    //int tgt_nb[100];
    //double angle[40];
    double E_br[Ndet], E_hit[Ndet], G_time[Ndet];
    int part_hit[Ndet];
    //double E_nai[40];
    //double E_tgt;
    //double angular[40];
    G4ThreeVector pos;
	int			fEvNo;		
	int			fRunNo;		
	int			fNPoints;	
	TArrayD			fPosX, fPosY, fPosZ; 
    TArrayD			fEdep;
    TArrayD         fKinEne;
    TArrayD			fGlobalTime;
	TArrayI			partType; 
	TArrayI			det_nb_step; 
    G4int			evNb;

// ROOT data saving stuff
    static TFile		*fRootFile;
    static TTree		*fTree;
	TString			fFileName;
};

#endif
