
#ifndef eliLaBr_TrackerSD_h
#define eliLaBr_TrackerSD_h 1

#include "globals.hh"
#include "G4VSensitiveDetector.hh"
#include "eliLaBr_TrackerHit.hh"
#include "G4ios.hh"

class G4Step;
class G4HCofThisEvent;


class eliLaBr_TrackerSD : public G4VSensitiveDetector
{
  public:
      eliLaBr_TrackerSD(G4String);
     ~eliLaBr_TrackerSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*, G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);

  private:
      eliLaBr_TrackerHitsCollection* trackerCollection;

};


#endif
