
#ifndef eliLaBr_PrimaryGeneratorAction_h
#define eliLaBr_PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "globals.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "G4RandomDirection.hh"
#include "Randomize.hh"
#include "eliLaBr_GammaSource.hh"
#include "eliLaBr_EventInformation.hh"

#include "TArrayD.h"
#include <TRandom2.h>
#include <TF1.h>
#include "TGraph2D.h"
#include "TH2D.h"
#include "TH1F.h"

#include <stdlib.h>
#include <cmath>

class eliLaBr_DetectorConstruction;
class G4ParticleGun;
class G4Event;


class eliLaBr_PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{

  public:
    eliLaBr_PrimaryGeneratorAction(eliLaBr_DetectorConstruction*);
    ~eliLaBr_PrimaryGeneratorAction();
    G4double ev_Weight, norm;
    G4int n_type;
    G4int NDet;



/*				
//   private: 
	G4int 			n_particle; 
	G4ParticleGun* 		particleGun; 
	G4ParticleTable* 	particleTable; 
	G4String 		particleName;
*/

private:
    eliLaBr_GammaSource *eliLaBr_Source;
    eliLaBr_DetectorConstruction* Detector;
    
  public: 
    void GeneratePrimaries	(G4Event* anEvent); 
	G4int			evNb; 

	TGraph2D*		neutron; 
	TGraph2D*		gamma; 
	TH2D*			th2neutron; 
	TH2D*			th2gamma; 
	TH1F*			gamma_en; 
	int*			tmp; 
	G4double		angle, 
				dummy, 
				energy; 
	G4double		Cos_th,
				theta, 
				phi; 
	G4double                PolAngleIn;
	G4ThreeVector		direction;
	G4ThreeVector		PolVector;
	G4bool			LCSource;

    G4int sign(G4double x);

};

#endif
