//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file eli/src/eliPhysListEmLivermorePolarized_usr.cc
/// \brief Implementation of the eliPhysListEmLivermorePolarized_usr class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "eliPhysListEmLivermorePolarized_usr.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

// gamma
#include "G4PolarizedPhotoElectricEffect.hh"
#include "G4LivermorePolarizedPhotoElectricModel.hh"

#include "G4PolarizedCompton.hh"
#include "G4LivermorePolarizedComptonModel.hh"

#include "G4PolarizedGammaConversion.hh"
#include "G4LivermorePolarizedGammaConversionModel.hh"

#include "G4RayleighScattering.hh"
#include "G4LivermorePolarizedRayleighModel.hh"

#include "G4PhotoNuclearProcess.hh"
#include "G4CascadeInterface.hh"

// e-
#include "G4eMultipleScattering.hh"
#include "G4UniversalFluctuation.hh"

#include "G4ePolarizedIonisation.hh"
#include "G4LivermoreIonisationModel.hh"

#include "G4ePolarizedBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh"

// e+
#include "G4eplusPolarizedAnnihilation.hh"

// mu

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

// hadrons

#include "G4hMultipleScattering.hh"
#include "G4MscStepLimitType.hh"

#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4alphaIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"

// msc models
#include "G4UrbanMscModel93.hh"
#include "G4WentzelVIModel.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4CoulombScattering.hh"

#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

eliPhysListEmLivermorePolarized_usr::eliPhysListEmLivermorePolarized_usr(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

eliPhysListEmLivermorePolarized_usr::~eliPhysListEmLivermorePolarized_usr()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eliPhysListEmLivermorePolarized_usr::ConstructProcess()
{
  // Add EM Processes from G4EmLivermorePolarizedPhysics builder

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    G4double LivermoreHighEnergyLimit = GeV;
     
    if (particleName == "gamma") {
     
      G4PolarizedPhotoElectricEffect* thePolarizedPhotoElectricEffect = new G4PolarizedPhotoElectricEffect();
      G4LivermorePolarizedPhotoElectricModel* theLivermorePolarizedPhotoElectricModel =
        new G4LivermorePolarizedPhotoElectricModel();
      theLivermorePolarizedPhotoElectricModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      thePolarizedPhotoElectricEffect->AddEmModel(0, theLivermorePolarizedPhotoElectricModel);
      pmanager->AddDiscreteProcess(thePolarizedPhotoElectricEffect);

      G4PolarizedCompton* thePolarizedCompton = new G4PolarizedCompton();
      G4LivermorePolarizedComptonModel* theLivermorePolarizedComptonModel =
        new G4LivermorePolarizedComptonModel();
      theLivermorePolarizedComptonModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      thePolarizedCompton->AddEmModel(0, theLivermorePolarizedComptonModel);
      pmanager->AddDiscreteProcess(thePolarizedCompton);

      G4PolarizedGammaConversion* thePolarizedGammaConversion = new G4PolarizedGammaConversion();
      G4LivermorePolarizedGammaConversionModel* theLivermorePolarizedGammaConversionModel =
        new G4LivermorePolarizedGammaConversionModel();
      theLivermorePolarizedGammaConversionModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      thePolarizedGammaConversion->AddEmModel(0, theLivermorePolarizedGammaConversionModel);
      pmanager->AddDiscreteProcess(thePolarizedGammaConversion);

      G4RayleighScattering* theRayleigh = new G4RayleighScattering();
      G4LivermorePolarizedRayleighModel* thePolarizedRayleighModel = new G4LivermorePolarizedRayleighModel();
      thePolarizedRayleighModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      theRayleigh->AddEmModel(0, thePolarizedRayleighModel);
      pmanager->AddDiscreteProcess(theRayleigh);
      
      /*G4PhotoNuclearProcess* process = new G4PhotoNuclearProcess();
      G4CascadeInterface* bertini = new G4CascadeInterface();
      bertini->SetMaxEnergy(10*GeV);
      process->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(process);*/
      
    } else if (particleName == "e-") {
  
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      //msc->AddEmModel(0, new G4UrbanMscModel93());
      msc->AddEmModel(0, new G4GoudsmitSaundersonMscModel());
      msc->SetStepLimitType(fUseDistanceToBoundary);
      pmanager->AddProcess(msc,                   -1, 1, 1);
      
      // Ionisation
      G4ePolarizedIonisation* ePolarizedIoni = new G4ePolarizedIonisation();
      G4LivermoreIonisationModel* theIoniLivermore = new
        G4LivermoreIonisationModel();
      theIoniLivermore->SetHighEnergyLimit(1*MeV); 
      ePolarizedIoni->AddEmModel(0, theIoniLivermore, new G4UniversalFluctuation() );
      ePolarizedIoni->SetStepFunction(0.2, 100*um); //
      pmanager->AddProcess(ePolarizedIoni,                 -1, 2, 2);
      
      // Bremsstrahlung
      G4ePolarizedBremsstrahlung* ePolarizedBrem = new G4ePolarizedBremsstrahlung();
      G4LivermoreBremsstrahlungModel* theBremLivermore = new
        G4LivermoreBremsstrahlungModel();
      theBremLivermore->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      ePolarizedBrem->AddEmModel(0, theBremLivermore);
      pmanager->AddProcess(ePolarizedBrem, -1,-3, 3);
            
    } else if (particleName == "e+") {

      // Identical to G4EmStandardPhysics_option3
      
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      //msc->AddEmModel(0, new G4UrbanMscModel93());
      msc->AddEmModel(0, new G4GoudsmitSaundersonMscModel());
      msc->SetStepLimitType(fUseDistanceToBoundary);
      pmanager->AddProcess(msc,                   -1, 1, 1);

      G4ePolarizedIonisation* ePolarizedIoni = new G4ePolarizedIonisation();
      ePolarizedIoni->SetStepFunction(0.2, 100*um);

      pmanager->AddProcess(ePolarizedIoni,                 -1, 2, 2);
      pmanager->AddProcess(new G4ePolarizedBremsstrahlung, -1,-3, 3);
      pmanager->AddProcess(new G4eplusPolarizedAnnihilation,0,-1, 4);

      
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {

      // Identical to G4EmStandardPhysics_option3
      
      G4MuMultipleScattering* msc = new G4MuMultipleScattering();
      msc->AddEmModel(0, new G4WentzelVIModel());
      pmanager->AddProcess(msc,                       -1, 1, 1);

      G4MuIonisation* muIoni = new G4MuIonisation();
      muIoni->SetStepFunction(0.2, 50*um);          

      pmanager->AddProcess(muIoni,                    -1, 2, 2);
      pmanager->AddProcess(new G4MuBremsstrahlung,    -1,-3, 3);
      pmanager->AddProcess(new G4MuPairProduction,    -1,-4, 4);
      pmanager->AddDiscreteProcess(new G4CoulombScattering());

    } else if (particleName == "GenericIon") {

      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);

      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetEmModel(new G4IonParametrisedLossModel());
      ionIoni->SetStepFunction(0.1, 10*um);
      pmanager->AddProcess(ionIoni,                   -1, 2, 2);
      pmanager->AddProcess(new G4NuclearStopping(),   -1, 3,-1);

    } else if (particleName == "alpha" ||
               particleName == "He3" ) {

      // Identical to G4EmStandardPhysics_option3
      
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);

      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetStepFunction(0.1, 20*um);
      pmanager->AddProcess(ionIoni,                   -1, 2, 2);
      pmanager->AddProcess(new G4NuclearStopping(),   -1, 3,-1);

    } else if (particleName == "pi+" ||
               particleName == "pi-" ||
               particleName == "kaon+" ||
               particleName == "kaon-" ||
               particleName == "proton" ) {

      // Identical to G4EmStandardPhysics_option3
      
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      
      G4hIonisation* hIoni = new G4hIonisation();
      hIoni->SetStepFunction(0.2, 50*um);

      pmanager->AddProcess(hIoni,                     -1, 2, 2);      
      pmanager->AddProcess(new G4hBremsstrahlung,     -1,-3, 3);
      pmanager->AddProcess(new G4hPairProduction,     -1,-4, 4);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

