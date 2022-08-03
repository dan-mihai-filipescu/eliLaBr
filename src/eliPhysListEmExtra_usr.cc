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
/// \file eli/src/eliPhysListEmExtra_usr.cc
/// \brief Implementation of the eliPhysListEmExtra_usr class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "eliPhysListEmExtra_usr.hh"

#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"

#include "G4PhotoNuclearProcess.hh"
#include "G4CascadeInterface.hh"
#include "G4GammaNuclearReaction.hh"
#include "G4ElectronNuclearProcess.hh"
#include "G4PositronNuclearProcess.hh"
#include "G4ElectroVDNuclearModel.hh"

// gamma

// e-

// e+

eliPhysListEmExtra_usr::eliPhysListEmExtra_usr(const G4String& name)
    :  G4VPhysicsConstructor(name) {
}

eliPhysListEmExtra_usr::~eliPhysListEmExtra_usr() {
}

/*void eliPhysListEmExtra_usr::ConstructParticle() {
	// In this method, static member functions should be called
	// for all particles which you want to use.
	// This ensures that objects of these particle types will be
	// created in the program.

	G4Geantino::GeantinoDefinition();
	G4Gamma::GammaDefinition();
	G4Electron::ElectronDefinition();
	G4Positron::PositronDefinition();

}*/

void eliPhysListEmExtra_usr::ConstructProcess() {

    G4ProcessManager * pManager = 0;

  //Add photonuclear reactions
    pManager = G4Gamma::Gamma()->GetProcessManager();
    G4PhotoNuclearProcess* process = new G4PhotoNuclearProcess();
  //  process->BiasCrossSectionByFactor(1000);
    G4CascadeInterface* bertini = new G4CascadeInterface();
    bertini->SetMaxEnergy(10*GeV);
    process->RegisterMe(bertini);
  //  G4GammaNuclearReaction* gammaNuclear = new G4GammaNuclearReaction();
  //  process->RegisterMe(gammaNuclear);
    pManager->AddDiscreteProcess(process);

  //Add electron nuclear reactions
    pManager = G4Electron::Electron()->GetProcessManager();
    G4ElectronNuclearProcess * theElectronNuclearProcess = new G4ElectronNuclearProcess;
    G4ElectroVDNuclearModel * theElectroReaction = new G4ElectroVDNuclearModel;
    theElectronNuclearProcess->RegisterMe(theElectroReaction);
    pManager->AddDiscreteProcess(theElectronNuclearProcess);

  //Add positron nuclear reactions
    pManager = G4Positron::Positron()->GetProcessManager();
    G4PositronNuclearProcess * thePositronNuclearProcess = new G4PositronNuclearProcess;
    thePositronNuclearProcess->RegisterMe(theElectroReaction);
    pManager->AddDiscreteProcess(thePositronNuclearProcess);

}
