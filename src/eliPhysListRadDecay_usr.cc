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
/// \file eli/src/eliPhysListRadDecay_usr.cc
/// \brief Implementation of the eliPhysListRadDecay_usr class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "eliPhysListRadDecay_usr.hh"

#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"

#include "G4PhysicsListHelper.hh"
#include "G4RadioactiveDecay.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"

eliPhysListRadDecay_usr::eliPhysListRadDecay_usr(const G4String& name)
    :  G4VPhysicsConstructor(name) {
}

eliPhysListRadDecay_usr::~eliPhysListRadDecay_usr() {
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

void eliPhysListRadDecay_usr::ConstructProcess() {

    G4RadioactiveDecay* radioactiveDecay = new G4RadioactiveDecay();
    radioactiveDecay->SetHLThreshold(-1.*s);
    radioactiveDecay->SetICM(true);                //Internal Conversion
    radioactiveDecay->SetARM(true);               //Atomic Rearangement

    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());

    // Deexcitation (in case of Atomic Rearangement)
    //
    G4UAtomicDeexcitation* de = new G4UAtomicDeexcitation();
    de->SetFluo(true);
    de->SetAuger(true);
    de->SetPIXE(true);
    G4LossTableManager::Instance()->SetAtomDeexcitation(de);

}
