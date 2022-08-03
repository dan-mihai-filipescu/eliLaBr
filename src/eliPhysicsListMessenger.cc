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
/// \file eli/src/eliPhysicsListMessenger.cc
/// \brief Implementation of the eliPhysicsListMessenger class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "eliPhysicsListMessenger.hh"

#include "eliPhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

eliPhysicsListMessenger::eliPhysicsListMessenger(eliPhysicsList* pPhys)
:fPPhysicsList(pPhys)
{   
  fPhysDir = new G4UIdirectory("/eli/phys/");
  fPhysDir->SetGuidance("physics control.");

  fGammaCutCmd = new G4UIcmdWithADoubleAndUnit("/eli/phys/setGammaCut",this);
  fGammaCutCmd->SetGuidance("Set gamma cut.");
  fGammaCutCmd->SetParameterName("Gcut",false);
  fGammaCutCmd->SetUnitCategory("Length");
  fGammaCutCmd->SetRange("Gcut>0.0");
  fGammaCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fElectCutCmd = new G4UIcmdWithADoubleAndUnit("/eli/phys/setElectronCut",this);
  fElectCutCmd->SetGuidance("Set electron cut.");
  fElectCutCmd->SetParameterName("Ecut",false);
  fElectCutCmd->SetUnitCategory("Length");
  fElectCutCmd->SetRange("Ecut>0.0");
  fElectCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fPosiCutCmd = new G4UIcmdWithADoubleAndUnit("/eli/phys/setPositronCut",this);
  fPosiCutCmd->SetGuidance("Set positron cut.");
  fPosiCutCmd->SetParameterName("Pcut",false);
  fPosiCutCmd->SetUnitCategory("Length");
  fPosiCutCmd->SetRange("Pcut>0.0");
  fPosiCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fProtoCutCmd = new G4UIcmdWithADoubleAndUnit("/eli/phys/setProtonCut",this);
  fProtoCutCmd->SetGuidance("Set proton cut.");
  fProtoCutCmd->SetParameterName("ProtCut",false);
  fProtoCutCmd->SetUnitCategory("Length");
  fProtoCutCmd->SetRange("ProtCut>0.0");
  fProtoCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fNeutCutCmd = new G4UIcmdWithADoubleAndUnit("/eli/phys/setNeutronCut",this);
  fNeutCutCmd->SetGuidance("Set neutron cut.");
  fNeutCutCmd->SetParameterName("NeutCut",false);
  fNeutCutCmd->SetUnitCategory("Length");
  fNeutCutCmd->SetRange("NeutCut>0.0");
  fNeutCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAllCutCmd = new G4UIcmdWithADoubleAndUnit("/eli/phys/setCuts",this);
  fAllCutCmd->SetGuidance("Set cut for all.");
  fAllCutCmd->SetParameterName("cut",false);
  fAllCutCmd->SetUnitCategory("Length");
  fAllCutCmd->SetRange("cut>0.0");
  fAllCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fECutCmd = new G4UIcmdWithADoubleAndUnit("/eli/phys/TargetCuts",this);
  fECutCmd->SetGuidance("Set cuts for the target");
  fECutCmd->SetParameterName("Ecut",false);
  fECutCmd->SetUnitCategory("Length");
  fECutCmd->SetRange("Ecut>0.0");
  fECutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fMCutCmd = new G4UIcmdWithADoubleAndUnit("/eli/phys/DetectorCuts",this);
  fMCutCmd->SetGuidance("Set cuts for the Detector");
  fMCutCmd->SetParameterName("Mcut",false);
  fMCutCmd->SetUnitCategory("Length");
  fMCutCmd->SetRange("Mcut>0.0");
  fMCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fPListCmd = new G4UIcmdWithAString("/eli/phys/SelectPhysics",this);
  fPListCmd->SetGuidance("Select modular physics list.");
  fPListCmd->SetParameterName("PList",false);
  fPListCmd->AvailableForStates(G4State_PreInit);

  listCmd = new G4UIcmdWithoutParameter("/eli/phys/ListPhysics",this);
  listCmd->SetGuidance("Available Physics Lists");
  listCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

eliPhysicsListMessenger::~eliPhysicsListMessenger()
{
  delete fPhysDir;
  delete fGammaCutCmd;
  delete fElectCutCmd;
  delete fPosiCutCmd;
  delete fProtoCutCmd;
  delete fNeutCutCmd;
  delete fAllCutCmd;
  delete fECutCmd;
  delete fMCutCmd;
  delete fPListCmd;
  delete listCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eliPhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{
  if( command == fGammaCutCmd )
   { fPPhysicsList->SetCutForGamma(fGammaCutCmd->GetNewDoubleValue(newValue));}

  if( command == fElectCutCmd )
   { fPPhysicsList->SetCutForElectron(fElectCutCmd->GetNewDoubleValue(newValue));}

  if( command == fPosiCutCmd )
   { fPPhysicsList->SetCutForPositron(fPosiCutCmd->GetNewDoubleValue(newValue));}

  if( command == fProtoCutCmd )
   { fPPhysicsList->SetCutForProton(fProtoCutCmd->GetNewDoubleValue(newValue));}

  if( command == fNeutCutCmd )
   { fPPhysicsList->SetCutForNeutron(fNeutCutCmd->GetNewDoubleValue(newValue));}

  if( command == fAllCutCmd )
    {
      G4double cut = fAllCutCmd->GetNewDoubleValue(newValue);
      fPPhysicsList->SetCutForGamma(cut);
      fPPhysicsList->SetCutForElectron(cut);
      fPPhysicsList->SetCutForPositron(cut);
      fPPhysicsList->SetCutForProton(cut);
      fPPhysicsList->SetCutForNeutron(cut);
    }

  if( command == fECutCmd )
   { fPPhysicsList->SetTargetCut(fECutCmd->GetNewDoubleValue(newValue));}

  if( command == fMCutCmd )
   { fPPhysicsList->SetDetectorCut(fMCutCmd->GetNewDoubleValue(newValue));}

  /*if( command == fPListCmd )
   { fPPhysicsList->SelectPhysicsList(newValue);}*/

  if( command == fPListCmd ) {
      if(fPPhysicsList) {
        G4String name = newValue;
        if(name == "PHYSLIST") {
      char* path = getenv(name);
      if (path) name = G4String(path);
      else {
        G4cout << "### PhysicsListMessenger WARNING: "
           << " environment variable PHYSLIST is not defined"
           << G4endl;
        return;
      }
        }
        fPPhysicsList->SelectPhysicsList(name);
      } else {
        G4cout << "### PhysicsListMessenger WARNING: "
           << " /eli/phys UI command is not available "
           << "for reference Physics List" << G4endl;
      }

    }

  if( command == listCmd ) {
      if(fPPhysicsList) {
        fPPhysicsList->List();
      } else {
        G4cout << "### PhysicsListMessenger WARNING: "
           << " /eli/phys UI command is not available "
           << "for reference Physics List" << G4endl;
      }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
