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
//
// $Id$
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "eliLaBr_GeometryMessenger.hh"

#include "eliLaBr_Geometry.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

eliLaBr_GeometryMessenger::eliLaBr_GeometryMessenger(eliLaBr_Geometry* Geo):Geometry(Geo)
{ 
  eliLaBrDir = new G4UIdirectory("/eli/");
  eliLaBrDir->SetGuidance("UI commands of this example");
  
  geometryDir = new G4UIdirectory("/eli/geometry/");
  geometryDir->SetGuidance("geometry control");

  SetRad1Cmd = new G4UIcmdWithADoubleAndUnit("/eli/geometry/SetRadius1",this);
  SetRad1Cmd->SetGuidance("Set Radius of the Inner Ring");
  SetRad1Cmd->SetParameterName("Size",false);
  SetRad1Cmd->SetRange("Size>=0.");
  SetRad1Cmd->SetUnitCategory("Length");
  SetRad1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetRad2Cmd = new G4UIcmdWithADoubleAndUnit("/eli/geometry/SetRadius2",this);
  SetRad2Cmd->SetGuidance("Set Radius of the Middle Ring");
  SetRad2Cmd->SetParameterName("Size",false);
  SetRad2Cmd->SetRange("Size>=0.");
  SetRad2Cmd->SetUnitCategory("Length");
  SetRad2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetRad3Cmd = new G4UIcmdWithADoubleAndUnit("/eli/geometry/SetRadius3",this);
  SetRad3Cmd->SetGuidance("Set Radius of the Outer Ring");
  SetRad3Cmd->SetParameterName("Size",false);
  SetRad3Cmd->SetRange("Size>=0.");
  SetRad3Cmd->SetUnitCategory("Length");
  SetRad3Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetNum1Cmd = new G4UIcmdWithAnInteger("/eli/geometry/SetNumber1",this);
  SetNum1Cmd->SetGuidance("Set number of Inner Detectors");
  SetNum1Cmd->SetParameterName("NbDet1",false);
  SetNum1Cmd->SetRange("NbDet1>=0");
  SetNum1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetNum2Cmd = new G4UIcmdWithAnInteger("/eli/geometry/SetNumber2",this);
  SetNum2Cmd->SetGuidance("Set number of Middle Detectors");
  SetNum2Cmd->SetParameterName("NbDet2",false);
  SetNum2Cmd->SetRange("NbDet2>=0");
  SetNum2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetNum3Cmd = new G4UIcmdWithAnInteger("/eli/geometry/SetNumber3",this);
  SetNum3Cmd->SetGuidance("Set number of Outer Detectors");
  SetNum3Cmd->SetParameterName("NbDet3",false);
  SetNum3Cmd->SetRange("NbDet3>=0");
  SetNum3Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetPhy1Cmd = new G4UIcmdWithADoubleAndUnit("/eli/geometry/SetPhase1",this);
  SetPhy1Cmd->SetGuidance("Set Angle Phase of the Inner Ring");
  SetPhy1Cmd->SetParameterName("Angle",false);
//  SetPhy1Cmd->SetRange("Angle>=0.");
  SetPhy1Cmd->SetUnitCategory("Angle");
  SetPhy1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetPhy2Cmd = new G4UIcmdWithADoubleAndUnit("/eli/geometry/SetPhase2",this);
  SetPhy2Cmd->SetGuidance("Set Angle Phase of the Middle Ring");
  SetPhy2Cmd->SetParameterName("Angle",false);
//  SetPhy2Cmd->SetRange("Angle>=0.");
  SetPhy2Cmd->SetUnitCategory("Angle");
  SetPhy2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetPhy3Cmd = new G4UIcmdWithADoubleAndUnit("/eli/geometry/SetPhase3",this);
  SetPhy3Cmd->SetGuidance("Set Angle Phase of the Outer Ring");
  SetPhy3Cmd->SetParameterName("Angle",false);
//  SetPhy3Cmd->SetRange("Angle>=0.");
  SetPhy3Cmd->SetUnitCategory("Angle");
  SetPhy3Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

eliLaBr_GeometryMessenger::~eliLaBr_GeometryMessenger()
{
  delete SetRad1Cmd;
  delete SetRad2Cmd;
  delete SetRad3Cmd;
  delete SetNum1Cmd;
  delete SetNum2Cmd;
  delete SetNum3Cmd;
  delete SetPhy1Cmd;
  delete SetPhy2Cmd;
  delete SetPhy3Cmd;
  delete geometryDir;
  delete eliLaBrDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eliLaBr_GeometryMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == SetRad1Cmd )
   { Geometry->SetRad1(SetRad1Cmd->GetNewDoubleValue(newValue));}

  if( command == SetRad2Cmd )
   { Geometry->SetRad2(SetRad2Cmd->GetNewDoubleValue(newValue));}

  if( command == SetRad3Cmd )
   { Geometry->SetRad3(SetRad3Cmd->GetNewDoubleValue(newValue));}

  if( command == SetNum1Cmd )
   { Geometry->SetNum1(SetNum1Cmd->GetNewIntValue(newValue));}

  if( command == SetNum2Cmd )
   { Geometry->SetNum2(SetNum2Cmd->GetNewIntValue(newValue));}

  if( command == SetNum3Cmd )
   { Geometry->SetNum3(SetNum3Cmd->GetNewIntValue(newValue));}

  if( command == SetPhy1Cmd )
   { Geometry->SetPhy1(SetPhy1Cmd->GetNewDoubleValue(newValue));}

  if( command == SetPhy2Cmd )
   { Geometry->SetPhy2(SetPhy2Cmd->GetNewDoubleValue(newValue));}

  if( command == SetPhy3Cmd )
   { Geometry->SetPhy3(SetPhy3Cmd->GetNewDoubleValue(newValue));}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
