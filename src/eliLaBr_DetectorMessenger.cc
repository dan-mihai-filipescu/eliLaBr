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

#include "eliLaBr_DetectorMessenger.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "eliLaBr_DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

eliLaBr_DetectorMessenger::eliLaBr_DetectorMessenger(eliLaBr_DetectorConstruction* Det):Detector(Det)
{ 
  G4UnitDefinition* PersonalUnits = new G4UnitDefinition ("inch" , "inch", "Length", 2.54*cm );
  PersonalUnits -> BuildUnitsTable();

  eliLaBrDir = new G4UIdirectory("/eli/");
  eliLaBrDir->SetGuidance("UI commands of this example");
  
  detDir = new G4UIdirectory("/eli/det/");
  detDir->SetGuidance("detector control");

  E8_FillingCmd = new G4UIcmdWithAString("/eli/det/setE8_Filling",this);
  E8_FillingCmd->SetGuidance("Set the Filling of E8 Hall (Vacuum or Air)");
  E8_FillingCmd->SetParameterName("choice",false);
  E8_FillingCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TGthicknessCmd = new G4UIcmdWithADoubleAndUnit("/eli/det/setTGThickness",this);
  TGthicknessCmd->SetGuidance("Set Thickness of the Target");
  TGthicknessCmd->SetParameterName("Size",false);
  TGthicknessCmd->SetRange("Size>=0.");
  TGthicknessCmd->SetUnitCategory("Length");
  TGthicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TGdiameterCmd = new G4UIcmdWithADoubleAndUnit("/eli/det/setTGDiameter",this);
  TGdiameterCmd->SetGuidance("Set Diameter of the Target");
  TGdiameterCmd->SetParameterName("Size",false);
  TGdiameterCmd->SetRange("Size>=0.");
  TGdiameterCmd->SetUnitCategory("Length");
  TGdiameterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  TGplacementCmd = new G4UIcmdWithAString("/eli/det/setTGPlacement",this);
  TGplacementCmd->SetGuidance("Place Target relative to the WORLD or DETECTOR");
  TGplacementCmd->SetParameterName("choice",false);
  TGplacementCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TGhousingCmd = new G4UIcmdWithABool("/eli/det/setTGHousing",this);
  TGhousingCmd->SetGuidance("USE or NOT Housing for the Target");
  TGhousingCmd->SetParameterName("choice",false);
  TGhousingCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TGpositionCmd = new G4UIcmdWithADoubleAndUnit("/eli/det/setTGPosition",this);
  TGpositionCmd->SetGuidance("Set Position of the Target Placement");
  TGpositionCmd->SetParameterName("Position",false);
//  TGpositionCmd->SetRange("Position>=0.");
  TGpositionCmd->SetUnitCategory("Length");
  TGpositionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TGmaterialCmd = new G4UIcmdWithAString("/eli/det/setTGMaterial",this);
  TGmaterialCmd->SetGuidance("Set Target Material");
  TGmaterialCmd->SetParameterName("choice",false);
  TGmaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ATthicknessCmd = new G4UIcmdWithADoubleAndUnit("/eli/det/setATThickness",this);
  ATthicknessCmd->SetGuidance("Set Thickness of the Attenuator");
  ATthicknessCmd->SetParameterName("Size",false);
  ATthicknessCmd->SetRange("Size>=0.");
  ATthicknessCmd->SetUnitCategory("Length");
  ATthicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ATholeCmd = new G4UIcmdWithADoubleAndUnit("/eli/det/setATHole",this);
  ATholeCmd->SetGuidance("Set hole size in the Attenuator");
  ATholeCmd->SetParameterName("Size",false);
  ATholeCmd->SetRange("(Size>=0.)&&(Size<=5.)");
  ATholeCmd->SetUnitCategory("Length");
  ATholeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ATpositionCmd = new G4UIcmdWithADoubleAndUnit("/eli/det/setATPosition",this);
  ATpositionCmd->SetGuidance("Set position of the Attenuator relative to the Target");
  ATpositionCmd->SetParameterName("Size",false);
  ATpositionCmd->SetRange("Size!=0.");
  ATpositionCmd->SetUnitCategory("Length");
  ATpositionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ATmaterialCmd = new G4UIcmdWithAString("/eli/det/setATMaterial",this);
  ATmaterialCmd->SetGuidance("Set Attenuator Material");
  ATmaterialCmd->SetParameterName("choice",false);
  ATmaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ATplacementCmd = new G4UIcmdWithABool("/eli/det/setATPlacement",this);
  ATplacementCmd->SetGuidance("PLACE or NOT the Attenuator");
  ATplacementCmd->SetParameterName("choice",false);
  ATplacementCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  PlaceWindowCmd = new G4UIcmdWithABool("/eli/det/PlaceWindow",this);
  PlaceWindowCmd->SetGuidance("Place or NOT Borosilicate vacuum window");
  PlaceWindowCmd->SetParameterName("choice",false);
  PlaceWindowCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  PlaceMirrorCmd = new G4UIcmdWithABool("/eli/det/PlaceMirror",this);
  PlaceMirrorCmd->SetGuidance("Place or NOT laser mirror");
  PlaceMirrorCmd->SetParameterName("choice",false);
  PlaceMirrorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  BADVacPresCmd = new G4UIcmdWithADoubleAndUnit("/eli/det/setBADVacuumPressure",this);
  BADVacPresCmd->SetGuidance("Set the Pressure of the Vacuum inside electron beamline");
  BADVacPresCmd->SetParameterName("Pressure",false);
  BADVacPresCmd->SetRange("Pressure>=0.");
  BADVacPresCmd->SetUnitCategory("Pressure");
  BADVacPresCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  SDpercentageCmd = new G4UIcmdWithADouble("/eli/det/setSDPercentage",this);
  SDpercentageCmd->SetGuidance("Set Percentage of the active isotope from the Sensitive Detector");
  SDpercentageCmd->SetParameterName("Percentage",false);
  SDpercentageCmd->SetRange("Percentage>=0. && Percentage<=100.");
  SDpercentageCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SDmaterialCmd = new G4UIcmdWithAString("/eli/det/setSDMaterial",this);
  SDmaterialCmd->SetGuidance("Set the Gas Type for the neutron counters");
  SDmaterialCmd->SetParameterName("choice",false);
  SDmaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SDpressureCmd = new G4UIcmdWithADoubleAndUnit("/eli/det/setSDPressure",this);
  SDpressureCmd->SetGuidance("Set the Pressure of the gas inside the neutron counters");
  SDpressureCmd->SetParameterName("Pressure",false);
  SDpressureCmd->SetRange("Pressure>=0.");
  SDpressureCmd->SetUnitCategory("Pressure");
  SDpressureCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  ModeratorSizeCmd = new G4UIcmdWithADoubleAndUnit("/eli/det/setModeratorSize",this);
  ModeratorSizeCmd->SetGuidance("Set the Size of the Moderator Section");
  ModeratorSizeCmd->SetParameterName("Dimension",false);
  ModeratorSizeCmd->SetRange("Dimension>=0.");
  ModeratorSizeCmd->SetUnitCategory("Length");
  ModeratorSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NumberOfRingsCmd = new G4UIcmdWithAnInteger("/eli/det/setNumberOfRings",this);
  NumberOfRingsCmd->SetGuidance("Set number of rings");
  NumberOfRingsCmd->SetParameterName("NbRings",false);
  NumberOfRingsCmd->SetRange("NbRings>0 && NbRings<4");
  NumberOfRingsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

/*
  UpdateCmd = new G4UIcmdWithoutParameter("/N03/det/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

eliLaBr_DetectorMessenger::~eliLaBr_DetectorMessenger()
{
/*
  delete NbLayersCmd;
  delete UpdateCmd;
*/
  delete E8_FillingCmd;
  delete TGthicknessCmd;
  delete TGdiameterCmd;
  delete TGplacementCmd;
  delete TGhousingCmd;
  delete TGpositionCmd;
  delete TGmaterialCmd;
  delete ATthicknessCmd;
  delete ATholeCmd;
  delete ATpositionCmd;
  delete ATmaterialCmd;
  delete ATplacementCmd;
  delete PlaceWindowCmd;
  delete PlaceMirrorCmd;
  delete BADVacPresCmd;
  delete SDpercentageCmd;
  delete SDmaterialCmd;
  delete SDpressureCmd;
  delete ModeratorSizeCmd;
  delete NumberOfRingsCmd;
  delete detDir;
  delete eliLaBrDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eliLaBr_DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == E8_FillingCmd )
   { Detector->SetE8_Filling(newValue);}

  if( command == TGthicknessCmd )
   { Detector->SetTarget_thickness(TGthicknessCmd
                                               ->GetNewDoubleValue(newValue));}

  if( command == TGdiameterCmd )
   { Detector->SetTarget_diameter(TGdiameterCmd
                                               ->GetNewDoubleValue(newValue));}

  if( command == TGplacementCmd )
   { Detector->SetTarget_placement(newValue);}

  if( command == TGhousingCmd )
   { Detector->SetTarget_housing(TGhousingCmd
                                               ->GetNewBoolValue(newValue));}

  if( command == TGpositionCmd )
   { Detector->SetTarget_position(TGpositionCmd
                                               ->GetNewDoubleValue(newValue));}

  if( command == TGmaterialCmd )
   { Detector->SetTarget_material(newValue);}

  if( command == ATthicknessCmd )
   { Detector->SetAttenuator_thickness(ATthicknessCmd
                                               ->GetNewDoubleValue(newValue));}

  if( command == ATholeCmd )
   { Detector->SetAttenuator_hole(ATholeCmd->GetNewDoubleValue(newValue));}

  if( command == ATpositionCmd )
   { Detector->SetAttenuator_position(ATpositionCmd
                                               ->GetNewDoubleValue(newValue));}

  if( command == ATmaterialCmd )
   { Detector->SetAttenuator_material(newValue);}

  if( command == ATplacementCmd )
   { Detector->SetAttenuator_placement(ATplacementCmd
                                               ->GetNewBoolValue(newValue));}

  if( command == PlaceWindowCmd )
   { Detector->SetPlace_window(PlaceWindowCmd->GetNewBoolValue(newValue));}

  if( command == PlaceMirrorCmd )
   { Detector->SetPlace_mirror(PlaceMirrorCmd->GetNewBoolValue(newValue));}

  if( command == BADVacPresCmd )
   { Detector->Set_BADVacPres(BADVacPresCmd
                                               ->GetNewDoubleValue(newValue));}

  if( command == SDpercentageCmd )
   { Detector->SetSD_percentage(SDpercentageCmd
                                               ->GetNewDoubleValue(newValue));}

  if( command == SDmaterialCmd )
   { Detector->SetSD_material(newValue);}

  if( command == SDpressureCmd )
   { Detector->SetSD_pressure(SDpressureCmd

                                               ->GetNewDoubleValue(newValue));}
  if( command == ModeratorSizeCmd )
   { Detector->SetModeratorSize(ModeratorSizeCmd
                                               ->GetNewDoubleValue(newValue));}

  if( command == NumberOfRingsCmd )
   { Detector->SetNumberOfRings(NumberOfRingsCmd->GetNewIntValue(newValue));}

/*   
  if( command == UpdateCmd )
   { Detector->UpdateGeometry(); }
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
