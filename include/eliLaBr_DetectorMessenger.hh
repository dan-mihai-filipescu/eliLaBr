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

#ifndef eliLaBr_DetectorMessenger_h
#define eliLaBr_DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class eliLaBr_DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class eliLaBr_DetectorMessenger: public G4UImessenger
{
  public:
    eliLaBr_DetectorMessenger(eliLaBr_DetectorConstruction* );
   ~eliLaBr_DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    eliLaBr_DetectorConstruction* Detector;
    
    G4UIdirectory*             eliLaBrDir;
    G4UIdirectory*             detDir;
    G4UIcmdWithAString*        E8_FillingCmd;
    G4UIcmdWithADoubleAndUnit* TGthicknessCmd;
    G4UIcmdWithADoubleAndUnit* TGdiameterCmd;
    G4UIcmdWithAString*        TGplacementCmd;
    G4UIcmdWithABool*          TGhousingCmd;
    G4UIcmdWithADoubleAndUnit* TGpositionCmd;
    G4UIcmdWithAString*        TGmaterialCmd;
    G4UIcmdWithADoubleAndUnit* ATthicknessCmd;
    G4UIcmdWithADoubleAndUnit* ATholeCmd;
    G4UIcmdWithADoubleAndUnit* ATpositionCmd;
    G4UIcmdWithAString*        ATmaterialCmd;
    G4UIcmdWithABool*          ATplacementCmd;
    G4UIcmdWithABool*          PlaceWindowCmd;
    G4UIcmdWithABool*          PlaceMirrorCmd;
    G4UIcmdWithADoubleAndUnit* BADVacPresCmd;
    G4UIcmdWithADouble*        SDpercentageCmd;
    G4UIcmdWithAString*        SDmaterialCmd;
    G4UIcmdWithADoubleAndUnit* SDpressureCmd;
    G4UIcmdWithADoubleAndUnit* ModeratorSizeCmd;
    G4UIcmdWithAnInteger*      NumberOfRingsCmd;
    
/*    
    G4UIcmdWithoutParameter*   UpdateCmd;
*/
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
