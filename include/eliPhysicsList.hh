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
/// \file eli/include/eliPhysicsList.hh
/// \brief Definition of the eliPhysicsList class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef eliPhysicsList_h
#define eliPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include <vector>

class G4VPhysicsConstructor;
class eliPhysicsListMessenger;
class G4ProductionCuts;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class eliPhysicsList: public G4VModularPhysicsList
{
public:
  eliPhysicsList();
  virtual ~eliPhysicsList();

  virtual void ConstructParticle();

  virtual void SetCuts();
  void SetCutForGamma(G4double);
  void SetCutForElectron(G4double);
  void SetCutForPositron(G4double);
  void SetCutForProton(G4double);
  void SetCutForNeutron(G4double);

  void SelectPhysicsList(const G4String& name);
  virtual void ConstructProcess();
  void AddStepMax();
  void List();

  void SetTargetCut(G4double val);
  void SetDetectorCut(G4double val);

private:

  void SetBuilderList ();
  void SetBuilderList0(G4bool flagHP = false);
  void SetBuilderList1(G4bool flagHP = false);
  void SetBuilderList2(G4bool addStopping = false);
  void SetBuilderList3();
  void SetBuilderList4();

  // hide assignment operator
  eliPhysicsList & operator=(const eliPhysicsList &right);
  eliPhysicsList(const eliPhysicsList&);

  G4double fCutForGamma;
  G4double fCutForElectron;
  G4double fCutForPositron;
  G4double fCutForNeutron;
  G4double fCutForProton;

  G4VPhysicsConstructor*  fEmPhysicsList;
  G4VPhysicsConstructor*  fRaddecayList;
  G4VPhysicsConstructor*  fParticleList;
  G4VPhysicsConstructor*  fEliParticleList;
  G4VPhysicsConstructor*  fHadPhysicsList;
  
  std::vector<G4VPhysicsConstructor*>  fHadronPhys;
  //G4int fNhadcomp;

  eliPhysicsListMessenger* fPMessenger;
  G4ProductionCuts* fDetectorCuts;
  G4ProductionCuts* fTargetCuts;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

