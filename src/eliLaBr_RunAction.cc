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
// $Id$
//
/// \file eliLaBr_RunAction.cc
/// \brief Implementation of the eliLaBr_RunAction class

#include "eliLaBr_RunAction.hh"
#include "eliLaBr_Analysis.hh"
#include "eliLaBr_HistoManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

eliLaBr_RunAction::eliLaBr_RunAction()
 : G4UserRunAction()
{ 
    fHistoManager = new eliLaBr_HistoManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

eliLaBr_RunAction::~eliLaBr_RunAction()
{
    delete fHistoManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eliLaBr_RunAction::BeginOfRunAction(const G4Run* run)
{ 
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;

  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  //histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->OpenFile();
  }

//  G4String fileName = "eliLaBr_OutputFile_hist";
//  analysisManager->OpenFile(fileName);



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eliLaBr_RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int nofEvents = aRun->GetNumberOfEvent();
  if ( nofEvents == 0 ) return;

  // print histogram statistics
  //
  //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

 /* if ( analysisManager->GetH1(1) ) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
       << " EAbs : mean = " << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy")
               << " rms = " << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") 
               << G4endl;
    G4cout                
       << " EGap : mean = " << G4BestUnit(analysisManager->GetH1(2)->mean(), "Energy") 
               << " rms = " << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Energy") 
               << G4endl;
  }*/

  // save histograms 
  //
  //analysisManager->Write();
  //analysisManager->CloseFile();
  
  // complete cleanup
  //
  //delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
