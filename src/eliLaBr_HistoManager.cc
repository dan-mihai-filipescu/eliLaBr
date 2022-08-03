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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "eliLaBr_HistoManager.hh"
#include "eliLaBr_Analysis.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

eliLaBr_HistoManager::eliLaBr_HistoManager()
  : fFileName("eliLaBr_HistFile_hist")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

eliLaBr_HistoManager::~eliLaBr_HistoManager()
{
    //We do not destruct here the G4AnalysisManager anymore
    //this done inside eliLaBr_Analysis class
    //delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eliLaBr_HistoManager::Book()
{
  // Book histograms, ntuple
  //

  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in eliLaBr_HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType()  << " analysis manager" << G4endl;

  // Create directories
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  // Open an output file
  //
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);     //enable inactivation of histograms
  analysisManager->SetFirstHistoId(1);
  
  // Define histograms start values
  //const G4int kMaxHisto = 9;
  G4int ih;

  // Create 1D histograms (activated)
  ih =analysisManager->CreateH1("H1Sngl","Single",     2000, 0., 20*MeV);
  analysisManager->SetActivation(G4VAnalysisManager::kH1, ih, true);
  ih =analysisManager->CreateH1("H1SnglN","SingleN",     2000, 0., 20*MeV);
  analysisManager->SetActivation(G4VAnalysisManager::kH1, ih, true);
  ih =analysisManager->CreateH1("H1Inc","Incident",   2000, 0., 20*MeV);
  analysisManager->SetActivation(G4VAnalysisManager::kH1, ih, true);
  ih =analysisManager->CreateH1("H1IncN","IncidentN",   2000, 0., 20*MeV);
  analysisManager->SetActivation(G4VAnalysisManager::kH1, ih, true);
  ih = analysisManager->CreateH1("Angle","PolAngle",   3600, 0., 360.);
  analysisManager->SetActivation(G4VAnalysisManager::kH1, ih, false);
  ih = analysisManager->CreateH1("AngleN","PolAngleN",   3600, 0., 360.);
  analysisManager->SetActivation(G4VAnalysisManager::kH1, ih, false);
  ih =analysisManager->CreateH1("H1Zint","Zint",     2020, -10.1*m, 10.1*m);
  analysisManager->SetActivation(G4VAnalysisManager::kH1, ih, true);
  ih =analysisManager->CreateH1("H1ZintN","ZintN",   2020, -10.1*m, 10.1*m);
  analysisManager->SetActivation(G4VAnalysisManager::kH1, ih, true);
  ih = analysisManager->CreateH1("AngleIni","PolAngleIni",   3600, 0., 360.);
  analysisManager->SetActivation(G4VAnalysisManager::kH1, ih, false);

  // Create 2D histograms (deactivated)
  //ih = analysisManager->CreateH1("Stokes","Stokes",   3600, -1., +1.);
  //analysisManager->SetActivation(G4VAnalysisManager::kH1, ih, false);
  ih = analysisManager->CreateH2("GammaBeam-profile","Gamma beam profile", 300, -2., 2., 300, -2., 2.,"mm","mm","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);
  ih = analysisManager->CreateH2("En-radius","En-radius", 300, 0., 20., 300, 0., 2.,"MeV","mm","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);
  ih = analysisManager->CreateH2("En-X","En-X", 300, -2., 2., 300, 0., 20., "mm", "MeV","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);
  ih = analysisManager->CreateH2("En-Y","En-Y", 300, -2., 2., 300, 0., 20., "mm", "MeV","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);

  ih = analysisManager->CreateH2("ele_pos","Electron space distribution", 300, -0.2, 0.2, 300, -0.2, 0.2,"mm","mm","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);
  ih = analysisManager->CreateH2("laser_pos","Interaction point XY distribution", 300, -0.2, 0.2, 300, -0.2, 0.2,"mm","mm","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);
  ih = analysisManager->CreateH2("horiz_xxp","Interaction HORIZONTAL spacephase distribution", 300, -0.2, 0.2, 300, -5.e-5, 5.e-5,"mm","radian","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);
  ih = analysisManager->CreateH2("vert_yyp","Interaction VERTICAL spacephase distribution", 300, -0.2, 0.2, 300, -5.e-5, 5.e-5,"mm","radian","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);
  ih = analysisManager->CreateH2("laser_imp","Laser Transversal Impulse distribution", 300, -0.3, 0.3, 300, -0.3, 0.3,"none","none","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);

  ih = analysisManager->CreateH2("stokes1","Stokes1", 1000, 0.0, 20., 1000, -1.01, 1.01, "MeV","none","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);
  ih = analysisManager->CreateH2("stokes2","Stokes2", 1000, 0.0, 20., 1000, -1.01, 1.01, "MeV","none","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);
  ih = analysisManager->CreateH2("Plin","Plin Sqrt(S1^2 + S2^2)", 1000, 0.0, 20., 1000, -1.01, 1.01, "MeV","none","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);
  ih = analysisManager->CreateH2("stokes3","Stokes3", 1000, 0.0, 20., 1000, -1.01, 1.01, "MeV","none","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);
  ih = analysisManager->CreateH2("PolAngle","PolAngle", 1000, 0.0, 20., 300, 0, 180., "MeV","deg","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);

  ih = analysisManager->CreateH2("source_aniso","source_aniso", 300, 0.0, 180., 300, 0.0, 360., "deg", "deg", "none", "none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);

  ih = analysisManager->CreateH2("target-ET","target-ET", 2000, 0.0, 20., 5000, 0.0, 500.);
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);

  ih = analysisManager->CreateH2("GammaBeam-EDep_profile","Gamma beam EneDep profile", 300, -2., 2., 300, -2., 2.,"mm","mm","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);
  ih = analysisManager->CreateH2("GammaBeam-EDepNorm_profile","Gamma beam EneDepNorm profile", 300, -2., 2., 300, -2., 2.,"mm","mm","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);

  ih = analysisManager->CreateH2("Perp_CS-Parallel_CS","Perp. CS / Parallel CS  polariz. scatt.", 300, 0., 180., 300, 0., 360.,"deg","deg","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);

  ih = analysisManager->CreateH2("stokes_surf1","StokesSurf1", 300, -2., 2., 300, -2., 2.,"mm","mm","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);
  ih = analysisManager->CreateH2("stokes_surf2","StokesSurf2", 300, -2., 2., 300, -2., 2.,"mm","mm","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);
  ih = analysisManager->CreateH2("stokes_surf3","StokesSurf3", 300, -2., 2., 300, -2., 2.,"mm","mm","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);

  ih = analysisManager->CreateH2("GammaBeam-AngleProfile","Gamma beam angle profile", 300, 0., 180., 300, 0., 360.,"deg","deg","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);
  ih = analysisManager->CreateH2("Polariz_Degree","Polarization Degree", 300, 0., 180., 300, 0., 360.,"deg","deg","none","none");
  analysisManager->SetActivation(G4VAnalysisManager::kH2, ih, false);

/*  const G4String id[] = {"0","1","2","3","4","5","6","7","8"};
  const G4String title[] = 
          { "dummy",                                    //0
            "energy spectrum (%): e+ e-",               //1
            "energy spectrum (%): nu_e anti_nu_e",      //2
            "energy spectrum (%): gamma",               //3                  
            "energy spectrum (%): alpha",               //4
            "energy spectrum (%): ions",                //5
            "total kinetic energy (Q)",                 //6                        
            "momentum balance",                         //7
            "total time of life of decay chain"         //8
          };  

  // Default values (to be reset via /analysis/h1/set command)               
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;


  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto; k++) {
    G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
    analysisManager->SetActivation(G4VAnalysisManager::kH1, ih, false);
  }*/

  /*  // Creating ntuple
    //
    analysisManager->CreateNtuple("B4", "Edep and TrackL");
    analysisManager->CreateNtupleDColumn("Eabs");
    analysisManager->CreateNtupleDColumn("Egap");
    analysisManager->CreateNtupleDColumn("Labs");
    analysisManager->CreateNtupleDColumn("Lgap");
    analysisManager->FinishNtuple();
  */

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
