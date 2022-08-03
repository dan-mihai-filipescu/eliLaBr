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
/// \file eli/src/PhysicsList.cc
/// \brief Implementation of the eliPhysicsList class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "eliPhysicsList.hh"
#include "eliPhysicsListMessenger.hh"

#include "eliPhysListParticles.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmLivermorePolarizedPhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "eliPhysListEmLivermore_usr.hh"
#include "eliPhysListEmLivermorePolarized_usr.hh"
#include "eliPhysListEmPenelope_usr.hh"
#include "eliPhysListHadron.hh"
#include "eliPhysListEmExtra_usr.hh"
#include "eliPhysListRadDecay_usr.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"

#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"

#include "HadronPhysicsFTFP_BERT.hh"
#include "HadronPhysicsFTF_BIC.hh"
#include "HadronPhysicsLHEP.hh"
#include "HadronPhysicsLHEP_EMV.hh"
#include "G4HadronInelasticQBBC.hh"
#include "HadronPhysicsQGSC_BERT.hh"
#include "HadronPhysicsQGSP.hh"
#include "HadronPhysicsQGSP_BERT.hh"
#include "HadronPhysicsQGSP_BERT_HP.hh"
#include "HadronPhysicsQGSP_BIC.hh"
#include "HadronPhysicsQGSP_BIC_HP.hh"
#include "HadronPhysicsQGSP_FTFP_BERT.hh"
#include "HadronPhysicsQGS_BIC.hh"

#include "G4EmExtraPhysics.hh"
#include "G4EmProcessOptions.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsXS.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronElasticPhysicsLHEP.hh"
#include "G4HadronQElasticPhysics.hh"
#include "G4ChargeExchangePhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4QStoppingPhysics.hh"
#include "G4LHEPStoppingPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4NeutronCrossSectionXS.hh"
#include "G4DecayPhysics.hh"

#include "G4LElastic.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

eliPhysicsList::eliPhysicsList() : G4VModularPhysicsList(),
 fHadPhysicsList(0),fDetectorCuts(0), fTargetCuts(0)
{
  G4LossTableManager::Instance();
  defaultCutValue =1.0*mm;
  fCutForGamma     = defaultCutValue;
  fCutForElectron  = defaultCutValue;
  fCutForPositron  = defaultCutValue;
  fCutForProton    = 0.01*mm;
  fCutForNeutron   = 0.00001*mm;

  //add new units for radioActive decays
  //
  const G4double minute = 60*second;
  const G4double hour   = 60*minute;
  const G4double day    = 24*hour;
  const G4double year   = 365*day;
  new G4UnitDefinition("minute", "min", "Time", minute);
  new G4UnitDefinition("hour",   "h",   "Time", hour);
  new G4UnitDefinition("day",    "d",   "Time", day);
  new G4UnitDefinition("year",   "y",   "Time", year);

  fPMessenger = new eliPhysicsListMessenger(this);

  SetVerboseLevel(1);
  
  //default physics
  fParticleList = new G4DecayPhysics();
  fEliParticleList = new eliPhysListParticles("particles");

  //default physics
  fRaddecayList = new G4RadioactiveDecayPhysics();

  // EM physics
  fEmPhysicsList = new G4EmStandardPhysics();
  //fEmPhysicsList = new G4EmLivermorePolarizedPhysics();
  //SelectPhysicsList("LowEnergy_EM");
  
  //SelectPhysicsList("QGSP_BIC_HP");
  //SelectPhysicsList("QGSP_BERT_HP");
  //SelectPhysicsList("Hadron_usr");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

eliPhysicsList::~eliPhysicsList()
{
  delete fPMessenger;
  delete fParticleList;
  delete fEliParticleList;
  delete fRaddecayList;
  delete fEmPhysicsList;
  if (fHadPhysicsList) delete fHadPhysicsList;
  if (fHadronPhys.size() > 0) {
    for(size_t i=0; i<fHadronPhys.size(); i++) {
      delete fHadronPhys[i];
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eliPhysicsList::ConstructParticle()
{
  //fParticleList->ConstructParticle();
  fEliParticleList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eliPhysicsList::ConstructProcess()
{
  AddTransportation();
  // em
  fEmPhysicsList->ConstructProcess();
  // decays
  fParticleList->ConstructProcess();
  fRaddecayList->ConstructProcess();
  // had
  if (fHadronPhys.size() > 0) {
    for(size_t i=0; i<fHadronPhys.size(); i++) {
      (fHadronPhys[i])->ConstructProcess();
    }
  }
  if (fHadPhysicsList) fHadPhysicsList->ConstructProcess();
  G4cout << "### eliPhysicsList::ConstructProcess is done" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eliPhysicsList::SelectPhysicsList(const G4String& name)
{
  if (verboseLevel>1) {
    G4cout << "eliPhysicsList::SelectPhysicsList: <" << name << ">" << G4endl;
  }

  if (name == "emstandard_opt0") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics();
  } else if (name == "emstandard_opt1") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option1();
  } else if (name == "emstandard_opt2") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option2();
  } else if (name == "emstandard_opt3") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option3();
  } else if (name == "EmLivermore") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLivermorePhysics();
  } else if (name == "EmLivermore_usr") {
    delete fEmPhysicsList;
    fEmPhysicsList = new eliPhysListEmLivermore_usr("Livermore_usr");
  } else if (name == "EmLivermore&SynchRad_usr") {
    delete fEmPhysicsList;
    fEmPhysicsList = new eliPhysListEmLivermore_usr("Livermore&SynchRad_usr");
  } else if (name == "EmLivermorePolarized") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLivermorePolarizedPhysics();
  } else if (name == "EmLivermorePolarized_usr") {
    delete fEmPhysicsList;
    fEmPhysicsList = new eliPhysListEmLivermorePolarized_usr("LivermorePolarized_usr");
  } else if (name == "EmPenelope") {
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmPenelopePhysics();
  } else if (name == "EmPenelope_usr") {
    delete fEmPhysicsList;
    fEmPhysicsList = new eliPhysListEmPenelope_usr("Penelope_usr");
  } else if (name == "FTFP_BERT" && !fHadPhysicsList) {
    SetBuilderList1();
    fHadronPhys.push_back( new HadronPhysicsFTFP_BERT());
  } else if (name == "FTFP_BERT_EMV") {
    SelectPhysicsList("emstandard_opt1");
    SelectPhysicsList("FTFP_BERT");
  } else if (name == "FTFP_BERT_EMX") {
    SelectPhysicsList("emstandard_opt2");
    SelectPhysicsList("FTFP_BERT");
  } else if (name == "FTF_BIC" && !fHadPhysicsList) {
    SetBuilderList0();
    fHadronPhys.push_back( new HadronPhysicsFTF_BIC());
    fHadronPhys.push_back( new G4NeutronCrossSectionXS(verboseLevel));
  } else if (name == "LHEP" && !fHadPhysicsList) {
    SetBuilderList2();
    fHadronPhys.push_back( new HadronPhysicsLHEP());
  } else if (name == "LHEP_EMV" && !fHadPhysicsList) {
    SelectPhysicsList("emstandard_opt1");
    SetBuilderList2(true);
    fHadronPhys.push_back( new HadronPhysicsLHEP_EMV());
  } else if (name == "QBBC" && !fHadPhysicsList) {
    SelectPhysicsList("emstandard_opt2");
    SetBuilderList3();
    fHadronPhys.push_back( new G4HadronInelasticQBBC());
  } else if (name == "QGSC_BERT" && !fHadPhysicsList) {
    SetBuilderList4();
    fHadronPhys.push_back( new HadronPhysicsQGSC_BERT());
  } else if (name == "QGS_BIC" && !fHadPhysicsList) {
    SetBuilderList0();
    fHadronPhys.push_back( new HadronPhysicsQGS_BIC());
    fHadronPhys.push_back( new G4NeutronCrossSectionXS(verboseLevel));
  } else if (name == "QGSP" && !fHadPhysicsList) {
    SetBuilderList1();
    fHadronPhys.push_back( new HadronPhysicsQGSP());
  } else if (name == "QGSP_BERT" && !fHadPhysicsList) {
    SetBuilderList1();
    fHadronPhys.push_back( new HadronPhysicsQGSP_BERT("std-hadron"));
  } else if (name == "QGSP_BERT_EMV") {
    SelectPhysicsList("emstandard_opt1");
    SelectPhysicsList("QGSP_BERT");
  } else if (name == "QGSP_BERT_EMX") {
    SelectPhysicsList("emstandard_opt2");
    SelectPhysicsList("QGSP_BERT");
  } else if (name == "QGSP_BERT_HP" && !fHadPhysicsList) {
    SetBuilderList1(true);
    fHadronPhys.push_back( new HadronPhysicsQGSP_BERT_HP("std-hadron"));
  } else if (name == "QGSP_FTFP_BERT" && !fHadPhysicsList) {
    SetBuilderList1();
    fHadronPhys.push_back( new HadronPhysicsQGSP_FTFP_BERT());
  } else if (name == "QGSP_BIC" && !fHadPhysicsList) {
    SetBuilderList0();
    fHadronPhys.push_back( new HadronPhysicsQGSP_BIC("std-hadron"));
  } else if (name == "QGSP_BIC_EMY" && !fHadPhysicsList) {
    SelectPhysicsList("emstandard_opt3");
    SetBuilderList0();
    fHadronPhys.push_back( new HadronPhysicsQGSP_BIC("std-hadron"));
  } else if (name == "QGSP_BIC_HP" && !fHadPhysicsList) {
    SetBuilderList0(true);
    fHadronPhys.push_back( new HadronPhysicsQGSP_BIC_HP("std-hadron"));
  } else if (name == "Hadron_usr" && !fHadPhysicsList) {
    SetBuilderList();
    fHadronPhys.push_back( new eliPhysListHadron("hadron"));
  } else if (name == "ExtraEmPhys") {
    fHadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
  } else if (name == "ExtraEmPhys_usr") {
    fHadronPhys.push_back( new eliPhysListEmExtra_usr("EmExtra_usr"));
  } else if (name == "RadDecay_usr") {
    delete fRaddecayList;
    fRaddecayList = new eliPhysListRadDecay_usr ("RadDecay_usr");
  } else {
      G4cout << "eliPhysicsList WARNING wrong or unkonwn <"
             << name << "> Physics " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eliPhysicsList::SetBuilderList()
{
  //fNhadcomp = 2;
  //fHadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
  //fHadronPhys.push_back( new G4HadronElasticPhysics("elastic",verboseLevel,flagHP));
  fHadronPhys.push_back( new G4StoppingPhysics("stopping",verboseLevel));
  //fHadronPhys.push_back( new G4IonBinaryCascadePhysics("ionBIC"));
  fHadronPhys.push_back( new G4NeutronTrackingCut("Neutron tracking cut",verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void eliPhysicsList::SetBuilderList0(G4bool flagHP)
{
  fHadronPhys.push_back( new G4EmExtraPhysics(verboseLevel));
  if(flagHP) {
    fHadronPhys.push_back( new G4HadronElasticPhysicsHP(verboseLevel) );
  } else {
    fHadronPhys.push_back( new G4HadronElasticPhysics(verboseLevel) );
  }
  fHadronPhys.push_back( new G4QStoppingPhysics(verboseLevel));
  fHadronPhys.push_back( new G4IonBinaryCascadePhysics(verboseLevel));
  fHadronPhys.push_back( new G4NeutronTrackingCut(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void eliPhysicsList::SetBuilderList1(G4bool flagHP)
{
  fHadronPhys.push_back( new G4EmExtraPhysics(verboseLevel));
  if(flagHP) {
    fHadronPhys.push_back( new G4HadronElasticPhysicsHP(verboseLevel) );
  } else {
    fHadronPhys.push_back( new G4HadronElasticPhysics(verboseLevel) );
  }
  fHadronPhys.push_back( new G4QStoppingPhysics(verboseLevel));
  fHadronPhys.push_back( new G4IonPhysics(verboseLevel));
  fHadronPhys.push_back( new G4NeutronTrackingCut(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void eliPhysicsList::SetBuilderList2(G4bool addStopping)
{
  fHadronPhys.push_back( new G4EmExtraPhysics(verboseLevel));
  fHadronPhys.push_back( new G4HadronElasticPhysicsLHEP(verboseLevel));
  if(addStopping) { fHadronPhys.push_back( new G4QStoppingPhysics(verboseLevel)); }
  fHadronPhys.push_back( new G4IonPhysics(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void eliPhysicsList::SetBuilderList3()
{
  fHadronPhys.push_back( new G4EmExtraPhysics(verboseLevel));
  RegisterPhysics( new G4HadronElasticPhysicsXS(verboseLevel) );
  fHadronPhys.push_back( new G4QStoppingPhysics(verboseLevel));
  fHadronPhys.push_back( new G4IonBinaryCascadePhysics(verboseLevel));
  fHadronPhys.push_back( new G4NeutronTrackingCut(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void eliPhysicsList::SetBuilderList4()
{
  fHadronPhys.push_back( new G4EmExtraPhysics(verboseLevel));
  fHadronPhys.push_back( new G4HadronQElasticPhysics(verboseLevel));
  fHadronPhys.push_back( new G4QStoppingPhysics(verboseLevel));
  fHadronPhys.push_back( new G4IonPhysics(verboseLevel));
  fHadronPhys.push_back( new G4NeutronTrackingCut(verboseLevel));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void eliPhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(fCutForGamma, "gamma");
  SetCutValue(fCutForElectron, "e-");
  SetCutValue(fCutForPositron, "e+");
  SetCutValue(fCutForProton, "proton");
  SetCutValue(fCutForNeutron, "neutron");
  G4cout << "world cuts are set" << G4endl;

  if( !fTargetCuts ) SetTargetCut(defaultCutValue);
  G4Region* region = (G4RegionStore::GetInstance())->GetRegion("Target");
  region->SetProductionCuts(fTargetCuts);
  G4cout << "Target cuts are set" << G4endl;

  if( !fDetectorCuts ) SetDetectorCut(defaultCutValue);
  region = (G4RegionStore::GetInstance())->GetRegion("Detector");
  region->SetProductionCuts(fDetectorCuts);
  G4cout << "Detector cuts are set" << G4endl;

  if (verboseLevel>0) {/*DumpList();*/ DumpCutValuesTable();}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eliPhysicsList::SetCutForGamma(G4double cut)
{
  fCutForGamma = cut;
  SetParticleCuts(fCutForGamma, G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eliPhysicsList::SetCutForElectron(G4double cut)
{
  fCutForElectron = cut;
  SetParticleCuts(fCutForElectron, G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eliPhysicsList::SetCutForPositron(G4double cut)
{
  fCutForPositron = cut;
  SetParticleCuts(fCutForPositron, G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eliPhysicsList::SetCutForProton(G4double cut)
{
  fCutForProton = cut;
  SetParticleCuts(fCutForProton, G4Proton::Proton());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eliPhysicsList::SetCutForNeutron(G4double cut)
{
  fCutForNeutron = cut;
  SetParticleCuts(fCutForNeutron, G4Neutron::Neutron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eliPhysicsList::SetTargetCut(G4double cut)
{
  if( !fTargetCuts ) fTargetCuts = new G4ProductionCuts();

  fTargetCuts->SetProductionCut(cut, idxG4GammaCut);
  fTargetCuts->SetProductionCut(cut, idxG4ElectronCut);
  fTargetCuts->SetProductionCut(cut, idxG4PositronCut);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eliPhysicsList::SetDetectorCut(G4double cut)
{
  if( !fDetectorCuts ) fDetectorCuts = new G4ProductionCuts();

  fDetectorCuts->SetProductionCut(cut, idxG4GammaCut);
  fDetectorCuts->SetProductionCut(cut, idxG4ElectronCut);
  fDetectorCuts->SetProductionCut(cut, idxG4PositronCut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// TAKEN FROM EXAMPLES/NOVICE/N02/........ TRYING TO APPLY CUTS ON NEUTRONS

#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

void eliPhysicsList::AddStepMax()
{

  // Step limitation seen as a process
  // G4StepLimiter* stepLimiter = new G4StepLimiter();
  G4UserSpecialCuts* userCuts = new G4UserSpecialCuts();

  theParticleIterator->reset();
  while ((*theParticleIterator)()){
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      pmanager ->AddDiscreteProcess(userCuts);
      // userCUTS (outside "IF")   10.000 events in 1 min 10 sec - SAME as inside "if"
/*
      if ( particle->GetPDGCharge() == 0.0 )
        {
51:25	  // pmanager ->AddDiscreteProcess(stepLimiter);
      pmanager ->AddDiscreteProcess(userCuts);

// using this "IF"........
// with STEP LIMITER:::   10.000 events in 2 min 30 sec
//	userCUTS:::::::   10.000 events in 1 min 10 sec

        }
*/
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void eliPhysicsList::List()
{
  G4cout << "### PhysicsLists available: " << G4endl;
  G4cout << "===> Electromagnetic: emstandard_opt0 emstandard_opt1 emstandard_opt2 emstandard_opt3" << G4endl;
  G4cout << "                      EmPenelope EmLivermore EmLivermorePolarized"                     << G4endl;
  G4cout << "                      EmPenelope_usr EmLivermore_usr EmLivermorePolarized_usr"         << G4endl;
  G4cout << "===> Hadron only:     FTFP_BERT FTF_BIC LHEP QGS_BIC QGSP QGSP_FTFP_BERT"              << G4endl;
  G4cout << "                      QGSC_BERT QGSP_BERT QGSP_BERT_HP QGSP_BIC QGSP_BIC_HP Hadron_usr"<< G4endl;
  G4cout << "===> Em+Hadron:       FTFP_BERT_EMV FTFP_BERT_EMX LHEP_EMV QBBC"                       << G4endl;
  G4cout << "                      QGSP_BERT_EMV QGSP_BERT_EMX QGSP_BIC_EMY"                        << G4endl;
  G4cout << "===> ExtraEM_Nucl:    ExtraEmPhys  ExtraEmPhys_usr RadDecay_usr"                       << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

