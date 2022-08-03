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
/// \file eli/src/eliPhysListHadron.cc
/// \brief Implementation of the eliPhysListHadron class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "eliPhysListHadron.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4Evaporation.hh"
#include "G4FermiBreakUp.hh"
#include "G4StatMF.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4PreCompoundModel.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

#include "G4LElastic.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"

//#include "G4CascadeInterface.hh"

#include "G4BinaryCascade.hh"
#include "G4BinaryLightIonReaction.hh"

//HPNeutron

#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPThermalScattering.hh"
#include "G4NeutronHPThermalScatteringData.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"
#include "G4NeutronHPFission.hh"
#include "G4NeutronHPFissionData.hh"

#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"

//#include "G4PhotoNuclearProcess.hh"
//#include "G4CascadeInterface.hh"
//#include "G4GammaNuclearReaction.hh"
//#include "G4ElectronNuclearProcess.hh"
//#include "G4PositronNuclearProcess.hh"
//#include "G4ElectroVDNuclearModel.hh"

//c-s
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"

#include "G4WilsonAblationModel.hh"
#include "G4WilsonAbrasionModel.hh"

// RadioactiveDecay
#include "G4RadioactiveDecay.hh"
#include "G4GenericIon.hh"

#include "G4SystemOfUnits.hh"


eliPhysListHadron::eliPhysListHadron(const G4String& name)
  :  G4VPhysicsConstructor(name),
    fTheNeutronElasticProcess(0),
   fTheFissionProcess(0),
   fTheCaptureProcess(0),fTheDeuteronInelasticProcess(0),
   fTheTritonInelasticProcess(0), fTheAlphaInelasticProcess(0),
   fTheIonInelasticProcess(0)
{}

eliPhysListHadron::~eliPhysListHadron()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void eliPhysListHadron::ConstructProcess()
{
   
  G4ProcessManager * pManager = 0;

//
/*  pManager = G4Gamma::Gamma()->GetProcessManager();
  G4PhotoNuclearProcess* process = new G4PhotoNuclearProcess();
//  process->BiasCrossSectionByFactor(1000);
  G4CascadeInterface* bertini = new G4CascadeInterface();
  bertini->SetMaxEnergy(10*GeV);
  process->RegisterMe(bertini);
//  G4GammaNuclearReaction* gammaNuclear = new G4GammaNuclearReaction();
//  process->RegisterMe(gammaNuclear);
  pManager->AddDiscreteProcess(process);

  pManager = G4Electron::Electron()->GetProcessManager();
  G4ElectronNuclearProcess * theElectronNuclearProcess = new G4ElectronNuclearProcess;
  G4ElectroVDNuclearModel * theElectroReaction = new G4ElectroVDNuclearModel;
  theElectronNuclearProcess->RegisterMe(theElectroReaction);
  pManager->AddDiscreteProcess(theElectronNuclearProcess);

  pManager = G4Positron::Positron()->GetProcessManager();
  G4PositronNuclearProcess * thePositronNuclearProcess = new G4PositronNuclearProcess;
  thePositronNuclearProcess->RegisterMe(theElectroReaction);
  pManager->AddDiscreteProcess(thePositronNuclearProcess);*/
  
  // this will be the model class for high energies
  G4TheoFSGenerator * theTheoModel = new G4TheoFSGenerator;
  // all models for treatment of thermal nucleus 
  G4Evaporation * theEvaporation = new G4Evaporation;
  G4FermiBreakUp * theFermiBreakUp = new G4FermiBreakUp;
  G4StatMF * theMF = new G4StatMF;
  // Evaporation logic
  G4ExcitationHandler * theHandler = new G4ExcitationHandler;
  theHandler->SetEvaporation(theEvaporation);
  theHandler->SetFermiModel(theFermiBreakUp);
  theHandler->SetMultiFragmentation(theMF);
  theHandler->SetMaxAandZForFermiBreakUp(12, 6);
  theHandler->SetMinEForMultiFrag(5*MeV);  

  // Pre equilibrium stage
  G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel(theHandler);
  // a no-cascade generator-precompound interaface
  G4GeneratorPrecompoundInterface * theCascade = new G4GeneratorPrecompoundInterface;
  theCascade->SetDeExcitation(thePreEquilib);  
        
  // here come the high energy parts
  // the string model; still not quite according to design
  // - Explicite use of the forseen interfaces
  G4VPartonStringModel * theStringModel;
  theStringModel = new G4QGSModel<G4QGSParticipants>;
  theTheoModel->SetTransport(theCascade);
  theTheoModel->SetHighEnergyGenerator(theStringModel);
  theTheoModel->SetMinEnergy(10*GeV);  // 15 GeV may be the right limit
  theTheoModel->SetMaxEnergy(100*TeV);
  
  G4VLongitudinalStringDecay * theFragmentation = new G4QGSMFragmentation;
  G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay(theFragmentation);
  theStringModel->SetFragmentationModel(theStringDecay);
  
  // Elastic Process
  G4LElastic* theElasticModel = new G4LElastic;
  theElasticModel->SetMaxEnergy(100000.*GeV);
  theElasticModel->SetMinEnergy(0.*MeV);
  fTheElasticProcess.RegisterMe(theElasticModel);

  // ---------------------------------------------------------------------------
  // Hadron elastic process
  // for all particles except neutrons

  theParticleIterator->reset();
  while( (*theParticleIterator)() ) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String particleName = particle->GetParticleName();
    if (particleName != "neutron")
        if ((particleName == "proton") || (particleName == "deuteron") || (particleName == "triton") || (particleName == "alpha")) {
      pManager = particle->GetProcessManager();
      if (particle->GetPDGMass() > 110.*MeV && fTheElasticProcess.IsApplicable(*particle)  && !particle->IsShortLived()) {
        pManager->AddDiscreteProcess(&fTheElasticProcess);
        //
                G4cout << "################################### Elastic model are registered for "
               << particle->GetParticleName()
               << G4endl;
      }
    }
  }
  //

  // ============================================= PROTON ==========================================================
  pManager = G4Proton::Proton()->GetProcessManager();
    // add inelastic process
  // Binary Cascade
  G4BinaryCascade * theBC = new G4BinaryCascade;
  theBC->SetMaxEnergy(10.*GeV);
  theBC->SetMinEnergy(0.*eV);
  fTheProtonInelastic.RegisterMe(theBC);
  // Higher energy
  fTheProtonInelastic.RegisterMe(theTheoModel);
  // now the cross-sections.
  G4ProtonInelasticCrossSection * theProtonData = new G4ProtonInelasticCrossSection;
  fTheProtonInelastic.AddDataSet(theProtonData);
  pManager->AddDiscreteProcess(&fTheProtonInelastic);
  // ================================================================================================================
  //
  // ============================================= NEUTRON ==========================================================
  pManager = G4Neutron::Neutron()->GetProcessManager();
  // add process
  // *** elastic scattering ***
  fTheNeutronElasticProcess = new G4HadronElasticProcess();

  G4LElastic* theElasticModel1 = new G4LElastic;
  theElasticModel1->SetMaxEnergy(100000.*GeV);
  theElasticModel1->SetMinEnergy(20.*MeV);
  fTheNeutronElasticProcess->RegisterMe(theElasticModel1);

  G4NeutronHPElastic * theElasticNeutron = new G4NeutronHPElastic;
  theElasticNeutron->SetMaxEnergy(20.*MeV);
  theElasticNeutron->SetMinEnergy ( 4.0*eV );
  fTheNeutronElasticProcess->RegisterMe(theElasticNeutron);
  G4NeutronHPElasticData * theNeutronData = new G4NeutronHPElasticData;
  fTheNeutronElasticProcess->AddDataSet(theNeutronData);

  G4NeutronHPThermalScattering* theNeutronThermalElasticModel = new G4NeutronHPThermalScattering();
  theNeutronThermalElasticModel->SetMaxEnergy ( 4.0*eV );
  theNeutronThermalElasticModel->SetMinEnergy ( 0.0*eV );
  fTheNeutronElasticProcess->RegisterMe(theNeutronThermalElasticModel);
  G4NeutronHPThermalScatteringData* theHPThermalScatteringData = new G4NeutronHPThermalScatteringData();
  fTheNeutronElasticProcess->AddDataSet( theHPThermalScatteringData );

  pManager->AddDiscreteProcess(fTheNeutronElasticProcess);

  // *** inelastic ***
  G4NeutronHPInelastic * theHPNeutronInelasticModel = new G4NeutronHPInelastic;
  theHPNeutronInelasticModel->SetMaxEnergy(20.*MeV);
  theHPNeutronInelasticModel->SetMinEnergy(0.*eV);
  fTheNeutronInelastic.RegisterMe(theHPNeutronInelasticModel);
  G4NeutronHPInelasticData * theNeutronData1 = new G4NeutronHPInelasticData;
  fTheNeutronInelastic.AddDataSet(theNeutronData1);
  // binary
  G4BinaryCascade * neutronBC = new G4BinaryCascade;
  neutronBC->SetMaxEnergy(10.*GeV);
  neutronBC->SetMinEnergy(20.*MeV);
  fTheNeutronInelastic.RegisterMe(neutronBC);
  G4NeutronInelasticCrossSection * theNeutronData2 = new G4NeutronInelasticCrossSection;
  fTheNeutronInelastic.AddDataSet(theNeutronData2);
  // higher energy
  fTheNeutronInelastic.RegisterMe(theTheoModel);  
  // now the cross-sections.
  pManager->AddDiscreteProcess(&fTheNeutronInelastic);

  // *** fission model ***
  fTheFissionProcess = new G4HadronFissionProcess;
  G4LFission* theFissionModel = new G4LFission;
  theFissionModel->SetMaxEnergy(100000.*GeV);
  theFissionModel->SetMinEnergy(20.*MeV);
  fTheFissionProcess->RegisterMe(theFissionModel);
  // fission HP
  G4NeutronHPFission* theHPFissionModel = new G4NeutronHPFission;
  theHPFissionModel->SetMaxEnergy(20.*MeV);
  theHPFissionModel->SetMinEnergy(0.*eV);
  fTheFissionProcess->RegisterMe(theHPFissionModel);
  G4NeutronHPFissionData* theNeutronDataf = new G4NeutronHPFissionData;
  fTheFissionProcess->AddDataSet(theNeutronDataf);

  pManager->AddDiscreteProcess(fTheFissionProcess);

  // *** capture ***
  fTheCaptureProcess = new G4HadronCaptureProcess;
  G4LCapture* theCaptureModel = new G4LCapture;
  theCaptureModel->SetMaxEnergy(100000.*GeV);
  theCaptureModel->SetMinEnergy(20.*MeV);
  fTheCaptureProcess->RegisterMe(theCaptureModel);
  //capture HP
  G4NeutronHPCapture * theHPNeutronCaptureModel = new G4NeutronHPCapture;
  theHPNeutronCaptureModel->SetMaxEnergy(20.*MeV);
  theHPNeutronCaptureModel->SetMinEnergy(0.*eV);
  fTheCaptureProcess->RegisterMe(theHPNeutronCaptureModel);
  G4NeutronHPCaptureData * theNeutronData3 = new G4NeutronHPCaptureData;
  fTheCaptureProcess->AddDataSet(theNeutronData3);

  pManager->AddDiscreteProcess(fTheCaptureProcess);
  // ================================================================================================================

  // now light ions
  // =========================================== Light Ion BINARY CASCADE ===========================================
  G4BinaryLightIonReaction * theIonBC= new G4BinaryLightIonReaction;
  theIonBC->SetMaxEnergy(10*GeV);
  theIonBC->SetMinEnergy(20*MeV);
  G4TripathiCrossSection * TripathiCrossSection= new G4TripathiCrossSection;
  G4IonsShenCrossSection * aShen = new G4IonsShenCrossSection;
    
  // ============================================= DEUTERON =========================================================
  pManager = G4Deuteron::Deuteron()->GetProcessManager();
  fTheDeuteronInelasticProcess =  new G4DeuteronInelasticProcess("inelastic");
  // higher energy
  fTheDeuteronInelasticProcess->RegisterMe(theTheoModel);
  fTheDeuteronInelasticProcess->RegisterMe(theIonBC);
  fTheDeuteronInelasticProcess->AddDataSet(TripathiCrossSection);
  fTheDeuteronInelasticProcess->AddDataSet(aShen);

  G4LEDeuteronInelastic* theDeuteronInelasticModel = new G4LEDeuteronInelastic;
  theDeuteronInelasticModel->SetMaxEnergy(20*MeV);
  theDeuteronInelasticModel->SetMinEnergy(0*MeV);
  fTheDeuteronInelasticProcess->RegisterMe(theDeuteronInelasticModel);

  pManager->AddDiscreteProcess(fTheDeuteronInelasticProcess);
  // ================================================================================================================

  // ============================================== TRITON ==========================================================
  pManager = G4Triton::Triton()->GetProcessManager();
  fTheTritonInelasticProcess =  new G4TritonInelasticProcess("inelastic");
  // higher energy
  fTheTritonInelasticProcess->RegisterMe(theTheoModel);
  fTheTritonInelasticProcess->RegisterMe(theIonBC);
  fTheTritonInelasticProcess->AddDataSet(TripathiCrossSection);
  fTheTritonInelasticProcess->AddDataSet(aShen);

  G4LETritonInelastic* theTritonInelasticModel = new G4LETritonInelastic;
  theTritonInelasticModel->SetMaxEnergy(20*MeV);
  theTritonInelasticModel->SetMinEnergy(0*MeV);
  fTheTritonInelasticProcess->RegisterMe(theTritonInelasticModel);

  pManager->AddDiscreteProcess(fTheTritonInelasticProcess);
  // ================================================================================================================

  // ============================================== ALPHA ===========================================================
  pManager = G4Alpha::Alpha()->GetProcessManager();
  fTheAlphaInelasticProcess = new G4AlphaInelasticProcess("inelastic");
  // higher energy
  fTheAlphaInelasticProcess->RegisterMe(theTheoModel);
  fTheAlphaInelasticProcess->RegisterMe(theIonBC);
  fTheAlphaInelasticProcess->AddDataSet(TripathiCrossSection);
  fTheAlphaInelasticProcess->AddDataSet(aShen);

  G4LEAlphaInelastic* theAlphaInelasticModel = new G4LEAlphaInelastic;
  theAlphaInelasticModel->SetMaxEnergy(20*MeV);
  theAlphaInelasticModel->SetMinEnergy(0*MeV);
  fTheAlphaInelasticProcess->RegisterMe(theAlphaInelasticModel);

  pManager->AddDiscreteProcess(fTheAlphaInelasticProcess);
  // ================================================================================================================

  // =========================================== GENERIC ION ========================================================
  pManager = G4GenericIon::GenericIon()->GetProcessManager();

  // need to add the elastic explicitly
  //pManager->AddDiscreteProcess(&fTheElasticProcess);

  fTheIonInelasticProcess = new G4IonInelasticProcess();

  fTheIonInelasticProcess->RegisterMe(theTheoModel);

  G4BinaryLightIonReaction * theGenIonBC= new G4BinaryLightIonReaction;
  theGenIonBC->SetMinEnergy(20.*MeV);
  theGenIonBC->SetMaxEnergy(70.*MeV);
  fTheIonInelasticProcess->RegisterMe(theGenIonBC);
  fTheIonInelasticProcess->AddDataSet(TripathiCrossSection);
  fTheIonInelasticProcess->AddDataSet(aShen);

  //Ablation & Abrasion
  G4WilsonAblationModel* ablModel = new G4WilsonAblationModel();
  G4WilsonAbrasionModel* ionModel = new G4WilsonAbrasionModel(ablModel);
  fTheIonInelasticProcess->RegisterMe(ionModel);

  pManager->AddDiscreteProcess(fTheIonInelasticProcess);
  // ================================================================================================================
}
