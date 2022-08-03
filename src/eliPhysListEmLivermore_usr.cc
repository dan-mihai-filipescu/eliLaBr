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
/// \file eli/src/eliPhysListEmLivermore_usr.cc
/// \brief Implementation of the eliPhysListEmLivermore_usr class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "eliPhysListEmLivermore_usr.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4EmProcessOptions.hh"

// gamma
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4ComptonScattering.hh"
#include "G4KleinNishinaModel.hh"
#include "G4LowEPComptonModel.hh"
#include "G4LivermoreComptonModel.hh"

#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4PenelopeGammaConversionModel.hh"

#include "G4RayleighScattering.hh" 
#include "G4LivermoreRayleighModel.hh"

#include "G4PhotoNuclearProcess.hh"
#include "G4CascadeInterface.hh"

// e-
#include "G4eMultipleScattering.hh"
#include "G4UniversalFluctuation.hh"

#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4PenelopeIonisationModel.hh"

#include "G4eBremsstrahlung.hh"
#include "G4Generator2BS.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4eBremsstrahlungRelModel.hh"
#include "G4LivermoreBremsstrahlungModel.hh"

// e+
#include "G4eplusAnnihilation.hh"

// mu

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

// hadrons

#include "G4hMultipleScattering.hh"
#include "G4MscStepLimitType.hh"

#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4alphaIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"

// msc models
#include "G4UrbanMscModel96.hh"
#include "G4WentzelVIModel.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"

#include "G4SynchrotronRadiation.hh"
#include "G4SynchrotronRadiationInMat.hh"

#include "G4PhysicsListHelper.hh"

#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

eliPhysListEmLivermore_usr::eliPhysListEmLivermore_usr(const G4String& name)
   :  G4VPhysicsConstructor(name)
{
  G4LossTableManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

eliPhysListEmLivermore_usr::~eliPhysListEmLivermore_usr()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eliPhysListEmLivermore_usr::ConstructProcess()
{
  // Add EM Processes from G4EmLivermorePhysics builder
  
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4NuclearStopping* pnuc = new G4NuclearStopping();
  
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    //G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    G4double HighEnergyLimit = GeV;
    G4double HigheEnergyLimit = 100.*MeV;
     
    if (particleName == "gamma") {
     
      // Photoelectric
      G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
      G4LivermorePhotoElectricModel* theLivermorePhotoElectricModel = 
        new G4LivermorePhotoElectricModel();
      theLivermorePhotoElectricModel->SetHighEnergyLimit(HighEnergyLimit);
      thePhotoElectricEffect->AddEmModel(0, theLivermorePhotoElectricModel);
      //pmanager->AddDiscreteProcess(thePhotoElectricEffect);
      ph->RegisterProcess(thePhotoElectricEffect, particle);

      // Compton scattering
      G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
      theComptonScattering->SetEmModel(new G4KleinNishinaModel(),1);
      //G4VEmModel* theLowEPComptonModel = new G4LivermoreComptonModel();
      G4VEmModel* theLowEPComptonModel = new G4LowEPComptonModel();
      theLowEPComptonModel->SetHighEnergyLimit(20.*MeV);
      theComptonScattering->AddEmModel(0, theLowEPComptonModel);
      //pmanager->AddDiscreteProcess(theComptonScattering);
      ph->RegisterProcess(theComptonScattering, particle);

      // Gamma conversion
      G4GammaConversion* theGammaConversion = new G4GammaConversion();
      //G4VEmModel* theGammaConversionModel = new G4LivermoreGammaConversionModel();
      G4VEmModel* theGammaConversionModel = new G4PenelopeGammaConversionModel();
      theGammaConversionModel->SetHighEnergyLimit(HighEnergyLimit);
      theGammaConversion->AddEmModel(0, theGammaConversionModel);
      //pmanager->AddDiscreteProcess(theGammaConversion);
      ph->RegisterProcess(theGammaConversion, particle);

      // Rayleigh scattering
      G4RayleighScattering* theRayleigh = new G4RayleighScattering();
      G4LivermoreRayleighModel* theRayleighModel = new G4LivermoreRayleighModel();
      theRayleighModel->SetHighEnergyLimit(HighEnergyLimit);
      theRayleigh->AddEmModel(0, theRayleighModel);
      //pmanager->AddDiscreteProcess(theRayleigh);
      ph->RegisterProcess(theRayleigh, particle);
      
      /*G4PhotoNuclearProcess* process = new G4PhotoNuclearProcess();
      G4CascadeInterface* bertini = new G4CascadeInterface();
      bertini->SetMaxEnergy(10*GeV);
      process->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(process);*/
      
    } else if (particleName == "e-") {
  
      // multiple scattering
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      G4GoudsmitSaundersonMscModel* msc1 = new G4GoudsmitSaundersonMscModel();
      //G4UrbanMscModel96* msc1 = new G4UrbanMscModel96();
      G4WentzelVIModel* msc2 = new G4WentzelVIModel();
      msc1->SetHighEnergyLimit(HigheEnergyLimit);
      msc2->SetLowEnergyLimit(HigheEnergyLimit);
      msc->SetStepLimitType(fUseDistanceToBoundary);
      msc->SetRangeFactor(0.01);
      msc->AddEmModel(0, msc1);
      msc->AddEmModel(0, msc2);
      //pmanager->AddProcess(msc,                   -1, 1, 1);
      ph->RegisterProcess(msc, particle);

      // Single Coulomb scattering
      G4CoulombScattering* ss = new G4CoulombScattering();
      G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel(); 
      ssm->SetLowEnergyLimit(HigheEnergyLimit);
      ssm->SetActivationLowEnergyLimit(HigheEnergyLimit);
      ss->SetMinKinEnergy(HigheEnergyLimit);
      ss->SetEmModel(ssm, 1); 
      ph->RegisterProcess(ss, particle);
      
      // Ionisation
      G4eIonisation* eIoni = new G4eIonisation();
      //G4VEmModel* theIoniModel = new G4LivermoreIonisationModel();
      G4VEmModel* theIoniModel = new G4PenelopeIonisationModel();
      theIoniModel->SetHighEnergyLimit(1*MeV); 
      eIoni->AddEmModel(0, theIoniModel, new G4UniversalFluctuation() );
      eIoni->SetStepFunction(0.2, 10*um); //     
      //pmanager->AddProcess(eIoni,                 -1, 2, 2);
      ph->RegisterProcess(eIoni, particle);
      
      // Bremsstrahlung
      G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
      //G4LivermoreBremsstrahlungModel* br1 = new G4LivermoreBremsstrahlungModel();
      G4SeltzerBergerModel* br1 = new G4SeltzerBergerModel();
      G4eBremsstrahlungRelModel* br2 = new G4eBremsstrahlungRelModel();
      br1->SetHighEnergyLimit(HighEnergyLimit);
      br1->SetAngularDistribution(new G4Generator2BS());
      br2->SetAngularDistribution(new G4Generator2BS());
      eBrem->AddEmModel(0, br1);
      eBrem->AddEmModel(0, br2);
      //pmanager->AddProcess(eBrem,                 -1, 3, 3);
      ph->RegisterProcess(eBrem, particle);
      
      // Synchrotron Radiation
      //pmanager->AddProcess(new G4SynchrotronRadiation,      -1,-1, 4);
      if(this->GetPhysicsName()=="Livermore&SynchRad_usr")
        {
	printf("\nEmLivermore_usr message: %s - include Synchrotron radiation physics\n",GetPhysicsName().data());
        ph->RegisterProcess(new G4SynchrotronRadiation(), particle);
	}
            
    } else if (particleName == "e+") {

      // Identical to G4EmStandardPhysics_option3
      
      // multiple scattering
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      G4GoudsmitSaundersonMscModel* msc1 = new G4GoudsmitSaundersonMscModel();
      //G4UrbanMscModel96* msc1 = new G4UrbanMscModel96();
      G4WentzelVIModel* msc2 = new G4WentzelVIModel();
      msc1->SetHighEnergyLimit(HigheEnergyLimit);
      msc2->SetLowEnergyLimit(HigheEnergyLimit);
      msc->SetStepLimitType(fUseDistanceToBoundary);
      msc->SetRangeFactor(0.01);
      msc->AddEmModel(0, msc1);
      msc->AddEmModel(0, msc2);      
      //pmanager->AddProcess(msc,                   -1, 1, 1);
      ph->RegisterProcess(msc, particle);

      // Single Coulomb scattering
      G4CoulombScattering* ss = new G4CoulombScattering();
      G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel(); 
      ssm->SetLowEnergyLimit(HigheEnergyLimit);
      ssm->SetActivationLowEnergyLimit(HigheEnergyLimit);
      ss->SetMinKinEnergy(HigheEnergyLimit);
      ss->SetEmModel(ssm, 1);
      ph->RegisterProcess(ss, particle); 

      // Ionisation
      G4eIonisation* eIoni = new G4eIonisation();
      //G4VEmModel* theIoniModel = new G4LivermoreIonisationModel();
      G4VEmModel* theIoniModel = new G4PenelopeIonisationModel();
      theIoniModel->SetHighEnergyLimit(1*MeV);
      eIoni->AddEmModel(0, theIoniModel, new G4UniversalFluctuation() );
      eIoni->SetStepFunction(0.2, 10*um);      
      //pmanager->AddProcess(eIoni,                 -1, 2, 2);
      ph->RegisterProcess(eIoni, particle);

      // Bremsstrahlung
      G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
      //G4LivermoreBremsstrahlungModel* br1 = new G4LivermoreBremsstrahlungModel();
      G4SeltzerBergerModel* br1 = new G4SeltzerBergerModel();
      G4eBremsstrahlungRelModel* br2 = new G4eBremsstrahlungRelModel();
      br1->SetHighEnergyLimit(HighEnergyLimit);
      br1->SetAngularDistribution(new G4Generator2BS());
      br2->SetAngularDistribution(new G4Generator2BS());
      eBrem->AddEmModel(0, br1);
      eBrem->AddEmModel(0, br2);
      //pmanager->AddProcess(eBrem,                -1, 3, 3); 
      ph->RegisterProcess(eBrem, particle);     

      // annihilation at rest and in flight
      //pmanager->AddProcess(new G4eplusAnnihilation,0,-1, 4);
      ph->RegisterProcess(new G4eplusAnnihilation(), particle);

      
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {

      // Identical to G4EmStandardPhysics_option3
      
      G4MuMultipleScattering* mumsc = new G4MuMultipleScattering();
      mumsc->AddEmModel(0, new G4WentzelVIModel());
      //pmanager->AddProcess(mumsc,                       -1, 1, 1);
      ph->RegisterProcess(mumsc, particle);

      G4MuIonisation* muIoni = new G4MuIonisation();
      muIoni->SetStepFunction(0.2, 50*um);          
      //pmanager->AddProcess(muIoni,                    -1, 2, 2);
      ph->RegisterProcess(muIoni, particle);

      //pmanager->AddProcess(new G4MuBremsstrahlung,    -1, 3, 3);
      ph->RegisterProcess(new G4MuBremsstrahlung(), particle);
      //pmanager->AddProcess(new G4MuPairProduction,    -1, 4, 4);
      ph->RegisterProcess(new G4MuPairProduction(), particle);
      //pmanager->AddDiscreteProcess(new G4CoulombScattering());
      ph->RegisterProcess(new G4CoulombScattering(), particle);

    } else if (particleName == "GenericIon") {

      //pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      ph->RegisterProcess(new G4hMultipleScattering("ionmsc"), particle);

      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetEmModel(new G4IonParametrisedLossModel());
      ionIoni->SetStepFunction(0.1, 1*um);
      //pmanager->AddProcess(ionIoni,                   -1, 2, 2);
      ph->RegisterProcess(ionIoni, particle);
      
      //pmanager->AddProcess(new G4NuclearStopping(),   -1, 3,-1);
      ph->RegisterProcess(pnuc, particle);

    } else if (particleName == "alpha" ||
               particleName == "He3" ) {

      // Identical to G4EmStandardPhysics_option3
      
      G4hMultipleScattering* msc = new G4hMultipleScattering();
      //pmanager->AddProcess(msc,                       -1, 1, 1);
      ph->RegisterProcess(msc, particle);

      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetStepFunction(0.1, 10*um);
      //pmanager->AddProcess(ionIoni,                   -1, 2, 2);
      ph->RegisterProcess(ionIoni, particle);
      //pmanager->AddProcess(new G4NuclearStopping(),   -1, 3,-1);
      ph->RegisterProcess(pnuc, particle);

    } else if (particleName == "pi+" ||
               particleName == "pi-" ||
               particleName == "kaon+" ||
               particleName == "kaon-" ||
               particleName == "proton" ) {

      // Identical to G4EmStandardPhysics_option3
      
      G4hMultipleScattering* pmsc = new G4hMultipleScattering();
      pmsc->SetEmModel(new G4WentzelVIModel());
      //pmanager->AddProcess(pmsc,                      -1, 1, 1);
      ph->RegisterProcess(pmsc, particle);
      
      G4hIonisation* hIoni = new G4hIonisation();
      hIoni->SetStepFunction(0.1, 10*um);
      //pmanager->AddProcess(hIoni,                     -1, 2, 2);
      ph->RegisterProcess(hIoni, particle);
            
      //pmanager->AddProcess(new G4hBremsstrahlung,     -1,-3, 3);
      ph->RegisterProcess(new G4hBremsstrahlung(), particle);
      //pmanager->AddProcess(new G4hPairProduction,     -1,-4, 4);
      ph->RegisterProcess(new G4hPairProduction(), particle);
      ph->RegisterProcess(new G4CoulombScattering(), particle);
      ph->RegisterProcess(pnuc, particle);
    } else if (particleName == "deuteron" ||
               particleName == "triton") {

      ph->RegisterProcess(new G4hMultipleScattering(), particle);
      ph->RegisterProcess(new G4hIonisation(), particle);
      ph->RegisterProcess(pnuc, particle);
    }
  }

  // Em options
  //      
  G4EmProcessOptions opt;
  //opt.SetVerbose(verbose);
  
  // Multiple Coulomb scattering
  //
  opt.SetPolarAngleLimit(CLHEP::pi);
    
  // Physics tables
  //
  opt.SetMinEnergy(100*eV);
  opt.SetMaxEnergy(10*TeV);
  opt.SetDEDXBinning(220);
  opt.SetLambdaBinning(220);

  // Nuclear stopping
  pnuc->SetMaxKinEnergy(MeV);

  // Deexcitation
  G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);
  de->SetFluo(true);
  de->SetAuger(true);
  de->SetPIXE(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

