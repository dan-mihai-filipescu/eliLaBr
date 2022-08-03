
#include "G4RunManager.hh"
//#include "G4UImanager.hh"

#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
//#include "PhysicsList.hh"
//#include "PhysicsListMessenger.hh"
#include "eliPhysicsList.hh"
#include "eliPhysicsListMessenger.hh"
#include "eliLaBr_DetectorConstruction.hh"
#include "eliLaBr_EventAction.hh"
#include "eliLaBr_PrimaryGeneratorAction.hh"
#include "eliLaBr_TrackerHit.hh"
#include "eliLaBr_TrackerSD.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif


#include "G4ios.hh"
using namespace std;

int main(int argc,char** argv)
{
  // Construct the default run manager
  //

    long seeds[2];
     time_t systime=time(NULL);
     seeds[0]=(long)systime;
     seeds[1]=(long)systime*G4UniformRand();
//     seeds[0]=(long)123.;
//     seeds[1]=(long)456.;
     CLHEP::HepRandom::setTheSeeds(seeds);
        CLHEP::HepRandom::showEngineStatus();
     cout<<"seeds are:"<<seeds[0]<<" and "<<seeds[1]<<endl;//print the seeds (to >garb) for debugging.

//TH1F *spectruSingle;
//spectruSingle  = new TH1F("spectruSingle", "spectruSingle", 8000, 0, 80);

	G4RunManager* runManager = new G4RunManager;

// set mandatory initialization classes


// -->> here's your detector......
//  G4VUserDetectorConstruction* detector = new eliLaBr_DetectorConstruction;
G4VUserDetectorConstruction* detector;
eliLaBr_DetectorConstruction* ELIdetector = new eliLaBr_DetectorConstruction;
detector = ELIdetector;
  runManager->SetUserInitialization(detector);
  G4cout << "\n\t Detectors DONE.....\n" << flush; 

// -->> here's your physics
  G4PhysListFactory factory;
  G4VModularPhysicsList* phys = 0;
  eliPhysicsListMessenger* mess = 0;
  G4String physName = "";

  // Physics List name defined via 2nd argument
  //if (argc==3) { physName = argv[2]; }

  // Physics List name defined via environment variable
  char* path = getenv("PHYSLIST");
  if (path) { physName = G4String(path); }

  // reference PhysicsList via its name
  if(factory.IsReferencePhysList(physName)) {
    phys = factory.GetReferencePhysList(physName);
    mess = new eliPhysicsListMessenger(0);
  }

  // local Physics List
  if(!phys) { phys = new eliPhysicsList(); mess = new eliPhysicsListMessenger(0);}

  // define physics
  runManager->SetUserInitialization(phys);
  //runManager->SetUserInitialization(new eliPhysicsList);
  G4cout << "\n\t\t Physics DONE.....\n" << flush; 

//  G4cout<<"\n\n HERE I AM !!!!!!!\n\n\n"<<flush;

  G4VUserPrimaryGeneratorAction* gen_action = new eliLaBr_PrimaryGeneratorAction(ELIdetector);

    runManager->SetUserAction(gen_action);

    runManager->SetUserAction(new eliLaBr_RunAction());

    runManager->SetUserAction(new eliLaBr_EventAction());

// Initialize G4 kernel
//    runManager->Initialize();
    G4cout << "\n\t\t\t RunManager launched.....\n" << flush; 


// get the pointer to the User Interface manager
  G4UImanager * UI = G4UImanager::GetUIpointer();  


#ifdef G4VIS_USE
// Visualization, if you choose to have it!
	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialize();
#endif


  if(argc==1)  // Define (G)UI terminal for interactive mode
  { 
#ifdef G4UI_USE_TCSH
    G4UIExecutive * ui = new G4UIExecutive(argc,argv);
#ifdef G4VIS_USE
//   UI->ApplyCommand("/control/execute vis1.mac");
//   UI->ApplyCommand("/control/execute vrml.mac");
#endif
    ui->SessionStart();
    delete ui;
#else
    G4UIsession * session = new G4UIterminal();
    session->SessionStart();
    delete session;
#endif
  }
  
  else   // Batch mode
  { 
    G4String fileName = argv[1];
    if(fileName == "dumb")
    {
      // G4UIterminal is a (dumb) terminal
      //
      G4UIsession * session = 0;
     #ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
     #ifdef G4VIS_USE
      //UI->ApplyCommand("/control/execute vis1.mac");
      //UI->ApplyCommand("/control/execute vrml.mac");
     #endif
     #else
      session = new G4UIterminal();
     #endif
      session->SessionStart();
      delete session;
    }
    else
    {
      G4String command = "/control/execute ";
      UI->ApplyCommand(command+fileName);
    }
  }
//  ____________________________________________________________________

#ifdef G4VIS_USE
	delete visManager;
#endif

  // job termination
  //
	delete runManager;
  return 0;
  
}
