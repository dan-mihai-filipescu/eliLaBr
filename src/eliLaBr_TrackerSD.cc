
#include "eliLaBr_TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"




eliLaBr_TrackerSD::eliLaBr_TrackerSD(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="trackerCollection");
}




eliLaBr_TrackerSD::~eliLaBr_TrackerSD(){ }




void eliLaBr_TrackerSD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new eliLaBr_TrackerHitsCollection
                          (SensitiveDetectorName,collectionName[0]); 
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, trackerCollection ); 
}




G4bool eliLaBr_TrackerSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{


  //G4double kineticEnergy = aStep->GetTrack()->GetKineticEnergy();
  G4int step_number = aStep->GetTrack()->GetCurrentStepNumber();
  G4double kineticEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();
  G4double local_time = aStep->GetTrack()->GetLocalTime();	// Time since the current track is created. 
  G4double global_time = aStep->GetTrack()->GetGlobalTime();	// Time since the event in which the track belongs is created. 
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4int trackID = aStep->GetTrack()->GetTrackID();
  G4int parentID = aStep->GetTrack()->GetParentID();
  G4String partName= aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
  G4String volName = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName(); 
  G4ThreeVector volPos = aStep->GetPreStepPoint()->GetTouchable()->GetTranslation(); 
  G4bool StStatIn, StStatOut;
  G4int replica_nb=0;

  if(aStep->GetPreStepPoint()->GetStepStatus()==fGeomBoundary) StStatIn = 1;
      else StStatIn = 0;
  if(aStep->GetPostStepPoint()->GetStepStatus()==fGeomBoundary) StStatOut = 1;
      else StStatOut = 0;

  G4int replica_nb0 = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(0); 
  G4int replica_nb1 = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(1);
//  G4int replica_nb2 = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(2); 
  if(volName=="he3_sensitive") replica_nb = replica_nb1;
  if(volName=="target") replica_nb = replica_nb0;

 // if ( edep == 0. ) return false;
/*if((volName=="target")&&(partName=="gamma")&&(trackID==1)&&(step_number==2)&&(StStatIn==1)&&(StStatOut==1))
      return false;*/

  eliLaBr_TrackerHit* newHit = new eliLaBr_TrackerHit();
  
//  newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchableHandle()
//                                               ->GetCopyNumber());
  newHit->SetKinEn    (kineticEnergy);
  newHit->SetGlobalTime(global_time);
  newHit->SetTrackID  (trackID);
  newHit->SetParentID  (parentID);
  newHit->SetEdep     (edep);
  newHit->SetPos      (aStep->GetPostStepPoint()->GetPosition());
  newHit->SetVname    (volName);
  newHit->SetPartName (partName);
  newHit->SetVolPos   (volPos);
  newHit->SetVolInStat (StStatIn);
  newHit->SetVolOutStat (StStatOut);
  newHit->SetRepNb   (replica_nb);
//  newHit->Set0RepNb   (replica_nb0);
//  newHit->Set1stRepNb (replica_nb1);
//  newHit->Set2ndRepNb (replica_nb2); 


/*  if (partName=="neutron"){
	  	  G4ThreeVector mom_pre=aStep->GetPreStepPoint()->GetMomentum();
		  G4ThreeVector mom_post=aStep->GetPostStepPoint()->GetMomentum();
		  G4ThreeVector mom_diff=mom_post-mom_pre;
		  G4double mass = mom_diff.mag2()/2./edep;
          G4cout<<"mass="<<mass<<G4endl;
		  if( 935.<mass && mass < 940.) partName="proton";
          if(11000.<mass && mass <13000.) partName="C12[0.0]";
  }*/

//  newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchableHandle()
//                                               ->GetCopyNumber());

/*  if ((partName=="neutron")  ) {
      G4cout<< "One neutron recoded!!!"<<G4endl;
     G4cout<< aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()<<G4endl;
  }*/
//if ((partName=="neutron")  )
  //if(aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="nCapture") G4cout<<"CAPTURE !!!"<<G4endl;
//  G4cout<< aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()<<G4endl;
    trackerCollection->insert( newHit );
//  G4cout<<"newhit:  KinEn="<<kineticEnergy<<"   edep="<<edep<<"    partName="<<partName<<"  trackID="<<trackID<<G4endl;
//  printf("Hit:%6s  KinEn=%9.6f\tEdep=%9.6f\tpartName=%10s\ttrackID=%3i\tparentID=%3i\tInVol=%1i\tOutVol=%1i\tDeHit=%3i\tStep#=%3i\n",volName.data(),kineticEnergy,edep,partName.data(),trackID,parentID,StStatIn,StStatOut,replica_nb,step_number);

  //newHit->Print();
  //newHit->Draw();

  return true;
}




void eliLaBr_TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel>0) { 
     G4int NbHits = trackerCollection->entries();
     G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
            << " hits in the tracker chambers: " << G4endl;
     for (G4int i=0;i<NbHits;i++) (*trackerCollection)[i]->Print();
    } 
}



