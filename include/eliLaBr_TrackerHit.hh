
#ifndef eliLaBr_TrackerHit_h
#define eliLaBr_TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ios.hh"

class eliLaBr_TrackerHit : public G4VHit
{
  public:

      eliLaBr_TrackerHit();
      ~eliLaBr_TrackerHit();
      eliLaBr_TrackerHit(const eliLaBr_TrackerHit &right);
      const eliLaBr_TrackerHit& operator=(const eliLaBr_TrackerHit &right);
      G4int operator==(const eliLaBr_TrackerHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      void Print();

  private:
      G4double KinEn; 
      G4double gTime;
      G4double edep; 
      G4ThreeVector pos; 
      const G4VPhysicalVolume* physVol; 
      const G4LogicalVolume* pLogV; 
      G4int trackID, parentID;
      G4String volName; 
      G4String particleName; 
      G4ThreeVector volPos; 
      G4int replica_nb; //replica_nb0, replica_nb1, replica_nb2;
      G4bool VolInStatus, VolOutStatus;


//	all this U define yourself in the ....SD.cc ( eliLaBr_TrackerSD.cc in this case...) 
//	and then attach all this functions to your own trackerCollection 

  public:
      inline void SetKinEn(G4double kineticEnergy)
      		{ KinEn = kineticEnergy; }
      inline G4double GetKinEn()
      		{ return KinEn; }
      inline void SetGlobalTime(G4double globalTime)
            { gTime = globalTime; }
      inline G4double GetGlobalTime()
            { return gTime; }
      inline void SetEdep(G4double de)
      		{ edep = de; }
      inline void AddEdep(G4double de) 
      		{ edep += de; }
      inline G4double GetEdep()
      		{ return edep; }
      inline void SetPos(G4ThreeVector xyz)
      		{ pos = xyz; }
      inline G4ThreeVector GetPos()
      		{ return pos; }
      inline const G4LogicalVolume * GetLogV()
      		{ return pLogV; }
      void SetVname(G4String name) 
      		{ volName = name; }
      inline G4String GetVname() 
      		{ return volName; }
	  void SetPartName(G4String name)
		    { particleName = name; }
	  inline G4String GetPartName()
	        { return particleName; }
      void SetVolPos(G4ThreeVector pos) 
      		{ volPos = pos; }
      inline G4ThreeVector GetVolPos() 
      		{ return volPos; }
      void SetTrackID(G4int value) 
      		{ trackID=value; }
      inline G4int GetTrackID() const
      		{ return trackID; }
      void SetParentID(G4int value)
            { parentID=value; }
      inline G4int GetParentID() const
            { return parentID; }


      void SetVolInStat(G4bool VolInStat)
            {VolInStatus=VolInStat; }
      inline G4bool GetVolInStat()
            {return VolInStatus; }
      void SetVolOutStat(G4bool VolOutStat)
            {VolOutStatus=VolOutStat; }
      inline G4bool GetVolOutStat()
            {return VolOutStatus; }

/*
      void SetReplicaNb(int nr) 
      		{ replica_nb = nr; }
      inline G4int GetReplicaNumber() 
      		{ return replica_nb; }
*/

      void SetRepNb(int nr)
            { replica_nb = nr; }
      inline G4int GetRepNb()
            { return replica_nb; }

/*      void Set0RepNb(int nr)
      		{ replica_nb0 = nr; } 
      inline G4int Get0RepNb() 
      		{ return replica_nb0; } 
      void Set1stRepNb(int nr)
      		{ replica_nb1 = nr; }
      inline G4int Get1stRepNb() 
      		{ return replica_nb1; }
      void Set2ndRepNb(int nr)
      		{ replica_nb2 = nr; }
      inline G4int Get2ndRepNb() 
      		{ return replica_nb2; }		*/

};

typedef G4THitsCollection<eliLaBr_TrackerHit> eliLaBr_TrackerHitsCollection;

extern G4Allocator<eliLaBr_TrackerHit> eliLaBr_TrackerHitAllocator;

inline void* eliLaBr_TrackerHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) eliLaBr_TrackerHitAllocator.MallocSingle();
  return aHit;
}

inline void eliLaBr_TrackerHit::operator delete(void *aHit)
{
  eliLaBr_TrackerHitAllocator.FreeSingle((eliLaBr_TrackerHit*) aHit);
}

#endif
