#ifndef eliLaBr_EventInformation_h
#define eliLaBr_EventInformation_h 1

#include "G4VUserEventInformation.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

#include <stdlib.h>
#include <cmath>

class eliLaBr_EventInformation : public G4VUserEventInformation
{
  public:
    eliLaBr_EventInformation();
    ~eliLaBr_EventInformation();
    inline void *operator new(size_t);
    inline void operator delete(void *aEventInfo);
    inline int operator ==(const eliLaBr_EventInformation& right) const
        {return (this==&right);}
    void Print() const;

  private:
    G4ThreeVector PrimaryPartPolarization;
    G4double eventWeight;
    G4double norm;
    G4double polAngle;
    G4int src_n_type;
    G4int NDet;

  public:
    inline G4ThreeVector Get_PrimaryPartPolarization() const { return PrimaryPartPolarization;}
    void Set_PrimaryPartPolarization(G4ThreeVector Polarization) {PrimaryPartPolarization = Polarization;}
    inline G4double Get_weight() const { return eventWeight;}
    void Set_weight(G4double weight) {eventWeight = weight;}
    inline G4double Get_norm() const { return norm;}
    void Set_norm(G4double norm_in) {norm = norm_in;}
    inline G4double Get_angle() const { return polAngle;}
    void Set_angle(G4double angle_in) {polAngle = angle_in;}
    inline G4int Get_src_n_type() const { return src_n_type;}
    void Set_src_n_type(G4int src_n_type_in) {src_n_type = src_n_type_in;}
    inline G4int Get_N_Det() const { return NDet;}
    void Set_N_Det(G4int NDet_in) {NDet = NDet_in;}
};

extern G4Allocator<eliLaBr_EventInformation> aEventInformationAllocator;

inline void* eliLaBr_EventInformation::operator new(size_t)
    { void* aEventInfo;
      aEventInfo = (void*)aEventInformationAllocator.MallocSingle();
      return aEventInfo;
    }

inline void eliLaBr_EventInformation::operator delete(void *aEventInfo)
    { aEventInformationAllocator.FreeSingle((eliLaBr_EventInformation*)aEventInfo);}


#endif
