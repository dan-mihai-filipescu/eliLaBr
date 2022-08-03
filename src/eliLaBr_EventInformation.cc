#include "globals.hh"
#include "eliLaBr_EventInformation.hh"

#include "G4ios.hh"
#include <iostream>
#include <fstream>
#include "Riostream.h"

using namespace std;

eliLaBr_EventInformation::eliLaBr_EventInformation() {eventWeight = 1.; norm = 1.; polAngle = 0.; src_n_type = 0; NDet = 0;}
eliLaBr_EventInformation::~eliLaBr_EventInformation() {;}
G4Allocator<eliLaBr_EventInformation> aEventInformationAllocator;
void eliLaBr_EventInformation::Print() const
  {
      G4cout << "Event WEIGHT: " << eventWeight  << G4endl;
      G4cout << "Total NORM: " << norm  << G4endl;
      G4cout << "Pol. Angle: " << polAngle << G4endl;
      G4cout << "Neutron Gun Status: " << G4endl;
      G4cout << "Number of Detectors: " << G4endl;

  }

