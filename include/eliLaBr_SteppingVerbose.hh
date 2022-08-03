class eliLaBr_SteppingVerbose;

#ifndef eliLaBr_SteppingVerbose_h
#define eliLaBr_SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"
#include "G4UnitsTable.hh"


class eliLaBr_SteppingVerbose : public G4SteppingVerbose
{
public:

	eliLaBr_SteppingVerbose();
	~eliLaBr_SteppingVerbose();

	void 	StepInfo();
	void 	TrackingStarted();
//	const G4String GetVolName();
	void	GetVolName();
	void	GetKineticEn();

};

#endif

