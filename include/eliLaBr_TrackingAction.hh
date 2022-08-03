#ifndef eliLaBr_TrackingAction_h
#define eliLaBr_TrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "globals.hh"


class eliLaBr_TrackingAction : public G4UserTrackingAction {

public:
	eliLaBr_TrackingAction();
	~eliLaBr_TrackingAction() {};

	void PreUserTrackingAction(const G4Track*);
	void PostUserTrackingAction(const G4Track*);

};

#endif
