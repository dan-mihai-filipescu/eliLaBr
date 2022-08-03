#ifndef eliLaBr_Counter_h
#define eliLaBr_Counter_h 1

#include "globals.hh"

class eliLaBr_Counter
{
public:
	eliLaBr_Counter();
	~eliLaBr_Counter();

    void InitialiseCounter(G4long value);
	void IncrementCounter();
    G4long  GetCounterValue();
	
    G4long MyCounter;

};

#endif
