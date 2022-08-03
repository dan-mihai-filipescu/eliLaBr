
#include "eliLaBr_Counter.hh"



eliLaBr_Counter::eliLaBr_Counter()
	{	MyCounter = 0;	}
	
eliLaBr_Counter::~eliLaBr_Counter()
	{	}



void eliLaBr_Counter::InitialiseCounter(G4long value)
	{	MyCounter = value;	}


void eliLaBr_Counter::IncrementCounter()
	{	MyCounter++;	}


G4long eliLaBr_Counter::GetCounterValue()
	{	return (MyCounter);	}

