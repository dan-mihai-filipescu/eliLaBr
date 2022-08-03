
#include "eliLaBr_TrackerHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"


G4Allocator<eliLaBr_TrackerHit> eliLaBr_TrackerHitAllocator;

eliLaBr_TrackerHit::eliLaBr_TrackerHit()
{;}

eliLaBr_TrackerHit::~eliLaBr_TrackerHit()
{;}

eliLaBr_TrackerHit::eliLaBr_TrackerHit(const eliLaBr_TrackerHit &right)
  : G4VHit()
{
  edep = right.edep;
  pos = right.pos;
}

const eliLaBr_TrackerHit& eliLaBr_TrackerHit::operator=(const eliLaBr_TrackerHit &right)
{
  edep = right.edep;
  pos = right.pos;
  return *this;
}

G4int eliLaBr_TrackerHit::operator==(const eliLaBr_TrackerHit &right) const
{
  return (this==&right) ? 1 : 0;
}

void eliLaBr_TrackerHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(pos);
    circle.SetScreenSize(0.04);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

void eliLaBr_TrackerHit::Print()
{;}


