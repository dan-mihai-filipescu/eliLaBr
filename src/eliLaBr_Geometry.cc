
#include "eliLaBr_Geometry.hh"
#include "eliLaBr_GeometryMessenger.hh"

eliLaBr_Geometry::eliLaBr_Geometry() 
   {
   Rad1 = 3.8 * cm;
   Rad2 = 7.0 * cm;
   Rad3 = 10. * cm;
   Num1 = 4;
   Num2 = 8;
   Num3 = 8;
   Phy1 = 3.*pi/8.;
   Phy2 = 0.;
   Phy3 = pi/8.;
   Rdist = 0.;
   ThetaAngle = 0.;
   
   // create commands for interactive definition of geometry
   geometryMessenger = new eliLaBr_GeometryMessenger(this);
   }

eliLaBr_Geometry::~eliLaBr_Geometry()
{
   delete geometryMessenger;
}

void eliLaBr_Geometry::GetDetPos(int nRing, int nTheta)
{

    G4int Nmax;
    if (nRing<1 || nRing>3 )
        G4cout<<"\n!!!! Ring number out of range, must be in [1,3] \n";
    if (nRing == 1) Nmax = Num1;
    else if (nRing == 2) Nmax = Num2;
    else Nmax = Num3;

    if (nTheta<1 || nTheta > Nmax)
        G4cout<<"\nDetector number out of range! (must be in [1,"<<Nmax<<"]\n";
    switch (nRing)
     {
      case 1:
         Rdist = Rad1;
	 step = 2.*pi/Num1;
         ThetaAngle = (nTheta-1)*step + Phy1;
         break;
      case 2:
        Rdist = Rad2;
	step = 2.*pi/Num2;
        ThetaAngle = (nTheta-1)*step + Phy2;
        break;
      case 3:
        Rdist = Rad3;
	step = 2.*pi/Num3;
        ThetaAngle = (nTheta-1)*step + Phy3;
        break;
     }
}

G4int eliLaBr_Geometry::GetNumDetectors(G4int nRing)
{
    switch (nRing)
     {
      case 1:
         return Num1;
         break;
      case 2:
        return Num2;
        break;
      case 3:
        return Num3;
        break;
     }
return 0;
}

G4double eliLaBr_Geometry::r() {return Rdist;}

G4double eliLaBr_Geometry::theta() {return ThetaAngle;}

void eliLaBr_Geometry::SetRad1(G4double val) {Rad1 = val;}
void eliLaBr_Geometry::SetRad2(G4double val) {Rad2 = val;}
void eliLaBr_Geometry::SetRad3(G4double val) {Rad3 = val;}
void eliLaBr_Geometry::SetNum1(G4int val) {Num1 = val;}
void eliLaBr_Geometry::SetNum2(G4int val) {Num2 = val;}
void eliLaBr_Geometry::SetNum3(G4int val) {Num3 = val;}
void eliLaBr_Geometry::SetPhy1(G4double val) {Phy1 = val;}
void eliLaBr_Geometry::SetPhy2(G4double val) {Phy2 = val;}
void eliLaBr_Geometry::SetPhy3(G4double val) {Phy3 = val;}
