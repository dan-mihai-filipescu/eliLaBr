
#ifndef eliLaBr_Geometry_hh
#define eliLaBr_Geometry_hh 1

#include <G4ThreeVector.hh>

class eliLaBr_GeometryMessenger;

class eliLaBr_Geometry {

private:
G4double Rdist,ThetaAngle, PhyAngle;
G4double step;
G4double Rad1, Rad2, Rad3;
G4int Num1, Num2, Num3;
G4double Phy1, Phy2, Phy3;

eliLaBr_GeometryMessenger* geometryMessenger;  //pointer to the Messenger

public:
eliLaBr_Geometry();
~eliLaBr_Geometry();
G4double r(),theta();
//    static G4ThreeVector GetClusterPos(int nr);
//    static G4double 	 GetClusterAngle(int nr);
void GetDetPos(int nRing, int nTheta);
G4int GetNumDetectors(G4int);

void SetRad1(G4double);
void SetRad2(G4double);
void SetRad3(G4double);
void SetNum1(G4int);
void SetNum2(G4int);
void SetNum3(G4int);
void SetPhy1(G4double);
void SetPhy2(G4double);
void SetPhy3(G4double);
};
#endif
