
#ifndef eliLaBr_DetectorConstruction_H
#define eliLaBr_DetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material; 
class G4Region;

#include "G4VUserDetectorConstruction.hh"
#include "eliLaBr_Geometry.hh"
#include "eliLaBr_TrackerHit.hh"
#include "eliLaBr_TrackerSD.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "G4TwoVector.hh"

#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4Polycone.hh"
#include "G4ExtrudedSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Orb.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "G4UImanager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"

#include "globals.hh"

class eliLaBr_DetectorMessenger;

class eliLaBr_DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    eliLaBr_DetectorConstruction();
    ~eliLaBr_DetectorConstruction();

    G4VPhysicalVolume* Construct();
  
  public:
    void SetE8_Filling(G4String);
  
    void SetTarget_thickness(G4double);
    void SetTarget_diameter(G4double);
    void SetTarget_placement(G4String);
    void SetTarget_housing(G4bool);
    void SetTarget_position(G4double);
    void SetTarget_material(G4String);
    
    void SetAttenuator_thickness(G4double);
    void SetAttenuator_hole(G4double);
    void SetAttenuator_position(G4double);
    void SetAttenuator_material(G4String);
    void SetAttenuator_placement(G4bool);
    
    void SetPlace_window(G4bool);
    void SetPlace_mirror(G4bool);
    
    void Set_BADVacPres(G4double);
    
    void SetSD_percentage(G4double);
    void SetSD_material (G4String);
    void SetSD_pressure (G4double);
    
    void SetModeratorSize(G4double);
    void SetNumberOfRings(G4int);
    G4int GetNumberOfDetectors();
    
    G4double expHall_x;
    G4double expHall_y;
    G4double expHall_z;
    G4double TargetPlacement_position;
    G4double GammaSourcePos;

  private:
    G4String E8_Filling;
    
    G4double Target_thickness;
    G4double Target_diameter;
    G4String Target_placement;
    G4bool   Target_housing;
    G4double Target_position;
    G4String Target_material;
    
    G4double Attenuator_thickness;
    G4double Attenuator_hole;
    G4double Attenuator_position;
    G4String Attenuator_material;
    G4bool Attenuator_placement;
    
    G4bool Place_window;
    G4bool Place_mirror;
    
    G4double BADVacPres;
    
    G4double pres,temp,SDdensity;
    G4double SD_percentage;
    G4String SD_material;
    
    G4double ModeratorSize;
    G4int NumberOfRings;
    G4int NumberOfDetectors;
    
    G4double DipoleRadius;
    G4double DipoleLength;
    G4double DipoleAngle;
    G4double EleEne;
    
    eliLaBr_DetectorMessenger* detectorMessenger;  //pointer to the Messenger
    eliLaBr_Geometry* SetupGeometry;
    
    
    // Logical volumes
    //
    G4LogicalVolume* experimentalHall_log;

    // Physical volumes
    //
    G4VPhysicalVolume* experimentalHall_phys;

    G4Region*   fTargetRegion;
    G4Region*   fDetectorRegion;

/*
    G4int    col_set_nb;
    G4double col_set_dist;
    G4int    col_set_slit_nb;
    G4double col_set_angle_span;
    G4double col_set_phi_zero[40];

    G4double col_slit_dist;
    G4double col_slit_z, col_slit_x, col_slit_y, col_slit_aperture;
    G4double col_slit_in_z, col_slit_in_x, col_slit_angle;


    G4double col_slit_max_aperture;
    G4double col_box_wall_width;
    G4int col_iter, col_aux_int;
    G4double col_min_dist;
*/

    G4UniformMagField* fMagField;      //pointer to the magnetic field


};

#endif

