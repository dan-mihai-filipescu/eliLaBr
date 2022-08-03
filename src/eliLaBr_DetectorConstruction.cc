
#include "eliLaBr_DetectorConstruction.hh"
#include "eliLaBr_DetectorMessenger.hh"
#include "G4NistManager.hh"

using namespace std;

eliLaBr_DetectorConstruction::eliLaBr_DetectorConstruction()
 :  experimentalHall_log(0), fMagField(0)
 
{
G4cout<<"\nDETECTOR#######################\n"<<flush;

expHall_x = 9.*m;
expHall_y = 5.*m;
expHall_z = 24.*m;

GammaSourcePos = -10.*m;

E8_Filling = "Air";
//E8_Filling = "Vacuum";

Target_thickness = 0.120*cm;
Target_diameter = 8.*mm;       //4Pi Neutron detector
//Target_diameter = 20.*mm;    //Monitor before collimator
Target_placement = "detector"; //Target is placed relative to the center of the 4Pi detector
//Target_placement = "world";    //Target is placed relative to the beggining of the WORLD
Target_housing = false;        //USE or NOT Target Housing
Target_position = 0.*mm;       //Target placement relative to the center of the Detector
//Target_position = 2.3*m;       //Target placement relative to the beggining of the WORLD
Target_material = "Lead";

Attenuator_thickness = 5.*cm;
Attenuator_hole = 0.*cm;
Attenuator_position = -2.*m;
Attenuator_material = "Lead";
Attenuator_placement = false;

Place_window = true;
Place_mirror = true;

BADVacPres = 1.e-6*atmosphere;

SD_percentage = 100.;         //Percentage of the active isotope from the Sensitive Detector
SD_material = "He3_gas";      //Set the Gas Type for the neutron counters

//SD_percentage = 96.;
//SD_material = "BF3_gas";

pres = STP_Pressure;
temp = STP_Temperature+20.*kelvin;

ModeratorSize = 36.*cm;
NumberOfRings = 3;
NumberOfDetectors = 0.;

DipoleRadius = 3.2168*m;
DipoleLength = 2.*m;
DipoleAngle = DipoleLength / DipoleRadius * radian;

// create commands for interactive definition of geometry
detectorMessenger = new eliLaBr_DetectorMessenger(this);
SetupGeometry = new eliLaBr_Geometry();

ifstream in;
/*
in.open("colimator.in",ios::in);
//Cate Seturi de colimatori?
in >> col_set_nb; cout<<"\n cate seturi: "<<col_set_nb;
//Ce distanta intre Seturile de colimatori?
in >> col_set_dist; col_set_dist = col_set_dist * cm;
        cout<<"\n dist: "<<col_set_dist<<" mm";
//Cate fante intr-un set?
in>>col_set_slit_nb; cout<<"\n cate fante: "<<col_set_slit_nb;
//Pe ce unghi se invart fantele dintr-un set?
in>>col_set_angle_span; col_set_angle_span=col_set_angle_span*deg;
        cout<<"\n atata span pe un set: "<<col_set_angle_span*180./pi<<" deg";
//Care e phi-ul fiecarui set de colimatori?
for (col_iter = 1; col_iter <= col_set_nb; col_iter ++)
{
    in >> col_set_phi_zero[col_iter];
    col_set_phi_zero[col_iter] = col_set_phi_zero[col_iter] * deg;
    cout<<"\n phi zero: "<<col_set_phi_zero[col_iter]*180./pi<<" deg";
}
//Ce distanta ai intre doua fante din acelasi set?
in>>col_slit_dist;   col_slit_dist = col_slit_dist * cm;
        cout<<"\n ce distanta intre fante: "<<col_slit_dist;
//Parametrii blocului de wolfram:
in>>col_slit_z>>col_slit_x>>col_slit_y>>col_slit_aperture;
col_slit_z = col_slit_z * cm;
col_slit_x = col_slit_x * cm;
col_slit_y = col_slit_y * cm;
col_slit_aperture = col_slit_aperture * cm;
cout<<"\n lead block Z: "<<col_slit_z<<" X: "<<col_slit_x<<" Y: "<<col_slit_y<<" slit aperture: "<<col_slit_aperture<<" in mm";
//Mancaturile din blocul de wolfram:
in>>col_slit_in_z>>col_slit_in_x>>col_slit_angle;
col_slit_in_z = col_slit_in_z * cm; //cat e mancat in grosime
col_slit_in_x = col_slit_in_x * cm; //cat e mancat in lungime
col_slit_angle = col_slit_angle * deg; //care e unghiul fantei
     cout<<"\n lead block mancat Z: "<<col_slit_in_z<<"in mm, X: "
        <<col_slit_in_x<<" in mm, unghiul in grade: "<<col_slit_angle*180./pi;
//Apertura maxima
in>>col_slit_max_aperture; col_slit_max_aperture = col_slit_max_aperture * cm;
     cout<<"\n Maximum aperture: "<<col_slit_max_aperture<< " mm";
//Grosimea peretilor cutiei in care sta fanta:
in>>col_box_wall_width; col_box_wall_width = col_box_wall_width * cm;
cout<<"\n Box walls width: "<<col_box_wall_width<< " mm";
//Distanta minima dintre doua suprafete:
in>>col_min_dist; col_min_dist = col_min_dist * cm;
cout<<"\n Minimum distance between two surfaces: "<<col_min_dist<< " mm";
in.close();
*/

G4double dummy;
in.open("incident_gamma.in",ios::in);
in >> dummy;
in >> dummy;
in >> EleEne;
in.close();
}




eliLaBr_DetectorConstruction::~eliLaBr_DetectorConstruction()
{
delete detectorMessenger;
}


G4VPhysicalVolume* eliLaBr_DetectorConstruction::Construct(){

//________MATERIALS DEFINITION:::___________________
//__________________________________________________
  G4String name,    // name
           symbol;  // symbol
  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;
  //G4int ncomponents;
  G4NistManager* man = G4NistManager::Instance();

//***Elements ______________________________________________________________________________
if ((SD_percentage<0.) || (SD_percentage>100.)) SD_percentage = 100.;

     G4Isotope* IHe3 = new G4Isotope("IHe3",2,3,3.016*g/mole);
     G4Isotope* IHe4 = new G4Isotope("IHe4",2,4,4.003*g/mole);
     G4Element* He_cust = new G4Element( "He_cust", "He_cust", 2 );
     He_cust->AddIsotope(IHe3,SD_percentage*perCent);
     He_cust->AddIsotope(IHe4,(100.-SD_percentage)*perCent);
     
     G4Isotope* IB10 = new G4Isotope("IB10",5,10,10.0129*g/mole);
     G4Isotope* IB11 = new G4Isotope("IB11",5,11,11.0093*g/mole);
     G4Element* B_cust = new G4Element( "B_cust", "B_cust", 2 );
     B_cust->AddIsotope(IB10,SD_percentage*perCent);
     B_cust->AddIsotope(IB11,(100.-SD_percentage)*perCent);

     G4Isotope* IC12 = new G4Isotope("IC12",6,12,12.0*g/mole);
     G4Element* C12 = new G4Element( "C12", "C12", 1 );
     C12->AddIsotope(IC12,1.);

     G4Isotope* IO16 = new G4Isotope("IO16",8,16,15.999*g/mole);
     G4Element* O16 = new G4Element( "O16", "O16", 1 );
     O16->AddIsotope(IO16,1.);

     //G4Element* H   = new G4Element( "H"  , "H"  , z=1. , a=  1.01*g/mole );
     //G4Element* H   = new G4Element(name="Hydrogen" , symbol="H"  , z=1. , a=  1.00794*g/mole);
     G4Element* HTSPol = new G4Element("TS_H_of_Polyethylene","H_Polyethylene", z=1., a=1.0079*g/mole);
     G4Element* Al  = new G4Element(name="Aluminum" , symbol="Al" , z=13., a= 26.98   *g/mole);
     G4Element* Fe  = new G4Element(name="Iron"     , symbol="Fe" , z=26., a= 55.845  *g/mole);
     G4Element* Cr  = new G4Element(name="Chromium" , symbol="Cr" , z=24., a= 51.996  *g/mole);
     G4Element* Ni  = new G4Element(name="Nickel"   , symbol="Ni" , z=28., a= 58.693  *g/mole);
     G4Element* Cu  = new G4Element(name="Copper"   , symbol="Cu" , z=29., a= 63.55   *g/mole);
     //G4Element* He3 = new G4Element( "He3", "He3", z=2. , a= 3.016*g/mole );
     //G4Element* C12 = new G4Element( "C12", "C12", z=6. , a= 12.00*g/mole );
     //G4Element* O16 = new G4Element( "O16", "O16", z=8. , a= 16.00*g/mole );
     G4Element* B   = new G4Element(name="Boron"    , symbol="B"  , z=5. , a= 10.81   *g/mole);
     G4Element* C   = new G4Element(name="Carbon"   , symbol="C"  , z=6. , a= 12.0107 *g/mole);
     G4Element* N   = new G4Element(name="Nitrogen" , symbol="N"  , z=7. , a= 14.01   *g/mole);
     G4Element* O   = new G4Element(name="Oxygen"   , symbol="O"  , z=8. , a= 15.999  *g/mole);
     G4Element* F   = new G4Element(name="Fluorine" , symbol="F"  , z=9. , a= 18.9984 *g/mole);
     G4Element* Na  = new G4Element(name="Sodium"   , symbol="Na" , z=11., a= 22.99   *g/mole);
     G4Element* Si  = new G4Element(name="Silicon"  , symbol="Si" , z=14., a= 28.08   *g/mole);
     G4Element* Ca  = new G4Element(name="Calcium"  , symbol="Ca" , z=20., a= 40.078  *g/mole);
     G4Element* Ti  = new G4Element(name="Titanium" , symbol="Ti" , z=22., a= 47.867  *g/mole);
     G4Element* Ge  = new G4Element(name="Germanium", symbol="Ge" , z=32., a= 72.630  *g/mole);
     G4Element* Cd  = new G4Element(name="Cadmium"  , symbol="Cd" , z=48., a=112.414  *g/mole);
     G4Element* La  = new G4Element(name="Lanthanum", symbol="La" , z=57., a=138.91   *g/mole);
     G4Element* Ce  = new G4Element(name="Cerium"   , symbol="Ce" , z=58., a=140.116  *g/mole);
     G4Element* Br  = new G4Element(name="Bromine"  , symbol="Br" , z=35., a= 79.90   *g/mole);
     G4Element* I   = new G4Element(name="Iodine"   , symbol="I"  , z=53., a=126.90   *g/mole);
     G4Element* W   = new G4Element(name="Tungsten" , symbol="W"  , z=74., a=183.84   *g/mole);
     G4Element* Au  = new G4Element(name="Gold"     , symbol="Au" , z=79., a=196.967  *g/mole);
     G4Element* Pb  = new G4Element(name="Lead"     , symbol="Pb" , z=82., a=207.2    *g/mole);
     G4Element* U   = new G4Element(name="Uranium"  , symbol="U"  , z=92., a=238.029  *g/mole);
     G4Element* Tm  = new G4Element(name="Thulium"  , symbol="Tm" , z=69., a=168.934  *g/mole);      
     G4Element* Tb  = new G4Element(name="Terbium"  , symbol="Tb" , z=65., a=158.925  *g/mole); 
     G4Element* Rh  = new G4Element(name="Rhodium"  , symbol="Rh" , z=45., a=102.905  *g/mole);
       //  G4Element* Li = new G4Element("Li"  , "Li" , z=3. , a=  6.94*g/mole );
       //  G4Element* Ge = new G4Element("Ge"  , "Ge" , z=32., a= 72.64*g/mole );
       //******************************************************************************************

//***Materials ///////////////////////////////////////////////////////////////////////////////////
     G4double p1, p2, p3, mass_numb;

     //Fe_Cr_Ni
     G4Material* Fe_Cr_Ni = new G4Material("Fe_Cr_Ni", density= 7.9*g/cm3, 3);
       Fe_Cr_Ni->AddElement(Fe, 72*perCent);
       Fe_Cr_Ni->AddElement(Cr, 18*perCent);
       Fe_Cr_Ni->AddElement(Ni, 10*perCent);

     //He3 + CO2 Gas Mixture
//     p1 = 87.2;
     p2 = 3.5;
//     p3 = 9.3;
     p3 = 2.*p2*G4Element::GetElement("Oxygen")->GetA()/G4Element::GetElement("Carbon")->GetA();
     p1 = 100. - (p2+p3);
     mass_numb = 1. / (p1*perCent / G4Element::GetElement("He_cust")->GetA() +
                 (p2+p3)*perCent / (G4Element::GetElement("Carbon")->GetA() + 2.* G4Element::GetElement("Oxygen")->GetA()) );
//     mass_numb = 1. / (p1*perCent / G4Element::GetElement("He_cust")->GetA() +
//                       p2*perCent / G4Element::GetElement("Carbon")->GetA() +
//		       p3*perCent / G4Element::GetElement("Oxygen")->GetA());

     density = mass_numb/(Avogadro*k_Boltzmann*temp/pres);
     if (SD_material == "He3_gas") SDdensity = density;
     G4Material* He3_gas = new G4Material("He3_gas", density, 3, kStateGas, temp, pres);
       He3_gas->AddElement(He_cust, p1*perCent);
       He3_gas->AddElement(C, p2*perCent);
       He3_gas->AddElement(O, p3*perCent);
     //G4Material* He3_gas = new G4Material ("U",z=92.,a=235.0439299*g/mole,density=19.1*g/cm3);

     //BF3
     p1 = 1.;
     p2 = 3.;
     mass_numb = p1 * G4Element::GetElement("B_cust")->GetA() +
                 p2 * G4Element::GetElement("Fluorine")->GetA();
     density = mass_numb/(Avogadro*k_Boltzmann*temp/pres);
     if (SD_material == "BF3_gas") SDdensity = density;
     G4Material* BF3_gas = new G4Material("BF3_gas", density, 2, kStateGas, temp, pres);
     BF3_gas->AddElement(B_cust, 1);
     BF3_gas->AddElement(F, 3);

     //Air (for low-pressure BAD Vacuum)
     p1 = 79.;
     p2 = 21.;
     mass_numb = 2. / (p1*perCent / G4Element::GetElement("Nitrogen")->GetA() + p2*perCent / G4Element::GetElement("Oxygen")->GetA() );
     density = mass_numb/(Avogadro*k_Boltzmann*temp/BADVacPres);
     G4Material* BAD_Vacuum = new G4Material("BAD_Vacuum", density, 2, kStateGas, temp, BADVacPres);
     BAD_Vacuum->AddElement(N, p1*perCent);
     BAD_Vacuum->AddElement(O, p2*perCent);
     

     //Paraffin
     G4Material* paraffin = man->FindOrBuildMaterial("G4_PARAFFIN");

     //Moderator
     //G4Material* polyethilene = man->FindOrBuildMaterial("G4_POLYETHYLENE");
     G4Material* polyethilene = new G4Material("polyethilene", density= 1.*g/cm3, 2);
       polyethilene->AddElement(HTSPol, 14.372*perCent);
       polyethilene->AddElement(C, 85.628*perCent);

       /*G4Material* polyethilene = new G4Material("Air", density= 1.29*mg/cm3, 2);
       polyethilene->AddElement(N, 70*perCent);
       polyethilene->AddElement(O, 30*perCent);*/

       
  //Calcium
  G4Material* Calcium = new G4Material ("Calcium", density = 1.55*g/cm3, 1);
  Calcium->AddElement(Ca, 100.*perCent);
  
  //Thulium
  G4Material* Thulium = new G4Material ("Thulium", density = 9.321*g/cm3, 1);
  Thulium->AddElement(Tm, 100.*perCent);

  //Terbium
  G4Material* Terbium = new G4Material ("Terbium", density = 8.22900*g/cm3, 1);
  Terbium->AddElement(Tb, 100.*perCent);  

  //Lanthanum
  G4Material* Lanthanum = new G4Material ("Lanthanum", density = 6.1540*g/cm3, 1);
  Lanthanum->AddElement(La, 100.*perCent);    
  
  //Rhodium
  G4Material* Rhodium = new G4Material ("Rhodium", density = 12.41*g/cm3, 1);
  Rhodium->AddElement(Rh, 100.*perCent);

  //Titanium
  G4Material* Titanium = new G4Material ("Titanium", density = 4.506*g/cm3, 1);
  Titanium->AddElement(Ti, 100.*perCent);

  //Cadmium
  G4Material* Cadmium  = new G4Material( "Cadmium" , density = 8.64 *g/cm3, 1);
  Cadmium->AddElement(Cd, 100.*perCent);

  //Aluminum
  G4Material* Aluminum = new G4Material ("Aluminum", density=2.7*g/cm3, 1);
  Aluminum->AddElement(Al, 100.*perCent);

  //Silicon
  G4Material* Silicon = new G4Material ("Silicon", density= 2.329*g/cm3 , 1);
  Silicon -> AddElement (Si, 100.*perCent);

  //Iron
  G4Material* Iron = new G4Material ("Iron", density = 7.86*g/cm3, 1);
  Iron->AddElement(Fe, 100.*perCent);

  //Nickel
  G4Material* Nickel = new G4Material ("Nickel", density= 8.908*g/cm3, 1);
  Nickel->AddElement(Ni, 100.*perCent);

  //Copper
  G4Material* Copper = new G4Material ("Copper",density=8.96*g/cm3, 1);
  Copper->AddElement(Cu, 100.*perCent);

  //Wolfram
  G4Material* Wolfram  = new G4Material( "Wolfram" , density = 19.25 *g/cm3, 1);
  Wolfram->AddElement(W, 100.*perCent);
  
  //Gold
  G4Material* Gold = new G4Material ("Gold", density= 19.3*g/cm3, 1);
  Gold->AddElement(Au, 100.*perCent);

  //Lead
  G4Material* Lead = new G4Material ("Lead", density= 11.34*g/cm3, 1);
  Lead->AddElement(Pb, 100.*perCent);

  //Uranium
  G4Material* Uranium = new G4Material ("Uranium", density= 19.1*g/cm3, 1);
  Uranium->AddElement(U, 100.*perCent);

  //Vacuum
  G4Material* Vacuum = new G4Material("Vacuum",z=1.,a=1.01*g/mole,
              density=universe_mean_density,kStateGas,0.1*kelvin,
              1.e-19*pascal);

  //Air
  G4Material* Air = new G4Material("Air", density= 1.29*mg/cm3, 2);
  Air->AddElement(N, 70*perCent);
  Air->AddElement(O, 30*perCent);

  //SiLi
  G4Material* SiLi = new G4Material ("SiLi", density= 2.330*g/cm3 , 1);
  SiLi -> AddElement (Si, 100.*perCent);


  //HPGe
  //G4Material* HPGe = new G4Material ("HPGe", z=32., a= 72.59*g/mole, density= 5.33*g/cm3);
  G4Material* HPGe = new G4Material ("HPGe", density= 5.323*g/cm3, 1);
  HPGe->AddElement(Ge, 100.*perCent);
  
  //LaBr3
  G4Material* LaBr3 = new G4Material ("LaBr3", density= 5.29*g/cm3 , 2);
  LaBr3->AddElement( La, 1);
  LaBr3->AddElement( Br, 3);
  
   //CeBr3
  G4Material* CeBr3 = new G4Material ("CeBr3", density= 5.1*g/cm3 , 2);
  CeBr3->AddElement( Ce, 1);
  CeBr3->AddElement( Br, 3);


   //NaI
  G4Material* NaI = new G4Material ("NaI", density= 3.67*g/cm3 , 2);
  NaI->AddElement( Na, 1);
  NaI->AddElement(  I, 3);

  //Fire bricks
  G4Material* Fire_brick = new G4Material ("Fire_brick", density= 2.00*g/cm3 , 3);
  Fire_brick->AddElement( Al, 12.*perCent);
  Fire_brick->AddElement(  O, 65.*perCent);
  Fire_brick->AddElement( Si, 23.*perCent);

  //Borosilicate glass: SiO2-80.6% B2O3-12.5% Na2O-4.2% Al2O3-2.2%
  G4Material* SiO2 = new G4Material ("SiO2", density= 2.196*g/cm3, 2);
  SiO2->AddElement( Si, 1);
  SiO2->AddElement( O, 2);
  G4Material* B2O3 = new G4Material ("B2O3", density= 3.13*g/cm3, 2);
  B2O3->AddElement( B, 2);
  B2O3->AddElement( O, 3);
  G4Material* Na2O = new G4Material ("Na2O", density= 2.27*g/cm3, 2);
  Na2O->AddElement( Na, 2);
  Na2O->AddElement( O, 1);
  G4Material* Al2O3 = new G4Material ("Al2O3", density= 3.987*g/cm3, 2);
  Al2O3->AddElement( Al, 2);
  Al2O3->AddElement( O, 3);
  G4Material* Borosilicate_glass = new G4Material ("Borosilicate_glass", density= 2.23*g/cm3 , 4);
  Borosilicate_glass->AddMaterial(SiO2 , 80.6 *perCent);
  //Borosilicate_glass->AddMaterial(B2O3 , 12.5 *perCent);
  Borosilicate_glass->AddMaterial(B2O3 , 13.0 *perCent);
  Borosilicate_glass->AddMaterial(Na2O ,  4.2 *perCent);
  Borosilicate_glass->AddMaterial(Al2O3,  2.2 *perCent);

  //Concrete
  /*G4Material* concrete = new G4Material(name="Concrete", density=2.3*g/cm3, ncomponents=6);
  concrete->AddElement(Si, 0.227915);
  concrete->AddElement(O, 0.60541);
  concrete->AddElement(H, 0.09972);
  concrete->AddElement(Ca, 0.04986);
  concrete->AddElement(Al, 0.014245);
  concrete->AddElement(Fe, 0.00285);*/
  G4Material* concrete = man->FindOrBuildMaterial("G4_CONCRETE");

  //Kapton
/*  G4Material* Kapton = new G4Material ("Kapton", density=1.42*g/cm3, 4);
  Kapton -> AddElement ( C, 22 );
  Kapton -> AddElement ( H, 10 );
  Kapton -> AddElement ( N,  2 );
  Kapton -> AddElement ( O,  4 );*/

  //Bakelite
/*  G4Material* Bakelite = new G4Material ("Bakelite", density=1.45*g/cm3, 3 );
  Bakelite -> AddElement ( H, 9 );
  Bakelite -> AddElement ( C, 9 );
  Bakelite -> AddElement ( O, 1 );*/

  //Glass
/*  G4Material* Glass = new G4Material("Glass", density=1.032*g/cm3,2);
  Glass->AddElement(C,91.533*perCent);
  Glass->AddElement(H,8.467*perCent);*/

//______ END OF MATERIALS DEFINITION _______________________________________
////////////////////////////////////////////////////////////////////////////
 
 
	 ////////////////////////////////////////////////////////////////////////////////
	 // ______________________________ COLOURS ______________________________________

		 G4Colour red		(1.,0.,0.);			// Define red color
		 G4Colour green		(0.,1.,0.);		// Define green color
		 G4Colour blue		(0.,0.,1.);		// Define blue color
		 G4Colour violet	(1.,0.,1.);
		 G4Colour turquoise	(0.,1.,1.);
		 G4Colour white		(1.,1.,1.);
		 G4Colour yellow	(1.,1.,0.);
		 G4Colour orange	(1.,0.5,0.);
		 G4Colour grey		(0.5,0.5,0.5);
		 G4VisAttributes attred(red);        	// Define a red visualization attribute
		 G4VisAttributes attgreen(green);
		 G4VisAttributes attblue(blue);
		 G4VisAttributes attviolet(violet);
		 G4VisAttributes attturq(turquoise);
		 G4VisAttributes* attturq_filled = new G4VisAttributes(G4Colour( 0., 1., 1., 1.));
		 		 attturq_filled->SetForceSolid(true);
		 G4VisAttributes attwhite(white);
		 G4VisAttributes attyellow(yellow);
		 G4VisAttributes attorange(orange);
		 G4VisAttributes attgrey(grey);
	 //////////////////////////////////////////////////////////////////////////////////



// ______________________________________________________________________________
// --------------------------------------- VOLUMES ------------------------------
//_______________________________________________________________________________



//_________________________________ EXPERIMENTAL HALL ______________________________


G4Box* experimentalHall_box = new G4Box	("expHall_box", 
  				      expHall_x,  expHall_y,  expHall_z); 
experimentalHall_log = new G4LogicalVolume	(experimentalHall_box, 
                                              Vacuum, "expHall_log", 0, 0, 0);
experimentalHall_log->SetVisAttributes(G4VisAttributes::Invisible);
experimentalHall_phys = new G4PVPlacement	(0,
                          G4ThreeVector(0.,0.,0.), experimentalHall_log,
					      "expHall", 0, false, 10);

// __________________________ END OF EXPERIMENTAL HALL ______________________________

//________________________________ TARGET GEOMETRY __________________________________
// Detector parameters
//G4double const inch_def = 25.4;
//G4double const det_radius = 1.5*inch_def*mm;
//G4double const det_length = 3.*inch_def*mm;
G4double target_x = Target_diameter;
//G4double target_y = 0.4*cm;
// *** TARGET THICKNESS ***
G4double target_z = Target_thickness;

//G4Box* target_box = new G4Box ("target_box",target_x/2.,target_y/2.,target_z/2.);
G4Tubs* target_tub = new G4Tubs("target_tub", 0.*mm , target_x/2. , target_z/2., 0.0, 360*degree);
//G4Tubs* target_tub = new G4Tubs("target_tub", 0.*mm , det_radius , det_length/2., 0.0, 360.*degree);

G4RotationMatrix rotmt  = G4RotationMatrix();
rotmt.rotateY(45.*deg);
rotmt.rotateZ(45.*deg);
G4Transform3D transformt = G4Transform3D(rotmt,G4ThreeVector(0.,0.,0.));

G4LogicalVolume* target_log;
G4Material* TargetMaterial;
G4Material* AttenuatorMaterial;
target_log = new G4LogicalVolume (target_tub,TargetMaterial->GetMaterial(Target_material,true),       "target_log");

/*     if (Target_material == "Air")        target_log = new G4LogicalVolume (target_tub,Air,       "target_log");
else if (Target_material == "Vacuum")     target_log = new G4LogicalVolume (target_tub,Vacuum,    "target_log");
else if (Target_material == "Titanium")   target_log = new G4LogicalVolume (target_tub,Titanium,  "target_log");
else if (Target_material == "Cadmium")    target_log = new G4LogicalVolume (target_tub,Cadmium,   "target_log");
else if (Target_material == "Aluminum")   target_log = new G4LogicalVolume (target_tub,Aluminum,  "target_log");
else if (Target_material == "Iron")       target_log = new G4LogicalVolume (target_tub,Iron,      "target_log");
else if (Target_material == "Nickel")     target_log = new G4LogicalVolume (target_tub,Nickel,    "target_log");
else if (Target_material == "Copper")     target_log = new G4LogicalVolume (target_tub,Copper,    "target_log");
else if (Target_material == "Wolfram")    target_log = new G4LogicalVolume (target_tub,Wolfram,   "target_log");
else if (Target_material == "Gold")       target_log = new G4LogicalVolume (target_tub,Gold,      "target_log");
else if (Target_material == "Lead")       target_log = new G4LogicalVolume (target_tub,Lead,      "target_log");
else if (Target_material == "Uranium")    target_log = new G4LogicalVolume (target_tub,Uranium,   "target_log");
else if (Target_material == "LaBr3")      target_log = new G4LogicalVolume (target_tub,LaBr3,     "target_log");
else if (Target_material == "CeBr3")      target_log = new G4LogicalVolume (target_tub,CeBr3,     "target_log");*/

fTargetRegion = new G4Region("Target");
fTargetRegion->AddRootLogicalVolume(target_log);
target_log->SetVisAttributes(G4Colour(1.,0.,0.));

/* ---------------- Placement of the TARGET directly into the WORLD -----------------*/
//new G4PVPlacement (transformt,
//new G4PVPlacement (0,G4ThreeVector(0.,0.,2992.*cm+det_length/2.-expHall_z+5.*m),
//new G4PVPlacement (0,G4ThreeVector(0.,0.,1000. *cm-expHall_z+10.*m),
//                                 target_log,"target",experimentalHall_log, false, 0);
/*------------------------------------------------------------------------------------*/

/* ---------------- TARGET Housing --------------------------------------------------*/
G4double TGhousing_thick = 1.*mm;
G4double TGhousing_space = 0.1*mm;
//Target housing
G4double zTGhousingCoord[4]     = {0.,
                                             TGhousing_thick,
                                             TGhousing_thick,
                                             TGhousing_thick+TGhousing_space+target_z};
G4double RinTGhousingCoord [4]  = {0., 0., TGhousing_space+target_x/2., TGhousing_space+target_x/2.};
G4double RoutTGhousingCoord [4] = {TGhousing_thick+TGhousing_space+target_x/2., TGhousing_thick+TGhousing_space+target_x/2.,
					     TGhousing_thick+TGhousing_space+target_x/2., TGhousing_thick+TGhousing_space+target_x/2.};
G4Polycone* TGhousing = new G4Polycone("TGhousing",0.*deg,360.*deg,4,
                                                  zTGhousingCoord,RinTGhousingCoord,RoutTGhousingCoord);
G4LogicalVolume* TGhousing_log = new G4LogicalVolume ( TGhousing, Aluminum, "TGhousing_log");
TGhousing_log->SetVisAttributes(G4Colour(1.,1.,1.));

//_____________________________END OF TARGET GEOMETRY _______________________________


//=====================LASER_Mirror & Borosilicate_window========================
// ___________________________LASER_Mirror_______________________________________
G4Tubs* LASER_Mirror_tub = new G4Tubs("LASER_Mirror_tub", 0.*mm , 4.*cm , 1.*mm, 0.0, 360*degree);
G4LogicalVolume* LASER_Mirror_log = new G4LogicalVolume ( LASER_Mirror_tub, Silicon, "LASER_Mirror_log");
G4RotationMatrix* LASER_Mirror_rot = new G4RotationMatrix();
LASER_Mirror_rot-> rotateY(-45.*deg);
if(Place_mirror) {
//G4VPhysicalVolume* LASER_Mirror = 
new G4PVPlacement ( LASER_Mirror_rot, G4ThreeVector(.0*cm, .0*cm, (1547.+5.-20.-100.) *cm + GammaSourcePos),
                          LASER_Mirror_log, "LASER_Mirror", experimentalHall_log, false, 0);
}

// ___________________________Borosilicate_window________________________________
G4Tubs* Borosilicate_window_tub = new G4Tubs("Borosilicate_window_tub", 0.*mm , 4.*cm , 2.25*mm, 0.0, 360*degree);
G4LogicalVolume* Borosilicate_window_log = new G4LogicalVolume ( Borosilicate_window_tub, Borosilicate_glass, "Borosilicate_window_log");
if(Place_window) {
//G4VPhysicalVolume* Borosilicate_window = 
new G4PVPlacement ( 0, G4ThreeVector(.0*cm, .0*cm, (1547.+5.-20.) *cm + GammaSourcePos),
                          Borosilicate_window_log, "Borosilicate_window", experimentalHall_log, false, 0);
}

//================== END OF LASER_Mirror & Borosilicate_window ==================


//_____________________________COLLIMATORS___________________________________________
// ___________________________dummy_____________________________________________

G4Tubs* dummy_tub = new G4Tubs("dummy_tub", 0.*mm , 150.*mm , 10./2.*mm, 0.0, 360*degree);
//G4LogicalVolume* dummy_log = 
new G4LogicalVolume ( dummy_tub, Lead, "dummy_log");

/*G4VPhysicalVolume* dummy1 = new G4PVPlacement ( 0, G4ThreeVector(.0, .0, -expHall_z + 10.*mm),
                          dummy_log, "Dummy1", experimentalHall_log, false, 0);

G4VPhysicalVolume* dummy2 = new G4PVPlacement ( 0, G4ThreeVector(.0, .0, expHall_z - 10.*mm),
                          dummy_log, "Dummy2", experimentalHall_log, false, 0);

G4VPhysicalVolume* dummy3 = new G4PVPlacement ( 0, G4ThreeVector(.0, .0, 4.5*m),
                          dummy_log, "Dummy3", experimentalHall_log, false, 0);*/

// ___________________________ABSorber___________________________________________
G4VSolid* ABS_box = new G4Box ( "ABS_Box", 30.*mm, 15.*mm, 10.*mm );
G4VSolid* ABS_cylinder = new G4Tubs ( "ABS_Cylinder", 0.*mm, 20.*mm, 40./2.*mm, 0.*deg, 360.*deg );
G4VSolid* ABS_subtraction = new G4SubtractionSolid ("ABS_subtraction", ABS_cylinder, ABS_box, 0, G4ThreeVector(0.,0.,0.));

G4LogicalVolume* ABS_subtraction_log = new G4LogicalVolume ( ABS_subtraction, Copper, "ABS_subtraction_log");

G4RotationMatrix* ABS_rot = new G4RotationMatrix();
ABS_rot-> rotateY(90.*deg);
ABS_rot-> rotateX(90.*deg);

//G4VPhysicalVolume* ABS = 
new G4PVPlacement ( ABS_rot, G4ThreeVector(.0*cm, .0*cm, 1300. *cm + GammaSourcePos),
                          ABS_subtraction_log, "ABS", experimentalHall_log, false, 0);



// ___________________________Beam Shutter___________________________________________
G4Tubs* beam_shutter_tub = new G4Tubs("beam_shutter_tub", 25.*mm , 150.*cm , 300./2.*mm, 0.0, 360*degree);
G4LogicalVolume* beam_shutter_log = new G4LogicalVolume ( beam_shutter_tub, Lead, "beam_shutter_log");

//G4VPhysicalVolume* beam_shutter = 
new G4PVPlacement ( 0, G4ThreeVector(.0*cm, .0*cm, (1707.+15.)*cm + GammaSourcePos),
                          beam_shutter_log, "BeamShutter", experimentalHall_log, false, 0);


// ___________________________Collimator_1___________________________________________

G4Tubs* fix_col_tub = new G4Tubs("fix_col_tub", 1.5*mm , 150.*mm , 100./2.*mm, 0.0, 360*degree);
G4LogicalVolume* fix_col_log = new G4LogicalVolume ( fix_col_tub, Lead, "fix_col_log");

//G4VPhysicalVolume* fix_col = 
new G4PVPlacement ( 0, G4ThreeVector(.0*cm, .0*cm, (1547.+5.)*cm + GammaSourcePos),
                          fix_col_log, "FixCollimator", experimentalHall_log, false, 0);

// ___________________________Collimator_2___________________________________________

G4Tubs* col_tub = new G4Tubs("col_tub", 1.0*mm , 150.*mm , 100./2.*mm, 0.0, 360*degree);
G4LogicalVolume* col_log = new G4LogicalVolume ( col_tub, Lead, "col_log");

//G4VPhysicalVolume* col = 
new G4PVPlacement ( 0, G4ThreeVector(.0*cm, .0*cm, (1847.+5.)*cm + GammaSourcePos),
                          col_log, "Collimator", experimentalHall_log, false, 0);   
//____________________________________________________________________________________

//------------------------------ ELI-NP COLLIMATORS ----------------------------------------------
// _____________________________ Define one slit ____________________________________
/*
    G4double col_set_mother_rad, col_set_mother_z;
    G4double col_slit_mother_x, col_slit_mother_y, col_slit_mother_z;
    G4double col_slit_case_x, col_slit_case_y, col_slit_case_z;
    G4double col_slit_case_x_subtract, col_slit_case_y_subtract, col_slit_case_z_subtract;
    col_slit_case_x = 2. *col_slit_x + col_slit_max_aperture
                    + 2. * col_box_wall_width + 3. * col_min_dist;
    col_slit_case_x_subtract = 2. *col_slit_x + col_slit_max_aperture
                    + 3. * col_min_dist;
    col_slit_case_y = col_slit_y + 2. * col_min_dist + 2. * col_box_wall_width;
    col_slit_case_y_subtract = col_slit_y + 2. * col_min_dist;
    col_slit_case_z = col_slit_z;
    col_slit_case_z_subtract = 1.1 * col_slit_z;
    cout<<"\n Dimensiunile cutiei:  "<<col_slit_case_x<<"  "<<col_slit_case_y<<"  "<<col_slit_case_z<<" in mm :)\n";

    G4Box* col_slit_case_big = new G4Box	("col_slit_case_big",
                          col_slit_case_x/2.,  col_slit_case_y/2.,  col_slit_case_z/2.);
    G4Box* col_slit_case_small = new G4Box	("col_slit_case_small",
                          col_slit_case_x_subtract/2.,  col_slit_case_y_subtract/2.,
                                             col_slit_case_z_subtract/2.);

    G4VSolid* col_slit_case_subtraction = new G4SubtractionSolid ("col_slit_case_subtraction",
                       col_slit_case_big, col_slit_case_small, 0, G4ThreeVector(0.,0.,0.));
    G4LogicalVolume* col_slit_case_log = new G4LogicalVolume (col_slit_case_subtraction,
                         Fe_Cr_Ni, "col_slit_case_log");
    col_slit_case_log->SetVisAttributes(G4Colour(0.9,0.9,0.9));
    col_slit_mother_x = col_slit_case_x + 2. * col_min_dist;
    col_slit_mother_y = col_slit_case_y + 2. * col_min_dist;
    col_slit_mother_z = col_slit_case_z + 2. * col_min_dist;
    G4Box* col_slit_mother_vol = new G4Box	("col_slit_mother_vol",
                          col_slit_mother_x/2.,  col_slit_mother_y/2.,  col_slit_mother_z/2.);
    G4LogicalVolume* col_slit_mother_log = new G4LogicalVolume (col_slit_mother_vol,
                         Vacuum, "col_slit_mother_log");
    col_slit_mother_log->SetVisAttributes(G4VisAttributes::Invisible);

    new G4PVPlacement ( 0, G4ThreeVector(0.,0.,0.),
            col_slit_case_log, "col_slit_case", col_slit_mother_log, false, 0);

    //Bucatile de plumb ale fantei
    G4Box* col_block_box = new G4Box	("col_block_box",
                          col_slit_x/2.,  col_slit_y/2.,  col_slit_z/2.);

    G4double col_aux;
    col_aux = col_slit_z - col_slit_in_z * 2.;
    col_aux = col_aux * tan(col_slit_angle);

    std::vector<G4TwoVector> polygon;
    polygon.push_back(G4TwoVector( 0.*mm, -col_slit_z/2.));
    polygon.push_back(G4TwoVector( 0.*mm,  col_slit_z/2.));
    polygon.push_back(G4TwoVector( col_slit_x-col_slit_in_x, col_slit_z/2.));
    polygon.push_back(G4TwoVector( col_slit_x-col_slit_in_x, col_slit_z/2.-col_slit_in_z));
    polygon.push_back(G4TwoVector( col_slit_x-col_aux, col_slit_z/2.-col_slit_in_z));
    polygon.push_back(G4TwoVector( col_slit_x,-(col_slit_z/2.-col_slit_in_z)));
    polygon.push_back(G4TwoVector( col_slit_x-col_slit_in_x,-(col_slit_z/2.-col_slit_in_z)));
    polygon.push_back(G4TwoVector( col_slit_x-col_slit_in_x,-col_slit_z/2.));

    std::vector<G4ExtrudedSolid::ZSection> zsections;
    zsections.push_back(G4ExtrudedSolid::ZSection(-col_slit_y/2, {0.,  0.}, 1.));
    zsections.push_back(G4ExtrudedSolid::ZSection( col_slit_y/2, {0.,  0.}, 1.));

    G4ExtrudedSolid* col_block_test = new G4ExtrudedSolid("col_extruded", polygon, zsections);
    G4LogicalVolume* col_block_test_log = new G4LogicalVolume (col_block_test,
                         Wolfram, "col_block_test_log");
    G4RotationMatrix* col_block_test_rot = new G4RotationMatrix();
    col_block_test_rot->rotateX(pi/2.);
    new G4PVPlacement ( col_block_test_rot, G4ThreeVector(-1./2.*col_slit_aperture-col_slit_x,0.,0.),
            col_block_test_log, "col_block", col_slit_mother_log, false, 0);
    G4RotationMatrix* wtf = new G4RotationMatrix();
    wtf->rotateX(pi/2.);
    wtf->rotateY(pi);
    new G4PVPlacement ( wtf, G4ThreeVector( 1./2.*col_slit_aperture+col_slit_x,0.,0.),
            col_block_test_log, "col_block", col_slit_mother_log, false, 0);



    G4LogicalVolume* col_block_log = new G4LogicalVolume (col_block_box,
                         Wolfram, "col_block_log");
    new G4PVPlacement ( 0, G4ThreeVector(-1./2.*col_slit_aperture-1./2.*col_slit_x,0.,0.),
            col_block_log, "col_block", col_slit_mother_log, false, 0);
    new G4PVPlacement ( 0, G4ThreeVector(1./2.*col_slit_aperture+1./2.*col_slit_x,0.,0.),
            col_block_log, "col_block", col_slit_mother_log, false, 0);

    col_set_mother_rad = col_slit_mother_x/2. * col_slit_mother_x/2.
                        + col_slit_mother_y/2. * col_slit_mother_y/2.;
    col_set_mother_rad = sqrt(col_set_mother_rad) + col_min_dist;
    col_set_mother_z = col_set_slit_nb * col_slit_mother_z +
              col_slit_dist * (col_set_slit_nb - 1) + 2. * col_min_dist;
    G4Tubs* col_set_mother_vol = new G4Tubs("col_set_mother_vol", 0.* cm ,
     col_set_mother_rad, col_set_mother_z/2., 0.0, 360*degree);
    G4LogicalVolume* col_set_mother_log = new G4LogicalVolume (col_set_mother_vol,
                         Vacuum, "col_set_mother_log");
    col_set_mother_log->SetVisAttributes(G4VisAttributes::Invisible);

    G4RotationMatrix* col_slit_rot[100];
    G4double col_slit_pos_z[100];
    for (col_iter = 1; col_iter <= col_set_slit_nb; col_iter ++)
        {
        col_slit_rot[col_iter] = new G4RotationMatrix();
        col_slit_rot[col_iter]->rotateZ(col_set_angle_span/col_set_slit_nb*(col_iter-1));
        col_slit_pos_z[col_iter] = -col_set_mother_z/2. + col_min_dist
            + (col_iter-1)* col_slit_dist + (col_iter-1./2.)*col_slit_mother_z;
        new G4PVPlacement ( col_slit_rot[col_iter], G4ThreeVector(0.,0.,col_slit_pos_z[col_iter]),
               col_slit_mother_log, "col_slit_mother", col_set_mother_log, false, 0);
        }


    G4RotationMatrix* col_set_rot[100];
    for (col_iter = 1; col_iter <= col_set_nb; col_iter ++)
        {
        col_set_rot[col_iter] = new G4RotationMatrix();
        col_set_rot[col_iter]->rotateZ(col_set_phi_zero[col_iter]);
        new G4PVPlacement ( col_set_rot[col_iter], G4ThreeVector(0.,0.,-expHall_z+889.7*cm+
                (col_iter-1./2.)*col_set_mother_z+(col_iter-1)*col_set_dist),
                col_set_mother_log, "col_set_mother", experimentalHall_log, false, 0);
        }
*/
//------------------------------ END OF ELI-NP COLLIMATORS ----------------------------------------------------------------
//============================== END OF COLLIMATORS =======================================================================

//============================== NewSUBARU DIPOLE Magnet ==================================================================
G4Tubs* DipoleMagnet_tub = new G4Tubs("DipoleMagnet_tub", DipoleRadius-200.*mm , DipoleRadius+200.*mm , 200.*mm, 180.*degree-DipoleAngle, DipoleAngle);
G4LogicalVolume* DipoleMagnet_log = new G4LogicalVolume ( DipoleMagnet_tub, Vacuum, "DipoleMagnet_log");
DipoleMagnet_log->SetVisAttributes(G4Colour(0.,0.,1.));

G4double fieldStrength = sqrt(2.*0.511*EleEne+EleEne*EleEne)/(299.792458*DipoleRadius/m);
G4cout<<"Dipole magnetic field: "<<fieldStrength<<" Tesla"<<G4endl;
fieldStrength = fieldStrength*tesla;
if( fMagField ) delete fMagField; //delete the existing mag field
fMagField = new G4UniformMagField(G4ThreeVector(0., fieldStrength, 0.));
G4FieldManager* fieldMgr = new G4FieldManager(fMagField);
fieldMgr->SetDetectorField(fMagField);
fieldMgr->CreateChordFinder(fMagField);
DipoleMagnet_log->SetFieldManager(fieldMgr, true);

G4RotationMatrix* DipoleMagnet_rot = new G4RotationMatrix();
DipoleMagnet_rot-> rotateX(-90.*deg);
//G4VPhysicalVolume* DipoleMagnet = 
new G4PVPlacement ( DipoleMagnet_rot, G4ThreeVector(DipoleRadius, .0*cm, (1847.+5.-841.)*cm + GammaSourcePos),
                          DipoleMagnet_log, "DipoleMagnet", experimentalHall_log, false, 0);
/*DipoleMagnet_rot-> rotateX(0.*deg);
G4VPhysicalVolume* DipoleMagnet = new G4PVPlacement ( DipoleMagnet_rot, G4ThreeVector(DipoleRadius, 0*cm, (1847.+5.-855.)*cm + GammaSourcePos),
                          DipoleMagnet_log, "DipoleMagnet", experimentalHall_log, false, 0);*/
//============================== END OF NewSUBARU DIPOLE Magnet ===========================================================

//================================= Electron beamline =====================================================================
G4Tubs* ElectronBeamline_tub = new G4Tubs("ElectronBeamline_tub", 0.*cm, 5.*cm, 1010.*cm, 0.*degree, 360.*degree);
G4LogicalVolume* ElectronBeamline_log = new G4LogicalVolume ( ElectronBeamline_tub, BAD_Vacuum, "ElectronBeamline_log");
ElectronBeamline_log->SetVisAttributes(G4Colour(0.,0.,0.5));
//G4VPhysicalVolume* ElectronBeamline = 
new G4PVPlacement ( 0, G4ThreeVector(0.*cm, 0.*cm, 0.*cm + GammaSourcePos),
                          ElectronBeamline_log, "ElectronBeamline", experimentalHall_log, false, 0);

G4Tubs* ElectronBeamlinePipe_tub = new G4Tubs("ElectronBeamlinePipe_tub", 5.01*cm, 5.16*cm, 1010.*cm, 0.*degree, 360.*degree);
G4LogicalVolume* ElectronBeamlinePipe_log = new G4LogicalVolume ( ElectronBeamlinePipe_tub, Fe_Cr_Ni, "ElectronBeamlinePipe_log");
ElectronBeamlinePipe_log->SetVisAttributes(G4Colour(0.5,0.5,0.5));
//G4VPhysicalVolume* ElectronBeamlinePipe = 
new G4PVPlacement ( 0, G4ThreeVector(0.*cm, 0.*cm, 0.*cm + GammaSourcePos),
                          ElectronBeamlinePipe_log, "ElectronBeamlinePipe", experimentalHall_log, false, 0);
//================================= END OF Electron beamline ==============================================================

//============================== ELI GEOMETRY =============================================================================

//Primul perete de caramida
G4Box* wall_Fire_brick_big_box   = new G4Box ("wall_Fire_brick_big_box",  2.*m/2.,2.*m/2.,150.*cm/2.);
G4Box* wall_Fire_brick_small_box = new G4Box ("wall_Fire_brick_small_box",2.*m/2.,2.*m/2.,100.*cm/2.);
G4Box* wall_concrete_big_box     = new G4Box ("wall_concrete_big_box",    2.*m/2.,2.*m/2.,200.*cm/2.);
G4Box* wall_concrete_small_box   = new G4Box ("wall_concrete_small_box",  2.*m/2.,2.*m/2.,100.*cm/2.);
G4Tubs* wall_Fire_brick_hole = new G4Tubs("wall_Fire_brick_hole", 0.*mm , 11.*mm/2. , 160./2.*cm, 0.0, 360*degree);
G4Tubs* wall_concrete_hole = new G4Tubs("wall_concrete_hole", 0.*mm , 7.*mm/2. , 210./2.*cm, 0.0, 360*degree);

G4VSolid* wall_Fire_brick_big_subtraction  = new G4SubtractionSolid ("wall_Fire_brick_big_subtraction",
                               wall_Fire_brick_big_box, wall_Fire_brick_hole, 0, G4ThreeVector(0.,0.,0.));
G4VSolid* wall_Fire_brick_small_subtraction  = new G4SubtractionSolid ("wall_Fire_brick_small_subtraction",
                               wall_Fire_brick_small_box, wall_Fire_brick_hole, 0, G4ThreeVector(0.,0.,0.));
G4VSolid* wall_concrete_big_subtraction  = new G4SubtractionSolid ("wall_concrete_big_subtraction",
                               wall_concrete_big_box, wall_concrete_hole, 0, G4ThreeVector(0.,0.,0.));
G4VSolid* wall_concrete_small_subtraction  = new G4SubtractionSolid ("wall_concrete_small_subtraction",
                               wall_concrete_small_box, wall_concrete_hole, 0, G4ThreeVector(0.,0.,0.));
G4LogicalVolume* wall_Fire_brick_big_log   = new G4LogicalVolume ( wall_Fire_brick_big_subtraction,   Fire_brick, "wall_Fire_brick_big_log");
G4LogicalVolume* wall_Fire_brick_small_log = new G4LogicalVolume ( wall_Fire_brick_small_subtraction, Fire_brick, "wall_Fire_brick_small_log");
G4LogicalVolume* wall_concrete_big_log     = new G4LogicalVolume ( wall_concrete_big_subtraction,     concrete,   "wall_concrete_big_log");
G4LogicalVolume* wall_concrete_small_log   = new G4LogicalVolume ( wall_concrete_small_subtraction,   concrete,   "wall_concrete_small_log");
wall_Fire_brick_big_log  ->SetVisAttributes(G4Colour(1.,1.,0.));
wall_Fire_brick_small_log->SetVisAttributes(G4Colour(1.,1.,0.));
wall_concrete_big_log    ->SetVisAttributes(G4Colour(0.,0.,1.));
wall_concrete_small_log  ->SetVisAttributes(G4Colour(0.,0.,1.));

/*new G4PVPlacement ( 0, G4ThreeVector(.0*cm, .0*cm, (258.7+150./2.)*cm-expHall_z),
                          wall_Fire_brick_big_log, "wall_Fire_brick_big", experimentalHall_log, false, 0);
new G4PVPlacement ( 0, G4ThreeVector(.0*cm, .0*cm, (968.35+100./2.)*cm-expHall_z),
                          wall_concrete_small_log, "wall_concrete_small", experimentalHall_log, false, 0);
new G4PVPlacement ( 0, G4ThreeVector(.0*cm, .0*cm, (1647.5+100./2.)*cm-expHall_z),
                          wall_Fire_brick_small_log, "wall_Fire_brick_small", experimentalHall_log, false, 0);*/
//new G4PVPlacement ( 0, G4ThreeVector(.0*cm, .0*cm, (2758.5+150./2.)*cm-expHall_z),
 //                         wall_Fire_brick_big_log, "wall_Fire_brick_big", experimentalHall_log, false, 0);
//============================== END OF ELI GEOMETRY ======================================================================


//============================== DETECTORS GEOMETRY =======================================================================

// _____________________________ DETECTOR POSITION __________________________________
    G4double x_cent, y_cent, z_cent;
    x_cent = 0. * cm;
    y_cent = 0. * cm;
    //z_cent = 3500.*cm+det_length/2.-expHall_z+1.*mm;
    z_cent = 4300.*cm-expHall_z+1.*mm;

// ______________________________ He3 DETECTORS _____________________________________
//Define sizes
G4double he3_sensitive_length, he3_front_dead_length, he3_rear_dead_length;
he3_sensitive_length = 45./2.* cm;
he3_front_dead_length = 2.5/2. * cm;
he3_rear_dead_length = 7.55/2. * cm;
G4double he3_length = he3_sensitive_length + he3_front_dead_length + he3_rear_dead_length;
G4double he3_radius = 1.2 * cm;
G4double he3_shell_length = he3_length; //carcasa det nu are capac
G4double he3_shell_irad   = he3_radius + 0.01 * mm;
G4double he3_shell_orad   = 1.25 * cm;
//Define positions of detectors
//HE3 MOTHER VOLUME
G4Tubs* he3_mother_vol = new G4Tubs("he3_mother_vol", 0.* cm , he3_shell_orad + 0.1*mm, he3_length + 0.1*mm, 0.0, 360*degree);
G4LogicalVolume* he3_mother_log = new G4LogicalVolume (he3_mother_vol, Vacuum, "he3_mother_log");
he3_mother_log->SetVisAttributes(G4VisAttributes::Invisible);
//HE3 Sensitive and NotSensitive Regions
G4Tubs* he3_sensitive_tube = new G4Tubs("he3_senzitive_tube", 0.* cm , he3_radius, he3_sensitive_length, 0.0, 360*degree);
G4Tubs* he3_rear_dead_tube = new G4Tubs("he3_rear_dead_tube", 0.* cm , he3_radius, he3_rear_dead_length, 0.0, 360*degree);
G4Tubs* he3_front_dead_tube = new G4Tubs("he3_front_dead_tube", 0.* cm , he3_radius, he3_front_dead_length, 0.0, 360*degree);
//HE3 INOX Shell
G4Tubs* he3_shell_tube = new G4Tubs("he3_shell_tube", he3_shell_irad, he3_shell_orad, he3_shell_length, 0.0, 360*degree);

G4Material* Sensitive_gas = He3_gas;  //Default He3_gas
     if (SD_material == "He3_gas") Sensitive_gas = He3_gas;
else if (SD_material == "BF3_gas") Sensitive_gas = BF3_gas;
G4cout<<"Counters Sensitive Gas Properties:"<<G4endl;
G4cout<<"Density: "<<SDdensity*cm3/mg<<" mg/cm3"<<G4endl;
G4cout<<"Temperature: "<< temp-STP_Temperature <<" Celsius"<<G4endl;
G4cout<<"Pressure: "<<pres/atmosphere<<" atm"<<G4endl;
G4cout<<Sensitive_gas<<G4endl;
G4LogicalVolume* he3_sensitive_log = new G4LogicalVolume ( he3_sensitive_tube, Sensitive_gas, "he3_sensitive_log");
he3_sensitive_log->SetVisAttributes(G4Colour(1.,0.,0.));
G4LogicalVolume* he3_rear_dead_log = new G4LogicalVolume ( he3_rear_dead_tube, Sensitive_gas, "he3_rear_dead_log");
he3_rear_dead_log->SetVisAttributes(G4Colour(1.,1.,1.));
G4LogicalVolume* he3_front_dead_log = new G4LogicalVolume ( he3_front_dead_tube, Sensitive_gas, "he3_front_dead_log");
he3_front_dead_log->SetVisAttributes(G4Colour(0.,0.,1.));
G4LogicalVolume* he3_shell_log = new G4LogicalVolume ( he3_shell_tube, Fe_Cr_Ni, "he3_shell_log");
he3_shell_log->SetVisAttributes(G4Colour(1.0,0.,0.));


//G4VPhysicalVolume* he3_front_dead = 
new G4PVPlacement ( 0,
 G4ThreeVector(0., 0., -he3_length+he3_front_dead_length), he3_front_dead_log,
      "he3_front_dead", he3_mother_log, false, 0);

//G4VPhysicalVolume* he3_sensitive = 
new G4PVPlacement ( 0,
 G4ThreeVector(0., 0., he3_front_dead_length-he3_rear_dead_length), he3_sensitive_log,
      "he3_sensitive", he3_mother_log, false, 0);

//G4VPhysicalVolume* he3_rear_dead = 
new G4PVPlacement ( 0,
 G4ThreeVector(0., 0., he3_front_dead_length+he3_sensitive_length), he3_rear_dead_log,
      "he3_rear_dead", he3_mother_log, false, 0);

//G4VPhysicalVolume* he3_shell = 
new G4PVPlacement ( 0,G4ThreeVector(0., 0., 0.),
              he3_shell_log, "he3_shell", he3_mother_log, false, 0);

//   G4VPhysicalVolume* he3_mother;


// ________________________________ MODERATOR _______________________________________
// __________________________________ SIZES _________________________________________
//Main body of moderator cube
G4double ModeratorLength = 50.*cm;
G4double ModeratorThick = 5.*cm;
G4double mod_body_x, mod_body_y, mod_body_z;
mod_body_x = ModeratorSize / 2.;
mod_body_y = ModeratorSize / 2.;
mod_body_z = ModeratorLength / 2.;
//Plates
G4double plate_thickness = ModeratorThick/2.;
//Front and read plates, square profile
G4double plate_front_side, plate_rear_side;
plate_front_side = (ModeratorSize + 0.1 * cm) / 2.;
plate_rear_side  = ModeratorSize / 2.;
//Lateral plates
G4double plate_lat_length, plate_lat_wide1, plate_lat_wide2;
plate_lat_length = (ModeratorLength + 2.*ModeratorThick + 0.1 * cm) / 2.;
plate_lat_wide1  = (ModeratorSize + 0.1 * cm) / 2.;
plate_lat_wide2  = (ModeratorSize + 2.*ModeratorThick + 0.1 * cm) / 2.;
//Cadmium plate thickness
G4double cd_thick = 0.05/2. * cm;
G4double cd_length = mod_body_z+plate_thickness+cd_thick;
//Gaps
G4double rad_gap_beam, rad_gap_detector, rad_gap_screw;
rad_gap_beam = 1.5 * cm;
rad_gap_detector = 1.35 * cm;
rad_gap_screw = 0.5 * cm;
//Screw
G4double rad_screw = 0.49 * cm;
//Spacing between the plates:
G4double plate_spacing = 0.1*mm;


// ___________________________ MODERATOR COMPONENTS __________________________________
//MOTHER VOLUME
G4Box* mother_vol = new G4Box("mother_vol", plate_lat_wide2+2*plate_spacing+2.*cd_thick+1.*cm,
                              plate_lat_wide2+2*plate_spacing+2.*cd_thick+1.*cm,
                              plate_lat_length+2*plate_spacing+2.*cd_thick+1.*cm);
//G4VSolid* mother_vol = new G4SubtractionSolid ("mother_vol", mother_vol_ini, target_tub, 0, G4ThreeVector(0.,0.,0.));
G4LogicalVolume* mother_log = new G4LogicalVolume (mother_vol, Vacuum, "mother_log");
fDetectorRegion   = new G4Region("Detector");
fDetectorRegion->AddRootLogicalVolume(mother_log);
mother_log->SetVisAttributes(G4VisAttributes::Invisible);

//CUBE AND PLATES
//Main body of moderator:
G4Box*  mod_cube = new G4Box ("mod_cube", mod_body_x,  mod_body_y,  mod_body_z);
//Rear plate (with gaps for detectors, beam and screws):
G4Box*  mod_plate_rear = new G4Box ("mod_plate_rear", plate_rear_side,  plate_rear_side,  plate_thickness);

//Front plate (with gap for beam):
G4Box*  mod_plate_front = new G4Box ("mod_plate_front", plate_front_side,  plate_front_side,  plate_thickness);

//Vertical plate:
G4Box*  mod_plate_lat1 = new G4Box ("mod_plate_lat1", plate_lat_wide1,  plate_thickness,  plate_lat_length);

//Horizontal plate:
G4Box*  mod_plate_lat2 = new G4Box ("mod_plate_lat2", plate_thickness, plate_lat_wide2,  plate_lat_length);

//Cadmium Rear plate (with gaps for detectors, beam and screws):
G4Box*  cd_plate_rear = new G4Box ("cd_plate_rear", plate_rear_side,  plate_rear_side,  cd_thick);

//Cadmium Front plate (with gap for beam):
G4Box*  cd_plate_front = new G4Box ("cd_plate_front", plate_front_side,  plate_front_side,  cd_thick);

//Cadmium Vertical plate:
G4Box*  cd_plate_lat1 = new G4Box ("cd_plate_lat1", cd_thick,  plate_rear_side,  cd_length);

//Cadmium Horizontal plate:
G4Box*  cd_plate_lat2 = new G4Box ("cd_plate_lat2", plate_front_side, cd_thick,  cd_length);

//Screws:
G4Tubs* screw_tube = new G4Tubs("screw_tube", 0.* cm , rad_screw, mod_body_z, 0.0, 360*degree);

//GAPS
G4Tubs* mod_beam_gap  = new G4Tubs("mod_beam_gap",  0.* cm , rad_gap_beam,     mod_body_z * 1.1, 0.0, 360*degree);
G4Tubs* mod_det_gap   = new G4Tubs("mod_det_gap",   0.* cm , rad_gap_detector, mod_body_z * 1.1, 0.0, 360*degree);
G4Tubs* mod_screw_gap = new G4Tubs("mod_screw_gap", 0.* cm , rad_gap_screw,    mod_body_z * 1.1, 0.0, 360*degree);

// ________________________ SUBTRACT FROM MODERATOR ____________________________________
// ____________BEAMP GAPS
G4VSolid* cube_subtraction           = new G4SubtractionSolid ("cube_subtraction",           mod_cube,        mod_beam_gap, 0, G4ThreeVector(0.,0.,0.));
G4VSolid* plate_rear_subtraction     = new G4SubtractionSolid ("plate_rear_subtraction",     mod_plate_rear,  mod_beam_gap, 0, G4ThreeVector(0.,0.,0.));
G4VSolid* plate_front_subtraction    = new G4SubtractionSolid ("plate_front_subtraction",    mod_plate_front, mod_beam_gap, 0, G4ThreeVector(0.,0.,0.));
G4VSolid* cd_plate_rear_subtraction  = new G4SubtractionSolid ("cd_plate_rear_subtraction",  cd_plate_rear,   mod_beam_gap, 0, G4ThreeVector(0.,0.,0.));
G4VSolid* cd_plate_front_subtraction = new G4SubtractionSolid ("cd_plate_front_subtraction", cd_plate_front,  mod_beam_gap, 0, G4ThreeVector(0.,0.,0.));
// ____________HE3 DETECTOR GAPS
G4ThreeVector uz, position;
G4double det_r, det_theta;
G4int iRing, iTheta, iThetaMax;
G4int iHeDet=0;
for (iRing = 1; iRing < (NumberOfRings+1); iRing ++)
  {
    iThetaMax = SetupGeometry->GetNumDetectors(iRing);
    for (iTheta = 1; iTheta < iThetaMax+1; iTheta ++)
      {
        SetupGeometry->GetDetPos(iRing, iTheta);
        det_r = SetupGeometry->r();
        det_theta = SetupGeometry->theta();
        G4cout<<"iRing = "<<iRing<<" iTheta = "<<iTheta;
        G4cout<<" r = "<<det_r<<" theta = "<<det_theta*180./pi<<G4endl;

        uz=G4ThreeVector(std::cos(det_theta), std::sin(det_theta), 0.);
        position = det_r * uz;
        cube_subtraction          = new G4SubtractionSolid ("cube_subtraction",          cube_subtraction,          mod_det_gap, 0, position);
        plate_rear_subtraction    = new G4SubtractionSolid ("plate_rear_subtraction",    plate_rear_subtraction,    mod_det_gap, 0, position);
        cd_plate_rear_subtraction = new G4SubtractionSolid ("cd_plate_rear_subtraction", cd_plate_rear_subtraction, mod_det_gap, 0, position);

        new G4PVPlacement ( 0, G4ThreeVector (det_r*std::cos(det_theta), det_r*std::sin(det_theta), he3_length-mod_body_z),
                   he3_mother_log, "he3_mother", mother_log, false, iHeDet);
        iHeDet++;
      }
  }
NumberOfDetectors = iHeDet;
// ____________SCREW GAPS
G4double DistScrewFromMargin, ScrewPosX, ScrewPosY;
DistScrewFromMargin = 3.*cm;
ScrewPosX = mod_body_x - DistScrewFromMargin;
ScrewPosY = mod_body_y - DistScrewFromMargin;
cube_subtraction = new G4SubtractionSolid ("cube_subtraction", cube_subtraction, mod_screw_gap, 0, G4ThreeVector( ScrewPosX, ScrewPosY,0.));
cube_subtraction = new G4SubtractionSolid ("cube_subtraction", cube_subtraction, mod_screw_gap, 0, G4ThreeVector(-ScrewPosX,-ScrewPosY,0.));
cube_subtraction = new G4SubtractionSolid ("cube_subtraction", cube_subtraction, mod_screw_gap, 0, G4ThreeVector( ScrewPosX,-ScrewPosY,0.));
cube_subtraction = new G4SubtractionSolid ("cube_subtraction", cube_subtraction, mod_screw_gap, 0, G4ThreeVector(-ScrewPosX, ScrewPosY,0.));

// ____________DECLARE LOGICAL VOLUMES + COLORS
//Moderator cube
G4LogicalVolume* cube_subtraction_log = new G4LogicalVolume ( cube_subtraction, polyethilene, "cube_subtraction_log");
cube_subtraction_log->SetVisAttributes(G4Colour(1.,1.,0.));

//Moderator Plates
G4LogicalVolume* plate_rear_subtraction_log = new G4LogicalVolume ( plate_rear_subtraction, polyethilene, "plate_rear_subtraction_log");
plate_rear_subtraction_log->SetVisAttributes(G4Colour(0.5,0.5,0.5));

G4LogicalVolume* plate_front_subtraction_log = new G4LogicalVolume ( plate_front_subtraction, polyethilene, "plate_front_subtraction_log");
plate_front_subtraction_log->SetVisAttributes(G4Colour(1.,0.5,0.));

G4LogicalVolume* mod_plate_lat1_log = new G4LogicalVolume ( mod_plate_lat1, polyethilene, "mod_plate_lat1_log");
mod_plate_lat1_log->SetVisAttributes(G4Colour(0.,0.,0.5));

G4LogicalVolume* mod_plate_lat2_log = new G4LogicalVolume ( mod_plate_lat2, polyethilene, "mod_plate_lat2_log");
mod_plate_lat2_log->SetVisAttributes(G4Colour(0.,0.5,0.));

//Cadmium Plates
G4LogicalVolume* cd_plate_rear_subtraction_log = new G4LogicalVolume ( cd_plate_rear_subtraction, Cadmium, "cd_plate_rear_subtraction_log");
cd_plate_rear_subtraction_log->SetVisAttributes(G4Colour(1.,0.,0.));

G4LogicalVolume* cd_plate_front_subtraction_log = new G4LogicalVolume ( cd_plate_front_subtraction, Cadmium, "cd_plate_front_subtraction_log");
cd_plate_front_subtraction_log->SetVisAttributes(G4Colour(1.,0.,0.));

G4LogicalVolume* cd_plate_lat1_log = new G4LogicalVolume ( cd_plate_lat1, Cadmium, "cd_plate_lat1_log");
cd_plate_lat1_log->SetVisAttributes(G4Colour(1.,0.,0.));

G4LogicalVolume* cd_plate_lat2_log = new G4LogicalVolume ( cd_plate_lat2, Cadmium, "cd_plate_lat2_log");
cd_plate_lat2_log->SetVisAttributes(G4Colour(1.,0.,0.));

G4LogicalVolume* screw_log = new G4LogicalVolume ( screw_tube, Fe_Cr_Ni, "screw_log");
screw_log->SetVisAttributes(G4Colour(1.,0.,0.));

// __________ PLACEMENTS _____________-

//G4VPhysicalVolume* cube = 
new G4PVPlacement ( 0, G4ThreeVector(0., 0., 0.), cube_subtraction_log, "cube", mother_log, false, 0);

//G4VPhysicalVolume* cd_plate_r = 
new G4PVPlacement ( 0, G4ThreeVector(0., 0., mod_body_z+plate_spacing+cd_thick),
        cd_plate_rear_subtraction_log,  "cd_plate_r", mother_log, false, 0);
//G4VPhysicalVolume* cd_plate_f = 
new G4PVPlacement ( 0, G4ThreeVector(0., 0.,-mod_body_z-plate_spacing-cd_thick),
        cd_plate_front_subtraction_log, "cd_plate_f", mother_log, false, 0);

//G4VPhysicalVolume* cd_plate_l1 = 
new G4PVPlacement ( 0, G4ThreeVector(-mod_body_y-plate_spacing-cd_thick, 0., plate_thickness),
                cd_plate_lat1_log, "cd_plate_l1", mother_log, false, 0);
                                 new G4PVPlacement ( 0, G4ThreeVector( mod_body_y+plate_spacing+cd_thick, 0., plate_thickness),
                cd_plate_lat1_log, "cd_plate_l1", mother_log, false, 1);
//G4VPhysicalVolume* cd_plate_l2 = 
new G4PVPlacement ( 0, G4ThreeVector(0., -mod_body_x-plate_spacing-cd_thick, plate_thickness),
                cd_plate_lat2_log, "cd_plate_l2", mother_log, false, 0);
                                 new G4PVPlacement ( 0, G4ThreeVector(0., +mod_body_x+plate_spacing+cd_thick, plate_thickness),
                cd_plate_lat2_log, "cd_plate_l2", mother_log, false, 1);

//G4VPhysicalVolume* plate_rear  = 
new G4PVPlacement ( 0, G4ThreeVector(0., 0.,  mod_body_z+2.*plate_spacing+2.*cd_thick+plate_thickness),
                   plate_rear_subtraction_log, "plate_rear", mother_log, false, 0);
//G4VPhysicalVolume* plate_front = 
new G4PVPlacement ( 0, G4ThreeVector(0., 0., -mod_body_z-2.*plate_spacing-2.*cd_thick-plate_thickness),
                   plate_front_subtraction_log, "plate_front", mother_log, false, 0);
//G4VPhysicalVolume* plate_lat1  = 
new G4PVPlacement ( 0, G4ThreeVector(0., -mod_body_y-2.*plate_spacing-2.*cd_thick-plate_thickness, 0.),
                   mod_plate_lat1_log, "plate_lat1", mother_log, false, 0);
                                 new G4PVPlacement ( 0, G4ThreeVector(0.,  mod_body_y+2.*plate_spacing+2.*cd_thick+plate_thickness, 0.),
                   mod_plate_lat1_log, "plate_lat1", mother_log, false, 1);
//G4VPhysicalVolume* plate_lat2 =  
new G4PVPlacement ( 0, G4ThreeVector(-mod_body_x-2.*plate_spacing-2.*cd_thick-plate_thickness, 0., 0.),
                   mod_plate_lat2_log, "plate_lat2", mother_log, false, 0); //AICI!!!
                                 new G4PVPlacement ( 0, G4ThreeVector( mod_body_x+2.*plate_spacing+2.*cd_thick+plate_thickness, 0., 0.),
                   mod_plate_lat2_log, "plate_lat2", mother_log, false, 1);

//G4VPhysicalVolume* screw = 
                           new G4PVPlacement ( 0, G4ThreeVector( ScrewPosX, ScrewPosY,0.), screw_log, "screw", mother_log, false, 0);
                           new G4PVPlacement ( 0, G4ThreeVector(-ScrewPosX, ScrewPosY,0.), screw_log, "screw", mother_log, false, 1);
                           new G4PVPlacement ( 0, G4ThreeVector(-ScrewPosX,-ScrewPosY,0.), screw_log, "screw", mother_log, false, 2);
                           new G4PVPlacement ( 0, G4ThreeVector( ScrewPosX,-ScrewPosY,0.), screw_log, "screw", mother_log, false, 3);

	//G4VPhysicalVolume* WTF = 
	new G4PVPlacement ( 0, G4ThreeVector(x_cent, y_cent, z_cent), mother_log, "mother", experimentalHall_log, false, 0);

//================================= DETECTORS GEOMETRY ======================================================================================


//================================= BEAM DUMP GEOMETRY ======================================================================================
// ____________________________________ SIZES ___________________________________________
G4double smallPb_radius = 10.*cm;
G4double bigPb_radius = 20.*cm;
G4double inPb_radius = (2.5/2.)*cm;
G4double bigParaffin_radius = 30.*cm;
G4double BD_PbShaftDepth = 60.*cm;
G4double BD_PbThickLength = 40.*cm;
G4double BD_PbThinLength = 20.*cm;
G4double BD_ParaffinOut = 20.*cm;
G4double BD_ParaffinLength = 79.3*cm;
G4double BD_CadmiumLength = 5*mm;
//____________DECLARE GEOMETRIC VOLUMES FOR BEAM DUMP
//Lead Beam Dump component
G4double zBeamDumpPbCoord[6]     = {0.,
                                    BD_PbShaftDepth,
                                    BD_PbShaftDepth,
                                    BD_PbShaftDepth+BD_PbThickLength,
                                    BD_PbShaftDepth+BD_PbThickLength,
                                    BD_PbShaftDepth+BD_PbThickLength+BD_PbThinLength};
G4double RinBeamDumpPbCoord [6]  = {inPb_radius, inPb_radius, 0., 0., 0., 0.};
G4double RoutBeamDumpPbCoord [6] = {bigPb_radius, bigPb_radius, bigPb_radius, bigPb_radius, smallPb_radius, smallPb_radius};
G4Polycone* BeamDump_PbCore = new G4Polycone("BeamDump_PbCore",0.*deg,360.*deg,6,zBeamDumpPbCoord,RinBeamDumpPbCoord,RoutBeamDumpPbCoord);
//Beam Dump Big Paraffin component
G4double zBeamDumpBigParaffinCoord[4]     = {0.,
                                             BD_ParaffinOut,
                                             BD_ParaffinOut,
                                             BD_ParaffinOut+BD_PbShaftDepth+BD_PbThickLength};
G4double RinBeamDumpBigParaffinCoord [4]  = {inPb_radius, inPb_radius, bigPb_radius+1.*mm, bigPb_radius+1.*mm};
G4double RoutBeamDumpBigParaffinCoord [4] = {bigParaffin_radius, bigParaffin_radius, bigParaffin_radius, bigParaffin_radius};
G4Polycone* BeamDump_BigParaffin = new G4Polycone("BeamDump_BigParaffin",0.*deg,360.*deg,4,
                                                  zBeamDumpBigParaffinCoord,RinBeamDumpBigParaffinCoord,RoutBeamDumpBigParaffinCoord);
//Beam Dump Paraffin tube component
G4Tubs* BeamDump_ParaffinTube = new G4Tubs("BeamDump_ParaffinTube",
                                              0.,
                                              smallPb_radius,
                                              BD_ParaffinLength/2.,
                                              0.0, 360*degree);
// ____________DECLARE LOGICAL VOLUMES + COLORS OF THE BEAM DUMP
//LOGICAL Lead Beam Dump component
G4LogicalVolume* BeamDump_PbCore_log = new G4LogicalVolume ( BeamDump_PbCore, Lead, "BeamDump_PbCore_log");
BeamDump_PbCore_log->SetVisAttributes(G4Colour(0.8,0.8,0.8));
//LOGICAL Beam Big Dump Paraffin component
G4LogicalVolume* BeamDump_BigParaffin_log = new G4LogicalVolume ( BeamDump_BigParaffin, paraffin, "BeamDump_BigParaffin_log");
BeamDump_BigParaffin_log->SetVisAttributes(G4Colour(1.,1.,0.));
//LOGICAL Beam Dump Paraffin tube component
G4LogicalVolume* BeamDump_ParaffinTube_log = new G4LogicalVolume ( BeamDump_ParaffinTube, paraffin, "BeamDump_ParaffinTube_log");
BeamDump_ParaffinTube_log->SetVisAttributes(G4Colour(1.,1.,0.));
//================================= BEAM DUMP GEOMETRY =========================================================================================


//================================= EXPERIMENTAL HALL GEOMETRY ==============================================================================
// ____________________________________ SIZES ___________________________________________
//Beam height
G4double beam_height = 1.5*m;
//Beam displacement
G4double beam_displacement = 1.5*m;
//Walls
G4double wall_thickness = 2. * m;
G4double sealing_thickness = 1.5 *m;
//Lateral walls
G4double lat_wall_length = 8. * m;
G4double lat_wall_height = 4. * m;
//Front and rear walls
G4double hall_wide = 8. * m;
G4double front_wall_width = hall_wide +2.*wall_thickness;
G4double front_wall_height =  lat_wall_height +2.*sealing_thickness;
//Floor and sealing
G4double floor_length  = lat_wall_length;
G4double floor_width = front_wall_width;
//Beam_hole
G4double rad_beam_hole = 6. * cm;
//G4double rad_beam_hole_vacuum = 2. *cm;
//Spacing between the walls
G4double wall_spacing = 0.1*mm;

//Placement of experimental hall
G4double Xhall_cent, Yhall_cent, Zhall_cent;
Xhall_cent = -beam_displacement;
Yhall_cent = lat_wall_height/2.-beam_height;
//Zhall_cent = 2758.5*cm+floor_width/2.-expHall_z; //z_cent;
Zhall_cent = 3470*cm+floor_width/2.-expHall_z; //z_cent;

//____________DECLARE GEOMETRIC VOLUMES FOR HALL WALLS
//MOTHER HALL VOLUME
G4Box* mother_hall_tot_vol = new G4Box("mother_hall_vol", floor_width/2. + 1.*cm,
                              front_wall_height/2. + 1.*cm,
                              floor_length/2. + wall_thickness + 1.*cm);
/*G4VSolid* mother_hall_vol1  = new G4SubtractionSolid ("detector_subtraction1", mother_hall_tot_vol,  mother_vol, 0,
                                                     G4ThreeVector(-Xhall_cent,-Yhall_cent,z_cent-Zhall_cent));*/
G4VSolid* mother_hall_vol  = new G4SubtractionSolid ("detector_subtraction1", mother_hall_tot_vol,  mother_vol, 0,
                                                     G4ThreeVector(-Xhall_cent,-Yhall_cent,z_cent-Zhall_cent));
/*G4Tubs* hall_beam_tube_vacuum = new G4Tubs("hall_beam_tube_vacuum", 0.* cm , rad_beam_hole_vacuum,
                        (lat_wall_length+wall_thickness-BD_ParaffinOut-2.*cm)/2., 0.0, 360*degree);
G4VSolid* mother_hall_vol  = new G4SubtractionSolid ("detector_subtraction", mother_hall_vol1,  hall_beam_tube_vacuum, 0,
                                                     G4ThreeVector(-Xhall_cent,-Yhall_cent,-(wall_thickness+BD_ParaffinOut)/2.));*/
G4LogicalVolume* mother_hall_log;
     if (E8_Filling=="Air")    mother_hall_log = new G4LogicalVolume (mother_hall_vol, Air,    "mother_hall_log");
else  mother_hall_log = new G4LogicalVolume (mother_hall_vol, Vacuum, "mother_hall_log");   // if (E8_Filling=="Vacuum")
mother_hall_log->SetVisAttributes(G4VisAttributes::Invisible);
//mother_hall_log->SetVisAttributes(G4VisAttributes(G4Colour(255.,0.,0.)));

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//G4LogicalVolume* dan_log = new G4LogicalVolume (hall_beam_tube_vacuum, Air, "dan_log");
//dan_log->SetVisAttributes(G4VisAttributes(G4Colour(255.,0.,0.)));
//new G4PVPlacement (0,G4ThreeVector(-Xhall_cent,-Yhall_cent,-(wall_thickness+BD_ParaffinOut)/2.), dan_log,"dan",mother_hall_log, false, 0);
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//Front wall (with hole for beam):
G4Box*  hall_wall_face = new G4Box ("hall_wall_face", front_wall_width/2.,  front_wall_height/2.,  wall_thickness/2.);
//Front beam hole
G4Tubs* hall_beam_tube = new G4Tubs("hall_beam_tube", 0.* cm , rad_beam_hole, wall_thickness/2.+1.*mm, 0.0, 360*degree);
G4VSolid* hall_wall_front  = new G4SubtractionSolid ("hall_wall_front", hall_wall_face, hall_beam_tube, 0,
                                                     G4ThreeVector(-Xhall_cent,-Yhall_cent,0.));
//Back beam dump
G4double zBeamDumpCoord[4]     = {-wall_thickness/2.-1.*mm, 0., 0., wall_thickness/2.+1.*mm};
G4double RinBeamDumpCoord [4]  = {0., 0., 0., 0.};
G4double RoutBeamDumpCoord [4] = {bigParaffin_radius+1.*mm, bigParaffin_radius+1.*mm, smallPb_radius+1.*mm, smallPb_radius+1.*mm};
G4Polycone* BeamDump_hole = new G4Polycone("BeamDump_hole",0.*deg,360.*deg,4,zBeamDumpCoord,RinBeamDumpCoord,RoutBeamDumpCoord);
G4VSolid* hall_wall_back  = new G4SubtractionSolid ("hall_wall_back", hall_wall_face, BeamDump_hole, 0,
                                                     G4ThreeVector(-Xhall_cent,-Yhall_cent, 0.));
//Lateral wall:
G4Box*  hall_wall_lat = new G4Box ("hall_wall_lat", wall_thickness/2., lat_wall_height/2., lat_wall_length/2. );
//Floor plate:
G4Box*  hall_plate_floor = new G4Box ("hall_plate_floor", floor_width/2.,  sealing_thickness/2., floor_length/2.);

// ____________DECLARE LOGICAL VOLUMES + COLORS OF EXPERIMENTAL HALL
//LOGICAL Front wall (with hole for beam):
G4LogicalVolume* hall_wall_front_log = new G4LogicalVolume ( hall_wall_front, concrete, "hall_wall_front_log");
hall_wall_front_log->SetVisAttributes(G4Colour(0.15,0.15,0.15));
G4LogicalVolume* hall_wall_back_log = new G4LogicalVolume ( hall_wall_back, concrete, "hall_wall_back_log");
hall_wall_back_log->SetVisAttributes(G4Colour(0.15,0.15,0.15));
//LOGICAL Lateral wall:
G4LogicalVolume* hall_wall_lat_log = new G4LogicalVolume ( hall_wall_lat, concrete, "hall_wall_lat_log");
hall_wall_lat_log->SetVisAttributes(G4Colour(0.3,0.3,0.3));
//LOGICAL Floor:
G4LogicalVolume* hall_plate_floor_log = new G4LogicalVolume ( hall_plate_floor, concrete, "hall_plate_floor_log");
hall_plate_floor_log->SetVisAttributes(G4Colour(0.2,0.2,0.2));
//hall_plate_floor_log->SetVisAttributes(G4VisAttributes::Invisible);

//______________HALL WALLS PLACEMENTS_______________
//Front & backward walls placement:
//G4VPhysicalVolume* hall_wall_front1 = 
new G4PVPlacement(0, G4ThreeVector( 0., 0., -(floor_length+wall_thickness)/2.-wall_spacing),
                                                        hall_wall_front_log, "hall_wall_front", mother_hall_log, false, 0);
//G4VPhysicalVolume* hall_wall_back1  = 
new G4PVPlacement(0, G4ThreeVector( 0., 0., +(floor_length+wall_thickness)/2.+wall_spacing),
                                                        hall_wall_back_log,  "hall_wall_back",  mother_hall_log, false, 0);
//Lateral walls placement:
//G4VPhysicalVolume* hall_wall_lat1 = 
new G4PVPlacement(0, G4ThreeVector( (front_wall_width-wall_thickness)/2.+wall_spacing, 0.,0.),
                                                        hall_wall_lat_log, "hall_wall_lat1", mother_hall_log, false, 0);
//G4VPhysicalVolume* hall_wall_lat2 = 
new G4PVPlacement(0, G4ThreeVector( -(front_wall_width-wall_thickness)/2.-wall_spacing, 0.,0.),
                                                        hall_wall_lat_log, "hall_wall_lat2", mother_hall_log, false, 1);
//Floor & sealing placement:
//G4VPhysicalVolume* hall_floor = 
new G4PVPlacement(0, G4ThreeVector( 0., (front_wall_height-sealing_thickness)/2.+wall_spacing,0.),
                                                        hall_plate_floor_log, "hall_plate_floor", mother_hall_log, false, 0);
//G4VPhysicalVolume* hall_sealing = 
new G4PVPlacement(0, G4ThreeVector( 0., -(front_wall_height-sealing_thickness)/2.-wall_spacing,0.),
                                                        hall_plate_floor_log, "hall_plate_sealing", mother_hall_log, false, 1);

/*
G4Orb* test_bila = new G4Orb ( "test_bila", 2.*cm);
G4LogicalVolume* test_bila_log = new G4LogicalVolume ( test_bila, concrete, "test_bila_log");
test_bila_log->SetVisAttributes(G4Colour(1.,0.,0.));
new G4PVPlacement ( 0, G4ThreeVector(0., 0., 2870*cm-expHall_z),
                           test_bila_log, "bila", experimentalHall_log, false, 0); */

/*----------------- TARGET PLACEMENT DIRECTLY INSIDE E8 EXPERIMENTAL AREA ----------------------------*/
//new G4PVPlacement (0,G4ThreeVector(-Xhall_cent,-Yhall_cent,3.*m), target_log,"target",mother_hall_log, false, 0);
/*----------------------------------------------------------------------------------------------------*/

//G4VPhysicalVolume* mother_hall = 
new G4PVPlacement ( 0, G4ThreeVector(Xhall_cent, Yhall_cent, Zhall_cent),
                                                     mother_hall_log, "mother_hall", experimentalHall_log, false, 0);

//=================================== EXPERIMENTAL HALL GEOMETRY ==========================================================================

//================================= BEAM DUMP PLACEMENTS ======================================================================================
//Lead Beam Dump component placement
//G4VPhysicalVolume* BeamDump_PbCore_Pys = 
new G4PVPlacement(0,
                                             G4ThreeVector( -Xhall_cent,-Yhall_cent, +(floor_length)/2.+wall_spacing),
                                             BeamDump_PbCore_log, "BeamDump_PbCore_Pys", mother_hall_log, false, 0);
//Beam Dump Big Paraffin component placement
//G4VPhysicalVolume* BeamDump_BigParaffin_Pys = 
new G4PVPlacement(0,
                                             G4ThreeVector( -Xhall_cent,
                                                            -Yhall_cent,
                                                            +(floor_length)/2.+wall_spacing-BD_ParaffinOut-1.*mm),
                                             BeamDump_BigParaffin_log, "BeamDump_BigParaffin_Pys", mother_hall_log, false, 0);
//Beam Dump Paraffin tube component placement
//G4VPhysicalVolume* BeamDump_ParaffinTube_Pys = 
new G4PVPlacement(0,
                                             G4ThreeVector( -Xhall_cent,
                                                            -Yhall_cent,
                                                            +(floor_length)/2.+wall_spacing+wall_thickness-BD_ParaffinLength/2.-BD_CadmiumLength-0.5*mm),
                                             BeamDump_ParaffinTube_log, "BeamDump_ParaffinTube_Pys", mother_hall_log, false, 0);
//================================= BEAM DUMP PLACEMENTS ======================================================================================


/*----------------- TARGET PLACEMENT INSIDE 4Pi DETECTOR ----------------------------*/
if (Target_placement=="detector") 
  {
  new G4PVPlacement (0,G4ThreeVector(0.,0.,Target_position), target_log,"target",mother_log, false, 0);
  TargetPlacement_position = z_cent + Target_position;
  if (Target_housing) 
  	new G4PVPlacement (0,G4ThreeVector(0.,0.,Target_position-(TGhousing_thick+TGhousing_space+target_z/2.)), TGhousing_log,"TargetHousing",mother_log, false, 0);
  }
/*----------------- TARGET PLACEMENT INSIDE E8 EXPERIMENTAL HALL ----------------------------*/
if (Target_placement=="E8") 
  {

/*  new G4PVPlacement (0,G4ThreeVector(-Xhall_cent,-Yhall_cent,Target_position-expHall_z-Zhall_cent), target_log,"target",mother_hall_log, false, 0);
  TargetPlacement_position = Target_position-expHall_z;
  if (Target_housing) 
  	new G4PVPlacement (0,G4ThreeVector(-Xhall_cent,-Yhall_cent,Target_position-expHall_z-Zhall_cent-(TGhousing_thick+TGhousing_space+target_z/2.)), TGhousing_log,"TargetHousing",mother_hall_log, false, 0);*/

  new G4PVPlacement (0,G4ThreeVector(-Xhall_cent,-Yhall_cent,Target_position + GammaSourcePos -Zhall_cent), target_log,"target",mother_hall_log, false, 0);
  TargetPlacement_position = Target_position + GammaSourcePos;
  if (Target_housing) 
  	new G4PVPlacement (0,G4ThreeVector(-Xhall_cent,-Yhall_cent,Target_position + GammaSourcePos -Zhall_cent-(TGhousing_thick+TGhousing_space+target_z/2.)), TGhousing_log,"TargetHousing",mother_hall_log, false, 0);

// ___________________________Brick_attenuator________________________________
  G4Tubs* Brick_attenuator_tub = new G4Tubs("Brick_attenuator_tub", Attenuator_hole/2. , 5.*cm , Attenuator_thickness/2., 0.0, 360*degree);
  G4LogicalVolume* Brick_attenuator_log = new G4LogicalVolume (Brick_attenuator_tub, AttenuatorMaterial->GetMaterial(Attenuator_material,true), "Brick_attenuator_log");
  if(Attenuator_placement)
    new G4PVPlacement ( 0, G4ThreeVector(-Xhall_cent, -Yhall_cent, Target_position + GammaSourcePos -Zhall_cent +Attenuator_position),
                          Brick_attenuator_log, "Brick_attenuator", mother_hall_log, false, 0);

  }
/* ---------------- PLACEMENT OF THE TARGET DIRECTLY INTO THE WORLD -----------------*/
else if (Target_placement=="world")
  {
  new G4PVPlacement (0,G4ThreeVector(0.,0.,Target_position-expHall_z), target_log,"target",experimentalHall_log, false, 0);
  TargetPlacement_position = Target_position-expHall_z;
  if (Target_housing) 
  	new G4PVPlacement (0,G4ThreeVector(0.,0.,Target_position-expHall_z-(TGhousing_thick+TGhousing_space+target_z/2.)), TGhousing_log,"TargetHousing",experimentalHall_log, false, 0);
  }
/*-----------------------------------------------------------------------------------*/



//__________________________________________________________________________________
  //------------------------------------------------ 
  // Sensitive detectors
  //------------------------------------------------ 

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String SDname = "eliLaBr_SD"; 
  eliLaBr_TrackerSD* aTrackerSD = new eliLaBr_TrackerSD( SDname );
  SDman->AddNewDetector( aTrackerSD ); 
  
  he3_sensitive_log -> SetSensitiveDetector( aTrackerSD );
  target_log -> SetSensitiveDetector( aTrackerSD );
  //dummy_log -> SetSensitiveDetector( aTrackerSD );
//  nai_log -> SetSensitiveDetector( aTrackerSD );
  
//  tg_log -> SetSensitiveDetector( aTrackerSD ); 
  
//__________________________________________________________________________________


  return experimentalHall_phys;
  
}

void eliLaBr_DetectorConstruction::SetE8_Filling(G4String val)      {E8_Filling = val;}
void eliLaBr_DetectorConstruction::SetTarget_thickness(G4double val){Target_thickness = val;}
void eliLaBr_DetectorConstruction::SetTarget_diameter(G4double val) {Target_diameter = val;}
void eliLaBr_DetectorConstruction::SetTarget_placement(G4String val){Target_placement = val;}
void eliLaBr_DetectorConstruction::SetTarget_housing(G4bool val)    {Target_housing = val;}
void eliLaBr_DetectorConstruction::SetTarget_position(G4double val) {Target_position = val;}
void eliLaBr_DetectorConstruction::SetTarget_material(G4String val) {Target_material = val;}
void eliLaBr_DetectorConstruction::SetAttenuator_thickness(G4double val){Attenuator_thickness = val;}
void eliLaBr_DetectorConstruction::SetAttenuator_hole(G4double val)     {Attenuator_hole = val;}
void eliLaBr_DetectorConstruction::SetAttenuator_position(G4double val) {Attenuator_position = val;}
void eliLaBr_DetectorConstruction::SetAttenuator_material(G4String val) {Attenuator_material = val;}
void eliLaBr_DetectorConstruction::SetAttenuator_placement(G4bool val)  {Attenuator_placement = val;}
void eliLaBr_DetectorConstruction::SetPlace_window(G4bool val)      {Place_window = val;}
void eliLaBr_DetectorConstruction::SetPlace_mirror(G4bool val)      {Place_mirror = val;}
void eliLaBr_DetectorConstruction::Set_BADVacPres(G4double val)     {BADVacPres = val;}
void eliLaBr_DetectorConstruction::SetSD_percentage(G4double val)   {SD_percentage = val;}
void eliLaBr_DetectorConstruction::SetSD_material(G4String val)     {SD_material = val;}
void eliLaBr_DetectorConstruction::SetSD_pressure(G4double val)     {pres = val;}
void eliLaBr_DetectorConstruction::SetModeratorSize(G4double val)   {ModeratorSize = val;}
void eliLaBr_DetectorConstruction::SetNumberOfRings(G4int val)      {NumberOfRings = val;}
G4int eliLaBr_DetectorConstruction::GetNumberOfDetectors()          {return NumberOfDetectors;}
