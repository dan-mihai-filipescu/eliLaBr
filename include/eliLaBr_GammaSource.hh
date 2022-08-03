/*
 * sbruGammaSource.hh
 *
 *  Created on: Mar 27, 2014
 *      Author: flippy
 */
//Class dealing with the interaction inside the gamma source.
//singleton
//holds all the information about the source description
//provides the Particle Gun with the photons
#ifndef ELILABRGAMMASOURCE_HH_
#define ELILABRGAMMASOURCE_HH_
#define VEC_DIM 2021
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4String.hh"
#include <vector>
#include "eliLaBr_Counter.hh"
#include "eliLaBr_Analysis.hh"
//#include "../include/sbruGlobalInfo.hh"
using namespace std;

class eliLaBr_DetectorConstruction;

class eliLaBr_GammaSource{

public:
    eliLaBr_GammaSource(eliLaBr_DetectorConstruction*);
    eliLaBr_GammaSource(G4double ene, G4double deene, G4double lwl);
    //void Set(double ene[10]);
    ~eliLaBr_GammaSource();
	void NextGamma();
    int HowMany();  //how many photons are to be shot by the source;

public:
 /* ==================== INPUT of GamaSource Class ==================== */
    G4double laser_power;  //Laser power in Watt
    G4double elec_intensity; //electron intensity nA
    G4double eene;   //MeV Electron beam kinetic energy in MeV
    G4double deene;  //Electron beam spread de/E in %
    G4int LaserType; // LaserType(G4int) can be 1(CO2), 2(INAZUMA), 3(TALON)
    G4double LaserToMidBL1, LaserToLens, LensToLWaist;
    G4double defaultLaserDiameter, LaserDiameter;
    G4double defaultSquareQF, SquareQF, defaultFocus, Focus, WaveLength;
    G4double RadiusBeforeFocus, RadiusAfterFocus;
    G4double RayleighBeforeFocus, RayleighAfterFocus;
    G4double lwl;    //laser beam wavelenght in nm;
    double mean;     //mean number of photons
    G4double xi, yi, zi, zmin, zmax;       //Absolute position of Electron beam waist and Gamma Source location along Z axis
    G4double sigx0, emmitx, sigy0, emmity; //Electron beam sizes at waist position and emmitances on transversal axes
    G4double laser_theta0, laser_phi0;     //Theta and Phi angles of incident laser beam
    G4double laser_phi0_ini;
    G4double ww, dispX, dispY, dispZ; //Laser beam waist and displacement of the focus
    G4double taui, taurot, tau;      //Laser polarization angle
    G4double Plin, Pcirc, PolLinValue;    //Probability of liniar and circular polarization of incident beam
    G4double PolDegree;
    G4int IDEAL_ELE_SW, IDEAL_LASER_SW, SOURCE_TYPE, Z_NORM_TYPE;
    G4int PLIN_MODEL_SW, VECTOR_MODEL_SW, STOKES_ELE_MODEL_SW, STOKES_LAB_MODEL_SW, STOKES_GRAPH_MODEL_SW;
    G4bool acceptZ, acceptXY, XfromFile, YfromFile;
    G4bool VectorModel, StokesEleModel, StokesLabModel, StokesLabGraph;
/* ===================================================================== */

 /* ==================== MAIN K-N generator METHODS ==================== */
 G4double PolAngleGenerator(G4double axis);
 void KNRand(G4double ep, G4double egmin, G4double egmax, G4double *eg, G4double *thrnd);
 //void Phi_gen(G4double t, G4double u, G4double tau, G4double Pt, G4double *Phi_ang );
 void Phi_gen(G4double t, G4double u, G4double Pt, G4double *Phi_ang );
 //G4ThreeVector SetNewPolarization(G4double epsilon, G4double ThetaKN, G4double phi, G4double *pBeta);
 void SetNewPolarization(G4bool VectorModel, G4double epsilon, G4double Pol, G4double ThetaKN, G4double phi, G4double *pBeta);
 /* ==================================================================== */

 /* === OUTPUT of GammaSource Class=== */
    G4double energy, ele_energy;
    G4ThreeVector direction, ele_direction, ex;
    G4ThreeVector position;
    G4ThreeVector Stokes, StokesLab;
    G4ThreeVector gammaPolarization1;
    G4double LCSweight;             //Weight of the generated event
    G4double norm;                  //Norm of the spectra
 /*====================================*/

private:
     eliLaBr_DetectorConstruction* Detector;
     //G4double test;
     G4double Z_position_read[VEC_DIM], betaX_read[VEC_DIM], betaY_read[VEC_DIM], alphaX_read[VEC_DIM], alphaY_read[VEC_DIM];
     G4double angleX_read[VEC_DIM], angleY_read[VEC_DIM], DispX[VEC_DIM], DispPX[VEC_DIM];
     G4double Z_position[VEC_DIM], betaX[VEC_DIM], betaY[VEC_DIM], alphaX[VEC_DIM], alphaY[VEC_DIM];
     G4double angleX[VEC_DIM], angleY[VEC_DIM], overlap_norm[VEC_DIM];
     G4double sigmaX2[VEC_DIM], sigmaY2[VEC_DIM], avgX[VEC_DIM], avgY[VEC_DIM];
     G4int index_laser[VEC_DIM];
     G4double *Z_laser, *laser_sig, *laser_sig2;
     G4double Z_max, Z_min, Z_pas, Z_pos, minnorm, maxnorm, sigma12, sigma22;
     G4double xint0, yint0, xang0, yang0, anglex, angley;
     G4double Total_integral, Partial_integral;
     G4int cnt, index_zmin, index_zmax, index_norm_min, index_norm_max, index_pas, index_laser_pas;
     G4String Obs_string;
     G4double eelec;                //Random energy of electron
     G4double ele_r2, ele_speed;    //Square electron radius and electron speed
     G4double Compton_CS;           //Total Comton CS
     G4double Gamma, Beta;          //Gamma and Beta values of electron
     G4double fgauge;               //Gauge constant
     G4double rnd, fact;
     G4double eph0, eph, deph;      //Laser photon energy
     G4double ele_theta0, ele_phi0; //Initial Theta and Phi angles of electron
     G4double theta, phi;           //Theta and Phi angles of the incident Laser photon in the Electron rest frame
     G4double x, y, t, u;           //Prerequisite variables for Phi generation from Klein-Nishina
     G4double egmin,egmax;          //Minimum and maximum energy of the scattered gamma rays in electron rest frame needed in K-N generator
     G4double egam, thrnd, Phi_ang; //Results of K-N generators (Energy, Theta, Phi)
     G4double Phi_ang_v;
     G4double CosTheta, CosThetaSq, SinThetaSq, SinPhi, SinPhiSq;
     G4double SinThSinPh2, PDterm1, PDterm2, PDterm3;
     G4LorentzVector eimp0;         //Lorentz vector of incident electron
     G4LorentzVector phimp0;        //Lorentz vector of incident Laser photon
     G4ThreeVector pol1dir;
     G4ThreeVector phot_dir0, phot_dir1;       //Direction of incident and scattered photon in the electron rest frame
     G4ThreeVector epsilon1, epsilon2, epsilonPerp, epsilonPar, versorPerp, versorPar;
     G4double Pol_eps1, Pol_eps2, termPerp, termPar, ratioPerpPar;
     G4ThreeVector eX, ey, ez, eZ, ez1, ez2;
     G4double slopeX, slopeY;
     G4double pBeta;
     G4double angle, Cosinus;
     G4LorentzVector Lepsilon1, Lepsilon2, LgammaPolarization1;
     G4LorentzVector pol0, pol1;    //Lorentz vector of incident Laser photon polarization, and after gauge transformation
     G4LorentzVector phimp1;        //Lorentz vector of outgoing Gamma photon
     G4LorentzVector eimp1;         //Lorentz vector of outgoing electron
     G4ThreeVector ele_dir0, laser_dir0, laser_dir; //3Vectors of incident Electrons, Laser photons (inside beam) and Direction of Laser beam
     G4double laser_theta, laser_phi;               //Theta and Phi angles of incident laser photon
     G4double xl,yl,zl,fl, laser_prob, l_prob, z_prob; //laser_strength   //position inside laser beam and probabilities
     //G4double sigx, sigy;                           //Electron beam size along Z axis
     G4double sw, betax0, betax, alphax, xang, betay0, betay, alphay, yang, kapa;  //Twiss parameters and angles of electron beam
     //G4double fx, fy, electron_prob, electron_strength;            //Electron probability
     G4double phix, ux, phiy, uy;                   //Variables needed to random generate electron beam space phase
     G4double xint, yint, zint;                     //position of interaction point inside electron beam
     G4double m2, wpi, wp, wp2, rayleigh;                //Laser parameters
     G4long Total_electron, Total_Zcycles, Total_XYcycles;           //Total number of cycles
     //G4double deltaZ_max;             //Maximum allowed space for interaction on Z axis
     G4double Stokes1, Stokes2, Stokes3; //Stokes parameters of incident laser photon
     G4double Stokes1f, Stokes2f, Stokes3f;  //Stokes parameters of outgoing gamma photon
     G4double Stokes1L, Stokes2L, Stokes3L;  //Stokes parameters of outgoing gamma photon in Lab. frame
     G4double F0, F0p, F1, F11, F22, F33, Fnum;         //Aditional Variables needed to compute Stokes parameters

     G4double dispXer, dispYer;
};


#endif /* ELILABRGAMMASOURCE_HH_ */
