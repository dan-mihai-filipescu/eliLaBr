#include "globals.hh"
#include "eliLaBr_PrimaryGeneratorAction.hh"
#include "eliLaBr_DetectorConstruction.hh"
#include "G4StokesVector.hh"
#include "G4ios.hh"
#include <iostream>
#include <fstream>
#include "Riostream.h"
#include "eliLaBr_GammaSource.hh"
#include "eliLaBr_Analysis.hh"
#include "TF1.h"
#include "TF2.h"

using namespace std;


  G4int n_particle = 1;
  G4double rand_prob;
  G4double prob_func;
  G4double gamma_energy_In, gamma_energy, neutron_energy;
  G4double rel_Zpos, emmiss_radius, emmiss_X, emmiss_Y, radius;
  G4ParticleGun *particleGun, *particleGun1;// = new G4ParticleGun(n_particle);
  G4ParticleTable* particleTable;// = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4double GammaProbability, NeutronProbability;
  G4int Pol_Model, PfromStokes, GS_izoAng, NS_izoAng;
  G4bool NS_spectrum;
  TF1* WeisskopfSpect = new TF1("Weisskopf","x*exp(-x/[0])/([0]*[0])",0.,20.);
  TF2* GammaAnizotropyE1 = new TF2("GammaAnizotropyE1","1.+((3.*x*x-1.)/2. - 3.*cos(2.*y)*(1.- x*x)/2.)/2.",-1.,1.,0.,2.*pi);
  TF2* GammaAnizotropyM1 = new TF2("GammaAnizotropyM1","1.+((3.*x*x-1.)/2. + 3.*cos(2.*y)*(1.- x*x)/2.)/2.",-1.,1.,0.,2.*pi);
  TF2* NeutronAnizotropyP = new TF2("NeutronAnizotropyP","(1.-x*x)*(1.+cos(2.*y))/2.",-1.,1.,0.,2.*pi);


  eliLaBr_PrimaryGeneratorAction::eliLaBr_PrimaryGeneratorAction(eliLaBr_DetectorConstruction* Det):Detector(Det)
{

  n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  particleGun1 = new G4ParticleGun(n_particle);
  particleTable = G4ParticleTable::GetParticleTable();
/*  	G4int seed=CLHEP::HepRandom::getTheSeed();
	cout<<"old seed="<<seed;
	seed=12345678;
	CLHEP::HepRandom::setTheSeed(seed);
	cout<<"                     new seed="<<seed<<endl;
*/

//_______________Gamma__Source___________

  eliLaBr_Source = new eliLaBr_GammaSource(Detector);
//  G4cout<<"\n\nElectron energy = "<<eliLaBr_Source->eene<<G4endl;

//  Probability = 1.;
  ifstream in;
  in.open("incident_energy.in",ios::in);
  in >> LCSource;
  in >> PolAngleIn;
  in >> Pol_Model >> PfromStokes;
  in >> GammaProbability;   if(GammaProbability<0.)   GammaProbability=0.;
  in >> NeutronProbability; if(NeutronProbability<0.) NeutronProbability=0.;  
  in >> gamma_energy_In;
  in >> neutron_energy; neutron_energy = neutron_energy * MeV;
  in >> GS_izoAng;
  in >> NS_izoAng;
  in >> NS_spectrum;
  in >> rel_Zpos; rel_Zpos = rel_Zpos * mm;
  in >> emmiss_radius; emmiss_radius = emmiss_radius * mm;

  if ((GammaProbability+NeutronProbability)>1.)
    {
    GammaProbability =   GammaProbability/(GammaProbability+NeutronProbability);
    NeutronProbability = NeutronProbability/(GammaProbability+NeutronProbability);
    }

  in.close();
  
  if(gamma_energy_In<0.)
    {
    char FileDan[100];     in.open("incident_fileName",ios::in);  in>>FileDan; in.close();        
    gamma_en = new TH1F( "gamma_en", "gamma_sp",5000,0,50.);
    //in.open("incident_EneSpec.in",ios::in);
    in.open(FileDan,ios::in);
    float tmp;
    for (int bb=0; bb<5000; bb++)
      {
      in >> tmp; 
      gamma_en->SetBinContent(bb,tmp);
      } 
    in.close();
    }
    else 
      gamma_energy = gamma_energy_In * MeV;
}

eliLaBr_PrimaryGeneratorAction::~eliLaBr_PrimaryGeneratorAction()
{

    delete eliLaBr_Source;
    delete particleGun;

}

void eliLaBr_PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

evNb = anEvent->GetEventID();

//{

//  particleGun->SetParticlePosition(G4ThreeVector( 0.0*cm, 1.*m, 1.*m ) ); // in front of 2R 
//  particleGun->SetParticlePosition(G4ThreeVector( 0.0*cm, -1.*m, .86*m ) ); // on the side of 2L 
  
//  particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="neutron"));
//  particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="gamma"));

//  double Eini=1.+(19. - 1.)*G4UniformRand();
//  double Eini=0.1 + (5.0 - 0.1)*G4UniformRand();
//  double Eini=0.7 + (19.0 - 0.7)*G4UniformRand();

//  particleGun->SetParticleEnergy(Eini*MeV); 

//  particleGun->SetParticleMomentumDirection( G4ThreeVector(0., 1., 0.) );	// towards 2R 
//  particleGun->SetParticleMomentumDirection( G4ThreeVector(0., -1., 0.) );	// towards 2L 

//}
/*
  angle  = 0.; 
  energy = 0.; 
	// n yield + g yield = 1.81 + 1.13 = 2.94 (Shen... d+Cu @ 33MeV )
  if ( gRandom->Rndm() * 2.94 > 1.81 )
	{ particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="gamma"));
	  th2gamma->GetRandom2(angle, energy); } 
  else	
	{ particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="neutron"));
	  th2neutron->GetRandom2(angle, energy); } 

// take care of the angle... 
  theta = angle; 
  phi = 2*pi*G4UniformRand(); 
  direction[0] = sin(theta)*cos(phi);
  direction[1] = sin(theta)*sin(phi);
  direction[2] = cos(theta);

  particleGun->SetParticleEnergy(energy*MeV); 
  particleGun->SetParticlePosition(G4ThreeVector( 0.0*cm, 0.*cm, 0.*cm ) ); 
  particleGun->SetParticleMomentumDirection( direction );
*/

//{
//  particleGun->SetParticlePosition(G4ThreeVector( 0.0*cm, 0.*cm, 0.*cm ) ); // 86.*cm

//  particleGun->SetParticlePosition(G4ThreeVector( (221.-150.)*mm, 221./2.*mm, 221./2.*mm ) ); // 15 cm in front of cluster 1 
  
//  particleGun->SetParticleMomentumDirection( G4ThreeVector() = G4RandomDirection() ); 
//  particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="gamma"));

//  particleGun->SetParticleEnergy(22.650*MeV); 	// 3.2 width state in 12C 
//  particleGun->GeneratePrimaryVertex(anEvent); 

/*  particleGun->SetParticleEnergy(18.212*MeV); 
  particleGun->GeneratePrimaryVertex(anEvent); 
  particleGun->SetParticleEnergy(4.438*MeV); 
  particleGun->GeneratePrimaryVertex(anEvent); */ 

//  th2neutron->GetRandom2( dummy, energy); 
//  energy = gamma_en->GetRandom();
//}

//  G4cout<<"\n\nElectron energy = "<<eliLaBr_Source->eene<<G4endl;

//G4bool LCSource=true;
G4double PolAngle,PolX,PolY,PolSelector;

//Default polarization
PolAngle = 0.;

//__________________________________GAMMA SOURCE_______________________________________
G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
if(LCSource) {

  eliLaBr_Source->NextGamma();

  if(eliLaBr_Source->SOURCE_TYPE!=0)
    {
    if((eliLaBr_Source->StokesEleModel || eliLaBr_Source->StokesLabModel) && (eliLaBr_Source->PLIN_MODEL_SW!=0)) PolSelector = 0.;
      else 
      PolSelector = G4UniformRand();
    if(!Pol_Model)
      {
      if(!PfromStokes)
        {
        PolVector = eliLaBr_Source->gammaPolarization1;
        if(PolSelector<eliLaBr_Source->PolDegree) PolAngle = PolVector.getPhi();
          else
          {
          PolAngle = 2.*pi * G4UniformRand() - pi;
          PolX = cos(PolAngle);
          PolY = sin(PolAngle);
          PolVector.set(PolX,PolY,0.);
          }
        }
        else
        {
        //PolVector = eliLaBr_Source->StokesLab;
        PolVector = eliLaBr_Source->Stokes;
        PolX = PolVector.getX();
        PolY = PolVector.getY();
        //if (PolSelector<(PolX*PolX+PolY*PolY))
        if (PolSelector<eliLaBr_Source->PolDegree)
          {
          if (PolX != 0.) PolAngle = atan(-PolY/PolX) / 2.;
          else if (PolY != 0.)
            {
            if(PolY > 0.)PolAngle = pi/4.;
            if(PolY < 0.) PolAngle = -pi/4.;
            } else PolAngle = pi * G4UniformRand();
//          if(sign(cos(2.*PolAngle)) != sign(PolX)) PolAngle = PolAngle + pi/2.;
	  if(sign(cos(2.*PolAngle)) != sign(PolX))
	    {
	    if(PolAngle>0.) PolAngle = PolAngle - pi/2.;
	      else PolAngle = PolAngle + pi/2.;
	    }
          //if(PolAngle<0.) PolAngle = PolAngle + pi;
	  //if(0.5<G4UniformRand())
	  //  {
	  //  if(PolAngle>0.) PolAngle = pi - PolAngle;
	  //    else PolAngle = -pi - PolAngle;
	  //  }
          }
          else PolAngle = pi * G4UniformRand();
        PolVector = (eliLaBr_Source->ex).rotate(eliLaBr_Source->direction,-PolAngle);
        PolAngle = PolVector.getPhi();
        //PolX = cos(PolAngle);
        //PolY = sin(PolAngle);
        //PolVector.set(PolX,PolY,0.);
        }
      }
      else
      {
      //PolVector = eliLaBr_Source->StokesLab;
      PolVector = eliLaBr_Source->Stokes;
      PolX = PolVector.getX();
      PolY = PolVector.getY();
      //if (PolSelector<(PolX*PolX+PolY*PolY))
      if (PolSelector<eliLaBr_Source->PolDegree)
        {
        if (PolX != 0.) PolAngle = atan(-PolY/PolX) / 2.;
        else if (PolY != 0.)
          {
          if(PolY > 0.)PolAngle = pi/4.;
          if(PolY < 0.) PolAngle = -pi/4.;
          } else PolAngle = pi * G4UniformRand();
//        if(sign(cos(2.*PolAngle)) != sign(PolX)) PolAngle = PolAngle + pi/2.;
	  if(sign(cos(2.*PolAngle)) != sign(PolX))
	    {
	    if(PolAngle>0.) PolAngle = PolAngle - pi/2.;
	      else PolAngle = PolAngle + pi/2.;
	    }
        //if(PolAngle<0.) PolAngle = PolAngle + pi;
        }
        else PolAngle = pi * G4UniformRand();
      PolVector = (eliLaBr_Source->ex).rotate(eliLaBr_Source->direction,-PolAngle);
      PolAngle = PolVector.getPhi();
      PolVector = eliLaBr_Source->Stokes;
      //PolVector = eliLaBr_Source->StokesLab;
      }

    analysisManager->FillH2(14, eliLaBr_Source->energy, PolAngle, 1.);

    particleGun->SetParticlePosition(eliLaBr_Source->position);
    particleGun->SetParticleMomentumDirection(eliLaBr_Source->direction);
    //particleGun->SetParticlePolarization( eliLaBr_Source->Stokes);
    //particleGun->SetParticlePolarization( G4ThreeVector(PolX,PolY,0.) );
    particleGun->SetParticlePolarization( PolVector );
    particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="gamma")); 
    particleGun->SetParticleEnergy( eliLaBr_Source->energy * MeV );
    }

  if(eliLaBr_Source->SOURCE_TYPE>=0)
    {
    particleGun1->SetParticlePosition(eliLaBr_Source->position);
    particleGun1->SetParticleMomentumDirection(eliLaBr_Source->ele_direction);
    particleGun1->SetParticleDefinition(particleTable->FindParticle(particleName="e-")); 
    particleGun1->SetParticleEnergy( eliLaBr_Source->ele_energy * MeV );
    }

  //particleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
  //particleGun->SetParticlePolarization( G4StokesVector(G4ThreeVector(1., 0., 0.)));
  //particleGun->SetParticleEnergy( 1.*MeV );
//  particleGun->SetParticleEnergy( energy*MeV );

//  G4cout<<"Gamma energy = "<<eliLaBr_Source->energy<<G4endl;

  if(eliLaBr_Source->SOURCE_TYPE!=0) particleGun->GeneratePrimaryVertex(anEvent);
  if(eliLaBr_Source->SOURCE_TYPE>=0) particleGun1->GeneratePrimaryVertex(anEvent);
  ev_Weight = eliLaBr_Source->LCSweight;
  norm = eliLaBr_Source->norm;
  n_type = 0;
}
//______________________________________________________________________________________

//_______________________________SIMPLE GAMMA SOURCE____________________________________
else {
  ev_Weight = 1.;
  norm = 1.;

rand_prob = G4UniformRand();
if(PolAngleIn>=0.) PolAngle = PolAngleIn*deg;
  else PolAngle = pi*G4UniformRand();
if (rand_prob < (1.-(GammaProbability+NeutronProbability))) 
  {
  if(!Pol_Model)
    {
    PolX = cos(PolAngle);
    PolY = sin(PolAngle);
    }
    else
    {
    PolX = cos(2.*PolAngle);
    PolY = -sin(2.*PolAngle);
    }

  particleGun->SetParticlePosition(G4ThreeVector( 0.0*cm, 0.*cm, -Detector->expHall_z+0.5*m) );      //3500.*cm-2200.*cm+1.*mm-15.*m) );
  particleGun->SetParticleMomentumDirection( G4ThreeVector(0., 0., 1.) );
  particleGun->SetParticlePolarization( G4ThreeVector(PolX,PolY,0.));
  particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="gamma"));
  if(gamma_energy_In<0.) gamma_energy = gamma_en->GetRandom() * MeV;
  particleGun->SetParticleEnergy( gamma_energy);
  //  particleGun->SetParticleEnergy( energy*MeV );
  particleGun->GeneratePrimaryVertex(anEvent);
  n_type = 0;
  }

 else if (rand_prob < (1.-NeutronProbability))
  {
  if(!emmiss_radius)
    {
    emmiss_X = 0.*cm;
    emmiss_Y = 0.*cm;
    }
    else
    {
    radius = emmiss_radius * sqrt(G4UniformRand());
    phi = 2.*pi*G4UniformRand();
    emmiss_X = radius*cos(phi);
    emmiss_Y = radius*cos(phi);
    }
  particleGun->SetParticlePosition(G4ThreeVector( emmiss_X, emmiss_Y, Detector->TargetPlacement_position + rel_Zpos) );
//  phi = 2*pi*G4UniformRand();
//  theta = acos(2.*G4UniformRand()-1.);
  if (GS_izoAng == 0) 
    {
    phi = 2.*pi*G4UniformRand();
    Cos_th = 2.*G4UniformRand()-1.;
    }
  else if (GS_izoAng == 1) GammaAnizotropyE1->GetRandom2(Cos_th,phi);
  else if (GS_izoAng == 2) GammaAnizotropyM1->GetRandom2(Cos_th,phi);
  theta=acos(Cos_th);
  direction[0] = sin(theta)*cos(phi);
  direction[1] = sin(theta)*sin(phi);
  direction[2] = cos(theta);
  if(PolAngle>0.) {direction.rotateZ(PolAngle); phi = phi+PolAngle; if(phi>(2.*pi)) phi = phi - 2.*pi;}
  particleGun->SetParticleMomentumDirection( direction );     
  particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="gamma"));
  if(gamma_energy_In<0.) gamma_energy = -gamma_energy_In * MeV;
  particleGun->SetParticleEnergy( gamma_energy );
  particleGun->GeneratePrimaryVertex(anEvent);
  //if ((GammaProbability+NeutronProbability)>=1.) 
  n_type = 1;
  //  else n_type = 0;
  }

 else
  {
  //  cout<<"Phys Event!!! "<< evNb<<"; Prob: "<<rand_prob<<G4endl;
  if(!emmiss_radius)
    {
    emmiss_X = 0.*cm;
    emmiss_Y = 0.*cm;
    }
    else
    {
    radius = emmiss_radius * sqrt(G4UniformRand());
    phi = 2.*pi*G4UniformRand();
    emmiss_X = radius*cos(phi);
    emmiss_Y = radius*cos(phi);
    }
  particleGun->SetParticlePosition(G4ThreeVector( emmiss_X, emmiss_Y, Detector->TargetPlacement_position + rel_Zpos) );
  //particleGun1->SetParticlePosition(G4ThreeVector( 0.0*cm, 0.*cm, Detector->TargetPlacement_position ) );
  //particleGun->SetParticlePosition(G4ThreeVector( 0.0*cm, 0.*cm, 3500.*cm-2200.*cm+1.*mm ) );
  /*do {phi = 2*pi*G4UniformRand();
      theta = acos(2.*G4UniformRand()-1.);
      prob_func = sin(theta)*sin(theta)*(1.+cos(2.*phi))/2.;
      if (NS_izoAng) rand_prob = G4UniformRand();
        else rand_prob = 0.;
      }
      while ( rand_prob >= prob_func );*/
  if (NS_izoAng == 0) 
    {
    phi = 2.*pi*G4UniformRand();
    Cos_th = 2.*G4UniformRand()-1.;
    }
  else if (NS_izoAng == 1) NeutronAnizotropyP->GetRandom2(Cos_th,phi);
  theta=acos(Cos_th);
  //analysisManager->FillH2(15, theta, phi, 1.);
  direction[0] = sin(theta)*cos(phi);
  direction[1] = sin(theta)*sin(phi);
  direction[2] = cos(theta);
  //direction = G4RandomDirection();
  if(PolAngle>0.) {direction.rotateZ(PolAngle); phi = phi+PolAngle; if(phi>(2.*pi)) phi = phi - 2.*pi;}
  analysisManager->FillH2(15, Cos_th, phi, 1.);
  particleGun->SetParticleMomentumDirection( direction );
  /*do {phi = 2*pi*G4UniformRand();
      theta = acos(2.*G4UniformRand()-1.);
      prob_func = sin(theta)*sin(theta)*(1.+cos(2.*phi))/2.;
      if (NS_izoAng) rand_prob = G4UniformRand();
        else rand_prob = 0.;
      }
      while ( rand_prob >= prob_func );*/
  /*if (NS_izoAng == 0) 
    {
    phi = 2.*pi*G4UniformRand();
    Cos_th = 2.*G4UniformRand()-1.;
    }
  else if (NS_izoAng == 1) NeutronAnizotropyP->GetRandom2(Cos_th,phi);
  theta=acos(Cos_th);*/
  //analysisManager->FillH2(15, theta, phi, 1.);
  /*direction[0] = sin(theta)*cos(phi);
  direction[1] = sin(theta)*sin(phi);
  direction[2] = cos(theta);*/
  //direction = G4RandomDirection();
  /*if(PolAngle>0.) {direction.rotateZ(PolAngle); phi = phi+PolAngle; if(phi>(2.*pi)) phi = phi - 2.*pi;}
  analysisManager->FillH2(15, Cos_th, phi, 1.);
  particleGun1->SetParticleMomentumDirection( direction );*/
  //particleGun->SetParticleMomentumDirection( G4ThreeVector() = G4RandomDirection() );
  particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="neutron"));
  //particleGun1->SetParticleDefinition(particleTable->FindParticle(particleName="neutron"));
  if (NS_spectrum) {
    //WeisskopfSpect->SetParameters(neutron_energy);
    WeisskopfSpect->SetParameter(0,neutron_energy);
    particleGun->SetParticleEnergy( WeisskopfSpect->GetRandom());
    //particleGun1->SetParticleEnergy( WeisskopfSpect->GetRandom());
    }
   else 
    {
    particleGun->SetParticleEnergy( neutron_energy);
    //particleGun1->SetParticleEnergy( neutron_energy);
    }
  particleGun->GeneratePrimaryVertex(anEvent);
  //particleGun1->GeneratePrimaryVertex(anEvent);
  //if ((GammaProbability+NeutronProbability)>=1.) 
  n_type = 1;
  //  else n_type = 0;
  }
  
}
//______________________________________________________________________________________

NDet = Detector->GetNumberOfDetectors();
eliLaBr_EventInformation* anInfo = new eliLaBr_EventInformation();
if(eliLaBr_Source->StokesLabGraph)
  anInfo->Set_PrimaryPartPolarization(eliLaBr_Source->StokesLab);
  else
  anInfo->Set_PrimaryPartPolarization(eliLaBr_Source->Stokes);
anInfo->Set_weight(ev_Weight);
anInfo->Set_norm(norm);
anInfo->Set_angle(PolAngle);
anInfo->Set_src_n_type(n_type);
anInfo->Set_N_Det(NDet);
G4Event* otherEvent = (G4Event*)anEvent;
otherEvent->SetUserInformation(anInfo);

}

G4int eliLaBr_PrimaryGeneratorAction::sign(G4double x)
{
    if(x<0.) return -1;
    if(x>0.) return 1;
    return 0;
}
