/*
 * eliLaBrGammaSource.cc
 *
 *  Created on: Mar 27, 2014
 *      Author: flippy
 */
#include "eliLaBr_GammaSource.hh"
#include "eliLaBr_DetectorConstruction.hh"
#include "Randomize.hh"

eliLaBr_Counter electron_counter, cycleZ_counter, cycleXY_counter;

eliLaBr_GammaSource::eliLaBr_GammaSource(eliLaBr_DetectorConstruction* Det):Detector(Det) {


    Z_pas = 10.*mm;
    Z_min = -10100.*mm;
    Z_max = 10100.*mm;
    Z_laser = new G4double[2*VEC_DIM];
    laser_sig = new G4double[2*VEC_DIM];
    laser_sig2 = new G4double[2*VEC_DIM];
    XfromFile = false;
    YfromFile = false;

    egam = 0.;
    thrnd = 0.;
    Phi_ang = 0.;
    eimp0 = G4LorentzVector();
    phimp0 = G4LorentzVector();
    phimp1 = G4LorentzVector();

    //deltaZ_max = 20. * m;
    //ww = 0.55* mm;      // Laser beam waist
    //m2 = 1.19;           // Laser beam quality factor, adimensional
    //deph = 0.05;
    //disp = 1830. * mm; // position of Laser beam waist (displacement from Electron waist)
//    disp = 0.;
//    zi = - 49.99 * m;
//    zi = -Detector->expHall_z+5.*cm;  // position of Electron beam waist
    zi = Detector->GammaSourcePos;

    ifstream in;
    in.open("incident_gamma.in",ios::in);

    in >> laser_power;      laser_power = laser_power * watt;
    in >> elec_intensity;   elec_intensity = elec_intensity * 1.e-3 *ampere;
    in >> eene;             eene = eene * MeV;
    in >> deene;            deene = deene*perCent;
    eelec = eene + electron_mass_c2;
    Gamma = eelec / electron_mass_c2;
    Beta = sqrt(1. - 1./(Gamma*Gamma));
    ele_speed = Beta * c_light;
    G4cout<<"Electron speed: "<<ele_speed*s/m<<" m/s"<<G4endl;

    in >> LaserType;  // LaserType(G4int) can be 1(CO2), 2(INAZUMA), 3(TALON)
    if(LaserType!=1)
      if(LaserType!=2)
	if(LaserType!=3)
          if(LaserType!=4) LaserType =2;  //The default is INAZUMA Laser

    if(LaserType==1)   // CO2 LASER parameters
      {
      LaserToMidBL1 = 16259.;      //mm.
      LaserToLens = 6536.;         //mm.
      defaultLaserDiameter = 2.4;  //mm.
      defaultSquareQF = 1.2;       // NO Dim!
      defaultFocus = 3250.;               //mm.
      WaveLength = 10500.;         //nm.
      }

    if(LaserType==2)   // INAZUMA LASER parameters
      {
      LaserToMidBL1 = 14982.;      //mm.
      LaserToLens = 8220.;         //mm.
      defaultLaserDiameter = 0.6;  //mm.
      defaultSquareQF = 1.2;       // NO Dim!
      //2022_01_04 changed FOCUS to 99% to reroduce dist=1830 mm
      defaultFocus = 5064.; //5115.09;             //mm.
      WaveLength = 1064.;          //nm.
      }

    if(LaserType==3)   // TALON LASER parameters
      {
      LaserToMidBL1 = 15259.;      //mm.
      LaserToLens = 6887.;         //mm.
      defaultLaserDiameter = 1.048;//mm.
      defaultSquareQF = 1.02;      // NO Dim!
      //2022_01_04 changed FOCUS to 99% similar to INAZUMA case (same lense is used for both lasers)
      defaultFocus = 4942.05; //4991.97;             //mm.
      WaveLength = 532.;           //nm.
      }

    if(LaserType==4)   // ELI-NP conditions (no lenses)
      {
      LaserToMidBL1 = 0.;          //mm.
      LaserToLens = 0.;            //mm.
      defaultLaserDiameter = 0.050;//mm.
      defaultSquareQF = 1.2;       // NO Dim!
      defaultFocus = 0.;           //mm.
      WaveLength = 532.;           //nm.
      //WaveLength = 0.05;
      }

    // Laser diameter at beam waist (from specifications)
    in >> LaserDiameter;  if(LaserDiameter<=0.) LaserDiameter = defaultLaserDiameter;
    // Laser beam quality factor, adimensional
    in >> SquareQF;  if(SquareQF<=0.) SquareQF = defaultSquareQF;
    // Laser focusing lens in mm
    in >> Focus;  if(Focus<=0.) Focus = defaultFocus;

    LaserToMidBL1 = LaserToMidBL1 * mm;
    LaserToLens = LaserToLens * mm;
    LaserDiameter = LaserDiameter * mm;
    Focus = Focus * mm;
    WaveLength = WaveLength * nm;
    lwl = WaveLength;

   RadiusBeforeFocus = LaserDiameter / 2.;
   RayleighBeforeFocus = pi * RadiusBeforeFocus * RadiusBeforeFocus / (SquareQF * WaveLength);
   if(LaserType!=4)
     {
// Gaussian optics computation for LASER focusing by lens
      LensToLWaist = 1. / (1./Focus - 
        1./(LaserToLens + RayleighBeforeFocus*RayleighBeforeFocus/(LaserToLens-Focus)));
      RayleighAfterFocus = sqrt((1. / (1./Focus - 1./LaserToLens) - LensToLWaist) * 
        (LensToLWaist - Focus));
      RadiusAfterFocus = sqrt(RayleighAfterFocus * SquareQF * WaveLength / pi);
// position of Laser beam waist (displacement from Electron waist)
      dispZ = LaserToMidBL1 - LensToLWaist;
      }
      else
      {
      RayleighAfterFocus = RayleighBeforeFocus;
      RadiusAfterFocus = RadiusBeforeFocus;
      dispZ = 0.*mm;
      }
    rayleigh = RayleighAfterFocus;

// Print LASER Parameters
    G4cout<<" === LASER PARAMETERS ==="<<G4endl;
    if(LaserType==1) G4cout<<"Laser Type: CO2"<<G4endl;
    if(LaserType==2) G4cout<<"Laser Type: INAZUMA"<<G4endl;
    if(LaserType==3) G4cout<<"Laser Type: TALON"<<G4endl;
    if(LaserType==4) G4cout<<"Laser Type: ELI-NP"<<G4endl;
    G4cout<<"Laser wavelength = "<<WaveLength/nm<<" nm."<<G4endl;
    G4cout<<"Laser quality factor ^2 = "<<SquareQF<<G4endl;
    G4cout<<"Laser diameter = "<<LaserDiameter/mm<<" mm."<<G4endl;
    if(LaserType!=4)
      {
      G4cout<<"Lens focal length = "<<Focus/mm<<" mm."<<G4endl;
      G4cout<<"Distance from Laser output to Electron waist = "<<LaserToMidBL1/mm<<" mm."<<G4endl;
      G4cout<<"Laser --> Lens Distance = "<<LaserToLens/mm<<" mm."<<G4endl;
      G4cout<<"Laser radius at waist before lens = "<<RadiusBeforeFocus/mm<<" mm."<<G4endl;
      G4cout<<"Raylegh before lens = "<<RayleighBeforeFocus/mm<<" mm."<<G4endl;
      G4cout<<"Lens --> Waist Distance = "<<LensToLWaist/mm<<" mm."<<G4endl;
      }
    G4cout<<"Laser radius at waist after lens = "<<RadiusAfterFocus/mm<<" mm."<<G4endl;
    G4cout<<"Raylegh after lens = "<<RayleighAfterFocus/mm<<" mm."<<G4endl;
    G4cout<<"Laser displacement = "<<dispZ/mm<<" mm."<<G4endl;

    ww = RadiusAfterFocus / 2.;  // Laser sigma at beam waist (according to ISO D=4sigma definition)
    
// OLD input code; Not needed any more. (08-Nov-2021)
/*  in >> ww; ww = ww * mm;
    in >> m2;
    in >> lwl; lwl = lwl * nm; */

    in >> deph;
    in >> dispX;            dispX = dispX * mm;
    in >> dispY;            dispY = dispY * mm;
// OLD input code; Not needed any more. (08-Nov-2021)
// There is NO need to input dispZ - it is computed in the Gaussian optics section
//  in >> dispZ;            dispZ = dispZ * mm;
    in >> xi;               xi = xi * mm;
    in >> yi;               yi = yi * mm;
    in >> zmin;             zmin = zmin * m;
    in >> zmax;             zmax = zmax * m;
    in >> sigx0;
    in >> emmitx;
    if(sigx0>0.) sigx0 = sigx0 * mm;
    if(emmitx>0.)
      {
      emmitx = emmitx * mm; 
      if(LaserType==4) emmitx /= (Beta*Gamma);
      }
    if((sigx0<=0.)||(emmitx<=0.)) XfromFile = true;
      else betax0 = sigx0 * sigx0 / emmitx;
    in >> sigy0;
    in >> emmity;
    if(sigy0>0.) sigy0 = sigy0 * mm;
    if(emmity>0.)
      {
      emmity = emmity * mm;
      if(LaserType==4) emmity /= (Beta*Gamma);
      }
    if((sigy0<=0.)||(emmity<=0.)) YfromFile = true;
      else betay0 = sigy0 * sigy0 / emmity;
    in >> laser_theta0;     laser_theta0 = laser_theta0 * deg;
    in >> laser_phi0_ini;    //laser_phi0 = laser_phi0 * deg;
    in >> taui;              taui = taui * deg;
    in >> Plin;
    in >> Pcirc;
    in >> PLIN_MODEL_SW >> VECTOR_MODEL_SW >> STOKES_ELE_MODEL_SW >> STOKES_LAB_MODEL_SW >> STOKES_GRAPH_MODEL_SW;
    if(VECTOR_MODEL_SW) VectorModel=true;
      else VectorModel=false;
    if(STOKES_ELE_MODEL_SW) StokesEleModel=true;
      else StokesEleModel=false;
    if(STOKES_LAB_MODEL_SW) StokesLabModel=true;
      else StokesLabModel=false;
    if(STOKES_GRAPH_MODEL_SW) StokesLabGraph=true;
      else StokesLabGraph=false;
    G4cout << "PlinModel = "<< PLIN_MODEL_SW << "; StokesEleModel = " << StokesEleModel << "; StokesLabModel = " << StokesLabModel <<G4endl;
    in >> IDEAL_ELE_SW;
    in >> IDEAL_LASER_SW;
    in >> SOURCE_TYPE;      //  -1 Just photons
                            //   0 Just electrons
                            //  +1 photons & electrons
    in >> Z_NORM_TYPE;      //   0 Accept any uniform random Z and multiply the matrices with Electron-Laser overlap integral
                            //   1 Accept Z values proportional with Electron-Laser overlap integral
    in.close();

    if (Plin>1.) Plin = 1.;
    //if ((Plin+Pcirc)>1.) Pcirc = 1. - Plin;
    if ((Plin*Plin+Pcirc*Pcirc)>1.) Pcirc = sqrt(1. - Plin*Plin);
    PolLinValue = Plin;
    
    if(PLIN_MODEL_SW>0) 
      {
      PolLinValue=1.;
      //Plin = 1.;
      }

// OLD code; Not needed any more. (08-Nov-2021)
// There is NO need to compute rayleigh here - it is computed in the Gaussian optics section
//  rayleigh = pi * ww * ww / (m2 * lwl);

 /*   sigx0 = 0.306 * mm;
    sigx = 0.985 * mm;
    sigy0 = 0.180 * mm;
    sigy = 0.240 * mm;
    sw = 7119. * mm;
    betax0 = sigx0*sw/sqrt(sigx*sigx-sigx0*sigx0);
    emmitx = sigx0*sigx0/betax0;
    betay0 = sigy0*sw/sqrt(sigy*sigy-sigy0*sigy0);
    emmity = sigy0*sigy0/betay0;*/

    G4cout<<"Electron radius = "<<classic_electr_radius/fermi<<" fermi"<<G4endl;
    ele_r2 = classic_electr_radius * classic_electr_radius;
    G4cout<<"Electron C.S. = "<<pi*ele_r2/barn<<" barn"<<G4endl;
    //Compute photon energy
    eph0 = hbarc * 2. * pi / lwl; G4cout<<"\n =============>>  ENERGIA FOTONULUI: "<<eph0/eV<<" eV"<<G4endl;

    //eph0 = 0.5*MeV;

    x = 2. * Gamma*eph0 * (1.+Beta) / electron_mass_c2;
    Compton_CS = 2.*pi*ele_r2 * ( (1.-4./x-8./(x*x))*log(1.+x) + 1./2. + 8./x - 1./(2.*(1.+x)*(1.+x)) ) / x;
    fact = 1. + Beta * cos(laser_theta0);

    char *line;
    line = (char *) calloc(160,sizeof(char));
    FILE *inputFile;
    inputFile = fopen("NewSUBARU_optics_BL01.txt","r");
    if(fgets(line,150,inputFile)==NULL) printf("\nERROR while reading Twiss parameters file\n");
    cnt = 0;
    while(fgets(line,150,inputFile))
      {
      //printf("%s\n",line);
      sscanf(line,"%lf%lf%lf%lf%lf%lf%lf%lf%lf", &Z_position_read[cnt],&betaX_read[cnt], &betaY_read[cnt], &DispX[cnt], &DispPX[cnt], &alphaX_read[cnt], &alphaY_read[cnt], &angleX_read[cnt], &angleY_read[cnt]);
      Z_position_read[cnt] *= meter;
      betaX_read[cnt] *= meter;
      betaY_read[cnt] *= meter;
      DispX[cnt] *= meter;
      //Multiply both  phases with Horizontal Betatron tune (6.30)
      //although vertical Betatron tune = 2.23
      angleX_read[cnt] *= 6.30;
      angleY_read[cnt] *= 6.30;
      cnt++;
      }
    fclose(inputFile);
    free(line);
    if(cnt!=VEC_DIM) G4cout<<"!!! WARNING !!! Dimension of Twiss parameters file is different from VEC_DIM"<<G4endl;
      else G4cout<<"Twiss parameters file has good size"<<G4endl;
    if(XfromFile)
      {
      betax0 = betaX_read[(VEC_DIM-1)/2];
      if(sigx0<=0.) sigx0 = sqrt(betax0*emmitx);
      if(emmitx<=0.) emmitx = sigx0*sigx0/betax0;
      }
    if(YfromFile)
      {
      betay0 = betaY_read[(VEC_DIM-1)/2];
      if(sigy0<=0.) sigy0 = sqrt(betay0*emmity);
      if(emmity<=0.) emmity = sigy0*sigy0/betay0;
      }

    G4cout<<"SigX0 = "<<sigx0/mm<<"mm; EmmitX = "<<emmitx/mm<<" mm"<<G4endl;
    G4cout<<"SigY0 = "<<sigy0/mm<<"mm; EmmitY = "<<emmity/mm<<" mm"<<G4endl;

    FILE *outputFile;
    outputFile = fopen("LCS_overlap.out","w");

    fprintf(outputFile,"SigX0 = %f mm; EmmitX = %e mm*rad;\n",sigx0/mm,emmitx/mm);
    fprintf(outputFile,"SigY0 = %f mm; EmmitY = %e mm*rad;\n",sigy0/mm,emmity/mm);

    zmin = floor((zmin-Z_min)/Z_pas)*Z_pas+Z_min;
    if(zmin<Z_min) zmin=Z_min;
    index_zmin = (G4int)round((zmin-Z_min)/Z_pas);
    zmax = ceil((zmax-Z_min)/Z_pas)*Z_pas+Z_min;
    if(zmax>Z_max) zmax=Z_max;
    index_zmax = (G4int)round((zmax-Z_min)/Z_pas);
    fprintf(outputFile,"Zmin = %f mm; Index_Zmin = %i;\n",zmin/mm,index_zmin);
    fprintf(outputFile,"Zmax = %f mm; Index_Zmax = %i;\n",zmax/mm,index_zmax);
    for(cnt=0;cnt<VEC_DIM;cnt++)
      {
      Z_pos = Z_min + Z_pas*cnt;
      Z_position[cnt] = Z_pos;
      if(XfromFile)
        {
        betaX[cnt] = betaX_read[cnt];
        alphaX[cnt] = alphaX_read[cnt];
        angleX[cnt] = angleX_read[cnt];
        }
        else
        {
        betaX[cnt] = betax0+Z_pos*Z_pos / betax0;
        alphaX[cnt] = -Z_pos / betax0;
        angleX[cnt] = atan(Z_pos/betax0);
        }
      if(YfromFile)
        {
        betaY[cnt] = betaY_read[cnt];
        alphaY[cnt] = alphaY_read[cnt];
        angleY[cnt] = angleY_read[cnt];
        }
        else
        {
        betaY[cnt] = betay0+Z_pos*Z_pos / betay0;
        alphaY[cnt] = -Z_pos / betay0;
        angleY[cnt] = atan(Z_pos/betay0);
        }
      }
    fprintf(outputFile,"\nComputed Laser sigma:\n");
    fprintf(outputFile,"Index  Z_pos (m)    Sigma (mm)\n");
    for(cnt=0;cnt<2*VEC_DIM;cnt++)
      {
      Z_pos = 2.*Z_min + Z_pas*cnt;
      Z_laser[cnt] = Z_pos;
      wpi = Z_pos / rayleigh;
      laser_sig[cnt] = ww * sqrt ( 1. + wpi*wpi ) ;
      laser_sig2[cnt] = ww * ww * ( 1. + wpi*wpi ) ;
      fprintf(outputFile,"%5i %9.3f %15.9f\n",cnt,Z_laser[cnt]/m, laser_sig[cnt]/mm);
      }
    //fprintf(outputFile,"\nElectron-Laser Overlap:\n");
    //fprintf(outputFile,"Index  Z_pos (m)    Sigma (mm)\n");
    //In case of ideal LASER compute the direction (slopes) of the laser
    //01.Dec.2021
    //Not only for IDEAL LASER. Compute laser direction in general
//    if(!IDEAL_LASER_SW)
//      {
      //Define both LASER versor ends
      ez1.set(0.,0.,0.);
      ez2.set(0.,0.,1.);
      //Perform Translation and Rotation transformations on Laser versor
      ez1.rotateY(laser_theta0);
      ez2.rotateY(laser_theta0);
      ez1.rotateZ(laser_phi0);
      ez2.rotateZ(laser_phi0);
      ez1.set(ez1.getX()+dispX, ez1.getY()+dispY, ez1.getZ()+dispZ);
      ez2.set(ez2.getX()+dispX, ez2.getY()+dispY, ez2.getZ()+dispZ);
      slopeX = (ez2.getX()-ez1.getX()) / (ez2.getZ()-ez1.getZ());
      slopeY = (ez2.getY()-ez1.getY()) / (ez2.getZ()-ez1.getZ());
//      }
    maxnorm = 0;
    Total_integral = 0.;
    Partial_integral = 0.;
    if((!IDEAL_ELE_SW)&&(!IDEAL_LASER_SW))
      {
      Total_integral = (Z_max-Z_min)/Compton_CS;
      Partial_integral = (zmax-zmin)/Compton_CS;
      minnorm = 1.;
      maxnorm = 1.;
      }
      else
      {
      for(cnt=0;cnt<VEC_DIM;cnt++)
        {
        Z_pos = Z_position[cnt];
        betax = betaX[cnt];
        betay = betaY[cnt];
        //01.Dec.2021
        //New algorithm of finding laser axis coordinates as a function of Z
        xl = ez1.getX() + slopeX*(Z_pos-ez1.getZ());
        yl = ez1.getY() + slopeY*(Z_pos-ez1.getZ());
	//2022_01_01 Correct the Z_pos shift 
        //zl = sqrt(xl*xl+yl*yl+Z_pos*Z_pos);
	zl = sqrt(xl*xl+yl*yl+(Z_pos-ez1.getZ())*(Z_pos-ez1.getZ()));
        //***** End of algorithm for laser axis coordinates calculation ****
        //Comment the OLD algorithm
        /*xl = - dispX;
        yl = - dispY;
        zl = Z_pos - dispZ;
//Align laser along z axis
        laser_dir.set(xl,yl,zl);
        laser_dir.rotateZ(-laser_phi0);
        laser_dir.rotateY(-laser_theta0);
        xl = laser_dir.getX();
        yl = laser_dir.getY();
        zl = laser_dir.getZ();*/
        index_laser_pas = (G4int)round((zl - 2.*Z_min)/Z_pas);
        index_laser[cnt] = index_laser_pas;
        wp2 = laser_sig2[index_laser_pas];
        if(!IDEAL_ELE_SW) overlap_norm[cnt] = exp(-xl*xl/wp2/2.)*exp(-yl*yl/wp2/2.)/(sqrt(wp2)*sqrt(wp2));
          else
          {
          sigma12 = betax*emmitx;
          sigma22 = betay*emmity;
          sigmaX2[cnt] = sigma12*wp2/(sigma12+wp2);
          sigmaY2[cnt] = sigma22*wp2/(sigma22+wp2);
          avgX[cnt] = xl*sigma12/(sigma12+wp2);
          avgY[cnt] = yl*sigma22/(sigma22+wp2);
          if(!IDEAL_LASER_SW) wp2 = 0.;
          overlap_norm[cnt] = exp(-xl*xl/(sigma12+wp2)/2.)*exp(-yl*yl/(sigma22+wp2)/2.)/(sqrt(sigma12+wp2)*sqrt(sigma22+wp2));
          }
        if((cnt==0) || (cnt==(VEC_DIM-1))) Total_integral += overlap_norm[cnt];
          else Total_integral += 2.*overlap_norm[cnt];
        //fprintf(outputFile,"%5i %9.3f %5i %9.3f %15.9f, %f\n",cnt,Z_position[cnt]/m, index_laser_pas, Z_laser[index_laser_pas]/m, laser_sig[index_laser_pas]/mm, overlap_norm[cnt]);
        if((index_zmin<=cnt)&&(cnt<=index_zmax))
          {
          if((cnt==index_zmin)||(cnt==index_zmax)) Partial_integral += overlap_norm[cnt];
            else Partial_integral += 2.*overlap_norm[cnt];
          if(cnt == index_zmin) {minnorm = overlap_norm[cnt]; index_norm_min = cnt;}
            else if(minnorm>overlap_norm[cnt]) {minnorm = overlap_norm[cnt]; index_norm_min = cnt;}
          if(maxnorm<overlap_norm[cnt]) {maxnorm = overlap_norm[cnt]; index_norm_max = cnt;}
          }
        }
      Total_integral *= Z_pas/(4.*pi);
      Partial_integral *= Z_pas/(4.*pi);
      }
    fprintf(outputFile,"\nIndex_MinNorm = %d; MinNorm = %e;\n", index_norm_min, minnorm);
    fprintf(outputFile,"Index_MaxNorm = %d; MaxNorm = %e;\n", index_norm_max, maxnorm);
    G4cout<<"Total_integral = "<<Total_integral*mm<<" 1/mm"<<G4endl;
    G4cout<<"Partial_integral = "<<Partial_integral*mm<<" 1/mm"<<G4endl;
    minnorm = minnorm/maxnorm;
    fprintf(outputFile,"\n                                   Computed Twisss parameters                                                     LASER Sigma                Electron-LASER Overlap\n");
    fprintf(outputFile,"Index  Z_pos (m)     BetaX (m)       BetaY (m)        AlphaX         AlphaY       PhiX (rad)     PhiY (rad)   LIndex L_pos (m)    LSigma (mm)        Norm            Obs.\n");
    for(cnt=0;cnt<VEC_DIM;cnt++)
      {
      overlap_norm[cnt] = overlap_norm[cnt]/maxnorm;
      index_laser_pas = index_laser[cnt];
      Obs_string = "";
      if((index_zmin<=cnt)&&(cnt<=index_zmax)) Obs_string += " *";
      if(cnt==index_zmin) Obs_string += " Z_min";
      if(cnt==index_zmax) Obs_string += " Z_max";
      if(cnt==index_norm_min) Obs_string += " N_min";
      if(cnt==index_norm_max) Obs_string += " N_max";
      fprintf(outputFile,"%5i %9.3f %15.9f %15.9f %14.9f %14.9f %14.9f %14.9f   %5i %9.3f %15.9f  %16.9e %s\n",cnt,Z_position[cnt]/m, betaX[cnt]/m, betaY[cnt]/m, alphaX[cnt], alphaY[cnt], angleX[cnt], angleY[cnt], index_laser_pas, Z_laser[index_laser_pas]/m, laser_sig[index_laser_pas]/mm, overlap_norm[cnt], Obs_string.data());
      }
    maxnorm = maxnorm/maxnorm;

    G4cout<<"Number of photons = "<<laser_power/eph0*second<<" ph/s"<<G4endl;
    G4cout<<"Number of electrons = "<<elec_intensity/eplus*second<<" e/s"<<G4endl;
    fprintf(outputFile,"\nNumber of photons = %14.6e ph/s\n",laser_power/eph0*second);
    fprintf(outputFile,"Number of electrons = %14.6e e/s\n",elec_intensity/eplus*second);
    G4cout<<"Total Compton CS: "<<Compton_CS/barn<<" barn"<<G4endl;
    fprintf(outputFile,"Total Compton CS: = %12.6f barn\n",Compton_CS/barn);
    fprintf(outputFile,"Total integral = %12.6f (1/mm)\n",Total_integral*mm);
    fprintf(outputFile,"Partial integral = %12.6f (1/mm)\n",Partial_integral*mm);
    G4cout<<"Fact = "<<fact<<G4endl;
    fprintf(outputFile,"Relativistic factor (1+Beta*Cos(theta)) = %12.6f\n",fact);
    fprintf(outputFile,"Total luminosity (per W & mA) = %14.6e (barn*s)^-1\n",fact*(watt/eph0)*(1.e-3*ampere/eplus)*Total_integral*barn*second/ele_speed);
    fprintf(outputFile,"Partial luminosity (per W & mA) = %14.6e (barn*s)^-1\n",fact*(watt/eph0)*(1.e-3*ampere/eplus)*Partial_integral*barn*second/ele_speed);
    fprintf(outputFile,"Total gamma rate = %14.6e gammas/s\n",fact*(laser_power/eph0)*(elec_intensity/eplus)*Total_integral*Compton_CS*second/ele_speed);
    //norm = fact*(laser_power/eph0)*(elec_intensity/e_SI)*Compton_CS*(zmax - zmin)/ele_speed;
    norm = fact*(laser_power/eph0)*(elec_intensity/eplus)*Compton_CS*Partial_integral*second/ele_speed;
    if(norm>(laser_power/eph0*second)) norm = laser_power/eph0*second;
    if(norm>(elec_intensity/eplus*second)) norm = elec_intensity/eplus*second;
    fprintf(outputFile,"Partial gamma rate (used for normalisation) = %14.6e gammas/s\n",norm);
    G4cout<<"NORM = "<<norm<<" gammas/s"<<G4endl;
    fclose(outputFile);

}

eliLaBr_GammaSource::eliLaBr_GammaSource(G4double ene, G4double dene, G4double plwl) {
        eene = ene; //MeV electron beam in MeV
        deene = dene; // beam spread de/E in %
        lwl = plwl; //laser beam wavelength in nm;
}


eliLaBr_GammaSource::~eliLaBr_GammaSource() {
    G4cout<<"\nLCS GammaSource:"<<G4endl;
    G4cout<<"Total number of Z shots: "<<Total_Zcycles<<G4endl;
    G4cout<<"Total number of XY shots: "<<Total_XYcycles<<G4endl;
    G4cout<<"Total number of generated electrons: "<<Total_electron<<G4endl;
}


void eliLaBr_GammaSource::NextGamma() {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

/*___________________Generate Laser rotation angle inside recirculator___*/
    if (laser_phi0_ini < 0.) laser_phi0 = floor(32.*G4UniformRand())*(360./32.)*deg;
      else laser_phi0 = laser_phi0_ini*deg;

/*___________________Generate Interaction Point__________________________*/
   // do
   // {
      electron_counter.IncrementCounter();
    do {
      zint = (zmax - zmin) * G4UniformRand() + zmin;
      if(zint<Z_min) index_pas = 0;
        else if(zint>Z_max) index_pas = VEC_DIM -1;
          else index_pas = (G4int)round((zint - Z_min)/Z_pas);
      //z_prob = (maxnorm-minnorm) * G4UniformRand() + minnorm;
      z_prob = maxnorm * G4UniformRand();
      cycleZ_counter.IncrementCounter();
      if((!Z_NORM_TYPE)||((!IDEAL_ELE_SW)&&(!IDEAL_LASER_SW))) acceptZ = false;
        else acceptZ = (z_prob>overlap_norm[index_pas]);
      } while (acceptZ);
    do {
      if(!IDEAL_ELE_SW)
        {
        xint=yint=0.; //=zint
        xang=yang=0.;
        }
        else
        {
        betax = betaX[index_pas];
        betay = betaY[index_pas];
        alphax = alphaX[index_pas];
        alphay = alphaY[index_pas];
        anglex = angleX[index_pas];
        angley = angleY[index_pas];
//Compute beta Twiss parameters along the line
        //betax = betax0+zint*zint / betax0;
        //betay = betay0+zint*zint / betay0;
//Generate phase space of electrons acording to C. Sun
        ux = G4RandExponential::shoot(1.);
        phix = 2. * pi * G4UniformRand();
        uy = G4RandExponential::shoot(1.);
        phiy = 2. * pi * G4UniformRand();
        if(!IDEAL_LASER_SW)
          {
          /*xl = - dispX;
          yl = - dispY;
          zl = zint - dispZ;
          //Align laser along z axis
          laser_dir.set(xl,yl,zl);
          //29.Nov.2021 Changed the order of rotations
          laser_dir.rotateY(-laser_theta0);
          laser_dir.rotateZ(-laser_phi0);
          xl = laser_dir.getX();
          yl = laser_dir.getY();
          zl = laser_dir.getZ();
          xint = xl;
          yint = yl;*/
          //Find the interaction position knowing Laser slope
          xint = ez1.getX() + slopeX*(zint-ez1.getZ());
          yint = ez1.getY() + slopeY*(zint-ez1.getZ());
          //Generate electron position in the focus
          xint0 = sqrt( 2. * ux * emmitx * betax0 ) * cos(phix);
          yint0 = sqrt( 2. * uy * emmity * betay0 ) * cos(phiy);
          //Solve the transfer equations in order to find electron momentum in the focus (knowing xint & yint)
          xang0 = xint - sqrt(betax/betax0)*cos(anglex)*xint0;
          if(XfromFile) xang0 -= DispX[index_pas]*deene;
          xang0 /= (sqrt(betax*betax0)*sin(anglex));
          yang0 = yint - sqrt(betay/betay0)*cos(angley)*yint0;
          yang0 /= (sqrt(betay*betay0)*sin(angley));
          //Apply transfer equations in order to get electron momentum in the interaction point
          xang = -(alphax*cos(anglex)+sin(anglex))*xint0/sqrt(betax*betax0) + (sqrt(betax0/betax)*(cos(anglex)-alphax*sin(anglex)))*xang0;
          if(XfromFile) xang += DispPX[index_pas]*deene;
          yang = -(alphay*cos(angley)+sin(angley))*yint0/sqrt(betay*betay0) + (sqrt(betay0/betay)*(cos(angley)-alphay*sin(angley)))*yang0;
          //On IDEAL LASER position of interaction is defined => deduce the angles, they are not random generated
          //xang = -sqrt(2.*ux*emmitx/betax)*(alphax*cos(phix)+sin(phix));
          //yang = -sqrt(2.*uy*emmity/betay)*(alphay*cos(phiy)+sin(phiy));
          }
          else
          {
          xint0 = sqrt( 2. * ux * emmitx * betax0 ) * cos(phix);
          xang0 = -sqrt(2.*ux*emmitx/betax0)*sin(phix);
          xint = sqrt(betax/betax0)*cos(anglex)*xint0 + sqrt(betax*betax0)*sin(anglex)*xang0;
          xang = -(alphax*cos(anglex)+sin(anglex))*xint0/sqrt(betax*betax0) + (sqrt(betax0/betax)*(cos(anglex)-alphax*sin(anglex)))*xang0;
          if(XfromFile)
            {
            xint += DispX[index_pas]*deene;
            xang += DispPX[index_pas]*deene;
            }
          yint0 = sqrt( 2. * uy * emmity * betay0 ) * cos(phiy);
          yang0 = -sqrt(2.*uy*emmity/betay0)*sin(phiy);
          yint = sqrt(betay/betay0)*cos(angley)*yint0 + sqrt(betay*betay0)*sin(angley)*yang0;
          yang = -(alphay*cos(angley)+sin(angley))*yint0/sqrt(betay*betay0) + (sqrt(betay0/betay)*(cos(angley)-alphay*sin(angley)))*yang0;
          }
        }

      if((!IDEAL_ELE_SW)||(!IDEAL_LASER_SW))
        {
        laser_prob = 1.;
        acceptXY = false;
        }
        else
        {
// Consider Z displacement of laser
        //dispXer = -1.5*mm + 3. * mm * G4UniformRand();
        dispXer = 0.;
        xl = xint - dispX - dispXer;
        yl = yint - dispY;
        zl = zint - dispZ;
//Align laser along z axis
        laser_dir.set(xl,yl,zl);
        laser_dir.rotateZ(-laser_phi0);
        laser_dir.rotateY(-laser_theta0);
        xl = laser_dir.getX();
        yl = laser_dir.getY();
        zl = laser_dir.getZ();
//Compute phasespace of laser
        index_laser_pas = (G4int)round((zl - 2.*Z_min)/Z_pas);
        //wpi = zl / rayleigh;
        //wp  = ww * sqrt ( 1. + wpi*wpi ) ;
        wp = laser_sig[index_laser_pas];
        wp2 = laser_sig2[index_laser_pas];
//Compute laser probability ************** !!!!!!!!!!!!!!!!!!!! **************
        //29.Nov.2021 for Z_NORM_TYPE=1 reject according to LASER prob
        if(Z_NORM_TYPE == 2) fl = (xint-avgX[index_pas])*(xint-avgX[index_pas])/sigmaX2[index_pas] + (yint-avgY[index_pas])*(yint-avgY[index_pas])/sigmaY2[index_pas];
          else fl = (xl*xl + yl*yl)/wp2;
        //laser_prob = exp(-2. * fl);
        laser_prob = exp(-fl/2.);
        l_prob = G4UniformRand();
        cycleXY_counter.IncrementCounter();
        if(Z_NORM_TYPE == 0) acceptXY = false;
          else acceptXY = (l_prob > laser_prob);
        //analysisManager->FillH2(5, xint, yint);
        }
      }  while (acceptXY);
        //} while (0);
        //}
      //cycle_counter.IncrementCounter();
      //laser_strength = 2. / (pi * wp * wp);
      //laser_strength = 1. / (2. * pi * wp2);
    if(Z_NORM_TYPE == 0) LCSweight = laser_prob*overlap_norm[index_pas];
      //29.Nov.2021 for Z_NORM_TYPE=1 reject according to LASER prob
      //else if(Z_NORM_TYPE == 1) LCSweight = laser_prob;
        else LCSweight = 1.;
    //LCSweight = laser_prob*overlap_norm[index_pas];
    //LCSweight = laser_strength * laser_prob;
    //if(!IDEAL_LASER_SW) LCSweight = 1.;
    analysisManager->FillH1(7, zint);
    analysisManager->FillH1(8, zint, LCSweight);
    analysisManager->FillH2(5, xint, yint);
    analysisManager->FillH2(6, xint, yint, LCSweight);
/*=========================================================================*/
    if(!SOURCE_TYPE) LCSweight = 1.;

//Compute electron probability
     /*sigx = sqrt(emmitx*betax);
     fx=xint/(sigx);
     sigy = sqrt(emmity*betay);
     fy=yint/(sigy);
     //Compute density of electrons
     electron_strength = 1. / (2.*pi*sigx*sigy);
     electron_prob =  electron_strength * exp(-fx*fx/2.)*exp(-fy*fy/2.);*/


/*______________Generate Electron and Photon Lorentz vectors_______________*/
//Compute Alpha Twiss parameters for electron beam
    //alphax = -zint / betax0;
    //alphay = -zint / betay0;
//Compute electron angles
    //xang = -sqrt(2.*ux*emmitx/betax)*(alphax*cos(phix)+sin(phix));
    //yang = -sqrt(2.*uy*emmity/betay)*(alphay*cos(phiy)+sin(phiy));
    analysisManager->FillH2(7, xint, xang, LCSweight);
    analysisManager->FillH2(8, yint, yang, LCSweight);
//Generate electron direction
    if(!IDEAL_ELE_SW)
      ele_theta0=ele_phi0=0.;
      else
      {
      //kapa = 1./sqrt(tan(xang)*tan(xang)+tan(yang)*tan(yang)+1.);
      //ele_dir0.set( tan(xang)*kapa, tan(yang)*kapa, kapa);
      kapa = 1./sqrt(xang*xang+yang*yang+1.);
      ele_dir0.set( xang*kapa, yang*kapa, kapa);
      ele_theta0 = ele_dir0.getTheta();
      ele_phi0 = ele_dir0.getPhi();
      if(ele_phi0<0.) ele_phi0=2.*pi+ele_phi0;
      //G4cout<<"X = "<<xang<<"; Y = "<<yang<<"; Phy = "<< ele_phi0*180/pi<<G4endl;
      }
    //analysisManager->FillH1(5,ele_phi0);
//Generate electron energy
    rnd = G4RandGauss::shoot(1.0,deene);
    eelec = eene * rnd + electron_mass_c2;
    Gamma = eelec / electron_mass_c2;
    Beta = sqrt(1. - 1./(Gamma*Gamma));
//Generate electron Lorentz vector
    eimp0.set(0.,0.,std::sqrt(eelec*eelec-electron_mass_c2*electron_mass_c2),eelec);

//Electron source ONLY! No LCS interaction
    if(!SOURCE_TYPE)
      {
      eimp1 = eimp0;
      eimp1.rotateY(ele_theta0);
      eimp1.rotateZ(ele_phi0);
      }
    else
      {
    taurot = 0;
    if(PLIN_MODEL_SW==1)
      {
      if (G4UniformRand()<Plin) taurot = 0.;
        else taurot = 2. * pi * G4UniformRand();
      }
    else if(PLIN_MODEL_SW==2)
      {
      if (G4UniformRand()<Plin) taurot = 0.;
        else taurot = PolAngleGenerator(1.-Plin);
      }
    //analysisManager->FillH1(9,taurot+taui);
        
//Generate photon direction and polarization Lorentz vector
    ey.set(0.,1.,0.);
    ey.rotateZ(taurot+taui);
    taurot = taui;
    if(!IDEAL_LASER_SW)
        {
        laser_dir0.set(0.,0.,1.);
        //pol0.set(laser_dir0.cross(ey),0.);
        //Rotate laser beam 180 deg against electron beam
        //laser_theta = 180.*deg - laser_theta0;
        //laser_phi = 180.*deg + laser_phi0;
        /*laser_dir0.rotateY(180.*deg - laser_theta0);
        laser_dir0.rotateZ(laser_phi0 + 180.*deg);
        laser_theta = laser_dir0.getTheta();
        laser_phi = laser_dir0.getPhi();*/
        }
    else
        {
        laser_dir0.set ( xl*zl/(rayleigh*rayleigh+zl*zl), yl*zl/(rayleigh*rayleigh+zl*zl), 1. );
        //G4cout<<xl*zl/(rayleigh*rayleigh+zl*zl)<<"; "<<yl*zl/(rayleigh*rayleigh+zl*zl)<<"; "<<1. <<G4endl;
        //G4double angle0 = laser_dir0.getPhi();
        //if(angle0<0.) angle0=2.*pi+angle0;
        //analysisManager->FillH1(5,angle0);
        //pol0.set(laser_dir0.cross(ey),0.);
        //pol0.set(ey.cross(laser_dir0),0.);
        //G4double angle0 = pol0.getV().getPhi();
        //if(angle0<0.) angle0=2.*pi+angle0;
        //analysisManager->FillH1(5,angle0);
        //Rotate laser beam 180 deg against electron beam
        /*laser_dir0.rotateY(180.*deg - laser_theta0);
        laser_dir0.rotateZ(laser_phi0 + 180.*deg);
        laser_theta = laser_dir0.getTheta();
        laser_phi = laser_dir0.getPhi();*/
        }
    laser_dir0.rotateY(180.*deg - laser_theta0);
    laser_dir0.rotateZ(laser_phi0 + 180.*deg);
    laser_theta = laser_dir0.getTheta();
    laser_phi = laser_dir0.getPhi();

//Generate photon Lorentz vector
    eph = eph0 * G4RandGauss::shoot(1.0,deph/100.);
    phimp0.set(0.,0.,eph,eph);
    //phimp0.set(0.,0.,eph0,eph0);
//Rotate photon Lorentz vector according to laser photon direction
    phimp0.rotateY(laser_theta);
    phimp0.rotateZ(laser_phi);
    analysisManager->FillH2(9, phimp0.getX()/eph, phimp0.getY()/eph);

//Rotate polarization Lorentz vector according to laser photon direction
//    pol0.rotateY(laser_theta0);
    //analysisManager->FillH1(5,pol0.getV().getPhi());
//    pol0.rotateZ(laser_phi0);
//    pol0.rotateY(laser_theta);
//    pol0.rotateZ(laser_phi);
    pol0.set(((phimp0.getV()).cross(ey)).unit(),0.);
    analysisManager->FillH1(9,(pol0.getV()).getPhi());

/*    taurot = 0;
    if(PLIN_MODEL_SW==1)
      {
      if (G4UniformRand()<Plin) taurot = 0.;
        else
	{
	taurot = 2. * pi * G4UniformRand();
        pol0.rotate(phimp0.getV(),taurot);
        }
      }
    else if(PLIN_MODEL_SW==2)
      {
      if (G4UniformRand()<Plin) taurot = 0.;
        //else taurot = PolAngleGenerator(std::sqrt(1.-Plin*Plin));
        else
	{
	taurot = PolAngleGenerator(1.-Plin);
        pol0.rotate(phimp0.getV(),taurot);
        }
      }
    analysisManager->FillH1(9,taurot+taui);
    taurot = taui;
    pol0.rotate(phimp0.getV(),-taurot);
*/    
    //G4double angle0 = pol0.getV().getPhi();
    //if(angle0<0.) angle0=2.*pi+angle0;
    //analysisManager->FillH1(5,taurot);
/*==========================================================================*/


/*___________Rotation with spherical angles of electron in order to have electron momentum along Z axis towards +____________*/
    phimp0.rotateZ(-ele_phi0);
    phimp0.rotateY(-ele_theta0);
    pol0.rotateZ(-ele_phi0);
    pol0.rotateY(-ele_theta0);
    //G4double angle0 = pol0.getV().getPhi();
    //if(angle0<0.) angle0=2.*pi+angle0;
    //analysisManager->FillH1(5,angle0);


/*___________Lorentz transformation into electron rest reference frame_______________________________________________________*/
    phimp0.boostZ(-Beta);
    pol0.boostZ(-Beta);
    eimp0.boostZ(-Beta);

/*___________Rotation with spherical angles of photon in order to have photon momentum along Z axis towards +________________*/
    theta = phimp0.theta();
    phi = phimp0.phi();
    phimp0.rotateZ(-phi);
    phimp0.rotateY(-theta);
    pol0.rotateZ(-phi);
    pol0.rotateY(-theta);
    //G4double angle0 = pol0.getV().getPhi();
    //if(angle0<0.) angle0=2.*pi+angle0;
    //analysisManager->FillH1(5,angle0);
/*___________Register this direction of incident photon______________________________________________________________________*/
    phot_dir0.set(phimp0.getX(),phimp0.getY(),phimp0.getZ());

/*___________Gauge transformation of polarization vector_____________________________________________________________________*/
    fgauge = pol0.getT() / phimp0.getT();
    pol1 = pol0 - fgauge*phimp0;
    pol1dir.set(pol1.getX(),pol1.getY(),pol1.getZ());
    pol1dir = pol1dir.unit();
    tau = pol1dir.getPhi();
    //tau=0;
    //if(tau<0.) tau=2.*pi+tau;
    //G4cout<<"Tau = "<<tau<<G4endl;
    //analysisManager->FillH1(5,tau);

//phimp0.setT(1.*electron_mass_c2);
/*___________Find minimum and maximum of the scattered gamma rays in electron rest frame_____________________________________*/
    egmin = phimp0.getT()/(1.+2.*phimp0.getT()/electron_mass_c2);
    egmax = phimp0.getT();

/*============================================= Klein - Nishina random generator ===========================================*/
/*___________Generate Theta and Energy of Gamma ray_________________________________________________________________________*/
    KNRand(phimp0.getT(), egmin, egmax, &egam, &thrnd);
    CosTheta = cos(thrnd);
    CosThetaSq = CosTheta*CosTheta;
    SinThetaSq = 1.-CosThetaSq;
    //thrnd = pi*G4UniformRand();
    //egam = phimp0.getT()/(1.+phimp0.getT()*(1.-std::cos(thrnd)));
    //G4cout<<phimp0.getT()/eV<<G4endl;
    //analysisManager->FillH1(5,thrnd);
/*___________Exit Lorentz vector of Gamma photon____________________________________________________________________________*/
    phimp1.set(0.,0.,egam,egam);
/*___________Generate Phi of Gamma ray photon_______________________________________________________________________________*/
    x = electron_mass_c2 / phimp0.getT();
    y = electron_mass_c2 / egam;
    t = x - y;
    t = t * t + 2. * t;
    G4double epsilon = egam/phimp0.getT();
    u = epsilon + 1./epsilon;

    //Phi_gen ( t, u, tau, sqrt(PolLinValue), &Phi_ang);
    Phi_gen ( t, u, PolLinValue, &Phi_ang);
    //Phi_ang = 2.*pi*G4UniformRand();
    //analysisManager->FillH1(5,Phi_ang);
    Phi_ang_v = Phi_ang - pi/2.;
    if(Phi_ang_v<0.) Phi_ang_v += 2.*pi;
    SinPhi = sin(Phi_ang);
    SinPhiSq = SinPhi*SinPhi;


//tau=0.;
//G4ThreeVector gammaPolarization1 = SetNewPolarization(epsilon, thrnd, tau-Phi_ang);
//gammaPolarization1 = SetNewPolarization(epsilon, thrnd, Phi_ang, &pBeta);
SetNewPolarization(VectorModel, epsilon, Plin, thrnd, Phi_ang_v, &pBeta);

/*___________Generate Stokes parameters ____________________________________________________________________________________*/
    //Stokes1 = -sqrt(Plin)*cos(2.*(tau-Phi_ang));
    //Stokes1 = -sqrt(Plin)*cos(2.*Phi_ang);
    Stokes1 = Plin*cos(2.*Phi_ang);
    //Stokes2 = sqrt(Plin)*sin(2.*(tau-Phi_ang));
    //Stokes2 = sqrt(Plin)*sin(2.*Phi_ang);
    Stokes2 = -Plin*sin(2.*Phi_ang);
    Stokes3 = Pcirc;

    F1 = SinThetaSq;
    F0 = u - F1;
    //F11 is different in M. Boca's report:
    //F11 = 1. + cos(thrnd);
    // This MUST be checked !!!
    //01.Dec.2021 Checked: Formula is good (taken from Landau), M. Boca formula is wrong
    //!!Landau convention:
    //Linear pol: eps3(Landau) = eps1(M.Boca); eps1(Landau) = eps2(M.Boca)
    //Circular pol: eps2(Landau)
    F11 = 1. + CosThetaSq;
    F22 = 2. * CosTheta;
    F33 = u * CosTheta;
    Fnum = F0 + F1*Stokes1;

/*==========================================================================================================================*/

/*___________Rotate exit Lorentz vector of Gamma photon according to Klein - Nishina results________________________________*/
    phimp1.rotateY(thrnd);
    phimp1.rotateZ(Phi_ang_v);

    gammaPolarization1.set(1.,0.,0.);
    versorPerp = (phimp1.getV()).cross(gammaPolarization1).unit();
    versorPar = (versorPerp.cross(phimp1.getV())).unit();
    //gammaPolarization1 = versorPar.rotate(phimp1.getV(),0.);
    gammaPolarization1 = versorPar.rotate(phimp1.getV(),pBeta);
    gammaPolarization1.rotateZ(tau);

    phimp1.rotateZ(tau);


/*___________Generate Lorentz vector of outgoing electron___________________________________________________________________*/
    eimp1 = eimp0 + phimp0 - phimp1;

 /*___________Register this direction of incident photon______________________________________________________________________*/
    phot_dir1.set(phimp1.getX(),phimp1.getY(),phimp1.getZ());

    eX.set(1.,0.,0.);
    //Compute the versor perpendicular on scattering plane (epsilon1)
    //epsilon1 = phot_dir0.cross(phot_dir1);
    //epsilon1 = (1./epsilon1.getR())*epsilon1;
    epsilon1 = (phot_dir0.cross(phot_dir1)).unit();  // this is eX versor in Landau
    //Compute the versor paralel with (contained into) scattering plane (epsilon2)
    //epsilon2 = phot_dir1.cross(epsilon1);
    //epsilon2 = (1./epsilon2.getR())*epsilon2;
    epsilon2 = (phot_dir1.cross(epsilon1)).unit();  //this is eY' versor in Landau

    LgammaPolarization1.set(gammaPolarization1,0.);
    Lepsilon1.set(epsilon1,0.);
    Lepsilon2.set(epsilon2,0.);

    if(!StokesEleModel)
      {
      Stokes1f = (F1 + F11*Stokes1) / Fnum;
      Stokes2f = F22*Stokes2 / Fnum;
      PolDegree = sqrt(Stokes1f*Stokes1f+Stokes2f*Stokes2f);
      ratioPerpPar = (u+2.*Stokes1)/(u-2.+2.*CosThetaSq*(1.-Stokes1));
      }
      else
      {
      SinThSinPh2 = SinThetaSq*SinPhiSq;
      PDterm1 = 4.*(Plin-SinThSinPh2)*(Plin-SinThSinPh2);
      PDterm2 = (1.-Plin)*SinThetaSq;
      PDterm3 = PDterm2 + 4.*Plin*(1.+SinThSinPh2)-4.*(Plin+1.)*SinThSinPh2*SinPhiSq;
      PolDegree = sqrt(PDterm1+PDterm2*PDterm3)/(u-PDterm2-2.*Plin*SinThSinPh2);
      //PolDegree = 1.;
      //Compute the final polarization component along epsilon1 versor
      Pol_eps1 = gammaPolarization1.dot(epsilon1);
      //Compute the final polarization component along epsilon2 versor
      Pol_eps2 = gammaPolarization1.dot(epsilon2);
      Stokes1f = PolDegree*(Pol_eps1*Pol_eps1 - Pol_eps2*Pol_eps2)/(Pol_eps1*Pol_eps1 + Pol_eps2*Pol_eps2);
      Stokes2f = 2.*PolDegree*Pol_eps1*Pol_eps2/(Pol_eps1*Pol_eps1 + Pol_eps2*Pol_eps2);
      epsilonPerp = Pol_eps1*epsilon1;  // this is polarization component on eX versor
      //Scalar product between initial polarization and final polarization perpendicular component
      termPerp = pol1dir.dot(epsilonPerp);
      //termPerp = eX.dot(epsilonPerp);
      termPerp = termPerp*termPerp;
      epsilonPar = Pol_eps2*epsilon2;  // this is polarization component on eY' versor
      //Scalar product between initial polarization and final polarization paralel component
      termPar = pol1dir.dot(epsilonPar);
      //termPar = eX.dot(epsilonPar);
      termPar = termPar*termPar;
      ratioPerpPar = (u-2.+4.*termPerp)/(u-2.+4.*termPar);
      }
    //Stokes1f = gammaPolarization1.getX()*gammaPolarization1.getX()-gammaPolarization1.getY()*gammaPolarization1.getY();
    //Stokes2f = 2.*gammaPolarization1.getX()*gammaPolarization1.getY();
    Stokes3f = F33*Stokes3 / Fnum;
    Stokes.set(Stokes1f,Stokes2f,Stokes3f);
    analysisManager->FillH2(19, thrnd, Phi_ang_v,ratioPerpPar);
    analysisManager->FillH2(23, thrnd, Phi_ang_v,LCSweight);
    analysisManager->FillH2(24, thrnd, Phi_ang_v,PolDegree);
    
/*___________Rotate back with spherical angles of photon done in order to have photon momentum along Z axis towards +_______*/
    phimp1.rotateY(theta);
    phimp1.rotateZ(phi);
    eimp1.rotateY(theta);
    eimp1.rotateZ(phi);
    LgammaPolarization1.rotateY(theta);
    LgammaPolarization1.rotateZ(phi);
    Lepsilon1.rotateY(theta);
    Lepsilon1.rotateZ(phi);
    Lepsilon2.rotateY(theta);
    Lepsilon2.rotateZ(phi);

/*___________Boost the Gamma ray photon along Z axis with Beta value of electrons___________________________________________*/
    phimp1.boostZ(Beta);
    eimp1.boostZ(Beta);
    LgammaPolarization1.boostZ(Beta);
    Lepsilon1.boostZ(Beta);
    Lepsilon2.boostZ(Beta);

/*___________Rotate final Lorentz vector of Gamma ray photon according to initial electron direction________________________*/
    phimp1.rotateY(ele_theta0);
    phimp1.rotateZ(ele_phi0);
    eimp1.rotateY(ele_theta0);
    eimp1.rotateZ(ele_phi0);
    LgammaPolarization1.rotateY(ele_theta0);
    LgammaPolarization1.rotateZ(ele_phi0);
    Lepsilon1.rotateY(ele_theta0);
    Lepsilon1.rotateZ(ele_phi0);
    Lepsilon2.rotateY(ele_theta0);
    Lepsilon2.rotateZ(ele_phi0);

    fgauge = LgammaPolarization1.getT() / phimp1.getT();
    LgammaPolarization1 = LgammaPolarization1 - fgauge*phimp1;
    gammaPolarization1 = (LgammaPolarization1.getV()).unit();
    //analysisManager->FillH2(13,phimp1.getT(), gammaPolarization1.getPhi(), 1.);
    fgauge = Lepsilon1.getT() / phimp1.getT();
    Lepsilon1 = Lepsilon1 - fgauge*phimp1;
    fgauge = Lepsilon2.getT() / phimp1.getT();
    Lepsilon2 = Lepsilon2 - fgauge*phimp1;
    
    ez = phimp1.getV().unit();
    eZ.set(0.,0.,1.);
    ex = (eZ.cross(ez)).unit();

    if(StokesLabModel)
      {
      //epsilon1 = (eZ.cross(ez)).unit();
      epsilon1 = ex;
      Pol_eps1 = gammaPolarization1.dot(epsilon1);
      epsilon2 = (ez.cross(epsilon1)).unit();
      Pol_eps2 = gammaPolarization1.dot(epsilon2);
      //Stokes1f = (Pol_eps1*Pol_eps1 - Pol_eps2*Pol_eps2)/(Pol_eps1*Pol_eps1 + Pol_eps2*Pol_eps2);
      Stokes1f = PolDegree*(Pol_eps1*Pol_eps1 - Pol_eps2*Pol_eps2);
      //Stokes2f = 2.*Pol_eps1*Pol_eps2/(Pol_eps1*Pol_eps1 + Pol_eps2*Pol_eps2);
      Stokes2f = 2.*PolDegree*Pol_eps1*Pol_eps2;
      Stokes.set(Stokes1f,Stokes2f,Stokes3f);
      
      epsilon1.set(1.,0.,0.);
      Pol_eps1 = gammaPolarization1.dot(epsilon1);
      epsilon2.set(0.,1.,0.);
      Pol_eps2 = gammaPolarization1.dot(epsilon2);
      //Stokes1L = (Pol_eps1*Pol_eps1 - Pol_eps2*Pol_eps2)/(Pol_eps1*Pol_eps1 + Pol_eps2*Pol_eps2);
      Stokes1L = PolDegree*(Pol_eps1*Pol_eps1 - Pol_eps2*Pol_eps2);
      //Stokes2L = 2.*Pol_eps1*Pol_eps2/(Pol_eps1*Pol_eps1 + Pol_eps2*Pol_eps2);
      Stokes2L = 2.*PolDegree*Pol_eps1*Pol_eps2;
      StokesLab.set(Stokes1L,Stokes2L,Stokes3L);
      }
      else
      {
      StokesLab = Stokes;
      //Cosinus = (ez.cross(ex)).dot((Lepsilon2.getV()).unit());
      Cosinus = ex.dot((Lepsilon1.getV()).unit());
      if(Cosinus<-1.) Cosinus = -1.;
      if(Cosinus>1.) Cosinus = 1.;
      angle = std::acos(Cosinus);
      //if(0.5<G4UniformRand()) angle = -angle;
      if((ex.cross(Lepsilon1.getV())).getZ() < 0.) angle = -angle;
      //angle = ex.angle((Lepsilon1.getV()).unit());
      //G4cout<<angle<<G4endl;
      //analysisManager->FillH1(5,angle);
      Stokes.rotateZ(2.*angle);
      //StokesLab.rotateZ(2.*angle);
      angle = Lepsilon1.phi();
      StokesLab.rotateZ(2.*angle);
      }

/*===========Set Direction Energy and Emission point position of the Laser-Compton Scattering Gamma Ray Source==============*/
    direction=phimp1.getV();
    energy = phimp1.getT();

    if(StokesLabGraph)
      {
      analysisManager->FillH2(10, energy, StokesLab.getX(), 1.);
      analysisManager->FillH2(11, energy, StokesLab.getY(), 1.);
      analysisManager->FillH2(12, energy, sqrt(StokesLab.getX()*StokesLab.getX()+StokesLab.getY()*StokesLab.getY()), 1.);
      //analysisManager->FillH2(12, energy, PolDegree, 1.);
      analysisManager->FillH2(13, energy, Stokes.getZ(), 1.);
      }
      else
      {
      analysisManager->FillH2(10, energy, Stokes.getX(), 1.);
      analysisManager->FillH2(11, energy, Stokes.getY(), 1.);
      analysisManager->FillH2(12, energy, sqrt(Stokes.getX()*Stokes.getX()+Stokes.getY()*Stokes.getY()), 1.);
      analysisManager->FillH2(13, energy, Stokes.getZ(), 1.);
      }

    //G4double angle = 2.*(-Lepsilon1.phi()+laser_dir0.getPhi());
//    G4double angle = Lepsilon1.phi() - laser_phi;
      
    //G4double angle = pi/4.;
    //G4cout<< (ez.cross(ey)).getX() << "; " << (ez.cross(ey)).getY() << "; " << (ez.cross(ey)).getZ() << "; " << sqrt((ez.cross(ey)).getX()*(ez.cross(ey)).getX() + (ez.cross(ey)).getY()*(ez.cross(ey)).getY() + (ez.cross(ey)).getZ()*(ez.cross(ey)).getZ()) << G4endl; 
      //G4cout << cos(Lepsilon1.phi()) << "; " << angle << G4endl;
/*    if(Lepsilon1.getY()>=0.)
    //if(ey.dot(Lepsilon2)>=0.)
        angle = -2.*(90.*deg-acos(angle));
        else
        angle = -2.*(90.*deg+acos(angle));*/
    //analysisManager->FillH1(5,angle*180./pi);
    //if(angle<0.) angle=2.*pi+angle;
    
    //analysisManager->FillH1(5,angle);
    //Stokes.rotateZ(2.*angle);
    }

    ele_direction = eimp1.getV();
    //12 Dec 2021 - Subtract the electron rest mass energy in ordet to pass only kinetic energy into PrimaryGenerator
    ele_energy = eimp1.getT() - electron_mass_c2;
    G4ThreeVector pos2(xint + xi + dispXer, yint + yi, zint + zi);
    position = pos2;
//Compute WEIGHT of the event and  transmi total nr. of cycles
    //LCSweight = electron_prob * laser_strength * laser_prob;
    //LCSweight = electron_prob * laser_strength;
    //LCSweight = electron_strength * laser_strength * laser_prob;
    //LCSweight = electron_strength * laser_strength;
    //LCSweight = laser_strength * laser_prob;
    Total_XYcycles = cycleXY_counter.GetCounterValue();
    Total_Zcycles = cycleZ_counter.GetCounterValue();
    Total_electron = electron_counter.GetCounterValue();
    //norm = ((G4double)Total_cycles / (G4double)Total_electron) *  ((zmax - zmin) / deltaZ_max) * 100.;

}

G4double eliLaBr_GammaSource::PolAngleGenerator(G4double axis)
{
  G4double rnum, angle, prob;
  G4double cosAngle, cosAngle2, sinAngle2;
  
  do
    {
    rnum = G4UniformRand();
    angle = pi*G4UniformRand();
    cosAngle = std::cos(angle);
    cosAngle2 = cosAngle*cosAngle;
    sinAngle2 = 1. - cosAngle2;
    prob = axis/std::sqrt(cosAngle2+axis*axis*sinAngle2);
    } while (rnum>prob);
  return angle-pi/2.;
}

void eliLaBr_GammaSource::KNRand(G4double ep, G4double egmin, G4double egmax, G4double *eg, G4double *thrnd)
{
    G4double rnum, x, y, t, u, feg;
    G4double kgmin, kgmax, kphot, kgam;
    G4double EneRatio, term, KNmax;
    
    kgmin = egmin/electron_mass_c2;
    kgmax = egmax/electron_mass_c2;
    kphot = ep / electron_mass_c2;
    x = 1. / kphot;
    term = 1. + 2.*kphot;
    KNmax = term + 1./term;
    // OLD
    //KNmax = 2. + 2./x;

//    do {
    do {
      rnum = G4UniformRand();
      //*eg = egmin + (egmax - egmin) * rnum;
      kgam = kgmin + (kgmax - kgmin) * rnum;

      //y = electron_mass_c2 / *eg;
      y = 1./kgam;
      t = x - y;
      EneRatio = x/y;
      //u = ep / *eg + *eg / ep;
      u = EneRatio + 1./EneRatio;
      //feg = (t*t + 2.*t + u)/(2. + 2./x);
      feg = (t*t + 2.*t + u)/KNmax;
      rnum = G4UniformRand();
      } while (rnum>feg);

    *eg = kgam*electron_mass_c2;
    //*thrnd = std::acos(1. - electron_mass_c2*(1./ *eg - 1./ep));
    *thrnd = std::acos(1. + t);
//    } while((89.<*thrnd/deg)&&(*thrnd/deg<91.));
}

void eliLaBr_GammaSource::Phi_gen(G4double t, G4double u, G4double Pt, G4double *Phi_ang)
//void eliLaBr_GammaSource::Phi_gen(G4double t, G4double u, G4double tau, G4double Pt, G4double *Phi_ang)
{
    G4double r5, rnum, g_phi;
    G4double g_phiMax;
    if(t<0) g_phiMax = (1. - Pt) *  t + u;
      else g_phiMax = (1. + Pt) *  t + u;
    //OLD
    //g_phiMax = 2.*pi*(t + u);
    do {
      r5   = G4UniformRand();
      rnum = G4UniformRand() * 2. * pi;
      //g_phi = 2. * cos ( tau - rnum ) * cos ( tau - rnum ) *  t + u;
      //g_phi = (1. + Pt * cos(2.*(tau - rnum)) ) *  t + u;
      g_phi = (1. - Pt * cos(2.*rnum) ) *  t + u;
      //g_phi = g_phi /( t + u );
      //g_phi = g_phi / pi / 2. ;
      g_phi /= g_phiMax;
      } while ( r5 >= g_phi );

    *Phi_ang = rnum;
}

//G4ThreeVector eliLaBr_GammaSource::SetNewPolarization(G4double epsilon, G4double ThetaKN, G4double phi, G4double *pBeta)
void eliLaBr_GammaSource::SetNewPolarization(G4bool VectorModel, G4double epsilon, G4double Pol, G4double ThetaKN, G4double phi, G4double *pBeta)
{
  //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4double rand1;
  G4double rand2;
  G4double cosPhi = std::cos(phi);
  //G4double sinPhi = std::sin(phi);
  G4double costheta = std::cos(ThetaKN);
  G4double sinSqrTh = 1. - costheta*costheta;
  //G4double sinTheta = std::sqrt(sinSqrTh);
  G4double cosSqrPhi = cosPhi*cosPhi;
  //  G4double cossqrth = 1.-sinSqrTh;
  //  G4double sinsqrphi = sinPhi*sinPhi;
  G4double normalisation = std::sqrt(1. - cosSqrPhi*sinSqrTh);
 

  // Determination of Theta 
  
  // ---- MGP ---- Commented out the following 3 lines to avoid compilation 
  // warnings (unused variables)
  G4double thetaProbability;
  G4double theta;
  G4double a, b;
  G4double cosTheta;
  //G4double cosBeta;
  //G4double sinBeta;
  a = 4.*normalisation*normalisation;
  b = epsilon + 1.0/epsilon;

  if(VectorModel)
    {
    //depaola method
    do
      {
      rand1 = G4UniformRand();
      rand2 = G4UniformRand();
      thetaProbability=0.;
      theta = twopi*rand1;
      cosTheta = std::cos(theta);
      thetaProbability = (b -2. + a*cosTheta*cosTheta)/(a+b-2.);
      //cosTheta = std::cos(theta);
      }
      while ( rand2 > thetaProbability );
    //*pBeta = theta;
    //cosBeta = cosTheta;
    }
  else
    // Dan Xu method (IEEE TNS, 52, 1160 (2005))
    {
    rand1 = G4UniformRand();
    rand2 = G4UniformRand();
    //if (rand1<(b-2.)/(2.0*b-4.0*sinSqrTh*cosSqrPhi))
    if (rand1<(b-2.)/(2.0*b-2.0*(1.-Pol*(1.-2.*cosSqrPhi))*sinSqrTh))
      {
      if (rand2<0.5) theta = pi/2.0;
      else theta = 3.0*pi/2.0;
      }
      /*{
      if (rand2<0.5) theta = 0;
      else theta = pi;
      }*/
    else
      {
      if (rand2<0.5) theta = 0;
      else theta = pi;
      }
      /*{
      if (rand2<0.5) theta = pi/2.0;
      else theta = 3.0*pi/2.0;
      }*/
    //cosBeta = std::cos(theta);
    }
  *pBeta = theta;
  //sinBeta = std::sqrt(1-cosBeta*cosBeta);
  
  /*G4ThreeVector gammaPolarization1;

  G4double xParallel = normalisation*cosBeta;
  G4double yParallel = -(sinSqrTh*cosPhi*sinPhi)*cosBeta/normalisation;
  G4double zParallel = -(costheta*sinTheta*cosPhi)*cosBeta/normalisation;
  G4double xPerpendicular = 0.;
  G4double yPerpendicular = (costheta)*sinBeta/normalisation;
  G4double zPerpendicular = -(sinTheta*sinPhi)*sinBeta/normalisation;

  G4double xTotal = (xParallel + xPerpendicular);
  G4double yTotal = (yParallel + yPerpendicular);
  G4double zTotal = (zParallel + zPerpendicular);
  
  gammaPolarization1.setX(xTotal);
  gammaPolarization1.setY(yTotal);
  gammaPolarization1.setZ(zTotal);*/
  
  //G4double WeightPolarization, WeightParallel, WeightPerpendicular, Limit;
  //Limit = 100.;
//  WeightParallel = std::fabs(xParallel+yParallel+zParallel);
  //WeightParallel = b + a*cosTheta*cosTheta;
//  WeightPerpendicular = std::fabs(xPerpendicular+yPerpendicular+zPerpendicular);
  //WeightPerpendicular = b;
  
  //if(WeightParallel==0.) WeightPolarization=Limit;
  //  else WeightPolarization = WeightPerpendicular/WeightParallel;
//    G4cout<<WeightPolarization<<G4endl;
//  WeightPolarization /= 1.e+20;
  //if(WeightPolarization>Limit) WeightPolarization=Limit;
  
  //if(phi<0.) phi = phi+2.*pi;
//  analysisManager->FillH2(18, ThetaKN, phi,WeightPolarization);
//  analysisManager->FillH2(18, ThetaKN, phi,1.);
  
  //return gammaPolarization1;
}


int eliLaBr_GammaSource::HowMany()
{
    if(mean>0) return CLHEP::RandPoisson::shoot(mean);
    else return (int)(abs(mean));
}
