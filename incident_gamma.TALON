14.
200.
982.428
0.04
1.058
1.02
532.
0.05
0.
0.
2548.
0.
0.
-10.10
10.10
-0.3
5.0e-05
-0.18
5.0e-06
0.
0.
0.
1.
0.
1
1
-1
2


    in >> laser_power;      laser_power = laser_power * watt;
    in >> elec_intensity;   elec_intensity = elec_intensity * 1.e-3;
    in >> eene;             eene = eene * MeV;
    in >> deene;
    eelec = eene + CLHEP::electron_mass_c2;
    Gamma = eelec / CLHEP::electron_mass_c2;
    Beta = sqrt(1. - 1./(Gamma*Gamma));
    ele_speed = Beta * CLHEP::c_light;
    G4cout<<"Electron speed: "<<ele_speed*s/m<<G4endl;

    in >> ww;               ww = ww * mm;   // Laser beam waist
    in >> m2;                               // Laser beam quality factor, adimensional
    in >> lwl;              lwl = lwl * nm;
    in >> deph;
    in >> dispX;            dispX = dispX * mm;
    in >> dispY;            dispY = dispY * mm;
    in >> dispZ;            dispZ = dispZ * mm;
    in >> xi;               xi = xi * mm;
    in >> yi;               yi = yi * mm;
    in >> zmin;             zmin = zmin * m;
    in >> zmax;             zmax = zmax * m;
    in >> sigx0;            sigx0 = sigx0 * mm;
    in >> emmitx;           emmitx = emmitx * mm / (Beta*Gamma);
    betax0 = sigx0 * sigx0 / emmitx;
    in >> sigy0;            sigy0 = sigy0 * mm;
    in >> emmity;           emmity = emmity * mm / (Beta*Gamma);
    betay0 = sigy0 * sigy0 / emmity;
    in >> laser_theta0;     laser_theta0 = laser_theta0 * deg;
    in >> laser_phi0_ini;    //laser_phi0 = laser_phi0 * deg;
    in >> taui;              taui = taui * deg;
    in >> Plin;
    in >> Pcirc;
    in >> IDEAL_ELE_SW;
    in >> IDEAL_LASER_SW;
    in >> SOURCE_TYPE;  -1 Just photons
                         0 Just electrons
			+1 photons & electrons
    in >> Z_NORM_TYPE  0 Accept any uniform random Z and multiply the matrices with Electron-Laser overlap integral
                       1 Accept/Reject Z values proportional with Electron-Laser overlap integral
		       2 Same rejection method as "1" used additional for selection of XY coord. of interaction point
		         In this latter case, LCSweight will be 1.

