# eliLaBr
GEANT4 LCS simulation on flat eff. n-detector based on 3He counters embedded into polyethylene moderator
========================================================================================================

## Developers

Dan Mihai FILIPESCU[^1] <dan.filipescu@nipne.ro>

Adriana Ioana GHEORGHE[^1] <ioana.gheorghe@nipne.ro>

[^1]: [Horia Hulubei - National Institute for Phisics and Nuclear Engineering, Bucharest-Magurele, ROMANIA](https://www.nipne.ro)

## 0 - Usage

Information regarding building & installation of code can be found in [INSTALL](INSTALL.md) file.

The code can be run both in 

  **batch mode** `eliLaBr commands`, or
  
  **interactive mode**
- `eliLaBr dumb` start a console interactive terminal;
- `eliLaBr` start a graphical Qt OpenGL interactive terminal.

In **interactive mode**, user should start with `control/execute commands`.

### *`commands`* macrofile

  This is the main macrofile of **eliLaBr** code.
  Depending on the purpose, it is recommended to perform the following changes to this macrofile:
  - when we are interested to launch **eliLaBr** in *batch mode* we should comment
  ```
  #/control/execute gui.mac
  #/control/execute vis1.mac
  ```
  and uncomment
  ```
  /run/beamOn #(int)NumberOfEvents
  ```
  - when we are interested to run **eliLaBr** in *interactive mode* we should do the other way arround:
  ```
  /control/execute gui.mac
  /control/execute vis1.mac
  #/run/beamOn #(int)NumberOfEvents
  ```
  
  ## 1 - GEOMETRY DEFINITION
  
  Several geometry macro-commands are implemented in `eliLaBr_DetectorMessenger` class.
  
  A list with these macro-commands can be found in **`SetGeometry.mac`** macrofile which should be loaded into **`commands`** macrofile.
  
  **TARGET** is made *sensitive* and can be used as a detector. Three example cases are listed into **`SetGeometry.mac`**,and the corresponding block of macro-commands can be un-commented depending on the needs:
  - *TARGET INSIDE 4Pi detector* - this is the normal case of a target placed in the middle of the neutron detector. Be carefull on the target dimensions not to exceed the inner diameter of the beamline pipe passing through the neutron detector central axis.
  - *TARGET INSIDE E8 Hall* - this is intended to monitor the beam spectrum and usually used to simulate a beam monitor detector. Do not forget to un-comment *TARGET HOUSING* section macro-command `/eli/det/setTGHousing true` in case you intend to.
  - *MONITOR FOR GAMMA SOURCE* - this is intended to be an imaginary plate placed in front of the gamma-beam prior to collimation, in order to histogram different characteristics of the gamma-beam like spatial, energy and/or polarization distributions.
  
  **Neutron counters** dimensions are defined by those of the detectors existing at NewSUBARU & ELI-NP facilities, namely: 1 inch diameter and 500 mm. length, while the detector housing is considered to be made from stainless-steel.
  
  Two types of sensitive gas are allowed:
  - *BF3 GAS*
  - *3He GAS*
  
  so please un-comment the corresponding block of macro-commands.
  
  There are allowed a maximum number of **three** concentric rings of neutron counters around the central beampipe passing through the axis of the polyehylene moderator cube.
  
  ## 2 - PHYSICS LIST
  
  Several physics macro-commands are implemented in `eliPhysicsListMessenger` class.
  
  A list with these macro-commands can be found in **`SetPhysics.mac`** macrofile which should be loaded into **`commands`** macrofile.
  
  A list with Electromagnetic interactions must always by used from the "*#===> Set Electromagnetic Interaction*" block of macro-commands. If multiple macro-commands from the Electromagnetic block are un-commented, only the last macro-commands will be active. The keywords ending with **_usr** correspond to the electromagnetic physics package classes build within **eliLaBr** simulation code.
  
  Specifying the list with hadronic interactions from "*#===> Set Hadronic Interaction Only*" block of macro-commands is optional. Only one macro-command from the Hadronic block can be un-commented.
  
  If a hadronic physics list is selected, in order to perform a proper simulation of the neutron transport, the following commands must be un-commented in the **`commands`** macrofile:
```
/physics_engine/neutron/energyLimit 0. GeV
/physics_engine/neutron/timeLimit 1500000. ns
```
  ## 3 - HISTOGRAMMING
  
  During the simulation run, several 1D & 2D histograms can be incremented.
  The lists with standard macro-commands implemented by `G4AnalysisManager` controlling activation of each histogram as well as fine tuning of each histogram properties are put together in two files with macro-commands:
  - **`histos.mac`** - this file is recommended to be loaded in **`commands`** macrofile in the normal case when the gamma beam is collimated and sent on target (we can monitor the beam spot on target for example);
  - **`histor_monitor.mac`** - this file is recommended to be loaded in **`commands`** macrofile when we are interested to monitor the gamma-beam properties prior to collimation.
  
  ## 4 - THE PRIMARY GENERATOR
  
  The properties of the source of the primary particles that initiate each simulation event can be specified in the file **`incident_energy.in`**. This file should always be present in the folder from where **eliLaBr** code was launched. Inside this file, the user will find the explanation of each parameter that can be set.
  
  The main source of primary generator consists from gamma rays generated along the straight synchrotron beamline outside the experimental hall. Thus, one important option to be set in **`incident_energy.in`** concerns the use of an ideal pencil-like gamma-beam or to perform a realistic simulation of Laser-Compton scattering (LCS) of eV photons on relativistic electrons considering a Gaussian model to describe the laser and Twiss parameter formalism in order to describe the electron beam. Laser polarization along with its influence is taken into account using the Stokes parameter formalism or the polarization vector formalism. Total and partial luminosity of the gamma-ray beam are also computed and used afterwards to normalize the spectra obtained in simulation. All the settings & parameters concerning the modeling of LCS interaction, along with their description, can be specified in the **`incident_gamma.in`** file. This file should also be present in the folder from where **eliLaBr** code was launched. Tabulated values of Twiss parameters for NewSUBARU synchrotron facility are put in **`NewSUBARU_optics_BL01.txt`** and this file should also be present in the folder where the simulation is performed.
  
  In case of the ideal pencil-like gamma-beam, the energy spectrum can be monochromatic or can be provided numerically into an ASCII file. The name of the ASCII file containing the gamma spectrum have to be specified into another file called **`incident_fileName`** and in such case these files must be present in the folder where the simulation is performed
  
  Besides considering a gamma-ray beam as primary source, the user can also specify in **`incident_energy.in`** file, probabilities of gamma and/or neutron emission for calibration purposes from inside/near detectors. If the sum of the neutron and gamma emission probabilities is less than 1, the remaining fraction up 100% is used to generate gamma-ray beam.
  Calibration gamma photons and/or neutrons can be generated from a planar disc with a specified radius within the target or we can specify a given displacement relative to the target of the planar emission disc.
  Calibration gamma ray spectrum is always monochromatic, and the angular emission distribution can be:
  - isotropic;
  - E1 distribution;
  - M1 distribution
  
  relative to the symmetry axis of the planar emission disc / neutron detector axis.
  
  The neutron spectrum can be monochromatic or according to Weiskopf-Ewing evaporation spectrum with a given temperature in MeV, while the angular emission distribution can be:
  - isotropic;
  - P wave neutrons (L=1).
  
  ## 5 - TRACKING & HITS
  
  In **eliLaBr** code the target and neutron counters are made *Sensitive Detectors*. Thus, for this detectors, a dedicated *Tracker* class is responsible to process each *Step* within an *Event* in order to extract information like:
  - Track ID;
  - Parent ID - the Track ID of the particle that generated the current particle on which *Track* we are currently;
  - Global Time - time since the *Event*, to which this *Track* belongs, is created;
  - Particle Kinetic Energy on current *Step*;
  - Total Energy Deposit - energy deposited along the current *Step*;
  - Particle name;
  - Volume name in which the current *Step* is located;
  - Position of current *Step*, taken from the *PreStepPoint*;
  - Boolean confirming (or not) that particle just Entered in the detector boundary - taken from the *PreStepPoint*;
  - Boolean confirming (or not) that the particle just Exit the detector boundary - taken from the *PostStepPoint*;
  - Replica number of the current detector in which the *Step* is located - this is useful to get the neutron counter index number.
  
  All these information collected from a *Step* residing in a *Sensitive Detector* is put together into a memory unit cell structure constituting the so-called *Hit*. Thus, for an *Event* we are putting together a *Collection (or Vector) of Hits* residing in the *Sensitive Detectors*.
  
  The information from *Hits Collection* is further analyzed and processed in *End Of Event Action* method from the **`EventAction`** class. Thus the energy deposition in time is analysed for each neutron counter separately and and also for target. The **ROOT** class **`TSpectrum`** is used to identify separate peak energy deposited over time. These chained vectors of pairs (DepositedEnergy, Time) are stored in a **`TVectors`** class, specially written for this task. The information written in such vectors is possible to be analyzed later, considering different time integration constants, in a similar way as the *Shaping time* is acting in real electronics, in order to be able to treat energy deposits with enough distance in time between them as different events. For example, if in the same GEANT event, after a first neutron is detected in one counter, if a second neutron is detected by the same counter at an enough time distance from the first neutron (due to longer thermalization time), we can consider two different neutron detection events for the same counter.
  
  In *End Of Event Action* method the energy deposited is integrated on each counter and in the target, the average time of energy deposited is computed (if we had EDep>0) for each counter and for target, and the multiplicity of detectors along with their index number is determined in order to put all these information into vectors.
  
  Also we check if any radiation passes the boundary of a detector/target (towards detector/target), without being interested if any energy was deposited or not, in that particular detector during the entire event, and we store the multiplicity of all detectors which were crossed by any radiation, the detector index number, what kind of radiation was, and with what incident energy.
  
  All information obtained in *End Of Event Action* method by processing the *Hits Collection* is written into a **ROOT** file, event by event, in order to be able to sort the events later, depending on the desired conditions, according to different imagined scenarios. Unfortunately, this event file is always named **`nai_9_MeV.root`**, but this limitation will be eliminated in a future release of **eliLaBr** code.
  
  ## Tvectors library
  
  **`Tvectors`** dynamic library provides the classes
  
  ## Papers
  
  Several papers made use of **eliLaBr** code over the time:
  - I. Gheorghe et al. *"Updated neutron-multiplicity sorting method for producing photoneutron average energies and resolving multiple firing events"*, NIMA **1019** 165867 (2021) [10.1016/j.nima.2021.165867](https://doi.org/10.1016/j.nima.2021.165867)
  - S. Belyshev et al. *"New Bi-209 photodisintegration data and physical criteria of data reliability"*, EPJ Web of Conf. **239** 01031 (2020) [10.1051/epjconf/202023901031](https://doi.org/10.1051/epjconf/202023901031)
  - H. Utsunomiya et al. *"GDR cross sections updated in the IAEA-CRP"*, EPJ Web of Conf. **239** 01002 (2020) [10.1051/epjconf/202023901002](https://doi.org/10.1051/epjconf/202023901002)
  - T. Kawano et al. *"IAEA Photonuclear Data Library 2019"*, Nuclear Data Sheets **163** 109-162 (2020) [10.1016/j.nds.2019.12.002](https://doi.org/10.1016/j.nds.2019.12.002)
  - M. Krzysiek et al. *"PHOTONEUTRON CROSS-SECTION MEASUREMENTS FOR Ho-165 BY THE DIRECT NEUTRON-MULTIPLICITY SORTING AT NEWSUBARU"*, ACTA PHYSICA POLONICA B **50** 487 (2019) [10.5506/APhysPolB.50.487](https://doi.org/10.5506/APhysPolB.50.487)
  - H. Utsunomiya et al. *"gamma-ray strength function for thallium isotopes relevant to the Pb-205-(TI)-T-205 chronometry"*, Phys. Rev. C **99** 024609 (2019) [10.1103/PhysRevC.99.024609](https://doi.org/10.1103/PhysRevC.99.024609)
  - M. Krzysiek et al. *"Simulation of the ELIGANT-GN array performances at ELI-NP for gamma beam energies larger than neutron threshold"*, NIMA **916**  257 (2019) [10.1016/j.nima.2018.11.058](https://doi.org/10.1016/j.nima.2018.11.058)
  - T. Renstrom et al. *"Verification of detailed balance for gamma absorption and emission in Dy isotopes"*, Phys. Rev. C **98**  054310 (2018) [10.1103/PhysRevC.98.054310](https://doi.org/10.1103/PhysRevC.98.054310)
  - H. Utsunomiya et al. *"Photon-flux determination by the Poisson-fitting technique with quenching corrections"*, NIMA **896** 103 (2018) [10.1016/j.nima.2018.04.021](https://doi.org/10.1016/j.nima.2018.04.021)
  - H. Utsunomiya et al. *"Photoneutron Reaction Data for Nuclear Physics and Astrophysics"*, EPJ Web of Conf. **178** 06003 (2018) [10.1051/epjconf/201817806003](https://doi.org/10.1051/epjconf/201817806003)
  - G. Gosta et al. *"Response function and linearity for high energy gamma-rays in large volume LaBr3:Ce detectors"*, NIMA **879** 92 (2018) [10.1016/j.nima.2017.10.018](https://doi.org/10.1016/j.nima.2017.10.018)
  - H. Utsunomiya et al. *"Direct neutron-multiplicity sorting with a flat-efficiency detector"*, NIMA **871** 135 (2017) [10.1016/j.nima.2017.08.001](https://doi.org/10.1016/j.nima.2017.08.001)
  - I. Gheorghe et al. *"Photoneutron cross-section measurements in the Bi-209(gamma, xn) reaction with a new method of direct neutron-multiplicity sorting"*, Phys. Rev. C **96** 044604 (2017) [10.1103/PhysRevC.96.044604](https://doi.org/10.1103/PhysRevC.96.044604)
  - I. Gheorghe et al. *"Partial photoneutron cross section measurements on Bi-209"*, EPJ Web of Conf. **146** 05011 (2017) [10.1051/epjconf/201714605011](https://doi.org/10.1051/epjconf/201714605011)
  - H. Utsunomiya et al. *"A unified understanding of (gamma,n) and (n,gamma) reactions and direct neutron-multiplicity sorting"*, EPJ Web of Conf. **146** 05002 (2017) [10.1051/epjconf/201714605002](https://doi.org/10.1051/epjconf/201714605002)
  - T. Renstrom et al. *"Low-energy enhancement in the gamma-ray strength functions of Ge-73,Ge-74"*, Phys. Rev. C **93** 064302 (2016) [10.1103/PhysRevC.93.064302](https://doi.org/10.1103/PhysRevC.93.064302)
  - W. Luo et al. *"Estimates for production of radioisotopes of medical interest at Extreme Light Infrastructure - Nuclear Physics facility"*, Applied Physics B **122** 8 (2016) [10.1007/s00340-015-6292-9](https://doi.org/10.1007/s00340-015-6292-9)
  
