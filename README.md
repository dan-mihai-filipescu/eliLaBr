# eliLaBr
GEANT4 LCS simulation on flat efficiency neutron detector based on <sup>3</sup>He counters embedded into polyethylene moderator
===================================================================================================================

[![DOI](https://zenodo.org/badge/520783465.svg)](https://zenodo.org/badge/latestdoi/520783465)
(DOI for the latest release of the `code`)

## Developers

Dan Mihai FILIPESCU[^1] <dan.filipescu@nipne.ro>

Adriana Ioana GHEORGHE[^1] <ioana.gheorghe@nipne.ro>

[^1]: [Horia Hulubei - National Institute for Physics and Nuclear Engineering, Bucharest-Magurele, ROMANIA](https://www.nipne.ro)

## How to cite

The code is described in the following articles: [JINST](https://doi.org/10.1088/1748-0221/17/11/P11006) ([arXiv](https://doi.org/10.48550/arXiv.2210.14669)) & [NIM](https://doi.org/10.1016/j.nima.2022.167885) ([arXiv](https://doi.org/10.48550/arXiv.2211.14650))
```
@article{Filipescu_2022,
doi = {10.1088/1748-0221/17/11/P11006},
url = {https://dx.doi.org/10.1088/1748-0221/17/11/P11006},
year = {2022},
month = {nov},
publisher = {IOP Publishing},
volume = {17},
number = {11},
pages = {P11006},
author = {Dan Filipescu},
title = {Monte Carlo simulation method of polarization effects in Laser Compton Scattering on relativistic electrons},
journal = {Journal of Instrumentation},
abstract = {Quasi-monochromatic, high energy and highly polarized γ-ray beam sources based on Compton scattering of laser photons (LCS) on relativistic electrons have developed for the last few decades as established instruments for nuclear physics studies. Following an extensive photoneutron experimental campaign at the LCS γ-ray beam line of the NewSUBARU synchrotron radiation facility at SPring8, Japan, a dedicated simulation code was developed for characterizing the incident γ-ray beams. The eliLaBr code is implemented using Geant4 and is available on the GitHub repository (github.com/dan-mihai-filipescu/eliLaBr). The present work describes step-by-step the Monte Carlo algorithm with focus on modeling the polarization properties of the scattered photon. The polarization is treated independently both in the Stokes parameters and in the polarization vector formalisms. An intervalidation between the two methods is given. Based on polarization state description requirements of different Geant4 physics classes, user recommendations are given on which of the two methods to be employed. The spatial and energy distributions for the LCS γ-ray beam and its Stokes parameters are obtained for head-on laser — relativistic electron collisions, where several incident laser polarization states were considered: linear, unpolarized, circular  and mixed linear and circular polarization. Results of previous investigations on the polarization of Compton scattered photons are reproduced. The influence of variable incident angle between photon and electron beam was also investigated. We show that the degree of polarization transfer from the incident photon to the scattered photon increases with the collision angle, where head-on is considered 0°. However, as the polarization transfer is strongly influenced by the incident photon energy, we show that, for γ-ray sources based on Compton scattering of laser photons on relativistic electrons, the polarization degree of the incident photon is almost completely transferred to the scattered photon for any incident angle.}
}
```

```
@article{FILIPESCU2023167885,
title = {Spectral distribution and flux of γ-ray beams produced through Compton scattering of unsynchronized laser and electron beams},
journal = {Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment},
volume = {1047},
pages = {167885},
year = {2023},
issn = {0168-9002},
doi = {https://doi.org/10.1016/j.nima.2022.167885},
url = {https://www.sciencedirect.com/science/article/pii/S0168900222011779},
author = {Dan Filipescu and Ioana Gheorghe and Konstantin Stopani and Sergey Belyshev and Satoshi Hashimoto and Shuji Miyamoto and Hiroaki Utsunomiya},
keywords = {Laser Compton scattering, Gamma ray source, Monte Carlo simulation, Photoneutron reactions, Electron synchrotron, LaBr:Ce detector},
abstract = {Intense, quasi-monochromatic, polarized γ-ray beams with high and tunable energy produced by Compton scattering of laser photons against relativistic electrons are used for fundamental studies and applications. Following a series of photoneutron cross section measurements in the Giant Dipole Resonance (GDR) energy region performed at the NewSUBARU synchrotron radiation facility, we have developed the eliLaBr Monte Carlo simulation code for characterization of the scattered γ-ray photon beams. The code is implemented using Geant4 and is available on the GitHub repository (https://github.com/dan-mihai-filipescu/eliLaBr). Here we report the validation of the eliLaBr code on NewSUBARU LCS γ-ray beam flux and spectral distribution data and two applications performed with it for asymmetric transverse emittance profiles electron beams, characteristic for synchrotrons. The first application is based on a systematic investigation of transverse collimator offsets relative to the laser and electron beam axis. We show that the maximum energy of the LCS γ-ray beam is altered by vertical collimator offsets, where the edge shifts towards lower energies with the increase in the offset. Secondly, using the eliLaBr code, we investigate the effect of the laser polarization plane orientation on the properties of the LCS γ-ray beams produced with asymmetric emittance electron beams. We show that: 1. The use of vertically polarized lasers contributes to the preservation of the LCS γ-ray beam maximum energy edge by increasing the precision in the vertical collimator alignment. 2. Under identical conditions for the electron and laser beams phase-space distributions, the energy spectrum of the scattered LCS γ-ray beam changes with the laser beam polarization plane orientation. More precisely, the use of vertically polarized laser beams slightly deteriorates the LCS γ-ray beam energy resolution.}
}
```

## Introduction

**eliLaBr** was created starting with 2013 during first photoneutron cross section measurements performed at [Laser Compton Scattering (LCS) Gamma-ray beam](https://www.lasti.u-hyogo.ac.jp/NS-en/facility/bl01/), beamline BL01/GACKO Hutch from [NewSUBARU](http://www.spring8.or.jp/en/about_us/whats_sp8/facilities/accelerators/new_subaru/) synchrotron, SPring8, JAPAN.

The code simulate Compton interaction between laser and relativistic electron beams. Laser is modelled using Gaussian beams, while the electron beam is described using Twiss parameter formalism. This code was the main tool used in conceiving ELI-GANT TDR (**E**xtreme **L**ight **I**nfrastructure - **G**amma **A**bove **N**eutron **T**hreshold - **T**echnical **D**esign **R**eport). More details are given in [ApendixC1](https://raw.githubusercontent.com/dan-mihai-filipescu/eliLaBr/main/doc/ELI_TDR3-GANT_ApendixC1_V1.9.4_April_2015.pdf) of ELI-GANT TDR. One of the main classes of the code, **`eliLaBr_GammaSource`**, provides the required parameters to the **`PrimaryGenerator`** class of **GEANT4** framework. Radiation is further transported through collimators in order to produce quasi-monochromatic gamma-ray beam which is incident on a target. The reaction products are detected with radiation detectors. Information regarding involved **GEANT4** physics classes are given in [ApendixC3](https://raw.githubusercontent.com/dan-mihai-filipescu/eliLaBr/main/doc/ELI_TDR3-GANT_ApendixC3_V1.9.4_April_2015.pdf) of GANT TDR.

This code was also the main tool used to design the flat efficiency neutron detector consisting from 31 <sup>3</sup>He counters (10 atm.) embedded into polyethylene moderator block. The present code was used in performing intensive simulations, testing many configurations, geometries, number of counters, reaching the most **Flat Efficiency** configuration which is currently used at the NewSUBARU Gamma-ray beam (gamma, xn) cross section measurements.

The same code was used to make the second design of flat efficiency neutron detector for ELI-NP, namely **GANT-TN** instrument which consists from 28 <sup>3</sup>He counters (12 atm.) which is now under operation at ELI-NP. More simulated results of **GANT-TN** instrument can be found in [ApendixC4](https://raw.githubusercontent.com/dan-mihai-filipescu/eliLaBr/main/doc/ELI_TDR3-GANT_ApendixC4_V1.9.4_April_2015.pdf) and [ApendixD](https://raw.githubusercontent.com/dan-mihai-filipescu/eliLaBr/main/doc/ELI_TDR3-GANT_ApendixD_V1.9.4_April_2015.pdf) of GANT TDR.

Another version of the same code (not uploaded yet on GitHub, to be uploaded in the future) was also used to design and make reaction rate estimates for the second GANT instrument: **GANT-GN**. [**GANT-GN**](http://www.eli-np.ro/rd4.php) (see ELIGANT tab) consists from 15 LaBr3(Ce) and 25 CeBr3 3x3 inch scintillator detectors for gamma-ray detection, and 35 liquid (20x5 cm) and 24 Li-glass (10x2 cm) scintillator detectors for fast-neutron detection with a variable flight-path between 0.5 m and 1.5 m. More details on: [ApendixE](https://raw.githubusercontent.com/dan-mihai-filipescu/eliLaBr/main/doc/ELI_TDR3-GANT_ApendixE_V1.9.4_April_2015.pdf) of GANT TDR.

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
  - when we are interested to run **eliLaBr** in *interactive mode* we should do the other way around:
  ```
  /control/execute gui.mac
  /control/execute vis1.mac
  #/run/beamOn #(int)NumberOfEvents
  ```
  
  **Important note**: the **eliLaBr** code generate the initial seeds of the Monte Carlo simulation based on the computer clock information collected when simulation starts, thus the simulation results cannot be identically reproduced, unless one input the same seeds at the beginning of the simulation (the seeds are displayed on screen at the beginning). This feature is important when using on multi-core computers - when launching **eliLaBr** multiple times using the same input, each process will run on a different core, and of course it will be launched at a slightly different time, thus the results of each process will be statistically different, allowing the user to add the simulation results of each simulation on different processor core into a single sum result.
  
  ## 1 - GEOMETRY DEFINITION
  
  The general layout of the experimental setup follows the NewSUBARU experimental facility and all characteristic sizes, distances and materials: beamline BL01, synchrotron dipole magnetic field, laser optical mirror used to insert the laser beam into the electron beamline, borosilicate vacuum window, copper gamma absorber installed into the electron beamline, gamma-beam shutter, primary and secondary gamma-beam collimators, and of course, flat efficiency neutron detector consisting from <sup>3</sup>He counters  embedded into the polyethylene moderator block. However, the experimental hall and the gamma-beam dump at the end of it, follows the ELI-NP E8 experimental hall design. This choice of experimental hall does not affect at all the simulation results in the case of <sup>3</sup>He counters 4Pi flat efficiency neutron detector. This have more influence in the case of fast-neutron detection array based on time of flight, due to the scattered neutrons on the experimental hall walls, and this justified the choice of sticking to the E8 geometry.

![NewSUBARU layout](/images/NewSUBARU.png)
  === NewSUBARU layout. ===
  
![4Pi Flat Efficiency neutron detector](/images/4Pi_Detector.png)
  === 4Pi Flat Efficiency neutron detector ===
  
  Several geometry macro-commands are implemented in `eliLaBr_DetectorMessenger` class.
  
  A list with these macro-commands can be found in **`SetGeometry.mac`** macrofile which should be loaded into **`commands`** macrofile.
  
  **TARGET** is made *sensitive* and can be used as a detector. Three example cases are listed into **`SetGeometry.mac`**,and the corresponding block of macro-commands can be un-commented depending on the needs:
  - *TARGET INSIDE 4Pi detector* - this is the normal case of a target placed in the middle of the neutron detector. Be careful on the target dimensions not to exceed the inner diameter of the beamline pipe passing through the neutron detector central axis.
  - *TARGET INSIDE E8 Hall* - this is intended to monitor the beam spectrum and usually used to simulate a beam monitor detector. Do not forget to un-comment *TARGET HOUSING* section macro-command `/eli/det/setTGHousing true` in case you intend to.
  - *MONITOR FOR GAMMA SOURCE* - this is intended to be an imaginary plate placed in front of the gamma-beam prior to collimation, in order to histogram different characteristics of the gamma-beam like spatial, energy and/or polarization distributions.
  
  **Neutron counters** dimensions are defined by those of the detectors existing at NewSUBARU & ELI-NP facilities, namely: 1 inch diameter and 500 mm. length, while the detector housing is considered to be made from stainless-steel.
  
  Two types of sensitive gas are allowed:
  - *BF3 GAS*
  - *3He GAS*
  
  so please un-comment the corresponding block of macro-commands.
  
  There are allowed a maximum number of **three** concentric rings of neutron counters around the central beampipe passing through the axis of the polyethylene moderator cube.
  
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
  
  In the root directory of the project, the user can find some example ROOT macros (read2.cc, read2_sum.cc, read_BeamSpot.cc, read_BeamSpot_sum.cc, etc.), that can be further extended and adapted in order to read and analyse the 1D & 2D histograms incremented during single- OR multiple-core simulations.
  
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
  
  The information from *Hits Collection* is further analyzed and processed in *End Of Event Action* method from the **`EventAction`** class. Thus the energy deposition in time is analyzed for each neutron counter separately and and also for target. The **ROOT** class **`TSpectrum`** is used to identify separate peak energy deposited over time. These chained vectors of pairs (DepositedEnergy, Time) are stored in a **`TVectors`** class, specially written for this task. The information written in such vectors is possible to be analyzed later, considering different time integration constants, in a similar way as the *Shaping time* is acting in real electronics, in order to be able to treat energy deposits with enough distance in time between them as different events. For example, if in the same GEANT event, after a first neutron is detected in one counter, if a second neutron is detected by the same counter at an enough time distance from the first neutron (due to longer thermalization time), we can consider two different neutron detection events for the same counter.
  
  In *End Of Event Action* method the energy deposited is integrated on each counter and in the target, the average time of energy deposited is computed (if we had EDep>0) for each counter and for target, and the multiplicity of detectors along with their index number is determined in order to put all these information into vectors.
  
  Also we check if any radiation passes the boundary of a detector/target (towards detector/target), without being interested if any energy was deposited or not, in that particular detector during the entire event, and we store the multiplicity of all detectors which were crossed by any radiation, the detector index number, what kind of radiation was, and with what incident energy.
  
  All information obtained in *End Of Event Action* method by processing the *Hits Collection* is written into a **ROOT** file, event by event, in order to be able to sort the events later, depending on the desired conditions, according to different imagined scenarios. Unfortunately, this event file is always named **`nai_9_MeV.root`**, but this limitation will be eliminated in a future release of **eliLaBr** code.
  
  ## Tvectors library
  
  **`Tvectors`** dynamic library provides the classes dealing with storage of vectors of pairs (DepositedEnergy, Time) for each detector, the vector of which elements consists from previous vectors for all detectors involved in the current event and also the necessary methods to access the vector element pairs, insert/append/delete/initialize/sort/contract elements, etc. The class is needed while performing the GEANT simulation, but also in the process of sorting/analyzing the simulation results.
  
  ## DataAnalysis package files
  
  This folder contains **`PileItUp`** code which is used to perform the sorting of events saved in the **`nai_9_MeV.root`** event file produced by the **eliLaBr** GEANT4 simulation code. The code relies on **`PileItUp.in`** input file which contains self-explanatory information.
  
  One important option specify how many consecutive events should be packed together in order to simulate the pile-up effect on the detector due to the multiple particle within one bunch. Thus, with one unique simulation performed, we can make multiple data analysis considering different scenarios: no pile-up, or pile-up considering different distributions, averages, etc. We can consider a fixed integer number of pile-up events, or a Poisson distribution of events, providing the average value of the Poisson distribution. **`Tvectors`** class is used in the pile-up process, appending (DepositedEnergy, Time) pair elements to the vectors of each detector, then sorting the elements and contracting them according to a *Shaping Time* packing provided at input in order to simulate the experimental effect of a shaper.
  
  The user can request analysis code to make output Energy/Time spectra files for individual <sup>3</sup>He counters AND/OR for 4Pi n-detector concentric rings AND totals if ANY/OR CERTAIN particle (neutron, gamma, proton, e+, e-, deuteron, alpha, ,<sup>12</sup>C, etc.) **deposited energy** OR **only crossed the counter border**. All these options better help understanding how different processes are taking place and unveil the details of the interaction.
  
  The most important task performed by the **`PileItUp`** code is the conversion of the simulated results into **GASP** list data format, which is the current list data format adopted for the acquired experimental data in the experiments performed at **RoSPHERE** array from 9MV Tandem, IFIN-HH, Romania and **GASP** array from INFN Legnaro, Italy laboratories. In this way, the user will be able to apply the same data analysis procedure to the simulated results as it would do for data acquired in a real experiment, using **GASPWare** data analysis package.
  
  In fact, **`PileItUp`** code converts the simulated list data results in two **GASP** files, using two configurations:
  1. In the first configuration it is considered that we have only one type of detector, namely the <sup>3</sup>He counter type, while the total number of such detectors is basically the total number of physical detectors (31 in the case of NewSUBARU flat efficiency detector). For each detector two parameters are stored (Energy and Time), at each event (piling-up or not) per target. Only one (the **First** one) Energy & Time per event and per detector are allowed. This is the common event list data adopted at **RoSPHERE** and **GASP** arrays. For this case, the user can check the demo GASP setup file **`Proje_indv.setup`**.
  2. In the second configuration it is considered that we have four types of detectors:
     - inner ring type of detectors (A);
     - middle ring type of detectors (B);
     - outer ring type of detectors (C);
     - all detectors (disregarding the ring);

  For each type of detector, the total number of detectors correspond to the total number of neutron hits per that type, considering that the vectors with (Energy, Time) pair hits corresponding to that detector was previously sorted according to time increase and contracted according to *Shaping Time*. In this second scenario, an event correspond to an entire time interval between two beam bunches, needed to allow thermalization and detection of all neutrons emitted in the reaction. This second configuration is the usual one used in the photo-neutron (g,xn) cross section measurements performed at the NewSUBARU gamma-beam facility and at the (particle, xn) cross section measurements performed at Tandem accelerators from IFIN-HH. For this case, the user can use the demo GASP setup file **`Proje.setup`**.
  
  User should provide to **`PileItUp`** code also the file **`files`** which contains the the total number of events that the user want to analyze and the names of the ROOT event files that have to be analyzed (this is useful if many such event files were produced in multiple-core simulations).
  
  ## Papers
  
  Several papers made use of **eliLaBr** code over the time:
  - D. Filipescu et al. *"Spectral distribution and flux of g-ray beams produced through Compton scattering of unsynchronized laser and electron beams"*, NIMA **1047** 167885 (2023) [10.1016/j.nima.2022.167885](https://doi.org/10.1016/j.nima.2022.167885) (arXiv: [arXiv:2211.14650]( 	
https://doi.org/10.48550/arXiv.2211.14650))
  - D. Filipescu *"Monte Carlo simulation method of polarization effects in Laser Compton Scattering on relativistic electrons"*, JINST **17** P11006 (2022) [10.1088/1748-0221/17/11/P11006](https://doi.org/10.1088/1748-0221/17/11/P11006) (arXiv: [arXiv:2210.14669](https://doi.org/10.48550/arXiv.2210.14669))
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
  - W. Luo et al. *"PRODUCTION OF RADIOISOTOPES OF MEDICAL INTEREST BY PHOTONUCLEAR REACTION USING ELI-NP gamma-RAY BEAM"*, ACTA PHYSICA POLONICA B **47** 763 (2016) [10.5506/APhysPolB.47.763](https://doi.org/10.5506/APhysPolB.47.763)
  - W. Luo et al. *"A data-based photonuclear simulation algorithm for determining specific activity of medical radioisotopes"*, Nuclear Science and Techniques **27** 113 (2016) [10.1007/s41365-016-0111-9](https://doi.org/10.1007/s41365-016-0111-9)
  - F. Camera et al. *"GAMMA ABOVE THE NEUTRON THRESHOLD EXPERIMENTS AT ELI-NP"*, Romanian Reports in Physics **68** S539 (2016) [RRP68_S539(2016)](http://www.rrp.infim.ro/2016_68_S.html)
  - D. Filipescu et al. *"Perspectives for photonuclear research at the Extreme Light Infrastructure - Nuclear Physics (ELI-NP) facility"*, EPJA **51** 185 (2015) [10.1140/epja/i2015-15185-9](https://doi.org/10.1140/epja/i2015-15185-9)
  - T. Renstrom et al. *"First evidence of low energy enhancement in Ge isotopes"*, EPJ Web of Conf. **93** 04003 (2015) [10.1051/epjconf/20159304003](https://doi.org/10.1051/epjconf/20159304003)
  - D. Filipescu et al. *"Photoneutron cross section measurements on Sm isotopes"*, EPJ Web of Conf. **93** 02006 (2015) [10.1051/epjconf/20159302006](https://doi.org/10.1051/epjconf/20159302006)
  - H. Utsunomiya et al. *"Photoneutron Reactions in Nuclear Astrophysics"*, Journal of Physics Conference Series **590** 012023 (2015) [10.1088/1742-6596/590/1/012023](https://doi.org/10.1088/1742-6596/590/1/012023)
  - H. Utsunomiya et al. *"Photodisintegration of 9Be through the 1/2+ state and cluster dipole resonance"*, **92** 064323 (2015) [10.1103/PhysRevC.92.064323](https://doi.org/10.1103/PhysRevC.92.064323)
  - H.T. Nyhus et al. *"Photoneutron cross sections for neodymium isotopes: Toward a unified understanding of (gamma, n) and (n, gamma) reactions in the rare earth region"*, Phys. Rev. C **91** 015808 (2015) [10.1103/PhysRevC.91.015808](https://doi.org/10.1103/PhysRevC.91.015808)
  - D. Filipescu et al. *"Geant4 simulations on Compton scattering of laser photons on relativistic electrons"*, AIP Conference Proceedings **1645** 322 (2015) [10.1063/1.4909594](https://doi.org/10.1063/1.4909594)
  - I. Gheorghe et al. *"Absolute photoneutron cross sections of Sm isotopes"*, AIP Conference Proceedings **1645** 327 (2015) [10.1063/1.4909595](https://doi.org/10.1063/1.4909595)
  - D. Filipescu et al. *"Photoneutron cross sections for samarium isotopes: Toward a unified understanding of (gamma, n) and (n,gamma) reactions in the rare earth region"*, Phys. Rev. C **90** 064616 (2014) [10.1103/PhysRevC.90.064616](https://doi.org/10.1103/PhysRevC.90.064616)
  - H. Utsunomiya et al. *"Energy Calibration of the NewSUBARU Storage Ring for Laser Compton-Scattering Gamma Rays and Applications"*, IEEE Transactions on Nuclear Science **61** 1252 (2014) [10.1109/TNS.2014.2312323](https://doi.org/10.1109/TNS.2014.2312323)
  
