# eliLaBr
GEANT4 LCS simulation on flat eff. n-detector based on 3He counters embedded into polyethylene moderator
========================================================================================================

## 0 - Usage

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
