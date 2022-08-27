**Download** the **eliLaBr** code from **GitHub**:
```
git clone --recurse-submodules https://github.com/dan-mihai-filipescu/eliLaBr.git eliLaBr
```

The Geant4 simulation package **eliLaBr** is intended to build executable and run it from 
`~/geant4_workdir` directory.

## Prerequisites

The package needs Geant4-9.6 version and ROOT-5.34/ROOT-6.12 (or newer).

In order to compile **`eliLaBr`** for ROOT-5.34, the user should un-comment
```
DEFINES =   -DOLD_ROOT
```
line in `GNUmakefile`.

We are currently working to adapt the package source code in order to be able
to compile also in Geant4-10 & Geant4-11. The `TSpectrum` ROOT class is also needed to be included and linked.

Geant4 package must include all optional physics databases (`G4NDL`, `G4EMLOW`, 
`Radioactive Decay`, etc.) but especially be careful to include the 
`G4NDL4.2/ThermalScattering` library in order to be able to perform the simulation
of neutron thermalization into polyethylene by elastic scattering of neutrons 
with energies below 4 eV using S formalism.
For visualization it is recommended to build Geant4 with OpenGL option and use 
the QT user interface.

Prior building the code, it is recommended that the user create `lib` folder
into his `/home`, if not existing already, and add the following command in 
`.bashrc`:
```
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/lib
```
In this folder will be placed the libraries `libTvectors.so` and `libPileItUp.so`
that are build during the build process and will be needed in order to run the 
package.

## Make

The code can be build by typing `make` in the main directory.
This is relying on `GNUmakefile` script and will make also `libTvectors.so` and 
`libPileItUp.so` libraries.

- `make clean` will clean the last build of `eliLaBr` code
- `make allclean` will clean the last build of `eliLaBr` code but also the 
**`Tvectors`** and the **`PileItUp`** builds.

## CMake

The **eliLaBr** code can also be build locally using **CMake** based on `CMakeLists.txt` file.
Go to `build` folder and read the instructions from `Readme.txt` file to build the executable and all needed libraries there.
