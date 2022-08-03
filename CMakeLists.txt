#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(Detector4Pi_new)
message("Adapted by Flippy for eliLaBr on 18.10.2020")

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find ROOT (optional package)
#
find_package(ROOT)
if(ROOT_FOUND)
  add_definitions(-DG4ANALYSIS_USE_ROOT)
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/Tvectors
                    ${Geant4_INCLUDE_DIR}
		    ${ROOT_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/Tvectors/*.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(eliLaBr eliLaBr.cc ${sources} ${headers})
target_link_libraries(eliLaBr ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Spectrum -L../Tvectors Tvectors)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build eliLaBr01. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(eliLaBr_SCRIPTS
    commands histos.mac SetPhysics.mac SetGeometry.mac incident_energy.in incident_gamma.in NewSUBARU_optics_BL01.txt vrml.mac vis.mac vis1.mac gui.mac
  )

foreach(_script ${eliLaBr_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS eliLaBr DESTINATION bin)

