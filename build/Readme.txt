
Instructions to build "eleLaBr" code locally.
=============================================

First build the needed dynamic libraries individually:
Go to "Tvectors" folder and type "make"
Go to "DataAnalysis" folder and type "make"

Go to "build" folder.
Type "cmake -DCMAKE_PREFIX_PATH=/home/geant/GEANT/geant4-install ../"
then type "make"

If everything works fine you should see the "eliLaBr" executable inside "build" folder.
Run the code within "build" folder.
