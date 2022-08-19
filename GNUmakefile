# $Id: GNUmakefile,v 1.1 2016/12/10 22:02:40 GANT Exp $
# --------------------------------------------------------------
# GNUmakefile for eliLaBr.  Dan Mihai FILIPESCU, 10/12/2016.
# --------------------------------------------------------------

name := eliLaBr
G4TARGET := $(name)
G4EXLIB := true

TVECTORSDIR	= ./Tvectors/
ANALYSISDIR	= ./DataAnalysis/

#DEFINES =	-DOLD_ROOT
DEFINES =

CPPFLAGS += $(DEFINES)
CPPFLAGS += $(shell root-config --cflags)
CPPFLAGS += -I$(TVECTORSDIR)
CPPFLAGS += -Wno-deprecated-declarations
EXTRALIBS += $(shell root-config --glibs)
EXTRALIBS += -lSpectrum -L$(TVECTORSDIR) -lTvectors

ifndef G4INSTALL
  G4INSTALL = /home/cern/geant/geant4.9.0
endif

.PHONY: elilabr
#all: lib bin;
all:	vectors analysis elilabr

elilabr:	lib bin

vectors:
	@cd $(TVECTORSDIR); $(MAKE)
	@echo $(G4BINDIR)/$(G4TARGET)
	@echo $(G4LIBDIR)
#	@cp $(TVECTORSDIR)libTvectors.so $(G4LIBDIR)/
#	@cp $(TVECTORSDIR)libTvectors.so $(ANALYSISDIR)/

analysis:
	@cd $(ANALYSISDIR); $(MAKE)

allclean: clean
	cd $(TVECTORSDIR); $(MAKE) distclean
	cd $(ANALYSISDIR); $(MAKE) distclean
	@rm ~/lib/libTvectors.so

include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/binmake.gmk
