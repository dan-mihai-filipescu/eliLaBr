# Makefile for the pileitup Data Analysis program using TSelector.
# This Makefile compile and link PileItUp.C, PileItUp.h & runPileit.cpp
# using the ROOT libraries on Linux platforms.
#
# Copyright (c) 2017 Dan Mihai FILIPESCU
#
# Author: Dan Mihai FILIPESCU, 08/02/2017

RC     := root-config
ifeq ($(shell which $(RC) 2>&1 | sed -ne "s@.*/$(RC)@$(RC)@p"),$(RC))
MKARCH := $(wildcard $(shell $(RC) --etcdir)/Makefile.arch)
RCONFIG := $(wildcard $(shell $(RC) --incdir)/RConfigure.h)
endif
ifneq ($(MKARCH),)
include $(MKARCH)
else
ifeq ($(ROOTSYS),)
ROOTSYS = ..
endif
include $(ROOTSYS)/etc/Makefile.arch
endif
-include ../MyConfig.mk

#------------------------------------------------------------------------------

TVECTORSDIR	= ../Tvectors/

NAME           = PileItUp
MAINNAME       = pileitup
TPILEITUPO     = $(NAME).$(ObjSuf) $(NAME)Dict.$(ObjSuf)
TPILEITUPS     = $(NAME).$(SrcSuf) $(NAME)Dict.$(SrcSuf)
TPILEITUPSO    = lib$(NAME).$(DllSuf)

PILEITUPO      = $(MAINNAME).$(ObjSuf)
PILEITUPS      = $(MAINNAME).$(SrcSuf)
PILEITUP       = $(MAINNAME)$(ExeSuf)

OBJS          = $(TPILEITUPO) $(PILEITUPO)

PROGRAMS      = $(TPILEITUPSO) $(PILEITUP)

#------------------------------------------------------------------------------

EXTRAINC = -I$(TVECTORSDIR)
LIBS += -L$(TVECTORSDIR) -lTvectors

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(ExeSuf)

all:            $(PROGRAMS)

$(PILEITUP):    $(PILEITUPO) $(TPILEITUPSO)
		$(LD) $(LDFLAGS) $(PILEITUPO) $(TPILEITUPO) $(LIBS) $(OutPutOpt)$@
		$(MT_EXE)
		@echo "$@ done"

$(TPILEITUPSO): $(TPILEITUPO)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(EXPLLINKLIBS) $(LIBS) $(OutPutOpt)$@
		$(MT_DLL)
		@echo "Generated dynamic library for " $(ARCH)

clean:
		@rm -f $(OBJS) core *Dict.*

distclean:      clean
		@rm -f $(PROGRAMS) *Dict.* *.def *.exp *.pcm\
		   *.ps *.so *.lib *.dll *.d *.log .def so_locations \
		   files/*
		@rm -rf cxx_repository

.SUFFIXES: .$(SrcSuf)

###
$(NAME).$(ObjSuf): $(NAME).h
$(NAME)Dict.$(SrcSuf): $(NAME).h $(NAME)LinkDef.h
	@echo "Generating dictionary $@..."
	rootcint -f $@ -c $(EXTRAINC) $^

.$(SrcSuf).$(ObjSuf):
	@echo "Compiling $@..."
	$(CXX) $(CXXFLAGS) $(EXTRAINC) -c $<
