# Makefile for the libTvectors.so dynamic library.
# This Makefile compile and link vectors.cxx & vectors.h
# using the ROOT libraries on Linux platforms.
#
# Copyright (c) 2016 Dan Mihai FILIPESCU
#
# Author: Dan Mihai FILIPESCU, 09/12/2016

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

NAME          = vectors
VECTORSO      = $(NAME).$(ObjSuf) $(NAME)Dict.$(ObjSuf)
VECTORSS      = $(NAME).$(SrcSuf) $(NAME)Dict.$(SrcSuf)
VECTORSSO     = libT$(NAME).$(DllSuf)

OBJS          = $(VECTORSO)

PROGRAMS      = $(VECTORSSO)

#------------------------------------------------------------------------------

#.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.SUFFIXES: .$(SrcSuf) .$(ObjSuf)

all:            $(PROGRAMS)

$(VECTORSSO):    $(VECTORSO)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(EXPLLINKLIBS) $(OutPutOpt)$@
		$(MT_DLL)
		@cp $@ ~/lib
		@echo "Generated dynamic library for " $(ARCH)

clean:
		@rm -f $(OBJS) core *Dict.*

distclean:      clean
		@rm -f $(PROGRAMS) *Dict.* *.def *.exp \
		   *.root *.ps *.so *.lib *.dll *.d *.log .def so_locations \
		   files/*
		@rm -rf cxx_repository

.SUFFIXES: .$(SrcSuf)

###
$(NAME)$(ObjSuf): $(NAME).h
$(NAME)Dict.$(SrcSuf): $(NAME).h $(NAME)LinkDef.h
	@echo "Generating dictionary $@..."
	rootcint -f $@ -c $^

.$(SrcSuf).$(ObjSuf):
	$(CXX)  $(CXXFLAGS) -c $<
