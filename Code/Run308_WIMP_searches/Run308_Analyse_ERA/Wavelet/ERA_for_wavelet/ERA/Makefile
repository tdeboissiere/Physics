
# Makefile for the EdwRootAna package
# PLATFORM DEPENDANT FLAGS : Lyon-SL5-64bit or Mac-Intel
# NB:  Former Lyon SL4 : compiler=linker="g++4"
OSNAME          = $(shell uname -s)
ifeq ($(OSNAME),Linux)
   CC = g++
   LINKER = g++
   LINKERFLAGS = -fPIC -shared
   MACFLAG = -fPIC -O2
else
ifeq ($(OSNAME),Darwin)
   CC = g++
   LINKER = libtool
#   LINKERFLAGS = -dynamic -undefined dynamic_lookup
   MACVERSION = $(shell sw_vers -productVersion)
   LINKERFLAGS = -dynamic -undefined dynamic_lookup -macosx_version_min $(MACVERSION)
   MACFLAG = -DDARWIN
endif
endif

# Compiler options, including ROOT
#DBG     = -g
DBG     =
#PROF	= -pg
PROF	=
INC     = -I/usr/include $(shell root-config --cflags)
# new INC     = -I/opt/root/include $(shell root-config --cflags)
#LDLIBS  = -L/usr/local/lib $(shell root-config --glibs)
LDLIBS  = $(shell root-config --glibs)
CFLAGS	= $(DBG) -c -Wall $(INC) $(MACFLAG) $(PROF)
LDFLAGS	= $(DBG) $(LDLIBS) $(MACFLAG) $(PROF)

# Class files and main file
UTILS	= EdwUtils.cc
CLASSES = EdwPulse.cc EdwEvent.cc EdwTemplate.cc FitPulse.cc NoiseSpectrum.cc
EXECS 	= SambaToRoot.exe \
	  BuildNoiseSpectra.exe \
	  SimpleProcess.exe \
	  WienerProcess.exe \
	  ClassifyEvent.exe \
	  get_trigged_trace_tree.exe

# Some variable substitutions
UTILS   := $(UTILS:.cc=.o)
UTILS   := $(patsubst %,lib/%,$(UTILS))
DICT	= $(CLASSES:.cc=Dict.cc)
DICT	:= $(patsubst %,lib/%,$(DICT))
OBJECTS = $(DICT:.cc=.o)

EXECS	:= $(patsubst %,bin/%,$(EXECS))
LIBSO	:= lib/EraLib.so

.PHONY: clean

all: $(UTILS) $(DICT) $(OBJECTS) $(EXECS) $(LIBSO)

# Utilities
$(UTILS): lib/%.o : src/%.cc
	$(CC) $(CFLAGS) $^ -o $@

# Classes: First (re)build the dictionary for ROOT function
$(DICT): lib/%Dict.cc: src/%.cc src/%.h
	@rm -f $(patsubst lib/%Dict.cc,lib/LinkDef_%.h,$@)
	@echo "#ifdef __CINT__" >> $(patsubst lib/%Dict.cc,lib/LinkDef_%.h,$@)
	@echo "#pragma link off all globals;" >> $(patsubst lib/%Dict.cc,lib/LinkDef_%.h,$@)
	@echo "#pragma link off all classes;" >> $(patsubst lib/%Dict.cc,lib/LinkDef_%.h,$@)
	@echo "#pragma link off all functions;" >> $(patsubst lib/%Dict.cc,lib/LinkDef_%.h,$@)
	@echo "#pragma link C++ class "$(patsubst lib/%Dict.cc,%,$@)"+;" >> $(patsubst lib/%Dict.cc,lib/LinkDef_%.h,$@)
	@echo "#endif" >> $(patsubst lib/%Dict.cc,lib/LinkDef_%.h,$@)
	@rm -f $@ $(@:.cc=.h)
	cd lib && rootcint $(patsubst lib/%,%,$@) -c $(patsubst src/%,../src/%,$<) LinkDef_$(patsubst lib/%Dict.cc,%,$@).h && cd ..

# Class libraries
$(OBJECTS): %.o: %.cc
	$(CC) $(CFLAGS) $^ -o $@

# Executables
$(EXECS): bin/%.exe: $(UTILS) $(OBJECTS) src/%.cc
	$(CC) $(CFLAGS) $(patsubst bin/%.exe,src/%.cc,$@) -o $(patsubst bin/%.exe,lib/%.o,$@) 
	$(CC) $(LDFLAGS) $(UTILS) $(OBJECTS) $(patsubst bin/%.exe,lib/%.o,$@) -o $@ `root-config --cflags --libs`

# External library for Root command line
$(LIBSO): $(OBJECTS)
	$(LINKER) $(LINKERFLAGS) -o $@ lib/EdwUtils.o lib/EdwPulseDict.o lib/EdwEventDict.o lib/FitPulseDict.o lib/EdwTemplateDict.o lib/NoiseSpectrumDict.o

clean:
	rm -f *~ lib/*.* bin/* python/AutoDict* python/*.pyc
