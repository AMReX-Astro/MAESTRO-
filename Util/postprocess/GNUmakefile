PRECISION = DOUBLE
PROFILE = FALSE
DEBUG = FALSE
DIM = 3

COMP = gnu

USE_MPI = TRUE
USE_OMP = FALSE

USE_REACT = TRUE

# programs to be compiled
ALL: radial_$(DIM)d.ex

Bpack   := ./Make.package
Blocs   := .

MAESTROEX_HOME := ../..

EOS_DIR := gamma_law_general
CONDUCTIVITY_DIR := constant
NETWORK_DIR := general_null
NETWORK_INPUTS := H_He.net

include $(MAESTROEX_HOME)/Exec/Make.Maestro

radial_$(DIM)d.ex: $(objForExecs)
	@echo Linking $@ ...
	$(SILENT) $(PRELINK) $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(libraries)
