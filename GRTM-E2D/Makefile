#===============================================================================
# Author: Tche L., USTC, seistche@gmail.com
# Created at: Wed 30 Dec 2020 09:32:34 PM CST
#-------------------------------------------------------------------------------

FC := gfortran
CLFS := -fopenmp

FFTWPATH := /opt/fftw3-gnu

CFLAGS :=
LFLAGS :=

SRCDIR := src
OBJDIR := obj
BINDIR := bin

STATIC :=
DEBUG :=
DOUBLETYPE := ON
FFTW := ON
IEEE := ON
PeiDHS := ON
SH := ON
PS := ON
COLORPRINT := ON
PROGBAR := ON

DFLAG_LIST := DEBUG DOUBLETYPE FFTW IEEE PeiDHS SH PS COLORPRINT PROGBAR

DFLAGS := $(foreach f, $(DFLAG_LIST), $(if $($(f)), -D$(f), ))
FCEXE := $(notdir $(FC))

CFLAGS += $(DFLAGS) -I$(OBJDIR)

ifeq ($(FCEXE), gfortran)
  CFLAGS += -J$(OBJDIR)
else ifeq ($(FCEXE), ifort)
  CFLAGS += -module $(OBJDIR)
else
  $(error "Unrecognized compiler <"$(FCEXE)">")
endif

ifdef DOUBLETYPE
  FFTWLIB = fftw3
else
  FFTWLIB = fftw3f
endif

ifdef STATIC
  CLFS += -fPIC
  ifeq ($(FCEXE), gfortran)
    LFLAGS += -static-libgfortran
  else ifeq ($(FCEXE), ifort)
    LFLAGS += -static-intel -qopenmp-link=static
  endif
endif

ifdef DEBUG
  CLFS += -O0
  CFLAGS += -g
  ifeq ($(FCEXE), gfortran)
	CLFS += -Wall -fcheck=all #-ffpe-trap=invalid,zero,overflow
  else ifeq ($(FCEXE), ifort)
	CLFS += -warn all -check all #-fpe:0
  endif
else
  CLFS += -O2
endif

ifdef FFTW
  CFLAGS += -I$(FFTWPATH)/include
  ifdef STATIC
    LFLAGS += $(FFTWPATH)/lib/lib$(FFTWLIB).a
  else
	LFLAGS += -L$(FFTWPATH)/lib -l$(FFTWLIB)
  endif
endif

CFLAGS := $(strip $(CFLAGS))
LFLAGS := $(strip $(LFLAGS))

SRC_MOD := common.F90 math.F90 parameter.F90 grtc.F90 dwim.F90
SRC_GRT := main.F90
EXE_GRT := grtcsgram

OBJ_MOD := $(foreach f, $(SRC_MOD), $(OBJDIR)/$(f:.F90=.o))
OBJ_GRT := $(OBJDIR)/$(SRC_GRT:.F90=.o)

all : prep link

prep:
	@mkdir -p $(BINDIR) $(OBJDIR)

link: $(BINDIR)/$(EXE_GRT)

$(BINDIR)/$(EXE_GRT): $(OBJ_MOD) $(OBJ_GRT)
	$(FC) $(CLFS) $^ $(LFLAGS) -o $@

clear: clean
	-rm -f $(BINDIR)/*

clean:
	-rm -f $(OBJDIR)/*

vpath %.F90 $(SRCDIR)
.SUFFIXES: .F90 .o

$(OBJDIR)/%.o : %.F90
	$(FC) $(CLFS) $(CFLAGS) $(SRCDIR)/$(<F) -c -o $@

# vim:ft=make noet
