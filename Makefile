#################################################################################
# Makefile THERMINATOR 2                                       ver. 2011-02-01  #
# Author: Mikolaj Chojnacki (Mikolaj.Chojnacki@ifj.edu.pl)                      #
#################################################################################
#
#helpUsage: (g)make [LABEL] [PREP_DEF]
#help  [LABEL]
#help    [all]     compile THERMINATOR2
#help    doc       generate documentation with Doxygen
#help    package   compress project to tar-ball file
#help    clean     remove object files and binaries
#help    cleandoc  remove Doxygen documentation
#help    help      this screen
#help
#help  [PREP_DEF]  list of preprocesor definitions
#help    DEBUG                     debugging level;                            default: DEBUG=0
#help    RESCALE_CHANNELS          rescale decay channels;                     default: not defined
#help    DISABLE_THREE_BODY_DECAYS disable three body decays;                  default: not defined
#help    DISABLE_TWO_BODY_DECAYS   disable two body decays;                    default: not defined
#help    BACK_FLOW                 particles emitted back to the hydro region; default: not defined
#

#################################################################################
# DEFINITIONS                                                                   #
#################################################################################

# version information form include/THGlobal.h
CALM_VERSION = $(shell grep "\#define _CALM_VERSION_" $(DIR_H)THGlobal.h | sed 's|[\#,_,2_,a-z,A-Z, ]*||')

# directory structure
DIR_MAIN    = ./
DIR_CXX     = $(DIR_MAIN)src/
DIR_H       = $(DIR_MAIN)include/
DIR_OBJ     = $(DIR_MAIN)obj/
DIR_EVENTS  = $(DIR_MAIN)events/

# search paths
vpath %.h   $(DIR_H)
vpath %.cxx $(DIR_CXX)
vpath %.o   $(DIR_OBJ)

# distribution files
F_INCLUDE   = $(DIR_H)*.h
F_SOURCE    = $(DIR_CXX)*.cxx
F_INI       = $(DIR_MAIN)*.ini
F_BASH      = $(DIR_MAIN)*.sh

# file lists
# CALM
BIN_EVENTS  = calm_events
HSRC_EVENTS = Parser.cxx Configurator.cxx ParticleDB.cxx ParticleType.cxx DecayTable.cxx DecayChannel.cxx EventGenerator.cxx Event.cxx Particle.cxx ParticleCoor.cxx ParticleDecayer.cxx Crc32.cxx Vector3D.cxx CALM.cxx
SRC_EVENTS  = $(HSRC_EVENTS:%=$(DIR_CXX)%) $(BIN_EVENTS:%=$(DIR_CXX)%.cxx)
OBJ_EVENTS  = $(SRC_EVENTS:$(DIR_CXX)%.cxx=$(DIR_OBJ)%.o)

# preprocessor
PREPROCESS  = -D_CXX_VER_="\"$(shell $(CXX) --version | grep $(CXX))\"" -D_ROOT_VER_="\"$(shell root-config --version)\""
ifdef DEBUG
  PREPROCESS  := $(PREPROCESS) -D_DEBUG_LEVEL_=$(DEBUG)
else
  PREPROCESS  := $(PREPROCESS) -D_DEBUG_LEVEL_=0
endif
ifdef BACK_FLOW
  PREPROCESS  := $(PREPROCESS) -D_MODEL_LHYQUID_ONLY_BACK_FLOW_=$(BACK_FLOW)
endif
ifdef RESCALE_CHANNELS
  PREPROCESS  := $(PREPROCESS) -D_PARTICLE_DECAYER_RESCALE_CHANNELS_=$(RESCALE_CHANNELS)
endif
ifdef DISABLE_THREE_BODY_DECAYS
  PREPROCESS  := $(PREPROCESS) -D_PARTICLE_DECAYER_DISABLE_THREE_BODY_DECAYS_=$(DISABLE_THREE_BODY_DECAYS)
endif
ifdef DISABLE_TWO_BODY_DECAYS
  PREPROCESS  := $(PREPROCESS) -D_PARTICLE_DECAYER_DISABLE_TWO_BODY_DECAYS_=$(DISABLE_TWO_BODY_DECAYS)
endif

# compilation
CXX         = g++
LD          = g++
CXXFLAGS    = -O0 -g -Wno-deprecated -I $(DIR_H) $(PREPROCESS) `root-config --cflags`
LFLAGS      = -lm -lgcc -g `root-config --libs`

#################################################################################
# RULES                                                                         #
#################################################################################
all: $(BIN_EVENTS:%=$(DIR_OBJ)%) $(BIN_FEMTO:%=$(DIR_OBJ)%) $(BIN_HBTFIT:%=$(DIR_OBJ)%)
	cp $^ $(DIR_MAIN)
	echo
	echo "Type \"./calm_events\" to generate events"
	echo

$(DIR_OBJ)calm_events: $(OBJ_EVENTS)
	echo "Linking:   $@ ($(LD))"
	$(LD) $^ -o $@  $(LFLAGS) 

$(DIR_OBJ)%.o: %.cxx
	@[ -d $(DIR_OBJ) ] || mkdir -p $(DIR_OBJ)
	echo "Compiling: $< ($(CXX))"
	$(CXX) $(CXXFLAGS) -c $< -o $@

package: $(F_INCLUDE) $(F_SOURCE) $(F_MACRO) $(F_FOMODEL) $(F_SHARE) $(F_DOXYGEN) $(F_INI) $(F_BASH) $(F_ADDONS) Makefile
	echo "$(CALM_VERSION)" > version
	tar zcvf $(F_PACK) $^ version
	echo "Package '$(F_PACK)' created."

clean:
	rm -rf $(DIR_OBJ)
	rm -f $(DIR_OBJ)$(BIN_EVENTS) $(DIR_MAIN)$(BIN_EVENTS)
	rm -rf *.log
	echo "*.o and binary files removed."

help:
	@grep -h "^#help" $(MAKEFILE_LIST) | sed 's|\#help||'

.SILENT :
