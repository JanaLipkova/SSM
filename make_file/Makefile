SHELL := /bin/bash

hn ?= $(shell hostname)

ifneq "$(findstring mpp,$(shell hostname))" ""
CPPFLAGS += -D LRZ_CLUSTER
include make.lrz
endif

ifneq "$(findstring garampat,$(shell whoami))" ""
CPPFLAGS += -D LRZ_GEORGE
include make.george
endif

ifeq "$(config)" "production"
CPPFLAGS += -DNDEBUG
CPPFLAGS += -D_FMMSILENT
endif

# Define Dimerization system
CPPFLAGS += -D Dimerization

LIBS += $(SBML_LIB_DIR)/libsbml.a 
LIBS += -L$(BLITZ_LIB_DIR) -lblitz -lexpat -lbz2 -lz 
LIBS += -L$(GSL_LIB_DIR) -lgsl -lgslcblas -lm

LD=$(CC)
.DEFAULT_GOAL := ssm

########################################################################
# List of all programs to be compiled
########################################################################
SRCS = $(wildcard  ../source/RNGLib/*.cpp)  \
        $(wildcard ../source/Timer/*.cpp)  \
        $(wildcard ../source/Methods/*.cpp)  \
        $(wildcard ../source/Jacobians/*.cpp)  \
        ../source/SBMLReaderAndParser.cpp \
        ../source/Simulation.cpp \
        ../source/SSMReaction.cpp \
        ../source/StochasticSimulationMethods.cpp

VPATH := ../source/ ../source/Timer/ ../source/Methods/ ../source/Jacobians/ ../source/RNGLib/

# object files, auto generated from source files
OBJECTS=$(SRCS:.cpp=.o)

# compile C source files
ssm: $(OBJECTS)
	$(CC) $^ -o $@ $(LIBS)
	@echo done


%.o: %.cpp
	$(CC) $(CPPFLAGS) -c $^ -o $@


clean:
	rm -f *.o
	rm -f ssm
	rm -f ../source/Timer/*.o
	rm -f ../source/*.o
	rm -r ../source/Methods/*.o
	rm -r ../source/Jacobians/*.o
	rm -r ../source/RNGLib/*.o	


