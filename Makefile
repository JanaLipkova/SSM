SBML_INC_DIR=/usr/local/include
SBML_LIB_DIR=/usr/local/lib

BLITZ_INC_DIR=/usr/local/include/
BLITZ_LIB_DIR=/usr/local/lib

BOOST_INC_DIR=/Users/garampat/Desktop/boost_1_65_1
BOOST_LIB_DIR=/Users/garampat/Desktop/boost_1_65_1/stage/lib

GSL_INC_DIR=/usr/local/include
GSL_LIB_DIR=/usr/local/lib

#========================================================================

CC=g++
CFLAGS = -w -c
LDFLAGS =


INC = -I$(SBML_INC_DIR) -I$(BLITZ_INC_DIR) -I$(BOOST_INC_DIR) -I$(GSL_INC_DIR)
LIB = -L$(BLITZ_LIB_DIR) -lblitz $(SBML_LIB_DIR)/libsbml.a  -lexpat -lbz2 -lz -L$(GSL_LIB_DIR) -lgsl -lgslcblas -lm


SRC = $(shell find ./source/RNGLib  -name *.cpp)  \
      $(shell find ./source/Timer   -name *.cpp)  \
      $(shell find ./source/Methods -name *.cpp)  \
      ./source/SBMLReaderAndParser.cpp \
	  ./source/Simulation.cpp \
	  ./source/SSMReaction.cpp \
	  ./source/StochasticSimulationMethods.cpp

OBJ = $(SRC:%.cpp=%.o)

BIN = StochasticSimulationMethods



all: $(SRC) $(BIN)

$(BIN): $(OBJ)
	$(CC)  $(LIB) $(OBJ) -o $@
	

.cpp.o:
	$(CC) $(CFLAGS) $(INC)  $< -o $@


clean:
	rm $(OBJ)
