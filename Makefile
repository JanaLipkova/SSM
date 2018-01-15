# output binary
BIN := ssm

# source files
#SRCS := $(shell find ./source/RNGLib  -name *.cpp)  \
#        $(shell find ./source/Timer   -name *.cpp)  \
#        $(shell find ./source/Methods -name *.cpp)  \
#        $(shell find ./source/Jacobians -name *.cpp)  \
#        ./source/SBMLReaderAndParser.cpp \
#	    ./source/Simulation.cpp \
#	    ./source/SSMReaction.cpp \
#	  	./source/StochasticSimulationMethods.cpp

SRCS := $(shell find ./source/RNGLib  -name *.cpp)  \
		$(shell find ./source/Timer   -name *.cpp)  \
        $(shell find ./source/Jacobians -name *.cpp)  \
		./source/Methods/SSA.cpp \
		./source/Methods/AdaptiveSLeaping.cpp \
		./source/Methods/SLeaping.cpp \
		./source/Methods/TauLeaping.cpp \
		./source/Methods/RLeapingJana.cpp \
		./source/Methods/RootFinderJacobian.cpp \
       	./source/SBMLReaderAndParser.cpp \
	    ./source/Simulation.cpp \
	    ./source/SSMReaction.cpp \
	  	./source/StochasticSimulationMethods.cpp
#		./source/Methods/AdaptiveTau.cpp \

SBML_INC_DIR=$(HOME)/usr/include
SBML_LIB_DIR=$(HOME)/usr/lib

BLITZ_INC_DIR=$(HOME)/usr/include
BLITZ_LIB_DIR=$(HOME)/usr/lib

BOOST_INC_DIR=$(BOOST_INCLUDEDIR)
BOOST_LIB_DIR=$(BOOST_LIBRARYDIR)

GSL_INC_DIR=
GSL_LIB_DIR=


# files included in the tarball generated by 'make dist' (e.g. add LICENSE file)
DISTFILES := $(BIN)

# filename of the tar archive generated by 'make dist'
DISTOUTPUT := $(BIN).tar.gz

# intermediate directory for generated object files
OBJDIR := .o
# intermediate directory for generated dependency files
DEPDIR := .d

# object files, auto generated from source files
OBJS := $(patsubst %,$(OBJDIR)/%.o,$(basename $(SRCS)))
# dependency files, auto generated from source files
DEPS := $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRCS)))

# compilers (at least gcc and clang) don't create the subdirectories automatically
$(shell mkdir -p $(dir $(OBJS)) >/dev/null)
$(shell mkdir -p $(dir $(DEPS)) >/dev/null)

# C compiler
CC := gcc
# C++ compiler
CXX := g++
# linker
LD := g++
# tar
TAR := tar

#INC = -I$(SBML_INC_DIR) -I$(BLITZ_INC_DIR) -I$(BOOST_INC_DIR) -I$(GSL_INC_DIR)
INC = -I$(SBML_INC_DIR) -I$(BLITZ_INC_DIR) -I$(BOOST_INC_DIR)
#LDLIBS = -L$(BLITZ_LIB_DIR) -lblitz $(SBML_LIB_DIR)/libsbml.a  -lexpat -lbz2 -lz -L$(GSL_LIB_DIR) -lgsl -lgslcblas -lm
LDLIBS = -L$(BLITZ_LIB_DIR) -lblitz -L$(SBML_LIB_DIR) -lsbml  -lexpat -lbz2 -lz -lgsl -lgslcblas -lm
WARNINGS = -Wno-absolute-value  -Wno-parentheses  -Wno-unused-variable  -Wno-unused-local-typedef \
		   -Wno-c++11-compat-deprecated-writable-strings -Wno-header-guard -Wno-format


# C flags
CFLAGS :=
# C++ flags
CXXFLAGS := 
# C/C++ flags
CPPFLAGS :=  $(INC)  $(WARNINGS) 
# linker flags
LDFLAGS :=  $(LIB)
# flags required for dependency generation; passed to compilers
DEPFLAGS = -MT $@ -MD -MP -MF $(DEPDIR)/$*.Td

# compile C source files
COMPILE.c = $(CC) $(DEPFLAGS) $(CFLAGS) $(CPPFLAGS) -c -o $@
# compile C++ source files
COMPILE.cc = $(CXX) $(DEPFLAGS) $(CXXFLAGS) $(CPPFLAGS) -c -o $@
# link object files to binary
LINK.o = $(LD) $(LDFLAGS) $(LDLIBS) -o $@
# precompile step
PRECOMPILE =
# postcompile step
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d


 # Define Dimerization system
CPPFLAGS += -DDimerization


all: $(BIN)

dist: $(DISTFILES)
	$(TAR) -cvzf $(DISTOUTPUT) $^

.PHONY: clean
clean:
	$(RM) -r $(OBJDIR) $(DEPDIR)

.PHONY: distclean
distclean: clean
	$(RM) $(BIN) $(DISTOUTPUT)

.PHONY: install
install:
	@echo no install tasks configured

.PHONY: uninstall
uninstall:
	@echo no uninstall tasks configured

.PHONY: check
check:
	@echo no tests configured

.PHONY: help
help:
	@echo available targets: all dist clean distclean install uninstall check


$(BIN): $(OBJS)
	$(LINK.o) $^
#	$(LD)  $(LDFLAGS) $(LDLIBS) $^  -o $@


$(OBJDIR)/%.o: %.c
$(OBJDIR)/%.o: %.c $(DEPDIR)/%.d
	$(PRECOMPILE)
	$(COMPILE.c) $<
	$(POSTCOMPILE)


$(OBJDIR)/%.o: %.cpp
$(OBJDIR)/%.o: %.cpp $(DEPDIR)/%.d
	$(PRECOMPILE)
	$(COMPILE.cc) $<
	$(POSTCOMPILE)


$(OBJDIR)/%.o: %.cc
$(OBJDIR)/%.o: %.cc $(DEPDIR)/%.d
	$(PRECOMPILE)
	$(COMPILE.cc) $<
	$(POSTCOMPILE)


$(OBJDIR)/%.o: %.cxx
$(OBJDIR)/%.o: %.cxx $(DEPDIR)/%.d
	$(PRECOMPILE)
	$(COMPILE.cc) $<
	$(POSTCOMPILE)

.PRECIOUS = $(DEPDIR)/%.d
$(DEPDIR)/%.d: ;


-include $(DEPS)




