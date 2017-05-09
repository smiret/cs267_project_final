#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

HOME = ..
VPATH = . $(HOME)/include/VisitWriter
CFLAGS = -I. -I$(HOME) -I$(HOME)/src -I$(HOME)/include/VisitWriter -std=c++11
EIGEN= -I /home/lbluque/SrcCode/eigen
OPENMP = -fopenmp
DIM=3
system := $(shell uname)
CXX=g++
#CXX=clang++
FFTWDIR = /usr
FFTINCFLAGS = -I$(FFTWDIR)/include
#CFLAGS += -g -Wall
CFLAGS += -O3
#CFLAGS += -std=c++11
ifeq ($(CXX) , clang++)
  CFLAGS += -stdlib=libc++
endif
osuffix:=$(DIM)D.o
dsuffix:=$(DIM)D.d
CPPFLAGS = -D DIM=$(DIM) 

#MYOBJS = FEPoissonOperator.o FEGrid.o Element.o Node.o SparseMatrix.o

#SRCFILES:= $(wildcard $(HOME)/src/*.cpp)
SRCFILES:= ../src/blmintegrals.cpp ../src/hsbounds.cpp ../src/mesher.cpp ../src/femfunctions.cpp ../src/GridWrite.cpp
SRCFILES += $(HOME)/include/VisitWriter/VisitWriter.cpp
OBJS:=$(patsubst %.cpp,%.o, $(SRCFILES)) 

all: project

%.o:%.cpp GNUmakefile
	$(CXX) -c $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

project:  blm_main.cpp $(HOME)/include/VisitWriter/VisitWriter.o $(OBJS)
	  $(CXX) $(CPPFLAGS) $(CFLAGS) $(OPENMP) $(EIGEN) -o parallel.exe blm_main.cpp $(OBJS)

clean:
	rm -r *.o *.d *.exe *.dSYM $(HOME)/src/*.o $(HOME)/src/*.d $(HOME)/src/*.exe $(HOME)/src/*.dSYM $(HOME)/include/VisitWriter/*.o $(HOME)/include/VisitWriter/*.d $(HOME)/include/VisitWriter/*.exe $(HOME)/include/VisitWriter/*.dSYM
-include $(OBJS:.o=.d)
