#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

HOME = ../..
#VPATH = . $(HOME)/../include/VisitWriter
#CFLAGS = -I. -I$(HOME) -I$(HOME)/../include/VisitWriter -std=c++11
CFLAGS = -I. -I$(HOME) -std=c++11
DIM=3
system := $(shell uname)
CXX=g++
#CXX=clang++
FFTWDIR = /usr
FFTINCFLAGS = -I$(FFTWDIR)/include
CFLAGS += -g -Wall
#CFLAGS += -O3
#CFLAGS += -std=c++11
ifeq ($(CXX) , clang++)
  CFLAGS += -stdlib=libc++
endif
osuffix:=$(DIM)D.o
dsuffix:=$(DIM)D.d
CPPFLAGS = -D DIM=$(DIM) 

#FEOBJS = FEPoissonOperator.o FEGrid.o Element.o Node.o SparseMatrix.o

SRCFILES:= $(HOME)/GridWrite.cpp $(HOME)/mesher.cpp
OBJS:=$(patsubst %.cpp,%.o, $(SRCFILES)) 

all: meshtest

%.o:%.cpp GNUmakefile
	$(CXX) -c $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

meshtest:  meshTest.cpp $(OBJS)
	  $(CXX) $(CPPFLAGS) $(CFLAGS) -o meshtest.exe meshTest.cpp $(OBJS)

clean:
	rm -r *.o *.d *.exe *.dSYM $(HOME)/*.o $(HOME)/*.d $(HOME)/*.exe $(HOME)/*.dSYM
-include $(OBJS:.o=.d)