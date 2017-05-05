#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

DIM = 2
HOME = ../..
VPATH = . $(HOME)/../include/VisitWriter
CFLAGS = -I. -I$(HOME) -I$(HOME)/../include/VisitWriter -std=c++11
DIM=2
system := $(shell uname)
CXX=g++
#CXX=clang++
FFTWDIR = /usr
FFTINCFLAGS = -I$(FFTWDIR)/include
#CFLAGS += -g -Wall
CFLAGS += -O3
CFLAGS += -std=c++11
ifeq ($(CXX) , clang++)
  CFLAGS += -stdlib=libc++
endif
osuffix:=$(DIM)D.o
dsuffix:=$(DIM)D.d
CPPFLAGS = -D DIM=$(DIM) 

#FEOBJS = FEPoissonOperator.o FEGrid.o Element.o Node.o SparseMatrix.o

SRCFILES:= $(wildcard $(HOME)/*.cpp)
SRCFILES += $(HOME)/../include/VisitWriter/VisitWriter.cpp
OBJS:=$(patsubst %.cpp,%.o, $(SRCFILES)) 

all: fetest

%.o:%.cpp GNUmakefile
	$(CXX) -c $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

fetest:  FEMain.cpp $(HOME)/../include/VisitWriter/VisitWriter.o $(OBJS)
	  $(CXX) $(CPPFLAGS) $(CFLAGS) -o fetest.exe FEMain.cpp $(OBJS)

clean:
	rm -r *.o *.d *.exe *.dSYM $(HOME)/*.o $(HOME)/*.d $(HOME)/*.exe $(HOME)/*.dSYM $(HOME)/../include/VisitWriter/*.o $(HOME)/../include/VisitWriter/*.d $(HOME)/../include/VisitWriter/*.exe $(HOME)/../include/VisitWriter/*.dSYM
-include $(OBJS:.o=.d)
