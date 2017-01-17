CXX=g++
OPTFLAGS=-O3 -fopenmp
LIBPATH=-Ilib/ -ldl 

all: CreateBloom FindExons

SRCS_B=CreateBloom.cpp lib/Uncompress.cpp lib/SignalHandler.cpp lib/Fcontrol.cpp 
SRCS_C=FindExons.cpp lib/FastaReader.cpp lib/Sequence.cpp lib/Options.cpp

CreateBloom: $(SRCS_B)
	$(CXX) $(OPTFLAGS) $(LIBPATH) -o $@ $^

FindExons: $(SRCS_C)
	$(CXX) $(OPTFLAGS) $(LIBPATH) -o $@ $^

clean:
	rm CreateBloom FindExons

