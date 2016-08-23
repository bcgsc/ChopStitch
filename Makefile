CXX=g++
OPTFLAGS=-O3 -fopenmp
LIBPATH=-Ilib -ldl

all: CreateBloom

SRCS=CreateBloom.cpp lib/Uncompress.cpp lib/SignalHandler.cpp lib/Fcontrol.cpp 

CreateBloom: $(SRCS)
	$(CXX) $(OPTFLAGS) $(LIBPATH) -o $@ $^

clean:
	rm CreateBloom 

