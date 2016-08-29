CXX=g++
OPTFLAGS=-O3 -fopenmp
LIBPATH=-Ilib/ -ldl 

all: CreateBloom ChopStitch

SRCS_B=CreateBloom.cpp lib/Uncompress.cpp lib/SignalHandler.cpp lib/Fcontrol.cpp 
SRCS_C=chopstitch.cpp

CreateBloom: $(SRCS_B)
	$(CXX) $(OPTFLAGS) $(LIBPATH) -o $@ $^

ChopStitch: $(SRCS_C)
	$(CXX) $(OPTFLAGS) -o $@ $^

clean:
	rm CreateBloom ChopStitch

