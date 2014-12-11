# modify the path of your Cantera.mak file
include  ($YOUR_CANTERA_INSTALLATION_FOLDER$)/build/platform/Cantera.mak

CC=gcc
CXX=g++
RM=rm -f
CCFLAGS=-g  -Wall -fpermissive #debug options
CPPFLAGS=$(CANTERA_INCLUDES)
LDFLAGS=
LDLIBS=$(CANTERA_LIBS) 

SRCS=flamelet.cpp
OBJS=  StFlow_2.o flamelet.o Lewis1Transport.o TransportFactory_2.o 

all: flamelet_code

flamelet_code: StFlow_2.o TransportFactory_2.o  Lewis1Transport.o flamelet.o flamelet 


StFlow_2.o: StFlow_2.cpp
	$(CXX) $(CCFLAGS) -c StFlow_2.cpp $(CPPFLAGS)

Lewis1Transport.o: Lewis1Transport.cpp 
	$(CXX) $(CCFLAGS) -c Lewis1Transport.cpp $(CPPFLAGS)

TransportFactory_2.o: TransportFactory_2.cpp
	$(CXX) $(CCFLAGS)-c TransportFactory_2.cpp $(CPPFLAGS)


flamelet.o: flamelet.cpp
	$(CXX)  -c  -g -Wall -fpermissive flamelet.cpp $(CPPFLAGS) 

flamelet: flamelet.o TransportFactory_2.o Lewis1Transport.o
	$(CXX)  flamelet.o StFlow_2.o TransportFactory_2.o Lewis1Transport.o -lf2c $(LDLIBS) -lblas -llapack  -o flamelet 
	cp flamelet ../ 

clean:
	$(RM) *.o 


