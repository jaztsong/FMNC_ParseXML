CC=g++
CFLAGS= -Wall -g 
LFLAGS=-fopenmp
LIBS=-lboost_system -lboost_filesystem 
SRCS= main.cc fmnc_parser.cc pugixml-1.7/src/pugixml.cpp

OBJS = $(SRCS:.cc=.o)

MAIN=parseXML

$(MAIN): $(OBJS) 
	        $(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

main.o: main.cc fmnc_parser.h
	    $(CC) $(CFLAGS) -fopenmp  -c main.cc

fmnc_parser.o: fmnc_parser.cc fmnc_parser.h
	    $(CC) $(CFLAGS) -c fmnc_parser.cc
pugixml.o: pugixml-1.7/src/pugixml.cpp pugixml-1.7/src/pugixml.hpp pugixml-1.7/src/pugiconfig.hpp
	    $(CC) $(CFLAGS) -c pugixml-1.7/src/pugixml.cpp
clean: 
		$(RM) $(MAIN) *.o *~
