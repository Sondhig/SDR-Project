# this makefile is intended for g++ on Linux

CC = g++
# CFLAGS = -c -Wall -Wpedantic
CFLAGS = -c -O3
LDFLAGS =
SOURCES = stereoThreadingMode0.cpp iofunc.cpp filter.cpp  fourier.cpp  genfunc.cpp  logfunc.cpp fmPll.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = stereoThreadingMode0

all: $(EXECUTABLE) clean

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -pthread -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	-rm $(OBJECTS)
