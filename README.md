## Real-time SDR for mono/stereo FM and RDS

The main objective of the project is to navigate a complex speciﬁcation and understand the challenges that must be addressed for a real-time implementation of a computing system that operates in a form factor-constrained environment. 

The project specification except the RDS part is in the [partial project document](doc/3dy4-project-partial.pdf). The RDS part will be released on or before March 8.

---
# How to Compile/Execute
## 1. Mono
- Edit `src/Makefile`
```Makefile
# this makefile is intended for g++ on Linux
CC = g++
# CFLAGS = -c -Wall -Wpedantic
CFLAGS = -c -O3
LDFLAGS =
SOURCES = mono.cpp mode0.cpp mode1.cpp iofunc.cpp filter.cpp  fourier.cpp  genfunc.cpp  logfunc.cpp fmPll.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = mono

all: $(EXECUTABLE) clean

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -pthread -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	-rm $(OBJECTS)
```
- Execute  

For Mode 0
```bash
make && cat ../data/samples.raw | ./mono | aplay -c 1 -f s16_le -r 48000
```
 or
 ```bash
make && cat ../data/samples.raw | ./mono 0 | aplay -c 1 -f s16_le -r 48000
```
and for Mode 1
```bash 
make && cat ../data/samples.raw | ./mono 1 | aplay -c 1 -f s16_le -r 48000
```

## 2. Stereo Mode 0
- Edit `src/Makefile`
```Makefile
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
```
- Execute  
```bash
make && cat ../data/samples.raw | ./stereoThreadingMode0 | aplay -c 2 -f s16_le -r 48000
```
## 3. Stereo Mode 1
- Edit `src/Makefile`
```Makefile
# this makefile is intended for g++ on Linux
CC = g++
# CFLAGS = -c -Wall -Wpedantic
CFLAGS = -c -O3
LDFLAGS =
SOURCES = stereoThreadingMode1.cpp iofunc.cpp filter.cpp  fourier.cpp  genfunc.cpp  logfunc.cpp fmPll.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = stereoThreadingMode1

all: $(EXECUTABLE) clean

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -pthread -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	-rm $(OBJECTS)
```
- Execute  
```bash
make && cat ../data/samples.raw | ./stereoThreadingMode1 | aplay -c 2 -f s16_le -r 48000
```
