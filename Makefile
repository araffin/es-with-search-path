## Makefile to build C example programs included with the COCO distribution
##
## NOTE: We have tried to make this Makefile as generic and portable
## as possible. However, there are many (incompatible) versions of
## make floating around. We regularly test using GNU make and BSD make
## from FreeBSD. If you have trouble compiling the examples, please
## try to use GNU make. 
##
## On Windows it is best to use either the included NMakefile by running
##
##   nmake -f NMakefile
##
## or installing Cygwin and running GNU make from within Cygwin.

LDFLAGS += -lm
#CCFLAGS = -g -ggdb -std=c89 -pedantic -Wall -Wextra -Wstrict-prototypes -Wshadow -Wno-sign-compare -Wconversion
CPPFLAGS = -std=c++11 -Wall -Wextra -Wstrict-prototypes -Wshadow -Wno-sign-compare -Wconversion -fpermissive
########################################################################
## Toplevel targets
all: experiment

clean:
	rm -f coco.o 
	rm -f experiment.o experiment 
	rm -f algo4.o

########################################################################
## Programs
experiment: experiment.o coco.o algo4.o
	g++ ${CPPFLAGS} -o experiment coco.o algo4.o experiment.o ${LDFLAGS}  

########################################################################
## Additional dependencies
algo4.o: algo4.cpp
	g++ -c ${CPPFLAGS} -o algo4.o algo4.cpp

coco.o: coco.h coco.c
	g++ -c ${CPPFLAGS} -o coco.o coco.c
experiment.o: coco.h coco.c experiment.c
	g++ -c ${CPPFLAGS} -o experiment.o experiment.c