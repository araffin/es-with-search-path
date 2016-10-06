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
CCFLAGS = -g -ggdb -std=c89 -pedantic -Wall -Wextra -Wstrict-prototypes -Wshadow -Wno-sign-compare -Wconversion

########################################################################
## Toplevel targets
all: algo_4

clean:
	rm -f coco.o 
	rm -f algo_4.o algo_4 

########################################################################
## Programs
algo_4: algo_4.o coco.o
	${CC} ${CCFLAGS} -o algo_4 coco.o algo_4.o ${LDFLAGS}  

########################################################################
## Additional dependencies
coco.o: coco.h coco.c
	${CC} -c ${CCFLAGS} -o coco.o coco.c
algo_4.o: coco.h coco.c algo_4.c
	${CC} -c ${CCFLAGS} -o algo_4.o algo_4.c