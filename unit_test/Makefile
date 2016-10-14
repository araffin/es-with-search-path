# Declaration of variables
CPPFLAGS = -std=c++11 -O3 -Wall -Wno-sign-compare -Wconversion -fpermissive

########################################################################
## Programs
algo4: algo4.o
	g++ ${CPPFLAGS} -o experiment algo4.o ${LDFLAGS}  

########################################################################
## Additional dependencies
algo4.o: test_algo.h test_algo.cpp
	g++ -c ${CPPFLAGS} -o algo4.o test_algo.cpp

# To remove generated files
clean:
	rm -f *.o
	rm -f *\~ *\#

all: algo4
