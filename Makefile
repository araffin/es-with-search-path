# Declaration of variables
CC = g++
CC_FLAGS = -Wall -std=c++11

# File names
EXEC = test.x
SOURCES = test.cpp

# Main target
$(EXEC): $(SOURCES)
	$(CC) $(SOURCES) $(CC_FLAGS) -o $(EXEC)

# To remove generated files
clean:
	rm -f $(EXEC)
	rm -f *\~ *\#

all: $(EXEC)
