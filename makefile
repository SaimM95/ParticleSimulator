CC=mpicc
CFLAGS=-c -g -Wall -std=c++11
LDFLAGS=-lm
OBJECTS=$(SOURCES:.cpp=.o)

SOURCES=main.cpp savebmp.cpp

EXECUTABLE=x.project

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC)  $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o x.*
