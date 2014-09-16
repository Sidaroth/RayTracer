# The compiler
CC=g++
# Flags for the compiler
CFLAGS=-c -Wall

# Name of the executable
EXECUTABLE=RayTracer.exe

# Output files from the executable
OUTPUT=

# Edit these to add additional libraries. 
# Libraries to include. (Prefix all additional libs with -L, like such: -L/example/dir/lib)
LIB=

# Library files to include (Prefix all addtional include dirs with -I, like such: -I/example/dir/include)
INC=

# Flags for the linker. (When adding new libraries, add another flag here, as such: -lexample-lib) PS: These are the libexample-lib.so files, with the lib prefix and .so extension removed. 
#CURLFLAGS=   -lcurl

LFLAGS= #$(CURLFLAGS)

# Source files (Add your files here for compilation. *Only the .cpp files*)
SOURCES= ./src/main.cpp


### !!!! DO NOT TOUCH ANYTHING BELOW THIS LINE EXCEPT CLEAN!!!! ###
# The objects
OBJECTS=$(SOURCES:.cpp=.o)

build: $(SOURCES) $(EXECUTABLE)

all: 
	build

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LIB) $(LFLAGS) -o $@

.cpp.o:
	$(CC) $(INC) $(CFLAGS) $< -o $@


# Run: make clean
# if you need to update object files / the executable and the make file does not want to update. This can sometimes happen if only header files have been edited. 
# Edit this (read: add to this) if you need to remove additional files when running "make clean". 
clean:
	rm -f *.o
	rm -f $(EXECUTABLE)
	rm -f $(OUTPUT)
