# Source, Executable, Includes, Library Defines
INCL   = def.h
SRC    = main.c func.c intensity.c neutron_star.c
OBJ    = $(SRC:.c=.o)
LIBS   = -lm -fopenmp
EXE    = disco

# Compiler, Linker Defines
CC      = gcc
#CFLAGS  = -ansi -pedantic -Wall -O2
CFLAGS  = -fopenmp
LIBPATH =
LDFLAGS = -o $(EXE) $(LIBPATH) $(LIBS)
#CFDEBUG = -ansi -pedantic -Wall -g -DDEBUG $(LDFLAGS)
CFDEBUG =
RM      = /bin/rm -f

# Compile and Assemble C Source Files into Object Files
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

# Link all Object Files with external Libraries into Binaries
$(EXE): $(OBJ)
	$(CC) $(LDFLAGS) $(OBJ)

# Objects depend on these Libraries
$(OBJ): $(INCL)

# Create a gdb/dbx Capable Executable with DEBUG flags turned on
# debug:
# 	$(CC) $(CFDEBUG) $(SRC)

# Clean Up Objects, Dumps out of source directory
clean:
	$(RM) $(OBJ) core a.out
