# Source, Executable, Includes, Library Defines
INCL   = def.h
#SRC    = main.c func.c intensity.c neutron_star.c operators.c star_shape.c disk_shape.c diagram_maker.c
SRC1    = main.c volume.c func.c intensity.c operators.c star_shape.c disk_shape.c sigma_phot.c
SRC2    = diagram_maker.c neutron_star.c operators.c
#SRC3    = volume.c func.c intensity.c operators.c star_shape.c disk_shape.c

OBJ1    = $(SRC1:.c=.o)
OBJ2    = $(SRC2:.c=.o)
#OBJ3    = $(SRC3:.c=.o)

LIBS   = -lm -fopenmp

EXE1	= disco
EXE2	= diagram
#EXE3    = volume



all: $(EXE1) $(EXE2) #$(EXE3)

# Compiler, Linker Defines
CC      = gcc
#CFLAGS  = -ansi -pedantic -Wall -O2
CFLAGS  = -fopenmp -g
LIBPATH =
LDFLAGS1 = -o $(EXE1) $(LIBPATH) $(LIBS)
LDFLAGS2 = -o $(EXE2) $(LIBPATH) $(LIBS)
#LDFLAGS3 = -o $(EXE3) $(LIBPATH) $(LIBS)

#CFDEBUG = -ansi -pedantic -Wall -g -DDEBUG $(LDFLAGS)
CFDEBUG =
RM      = /bin/rm -f

# Compile and Assemble C Source Files into Object Files
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c	

# Link all Object Files with external Libraries into Binaries
$(EXE1): $(OBJ1)
	$(CC) $(OBJ1) $(LDFLAGS1)

$(EXE2): $(OBJ2)
	$(CC) $(OBJ2) $(LDFLAGS2)

#$(EXE3): $(OBJ3)
#	$(CC) $(OBJ3) $(LDFLAGS3)



# Objects depend on these Libraries
$(OBJ1): $(INCL)
$(OBJ2): $(INCL)
#$(OBJ3): $(INCL)


# Create a gdb/dbx Capable Executable with DEBUG flags turned on
# debug:
# 	$(CC) $(CFDEBUG) $(SRC)

# Clean Up Objects, Dumps out of source directory
clean:
	$(RM) $(OBJ1) $(OBJ2) core a.out
