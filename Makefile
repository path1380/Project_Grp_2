# Makefile for project

FC = ifort
LD = ifort
SRC = src
BIN = bin
#F90FLAGS = -fbounds-check -Wall -fbacktrace
F90FLAGS = -fopenmp
FFLAGS = -O3 -fopenmp
LDFLAGS = -mkl -qopt-matmul 
#LDFLAGS= -lblas -llapack
_OBJECTS = type_defs.o quad_element.o legendre_module.o problemsetup.o weights.o main.o
OBJECTS = $(patsubst %,$(BIN)/%,$(_OBJECTS))
EXECUTABLE = main.x

.PHONY: clean
$(EXECUTABLE): $(OBJECTS)
	       $(LD) -o $(EXECUTABLE) $(OBJECTS) -I$(BIN)/ $(LDFLAGS)
$(BIN)/%.o : $(SRC)/%.f90
	   $(FC) $(F90FLAGS) -c -o $@ $< -L$(BIN)/
$(BIN)/%.o : $(SRC)/%.f
	   $(FC) $(FFLAGS) -c -o $@ $<
clean:
	rm -f $(OBJECTS) $(EXECUTABLE) $(BIN)/*.mod
