# Makefile for Homework 5 main file

FC = gfortran
LD = gfortran
SRC = src
BIN = bin
F90FLAGS = -fbounds-check -Wall -fbacktrace -g
FFLAGS = -O3
_SRCS = type_defs.f90 quad_element.f90 legendre_module.f90 problemsetup.f90 assemble.f90 coeff2d.f90 main.f90
SRCS = $(patsubst %,$(SRC)/%,$(_SRCS))
_OBJS = $(_SRCS:.f90=.o)
OBJS = $(patsubst %,$(BIN)/%,$(_OBJS))
EXECUTABLE = test.x
# EXECUTABLE = test_pis.x

.PHONY: clean compile build run

compile:
	$(FC) $(FFLAGS) -c $(SRCS)
	mv *.mod $(BIN)
	mv *.o $(BIN)

build: $(OBJS)
	$(LD) -o $(EXECUTABLE) $(OBJS) -llapack

run: $(EXECUTABLE)
	./$(EXECUTABLE)
	
clean:
	rm -f $(OBJS) $(EXECUTABLE) $(BIN)/*.mod
	rm -f *.o *.mod
