# Makefile for project

FC = gfortran
LD = gfortran
SRC = src
BIN = bin
F90FLAGS = -fbounds-check -Wall -fbacktrace -g
FFLAGS = -O3
<<<<<<< HEAD
_SRCS = type_defs.f90 quad_element.f90 legendre_module.f90 problemsetup.f90 solve_quad_flux.f90 main.f90
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
	
=======
LDFLAGS= -llapack -lblas
_OBJECTS = type_defs.o quad_element.o legendre_module.o problemsetup.o weights.o main.o
OBJECTS = $(patsubst %,$(BIN)/%,$(_OBJECTS))
EXECUTABLE = main.x

.PHONY: clean
$(EXECUTABLE): $(OBJECTS)
	       $(LD) -o $(EXECUTABLE) $(OBJECTS) -I$(BIN)/ $(LDFLAGS)
$(BIN)/%.o : $(SRC)/%.f90
	   $(FC) $(F90FLAGS) -c -o $@ $< -J$(BIN)/
$(BIN)/%.o : $(SRC)/%.f
	   $(FC) $(FFLAGS) -c -o $@ $<
>>>>>>> f994267bfc8ba6e0efb47082e09d62d24017ea2f
clean:
	rm -f $(OBJECTS) $(EXECUTABLE) $(BIN)/*.mod
