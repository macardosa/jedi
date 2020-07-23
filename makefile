FC      = gfortran
FFLAGS  = -O3
OBJ     = chem.o constants.o jedi.o
TARGET  = jedi

.PHONY:	clean all test

all: $(TARGET) clean test

%.o: %.f90
		$(FC) -c $< $(FFLAGS)

$(TARGET): $(OBJ)
		$(FC) -o $@ $^
test:
		./$(TARGET) geom_eq.xyz geom_def.xyz redundant.dat hessianINT.dat

clean:
		rm *.o *.mod
cpbin:
		cp $(TARGET) $(BINPATH)
