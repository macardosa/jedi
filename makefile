FC      = gfortran
FFLAGS  = -O3
OBJ     = chem.o constants.o jedi.o
TARGET  = jedi

all: $(TARGET) clean

%.o: %.f90
		$(FC) -c $< $(FFLAGS)

$(TARGET): $(OBJ)
		$(FC) -o $@ $^

.PHONY:	clean all

clean:
		rm *.o *.mod
