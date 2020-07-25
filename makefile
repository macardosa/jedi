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
		./$(TARGET) test/1.xyz test/5.xyz test/redundant.dat test/hessianINT.dat
		mv strain.dat test

clean:
		rm *.o *.mod
install:
		cp $(TARGET) $(BINPATH)
