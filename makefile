FC      := gfortran
FFLAGS  := -O3
OBJ     := chem.o constants.o jedi.o
TARGET  := jedi
BINPATH := $(HOME)/.local/bin

.PHONY:	clean all test

all: $(TARGET) clean test

%.o: %.f90
		@$(FC) -c $< $(FFLAGS)

$(TARGET): $(OBJ)
		@$(FC) -o $@ $^

test:
		@./$(TARGET) test/1.xyz test/5.xyz test/redundant.dat test/hessianINT.dat
		@if [ -s ./strain.dat ]; then \
			echo "The program passed the test"; \
		else \
			echo "Test failed"; \
		fi
clean:
		@rm *.o *.mod
install:
		@cp $(TARGET) $(BINPATH)
		@echo "Jedi installed in $(HOME)/.local/bin"
