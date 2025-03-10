# Makefile for compiling and running Kobayashi.f90

# Set the Fortran compiler to use (ifort or gfortran)
FC = gfortran

# Compiler flags
FLAGS = -O2 -Wall

# Source file
SRC = Kobayashi.f90

# Output executable
OUT = kobayashi

.PHONY: all clean run ifort gfortran

# Default target
all: $(OUT)

# Compile the source file into the executable
$(OUT): $(SRC)
ifeq ($(FC), ifort)
	$(FC) $(FLAGS) -module $(MOD_DIR) -o $(OUT) $(SRC)
else
	$(FC) $(FLAGS) -o $(OUT) $(SRC)
endif

# Run the executable
run: $(OUT)
	./$(OUT)

# Clean the build files
clean:
	rm -f $(OUT) *.o *.mod *.txt *.out *.mp4

# Switch to ifort and compile
ifort:
	$(MAKE) clean
	$(MAKE) FC=ifort FLAGS="-O2 -warn all -I/path/to/intel/mod"

# Switch to gfortran and compile
gfortran:
	$(MAKE) clean
	$(MAKE) FC=gfortran FLAGS="-O2 -Wall"
