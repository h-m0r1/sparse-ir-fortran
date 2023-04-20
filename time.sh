#!/bin/sh

rm Makefile
cat <<EOF > Makefile
FC = ifort

#LDLIBS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
LDLIBS = -mkl=sequential
LDFLAGS = -static-intel
F90FLAGS = -O2 -g -traceback

TARGET = time_measurement.x
OBJS = time_measurement.o sparse_ir.o sparse_ir_io.o sparse_ir_preset.o

.SUFFIXES: .f90

%.o: %.f90
	$(COMPILE.f) $(OUTPUT_OPTION) $< $(F90FLAGS)

%.mod: %.f90 %.o
	@:

$(TARGET): $(OBJS)
	$(LINK.f) $^ $(LDLIBS) $(LDFLAGS) -o $@

time_measurement.o: sparse_ir.mod sparse_ir_io.mod sparse_ir_preset.mod

.PHONY: test
test: $(TARGET)

.PHONY: clean
clean:
	rm -f *.o *.mod
EOF
make clean
make test
./time_measurement.x

rm Makefile
cat <<EOF > Makefile
FC = ifort

#LDLIBS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
LDLIBS = -mkl=sequential
LDFLAGS = -static-intel
F90FLAGS = -O2 -g -traceback

TARGET = time_measurement.x
OBJS = time_measurement.o sparse_ir1.o sparse_ir_io.o sparse_ir_preset.o

.SUFFIXES: .f90

%.o: %.f90
	$(COMPILE.f) $(OUTPUT_OPTION) $< $(F90FLAGS)

%.mod: %.f90 %.o
	@:

$(TARGET): $(OBJS)
	$(LINK.f) $^ $(LDLIBS) $(LDFLAGS) -o $@

time_measurement.o: sparse_ir1.mod sparse_ir_io.mod sparse_ir_preset.mod

.PHONY: test
test: $(TARGET)

.PHONY: clean
clean:
	rm -f *.o *.mod
EOF
make clean
make test
./time_measurement.x

rm Makefile
cat <<EOF > Makefile
FC = ifort

#LDLIBS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
LDLIBS = -mkl=sequential
LDFLAGS = -static-intel
F90FLAGS = -O2 -g -traceback

TARGET = time_measurement.x
OBJS = time_measurement.o sparse_ir2.o sparse_ir_io.o sparse_ir_preset.o

.SUFFIXES: .f90

%.o: %.f90
	$(COMPILE.f) $(OUTPUT_OPTION) $< $(F90FLAGS)

%.mod: %.f90 %.o
	@:

$(TARGET): $(OBJS)
	$(LINK.f) $^ $(LDLIBS) $(LDFLAGS) -o $@

time_measurement.o: sparse_ir2.mod sparse_ir_io.mod sparse_ir_preset.mod

.PHONY: test
test: $(TARGET)

.PHONY: clean
clean:
	rm -f *.o *.mod
EOF
make clean
make test
./time_measurement.x