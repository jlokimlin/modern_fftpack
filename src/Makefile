
include ../make.inc

LINKS = \
-I$(OBJ_DIRECTORY) \
-I$(LIB_DIRECTORY)

EXECUTABLE = ./build_lib.exe

type_FFTpack.o: \
fftpack_routines.o \
type_FFTpack.f90
	$(FC) $(FFLAGS) -c type_FFTpack.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

fftpack_precision.o: \
fftpack_precision.f90
	$(FC) $(FFLAGS) -c fftpack_precision.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

auxiliary_routines.o: \
fftpack_precision.o \
auxiliary_routines.f90
	$(FC) $(FFLAGS) -c auxiliary_routines.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

error_handlers.o: \
auxiliary_routines.o \
error_handlers.f90
	$(FC) $(FFLAGS) -c error_handlers.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

complex_transform_routines.o: \
auxiliary_routines.o \
complex_transform_routines.f90
	$(FC) $(FFLAGS) -c complex_transform_routines.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

complex_initialization_routines.o: \
complex_transform_routines.o \
complex_initialization_routines.f90
	$(FC) $(FFLAGS) -c complex_initialization_routines.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

complex_backward_1d.o: \
complex_transform_routines.o \
complex_backward_1d.f90
	$(FC) $(FFLAGS) -c complex_backward_1d.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

complex_forward_1d.o: \
complex_transform_routines.o \
complex_forward_1d.f90
	$(FC) $(FFLAGS) -c complex_forward_1d.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

multiple_complex_forward.o: \
complex_transform_routines.o \
multiple_complex_forward.f90
	$(FC) $(FFLAGS) -c multiple_complex_forward.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

complex_forward_2d.o: \
multiple_complex_forward.o \
complex_forward_2d.f90
	$(FC) $(FFLAGS) -c complex_forward_2d.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

multiple_complex_backward.o: \
complex_transform_routines.o \
multiple_complex_backward.f90
	$(FC) $(FFLAGS) -c multiple_complex_backward.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

complex_backward_2d.o: \
multiple_complex_backward.o \
complex_backward_2d.f90
	$(FC) $(FFLAGS) -c complex_backward_2d.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

real_transform_routines.o: \
auxiliary_routines.o \
complex_transform_routines.o \
real_transform_routines.f90
	$(FC) $(FFLAGS) -c real_transform_routines.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

real_initialization_routines.o: \
real_transform_routines.o \
real_initialization_routines.f90
	$(FC) $(FFLAGS) -c real_initialization_routines.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

real_backward_1d.o: \
real_transform_routines.o \
real_backward_1d.f90
	$(FC) $(FFLAGS) -c real_backward_1d.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

real_forward_1d.o: \
real_transform_routines.o \
real_forward_1d.f90
	$(FC) $(FFLAGS) -c real_forward_1d.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

fftpack_routines.o: \
auxiliary_routines.o \
complex_transform_routines.o \
real_transform_routines.o \
fftpack_routines.f90
	$(FC) $(FFLAGS) -c fftpack_routines.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

main.o: \
type_FFTpack.o \
main.f90
	$(FC) $(FFLAGS) -c main.f90 $(LINKS) \
	-o $(OBJ_DIRECTORY)/$@

exec: $(OBJ)
	$(FC) $(FFLAGS) -o $(EXECUTABLE) $(OBJS) $(LINKS)

lib: 
	$(AR) rc lib$(LIB_NAME).a $(OBJ_DIRECTORY)/*.o
	mv lib$(LIB_NAME).a $(LIB_DIRECTORY)

all:
	make clean
	make exec; $(EXECUTABLE)
	make lib

clean:
	rm -f $(OBJS) $(LIB) $(LIB_DIRECTORY)/*.mod $(LIB_DIRECTORY)/*.smod
	rm -f $(EXECUTABLE)