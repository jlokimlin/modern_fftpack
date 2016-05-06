
include ../make.inc

SRC=type_FFTpack.f90

OBJ=$(subst .f90,.o,$(SRC))
OBJS=$(addprefix ../objs/,$(OBJ))

$(LIB) : $(OBJS)
	$(AR) -rv $@ $? 

../objs/%.o : %.f90
	$(FC) -c $< -o $@

clean:
	rm -f $(LIB) $(OBJS) 
