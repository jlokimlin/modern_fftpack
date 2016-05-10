
include ../make.inc

SRC=$(subst .o,.f90, $(OBJ))

$(LIB): $(OBJS)
	$(AR) -rv $@ $? 

../objs/%.o : %.f90
	$(FC) -c $< -o $@

clean:
	rm -f $(LIB) $(OBJS) 
