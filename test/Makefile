
include ../make.inc

SRC=$(subst .o,.f90, $(OBJ))

EXES=$(subst .o,.exe, $(OBJ))

FFTPACK_MODULE_PATH = $(LIB_DIRECTORY)

LINKS = -L$(LIB_DIRECTORY) -l$(LIB_NAME) -I$(FFTPACK_MODULE_PATH) -I$(OBJ_DIRECTORY)

all: $(EXES)

$(EXES): $(SRC)

%.exe : %.f90
	rm -f $@
	$(FC) $< -o $@ $(LINKS)
	./$@


clean :
	rm -f $(OBJS) $(EXES)