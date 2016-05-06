
include ../make.inc

SRC=tcfft1.f90 trfft1.f90 tcosq1.f90 tcost1.f90 tsinq1.f90 tsint1.f90 \
    tcfft2.f90 trfft2.f90

EXES=$(subst .f90,.exe, $(SRC))

all : $(EXES)

$(EXES) : $(SRC)

%.exe : %.f90
	rm -f $@
	$(FC) $< -o $@ -L../lib -l$(LIB_NAME)
	./$@


clean :
	rm -f $(EXES)
