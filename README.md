# **modern\_fftpack - a Fortran library of fast Fourier transforms**

An object-oriented modernization of NCAR's FFTPACK5.1.


The library is fully Fortran 2008 (ISO/IEC 1539-1:2010) compliant, entirely free of **go to**, and employs the instrinsic module *iso\_fortran\_env* to maintain 64-bit floating-point precision. 


Test programs are provided for the transforms. Each serves two purposes: as a template to guide you in writing your own codes utilizing the modern\_fttpack library, and as a demonstration that you can correctly produce the executables. 

-----------------------------------------------------------------------------

## Requirements

* The GNU Make tool https://www.gnu.org/software/make/
* The GNU gfortran compiler https://gcc.gnu.org/wiki/GFortran

-----------------------------------------------------------------------------


## To build the project

Type the following command line arguments

```bash

	git clone https://github.com/jlokimlin/modern_fftpack.git
	
	cd modern_fftpack; make all
```

-----------------------------------------------------------------------------

## Contributing

This project is still a work in progress and anyone is free to contribute under the proviso that they abstain from using the dreaded **go to**. 

-----------------------------------------------------------------------------

## TODO
* Hide dynamically allocated workspace arrays from the end-user 
* Ideally, **type**(FFTpack)'s functionality and usage will evolve into something similar to scipy.fftpack
* Implement generic type-bound procedures %fft and %ifft to handle both complex and real transforms
* Add differential and psuedo-differential operators

-----------------------------------------------------------------------------


## Bibliography

[1] Swarztrauber, Paul N. "Symmetric FFTs." Mathematics of Computation 47.175 (1986): 323-346.

[2] Swarztrauber, Paul N. "FFT algorithms for vector computers." Parallel Computing 1.1 (1984): 45-63.

-----------------------------------------------------------------------------

## Result

	```

	program tcfft1 and related messages:
	cfft1 backward-forward max error = 0.794798730345622E-015
	program tcfft1 and related messages:
	cfft1 forward-backward max error = 0.867111901826273E-015
	end program tcfft1 and related messages
	
	program trfft1 and related messages:
	rfft1 forward-backward max error = 0.555111512312578E-015
	rfft1 backward-forward max error = 0.777156117237610E-015
	end program trfft1 and related messages
	
	program tcosq1 and related messages:
	cosq1 forward-backward max error = 0.777156117237610E-015
	cosq1 backward-forward max error = 0.666133814775094E-015
	end program tcosq1 and related messages

	program tcost1 and related messages:
	cost1 forward-backward max error = 0.233146835171283E-013
	cost1 backward-forward max error = 0.224265050974282E-013
	end program tcost1 and related messages
	
	program tsinq1 and related messages:
	sinq1 forward-backward max error = 0.777156117237610E-015
	sinq1 backward-forward max error = 0.777156117237610E-015
	end program tsinq1 and related messages
	

	program tsint1 and related messages:
	sint1 forward-backward max error = 0.455191440096314E-013
	sint1 backward-forward max error = 0.395239396766556E-013
	end program tsint1 and related messages
	

	program tcfft2 and related messages:
	cfft2 forward-backward max error = 0.948574968053509E-015
	cfft2 backward-forward max error = 0.915513359704447E-015
	end program tcfft2 and related messages
	
	program trfft2 and related messages:
	rfft2 forward-backward max error = 0.777156117237610E-015
	rfft2 backward-forward max error = 0.832667268468867E-015
	end program trfft2 and related messages
	
	 
	This result was compiled by GCC version 5.3.1 20160409 using the options -I ../lib -mtune=generic -march=x86-64 -O3 -Wall -J ../lib
	
	```

