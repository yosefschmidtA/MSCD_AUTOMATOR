FFLAGS = -O3
gfortranFLAGS=$(FFLAGS)
gfortran=$(FC)
.o.a:

prog:	phsh0

run:	phsh0
	phsh0

phsh0= abinitio.o atsolve.o augment.o clebschgordan.o elsolve.o \
exchcorr.o fitx0.o fourier.o getillls.o getpot.o hb.o hfdisk.o \
initiali.o integ.o parabreg.o phsh0.o pseudize.o pseudo.o \
setgrid.o setqmm.o

phsh0: $(phsh0)
	 $(FC) $(FFLAGS) $(phsh0) -lm -o phsh0
