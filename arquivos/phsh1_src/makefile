FFLAGS = -g
FORT77FLAGS=$(FFLAGS)
FORT77=$(FC)
.o.a:

prog:	phsh1

run:	phsh1
	phsh1

phsh1= chgrid.o clemin.o hsin.o mad.o mtz.o mtzm.o nbr.o phsh1_mscd.o \
poison.o rela.o sumax.o

phsh1: $(phsh1)
	 $(FC) $(FFLAGS) $(phsh1) -lm -o phsh1
