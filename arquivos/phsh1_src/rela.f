      SUBROUTINE RELA(RHO,RX,NX,NGRID)
      DIMENSION RHO(NGRID),RX(NGRID)
      REAL NAME(4),SUM
      COMMON /WK/ RR(2000),RS(2000)
C  ROUTINE FOR INPUT OF CHARGE DENSITY FROM RELATIVISTIC ORBITALS
C  (ERIC SHIRLEY PROGRAM), AND CALCULATION OF CHARGE DENSITY ON 
C  THE RADIAL MESH RX
C    RHO = 4*PI*SUM OVER STATES OF (MODULUS(WAVE FN)**2) *
C          RADIUS**2
C    RMIN= minimum radial coordinate defining the logarithmic mesh used
C          in relativistic calculation
C    RMAX= maximum radial coordinate defining the logarithmic mesh used
C          in relativistic calculation
C    NR  = number of points in the mesh
C the mesh is defined as r(i)=rmin*(rmax/rmin)**(dfloat(i)/dfloat(nr))
C FOR EACH ATOMIC STATE I?
      READ(4,100)NAME,IPRINT
      read(4,54) rmin,rmax,nr,z
54    format (d15.8,d15.8,i5,f5.2)
c
c initialization of logarithmic grid
c
	do 5 i=1,nr
	   rr(i)=rmin*(rmax/rmin)**(dfloat(i)/dfloat(nr))
5       continue
      NS=nr
C  read in charge density
	read(4,56) (rs(j),j=1,nr)
56      format (f15.10)
c
c  INTERPOLATION TO GRID RX
c
      NX=NGRID
      CALL CHGRID(RS,RR,NS,RHO,RX,NX)
      IF(IPRINT.EQ.0)RETURN
      WRITE(11,200)NAME,(RR(IX),IX=1,NS)
C      write(11,200)(RR(IX),IX=1,NS)
      DO 7 IX=1,NX
      IF(RHO(IX).LT.1.0E-9)GOTO 8
7     CONTINUE
8     NX=IX
      WRITE(11,202)RX(NX),NX,(RHO(IX),IX=1,NX)
      RETURN
100   FORMAT(4A4/I4)
102   FORMAT(5F9.4)
200   FORMAT(1H1,4A4,36H RELAT. WAVEFUNCTIONS (ERIC SHIRLEY),
     + 9H R RADIUS,17H LOGARITHMIC MESH,/(10F12.5/))
201   FORMAT(3H0L?,I3//5(10F11.5/))
202   FORMAT(29H0CHARGE DENSITY OUT TO RADIUS,F12.5,10X,
     + 3HNX?,I4//5(10E12.4/))
      END
