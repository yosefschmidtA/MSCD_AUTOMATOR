C
C---------------------------------------------------------------------
C  program PHSH1.FOR
C---------------------------------------------------------------------
C
C  adated from CAVPOT 
C
      PROGRAM CAVPOT
C  PHASE SHIFT PROGRAM FROM CAVLEED PACKAGE (PENDRY-TITTERINGTON).
C
C  MAIN PROGRAM FOR COMPUTATION OF MUFFIN-TIN POTENTIAL
C  AFTER MATTHEISS.   REF - LOUCKS, 'APW METHOD', 1967
C  INCLUDING? HARTREE POTENTIAL, SLATER-TYPE STATISTICAL
C  EXCHANGE TERM, AND MADELUNG CORRECTION FOR IONIC
C  MATERIALS
C  DIMENSIONED FOR UP TO NTOT (24) ATOMS IN THE UNIT CELL, OF NIEQ
C  INEQUIVALENT TYPES (BUT NIEQ<14)
      PARAMETER (NIEQ=18,NTOT=40)
      DIMENSION SIG(250,NIEQ),RHO(250,NIEQ),VH(250,NIEQ),
     + VS(250,NIEQ),VMAD(550,NIEQ),RX(550),RS(550),POT(550),rmscd(550)
      DIMENSION RC(3,3),RK(3,NTOT),ZM(NTOT),Z(NIEQ),ZC(NIEQ),
     + RMT(NIEQ),JRMT(NIEQ),JRMT2(NIEQ),NRR(NIEQ),NCON(NIEQ),NX(NIEQ)
      DIMENSION IA(30,NIEQ),NA(30,NIEQ),AD(30,NIEQ)
c      REAL TITLE(20),NAME(4,NIEQ),WFN,WFN0,WFN1,WFN2,WFN3,SUM
      character*10 title(20)
      character name(4,nieq)
      character*4 wfn,wfn0,wfn1,wfn2,wfn3
      double precision sum
      COMMON /WK/ WK1(250),WK2(250)
      COMMON /WF/ WF2(250,14),WC(14),LC(14)
      DATA NGRID,MC,PI/250,30,3.1415926536/,WFN0,WFN1,WFN2,WFN3/
     +4HRELA,4HHERM,4HCLEM,4HPOTE/
      INDEX(X)=20.0*(ALOG(X)+8.8)+2.0
C
C First input channels
C
C     OPEN (UNIT=35,FILE='input_0.dat',STATUS='OLD')
      OPEN (UNIT=4,FILE='atomic.i',STATUS='OLD')
      OPEN (UNIT=7,FILE='cluster.i',STATUS='OLD')
C
C Now output channels
C
      OPEN (UNIT=11,FILE='check.o',STATUS='UNKNOWN')
      OPEN (UNIT=9,FILE='mufftin.d',STATUS='UNKNOWN')
      OPEN (UNIT=35,FILE='input.dat')
C
C  INITIALISATION OF LOUCKS' EXPONENTIAL MESH
      X=-8.8
      DO 1 IX=1,NGRID
      RX(IX)=EXP(X)
1     X=X+0.05
C
      READ(7,100)TITLE
      WRITE(11,200)TITLE
C
C  INPUT OF CRYSTALLOGRAPHIC DATA
C    SPA = LATTICE CONSTANT IN A.U.
C    RC(I,J) = I'TH COORDINATE OF THE J'TH AXIS OF UNIT CELL,
C    IN UNITS OF SPA
C    RK(I,J) = I'TH COORDINATE OF THE J'TH ATOM IN UNIT CELL,
C    IN UNITS OF SPA
C    NR = NUMBER OF INEQUIVALENT ATOMS IN UNIT CELL
C  FOR AN ATOM OF TYPE IR?
C    NRR(IR) = NUMBER IN UNIT CELL
C    Z(IR) = ATOMIC NUMBER
C    ZC(IR) = VALENCE CHARGE
C    RMT(IR) = MUFFIN-TIN RADIUS
      READ(7,101)SPA
      READ(7,101)((RC(I,J),I=1,3),J=1,3)
      DO 2 I=1,3
      DO 2 J=1,3
2     RC(I,J)=SPA*RC(I,J)
      READ(7,102)NR
      DO 3 IR=1,NR
      DO 3 I=1,NGRID
      VH(I,IR)=0.0
      VS(I,IR)=0.0
      VMAD(I,IR)=0.0
      SIG(I,IR)=0.0
3     RHO(I,IR)=0.0
      VHAR=0.0
      VEX=0.0
C
      JJ=0
      ZZ=0.0
      DO 4 IR=1,NR
      READ(7,100)(NAME(I,IR),I=1,4)
      READ(7,103)NRR(IR),Z(IR),ZC(IR),RMT(IR)
      ZZ=ZZ+ABS(ZC(IR))
      JRMT(IR)=INDEX(RMT(IR))
      N=NRR(IR)
      DO 4 J=1,N
      JJ=JJ+1
      ZM(JJ)=ZC(IR)
      READ(7,101)(RK(I,JJ),I=1,3)
      DO 4 I=1,3
4     RK(I,JJ)=SPA*RK(I,JJ)
C    N = TOTAL NUMBER OF ATOMS IN UNIT CELL
C    AV = TOTAL VOLUME OF UNIT CELL
C    OMA = ATOMIC VOLUME
C    RWS = WIGNER-SEITZ RADIUS
      N=JJ
      RCC1=RC(2,2)*RC(3,3)-RC(3,2)*RC(2,3)
      RCC2=RC(3,2)*RC(1,3)-RC(1,2)*RC(3,3)
      RCC3=RC(1,2)*RC(2,3)-RC(2,2)*RC(1,3)
      AV=ABS(RC(1,1)*RCC1+RC(2,1)*RCC2+RC(3,1)*RCC3)
      OMA=AV/FLOAT(N)
      RWS=(0.75*OMA/PI)**(1.0/3.0)
      JRWS=INDEX(RWS)
      WRITE(11,201)((RC(I,J),I=1,3),J=1,3)
      WRITE(11,202)AV,OMA,RWS
      JJ=0
      DO 6 IR=1,NR
      WRITE(11,203)IR,(NAME(I,IR),I=1,4),NRR(IR)
      INR=NRR(IR)
      DO 5 IIR=1,INR
      JJ=JJ+1
5     WRITE(11,204)(RK(I,JJ),I=1,3)
6     WRITE(11,205)Z(IR),ZC(IR),RMT(IR)
      WRITE(11,216)(RX(IX),IX=1,NGRID)
C
C  FOR EACH ATOMIC TYPE, READ IN ATOMIC WAVEFUNCTIONS FOR NEUTRAL
C  ATOM, IN EITHER THE HERMAN-SKILLMAN OR CLEMENTI FORM, PRODUCING?
C    RHO = 4*PI*CHARGE DENSITY * RADIUS**2
      MIX=0
      DO 11 IR=1,NR
      READ(4,100)WFN
C  OPTION 0)  RELATIVISTIC CHARGE DENSITY INPUT
      IF(WFN.EQ.WFN0)CALL RELA(RHO(1,IR),RX,NX(IR),NGRID)
C  OPTION 1)  HERMAN-SKILLMAN INPUT
      IF(WFN.EQ.WFN1)CALL HSIN(RHO(1,IR),RX,NX(IR),NGRID)
C  OPTION 2)  CLEMENTI INPUT
      IF(WFN.EQ.WFN2)CALL CLEMIN(RHO(1,IR),RX,NX(IR),NGRID)
C  OPTION 3)  POTENTIAL INPUT
      IF(WFN.EQ.WFN3)GOTO 14
C  RHO IS NORMALISED USING TOTAL ELECTRONIC CHARGE ON THE ATOM
C  CALCULATED BY THE TRAPEZOIDAL RULE
7     NIX=NX(IR)
      MIX=MAX0(NIX,MIX)
      SUM=0.0D0
      W1=0.025*RHO(1,IR)*RX(1)
C      JRXX=JRMT(IR)
      DO 8 IX=2,NIX
      W2=0.025*RHO(IX,IR)*RX(IX)
      SUM=SUM+W1+W2
8     W1=W2
      ZE=SUM
C      ANORM=Z(IR)/ZE
C      DO 9 IX=1,NIX
C9     RHO(IX,IR)=RHO(IX,IR)*ANORM
C  SOLVE POISSON'S EQUATION
C    SIG = COULOMB POTENTIAL
C    RHO = 4*PI*CHARGE DENSITY*RADIUS SQUARED
      CALL POISON(RHO(1,IR),Z(IR),NIX,SIG(1,IR))
      X=-8.8
      DO 10 IX=1,NIX
      CE=EXP(-0.5*X)
      SIG(IX,IR)=CE*(-2.0*Z(IR)*CE+SIG(IX,IR))
      RHO(IX,IR)=RHO(IX,IR)/(RX(IX)**2)
10    X=X+0.05
      WRITE(11,206)(NAME(I,IR),I=1,4),ZE,RX(NIX),NIX
      WRITE(11,207)(SIG(IX,IR),IX=1,NIX)
11    CONTINUE
C
C  DETAILS OF NEIGHBOURING SHELLS FOR EACH ATOMIC TYPE IR?
C    NCON(IR) = NUMBER OF SHELLS INCLUDED
C    IA(J,IR) = ATOMIC TYPE IN J'TH SHELL
C    NA(J,IR) = NUMBER OF ATOMS IN J'TH SHELL
C    AD(J,IR) = DISTANCE TO J'TH SHELL
      RMAX=RX(MIX)
      CALL NBR(IA,NA,AD,NCON,NRR,NR,RC,RK,N,RMAX,MC)
      WRITE(11,208)
      DO 12 IR=1,NR
      WRITE(11,209)IR
      NC=NCON(IR)
      IC=(NC-1)/12+1
      KC=0
      DO 12 I=1,IC
      JC=KC+1
      KC=MIN0(NC,KC+12)
      WRITE(11,210)(AD(J,IR),J=JC,KC)
      WRITE(11,211)(NA(J,IR),J=JC,KC)
12    WRITE(11,212)(IA(J,IR),J=JC,KC)
      READ(7,102) nform
C
C  CALCULATION OF THE MUFFIN-TIN POTENTIAL FOR EACH NEUTRAL
C  ATOM, FOLLOWING THE MATTHEISS PRESCRIPTION
C  READ IN ALPHA FOR THE SLATER EXCHANGE TERM
      READ(7,101)ALPHA
      WRITE(11,215)ALPHA
      PD=6.0/(PI*PI)
      DO 13 IR=1,NR
      JRX=MAX0(JRWS,JRMT(IR))
C  SUMMING THE POTENTIALS FROM NEUTRAL ATOMS
C    VH = HARTREE POTENTIAL
      CALL SUMAX(VH(1,IR),SIG,RX,NX,NCON(IR),IA(1,IR),NA(1,IR),
     + AD(1,IR),JRX,NGRID,NR)
C  SUMMING THE CHARGE DENSITY ABOUT EACH ATOMIC TYPE
C    VS = TOTAL CHARGE DENSITY, THEN SLATER EXCHANGE TERM
      CALL SUMAX(VS(1,IR),RHO,RX,NX,NCON(IR),IA(1,IR),NA(1,IR),
     + AD(1,IR),JRX,NGRID,NR)
      DO 13 IX=1,JRX
13    VS(IX,IR)=-1.5*ALPHA*(PD*VS(IX,IR))**(1.0/3.0)
C
C  CALCULATE THE MUFFIN-TIN ZERO
      VINT=0.
      READ(7,102)NH
      IF(NH.EQ.0.AND.NR.EQ.1)CALL MTZM(VH(1,1),VS(1,1),RX,NGRID,
     + RMT(1),RWS,JRMT(1),JRWS,VHAR,VEX)
      IF(NH.NE.0)CALL MTZ(SIG,RHO,RX,NGRID,RMT,NRR,NX,NR,RC,RK,N,
     + VHAR,VEX,ALPHA,AV,NH)
C      write(*,*) 'Slab or Bulk calculation?'
C      write(*,*) 'input 1 for Slab or 0 for Bulk' 
       nbulk=0
C      read(*,*) nbulk 
C      read(33,*) nbulk
      if(nbulk.eq.1) then
C	 write(*,*) 'Input the MTZ value from the substrate calculation'
C	 read(*,*) esht
         read(33,*)esht
	 esh=esht-(vhar+vex)
      else
C	 write(*,*) 'If you are interested in adatoms on this substrate'
C	 write(*,*) 'rerun a slab calculation with the adatoms'
C	 write(*,*) 'and use this MTZ value as input when asked '
C        write(35,*) '1'
	 write(35,*) vhar+vex
      endif
      GOTO 16
C
C  OPTION 3)  READ IN POTENTIAL OF NEUTRAL ATOM, VH, ON RADIAL
C  GRID, RX, FOR CORRECTION BY MADELUNG SUMMATION
14    READ(4,104)NGRID,(RX(IX),IX=1,NGRID)
      DO 15 IR=1,NR
      READ(4,104)JRX,(VH(IX,IR),IX=1,JRX)
15    JRMT(IR)=JRX
C
C  THE MADELUNG CORRECTION FOR IONIC MATERIALS.   SUBROUTINE MAD
C  COMPUTES THE SPHERICALLY AND SPATIALLY AVERAGED FIELDS FOR
C  THE LATTICE OF POINT CHARGES ABOUT EACH ATOMIC TYPE
16    IF(ZZ.NE.0)CALL MAD(VMAD,RX,NGRID,RMT,NRR,JRMT,NR,
     + RC,RK,ZM,N,AV)
C
C  THE TOTAL MUFFIN-TIN POTENTIAL IS ACCUMULATED INTO SIG,
C  REFERRED TO THE MUFFIN-TIN ZERO
      VINT=VHAR+VEX
      if (nform.eq.0)write(9,102)NR
      DO 17 IR=1,NR
      WRITE(11,213)(NAME(I,IR),I=1,4),VINT,RMT(IR)
      JRX=JRMT(IR)
      DO 17 IX=1,JRX
      VH(IX,IR)=VH(IX,IR)-VHAR
      VS(IX,IR)=VS(IX,IR)-VEX
      SIG(IX,IR)=VH(IX,IR)+VS(IX,IR)+VMAD(IX,IR)
17    WRITE(11,214)RX(IX),VH(IX,IR),VS(IX,IR),VMAD(IX,IR),SIG(IX,IR)
C
C     WRITE(9,219)NGRID,(RX(IX),IX=1,NGRID)
C write output in a format to be read by WILLIAMS phase shift program (NFORM=1)
C by CAVLEED phase shift program (NFORM=0), or by the relativistic phase
C shift program (NFORM=2)
C      
C Also prepare to shift the potential by an amount of the order 
C of the bulk muffintin zero.
C This is needed only if the cluster.i file correspond to a surface adsorbate
C      esh=SIG(JRX,IR)
C      esh=-1.07
      if (nform.eq.1) write(9,220) NR
      if (nform.eq.2) then
c
c define german grid RX and save old grid in RS
c
	 RM=60.0
	 DX=0.03125
	 NMX=421
	 RS(1)=RX(1)
	 RX(1)=RM*EXP(DX*(1-NMX))
	 J=1
	 RM= EXP(DX)
  110    K=J+1
	 RS(K)=RX(K)
	 RX(K)=RM*RX(J)
	 J=K
	 IF (J.LT.NMX)  GO TO 110
      endif
      DO 18 IR=1,NR
      JRX=JRMT(IR)
      if (nform.eq.0) then
	WRITE(9,217)(NAME(I,IR),I=1,4)
	WRITE(9,218)Z(IR),RMT(IR),VINT
      elseif(nform.eq.1)then
	WRITE(9,221)Z(IR),RMT(IR) 
      elseif(nform.eq.2) then
c
c es=Emin for phase shift calculation (ev)
c de=delta E for phase shift calculation (ev)
c ue=Emax for phase shift calculation (ev)
c lsm=maximum number of phase shifts desired
	es=40.
	de=5.
	ue=480.
	lsm=12
	WRITE(9,217)(NAME(I,IR),I=1,4)
	WRITE(9,111)ES,DE,UE,LSM,VINT                                    C
111     FORMAT (3D12.4,4X,I3,4X,D12.4)
c
c  INTERPOLATION TO GRID RX
c
	 do 188 k=1,jrx
188      sig(k,IR)=(sig(k,IR)-esh)*rs(k)
	 NMXX=NMX
	 CALL CHGRID(SIG(1,IR),RS,JRX,POT,RX,NMXX)
	 IZ=Z(IR)
	 WRITE(9,105)IZ,RMT(IR),NMXX
105      FORMAT(I4,F10.6,I4)
	 JRX=NMXX
      endif
      if (nform.eq.0)write(9,102)JRX
      if(nform.eq.1)then
	DO 19 IX=1,JRX
19      WRITE(9,219)RX(IX),RX(IX)*(SIG(IX,IR)-esh)
	rneg=-1.
	WRITE(9,219) rneg
C        if (nform.eq.1) WRITE(9,219) rneg
      elseif(nform.eq.0) then
	DO 199 IX=1,JRX
199      WRITE(9,219)RX(IX),(SIG(IX,IR)-esh)

      elseif (nform.eq.2) then
	WRITE(9,106)(POT(IX),IX=1,JRX)
106     FORMAT(5E14.7)
      endif
      if (nform.eq.3) then
c         zero=log(rx(1))+0.00001
          write(9,220) NR
        do icounter=1,jrx
c         print *, rx(icounter), '   ', log(rx(icounter))
c         rmscd(icounter)=log(rx(icounter))-zero
         rmscd(icounter)=rx(icounter)
         write(9,219) rmscd(icounter),-1.0*rmscd(icounter)
     +                *(sig(icounter,ir)-esh)
        enddo
      endif
18    CONTINUE
C
      STOP
C
100   FORMAT(20A4)
101   FORMAT(3F8.4)
102   FORMAT(I4)
103   FORMAT(I4,3F8.4)
104   FORMAT(I4/(5E14.5))
200   FORMAT(30H1MUFFIN-TIN POTENTIAL PROGRAM?,5X,20A4)
201   FORMAT(///18H AXES OF UNIT CELL/(6X,3F8.4))
202   FORMAT(18H0UNIT CELL VOLUME?,F15.4/
     + 15H ATOMIC VOLUME?,F18.4/
     + 21H WIGNER-SEITZ RADIUS?,F12.4)
203   FORMAT(///5H TYPE,I2,6H ATOM?,2X,4A4/
     + I4,19H ATOMS IN UNIT CELL)
204   FORMAT(6X,3F8.4)
205   FORMAT(15H0ATOMIC NUMBER?,F15.1/9H VALENCE?,F21.1/
     + 19H MUFFIN-TIN RADIUS?,F14.4)
206   FORMAT(///1H ,4A4,19H ELECTRONIC CHARGE?,F12.5/
     + 51H0COULOMB POTENTIAL FOR ISOLATED ATOM, OUT TO RADIUS,
     + F12.5,10X,3HNX?,I4/)
207   FORMAT(5(10E12.4/))
208   FORMAT(1H1)
209   FORMAT(//34H0NEAREST NEIGHBOUR SHELLS FOR TYPE,I2,5H ATOM)
210   FORMAT(9H DISTANCE,1X,15(F8.4))
211   FORMAT(7H NUMBER,3X,15(I5,3X))
212   FORMAT(5H TYPE,5X,15(I5,3X))
213   FORMAT(1H1,4A4,5X,33HPOTENTIALS IN RYDBERGS CORRECT TO,
     + 17H MUFFIN-TIN ZERO?,F8.4/19H0MUFFIN-TIN RADIUS?,F8.4//
     + 5X,6HRADIUS,5X,17HHARTREE POTENTIAL,9X,
     + 8HEXCHANGE,4X,19HMADELUNG CORRECTION,5X,
     + 15HTOTAL POTENTIAL)
214   FORMAT(F12.5,4E20.6)
215   FORMAT(///39H0STATISTICAL EXCHANGE PARAMETER, ALPHA?,F10.4)
216   FORMAT(///20H0LOUCKS' RADIAL MESH//5(10F11.5/))
217   FORMAT(4A4)
218   FORMAT(3F8.4)
219   FORMAT(2E14.5)
220   FORMAT(10H &NL2 NRR=,i2,5H &END)
221   FORMAT(9H &NL16 Z=,f7.4,4H,RT=,f7.4,5H &END)                                             C
      END
