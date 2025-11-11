      SUBROUTINE POISON(PSQ,Z,J,W)
      DIMENSION PSQ(J),W(J)
      DOUBLE PRECISION E(250),F(250),ACC,A,B,C,D,C2
C TAKEN FROM LOUCKS' BOOK, APPENDIX 1
      A=1.0D0-0.0025D0/48.0D0
C EQ. A1.11
      B=-2.0D0-0.025D0/48.0D0
C EQ. A1.12
      C=0.0025D0/6.0D0
      D=DEXP(0.025D0)
      C2=-B/A
      E(1)=0.0D0
C EQ. A1.29
      F(1)=D
C EQ.A1.30
      X=-8.75
      J1=J-1
      DO 1 I=2,J1
      ACC=C*EXP(0.5*X)*(D*PSQ(I+1)+10.0*PSQ(I)+PSQ(I-1)/D)
C EQS. A1.13, A1.6
      F(I)=C2-1.0/F(I-1)
C EQ. A1.20
      E(I)=(ACC/A+E(I-1))/F(I)
C EQ. A1.21
1     X=X+0.05
      W(J)=2.0*Z*EXP(-0.5*X)
      ACC=W(J)
      DO 2 I=1,J1
      JC=J-I
      ACC=E(JC)+ACC/F(JC)
2     W(JC)=ACC
C EQ.A1.15
      RETURN
      END
