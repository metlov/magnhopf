C     THIS PROGRAM COMPUTES THE SOLUTION OF THE EXAMPLE PROBLEM GIVEN
C     IN THE DOCUMENTATION FOR SUBROUTINE MUSL.
C
      DOUBLE PRECISION A,B,MA(3,3),MB(3,3),BCV(3),AMP,ER(5),TI(15),
     1 X(3,15),U(6,15),Q(3,3,15),D(3,15),PHIREC(6,15),W(42),
     2 EXSOL,AE
      INTEGER IW(9)
      EXTERNAL FLIN,FDIF
C
C     SETTING OF THE INPUT PARAMETERS
C
      N = 3
      IHOM = 1
      ER(1) = 1.D-11
      ER(2) = 1.D-6
      ER(3) = 1.1D-16
      NRTI = 10
      NTI = 15
      NU = 6
      LW = 42
      LIW = 9
      A = 0.D0
      B = 6.D0
      AMP = 0.D0
C
C     SETTING THE BC MATRICES MA AND MB
C
      DO 1100 I = 1 , N
        DO 1000 J = 1 , N
          MA(I,J) = 0.D0
          MB(I,J) = 0.D0
 1000   CONTINUE
        MA(I,I) = 1.D0
        MB(I,I) = 1.D0
 1100 CONTINUE
C
C     SETTING THE BC VECTOR BCV
C
      BCV(1) = 1.D0 + DEXP(6.D0)
      BCV(2) = BCV(1)
      BCV(3) = BCV(1)
C
C     CALL MUSL
C
      CALL MUSL(FLIN,FDIF,N,IHOM,A,B,MA,MB,BCV,AMP,ER,NRTI,TI,NTI,
     1          X,U,NU,Q,D,KPART,PHIREC,W,LW,IW,LIW,IERROR)
      IF ((IERROR.NE.0).AND.(IERROR.NE.200).AND.(IERROR.NE.213).AND.
     1 (IERROR.NE.240)) GOTO 5000
C
C     COMPUTATION OF THE ABSOLUTE ERROR IN THE SOLUTION AND WRITING
C     OF THE SOLUTION AT THE OUTPUTPOINTS
C
      WRITE(*,200)
      WRITE(*,190) ER(4),ER(5)
      WRITE(*,210)
      WRITE(*,200)
      DO 1500 K = 1 , NRTI
        EXSOL = DEXP(TI(K))
        AE = EXSOL - X(1,K)
        WRITE(*,220) K,TI(K),X(1,K),EXSOL,AE
        DO 1300 I = 2 , N
          AE = EXSOL - X(I,K)
          WRITE(*,230) X(I,K),EXSOL,AE
 1300   CONTINUE
 1500 CONTINUE
      STOP
 5000 WRITE(*,300) IERROR
      STOP
C
  190 FORMAT(' CONDITION NUMBER     = ',D10.3,/,
     1       ' AMPLIFICATION FACTOR = ',D10.3,/)
  200 FORMAT(' ')
  210 FORMAT('  I ',6X,'T',8X,'APPROX. SOL.',9X,'EXACT SOL.',8X,
     1 'ABS. ERROR')
  220 FORMAT(' ',I3,3X,F7.4,3(3X,D16.9))
  230 FORMAT(' ',13X,3(3X,D16.9))
  300 FORMAT(' TERMINAL ERROR IN MUSL: IERROR = ',I4)
C
      END
C
      SUBROUTINE FLIN(T,Y,F)
C     ----------------------
C
      DOUBLE PRECISION T,Y(3),F(3)
      DOUBLE PRECISION TI,SI,CO
C
      TI = 2.D0 * T
      SI = 2.D0 * DSIN(TI)
      CO = 2.D0 * DCOS(TI)
      F(1) = (1.D0 - CO) * Y(1) + (1.D0 + SI) * Y(3)
      F(2) = 2.D0 * Y(2)
      F(3) = (-1.D0 + SI) * Y(1) + (1.D0 + CO) * Y(3)
C
      RETURN
C     END OF FLIN
      END
C
      SUBROUTINE FDIF(T,Y,F)
C     ----------------------
C
      DOUBLE PRECISION T,Y(3),F(3)
      DOUBLE PRECISION TI,SI,CO
C
      CALL FLIN(T,Y,F)
      TI = 2.D0 * T
      SI = 2.D0 * DSIN(TI)
      CO = 2.D0 * DCOS(TI)
      TI = DEXP(T)
      F(1) = F(1) + (-1.D0 + CO - SI)*TI
      F(2) = F(2) - TI
      F(3) = F(3) + (1.D0 - CO - SI)*TI
C
      RETURN
C     END OF FDIF
      END
