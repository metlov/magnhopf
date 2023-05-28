C     TEST PROGRAMMA VOOR NIET LINEAIRE RANDWAARDE PROBLEEM
C
C          EXAMPLE  3 FROM SIMINARBERICHT 67
C             HANKE, LAMOUR, WINKLER
C          SEKTION MATHEMATIK, HUMBOLDT-UNIVERSITAT, BERLIN
C
C
C          X1' = A.X1/X2.(X3-X1)
C          X2' = -A.(X3-X1)
C          X3' = 1/X4.(B-C.(X3-X5)-A.X3(X3-X1))
C          X4' = A.(X3-X1)
C          X5' = -C/D.(X5-X3)
C
C          X1(0) = X2(0) = X3(0) = 1 ; X4(0) = -10
C          X3(1) = X5(1)
C
C          A = 0.5 , B = 0.9 , C = 1000 , D = 10.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ER(5),TI(12),X(5,12),Q(5,5,12),U(15,12),D(5,12),
     1          PHIREC(15,12),W(315),WGR(20)
      INTEGER IW(27)
      EXTERNAL FDIF,X0T,G
C
      N = 5
      NU = 15
      NTI = 12
      LW = 315
      LIW = 27
      ER(3) = 1.1D-15
      LWG = 20
      ER(1) = 1.d-6
      ER(2) = 1.D-2
      A = 0.D0
      B = 1.D0
      NRTI=10
      AMP=100
      ITLIM=20
      IERROR=1
      CALL MUSN(FDIF,X0T,G,N,A,B,ER,TI,NTI,NRTI,AMP,ITLIM,X,Q,U,NU,D,
     1          PHIREC,KPART,W,LW,IW,LIW,WGR,LWG,IERROR)
      WRITE(*,*) ' MUSN: IERROR =',IERROR
      WRITE(*,200) A,B,ER(1),ER(2),ER(4),ER(5),KPART
  200 FORMAT(' A = ',F8.4,3X,'B = ',F8.4,/,' REQUIRED TOLERANCE = ',1P,
     1 D12.5,3X,'START TOLERANCE  = ',D12.5,/,
     2 ' CONDITION NUMBER = ',D12.5,3X,
     3 'AMPLIFICATION FACTOR = ',D12.5,/,' K-PARTITIONING =',I2,/)
      IF (IERROR.NE.0) GOTO 3000
      WRITE(*,215)
  215 FORMAT(' I ',4X,'T',9X,'Y1',12X,'Y2',12X,'Y3',12X,'Y4',12X,'Y5',
     1 /)
      DO 2200 K = 1 , NRTI
        WRITE(*,220) K,TI(K),(X(J,K),J=1,N)
 2200 CONTINUE
  220 FORMAT(' ',I2,1X,F6.4,1P,5(2X,D12.5))
 3000 CONTINUE
      STOP
      END
      SUBROUTINE FDIF(T,Y,F)
C     ----------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(5),F(5)
C
      Y3MY1 = Y(3) - Y(1)
      Y3MY5 = Y(3) - Y(5)
      F(1) = 0.5D0 * Y(1) * Y3MY1 / Y(2)
      F(2) = - 0.5D0 * Y3MY1
      F(3) = (0.9D0 - 1.D3 * Y3MY5 - 0.5D0 * Y(3) * Y3MY1) / Y(4)
      F(4) = 0.5D0 * Y3MY1
      F(5) = 1.D2 * Y3MY5
      RETURN
C     END OF FDIF
      END
      SUBROUTINE X0T(T,X)
C     -------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(5)
C
      X(1) = 1.D0
      X(2) = 1.D0
      X(4) = -10.D0
      X(3) = - 4.5D0*T*T + 8.91D0 * T + 1.D0
      X(5) = - 4.5D0*T*T + 9.D0 * T + 0.91D0
      RETURN
C     END OF X0T
      END
      SUBROUTINE G(N,XA,XB,FG,DGA,DGB)
C     --------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XA(N),XB(N),FG(N),DGA(N,N),DGB(N,N)
C
      DO 1100 I = 1 , N
      DO 1100 J = 1 , N
        DGA(I,J) = 0.D0
        DGB(I,J) = 0.D0
 1100 CONTINUE
      DGA(1,1) = 1.D0
      DGA(2,2) = 1.D0
      DGA(3,3) = 1.D0
      DGA(4,4) = 1.D0
      DGB(5,3) = 1.D0
      DGB(5,5) = -1.D0
      FG(1) = XA(1) - 1.D0
      FG(2) = XA(2) - 1.D0
      FG(3) = XA(3) - 1.D0
      FG(4) = XA(4) + 10.D0
      FG(5) = XB(3) - XB(5)
      RETURN
C     END OF G
      END
