      SUBROUTINE DBCMAV(N,MA,MB,BCV,IHOM,NRTI,Q,NTI,PHIREC,NU,Z,KKK,
     1                  R,BB,WS)
C     --------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MA(N,N),MB(N,N)
      DIMENSION BCV(N),Q(N,N,NTI),PHIREC(NU,NTI),Z(N,NTI),
     1          R(N,N),BB(N),WS(N,N)
      INTEGER KKK(N)
      LOGICAL HOM0,HOM1
C
      HOM0 = IHOM.EQ.0
      HOM1 = IHOM.NE.0
      DO 1000 I = 1 , N
        KKK(I) = I*(I+1)/2
 1000 CONTINUE
      DO 1010 I = 1 , N
        BB(I) = BCV(I)
 1010 CONTINUE
C     COMPUTATION OF MA.Q(1) STORED IN WS
      DO 1120 K = 1 , N
        CALL DMATVC(MA,N,N,N,N,Q(1,K,1),WS(1,K))
 1120 CONTINUE
      IF (HOM1) THEN
C     COMPUTATION OF MA.Q(1).Z(1) AND BCV-MA.Q(1).Z(1) STORED IN BB.
        DO 1210 I = 1 , N
          SOM = 0.D0
          DO 1220 J = 1 , N
            SOM = SOM + WS(I,J) * Z(J,1)
 1220     CONTINUE
          BB(I) = BB(I) - SOM
 1210   CONTINUE
      ENDIF
C     COMPUTATION OF MA.Q(1).PHIREC(1) STORED IN R
      DO 1420 K = 1 , N
        NK = KKK(K) - K
        DO 1410 I = 1 , N
          SOM = 0.D+0
          DO 1400 J = 1 , K
            NJ = NK + J
            SOM = SOM + WS(I,J) * PHIREC(NJ,1)
 1400     CONTINUE
          R(I,K) = SOM
 1410   CONTINUE
 1420 CONTINUE
C     COMPUTATION OF MB.Q(N) STORED IN WS
      DO 1520 K = 1 , N
        CALL DMATVC(MB,N,N,N,N,Q(1,K,NRTI),WS(1,K))
 1520 CONTINUE
      IF (HOM1) THEN
C     COMPUTATION OF MB.Q(N).Z(N) AND BCV-MB.Q(N).Z(N) STORED IN BB
        DO 1610 I = 1 , N
          SOM = 0.D0
          DO 1600 J = 1 , N
            SOM = SOM + WS(I,J) * Z(J,NRTI)
 1600     CONTINUE
          BB(I) = BB(I) - SOM
 1610   CONTINUE
      ENDIF
C     COMPUTATION OF MB.Q(N).PHIREC(N) AND
C     MA.Q(1).PHIREC(1) + MB.Q(N).PHIREC(N) STORED IN R
      DO 1820 K = 1 , N
        NK = KKK(K) - K
        DO 1810 I = 1 , N
          SOM = 0.D+0
          DO 1800 J = 1 , K
            NJ = NK + J
            SOM = SOM + WS(I,J) * PHIREC(NJ,NRTI)
 1800     CONTINUE
        R(I,K) = R(I,K) + SOM
 1810   CONTINUE
 1820 CONTINUE
      RETURN
C     END OF DBCMAV
      END
      SUBROUTINE DFQUS(PHI,N,IRM,IRN,DIAG,BET,K,KKK,U,NU,NTI,Q,
     1                 WPHI0,WSUM)
C     ----------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PHI(IRM,IRN),DIAG(IRN),BET(IRN),U(NU,NTI),
     1          Q(IRM,IRM,NTI),WPHI0(IRM,IRM),WSUM(IRN)
      INTEGER KKK(N)
C
      CALL DQEVAL(N,N,PHI,IRM,IRN,BET,WPHI0,WSUM)
      DO 1000 I = 1 , N
      DO 1000 J = 1 , N
      Q(I,J,K) = WPHI0(I,J)
 1000 CONTINUE
      U(1,K) = DIAG(1)
      IF (N.LT.2) GOTO 2000
      DO 1200 I = 2 , N
      NI = KKK(I)
      IMIN1 = I - 1
      NJ = NI - I
      U(NI,K) = DIAG(I)
      DO 1100 J = 1 , IMIN1
      NK = NJ + J
      U(NK,K) = PHI(J,I)
 1100 CONTINUE
 1200 CONTINUE
 2000 DO 2300 I = 1 , N
        NI = KKK(I)
        IF (U(NI,K).GE.0.D0) GOTO 2300
        DO 2100 J = 1 , N
          Q(J,I,K) = -Q(J,I,K)
 2100   CONTINUE
        DO 2200 J = I , N
          NJ = KKK(J) - J + I
          U(NJ,K) = -U(NJ,K)
 2200   CONTINUE
 2300 CONTINUE
      RETURN
C     END OF DFQUS
      END
      SUBROUTINE DFUNRC(N,KPART,NRTI,U,NU,NTI,KKK,PHIREC,W1,W2,
     1                  IERROR)
C     ----------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PHIREC(NU,NTI),U(NU,NTI),W1(N),W2(N)
      INTEGER KKK(N)
C
      IERROR = 0
      IF (KPART.LT.0) KPART=0
      IF (KPART.GT.N) KPART=N
      KPPL1 = KPART + 1
C     COMPUTATION OF THE PHI1'S FOR I <= KPART
      IF (KPART.GT.0) THEN
      DO 1400 K = 1 , KPART
      KK = KKK(K)
      NK = KK - K
      NK1 = NK + 1
      DO 1150 I = 1 , K
      PHIREC(NK + I,NRTI) = 0.0D+0
 1150 CONTINUE
      PHIREC(KK,NRTI) = 1.0D+0
      DO 1350 I = 2 , NRTI
      IN = NRTI + 1 - I
      IN1 = IN + 1
      DO 1250 IJ = 1 , K
      W2(IJ) = PHIREC(NK + IJ,IN1)
 1250 CONTINUE
      CALL DSOLUP(K,U(1,IN1),W2,W1,IERROR)
      IF (IERROR.GT.0) GOTO 3000
      DO 1300 IJ = 1 , K
      PHIREC(NK + IJ,IN) = W1(IJ)
 1300 CONTINUE
 1350 CONTINUE
 1400 CONTINUE
      ENDIF
C     THE REQUIRED PART OF PHIREC1 IS NOW COMPUTED.
C     THE REQUIRED PART OF PHIREC2 WILL BE COMPUTED NOW.
      IF (KPPL1.LE.N) THEN
C     COMPUTATION OF THE PHI2'S FOR I > KPART
      DO 1700 K = KPPL1 , N
      KK = KKK(K)
      NK = KK - K
      NK1 = NK + KPPL1
      DO 1500 I = NK1 , KK
      PHIREC(I,1) = 0.0D+0
 1500 CONTINUE
      PHIREC(KK,1) = 1.0D+0
      DO 1650 I = 2 , NRTI
      IN1 = I - 1
      DO 1600 IJ = KPPL1 , K
      SOM = 0.0D+0
      DO 1550 IN = IJ , K
      MK = KKK(IN) - IN
      SOM = SOM + PHIREC(NK + IN , IN1) * U(MK + IJ , I)
 1550 CONTINUE
      PHIREC(NK + IJ , I) = SOM
 1600 CONTINUE
 1650 CONTINUE
 1700 CONTINUE
C     COMPUTATION OF THE PHI1'S FOR I > KPART, IF KPART > 0
      IF (KPART.GT.0) THEN
      DO 2050 K = KPPL1 , N
      KK = KKK(K)
      NK = KK - K
      DO 1750 I = 1, KPART
      PHIREC(NK + I, NRTI) = 0.0D+0
 1750 CONTINUE
      DO 2000 I = 2 , NRTI
      IN = NRTI + 1 - I
      IN1 = IN + 1
      DO 1850 IJ = 1 , KPART
      SOM = 0.0D+0
      DO 1800 NK1 = KPPL1 , K
      MK = KKK(NK1 - 1)
      SOM = SOM + U(MK + IJ,IN1) * PHIREC(NK + NK1 ,IN)
 1800 CONTINUE
      W2(IJ) = - SOM + PHIREC(NK + IJ , IN1)
 1850 CONTINUE
      CALL DSOLUP(KPART,U(1,IN1),W2,W1,IERROR)
      IF (IERROR.GT.0) GOTO 3000
      DO 1950 IJ = 1 , KPART
      PHIREC(NK + IJ , IN) = W1(IJ)
 1950 CONTINUE
 2000 CONTINUE
 2050 CONTINUE
      ENDIF
      ENDIF
      RETURN
 3000 IERROR = 250
      RETURN
C     END OF DFUNRC
      END
      SUBROUTINE DGTUR(N,IHOM,NRTI,U,KU,NTI,Q,D,ER,KKK,KPART,IP1,IP2,
     1                 WDIA,WBET,WSUM,WPHI,WPHI0,IERROR)
C     ---------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(KU,NTI),Q(N,N,NTI),D(N,NTI),ER(5),WDIA(N),
     1          WBET(N),WSUM(N),WPHI(N,N),WPHI0(N,N)
      INTEGER KKK(N),IP1(N),IP2(N)
      LOGICAL NOSORT
C
      IERROR = 0
      ISORT = IP1(1)
      NKP = N
      IF (ISORT.NE.0) NKP = IP1(2)
      IF (NKP.LE.0) THEN
        KPART = 0
        RETURN
      ENDIF
      NTI1 = NTI - 1
      NOSORT = .FALSE.
C.....CHECK WHETHER TRANSFORMATION IS ALLOWED.
      TOL = DMAX1(ER(1),ER(2))
      TOL1 = TOL / ER(3)
      DO 1100 I = 1 , N
        IP1(I) = I
 1100 CONTINUE
      DO 1300 I = 1 , NKP
        IK = KKK(I)
        WDIA(I) = 0.0D+0
        DO 1200 J = 2 , NRTI
          SOM = DABS(U(IK,J))
          IF (SOM.GE.TOL1) NOSORT = .TRUE.
          WDIA(I) = WDIA(I) + DLOG10(SOM)
 1200 CONTINUE
 1300 CONTINUE
      KPART = 0
      DO 1400 I = 1 , NKP
        IF (WDIA(I).GT.0) KPART = KPART + 1
 1400 CONTINUE
C     PARTITIONING IS DETERMINED
      IF ((KPART.EQ.0).OR.(KPART.EQ.NKP)) RETURN
      IF (NOSORT) RETURN
      IFLAG = 0
      K = 0
      N1 = NKP - 1
C     DETERMINATION OF THE PERMUTATION WHICH GLOBALLY ORDERS THE MODES
      DO 1700 I = 1 , N1
        SORT = WDIA(I)
        IF (SORT.LE.0.D0) THEN
          IPL1 = I + 1
          DO 1500 J = IPL1 , NKP
            SOM = WDIA(J)
            IF ((SOM.GT.0.D0).AND.(SOM-SORT.GT.1.D0)) THEN
              IFLAG = 1
              WDIA(J) = SORT
              WDIA(I) = SOM
              JM = IP1(I)
              IP1(I) = IP1(J)
              IP1(J) = JM
              GOTO 1600
            ENDIF
 1500     CONTINUE
 1600     CONTINUE
        ENDIF
 1700 CONTINUE
      IF (IFLAG.EQ.0) RETURN
C.....THE PERMUTATION P IS NOW DETERMINED
C.....COMPUTATION OF Q1P,U2P,INV(P)D1,INV(P)W1
      IERROR = 1
      DO 2100 I = 1 , N
      DO 2100 J = 1 , N
        Q(I,J,NTI) = Q(I,J,1)
        WPHI(I,J) = 0.0D+0
 2100 CONTINUE
      DO 2500 I = 1 , N
        IK = IP1(I)
        IP1(I) = IP2(IK)
        IK1 = KKK(IK) - IK
        DO 2300 J = 1 , N
           Q(J,I,1) = Q(J,IK,NTI)
 2300   CONTINUE
        DO 2400 J = 1 , IK
          WPHI(J,I) = U(IK1+J,2)
 2400   CONTINUE
 2500 CONTINUE
      DO 2600 I = 1 , N
        IP2(I) = IP1(I)
 2600 CONTINUE
C.....TRANSFORMATION OF THE U'S,Q'S,D'S AND WI'S
      DO 3400 I = 2 , NRTI
        CALL DQUDEC(WPHI,N,N,N,N,WDIA,WBET,WSUM)
        CALL DFQUS(WPHI,N,N,N,WDIA,WBET,NTI,KKK,U,KU,NTI,Q,
     1             WPHI0,WSUM)
        DO 2700 J = 1 , KU
          U(J,I) = U(J,NTI)
 2700   CONTINUE
        IF (IHOM.EQ.0) GOTO 2950
        CALL DTAMVC(Q(1,1,NTI),N,N,N,N,D(1,I),WBET)
        DO 2900 J = 1 , N
          D(J,I) = WBET(J)
 2900   CONTINUE
 2950   DO 3000 J = 1 , N
          CALL DMATVC(Q(1,1,I),N,N,N,N,Q(1,J,NTI),WPHI0(1,J))
 3000   CONTINUE
        DO 3100 J = 1 , N
        DO 3100 L = 1 , N
          Q(J,L,I) = WPHI0(J,L)
 3100   CONTINUE
        I1 = I + 1
        IF (I1.GT.NRTI) GOTO 3400
C.....COMPUTATION OF UI+1.OI
        DO 3300 J = 1 , N
        DO 3300 L = 1 , N
          SOM = 0.0D+0
          DO 3200 K = L , N
            IK = KKK(K) - K
            SOM = SOM + U(IK+L,I1) * Q(K,J,NTI)
 3200     CONTINUE
          WPHI(L,J) = SOM
 3300  CONTINUE
 3400 CONTINUE
      RETURN
C     END OF DGTUR
      END
       SUBROUTINE DINPRO(A,NA,B,NB,N,C)
C      --------------------------------
C
       DOUBLE PRECISION A(*),B(*),C
C
C      LOCAL VARIABLES
       INTEGER NN,JL,JJ
       C = 0.0D+0
       IF (N.LE.0) GOTO 2000
       NN = 1 - NB
       JL = NA * N
       DO 1000 JJ = 1 , JL , NA
        NN = NN + NB
        C = C + A(JJ) * B(NN)
 1000  CONTINUE
 2000  CONTINUE
       RETURN
C      END OF DINPRO
C
C      DINPRO IS THE DOUBLE PRECISION VERSION OF SUBROUTINE
C      VVIPP FROM THE ACCU LIBRARY
C
      END
      SUBROUTINE DKPCH(N,U,NU,NTI,NRTI,KPART,KKK,ER,IER)
C     --------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(NU,NTI),ER(5)
      INTEGER KKK(N)
C
      IER = 0
      PSI = DMAX1(ER(1),ER(2)) / ER(3)
      ER(5) = 0.0D+0
      BLOWUQ = 0.0D+0
      KPPL1 = KPART + 1
      IF (KPART.GE.N) GOTO 1300
      DO 1200 K = KPPL1 , N
        PROD = 1.0D+0
        KK = KKK(K)
        DO 1100 I = 2 , NRTI
          PROD = 1.0D+0 + PROD * DABS(U(KK,I))
          IF (BLOWUQ.LT.PROD) BLOWUQ = PROD
 1100   CONTINUE
 1200 CONTINUE
 1300 IF (KPART.LT.1) GOTO 2000
      DO 1500 K = 1 , KPART
        PROD = 1.0D+0
        KK = KKK(K)
        DO 1400 I = 2 , NRTI
          IN = NRTI + 2 - I
          PROD = 1.0D+0 + PROD / DABS(U(KK,IN))
          IF (ER(5).LT.PROD) ER(5) = PROD
 1400   CONTINUE
 1500 CONTINUE
      IF (KPART.GE.N) GOTO 2000
      BLOW1 = 0.0D+0
      BLOW2 = 0.0D+0
      DO 1700 K = 1 , KPART
        KK = KKK(K)
        AMX = 1.0D+0 / DABS(U(KK,NRTI))
        AMX2 = AMX
        IF (NRTI.LT.3) GOTO 1650
        DO 1600 I = 3 , NRTI
        IN = NRTI + 2 - I
        PROD = 1.0D+0 / DABS(U(KK,IN))
        AMX1 = PROD * AMX2
        AMX2 = DMAX1(AMX1,PROD)
        AMX = DMAX1(AMX,AMX2)
 1600   CONTINUE
 1650   IF (AMX.GT.BLOW1) BLOW1 = AMX
 1700 CONTINUE
      DO 1900 K = KPPL1 , N
        KK = KKK(K)
        AMX = DABS(U(KK,2))
        AMX2 = AMX
        DO 1800 I = 3 , NRTI
        PROD = DABS(U(KK,I))
        AMX1 = PROD * AMX2
        AMX2 = DMAX1(PROD,AMX1)
        AMX = DMAX1(AMX,AMX2)
 1800   CONTINUE
        IF (AMX.GT.BLOW2) BLOW2 = AMX
 1900 CONTINUE
      BLOW = BLOW1 * BLOW2
      ER(5) = ER(5) + BLOW
 2000 ER(5) = DMAX1(ER(5),BLOWUQ)
      IF (ER(5).LT.PSI) RETURN
      IER = 240
      RETURN
C     END OF DKPCH
      END
      SUBROUTINE DMATVC(A,M,N,IRM,IRN,Y,X)
C     ------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IRM,IRN),Y(IRN),X(IRM)
C
      IF (M.LE.0.OR.N.LE.0) GOTO 2000
      DO 1100 I=1,M
      S = 0.0D+0
      DO 1000 J=1,N
      S = S + A(I,J) * Y(J)
 1000 CONTINUE
      X(I) = S
 1100 CONTINUE
 2000 CONTINUE
      RETURN
C     END OF DMATVC
      END
      SUBROUTINE DPSR(N,KPART,NRTI,U,NU,NTI,D,KKK,V,W1,W2,IERROR)
C     -----------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(NU,NTI),D(N,NTI),V(N,NTI),W1(N),W2(N)
      INTEGER KKK(N)
C
      IERROR = 0
      KPPL1 = KPART + 1
C     COMPUTATION OF VI2
      IF (KPART.GE.N) GOTO 1500
      DO 1300 I = KPPL1 , N
      V(I,1) = 0.0D+0
 1300 CONTINUE
      DO 1450 I = 2 , NRTI
      IN = I - 1
      DO 1400 K = KPPL1 , N
      SOM = 0.0D+0
      DO 1350 L = K , N
      MK = KKK(L) - L
      SOM = SOM + U(MK+K,I) * V(L,IN)
      MK1 = MK + K
 1350 CONTINUE
      V(K,I) = SOM + D(K,I)
 1400 CONTINUE
 1450 CONTINUE
C     COMPUTATION OF VI1
 1500 IF (KPART.LT.1) RETURN
      DO 1600 I = 1 , KPART
      V(I,NRTI) = 0.0D+0
 1600 CONTINUE
      DO 1850 I = 2 , NRTI
      IN = NRTI + 1 - I
      IN1 = IN + 1
      DO 1700 K = 1 , KPART
      SOM = 0.0D+0
      IF (KPPL1.GT.N) GOTO 1660
      DO 1650 L = KPPL1 , N
      MK = KKK(L) - L
      SOM = SOM + U(MK+K,IN1) * V(L,IN)
 1650 CONTINUE
 1660 W2(K) = - SOM - D(K,IN1) + V(K,IN1)
 1700 CONTINUE
      CALL DSOLUP(KPART,U(1,IN1),W2,W1,IERROR)
      IF (IERROR.NE.0) GOTO 3000
      DO 1800 K = 1 , KPART
      V(K,IN) = W1(K)
 1800 CONTINUE
 1850 CONTINUE
      RETURN
 3000 IERROR = 250
      RETURN
C     END OF DPSR
      END
      SUBROUTINE DQEVAL(M,N,QR,IRM,IRN,BET,QREV,WD)
C     --------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QR(IRM,IRN),QREV(IRM,IRM),BET(IRN),WD(IRM)
C
      DO 1100 I=N,M
      DO 1050 J=I,M
        QREV(I,J) = -QR(I,N) * QR(J,N) * BET(N)
        QREV(J,I) = QREV(I,J)
 1050 CONTINUE
      QREV(I,I) = QREV(I,I) + 1.0D+0
 1100 CONTINUE
      IF (N.LT.2) RETURN
      DO 1300 IK = 2 , N
      I = N + 2 - IK
C     THE NECESARY INNERPRODUCTS OF THE UK'S ARE STORED IN THE
C     ARRAY WD(1:IRM)
C
      IMIN1 = I - 1
      WD(IMIN1) = QR(IMIN1,IMIN1)
      IM = M - I + 1
      DO 1150 J = I , M
        CALL DINPRO(QR(I,IMIN1),1,QREV(I,J),1,IM,WD(J))
C
C     COMPUTATION OF (I - BET(I-1).U(I-1)') * PRODUCT
C     ( ( -BET(K).U(K).U(K)') K= I TO N)
C
 1150 CONTINUE
      DO 1200 K = I , M
      DO 1200 L = I , M
        QREV(K,L) = QREV(K,L) - BET(IMIN1) * QR(K,IMIN1) * WD(L)
 1200 CONTINUE
      DO 1250 J = I , M
        QREV(J,IMIN1) = - BET(IMIN1) * QR(J,IMIN1) * WD(IMIN1)
        QREV(IMIN1,J) = - BET(IMIN1) * QR(IMIN1,IMIN1) * WD(J)
 1250 CONTINUE
      QREV(IMIN1,IMIN1) = 1.0D+0 - BET(IMIN1)*WD(IMIN1)*WD(IMIN1)
 1300 CONTINUE
      RETURN
C     END OF DQEVAL
      END
      SUBROUTINE DQUDEC(A,M,N,IRM,IRN,DIAG,BET,WY)
C     ----------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IRM,IRN),DIAG(IRN),BET(IRN),WY(IRN)
C
      IF (M.GT.IRM) GOTO 3000
      IF (N.GT.IRN) GOTO 3010
C
      DO 2000 K = 1 , N
       K1 = K + 1
       IK = M - K + 1
       CALL DINPRO(A(K,K),1,A(K,K),1,IK,SIGMA)
       IF (SIGMA.EQ.0.D0) GOTO 1400
       QRKK = A(K,K)
       ALPHAK = DSQRT(SIGMA)
       IF (QRKK.GE.0.0D+0) ALPHAK = - ALPHAK
       DIAG(K) = ALPHAK
       BETA = 1.0D+0 / (SIGMA - QRKK * ALPHAK)
       BET(K) = BETA
       A(K,K) = QRKK - ALPHAK
       IF (K1.GT.N) GOTO 2000
       DO 1200 J = K1 , N
        CALL DINPRO(A(K,K),1,A(K,J),1,IK,WY(J))
        WY(J) = WY(J) * BETA
 1200  CONTINUE
       DO 1300 J = K1 , N
       DO 1300 I = K  , M
        A(I,J) = A(I,J) - A(I,K) * WY(J)
 1300  CONTINUE
       GOTO 2000
 1400  BET(K) = 0.D0
       DIAG(K) = 0.D0
 2000 CONTINUE
      RETURN
 3000 WRITE(*,100)
      RETURN
 3010 WRITE(*,110)
      RETURN
C
  100 FORMAT(' ERROR: M > IRM ; NO DECOMPOSITION DONE')
  110 FORMAT(' ERROR: N > IRN ; NO DECOMPOSITION DONE')
C
C     END OF DQUDEC
      END
      SUBROUTINE DSBVP(N,IHOM,M0,MN,NRTI,Q,NTI,BCV,PHIREC,NU,V,
     1                 WI,EPS,COND,KKK,IP,WQ,WS,W1,W2,IER)
C    ------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION M0(N,N),MN(N,N)
      DIMENSION Q(N,N,NTI),BCV(N),PHIREC(NU,NTI),V(N,NTI),
     1          WI(N,NTI),WQ(N,N),WS(N,N),W1(N),W2(N)
      INTEGER KKK(N),IP(N)
C
      CALL DBCMAV(N,M0,MN,BCV,IHOM,NRTI,Q,NTI,PHIREC,NU,WI,KKK,WQ,
     1            W1,WS)
      CALL DCROUT(WQ,N,N,W1,EPS,IP,W2,IER)
      IF (IER.NE.0) THEN
       IER = 260
       RETURN
      ENDIF
C     COMPUTATION OF THE CONDITION NUMBER
      COND = 0.0D+0
      DO 1300 I = 1 , N
        DO 1100 J = 1 , N
          W2(J) = 0.0D+0
 1100   CONTINUE
        W2(I) = 1.0D+0
        CALL DSOLDE(WQ,N,N,IP,W2)
        SOMIN = 0.0D+0
        DO 1200 J = 1 , N
          SOMIN = SOMIN + DABS(W2(J))
 1200   CONTINUE
        IF (SOMIN.GT.COND) COND = SOMIN
 1300 CONTINUE
C     COMPUTATION OF EI'=PHIREC(I)Y' AND EI=Q(I)EI'
      DO 1500 K=1,NRTI
      DO 1500 I=1,N
        SOMIN=0.0D+0
        DO 1400 J=I,N
          NJ=KKK(J) - J
          NI= NJ + I
          SOMIN = SOMIN + PHIREC(NI,K) * W1(J)
 1400   CONTINUE
        IF (IHOM.EQ.0) THEN
          WI(I,K) = SOMIN
        ELSE
          WI(I,K) = WI(I,K) + SOMIN
        ENDIF
 1500 CONTINUE
C     COMPUTATION OF Q(I)XI'
      DO 2000 K=1,NRTI
        CALL DMATVC(Q(1,1,K),N,N,N,N,WI(1,K),W2)
        DO 1800 J=1,N
          V(J,K) = W2(J)
 1800   CONTINUE
 2000 CONTINUE
      RETURN
C     END OF DSBVP
      END
      SUBROUTINE DSOLUP(M,U,B,X,IER)
C     ------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(*),X(*),B(*)
C
      IER = 0
      IF (M.LE.0) RETURN
      M1 = M * (M+1) / 2
      IF (U(M1).EQ.0.0D+0) GOTO 1200
      X(M) = B(M) / U(M1)
      IF (M.EQ.1) RETURN
      MMIN1 = M - 1
      DO 1100 I = 1 , MMIN1
      MI = M - I
      SOMIN = 0.0D+0
      MJ = (MI + 1) * MI / 2
      MIJ = MJ
      DO 1000 J = 1 , I
      MIJ = MIJ + MI - 1 + J
      SOMIN = SOMIN + U(MIJ) * X(MI+J)
 1000 CONTINUE
      IF (U(MJ).EQ.0.0D+0) GOTO 1200
      X(MI) = (B(MI) - SOMIN) / U(MJ)
 1100 CONTINUE
      RETURN
 1200 IER = 1
      RETURN
C     END OF DSOLUP
      END
      SUBROUTINE DUPUP(N,U1,U2,U3)
C     ----------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U1(*),U2(*),U3(*)
C
      DO 2000 I = 1 , N
      IK = I * (I+1) / 2 - I
        DO 2000 J = 1 , I
          SOM = 0.0D+0
          DO 1000 L = J , I
            IL = L * (L+1) / 2 - L
            SOM = SOM + U1(IL+J) * U2(IK+L)
 1000     CONTINUE
          U3(IK+J) = SOM
 2000 CONTINUE
      RETURN
C     END OF DUPUP
      END
C
      SUBROUTINE DAMTES(AMP,WDIA,N,IFLAG)
C     -----------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WDIA(*)
C
      AMAXU = 0.0D+0
      DO 1000 I = 1 , N
      ABSDIA = DABS(WDIA(I))
      IF (ABSDIA.GT.AMAXU) AMAXU = ABSDIA
 1000 CONTINUE
      IF (AMAXU.GT.AMP) GOTO 1020
      IFLAG = 0
      RETURN
 1020 IFLAG = 1
      RETURN
      END
      SUBROUTINE DTAMVC(A,M,N,IM,IN,Y,X)
C     ----------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(IM,IN),Y(IM),X(IN)
C
      IF (M.LE.0 .OR. N.LE.0) RETURN
      DO 1100 I = 1 , N
        SOM = 0.D0
        DO 1000 J = 1 , M
          SOM = SOM + A(J,I) * Y(J)
 1000   CONTINUE
        X(I) = SOM
 1100 CONTINUE
      RETURN
C     END OF DTAMVC
      END
      SUBROUTINE DCDI(D,N,NTI,Q,K,WSUM)
C     ---------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION D(N,NTI),Q(N,N,NTI),WSUM(N)
C
      CALL DTAMVC(Q(1,1,K),N,N,N,N,D(1,K),WSUM)
      DO 1000 I = 1 , N
        D(I,K) = WSUM(I)
 1000 CONTINUE
      RETURN
      END
      SUBROUTINE DCNRHS(DIAG,N,NRHS)
C     ------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DIAG(N)
C
      A = DABS(DIAG(1))
      NRHS = 1
      DO 1000 I = 2 , N
        B = DABS(DIAG(I))
        IF (B.LT.1.D0) B = 1.D0 / B
        IF (B.LE.A) GOTO 1000
        A = B
        NRHS = I
 1000 CONTINUE
      RETURN
C     END OF DCNRHS
      END
      SUBROUTINE DCPHIS(FLIN,FDIF,N,A,B,ER,K,IHOM,NRHS,WY,WI,Q,NTI,W,
     1                  WPHI,WPHI0,IER)
C     ---------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ER(5),WY(N),WI(N,NTI),Q(N,N,NTI),
     1 W(N,7),WPHI(N,N),WPHI0(N,N)
      EXTERNAL FLIN,FDIF
C     LOCAL ARRAYS
      DIMENSION HI(5),HI2(5)
C
C     COMPUTATION OF THE FUNDAMENTAL SOLUTION
      KMIN1 = K - 1
      NTI1 = NTI - 1
      IF ((NRHS.NE.0).OR.(IHOM.NE.0)) GOTO 1600
C     NRHS=0 AND IHOM=0 , DETERMONATION OF NRHS ON THE FIRST MINOR
C     SHOOTING INTERVAL
      DO 1500 L = 1 , N
        T1 = A
        T2 = B
        DO 1100 I = 1 , N
         WY(I) = Q(I,L,1)
 1100   CONTINUE
        IER = 1
        CALL DRKFSM(FLIN,N,WY,T1,T2,ER(1),ER(2),ER(3),HI2,NHI2,W,IER)
        IF (IER.GT.3) RETURN
        DEL2 = (T2 - A) / NHI2
        IF (L.EQ.1) GOTO 1200
        IF (DABS(DEL2).GT.DABS(DEL1)) GOTO 1500
 1200   DEL1 = DEL2
        NRHS = L
        NHI = NHI2
        DO 1300 I = 1 , NHI
          HI(I) = HI2(I)
 1300   CONTINUE
        T3 = T2
        DO 1400 I = 1 , N
          WPHI(I,NRHS) = WY(I)
          WPHI0(I,NRHS)= WY(I)
 1400   CONTINUE
 1500 CONTINUE
      T1 = A
      T2 = B
      B = T3
      GOTO 2400
 1600 IF (IHOM.EQ.0) GOTO 1900
C     INHOMOGENEOUS SYSTEM, COMPUTATION OF THE PARTICULAR SOLUTION
      DO 1700 I = 1 , N
        WY(I) = 0.D0
 1700 CONTINUE
      IER = 1
      T1 = A
      CALL DRKFSM(FDIF,N,WY,T1,B,ER(1),ER(2),ER(3),HI,NHI,W,IER)
      IF (IER.GT.3) RETURN
      DO 1800 I = 1 , N
        WI(I,K) = WY(I)
 1800 CONTINUE
      NRHS = 0
      GOTO 2400
 1900 DO 2000 I = 1 , N
        WY(I) = Q(I,NRHS,KMIN1)
 2000 CONTINUE
      IER = 1
      T1 = A
      CALL DRKFSM(FLIN,N,WY,T1,B,ER(1),ER(2),ER(3),HI,NHI,W,IER)
      IF (IER.GT.3) RETURN
      IF (KMIN1.GT.1) GOTO 2200
      DO 2100 I = 1 , N
        WPHI(I,NRHS) = WY(I)
        WPHI0(I,NRHS)= WY(I)
 2100 CONTINUE
      GOTO 2400
 2200 DO 2300 I = 1 , N
        WPHI(I,NRHS) = WY(I)
 2300 CONTINUE
      GOTO 2400
 2400 CONTINUE
C     COMPUTATION OF THE BASIC HOMOGENEOUS SOLUTIONS, WHICH HAVE NOT
C     BEEN COMPUTED YET
      DO 2900 I = 1 , N
        IF (I.EQ.NRHS) GOTO 2900
        DO 2500 J = 1 , N
          WY(J) = Q(J,I,KMIN1)
 2500   CONTINUE
        T1 = A
        CALL DRKFGS(FLIN,N,WY,T1,HI,NHI,W)
        IF (KMIN1.GT.1) GOTO 2700
        DO 2600 J = 1 , N
          WPHI(J,I) = WY(J)
          WPHI0(J,I) = WY(J)
 2600   CONTINUE
        GOTO 2900
 2700   DO 2800 J = 1 , N
          WPHI(J,I) = WY(J)
 2800   CONTINUE
 2900 CONTINUE
      CALL DQUDEC(WPHI,N,N,N,N,W(1,1),W(1,2),WY)
      RETURN
C     END OF DCPHIS
      END
      SUBROUTINE DCROUT(A,N,M,B,EPS,P,V,IER)
C     --------------------------------------
C
      DOUBLE PRECISION A(M,*),B(*),V(*),EPS
      INTEGER P(*)
C
      IER = 0
      CALL DLUDEC(A,N,M,EPS,P,V,IER)
      IF (IER.NE.0) RETURN
      CALL DSOLDE(A,N,M,P,B)
      RETURN
C     END OF DCROUT
C
C     DCROUT IS AN DOUBLE PRECISION FORTRAN VERSION OF
C     T.DEKKER'S SOLDEC
      END
      SUBROUTINE DINTCH(A,NA,B,NB,N)
C     ------------------------------
C
      DOUBLE PRECISION A(*),B(*),C
C
      IF (N.LE.0) RETURN
      NN = 1 - NB
      JL = NA * N
      DO 1000 JJ = 1 , JL , NA
        NN = NN + NB
        C = B(NN)
        B(NN) = A(JJ)
        A(JJ) = C
 1000 CONTINUE
      RETURN
C     END OF DINTCH
C
C     DINTCH IS A DOUBLE PRECISION VERSION OF
C     J.M. VAN KATS SUBROUTINE INTCHA
      END
      SUBROUTINE DLUDEC(A,N,M,EPS,P,V,IER)
C     ------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,*),V(*)
      INTEGER P(*),PK
C
      IER = 0
      R = -1.0D+0
      Z = 1.0D+0
      DO 1000 I = 1 , N
        CALL DINPRO(A(I,1),M,A(I,1),M,N,S)
        S = DSQRT(S)
        IF (S.GT.R) R = S
        V(I) = 1.0D+0 / S
 1000 CONTINUE
      E = EPS * R
      DO 5000 K = 1 , N
        R = -1.0D+0
        K1 = K - 1
        DO 2000 I = K , N
          CALL DINPRO(A(I,1),M,A(1,K),1,K1,S)
          A(I,K) = A(I,K) - S
          S = DABS(A(I,K)) * V(I)
          IF (S.LE.R) GOTO 2000
          R = S
          PK = I
 2000   CONTINUE
        P(K) = PK
        V(PK) = V(K)
        Z = A(PK,K)
        IF (DABS(Z).LT.E) GOTO 7000
        IF (PK.EQ.K) GOTO 3000
        CALL DINTCH(A(K,1),M,A(PK,1),M,N)
 3000   KK = K + 1
        IF (K.EQ.N) RETURN
        DO 4000 I = KK , N
          CALL DINPRO(A(K,1),M,A(1,I),1,K1,S)
          A(K,I) = (A(K,I) - S) / Z
 4000   CONTINUE
 5000 CONTINUE
 7000 IER = 1
      RETURN
C     END OF DLUDEC
C
C     DLUDEC IS AN FORTRAN VERSION OF T. DEKKER'S DEC
      END
      SUBROUTINE DSOLDE(A,N,M,P,B)
C     ----------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,*),B(*)
      INTEGER P(*),PK
C
      DO 1000 K = 1 , N
        R = B(K)
        PK = P(K)
        CALL DINPRO(A(K,1),M,B,1,K-1,C)
        B(K) = (B(PK) - C) / A(K,K)
        IF (PK.NE.K) B(PK) = R
 1000 CONTINUE
      IF (N.EQ.1) RETURN
      K = N
 2000 K = K - 1
      CALL DINPRO(A(K,K+1),M,B(K+1),1,N-K,C)
      B(K) = B(K) - C
      IF (K.GT.1) GOTO 2000
      RETURN
C     END OF DSOLDE
C
C     DSOLDE IS DOUBLE PRECISION VERSION OF H.ZWEERUS-VINK'S
C     AND J.M. VAN KAT'S SUBROUTINE SOLDEC
      END
