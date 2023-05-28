C                         ******************
C**************************       MUS      *****************************
C                         ******************
C
C                      DOUBLE PRECISION VERSION
C
C
C     AUTHORS : R.M.M. MATTHEIJ , G.W.M. STAARINK.
C
C     RELEASE DATE: FEB. 15 1988
C
C     LAST UPDATE: JUNE 1992
C
C     ADDRESS : G.W.M. STAARINK
C               SINGEL 90
C               6584 BK MOLENHOEK (LB)
C               THE NETHERLANDS
C
C               R.M.M. MATTHEIJ
C               DEPT. OF MATHEMATICS AND COMPUTING SCIENCE
C               EINDHOVEN UNIVERSITY OF TECHNOLOGY
C               P.O. BOX 513
C               5600 MB EINDHOVEN
C               THE NETHERLANDS
C
C     REFERENCES: R.M.M. MATTHEIJ, G.W.M. STAARINK, BOUNDPAK, NUMERICAL
C                 SOFTWARE FOR LINEAR BVP, EUT REPORT 92-WSK-01 (1992).
C                 R.M.M. MATTHEIJ, G.W.M. STAARINK, IMPLEMENTING
C                 MULTIPLE SHOOTING FOR NONLINEAR BVP, EUT, RANA 87-14,
C                 (1987).
C
C***********************************************************************
C
      SUBROUTINE MUSL(FLIN,FDIF,N,IHOM,A,B,MA,MB,BCV,AMP,ER,NRTI,TI,
     1                NTI,X,U,NU,Q,D,KPART,PHIREC,W,LW,IW,LIW,IERROR)
C     -----------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MA(N,N),MB(N,N)
      DIMENSION BCV(N),ER(5),TI(NTI),X(N,NTI),U(NU,NTI),Q(N,N,NTI),
     1          D(N,NTI),PHIREC(NU,NTI),W(LW)
C
      INTEGER IW(LIW)
C
      EXTERNAL FLIN,FDIF
C
      I2 = 1 + N
      I3 = I2 + N
      J2 = 1 + N
      J3 = J2 + 7*N
      J4 = J3 + N*N
      IW(I2) = 0
      IERROR = 0
      I = N * (N+1) / 2
      IF ((N.LT.1).OR.(IHOM.LT.0).OR.(NRTI.LT.0).OR.(NTI.LT.5).OR.
     1 (NU.LT.I).OR.(A.EQ.B)) THEN
        IERROR = 100
        GOTO 2000
      ENDIF
      IF ((ER(1).LT.0.D0).OR.(ER(2).LT.0.D0).OR.(ER(3).LT.0.D0)) THEN
        IERROR = 101
        GOTO 2000
      ENDIF
      IF ((LW.LT.8*N+2*N*N).OR.(LIW.LT.3*N)) THEN
        IERROR = 103
        GOTO 2000
      ENDIF
C     SETTING Q(I) = I.
      DO 1050 I = 1 , N
      DO 1000 J = 1 , N
        Q(I,J,1) = 0.0D+0
 1000 CONTINUE
        Q(I,I,1) = 1.D+0
 1050 CONTINUE
      CALL DDUR(FLIN,FDIF,N,IHOM,A,B,NRTI,AMP,TI,NTI,ER,Q,U,NU,D,D,
     1          KPART,W(1),W(J2),W(J3),W(J4),IW(1),IW(I2),IW(I3),IERROR)
      IF ((IERROR.NE.0).AND.(IERROR.NE.200).AND.(IERROR.NE.213))
     1  GOTO 2000
      IF (IERROR.NE.0) THEN
        WRITE(*,110) IERROR
        CALL ERRHAN(IERROR,ER,0)
      ENDIF
      IW(I2) = 0
      J5 = J2+N
      CALL DGTUR(N,IHOM,NRTI,U,NU,NTI,Q,D,ER,IW(1),KPART,IW(I2),
     1           IW(I3),W(1),W(J2),W(J5),W(J3),W(J4),IERROR)
      CALL DKPCH(N,U,NU,NTI,NRTI,KPART,IW(1),ER,IERROR)
      IF (IERROR.NE.0) THEN
        WRITE(*,110) IERROR
        CALL ERRHAN(IERROR,ER,0)
      ENDIF
      CALL DFUNRC(N,KPART,NRTI,U,NU,NTI,IW(1),PHIREC,W(1),W(J2),IERROR)
      IF (IERROR.NE.0) GOTO 2000
      IF (IHOM.NE.0) THEN
      CALL DPSR(N,KPART,NRTI,U,NU,NTI,D,IW(1),X,W(1),W(J2),IERROR)
      IF (IERROR.NE.0) GOTO 2000
      ENDIF
      CALL DSBVP(N,IHOM,MA,MB,NRTI,Q,NTI,BCV,PHIREC,NU,X,X,ER(3),
     1           ER(4),IW(1),IW(I2),W(J3),W(J4),W(1),W(J2),IERROR)
      IF (IERROR.NE.0) GOTO 2000
      RETURN
 2000 WRITE(*,100) IERROR
      CALL ERRHAN(IERROR,ER,0)
      RETURN
  100 FORMAT(' TERMINAL ERROR IN MUSL : IERROR =',I4)
  110 FORMAT(' WARNING ERROR IN MUSL : IERROR =',I4)
C     END OF MUSL
      END
      SUBROUTINE MUSN(FDIF,X0T,G,N,A,B,ER,TI,NTI,NRTI,AMP,ITLIM,X,Q,
     1                U,NU,D,PHI,KP,W,LW,IW,LIW,WGR,LWG,IERROR)
C     ---------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ER(5),TI(NTI),X(N,NTI),Q(N,N,NTI),U(NU,NTI),D(N,NTI),
     1          PHI(NU,NTI),W(LW),ER1(5),WGR(LWG)
      INTEGER IW(LIW)
      EXTERNAL FDIF,X0T,G
      DIMENSION ALINC(4)
      LOGICAL DIAGNO
C
      DIAGNO = IERROR.EQ.1
      IERROR = 0
C     CHECK INPUT PARAMETERS
C     INPUT ERROR 105
      IF ((N.LT.1).OR.(NRTI.LT.0).OR.(NTI.LT.3).OR.(NU.LT.N*(N+1)/2)
     1 .OR.(A.EQ.B)) GOTO 5000
C     INPUT ERROR 101
      IF ((ER(1).LT.0.D0).OR.(ER(2).LT.0.D0).OR.(ER(3).LT.0.D0))
     1 GOTO 5010
C     INPUT ERROR 106
      K = N * NTI
      L = N * N
      IF ((LW.LT.(7*N+3*K+4*L)).OR.(LIW.LT.3*N+NTI)) GOTO 5020
C     SETTING POINTERS FOR THE WORK ARRAYS W AND IW
      I1 = 7 * N + 1
      I2 = I1 + K
      I3 = I2 + K
      I4 = I3 + K
      I5 = I4 + L
      I6 = I5 + L
      I7 = I6 + L
      J1 = N + 1
      J2 = J1 + N
      J3 = J2 + N
      TOL = ER(1)
      ER1(3) = ER(3)
      ER1(1) = 1.1D-12 + 3.D0 * ER(3)
      TOL1 = DMAX1(ER(1),DMIN1(ER(2),1.D-2))
      ER1(2) = TOL1
      EPS = DSQRT(ER(3))
      ALINCR = DSQRT(ER(1)/ER(3)) / 4.D0
      IF ((AMP.LE.1.D0).OR.(AMP.GT.ALINCR)) THEN
        ALINC(1) = ALINCR
      ELSE
        ALINC(1) = AMP
      ENDIF
      ALINC(2) = ALINC(1) * ALINC(1)
      ALINC(3) = ALINC(1) * EPS
      ALINC(4) = ALINC(3) * ALINC(1)
      IF (DIAGNO) WRITE(*,200) TOL1
      ITER = 0
      ITTOT = 0
      IHL = 0
      ALAMIN = 0.007D0
 1000 CALL DCSHPO(FDIF,X0T,N,A,B,ER1,TI,NTI,NRTI,ALINC,X,W(I1),Q,U,NU,
     1            W(1),W(I2),W(I4),W(I5),IW(1),IW(J1),WGR,LWG,IERROR)
      NEG = IW(J1)
      IF ((IERROR.NE.0).AND.(IERROR.NE.213)) GOTO 6000
      IF (IERROR.NE.0) THEN
        WRITE(*,110) IERROR
        CALL ERRHAN(IERROR,ER,NEG)
      ENDIF
      NOSH = NRTI
      JAC = 1
      NG = 0
      CALL DCSAOJ(FDIF,N,ER1,TI,NTI,NRTI,X,JAC,NG,IW(1),W(I1),Q,U,NU,
     1            W(1),W(I4),W(I5),IW(J1),WGR,LWG,IERROR)
      NEG = IW(J1)
      IF ((IERROR.NE.0).AND.(IERROR.NE.213)) GOTO 6000
      IF (DIAGNO) THEN
        WRITE(*,210)
        DO 1050 I = 1 , NRTI , 5
          II = MIN0(NRTI,I+4)
          WRITE(*,220) (TI(J),J=I,II)
 1050   CONTINUE
      ENDIF
      JN = I1 + N * NRTI -1
 1100 CALL DJINGX(G,N,NRTI,X,W(I1),NTI,Q,U,NU,ER,IW(1),JAC,D,KP,PHI,
     1            W(1),W(I4),W(I5),W(I6),W(I7),IW(J1),IW(J2),IERROR)
      IF ((IERROR.NE.0).AND.(IERROR.NE.240)) GOTO 6000
      IF (IERROR.NE.0) THEN
        WRITE(*,110) IERROR
        CALL ERRHAN(IERROR,ER,NEG)
      ENDIF
C     XI IN W(I1)
      NM = N * NRTI
      CALL DINPRO(W(I1),1,W(I1),1,NM,RKSI)
      RKSI = DSQRT(RKSI)
      IF (DIAGNO) THEN
        WRITE(*,230) RKSI
        WRITE(*,240) ER(4),ER(5)
      ENDIF
      IF (RKSI.LT.TOL1) GOTO 2500
      ALAM = 1.D0
      ALM1 = 1.D0
C     COMPUTE NEXT ITERATION
 1400 DO 1500 I = 1 , NRTI
        L1 = (I-1) * N + I3 - 1
        L2 = (I-1) * N + I1 - 1
        DO 1500 J = 1 , N
          W(L1+J) = X(J,I) + ALAM * W(L2+J)
 1500 CONTINUE
C     SI+ STORED IN W(I3)
      IF (DIAGNO) THEN
        WRITE(*,250) ALAM
      ENDIF
      NG = 0
      CALL DCSAOJ(FDIF,N,ER1,TI,NTI,NRTI,W(I3),JAC,NG,IW(1),W(I2),Q,U,
     1            NU,W(1),W(I4),W(I5),IW(J1),WGR,LWG,IERROR)
      NEG = IW(J1)
      IF ((IERROR.NE.0).AND.(IERROR.NE.213)) GOTO 6000
      IF (IERROR.NE.0) THEN
        WRITE(*,110) IERROR
        CALL ERRHAN(IERROR,ER,NEG)
      ENDIF
      IF (JAC.NE.0) CALL DCHINC(W(I3),W(I2),N,NTI,U,NU,NRTI,ALINC(2),
     1                          IW(1),IW(J3))
      CALL DJINGX(G,N,NRTI,W(I3),W(I2),NTI,Q,U,NU,ER,IW(1),JAC,D,KP,PHI,
     1            W(1),W(I4),W(I5),W(I6),W(I7),IW(J1),IW(J2),IERROR)
C     XI+ STORED IN W(I2)
      IF ((IERROR.EQ.260).AND.(JAC.NE.0)) THEN
        ALAM = ALAM / 2.D0
        IF (DIAGNO) THEN
          WRITE(*,260) ALAM
        ENDIF
        IF (ALAM.LT.ALAMIN) GOTO 5100
        CALL DCSAOJ(FDIF,N,ER1,TI,NTI,NRTI,X,JAC,NG,IW(1),W(I1),Q,U,NU,
     1              W(1),W(I4),W(I5),IW(J1),WGR,LWG,IERROR)
        NEG = IW(J1)
        CALL DJINGX(G,N,NRTI,X,W(I1),NTI,Q,U,NU,ER,IW(1),JAC,D,KP,PHI,
     1              W(1),W(I4),W(I5),W(I6),W(I7),IW(J1),IW(J2),IERROR)
        JAC = 0
        IW(J3) = 0
        GOTO 1400
      ENDIF
      IF ((IERROR.NE.0).AND.(IERROR.NE.240)) GOTO 6000
      IF (IERROR.NE.0) THEN
        WRITE(*,110) IERROR
        CALL ERRHAN(IERROR,ER,NEG)
      ENDIF
      CALL DINPRO(W(I2),1,W(I2),1,NM,RKSI0)
      RKSI0 = DSQRT(RKSI0)
      IF (DIAGNO) THEN
        WRITE(*,230) RKSI0
        WRITE(*,240) ER(4),ER(5)
      ENDIF
      IF (RKSI0.GT.RKSI) THEN
        ALAM = ALAM / 2.D0
        IHL = 1
        IF (ALAM.LT.ALAMIN) GOTO 3000
        IF (JAC.EQ.1) THEN
        CALL DCSAOJ(FDIF,N,ER1,TI,NTI,NRTI,X,JAC,NG,IW(1),W(I1),Q,U,
     1              NU,W(1),W(I4),W(I5),IW(J1),WGR,LWG,IERROR)
        NEG = IW(J1)
        CALL DJINGX(G,N,NRTI,X,W(I1),NTI,Q,U,NU,ER,IW(1),JAC,D,KP,PHI,
     1              W(1),W(I4),W(I5),W(I6),W(I7),IW(J1),IW(J2),IERROR)
        ENDIF
        JAC = 0
        IW(J3) = 0
        GOTO 1400
      ELSE
        JAC1 = JAC
        IF ((RKSI0.LT.0.6D0*RKSI).AND.(ALAM.EQ.1.D0)) JAC = 1
        IF (JAC1.NE.JAC) THEN
          DO 1650 I = J3 , J3+NRTI-1
             IW(I) = 0
 1650     CONTINUE
          IW(J3) = -1
        ENDIF
        RKSI = RKSI0
      ENDIF
C     LAMBDA IS ACCEPTED
C     COMPUTING NEXT PREDICTED LAMBDA
      IF ((ALAM.EQ.ALM1).AND.(IHL.EQ.0)) THEN
        ALAM = DMIN1(2.D0 * ALAM,1.D0)
      ELSE
        IHL = 0
        ALM1 = ALAM
      ENDIF
      DO 1700 I = 1 , NRTI
        L1 = (I-1) * N + I3 - 1
        DO 1700 J = 1 , N
          X(J,I) = W(L1+J)
 1700 CONTINUE
C     X := SI+ STORED IN X
      DO 1800 I = 1 , NM
        W(I1-1+I) = W(I2-1+I)
 1800 CONTINUE
C     XI := XI+  STORED IN W(I1)
      ITER = ITER + 1
      ITTOT = ITTOT + 1
      IF (DIAGNO) THEN
        IF (ALM1.EQ.1.D0) THEN
          WRITE(*,280) ITTOT
        ELSE
          WRITE(*,290) ITTOT,ALM1
        ENDIF
      ENDIF
      IF (RKSI.LT.TOL1) GOTO 2500
      IF (ITER.GT.ITLIM) GOTO 5200
      IF ((JAC.EQ.0).OR.(IW(J3).EQ.0)) GOTO 1400
 2400 CALL DNEWPO(FDIF,N,TOL1,ER(3),TI,NTI,NRTI,ALINC,X,W(I1),IW(1),Q,
     1            U,NU,IW(J3),WGR,LWG,W(1),W(I4),W(I5),IERROR)
      NEG = IW(J3)
      IF (IERROR.NE.0) GOTO 6000
      IF (DIAGNO.AND.NRTI.NE.NOSH) THEN
        WRITE(*,300) NRTI-NOSH
        NOSH = NRTI
      ENDIF
      GOTO 1100
 2500 IF (TOL1.LE.TOL) THEN
       WRITE(*,*) 'RETURN FROM MUSN : NO ITERATIONS =',ITTOT
       RETURN
      ENDIF
      TOL1 = DMAX1(TOL1*TOL1,TOL)
      IF (DABS(TOL-TOL1).LT.100.D0*ER(3)) TOL1 = TOL
      ER1(1) = 1.1D-12 + 3.D0 * ER(3)
      ER1(2) = TOL1
      ITER = 0
      IF (DIAGNO) THEN
        WRITE(*,310) TOL1
      ENDIF
      JAC = 0
      NG = 1
      CALL DCSAOJ(FDIF,N,ER1,TI,NTI,NRTI,X,JAC,NG,IW(1),W(I1),Q,U,NU,
     1            W(1),W(I4),W(I5),IW(J1),WGR,LWG,IERROR)
      NEG = IW(J1)
      IF ((IERROR.NE.0).AND.(IERROR.NE.213)) GOTO 6000
      IW(J3) = -1
      DO 2600 I = J3+1 , J3+NRTI-1
        IW(I) = 0
 2600 CONTINUE
      JAC = 1
      GOTO 2400
 3000 TOL1 = 1.D-1 * TOL1
      IF (DABS(TOL-TOL1).LT.100.D0*ER(3)) TOL1 = TOL
      ER1(1) = 1.1D-12 + 3.D0 * ER(3)
      ER1(2) = TOL1
      IF (TOL1.LT.TOL) GOTO 5100
      IF (DIAGNO) THEN
        WRITE(*,320) TOL1
      ENDIF
      NRTI = 1
      GOTO 1000
 5000 IERROR = 105
      GOTO 6000
 5010 IERROR = 101
      GOTO 6000
 5020 IERROR = 106
      GOTO 6000
 5100 IERROR = 230
      GOTO 6000
 5200 IERROR = 231
 6000 WRITE(*,100) IERROR
      CALL ERRHAN(IERROR,ER,NEG)
      RETURN
  100 FORMAT(' TERMINAL ERROR IN MUSN : IERROR = ',I4,/)
  110 FORMAT(' WARNING ERROR IN MUSN : IERROR = ',I4,/)
  200 FORMAT(' ',4('*****'),' START TOLERANCE : ',8X,1P,D12.5)
  210 FORMAT(5X,'SHOOTING POINTS :')
  220 FORMAT(5X,5(D12.5,3X))
  230 FORMAT(5X,'NORMALIZED RESIDUE : ',22X,1P,D16.9)
  240 FORMAT(5X,'COND.NUMBER = ',1P,D12.5,2X,'AMPLI.FACTOR = ',D12.5)
  250 FORMAT(5X,'PREDICTED DAMPING FACTOR : ',16X,F10.7)
  260 FORMAT(' IERROR = 260 , NEWTON UPDATE CAUSES ILL-CONDITIONING',/,
     1 ' WITH RESPECT TO THE BC; WILL TRY A DAMPING FACTOR=',F10.7)
  280 FORMAT(' ',10('******'),/,' ',4('****'),2X,I4,
     1 ' ITERATION COMPLETED',2X,4('****'),/,' ',4('*****'),2X,
     2 'FULL NEWTON STEP',2X,4('*****'),/,' ',10('******'),/)
  290 FORMAT(' ',10('******'),/,' ',4('****'),2X,I4,
     1 ' ITERATION COMPLETED',2X,4('****'),/,' ***** ',
     2 'DAMPED NEWTON STEP ; DAMPING FACTOR = ',F10.7,' *****',/,' ',
     3 10('******'),/)
  300 FORMAT(' ',I3,' NEW SHOOTING POINTS INSERTED.',/)
  310 FORMAT(' ',4('*****'),' NEW TOLERANCE : ',10X,1P,D12.5)
  320 FORMAT(' ***** NEWTON FAILED WITH THIS TOLERANCE. ***** ',/,
     1 ' ***** WILL TRY NEW TOLERANCE : ',10X,'TOL = ',1P,D12.5)
C     END OF MUSN
      END
      SUBROUTINE DCSHPO(FDIF,X0T,N,A,B,ER,TI,NTI,NRTI,ALINC,X,S,Q,U,NU,
     1                  W,WTI,WF,WF0,KKK,IP,WGR,LWGR,IERROR)
C     ----------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ER(5),TI(NTI),ALINC(4),X(N,NTI),S(N,NTI),Q(N,N,NTI),
     1          U(NU,NTI),W(N,7),WTI(NTI),WF(N,N),WF0(N,N),WGR(LWGR)
      INTEGER KKK(N),IP(N)
      EXTERNAL FDIF,X0T
C
      IERROR = 0
      EPS = DSQRT(ER(3))
      TI(1) = A
      WTI(1) = A
      NOTI = 2
      IF (NRTI.EQ.0) TI(2) = B
      IF (NRTI.EQ.1) THEN
        IF (A.LE.B) THEN
          DO 1100 I = 2 , NTI
            IF (TI(I-1).GE.TI(I)) GOTO 5200
            IF (TI(I).GE.B) GOTO 1200
 1100     CONTINUE
 1200     NOTI=I
        ELSE
          DO 1300 I = 2 , NTI
            IF (TI(I-1).LE.TI(I)) GOTO 5200
            IF (TI(I).LE.B) GOTO 1400
 1300     CONTINUE
 1400     NOTI=I
        ENDIF
        IF (TI(NOTI).NE.B) GOTO 5300
      ENDIF
      IF (NRTI.GT.1) THEN
        NOTI = NRTI + 1
        ST = (B-A) / NRTI
        TI(NRTI+1) = B
        DO 1500 I = 2 , NRTI
          TI(I) = TI(I-1) + ST
 1500   CONTINUE
      ENDIF
      IF (LWGR.LT.NOTI) GOTO 5550
C     DETERMINATION OF THE SHOOTING POINTS. GIVEN OUTPUTPOINTS WILL
C     BE SHOOTING POINTS. IF THE INCREMENT BETWEEN TWO OUTPUTPOINTS
C     IS GREATER THAN ALINC(1), A NEW SHOOTING POINT IS INSERTED.
C     DURING THE DETERMINATION THE SHOOTING POINTS WILL BE STORED IN
C     THE ARRAY WTI.
      DO 1600 I = 1 , N
        KKK(I) = I*(I+1) / 2
      DO 1600 J = 1 , N
        IF (I.EQ.J) THEN
          Q(I,I,1) = 1.D0
        ELSE
          Q(I,J,1) = 0.D0
        ENDIF
        Q(I,J,NTI) = Q(I,J,1)
 1600 CONTINUE
      WTI(1) = TI(1)
      JTI = 2
      NHI = 5
      NRTI = 2
      MWGR = 0
      DO 2600 K = 2 , NTI
        KM1 = K - 1
        T1 = WTI(KM1)
        T2 = TI(JTI)
 1650   CALL X0T(T1,X(1,KM1))
        IF (T1.EQ.B) GOTO 2700
        X1 = 0.D0
        DO 1700 I = 1 , N
          ST = DABS(X(I,KM1))
          S(I,K) = X(I,KM1)
          IF (ST.GT.X1) X1 = ST
          DO 1700 J = 1 , N
            WF(J,I) = X(J,KM1) + EPS * Q(J,I,KM1)
 1700   CONTINUE
        IF (X1.EQ.0.D0) X1 = 1.D0
        IWGR = 0
        IER = 1
        IF (K.GT.2) IER=-1
 1900   CALL DRKFMS(FDIF,N,S(1,K),T1,T2,ER(1),ER(2),ER(3),WF,N,HI,NHI,
     1              W,IER)
        IF (IER.GT.3) GOTO 5400
C     CHECK WHETHER T2 SHOULD BE A SHOOTING POINT
        X2 = 0.D0
        X3 = 0.D0
        DO 2000 I = 1 , N
          ST = DABS(S(I,K))
          IF (ST.GT.X2) X2 = ST
          DO 2000 J = 1 , N
            ST = DABS(WF(J,I)-S(J,K))
            IF (ST.GT.X3) X3 = ST
 2000   CONTINUE
        IF ((X2/X1.LT.ALINC(1)).AND.(X3.LT.ALINC(3))) THEN
          IWGR = IWGR + 1
          MWGR = MWGR + 1
          IF (MWGR.GT.LWGR) GOTO 5600
          T1 = T2
          WGR(MWGR) = T2
          IER = -1
          IF (T2.EQ.TI(JTI)) GOTO 2040
          T2 = TI(JTI)
          GOTO 1900
        ENDIF
        IF ((X2/X1.LE.ALINC(2)).AND.(X3.LE.ALINC(4))) THEN
          MWGR = MWGR + 1
          IF (MWGR.GT.LWGR) GOTO 5600
          IWGR = IWGR + 1
          WGR(MWGR) = T2
          IF (T2.EQ.TI(JTI)) GOTO 2040
          WTI(K) = T2
          GOTO 2050
        ENDIF
C     INCREMENT GREATER THAN ALINC(2)
        IF (IWGR.EQ.0) THEN
          T1 = WTI(KM1)
          T2 = (T2 + T1) / 2.D0
          GOTO 1650
        ELSE
          WTI(K) = WGR(MWGR)
          GOTO 2050
        ENDIF
C     A SHOOTING POINT IS DETERMINED.
C     COMPUTATION OF JK, QK AND UK
 2040   JTI = JTI + 1
        WTI(K) = T2
 2050   DO 2100 I = 1 , N
        DO 2100 J = 1 , N
          WF(J,I) = (WF(J,I) - S(J,K)) / EPS
          IF (K.EQ.2) WF0(J,I) = WF(J,I)
 2100   CONTINUE
        CALL DQUDEC(WF,N,N,N,N,W(1,1),W(1,2),W(1,3))
        IF (K.EQ.2) THEN
          ISORT = 0
          DO 2400 ITEL = 1 , N
            CALL DSORTD(W(1,1),N,N,KP,IP,ISORT,IER)
            IF (IER.EQ.0) GOTO 2500
            DO 2200 I = 1 , N
              IK = IP(I)
            DO 2200 J = 1 , N
              WF(J,I) = WF0(J,IK)
              Q(J,I,1) = Q(J,IK,NTI)
 2200       CONTINUE
            DO 2300 I = 1 , N
            DO 2300 J = 1 , N
              WF0(I,J) = WF(I,J)
              Q(I,J,NTI) = Q(I,J,1)
 2300       CONTINUE
            CALL DQUDEC(WF,N,N,N,N,W(1,1),W(1,2),W(1,3))
 2400     CONTINUE
        ENDIF
 2500   CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),K,KKK,U,NU,NTI,Q,WF0,W(1,3))
        NRTI = K
 2600 CONTINUE
 2700 DO 2800 I = 1 , NRTI
        TI(I) = WTI(I)
 2800 CONTINUE
      IF (TI(NRTI).NE.B) GOTO 5500
      RETURN
 5200 IERROR = 120
      RETURN
 5300 IERROR = 121
      RETURN
 5400 IERROR = 210 + IER
      RETURN
 5500 IERROR = 122
      RETURN
 5550 IERROR = 123
      IP(1) = NOTI
      RETURN
 5600 IERROR = 219
      IP(1) = (B - A) / (T1 - A) * (LWGR+1) + 1
      RETURN
C     END OF DCSHPO
      END
      SUBROUTINE DCSAOJ(FDIF,N,ER,TI,NTI,NRTI,X,JA,NG,KKK,S,Q,U,NU,W,
     1                  WF,WF0,IP,WGR,LWGR,IERROR)
C     ----------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ER(5),TI(NTI),X(N,NTI),S(N,NTI),Q(N,N,NTI),U(NU,NTI),
     1          W(N,7),WF(N,N),WF0(N,N),WGR(LWGR)
      INTEGER KKK(N),IP(N)
      DIMENSION HS(5)
      EXTERNAL FDIF
C
C     INPUT ERROR 122
      IF (NTI-NRTI.LT.1) GOTO 5000
      IERROR = 0
      NHI = 0
      EPS = DSQRT(ER(3))
      IF (JA.EQ.1) THEN
        DO 1100 I = 1 , N
        DO 1100 J = 1 , N
          Q(I,J,NTI) = Q(I,J,1)
 1100   CONTINUE
      ENDIF
      MWGR = 1
      IF (NG.EQ.1) MWGR = 0
      DO 1900 K = 2 , NRTI
        KM1 = K - 1
        T1 = TI(KM1)
        T2 = TI(K)
        DO 1200 I = 1 , N
          S(I,K) = X(I,KM1)
 1200   CONTINUE
        IF (JA.EQ.1) THEN
          DO 1300 I = 1 , N
          DO 1300 J = 1 , N
            WF(J,I) = X(J,KM1) + EPS * Q(J,I,KM1)
 1300     CONTINUE
          NF = N
        ELSE
          NF = 0
        ENDIF
        IF (NG.EQ.1) THEN
         NHI = 5
         IER = -1
         IF (K.EQ.2) IER = 1
 1320    CALL DRKFMS(FDIF,N,S(1,K),T1,T2,ER(1),ER(2),ER(3),WF,NF,HI,NHI,
     1               W,IER)
         IF (IER.GT.3) GOTO 5100
         MWGR = MWGR + 1
         IF (MWGR.GT.LWGR) GOTO 5200
         WGR(MWGR) = T2
         IF (T2.NE.TI(K)) THEN
           T1 = T2
           T2 = TI(K)
           IER = -1
           GOTO 1320
         ENDIF
        ELSE
          T0 = T1
 1340     ST = (WGR(MWGR) - T0) / 5.D0
          DO 1360 I = 1 , 5
            HS(I) = ST
 1360     CONTINUE
          CALL DRKFGS(FDIF,N,S(1,K),T1,HS,5,W)
          IF (JA.EQ.1) THEN
            DO 1380 I = 1 , N
              T1 = T0
              CALL DRKFGS(FDIF,N,WF(1,I),T1,HS,5,W)
 1380       CONTINUE
          ENDIF
            T0 = WGR(MWGR)
            MWGR = MWGR + 1
            IF (T0.NE.T2) GOTO 1340
        ENDIF
C     IF JA=1 COMPUTATION OF JK, QK AND UK
        IF (JA.EQ.1) THEN
          DO 1400 I = 1 , N
          DO 1400 J = 1 , N
            WF(J,I) = (WF(J,I) - S(J,K)) / EPS
            IF (K.EQ.2) WF0(J,I) = WF(J,I)
 1400     CONTINUE
          CALL DQUDEC(WF,N,N,N,N,W(1,1),W(1,2),W(1,3))
          IF (K.EQ.2) THEN
            ISORT = 0
            DO 1700 ITEL = 1 , N
              CALL DSORTD(W(1,1),N,N,KP,IP,ISORT,IER)
              IF (IER.EQ.0) GOTO 1800
              DO 1500 I = 1 , N
                IK = IP(I)
              DO 1500 J = 1 , N
                WF(J,I) = WF0(J,IK)
                Q(J,I,1) = Q(J,IK,NTI)
 1500         CONTINUE
              DO 1600 I = 1 , N
              DO 1600 J = 1 , N
                WF0(I,J) = WF(I,J)
                 Q(I,J,NTI) = Q(I,J,1)
 1600         CONTINUE
              CALL DQUDEC(WF,N,N,N,N,W(1,1),W(1,2),W(1,3))
 1700       CONTINUE
          ENDIF
 1800     CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),K,KKK,U,NU,NTI,Q,WF0,W(1,3))
        ENDIF
 1900 CONTINUE
      RETURN
 5000 IERROR = 122
      RETURN
 5100 IERROR = 210 + IER
      RETURN
 5200 IERROR = 219
      IP(1) = (TI(NRTI) - TI(1)) / (T1 - TI(1)) * (LWGR+1) + 1
      RETURN
C     END OF DCSAOJ
      END
      SUBROUTINE DJINGX(G,N,NRTI,X,S,NTI,Q,U,NU,ER,KKK,JAC,D,KPART,
     1                  PHIREC,W,WF,WF0,WBMA,WBMB,IP1,IP2,IERROR)
C     ---------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N,NTI),S(N,NTI),Q(N,N,NTI),U(NU,NTI),ER(5),D(N,NTI),
     1          PHIREC(NU,NTI),W(N,3),WF(N,N),WF0(N,N),WBMA(N,N),
     2          WBMB(N,N)
      INTEGER KKK(N),IP1(N),IP2(N)
      EXTERNAL G
C
      IERROR = 0
C     COMPUTATION OF THE INHOMOGENEOUS TERMS OF THE RECURSION.
      DO 1400 I = 2 , NRTI
        DO 1100 J = 1 , N
           W(J,1) = S(J,I) - X(J,I)
 1100   CONTINUE
        DO 1300 J = 1 , N
          SOM = 0.D0
          DO 1200 K = 1 , N
            SOM = SOM + Q(K,J,I) * W(K,1)
 1200     CONTINUE
          D(J,I) = SOM
 1300   CONTINUE
 1400 CONTINUE
      IF (JAC.NE.0) THEN
      IP1(1) = 0
      CALL DGTUR(N,1,NRTI,U,NU,NTI,Q,D,ER,KKK,KPART,IP1,IP2,W(1,1),
     1           W(1,2),W(1,3),WF,WF0,IERROR)
      CALL DKPCH(N,U,NU,NTI,NRTI,KPART,KKK,ER,IERROR)
      CALL DFUNRC(N,KPART,NRTI,U,NU,NTI,KKK,PHIREC,W(1,1),W(1,2),IERROR)
      IF (IERROR.NE.0) RETURN
      ENDIF
      CALL DPSR(N,KPART,NRTI,U,NU,NTI,D,KKK,S,W(1,1),W(1,2),IERROR)
      IF (IERROR.NE.0) RETURN
C     COMPUTATION OF THE BOUNDARY CONDITIONS.
      IF (JAC.EQ.0) THEN
        CALL G(N,X(1,1),X(1,NRTI),W(1,1),WF,WF0)
      ELSE
        CALL G(N,X(1,1),X(1,NRTI),W(1,1),WBMA,WBMB)
      ENDIF
      DO 1500 I = 1 , N
        W(I,1) = -W(I,1)
 1500 CONTINUE
      CALL DSBVP(N,1,WBMA,WBMB,NRTI,Q,NTI,W(1,1),PHIREC,NU,S,S,
     1           ER(3),ER(4),KKK,IP1,WF,WF0,W(1,2),W(1,3),IERROR)
      RETURN
C     END OF DJINGX
      END
      SUBROUTINE DCHINC(X,S,N,NTI,U,NU,NRTI,ALINCR,KKK,IPG)
C     -----------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N,NTI),S(N,NTI),U(NU,NTI)
      INTEGER KKK(N),IPG(NTI)
C
      IPG(1) = 0
      DO 1500 K = 2 , NRTI
        KM1 = K - 1
        E1 = 0.D0
        E2 = 0.D0
        DO 1000 J = 1 , N
          E3 = DABS(X(J,KM1))
          IF (E3.GT.E1) E1 = E3
          E3 = DABS(S(J,K))
          IF (E3.GT.E2) E2 = E3
 1000   CONTINUE
        E3 = E2 / E1
        IF (E3.GT.ALINCR) THEN
          IPG(1) = 1
          IPG(K) = 1
        ELSE
          IPG(K) = 0
        ENDIF
 1500 CONTINUE
      DO 2500 K = 2 , NRTI
        E1 = 0.D0
        DO 2000 J = 1 , N
          E2 = U(KKK(J),K)
          IF (E2.GT.E1) E1 = E2
 2000   CONTINUE
        IF (E1.GT.ALINCR) THEN
          IPG(1) = 1
          IPG(K) = 1
        ENDIF
 2500 CONTINUE
      RETURN
C     END OF DCHINC
      END
      SUBROUTINE DNEWPO(FDIF,N,TOL,EPSMA,TI,NTI,NRTI,ALINC,X,S,KKK,Q,U,
     1                  NU,IPG,WG,LWG,W,WF,WF0,IERROR)
C     -----------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TI(NTI),ALINC(4),X(N,NTI),S(N,NTI),Q(N,N,NTI),
     1          U(NU,NTI),WG(LWG),W(N,7),WF(N,N),WF0(N,N)
      INTEGER KKK(N),IPG(NTI)
      DIMENSION HS(10)
      EXTERNAL FDIF
C
      IERROR = 0
      IF (IPG(1).EQ.0) RETURN
C     INPUT ERROR 122
      IF (NTI-NRTI.LT.1) GOTO 5000
      LHI = 0
      NF = 0
      JAC = 0
      IF (IPG(1).EQ.-1) JAC = 1
      B = TI(NRTI)
      EPS = DSQRT(EPSMA)
      MWGR = 1
      K = 2
C     COMPUTING TOTAL NUMBER OF CURRENT GRID POINTS AND EXPECTED NEEDED
C     NUMBER OF GRIDPOINTS
      NEG = 0
      IK = 0
      IJ = 1
      DO 1300 I = 2 , NRTI
        IM1 = I - 1
        DO 1100 J = IJ , LWG
          IF (WG(J).EQ.TI(IM1)) IK = J
          IF (WG(J).EQ.TI(I)) GOTO 1200
 1100   CONTINUE
 1200   IJ = J
        IF (IPG(I).NE.0) NEG = NEG + IJ - IK
 1300 CONTINUE
      NTOTGR = IJ
      NEG = NEG + NTOTGR
 1500 IF ((K.GT.NTI).OR.(NRTI.GE.NTI)) GOTO 5000
      NSHPO = 0
      KM1 = K - 1
      KPL1 = K + 1
      T1 = TI(KM1)
      IF (T1.GE.B) GOTO 4900
      T2 = TI(K)
C     IF IPG(K)=0 NO SHOOTING POINT HAS TO BE INSERTED IN THE INTERVAL
C     (TI(K-1),TI(K)).
C     IF JAC=0 NO UPDATE FOR THE JACOBIAN IS NEEDED.
      IF (IPG(K).EQ.0) THEN
        DO 1600 I = 1 , N
          S(I,K) = X(I,KM1)
          IF (JAC.NE.0) THEN
            DO 1550 J = 1 , N
              WF(J,I) = X(J,KM1) + EPS * Q(J,I,KM1)
 1550       CONTINUE
          ENDIF
 1600   CONTINUE
        T0 = T1
 1700   ST = (WG(MWGR) - T0) / 5.D0
        DO 1800 I = 1 , 5
          HS(I) = ST
 1800   CONTINUE
        CALL DRKFGS(FDIF,N,S(1,K),T1,HS,5,W)
        IF (JAC.NE.0) THEN
          DO 1900 I = 1 , N
            T1 = T0
            CALL DRKFGS(FDIF,N,WF(1,I),T1,HS,5,W)
 1900     CONTINUE
        ENDIF
        T0 = WG(MWGR)
        T1 = T0
        MWGR = MWGR + 1
        IF (T1.NE.T2) GOTO 1700
        IF (JAC.NE.0) THEN
          DO 1950 I = 1 , N
          DO 1950 J = 1 , N
            WF(J,I) = (WF(J,I) - S(J,K)) / EPS
 1950     CONTINUE
          CALL DQUDEC(WF,N,N,N,N,W(1,1),W(1,2),W(1,3))
          CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),K,KKK,U,NU,NTI,Q,WF0,W(1,3))
        ENDIF
        K = K + 1
        GOTO 1500
      ENDIF
C     A NEW SHOOTING POINT IS NEEDED AND FROM NOW ON JACOBIAN UPDATES
C     ARE NECESSARY.
      JAC = 1
      NSHPO = 0
      KWGR = MWGR
      KWG = MWGR
      X1 = 0.D0
      DO 2000 I = 1 , N
        X2 = DABS(X(I,KM1))
        IF (X2.GT.X1) X1 = X2
        S(I,K) = X(I,KM1)
        S(I,KPL1) = X(I,KM1)
        Q(I,1,NTI) = 1.D0
        DO 2000 J = 1 , N
          WF(J,I) = X(J,KM1) + EPS * Q(J,I,KM1)
          WF0(J,I) = WF(J,I)
 2000 CONTINUE
      IWGR = 0
      IWG = 0
 2100 CALL DRKFGG(FDIF,N,T1,WG(MWGR),S(1,K),WF,HS,5,W)
C     CHECK INCREMENT
      X2 = 0.D0
      X3 = 0.D0
      DO 2300 I = 1 , N
        ST = DABS(S(I,K))
        IF (ST.GT.X2) X2 = ST
        X4 = 0.D0
        DO 2200 J = 1 , N
          ST = DABS((WF(J,I)-S(I,K))/Q(I,1,NTI))
          IF (ST.GT.X3) X3 = ST
          IF (ST.GT.X4) X4 = ST
 2200   CONTINUE
        Q(I,2,NTI) = X4
 2300 CONTINUE
      IF ((X2/X1.LT.ALINC(1)).AND.(X3.LT.ALINC(3))) THEN
        IWGR = IWGR + 1
        T1 = WG(MWGR)
        MWGR = MWGR + 1
        IF (T1.EQ.T2) GOTO 3500
        GOTO 2100
      ENDIF
      NSHPO = 1
      IF ((IWGR.EQ.0).AND.(WG(MWGR).EQ.T2)) GOTO 2800
      IF (KWGR.EQ.1) THEN
        T1 = TI(KM1)
      ELSE
        T1 = WG(KWGR-1)
      ENDIF
      DO 2400 I = KWGR , KWGR+IWGR
        CALL DRKFGG(FDIF,N,T1,WG(I),S(1,KPL1),WF0,HS,10,W)
        T1 = WG(I)
 2400 CONTINUE
      X4 = 0.D0
      DO 2600 I = 1 , N
        ST = DABS(S(I,KPL1))
        IF (ST.GT.X4) X4 = ST
        X5 = 0.D0
        DO 2500 J = 1 , N
          ST = DABS((WF0(J,I)-S(J,KPL1))/Q(I,1,NTI))
          IF (ST.GT.X5) X5 = ST
 2500   CONTINUE
        Q(I,2,NTI) = Q(I,2,NTI) / X5
 2600 CONTINUE
      X5 = X2 / X4
      DO 2700 I = 1 , N
        IF (Q(I,2,NTI).GT.X5) X5 = Q(I,2,NTI)
 2700 CONTINUE
      IF ((X5.GE.0.9D0).AND.(X5.LE.1.12D0)) GOTO 3200
C     NEW GRID NEEDED.
 2800 JWGR = IWGR + 1
      IF (NTOTGR+JWGR.GT.LWG) GOTO 5200
      DO 2900 I = NTOTGR , MWGR , -1
        WG(I+JWGR) = WG(I)
 2900 CONTINUE
      J = KWGR - 1
      DO 3000 I = JWGR , 1 , -1
        WG(J+2*I) = WG(J+I)
 3000 CONTINUE
      J = KWGR - 2
      IF (KWGR.EQ.1) THEN
        WG(KWGR) = (WG(KWGR+1) + TI(1)) / 2.D0
        IK = 2
      ELSE
        IK = 1
      ENDIF
      DO 3100 I = IK , JWGR
        I1 = J + 2*I
        WG(I1) = (WG(I1+1) + WG(I1-1)) / 2.D0
 3100 CONTINUE
      MWGR = MWGR + JWGR
 3200 IF (WG(MWGR).NE.T2) THEN
        MWGR = MWGR + 1
        KWGR = MWGR
        IWGR = 0
        X1 = 0.D0
        DO 3400 I = 1 , N
          ST = DABS(S(I,K))
          IF (ST.GT.X1) X1 = ST
          S(I,KPL1) = S(I,K)
          X4 = 0.D0
          DO 3300 J = 1 , N
            ST = DABS(WF(J,I))
            IF (ST.GT.X4) X4 = ST
            WF0(J,I) = WF(J,I)
 3300     CONTINUE
          Q(I,1,NTI) = X4
 3400   CONTINUE
        GOTO 2100
      ELSE
        MWGR = MWGR + 1
      ENDIF
 3500 IF (NSHPO.EQ.0) GOTO 4100
C     DETERMINATION OF THE NEW SHOOTING POINT.
      T1 = TI(KM1)
      T2 = TI(K)
      TP = (T2 + T1) / 2.D0
      E1 = DABS(T2-T1)
      DO 3600 I = KWG , LWG
        ST = DABS(TP-WG(I))
        IF (ST.LT.E1) THEN
          E1 = ST
          JP = I
        ELSE
          GOTO 3700
        ENDIF
 3600 CONTINUE
 3700 TP = WG(JP)
C     INSERT NEW SHOOTING POINT
      NRTI = NRTI + 1
      DO 3800 I = NRTI , KPL1 , -1
        IM1 = I-1
        TI(I) = TI(IM1)
        IPG(I) = IPG(IM1)
        DO 3800 J = 1 , N
          X(J,I) = X(J,IM1)
 3800 CONTINUE
      TI(K) = TP
C     COMPUTATION OF S(.,K) AND WF
      DO 3900 I = 1 , N
        S(I,K) = X(I,KM1)
        DO 3900 J = 1 , N
          WF(J,I) = X(J,KM1) + EPS * Q(J,I,KM1)
 3900 CONTINUE
      DO 4000 I = KWG , JP
        CALL DRKFGG(FDIF,N,T1,WG(I),S(1,K),WF,HS,5,W)
        T1 = WG(I)
 4000 CONTINUE
 4100 DO 4150 I = 1 , N
      DO 4150 J = 1 , N
        WF(J,I) = (WF(J,I) - S(J,K)) / EPS
 4150 CONTINUE
      CALL DQUDEC(WF,N,N,N,N,W(1,1),W(1,2),W(1,3))
      CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),K,KKK,U,NU,NTI,Q,WF0,W(1,3))
      IF (NSHPO.EQ.0) THEN
        K = K + 1
        GOTO 1500
      ENDIF
C     COMPUTATION NEW STARTING VALUE FOR X(.,K).
      X2 = U(1,K)
      DO 4200 I = 2 , N
        IK = KKK(I)
        IF (U(IK,K).LT.X2) X2 = U(IK,K)
 4200 CONTINUE
      ST = DLOG(X2) * (TI(K) - TI(KPL1)) / (TI(K) - TI(KM1))
      ST = DEXP(ST)
      AE = DMIN1(TOL,TOL/ST)
      RE=1.1D-12 + 3.D0 * EPSMA
      T1 = TI(KPL1)
      T2 = TI(K)
      IER = 1
      DO 4300 I = 1 , N
        X(I,K) = X(I,KPL1)
 4300 CONTINUE
      CALL DRKFMS(FDIF,N,X(1,K),T1,T2,RE,AE,EPSMA,WF,NF,HI,LHI,W,IER)
      IF (IER.GT.3) GOTO 5100
      DO 4400 I = 1 , N
        S(I,KPL1) = X(I,K)
        DO 4400 J = 1 , N
          WF(J,I) = X(J,K) + EPS * Q(J,I,K)
 4400 CONTINUE
      T1 = TP
      DO 4500 I = JP+1 , MWGR-1
        CALL DRKFGG(FDIF,N,T1,WG(I),S(1,KPL1),WF,HS,5,W)
        T1 = WG(I)
 4500 CONTINUE
      DO 4600 I = 1 , N
      DO 4600 J = 1 , N
        WF(J,I) = (WF(J,I) - S(J,KPL1)) / EPS
 4600 CONTINUE
      CALL DQUDEC(WF,N,N,N,N,W(1,1),W(1,2),W(1,3))
      CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),KPL1,KKK,U,NU,NTI,Q,WF0,W(1,3))
      K = K + 2
      GOTO 1500
 4900 IF (TI(NRTI).NE.B) GOTO 5000
      RETURN
 5000 IERROR = 122
      RETURN
 5100 IERROR = IER + 210
      RETURN
 5200 IERROR = 219
      IPG(1) = NEG
      RETURN
C     END OF DNEWPO
      END
