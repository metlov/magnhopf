      SUBROUTINE DDUR(FLIN,FDIF,N,IHOM,A,B,NRTI,AMP,TI,NTI,ER,Q,U,KU,WI,
     1                D,KPART,WY,W,WF,WF0,KKK,IP1,IP2,IERROR)
C     --------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TI(NTI),ER(5),Q(N,N,NTI),U(KU,NTI),WI(N,NTI),
     1 D(N,NTI),WY(N),W(N,7),WF(N,N),WF0(N,N)
C
      INTEGER KKK(N),IP1(N),IP2(N)
C
      LOGICAL HOM0,HOM1,SORT,INC
C
      EXTERNAL FLIN,FDIF
C
C.....SETTING STARTING PARAMETERS
      IERROR = 0
      JFLAG = 0
      ISORT = IP1(1)
      NKP = N
      IF (ISORT.EQ.1) NKP = IP1(2)
      JSORT = ISORT + 2
      HOM0 = IHOM.EQ.0
      HOM1 = IHOM.NE.0
      SORT = .TRUE.
      ITH = IHOM
      TI(1) = A
      NTI2 = NTI - 3
      NTI1 = NTI - 1
      IF (NRTI.GT.NTI2) GOTO 7200
      TOL = DMAX1(ER(1),ER(2))
      AMP3 = TOL / ER(3)
C.....DETERMINE OPTION FOR THE OUTPUT POINTS
      IF (NRTI.EQ.0) THEN
C.....NRTI = 0
        IF (AMP.LE.1.D0) AMP = TOL / (10.D0*ER(3))
        NOTI = 2
        TI(2) = B
      ENDIF
      IF (NRTI.EQ.1) THEN
C.....NRTI = 1
        IF (A.LT.B) THEN
          DO 1100 I = 2 , NTI
            IF (TI(I-1).GE.TI(I)) GOTO 7000
            IF (TI(I).GE.B) GOTO 1300
 1100     CONTINUE
        ELSE
          DO 1200 I = 2 , NTI
            IF (TI(I-1).LE.TI(I)) GOTO 7000
            IF (TI(I).LE.B) GOTO  1300
 1200     CONTINUE
        ENDIF
 1300   IF (TI(I).NE.B) GOTO 7100
        NOTI = I
      ENDIF
      IF (NRTI.GT.1) THEN
C.....NRTI > 1
        SOM = (B-A) / NRTI
        NOTI = NRTI + 1
        TI(NOTI) = B
        DO 1400 I = 2 , NRTI
          TI(I) = TI(I-1) + SOM
 1400   CONTINUE
      ENDIF
      IF (AMP.GT.1.D0) THEN
C.....AMP > 1
        INC = .TRUE.
        AMP1 = AMP
        AMP2 = 2.D0 * AMP
      ELSE
        INC = .FALSE.
      ENDIF
      IF (NOTI.GT.NTI2) GOTO 7200
      DO 1500 I = 1 , N
        KKK(I) = I * (I+1) / 2
        IP2(I) = I
      DO 1500 J = 1 , N
        Q(J,I,NTI) = Q(J,I,1)
 1500 CONTINUE
      IU = 0
C.....Q1 = I , QNTI = Q1 FOR MAKING IT POSSIBLE TO PERMUTATE Q1
C.....LATERON
C.....COMPUTATION OF Q2 AND U2
      T1 = TI(1)
      T2 = TI(2)
      NRHS = 1
      IF (HOM0) NRHS = 0
      CALL DCPHIS(FLIN,FDIF,N,T1,T2,ER,2,IHOM,NRHS,WY,WI,Q,NTI,W,
     1            WF,WF0,IER)
      IF (IER.GT.3) GOTO 7300
      IF (IER.EQ.3) IERROR = 213
      CALL DAMTES(AMP3,W(1,1),N,IFLAG)
      IF (IFLAG.NE.0) SORT = .FALSE.
C.....ORDENING OF U2
      IF (SORT) THEN
        CALL DSORTD(W(1,1),N,NKP,KPART,IP1,ISORT,IFLAG)
        IF (IFLAG.NE.0) THEN
          JFLAG = 1
          DO 1600 I = 1 , N
            IK = IP1(I)
            IP1(I) = IP2(IK)
          DO 1600 J = 1 , N
            WF(J,I) = WF0(J,IK)
 1600     CONTINUE
          DO 1700 I = 1 , N
            IP2(I) = IP1(I)
 1700     CONTINUE
          CALL DQUDEC(WF,N,N,N,N,W(1,1),W(1,2),WY)
        ENDIF
      ENDIF
      IF (HOM0) CALL DCNRHS(W(1,1),N,NRHS)
      CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),2,KKK,U,KU,NTI,Q,
     1           WF0,W(1,3))
      IF (HOM1) CALL DCDI(D,N,NTI,Q,2,W(1,3))
      IFLAG = 0
      IF (INC) THEN
        CALL DAMTES(AMP2,W(1,1),N,IFLAG)
        IF (IFLAG.EQ.1) IERROR = 200
      ENDIF
      T1 = T2
C.....COMPUTATION OF V2
 1800 IF (T1.EQ.TI(2)) GOTO 3500
      IFLAG = 0
      IF (INC) CALL DAMTES(AMP1,W(1,1),N,IFLAG)
      T2 = TI(2)
      CALL DCPHIS(FLIN,FDIF,N,T1,T2,ER,3,IHOM,NRHS,WY,WI,Q,NTI,W,WF,
     1            WF0,IER)
      IF (IER.GT.3) GOTO 7300
      IF (HOM0) CALL DCNRHS(W(1,1),N,NRHS)
      CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),3,KKK,U,KU,NTI,Q,WF0,W(1,3))
 1850 IF (IFLAG.EQ.1) THEN
        NOTI = NOTI + 1
        IF (NOTI.GT.NTI2) GOTO 7200
        DO 1900 I = NOTI , 3 , -1
          TI(I) = TI(I-1)
 1900   CONTINUE
        TI(2) = T1
        IU = 1
        T1 = T2
        GOTO 3500
      ENDIF
      CALL DUPUP(N,U(1,3),U(1,2),U(1,4))
      DO 2000 I = 1 , N
        IK = KKK(I)
        W(I,4) = U(IK,4)
 2000 CONTINUE
      IF (INC) THEN
        CALL DAMTES(AMP2,W(1,4),N,IFLAG)
        IF (IFLAG.EQ.1) GOTO 1850
      ENDIF
      IF (HOM1) THEN
C.....COMPUTATION OF U3.D2
      DO 2200 I = 1 , N
        SOM = 0.0D+0
        DO 2100 J = I , N
          IK = KKK(J) - J
          SOM = SOM + U(IK+I,3) * D(J,2)
 2100   CONTINUE
        D(I,2) = SOM
 2200 CONTINUE
C.....COMPUTATION OF D'I+1, STORED IN D3
      CALL DCDI(D,N,NTI,Q,3,W(1,3))
      DO 2300 I = 1 , N
        D(I,2) = D(I,2) + D(I,3)
 2300 CONTINUE
      ENDIF
C.....CHECK ORDERING OF U4, IF NOT ORDERED U4.P = Q4.U4' IS
C.....COMPUTED, STORED IN Q4 AND U4
      IFLAG = 0
      IF (SORT) THEN
        CALL DAMTES(AMP3,W(1,4),N,IK)
        IF (IK.NE.0) THEN
          SORT = .FALSE.
          GOTO  3100
        ENDIF
        DO 2400 I = 1 , N
          IK = KKK(I) - I
        DO 2400 J = 1 , N
          IF (J.LE.I) THEN
            WF0(J,I) = U(IK+J,4)
          ELSE
            WF0(J,I) = 0.D0
          ENDIF
 2400   CONTINUE
        CALL DSORTD(W(1,4),N,NKP,KPART,IP1,JSORT,IFLAG)
        IF (IFLAG.NE.0) THEN
          JFLAG = 1
          DO 2500 I = 1 , N
            IK = IP1(I)
            IP1(I) = IP2(IK)
          DO 2500 J = 1 , N
            WF(J,I) = WF0(J,IK)
 2500     CONTINUE
          DO 2600 I = 1 , N
            IP2(I) = IP1(I)
 2600     CONTINUE
          CALL DQUDEC(WF,N,N,N,N,W(1,1),W(1,2),WY)
C.....U4 WAS NOT ORDERED, A NEW U4' IS COMPUTED, STORED IN U4
          CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),4,KKK,U,KU,NTI,Q,WF0,W(1,3))
C.....COMPUTATION OF NEW Q3' = Q3.Q4, STORED IN Q2
          DO 2700 J = 1 , N
            CALL DMATVC(Q(1,1,3),N,N,N,N,Q(1,J,4),Q(1,J,2))
 2700     CONTINUE
          IF (HOM0) CALL DCNRHS(W(1,1),N,NRHS)
          IF (HOM1) THEN
C.....COMPUTATION OF NEW D2' = INV(Q4).D2, STORED IN D2
            CALL DTAMVC(Q(1,1,4),N,N,N,N,D(1,2),W(1,3))
            DO 3000 I = 1 , N
              D(I,2) = W(I,3)
 3000       CONTINUE
          ENDIF
          GOTO 3250
        ENDIF
      ENDIF
C.....U4 WAS ALREADY ORDERED SO Q2 = Q3
 3100 DO 3200 J = 1 , N
      DO 3200 I = 1 , N
        Q(J,I,2) = Q(J,I,3)
 3200 CONTINUE
C.....U2 := U4
 3250 DO 3300 I = 1 , KU
       U(I,2) = U(I,4)
 3300 CONTINUE
      T1 = T2
      GOTO 1800
 3500 DO 3600 I = 1 , N
        IK = IP2(I)
        DO 3600 J = 1 , N
          Q(J,I,1) = Q(J,IK,NTI)
 3600 CONTINUE
      IF (TI(2).EQ.B) GOTO 6000
      DO 5000 K = 3 , NTI2
        KMIN1 = K - 1
        KPL1 = K + 1
        KPL2 = K + 2
        IF (TI(KMIN1).EQ.B) GOTO 6000
        IF (IU.EQ.1) THEN
          CALL DAMTES(AMP2,W(1,1),N,IFLAG)
          IF (IFLAG.EQ.1) IERROR = 200
          IU = 0
          T2 = TI(K)
        ELSE
          T1 = TI(KMIN1)
          T2 = TI(K)
          CALL DCPHIS(FLIN,FDIF,N,T1,T2,ER,K,IHOM,NRHS,WY,WI,Q,NTI,W,
     1                WF,WF0,IER)
          IF (IER.GT.3) GOTO 7300
          IF (HOM0) CALL DCNRHS(W(1,1),N,NRHS)
          CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),K,KKK,U,KU,NTI,Q,WF0,W(1,3))
          T1 = T2
          T2 = TI(K)
        ENDIF
        IF (HOM1) CALL DCDI(D,N,NTI,Q,K,W(1,3))
 4000   IF (T1.EQ.TI(K)) GOTO 5000
        IFLAG = 0
        IF (INC) CALL DAMTES(AMP,W(1,1),N,IFLAG)
        CALL DCPHIS(FLIN,FDIF,N,T1,T2,ER,KPL1,IHOM,NRHS,WY,WI,Q,NTI,W,
     1              WF,WF0,IER)
        IF (IER.GT.3) GOTO 7300
        IF (HOM0) CALL DCNRHS(W(1,1),N,NRHS)
        CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),KPL1,KKK,U,KU,NTI,Q,
     1             WF0,W(1,3))
 4100   IF (IFLAG.EQ.1) THEN
          NOTI = NOTI + 1
          IF (NOTI.GT.NTI2) GOTO 7200
          DO 4200 I = NOTI , KPL1 , -1
            TI(I) = TI(I-1)
 4200     CONTINUE
          TI(K) = T1
          IU = 1
          T1 = T2
          GOTO 5000
        ENDIF
        CALL DUPUP(N,U(1,KPL1),U(1,K),U(1,KPL2))
        DO 4300 I = 1 , N
          IK = KKK(I)
          W(I,4) = U(IK,KPL2)
 4300   CONTINUE
        IFLAG = 0
        IF (INC) THEN
          CALL DAMTES(AMP2,W(1,4),N,IFLAG)
          IF (IFLAG.EQ.1) GOTO 4100
        ENDIF
        IF (HOM1) THEN
C.....COMPUTATION OF UK+1.DK
          DO 4500 I = 1 , N
            SOM = 0.0D+0
            DO 4400 J = I , N
              IK = KKK(J) - J
              SOM = SOM + U(IK+I,KPL1) * D(J,K)
 4400       CONTINUE
           D(I,K) = SOM
 4500     CONTINUE
          CALL DCDI(D,N,NTI,Q,KPL1,W(1,3))
          DO 4600 I = 1 , N
            D(I,K) = D(I,K) + D(I,KPL1)
 4600     CONTINUE
        ENDIF
        DO 4700 I = 1 , N
        DO 4700 J = 1 , N
          Q(I,J,K) =Q(I,J,KPL1)
 4700   CONTINUE
        DO 4800 I = 1 , KU
          U(I,K) = U(I,KPL2)
 4800   CONTINUE
        T1 = T2
        T2 = TI(K)
        GOTO 4000
 5000 CONTINUE
 6000 NRTI = NOTI
      IF (TI(NRTI).NE.B) GOTO 7200
      RETURN
 7000 IERROR = 120
      RETURN
 7100 IERROR = 121
      RETURN
 7200 IERROR = 122
      RETURN
 7300 IERROR = IER + 210
      RETURN
C     END OF DDUR
      END
      SUBROUTINE DSORTD(DIAG,N,NKP,KPART,IP,ISORT,IFLAG)
C     ---------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DIAG(N)
      INTEGER IP(N)
C
      IFLAG = 0
      DO 1000 I = 1 , N
        IP(I) = I
 1000 CONTINUE
      NN = N
      IF (ISORT.EQ.1) THEN
C     ISORT = 1 ; CHECK WHETHER PERMUTATION IS NECESSARY.
        J = 0
        DO 1100 I = 1 , NKP
          IF (DABS(DIAG(I)).GT.1.D0) J = J + 1
 1100   CONTINUE
        IF (J.EQ.NKP) THEN
          KPART = J
          RETURN
        ENDIF
        NN = NKP
      ENDIF
C     ENDIF ISORT = 1
      IF (ISORT.LT.2) THEN
C     MONOTONIC ORDENING
        K = 0
        N1 = NN - 1
        DO 1400 I = 1 , N1
          SORT = 0.D0
          DO 1300 J = I , NN
            A = DABS(DIAG(J))
            IF (A.GT.SORT) THEN
              SORT = A
              JM = J
            ENDIF
 1300     CONTINUE
          IF (SORT.GT.1.D0) K = K+1
          IF (JM.NE.I) THEN
            IFLAG = IFLAG + 1
            A = DIAG(JM)
            DIAG(JM) = DIAG(I)
            DIAG(I) = A
            J = IP(JM)
            IP(JM) = IP(I)
            IP(I) = J
          ENDIF
 1400   CONTINUE
        KPART = K
        IF (DABS(DIAG(NN)).GT.1.D0) KPART = KPART + 1
        RETURN
      ENDIF
C     ENDIF ISORT < 2
C     ISORT = 2 OR 3; NO MONOTONIC ORDENING
      IF (ISORT.EQ.3) NN = NKP
      N1 = NN - 1
      DO 2200 I = 1 , N1
        A = DABS(DIAG(I))
        IF (A.LT.1.D0) THEN
          DO 2000 J = I+1 , NN
            B = DABS(DIAG(J))
            IF ((B.GT.1.D0).AND.(B/A.GE.10.D0)) THEN
              IFLAG = IFLAG + 1
              SORT = DIAG(J)
              DIAG(J) = DIAG(I)
              DIAG(I) = SORT
              JM = IP(J)
              IP(J) = IP(I)
              IP(I) = JM
              GOTO 2100
            ENDIF
 2000     CONTINUE
 2100     CONTINUE
        ENDIF
 2200 CONTINUE
      KPART = 0
      DO 2300 I = 1 , NN
        IF (DABS(DIAG(I)).GT.1.D0) KPART = KPART + 1
 2300 CONTINUE
      RETURN
C     END OF DSORTD
      END
      SUBROUTINE DRKFGS(F,NEQN,Y,T,HI,NHI,W)
C     --------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(NEQN),HI(5),W(NEQN,7)
      EXTERNAL F
C
      CALL F(T,Y,W(1,1))
      DO 2000 I = 1 , NHI
        CALL DRKF1S(F,NEQN,Y,T,HI(I),W)
        DO 1000 J = 1 , NEQN
          Y(J) = W(J,7)
 1000   CONTINUE
        T = T + HI(I)
        IF (I.LT.NHI) CALL F(T,Y,W(1,1))
 2000 CONTINUE
      RETURN
C     END OF DRKFGS
      END
      SUBROUTINE DRKFSM(F,NEQN,Y,T,TOUT,RE,AE,EPS,HI,NHI,W,IFLAG)
C     -----------------------------------------------------------
C
C     FEHLBERG FOURTH-FIFTH ORDER RUNGE-KUTTA METHOD
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(NEQN),HI(5),W(NEQN,7)
      LOGICAL HFAILD,OUTPUT
      EXTERNAL F
C
C     REMIN IS THE MINIMUM ACCEPTABLE VALUE OF RELERR.  ATTEMPTS
C     TO OBTAIN HIGHER ACCURACY WITH THIS SUBROUTINE ARE USUALLY
C     VERY EXPENSIVE AND OFTEN UNSUCCESSFUL.
C
      DATA REMIN/1.D-12/
C
      NHI = 0
      KFLAG = 0
C
C
C     CHECK INPUT PARAMETERS
C
C
      IF (NEQN .LT. 1) GO TO 10
      IF ((RE .LT. 0.0D0)  .OR.  (AE .LT. 0.0D0)) GO TO 10
      MFLAG=IABS(IFLAG)
      IF (IFLAG.NE.1) GOTO 10
C
      U26 = 26.0D0*EPS
      GO TO 50
C
C     INVALID INPUT
   10 IFLAG=8
      RETURN
C
   50 RER=2.0D0*EPS+REMIN
      IF (RE .GE. RER) GO TO 55
C
C     RELATIVE ERROR TOLERANCE TOO SMALL
      RE=RER
      KFLAG=3
C
   55 DT=TOUT-T
C
      A=T
      CALL F(A,Y,W(1,1))
      IF (T .NE. TOUT) GO TO 65
      IFLAG=2
      IF (KFLAG.EQ.3) IFLAG = 3
      RETURN
C
C
   65 H=DABS(DT)
      TOLN=0.
      DO 70 K=1,NEQN
        S = RE * DABS(Y(K))
        TOL = S + DMAX1(RE,AE)
        IF (TOL .LE. 0.D0) GO TO 70
        TOLN=TOL
        YPK=DABS(W(K,1))
        IF (YPK*H**5 .GT. TOL) H=(TOL/YPK)**0.2D0
   70 CONTINUE
      IF (TOLN .LE. 0.0D0) H=0.0D0
      H=DMAX1(H,U26*DMAX1(DABS(T),DABS(DT)))
C
C
C     SET STEPSIZE FOR INTEGRATION IN THE DIRECTION FROM T TO TOUT
C
      H=DSIGN(H,DT)
C
C     INITIALIZE OUTPUT POINT INDICATOR
C
      OUTPUT= .FALSE.
C
C     TO AVOID PREMATURE UNDERFLOW IN THE ERROR TOLERANCE FUNCTION,
C     SCALE THE ERROR TOLERANCES
C
      SCALE=2.0D0/RE
      AAE=SCALE*AE
C
C
C     STEP BY STEP INTEGRATION
C
  100 HFAILD= .FALSE.
C
C     SET SMALLEST ALLOWABLE STEPSIZE
C
      HMIN=U26*DABS(T)
C
C     ADJUST STEPSIZE IF NECESSARY TO HIT THE OUTPUT POINT.
C     LOOK AHEAD TWO STEPS TO AVOID DRASTIC CHANGES IN THE STEPSIZE AND
C     THUS LESSEN THE IMPACT OF OUTPUT POINTS ON THE CODE.
C
      DT=TOUT-T
      IF (DABS(DT) .GE. 2.0D0*DABS(H)) GO TO 200
      IF (DABS(DT) .GT. DABS(H)) GO TO 150
C
C     THE NEXT SUCCESSFUL STEP WILL COMPLETE THE INTEGRATION TO THE
C     OUTPUT POINT
C
      OUTPUT= .TRUE.
      H=DT
      GO TO 200
C
  150 H=0.5D0*DT
C
C     ADVANCE AN APPROXIMATE SOLUTION OVER ONE STEP OF LENGTH H
C
  200 CALL DRKF1S(F,NEQN,Y,T,H,W)
C
C     COMPUTE AND TEST ALLOWABLE TOLERANCES VERSUS LOCAL ERROR ESTIMATES
C     AND REMOVE SCALING OF TOLERANCES. NOTE THAT RELATIVE ERROR IS
C     MEASURED WITH RESPECT TO THE AVERAGE OF THE MAGNITUDES OF THE
C     SOLUTION AT THE BEGINNING AND END OF THE STEP.
C
      EEOET=0.0D0
      DO 250 K=1,NEQN
        ET=DABS(Y(K))+DABS(W(K,7))+AAE
        IF (ET .GT. 0.0D0) GO TO 240
C
C       INAPPROPRIATE ERROR TOLERANCE
        IFLAG=5
        RETURN
C
  240   EE=DABS((-2090.0D0*W(K,1)+(21970.0D0*W(K,4)-15048.0D0*W(K,5)))+
     1                        (22528.0D0*W(K,3)-27360.0D0*W(K,6)))
  250   EEOET=DMAX1(EEOET,EE/ET)
C
      ESTTOL=DABS(H)*EEOET*SCALE/752400.0D0
C
      IF (ESTTOL .LE. 1.0D0) GO TO 260
C
C
C     UNSUCCESSFUL STEP
C                       REDUCE THE STEPSIZE , TRY AGAIN
C                       THE DECREASE IS LIMITED TO A FACTOR OF 1/10
C
      HFAILD= .TRUE.
      OUTPUT= .FALSE.
      S=0.1D0
      IF (ESTTOL .LT. 59049.0D0) S=0.9D0/ESTTOL**0.2D0
      H=S*H
      IF (DABS(H) .GT. HMIN) GO TO 200
C
C     REQUESTED ERROR UNATTAINABLE AT SMALLEST ALLOWABLE STEPSIZE
      IFLAG=6
      KFLAG=6
      RETURN
C
C
C     SUCCESSFUL STEP
C                        STORE SOLUTION AT T+H
C                        AND EVALUATE DERIVATIVES THERE
C
  260 T=T+H
      NHI = NHI + 1
      HI(NHI) = H
      DO 270 K=1,NEQN
  270   Y(K)=W(K,7)
      IF (NHI.EQ.5) GOTO 320
      A=T
      CALL F(A,Y,W(1,1))
C
C
C                       CHOOSE NEXT STEPSIZE
C                       THE INCREASE IS LIMITED TO A FACTOR OF 5
C                       IF STEP FAILURE HAS JUST OCCURRED, NEXT
C                          STEPSIZE IS NOT ALLOWED TO INCREASE
C
      S=5.0D0
      IF (ESTTOL .GT. 1.889568D-4) S=0.9D0/ESTTOL**0.2D0
      IF (HFAILD) S=DMIN1(S,1.0D0)
      H=DSIGN(DMAX1(S*DABS(H),HMIN),H)
C
C     END OF CORE INTEGRATOR
C
C
C     SHOULD WE TAKE ANOTHER STEP
C
      IF (OUTPUT) GO TO 300
      IF (IFLAG .GT. 0) GO TO 100
C
C
C     INTEGRATION SUCCESSFULLY COMPLETED
C
C     INTERVAL MODE
  300 T=TOUT
      IFLAG=2
      IF (KFLAG.EQ.3) IFLAG = 3
      RETURN
C
  320 TOUT = T
      IFLAG = 2
      IF (KFLAG.EQ.3) IFLAG = 3
      RETURN
C
C     END OF DRKFSM
      END
      SUBROUTINE DRKF1S(F,NEQN,Y,T,H,W)
C     ---------------------------------
C
C     FEHLBERG FOURTH-FIFTH ORDER RUNGE-KUTTA METHOD
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(NEQN),W(NEQN,7)
      EXTERNAL F
C
      CH=H/4.0D0
      DO 221 K=1,NEQN
  221   W(K,6)=Y(K)+CH*W(K,1)
      CALL F(T+CH,W(1,6),W(1,2))
C
      CH=3.0D0*H/32.0D0
      DO 222 K=1,NEQN
  222   W(K,6)=Y(K)+CH*(W(K,1)+3.0D0*W(K,2))
      CALL F(T+3.0D0*H/8.0D0,W(1,6),W(1,3))
C
      CH=H/2197.0D0
      DO 223 K=1,NEQN
  223   W(K,6) = Y(K) + CH * (1932.0D0 * W(K,1)
     1                + (7296.0D0 * W(K,3) - 7200.0D0 * W(K,2)))
      CALL F(T+12.0D0*H/13.0D0,W(1,6),W(1,4))
C
      CH=H/4104.0D0
      DO 224 K=1,NEQN
  224   W(K,6)=Y(K)+CH*((8341.0D0*W(K,1)-845.0D0*W(K,4))+
     1                    (29440.0D0*W(K,3)-32832.0D0*W(K,2)))
      CALL F(T+H,W(1,6),W(1,5))
C
      CH=H/20520.0D0
      DO 225 K=1,NEQN
  225   W(K,2)=Y(K)+CH*((-6080.0D0*W(K,1)+(9295.0D0*W(K,4)-
     1         5643.0D0*W(K,5)))+(41040.0D0*W(K,2)-28352.0D0*W(K,3)))
      CALL F(T+H/2.0D0,W(1,2),W(1,6))
C
C     COMPUTE APPROXIMATE SOLUTION AT T+H
C
      CH=H/7618050.0D0
      DO 230 K=1,NEQN
  230   W(K,7)=Y(K)+CH*((902880.0D0*W(K,1)+(3855735.0D0*W(K,4)-
     1        1371249.0D0*W(K,5)))+(3953664.0D0*W(K,3)+
     2        277020.0D0*W(K,6)))
C
      RETURN
C     END OF DRKF1S
      END
      SUBROUTINE DRKFMS(F,NEQN,Y,T,TOUT,RE,AE,EPS,PHI,NPHI,
     1                  HI,NHI,W,IFLAG)
C     ----------------------------------------------------------------
C
C     FEHLBERG FOURTH-FIFTH ORDER RUNGE-KUTTA METHOD
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(NEQN),PHI(NEQN,NEQN),W(NEQN,7)
      LOGICAL HFAILD,OUTPUT
      EXTERNAL F
C
C     REMIN IS THE MINIMUM ACCEPTABLE VALUE OF RELERR.  ATTEMPTS
C     TO OBTAIN HIGHER ACCURACY WITH THIS SUBROUTINE ARE USUALLY
C     VERY EXPENSIVE AND OFTEN UNSUCCESSFUL.
C
      DATA REMIN/1.D-12/
C
      LHI = 0
      IST = 0
      IUST = 0
      KFLAG = 0
C
C
C     CHECK INPUT PARAMETERS
C
C
      IF (NEQN .LT. 1) GO TO 10
      IF ((RE .LT. 0.0D0)  .OR.  (AE .LT. 0.0D0)) GO TO 10
      MFLAG=IABS(IFLAG)
      IF (MFLAG.NE.1) GOTO 10
C
      U26 = 26.0D0*EPS
      GO TO 50
C
C     INVALID INPUT
   10 IFLAG=8
      RETURN
C
   50 RER=2.0D0*EPS+REMIN
      IF (RE .GE. RER) GO TO 55
C
C     RELATIVE ERROR TOLERANCE TOO SMALL
      RE=RER
      KFLAG=3
C
   55 DT=TOUT-T
C
      A=T
      CALL F(A,Y,W(1,1))
      IF (T .NE. TOUT) GO TO 65
      IFLAG=2
      IF (KFLAG.EQ.3) IFLAG = 3
      RETURN
C
C
   65 IF (IFLAG.EQ.1) THEN
      H=DABS(DT)
      TOLN=0.
      DO 70 K=1,NEQN
        S = RE * DABS(Y(K))
        TOL = S + DMAX1(RE,AE)
        IF (TOL .LE. 0.D0) GO TO 70
        TOLN=TOL
        YPK=DABS(W(K,1))
        IF (YPK*H**5 .GT. TOL) H=(TOL/YPK)**0.2D0
   70 CONTINUE
      IF (TOLN .LE. 0.0D0) H=0.0D0
      H=DMAX1(H,U26*DMAX1(DABS(T),DABS(DT)))
      HI = H
      ELSE
      H = DMIN1(DABS(HI),DABS(DT))
      ENDIF
C
C     SET STEPSIZE FOR INTEGRATION IN THE DIRECTION FROM T TO TOUT
C
      H=DSIGN(H,DT)
C
C     INITIALIZE OUTPUT POINT INDICATOR
C
      OUTPUT= .FALSE.
C
C     TO AVOID PREMATURE UNDERFLOW IN THE ERROR TOLERANCE FUNCTION,
C     SCALE THE ERROR TOLERANCES
C
      SCALE=2.0D0/RE
      AAE=SCALE*AE
C
C
C     STEP BY STEP INTEGRATION
C
  100 HFAILD= .FALSE.
C
C     SET SMALLEST ALLOWABLE STEPSIZE
C
      HMIN=U26*DABS(T)
C
C     ADJUST STEPSIZE IF NECESSARY TO HIT THE OUTPUT POINT.
C     LOOK AHEAD TWO STEPS TO AVOID DRASTIC CHANGES IN THE STEPSIZE AND
C     THUS LESSEN THE IMPACT OF OUTPUT POINTS ON THE CODE.
C
      DT=TOUT-T
      IF (DABS(DT) .GE. 2.0D0*DABS(H)) GO TO 200
      IF (DABS(DT) .GT. DABS(H)) GO TO 150
C
C     THE NEXT SUCCESSFUL STEP WILL COMPLETE THE INTEGRATION TO THE
C     OUTPUT POINT
C
      OUTPUT= .TRUE.
      H=DT
      GO TO 200
C
  150 H=0.5D0*DT
C
C     ADVANCE AN APPROXIMATE SOLUTION OVER ONE STEP OF LENGTH H
C
  200 CALL DRKF1S(F,NEQN,Y,T,H,W)
      IST = IST + 1
C
C     COMPUTE AND TEST ALLOWABLE TOLERANCES VERSUS LOCAL ERROR ESTIMATES
C     AND REMOVE SCALING OF TOLERANCES. NOTE THAT RELATIVE ERROR IS
C     MEASURED WITH RESPECT TO THE AVERAGE OF THE MAGNITUDES OF THE
C     SOLUTION AT THE BEGINNING AND END OF THE STEP.
C
      EEOET=0.0D0
      DO 250 K=1,NEQN
        ET=DABS(Y(K))+DABS(W(K,7))+AAE
        IF (ET .GT. 0.0D0) GO TO 240
C
C       INAPPROPRIATE ERROR TOLERANCE
        IFLAG=5
        RETURN
C
  240   EE=DABS((-2090.0D0*W(K,1)+(21970.0D0*W(K,4)-15048.0D0*W(K,5)))+
     1                        (22528.0D0*W(K,3)-27360.0D0*W(K,6)))
  250   EEOET=DMAX1(EEOET,EE/ET)
C
      ESTTOL=DABS(H)*EEOET*SCALE/752400.0D0
C
      IF (ESTTOL .LE. 1.0D0) GO TO 260
C
C
C     UNSUCCESSFUL STEP
C                       REDUCE THE STEPSIZE , TRY AGAIN
C                       THE DECREASE IS LIMITED TO A FACTOR OF 1/10
C
      IUST = IUST + 1
      HFAILD= .TRUE.
      OUTPUT= .FALSE.
      S=0.1D0
      IF (ESTTOL .LT. 59049.0D0) S=0.9D0/ESTTOL**0.2D0
      H=S*H
      IF (DABS(H) .GT. HMIN) GO TO 200
C
C     REQUESTED ERROR UNATTAINABLE AT SMALLEST ALLOWABLE STEPSIZE
      IFLAG=6
      KFLAG=6
      RETURN
C
C
C     SUCCESSFUL STEP
C                        STORE SOLUTION AT T+H
C                        AND EVALUATE DERIVATIVES THERE
C
  260 LHI = LHI + 1
      DO 270 K = 1 , NEQN
        Y(K) = W(K,7)
  270 CONTINUE
      IF (NPHI.GT.0) THEN
      DO 290 K = 1 , NPHI
        CALL F(T,PHI(1,K),W(1,1))
        CALL DRKF1S(F,NEQN,PHI(1,K),T,H,W)
        DO 280 L = 1 , NEQN
          PHI(L,K) = W(L,7)
  280   CONTINUE
  290 CONTINUE
      ENDIF
      T = T + H
      IF (LHI.EQ.NHI) GOTO 320
      A=T
      CALL F(A,Y,W(1,1))
C
C
C                       CHOOSE NEXT STEPSIZE
C                       THE INCREASE IS LIMITED TO A FACTOR OF 5
C                       IF STEP FAILURE HAS JUST OCCURRED, NEXT
C                          STEPSIZE IS NOT ALLOWED TO INCREASE
C
      S=5.0D0
      IF (ESTTOL .GT. 1.889568D-4) S=0.9D0/ESTTOL**0.2D0
      IF (HFAILD) S=DMIN1(S,1.0D0)
      H=DSIGN(DMAX1(S*DABS(H),HMIN),H)
      HI = H
C
C     END OF CORE INTEGRATOR
C
C     SHOULD WE TAKE ANOTHER STEP
C
      IF (OUTPUT) GO TO 300
      GO TO 100
C
C
C     INTEGRATION SUCCESSFULLY COMPLETED
C
C     INTERVAL MODE
  300 T=TOUT
      IFLAG=2
      IF (KFLAG.EQ.3) IFLAG = 3
      RETURN
C
  320 TOUT = T
      IFLAG = 2
      IF (KFLAG.EQ.3) IFLAG = 3
      RETURN
C
C     END OF DRKFSM
      END
      SUBROUTINE DRKFGG(FDIF,N,T1,T2,X,PHI,HS,NHS,W)
C     ----------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),PHI(N,N),HS(NHS),W(N,7)
      EXTERNAL FDIF
C
      T0 = T1
      ST = (T2-T1) / NHS
      DO 1000 I = 1 , NHS
        HS(I) = ST
 1000 CONTINUE
      CALL DRKFGS(FDIF,N,X,T1,HS,NHS,W)
      DO 1100 I = 1 , N
        T1 = T0
        CALL DRKFGS(FDIF,N,PHI(1,I),T1,HS,NHS,W)
 1100 CONTINUE
      RETURN
C     END OF DRKFGG
      END
      SUBROUTINE ERRHAN(IERROR,ER,NEG)
C     --------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ER(5)
C
      IF (IERROR.EQ.100) WRITE(*,100)
      IF (IERROR.EQ.101) WRITE(*,101)
      IF (IERROR.EQ.103) WRITE(*,103)
      IF (IERROR.EQ.105) WRITE(*,105)
      IF (IERROR.EQ.106) WRITE(*,106)
      IF (IERROR.EQ.120) WRITE(*,120)
      IF (IERROR.EQ.121) WRITE(*,121)
      IF (IERROR.EQ.122) WRITE(*,122)
      IF (IERROR.EQ.123) WRITE(*,123) NEG
      IF (IERROR.EQ.200) WRITE(*,200)
      IF (IERROR.EQ.213) WRITE(*,213) ER(1)
      IF (IERROR.EQ.215) WRITE(*,215)
      IF (IERROR.EQ.218) WRITE(*,218)
      IF (IERROR.EQ.219) WRITE(*,219) NEG
      IF (IERROR.EQ.230) WRITE(*,230)
      IF (IERROR.EQ.231) WRITE(*,231)
      IF (IERROR.EQ.240) WRITE(*,240) ER(5)
      IF (IERROR.EQ.250) WRITE(*,250)
      IF (IERROR.EQ.260) WRITE(*,260)
      RETURN
  100 FORMAT(' INPUT ERROR: N<1 OR IHOM<0 OR NRTI<0 OR NTI<5 OR A=B OR',
     1 ' NU<N*(N+1)/2')
  101 FORMAT(' INPUT ERROR: ER(1), ER(2) OR ER(3) IS NEGATIVE')
  103 FORMAT(' INPUT ERROR: LW<8*N+2*N*N OR LIW<3*N')
  105 FORMAT(' INPUT ERROR: N<1 OR NRTI<0 OR NTI<3 OR A=B OR',
     1 ' NU<N*(N+1)/2')
  106 FORMAT(' INPUT ERROR: LW<7*N+3*N*NTI+4*N*N OR LIW<3*N+NTI')
  120 FORMAT(' INPUT ERROR: NRTI=1, BUT SUPPLIED OUTPUT POINTS IN TI ',
     1 /,' ARE NOT IN STRICT MONOTONE ORDER')
  121 FORMAT(' INPUT ERROR: NRTI=1, BUT A OR B ARE NOT INCLUDED IN THE',
     1 /,' GIVEN OUTPUTPOINTS IN TI')
  122 FORMAT(' INPUT ERROR: NTI TOO SMALL !')
  123 FORMAT(' INPUT ERROR: LWG IS LESS THAN THE NUMBER OF OUTPUT  ',
     1 'POINTS.',/,' LWG MUST BE GREATER THAN:',I10)
  200 FORMAT(' INCREMENT ON MINOR SHOOTING INTERVAL IS',
     1 ' GREATER THAN AMP')
  213 FORMAT(' ER(1) TOO SMALL, CHANGED INTO:',D16.9)
  215 FORMAT(' PURE RELATIVE ERROR TEST IMPOSSIBLE; MUST USE NON-ZERO',
     1 /,' ABSOLUTE TOLERANCE')
  218 FORMAT(' N, ER(1) OR ER(2) IS NEGATIVE')
  219 FORMAT(' ARRAY WG IS TOO SMALL. ESTIMATED VALUE FOR LWG:',I10)
  230 FORMAT(' NECESSARY DAMPING FACTOR IS TOO SMALL; YOU MAY TRY MORE',
     1 /,' OUTPUT POINTS OR START WITH A BETTER APPROXIMATION OF',/,
     2 ' THE SOLUTION.')
  231 FORMAT(' NUMBER OF ITERATIONS > ITLIM')
  240 FORMAT(' BAD DICHOTOMY; AMPLIFICATION FACTOR = ',D12.5)
  250 FORMAT(' ONE OF THE UPPER TRIANGULAR MATRIX OF THE MULTIPLE',
     1 /,' SHOOTING RECURSION IS SINGULAR')
  260 FORMAT(' PROBLEM TOO ILL-CONDITIONED WITH RESPECT TO THE BC')
C     END OF ERRHAN
      END
