      SUBROUTINE MBHINIT
*
*
*       Static massive black hole at the origin (KZ(44) = 1).
*       -----------------------------------------------------
*
*       Reads &INMBH after the IMBH input of START: MBH [Msun], EPSPC
*       [pc], RTMODE, RACCRPC [pc], FACC, BETAF (mbh.h). The black hole
*       is a fixed point mass at the coordinate origin acting on the
*       irregular force of every body (MBHFORCE), on KS perturbations as
*       a tidal difference and in the initial polynomials (MBHFPOLY); its
*       potential energy enters EPL. A single star whose orbit about it
*       reaches r_t/BETAF, or 4 G M/c^2, is removed at its pericentre
*       passage (MBHCHECK, MBHTDE). Requires KZ(14) = 0.
*
      INCLUDE 'common6.h'
      INCLUDE 'mbh.h'
      REAL*8  MBH,EPSPC,RACCRPC
      NAMELIST /INMBH/ MBH,EPSPC,RTMODE,RACCRPC,FACC,BETAF
*
      MBH = 0.0D0
      EPSPC = 0.0D0
      RTMODE = 2
      RACCRPC = 0.0D0
      FACC = 0.5D0
      BETAF = 1.0D0
      IF (rank.eq.0) THEN
          READ (NML=INMBH, IOSTAT=IIC, UNIT=5)
          IF (IIC.NE.0) THEN
              WRITE (6,'(A)') ' MBHINIT: invalid &INMBH namelist'
              STOP
          END IF
      END IF
#if MPIINIT
      CALL MPI_BCAST(MBH,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(EPSPC,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RTMODE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(RACCRPC,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(FACC,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(BETAF,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif
      IF (MBH.LE.0.0D0.OR.KZ(14).NE.0.OR.EPSPC.LT.0.0D0.OR.
     &    BETAF.LE.0.0D0.OR.FACC.LT.0.0D0.OR.FACC.GT.1.0D0.OR.
     &    (RTMODE.EQ.1.AND.KZ(19).LT.3).OR.
     &    (RTMODE.EQ.2.AND.RACCRPC.LE.0.0D0).OR.
     &    RTMODE.LT.1.OR.RTMODE.GT.2) THEN
          IF (rank.eq.0) WRITE (6,10)  MBH, EPSPC, RTMODE, RACCRPC,
     &                                 FACC, BETAF, KZ(14), KZ(19)
   10     FORMAT (' MBHINIT: invalid MBH',1P,2E12.4,I3,3E12.4,
     &            ' KZ(14), KZ(19)',2I4)
          STOP
      END IF
      CMBH = MBH/ZMBAR
      EPSBH = EPSPC/RBAR
      RACCR = RACCRPC/RBAR
      NTDE = 0
      NPLUNGE = 0
      DMTDE = 0.0D0
      EMBHC = 0.0D0
      IF (rank.eq.0) WRITE (6,20)  MBH, CMBH, EPSPC, RTMODE, RACCRPC,
     &                             FACC, BETAF
   20 FORMAT (/,' MBH STATIC: M[Msun] M[NB] EPS[pc] RTMODE RACCR[pc] ',
     &        'FACC BETAF ',1P,3E13.5,I3,3E13.5)
*
      RETURN
*
      END
*
*
      SUBROUTINE MBHFORCE(XI,VI,FF,FFD)
*
*
*       Force and first derivative of the static massive black hole.
*       ------------------------------------------------------------
*
*       Added to FF and FFD for the body at XI with velocity VI.
*
      INCLUDE 'common6.h'
      INCLUDE 'mbh.h'
      REAL*8  XI(3),VI(3),FF(3),FFD(3)
*
      R2 = XI(1)**2 + XI(2)**2 + XI(3)**2 + EPSBH**2
      RV = XI(1)*VI(1) + XI(2)*VI(2) + XI(3)*VI(3)
      A3 = CMBH/(R2*SQRT(R2))
      DO 10 K = 1,3
          FF(K) = FF(K) - A3*XI(K)
          FFD(K) = FFD(K) - A3*(VI(K) - 3.0D0*RV*XI(K)/R2)
   10 CONTINUE
*
      RETURN
*
      END
*
*
      SUBROUTINE MBHFPOLY(I1,I2)
*
*
*       Black-hole force in the initial irregular polynomials.
*       ------------------------------------------------------
*
*       Adds the force and derivative of the static black hole to FI and
*       D1 of the bodies I1..I2 (with ghosts left out), as XTRNLD does
*       for an external field.
*
      INCLUDE 'common6.h'
      INCLUDE 'mbh.h'
*
      DO 10 I = I1,I2
          IF (BODY(I).GT.0.0D0) THEN
              CALL MBHFORCE(X(1,I),XDOT(1,I),FI(1,I),D1(1,I))
          END IF
   10 CONTINUE
*
      RETURN
*
      END
*
*
      REAL*8 FUNCTION MBHPOTE()
*
*
*       Potential energy of the bodies in the black-hole field.
*       -------------------------------------------------------
*
*       Sum of -m M/(r^2 + eps^2)^(1/2) over single bodies and KS centres
*       of mass (the tidal part of a pair's energy is in its h).
*
      INCLUDE 'common6.h'
      INCLUDE 'mbh.h'
*
      EP = 0.0D0
      DO 10 I = IFIRST,NTOT
          IF (BODY(I).GT.0.0D0) THEN
              R2 = X(1,I)**2 + X(2,I)**2 + X(3,I)**2 + EPSBH**2
              EP = EP - BODY(I)*CMBH/SQRT(R2)
          END IF
   10 CONTINUE
      MBHPOTE = EP
*
      RETURN
*
      END
*
*
      SUBROUTINE MBHCHECK(IACC)
*
*
*       Pericentre passages inside the black hole's capture radius.
*       -----------------------------------------------------------
*
*       For every single star of the current block, the two-body orbit
*       about the black hole (corrected X0, X0DOT) gives the pericentre q
*       and the time to (or since) it. The star is removed (MBHTDE) when
*       q is below the capture radius, max(r_t/BETAF, 4 G M/c^2), and the
*       pericentre lies within its next step ahead or 2 steps behind, so
*       that a passage between two steps is not missed (a test of the
*       sampled position alone would be). r_t = R* (M/m)^(1/3) under
*       RTMODE 1, RACCR under RTMODE 2. IACC = number of removals.
*
      INCLUDE 'common6.h'
      INCLUDE 'tlist.h'
      INCLUDE 'mbh.h'
      INTEGER  LCAND(128),KCAND(128)
      REAL*8  QCAND(128),RCAND(128)
*
      IACC = 0
      NCAND = 0
      RPL = 4.0D0*CMBH/CLIGHT**2
      DO 20 L = 1,NXTLEN
          I = NXTLST(L)
          IF (I.GT.N.OR.I.LT.IFIRST) GO TO 20
          IF (BODY(I).LE.0.0D0) GO TO 20
          ZMU = CMBH + BODY(I)
          RR = SQRT(X0(1,I)**2 + X0(2,I)**2 + X0(3,I)**2)
          V2 = X0DOT(1,I)**2 + X0DOT(2,I)**2 + X0DOT(3,I)**2
          RD = (X0(1,I)*X0DOT(1,I) + X0(2,I)*X0DOT(2,I) +
     &          X0(3,I)*X0DOT(3,I))/RR
          HX = X0(2,I)*X0DOT(3,I) - X0(3,I)*X0DOT(2,I)
          HY = X0(3,I)*X0DOT(1,I) - X0(1,I)*X0DOT(3,I)
          HZ = X0(1,I)*X0DOT(2,I) - X0(2,I)*X0DOT(1,I)
          H2 = HX**2 + HY**2 + HZ**2
          EPS = 0.5D0*V2 - ZMU/RR
          ECC = SQRT(MAX(1.0D0 + 2.0D0*EPS*H2/ZMU**2,0.0D0))
          QP = H2/(ZMU*(1.0D0 + ECC))
          IF (RTMODE.EQ.1) THEN
              RT = RADIUS(I)*(CMBH/BODY(I))**(1.0D0/3.0D0)
          ELSE
              RT = RACCR
          END IF
          RCAP = MAX(RT/BETAF,RPL)
          IF (QP.GE.RCAP) GO TO 20
*       Time between the present position and the pericentre.
          TP = PERITIME(RR,QP,ECC,EPS,ZMU)
          IF ((RD.LT.0.0D0.AND.TP.LE.STEP(I)).OR.
     &        (RD.GE.0.0D0.AND.TP.LE.2.0D0*STEP(I))) THEN
              IF (NCAND.LT.128) THEN
                  NCAND = NCAND + 1
                  LCAND(NCAND) = I
                  QCAND(NCAND) = QP
                  RCAND(NCAND) = RT
                  KCAND(NCAND) = 1
                  IF (RT.LT.RPL*BETAF) KCAND(NCAND) = 2
              END IF
          END IF
   20 CONTINUE
*
*       Removals after the scan (they modify the block list).
      DO 30 L = 1,NCAND
          CALL MBHTDE(LCAND(L),QCAND(L),RCAND(L),KCAND(L))
          IACC = IACC + 1
   30 CONTINUE
*
      RETURN
*
      END
*
*
      REAL*8 FUNCTION PERITIME(RR,QP,ECC,EPS,ZMU)
*
*
*       Two-body time between the distance RR and the pericentre QP.
*       ------------------------------------------------------------
*
*       Elliptic: Kepler's equation in E; hyperbolic: in the hyperbolic
*       anomaly; |e - 1| < 1e-6: Barker's equation.
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
      IF (ABS(ECC - 1.0D0).LT.1.0D-06) THEN
          D = SQRT(MAX(RR/QP - 1.0D0,0.0D0))
          PERITIME = SQRT(2.0D0*QP**3/ZMU)*(D + D**3/3.0D0)
      ELSE IF (ECC.LT.1.0D0) THEN
          SEMI = -0.5D0*ZMU/EPS
          CE = MIN(MAX((1.0D0 - RR/SEMI)/ECC,-1.0D0),1.0D0)
          EA = ACOS(CE)
          PERITIME = (EA - ECC*SIN(EA))*SQRT(SEMI**3/ZMU)
      ELSE
          SEMI = 0.5D0*ZMU/EPS
          CH = MAX((1.0D0 + RR/SEMI)/ECC,1.0D0)
          FA = LOG(CH + SQRT(CH**2 - 1.0D0))
          PERITIME = (ECC*SINH(FA) - FA)*SQRT(SEMI**3/ZMU)
      END IF
*
      RETURN
*
      END
*
*
      SUBROUTINE MBHTDE(I,QP,RT,KIND)
*
*
*       Removal of a star at the black hole (disruption or plunge).
*       -----------------------------------------------------------
*
*       The star's energy in the cluster and black-hole potentials,
*       E = m v^2/2 - m sum m_j/r_ij - m M/r (the other bodies predicted
*       to TIME), is booked into ECOLL, so that the total energy stays
*       continuous when its kinetic and potential terms vanish; its mass
*       leaves ZMASS. The star becomes a ghost as in COAL: zero mass, far
*       from the cluster, in the ghost part of the block list. Stamped
*       line 'TDE' (KIND 1, disruption) or 'TDE PLUNGE' (KIND 2,
*       q < 4 G M/c^2 above r_t): time, name, type, mass, radius, q,
*       r_t, beta, booked energy, black-hole mass.
*
      INCLUDE 'common6.h'
      INCLUDE 'mbh.h'
      REAL*8  XJ(3)
      CHARACTER*10  WHICH
*
      ZM = BODY(I)
      POTI = 0.0D0
      DO 10 J = IFIRST,NTOT
          IF (J.EQ.I.OR.BODY(J).LE.0.0D0) GO TO 10
          DT = TIME - T0(J)
          DO 5 K = 1,3
              XJ(K) = ((FDOT(K,J)*DT + F(K,J))*DT + X0DOT(K,J))*DT +
     &                X0(K,J)
    5     CONTINUE
          RIJ = SQRT((XJ(1) - X0(1,I))**2 + (XJ(2) - X0(2,I))**2 +
     &               (XJ(3) - X0(3,I))**2)
          POTI = POTI - BODY(J)/RIJ
   10 CONTINUE
      R2 = X0(1,I)**2 + X0(2,I)**2 + X0(3,I)**2 + EPSBH**2
      V2 = X0DOT(1,I)**2 + X0DOT(2,I)**2 + X0DOT(3,I)**2
      EI = ZM*(0.5D0*V2 + POTI - CMBH/SQRT(R2))
      ECOLL = ECOLL + EI
      EMBHC = EMBHC + EI
      ZMASS = ZMASS - ZM
      DMTDE = DMTDE + ZM
      IF (KIND.EQ.1) THEN
          NTDE = NTDE + 1
          WHICH = 'TDE       '
      ELSE
          NPLUNGE = NPLUNGE + 1
          WHICH = 'TDE PLUNGE'
      END IF
      BETA = RT/QP
      IF (rank.eq.0) WRITE (6,20)  WHICH, TIME+TOFF, NAME(I), KSTAR(I),
     &    ZM*ZMBAR, RADIUS(I)*SU, QP*RBAR, RT*RBAR, BETA, EI, CMBH*ZMBAR
   20 FORMAT (1X,A10,' T[NB] NAME K* M[M*] R*[R*] Q[PC] RT[PC] BETA ',
     &        'E[NB] MBH[M*]',1P,E17.10,I10,I4,4E13.5,E11.3,2E13.5)
*
*       Ghost as in COAL: zero mass, far away, moved to the ghost list.
      RI = SQRT(R2)
      VI = SQRT(V2)
      BODY(I) = 0.0D0
      RADIUS(I) = 0.0D0
      TEV(I) = 1.0D+10
      T0(I) = TADJ + DTADJ
      CALL DELAY_REMOVE_TLIST(I,STEP,DTK)
      STEP(I) = 2.0D0*DTK(1)
      CALL ADD_TLIST(I,STEP,DTK)
      DO 30 K = 1,3
          X0(K,I) = MIN(1.0D+04 + X0(K,I),1000.0D0*RSCALE*X0(K,I)/RI)
          X(K,I) = X0(K,I)
          X0DOT(K,I) = SQRT(0.004D0*ZMASS/RSCALE)*X0DOT(K,I)/
     &                 MAX(VI,1.0D-20)
          XDOT(K,I) = X0DOT(K,I)
          F(K,I) = 0.0D0
          FDOT(K,I) = 0.0D0
          FI(K,I) = 0.0D0
          FIDOT(K,I) = 0.0D0
          D0(K,I) = 0.0D0
          D1(K,I) = 0.0D0
          D2(K,I) = 0.0D0
          D3(K,I) = 0.0D0
          D0R(K,I) = 0.0D0
          D1R(K,I) = 0.0D0
          D2R(K,I) = 0.0D0
          D3R(K,I) = 0.0D0
   30 CONTINUE
*
      RETURN
*
      END
