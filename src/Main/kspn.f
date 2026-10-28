      SUBROUTINE KSPNACC(XR,VR,ZM1,ZM2,FEXT,IPN,AP,APD,AR,ARD)
*
*
*       Post-Newtonian relative acceleration of a two-body system.
*       ----------------------------------------------------------
*
*       Harmonic coordinates, centre-of-mass frame, order-reduced
*       (Blanchet 2024, Living Rev. Relativ. 27, 4, Eqs. for the
*       coefficients A and B; 10.1007/s41114-024-00050-z):
*           a = -(m/r^2) [(1 + A) n + B v],
*       x = x1 - x2, v = v1 - v2, n = x/r, RD = n.v, m = m1 + m2,
*       nu = m1 m2/m^2, G = 1, c = CLIGHT. IPN selects the orders by
*       bits: 1 = 1PN, 2 = 2PN, 4 = 2.5PN (radiation reaction).
*       AP is the PN part of the acceleration, APD its time derivative
*       along the motion with the total acceleration -m x/r^3 + FEXT
*       (FEXT the relative perturbation); AR, ARD the 2.5PN part alone.
*
      INCLUDE 'common6.h'
      REAL*8  XR(3),VR(3),FEXT(3),AP(3),APD(3),AR(3),ARD(3),AT(3),
     &        UN(3),UND(3)
      INTEGER  IPN
*
      ZM = ZM1 + ZM2
      ZNU = ZM1*ZM2/ZM**2
      R2 = XR(1)**2 + XR(2)**2 + XR(3)**2
      RR = SQRT(R2)
      RIN = 1.0D0/RR
      RD = (XR(1)*VR(1) + XR(2)*VR(2) + XR(3)*VR(3))*RIN
      W = VR(1)**2 + VR(2)**2 + VR(3)**2
      ZU = ZM*RIN
      CI2 = 1.0D0/CLIGHT**2
      CI4 = CI2*CI2
      CI5 = CI4/CLIGHT
*
*       Total relative acceleration and the derivatives of r, RD, U, W.
      DO 5 K = 1,3
          AT(K) = -ZU*RIN*RIN*XR(K) + FEXT(K)
          UN(K) = XR(K)*RIN
          UND(K) = (VR(K) - RD*UN(K))*RIN
    5 CONTINUE
      RDD = (W - RD**2)*RIN + UN(1)*AT(1) + UN(2)*AT(2) + UN(3)*AT(3)
      UD = -ZU*RD*RIN
      WD = 2.0D0*(VR(1)*AT(1) + VR(2)*AT(2) + VR(3)*AT(3))
*
*       Coefficients and their partial derivatives in U, RD and W.
      PA = 0.0D0
      PAU = 0.0D0
      PAR = 0.0D0
      PAW = 0.0D0
      PB = 0.0D0
      PBU = 0.0D0
      PBR = 0.0D0
      PBW = 0.0D0
      IF (MOD(IPN,2).EQ.1) THEN
          PA = PA + CI2*(-1.5D0*ZNU*RD**2 + (1.0D0 + 3.0D0*ZNU)*W
     &                 - (4.0D0 + 2.0D0*ZNU)*ZU)
          PAU = PAU - CI2*(4.0D0 + 2.0D0*ZNU)
          PAR = PAR - CI2*3.0D0*ZNU*RD
          PAW = PAW + CI2*(1.0D0 + 3.0D0*ZNU)
          PB = PB + CI2*(-4.0D0 + 2.0D0*ZNU)*RD
          PBR = PBR + CI2*(-4.0D0 + 2.0D0*ZNU)
      END IF
      IF (MOD(IPN/2,2).EQ.1) THEN
          C1 = 1.875D0*ZNU - 5.625D0*ZNU**2
          C2 = -4.5D0*ZNU + 6.0D0*ZNU**2
          C3 = 3.0D0*ZNU - 4.0D0*ZNU**2
          C4 = -2.0D0 - 25.0D0*ZNU - 2.0D0*ZNU**2
          C5 = -6.5D0*ZNU + 2.0D0*ZNU**2
          C6 = 9.0D0 + 21.75D0*ZNU
          PA = PA + CI4*(C1*RD**4 + C2*RD**2*W + C3*W**2
     &                 + ZU*(C4*RD**2 + C5*W) + C6*ZU**2)
          PAU = PAU + CI4*(C4*RD**2 + C5*W + 2.0D0*C6*ZU)
          PAR = PAR + CI4*(4.0D0*C1*RD**3 + 2.0D0*C2*RD*W
     &                     + 2.0D0*C4*ZU*RD)
          PAW = PAW + CI4*(C2*RD**2 + 2.0D0*C3*W + C5*ZU)
          DD1 = 4.5D0*ZNU + 3.0D0*ZNU**2
          DD2 = -7.5D0*ZNU - 2.0D0*ZNU**2
          DD3 = 2.0D0 + 20.5D0*ZNU + 4.0D0*ZNU**2
          PB = PB + CI4*(DD1*RD**3 + DD2*RD*W + DD3*ZU*RD)
          PBU = PBU + CI4*DD3*RD
          PBR = PBR + CI4*(3.0D0*DD1*RD**2 + DD2*W + DD3*ZU)
          PBW = PBW + CI4*DD2*RD
      END IF
*       Radiation reaction (2.5PN) kept apart for the energy booking.
      ARR = 0.0D0
      ARRU = 0.0D0
      ARRR = 0.0D0
      ARRW = 0.0D0
      BRR = 0.0D0
      BRRU = 0.0D0
      BRRW = 0.0D0
      IF (MOD(IPN/4,2).EQ.1) THEN
          FRAD = 1.6D0*ZNU*CI5
          ARR = -FRAD*(17.0D0/3.0D0*ZU**2*RD + 3.0D0*ZU*RD*W)
          ARRU = -FRAD*(34.0D0/3.0D0*ZU*RD + 3.0D0*RD*W)
          ARRR = -FRAD*(17.0D0/3.0D0*ZU**2 + 3.0D0*ZU*W)
          ARRW = -FRAD*3.0D0*ZU*RD
          BRR = FRAD*(3.0D0*ZU**2 + ZU*W)
          BRRU = FRAD*(6.0D0*ZU + W)
          BRRW = FRAD*ZU
      END IF
*
*       a = -(U/r) [A n + B v] and its derivative
*       -(U/r) [(A' - 2 RD A/r) n + A n' + (B' - 2 RD B/r) v + B a].
      Q = -ZU*RIN
      CALL KSPNTERM(Q,RIN,RD,UD,RDD,WD,PA,PAU,PAR,PAW,PB,PBU,PBR,PBW,
     &              UN,UND,VR,AT,AP,APD)
      CALL KSPNTERM(Q,RIN,RD,UD,RDD,WD,ARR,ARRU,ARRR,ARRW,BRR,BRRU,
     &              0.0D0,BRRW,UN,UND,VR,AT,AR,ARD)
      DO 10 K = 1,3
          AP(K) = AP(K) + AR(K)
          APD(K) = APD(K) + ARD(K)
   10 CONTINUE
*
      RETURN
*
      END
*
*
      SUBROUTINE KSPNTERM(Q,RIN,RD,UD,RDD,WD,A,AU,AR1,AW,B,BU,BR1,BW,
     &                    UN,UND,VR,AT,AP,APD)
*
*
*       One PN term Q [A n + B v] and its time derivative.
*       --------------------------------------------------
*
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  Q,RIN,RD,UD,RDD,WD,A,AU,AR1,AW,B,BU,BR1,BW,UN(3),UND(3),
     &        VR(3),AT(3),AP(3),APD(3)
*
      AD = AU*UD + AR1*RDD + AW*WD - 2.0D0*RD*A*RIN
      BD = BU*UD + BR1*RDD + BW*WD - 2.0D0*RD*B*RIN
      DO 10 K = 1,3
          AP(K) = Q*(A*UN(K) + B*VR(K))
          APD(K) = Q*(AD*UN(K) + A*UND(K) + BD*VR(K) + B*AT(K))
   10 CONTINUE
*
      RETURN
*
      END
*
*
      SUBROUTINE KSPNADD(IPAIR,I1,XI,VI,FP,FD)
*
*
*       Post-Newtonian terms added to the perturbation of a KS pair.
*       ------------------------------------------------------------
*
*       Under KZ(43) > 0 a bound pair of two neutron stars or black
*       holes (KSTAR 13, 14; any bound pair with bit 8, a test mode)
*       whose pericentre compactness m/(c^2 q) is at least EPSPN gets
*       the selected PN orders (KSPNACC) added to its relative
*       perturbation FP and derivative FD. Called after GAMMA has been
*       set, so that termination, slow-down and perturber selection
*       stay Newtonian. KSINT then skips the orbit-averaged Peters decay
*       of the step (GWSKIP): the radiation reaction is the 2.5PN term. With bit 16 every 64th PN step writes a line
*       with the corrected elements and the booked work (KSPNBK).
*
      INCLUDE 'common6.h'
      INCLUDE 'kspn.h'
      REAL*8  XI(6),VI(6),FP(6),FD(6),XR(3),VR(3),FEXT(3),CX(3)
*
      LPNSTEP = 0
      IPN = MOD(KZ(43),8)
      IF (IPN.EQ.0) RETURN
      I2 = I1 + 1
      IF (MOD(KZ(43)/8,2).EQ.0) THEN
          IF (KSTAR(I1).LT.13.OR.KSTAR(I1).GT.14.OR.
     &        KSTAR(I2).LT.13.OR.KSTAR(I2).GT.14) RETURN
      END IF
      ZM1 = BODY(I1)
      ZM2 = BODY(I2)
      IF (ZM1.LE.0.0D0.OR.ZM2.LE.0.0D0) RETURN
      ZM = ZM1 + ZM2
      DO 5 K = 1,3
          XR(K) = XI(K) - XI(K+3)
          VR(K) = VI(K) - VI(K+3)
          FEXT(K) = FP(K)
    5 CONTINUE
      RR = SQRT(XR(1)**2 + XR(2)**2 + XR(3)**2)
      V2 = VR(1)**2 + VR(2)**2 + VR(3)**2
      HB = 0.5D0*V2 - ZM/RR
      IF (HB.GE.0.0D0) RETURN
      SEMI = -0.5D0*ZM/HB
      CX(1) = XR(2)*VR(3) - XR(3)*VR(2)
      CX(2) = XR(3)*VR(1) - XR(1)*VR(3)
      CX(3) = XR(1)*VR(2) - XR(2)*VR(1)
      ZL2 = CX(1)**2 + CX(2)**2 + CX(3)**2
      ECC = SQRT(MAX(1.0D0 - ZL2/(ZM*SEMI),0.0D0))
      QP = SEMI*(1.0D0 - ECC)
      IF (ZM.LT.EPSPN*CLIGHT**2*QP) RETURN
*
      CALL KSPNACC(XR,VR,ZM1,ZM2,FEXT,IPN,FPN,FDPN,FRR,FDRR)
      DO 10 K = 1,3
          FP(K) = FP(K) + FPN(K)
          FD(K) = FD(K) + FDPN(K)
   10 CONTINUE
      LPNSTEP = 1
      NPNSTP = NPNSTP + 1
*
*
      RETURN
*
      END
*
*
      SUBROUTINE KSPNUPD(IPAIR,UI,UIDOT,RI,FP,FD1,FD)
*
*
*       PN terms of a KS step at the corrected state.
*       ---------------------------------------------
*
*       Called by KSCORR between its iterations: the PN part of FP and
*       FD1 (physical derivative) is replaced by the terms at the
*       relative motion of UI, UIDOT, FD = RI*FD1 (regularized), and
*       FPN, FDPN, FRR, FDRR are updated for the energy booking. The
*       external perturbation enters the derivative as before.
*
      INCLUDE 'common6.h'
      INCLUDE 'kspn.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
      REAL*8  UI(4),UIDOT(4),FP(6),FD1(6),FD(6),A1(3,4),XR(3),VR(3),
     &        FEXT(3),PNN(3),PNDN(3),RRN(3),RRDN(3)
*
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      ZMOD = 1.0D0
      IF (KZ(26).GT.0.AND.KSLOW(IPAIR).GT.1) THEN
          ZMOD = FLOAT(ISLOW(KSLOW(IPAIR)))
      END IF
      XR(1) = UI(1)**2 - UI(2)**2 - UI(3)**2 + UI(4)**2
      XR(2) = 2.0D0*(UI(1)*UI(2) - UI(3)*UI(4))
      XR(3) = 2.0D0*(UI(1)*UI(3) + UI(2)*UI(4))
      CALL MATRIX(UI,A1)
      DO 5 L = 1,3
          VR(L) = 2.0D0*(A1(L,1)*UIDOT(1) + A1(L,2)*UIDOT(2) +
     &                   A1(L,3)*UIDOT(3) + A1(L,4)*UIDOT(4))/RI
          FEXT(L) = FP(L)/ZMOD - FPN(L)
    5 CONTINUE
      CALL KSPNACC(XR,VR,BODY(I1),BODY(I2),FEXT,MOD(KZ(43),8),
     &             PNN,PNDN,RRN,RRDN)
      DO 10 K = 1,3
          FP(K) = FP(K) + ZMOD*(PNN(K) - FPN(K))
          FD1(K) = FD1(K) + ZMOD*(PNDN(K) - FDPN(K))
          FD(K) = RI*FD1(K)
          FPN(K) = PNN(K)
          FDPN(K) = PNDN(K)
          FRR(K) = RRN(K)
          FDRR(K) = RRDN(K)
   10 CONTINUE
*
      RETURN
*
      END
*
*
      SUBROUTINE KSPNBK(IPAIR,UI,UIDOT,FP,RI,HNEW,DTU)
*
*
*       Work of the post-Newtonian terms over one KS step.
*       --------------------------------------------------
*
*       The 1PN and 2PN terms conserve E_N + E_PN (PNENER) to O(c^-6):
*       their work is booked as mu times the change of E_PN between the
*       ends of the step, which sums exactly over steps; a quadrature of
*       their work, oscillating by ~ m/(c^2 r) E_N per orbit, would not.
*       The 2.5PN work enters by the two-point Hermite quadrature of
*       KSCORR, DH = (HD0 + HD1) dtau/2 + (HD20 - HD21) dtau^2/12, with
*       HD = 2 u'.L(u)^T P and HD2 = 2 u''.L(u)^T P + 2 u'.(L(u)^T P'
*       + L(u')^T P) at the corrected end (terms of the step, total FP)
*       and at the start (U0, UDOT, Newtonian acceleration). The sum
*       enters ECOLL and E(10), the 2.5PN part also EGRAV, so that the
*       total energy stays conserved with h Newtonian.
*
      INCLUDE 'common6.h'
      INCLUDE 'kspn.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
      REAL*8  UI(4),UIDOT(4),FP(6),A1(3,4),U2E(4),U2S(4),UZ(4),UZD(4),
     &        XR(3),VR(3),ZERO3(3),P0(3),PD0(3),QR0(3),QRD0(3),PNENER,
     &        GV(3)
*
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      ZM1 = BODY(I1)
      ZM2 = BODY(I2)
      ZMU = ZM1*ZM2/(ZM1 + ZM2)
      ZMOD = 1.0D0
      IF (KZ(26).GT.0.AND.KSLOW(IPAIR).GT.1) THEN
          ZMOD = FLOAT(ISLOW(KSLOW(IPAIR)))
      END IF
*
*       End of the step: 2u'' = h u + r L^T P (P total, slow-down incl.).
      CALL MATRIX(UI,A1)
      DO 5 K = 1,4
          U2E(K) = HNEW*UI(K) + RI*(A1(1,K)*FP(1) + A1(2,K)*FP(2) +
     &                              A1(3,K)*FP(3))
    5 CONTINUE
      CALL KSPNHD(UI,UIDOT,U2E,R(IPAIR),ZMOD,FRR,FDRR,GE,GE2)
*
*       Start of the step: relative motion from U0, UDOT.
      DO 10 K = 1,4
          UZ(K) = U0(K,IPAIR)
          UZD(K) = UDOT(K,IPAIR)
          U2S(K) = H0(IPAIR)*UZ(K) + 2.0D0*FP0(K,IPAIR)
   10 CONTINUE
      R0S = UZ(1)**2 + UZ(2)**2 + UZ(3)**2 + UZ(4)**2
      XR(1) = UZ(1)**2 - UZ(2)**2 - UZ(3)**2 + UZ(4)**2
      XR(2) = 2.0D0*(UZ(1)*UZ(2) - UZ(3)*UZ(4))
      XR(3) = 2.0D0*(UZ(1)*UZ(3) + UZ(2)*UZ(4))
      CALL MATRIX(UZ,A1)
      DO 15 L = 1,3
          VR(L) = 2.0D0*(A1(L,1)*UZD(1) + A1(L,2)*UZD(2) +
     &                   A1(L,3)*UZD(3) + A1(L,4)*UZD(4))/R0S
          ZERO3(L) = 0.0D0
   15 CONTINUE
      CALL KSPNACC(XR,VR,ZM1,ZM2,ZERO3,MOD(KZ(43),8),P0,PD0,QR0,QRD0)
      CALL KSPNHD(UZ,UZD,U2S,R0S,ZMOD,QR0,QRD0,GS,GS2)
      ES = PNENER(XR,VR,ZM1,ZM2,MOD(KZ(43),8))
      CALL PNEGRAD(XR,VR,ZM1,ZM2,MOD(KZ(43),8),GV)
      CRS = R0S*(GV(1)*QR0(1) + GV(2)*QR0(2) + GV(3)*QR0(3))
*
*       End state for the conservative energy.
      XR(1) = UI(1)**2 - UI(2)**2 - UI(3)**2 + UI(4)**2
      XR(2) = 2.0D0*(UI(1)*UI(2) - UI(3)*UI(4))
      XR(3) = 2.0D0*(UI(1)*UI(3) + UI(2)*UI(4))
      CALL MATRIX(UI,A1)
      DO 20 L = 1,3
          VR(L) = 2.0D0*(A1(L,1)*UIDOT(1) + A1(L,2)*UIDOT(2) +
     &                   A1(L,3)*UIDOT(3) + A1(L,4)*UIDOT(4))/RI
   20 CONTINUE
      EE = PNENER(XR,VR,ZM1,ZM2,MOD(KZ(43),8))
      CALL PNEGRAD(XR,VR,ZM1,ZM2,MOD(KZ(43),8),GV)
      CRE = RI*(GV(1)*FRR(1) + GV(2)*FRR(2) + GV(3)*FRR(3))
*
*       Conservative orders: the change of their energy (exact sum over
*       steps) less its part due to the 2.5PN force, dE_PN/dv . a_RR
*       (an O(c^-7) cross term, trapezoidal rule in tau); 2.5PN: the
*       Hermite quadrature of its work.
      DHRR = 0.5D0*(GS + GE)*DTU + (GS2 - GE2)*DTU**2/12.0D0
      DERR = -ZMU*DHRR
      DEPN = ZMU*(EE - ES - 0.5D0*(CRS + CRE)*DTU) + DERR
      ECOLL = ECOLL + DEPN
      E(10) = E(10) + DEPN
      EGRAV = EGRAV + DERR
      EPNW = EPNW + DEPN
      EPNRR = EPNRR + DERR
*
*       Diagnostic line at the corrected end of the step (KZ(43) bit 16):
*       elements, eccentricity vector, distance, h, booked work.
      IF (MOD(KZ(43)/16,2).EQ.1.AND.MOD(NPNSTP,64).EQ.1.AND.
     &    rank.eq.0) THEN
          XR(1) = UI(1)**2 - UI(2)**2 - UI(3)**2 + UI(4)**2
          XR(2) = 2.0D0*(UI(1)*UI(2) - UI(3)*UI(4))
          XR(3) = 2.0D0*(UI(1)*UI(3) + UI(2)*UI(4))
          CALL MATRIX(UI,A1)
          DO 25 L = 1,3
              VR(L) = 2.0D0*(A1(L,1)*UIDOT(1) + A1(L,2)*UIDOT(2) +
     &                       A1(L,3)*UIDOT(3) + A1(L,4)*UIDOT(4))/RI
   25     CONTINUE
          ZM = ZM1 + ZM2
          RR = SQRT(XR(1)**2 + XR(2)**2 + XR(3)**2)
          HB = 0.5D0*(VR(1)**2 + VR(2)**2 + VR(3)**2) - ZM/RR
          SEMI = -0.5D0*ZM/HB
          CX1 = XR(2)*VR(3) - XR(3)*VR(2)
          CX2 = XR(3)*VR(1) - XR(1)*VR(3)
          CX3 = XR(1)*VR(2) - XR(2)*VR(1)
          ECC = SQRT(MAX(1.0D0 - (CX1**2 + CX2**2 + CX3**2)/
     &                   (ZM*SEMI),0.0D0))
          EV1 = (VR(2)*CX3 - VR(3)*CX2)/ZM - XR(1)/RR
          EV2 = (VR(3)*CX1 - VR(1)*CX3)/ZM - XR(2)/RR
          EV3 = (VR(1)*CX2 - VR(2)*CX1)/ZM - XR(3)/RR
          WRITE (6,30)  TIME+TOFF, IPAIR, NAME(I1), NAME(I2), SEMI,
     &                  ECC, EV1, EV2, EV3, RR, HB, EPNW, EPNRR
   30     FORMAT (' PN KSINT T IP NAM A E EV R H EPN ERR ',1P,E24.16,
     &            I7,2I10,E24.16,E20.12,3E20.12,E15.7,E20.12,2E15.7)
      END IF
*
      RETURN
*
      END
*
*
      SUBROUTINE KSPNHD(UI,UIDOT,U2X2,RF,ZMOD,P,PD,HD,HD2)
*
*
*       PN part of dh/dtau and of its derivative at one KS state.
*       ---------------------------------------------------------
*
*       U2X2 = 2u'', RF the distance converting dP/dt to dP/dtau, ZMOD
*       the slow-down factor applied to the perturbation.
*
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  UI(4),UIDOT(4),U2X2(4),P(3),PD(3),A1(3,4),A2(3,4)
*
      CALL MATRIX(UI,A1)
      CALL MATRIX(UIDOT,A2)
      HD = 0.0D0
      HD2 = 0.0D0
      DO 10 K = 1,4
          AK = ZMOD*(A1(1,K)*P(1) + A1(2,K)*P(2) + A1(3,K)*P(3))
          AK4 = ZMOD*(RF*(A1(1,K)*PD(1) + A1(2,K)*PD(2) +
     &                    A1(3,K)*PD(3)) +
     &                A2(1,K)*P(1) + A2(2,K)*P(2) + A2(3,K)*P(3))
          HD = HD + 2.0D0*UIDOT(K)*AK
          HD2 = HD2 + U2X2(K)*AK + 2.0D0*UIDOT(K)*AK4
   10 CONTINUE
*
      RETURN
*
      END
*
*
      REAL*8 FUNCTION PNENER(XR,VR,ZM1,ZM2,IPN)
*
*
*       Conservative PN energy of a two-body system per reduced mass.
*       -------------------------------------------------------------
*
*       Centre-of-mass frame, harmonic coordinates, the 1PN (bit 1) and
*       2PN (bit 2) parts of E/mu (Kidder 1995, Phys. Rev. D 52, 821,
*       Eq. 2.7; 10.1103/PhysRevD.52.821), G = 1, c = CLIGHT; with the
*       accelerations of KSPNACC, E_N + E_PN is conserved to O(c^-6).
*
      INCLUDE 'common6.h'
      REAL*8  XR(3),VR(3)
      INTEGER  IPN
*
      ZM = ZM1 + ZM2
      ZNU = ZM1*ZM2/ZM**2
      RR = SQRT(XR(1)**2 + XR(2)**2 + XR(3)**2)
      ZU = ZM/RR
      RD = (XR(1)*VR(1) + XR(2)*VR(2) + XR(3)*VR(3))/RR
      W = VR(1)**2 + VR(2)**2 + VR(3)**2
      CI2 = 1.0D0/CLIGHT**2
      EPN = 0.0D0
      IF (MOD(IPN,2).EQ.1) THEN
          EPN = EPN + CI2*(0.375D0*(1.0D0 - 3.0D0*ZNU)*W**2 +
     &          0.5D0*ZU*((3.0D0 + ZNU)*W + ZNU*RD**2 + ZU))
      END IF
      IF (MOD(IPN/2,2).EQ.1) THEN
          EPN = EPN + CI2**2*(
     &        0.3125D0*(1.0D0 - 7.0D0*ZNU + 13.0D0*ZNU**2)*W**3 +
     &        0.125D0*ZU*((21.0D0 - 23.0D0*ZNU - 27.0D0*ZNU**2)*W**2
     &          + 2.0D0*ZNU*(1.0D0 - 15.0D0*ZNU)*W*RD**2
     &          - 3.0D0*ZNU*(1.0D0 - 3.0D0*ZNU)*RD**4) +
     &        0.125D0*ZU**2*((14.0D0 - 55.0D0*ZNU + 4.0D0*ZNU**2)*W
     &          + (4.0D0 + 69.0D0*ZNU + 12.0D0*ZNU**2)*RD**2) -
     &        0.25D0*ZU**3*(2.0D0 + 15.0D0*ZNU))
      END IF
      PNENER = EPN
*
      RETURN
*
      END
*
*
      SUBROUTINE PNEGRAD(XR,VR,ZM1,ZM2,IPN,GV)
*
*
*       Velocity gradient of the conservative PN energy (PNENER).
*       ---------------------------------------------------------
*
*       dE_PN/dv = 2 (dE/dW) v + (dE/dRD) n with W = v.v, RD = n.v.
*
      INCLUDE 'common6.h'
      REAL*8  XR(3),VR(3),GV(3)
      INTEGER  IPN
*
      ZM = ZM1 + ZM2
      ZNU = ZM1*ZM2/ZM**2
      RR = SQRT(XR(1)**2 + XR(2)**2 + XR(3)**2)
      ZU = ZM/RR
      RD = (XR(1)*VR(1) + XR(2)*VR(2) + XR(3)*VR(3))/RR
      W = VR(1)**2 + VR(2)**2 + VR(3)**2
      CI2 = 1.0D0/CLIGHT**2
      EW = 0.0D0
      ER = 0.0D0
      IF (MOD(IPN,2).EQ.1) THEN
          EW = EW + CI2*(0.75D0*(1.0D0 - 3.0D0*ZNU)*W +
     &         0.5D0*ZU*(3.0D0 + ZNU))
          ER = ER + CI2*ZU*ZNU*RD
      END IF
      IF (MOD(IPN/2,2).EQ.1) THEN
          EW = EW + CI2**2*(
     &        0.9375D0*(1.0D0 - 7.0D0*ZNU + 13.0D0*ZNU**2)*W**2 +
     &        0.125D0*ZU*(2.0D0*(21.0D0 - 23.0D0*ZNU - 27.0D0*ZNU**2)*W
     &          + 2.0D0*ZNU*(1.0D0 - 15.0D0*ZNU)*RD**2) +
     &        0.125D0*ZU**2*(14.0D0 - 55.0D0*ZNU + 4.0D0*ZNU**2))
          ER = ER + CI2**2*(
     &        0.125D0*ZU*(4.0D0*ZNU*(1.0D0 - 15.0D0*ZNU)*W*RD
     &          - 12.0D0*ZNU*(1.0D0 - 3.0D0*ZNU)*RD**3) +
     &        0.25D0*ZU**2*(4.0D0 + 69.0D0*ZNU + 12.0D0*ZNU**2)*RD)
      END IF
      DO 10 K = 1,3
          GV(K) = 2.0D0*EW*VR(K) + ER*XR(K)/RR
   10 CONTINUE
*
      RETURN
*
      END
