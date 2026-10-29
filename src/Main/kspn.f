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
      SUBROUTINE KSPNBK(IPAIR,UI,UIDOT,RI,DTU,HP1,H2P1,HR1,H2R1,
     &                  HD1,HD21)
*
*
*       Work of the post-Newtonian terms over one KS step.
*       --------------------------------------------------
*
*       KSCORR advances h by the two-point Hermite rule
*       DH = (HD0 + HD1) dtau/2 + (HD20 - HD21) dtau^2/12, linear in the
*       four derivatives. Its part from the PN terms (start values HPNS,
*       end values HP1, H2P1 of all orders; HRRS, HR1, H2R1 of 2.5PN) is
*       booked as their work, into ECOLL and E(10), the 2.5PN part also
*       into EGRAV, so that the energy check sees the h change of the
*       Newtonian perturbation alone, exactly. The start values hold
*       while HDOT and HDOT2 are those this routine saw at the end of
*       the last step (HPKEY); KSPOLY recomputes them without the PN
*       terms. HD1, HD21: the new HDOT, HDOT2. RESOLV books the PN part
*       of a prediction into the step (DHPNP) when it ends the step early,
*       and ENERGY leaves it out of EBIN. With KZ(43) bit 16 every 64th PN
*       step writes the corrected elements and the booked work.
*
      INCLUDE 'common6.h'
      INCLUDE 'kspn.h'
      REAL*8  UI(4),UIDOT(4),A1(3,4),XR(3),VR(3)
*
      HP0 = 0.0D0
      H2P0 = 0.0D0
      HR0 = 0.0D0
      H2R0 = 0.0D0
      IF (HPKEY(1,IPAIR).EQ.HDOT(IPAIR).AND.
     &    HPKEY(2,IPAIR).EQ.HDOT2(IPAIR)) THEN
          HP0 = HPNS(1,IPAIR)
          H2P0 = HPNS(2,IPAIR)
          HR0 = HRRS(1,IPAIR)
          H2R0 = HRRS(2,IPAIR)
      END IF
*       New values and the PN parts of HDOT3, HDOT4 as KSCORR forms them.
      HPNS(1,IPAIR) = HP1
      HPNS(2,IPAIR) = H2P1
      HRRS(1,IPAIR) = HR1
      HRRS(2,IPAIR) = H2R1
      CALL KSPNH34(HP0,HP1,H2P0,H2P1,DTU,HPNS(3,IPAIR),HPNS(4,IPAIR))
      CALL KSPNH34(HR0,HR1,H2R0,H2R1,DTU,HRRS(3,IPAIR),HRRS(4,IPAIR))
      HPKEY(1,IPAIR) = HD1
      HPKEY(2,IPAIR) = HD21
      IF (HP0.EQ.0.0D0.AND.H2P0.EQ.0.0D0.AND.HP1.EQ.0.0D0.AND.
     &    H2P1.EQ.0.0D0) RETURN
*
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      ZM1 = BODY(I1)
      ZM2 = BODY(I2)
      ZMU = ZM1*ZM2/(ZM1 + ZM2)
      DHPN = 0.5D0*(HP0 + HP1)*DTU + (H2P0 - H2P1)*DTU**2/12.0D0
      DHRR = 0.5D0*(HR0 + HR1)*DTU + (H2R0 - H2R1)*DTU**2/12.0D0
      DEPN = -ZMU*DHPN
      DERR = -ZMU*DHRR
      ECOLL = ECOLL + DEPN
      E(10) = E(10) + DEPN
      EGRAV = EGRAV + DERR
      EPNW = EPNW + DEPN
      EPNRR = EPNRR + DERR
*
*       Diagnostic line at the corrected end of the step (KZ(43) bit 16):
*       elements, eccentricity vector, distance, h, booked work.
      IF (MOD(KZ(43)/16,2).EQ.1.AND.MOD(NPNSTP,64).EQ.1.AND.
     &    LPNSTEP.EQ.1.AND.rank.eq.0) THEN
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
      SUBROUTINE KSPNH34(HD0,HD1,HD20,HD21,DTU,HD3,HD4)
*
*
*       Third and fourth derivatives of h as KSCORR forms them.
*       -------------------------------------------------------
*
*       Linear in the four end values, so applied to their PN parts it
*       gives the PN parts of HDOT3 and HDOT4.
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
      DHD = HD0 - HD1
      SUM = HD20 + HD21
      BT2 = -3.0D0*DHD - (SUM + HD20)*DTU
      AT3 = 2.0D0*DHD + SUM*DTU
      HD3 = (3.0D0*AT3 + BT2)*2.0D0/DTU**2
      HD4 = AT3*6.0D0/DTU**3
*
      RETURN
*
      END
*
*
      SUBROUTINE KSPNPRD(IPAIR,KDUM,DTU)
*
*
*       PN part of the h prediction of RESOLV.
*       --------------------------------------
*
*       DHPNP, DHRRP: the parts of HT = H + Taylor series in HDOT..HDOT4
*       due to the PN terms (stored by KSPNBK, valid while HPKEY holds).
*       With KDUM > 1 the prediction replaces the step (termination, or
*       H = HT in MDOT): the part is booked like a step and the stored
*       values cleared, so that it is booked once.
*
      INCLUDE 'common6.h'
      INCLUDE 'kspn.h'
*
      DHPNP = 0.0D0
      DHRRP = 0.0D0
      IF (HPKEY(1,IPAIR).NE.HDOT(IPAIR).OR.
     &    HPKEY(2,IPAIR).NE.HDOT2(IPAIR)) RETURN
      DTU2 = DTU/24.0D0
      DHPNP = (((HPNS(4,IPAIR)*DTU2 + ONE6*HPNS(3,IPAIR))*DTU +
     &         0.5D0*HPNS(2,IPAIR))*DTU + HPNS(1,IPAIR))*DTU
      DHRRP = (((HRRS(4,IPAIR)*DTU2 + ONE6*HRRS(3,IPAIR))*DTU +
     &         0.5D0*HRRS(2,IPAIR))*DTU + HRRS(1,IPAIR))*DTU
      IF (KDUM.LT.2.OR.DHPNP.EQ.0.0D0) RETURN
*
      I1 = 2*IPAIR - 1
      ZMU = BODY(I1)*BODY(I1+1)/(BODY(I1) + BODY(I1+1))
      DEPN = -ZMU*DHPNP
      DERR = -ZMU*DHRRP
      ECOLL = ECOLL + DEPN
      E(10) = E(10) + DEPN
      EGRAV = EGRAV + DERR
      EPNW = EPNW + DEPN
      EPNRR = EPNRR + DERR
      DO 10 K = 1,4
          HPNS(K,IPAIR) = 0.0D0
          HRRS(K,IPAIR) = 0.0D0
   10 CONTINUE
*
      RETURN
*
      END
