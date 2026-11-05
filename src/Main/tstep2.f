      REAL*8 FUNCTION TSTEP2(XDOT,F,FDOT,ETA)
*
*
*       General time-step criterion by xdot, f and fdot.
*       ----------------------------
*
      IMPLICIT NONE
      REAL*8  XDOT(3),F(3),FDOT(3)
      REAL*8 F2,FDOT2,V2,ETA,DT2
*
*       Obtain new integration step using composite expression.
*       STEP = (ETA*(F*F2DOT + FDOT**2)/(FDOT*F3DOT + F2DOT**2))**0.5.
*
      V2 = XDOT(1)**2 + XDOT(2)**2 + XDOT(3)**2
      F2 = F(1)**2 + F(2)**2 + F(3)**2
      FDOT2 = FDOT(1)**2 + FDOT(2)**2 + FDOT(3)**2
*
*       Squared time scales: (v/F)**2 and (F/FDOT)**2. A body at (or near)
*       rest has no velocity time scale: v/F below 1e-2 F/FDOT leaves the
*       force criterion alone (a zero velocity gave a zero step). Without
*       a jerk only the velocity scale remains; without both (a body at
*       rest under a constant force) no scale exists here and the caller
*       bounds the step (STEPS2).
      DT2 = 0.0D0
      IF (F2.GT.0.0D0) THEN
          IF (FDOT2.GT.0.0D0) THEN
              DT2 = F2/FDOT2
              IF (V2*FDOT2.GE.1.0D-04*F2**2) DT2 = MIN(V2/F2,DT2)
          ELSE IF (V2.GT.0.0D0) THEN
              DT2 = V2/F2
          END IF
      END IF
      IF (DT2.GT.0.0D0) THEN
          TSTEP2 = 0.25*ETA*SQRT(DT2)
      ELSE
          TSTEP2 = 1.0D10
      END IF
*
      RETURN
*
      END
