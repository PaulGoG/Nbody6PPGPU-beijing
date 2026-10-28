      REAL*8 FUNCTION TSTEP2(XDOT,F,FDOT,ETA)
*
*
*       General time-step criterion by xdot, f and fdot.
*       ----------------------------
*
      IMPLICIT NONE
      REAL*8  XDOT(3),F(3),FDOT(3)
      REAL*8 F2,FDOT2,V2,ETA
*
*       Obtain new integration step using composite expression.
*       STEP = (ETA*(F*F2DOT + FDOT**2)/(FDOT*F3DOT + F2DOT**2))**0.5.
*
      V2 = XDOT(1)**2 + XDOT(2)**2 + XDOT(3)**2
      F2 = F(1)**2 + F(2)**2 + F(3)**2
      FDOT2 = FDOT(1)**2 + FDOT(2)**2 + FDOT(3)**2
*
*       A body at (or near) rest has no velocity time scale: v/F below
*       1e-2 F/FDOT leaves the force criterion alone (a zero velocity
*       gave a zero step).
      IF (V2*FDOT2.LT.1.0D-04*F2**2) THEN
          TSTEP2 = 0.25*ETA*SQRT(F2/FDOT2)
      ELSE
          TSTEP2 = 0.25*ETA*SQRT(MIN(V2/F2,F2/FDOT2))
      END IF
*
      RETURN
*
      END
