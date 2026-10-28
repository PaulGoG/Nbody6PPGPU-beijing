*       Post-Newtonian terms of KS pairs (KZ(43); routines in kspn.f).
*       FPN, FDPN: PN relative acceleration and its time derivative of
*       the current step (all selected orders); FRR, FDRR: their 2.5PN
*       part. LPNSTEP = 1 while the current KSINT step carries them.
*       EPNW, EPNRR: work booked since the start of this process (all
*       orders, 2.5PN part). NPNSTP counts PN steps. EPSPN: smallest
*       pericentre compactness m/(c^2 q) at which the terms act.
      REAL*8  FPN,FDPN,FRR,FDRR,EPNW,EPNRR,EPSPN
      INTEGER  LPNSTEP,NPNSTP
      PARAMETER (EPSPN = 1.0D-06)
      COMMON/KSPNC/  FPN(3),FDPN(3),FRR(3),FDRR(3),EPNW,EPNRR,
     &               LPNSTEP,NPNSTP
