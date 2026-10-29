*       Post-Newtonian terms of KS pairs (KZ(43); routines in kspn.f).
*       FPN, FDPN: PN relative acceleration and its time derivative of
*       the current step (all selected orders); FRR, FDRR: their 2.5PN
*       part. LPNSTEP = 1 while the current KSINT step carries them.
*       EPNW, EPNRR: work booked since the start of the run (all
*       orders, 2.5PN part). NPNSTP counts PN steps. EPSPN: smallest
*       pericentre compactness m/(c^2 q) at which the terms act.
*       HPNS(1:4,IPAIR), HRRS(1:4,IPAIR): PN and 2.5PN parts of HDOT,
*       HDOT2, HDOT3 and HDOT4 at the end of the pair's last corrector
*       step, valid while HPKEY(1:2,IPAIR) equal HDOT and HDOT2 (KSPOLY
*       recomputes these without the PN terms). DHPNP, DHRRP: their part
*       of the h prediction of the last RESOLV call.
      REAL*8  FPN,FDPN,FRR,FDRR,EPNW,EPNRR,EPSPN,HPNS,HRRS,HPKEY,
     &        DHPNP,DHRRP
      INTEGER  LPNSTEP,NPNSTP
      PARAMETER (EPSPN = 1.0D-06)
      COMMON/KSPNC/  FPN(3),FDPN(3),FRR(3),FDRR(3),EPNW,EPNRR,
     &               LPNSTEP,NPNSTP
      COMMON/KSPNS/  HPNS(4,KMAX),HRRS(4,KMAX),HPKEY(2,KMAX),DHPNP,
     &               DHRRP
