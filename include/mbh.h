*       Static massive black hole at the origin (KZ(44) = 1; mbh.f).
*       CMBH mass, EPSBH Plummer softening, RACCR fixed disruption
*       radius (N-body units); RTMODE 1: r_t = R* (M/m)^(1/3) with the
*       stellar radius (KZ(19) >= 3), 2: r_t = RACCR for every star.
*       BETAF: disruption when r_t/q >= BETAF. FACC: fraction of the
*       disrupted mass counted as accreted (the mass itself is static).
*       NTDE, NPLUNGE: events; DMTDE: mass removed; EMBHC: energy booked.
      REAL*8  CMBH,EPSBH,RACCR,FACC,BETAF,DMTDE,EMBHC
      INTEGER  RTMODE,NTDE,NPLUNGE
      COMMON/MBHC/  CMBH,EPSBH,RACCR,FACC,BETAF,DMTDE,EMBHC,
     &              RTMODE,NTDE,NPLUNGE
