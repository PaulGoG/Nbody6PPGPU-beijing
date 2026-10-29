*       Static massive black hole at the origin (KZ(44) = 1; mbh.F).
*       CMBH mass, EPSBH Plummer softening, RACCR fixed disruption
*       radius (N-body units); RTMODE 1: r_t = R* (M/m)^(1/3) with the
*       stellar radius (KZ(19) >= 3), 2: r_t = RACCR for every star.
*       IPART 0: full disruption when r_t/q >= BETAF; 1, 2: partial and
*       full disruption by the stripped fraction of Guillochon & Ramirez-
*       Ruiz (2013) for gamma = 4/3, 5/3 (BETAF unused). FACC: fraction of
*       the disrupted mass the black hole would gain (growth not yet
*       implemented). NTDE, NPART, NPLUNGE: events; DMTDE: mass removed;
*       EMBHC: energy booked; QPART(NAME): pericentre of a star's last
*       partial event, cleared once the star is beyond 3 QPART again.
      REAL*8  CMBH,EPSBH,RACCR,FACC,BETAF,DMTDE,EMBHC,QPART
      INTEGER  RTMODE,IPART,NTDE,NPART,NPLUNGE
      COMMON/MBHC/  CMBH,EPSBH,RACCR,FACC,BETAF,DMTDE,EMBHC,
     &              RTMODE,IPART,NTDE,NPART,NPLUNGE
      COMMON/MBHT/  QPART(NMAX)
