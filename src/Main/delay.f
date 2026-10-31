      SUBROUTINE DELAY(KCH,KS)
*
*
*       Delay of multiple regularization.
*       ---------------------------------
*
      INCLUDE 'common6.h'
      COMMON/SAVEIT/  IPH0,KS0,KCH0,KS20,JC0,JCL0,PCR0
*       Names of the pair and intruder of a deferred request and the pair
*       count at its time (verified before the request is carried out).
      COMMON/SAVEIT2/ NMA0,NMB0,NMJ0,NMJ20,NP0
*
*
*       Check for saving KS pointers from IMPACT or termination at TBLOCK.
      IF (KS.GE.0) THEN
          IPH0 = IPHASE
          KS0 = KSPAIR
          KCH0 = KCH
          KS20 = KS
          JC0 = JCOMP
          JCL0 = JCLOSE
          PCR0 = PCRIT
          NMA0 = NAME(2*KSPAIR-1)
          NMB0 = NAME(2*KSPAIR)
          NMJ0 = 0
          NMJ20 = 0
          IF (JCOMP.GT.N) THEN
              NMJ0 = NAME(2*(JCOMP-N)-1)
              NMJ20 = NAME(2*(JCOMP-N))
          ELSE IF (JCOMP.GT.0) THEN
              NMJ0 = NAME(JCOMP)
          END IF
          NP0 = NPAIRS
      ELSE
*       Copy relevant indices for delayed KS termination.
          IPHASE = IPH0
          KSPAIR = KS0
          KCHAIN = KCH0
          KS2 = KS20
          JCOMP = JC0
          JCLOSE = JCL0
          PCRIT = PCR0
*      Preserve contents of KSAVE during chain regularization.
          IF (NCH.EQ.0) THEN
              KSAVE(1) = 0
              KSAVE(2) = 0
          END IF
          KSKIP = 0
*
*      Exit in case of new merger or merger termination.
          IF (IPHASE.EQ.6.OR.IPHASE.EQ.7) THEN
              GO TO 10
          END IF
*       Abandon a request whose pair or intruder no longer sits at the
*       saved index, or whose pair count changed since (a termination or
*       regularization in between renumbers the pairs): a chain of other
*       bodies was started in that case, with the c.m. of the intended
*       pair, and the run halted. IMPACT raises the request again.
          IFAIL = 0
          IF (KSPAIR.LT.1.OR.KSPAIR.GT.NPAIRS) THEN
              IFAIL = 1
          ELSE IF (NAME(2*KSPAIR-1).NE.NMA0.OR.NAME(2*KSPAIR).NE.NMB0)
     &        THEN
              IFAIL = 1
          END IF
          IF (JCOMP.GT.N) THEN
              IF (JCOMP-N.GT.NPAIRS) THEN
                  IFAIL = 2
              ELSE IF (NAME(2*(JCOMP-N)-1).NE.NMJ0.OR.
     &                 NAME(2*(JCOMP-N)).NE.NMJ20) THEN
                  IFAIL = 2
              END IF
          ELSE IF (JCOMP.GT.0) THEN
              IF (JCOMP.GT.NTOT) THEN
                  IFAIL = 3
              ELSE IF (NAME(JCOMP).NE.NMJ0) THEN
                  IFAIL = 3
              END IF
          END IF
          IF (NPAIRS.NE.NP0) IFAIL = 4
          IF (IFAIL.GT.0) THEN
              if(rank.eq.0) WRITE (6,2) IFAIL, IPHASE, NMA0, NMB0, NMJ0,
     &                                 NMJ20, NP0, NPAIRS
    2         FORMAT (' DELAY: REQUEST ABANDONED (RENUMBERED)  CASE',
     &                ' IPHASE NAMES NP0 NPAIRS ',2I3,4I8,2I7)
              IPHASE = 0
              GO TO 10
          END IF
*
*       Include the case of two interacting KS pairs.
          IF (JCOMP.GT.N) THEN
              IF (KCHAIN.GT.0.AND.KSTAR(N+KSPAIR).NE.0) THEN
                  KSAVE(1) = KSTAR(N+KSPAIR)
                  KSAVE(2) = NAME(2*KSPAIR-1) + NAME(2*KSPAIR)
                  KSKIP = 1
              END IF
*       Terminate smallest pair first and copy second pair index.
              CALL KSTERM
              KSPAIR = KS2
*       Specify JCOMP < 0 to prevent spurious prediction second KSTERM call.
              JCOMP = -1
              IF(KZ(26).LT.2) JCLOSE = 0
          END IF
*
*       Save KSTAR (> 0) and sum of component names (for chain termination).
      IF (KCHAIN.GT.0.AND.KSTAR(N+KSPAIR).NE.0.AND.KSKIP.EQ.0) THEN
          KSAVE(1) = KSTAR(N+KSPAIR)
          KSAVE(2) = NAME(2*KSPAIR-1) + NAME(2*KSPAIR)
      END IF
*
*       Terminate binary in triple or widest binary-binary collision pair.
          CALL KSTERM
*
*       See whether chain regularization indicator should be switched on.
          IF (KCHAIN.GT.0) THEN
              IPHASE = 8
          END IF
      END IF
*
   10 RETURN
*
      END
