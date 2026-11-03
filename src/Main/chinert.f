      LOGICAL FUNCTION CHINERT(IPAIR)
*
*
*       Chain-inert eligibility of a KS pair.
*       -------------------------------------
*
*       A bound pair whose semi-major axis is below CHINRT*RMIN cannot be
*       resolved consistently by the chain integrator next to members at
*       the chain scale (geometric growth of the chain inconsistency,
*       lifetime cuts); such a pair enters a chain as a dormant c.m. body.
*
      INCLUDE 'common6.h'
      INCLUDE 'chinert.h'
*
*
      CHINERT = .FALSE.
      IF (IPAIR.LT.1.OR.IPAIR.GT.NPAIRS) RETURN
      IF (H(IPAIR).GE.0.0D0.OR.BODY(N+IPAIR).LE.0.0D0) RETURN
      SEMI = -0.5D0*BODY(N+IPAIR)/H(IPAIR)
      IF (SEMI.LT.CHINRT*RMIN) CHINERT = .TRUE.
*
      RETURN
*
      END


      SUBROUTINE CHDORMA(NAMEI,B1,B2,SEMI)
*
*
*       Registration of a dormant pair.
*       -------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'chinert.h'
*
*
*       Replace an existing entry of the same c.m. name.
      DO 1 K = 1,NDORM
          IF (NMDORM(K).EQ.NAMEI) THEN
              BDORM(1,K) = B1
              BDORM(2,K) = B2
              ADORM(K) = SEMI
              RETURN
          END IF
    1 CONTINUE
      IF (NDORM.GE.2) THEN
          if(rank.eq.0) WRITE (6,2)  NAMEI, NMDORM(1), NMDORM(2)
    2     FORMAT (' CHAIN INERT: REGISTRY FULL    NM NMDORM ',3I8)
          RETURN
      END IF
      NDORM = NDORM + 1
      NMDORM(NDORM) = NAMEI
      BDORM(1,NDORM) = B1
      BDORM(2,NDORM) = B2
      ADORM(NDORM) = SEMI
*
      RETURN
*
      END


      SUBROUTINE CHDORMR(NAMEI,B1,B2,IFOUND)
*
*
*       Removal of a dormant pair from the registry.
*       --------------------------------------------
*
      INCLUDE 'common6.h'
      INCLUDE 'chinert.h'
*
*
      IFOUND = 0
      DO 1 K = 1,NDORM
          IF (NMDORM(K).EQ.NAMEI) THEN
              B1 = BDORM(1,K)
              B2 = BDORM(2,K)
              IFOUND = K
          END IF
    1 CONTINUE
      IF (IFOUND.EQ.0) RETURN
      DO 2 K = IFOUND,NDORM-1
          NMDORM(K) = NMDORM(K+1)
          ADORM(K) = ADORM(K+1)
          BDORM(1,K) = BDORM(1,K+1)
          BDORM(2,K) = BDORM(2,K+1)
    2 CONTINUE
      NDORM = NDORM - 1
*
      RETURN
*
      END
