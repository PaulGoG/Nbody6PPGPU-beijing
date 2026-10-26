      SUBROUTINE IMBHINIT
*
*      
*     Initialization of massive black hole
*
      INCLUDE 'common6.h'
      REAL*8 XIMBH(3),VIMBH(3)
      
*     Read massive black hole mass and output time interval
      WRITE(6,53)
   53 FORMAT(//,'**** imbhinit KZ(24).EQ.1: rewrite input!! ****',//)
      if(rank.eq.0) READ (5,*) BIMBH, XIMBH(1:3),VIMBH(1:3), DTBH

#ifdef PARALLEL
      CALL MPI_BCAST(BIMBH,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(XIMBH,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(VIMBH,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(DTBH,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif      

*     Reset mass of particle N, treat N as IMBH
      BODY(N) = BIMBH/ZMBAR
      KSTAR(N) = 14
      TEV(N) = MAX(2*TCRIT,1000000.0)
      DO K = 1,3
         X(K,N) = XIMBH(K)/RBAR
         XDOT(K,N) = VIMBH(K)/VSTAR
      END DO
*     Store the name of massive black hole (and the list of the flag route)
      NIMBH = N
      NIMBHS = 1
      NAMEBH(1) = N

*     Printing initial massive black hole data
      IF(rank.eq.0) write(6,1) NIMBH,BIMBH,XIMBH(1:3),
     &     VIMBH(1:3)
 1    format('MBH Initialization:  NAME',I12,' M[M*]',E15.5,
     &     ' X(1:3)[PC]',3E15.5,' V(1:3)[km/s]',3E15.5)
      
      RETURN
      
      END

***** Note:------------------------------------------------------------**

      SUBROUTINE IMBHFLAG
*
*
*     Flag route of the massive black holes (KZ(24) = 3): up to 16
*     bodies of dat.10, named by their index, keep their mass, position
*     and velocity and become stellar type 14. Input after &INSCALE
*     (and &INXTRNL0): NIMBHS DTBH on one line, the names on the next.
*
      INCLUDE 'common6.h'
*
      WRITE (6,53)
   53 FORMAT (//,'**** imbhflag KZ(24).EQ.3: NIMBHS DTBH / NAMES ****',
     &     //)
      if(rank.eq.0) READ (5,*) NIMBHS, DTBH
#ifdef PARALLEL
      CALL MPI_BCAST(NIMBHS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(DTBH,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif
      IF (NIMBHS.LT.1.OR.NIMBHS.GT.16) THEN
         if(rank.eq.0) WRITE (6,*) 'imbhflag: NIMBHS outside 1..16: ',
     &        NIMBHS
         call flush(6)
         call abort
      END IF
      if(rank.eq.0) READ (5,*) (NAMEBH(K),K=1,NIMBHS)
#ifdef PARALLEL
      CALL MPI_BCAST(NAMEBH,16,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
*
*     The bodies keep the dat.10 state; type, evolution time and the
*     total mass BIMBH (used by the Lagrangian radii) are set here.
      BIMBH = 0.0D0
      DO K = 1,NIMBHS
         I = NAMEBH(K)
         IF (I.LT.1.OR.I.GT.N) THEN
            if(rank.eq.0) WRITE (6,*) 'imbhflag: name outside 1..N: ',I
            call flush(6)
            call abort
         END IF
         DO L = 1,K-1
            IF (NAMEBH(L).EQ.I) THEN
               if(rank.eq.0) WRITE (6,*) 'imbhflag: repeated name ',I
               call flush(6)
               call abort
            END IF
         END DO
         KSTAR(I) = 14
         TEV(I) = MAX(2*TCRIT,1000000.0)
         BIMBH = BIMBH + BODY(I)*ZMBAR
         IF (rank.eq.0) WRITE (6,1) I,BODY(I)*ZMBAR,
     &        (X(KK,I)*RBAR,KK=1,3),(XDOT(KK,I)*VSTAR,KK=1,3)
      END DO
    1 FORMAT ('MBH flag:  NAME',I12,' M[M*]',E15.5,
     &     ' X(1:3)[PC]',3E15.5,' V(1:3)[km/s]',3E15.5)
      NIMBH = NAMEBH(1)
      call flush(6)
*
      RETURN
*
      END

***** Note:------------------------------------------------------------**

      LOGICAL FUNCTION IMBHNAME(NM)
*
*
*     True when NM is the name of a massive black hole of the list.
*
      INCLUDE 'common6.h'
*
      IMBHNAME = .FALSE.
      DO K = 1,NIMBHS
         IF (NAMEBH(K).EQ.NM) THEN
            IMBHNAME = .TRUE.
            RETURN
         END IF
      END DO
*
      RETURN
*
      END

***** Note:------------------------------------------------------------**

      SUBROUTINE IMBHRESTART
*
*
*     Massive black holes after a restart: every name of the list is
*     looked up; a missing one (escaped and removed, merged) is dropped
*     with a warning instead of ending the run. Dumps written before the
*     list existed carry NIMBHS = 0 and rebuild it from NIMBH.
*
      INCLUDE 'common6.h'
*
      IF (NIMBHS.EQ.0.AND.KZ(24).EQ.1) THEN
         NIMBHS = 1
         NAMEBH(1) = NIMBH
      END IF
      BIMBH = 0.0D0
      KEEP = 0
      DO K = 1,NIMBHS
         IBH = 0
         DO L = N,1,-1
            IF (NAME(L).EQ.NAMEBH(K)) THEN
               IBH = L
               GO TO 10
            END IF
         END DO
 10      IF (IBH.EQ.0) THEN
            if(rank.eq.0) WRITE (6,*) 'MBH not found after restart,',
     &           ' dropped from the list: NAME',NAMEBH(K)
         ELSE
            KEEP = KEEP + 1
            NAMEBH(KEEP) = NAMEBH(K)
            BIMBH = BIMBH + BODY(IBH)*ZMBAR
            IF(rank.eq.0) write(6,1) IBH,NAMEBH(KEEP),BODY(IBH)*ZMBAR,
     &           X(1:3,IBH)*RBAR,XDOT(1:3,IBH)*VSTAR
         END IF
      END DO
      NIMBHS = KEEP
      IF (NIMBHS.GT.0) NIMBH = NAMEBH(1)
      IF (NIMBHS.EQ.0.and.rank.eq.0) WRITE (6,*)
     &     'MBH: no massive black hole left after the restart'
 1    format('MBH: IBH',I12,' NAME',I12,' M[M*]',E15.5,
     &     ' X(1:3)[PC]',3E15.5,' V(1:3)[km/s]',3E15.5)
      call flush(6)
*
      RETURN
*
      END
      
