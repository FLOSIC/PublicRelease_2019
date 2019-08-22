C UTEP Electronic Structure Lab (2019)
C
      SUBROUTINE WFWIND(EMIN,EMAX,SPN1,SPN2,IFAIL)
C
C     WFWIND VERSION DIRK POREZAG NOVEMBER 1996
C
C     ------------------------------------------------------------------
C
       use for_diag1
       use hstor1,only : hstor
       use debug1
       use common2,only : ISPN, NSPN
       use common5,only : PSI_COEF, OCCUPANCY, N_OCC, PSI,
     &   NWF, NWFS, EVLOCC
       use common8,only : REP, N_REP, NDMREP, IGEN, NS_TOT
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:07 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IFAIL, I, IB, IBAS, IDEG, IEIG, IERR, IND1, IND2,
     & IOFS, IREP, ISPBEG, ISPEND, JBAS, JSPN, KBAS, NBAS, NEIG, NSYMM,
     & NUSYM
       REAL*8 :: EMIN , EMAX, TIME1, TIME2
       SAVE
       LOGICAL SPN1,SPN2,DOIT
C
C     ------------------------------------------------------------------
C
C     --- CHECKING AND SETTING UP SOME STUFF ---
C
       IFAIL=0
       CALL GTTIME(TIME1)
       IF (DEBUG) PRINT '(A,2(1X,F15.5))',
     &            ' IN WFWIND: EMIN, EMAX:',EMIN,EMAX
       IF(N_REP.GT.MAX_REP)THEN
         PRINT *,'WFWIND: MAX_REP MUST BE AT LEAST: ',N_REP
         CALL STOPIT
       END IF
C
C     --- LOOP OVER SPIN ---
C
      NSYMM=0
      NUSYM=0
      ISPBEG=NSPN
      ISPEND=1
      IF(SPN1)ISPBEG=1
      IF(SPN2)ISPEND=NSPN
C
      NWF=0
      NWFS   (1)=0
      NWFS(NSPN)=0
      DO JSPN=1,NSPN
        DO IREP=1,N_REP
          N_OCC(IREP,JSPN)=0
        END DO
      END DO
C
      DO 240 ISPN=ISPBEG,ISPEND
        DOIT=.FALSE.
        IF (ISPN.EQ.1.AND.SPN1) DOIT=.TRUE.
        IF (ISPN.EQ.2.AND.SPN2) DOIT=.TRUE.
C
        IF(DOIT)THEN
          IF (DEBUG) PRINT *,'WFWIND CALLS OVERLAP MODE: 1'
          CALL OVERLAP(1)
          IF (DEBUG) PRINT *,'WFWIND CALLS OVERLAP MODE: 2'
          CALL OVERLAP(2)
C
C     --- LOOP OVER REPRESENTATIONS ---
C     --- GET MATRIX ELEMENTS       ---
C
          KBAS=0
          DO 130 IREP=1,N_REP
            NBAS=NS_TOT(IREP)
            IF(NBAS.GT.NDH)THEN
              PRINT *,'WFWIND: NDH MUST BE AT LEAST: ',NBAS
              CALL STOPIT
            END IF
C
C     --- ALLOCATE LOCAL FIELDS ---
C
          ALLOCATE(AHAM(NBAS,NBAS),STAT=IERR)
          IF(IERR.NE.0)THEN
            WRITE(6,*)'wfwind:Error allocating Ham'
          ENDIF 
          ALLOCATE(AOVER(NBAS,NBAS),STAT=IERR)
          IF(IERR.NE.0)THEN
            WRITE(6,*)'wfwind:Error allocating Overlap'
          ENDIF 
          ALLOCATE(AEVAL(NBAS),STAT=IERR)
          IF(IERR.NE.0)THEN
            WRITE(6,*)'wfwind:Error allocating Eval'
          ENDIF 
          ALLOCATE(ASC1(NBAS),STAT=IERR)
          IF(IERR.NE.0)THEN
            WRITE(6,*)'wfwind:Error allocating Sc1'
          ENDIF 
C
          DO 80 IBAS=1,NBAS
            DO 70 JBAS=IBAS,NBAS
              KBAS=KBAS+1
              AOVER(JBAS,IBAS)=HSTOR(KBAS,1)
              AHAM (JBAS,IBAS)=HSTOR(KBAS,2)
   70       CONTINUE
   80     CONTINUE
C
CJK01/2001
c         CALL SCISSOR(IREP)
CJK01/2001       
C
C     --- GET EIGENVECTORS AND EIGENVALUES   ---
C     --- (WILL BE RETURNED IN HAM AND EVAL) ---
C
          IF(NBAS.GT.0)THEN
            CALL DIAGGE(NBAS,NBAS,AHAM,AOVER,AEVAL,ASC1,1)
          ENDIF
C
C     --- DETERMINE WHICH EIGENSTATES ARE NEEDED ---
C
          IND1=0
          IND2=0
          DO I=1,NBAS
            IF (IND1.EQ.0) THEN
              IF (AEVAL(I) .GE. EMIN) IND1=I
            END IF
            IF (IND2.EQ.0) THEN
              IF (AEVAL(I) .GT. EMAX) IND2=I
            END IF
          END DO
C
          IF (IND1.EQ.0) THEN
            IOFS=0
            NEIG=0
          ELSE IF (IND2.EQ.0) THEN
            IOFS=IND1-1
            NEIG=NBAS-IOFS
          ELSE
            IOFS=IND1-1
            NEIG=IND2-IND1
          END IF
C
          IF (NEIG.GT.MAX_VIRT_PER_SYM) THEN
            PRINT *,'WFWIND: MAX_VIRT_PER_SYM MUST BE AT LEAST: ',NEIG
            IFAIL=1
            RETURN
          END IF
C
C     --- STORE WAVEFUNCTIONS AND EIGENVALUES ---
C
          N_OCC(IREP,ISPN)=NEIG
          NWFS(ISPN)=NWFS(ISPN)+NEIG*NDMREP(IREP)
          NWF=NWF+NEIG*NDMREP(IREP)
          DO IEIG=1,NEIG
            NSYMM=NSYMM+1
            OCCUPANCY(NSYMM)=1.0D0
            DO IB=1,NBAS
              PSI_COEF(IB,IEIG,IREP,ISPN)=AHAM(IB,IEIG+IOFS)
            END DO
C
            DO IDEG=1,NDMREP(IREP)
              NUSYM=NUSYM+1
              IF (NUSYM.GT.MAX_OCC) THEN
                PRINT *,'WFWIND: MAX_OCC MUST BE AT LEAST ',NUSYM
                IFAIL=1
                RETURN
              END IF
              EVLOCC(NUSYM)=AEVAL(IEIG+IOFS)
            END DO
          END DO
C
C     --- DEALLOCATE LOCAL FIELDS ---
C
          DEALLOCATE(AHAM,STAT=IERR)
          IF(IERR.NE.0)THEN
            WRITE(6,*)'wfwind:Error deallocating Ham'
          END IF 
          DEALLOCATE(AOVER,STAT=IERR)
          IF(IERR.NE.0)THEN
            WRITE(6,*)'wfwind:Error deallocating Overlap'
          ENDIF 
          DEALLOCATE(AEVAL,STAT=IERR)
          IF(IERR.NE.0)THEN
            WRITE(6,*)'wfwind:Error deallocating Eval'
          ENDIF 
          DEALLOCATE(ASC1,STAT=IERR)
          IF(IERR.NE.0)THEN
            WRITE(6,*)'wfwind:Error deallocating Sc1'
          ENDIF 
C
  130     CONTINUE
C
        END IF
  240 CONTINUE
C
       CALL GTTIME(TIME2)
       CALL TIMOUT('WAVEFUNCTIONS IN ENERGY WINDOW:    ',TIME2-TIME1)
C
       RETURN
C
C     ------------------------------------------------------------------
C
       END
