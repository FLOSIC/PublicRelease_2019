C UTEP Electronic Structure Lab (2019)
C
C ****************************************************************
C
       SUBROUTINE ORBSIC(LSIC)
C       INCLUDE 'PARAMS'
C       INCLUDE 'commons.inc'
C
C
       use debug1
       use common2,only : N_CON, LSYMMAX, ISPN, NSPN, WFFILE
       use common5,only : PSI_COEF, OCCUPANCY, N_OCC, PSI, NWF, NWFS,
     & EVLOCC
       use common7,only : T1UNRV, T2UNRV
       use for_diag1,only : HAM, OVER, filo, SC1, eval
       use hstor1,only : HSTOR
       use common8,only : REP, N_REP, NDMREP, U_MAT, N_SALC, INDBEG
     &                   ,NS_TOT
!SIC module
       use LOCORB,only : TMAT,MORB,ZSIC,IRBSIC
       use MOCORB,only : SLAT,NFRM,ZTZL,JJJJJJ
       use ORBENG,only : EVALOCC
       INCLUDE  'PARAMA2'
       INTEGER,PARAMETER :: MAX_TOT=NDH*MAX_REP
       LOGICAL FAILED
       LOGICAL EXIST, LSIC
       CHARACTER*12 EVALSTR
!      REAL*8, DIMENSION (NDH):: EVALSAV(NDH), OCCTMP(MAX_TOT*MXSPN)
!      REAL*8, DIMENSION (NDH):: EVL(NDH)
       REAL*8, allocatable :: EVALSAV(:),OCCTMP(:),EVL(:)
C      COMMON/ORBENG/EVALOCC(MAX_OCC)
C      COMMON/EVLTMP/EVALSAV(MAX_TOT*MXSPN),OCCTMP(MAX_TOT*MXSPN)
C      INCLUDE 'commons.inc'
       INTEGER :: MD, KSPN, KORB, IFNCT, ISHELLA, I_SITE, N_NUC, ILOC,
     & I_LOCAL, I_SALC, IBASE, ICALL1, ICALL2, IL, IMS, IND_SALC,
     & INDEX, IOCC, IORB, IORBB, IORBE, IQ, IQ_BEG, IROW, ISHELL,
     & ISHELLV, ITOT, IWF, J_LOCAL, JORB, JWF, K_REP, K_ROW, KSALC, LI,
     & LSPN, LSPNB, LSPNE, MLOCAL, MO, MSTOT, MU, NDEG
       REAL*8 :: RVEC , RVECI, ADD1, ADD2, ADD3, ADDTMAT, PHILOC,
     & TIMER1, TIMER2, TMATMAX
       INTEGER :: MSPN, I, MREP, I_REP,I_OCC,NBASF, KB, IB, JB, IREP
       INTEGER :: ITER, NBAS, J, LBAS, KBAS, IBAS, JBAS,K
       REAL*8 :: ERR, DOT,ER1, ERROR, EVV

       ALLOCATE(EVALSAV(MAX_OCC))  !Bug: changed NDH->MAX_OCC
       ALLOCATE(EVL(NDH))
C
!      write(6,*)'from orbsic'
       DO 250 ISPN=1,NSPN
        DO I=1,NWFS(ISPN)
         J=I+(ISPN-1)*NWFS(1)
         EVL(I)=EVLOCC(J)
        END DO
   
!       WRITE(6,*)'EVAL BEFORE DIAG'
!       DO I=1,NS_TOT(1)
!        WRITE(6,*) EVL(I)
!       END DO

        CALL OVERLAP(1)
        KB=0
        DO IREP=1,N_REP
         NBAS=NS_TOT(IREP)
         DO IB=1,NBAS
          DO JB=1,NBAS
           HAM(IB,JB)=0.0D0
           OVER(IB,JB)=0.0D0
           FILO(IB,JB)=0.0D0
          END DO
         END DO


         DO 160 IB=1 ,NBAS
         DO 160 JB=IB,NBAS
          KB=KB+1
          OVER(JB,IB)=HSTOR(KB,1)
          OVER(IB,JB)=HSTOR(KB,1)
          HAM(JB,IB)=0.0D0
          HAM(IB,JB)=0.0D0
          FILO(IB,JB)=0.0D0
          FILO(JB,IB)=0.0D0
  160    CONTINUE
         DO KBAS=1,NBAS
          DO LBAS=1,NBAS
           DO IWF=1,NWFS(ISPN)
            HAM(KBAS,LBAS)=HAM(KBAS,LBAS)+
     &      EVL(IWF)*PSI_COEF(KBAS,IWF,IREP,ISPN)
     &              *PSI_COEF(LBAS,IWF,IREP,ISPN)
           END DO
          END DO
         END DO
         ERROR=0.0D0
         DO KBAS=1,NBAS
          DO LBAS=KBAS,NBAS
           ERROR=ERROR+ABS(HAM(KBAS,LBAS)-HAM(LBAS,KBAS))
          END DO
         END DO
         WRITE(6,*)'ERROR=',ERROR
        
         DO KBAS=1,NBAS
          DO JBAS=1,NBAS
           DO LBAS=1,NBAS
            FILO(KBAS,JBAS)=FILO(KBAS,JBAS)
     &                     +HAM(KBAS,LBAS)*OVER(LBAS,JBAS) 
           END DO
          END DO
         END DO
             
         DO IBAS=1,NBAS
          DO JBAS=1,NBAS
           HAM(IBAS,JBAS)=0.0D0
          END DO
         END DO
         DO IBAS=1,NBAS
          DO JBAS=1,NBAS
           DO KBAS=1,NBAS
            HAM(IBAS,JBAS)=HAM(IBAS,JBAS)
     &                    +OVER(IBAS,KBAS)*FILO(KBAS,JBAS) 
           END DO
          END DO
         END DO
         
         CALL DIAGGE(NDH,NBAS,HAM,OVER,EVAL,SC1,1)
!        WRITE(6,*)'AFTER DIAG'
C        DO I=1,NBAS
C         WRITE(6,*) EVAL(I)
          EVV=EVAL(NWFS(ISPN))
C        END DO
         DO I=1,NWFS(ISPN)
          J=I+(ISPN-1)*NS_TOT(1)
          EVALSAV(J)=EVAL(I)
         END DO
         DO I=NWFS(ISPN)+1,NS_TOT(IREP)
          J=I+(ISPN-1)*NS_TOT(1)
          EVALSAV(J)=EVV+0.1*I
         END DO
    
         DO I=1,NS_TOT(IREP)
!         WRITE(6,*) EVALSAV(I)
          DO IBAS=1,NS_TOT(IREP)
           PSI_COEF(IBAS,I,IREP,ISPN)=HAM(IBAS,I)
          END DO
         END DO
        END DO
  250  CONTINUE

       DEALLOCATE(EVL)

       ALLOCATE(OCCTMP(MAX_TOT*MXSPN))

       DO I=1,MAX_TOT*MXSPN
        OCCTMP(I)=0.0D0
       END DO
       J=0
       DO ISPN=1,NSPN
        K=(ISPN-1)*NWFS(1)
        DO IREP=1,N_REP
         DO IBAS=1,N_OCC(IREP,ISPN)
          J=J+1
          K=K+1
          OCCTMP(J)=OCCUPANCY(K)
         END DO
         DO IBAS=N_OCC(IREP,ISPN)+1,NS_TOT(IREP)
          J=J+1
          OCCTMP(J)=0.0D0
         END DO
        END DO
       END DO
       WRITE(6,*) 'ISPN, IREP, J, OCCUPANCY(J)'
       J=0
       NWF=0
       DO ISPN=1,NSPN
        NWFS(ISPN)=NS_TOT(1)
        NWF=NWF+NWFS(ISPN)
        DO IREP=1,N_REP
         N_OCC(IREP,ISPN)=NWFS(ISPN)
         DO IBAS=1,NS_TOT(IREP)
          J=J+1
          OCCUPANCY(J)=OCCTMP(J)
          EVALOCC(J)=EVALSAV(J)
          WRITE(6,*) ISPN, IREP, J, OCCUPANCY(J),EVALOCC(J)
 600      FORMAT(I3,I3,I5,F8.4,F12.8)
         END DO
        END DO
       END DO

       DEALLOCATE(OCCTMP)
       DEALLOCATE(EVALSAV)

       WRITE(6,*) 'ORBSIC:', NWF, NWFS(1),NWFS(2)

         


C        ERR=0.0D0
C        DO 175 IWF=1,N_OCC(IREP,ISPN)
C        DO 175 JWF=1,IWF
C         DOT=0.0D0
C         IF (JWF.EQ.IWF) THEN
C          DO IB=1,NBAS
C           DO JB=1,NBAS
C            DOT=DOT+OVER(JB,IB)*PSI_COEF(JB,JWF,IREP,ISPN)
C     &                         *PSI_COEF(IB,IWF,IREP,ISPN)
C           END DO
C          END DO
C          ERR=ERR+ABS(DOT-1.0D0)
C          DO IB=1,NBAS
C           PSI_COEF(IB,IWF,IREP,ISPN)=
C     &     PSI_COEF(IB,IWF,IREP,ISPN)/SQRT(DOT)
C          END DO
C          DO IB=1,NBAS
C           HAM(IB,IWF)=0.0D0
C           DO JB=1,NBAS
C            HAM(IB,IWF)=HAM(IB,IWF)
C     &                 +PSI_COEF(JB,IWF,IREP,ISPN)*OVER(JB,IB)
C           END DO
C          END DO
C         ELSE
C          DOT=0.0D0
C          DO IB=1,NBAS
C           DOT=DOT+PSI_COEF(IB,IWF,IREP,ISPN)*HAM(IB,JWF)
C          END DO
C          ERR=ERR+ABS(DOT      )
C          DO IB=1,NBAS
C           PSI_COEF(IB,IWF,IREP,ISPN)=PSI_COEF(IB,IWF,IREP,ISPN)
C     &    -PSI_COEF(IB,JWF,IREP,ISPN)*DOT
C          END DO
C         END IF
C  175   CONTINUE
C        ER1=ERR
CC
CC TEST
C
C        IF(DEBUG) THEN
C        ERR=0.0D0
C        DO 185 IWF=1,N_OCC(IREP,ISPN)
C        DO 185 JWF=1,IWF
C         DOT=0.0D0
C         DO IB=1,NBAS
C          DO JB=1,NBAS
C           DOT=DOT+OVER(JB,IB)*PSI_COEF(JB,JWF,IREP,ISPN)
C     &                        *PSI_COEF(IB,IWF,IREP,ISPN)
C          END DO
C         END DO
C         IF (IWF.EQ.JWF) THEN
C          ERR=ERR+ABS(DOT-1.0D0)
C         ELSE
C          ERR=ERR+ABS(DOT)
C         END IF
C  185   CONTINUE
C        PRINT *,'REORTHOGONALIZATION ERROR:',ERR,ER1
C        ENDIF
C  240  CONTINUE
C
C REMOVE OLD EVALUE FILES AND CREATE DUMMY FILE FOR FIRST ITERATION
C
       ITER=0
  300  ITER=ITER+1
       WRITE(EVALSTR,'(A,I3.3)')'EVAL',ITER
       INQUIRE(FILE=EVALSTR,EXIST=EXIST)
       IF (EXIST) THEN
        OPEN(98,FILE=EVALSTR,FORM='FORMATTED',STATUS='OLD')
        CLOSE(98,STATUS='DELETE')
        GOTO 300
       END IF
       CONTINUE
       OPEN(98,FILE='EVAL001',FORM='FORMATTED',STATUS='UNKNOWN')
       REWIND(98)
       WRITE(98,*) 'NO EIGENVALUES FOR READ WAVEFUNCTIONS'
       CLOSE(98)
       RETURN
C
C FAILURE
C
  900  PRINT '(A)','UNABLE TO READ OR PROCESS OLD WAVEFUNCTIONS'
       PRINT '(A)','PROCEEDING WITH DEFAULT STARTING POINT'
       FAILED= .TRUE.
       RETURN
       END
C
