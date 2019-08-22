C UTEP Electronic Structure Lab (2019)
C
      SUBROUTINE BFGS(DE1,DE2,DR1,DDEG,DUPD,UPDATED)
C
C     BFGS HESSIAN UPDATE
C
C     REFERENCES: C.G. Broyden, J. Inst. Maths. Applns. 6, 76 (1970)
C                 R. Fletcher, Computer J. 13, 317 (1970)
C                 D. Goldfarb, Maths. Comp. 24, 23 (1970)
C                 D.F. Shanno, Maths. Comp. 24, 647 (1970)
C
C     BY ULISES REVELES, JULY 2013.
C
C     ------------------------------------------------------------------
C
C     LOCAL DIMENSIONS:
C
C     DDEG: Dimension of degrees of freedom.
C     DUPD: Dimension of steps in Hessian update.
C
C     LOCAL VARIABLES:
C
C     DE1: Energy gradients of optimization steps.
C     DE2: Hessian matrix -> updated Hessian matrix.
C     DR1: Coordinates of optimization steps.
C
C     LOCAL DYNAMICAL FIELDS:
C
C     DELTA: Geometry difference vector.
C     GAMMA: Gradient difference vector.
C     WORK : Work vector.
C
C     ------------------------------------------------------------------
C
      REAL*8 EPS,FCLEAN
      PARAMETER (EPS = 1.E-8, FCLEAN = 1.E7)
C
      INTEGER DDEG,DUPD
      REAL*8 DE1(DDEG,DUPD),DE2(DDEG,DDEG),DR1(DDEG,DUPD)
C
      CHARACTER*80 PRTSTR,STRCOMP
      LOGICAL UPDATED
C
      INTEGER I,J
      REAL*8 NUMRND
      REAL*8 FACTOR(2)
C
      INTEGER ALLOCATION
      REAL*8,ALLOCATABLE :: DELTA(:),GAMMA(:),WORK(:)
C
C     ------------------------------------------------------------------
C
C     --- ALLOCATE LOCAL FIELDS ---
C
      ALLOCATE(DELTA(DDEG),GAMMA(DDEG),WORK(DDEG),STAT=ALLOCATION)
      IF (ALLOCATION.GT.0) THEN
        WRITE(*,*)'ERROR IN BFGS: ALLOCATION FAILED'
      END IF
C
C     --- INITIALIZATION ---
C
      DELTA(:) = 0.0
      GAMMA(:) = 0.0
      WORK(:) = 0.0
C
C     --- CHECK IF PREVIOS CYCLE INFORMATION EXIST ---
C
      DELTA(:) = DR1(:,2)
      GAMMA(:) = DE1(:,2)
C
      FACTOR(1) = DOT_PRODUCT(DELTA,DELTA)
      FACTOR(2) = DOT_PRODUCT(GAMMA,GAMMA)
      IF ((FACTOR(1).LT.EPS).OR.(FACTOR(2).LT.EPS)) THEN
        WRITE(*,*)'HESSIAN PRESERVED AND NOT UPDATED'
        GO TO 1000
      END IF
C
C     -- NOW TRY TO UPDATE HESSIAN MATRIX ---
C
      UPDATED = .TRUE.
C
C     --- CALCULATE GRADIENT AND GEOMETRY DIFFERENCE VECTORS ---
C
      DELTA(:) = 0.0
      GAMMA(:) = 0.0
      DELTA(:) = DR1(:,1) - DR1(:,2)
      GAMMA(:) = DE1(:,1) - DE1(:,2)
C
C     --- CLEAN NUMERICAL NOISE IN STEPS AND GRADIENTS ---
C
      DO I=1,DDEG
        DELTA(I) = NUMRND(DELTA(I),FCLEAN)
        GAMMA(I) = NUMRND(GAMMA(I),FCLEAN)
      END DO
C
      FACTOR(1) = 1.0/DOT_PRODUCT(GAMMA,DELTA)
C
      WORK(1:DDEG) = MATMUL(DE2(1:DDEG,1:DDEG),DELTA(1:DDEG))
C
      FACTOR(2) = 1.0/DOT_PRODUCT(DELTA,WORK)
C
C     --- BFGS UPDATE PRESERVING POSITIVE DEFINITNESS ---
C
      WRITE (PRTSTR,5000) 1.0/FACTOR(1)
 5000 FORMAT ('BFGS DELTA*GAMMA: ',G10.3)
      WRITE(*,*)STRCOMP(PRTSTR)
      WRITE (PRTSTR,5010) 1.0/FACTOR(2)
 5010 FORMAT ('BFGS DELTA*DE2*DELTA: ',G10.3)
      WRITE(*,*)STRCOMP(PRTSTR)
C
      IF (1.0/FACTOR(1).LT.EPS) THEN
C       IF (REDOSTEP) OPTRST = .TRUE.
        WRITE(*,*)'HESSIAN PRESERVED AND NOT UPDATED'
        UPDATED = .FALSE.
        GO TO 1000
      END IF
C
C     --- UPDATE HESSIAN MATRIX ---
C
      DO I=1,DDEG
        DO J=1,DDEG
          DE2(I,J) = DE2(I,J) + FACTOR(1)*GAMMA(I)*GAMMA(J) -
     $               FACTOR(2)*WORK(I)*WORK(J)
        END DO
      END DO
C
 1000 CONTINUE
C
C     --- DEALLOCATE LOCAL VECTORS ---  
C
      DEALLOCATE(DELTA,GAMMA,WORK,STAT=ALLOCATION)
      IF (ALLOCATION.GT.0) THEN
        WRITE(*,*)'ERROR IN BFGS: ALLOCATION FAILED'
      END IF
C
C     ------------------------------------------------------------------
C
      END
