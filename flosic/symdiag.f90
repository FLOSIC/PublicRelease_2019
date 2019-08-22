! UTEP Electronic Structure Lab (2019)

      SUBROUTINE SYMDIAG(MATRIX,DMAT,N,EIGVAL)
!
!     Purpose: DIAMAT DIAgonalizes the real symmetric MATrix 
!              using LAPACK.
!
!     BY LUIS BASURTO, SEPTEMBER 2013
!
!     ------------------------------------------------------------------
!
!     INPUT:
!
!     DMAT  : Dimension of matrix.
!     MATRIX: Matrix to be diagonalized (will be overwritten).
!     N     : Number of rows or columns of matrix.
!
!     OUTPUT:
!
!     EIGVAL: Eigenvalues.
!     MATRIX: Eigenvectors (transformation matrix).
!
!     LOCAL VARIABLES:
!
!     IERR   : INTEGER ERRO CODE.
!     LWORK  : INTEGER WORK ARRAY SIZE.
!     LIWORK : INTEGER IWORK ARRAY SIZE.
!     ALLOCATION : INTEGER ALLOCATION RESULT CODE.
!
!     LOCAL DYNAMICAL FIELDS:
!
!     LWORK: INTEGER SIZE OF WORK ARRAY.
!     WORK: REAL WORK VECTOR.
!
!     ------------------------------------------------------------------
!

      IMPLICIT NONE
!
      INTEGER,INTENT(IN)   :: DMAT,N
      REAL*8,INTENT(INOUT) :: MATRIX(DMAT,DMAT)
      REAL*8,INTENT(OUT) :: EIGVAL(DMAT)
!
      INTEGER :: LWORK,LIWORK,IERR,ALLOCATION
!
      REAL*8,ALLOCATABLE  :: WORK(:)
      INTEGER,ALLOCATABLE :: IWORK(:)
!
!     ------------------------------------------------------------------
!
!     --- ALLOCATE INITIAL WORK ARRAY ---
!
      ALLOCATE(WORK(1),IWORK(1),STAT=ALLOCATION)
      IF (ALLOCATION.GT.0) THEN
        WRITE(*,*) 'ERROR IN SYMDIAG: INITIAL ALLOCATION FAILED'
      END IF
!
!     -- REQUEST WORK ARRAY SIZE
!
      LWORK=-1
      LIWORK=-1
      CALL DSYEVD('V','U',N,MATRIX,DMAT,EIGVAL,WORK,LWORK,IWORK,LIWORK,IERR)
!
!     --- REALLOCATE TO NEW WORK ARRAY SIZE
!
      LWORK=WORK(1)
      LIWORK=IWORK(1)
      DEALLOCATE(WORK,IWORK)
      ALLOCATE(WORK(LWORK),IWORK(LIWORK),STAT=ALLOCATION)
      IF (ALLOCATION.GT.0) THEN
        WRITE(*,*) 'ERROR IN SYMDIAG: SECOND ALLOCATION FAILED'
      END IF
!
!     --- CALL EIGENSOLVER ---
!
      CALL DSYEVD('V','U',N,MATRIX,DMAT,EIGVAL,WORK,LWORK,IWORK,LIWORK,IERR)
      IF (IERR.NE.0) THEN
        WRITE(6,*) 'ERROR IN SYMDIAG, ERROR:',IERR
      END IF
!
!     --- DEALLOCATE LOCAL FIELDS ---
!
      DEALLOCATE(WORK,IWORK,STAT=ALLOCATION)
      IF (ALLOCATION.GT.0) THEN
        WRITE(*,*) 'ERROR IN DIAMAT: DEALLOCATION FAILED'
      END IF
!
!     --- ZERO OUT UNUSED MATRIX SPACE ---
!
      IF (N.LT.DMAT) THEN
        MATRIX(N+1:DMAT,1:DMAT) = 0.0
        MATRIX(1:DMAT,N+1:DMAT) = 0.0
      END IF
!
!     ------------------------------------------------------------------
!
      END
