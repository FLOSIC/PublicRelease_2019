C UTEP Electronic Structure Lab (2019)
C
      REAL*8 FUNCTION NUMRND(VALUE,ROUND)
C
C     ROUND REAL NUMBER TO REQUESTED LIMIT
C
C     BY ULISES REVELES, JUNE 2013.
C
C     ******************************************************************
C
C     List of variables:
C
C     VALUE: Value to round.
C     ROUND: Rounding limit.
C
C     ------------------------------------------------------------------
C
      IMPLICIT NONE
C
      REAL*8 VALUE,ROUND
C
      REAL*8 NK
      INTEGER VK
C
C     ------------------------------------------------------------------
C
C     --- GET INTEGER VALUE AND DIFFERENCE ---
C
      VK = INT(VALUE)
      NK = VALUE - VK
C
C     --- NOW ROUND ---
C
      IF (VALUE.GT.0.0) THEN
        NK = AINT(NK*ROUND + 0.5)/ROUND
      ELSE
        NK = AINT(NK*ROUND - 0.5)/ROUND
      END IF
      NUMRND = VK + NK
C
C     ------------------------------------------------------------------
C
      END
