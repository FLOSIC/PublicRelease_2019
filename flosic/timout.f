C UTEP Electronic Structure Lab (2019)
C
C
C *********************************************************************
C
        SUBROUTINE TIMOUT(STRING,TIME)
         use mpidat1,only  : IRANK
         CHARACTER(LEN=*) :: STRING
         REAL*8 :: TIME
         WRITE (6+IRANK, 1000)STRING,TIME
 1000    FORMAT('TIME FOR ',A,' ',F12.3)
         RETURN
        END
