C UTEP Electronic Structure Lab (2019)
C
      SUBROUTINE RESTCST
C
C     RESTORE INTERNAL CONSTRAINS.
C
C     ------------------------------------------------------------------
C
C     GLOBAL DIMENSIONS:
C
C     NPRI: Dimension of primitive internal coordinates.
C
C     ------------------------------------------------------------------
C
      use common2,only : CPRIMS,QPRIMS,NPRI
C
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:01 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
C
      INTEGER IPRI
C
C     ------------------------------------------------------------------
C
C     --- RESTORE INTERNAL CONSTRAINS ---
C
      DO IPRI=1,NPRI
        IF ((CPRIMS(IPRI,5).EQ.-99999).AND.
     $      (CPRIMS(IPRI,6).EQ.-99999)) THEN
          QPRIMS(IPRI,1) = QPRIMS(IPRI,2)
        END IF
      END DO
C
C     ------------------------------------------------------------------
C
      END
