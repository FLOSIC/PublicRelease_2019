C UTEP Electronic Structure Lab (2019)
C
C    ------------------------------------------------------------------
C
C     **************************
C     LINE SEARCH ROUTINE MCSRCH
C     **************************
C
      SUBROUTINE MCSRCH(N,X,F,G,S,STP,FLTOL,XTOL,MAXFEV,INFO,NFEV,WA)
      INTEGER N,MAXFEV,INFO,NFEV
      DOUBLE PRECISION F,STP,FLTOL,GLTOL,XTOL,STPMIN,STPMAX
      DOUBLE PRECISION X(N),G(N),S(N),WA(N)
      COMMON /LB3/MP,LP,GLTOL,STPMIN,STPMAX
      SAVE
C
C                     SUBROUTINE MCSRCH
C                
C     A slight modification of the subroutine CSRCH of More' and Thuente.
C     The changes are to allow reverse communication, and do not affect
C     the performance of the routine. 
C
C     THE PURPOSE OF MCSRCH IS TO FIND A STEP WHICH SATISFIES
C     A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION.
C
C     AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF
C     UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF
C     UNCERTAINTY IS INITIALLY CHOSEN SO THAT IT CONTAINS A
C     MINIMIZER OF THE MODIFIED FUNCTION
C
C          F(X+STP*S) - F(X) - FLTOL*STP*(GRADF(X)'S).
C
C     IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION
C     HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE,
C     THEN THE INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT
C     CONTAINS A MINIMIZER OF F(X+STP*S).
C
C     THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES
C     THE SUFFICIENT DECREASE CONDITION
C
C           F(X+STP*S) .LE. F(X) + FLTOL*STP*(GRADF(X)'S),
C
C     AND THE CURVATURE CONDITION
C
C           ABS(GRADF(X+STP*S)'S)) .LE. GLTOL*ABS(GRADF(X)'S).
C
C     IF FLTOL IS LESS THAN GLTOL AND IF, FOR EXAMPLE, THE FUNCTION
C     IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES
C     BOTH CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH
C     CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING
C     ERRORS PREVENT FURTHER PROGRESS. IN THIS CASE STP ONLY
C     SATISFIES THE SUFFICIENT DECREASE CONDITION.
C
C     THE SUBROUTINE STATEMENT IS
C
C        SUBROUTINE MCSRCH(N,X,F,G,S,STP,FLTOL,XTOL, MAXFEV,INFO,NFEV,WA)
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF VARIABLES.
C
C       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
C         BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS
C         X + STP*S.
C
C       F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F
C         AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S.
C
C       G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
C         GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT
C         OF F AT X + STP*S.
C
C       S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE
C         SEARCH DIRECTION.
C
C       STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN
C         INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT
C         STP CONTAINS THE FINAL ESTIMATE.
C
C       FLTOL AND GLTOL ARE NONNEGATIVE INPUT VARIABLES. (In this reverse
C         communication implementation GLTOL is defined in a COMMON
C         statement.) TERMINATION OCCURS WHEN THE SUFFICIENT DECREASE
C         CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
C         SATISFIED.
C
C       XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS
C         WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
C         IS AT MOST XTOL.
C
C       STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH
C         SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP. (In this reverse
C         communication implementatin they are defined in a COMMON
C         statement).
C
C       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
C         OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST
C         MAXFEV BY THE END OF AN ITERATION.
C
C       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
C
C         INFO = 0  IMPROPER INPUT PARAMETERS.
C
C         INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT.
C
C         INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
C                   DIRECTIONAL DERIVATIVE CONDITION HOLD.
C
C         INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
C                   IS AT MOST XTOL.
C
C         INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.
C
C         INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.
C
C         INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.
C
C         INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
C                   THERE MAY NOT BE A STEP WHICH SATISFIES THE
C                   SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
C                   TOLERANCES MAY BE TOO SMALL.
C
C       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF
C         CALLS TO FCN.
C
C       WA IS A WORK ARRAY OF LENGTH N.
C
C     SUBPROGRAMS CALLED
C
C       MCSTEP
C
C       FORTRAN-SUPPLIED...ABS,MAX,MIN
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
C     JORGE J. MORE', DAVID J. THUENTE
C
C     **********
      INTEGER INFOC,J
      LOGICAL BRACKT,STAGE1,EXIST
      DOUBLE PRECISION DG,DGM,DGINIT,DGTEST,DGX,DGXM,DGY,DGYM,
     &       FINIT,FTEST1,FM,FX,FXM,FY,FYM,P5,P66,STX,STY,
     &       STMIN,STMAX,WIDTH,WIDTH1,XTRAPF,ZERO
      DATA P5,P66,XTRAPF,ZERO /0.5D0,0.66D0,4.0D0,0.0D0/

*********************CHANGES MADE HERE*****************************
        INQUIRE(FILE="STEP.LBF",EXIST=EXIST)
        IF (EXIST) THEN
        OPEN(UNIT=8,FILE="STEP.LBF",STATUS="UNKNOWN")
        READ(8,*)BRACKT
        READ(8,*)STAGE1
        READ(8,*)FINIT
        READ(8,*)DGTEST
        READ(8,*)WIDTH
        READ(8,*)WIDTH1
        READ(8,*)DGINIT
        READ(8,*)INFOC
        READ(8,*)STX
        READ(8,*)FX
        READ(8,*)DGX
        READ(8,*)STY
        READ(8,*)FY
        READ(8,*)DGY
        READ(8,*)STMIN
        READ(8,*)STMAX
        CLOSE(UNIT=8)
        ENDIF
************************TILL HERE***********************************

      IF(INFO.EQ.-1) GO TO 45
      INFOC = 1
C
C     CHECK THE INPUT PARAMETERS FOR ERRORS.
C
      IF (N .LE. 0 .OR. STP .LE. ZERO .OR. FLTOL .LT. ZERO .OR.
     &    GLTOL .LT. ZERO .OR. XTOL .LT. ZERO .OR. STPMIN .LT. ZERO
     &    .OR. STPMAX .LT. STPMIN .OR. MAXFEV .LE. 0) RETURN
C
C     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
C     AND CHECK THAT S IS A DESCENT DIRECTION.
C
      DGINIT = ZERO
      DO 10 J = 1, N
         DGINIT = DGINIT + G(J)*S(J)
   10    CONTINUE
      IF (DGINIT .GE. ZERO) then
         write(LP,15)
   15    FORMAT(/'  THE SEARCH DIRECTION IS NOT A DESCENT DIRECTION')
         RETURN
         ENDIF
C
C     INITIALIZE LOCAL VARIABLES.
C
      BRACKT = .FALSE.
      STAGE1 = .TRUE.
      NFEV = 0
      FINIT = F
      DGTEST = FLTOL*DGINIT
      WIDTH = STPMAX - STPMIN
      WIDTH1 = WIDTH/P5

	DO 20 J = 1, N
         WA(J) = X(J)
   20    CONTINUE
C
C     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
C     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
C     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
C     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
C     THE INTERVAL OF UNCERTAINTY.
C     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
C     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
C
      STX = ZERO
      FX = FINIT
      DGX = DGINIT
      STY = ZERO
      FY = FINIT
      DGY = DGINIT
C
C     START OF ITERATION.
C
   30 CONTINUE
C
C        SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
C        TO THE PRESENT INTERVAL OF UNCERTAINTY.
C
         IF (BRACKT) THEN
            STMIN = MIN(STX,STY)
            STMAX = MAX(STX,STY)
         ELSE
            STMIN = STX
            STMAX = STP + XTRAPF*(STP - STX)
            END IF
C
C        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
C
         STP = MAX(STP,STPMIN)
         STP = MIN(STP,STPMAX)
C
C        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
C        STP BE THE LOWEST POINT OBTAINED SO FAR.
C
         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))
     &      .OR. NFEV .GE. MAXFEV-1 .OR. INFOC .EQ. 0
     &      .OR. (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX)) STP = STX

*********************CHANGES MADE HERE*****************************
        OPEN(UNIT=8,FILE="STEP.LBF",STATUS="UNKNOWN")
        WRITE(8,*)BRACKT," BRACKT"
        WRITE(8,*)STAGE1," STAGE1"
        WRITE(8,*)FINIT," FINIT"
        WRITE(8,*)DGTEST," DGTEST"
        WRITE(8,*)WIDTH," WIDTH"
        WRITE(8,*)WIDTH1," WIDTH1"
        WRITE(8,*)DGINIT," DGINIT"
        WRITE(8,*)INFOC, " INFOC"
	WRITE(8,*)STX, " STX"
	WRITE(8,*)FX, " FX"
	WRITE(8,*)DGX," DGX"
	WRITE(8,*)STY," STY"
	WRITE(8,*)FY, " FY"
	WRITE(8,*)DGY," DGY"
	WRITE(8,*)STMIN, " STMIN"
	WRITE(8,*)STMAX, " STMAX"
        CLOSE(UNIT=8)
************************TILL HERE***********************************

C
C        EVALUATE THE FUNCTION AND GRADIENT AT STP
C        AND COMPUTE THE DIRECTIONAL DERIVATIVE.
C        We return to main program to obtain F and G.
C
         DO 40 J = 1, N
            X(J) = WA(J) + STP*S(J)
   40       CONTINUE
         INFO=-1
         RETURN
C
   45    INFO=0
         NFEV = NFEV + 1
         DG = ZERO
         DO 50 J = 1, N
            DG = DG + G(J)*S(J)
   50       CONTINUE
         FTEST1 = FINIT + STP*DGTEST
C
C        TEST FOR CONVERGENCE.
C
         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))
     &      .OR. INFOC .EQ. 0) INFO = 6
         IF (STP .EQ. STPMAX .AND.
     &       F .LE. FTEST1 .AND. DG .LE. DGTEST) INFO = 5
         IF (STP .EQ. STPMIN .AND.
     &       (F .GT. FTEST1 .OR. DG .GE. DGTEST)) INFO = 4
         IF (NFEV .GE. MAXFEV) INFO = 3
         IF (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX) INFO = 2
         IF (F .LE. FTEST1 .AND. ABS(DG) .LE. GLTOL*(-DGINIT)) INFO = 1
C
C        CHECK FOR TERMINATION.
C
         IF (INFO .NE. 0) RETURN
C
C        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
C        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
C
         IF (STAGE1 .AND. F .LE. FTEST1 .AND.
     &       DG .GE. MIN(FLTOL,GLTOL)*DGINIT) STAGE1 = .FALSE.
C
C        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
C        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
C        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
C        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
C        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
C
         IF (STAGE1 .AND. F .LE. FX .AND. F .GT. FTEST1) THEN
C
C           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
C
            FM = F - STP*DGTEST
            FXM = FX - STX*DGTEST
            FYM = FY - STY*DGTEST

           DGM = DG - DGTEST
            DGXM = DGX - DGTEST
            DGYM = DGY - DGTEST
C
C           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
C           AND TO COMPUTE THE NEW STEP.
C
            CALL MCSTEP(STX,FXM,DGXM,STY,FYM,DGYM,STP,FM,DGM,
     &                 BRACKT,STMIN,STMAX,INFOC)
C
C           RESET THE FUNCTION AND GRADIENT VALUES FOR F.
C
            FX = FXM + STX*DGTEST
            FY = FYM + STY*DGTEST
            DGX = DGXM + DGTEST
            DGY = DGYM + DGTEST
         ELSE
C
C           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
C           AND TO COMPUTE THE NEW STEP.
C
            CALL MCSTEP(STX,FX,DGX,STY,FY,DGY,STP,F,DG,
     &                 BRACKT,STMIN,STMAX,INFOC)
            END IF
C
C        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
C        INTERVAL OF UNCERTAINTY.
C
         IF (BRACKT) THEN
            IF (ABS(STY-STX) .GE. P66*WIDTH1)
     &         STP = STX + P5*(STY - STX)
            WIDTH1 = WIDTH
            WIDTH = ABS(STY-STX)
            END IF
C
C        END OF ITERATION.
C
         GO TO 30
C
C     LAST LINE OF SUBROUTINE MCSRCH.
C
      END
