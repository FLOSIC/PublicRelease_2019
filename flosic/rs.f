C UTEP Electronic Structure Lab (2019)
C
      SUBROUTINE RS (NM,N,A,W,MATZ,Z,FV1,FV2,IERR)
C
      INTEGER N,NM,IERR,MATZ
      REAL*8 A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)
C
C     this subroutine calls the recommended sequence of
C     subroutines from the eigensystem subroutine package (eispack)
C     to find the eigenvalues and eigenvectors (if desired)
C     of a real symmetric matrix.
C
C     on input
C
C     nm  must be set to the row dimension of the two-dimensional
C     array parameters as declared in the calling program
C     dimension statement.
C
C     n  is the order of the matrix  a.
C
C     a  contains the real symmetric matrix.
C
C     matz  is an integer variable set equal to zero if
C     only eigenvalues are desired.  otherwise it is set to
C     any non-zero integer for both eigenvalues and eigenvectors.
C
C     on output
C
C     w  contains the eigenvalues in ascending order.
C
C     z  contains the eigenvectors if matz is not zero.
C
C     ierr  is an integer output variable set equal to an error
C     completion code described in the documentation for tqlrat
C     and tql2.  the normal completion code is zero.
C
C     fv1  and  fv2  are temporary storage arrays.
C
C     questions and comments should be directed to burton s. garbow,
C     mathematics and computer science div, argonne national laboratory
C
C     this version dated august 1983.
C
C     ------------------------------------------------------------------
C
      IF (N.LE.NM) GO TO 10
      IERR = 10*N
      GO TO 30
C
   10 CONTINUE
      IF (MATZ.NE.0) GO TO 20
C
C     .......... find eigenvalues only ..........
C
      CALL TRED1 (NM,N,A,W,FV1,FV2)
C
C     tqlrat encounters catastrophic underflow on the Vax
C     call tqlrat (n,w,fv2,ierr)
C
      CALL TQL1 (N,W,FV1,IERR)
      GO TO 30
C
C     .......... find both eigenvalues and eigenvectors ..........
C
   20 CONTINUE
      CALL TRED2 (NM,N,A,W,FV1,Z)
      CALL TQL2 (NM,N,W,FV1,Z,IERR)
   30 CONTINUE
      END
      SUBROUTINE TQL1 (N,D,E,IERR)
C
      INTEGER I,J,L,M,N,II,L1,L2,MML,IERR
      REAL*8 D(N),E(N)
      REAL*8 C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2,TST1,TST2,PYTHAG,RONE
C
C     this subroutine is a translation of the algol procedure tql1,
C     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
C     wilkinson.
C     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
C
C     this subroutine finds the eigenvalues of a symmetric
C     tridiagonal matrix by the ql method.
C
C     on input
C
C     n is the order of the matrix.
C
C     d contains the diagonal elements of the input matrix.
C
C     e contains the subdiagonal elements of the input matrix
C     in its last n-1 positions.  e(1) is arbitrary.
C
C     on output
C
C     d contains the eigenvalues in ascending order.  if an
C     error exit is made, the eigenvalues are correct and
C     ordered for indices 1,2,...ierr-1, but may not be
C     the smallest eigenvalues.
C
C     e has been destroyed.
C
C     ierr is set to
C     zero       for normal return,
C     j          if the j-th eigenvalue has not been
C     determined after 30 iterations.
C
C     calls pythag for  dsqrt(a*a + b*b) .
C
C     questions and comments should be directed to burton s. garbow,
C     mathematics and computer science div, argonne national laboratory
C
C     this version dated august 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N.EQ.1) GO TO 140
C
      DO 10 I = 2,N
        E(I - 1) = E(I)
   10 CONTINUE
C
      F = 0.0
      TST1 = 0.0
      E(N) = 0.0
C
      DO 120 L = 1,N
        J = 0
        H = ABS(D(L)) + ABS(E(L))
        IF (TST1.LT.H) TST1 = H
C
C       .......... look for small sub-diagonal element ..........
C
        DO 20 M = L,N
          TST2 = TST1 + ABS(E(M))
          IF (TST2.EQ.TST1) GO TO 30
C
C         .......... e(n) is always zero, so there is no exit
C         through the bottom of the loop ..........
C
   20   CONTINUE
C
   30   CONTINUE
        IF (M.EQ.L) GO TO 80
   40   CONTINUE
        IF (J.EQ.30) GO TO 130
        J = J + 1
C
C       .......... form shift ..........
C
        L1 = L + 1
        L2 = L1 + 1
        G = D(L)
        P = (D(L1) - G)/ (2.0*E(L))
        RONE = 1.0
        R = PYTHAG(P,RONE)
        D(L) = E(L)/ (P + SIGN(R,P))
        D(L1) = E(L)* (P + SIGN(R,P))
        DL1 = D(L1)
        H = G - D(L)
        IF (L2.GT.N) GO TO 60
C
        DO 50 I = L2,N
          D(I) = D(I) - H
   50   CONTINUE
C
   60   CONTINUE
        F = F + H
C
C       .......... ql transformation ..........
C
        P = D(M)
        C = 1.0
        C2 = C
        EL1 = E(L1)
        S = 0.0
        MML = M - L
C
C       .......... for i=m-1 step -1 until l do -- ..........
C
        DO 70 II = 1,MML
          C3 = C2
          C2 = C
          S2 = S
          I = M - II
          G = C*E(I)
          H = C*P
          R = PYTHAG(P,E(I))
          E(I + 1) = S*R
          S = E(I)/R
          C = P/R
          P = C*D(I) - S*G
          D(I + 1) = H + S* (C*G + S*D(I))
   70   CONTINUE
C
        P = -S*S2*C3*EL1*E(L)/DL1
        E(L) = S*P
        D(L) = C*P
        TST2 = TST1 + ABS(E(L))
        IF (TST2.GT.TST1) GO TO 40
   80   CONTINUE
        P = D(L) + F
C
C       .......... order eigenvalues ..........
C
        IF (L.EQ.1) GO TO 100
C
C       .......... for i=l step -1 until 2 do -- ..........
C
        DO 90 II = 2,L
          I = L + 2 - II
          IF (P.GE.D(I - 1)) GO TO 110
          D(I) = D(I - 1)
   90   CONTINUE
C
  100   CONTINUE
        I = 1
  110   CONTINUE
        D(I) = P
  120 CONTINUE
C
      GO TO 140
C
C     .......... set error -- no convergence to an
C     eigenvalue after 30 iterations ..........
C
  130 CONTINUE
      IERR = L
  140 CONTINUE
      END
      SUBROUTINE TQL2 (NM,N,D,E,Z,IERR)
C
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      REAL*8 D(N),E(N),Z(NM,N)
      REAL*8 C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2,TST1,TST2,PYTHAG,RONE
C
C     this subroutine is a translation of the algol procedure tql2,
C     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
C     wilkinson.
C     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
C
C     this subroutine finds the eigenvalues and eigenvectors
C     of a symmetric tridiagonal matrix by the ql method.
C     the eigenvectors of a full symmetric matrix can also
C     be found if  tred2  has been used to reduce this
C     full matrix to tridiagonal form.
C
C     on input
C
C     nm must be set to the row dimension of two-dimensional
C     array parameters as declared in the calling program
C     dimension statement.
C
C     n is the order of the matrix.
C
C     d contains the diagonal elements of the input matrix.
C
C     e contains the subdiagonal elements of the input matrix
C     in its last n-1 positions.  e(1) is arbitrary.
C
C     z contains the transformation matrix produced in the
C     reduction by  tred2, if performed.  if the eigenvectors
C     of the tridiagonal matrix are desired, z must contain
C     the identity matrix.
C
C     on output
C
C     d contains the eigenvalues in ascending order.  if an
C     error exit is made, the eigenvalues are correct but
C     unordered for indices 1,2,...,ierr-1.
C
C     e has been destroyed.
C
C     z contains orthonormal eigenvectors of the symmetric
C     tridiagonal (or full) matrix.  if an error exit is made,
C     z contains the eigenvectors associated with the stored
C     eigenvalues.
C
C     ierr is set to
C     zero       for normal return,
C     j          if the j-th eigenvalue has not been
C     determined after 30 iterations.
C
C     calls pythag for  dsqrt(a*a + b*b) .
C
C     questions and comments should be directed to burton s. garbow,
C     mathematics and computer science div, argonne national laboratory
C
C     this version dated august 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N.EQ.1) GO TO 150
C
      DO 10 I = 2,N
        E(I - 1) = E(I)
   10 CONTINUE
C
      F = 0.0
      TST1 = 0.0
      E(N) = 0.0
C
      DO 100 L = 1,N
        J = 0
        H = ABS(D(L)) + ABS(E(L))
        IF (TST1.LT.H) TST1 = H
C
C       .......... look for small sub-diagonal element ..........
C
        DO 20 M = L,N
          TST2 = TST1 + ABS(E(M))
          IF (TST2.EQ.TST1) GO TO 30
C
C         .......... e(n) is always zero, so there is no exit
C         through the bottom of the loop ..........
C
   20   CONTINUE
C
   30   CONTINUE
        IF (M.EQ.L) GO TO 90
   40   CONTINUE
        IF (J.EQ.30) GO TO 140
        J = J + 1
C
C       .......... form shift ..........
C
        L1 = L + 1
        L2 = L1 + 1
        G = D(L)
        P = (D(L1) - G)/ (2.0*E(L))
        RONE = 1.0
        R = PYTHAG(P,RONE)
        D(L) = E(L)/ (P + SIGN(R,P))
        D(L1) = E(L)* (P + SIGN(R,P))
        DL1 = D(L1)
        H = G - D(L)
        IF (L2.GT.N) GO TO 60
C
        DO 50 I = L2,N
          D(I) = D(I) - H
   50   CONTINUE
C
   60   CONTINUE
        F = F + H
C
C       .......... ql transformation ..........
C
        P = D(M)
        C = 1.0
        C2 = C
        EL1 = E(L1)
        S = 0.0
        MML = M - L
C
C       .......... for i=m-1 step -1 until l do -- ..........
C
        DO 80 II = 1,MML
          C3 = C2
          C2 = C
          S2 = S
          I = M - II
          G = C*E(I)
          H = C*P
          R = PYTHAG(P,E(I))
          E(I + 1) = S*R
          S = E(I)/R
          C = P/R
          P = C*D(I) - S*G
          D(I + 1) = H + S* (C*G + S*D(I))
C
C         .......... form vector ..........
C
          DO 70 K = 1,N
            H = Z(K,I + 1)
            Z(K,I + 1) = S*Z(K,I) + C*H
            Z(K,I) = C*Z(K,I) - S*H
   70     CONTINUE
C
   80   CONTINUE
C
        P = -S*S2*C3*EL1*E(L)/DL1
        E(L) = S*P
        D(L) = C*P
        TST2 = TST1 + ABS(E(L))
        IF (TST2.GT.TST1) GO TO 40
   90   CONTINUE
        D(L) = D(L) + F
  100 CONTINUE
C
C     .......... order eigenvalues and eigenvectors ..........
C
      DO 130 II = 2,N
        I = II - 1
        K = I
        P = D(I)
C
        DO 110 J = II,N
          IF (D(J).GE.P) GO TO 110
          K = J
          P = D(J)
  110   CONTINUE
C
        IF (K.EQ.I) GO TO 130
        D(K) = D(I)
        D(I) = P
C
        DO 120 J = 1,N
          P = Z(J,I)
          Z(J,I) = Z(J,K)
          Z(J,K) = P
  120   CONTINUE
C
  130 CONTINUE
C
      GO TO 150
C
C     .......... set error -- no convergence to an
C     eigenvalue after 30 iterations ..........
C
  140 CONTINUE
      IERR = L
  150 CONTINUE
      END
      SUBROUTINE TRED1 (NM,N,A,D,E,E2)
C
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL*8 A(NM,N),D(N),E(N),E2(N)
      REAL*8 F,G,H,SCALE
C
C     this subroutine is a translation of the algol procedure tred1,
C     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
C     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
C
C     this subroutine reduces a real symmetric matrix
C     to a symmetric tridiagonal matrix using
C     orthogonal similarity transformations.
C
C     on input
C
C     nm must be set to the row dimension of two-dimensional
C     array parameters as declared in the calling program
C     dimension statement.
C
C     n is the order of the matrix.
C
C     a contains the real symmetric input matrix.  only the
C     lower triangle of the matrix need be supplied.
C
C     on output
C
C     a contains information about the orthogonal trans-
C     formations used in the reduction in its strict lower
C     triangle.  the full upper triangle of a is unaltered.
C
C     d contains the diagonal elements of the tridiagonal matrix.
C
C     e contains the subdiagonal elements of the tridiagonal
C     matrix in its last n-1 positions.  e(1) is set to zero.
C
C     e2 contains the squares of the corresponding elements of e.
C     e2 may coincide with e if the squares are not needed.
C
C     questions and comments should be directed to burton s. garbow,
C     mathematics and computer science div, argonne national laboratory
C
C     this version dated august 1983.
C
C     ------------------------------------------------------------------
C
      DO 10 I = 1,N
        D(I) = A(N,I)
        A(N,I) = A(I,I)
   10 CONTINUE
C
C     .......... for i=n step -1 until 1 do -- ..........
C
      DO 170 II = 1,N
        I = N + 1 - II
        L = I - 1
        H = 0.0
        SCALE = 0.0
        IF (L.LT.1) GO TO 40
C
C       .......... scale row (algol tol then not needed) ..........
C
        DO 20 K = 1,L
          SCALE = SCALE + ABS(D(K))
   20   CONTINUE
C
        IF (SCALE.NE.0.0) GO TO 50
C
        DO 30 J = 1,L
          D(J) = A(L,J)
          A(L,J) = A(I,J)
          A(I,J) = 0.0
   30   CONTINUE
C
   40   CONTINUE
        E(I) = 0.0
        E2(I) = 0.0
        GO TO 170
C
   50   CONTINUE
        DO 60 K = 1,L
          D(K) = D(K)/SCALE
          H = H + D(K)*D(K)
   60   CONTINUE
C
        E2(I) = SCALE*SCALE*H
        F = D(L)
        G = -SIGN(SQRT(H),F)
        E(I) = SCALE*G
        H = H - F*G
        D(L) = F - G
        IF (L.EQ.1) GO TO 150
C
C       .......... form a*u ..........
C
        DO 70 J = 1,L
          E(J) = 0.0
   70   CONTINUE
C
        DO 100 J = 1,L
          F = D(J)
          G = E(J) + A(J,J)*F
          JP1 = J + 1
          IF (L.LT.JP1) GO TO 90
C
          DO 80 K = JP1,L
            G = G + A(K,J)*D(K)
            E(K) = E(K) + A(K,J)*F
   80     CONTINUE
C
   90     CONTINUE
          E(J) = G
  100   CONTINUE
C
C       .......... form p ..........
C
        F = 0.0
C
        DO 110 J = 1,L
          E(J) = E(J)/H
          F = F + E(J)*D(J)
  110   CONTINUE
C

        H = F/ (H + H)
C       .......... form q ..........

        DO 120 J = 1,L
          E(J) = E(J) - H*D(J)
  120   CONTINUE
C
C       .......... form reduced a ..........
C
        DO 140 J = 1,L
          F = D(J)
          G = E(J)
C
          DO 130 K = J,L
            A(K,J) = A(K,J) - F*E(K) - G*D(K)
  130     CONTINUE
C
  140   CONTINUE
C
  150   CONTINUE
        DO 160 J = 1,L
          F = D(J)
          D(J) = A(L,J)
          A(L,J) = A(I,J)
          A(I,J) = F*SCALE
  160   CONTINUE
C
  170 CONTINUE
C
      END
      SUBROUTINE TRED2 (NM,N,A,D,E,Z)
C
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL*8 A(NM,N),D(N),E(N),Z(NM,N)
      REAL*8 F,G,H,HH,SCALE
C
C     this subroutine is a translation of the algol procedure tred2,
C     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
C     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
C
C     this subroutine reduces a real symmetric matrix to a
C     symmetric tridiagonal matrix using and accumulating
C     orthogonal similarity transformations.
C
C     on input
C
C     nm must be set to the row dimension of two-dimensional
C     array parameters as declared in the calling program
C     dimension statement.
C
C     n is the order of the matrix.
C
C     a contains the real symmetric input matrix.  only the
C     lower triangle of the matrix need be supplied.
C
C     on output
C
C     d contains the diagonal elements of the tridiagonal matrix.
C
C     e contains the subdiagonal elements of the tridiagonal
C     matrix in its last n-1 positions.  e(1) is set to zero.
C
C     z contains the orthogonal transformation matrix
C     produced in the reduction.
C
C     a and z may coincide.  if distinct, a is unaltered.
C
C     questions and comments should be directed to burton s. garbow,
C     mathematics and computer science div, argonne national laboratory
C
C     this version dated august 1983.
C
C     ------------------------------------------------------------------
C
      DO 20 I = 1,N
C
        DO 10 J = I,N
          Z(J,I) = A(J,I)
   10   CONTINUE
C
        D(I) = A(N,I)
   20 CONTINUE
C
      IF (N.EQ.1) GO TO 250
C
C     .......... for i=n step -1 until 2 do -- ..........
C
      DO 170 II = 2,N
        I = N + 2 - II
        L = I - 1
        H = 0.0
        SCALE = 0.0
        IF (L.LT.2) GO TO 40
C
C       .......... scale row (algol tol then not needed) ..........
C
        DO 30 K = 1,L
          SCALE = SCALE + ABS(D(K))
   30   CONTINUE
C
        IF (SCALE.NE.0.0) GO TO 60
   40   CONTINUE
        E(I) = D(L)
C
        DO 50 J = 1,L
          D(J) = Z(L,J)
          Z(I,J) = 0.0
          Z(J,I) = 0.0
   50   CONTINUE
C
        GO TO 160
C
   60   CONTINUE
        DO 70 K = 1,L
          D(K) = D(K)/SCALE
          H = H + D(K)*D(K)
   70   CONTINUE
C
        F = D(L)
        G = -SIGN(SQRT(H),F)
        E(I) = SCALE*G
        H = H - F*G
        D(L) = F - G
C
C       .......... form a*u ..........
C
        DO 80 J = 1,L
          E(J) = 0.0
   80   CONTINUE
C
        DO 110 J = 1,L
          F = D(J)
          Z(J,I) = F
          G = E(J) + Z(J,J)*F
          JP1 = J + 1
          IF (L.LT.JP1) GO TO 100
C
          DO 90 K = JP1,L
            G = G + Z(K,J)*D(K)
            E(K) = E(K) + Z(K,J)*F
   90     CONTINUE
C
  100     CONTINUE
          E(J) = G
  110   CONTINUE
C
C       .......... form p ..........
C
        F = 0.0
C
        DO 120 J = 1,L
          E(J) = E(J)/H
          F = F + E(J)*D(J)
  120   CONTINUE
C
        HH = F/ (H + H)
C
C       .......... form q ..........
C
        DO 130 J = 1,L
          E(J) = E(J) - HH*D(J)
  130   CONTINUE
C
C       .......... form reduced a ..........
C
        DO 150 J = 1,L
          F = D(J)
          G = E(J)
C
          DO 140 K = J,L
            Z(K,J) = Z(K,J) - F*E(K) - G*D(K)
  140     CONTINUE
C
          D(J) = Z(L,J)
          Z(I,J) = 0.0
  150   CONTINUE
C
  160   CONTINUE
        D(I) = H
  170 CONTINUE
C
C     .......... accumulation of transformation matrices ..........
C
      DO 240 I = 2,N
        L = I - 1
        Z(N,L) = Z(L,L)
        Z(L,L) = 1.0
        H = D(I)
        IF (H.EQ.0.0) GO TO 220
C
        DO 180 K = 1,L
          D(K) = Z(K,I)/H
  180   CONTINUE
C
        DO 210 J = 1,L
          G = 0.0
C
          DO 190 K = 1,L
            G = G + Z(K,I)*Z(K,J)
  190     CONTINUE
C
          DO 200 K = 1,L
            Z(K,J) = Z(K,J) - G*D(K)
  200     CONTINUE
  210   CONTINUE
C
  220   CONTINUE
        DO 230 K = 1,L
          Z(K,I) = 0.0
  230   CONTINUE
C
  240 CONTINUE
C
  250 CONTINUE
      DO 260 I = 1,N
        D(I) = Z(N,I)
        Z(N,I) = 0.0
  260 CONTINUE
C
      Z(N,N) = 1.0
      E(1) = 0.0
      END
