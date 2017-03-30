module ode_lsoda_aux2

  implicit none
  public

contains
  !*==DGEFA.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
  !DECK DGEFA
        SUBROUTINE DGEFA(A,Lda,N,Ipvt,Info)
        IMPLICIT NONE
  !*--DGEFA5
  !***BEGIN PROLOGUE  DGEFA
  !***PURPOSE  Factor a matrix using Gaussian elimination.
  !***CATEGORY  D2A1
  !***TYPE      DOUBLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)
  !***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
  !             MATRIX FACTORIZATION
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     DGEFA factors a double precision matrix by Gaussian elimination.
  !
  !     DGEFA is usually called by DGECO, but it can be called
  !     directly with a saving in time if  RCOND  is not needed.
  !     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
  !
  !     On Entry
  !
  !        A       DOUBLE PRECISION(LDA, N)
  !                the matrix to be factored.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  A .
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !     On Return
  !
  !        A       an upper triangular matrix and the multipliers
  !                which were used to obtain it.
  !                The factorization can be written  A = L*U  where
  !                L  is a product of permutation and unit lower
  !                triangular matrices and  U  is upper triangular.
  !
  !        IPVT    INTEGER(N)
  !                an integer vector of pivot indices.
  !
  !        INFO    INTEGER
  !                = 0  normal value.
  !                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
  !                     condition for this subroutine, but it does
  !                     indicate that DGESL or DGEDI will divide by zero
  !                     if called.  Use  RCOND  in DGECO for a reliable
  !                     indication of singularity.
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DGEFA
        INTEGER Lda , N , Ipvt(*) , Info
        DOUBLE PRECISION A(Lda,*)
  !
        DOUBLE PRECISION t
        INTEGER j , k , kp1 , l , nm1
  !
  !     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
  !
  !***FIRST EXECUTABLE STATEMENT  DGEFA
        Info = 0
        nm1 = N - 1
        IF ( nm1.GE.1 ) THEN
           DO k = 1 , nm1
              kp1 = k + 1
  !
  !        FIND L = PIVOT INDEX
  !
              l = IDAMAX(N-k+1,A(k,k),1) + k - 1
              Ipvt(k) = l
  !
  !        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
  !
              IF ( A(l,k).EQ.0.0D0 ) THEN
                 Info = k
              ELSE
  !
  !           INTERCHANGE IF NECESSARY
  !
                 IF ( l.NE.k ) THEN
                    t = A(l,k)
                    A(l,k) = A(k,k)
                    A(k,k) = t
                 ENDIF
  !
  !           COMPUTE MULTIPLIERS
  !
                 t = -1.0D0/A(k,k)
                 CALL DSCAL(N-k,t,A(k+1,k),1)
  !
  !           ROW ELIMINATION WITH COLUMN INDEXING
  !
                 DO j = kp1 , N
                    t = A(l,j)
                    IF ( l.NE.k ) THEN
                       A(l,j) = A(k,j)
                       A(k,j) = t
                    ENDIF
                    CALL DAXPY(N-k,t,A(k+1,k),1,A(k+1,j),1)
                 ENDDO
              ENDIF
           ENDDO
        ENDIF
        Ipvt(N) = N
        IF ( A(N,N).EQ.0.0D0 ) Info = N
        END
  !*==DGESL.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
  !DECK DGESL
        SUBROUTINE DGESL(A,Lda,N,Ipvt,B,Job)
        IMPLICIT NONE
  !*--DGESL122
  !***BEGIN PROLOGUE  DGESL
  !***PURPOSE  Solve the real system A*X=B or TRANS(A)*X=B using the
  !            factors computed by DGECO or DGEFA.
  !***CATEGORY  D2A1
  !***TYPE      DOUBLE PRECISION (SGESL-S, DGESL-D, CGESL-C)
  !***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     DGESL solves the double precision system
  !     A * X = B  or  TRANS(A) * X = B
  !     using the factors computed by DGECO or DGEFA.
  !
  !     On Entry
  !
  !        A       DOUBLE PRECISION(LDA, N)
  !                the output from DGECO or DGEFA.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  A .
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !        IPVT    INTEGER(N)
  !                the pivot vector from DGECO or DGEFA.
  !
  !        B       DOUBLE PRECISION(N)
  !                the right hand side vector.
  !
  !        JOB     INTEGER
  !                = 0         to solve  A*X = B ,
  !                = nonzero   to solve  TRANS(A)*X = B  where
  !                            TRANS(A)  is the transpose.
  !
  !     On Return
  !
  !        B       the solution vector  X .
  !
  !     Error Condition
  !
  !        A division by zero will occur if the input factor contains a
  !        zero on the diagonal.  Technically this indicates singularity
  !        but it is often caused by improper arguments or improper
  !        setting of LDA .  It will not occur if the subroutines are
  !        called correctly and if DGECO has set RCOND .GT. 0.0
  !        or DGEFA has set INFO .EQ. 0 .
  !
  !     To compute  INVERSE(A) * C  where  C  is a matrix
  !     with  P  columns
  !           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
  !           IF (RCOND is too small) GO TO ...
  !           DO 10 J = 1, P
  !              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
  !        10 CONTINUE
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  DAXPY, DDOT
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DGESL
        INTEGER Lda , N , Ipvt(*) , Job
        DOUBLE PRECISION A(Lda,*) , B(*)
  !
        DOUBLE PRECISION t
        INTEGER k , kb , l , nm1
  !***FIRST EXECUTABLE STATEMENT  DGESL
        nm1 = N - 1
        IF ( Job.NE.0 ) THEN
  !
  !        JOB = NONZERO, SOLVE  TRANS(A) * X = B
  !        FIRST SOLVE  TRANS(U)*Y = B
  !
           DO k = 1 , N
              t = DDOT(k-1,A(1,k),1,B(1),1)
              B(k) = (B(k)-t)/A(k,k)
           ENDDO
  !
  !        NOW SOLVE TRANS(L)*X = Y
  !
           IF ( nm1.GE.1 ) THEN
              DO kb = 1 , nm1
                 k = N - kb
                 B(k) = B(k) + DDOT(N-k,A(k+1,k),1,B(k+1),1)
                 l = Ipvt(k)
                 IF ( l.NE.k ) THEN
                    t = B(l)
                    B(l) = B(k)
                    B(k) = t
                 ENDIF
              ENDDO
           ENDIF
        ELSE
  !
  !        JOB = 0 , SOLVE  A * X = B
  !        FIRST SOLVE  L*Y = B
  !
           IF ( nm1.GE.1 ) THEN
              DO k = 1 , nm1
                 l = Ipvt(k)
                 t = B(l)
                 IF ( l.NE.k ) THEN
                    B(l) = B(k)
                    B(k) = t
                 ENDIF
                 CALL DAXPY(N-k,t,A(k+1,k),1,B(k+1),1)
              ENDDO
           ENDIF
  !
  !        NOW SOLVE  U*X = Y
  !
           DO kb = 1 , N
              k = N + 1 - kb
              B(k) = B(k)/A(k,k)
              t = -B(k)
              CALL DAXPY(k-1,t,A(1,k),1,B(1),1)
           ENDDO
        ENDIF
        END
  !*==DGBFA.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
  !DECK DGBFA
        SUBROUTINE DGBFA(Abd,Lda,N,Ml,Mu,Ipvt,Info)
        IMPLICIT NONE
  !*--DGBFA253
  !***BEGIN PROLOGUE  DGBFA
  !***PURPOSE  Factor a band matrix using Gaussian elimination.
  !***CATEGORY  D2A2
  !***TYPE      DOUBLE PRECISION (SGBFA-S, DGBFA-D, CGBFA-C)
  !***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     DGBFA factors a double precision band matrix by elimination.
  !
  !     DGBFA is usually called by DGBCO, but it can be called
  !     directly with a saving in time if  RCOND  is not needed.
  !
  !     On Entry
  !
  !        ABD     DOUBLE PRECISION(LDA, N)
  !                contains the matrix in band storage.  The columns
  !                of the matrix are stored in the columns of  ABD  and
  !                the diagonals of the matrix are stored in rows
  !                ML+1 through 2*ML+MU+1 of  ABD .
  !                See the comments below for details.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  ABD .
  !                LDA must be .GE. 2*ML + MU + 1 .
  !
  !        N       INTEGER
  !                the order of the original matrix.
  !
  !        ML      INTEGER
  !                number of diagonals below the main diagonal.
  !                0 .LE. ML .LT.  N .
  !
  !        MU      INTEGER
  !                number of diagonals above the main diagonal.
  !                0 .LE. MU .LT.  N .
  !                More efficient if  ML .LE. MU .
  !     On Return
  !
  !        ABD     an upper triangular matrix in band storage and
  !                the multipliers which were used to obtain it.
  !                The factorization can be written  A = L*U  where
  !                L  is a product of permutation and unit lower
  !                triangular matrices and  U  is upper triangular.
  !
  !        IPVT    INTEGER(N)
  !                an integer vector of pivot indices.
  !
  !        INFO    INTEGER
  !                = 0  normal value.
  !                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
  !                     condition for this subroutine, but it does
  !                     indicate that DGBSL will divide by zero if
  !                     called.  Use  RCOND  in DGBCO for a reliable
  !                     indication of singularity.
  !
  !     Band Storage
  !
  !           If  A  is a band matrix, the following program segment
  !           will set up the input.
  !
  !                   ML = (band width below the diagonal)
  !                   MU = (band width above the diagonal)
  !                   M = ML + MU + 1
  !                   DO 20 J = 1, N
  !                      I1 = MAX(1, J-MU)
  !                      I2 = MIN(N, J+ML)
  !                      DO 10 I = I1, I2
  !                         K = I - J + M
  !                         ABD(K,J) = A(I,J)
  !                10    CONTINUE
  !                20 CONTINUE
  !
  !           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
  !           In addition, the first  ML  rows in  ABD  are used for
  !           elements generated during the triangularization.
  !           The total number of rows needed in  ABD  is  2*ML+MU+1 .
  !           The  ML+MU by ML+MU  upper left triangle and the
  !           ML by ML  lower right triangle are not referenced.
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DGBFA
        INTEGER Lda , N , Ml , Mu , Ipvt(*) , Info
        DOUBLE PRECISION Abd(Lda,*)
  !
        DOUBLE PRECISION t
        INTEGER i , i0 , j , ju , jz , j0 , j1 , k , kp1 , l ,   &
       &        lm , m , mm , nm1
  !
  !***FIRST EXECUTABLE STATEMENT  DGBFA
        m = Ml + Mu + 1
        Info = 0
  !
  !     ZERO INITIAL FILL-IN COLUMNS
  !
        j0 = Mu + 2
        j1 = MIN(N,m) - 1
        IF ( j1.GE.j0 ) THEN
           DO jz = j0 , j1
              i0 = m + 1 - jz
              DO i = i0 , Ml
                 Abd(i,jz) = 0.0D0
              ENDDO
           ENDDO
        ENDIF
        jz = j1
        ju = 0
  !
  !     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
  !
        nm1 = N - 1
        IF ( nm1.GE.1 ) THEN
           DO k = 1 , nm1
              kp1 = k + 1
  !
  !        ZERO NEXT FILL-IN COLUMN
  !
              jz = jz + 1
              IF ( jz.LE.N ) THEN
                 IF ( Ml.GE.1 ) THEN
                    DO i = 1 , Ml
                       Abd(i,jz) = 0.0D0
                    ENDDO
                 ENDIF
              ENDIF
  !
  !        FIND L = PIVOT INDEX
  !
              lm = MIN(Ml,N-k)
              l = IDAMAX(lm+1,Abd(m,k),1) + m - 1
              Ipvt(k) = l + k - m
  !
  !        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
  !
              IF ( Abd(l,k).EQ.0.0D0 ) THEN
                 Info = k
              ELSE
  !
  !           INTERCHANGE IF NECESSARY
  !
                 IF ( l.NE.m ) THEN
                    t = Abd(l,k)
                    Abd(l,k) = Abd(m,k)
                    Abd(m,k) = t
                 ENDIF
  !
  !           COMPUTE MULTIPLIERS
  !
                 t = -1.0D0/Abd(m,k)
                 CALL DSCAL(lm,t,Abd(m+1,k),1)
  !
  !           ROW ELIMINATION WITH COLUMN INDEXING
  !
                 ju = MIN(MAX(ju,Mu+Ipvt(k)),N)
                 mm = m
                 IF ( ju.GE.kp1 ) THEN
                    DO j = kp1 , ju
                       l = l - 1
                       mm = mm - 1
                       t = Abd(l,j)
                       IF ( l.NE.mm ) THEN
                          Abd(l,j) = Abd(mm,j)
                          Abd(mm,j) = t
                       ENDIF
                       CALL DAXPY(lm,t,Abd(m+1,k),1,Abd(mm+1,j),1)
                    ENDDO
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
        Ipvt(N) = N
        IF ( Abd(m,N).EQ.0.0D0 ) Info = N
        END
  !*==DGBSL.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
  !DECK DGBSL
        SUBROUTINE DGBSL(Abd,Lda,N,Ml,Mu,Ipvt,B,Job)
        IMPLICIT NONE
  !*--DGBSL442
  !***BEGIN PROLOGUE  DGBSL
  !***PURPOSE  Solve the real band system A*X=B or TRANS(A)*X=B using
  !            the factors computed by DGBCO or DGBFA.
  !***CATEGORY  D2A2
  !***TYPE      DOUBLE PRECISION (SGBSL-S, DGBSL-D, CGBSL-C)
  !***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     DGBSL solves the double precision band system
  !     A * X = B  or  TRANS(A) * X = B
  !     using the factors computed by DGBCO or DGBFA.
  !
  !     On Entry
  !
  !        ABD     DOUBLE PRECISION(LDA, N)
  !                the output from DGBCO or DGBFA.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  ABD .
  !
  !        N       INTEGER
  !                the order of the original matrix.
  !
  !        ML      INTEGER
  !                number of diagonals below the main diagonal.
  !
  !        MU      INTEGER
  !                number of diagonals above the main diagonal.
  !
  !        IPVT    INTEGER(N)
  !                the pivot vector from DGBCO or DGBFA.
  !
  !        B       DOUBLE PRECISION(N)
  !                the right hand side vector.
  !
  !        JOB     INTEGER
  !                = 0         to solve  A*X = B ,
  !                = nonzero   to solve  TRANS(A)*X = B , where
  !                            TRANS(A)  is the transpose.
  !
  !     On Return
  !
  !        B       the solution vector  X .
  !
  !     Error Condition
  !
  !        A division by zero will occur if the input factor contains a
  !        zero on the diagonal.  Technically this indicates singularity
  !        but it is often caused by improper arguments or improper
  !        setting of LDA .  It will not occur if the subroutines are
  !        called correctly and if DGBCO has set RCOND .GT. 0.0
  !        or DGBFA has set INFO .EQ. 0 .
  !
  !     To compute  INVERSE(A) * C  where  C  is a matrix
  !     with  P  columns
  !           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
  !           IF (RCOND is too small) GO TO ...
  !           DO 10 J = 1, P
  !              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
  !        10 CONTINUE
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  DAXPY, DDOT
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DGBSL
        INTEGER Lda , N , Ml , Mu , Ipvt(*) , Job
        DOUBLE PRECISION Abd(Lda,*) , B(*)
  !
        DOUBLE PRECISION t
        INTEGER k , kb , l , la , lb , lm , m , nm1
  !***FIRST EXECUTABLE STATEMENT  DGBSL
        m = Mu + Ml + 1
        nm1 = N - 1
        IF ( Job.NE.0 ) THEN
  !
  !        JOB = NONZERO, SOLVE  TRANS(A) * X = B
  !        FIRST SOLVE  TRANS(U)*Y = B
  !
           DO k = 1 , N
              lm = MIN(k,m) - 1
              la = m - lm
              lb = k - lm
              t = DDOT(lm,Abd(la,k),1,B(lb),1)
              B(k) = (B(k)-t)/Abd(m,k)
           ENDDO
  !
  !        NOW SOLVE TRANS(L)*X = Y
  !
           IF ( Ml.NE.0 ) THEN
              IF ( nm1.GE.1 ) THEN
                 DO kb = 1 , nm1
                    k = N - kb
                    lm = MIN(Ml,N-k)
                    B(k) = B(k) + DDOT(lm,Abd(m+1,k),1,B(k+1),1)
                    l = Ipvt(k)
                    IF ( l.NE.k ) THEN
                       t = B(l)
                       B(l) = B(k)
                       B(k) = t
                    ENDIF
                 ENDDO
              ENDIF
           ENDIF
        ELSE
  !
  !        JOB = 0 , SOLVE  A * X = B
  !        FIRST SOLVE L*Y = B
  !
           IF ( Ml.NE.0 ) THEN
              IF ( nm1.GE.1 ) THEN
                 DO k = 1 , nm1
                    lm = MIN(Ml,N-k)
                    l = Ipvt(k)
                    t = B(l)
                    IF ( l.NE.k ) THEN
                       B(l) = B(k)
                       B(k) = t
                    ENDIF
                    CALL DAXPY(lm,t,Abd(m+1,k),1,B(k+1),1)
                 ENDDO
              ENDIF
           ENDIF
  !
  !        NOW SOLVE  U*X = Y
  !
           DO kb = 1 , N
              k = N + 1 - kb
              B(k) = B(k)/Abd(m,k)
              lm = MIN(k,m) - 1
              la = m - lm
              lb = k - lm
              t = -B(k)
              CALL DAXPY(lm,t,Abd(la,k),1,B(lb),1)
           ENDDO
        ENDIF
        END
  !*==DAXPY.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
  !DECK DAXPY
        SUBROUTINE DAXPY(N,Da,Dx,Incx,Dy,Incy)
        IMPLICIT NONE
  !*--DAXPY593
  !*** Start of declarations inserted by SPAG
        INTEGER i , Incx , Incy , ix , iy , m , mp1 , N , ns
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DAXPY
  !***PURPOSE  Compute a constant times a vector plus a vector.
  !***CATEGORY  D1A7
  !***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
  !***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
  !***AUTHOR  Lawson, C. L., (JPL)
  !           Hanson, R. J., (SNLA)
  !           Kincaid, D. R., (U. of Texas)
  !           Krogh, F. T., (JPL)
  !***DESCRIPTION
  !
  !                B L A S  Subprogram
  !    Description of Parameters
  !
  !     --Input--
  !        N  number of elements in input vector(s)
  !       DA  double precision scalar multiplier
  !       DX  double precision vector with N elements
  !     INCX  storage spacing between elements of DX
  !       DY  double precision vector with N elements
  !     INCY  storage spacing between elements of DY
  !
  !     --Output--
  !       DY  double precision result (unchanged if N .LE. 0)
  !
  !     Overwrite double precision DY with double precision DA*DX + DY.
  !     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
  !       DY(LY+I*INCY),
  !     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
  !     defined in a similar way using INCY.
  !
  !***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !                 Krogh, Basic linear algebra subprograms for Fortran
  !                 usage, Algorithm No. 539, Transactions on Mathematical
  !                 Software 5, 3 (September 1979), pp. 308-323.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DAXPY
        DOUBLE PRECISION Dx(*) , Dy(*) , Da
  !***FIRST EXECUTABLE STATEMENT  DAXPY
        IF ( N.LE.0 .OR. Da.EQ.0.0D0 ) RETURN
        IF ( Incx.EQ.Incy ) THEN
           IF ( Incx.LT.1 ) THEN
           ELSEIF ( Incx.EQ.1 ) THEN
  !
  !     Code for both increments equal to 1.
  !
  !     Clean-up loop so remaining vector length is a multiple of 4.
  !
              m = MOD(N,4)
              IF ( m.NE.0 ) THEN
                 DO i = 1 , m
                    Dy(i) = Dy(i) + Da*Dx(i)
                 ENDDO
                 IF ( N.LT.4 ) RETURN
              ENDIF
              GOTO 100
           ELSE
  !
  !     Code for equal, positive, non-unit increments.
  !
              ns = N*Incx
              DO i = 1 , ns , Incx
                 Dy(i) = Da*Dx(i) + Dy(i)
              ENDDO
              GOTO 99999
           ENDIF
        ENDIF
  !
  !     Code for unequal or nonpositive increments.
  !
        ix = 1
        iy = 1
        IF ( Incx.LT.0 ) ix = (-N+1)*Incx + 1
        IF ( Incy.LT.0 ) iy = (-N+1)*Incy + 1
        DO i = 1 , N
           Dy(iy) = Dy(iy) + Da*Dx(ix)
           ix = ix + Incx
           iy = iy + Incy
        ENDDO
        RETURN
   100  mp1 = m + 1
        DO i = mp1 , N , 4
           Dy(i) = Dy(i) + Da*Dx(i)
           Dy(i+1) = Dy(i+1) + Da*Dx(i+1)
           Dy(i+2) = Dy(i+2) + Da*Dx(i+2)
           Dy(i+3) = Dy(i+3) + Da*Dx(i+3)
        ENDDO
        RETURN
  99999 END
  !*==DCOPY.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
  !DECK DCOPY
        SUBROUTINE DCOPY(N,Dx,Incx,Dy,Incy)
        IMPLICIT NONE
  !*--DCOPY697
  !*** Start of declarations inserted by SPAG
        INTEGER i , Incx , Incy , ix , iy , m , mp1 , N , ns
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DCOPY
  !***PURPOSE  Copy a vector.
  !***CATEGORY  D1A5
  !***TYPE      DOUBLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
  !***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
  !***AUTHOR  Lawson, C. L., (JPL)
  !           Hanson, R. J., (SNLA)
  !           Kincaid, D. R., (U. of Texas)
  !           Krogh, F. T., (JPL)
  !***DESCRIPTION
  !
  !                B L A S  Subprogram
  !    Description of Parameters
  !
  !     --Input--
  !        N  number of elements in input vector(s)
  !       DX  double precision vector with N elements
  !     INCX  storage spacing between elements of DX
  !       DY  double precision vector with N elements
  !     INCY  storage spacing between elements of DY
  !
  !     --Output--
  !       DY  copy of vector DX (unchanged if N .LE. 0)
  !
  !     Copy double precision DX to double precision DY.
  !     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
  !     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
  !     defined in a similar way using INCY.
  !
  !***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !                 Krogh, Basic linear algebra subprograms for Fortran
  !                 usage, Algorithm No. 539, Transactions on Mathematical
  !                 Software 5, 3 (September 1979), pp. 308-323.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DCOPY
        DOUBLE PRECISION Dx(*) , Dy(*)
  !***FIRST EXECUTABLE STATEMENT  DCOPY
        IF ( N.LE.0 ) RETURN
        IF ( Incx.EQ.Incy ) THEN
           IF ( Incx.LT.1 ) THEN
           ELSEIF ( Incx.EQ.1 ) THEN
  !
  !     Code for both increments equal to 1.
  !
  !     Clean-up loop so remaining vector length is a multiple of 7.
  !
              m = MOD(N,7)
              IF ( m.NE.0 ) THEN
                 DO i = 1 , m
                    Dy(i) = Dx(i)
                 ENDDO
                 IF ( N.LT.7 ) RETURN
              ENDIF
              GOTO 100
           ELSE
  !
  !     Code for equal, positive, non-unit increments.
  !
              ns = N*Incx
              DO i = 1 , ns , Incx
                 Dy(i) = Dx(i)
              ENDDO
              GOTO 99999
           ENDIF
        ENDIF
  !
  !     Code for unequal or nonpositive increments.
  !
        ix = 1
        iy = 1
        IF ( Incx.LT.0 ) ix = (-N+1)*Incx + 1
        IF ( Incy.LT.0 ) iy = (-N+1)*Incy + 1
        DO i = 1 , N
           Dy(iy) = Dx(ix)
           ix = ix + Incx
           iy = iy + Incy
        ENDDO
        RETURN
   100  mp1 = m + 1
        DO i = mp1 , N , 7
           Dy(i) = Dx(i)
           Dy(i+1) = Dx(i+1)
           Dy(i+2) = Dx(i+2)
           Dy(i+3) = Dx(i+3)
           Dy(i+4) = Dx(i+4)
           Dy(i+5) = Dx(i+5)
           Dy(i+6) = Dx(i+6)
        ENDDO
        RETURN
  99999 END
  !*==DDOT.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
  !DECK DDOT
        DOUBLE PRECISION FUNCTION DDOT(N,Dx,Incx,Dy,Incy)
        IMPLICIT NONE
  !*--DDOT802
  !*** Start of declarations inserted by SPAG
        INTEGER i , Incx , Incy , ix , iy , m , mp1 , N , ns
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DDOT
  !***PURPOSE  Compute the inner product of two vectors.
  !***CATEGORY  D1A4
  !***TYPE      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
  !***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
  !***AUTHOR  Lawson, C. L., (JPL)
  !           Hanson, R. J., (SNLA)
  !           Kincaid, D. R., (U. of Texas)
  !           Krogh, F. T., (JPL)
  !***DESCRIPTION
  !
  !                B L A S  Subprogram
  !    Description of Parameters
  !
  !     --Input--
  !        N  number of elements in input vector(s)
  !       DX  double precision vector with N elements
  !     INCX  storage spacing between elements of DX
  !       DY  double precision vector with N elements
  !     INCY  storage spacing between elements of DY
  !
  !     --Output--
  !     DDOT  double precision dot product (zero if N .LE. 0)
  !
  !     Returns the dot product of double precision DX and DY.
  !     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY),
  !     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
  !     defined in a similar way using INCY.
  !
  !***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !                 Krogh, Basic linear algebra subprograms for Fortran
  !                 usage, Algorithm No. 539, Transactions on Mathematical
  !                 Software 5, 3 (September 1979), pp. 308-323.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DDOT
        DOUBLE PRECISION Dx(*) , Dy(*)
  !***FIRST EXECUTABLE STATEMENT  DDOT
        DDOT = 0.0D0
        IF ( N.LE.0 ) RETURN
        IF ( Incx.EQ.Incy ) THEN
           IF ( Incx.LT.1 ) THEN
           ELSEIF ( Incx.EQ.1 ) THEN
  !
  !     Code for both increments equal to 1.
  !
  !     Clean-up loop so remaining vector length is a multiple of 5.
  !
              m = MOD(N,5)
              IF ( m.NE.0 ) THEN
                 DO i = 1 , m
                    DDOT = DDOT + Dx(i)*Dy(i)
                 ENDDO
                 IF ( N.LT.5 ) RETURN
              ENDIF
              GOTO 100
           ELSE
  !
  !     Code for equal, positive, non-unit increments.
  !
              ns = N*Incx
              DO i = 1 , ns , Incx
                 DDOT = DDOT + Dx(i)*Dy(i)
              ENDDO
              GOTO 99999
           ENDIF
        ENDIF
  !
  !     Code for unequal or nonpositive increments.
  !
        ix = 1
        iy = 1
        IF ( Incx.LT.0 ) ix = (-N+1)*Incx + 1
        IF ( Incy.LT.0 ) iy = (-N+1)*Incy + 1
        DO i = 1 , N
           DDOT = DDOT + Dx(ix)*Dy(iy)
           ix = ix + Incx
           iy = iy + Incy
        ENDDO
        RETURN
   100  mp1 = m + 1
        DO i = mp1 , N , 5
           DDOT = DDOT + Dx(i)*Dy(i) + Dx(i+1)*Dy(i+1) + Dx(i+2)*Dy(i+2)  &
       &          + Dx(i+3)*Dy(i+3) + Dx(i+4)*Dy(i+4)
        ENDDO
        RETURN
  99999 END
  !*==DNRM2.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
  !DECK DNRM2
        DOUBLE PRECISION FUNCTION DNRM2(N,Dx,Incx)
        IMPLICIT NONE
  !*--DNRM2903
  !*** Start of declarations inserted by SPAG
        INTEGER i , Incx , j , N , nn
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DNRM2
  !***PURPOSE  Compute the Euclidean length (L2 norm) of a vector.
  !***CATEGORY  D1A3B
  !***TYPE      DOUBLE PRECISION (SNRM2-S, DNRM2-D, SCNRM2-C)
  !***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,
  !             LINEAR ALGEBRA, UNITARY, VECTOR
  !***AUTHOR  Lawson, C. L., (JPL)
  !           Hanson, R. J., (SNLA)
  !           Kincaid, D. R., (U. of Texas)
  !           Krogh, F. T., (JPL)
  !***DESCRIPTION
  !
  !                B L A S  Subprogram
  !    Description of parameters
  !
  !     --Input--
  !        N  number of elements in input vector(s)
  !       DX  double precision vector with N elements
  !     INCX  storage spacing between elements of DX
  !
  !     --Output--
  !    DNRM2  double precision result (zero if N .LE. 0)
  !
  !     Euclidean norm of the N-vector stored in DX with storage
  !     increment INCX.
  !     If N .LE. 0, return with result = 0.
  !     If N .GE. 1, then INCX must be .GE. 1
  !
  !     Four phase method using two built-in constants that are
  !     hopefully applicable to all machines.
  !         CUTLO = maximum of  SQRT(U/EPS)  over all known machines.
  !         CUTHI = minimum of  SQRT(V)      over all known machines.
  !     where
  !         EPS = smallest no. such that EPS + 1. .GT. 1.
  !         U   = smallest positive no.   (underflow limit)
  !         V   = largest  no.            (overflow  limit)
  !
  !     Brief outline of algorithm.
  !
  !     Phase 1 scans zero components.
  !     move to phase 2 when a component is nonzero and .LE. CUTLO
  !     move to phase 3 when a component is .GT. CUTLO
  !     move to phase 4 when a component is .GE. CUTHI/M
  !     where M = N for X() real and M = 2*N for complex.
  !
  !     Values for CUTLO and CUTHI.
  !     From the environmental parameters listed in the IMSL converter
  !     document the limiting values are as follows:
  !     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
  !                   Univac and DEC at 2**(-103)
  !                   Thus CUTLO = 2**(-51) = 4.44089E-16
  !     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
  !                   Thus CUTHI = 2**(63.5) = 1.30438E19
  !     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
  !                   Thus CUTLO = 2**(-33.5) = 8.23181D-11
  !     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
  !     DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
  !     DATA CUTLO, CUTHI /4.441E-16,  1.304E19/
  !
  !***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !                 Krogh, Basic linear algebra subprograms for Fortran
  !                 usage, Algorithm No. 539, Transactions on Mathematical
  !                 Software 5, 3 (September 1979), pp. 308-323.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DNRM2
        INTEGER next
        DOUBLE PRECISION Dx(*) , cutlo , cuthi , hitest , sum , xmax ,    &
       &                 zero , one
        SAVE cutlo , cuthi , zero , one
        DATA zero , one/0.0D0 , 1.0D0/
  !
        DATA cutlo , cuthi/8.232D-11 , 1.304D19/
  !***FIRST EXECUTABLE STATEMENT  DNRM2
        IF ( N.GT.0 ) THEN
  !
           ASSIGN 200 TO next
           sum = zero
           nn = N*Incx
  !
  !                                                 BEGIN MAIN LOOP
  !
           i = 1
        ELSE
           DNRM2 = zero
           GOTO 99999
        ENDIF
   100  GOTO next
   200  IF ( ABS(Dx(i)).GT.cutlo ) GOTO 800
        ASSIGN 300 TO next
        xmax = zero
  !
  !                        PHASE 1.  SUM IS ZERO
  !
   300  IF ( Dx(i).EQ.zero ) GOTO 900
        IF ( ABS(Dx(i)).GT.cutlo ) GOTO 800
  !
  !                                PREPARE FOR PHASE 2.
  !
        ASSIGN 600 TO next
        GOTO 500
  !
  !                                PREPARE FOR PHASE 4.
  !
   400  i = j
        ASSIGN 700 TO next
        sum = (sum/Dx(i))/Dx(i)
   500  xmax = ABS(Dx(i))
  !
        sum = sum + (Dx(i)/xmax)**2
        GOTO 900
  !
  !                   PHASE 2.  SUM IS SMALL.
  !                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
  !
   600  IF ( ABS(Dx(i)).GT.cutlo ) THEN
  !
  !                  PREPARE FOR PHASE 3.
  !
           sum = (sum*xmax)*xmax
           GOTO 800
        ENDIF
  !
  !                     COMMON CODE FOR PHASES 2 AND 4.
  !                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
  !
   700  IF ( ABS(Dx(i)).LE.xmax ) THEN
           sum = sum + (Dx(i)/xmax)**2
        ELSE
           sum = one + sum*(xmax/Dx(i))**2
           xmax = ABS(Dx(i))
        ENDIF
        GOTO 900
  !
  !     FOR REAL OR D.P. SET HITEST = CUTHI/N
  !     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
  !
   800  hitest = cuthi/N
  !
  !                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
  !
        DO j = i , nn , Incx
           IF ( ABS(Dx(j)).GE.hitest ) GOTO 400
           sum = sum + Dx(j)**2
        ENDDO
        DNRM2 = SQRT(sum)
        GOTO 99999
  !
   900  i = i + Incx
        IF ( i.LE.nn ) GOTO 100
  !
  !              END OF MAIN LOOP.
  !
  !              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
  !
        DNRM2 = xmax*SQRT(sum)
  99999 END
  !*==DSCAL.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
  !DECK DSCAL
        SUBROUTINE DSCAL(N,Da,Dx,Incx)
        IMPLICIT NONE
  !*--DSCAL1074
  !***BEGIN PROLOGUE  DSCAL
  !***PURPOSE  Multiply a vector by a constant.
  !***CATEGORY  D1A6
  !***TYPE      DOUBLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)
  !***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
  !***AUTHOR  Lawson, C. L., (JPL)
  !           Hanson, R. J., (SNLA)
  !           Kincaid, D. R., (U. of Texas)
  !           Krogh, F. T., (JPL)
  !***DESCRIPTION
  !
  !                B L A S  Subprogram
  !    Description of Parameters
  !
  !     --Input--
  !        N  number of elements in input vector(s)
  !       DA  double precision scale factor
  !       DX  double precision vector with N elements
  !     INCX  storage spacing between elements of DX
  !
  !     --Output--
  !       DX  double precision result (unchanged if N.LE.0)
  !
  !     Replace double precision DX by double precision DA*DX.
  !     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),
  !     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
  !
  !***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !                 Krogh, Basic linear algebra subprograms for Fortran
  !                 usage, Algorithm No. 539, Transactions on Mathematical
  !                 Software 5, 3 (September 1979), pp. 308-323.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900821  Modified to correct problem with a negative increment.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DSCAL
        DOUBLE PRECISION Da , Dx(*)
        INTEGER i , Incx , ix , m , mp1 , N
  !***FIRST EXECUTABLE STATEMENT  DSCAL
        IF ( N.LE.0 ) RETURN
        IF ( Incx.EQ.1 ) THEN
  !
  !     Code for increment equal to 1.
  !
  !     Clean-up loop so remaining vector length is a multiple of 5.
  !
           m = MOD(N,5)
           IF ( m.NE.0 ) THEN
              DO i = 1 , m
                 Dx(i) = Da*Dx(i)
              ENDDO
              IF ( N.LT.5 ) RETURN
           ENDIF
           mp1 = m + 1
           DO i = mp1 , N , 5
              Dx(i) = Da*Dx(i)
              Dx(i+1) = Da*Dx(i+1)
              Dx(i+2) = Da*Dx(i+2)
              Dx(i+3) = Da*Dx(i+3)
              Dx(i+4) = Da*Dx(i+4)
           ENDDO
        ELSE
  !
  !     Code for increment not equal to 1.
  !
           ix = 1
           IF ( Incx.LT.0 ) ix = (-N+1)*Incx + 1
           DO i = 1 , N
              Dx(ix) = Da*Dx(ix)
              ix = ix + Incx
           ENDDO
           RETURN
        ENDIF
        END
  !*==IDAMAX.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
  !DECK IDAMAX
        INTEGER FUNCTION IDAMAX(N,Dx,Incx)
        IMPLICIT NONE
  !*--IDAMAX1158
  !***BEGIN PROLOGUE  IDAMAX
  !***PURPOSE  Find the smallest index of that component of a vector
  !            having the maximum magnitude.
  !***CATEGORY  D1A2
  !***TYPE      DOUBLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C)
  !***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR
  !***AUTHOR  Lawson, C. L., (JPL)
  !           Hanson, R. J., (SNLA)
  !           Kincaid, D. R., (U. of Texas)
  !           Krogh, F. T., (JPL)
  !***DESCRIPTION
  !
  !                B L A S  Subprogram
  !    Description of Parameters
  !
  !     --Input--
  !        N  number of elements in input vector(s)
  !       DX  double precision vector with N elements
  !     INCX  storage spacing between elements of DX
  !
  !     --Output--
  !   IDAMAX  smallest index (zero if N .LE. 0)
  !
  !     Find smallest index of maximum magnitude of double precision DX.
  !     IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX+(I-1)*INCX)),
  !     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
  !
  !***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
  !                 Krogh, Basic linear algebra subprograms for Fortran
  !                 usage, Algorithm No. 539, Transactions on Mathematical
  !                 Software 5, 3 (September 1979), pp. 308-323.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900821  Modified to correct problem with a negative increment.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  IDAMAX
        DOUBLE PRECISION Dx(*) , dmax , xmag
        INTEGER i , Incx , ix , N
  !***FIRST EXECUTABLE STATEMENT  IDAMAX
        IDAMAX = 0
        IF ( N.LE.0 ) RETURN
        IDAMAX = 1
        IF ( N.EQ.1 ) RETURN
  !
        IF ( Incx.EQ.1 ) THEN
  !
  !     Code for increments equal to 1.
  !
           dmax = ABS(Dx(1))
           DO i = 2 , N
              xmag = ABS(Dx(i))
              IF ( xmag.GT.dmax ) THEN
                 IDAMAX = i
                 dmax = xmag
              ENDIF
           ENDDO
           GOTO 99999
        ENDIF
  !
  !     Code for increments not equal to 1.
  !
        ix = 1
        IF ( Incx.LT.0 ) ix = (-N+1)*Incx + 1
        dmax = ABS(Dx(ix))
        ix = ix + Incx
        DO i = 2 , N
           xmag = ABS(Dx(ix))
           IF ( xmag.GT.dmax ) THEN
              IDAMAX = i
              dmax = xmag
           ENDIF
           ix = ix + Incx
        ENDDO
        RETURN
  99999 END
  !*==XERRWD.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
  !DECK XERRWD
        SUBROUTINE XERRWD(Msg,Nmes,Nerr,Level,Ni,I1,I2,Nr,R1,R2)
        IMPLICIT NONE
  !*--XERRWD1243
  !***BEGIN PROLOGUE  XERRWD
  !***SUBSIDIARY
  !***PURPOSE  Write error message with values.
  !***CATEGORY  R3C
  !***TYPE      DOUBLE PRECISION (XERRWV-S, XERRWD-D)
  !***AUTHOR  Hindmarsh, Alan C., (LLNL)
  !***DESCRIPTION
  !
  !  Subroutines XERRWD, XSETF, XSETUN, and the function routine IXSAV,
  !  as given here, constitute a simplified version of the SLATEC error
  !  handling package.
  !
  !  All arguments are input arguments.
  !
  !  MSG    = The message (character array).
  !  NMES   = The length of MSG (number of characters).
  !  NERR   = The error number (not used).
  !  LEVEL  = The error level..
  !           0 or 1 means recoverable (control returns to caller).
  !           2 means fatal (run is aborted--see note below).
  !  NI     = Number of integers (0, 1, or 2) to be printed with message.
  !  I1,I2  = Integers to be printed, depending on NI.
  !  NR     = Number of reals (0, 1, or 2) to be printed with message.
  !  R1,R2  = Reals to be printed, depending on NR.
  !
  !  Note..  this routine is machine-dependent and specialized for use
  !  in limited context, in the following ways..
  !  1. The argument MSG is assumed to be of type CHARACTER, and
  !     the message is printed with a format of (1X,A).
  !  2. The message is assumed to take only one line.
  !     Multi-line messages are generated by repeated calls.
  !  3. If LEVEL = 2, control passes to the statement   STOP
  !     to abort the run.  This statement may be machine-dependent.
  !  4. R1 and R2 are assumed to be in double precision and are printed
  !     in D21.13 format.
  !
  !***ROUTINES CALLED  IXSAV
  !***REVISION HISTORY  (YYMMDD)
  !   920831  DATE WRITTEN
  !   921118  Replaced MFLGSV/LUNSAV by IXSAV. (ACH)
  !   930329  Modified prologue to SLATEC format. (FNF)
  !   930407  Changed MSG from CHARACTER*1 array to variable. (FNF)
  !   930922  Minor cosmetic change. (FNF)
  !***END PROLOGUE  XERRWD
  !
  !*Internal Notes:
  !
  ! For a different default logical unit number, IXSAV (or a subsidiary
  ! routine that it calls) will need to be modified.
  ! For a different run-abort command, change the statement following
  ! statement 100 at the end.
  !-----------------------------------------------------------------------
  ! Subroutines called by XERRWD.. None
  ! Function routine called by XERRWD.. IXSAV
  !-----------------------------------------------------------------------
  !**End
  !
  !  Declare arguments.
  !
        DOUBLE PRECISION R1 , R2
        INTEGER Nmes , Nerr , Level , Ni , I1 , I2 , Nr
        CHARACTER*(*) Msg
  !
  !  Declare local variables.
  !
        INTEGER lunit , mesflg
  !
  !  Get logical unit number and message print flag.
  !
  !***FIRST EXECUTABLE STATEMENT  XERRWD
        lunit = IXSAV(1,0,.FALSE.)
        mesflg = IXSAV(2,0,.FALSE.)
        IF ( mesflg.NE.0 ) THEN
  !
  !  Write the message.
  !
           WRITE (lunit,99001) Msg
  99001    FORMAT (1X,A)
           IF ( Ni.EQ.1 ) WRITE (lunit,99002) I1
  99002    FORMAT (6X,'In above message,  I1 =',I10)
           IF ( Ni.EQ.2 ) WRITE (lunit,99003) I1 , I2
  99003    FORMAT (6X,'In above message,  I1 =',I10,3X,'I2 =',I10)
           IF ( Nr.EQ.1 ) WRITE (lunit,99004) R1
  99004    FORMAT (6X,'In above message,  R1 =',D21.13)
           IF ( Nr.EQ.2 ) WRITE (lunit,99005) R1 , R2
  99005    FORMAT (6X,'In above,  R1 =',D21.13,3X,'R2 =',D21.13)
        ENDIF
  !
  !  Abort the run if LEVEL = 2.
  !
        IF ( Level.NE.2 ) RETURN
        STOP
  !----------------------- End of Subroutine XERRWD ----------------------
        END
  !*==XSETF.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
  !DECK XSETF
        SUBROUTINE XSETF(Mflag)
        IMPLICIT NONE
  !*--XSETF1342
  !***BEGIN PROLOGUE  XSETF
  !***PURPOSE  Reset the error print control flag.
  !***CATEGORY  R3A
  !***TYPE      ALL (XSETF-A)
  !***KEYWORDS  ERROR CONTROL
  !***AUTHOR  Hindmarsh, Alan C., (LLNL)
  !***DESCRIPTION
  !
  !   XSETF sets the error print control flag to MFLAG:
  !      MFLAG=1 means print all messages (the default).
  !      MFLAG=0 means no printing.
  !
  !***SEE ALSO  XERRWD, XERRWV
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  IXSAV
  !***REVISION HISTORY  (YYMMDD)
  !   921118  DATE WRITTEN
  !   930329  Added SLATEC format prologue. (FNF)
  !   930407  Corrected SEE ALSO section. (FNF)
  !   930922  Made user-callable, and other cosmetic changes. (FNF)
  !***END PROLOGUE  XSETF
  !
  ! Subroutines called by XSETF.. None
  ! Function routine called by XSETF.. IXSAV
  !-----------------------------------------------------------------------
  !**End
        INTEGER Mflag , junk
  !
  !***FIRST EXECUTABLE STATEMENT  XSETF
        IF ( Mflag.EQ.0 .OR. Mflag.EQ.1 ) junk = IXSAV(2,Mflag,.TRUE.)
  !----------------------- End of Subroutine XSETF -----------------------
        END
  !*==XSETUN.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
  !DECK XSETUN
        SUBROUTINE XSETUN(Lun)
        IMPLICIT NONE
  !*--XSETUN1379
  !***BEGIN PROLOGUE  XSETUN
  !***PURPOSE  Reset the logical unit number for error messages.
  !***CATEGORY  R3B
  !***TYPE      ALL (XSETUN-A)
  !***KEYWORDS  ERROR CONTROL
  !***DESCRIPTION
  !
  !   XSETUN sets the logical unit number for error messages to LUN.
  !
  !***AUTHOR  Hindmarsh, Alan C., (LLNL)
  !***SEE ALSO  XERRWD, XERRWV
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  IXSAV
  !***REVISION HISTORY  (YYMMDD)
  !   921118  DATE WRITTEN
  !   930329  Added SLATEC format prologue. (FNF)
  !   930407  Corrected SEE ALSO section. (FNF)
  !   930922  Made user-callable, and other cosmetic changes. (FNF)
  !***END PROLOGUE  XSETUN
  !
  ! Subroutines called by XSETUN.. None
  ! Function routine called by XSETUN.. IXSAV
  !-----------------------------------------------------------------------
  !**End
        INTEGER Lun , junk
  !
  !***FIRST EXECUTABLE STATEMENT  XSETUN
        IF ( Lun.GT.0 ) junk = IXSAV(1,Lun,.TRUE.)
  !----------------------- End of Subroutine XSETUN ----------------------
        END
  !*==IXSAV.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
  !DECK IXSAV
        INTEGER FUNCTION IXSAV(Ipar,Ivalue,Iset)
        IMPLICIT NONE
  !*--IXSAV1414
  !***BEGIN PROLOGUE  IXSAV
  !***SUBSIDIARY
  !***PURPOSE  Save and recall error message control parameters.
  !***CATEGORY  R3C
  !***TYPE      ALL (IXSAV-A)
  !***AUTHOR  Hindmarsh, Alan C., (LLNL)
  !***DESCRIPTION
  !
  !  IXSAV saves and recalls one of two error message parameters:
  !    LUNIT, the logical unit number to which messages are printed, and
  !    MESFLG, the message print flag.
  !  This is a modification of the SLATEC library routine J4SAVE.
  !
  !  Saved local variables..
  !   LUNIT  = Logical unit number for messages.  The default is obtained
  !            by a call to IUMACH (may be machine-dependent).
  !   MESFLG = Print control flag..
  !            1 means print all messages (the default).
  !            0 means no printing.
  !
  !  On input..
  !    IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG).
  !    IVALUE = The value to be set for the parameter, if ISET = .TRUE.
  !    ISET   = Logical flag to indicate whether to read or write.
  !             If ISET = .TRUE., the parameter will be given
  !             the value IVALUE.  If ISET = .FALSE., the parameter
  !             will be unchanged, and IVALUE is a dummy argument.
  !
  !  On return..
  !    IXSAV = The (old) value of the parameter.
  !
  !***SEE ALSO  XERRWD, XERRWV
  !***ROUTINES CALLED  IUMACH
  !***REVISION HISTORY  (YYMMDD)
  !   921118  DATE WRITTEN
  !   930329  Modified prologue to SLATEC format. (FNF)
  !   930915  Added IUMACH call to get default output unit.  (ACH)
  !   930922  Minor cosmetic changes. (FNF)
  !   010425  Type declaration for IUMACH added. (ACH)
  !***END PROLOGUE  IXSAV
  !
  ! Subroutines called by IXSAV.. None
  ! Function routine called by IXSAV.. IUMACH
  !-----------------------------------------------------------------------
  !**End
        LOGICAL Iset
        INTEGER Ipar , Ivalue
  !-----------------------------------------------------------------------
        INTEGER lunit , mesflg
  !-----------------------------------------------------------------------
  ! The following Fortran-77 declaration is to cause the values of the
  ! listed (local) variables to be saved between calls to this routine.
  !-----------------------------------------------------------------------
        SAVE lunit , mesflg
        DATA lunit/ - 1/ , mesflg/1/
  !
  !***FIRST EXECUTABLE STATEMENT  IXSAV
        IF ( Ipar.EQ.1 ) THEN
           IF ( lunit.EQ.-1 ) lunit = IUMACH()
           IXSAV = lunit
           IF ( Iset ) lunit = Ivalue
        ENDIF
  !
        IF ( Ipar.EQ.2 ) THEN
           IXSAV = mesflg
           IF ( Iset ) mesflg = Ivalue
        ENDIF
  !
  !----------------------- End of Function IXSAV -------------------------
        END
  !*==IUMACH.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
  !DECK IUMACH
        INTEGER FUNCTION IUMACH()
        IMPLICIT NONE
  !*--IUMACH1489
  !***BEGIN PROLOGUE  IUMACH
  !***PURPOSE  Provide standard output unit number.
  !***CATEGORY  R1
  !***TYPE      INTEGER (IUMACH-I)
  !***KEYWORDS  MACHINE CONSTANTS
  !***AUTHOR  Hindmarsh, Alan C., (LLNL)
  !***DESCRIPTION
  ! *Usage:
  !        INTEGER  LOUT, IUMACH
  !        LOUT = IUMACH()
  !
  ! *Function Return Values:
  !     LOUT : the standard logical unit for Fortran output.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   930915  DATE WRITTEN
  !   930922  Made user-callable, and other cosmetic changes. (FNF)
  !***END PROLOGUE  IUMACH
  !
  !*Internal Notes:
  !  The built-in value of 6 is standard on a wide range of Fortran
  !  systems.  This may be machine-dependent.
  !**End
  !***FIRST EXECUTABLE STATEMENT  IUMACH
        IUMACH = 6
  !
  !----------------------- End of Function IUMACH ------------------------
        END

end module ode_lsoda_aux2
