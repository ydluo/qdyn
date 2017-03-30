module ode_lsoda_aux1

  use ode_lsoda_aux2

  implicit none
  public

contains
  !*==DUMACH.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DUMACH
        DOUBLE PRECISION FUNCTION DUMACH()
        IMPLICIT NONE
  !*--DUMACH5
  !***BEGIN PROLOGUE  DUMACH
  !***PURPOSE  Compute the unit roundoff of the machine.
  !***CATEGORY  R1
  !***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D)
  !***KEYWORDS  MACHINE CONSTANTS
  !***AUTHOR  Hindmarsh, Alan C., (LLNL)
  !***DESCRIPTION
  ! *Usage:
  !        DOUBLE PRECISION  A, DUMACH
  !        A = DUMACH()
  !
  ! *Function Return Values:
  !     A : the unit roundoff of the machine.
  !
  ! *Description:
  !     The unit roundoff is defined as the smallest positive machine
  !     number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
  !     in a machine-independent manner.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  DUMSUM
  !***REVISION HISTORY  (YYYYMMDD)
  !   19930216  DATE WRITTEN
  !   19930818  Added SLATEC-format prologue.  (FNF)
  !   20030707  Added DUMSUM to force normal storage of COMP.  (ACH)
  !***END PROLOGUE  DUMACH
  !
        DOUBLE PRECISION u , comp
  !***FIRST EXECUTABLE STATEMENT  DUMACH
        u = 1.0D0
   100  u = u*0.5D0
        CALL DUMSUM(1.0D0,u,comp)
        IF ( comp.NE.1.0D0 ) GOTO 100
        DUMACH = u*2.0D0
  !----------------------- End of Function DUMACH ------------------------
        END
  !*==DUMSUM.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
        SUBROUTINE DUMSUM(A,B,C)
        IMPLICIT NONE
  !*--DUMSUM45
  !     Routine to force normal storing of A + B, for DUMACH.
        DOUBLE PRECISION A , B , C
        C = A + B
        END
  !*==DCFODE.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DCFODE
        SUBROUTINE DCFODE(Meth,Elco,Tesco)
        IMPLICIT NONE
  !*--DCFODE54
  !***BEGIN PROLOGUE  DCFODE
  !***SUBSIDIARY
  !***PURPOSE  Set ODE integrator coefficients.
  !***TYPE      DOUBLE PRECISION (SCFODE-S, DCFODE-D)
  !***AUTHOR  Hindmarsh, Alan C., (LLNL)
  !***DESCRIPTION
  !
  !  DCFODE is called by the integrator routine to set coefficients
  !  needed there.  The coefficients for the current method, as
  !  given by the value of METH, are set for all orders and saved.
  !  The maximum order assumed here is 12 if METH = 1 and 5 if METH = 2.
  !  (A smaller value of the maximum order is also allowed.)
  !  DCFODE is called once at the beginning of the problem,
  !  and is not called again unless and until METH is changed.
  !
  !  The ELCO array contains the basic method coefficients.
  !  The coefficients el(i), 1 .le. i .le. nq+1, for the method of
  !  order nq are stored in ELCO(i,nq).  They are given by a genetrating
  !  polynomial, i.e.,
  !      l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
  !  For the implicit Adams methods, l(x) is given by
  !      dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
  !  For the BDF methods, l(x) is given by
  !      l(x) = (x+1)*(x+2)* ... *(x+nq)/K,
  !  where         K = factorial(nq)*(1 + 1/2 + ... + 1/nq).
  !
  !  The TESCO array contains test constants used for the
  !  local error test and the selection of step size and/or order.
  !  At order nq, TESCO(k,nq) is used for the selection of step
  !  size at order nq - 1 if k = 1, at order nq if k = 2, and at order
  !  nq + 1 if k = 3.
  !
  !***SEE ALSO  DLSODE
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   791129  DATE WRITTEN
  !   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
  !   890503  Minor cosmetic changes.  (FNF)
  !   930809  Renamed to allow single/double precision versions. (ACH)
  !***END PROLOGUE  DCFODE
  !**End
        INTEGER Meth
        INTEGER i , ib , nq , nqm1 , nqp1
        DOUBLE PRECISION Elco , Tesco
        DOUBLE PRECISION agamq , fnq , fnqm1 , pc , pint , ragq , rqfac , &
       &                 rq1fac , tsign , xpin
        DIMENSION Elco(13,12) , Tesco(3,12)
        DIMENSION pc(12)
  !
  !***FIRST EXECUTABLE STATEMENT  DCFODE
        GOTO (100,200) , Meth
  !
   100  Elco(1,1) = 1.0D0
        Elco(2,1) = 1.0D0
        Tesco(1,1) = 0.0D0
        Tesco(2,1) = 2.0D0
        Tesco(1,2) = 1.0D0
        Tesco(3,12) = 0.0D0
        pc(1) = 1.0D0
        rqfac = 1.0D0
        DO nq = 2 , 12
  !-----------------------------------------------------------------------
  ! The PC array will contain the coefficients of the polynomial
  !     p(x) = (x+1)*(x+2)*...*(x+nq-1).
  ! Initially, p(x) = 1.
  !-----------------------------------------------------------------------
           rq1fac = rqfac
           rqfac = rqfac/nq
           nqm1 = nq - 1
           fnqm1 = nqm1
           nqp1 = nq + 1
  ! Form coefficients of p(x)*(x+nq-1). ----------------------------------
           pc(nq) = 0.0D0
           DO ib = 1 , nqm1
              i = nqp1 - ib
              pc(i) = pc(i-1) + fnqm1*pc(i)
           ENDDO
           pc(1) = fnqm1*pc(1)
  ! Compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
           pint = pc(1)
           xpin = pc(1)/2.0D0
           tsign = 1.0D0
           DO i = 2 , nq
              tsign = -tsign
              pint = pint + tsign*pc(i)/i
              xpin = xpin + tsign*pc(i)/(i+1)
           ENDDO
  ! Store coefficients in ELCO and TESCO. --------------------------------
           Elco(1,nq) = pint*rq1fac
           Elco(2,nq) = 1.0D0
           DO i = 2 , nq
              Elco(i+1,nq) = rq1fac*pc(i)/i
           ENDDO
           agamq = rqfac*xpin
           ragq = 1.0D0/agamq
           Tesco(2,nq) = ragq
           IF ( nq.LT.12 ) Tesco(1,nqp1) = ragq*rqfac/nqp1
           Tesco(3,nqm1) = ragq
        ENDDO
        RETURN
  !
   200  pc(1) = 1.0D0
        rq1fac = 1.0D0
        DO nq = 1 , 5
  !-----------------------------------------------------------------------
  ! The PC array will contain the coefficients of the polynomial
  !     p(x) = (x+1)*(x+2)*...*(x+nq).
  ! Initially, p(x) = 1.
  !-----------------------------------------------------------------------
           fnq = nq
           nqp1 = nq + 1
  ! Form coefficients of p(x)*(x+nq). ------------------------------------
           pc(nqp1) = 0.0D0
           DO ib = 1 , nq
              i = nq + 2 - ib
              pc(i) = pc(i-1) + fnq*pc(i)
           ENDDO
           pc(1) = fnq*pc(1)
  ! Store coefficients in ELCO and TESCO. --------------------------------
           DO i = 1 , nqp1
              Elco(i,nq) = pc(i)/pc(2)
           ENDDO
           Elco(2,nq) = 1.0D0
           Tesco(1,nq) = rq1fac
           Tesco(2,nq) = nqp1/Elco(1,nq)
           Tesco(3,nq) = (nq+2)/Elco(1,nq)
           rq1fac = rq1fac/fnq
        ENDDO
  !----------------------- END OF SUBROUTINE DCFODE ----------------------
        END
  !*==DINTDY.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DINTDY
        SUBROUTINE DINTDY(T,K,Yh,Nyh,Dky,Iflag)
        IMPLICIT NONE
  !*--DINTDY189
  !***BEGIN PROLOGUE  DINTDY
  !***SUBSIDIARY
  !***PURPOSE  Interpolate solution derivatives.
  !***TYPE      DOUBLE PRECISION (SINTDY-S, DINTDY-D)
  !***AUTHOR  Hindmarsh, Alan C., (LLNL)
  !***DESCRIPTION
  !
  !  DINTDY computes interpolated values of the K-th derivative of the
  !  dependent variable vector y, and stores it in DKY.  This routine
  !  is called within the package with K = 0 and T = TOUT, but may
  !  also be called by the user for any K up to the current order.
  !  (See detailed instructions in the usage documentation.)
  !
  !  The computed values in DKY are gotten by interpolation using the
  !  Nordsieck history array YH.  This array corresponds uniquely to a
  !  vector-valued polynomial of degree NQCUR or less, and DKY is set
  !  to the K-th derivative of this polynomial at T.
  !  The formula for DKY is:
  !               q
  !   DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)
  !              j=K
  !  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR.
  !  The quantities  nq = NQCUR, l = nq+1, N = NEQ, tn, and h are
  !  communicated by COMMON.  The above sum is done in reverse order.
  !  IFLAG is returned negative if either K or T is out of bounds.
  !
  !***SEE ALSO  DLSODE
  !***ROUTINES CALLED  XERRWD
  !***COMMON BLOCKS    DLS001
  !***REVISION HISTORY  (YYMMDD)
  !   791129  DATE WRITTEN
  !   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
  !   890503  Minor cosmetic changes.  (FNF)
  !   930809  Renamed to allow single/double precision versions. (ACH)
  !   010418  Reduced size of Common block /DLS001/. (ACH)
  !   031105  Restored 'own' variables to Common block /DLS001/, to
  !           enable interrupt/restart feature. (ACH)
  !   050427  Corrected roundoff decrement in TP. (ACH)
  !***END PROLOGUE  DINTDY
  !**End
        INTEGER K , Nyh , Iflag
        DOUBLE PRECISION T , Yh , Dky
        DIMENSION Yh(Nyh,*) , Dky(*)
        INTEGER IOWnd , IOWns , ICF , IERpj , IERsl , JCUr , JSTart ,     &
       &        KFLag , L , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,     &
       &        METh , MITer , MAXord , MAXcor , MSBp , MXNcf , N , NQ ,  &
       &        NST , NFE , NJE , NQU
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , IOWnd(6) , IOWns(6) , ICF ,    &
       &                IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , &
       &                LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,&
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        INTEGER i , ic , j , jb , jb2 , jj , jj1 , jp1
        DOUBLE PRECISION c , r , s , tp
        CHARACTER*80 msg
  !
  !***FIRST EXECUTABLE STATEMENT  DINTDY
        Iflag = 0
        IF ( K.LT.0 .OR. K.GT.NQ ) THEN
  !
           msg = 'DINTDY-  K (=I1) illegal      '
           CALL XERRWD(msg,30,51,0,1,K,0,0,0.0D0,0.0D0)
           Iflag = -1
           RETURN
        ELSE
           tp = TN - HU - 100.0D0*UROund*SIGN(ABS(TN)+ABS(HU),HU)
           IF ( (T-tp)*(T-TN).GT.0.0D0 ) THEN
              msg = 'DINTDY-  T (=R1) illegal      '
              CALL XERRWD(msg,30,52,0,0,0,0,1,T,0.0D0)
              msg =                                                       &
       &    '      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
              CALL XERRWD(msg,60,52,0,0,0,0,2,tp,TN)
              Iflag = -2
              GOTO 99999
           ELSE
  !
              s = (T-TN)/H
              ic = 1
              IF ( K.NE.0 ) THEN
                 jj1 = L - K
                 DO jj = jj1 , NQ
                    ic = ic*jj
                 ENDDO
              ENDIF
              c = ic
              DO i = 1 , N
                 Dky(i) = c*Yh(i,L)
              ENDDO
              IF ( K.NE.NQ ) THEN
                 jb2 = NQ - K
                 DO jb = 1 , jb2
                    j = NQ - jb
                    jp1 = j + 1
                    ic = 1
                    IF ( K.NE.0 ) THEN
                       jj1 = jp1 - K
                       DO jj = jj1 , j
                          ic = ic*jj
                       ENDDO
                    ENDIF
                    c = ic
                    DO i = 1 , N
                       Dky(i) = c*Yh(i,jp1) + s*Dky(i)
                    ENDDO
                 ENDDO
                 IF ( K.EQ.0 ) RETURN
              ENDIF
           ENDIF
        ENDIF
        r = H**(-K)
        DO i = 1 , N
           Dky(i) = r*Dky(i)
        ENDDO
        RETURN
  !----------------------- END OF SUBROUTINE DINTDY ----------------------
  99999 END
  !*==DPREPJ.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DPREPJ
        SUBROUTINE DPREPJ(Neq,Y,Yh,Nyh,Ewt,Ftem,Savf,Wm,Iwm,F,JAC)
        IMPLICIT NONE
  !*--DPREPJ313
  !***BEGIN PROLOGUE  DPREPJ
  !***SUBSIDIARY
  !***PURPOSE  Compute and process Newton iteration matrix.
  !***TYPE      DOUBLE PRECISION (SPREPJ-S, DPREPJ-D)
  !***AUTHOR  Hindmarsh, Alan C., (LLNL)
  !***DESCRIPTION
  !
  !  DPREPJ is called by DSTODE to compute and process the matrix
  !  P = I - h*el(1)*J , where J is an approximation to the Jacobian.
  !  Here J is computed by the user-supplied routine JAC if
  !  MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5.
  !  If MITER = 3, a diagonal approximation to J is used.
  !  J is stored in WM and replaced by P.  If MITER .ne. 3, P is then
  !  subjected to LU decomposition in preparation for later solution
  !  of linear systems with P as coefficient matrix.  This is done
  !  by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
  !
  !  In addition to variables described in DSTODE and DLSODE prologues,
  !  communication with DPREPJ uses the following:
  !  Y     = array containing predicted values on entry.
  !  FTEM  = work array of length N (ACOR in DSTODE).
  !  SAVF  = array containing f evaluated at predicted y.
  !  WM    = real work space for matrices.  On output it contains the
  !          inverse diagonal matrix if MITER = 3 and the LU decomposition
  !          of P if MITER is 1, 2 , 4, or 5.
  !          Storage of matrix elements starts at WM(3).
  !          WM also contains the following matrix-related data:
  !          WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
  !          WM(2) = H*EL0, saved for later use if MITER = 3.
  !  IWM   = integer work space containing pivot information, starting at
  !          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
  !          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
  !  EL0   = EL(1) (input).
  !  IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
  !          P matrix found to be singular.
  !  JCUR  = output flag = 1 to indicate that the Jacobian matrix
  !          (or approximation) is now current.
  !  This routine also uses the COMMON variables EL0, H, TN, UROUND,
  !  MITER, N, NFE, and NJE.
  !
  !***SEE ALSO  DLSODE
  !***ROUTINES CALLED  DGBFA, DGEFA, DVNORM
  !***COMMON BLOCKS    DLS001
  !***REVISION HISTORY  (YYMMDD)
  !   791129  DATE WRITTEN
  !   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
  !   890504  Minor cosmetic changes.  (FNF)
  !   930809  Renamed to allow single/double precision versions. (ACH)
  !   010418  Reduced size of Common block /DLS001/. (ACH)
  !   031105  Restored 'own' variables to Common block /DLS001/, to
  !           enable interrupt/restart feature. (ACH)
  !***END PROLOGUE  DPREPJ
  !**End
        EXTERNAL F , JAC
        INTEGER Neq , Nyh , Iwm
        DOUBLE PRECISION Y , Yh , Ewt , Ftem , Savf , Wm
        DIMENSION Neq(*) , Y(*) , Yh(Nyh,*) , Ewt(*) , Ftem(*) , Savf(*) ,&
       &          Wm(*) , Iwm(*)
        INTEGER IOWnd , IOWns , ICF , IERpj , IERsl , JCUr , JSTart ,     &
       &        KFLag , L , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,     &
       &        METh , MITer , MAXord , MAXcor , MSBp , MXNcf , N , NQ ,  &
       &        NST , NFE , NJE , NQU
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , IOWnd(6) , IOWns(6) , ICF ,    &
       &                IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , &
       &                LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,&
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        INTEGER i , i1 , i2 , ier , ii , j , j1 , jj , lenp , mba ,       &
       &        mband , meb1 , meband , ml , ml3 , mu , np1
        DOUBLE PRECISION con , di , fac , hl0 , r , r0 , srur , yi , yj , &
       &                 yjj
  !
  !***FIRST EXECUTABLE STATEMENT  DPREPJ
        NJE = NJE + 1
        IERpj = 0
        JCUr = 1
        hl0 = H*EL0
        GOTO (100,200,400,600,700) , MITer
  ! If MITER = 1, call JAC and multiply by scalar. -----------------------
   100  lenp = N*N
        DO i = 1 , lenp
           Wm(i+2) = 0.0D0
        ENDDO
        CALL JAC(Neq,TN,Y,0,0,Wm(3),N)
        con = -hl0
        DO i = 1 , lenp
           Wm(i+2) = Wm(i+2)*con
        ENDDO
        GOTO 300
  ! If MITER = 2, make N calls to F to approximate J. --------------------
   200  fac = DVNORM(N,Savf,Ewt)
        r0 = 1000.0D0*ABS(H)*UROund*N*fac
        IF ( r0.EQ.0.0D0 ) r0 = 1.0D0
        srur = Wm(1)
        j1 = 2
        DO j = 1 , N
           yj = Y(j)
           r = MAX(srur*ABS(yj),r0/Ewt(j))
           Y(j) = Y(j) + r
           fac = -hl0/r
           CALL F(Neq,TN,Y,Ftem)
           DO i = 1 , N
              Wm(i+j1) = (Ftem(i)-Savf(i))*fac
           ENDDO
           Y(j) = yj
           j1 = j1 + N
        ENDDO
        NFE = NFE + N
  ! Add identity matrix. -------------------------------------------------
   300  j = 3
        np1 = N + 1
        DO i = 1 , N
           Wm(j) = Wm(j) + 1.0D0
           j = j + np1
        ENDDO
  ! Do LU decomposition on P. --------------------------------------------
        CALL DGEFA(Wm(3),N,N,Iwm(21),ier)
        IF ( ier.NE.0 ) IERpj = 1
        RETURN
  ! If MITER = 3, construct a diagonal approximation to J and P. ---------
   400  Wm(2) = hl0
        r = EL0*0.1D0
        DO i = 1 , N
           Y(i) = Y(i) + r*(H*Savf(i)-Yh(i,2))
        ENDDO
        CALL F(Neq,TN,Y,Wm(3))
        NFE = NFE + 1
        DO i = 1 , N
           r0 = H*Savf(i) - Yh(i,2)
           di = 0.1D0*r0 - H*(Wm(i+2)-Savf(i))
           Wm(i+2) = 1.0D0
           IF ( ABS(r0).GE.UROund/Ewt(i) ) THEN
              IF ( ABS(di).EQ.0.0D0 ) GOTO 500
              Wm(i+2) = 0.1D0*r0/di
           ENDIF
        ENDDO
        RETURN
   500  IERpj = 1
        RETURN
  ! If MITER = 4, call JAC and multiply by scalar. -----------------------
   600  ml = Iwm(1)
        mu = Iwm(2)
        ml3 = ml + 3
        mband = ml + mu + 1
        meband = mband + ml
        lenp = meband*N
        DO i = 1 , lenp
           Wm(i+2) = 0.0D0
        ENDDO
        CALL JAC(Neq,TN,Y,ml,mu,Wm(ml3),meband)
        con = -hl0
        DO i = 1 , lenp
           Wm(i+2) = Wm(i+2)*con
        ENDDO
        GOTO 800
  ! If MITER = 5, make MBAND calls to F to approximate J. ----------------
   700  ml = Iwm(1)
        mu = Iwm(2)
        mband = ml + mu + 1
        mba = MIN(mband,N)
        meband = mband + ml
        meb1 = meband - 1
        srur = Wm(1)
        fac = DVNORM(N,Savf,Ewt)
        r0 = 1000.0D0*ABS(H)*UROund*N*fac
        IF ( r0.EQ.0.0D0 ) r0 = 1.0D0
        DO j = 1 , mba
           DO i = j , N , mband
              yi = Y(i)
              r = MAX(srur*ABS(yi),r0/Ewt(i))
              Y(i) = Y(i) + r
           ENDDO
           CALL F(Neq,TN,Y,Ftem)
           DO jj = j , N , mband
              Y(jj) = Yh(jj,1)
              yjj = Y(jj)
              r = MAX(srur*ABS(yjj),r0/Ewt(jj))
              fac = -hl0/r
              i1 = MAX(jj-mu,1)
              i2 = MIN(jj+ml,N)
              ii = jj*meb1 - ml + 2
              DO i = i1 , i2
                 Wm(ii+i) = (Ftem(i)-Savf(i))*fac
              ENDDO
           ENDDO
        ENDDO
        NFE = NFE + mba
  ! Add identity matrix. -------------------------------------------------
   800  ii = mband + 2
        DO i = 1 , N
           Wm(ii) = Wm(ii) + 1.0D0
           ii = ii + meband
        ENDDO
  ! Do LU decomposition of P. --------------------------------------------
        CALL DGBFA(Wm(3),meband,N,ml,mu,Iwm(21),ier)
        IF ( ier.NE.0 ) IERpj = 1
  !----------------------- END OF SUBROUTINE DPREPJ ----------------------
        END
  !*==DSOLSY.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DSOLSY
        SUBROUTINE DSOLSY(Wm,Iwm,X,Tem)
        IMPLICIT NONE
  !*--DSOLSY519
  !***BEGIN PROLOGUE  DSOLSY
  !***SUBSIDIARY
  !***PURPOSE  ODEPACK linear system solver.
  !***TYPE      DOUBLE PRECISION (SSOLSY-S, DSOLSY-D)
  !***AUTHOR  Hindmarsh, Alan C., (LLNL)
  !***DESCRIPTION
  !
  !  This routine manages the solution of the linear system arising from
  !  a chord iteration.  It is called if MITER .ne. 0.
  !  If MITER is 1 or 2, it calls DGESL to accomplish this.
  !  If MITER = 3 it updates the coefficient h*EL0 in the diagonal
  !  matrix, and then computes the solution.
  !  If MITER is 4 or 5, it calls DGBSL.
  !  Communication with DSOLSY uses the following variables:
  !  WM    = real work space containing the inverse diagonal matrix if
  !          MITER = 3 and the LU decomposition of the matrix otherwise.
  !          Storage of matrix elements starts at WM(3).
  !          WM also contains the following matrix-related data:
  !          WM(1) = SQRT(UROUND) (not used here),
  !          WM(2) = HL0, the previous value of h*EL0, used if MITER = 3.
  !  IWM   = integer work space containing pivot information, starting at
  !          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
  !          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
  !  X     = the right-hand side vector on input, and the solution vector
  !          on output, of length N.
  !  TEM   = vector of work space of length N, not used in this version.
  !  IERSL = output flag (in COMMON).  IERSL = 0 if no trouble occurred.
  !          IERSL = 1 if a singular matrix arose with MITER = 3.
  !  This routine also uses the COMMON variables EL0, H, MITER, and N.
  !
  !***SEE ALSO  DLSODE
  !***ROUTINES CALLED  DGBSL, DGESL
  !***COMMON BLOCKS    DLS001
  !***REVISION HISTORY  (YYMMDD)
  !   791129  DATE WRITTEN
  !   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
  !   890503  Minor cosmetic changes.  (FNF)
  !   930809  Renamed to allow single/double precision versions. (ACH)
  !   010418  Reduced size of Common block /DLS001/. (ACH)
  !   031105  Restored 'own' variables to Common block /DLS001/, to
  !           enable interrupt/restart feature. (ACH)
  !***END PROLOGUE  DSOLSY
  !**End
        INTEGER Iwm
        DOUBLE PRECISION Wm , X , Tem
        DIMENSION Wm(*) , Iwm(*) , X(*) , Tem(*)
        INTEGER IOWnd , IOWns , ICF , IERpj , IERsl , JCUr , JSTart ,     &
       &        KFLag , L , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,     &
       &        METh , MITer , MAXord , MAXcor , MSBp , MXNcf , N , NQ ,  &
       &        NST , NFE , NJE , NQU
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , IOWnd(6) , IOWns(6) , ICF ,    &
       &                IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , &
       &                LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,&
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        INTEGER i , meband , ml , mu
        DOUBLE PRECISION di , hl0 , phl0 , r
  !
  !***FIRST EXECUTABLE STATEMENT  DSOLSY
        IERsl = 0
        GOTO (100,100,200,400,400) , MITer
   100  CALL DGESL(Wm(3),N,N,Iwm(21),X,0)
        RETURN
  !
   200  phl0 = Wm(2)
        hl0 = H*EL0
        Wm(2) = hl0
        IF ( hl0.NE.phl0 ) THEN
           r = hl0/phl0
           DO i = 1 , N
              di = 1.0D0 - r*(1.0D0-1.0D0/Wm(i+2))
              IF ( ABS(di).EQ.0.0D0 ) GOTO 300
              Wm(i+2) = 1.0D0/di
           ENDDO
        ENDIF
        DO i = 1 , N
           X(i) = Wm(i+2)*X(i)
        ENDDO
        RETURN
   300  IERsl = 1
        RETURN
  !
   400  ml = Iwm(1)
        mu = Iwm(2)
        meband = 2*ml + mu + 1
        CALL DGBSL(Wm(3),meband,N,ml,mu,Iwm(21),X,0)
  !----------------------- END OF SUBROUTINE DSOLSY ----------------------
        END
  !*==DSRCOM.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DSRCOM
        SUBROUTINE DSRCOM(Rsav,Isav,Job)
        IMPLICIT NONE
  !*--DSRCOM615
  !***BEGIN PROLOGUE  DSRCOM
  !***SUBSIDIARY
  !***PURPOSE  Save/restore ODEPACK COMMON blocks.
  !***TYPE      DOUBLE PRECISION (SSRCOM-S, DSRCOM-D)
  !***AUTHOR  Hindmarsh, Alan C., (LLNL)
  !***DESCRIPTION
  !
  !  This routine saves or restores (depending on JOB) the contents of
  !  the COMMON block DLS001, which is used internally
  !  by one or more ODEPACK solvers.
  !
  !  RSAV = real array of length 218 or more.
  !  ISAV = integer array of length 37 or more.
  !  JOB  = flag indicating to save or restore the COMMON blocks:
  !         JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV)
  !         JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV)
  !         A call with JOB = 2 presumes a prior call with JOB = 1.
  !
  !***SEE ALSO  DLSODE
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    DLS001
  !***REVISION HISTORY  (YYMMDD)
  !   791129  DATE WRITTEN
  !   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
  !   890503  Minor cosmetic changes.  (FNF)
  !   921116  Deleted treatment of block /EH0001/.  (ACH)
  !   930801  Reduced Common block length by 2.  (ACH)
  !   930809  Renamed to allow single/double precision versions. (ACH)
  !   010418  Reduced Common block length by 209+12. (ACH)
  !   031105  Restored 'own' variables to Common block /DLS001/, to
  !           enable interrupt/restart feature. (ACH)
  !   031112  Added SAVE statement for data-loaded constants.
  !***END PROLOGUE  DSRCOM
  !**End
        INTEGER Isav , Job
        INTEGER ILS
        INTEGER i , lenils , lenrls
        DOUBLE PRECISION Rsav , RLS
        DIMENSION Rsav(*) , Isav(*)
        SAVE lenrls , lenils
        COMMON /DLS001/ RLS(218) , ILS(37)
        DATA lenrls/218/ , lenils/37/
  !
  !***FIRST EXECUTABLE STATEMENT  DSRCOM
        IF ( Job.EQ.2 ) THEN
  !
           DO i = 1 , lenrls
              RLS(i) = Rsav(i)
           ENDDO
           DO i = 1 , lenils
              ILS(i) = Isav(i)
           ENDDO
           GOTO 99999
        ENDIF
  !
        DO i = 1 , lenrls
           Rsav(i) = RLS(i)
        ENDDO
        DO i = 1 , lenils
           Isav(i) = ILS(i)
        ENDDO
        RETURN
  !----------------------- END OF SUBROUTINE DSRCOM ----------------------
  99999 END
  !*==DSTODE.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DSTODE
        SUBROUTINE DSTODE(Neq,Y,Yh,Nyh,Yh1,Ewt,Savf,Acor,Wm,Iwm,F,JAC,    &
       &                  PJAC,SLVS)
        IMPLICIT NONE
  !*--DSTODE685
  !*** Start of declarations inserted by SPAG
        !INTEGER JAC
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DSTODE
  !***SUBSIDIARY
  !***PURPOSE  Performs one step of an ODEPACK integration.
  !***TYPE      DOUBLE PRECISION (SSTODE-S, DSTODE-D)
  !***AUTHOR  Hindmarsh, Alan C., (LLNL)
  !***DESCRIPTION
  !
  !  DSTODE performs one step of the integration of an initial value
  !  problem for a system of ordinary differential equations.
  !  Note:  DSTODE is independent of the value of the iteration method
  !  indicator MITER, when this is .ne. 0, and hence is independent
  !  of the type of chord method used, or the Jacobian structure.
  !  Communication with DSTODE is done with the following variables:
  !
  !  NEQ    = integer array containing problem size in NEQ(1), and
  !           passed as the NEQ argument in all calls to F and JAC.
  !  Y      = an array of length .ge. N used as the Y argument in
  !           all calls to F and JAC.
  !  YH     = an NYH by LMAX array containing the dependent variables
  !           and their approximate scaled derivatives, where
  !           LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
  !           j-th derivative of y(i), scaled by h**j/factorial(j)
  !           (j = 0,1,...,NQ).  on entry for the first step, the first
  !           two columns of YH must be set from the initial values.
  !  NYH    = a constant integer .ge. N, the first dimension of YH.
  !  YH1    = a one-dimensional array occupying the same space as YH.
  !  EWT    = an array of length N containing multiplicative weights
  !           for local error measurements.  Local errors in Y(i) are
  !           compared to 1.0/EWT(i) in various error tests.
  !  SAVF   = an array of working storage, of length N.
  !           Also used for input of YH(*,MAXORD+2) when JSTART = -1
  !           and MAXORD .lt. the current order NQ.
  !  ACOR   = a work array of length N, used for the accumulated
  !           corrections.  On a successful return, ACOR(i) contains
  !           the estimated one-step local error in Y(i).
  !  WM,IWM = real and integer work arrays associated with matrix
  !           operations in chord iteration (MITER .ne. 0).
  !  PJAC   = name of routine to evaluate and preprocess Jacobian matrix
  !           and P = I - h*el0*JAC, if a chord method is being used.
  !  SLVS   = name of routine to solve linear system in chord iteration.
  !  CCMAX  = maximum relative change in h*el0 before PJAC is called.
  !  H      = the step size to be attempted on the next step.
  !           H is altered by the error control algorithm during the
  !           problem.  H can be either positive or negative, but its
  !           sign must remain constant throughout the problem.
  !  HMIN   = the minimum absolute value of the step size h to be used.
  !  HMXI   = inverse of the maximum absolute value of h to be used.
  !           HMXI = 0.0 is allowed and corresponds to an infinite hmax.
  !           HMIN and HMXI may be changed at any time, but will not
  !           take effect until the next change of h is considered.
  !  TN     = the independent variable. TN is updated on each step taken.
  !  JSTART = an integer used for input only, with the following
  !           values and meanings:
  !                0  perform the first step.
  !            .gt.0  take a new step continuing from the last.
  !               -1  take the next step with a new value of H, MAXORD,
  !                     N, METH, MITER, and/or matrix parameters.
  !               -2  take the next step with a new value of H,
  !                     but with other inputs unchanged.
  !           On return, JSTART is set to 1 to facilitate continuation.
  !  KFLAG  = a completion code with the following meanings:
  !                0  the step was succesful.
  !               -1  the requested error could not be achieved.
  !               -2  corrector convergence could not be achieved.
  !               -3  fatal error in PJAC or SLVS.
  !           A return with KFLAG = -1 or -2 means either
  !           abs(H) = HMIN or 10 consecutive failures occurred.
  !           On a return with KFLAG negative, the values of TN and
  !           the YH array are as of the beginning of the last
  !           step, and H is the last step size attempted.
  !  MAXORD = the maximum order of integration method to be allowed.
  !  MAXCOR = the maximum number of corrector iterations allowed.
  !  MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0).
  !  MXNCF  = maximum number of convergence failures allowed.
  !  METH/MITER = the method flags.  See description in driver.
  !  N      = the number of first-order differential equations.
  !  The values of CCMAX, H, HMIN, HMXI, TN, JSTART, KFLAG, MAXORD,
  !  MAXCOR, MSBP, MXNCF, METH, MITER, and N are communicated via COMMON.
  !
  !***SEE ALSO  DLSODE
  !***ROUTINES CALLED  DCFODE, DVNORM
  !***COMMON BLOCKS    DLS001
  !***REVISION HISTORY  (YYMMDD)
  !   791129  DATE WRITTEN
  !   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
  !   890503  Minor cosmetic changes.  (FNF)
  !   930809  Renamed to allow single/double precision versions. (ACH)
  !   010418  Reduced size of Common block /DLS001/. (ACH)
  !   031105  Restored 'own' variables to Common block /DLS001/, to
  !           enable interrupt/restart feature. (ACH)
  !***END PROLOGUE  DSTODE
  !**End
        EXTERNAL F , JAC , PJAC , SLVS
        INTEGER Neq , Nyh , Iwm
        DOUBLE PRECISION Y , Yh , Yh1 , Ewt , Savf , Acor , Wm
        DIMENSION Neq(*) , Y(*) , Yh(Nyh,*) , Yh1(*) , Ewt(*) , Savf(*) , &
       &          Acor(*) , Wm(*) , Iwm(*)
        INTEGER IOWnd , IALth , IPUp , LMAx , MEO , NQNyh , NSLp , ICF ,  &
       &        IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , LEWt ,  &
       &        LACor , LSAvf , LWM , LIWm , METh , MITer , MAXord ,      &
       &        MAXcor , MSBp , MXNcf , N , NQ , NST , NFE , NJE , NQU
        INTEGER i , i1 , iredo , iret , j , jb , m , ncf , newq
        DOUBLE PRECISION CONit , CRAte , EL , ELCo , HOLd , RMAx , TESco ,&
       &                 CCMax , EL0 , H , HMIn , HMXi , HU , RC , TN ,   &
       &                 UROund
        DOUBLE PRECISION dcon , ddn , del , delp , dsm , dup , exdn ,     &
       &                 exsm , exup , r , rh , rhdn , rhsm , rhup , told
        COMMON /DLS001/ CONit , CRAte , EL(13) , ELCo(13,12) , HOLd ,     &
       &                RMAx , TESco(3,12) , CCMax , EL0 , H , HMIn ,     &
       &                HMXi , HU , RC , TN , UROund , IOWnd(6) , IALth , &
       &                IPUp , LMAx , MEO , NQNyh , NSLp , ICF , IERpj ,  &
       &                IERsl , JCUr , JSTart , KFLag , L , LYH , LEWt ,  &
       &                LACor , LSAvf , LWM , LIWm , METh , MITer ,       &
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
  !
  !***FIRST EXECUTABLE STATEMENT  DSTODE
        KFLag = 0
        told = TN
        ncf = 0
        IERpj = 0
        IERsl = 0
        JCUr = 0
        ICF = 0
        delp = 0.0D0
        IF ( JSTart.GT.0 ) GOTO 500
        IF ( JSTart.EQ.-1 ) THEN
  !-----------------------------------------------------------------------
  ! The following block handles preliminaries needed when JSTART = -1.
  ! IPUP is set to MITER to force a matrix update.
  ! If an order increase is about to be considered (IALTH = 1),
  ! IALTH is reset to 2 to postpone consideration one more step.
  ! If the caller has changed METH, DCFODE is called to reset
  ! the coefficients of the method.
  ! If the caller has changed MAXORD to a value less than the current
  ! order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
  ! If H is to be changed, YH must be rescaled.
  ! If H or METH is being changed, IALTH is reset to L = NQ + 1
  ! to prevent further changes in H for that many steps.
  !-----------------------------------------------------------------------
           IPUp = MITer
           LMAx = MAXord + 1
           IF ( IALth.EQ.1 ) IALth = 2
           IF ( METh.NE.MEO ) THEN
              CALL DCFODE(METh,ELCo,TESco)
              MEO = METh
              IF ( NQ.LE.MAXord ) THEN
                 IALth = L
                 iret = 1
                 GOTO 100
              ENDIF
           ELSEIF ( NQ.LE.MAXord ) THEN
              GOTO 200
           ENDIF
           NQ = MAXord
           L = LMAx
           DO i = 1 , L
              EL(i) = ELCo(i,NQ)
           ENDDO
           NQNyh = NQ*Nyh
           RC = RC*EL(1)/EL0
           EL0 = EL(1)
           CONit = 0.5D0/(NQ+2)
           ddn = DVNORM(N,Savf,Ewt)/TESco(1,L)
           exdn = 1.0D0/L
           rhdn = 1.0D0/(1.3D0*ddn**exdn+0.0000013D0)
           rh = MIN(rhdn,1.0D0)
           iredo = 3
           IF ( H.EQ.HOLd ) GOTO 300
           rh = MIN(rh,ABS(H/HOLd))
           H = HOLd
           GOTO 400
        ELSE
           IF ( JSTart.EQ.-2 ) GOTO 200
  !-----------------------------------------------------------------------
  ! On the first call, the order is set to 1, and other variables are
  ! initialized.  RMAX is the maximum ratio by which H can be increased
  ! in a single step.  It is initially 1.E4 to compensate for the small
  ! initial H, but then is normally equal to 10.  If a failure
  ! occurs (in corrector convergence or error test), RMAX is set to 2
  ! for the next increase.
  !-----------------------------------------------------------------------
           LMAx = MAXord + 1
           NQ = 1
           L = 2
           IALth = 2
           RMAx = 10000.0D0
           RC = 0.0D0
           EL0 = 1.0D0
           CRAte = 0.7D0
           HOLd = H
           MEO = METh
           NSLp = 0
           IPUp = MITer
           iret = 3
  !-----------------------------------------------------------------------
  ! DCFODE is called to get all the integration coefficients for the
  ! current METH.  Then the EL vector and related constants are reset
  ! whenever the order NQ is changed, or at the start of the problem.
  !-----------------------------------------------------------------------
           CALL DCFODE(METh,ELCo,TESco)
        ENDIF
   100  DO i = 1 , L
           EL(i) = ELCo(i,NQ)
        ENDDO
        NQNyh = NQ*Nyh
        RC = RC*EL(1)/EL0
        EL0 = EL(1)
        CONit = 0.5D0/(NQ+2)
        GOTO (200,300,500) , iret
  !-----------------------------------------------------------------------
  ! If H is being changed, the H ratio RH is checked against
  ! RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
  ! L = NQ + 1 to prevent a change of H for that many steps, unless
  ! forced by a convergence or error test failure.
  !-----------------------------------------------------------------------
   200  IF ( H.EQ.HOLd ) GOTO 500
        rh = H/HOLd
        H = HOLd
        iredo = 3
        GOTO 400
   300  rh = MAX(rh,HMIn/ABS(H))
   400  rh = MIN(rh,RMAx)
        rh = rh/MAX(1.0D0,ABS(H)*HMXi*rh)
        r = 1.0D0
        DO j = 2 , L
           r = r*rh
           DO i = 1 , N
              Yh(i,j) = Yh(i,j)*r
           ENDDO
        ENDDO
        H = H*rh
        RC = RC*rh
        IALth = L
        IF ( iredo.EQ.0 ) THEN
           RMAx = 10.0D0
           GOTO 1300
        ENDIF
  !-----------------------------------------------------------------------
  ! This section computes the predicted values by effectively
  ! multiplying the YH array by the Pascal Triangle matrix.
  ! RC is the ratio of new to old values of the coefficient  H*EL(1).
  ! When RC differs from 1 by more than CCMAX, IPUP is set to MITER
  ! to force PJAC to be called, if a Jacobian is involved.
  ! In any case, PJAC is called at least every MSBP steps.
  !-----------------------------------------------------------------------
   500  IF ( ABS(RC-1.0D0).GT.CCMax ) IPUp = MITer
        IF ( NST.GE.NSLp+MSBp ) IPUp = MITer
        TN = TN + H
        i1 = NQNyh + 1
        DO jb = 1 , NQ
           i1 = i1 - Nyh
  !dir$ ivdep
           DO i = i1 , NQNyh
              Yh1(i) = Yh1(i) + Yh1(i+Nyh)
           ENDDO
        ENDDO
  !-----------------------------------------------------------------------
  ! Up to MAXCOR corrector iterations are taken.  A convergence test is
  ! made on the R.M.S. norm of each correction, weighted by the error
  ! weight vector EWT.  The sum of the corrections is accumulated in the
  ! vector ACOR(i).  The YH array is not altered in the corrector loop.
  !-----------------------------------------------------------------------
   600  m = 0
        DO i = 1 , N
           Y(i) = Yh(i,1)
        ENDDO
        CALL F(Neq,TN,Y,Savf)
        NFE = NFE + 1
        IF ( IPUp.GT.0 ) THEN
  !-----------------------------------------------------------------------
  ! If indicated, the matrix P = I - h*el(1)*J is reevaluated and
  ! preprocessed before starting the corrector iteration.  IPUP is set
  ! to 0 as an indicator that this has been done.
  !-----------------------------------------------------------------------
           CALL PJAC(Neq,Y,Yh,Nyh,Ewt,Acor,Savf,Wm,Iwm,F,JAC)
           IPUp = 0
           RC = 1.0D0
           NSLp = NST
           CRAte = 0.7D0
           IF ( IERpj.NE.0 ) GOTO 900
        ENDIF
        DO i = 1 , N
           Acor(i) = 0.0D0
        ENDDO
   700  IF ( MITer.NE.0 ) THEN
  !-----------------------------------------------------------------------
  ! In the case of the chord method, compute the corrector error,
  ! and solve the linear system with that as right-hand side and
  ! P as coefficient matrix.
  !-----------------------------------------------------------------------
           DO i = 1 , N
              Y(i) = H*Savf(i) - (Yh(i,2)+Acor(i))
           ENDDO
           CALL SLVS(Wm,Iwm,Y,Savf)
           IF ( IERsl.LT.0 ) GOTO 900
           IF ( IERsl.GT.0 ) GOTO 800
           del = DVNORM(N,Y,Ewt)
           DO i = 1 , N
              Acor(i) = Acor(i) + Y(i)
              Y(i) = Yh(i,1) + EL(1)*Acor(i)
           ENDDO
        ELSE
  !-----------------------------------------------------------------------
  ! In the case of functional iteration, update Y directly from
  ! the result of the last function evaluation.
  !-----------------------------------------------------------------------
           DO i = 1 , N
              Savf(i) = H*Savf(i) - Yh(i,2)
              Y(i) = Savf(i) - Acor(i)
           ENDDO
           del = DVNORM(N,Y,Ewt)
           DO i = 1 , N
              Y(i) = Yh(i,1) + EL(1)*Savf(i)
              Acor(i) = Savf(i)
           ENDDO
        ENDIF
  !-----------------------------------------------------------------------
  ! Test for convergence.  If M.gt.0, an estimate of the convergence
  ! rate constant is stored in CRATE, and this is used in the test.
  !-----------------------------------------------------------------------
        IF ( m.NE.0 ) CRAte = MAX(0.2D0*CRAte,del/delp)
        dcon = del*MIN(1.0D0,1.5D0*CRAte)/(TESco(2,NQ)*CONit)
        IF ( dcon.LE.1.0D0 ) THEN
  !-----------------------------------------------------------------------
  ! The corrector has converged.  JCUR is set to 0
  ! to signal that the Jacobian involved may need updating later.
  ! The local error test is made and control passes to statement 500
  ! if it fails.
  !-----------------------------------------------------------------------
           JCUr = 0
           IF ( m.EQ.0 ) dsm = del/TESco(2,NQ)
           IF ( m.GT.0 ) dsm = DVNORM(N,Acor,Ewt)/TESco(2,NQ)
           IF ( dsm.GT.1.0D0 ) THEN
  !-----------------------------------------------------------------------
  ! The error test failed.  KFLAG keeps track of multiple failures.
  ! Restore TN and the YH array to their previous values, and prepare
  ! to try the step again.  Compute the optimum step size for this or
  ! one lower order.  After 2 or more failures, H is forced to decrease
  ! by a factor of 0.2 or less.
  !-----------------------------------------------------------------------
              KFLag = KFLag - 1
              TN = told
              i1 = NQNyh + 1
              DO jb = 1 , NQ
                 i1 = i1 - Nyh
  !dir$ ivdep
                 DO i = i1 , NQNyh
                    Yh1(i) = Yh1(i) - Yh1(i+Nyh)
                 ENDDO
              ENDDO
              RMAx = 2.0D0
              IF ( ABS(H).LE.HMIn*1.00001D0 ) THEN
  !-----------------------------------------------------------------------
  ! All returns are made through this section.  H is saved in HOLD
  ! to allow the caller to change H on the next step.
  !-----------------------------------------------------------------------
                 KFLag = -1
                 GOTO 1400
              ELSEIF ( KFLag.LE.-3 ) THEN
  !-----------------------------------------------------------------------
  ! Control reaches this section if 3 or more failures have occured.
  ! If 10 failures have occurred, exit with KFLAG = -1.
  ! It is assumed that the derivatives that have accumulated in the
  ! YH array have errors of the wrong order.  Hence the first
  ! derivative is recomputed, and the order is set to 1.  Then
  ! H is reduced by a factor of 10, and the step is retried,
  ! until it succeeds or H reaches HMIN.
  !-----------------------------------------------------------------------
                 IF ( KFLag.EQ.-10 ) THEN
                    KFLag = -1
                    GOTO 1400
                 ELSE
                    rh = 0.1D0
                    rh = MAX(HMIn/ABS(H),rh)
                    H = H*rh
                    DO i = 1 , N
                       Y(i) = Yh(i,1)
                    ENDDO
                    CALL F(Neq,TN,Y,Savf)
                    NFE = NFE + 1
                    DO i = 1 , N
                       Yh(i,2) = H*Savf(i)
                    ENDDO
                    IPUp = MITer
                    IALth = 5
                    IF ( NQ.EQ.1 ) GOTO 500
                    NQ = 1
                    L = 2
                    iret = 3
                    GOTO 100
                 ENDIF
              ELSE
                 iredo = 2
                 rhup = 0.0D0
                 GOTO 1000
              ENDIF
           ELSE
  !-----------------------------------------------------------------------
  ! After a successful step, update the YH array.
  ! Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
  ! If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
  ! use in a possible order increase on the next step.
  ! If a change in H is considered, an increase or decrease in order
  ! by one is considered also.  A change in H is made only if it is by a
  ! factor of at least 1.1.  If not, IALTH is set to 3 to prevent
  ! testing for that many steps.
  !-----------------------------------------------------------------------
              KFLag = 0
              iredo = 0
              NST = NST + 1
              HU = H
              NQU = NQ
              DO j = 1 , L
                 DO i = 1 , N
                    Yh(i,j) = Yh(i,j) + EL(j)*Acor(i)
                 ENDDO
              ENDDO
              IALth = IALth - 1
              IF ( IALth.EQ.0 ) THEN
  !-----------------------------------------------------------------------
  ! Regardless of the success or failure of the step, factors
  ! RHDN, RHSM, and RHUP are computed, by which H could be multiplied
  ! at order NQ - 1, order NQ, or order NQ + 1, respectively.
  ! In the case of failure, RHUP = 0.0 to avoid an order increase.
  ! The largest of these is determined and the new order chosen
  ! accordingly.  If the order is to be increased, we compute one
  ! additional scaled derivative.
  !-----------------------------------------------------------------------
                 rhup = 0.0D0
                 IF ( L.NE.LMAx ) THEN
                    DO i = 1 , N
                       Savf(i) = Acor(i) - Yh(i,LMAx)
                    ENDDO
                    dup = DVNORM(N,Savf,Ewt)/TESco(3,NQ)
                    exup = 1.0D0/(L+1)
                    rhup = 1.0D0/(1.4D0*dup**exup+0.0000014D0)
                 ENDIF
                 GOTO 1000
              ELSE
                 IF ( IALth.LE.1 ) THEN
                    IF ( L.NE.LMAx ) THEN
                       DO i = 1 , N
                          Yh(i,LMAx) = Acor(i)
                       ENDDO
                    ENDIF
                 ENDIF
                 GOTO 1300
              ENDIF
           ENDIF
        ELSE
           m = m + 1
           IF ( m.NE.MAXcor ) THEN
              IF ( m.LT.2 .OR. del.LE.2.0D0*delp ) THEN
                 delp = del
                 CALL F(Neq,TN,Y,Savf)
                 NFE = NFE + 1
                 GOTO 700
              ENDIF
           ENDIF
        ENDIF
  !-----------------------------------------------------------------------
  ! The corrector iteration failed to converge.
  ! If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
  ! the next try.  Otherwise the YH array is retracted to its values
  ! before prediction, and H is reduced, if possible.  If H cannot be
  ! reduced or MXNCF failures have occurred, exit with KFLAG = -2.
  !-----------------------------------------------------------------------
   800  IF ( MITer.NE.0 .AND. JCUr.NE.1 ) THEN
           ICF = 1
           IPUp = MITer
           GOTO 600
        ENDIF
   900  ICF = 2
        ncf = ncf + 1
        RMAx = 2.0D0
        TN = told
        i1 = NQNyh + 1
        DO jb = 1 , NQ
           i1 = i1 - Nyh
  !dir$ ivdep
           DO i = i1 , NQNyh
              Yh1(i) = Yh1(i) - Yh1(i+Nyh)
           ENDDO
        ENDDO
        IF ( IERpj.LT.0 .OR. IERsl.LT.0 ) THEN
           KFLag = -3
           GOTO 1400
        ELSEIF ( ABS(H).LE.HMIn*1.00001D0 ) THEN
           KFLag = -2
           GOTO 1400
        ELSEIF ( ncf.EQ.MXNcf ) THEN
           KFLag = -2
           GOTO 1400
        ELSE
           rh = 0.25D0
           IPUp = MITer
           iredo = 1
           GOTO 300
        ENDIF
   1000 exsm = 1.0D0/L
        rhsm = 1.0D0/(1.2D0*dsm**exsm+0.0000012D0)
        rhdn = 0.0D0
        IF ( NQ.NE.1 ) THEN
           ddn = DVNORM(N,Yh(1,L),Ewt)/TESco(1,NQ)
           exdn = 1.0D0/NQ
           rhdn = 1.0D0/(1.3D0*ddn**exdn+0.0000013D0)
        ENDIF
        IF ( rhsm.GE.rhup ) THEN
           IF ( rhsm.GE.rhdn ) THEN
              newq = NQ
              rh = rhsm
              GOTO 1100
           ENDIF
        ELSEIF ( rhup.GT.rhdn ) THEN
           newq = L
           rh = rhup
           IF ( rh.LT.1.1D0 ) THEN
              IALth = 3
              GOTO 1300
           ELSE
              r = EL(L)/L
              DO i = 1 , N
                 Yh(i,newq+1) = Acor(i)*r
              ENDDO
              GOTO 1200
           ENDIF
        ENDIF
        newq = NQ - 1
        rh = rhdn
        IF ( KFLag.LT.0 .AND. rh.GT.1.0D0 ) rh = 1.0D0
   1100 IF ( (KFLag.EQ.0) .AND. (rh.LT.1.1D0) ) THEN
           IALth = 3
           GOTO 1300
        ELSE
           IF ( KFLag.LE.-2 ) rh = MIN(rh,0.2D0)
  !-----------------------------------------------------------------------
  ! If there is a change of order, reset NQ, l, and the coefficients.
  ! In any case H is reset according to RH and the YH array is rescaled.
  ! Then exit from 690 if the step was OK, or redo the step otherwise.
  !-----------------------------------------------------------------------
           IF ( newq.EQ.NQ ) GOTO 300
        ENDIF
   1200 NQ = newq
        L = NQ + 1
        iret = 2
        GOTO 100
   1300 r = 1.0D0/TESco(2,NQU)
        DO i = 1 , N
           Acor(i) = Acor(i)*r
        ENDDO
   1400 HOLd = H
        JSTart = 1
  !----------------------- END OF SUBROUTINE DSTODE ----------------------
        END
  !*==DEWSET.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DEWSET
        SUBROUTINE DEWSET(N,Itol,Rtol,Atol,Ycur,Ewt)
        IMPLICIT NONE
  !*--DEWSET1249
  !***BEGIN PROLOGUE  DEWSET
  !***SUBSIDIARY
  !***PURPOSE  Set error weight vector.
  !***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D)
  !***AUTHOR  Hindmarsh, Alan C., (LLNL)
  !***DESCRIPTION
  !
  !  This subroutine sets the error weight vector EWT according to
  !      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
  !  with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
  !  depending on the value of ITOL.
  !
  !***SEE ALSO  DLSODE
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   791129  DATE WRITTEN
  !   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
  !   890503  Minor cosmetic changes.  (FNF)
  !   930809  Renamed to allow single/double precision versions. (ACH)
  !***END PROLOGUE  DEWSET
  !**End
        INTEGER N , Itol
        INTEGER i
        DOUBLE PRECISION Rtol , Atol , Ycur , Ewt
        DIMENSION Rtol(*) , Atol(*) , Ycur(N) , Ewt(N)
  !
  !***FIRST EXECUTABLE STATEMENT  DEWSET
        GOTO (100,200,300,400) , Itol
   100  DO i = 1 , N
           Ewt(i) = Rtol(1)*ABS(Ycur(i)) + Atol(1)
        ENDDO
        RETURN
   200  DO i = 1 , N
           Ewt(i) = Rtol(1)*ABS(Ycur(i)) + Atol(i)
        ENDDO
        RETURN
   300  DO i = 1 , N
           Ewt(i) = Rtol(i)*ABS(Ycur(i)) + Atol(1)
        ENDDO
        RETURN
   400  DO i = 1 , N
           Ewt(i) = Rtol(i)*ABS(Ycur(i)) + Atol(i)
        ENDDO
  !----------------------- END OF SUBROUTINE DEWSET ----------------------
        END
  !*==DVNORM.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DVNORM
        DOUBLE PRECISION FUNCTION DVNORM(N,V,W)
        IMPLICIT NONE
  !*--DVNORM1299
  !***BEGIN PROLOGUE  DVNORM
  !***SUBSIDIARY
  !***PURPOSE  Weighted root-mean-square vector norm.
  !***TYPE      DOUBLE PRECISION (SVNORM-S, DVNORM-D)
  !***AUTHOR  Hindmarsh, Alan C., (LLNL)
  !***DESCRIPTION
  !
  !  This function routine computes the weighted root-mean-square norm
  !  of the vector of length N contained in the array V, with weights
  !  contained in the array W of length N:
  !    DVNORM = SQRT( (1/N) * SUM( V(i)*W(i) )**2 )
  !
  !***SEE ALSO  DLSODE
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   791129  DATE WRITTEN
  !   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
  !   890503  Minor cosmetic changes.  (FNF)
  !   930809  Renamed to allow single/double precision versions. (ACH)
  !***END PROLOGUE  DVNORM
  !**End
        INTEGER N , i
        DOUBLE PRECISION V , W , sum
        DIMENSION V(N) , W(N)
  !
  !***FIRST EXECUTABLE STATEMENT  DVNORM
        sum = 0.0D0
        DO i = 1 , N
           sum = sum + (V(i)*W(i))**2
        ENDDO
        DVNORM = SQRT(sum/N)
  !----------------------- END OF FUNCTION DVNORM ------------------------
        END
  !*==DIPREP.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DIPREP
        SUBROUTINE DIPREP(Neq,Y,Rwork,Ia,Ja,Ipflag,F,JAC)
        IMPLICIT NONE
  !*--DIPREP1337
  !*** Start of declarations inserted by SPAG
        !REAL F
        !INTEGER JAC
  !*** End of declarations inserted by SPAG
        EXTERNAL F , JAC
        INTEGER Neq , Ia , Ja , Ipflag
        DOUBLE PRECISION Y , Rwork
        DIMENSION Neq(*) , Y(*) , Rwork(*) , Ia(*) , Ja(*)
        INTEGER IOWnd , IOWns , ICF , IERpj , IERsl , JCUr , JSTart ,     &
       &        KFLag , L , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,     &
       &        METh , MITer , MAXord , MAXcor , MSBp , MXNcf , N , NQ ,  &
       &        NST , NFE , NJE , NQU
        INTEGER IPLost , IESp , ISTatc , IYS , IBA , IBIan , IBJan ,      &
       &        IBJgp , IPIan , IPJan , IPJgp , IPIgp , IPR , IPC , IPIc ,&
       &        IPIsp , IPRsp , IPA , LENyh , LENyhm , LENwk , LREq ,     &
       &        LRAt , LREst , LWMin , MOSs , MSBj , NSLj , NGP , NLU ,   &
       &        NNZ , NSP , NZL , NZU
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        DOUBLE PRECISION RLSs
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , IOWnd(6) , IOWns(6) , ICF ,    &
       &                IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , &
       &                LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,&
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        COMMON /DLSS01/ RLSs(6) , IPLost , IESp , ISTatc , IYS , IBA ,    &
       &                IBIan , IBJan , IBJgp , IPIan , IPJan , IPJgp ,   &
       &                IPIgp , IPR , IPC , IPIc , IPIsp , IPRsp , IPA ,  &
       &                LENyh , LENyhm , LENwk , LREq , LRAt , LREst ,    &
       &                LWMin , MOSs , MSBj , NSLj , NGP , NLU , NNZ ,    &
       &                NSP , NZL , NZU
        INTEGER i , imax , lewtn , lyhd , lyhn
        INTEGER, DIMENSION(LENwk) :: Rwork_i
  !-----------------------------------------------------------------------
  ! This routine serves as an interface between the driver and
  ! Subroutine DPREP.  It is called only if MITER is 1 or 2.
  ! Tasks performed here are:
  !  * call DPREP,
  !  * reset the required WM segment length LENWK,
  !  * move YH back to its final location (following WM in RWORK),
  !  * reset pointers for YH, SAVF, EWT, and ACOR, and
  !  * move EWT to its new position if ISTATE = 1.
  ! IPFLAG is an output error indication flag.  IPFLAG = 0 if there was
  ! no trouble, and IPFLAG is the value of the DPREP error flag IPPER
  ! if there was trouble in Subroutine DPREP.
  !-----------------------------------------------------------------------
        Ipflag = 0
  ! Call DPREP to do matrix preprocessing operations. --------------------

        Rwork_i = Rwork(LWM)
        CALL DPREP(Neq,Y,Rwork(LYH),Rwork(LSAvf),Rwork(LEWt),Rwork(LACor),&
       &           Ia,Ja,Rwork(LWM),Rwork_i,Ipflag,F,JAC)
        LENwk = MAX(LREq,LWMin)
        IF ( Ipflag.LT.0 ) RETURN
  ! If DPREP was successful, move YH to end of required space for WM. ----
        lyhn = LWM + LENwk
        IF ( lyhn.GT.LYH ) RETURN
        lyhd = LYH - lyhn
        IF ( lyhd.NE.0 ) THEN
           imax = lyhn - 1 + LENyhm
           DO i = lyhn , imax
              Rwork(i) = Rwork(i+lyhd)
           ENDDO
           LYH = lyhn
        ENDIF
  ! Reset pointers for SAVF, EWT, and ACOR. ------------------------------
        LSAvf = LYH + LENyh
        lewtn = LSAvf + N
        LACor = lewtn + N
        IF ( ISTatc.NE.3 ) THEN
  ! If ISTATE = 1, move EWT (left) to its new position. ------------------
           IF ( lewtn.GT.LEWt ) RETURN
           DO i = 1 , N
              Rwork(i+lewtn-1) = Rwork(i+LEWt-1)
           ENDDO
        ENDIF
        LEWt = lewtn
  !----------------------- End of Subroutine DIPREP ----------------------
        END
  !*==DPREP.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DPREP
        SUBROUTINE DPREP(Neq,Y,Yh,Savf,Ewt,Ftem,Ia,Ja,Wk,Iwk,Ipper,F,JAC)
        IMPLICIT NONE
  !*--DPREP1419
        EXTERNAL F , JAC
        INTEGER Neq , Ia , Ja , Iwk , Ipper
        DOUBLE PRECISION Y , Yh , Savf , Ewt , Ftem , Wk
        DIMENSION Neq(*) , Y(*) , Yh(*) , Savf(*) , Ewt(*) , Ftem(*) ,    &
       &          Ia(*) , Ja(*) , Wk(*) , Iwk(*)
        INTEGER IOWnd , IOWns , ICF , IERpj , IERsl , JCUr , JSTart ,     &
       &        KFLag , L , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,     &
       &        METh , MITer , MAXord , MAXcor , MSBp , MXNcf , N , NQ ,  &
       &        NST , NFE , NJE , NQU
        INTEGER IPLost , IESp , ISTatc , IYS , IBA , IBIan , IBJan ,      &
       &        IBJgp , IPIan , IPJan , IPJgp , IPIgp , IPR , IPC , IPIc ,&
       &        IPIsp , IPRsp , IPA , LENyh , LENyhm , LENwk , LREq ,     &
       &        LRAt , LREst , LWMin , MOSs , MSBj , NSLj , NGP , NLU ,   &
       &        NNZ , NSP , NZL , NZU
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        DOUBLE PRECISION CON0 , CONmin , CCMxj , PSMall , RBIg , SETh
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , IOWnd(6) , IOWns(6) , ICF ,    &
       &                IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , &
       &                LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,&
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        COMMON /DLSS01/ CON0 , CONmin , CCMxj , PSMall , RBIg , SETh ,    &
       &                IPLost , IESp , ISTatc , IYS , IBA , IBIan ,      &
       &                IBJan , IBJgp , IPIan , IPJan , IPJgp , IPIgp ,   &
       &                IPR , IPC , IPIc , IPIsp , IPRsp , IPA , LENyh ,  &
       &                LENyhm , LENwk , LREq , LRAt , LREst , LWMin ,    &
       &                MOSs , MSBj , NSLj , NGP , NLU , NNZ , NSP , NZL ,&
       &                NZU
        INTEGER i , ibr , ier , ipil , ipiu , iptt1 , iptt2 , j , jfound ,&
       &        k , knew , kmax , kmin , ldif , lenigp , liwk , maxg ,    &
       &        np1 , nzsut
        DOUBLE PRECISION dq , dyj , erwt , fac , yj
  !-----------------------------------------------------------------------
  ! This routine performs preprocessing related to the sparse linear
  ! systems that must be solved if MITER = 1 or 2.
  ! The operations that are performed here are:
  !  * compute sparseness structure of Jacobian according to MOSS,
  !  * compute grouping of column indices (MITER = 2),
  !  * compute a new ordering of rows and columns of the matrix,
  !  * reorder JA corresponding to the new ordering,
  !  * perform a symbolic LU factorization of the matrix, and
  !  * set pointers for segments of the IWK/WK array.
  ! In addition to variables described previously, DPREP uses the
  ! following for communication:
  ! YH     = the history array.  Only the first column, containing the
  !          current Y vector, is used.  Used only if MOSS .ne. 0.
  ! SAVF   = a work array of length NEQ, used only if MOSS .ne. 0.
  ! EWT    = array of length NEQ containing (inverted) error weights.
  !          Used only if MOSS = 2 or if ISTATE = MOSS = 1.
  ! FTEM   = a work array of length NEQ, identical to ACOR in the driver,
  !          used only if MOSS = 2.
  ! WK     = a real work array of length LENWK, identical to WM in
  !          the driver.
  ! IWK    = integer work array, assumed to occupy the same space as WK.
  ! LENWK  = the length of the work arrays WK and IWK.
  ! ISTATC = a copy of the driver input argument ISTATE (= 1 on the
  !          first call, = 3 on a continuation call).
  ! IYS    = flag value from ODRV or CDRV.
  ! IPPER  = output error flag with the following values and meanings:
  !          0  no error.
  !         -1  insufficient storage for internal structure pointers.
  !         -2  insufficient storage for JGROUP.
  !         -3  insufficient storage for ODRV.
  !         -4  other error flag from ODRV (should never occur).
  !         -5  insufficient storage for CDRV.
  !         -6  other error flag from CDRV.
  !-----------------------------------------------------------------------
        IBIan = LRAt*2
        IPIan = IBIan + 1
        np1 = N + 1
        IPJan = IPIan + np1
        IBJan = IPJan - 1
        liwk = LENwk*LRAt
        IF ( IPJan+N-1.GT.liwk ) GOTO 400
        IF ( MOSs.NE.0 ) THEN
  !
           IF ( ISTatc.NE.3 ) THEN
  ! ISTATE = 1 and MOSS .ne. 0.  Perturb Y for structure determination. --
              DO i = 1 , N
                 erwt = 1.0D0/Ewt(i)
                 fac = 1.0D0 + 1.0D0/(i+1.0D0)
                 Y(i) = Y(i) + fac*SIGN(erwt,Y(i))
              ENDDO
              GOTO (100,200) , MOSs
           ENDIF
  !
  ! ISTATE = 3 and MOSS .ne. 0.  Load Y from YH(*,1). --------------------
           DO i = 1 , N
              Y(i) = Yh(i)
           ENDDO
           GOTO (100,200) , MOSs
        ENDIF
  !
  ! MOSS = 0.  Process user's IA,JA.  Add diagonal entries if necessary. -
        knew = IPJan
        kmin = Ia(1)
        Iwk(IPIan) = 1
        DO j = 1 , N
           jfound = 0
           kmax = Ia(j+1) - 1
           IF ( kmin.LE.kmax ) THEN
              DO k = kmin , kmax
                 i = Ja(k)
                 IF ( i.EQ.j ) jfound = 1
                 IF ( knew.GT.liwk ) GOTO 400
                 Iwk(knew) = i
                 knew = knew + 1
              ENDDO
              IF ( jfound.EQ.1 ) GOTO 50
           ENDIF
           IF ( knew.GT.liwk ) GOTO 400
           Iwk(knew) = j
           knew = knew + 1
   50      Iwk(IPIan+j) = knew + 1 - IPJan
           kmin = kmax + 1
        ENDDO
        GOTO 300
  !
  ! MOSS = 1.  Compute structure from user-supplied Jacobian routine JAC.
  ! A dummy call to F allows user to create temporaries for use in JAC. --
   100  CALL F(Neq,TN,Y,Savf)
        k = IPJan
        Iwk(IPIan) = 1
        DO j = 1 , N
           IF ( k.GT.liwk ) GOTO 400
           Iwk(k) = j
           k = k + 1
           DO i = 1 , N
              Savf(i) = 0.0D0
           ENDDO
           CALL JAC(Neq,TN,Y,j,Iwk(IPIan),Iwk(IPJan),Savf)
           DO i = 1 , N
              IF ( ABS(Savf(i)).GT.SETh ) THEN
                 IF ( i.NE.j ) THEN
                    IF ( k.GT.liwk ) GOTO 400
                    Iwk(k) = i
                    k = k + 1
                 ENDIF
              ENDIF
           ENDDO
           Iwk(IPIan+j) = k + 1 - IPJan
        ENDDO
        GOTO 300
  !
  ! MOSS = 2.  Compute structure from results of N + 1 calls to F. -------
   200  k = IPJan
        Iwk(IPIan) = 1
        CALL F(Neq,TN,Y,Savf)
        DO j = 1 , N
           IF ( k.GT.liwk ) GOTO 400
           Iwk(k) = j
           k = k + 1
           yj = Y(j)
           erwt = 1.0D0/Ewt(j)
           dyj = SIGN(erwt,yj)
           Y(j) = yj + dyj
           CALL F(Neq,TN,Y,Ftem)
           Y(j) = yj
           DO i = 1 , N
              dq = (Ftem(i)-Savf(i))/dyj
              IF ( ABS(dq).GT.SETh ) THEN
                 IF ( i.NE.j ) THEN
                    IF ( k.GT.liwk ) GOTO 400
                    Iwk(k) = i
                    k = k + 1
                 ENDIF
              ENDIF
           ENDDO
           Iwk(IPIan+j) = k + 1 - IPJan
        ENDDO
  !
   300  IF ( MOSs.NE.0 .AND. ISTatc.EQ.1 ) THEN
  ! If ISTATE = 1 and MOSS .ne. 0, restore Y from YH. --------------------
           DO i = 1 , N
              Y(i) = Yh(i)
           ENDDO
        ENDIF
        NNZ = Iwk(IPIan+N) - 1
        lenigp = 0
        IPIgp = IPJan + NNZ
        IF ( MITer.EQ.2 ) THEN
  !
  ! Compute grouping of column indices (MITER = 2). ----------------------
           maxg = np1
           IPJgp = IPJan + NNZ
           IBJgp = IPJgp - 1
           IPIgp = IPJgp + N
           iptt1 = IPIgp + np1
           iptt2 = iptt1 + N
           LREq = iptt2 + N - 1
           IF ( LREq.GT.liwk ) GOTO 500
           CALL JGROUP(N,Iwk(IPIan),Iwk(IPJan),maxg,NGP,Iwk(IPIgp),       &
       &               Iwk(IPJgp),Iwk(iptt1),Iwk(iptt2),ier)
           IF ( ier.NE.0 ) GOTO 500
           lenigp = NGP + 1
        ENDIF
  !
  ! Compute new ordering of rows/columns of Jacobian. --------------------
        IPR = IPIgp + lenigp
        IPC = IPR
        IPIc = IPC + N
        IPIsp = IPIc + N
        IPRsp = (IPIsp-2)/LRAt + 2
        IESp = LENwk + 1 - IPRsp
        IF ( IESp.LT.0 ) GOTO 600
        ibr = IPR - 1
        DO i = 1 , N
           Iwk(ibr+i) = i
        ENDDO
        NSP = liwk + 1 - IPIsp
        CALL ODRV(N,Iwk(IPIan),Iwk(IPJan),Wk,Iwk(IPR),Iwk(IPIc),NSP,      &
       &          Iwk(IPIsp),1,IYS)
        IF ( IYS.EQ.11*N+1 ) THEN
  !
           Ipper = -4
           RETURN
        ELSE
           IF ( IYS.NE.0 ) GOTO 600
  !
  ! Reorder JAN and do symbolic LU factorization of matrix. --------------
           IPA = LENwk + 1 - NNZ
           NSP = IPA - IPRsp
           LREq = MAX(12*N/LRAt,6*N/LRAt+2*N+NNZ) + 3
           LREq = LREq + IPRsp - 1 + NNZ
           IF ( LREq.GT.LENwk ) GOTO 700
           IBA = IPA - 1
           DO i = 1 , NNZ
              Wk(IBA+i) = 0.0D0
           ENDDO
           IPIsp = LRAt*(IPRsp-1) + 1
           CALL CDRV(N,Iwk(IPR),Iwk(IPC),Iwk(IPIc),Iwk(IPIan),Iwk(IPJan), &
       &             Wk(IPA),Wk(IPA),Wk(IPA),NSP,Iwk(IPIsp),Wk(IPRsp),    &
       &             IESp,5,IYS)
           LREq = LENwk - IESp
           IF ( IYS.EQ.10*N+1 ) GOTO 700
           IF ( IYS.NE.0 ) THEN
  !
              Ipper = -6
              LREq = LENwk
              GOTO 99999
           ELSE
              ipil = IPIsp
              ipiu = ipil + 2*N + 1
              NZU = Iwk(ipil+N) - Iwk(ipil)
              NZL = Iwk(ipiu+N) - Iwk(ipiu)
              IF ( LRAt.LE.1 ) THEN
                 CALL ADJLR(N,Iwk(IPIsp),ldif)
                 LREq = LREq + ldif
              ENDIF
              IF ( LRAt.EQ.2 .AND. NNZ.EQ.N ) LREq = LREq + 1
              NSP = NSP + LREq - LENwk
              IPA = LREq + 1 - NNZ
              IBA = IPA - 1
              Ipper = 0
              RETURN
           ENDIF
        ENDIF
  !
   400  Ipper = -1
        LREq = 2 + (2*N+1)/LRAt
        LREq = MAX(LENwk+1,LREq)
        RETURN
  !
   500  Ipper = -2
        LREq = (LREq-1)/LRAt + 1
        RETURN
  !
   600  Ipper = -3
        CALL CNTNZU(N,Iwk(IPIan),Iwk(IPJan),nzsut)
        LREq = LENwk - IESp + (3*N+4*nzsut-1)/LRAt + 1
        RETURN
  !
   700  Ipper = -5
        RETURN
  !----------------------- End of Subroutine DPREP -----------------------
  99999 END
  !*==JGROUP.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK JGROUP
        SUBROUTINE JGROUP(N,Ia,Ja,Maxg,Ngrp,Igp,Jgp,Incl,Jdone,Ier)
        IMPLICIT NONE
  !*--JGROUP1702
        INTEGER N , Ia , Ja , Maxg , Ngrp , Igp , Jgp , Incl , Jdone , Ier
        DIMENSION Ia(*) , Ja(*) , Igp(*) , Jgp(*) , Incl(*) , Jdone(*)
  !-----------------------------------------------------------------------
  ! This subroutine constructs groupings of the column indices of
  ! the Jacobian matrix, used in the numerical evaluation of the
  ! Jacobian by finite differences.
  !
  ! Input:
  ! N      = the order of the matrix.
  ! IA,JA  = sparse structure descriptors of the matrix by rows.
  ! MAXG   = length of available storage in the IGP array.
  !
  ! Output:
  ! NGRP   = number of groups.
  ! JGP    = array of length N containing the column indices by groups.
  ! IGP    = pointer array of length NGRP + 1 to the locations in JGP
  !          of the beginning of each group.
  ! IER    = error indicator.  IER = 0 if no error occurred, or 1 if
  !          MAXG was insufficient.
  !
  ! INCL and JDONE are working arrays of length N.
  !-----------------------------------------------------------------------
        INTEGER i , j , k , kmin , kmax , ncol , ng
  !
        Ier = 0
        DO j = 1 , N
           Jdone(j) = 0
        ENDDO
        ncol = 1
        DO ng = 1 , Maxg
           Igp(ng) = ncol
           DO i = 1 , N
              Incl(i) = 0
           ENDDO
           DO j = 1 , N
  ! Reject column J if it is already in a group.--------------------------
              IF ( Jdone(j).NE.1 ) THEN
                 kmin = Ia(j)
                 kmax = Ia(j+1) - 1
                 DO k = kmin , kmax
  ! Reject column J if it overlaps any column already in this group.------
                    i = Ja(k)
                    IF ( Incl(i).EQ.1 ) GOTO 50
                 ENDDO
  ! Accept column J into group NG.----------------------------------------
                 Jgp(ncol) = j
                 ncol = ncol + 1
                 Jdone(j) = 1
                 DO k = kmin , kmax
                    i = Ja(k)
                    Incl(i) = 1
                 ENDDO
              ENDIF
   50      ENDDO
  ! Stop if this group is empty (grouping is complete).-------------------
           IF ( ncol.EQ.Igp(ng) ) GOTO 100
        ENDDO
  ! Error return if not all columns were chosen (MAXG too small).---------
        IF ( ncol.LE.N ) THEN
           Ier = 1
           GOTO 99999
        ELSE
           ng = Maxg
        ENDIF
   100  Ngrp = ng - 1
        RETURN
  !----------------------- End of Subroutine JGROUP ----------------------
  99999 END
  !*==ADJLR.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK ADJLR
        SUBROUTINE ADJLR(N,Isp,Ldif)
        IMPLICIT NONE
  !*--ADJLR1775
        INTEGER N , Isp , Ldif
        DIMENSION Isp(*)
  !-----------------------------------------------------------------------
  ! This routine computes an adjustment, LDIF, to the required
  ! integer storage space in IWK (sparse matrix work space).
  ! It is called only if the word length ratio is LRAT = 1.
  ! This is to account for the possibility that the symbolic LU phase
  ! may require more storage than the numerical LU and solution phases.
  !-----------------------------------------------------------------------
        INTEGER ip , jlmax , jumax , lnfc , lsfc , nzlu
  !
        ip = 2*N + 1
  ! Get JLMAX = IJL(N) and JUMAX = IJU(N) (sizes of JL and JU). ----------
        jlmax = Isp(ip)
        jumax = Isp(ip+ip)
  ! NZLU = (size of L) + (size of U) = (IL(N+1)-IL(1)) + (IU(N+1)-IU(1)).
        nzlu = Isp(N+1) - Isp(1) + Isp(ip+N+1) - Isp(ip+1)
        lsfc = 12*N + 3 + 2*MAX(jlmax,jumax)
        lnfc = 9*N + 2 + jlmax + jumax + nzlu
        Ldif = MAX(0,lsfc-lnfc)
  !----------------------- End of Subroutine ADJLR -----------------------
        END
  !*==CNTNZU.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK CNTNZU
        SUBROUTINE CNTNZU(N,Ia,Ja,Nzsut)
        IMPLICIT NONE
  !*--CNTNZU1802
        INTEGER N , Ia , Ja , Nzsut
        DIMENSION Ia(*) , Ja(*)
  !-----------------------------------------------------------------------
  ! This routine counts the number of nonzero elements in the strict
  ! upper triangle of the matrix M + M(transpose), where the sparsity
  ! structure of M is given by pointer arrays IA and JA.
  ! This is needed to compute the storage requirements for the
  ! sparse matrix reordering operation in ODRV.
  !-----------------------------------------------------------------------
        INTEGER ii , jj , j , jmin , jmax , k , kmin , kmax , num
  !
        num = 0
        DO ii = 1 , N
           jmin = Ia(ii)
           jmax = Ia(ii+1) - 1
           IF ( jmin.LE.jmax ) THEN
              DO j = jmin , jmax
                 IF ( Ja(j).LT.ii ) THEN
                    jj = Ja(j)
                    kmin = Ia(jj)
                    kmax = Ia(jj+1) - 1
                    IF ( kmin.LE.kmax ) THEN
                       DO k = kmin , kmax
                          IF ( Ja(k).EQ.ii ) GOTO 20
                       ENDDO
                    ENDIF
                 ELSEIF ( Ja(j).EQ.ii ) THEN
                    GOTO 20
                 ENDIF
                 num = num + 1
   20         ENDDO
           ENDIF
        ENDDO
        Nzsut = num
  !----------------------- End of Subroutine CNTNZU ----------------------
        END
  !*==DPRJS.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DPRJS
        SUBROUTINE DPRJS(Neq,Y,Yh,Nyh,Ewt,Ftem,Savf,Wk,Iwk,F,JAC)
        IMPLICIT NONE
  !*--DPRJS1843
        EXTERNAL F , JAC
        INTEGER Neq , Nyh , Iwk
        DOUBLE PRECISION Y , Yh , Ewt , Ftem , Savf , Wk
        DIMENSION Neq(*) , Y(*) , Yh(Nyh,*) , Ewt(*) , Ftem(*) , Savf(*) ,&
       &          Wk(*) , Iwk(*)
        INTEGER IOWnd , IOWns , ICF , IERpj , IERsl , JCUr , JSTart ,     &
       &        KFLag , L , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,     &
       &        METh , MITer , MAXord , MAXcor , MSBp , MXNcf , N , NQ ,  &
       &        NST , NFE , NJE , NQU
        INTEGER IPLost , IESp , ISTatc , IYS , IBA , IBIan , IBJan ,      &
       &        IBJgp , IPIan , IPJan , IPJgp , IPIgp , IPR , IPC , IPIc ,&
       &        IPIsp , IPRsp , IPA , LENyh , LENyhm , LENwk , LREq ,     &
       &        LRAt , LREst , LWMin , MOSs , MSBj , NSLj , NGP , NLU ,   &
       &        NNZ , NSP , NZL , NZU
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        DOUBLE PRECISION CON0 , CONmin , CCMxj , PSMall , RBIg , SETh
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , IOWnd(6) , IOWns(6) , ICF ,    &
       &                IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , &
       &                LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,&
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        COMMON /DLSS01/ CON0 , CONmin , CCMxj , PSMall , RBIg , SETh ,    &
       &                IPLost , IESp , ISTatc , IYS , IBA , IBIan ,      &
       &                IBJan , IBJgp , IPIan , IPJan , IPJgp , IPIgp ,   &
       &                IPR , IPC , IPIc , IPIsp , IPRsp , IPA , LENyh ,  &
       &                LENyhm , LENwk , LREq , LRAt , LREst , LWMin ,    &
       &                MOSs , MSBj , NSLj , NGP , NLU , NNZ , NSP , NZL ,&
       &                NZU
        INTEGER i , imul , j , jj , jok , jmax , jmin , k , kmax , kmin , &
       &        ng
        DOUBLE PRECISION con , di , fac , hl0 , pij , r , r0 , rcon ,     &
       &                 rcont , srur
  !-----------------------------------------------------------------------
  ! DPRJS is called to compute and process the matrix
  ! P = I - H*EL(1)*J , where J is an approximation to the Jacobian.
  ! J is computed by columns, either by the user-supplied routine JAC
  ! if MITER = 1, or by finite differencing if MITER = 2.
  ! if MITER = 3, a diagonal approximation to J is used.
  ! if MITER = 1 or 2, and if the existing value of the Jacobian
  ! (as contained in P) is considered acceptable, then a new value of
  ! P is reconstructed from the old value.  In any case, when MITER
  ! is 1 or 2, the P matrix is subjected to LU decomposition in CDRV.
  ! P and its LU decomposition are stored (separately) in WK.
  !
  ! In addition to variables described previously, communication
  ! with DPRJS uses the following:
  ! Y     = array containing predicted values on entry.
  ! FTEM  = work array of length N (ACOR in DSTODE).
  ! SAVF  = array containing f evaluated at predicted y.
  ! WK    = real work space for matrices.  On output it contains the
  !         inverse diagonal matrix if MITER = 3, and P and its sparse
  !         LU decomposition if MITER is 1 or 2.
  !         Storage of matrix elements starts at WK(3).
  !         WK also contains the following matrix-related data:
  !         WK(1) = SQRT(UROUND), used in numerical Jacobian increments.
  !         WK(2) = H*EL0, saved for later use if MITER = 3.
  ! IWK   = integer work space for matrix-related data, assumed to
  !         be equivalenced to WK.  In addition, WK(IPRSP) and IWK(IPISP)
  !         are assumed to have identical locations.
  ! EL0   = EL(1) (input).
  ! IERPJ = output error flag (in Common).
  !       = 0 if no error.
  !       = 1  if zero pivot found in CDRV.
  !       = 2  if a singular matrix arose with MITER = 3.
  !       = -1 if insufficient storage for CDRV (should not occur here).
  !       = -2 if other error found in CDRV (should not occur here).
  ! JCUR  = output flag showing status of (approximate) Jacobian matrix:
  !          = 1 to indicate that the Jacobian is now current, or
  !          = 0 to indicate that a saved value was used.
  ! This routine also uses other variables in Common.
  !-----------------------------------------------------------------------
        hl0 = H*EL0
        con = -hl0
        IF ( MITer.EQ.3 ) THEN
  !
  ! If MITER = 3, construct a diagonal approximation to J and P. ---------
           JCUr = 1
           NJE = NJE + 1
           Wk(2) = hl0
           IERpj = 0
           r = EL0*0.1D0
           DO i = 1 , N
              Y(i) = Y(i) + r*(H*Savf(i)-Yh(i,2))
           ENDDO
           CALL F(Neq,TN,Y,Wk(3))
           NFE = NFE + 1
           DO i = 1 , N
              r0 = H*Savf(i) - Yh(i,2)
              di = 0.1D0*r0 - H*(Wk(i+2)-Savf(i))
              Wk(i+2) = 1.0D0
              IF ( ABS(r0).GE.UROund/Ewt(i) ) THEN
                 IF ( ABS(di).EQ.0.0D0 ) GOTO 400
                 Wk(i+2) = 0.1D0*r0/di
              ENDIF
           ENDDO
           RETURN
        ELSE
  ! See whether J should be reevaluated (JOK = 0) or not (JOK = 1). ------
           jok = 1
           IF ( NST.EQ.0 .OR. NST.GE.NSLj+MSBj ) jok = 0
           IF ( ICF.EQ.1 .AND. ABS(RC-1.0D0).LT.CCMxj ) jok = 0
           IF ( ICF.EQ.2 ) jok = 0
           IF ( jok.EQ.1 ) THEN
  !
  ! If JOK = 1, reconstruct new P from old P. ----------------------------
              JCUr = 0
              rcon = con/CON0
              rcont = ABS(con)/CONmin
              IF ( rcont.LE.RBIg .OR. IPLost.NE.1 ) THEN
                 kmin = Iwk(IPIan)
                 DO j = 1 , N
                    kmax = Iwk(IPIan+j) - 1
                    DO k = kmin , kmax
                       i = Iwk(IBJan+k)
                       pij = Wk(IBA+k)
                       IF ( i.EQ.j ) THEN
                          pij = pij - 1.0D0
                          IF ( ABS(pij).LT.PSMall ) THEN
                             IPLost = 1
                             CONmin = MIN(ABS(CON0),CONmin)
                          ENDIF
                       ENDIF
                       pij = pij*rcon
                       IF ( i.EQ.j ) pij = pij + 1.0D0
                       Wk(IBA+k) = pij
                    ENDDO
                    kmin = kmax + 1
                 ENDDO
                 GOTO 300
              ENDIF
           ENDIF
  !
  ! MITER = 1 or 2, and the Jacobian is to be reevaluated. ---------------
           JCUr = 1
           NJE = NJE + 1
           NSLj = NST
           IPLost = 0
           CONmin = ABS(con)
           GOTO (100,200) , MITer
        ENDIF
  !
  ! If MITER = 1, call JAC, multiply by scalar, and add identity. --------
   100  kmin = Iwk(IPIan)
        DO j = 1 , N
           kmax = Iwk(IPIan+j) - 1
           DO i = 1 , N
              Ftem(i) = 0.0D0
           ENDDO
           CALL JAC(Neq,TN,Y,j,Iwk(IPIan),Iwk(IPJan),Ftem)
           DO k = kmin , kmax
              i = Iwk(IBJan+k)
              Wk(IBA+k) = Ftem(i)*con
              IF ( i.EQ.j ) Wk(IBA+k) = Wk(IBA+k) + 1.0D0
           ENDDO
           kmin = kmax + 1
        ENDDO
        GOTO 300
  !
  ! If MITER = 2, make NGP calls to F to approximate J and P. ------------
   200  fac = DVNORM(N,Savf,Ewt)
        r0 = 1000.0D0*ABS(H)*UROund*N*fac
        IF ( r0.EQ.0.0D0 ) r0 = 1.0D0
        srur = Wk(1)
        jmin = Iwk(IPIgp)
        DO ng = 1 , NGP
           jmax = Iwk(IPIgp+ng) - 1
           DO j = jmin , jmax
              jj = Iwk(IBJgp+j)
              r = MAX(srur*ABS(Y(jj)),r0/Ewt(jj))
              Y(jj) = Y(jj) + r
           ENDDO
           CALL F(Neq,TN,Y,Ftem)
           DO j = jmin , jmax
              jj = Iwk(IBJgp+j)
              Y(jj) = Yh(jj,1)
              r = MAX(srur*ABS(Y(jj)),r0/Ewt(jj))
              fac = -hl0/r
              kmin = Iwk(IBIan+jj)
              kmax = Iwk(IBIan+jj+1) - 1
              DO k = kmin , kmax
                 i = Iwk(IBJan+k)
                 Wk(IBA+k) = (Ftem(i)-Savf(i))*fac
                 IF ( i.EQ.jj ) Wk(IBA+k) = Wk(IBA+k) + 1.0D0
              ENDDO
           ENDDO
           jmin = jmax + 1
        ENDDO
        NFE = NFE + NGP
  !
  ! Do numerical factorization of P matrix. ------------------------------
   300  NLU = NLU + 1
        CON0 = con
        IERpj = 0
        DO i = 1 , N
           Ftem(i) = 0.0D0
        ENDDO
        CALL CDRV(N,Iwk(IPR),Iwk(IPC),Iwk(IPIc),Iwk(IPIan),Iwk(IPJan),    &
       &          Wk(IPA),Ftem,Ftem,NSP,Iwk(IPIsp),Wk(IPRsp),IESp,2,IYS)
        IF ( IYS.EQ.0 ) RETURN
        imul = (IYS-1)/N
        IERpj = -2
        IF ( imul.EQ.8 ) IERpj = 1
        IF ( imul.EQ.10 ) IERpj = -1
        RETURN
   400  IERpj = 2
  !----------------------- End of Subroutine DPRJS -----------------------
        END
  !*==DSOLSS.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DSOLSS
        SUBROUTINE DSOLSS(Wk,Iwk,X,Tem)
        IMPLICIT NONE
  !*--DSOLSS2057
        INTEGER Iwk
        DOUBLE PRECISION Wk , X , Tem
        DIMENSION Wk(*) , Iwk(*) , X(*) , Tem(*)
        INTEGER IOWnd , IOWns , ICF , IERpj , IERsl , JCUr , JSTart ,     &
       &        KFLag , L , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,     &
       &        METh , MITer , MAXord , MAXcor , MSBp , MXNcf , N , NQ ,  &
       &        NST , NFE , NJE , NQU
        INTEGER IPLost , IESp , ISTatc , IYS , IBA , IBIan , IBJan ,      &
       &        IBJgp , IPIan , IPJan , IPJgp , IPIgp , IPR , IPC , IPIc ,&
       &        IPIsp , IPRsp , IPA , LENyh , LENyhm , LENwk , LREq ,     &
       &        LRAt , LREst , LWMin , MOSs , MSBj , NSLj , NGP , NLU ,   &
       &        NNZ , NSP , NZL , NZU
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        DOUBLE PRECISION RLSs
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , IOWnd(6) , IOWns(6) , ICF ,    &
       &                IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , &
       &                LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,&
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        COMMON /DLSS01/ RLSs(6) , IPLost , IESp , ISTatc , IYS , IBA ,    &
       &                IBIan , IBJan , IBJgp , IPIan , IPJan , IPJgp ,   &
       &                IPIgp , IPR , IPC , IPIc , IPIsp , IPRsp , IPA ,  &
       &                LENyh , LENyhm , LENwk , LREq , LRAt , LREst ,    &
       &                LWMin , MOSs , MSBj , NSLj , NGP , NLU , NNZ ,    &
       &                NSP , NZL , NZU
        INTEGER i
        DOUBLE PRECISION di , hl0 , phl0 , r
  !-----------------------------------------------------------------------
  ! This routine manages the solution of the linear system arising from
  ! a chord iteration.  It is called if MITER .ne. 0.
  ! If MITER is 1 or 2, it calls CDRV to accomplish this.
  ! If MITER = 3 it updates the coefficient H*EL0 in the diagonal
  ! matrix, and then computes the solution.
  ! communication with DSOLSS uses the following variables:
  ! WK    = real work space containing the inverse diagonal matrix if
  !         MITER = 3 and the LU decomposition of the matrix otherwise.
  !         Storage of matrix elements starts at WK(3).
  !         WK also contains the following matrix-related data:
  !         WK(1) = SQRT(UROUND) (not used here),
  !         WK(2) = HL0, the previous value of H*EL0, used if MITER = 3.
  ! IWK   = integer work space for matrix-related data, assumed to
  !         be equivalenced to WK.  In addition, WK(IPRSP) and IWK(IPISP)
  !         are assumed to have identical locations.
  ! X     = the right-hand side vector on input, and the solution vector
  !         on output, of length N.
  ! TEM   = vector of work space of length N, not used in this version.
  ! IERSL = output flag (in Common).
  !         IERSL = 0  if no trouble occurred.
  !         IERSL = -1 if CDRV returned an error flag (MITER = 1 or 2).
  !                    This should never occur and is considered fatal.
  !         IERSL = 1  if a singular matrix arose with MITER = 3.
  ! This routine also uses other variables in Common.
  !-----------------------------------------------------------------------
        IERsl = 0
        GOTO (100,100,200) , MITer
   100  CALL CDRV(N,Iwk(IPR),Iwk(IPC),Iwk(IPIc),Iwk(IPIan),Iwk(IPJan),    &
       &          Wk(IPA),X,X,NSP,Iwk(IPIsp),Wk(IPRsp),IESp,4,IERsl)
        IF ( IERsl.NE.0 ) IERsl = -1
        RETURN
  !
   200  phl0 = Wk(2)
        hl0 = H*EL0
        Wk(2) = hl0
        IF ( hl0.NE.phl0 ) THEN
           r = hl0/phl0
           DO i = 1 , N
              di = 1.0D0 - r*(1.0D0-1.0D0/Wk(i+2))
              IF ( ABS(di).EQ.0.0D0 ) GOTO 300
              Wk(i+2) = 1.0D0/di
           ENDDO
        ENDIF
        DO i = 1 , N
           X(i) = Wk(i+2)*X(i)
        ENDDO
        RETURN
   300  IERsl = 1
  !
  !----------------------- End of Subroutine DSOLSS ----------------------
        END
  !*==DSRCMS.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DSRCMS
        SUBROUTINE DSRCMS(Rsav,Isav,Job)
        IMPLICIT NONE
  !*--DSRCMS2143
  !-----------------------------------------------------------------------
  ! This routine saves or restores (depending on JOB) the contents of
  ! the Common blocks DLS001, DLSS01, which are used
  ! internally by one or more ODEPACK solvers.
  !
  ! RSAV = real array of length 224 or more.
  ! ISAV = integer array of length 71 or more.
  ! JOB  = flag indicating to save or restore the Common blocks:
  !        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
  !        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
  !        A call with JOB = 2 presumes a prior call with JOB = 1.
  !-----------------------------------------------------------------------
        INTEGER Isav , Job
        INTEGER ILS , ILSs
        INTEGER i , lenils , leniss , lenrls , lenrss
        DOUBLE PRECISION Rsav , RLS , RLSs
        DIMENSION Rsav(*) , Isav(*)
        SAVE lenrls , lenils , lenrss , leniss
        COMMON /DLS001/ RLS(218) , ILS(37)
        COMMON /DLSS01/ RLSs(6) , ILSs(34)
        DATA lenrls/218/ , lenils/37/ , lenrss/6/ , leniss/34/
  !
        IF ( Job.EQ.2 ) THEN
  !
           DO i = 1 , lenrls
              RLS(i) = Rsav(i)
           ENDDO
           DO i = 1 , lenrss
              RLSs(i) = Rsav(lenrls+i)
           ENDDO
  !
           DO i = 1 , lenils
              ILS(i) = Isav(i)
           ENDDO
           DO i = 1 , leniss
              ILSs(i) = Isav(lenils+i)
           ENDDO
           GOTO 99999
        ENDIF
        DO i = 1 , lenrls
           Rsav(i) = RLS(i)
        ENDDO
        DO i = 1 , lenrss
           Rsav(lenrls+i) = RLSs(i)
        ENDDO
  !
        DO i = 1 , lenils
           Isav(i) = ILS(i)
        ENDDO
        DO i = 1 , leniss
           Isav(lenils+i) = ILSs(i)
        ENDDO
  !
        RETURN
  !
  !----------------------- End of Subroutine DSRCMS ----------------------
  99999 END
  !*==ODRV.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK ODRV
        SUBROUTINE ODRV(N,Ia,Ja,A,P,Ip,Nsp,Isp,Path,Flag)
        IMPLICIT NONE
  !*--ODRV2205
  !*** Start of declarations inserted by SPAG
        INTEGER max , N , next , Nsp
  !*** End of declarations inserted by SPAG
  !                                                                 5/2/83
  !***********************************************************************
  !  odrv -- driver for sparse matrix reordering routines
  !***********************************************************************
  !
  !  description
  !
  !    odrv finds a minimum degree ordering of the rows and columns
  !    of a matrix m stored in (ia,ja,a) format (see below).  for the
  !    reordered matrix, the work and storage required to perform
  !    gaussian elimination is (usually) significantly less.
  !
  !    note.. odrv and its subordinate routines have been modified to
  !    compute orderings for general matrices, not necessarily having any
  !    symmetry.  the miminum degree ordering is computed for the
  !    structure of the symmetric matrix  m + m-transpose.
  !    modifications to the original odrv module have been made in
  !    the coding in subroutine mdi, and in the initial comments in
  !    subroutines odrv and md.
  !
  !    if only the nonzero entries in the upper triangle of m are being
  !    stored, then odrv symmetrically reorders (ia,ja,a), (optionally)
  !    with the diagonal entries placed first in each row.  this is to
  !    ensure that if m(i,j) will be in the upper triangle of m with
  !    respect to the new ordering, then m(i,j) is stored in row i (and
  !    thus m(j,i) is not stored),  whereas if m(i,j) will be in the
  !    strict lower triangle of m, then m(j,i) is stored in row j (and
  !    thus m(i,j) is not stored).
  !
  !
  !  storage of sparse matrices
  !
  !    the nonzero entries of the matrix m are stored row-by-row in the
  !    array a.  to identify the individual nonzero entries in each row,
  !    we need to know in which column each entry lies.  these column
  !    indices are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
  !    ja(k) = j.  to identify the individual rows, we need to know where
  !    each row starts.  these row pointers are stored in the array ia.
  !    i.e., if m(i,j) is the first nonzero entry (stored) in the i-th row
  !    and  a(k) = m(i,j),  then  ia(i) = k.  moreover, ia(n+1) points to
  !    the first location following the last element in the last row.
  !    thus, the number of entries in the i-th row is  ia(i+1) - ia(i),
  !    the nonzero entries in the i-th row are stored consecutively in
  !
  !            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
  !
  !    and the corresponding column indices are stored consecutively in
  !
  !            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
  !
  !    when the coefficient matrix is symmetric, only the nonzero entries
  !    in the upper triangle need be stored.  for example, the matrix
  !
  !             ( 1  0  2  3  0 )
  !             ( 0  4  0  0  0 )
  !         m = ( 2  0  5  6  0 )
  !             ( 3  0  6  7  8 )
  !             ( 0  0  0  8  9 )
  !
  !    could be stored as
  !
  !            - 1  2  3  4  5  6  7  8  9 10 11 12 13
  !         ---+--------------------------------------
  !         ia - 1  4  5  8 12 14
  !         ja - 1  3  4  2  1  3  4  1  3  4  5  4  5
  !          a - 1  2  3  4  2  5  6  3  6  7  8  8  9
  !
  !    or (symmetrically) as
  !
  !            - 1  2  3  4  5  6  7  8  9
  !         ---+--------------------------
  !         ia - 1  4  5  7  9 10
  !         ja - 1  3  4  2  3  4  4  5  5
  !          a - 1  2  3  4  5  6  7  8  9          .
  !
  !
  !  parameters
  !
  !    n    - order of the matrix
  !
  !    ia   - integer one-dimensional array containing pointers to delimit
  !           rows in ja and a.  dimension = n+1
  !
  !    ja   - integer one-dimensional array containing the column indices
  !           corresponding to the elements of a.  dimension = number of
  !           nonzero entries in (the upper triangle of) m
  !
  !    a    - real one-dimensional array containing the nonzero entries in
  !           (the upper triangle of) m, stored by rows.  dimension =
  !           number of nonzero entries in (the upper triangle of) m
  !
  !    p    - integer one-dimensional array used to return the permutation
  !           of the rows and columns of m corresponding to the minimum
  !           degree ordering.  dimension = n
  !
  !    ip   - integer one-dimensional array used to return the inverse of
  !           the permutation returned in p.  dimension = n
  !
  !    nsp  - declared dimension of the one-dimensional array isp.  nsp
  !           must be at least  3n+4k,  where k is the number of nonzeroes
  !           in the strict upper triangle of m
  !
  !    isp  - integer one-dimensional array used for working storage.
  !           dimension = nsp
  !
  !    path - integer path specification.  values and their meanings are -
  !             1  find minimum degree ordering only
  !             2  find minimum degree ordering and reorder symmetrically
  !                  stored matrix (used when only the nonzero entries in
  !                  the upper triangle of m are being stored)
  !             3  reorder symmetrically stored matrix as specified by
  !                  input permutation (used when an ordering has already
  !                  been determined and only the nonzero entries in the
  !                  upper triangle of m are being stored)
  !             4  same as 2 but put diagonal entries at start of each row
  !             5  same as 3 but put diagonal entries at start of each row
  !
  !    flag - integer error flag.  values and their meanings are -
  !               0    no errors detected
  !              9n+k  insufficient storage in md
  !             10n+1  insufficient storage in odrv
  !             11n+1  illegal path specification
  !
  !
  !  conversion from real to double precision
  !
  !    change the real declarations in odrv and sro to double precision
  !    declarations.
  !
  !-----------------------------------------------------------------------
  !
        INTEGER Ia(*) , Ja(*) , P(*) , Ip(*) , Isp(*) , Path , Flag , v , &
       &        l , head , tmp , q
  !...  real  a(*)
        DOUBLE PRECISION A(*)
        LOGICAL dflag
  !
  !----initialize error flag and validate path specification
        Flag = 0
        IF ( Path.LT.1 .OR. 5.LT.Path ) THEN
  ! ** error -- illegal path specified
           Flag = 11*N + 1
           GOTO 99999
        ELSE
  !
  !----allocate storage and find minimum degree ordering
           IF ( (Path-1)*(Path-2)*(Path-4).EQ.0 ) THEN
              max = (Nsp-N)/2
              v = 1
              l = v + max
              head = l + max
              next = head + N
              IF ( max.LT.N ) GOTO 100
  !
              CALL MD(N,Ia,Ja,max,Isp(v),Isp(l),Isp(head),P,Ip,Isp(v),    &
       &              Flag)
  !
  ! ** error -- error detected in md
              IF ( Flag.NE.0 ) RETURN
           ENDIF
  !
  !----allocate storage and symmetrically reorder matrix
           IF ( (Path-2)*(Path-3)*(Path-4)*(Path-5).EQ.0 ) THEN
              tmp = (Nsp+1) - N
              q = tmp - (Ia(N+1)-1)
              IF ( q.LT.1 ) GOTO 100
  !
              dflag = Path.EQ.4 .OR. Path.EQ.5
              CALL SRO(N,Ip,Ia,Ja,A,Isp(tmp),Isp(q),dflag)
           ENDIF
  !
           RETURN
        ENDIF
  ! ** error -- insufficient storage
   100  Flag = 10*N + 1
        RETURN
  99999 END
  !*==MD.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
        SUBROUTINE MD(N,Ia,Ja,Max,V,L,Head,Last,Next,Mark,Flag)
        IMPLICIT NONE
  !*--MD2389
  !*** Start of declarations inserted by SPAG
        INTEGER k , Max , N
  !*** End of declarations inserted by SPAG
  !***********************************************************************
  !  md -- minimum degree algorithm (based on element model)
  !***********************************************************************
  !
  !  description
  !
  !    md finds a minimum degree ordering of the rows and columns of a
  !    general sparse matrix m stored in (ia,ja,a) format.
  !    when the structure of m is nonsymmetric, the ordering is that
  !    obtained for the symmetric matrix  m + m-transpose.
  !
  !
  !  additional parameters
  !
  !    max  - declared dimension of the one-dimensional arrays v and l.
  !           max must be at least  n+2k,  where k is the number of
  !           nonzeroes in the strict upper triangle of m + m-transpose
  !
  !    v    - integer one-dimensional work array.  dimension = max
  !
  !    l    - integer one-dimensional work array.  dimension = max
  !
  !    head - integer one-dimensional work array.  dimension = n
  !
  !    last - integer one-dimensional array used to return the permutation
  !           of the rows and columns of m corresponding to the minimum
  !           degree ordering.  dimension = n
  !
  !    next - integer one-dimensional array used to return the inverse of
  !           the permutation returned in last.  dimension = n
  !
  !    mark - integer one-dimensional work array (may be the same as v).
  !           dimension = n
  !
  !    flag - integer error flag.  values and their meanings are -
  !             0     no errors detected
  !             9n+k  insufficient storage in md
  !
  !
  !  definitions of internal parameters
  !
  !    ---------+---------------------------------------------------------
  !    v(s)     - value field of list entry
  !    ---------+---------------------------------------------------------
  !    l(s)     - link field of list entry  (0 =) end of list)
  !    ---------+---------------------------------------------------------
  !    l(vi)    - pointer to element list of uneliminated vertex vi
  !    ---------+---------------------------------------------------------
  !    l(ej)    - pointer to boundary list of active element ej
  !    ---------+---------------------------------------------------------
  !    head(d)  - vj =) vj head of d-list d
  !             -  0 =) no vertex in d-list d
  !
  !
  !             -                  vi uneliminated vertex
  !             -          vi in ek           -       vi not in ek
  !    ---------+-----------------------------+---------------------------
  !    next(vi) - undefined but nonnegative   - vj =) vj next in d-list
  !             -                             -  0 =) vi tail of d-list
  !    ---------+-----------------------------+---------------------------
  !    last(vi) - (not set until mdp)         - -d =) vi head of d-list d
  !             --vk =) compute degree        - vj =) vj last in d-list
  !             - ej =) vi prototype of ej    -  0 =) vi not in any d-list
  !             -  0 =) do not compute degree -
  !    ---------+-----------------------------+---------------------------
  !    mark(vi) - mark(vk)                    - nonneg. tag .lt. mark(vk)
  !
  !
  !             -                   vi eliminated vertex
  !             -      ei active element      -           otherwise
  !    ---------+-----------------------------+---------------------------
  !    next(vi) - -j =) vi was j-th vertex    - -j =) vi was j-th vertex
  !             -       to be eliminated      -       to be eliminated
  !    ---------+-----------------------------+---------------------------
  !    last(vi) -  m =) size of ei = m        - undefined
  !    ---------+-----------------------------+---------------------------
  !    mark(vi) - -m =) overlap count of ei   - undefined
  !             -       with ek = m           -
  !             - otherwise nonnegative tag   -
  !             -       .lt. mark(vk)         -
  !
  !-----------------------------------------------------------------------
  !
        INTEGER Ia(*) , Ja(*) , V(*) , L(*) , Head(*) , Last(*) , Next(*) &
       &        , Mark(*) , Flag , tag , dmin , vk , ek , tail
        EQUIVALENCE (vk,ek)
  !
  !----initialization
        tag = 0
        CALL MDI(N,Ia,Ja,Max,V,L,Head,Last,Next,Mark,tag,Flag)
        IF ( Flag.NE.0 ) RETURN
  !
        k = 0
        dmin = 1
  !
  !----while  k .lt. n  do
   100  IF ( k.GE.N ) THEN
  !
  !----generate inverse permutation from permutation
           DO k = 1 , N
              Next(k) = -Next(k)
              Last(Next(k)) = k
           ENDDO
        ELSE
  !
  !------search for vertex of minimum degree
   150     IF ( Head(dmin).GT.0 ) THEN
  !
  !------remove vertex vk of minimum degree from degree list
              vk = Head(dmin)
              Head(dmin) = Next(vk)
              IF ( Head(dmin).GT.0 ) Last(Head(dmin)) = -dmin
  !
  !------number vertex vk, adjust tag, and tag vk
              k = k + 1
              Next(vk) = -k
              Last(ek) = dmin - 1
              tag = tag + Last(ek)
              Mark(vk) = tag
  !
  !------form element ek from uneliminated neighbors of vk
              CALL MDM(vk,tail,V,L,Last,Next,Mark)
  !
  !------purge inactive elements and do mass elimination
              CALL MDP(k,ek,tail,V,L,Head,Last,Next,Mark)
  !
  !------update degrees of uneliminated vertices in ek
              CALL MDU(ek,dmin,V,L,Head,Last,Next,Mark)
  !
              GOTO 100
           ELSE
              dmin = dmin + 1
              GOTO 150
           ENDIF
        ENDIF
  !
        END
  !*==MDI.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
        SUBROUTINE MDI(N,Ia,Ja,Max,V,L,Head,Last,Next,Mark,Tag,Flag)
        IMPLICIT NONE
  !*--MDI2533
  !*** Start of declarations inserted by SPAG
        INTEGER j , jmax , jmin , k , kmax , lvk , Max , N , nextvi
  !*** End of declarations inserted by SPAG
  !***********************************************************************
  !  mdi -- initialization
  !***********************************************************************
        INTEGER Ia(*) , Ja(*) , V(*) , L(*) , Head(*) , Last(*) , Next(*) &
       &        , Mark(*) , Tag , Flag , sfs , vi , dvi , vj
  !
  !----initialize degrees, element lists, and degree lists
        DO vi = 1 , N
           Mark(vi) = 1
           L(vi) = 0
           Head(vi) = 0
        ENDDO
        sfs = N + 1
  !
  !----create nonzero structure
  !----for each nonzero entry a(vi,vj)
        DO vi = 1 , N
           jmin = Ia(vi)
           jmax = Ia(vi+1) - 1
           IF ( jmin.LE.jmax ) THEN
              DO j = jmin , jmax
                 vj = Ja(j)
                 IF ( vj.LT.vi ) THEN
  !
  !------if a(vi,vj) is in strict lower triangle
  !------check for previous occurrence of a(vj,vi)
                    lvk = vi
                    kmax = Mark(vi) - 1
                    IF ( kmax.NE.0 ) THEN
                       DO k = 1 , kmax
                          lvk = L(lvk)
                          IF ( V(lvk).EQ.vj ) GOTO 20
                       ENDDO
                    ENDIF
                 ELSEIF ( vj.EQ.vi ) THEN
                    GOTO 20
                 ENDIF
  !----for unentered entries a(vi,vj)
                 IF ( sfs.GE.Max ) GOTO 100
  !
  !------enter vj in element list for vi
                 Mark(vi) = Mark(vi) + 1
                 V(sfs) = vj
                 L(sfs) = L(vi)
                 L(vi) = sfs
                 sfs = sfs + 1
  !
  !------enter vi in element list for vj
                 Mark(vj) = Mark(vj) + 1
                 V(sfs) = vi
                 L(sfs) = L(vj)
                 L(vj) = sfs
                 sfs = sfs + 1
   20         ENDDO
           ENDIF
        ENDDO
  !
  !----create degree lists and initialize mark vector
        DO vi = 1 , N
           dvi = Mark(vi)
           Next(vi) = Head(dvi)
           Head(dvi) = vi
           Last(vi) = -dvi
           nextvi = Next(vi)
           IF ( nextvi.GT.0 ) Last(nextvi) = vi
           Mark(vi) = Tag
        ENDDO
  !
        RETURN
  !
  ! ** error-  insufficient storage
   100  Flag = 9*N + vi
        END
  !*==MDM.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
        SUBROUTINE MDM(Vk,Tail,V,L,Last,Next,Mark)
        IMPLICIT NONE
  !*--MDM2613
  !***********************************************************************
  !  mdm -- form element from uneliminated neighbors of vk
  !***********************************************************************
        INTEGER Vk , Tail , V(*) , L(*) , Last(*) , Next(*) , Mark(*) ,   &
       &        tag , s , ls , vs , es , b , lb , vb , blp , blpmax
        EQUIVALENCE (vs,es)
  !
  !----initialize tag and list of uneliminated neighbors
        tag = Mark(Vk)
        Tail = Vk
  !
  !----for each vertex/element vs/es in element list of vk
        ls = L(Vk)
   100  s = ls
        IF ( s.EQ.0 ) THEN
  !
  !----terminate list of uneliminated neighbors
           L(Tail) = 0
        ELSE
           ls = L(s)
           vs = V(s)
           IF ( Next(vs).LT.0 ) THEN
  !
  !------if es is active element, then ...
  !--------for each vertex vb in boundary list of element es
              lb = L(es)
              blpmax = Last(es)
              DO blp = 1 , blpmax
                 b = lb
                 lb = L(b)
                 vb = V(b)
  !
  !----------if vb is untagged vertex, then tag and append to list of
  !----------uneliminated neighbors
                 IF ( Mark(vb).LT.tag ) THEN
                    Mark(vb) = tag
                    L(Tail) = b
                    Tail = b
                 ENDIF
              ENDDO
  !
  !--------mark es inactive
  !
              Mark(es) = tag
           ELSE
  !
  !------if vs is uneliminated vertex, then tag and append to list of
  !------uneliminated neighbors
              Mark(vs) = tag
              L(Tail) = s
              Tail = s
           ENDIF
           GOTO 100
        ENDIF
  !
        END
  !*==MDP.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
        SUBROUTINE MDP(K,Ek,Tail,V,L,Head,Last,Next,Mark)
        IMPLICIT NONE
  !*--MDP2673
  !*** Start of declarations inserted by SPAG
        INTEGER i , K
  !*** End of declarations inserted by SPAG
  !***********************************************************************
  !  mdp -- purge inactive elements and do mass elimination
  !***********************************************************************
        INTEGER Ek , Tail , V(*) , L(*) , Head(*) , Last(*) , Next(*) ,   &
       &        Mark(*) , tag , free , li , vi , lvi , evi , s , ls , es ,&
       &        ilp , ilpmax
  !
  !----initialize tag
        tag = Mark(Ek)
  !
  !----for each vertex vi in ek
        li = Ek
        ilpmax = Last(Ek)
        IF ( ilpmax.GT.0 ) THEN
           DO ilp = 1 , ilpmax
              i = li
              li = L(i)
              vi = V(li)
  !
  !------remove vi from degree list
              IF ( Last(vi).NE.0 ) THEN
                 IF ( Last(vi).GT.0 ) THEN
                    Next(Last(vi)) = Next(vi)
                 ELSE
                    Head(-Last(vi)) = Next(vi)
                 ENDIF
                 IF ( Next(vi).GT.0 ) Last(Next(vi)) = Last(vi)
              ENDIF
  !
  !------remove inactive items from element list of vi
              ls = vi
   20         s = ls
              ls = L(s)
              IF ( ls.EQ.0 ) THEN
  !
  !------if vi is interior vertex, then remove from list and eliminate
                 lvi = L(vi)
                 IF ( lvi.NE.0 ) THEN
  !
  !------else ...
  !--------classify vertex vi
                    IF ( L(lvi).NE.0 ) THEN
  !
  !----------else mark vi to compute degree
                       Last(vi) = -Ek
                    ELSE
                       evi = V(lvi)
                       IF ( Next(evi).GE.0 ) THEN
                          Last(vi) = -Ek
                       ELSEIF ( Mark(evi).LT.0 ) THEN
  !
  !----------else if vi is duplicate vertex, then mark as such and adjust
  !----------overlap count for corresponding element
                          Last(vi) = 0
                          Mark(evi) = Mark(evi) - 1
                       ELSE
  !
  !----------if vi is prototype vertex, then mark as such, initialize
  !----------overlap count for corresponding element, and move vi to end
  !----------of boundary list
                          Last(vi) = evi
                          Mark(evi) = -1
                          L(Tail) = li
                          Tail = li
                          L(i) = L(li)
                          li = i
                       ENDIF
                    ENDIF
  !
  !--------insert ek in element list of vi
                    V(free) = Ek
                    L(free) = L(vi)
                    L(vi) = free
                 ELSE
                    L(i) = L(li)
                    li = i
  !
                    K = K + 1
                    Next(vi) = -K
                    Last(Ek) = Last(Ek) - 1
                 ENDIF
              ELSE
                 es = V(ls)
                 IF ( Mark(es).GE.tag ) THEN
                    free = ls
                    L(s) = L(ls)
                    ls = s
                 ENDIF
                 GOTO 20
              ENDIF
           ENDDO
        ENDIF
  !
  !----terminate boundary list
        L(Tail) = 0
  !
        END
  !*==MDU.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
        SUBROUTINE MDU(Ek,Dmin,V,L,Head,Last,Next,Mark)
        IMPLICIT NONE
  !*--MDU2777
  !*** Start of declarations inserted by SPAG
        INTEGER i
  !*** End of declarations inserted by SPAG
  !***********************************************************************
  !  mdu -- update degrees of uneliminated vertices in ek
  !***********************************************************************
        INTEGER Ek , Dmin , V(*) , L(*) , Head(*) , Last(*) , Next(*) ,   &
       &        Mark(*) , tag , vi , evi , dvi , s , vs , es , b , vb ,   &
       &        ilp , ilpmax , blp , blpmax
        EQUIVALENCE (vs,es)
  !
  !----initialize tag
        tag = Mark(Ek) - Last(Ek)
  !
  !----for each vertex vi in ek
        i = Ek
        ilpmax = Last(Ek)
        IF ( ilpmax.GT.0 ) THEN
           DO ilp = 1 , ilpmax
              i = L(i)
              vi = V(i)
              IF ( Last(vi).LT.0 ) THEN
  !
  !------if vi neither prototype nor duplicate vertex, then merge elements
  !------to compute degree
                 tag = tag + 1
                 dvi = Last(Ek)
  !
  !--------for each vertex/element vs/es in element list of vi
                 s = L(vi)
   10            s = L(s)
                 IF ( s.NE.0 ) THEN
                    vs = V(s)
                    IF ( Next(vs).GE.0 ) THEN
  !
  !----------if vs is uneliminated vertex, then tag and adjust degree
                       Mark(vs) = tag
                       dvi = dvi + 1
                       GOTO 10
  !
  !----------if es is active element, then expand
  !------------check for outmatched vertex
                    ELSEIF ( Mark(es).LT.0 ) THEN
  !
  !------else if vi is outmatched vertex, then adjust overlaps but do not
  !------compute degree
                       Last(vi) = 0
                       Mark(es) = Mark(es) - 1
   12                  s = L(s)
                       IF ( s.EQ.0 ) GOTO 50
                       es = V(s)
                       IF ( Mark(es).LT.0 ) Mark(es) = Mark(es) - 1
                       GOTO 12
                    ELSE
  !
  !------------for each vertex vb in es
                       b = es
                       blpmax = Last(es)
                       DO blp = 1 , blpmax
                          b = L(b)
                          vb = V(b)
  !
  !--------------if vb is untagged, then tag and adjust degree
                          IF ( Mark(vb).LT.tag ) THEN
                             Mark(vb) = tag
                             dvi = dvi + 1
                          ENDIF
                       ENDDO
  !
                       GOTO 10
                    ENDIF
                 ENDIF
              ELSEIF ( Last(vi).EQ.0 ) THEN
                 GOTO 50
              ELSE
  !
  !------else if vi is prototype vertex, then calculate degree by
  !------inclusion/exclusion and reset overlap count
                 evi = Last(vi)
                 dvi = Last(Ek) + Last(evi) + Mark(evi)
                 Mark(evi) = 0
              ENDIF
  !
  !------insert vi in appropriate degree list
              Next(vi) = Head(dvi)
              Head(dvi) = vi
              Last(vi) = -dvi
              IF ( Next(vi).GT.0 ) Last(Next(vi)) = vi
              IF ( dvi.LT.Dmin ) Dmin = dvi
  !
   50      ENDDO
        ENDIF
  !
        END
  !*==SRO.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
        SUBROUTINE SRO(N,Ip,Ia,Ja,A,Q,R,Dflag)
        IMPLICIT NONE
  !*--SRO2875
  !*** Start of declarations inserted by SPAG
        INTEGER i , ilast , j , jak , jdummy , jmax , jmin , k , N
  !*** End of declarations inserted by SPAG
  !***********************************************************************
  !  sro -- symmetric reordering of sparse symmetric matrix
  !***********************************************************************
  !
  !  description
  !
  !    the nonzero entries of the matrix m are assumed to be stored
  !    symmetrically in (ia,ja,a) format (i.e., not both m(i,j) and m(j,i)
  !    are stored if i ne j).
  !
  !    sro does not rearrange the order of the rows, but does move
  !    nonzeroes from one row to another to ensure that if m(i,j) will be
  !    in the upper triangle of m with respect to the new ordering, then
  !    m(i,j) is stored in row i (and thus m(j,i) is not stored),  whereas
  !    if m(i,j) will be in the strict lower triangle of m, then m(j,i) is
  !    stored in row j (and thus m(i,j) is not stored).
  !
  !
  !  additional parameters
  !
  !    q     - integer one-dimensional work array.  dimension = n
  !
  !    r     - integer one-dimensional work array.  dimension = number of
  !            nonzero entries in the upper triangle of m
  !
  !    dflag - logical variable.  if dflag = .true., then store nonzero
  !            diagonal elements at the beginning of the row
  !
  !-----------------------------------------------------------------------
  !
        INTEGER Ip(*) , Ia(*) , Ja(*) , Q(*) , R(*)
  !...  real  a(*),  ak
        DOUBLE PRECISION A(*) , ak
        LOGICAL Dflag
  !
  !
  !--phase 1 -- find row in which to store each nonzero
  !----initialize count of nonzeroes to be stored in each row
        DO i = 1 , N
           Q(i) = 0
        ENDDO
  !
  !----for each nonzero element a(j)
        DO i = 1 , N
           jmin = Ia(i)
           jmax = Ia(i+1) - 1
           IF ( jmin.LE.jmax ) THEN
              DO j = jmin , jmax
  !
  !--------find row (=r(j)) and column (=ja(j)) in which to store a(j) ...
                 k = Ja(j)
                 IF ( Ip(k).LT.Ip(i) ) Ja(j) = i
                 IF ( Ip(k).GE.Ip(i) ) k = i
                 R(j) = k
  !
  !--------... and increment count of nonzeroes (=q(r(j)) in that row
                 Q(k) = Q(k) + 1
              ENDDO
           ENDIF
        ENDDO
  !
  !
  !--phase 2 -- find new ia and permutation to apply to (ja,a)
  !----determine pointers to delimit rows in permuted (ja,a)
        DO i = 1 , N
           Ia(i+1) = Ia(i) + Q(i)
           Q(i) = Ia(i+1)
        ENDDO
  !
  !----determine where each (ja(j),a(j)) is stored in permuted (ja,a)
  !----for each nonzero element (in reverse order)
        ilast = 0
        jmin = Ia(1)
        jmax = Ia(N+1) - 1
        j = jmax
        DO jdummy = jmin , jmax
           i = R(j)
           IF ( .NOT.Dflag .OR. Ja(j).NE.i .OR. i.EQ.ilast ) THEN
  !
  !------put (off-diagonal) nonzero in last unused location in row
              Q(i) = Q(i) - 1
              R(j) = Q(i)
           ELSE
  !
  !------if dflag, then put diagonal nonzero at beginning of row
              R(j) = Ia(i)
              ilast = i
           ENDIF
  !
           j = j - 1
        ENDDO
  !
  !
  !--phase 3 -- permute (ja,a) to upper triangular form (wrt new ordering)
        DO j = jmin , jmax
   50      IF ( R(j).NE.j ) THEN
              k = R(j)
              R(j) = R(k)
              R(k) = k
              jak = Ja(k)
              Ja(k) = Ja(j)
              Ja(j) = jak
              ak = A(k)
              A(k) = A(j)
              A(j) = ak
              GOTO 50
           ENDIF
        ENDDO
  !
        END
  !*==CDRV.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK CDRV
        SUBROUTINE CDRV(N,R,C,Ic,Ia,Ja,A,B,Z,Nsp,Isp,Rsp,Esp,Path,Flag)
        IMPLICIT NONE
  !*--CDRV2993
  !*** Start of declarations inserted by SPAG
        INTEGER i , ijl , iju , il , ira , irac , irl , iru , iu , j ,    &
       &        jl , jlmax , jra , jrl , jru , ju , jumax , jutmp , l ,   &
       &        lmax
        INTEGER lratio , max , N , Nsp
  !*** End of declarations inserted by SPAG
  !*** subroutine cdrv
  !*** driver for subroutines for solving sparse nonsymmetric systems of
  !       linear equations (compressed pointer storage)
  !
  !
  !    parameters
  !    class abbreviations are--
  !       n - integer variable
  !       f - real variable
  !       v - supplies a value to the driver
  !       r - returns a result from the driver
  !       i - used internally by the driver
  !       a - array
  !
  ! class - parameter
  ! ------+----------
  !       -
  !         the nonzero entries of the coefficient matrix m are stored
  !    row-by-row in the array a.  to identify the individual nonzero
  !    entries in each row, we need to know in which column each entry
  !    lies.  the column indices which correspond to the nonzero entries
  !    of m are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
  !    ja(k) = j.  in addition, we need to know where each row starts and
  !    how long it is.  the index positions in ja and a where the rows of
  !    m begin are stored in the array ia.  i.e., if m(i,j) is the first
  !    nonzero entry (stored) in the i-th row and a(k) = m(i,j),  then
  !    ia(i) = k.  moreover, the index in ja and a of the first location
  !    following the last element in the last row is stored in ia(n+1).
  !    thus, the number of entries in the i-th row is given by
  !    ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
  !    consecutively in
  !            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
  !    and the corresponding column indices are stored consecutively in
  !            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
  !    for example, the 5 by 5 matrix
  !                ( 1. 0. 2. 0. 0.)
  !                ( 0. 3. 0. 0. 0.)
  !            m = ( 0. 4. 5. 6. 0.)
  !                ( 0. 0. 0. 7. 0.)
  !                ( 0. 0. 0. 8. 9.)
  !    would be stored as
  !               - 1  2  3  4  5  6  7  8  9
  !            ---+--------------------------
  !            ia - 1  3  4  7  8 10
  !            ja - 1  3  2  2  3  4  4  4  5
  !             a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
  !
  ! nv    - n     - number of variables/equations.
  ! fva   - a     - nonzero entries of the coefficient matrix m, stored
  !       -           by rows.
  !       -           size = number of nonzero entries in m.
  ! nva   - ia    - pointers to delimit the rows in a.
  !       -           size = n+1.
  ! nva   - ja    - column numbers corresponding to the elements of a.
  !       -           size = size of a.
  ! fva   - b     - right-hand side b.  b and z can the same array.
  !       -           size = n.
  ! fra   - z     - solution x.  b and z can be the same array.
  !       -           size = n.
  !
  !         the rows and columns of the original matrix m can be
  !    reordered (e.g., to reduce fillin or ensure numerical stability)
  !    before calling the driver.  if no reordering is done, then set
  !    r(i) = c(i) = ic(i) = i  for i=1,...,n.  the solution z is returned
  !    in the original order.
  !         if the columns have been reordered (i.e.,  c(i).ne.i  for some
  !    i), then the driver will call a subroutine (nroc) which rearranges
  !    each row of ja and a, leaving the rows in the original order, but
  !    placing the elements of each row in increasing order with respect
  !    to the new ordering.  if  path.ne.1,  then nroc is assumed to have
  !    been called already.
  !
  ! nva   - r     - ordering of the rows of m.
  !       -           size = n.
  ! nva   - c     - ordering of the columns of m.
  !       -           size = n.
  ! nva   - ic    - inverse of the ordering of the columns of m.  i.e.,
  !       -           ic(c(i)) = i  for i=1,...,n.
  !       -           size = n.
  !
  !         the solution of the system of linear equations is divided into
  !    three stages --
  !      nsfc -- the matrix m is processed symbolically to determine where
  !               fillin will occur during the numeric factorization.
  !      nnfc -- the matrix m is factored numerically into the product ldu
  !               of a unit lower triangular matrix l, a diagonal matrix
  !               d, and a unit upper triangular matrix u, and the system
  !               mx = b  is solved.
  !      nnsc -- the linear system  mx = b  is solved using the ldu
  !  or           factorization from nnfc.
  !      nntc -- the transposed linear system  mt x = b  is solved using
  !               the ldu factorization from nnf.
  !    for several systems whose coefficient matrices have the same
  !    nonzero structure, nsfc need be done only once (for the first
  !    system).  then nnfc is done once for each additional system.  for
  !    several systems with the same coefficient matrix, nsfc and nnfc
  !    need be done only once (for the first system).  then nnsc or nntc
  !    is done once for each additional right-hand side.
  !
  ! nv    - path  - path specification.  values and their meanings are --
  !       -           1  perform nroc, nsfc, and nnfc.
  !       -           2  perform nnfc only  (nsfc is assumed to have been
  !       -               done in a manner compatible with the storage
  !       -               allocation used in the driver).
  !       -           3  perform nnsc only  (nsfc and nnfc are assumed to
  !       -               have been done in a manner compatible with the
  !       -               storage allocation used in the driver).
  !       -           4  perform nntc only  (nsfc and nnfc are assumed to
  !       -               have been done in a manner compatible with the
  !       -               storage allocation used in the driver).
  !       -           5  perform nroc and nsfc.
  !
  !         various errors are detected by the driver and the individual
  !    subroutines.
  !
  ! nr    - flag  - error flag.  values and their meanings are --
  !       -             0     no errors detected
  !       -             n+k   null row in a  --  row = k
  !       -            2n+k   duplicate entry in a  --  row = k
  !       -            3n+k   insufficient storage in nsfc  --  row = k
  !       -            4n+1   insufficient storage in nnfc
  !       -            5n+k   null pivot  --  row = k
  !       -            6n+k   insufficient storage in nsfc  --  row = k
  !       -            7n+1   insufficient storage in nnfc
  !       -            8n+k   zero pivot  --  row = k
  !       -           10n+1   insufficient storage in cdrv
  !       -           11n+1   illegal path specification
  !
  !         working storage is needed for the factored form of the matrix
  !    m plus various temporary vectors.  the arrays isp and rsp should be
  !    equivalenced.  integer storage is allocated from the beginning of
  !    isp and real storage from the end of rsp.
  !
  ! nv    - nsp   - declared dimension of rsp.  nsp generally must
  !       -           be larger than  8n+2 + 2k  (where  k = (number of
  !       -           nonzero entries in m)).
  ! nvira - isp   - integer working storage divided up into various arrays
  !       -           needed by the subroutines.  isp and rsp should be
  !       -           equivalenced.
  !       -           size = lratio*nsp.
  ! fvira - rsp   - real working storage divided up into various arrays
  !       -           needed by the subroutines.  isp and rsp should be
  !       -           equivalenced.
  !       -           size = nsp.
  ! nr    - esp   - if sufficient storage was available to perform the
  !       -           symbolic factorization (nsfc), then esp is set to
  !       -           the amount of excess storage provided (negative if
  !       -           insufficient storage was available to perform the
  !       -           numeric factorization (nnfc)).
  !
  !
  !  conversion to double precision
  !
  !    to convert these routines for double precision arrays..
  !    (1) use the double precision declarations in place of the real
  !    declarations in each subprogram, as given in comment cards.
  !    (2) change the data-loaded value of the integer  lratio
  !    in subroutine cdrv, as indicated below.
  !    (3) change e0 to d0 in the constants in statement number 10
  !    in subroutine nnfc and the line following that.
  !
        INTEGER R(*) , C(*) , Ic(*) , Ia(*) , Ja(*) , Isp(*) , Esp ,      &
       &        Path , Flag , d , u , q , row , tmp , ar , umax
  !     real  a(*), b(*), z(*), rsp(*)
        DOUBLE PRECISION A(*) , B(*) , Z(*) , Rsp(*)
  !
  !  set lratio equal to the ratio between the length of floating point
  !  and integer array data.  e. g., lratio = 1 for (real, integer),
  !  lratio = 2 for (double precision, integer)
  !
        DATA lratio/2/
  !
        IF ( Path.LT.1 .OR. 5.LT.Path ) THEN
  ! ** error.. illegal path specification
           Flag = 11*N + 1
           GOTO 99999
        ELSE
  !******initialize and divide up temporary storage  *******************
           il = 1
           ijl = il + (N+1)
           iu = ijl + N
           iju = iu + (N+1)
           irl = iju + N
           jrl = irl + N
           jl = jrl + N
  !
  !  ******  reorder a if necessary, call nsfc if flag is set  ***********
           IF ( (Path-1)*(Path-5).NE.0 ) GOTO 200
           max = (lratio*Nsp+1-jl) - (N+1) - 5*N
           jlmax = max/2
           q = jl + jlmax
           ira = q + (N+1)
           jra = ira + N
           irac = jra + N
           iru = irac + N
           jru = iru + N
           jutmp = jru + N
           jumax = lratio*Nsp + 1 - jutmp
           Esp = max/lratio
           IF ( jlmax.LE.0 .OR. jumax.LE.0 ) GOTO 400
  !
           DO i = 1 , N
              IF ( C(i).NE.i ) GOTO 50
           ENDDO
           GOTO 100
   50      ar = Nsp + 1 - N
           CALL NROC(N,Ic,Ia,Ja,A,Isp(il),Rsp(ar),Isp(iu),Flag)
           IF ( Flag.NE.0 ) GOTO 300
  !
   100     CALL NSFC(N,R,Ic,Ia,Ja,jlmax,Isp(il),Isp(jl),Isp(ijl),jumax,   &
       &             Isp(iu),Isp(jutmp),Isp(iju),Isp(q),Isp(ira),Isp(jra),&
       &             Isp(irac),Isp(irl),Isp(jrl),Isp(iru),Isp(jru),Flag)
           IF ( Flag.NE.0 ) GOTO 300
  !  ******  move ju next to jl  *****************************************
           jlmax = Isp(ijl+N-1)
           ju = jl + jlmax
           jumax = Isp(iju+N-1)
           IF ( jumax.GT.0 ) THEN
              DO j = 1 , jumax
                 Isp(ju+j-1) = Isp(jutmp+j-1)
              ENDDO
           ENDIF
        ENDIF
  !
  !  ******  call remaining subroutines  *********************************
   200  jlmax = Isp(ijl+N-1)
        ju = jl + jlmax
        jumax = Isp(iju+N-1)
        l = (ju+jumax-2+lratio)/lratio + 1
        lmax = Isp(il+N) - 1
        d = l + lmax
        u = d + N
        row = Nsp + 1 - N
        tmp = row - N
        umax = tmp - u
        Esp = umax - (Isp(iu+N)-1)
  !
        IF ( (Path-1)*(Path-2).EQ.0 ) THEN
           IF ( umax.LT.0 ) GOTO 400
           CALL NNFC(N,R,C,Ic,Ia,Ja,A,Z,B,lmax,Isp(il),Isp(jl),Isp(ijl),  &
       &             Rsp(l),Rsp(d),umax,Isp(iu),Isp(ju),Isp(iju),Rsp(u),  &
       &             Rsp(row),Rsp(tmp),Isp(irl),Isp(jrl),Flag)
           IF ( Flag.NE.0 ) GOTO 300
        ENDIF
  !
        IF ( (Path-3).EQ.0 ) CALL NNSC(N,R,C,Isp(il),Isp(jl),Isp(ijl),    &
       &                               Rsp(l),Rsp(d),Isp(iu),Isp(ju),     &
       &                               Isp(iju),Rsp(u),Z,B,Rsp(tmp))
  !
        IF ( (Path-4).EQ.0 ) CALL NNTC(N,R,C,Isp(il),Isp(jl),Isp(ijl),    &
       &                               Rsp(l),Rsp(d),Isp(iu),Isp(ju),     &
       &                               Isp(iju),Rsp(u),Z,B,Rsp(tmp))
        RETURN
  !
  ! ** error.. error detected in nroc, nsfc, nnfc, or nnsc
   300  RETURN
  ! ** error.. insufficient storage
   400  Flag = 10*N + 1
        RETURN
  99999 END
  !*==NROC.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
        SUBROUTINE NROC(N,Ic,Ia,Ja,A,Jar,Ar,P,Flag)
        IMPLICIT NONE
  !*--NROC3263
  !*** Start of declarations inserted by SPAG
        INTEGER i , j , jmax , jmin , k , N , newj
  !*** End of declarations inserted by SPAG
  !
  !       ----------------------------------------------------------------
  !
  !               yale sparse matrix package - nonsymmetric codes
  !                    solving the system of equations mx = b
  !
  !    i.   calling sequences
  !         the coefficient matrix can be processed by an ordering routine
  !    (e.g., to reduce fillin or ensure numerical stability) before using
  !    the remaining subroutines.  if no reordering is done, then set
  !    r(i) = c(i) = ic(i) = i  for i=1,...,n.  if an ordering subroutine
  !    is used, then nroc should be used to reorder the coefficient matrix
  !    the calling sequence is --
  !        (       (matrix ordering))
  !        (nroc   (matrix reordering))
  !         nsfc   (symbolic factorization to determine where fillin will
  !                  occur during numeric factorization)
  !         nnfc   (numeric factorization into product ldu of unit lower
  !                  triangular matrix l, diagonal matrix d, and unit
  !                  upper triangular matrix u, and solution of linear
  !                  system)
  !         nnsc   (solution of linear system for additional right-hand
  !                  side using ldu factorization from nnfc)
  !    (if only one system of equations is to be solved, then the
  !    subroutine trk should be used.)
  !
  !    ii.  storage of sparse matrices
  !         the nonzero entries of the coefficient matrix m are stored
  !    row-by-row in the array a.  to identify the individual nonzero
  !    entries in each row, we need to know in which column each entry
  !    lies.  the column indices which correspond to the nonzero entries
  !    of m are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
  !    ja(k) = j.  in addition, we need to know where each row starts and
  !    how long it is.  the index positions in ja and a where the rows of
  !    m begin are stored in the array ia.  i.e., if m(i,j) is the first
  !    (leftmost) entry in the i-th row and  a(k) = m(i,j),  then
  !    ia(i) = k.  moreover, the index in ja and a of the first location
  !    following the last element in the last row is stored in ia(n+1).
  !    thus, the number of entries in the i-th row is given by
  !    ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
  !    consecutively in
  !            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
  !    and the corresponding column indices are stored consecutively in
  !            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
  !    for example, the 5 by 5 matrix
  !                ( 1. 0. 2. 0. 0.)
  !                ( 0. 3. 0. 0. 0.)
  !            m = ( 0. 4. 5. 6. 0.)
  !                ( 0. 0. 0. 7. 0.)
  !                ( 0. 0. 0. 8. 9.)
  !    would be stored as
  !               - 1  2  3  4  5  6  7  8  9
  !            ---+--------------------------
  !            ia - 1  3  4  7  8 10
  !            ja - 1  3  2  2  3  4  4  4  5
  !             a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
  !
  !         the strict upper (lower) triangular portion of the matrix
  !    u (l) is stored in a similar fashion using the arrays  iu, ju, u
  !    (il, jl, l)  except that an additional array iju (ijl) is used to
  !    compress storage of ju (jl) by allowing some sequences of column
  !    (row) indices to used for more than one row (column)  (n.b., l is
  !    stored by columns).  iju(k) (ijl(k)) points to the starting
  !    location in ju (jl) of entries for the kth row (column).
  !    compression in ju (jl) occurs in two ways.  first, if a row
  !    (column) i was merged into the current row (column) k, and the
  !    number of elements merged in from (the tail portion of) row
  !    (column) i is the same as the final length of row (column) k, then
  !    the kth row (column) and the tail of row (column) i are identical
  !    and iju(k) (ijl(k)) points to the start of the tail.  second, if
  !    some tail portion of the (k-1)st row (column) is identical to the
  !    head of the kth row (column), then iju(k) (ijl(k)) points to the
  !    start of that tail portion.  for example, the nonzero structure of
  !    the strict upper triangular part of the matrix
  !            d 0 x x x
  !            0 d 0 x x
  !            0 0 d x 0
  !            0 0 0 d x
  !            0 0 0 0 d
  !    would be represented as
  !                - 1 2 3 4 5 6
  !            ----+------------
  !             iu - 1 4 6 7 8 8
  !             ju - 3 4 5 4
  !            iju - 1 2 4 3           .
  !    the diagonal entries of l and u are assumed to be equal to one and
  !    are not stored.  the array d contains the reciprocals of the
  !    diagonal entries of the matrix d.
  !
  !    iii. additional storage savings
  !         in nsfc, r and ic can be the same array in the calling
  !    sequence if no reordering of the coefficient matrix has been done.
  !         in nnfc, r, c, and ic can all be the same array if no
  !    reordering has been done.  if only the rows have been reordered,
  !    then c and ic can be the same array.  if the row and column
  !    orderings are the same, then r and c can be the same array.  z and
  !    row can be the same array.
  !         in nnsc or nntc, r and c can be the same array if no
  !    reordering has been done or if the row and column orderings are the
  !    same.  z and b can be the same array.  however, then b will be
  !    destroyed.
  !
  !    iv.  parameters
  !         following is a list of parameters to the programs.  names are
  !    uniform among the various subroutines.  class abbreviations are --
  !       n - integer variable
  !       f - real variable
  !       v - supplies a value to a subroutine
  !       r - returns a result from a subroutine
  !       i - used internally by a subroutine
  !       a - array
  !
  ! class - parameter
  ! ------+----------
  ! fva   - a     - nonzero entries of the coefficient matrix m, stored
  !       -           by rows.
  !       -           size = number of nonzero entries in m.
  ! fva   - b     - right-hand side b.
  !       -           size = n.
  ! nva   - c     - ordering of the columns of m.
  !       -           size = n.
  ! fvra  - d     - reciprocals of the diagonal entries of the matrix d.
  !       -           size = n.
  ! nr    - flag  - error flag.  values and their meanings are --
  !       -            0     no errors detected
  !       -            n+k   null row in a  --  row = k
  !       -           2n+k   duplicate entry in a  --  row = k
  !       -           3n+k   insufficient storage for jl  --  row = k
  !       -           4n+1   insufficient storage for l
  !       -           5n+k   null pivot  --  row = k
  !       -           6n+k   insufficient storage for ju  --  row = k
  !       -           7n+1   insufficient storage for u
  !       -           8n+k   zero pivot  --  row = k
  ! nva   - ia    - pointers to delimit the rows of a.
  !       -           size = n+1.
  ! nvra  - ijl   - pointers to the first element in each column in jl,
  !       -           used to compress storage in jl.
  !       -           size = n.
  ! nvra  - iju   - pointers to the first element in each row in ju, used
  !       -           to compress storage in ju.
  !       -           size = n.
  ! nvra  - il    - pointers to delimit the columns of l.
  !       -           size = n+1.
  ! nvra  - iu    - pointers to delimit the rows of u.
  !       -           size = n+1.
  ! nva   - ja    - column numbers corresponding to the elements of a.
  !       -           size = size of a.
  ! nvra  - jl    - row numbers corresponding to the elements of l.
  !       -           size = jlmax.
  ! nv    - jlmax - declared dimension of jl.  jlmax must be larger than
  !       -           the number of nonzeros in the strict lower triangle
  !       -           of m plus fillin minus compression.
  ! nvra  - ju    - column numbers corresponding to the elements of u.
  !       -           size = jumax.
  ! nv    - jumax - declared dimension of ju.  jumax must be larger than
  !       -           the number of nonzeros in the strict upper triangle
  !       -           of m plus fillin minus compression.
  ! fvra  - l     - nonzero entries in the strict lower triangular portion
  !       -           of the matrix l, stored by columns.
  !       -           size = lmax.
  ! nv    - lmax  - declared dimension of l.  lmax must be larger than
  !       -           the number of nonzeros in the strict lower triangle
  !       -           of m plus fillin  (il(n+1)-1 after nsfc).
  ! nv    - n     - number of variables/equations.
  ! nva   - r     - ordering of the rows of m.
  !       -           size = n.
  ! fvra  - u     - nonzero entries in the strict upper triangular portion
  !       -           of the matrix u, stored by rows.
  !       -           size = umax.
  ! nv    - umax  - declared dimension of u.  umax must be larger than
  !       -           the number of nonzeros in the strict upper triangle
  !       -           of m plus fillin  (iu(n+1)-1 after nsfc).
  ! fra   - z     - solution x.
  !       -           size = n.
  !
  !       ----------------------------------------------------------------
  !
  !*** subroutine nroc
  !*** reorders rows of a, leaving row order unchanged
  !
  !
  !       input parameters.. n, ic, ia, ja, a
  !       output parameters.. ja, a, flag
  !
  !       parameters used internally..
  ! nia   - p     - at the kth step, p is a linked list of the reordered
  !       -           column indices of the kth row of a.  p(n+1) points
  !       -           to the first entry in the list.
  !       -           size = n+1.
  ! nia   - jar   - at the kth step,jar contains the elements of the
  !       -           reordered column indices of a.
  !       -           size = n.
  ! fia   - ar    - at the kth step, ar contains the elements of the
  !       -           reordered row of a.
  !       -           size = n.
  !
        INTEGER Ic(*) , Ia(*) , Ja(*) , Jar(*) , P(*) , Flag
  !     real  a(*), ar(*)
        DOUBLE PRECISION A(*) , Ar(*)
  !
  !  ******  for each nonempty row  *******************************
        DO k = 1 , N
           jmin = Ia(k)
           jmax = Ia(k+1) - 1
           IF ( jmin.LE.jmax ) THEN
              P(N+1) = N + 1
  !  ******  insert each element in the list  *********************
              DO j = jmin , jmax
                 newj = Ic(Ja(j))
                 i = N + 1
   10            IF ( P(i).GE.newj ) THEN
                    IF ( P(i).EQ.newj ) GOTO 100
                    P(newj) = P(i)
                    P(i) = newj
                    Jar(newj) = Ja(j)
                    Ar(newj) = A(j)
                 ELSE
                    i = P(i)
                    GOTO 10
                 ENDIF
              ENDDO
  !  ******  replace old row in ja and a  *************************
              i = N + 1
              DO j = jmin , jmax
                 i = P(i)
                 Ja(j) = Jar(i)
                 A(j) = Ar(i)
              ENDDO
           ENDIF
        ENDDO
        Flag = 0
        RETURN
  !
  ! ** error.. duplicate entry in a
   100  Flag = N + k
        END
  !*==NSFC.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
        SUBROUTINE NSFC(N,R,Ic,Ia,Ja,Jlmax,Il,Jl,Ijl,Jumax,Iu,Ju,Iju,Q,   &
       &                Ira,Jra,Irac,Irl,Jrl,Iru,Jru,Flag)
        IMPLICIT NONE
  !*--NSFC3507
  !*** Start of declarations inserted by SPAG
        INTEGER i , i1 , iak , irai , irll , irul , j , jaiak , jairai ,  &
       &        Jlmax , jlmin , jlptr , jmax , jmin , jtmp , Jumax ,      &
       &        jumin , juptr , k , lasti
        INTEGER lastid , long , luk , m , N , np1
  !*** End of declarations inserted by SPAG
  !*** subroutine nsfc
  !*** symbolic ldu-factorization of nonsymmetric sparse matrix
  !      (compressed pointer storage)
  !
  !
  !       input variables.. n, r, ic, ia, ja, jlmax, jumax.
  !       output variables.. il, jl, ijl, iu, ju, iju, flag.
  !
  !       parameters used internally..
  ! nia   - q     - suppose  m*  is the result of reordering  m.  if
  !       -           processing of the ith row of  m*  (hence the ith
  !       -           row of  u) is being done,  q(j)  is initially
  !       -           nonzero if  m*(i,j) is nonzero (j.ge.i).  since
  !       -           values need not be stored, each entry points to the
  !       -           next nonzero and  q(n+1)  points to the first.  n+1
  !       -           indicates the end of the list.  for example, if n=9
  !       -           and the 5th row of  m*  is
  !       -              0 x x 0 x 0 0 x 0
  !       -           then  q  will initially be
  !       -              a a a a 8 a a 10 5           (a - arbitrary).
  !       -           as the algorithm proceeds, other elements of  q
  !       -           are inserted in the list because of fillin.
  !       -           q  is used in an analogous manner to compute the
  !       -           ith column of  l.
  !       -           size = n+1.
  ! nia   - ira,  - vectors used to find the columns of  m.  at the kth
  ! nia   - jra,      step of the factorization,  irac(k)  points to the
  ! nia   - irac      head of a linked list in  jra  of row indices i
  !       -           such that i .ge. k and  m(i,k)  is nonzero.  zero
  !       -           indicates the end of the list.  ira(i)  (i.ge.k)
  !       -           points to the smallest j such that j .ge. k and
  !       -           m(i,j)  is nonzero.
  !       -           size of each = n.
  ! nia   - irl,  - vectors used to find the rows of  l.  at the kth step
  ! nia   - jrl       of the factorization,  jrl(k)  points to the head
  !       -           of a linked list in  jrl  of column indices j
  !       -           such j .lt. k and  l(k,j)  is nonzero.  zero
  !       -           indicates the end of the list.  irl(j)  (j.lt.k)
  !       -           points to the smallest i such that i .ge. k and
  !       -           l(i,j)  is nonzero.
  !       -           size of each = n.
  ! nia   - iru,  - vectors used in a manner analogous to  irl and jrl
  ! nia   - jru       to find the columns of  u.
  !       -           size of each = n.
  !
  !  internal variables..
  !    jlptr - points to the last position used in  jl.
  !    juptr - points to the last position used in  ju.
  !    jmin,jmax - are the indices in  a or u  of the first and last
  !                elements to be examined in a given row.
  !                for example,  jmin=ia(k), jmax=ia(k+1)-1.
  !
        INTEGER cend , qm , rend , rk , vj
        INTEGER Ia(*) , Ja(*) , Ira(*) , Jra(*) , Il(*) , Jl(*) , Ijl(*)
        INTEGER Iu(*) , Ju(*) , Iju(*) , Irl(*) , Jrl(*) , Iru(*) , Jru(*)
        INTEGER R(*) , Ic(*) , Q(*) , Irac(*) , Flag
  !
  !  ******  initialize pointers  ****************************************
        np1 = N + 1
        jlmin = 1
        jlptr = 0
        Il(1) = 1
        jumin = 1
        juptr = 0
        Iu(1) = 1
        DO k = 1 , N
           Irac(k) = 0
           Jra(k) = 0
           Jrl(k) = 0
           Jru(k) = 0
        ENDDO
  !  ******  initialize column pointers for a  ***************************
        DO k = 1 , N
           rk = R(k)
           iak = Ia(rk)
           IF ( iak.GE.Ia(rk+1) ) GOTO 400
           jaiak = Ic(Ja(iak))
           IF ( jaiak.GT.k ) GOTO 700
           Jra(k) = Irac(jaiak)
           Irac(jaiak) = k
           Ira(k) = iak
        ENDDO
  !
  !  ******  for each column of l and row of u  **************************
        DO k = 1 , N
  !
  !  ******  initialize q for computing kth column of l  *****************
           Q(np1) = np1
           luk = -1
  !  ******  by filling in kth column of a  ******************************
           vj = Irac(k)
           IF ( vj.NE.0 ) THEN
              qm = np1
   20         m = qm
              qm = Q(m)
              IF ( qm.LT.vj ) GOTO 20
              IF ( qm.EQ.vj ) GOTO 500
              luk = luk + 1
              Q(m) = vj
              Q(vj) = qm
              vj = Jra(vj)
              IF ( vj.NE.0 ) THEN
                 qm = np1
                 GOTO 20
              ENDIF
           ENDIF
  !  ******  link through jru  *******************************************
           lastid = 0
           lasti = 0
           Ijl(k) = jlptr
           i = k
   50      i = Jru(i)
           IF ( i.EQ.0 ) THEN
  !  ******  lasti is the longest column merged into the kth  ************
  !  ******  see if it equals the entire kth column  *********************
              qm = Q(np1)
              IF ( qm.NE.k ) GOTO 700
              IF ( luk.EQ.0 ) GOTO 100
              IF ( lastid.NE.luk ) THEN
  !  ******  if not, see if kth column can overlap the previous one  *****
                 IF ( jlmin.LE.jlptr ) THEN
                    qm = Q(qm)
                    DO j = jlmin , jlptr
                       IF ( Jl(j).LT.qm ) THEN
                       ELSEIF ( Jl(j).EQ.qm ) THEN
                          GOTO 60
                       ELSE
                          GOTO 80
                       ENDIF
                    ENDDO
                 ENDIF
                 GOTO 80
              ELSE
  !  ******  if so, jl can be compressed  ********************************
                 irll = Irl(lasti)
                 Ijl(k) = irll + 1
                 IF ( Jl(irll).NE.k ) Ijl(k) = Ijl(k) - 1
                 GOTO 100
              ENDIF
   60         Ijl(k) = j
              DO i = j , jlptr
                 IF ( Jl(i).NE.qm ) GOTO 80
                 qm = Q(qm)
                 IF ( qm.GT.N ) GOTO 100
              ENDDO
              jlptr = j - 1
  !  ******  move column indices from q to jl, update vectors  ***********
   80         jlmin = jlptr + 1
              Ijl(k) = jlmin
              IF ( luk.NE.0 ) THEN
                 jlptr = jlptr + luk
                 IF ( jlptr.GT.Jlmax ) GOTO 600
                 qm = Q(np1)
                 DO j = jlmin , jlptr
                    qm = Q(qm)
                    Jl(j) = qm
                 ENDDO
              ENDIF
   100        Irl(k) = Ijl(k)
              Il(k+1) = Il(k) + luk
  !
  !  ******  initialize q for computing kth row of u  ********************
              Q(np1) = np1
              luk = -1
  !  ******  by filling in kth row of reordered a  ***********************
              rk = R(k)
              jmin = Ira(k)
              jmax = Ia(rk+1) - 1
              IF ( jmin.LE.jmax ) THEN
                 DO j = jmin , jmax
                    vj = Ic(Ja(j))
                    qm = np1
   105              m = qm
                    qm = Q(m)
                    IF ( qm.LT.vj ) GOTO 105
                    IF ( qm.EQ.vj ) GOTO 500
                    luk = luk + 1
                    Q(m) = vj
                    Q(vj) = qm
                 ENDDO
              ENDIF
  !  ******  link through jrl,  ******************************************
              lastid = 0
              lasti = 0
              Iju(k) = juptr
              i = k
              i1 = Jrl(k)
   120        i = i1
              IF ( i.EQ.0 ) THEN
  !  ******  update jrl(k) and irl(k)  ***********************************
                 IF ( Il(k+1).GT.Il(k) ) THEN
                    j = Jl(Irl(k))
                    Jrl(k) = Jrl(j)
                    Jrl(j) = k
                 ENDIF
  !  ******  lasti is the longest row merged into the kth  ***************
  !  ******  see if it equals the entire kth row  ************************
                 qm = Q(np1)
                 IF ( qm.NE.k ) GOTO 700
                 IF ( luk.EQ.0 ) GOTO 150
                 IF ( lastid.NE.luk ) THEN
  !  ******  if not, see if kth row can overlap the previous one  ********
                    IF ( jumin.LE.juptr ) THEN
                       qm = Q(qm)
                       DO j = jumin , juptr
                          IF ( Ju(j).LT.qm ) THEN
                          ELSEIF ( Ju(j).EQ.qm ) THEN
                             GOTO 130
                          ELSE
                             GOTO 140
                          ENDIF
                       ENDDO
                    ENDIF
                    GOTO 140
                 ELSE
  !  ******  if so, ju can be compressed  ********************************
                    irul = Iru(lasti)
                    Iju(k) = irul + 1
                    IF ( Ju(irul).NE.k ) Iju(k) = Iju(k) - 1
                    GOTO 150
                 ENDIF
   130           Iju(k) = j
                 DO i = j , juptr
                    IF ( Ju(i).NE.qm ) GOTO 140
                    qm = Q(qm)
                    IF ( qm.GT.N ) GOTO 150
                 ENDDO
                 juptr = j - 1
  !  ******  move row indices from q to ju, update vectors  **************
   140           jumin = juptr + 1
                 Iju(k) = jumin
                 IF ( luk.NE.0 ) THEN
                    juptr = juptr + luk
                    IF ( juptr.GT.Jumax ) GOTO 800
                    qm = Q(np1)
                    DO j = jumin , juptr
                       qm = Q(qm)
                       Ju(j) = qm
                    ENDDO
                 ENDIF
   150           Iru(k) = Iju(k)
                 Iu(k+1) = Iu(k) + luk
  !
  !  ******  update iru, jru  ********************************************
                 i = k
              ELSE
                 i1 = Jrl(i)
                 qm = np1
                 jmin = Iru(i)
                 jmax = Iju(i) + Iu(i+1) - Iu(i) - 1
                 long = jmax - jmin
                 IF ( long.GE.0 ) THEN
                    jtmp = Ju(jmin)
                    IF ( jtmp.NE.k ) THEN
  !  ******  update irl and jrl, *****************************************
                       long = long + 1
                       cend = Ijl(i) + Il(i+1) - Il(i)
                       Irl(i) = Irl(i) + 1
                       IF ( Irl(i).LT.cend ) THEN
                          j = Jl(Irl(i))
                          Jrl(i) = Jrl(j)
                          Jrl(j) = i
                       ENDIF
                    ENDIF
                    IF ( lastid.LT.long ) THEN
                       lasti = i
                       lastid = long
                    ENDIF
  !  ******  and merge the corresponding rows into the kth row  **********
                    DO j = jmin , jmax
                       vj = Ju(j)
   152                 m = qm
                       qm = Q(m)
                       IF ( qm.LT.vj ) GOTO 152
                       IF ( qm.NE.vj ) THEN
                          luk = luk + 1
                          Q(m) = vj
                          Q(vj) = qm
                          qm = vj
                       ENDIF
                    ENDDO
                 ENDIF
                 GOTO 120
              ENDIF
           ELSE
              qm = np1
              jmin = Irl(i)
              jmax = Ijl(i) + Il(i+1) - Il(i) - 1
              long = jmax - jmin
              IF ( long.GE.0 ) THEN
                 jtmp = Jl(jmin)
                 IF ( jtmp.NE.k ) long = long + 1
                 IF ( jtmp.EQ.k ) R(i) = -R(i)
                 IF ( lastid.LT.long ) THEN
                    lasti = i
                    lastid = long
                 ENDIF
  !  ******  and merge the corresponding columns into the kth column  ****
                 DO j = jmin , jmax
                    vj = Jl(j)
   155              m = qm
                    qm = Q(m)
                    IF ( qm.LT.vj ) GOTO 155
                    IF ( qm.NE.vj ) THEN
                       luk = luk + 1
                       Q(m) = vj
                       Q(vj) = qm
                       qm = vj
                    ENDIF
                 ENDDO
              ENDIF
              GOTO 50
           ENDIF
   200     i1 = Jru(i)
           IF ( R(i).LT.0 ) THEN
              R(i) = -R(i)
           ELSE
              rend = Iju(i) + Iu(i+1) - Iu(i)
              IF ( Iru(i).LT.rend ) THEN
                 j = Ju(Iru(i))
                 Jru(i) = Jru(j)
                 Jru(j) = i
              ENDIF
           ENDIF
           i = i1
           IF ( i.EQ.0 ) THEN
  !
  !  ******  update ira, jra, irac  **************************************
              i = Irac(k)
              IF ( i.EQ.0 ) GOTO 300
           ELSE
              Iru(i) = Iru(i) + 1
              GOTO 200
           ENDIF
   250     i1 = Jra(i)
           Ira(i) = Ira(i) + 1
           IF ( Ira(i).LT.Ia(R(i)+1) ) THEN
              irai = Ira(i)
              jairai = Ic(Ja(irai))
              IF ( jairai.LE.i ) THEN
                 Jra(i) = Irac(jairai)
                 Irac(jairai) = i
              ENDIF
           ENDIF
           i = i1
           IF ( i.NE.0 ) GOTO 250
   300  ENDDO
  !
        Ijl(N) = jlptr
        Iju(N) = juptr
        Flag = 0
        RETURN
  !
  ! ** error.. null row in a
   400  Flag = N + rk
        RETURN
  ! ** error.. duplicate entry in a
   500  Flag = 2*N + rk
        RETURN
  ! ** error.. insufficient storage for jl
   600  Flag = 3*N + k
        RETURN
  ! ** error.. null pivot
   700  Flag = 5*N + k
        RETURN
  ! ** error.. insufficient storage for ju
   800  Flag = 6*N + k
        END
  !*==NNFC.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
        SUBROUTINE NNFC(N,R,C,Ic,Ia,Ja,A,Z,B,Lmax,Il,Jl,Ijl,L,D,Umax,Iu,  &
       &                Ju,Iju,U,Row,Tmp,Irl,Jrl,Flag)
        IMPLICIT NONE
  !*--NNFC3886
  !*** Start of declarations inserted by SPAG
        INTEGER i , i1 , i2 , ijlb , j , jmax , jmin , k , Lmax , mu , N
  !*** End of declarations inserted by SPAG
  !*** subroutine nnfc
  !*** numerical ldu-factorization of sparse nonsymmetric matrix and
  !      solution of system of linear equations (compressed pointer
  !      storage)
  !
  !
  !       input variables..  n, r, c, ic, ia, ja, a, b,
  !                          il, jl, ijl, lmax, iu, ju, iju, umax
  !       output variables.. z, l, d, u, flag
  !
  !       parameters used internally..
  ! nia   - irl,  - vectors used to find the rows of  l.  at the kth step
  ! nia   - jrl       of the factorization,  jrl(k)  points to the head
  !       -           of a linked list in  jrl  of column indices j
  !       -           such j .lt. k and  l(k,j)  is nonzero.  zero
  !       -           indicates the end of the list.  irl(j)  (j.lt.k)
  !       -           points to the smallest i such that i .ge. k and
  !       -           l(i,j)  is nonzero.
  !       -           size of each = n.
  ! fia   - row   - holds intermediate values in calculation of  u and l.
  !       -           size = n.
  ! fia   - tmp   - holds new right-hand side  b*  for solution of the
  !       -           equation ux = b*.
  !       -           size = n.
  !
  !  internal variables..
  !    jmin, jmax - indices of the first and last positions in a row to
  !      be examined.
  !    sum - used in calculating  tmp.
  !
        INTEGER rk , Umax
        INTEGER R(*) , C(*) , Ic(*) , Ia(*) , Ja(*) , Il(*) , Jl(*) ,     &
       &        Ijl(*)
        INTEGER Iu(*) , Ju(*) , Iju(*) , Irl(*) , Jrl(*) , Flag
  !     real  a(*), l(*), d(*), u(*), z(*), b(*), row(*)
  !     real tmp(*), lki, sum, dk
        DOUBLE PRECISION A(*) , L(*) , D(*) , U(*) , Z(*) , B(*) , Row(*)
        DOUBLE PRECISION Tmp(*) , lki , sum , dk
  !
  !  ******  initialize pointers and test storage  ***********************
        IF ( Il(N+1)-1.GT.Lmax ) THEN
  !
  ! ** error.. insufficient storage for l
           Flag = 4*N + 1
           RETURN
        ELSE
           IF ( Iu(N+1)-1.GT.Umax ) THEN
  ! ** error.. insufficient storage for u
              Flag = 7*N + 1
              RETURN
           ELSE
              DO k = 1 , N
                 Irl(k) = Il(k)
                 Jrl(k) = 0
              ENDDO
  !
  !  ******  for each row  ***********************************************
              DO k = 1 , N
  !  ******  reverse jrl and zero row where kth row of l will fill in  ***
                 Row(k) = 0
                 i1 = 0
                 IF ( Jrl(k).NE.0 ) THEN
                    i = Jrl(k)
   5                i2 = Jrl(i)
                    Jrl(i) = i1
                    i1 = i
                    Row(i) = 0
                    i = i2
                    IF ( i.NE.0 ) GOTO 5
                 ENDIF
  !  ******  set row to zero where u will fill in  ***********************
                 jmin = Iju(k)
                 jmax = jmin + Iu(k+1) - Iu(k) - 1
                 IF ( jmin.LE.jmax ) THEN
                    DO j = jmin , jmax
                       Row(Ju(j)) = 0
                    ENDDO
                 ENDIF
  !  ******  place kth row of a in row  **********************************
                 rk = R(k)
                 jmin = Ia(rk)
                 jmax = Ia(rk+1) - 1
                 DO j = jmin , jmax
                    Row(Ic(Ja(j))) = A(j)
                 ENDDO
  !  ******  initialize sum, and link through jrl  ***********************
                 sum = B(rk)
                 i = i1
                 IF ( i.EQ.0 ) GOTO 20
  !  ******  assign the kth row of l and adjust row, sum  ****************
   10            lki = -Row(i)
  !  ******  if l is not required, then comment out the following line  **
                 L(Irl(i)) = -lki
                 sum = sum + lki*Tmp(i)
                 jmin = Iu(i)
                 jmax = Iu(i+1) - 1
                 IF ( jmin.LE.jmax ) THEN
                    mu = Iju(i) - jmin
                    DO j = jmin , jmax
                       Row(Ju(mu+j)) = Row(Ju(mu+j)) + lki*U(j)
                    ENDDO
                 ENDIF
                 i = Jrl(i)
                 IF ( i.NE.0 ) GOTO 10
  !
  !  ******  assign kth row of u and diagonal d, set tmp(k)  *************
   20            IF ( Row(k).EQ.0.0D0 ) GOTO 100
                 dk = 1.0D0/Row(k)
                 D(k) = dk
                 Tmp(k) = sum*dk
                 IF ( k.EQ.N ) GOTO 60
                 jmin = Iu(k)
                 jmax = Iu(k+1) - 1
                 IF ( jmin.LE.jmax ) THEN
                    mu = Iju(k) - jmin
                    DO j = jmin , jmax
                       U(j) = Row(Ju(mu+j))*dk
                    ENDDO
                 ENDIF
  !
  !  ******  update irl and jrl, keeping jrl in decreasing order  ********
                 i = i1
                 IF ( i.EQ.0 ) GOTO 40
   30            Irl(i) = Irl(i) + 1
                 i1 = Jrl(i)
                 IF ( Irl(i).LT.Il(i+1) ) THEN
                    ijlb = Irl(i) - Il(i) + Ijl(i)
                    j = Jl(ijlb)
   35               IF ( i.GT.Jrl(j) ) THEN
                       Jrl(i) = Jrl(j)
                       Jrl(j) = i
                    ELSE
                       j = Jrl(j)
                       GOTO 35
                    ENDIF
                 ENDIF
                 i = i1
                 IF ( i.NE.0 ) GOTO 30
   40            IF ( Irl(k).LT.Il(k+1) ) THEN
                    j = Jl(Ijl(k))
                    Jrl(k) = Jrl(j)
                    Jrl(j) = k
                 ENDIF
   60         ENDDO
  !
  !  ******  solve  ux = tmp  by back substitution  **********************
              k = N
              DO i = 1 , N
                 sum = Tmp(k)
                 jmin = Iu(k)
                 jmax = Iu(k+1) - 1
                 IF ( jmin.LE.jmax ) THEN
                    mu = Iju(k) - jmin
                    DO j = jmin , jmax
                       sum = sum - U(j)*Tmp(Ju(mu+j))
                    ENDDO
                 ENDIF
                 Tmp(k) = sum
                 Z(C(k)) = sum
                 k = k - 1
              ENDDO
              Flag = 0
              RETURN
           ENDIF
  ! ** error.. zero pivot
   100     Flag = 8*N + k
        ENDIF
        END
  !*==NNSC.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
        SUBROUTINE NNSC(N,R,C,Il,Jl,Ijl,L,D,Iu,Ju,Iju,U,Z,B,Tmp)
        IMPLICIT NONE
  !*--NNSC4061
  !*** Start of declarations inserted by SPAG
        INTEGER i , j , jmax , jmin , k , ml , mu , N
  !*** End of declarations inserted by SPAG
  !*** subroutine nnsc
  !*** numerical solution of sparse nonsymmetric system of linear
  !      equations given ldu-factorization (compressed pointer storage)
  !
  !
  !       input variables..  n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, b
  !       output variables.. z
  !
  !       parameters used internally..
  ! fia   - tmp   - temporary vector which gets result of solving  ly = b.
  !       -           size = n.
  !
  !  internal variables..
  !    jmin, jmax - indices of the first and last positions in a row of
  !      u or l  to be used.
  !
        INTEGER R(*) , C(*) , Il(*) , Jl(*) , Ijl(*) , Iu(*) , Ju(*) ,    &
       &        Iju(*)
  !     real l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk, sum
        DOUBLE PRECISION L(*) , D(*) , U(*) , B(*) , Z(*) , Tmp(*) ,      &
       &                 tmpk , sum
  !
  !  ******  set tmp to reordered b  *************************************
        DO k = 1 , N
           Tmp(k) = B(R(k))
        ENDDO
  !  ******  solve  ly = b  by forward substitution  *********************
        DO k = 1 , N
           jmin = Il(k)
           jmax = Il(k+1) - 1
           tmpk = -D(k)*Tmp(k)
           Tmp(k) = -tmpk
           IF ( jmin.LE.jmax ) THEN
              ml = Ijl(k) - jmin
              DO j = jmin , jmax
                 Tmp(Jl(ml+j)) = Tmp(Jl(ml+j)) + tmpk*L(j)
              ENDDO
           ENDIF
        ENDDO
  !  ******  solve  ux = y  by back substitution  ************************
        k = N
        DO i = 1 , N
           sum = -Tmp(k)
           jmin = Iu(k)
           jmax = Iu(k+1) - 1
           IF ( jmin.LE.jmax ) THEN
              mu = Iju(k) - jmin
              DO j = jmin , jmax
                 sum = sum + U(j)*Tmp(Ju(mu+j))
              ENDDO
           ENDIF
           Tmp(k) = -sum
           Z(C(k)) = -sum
           k = k - 1
        ENDDO
        END
  !*==NNTC.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
        SUBROUTINE NNTC(N,R,C,Il,Jl,Ijl,L,D,Iu,Ju,Iju,U,Z,B,Tmp)
        IMPLICIT NONE
  !*--NNTC4124
  !*** Start of declarations inserted by SPAG
        INTEGER i , j , jmax , jmin , k , ml , mu , N
  !*** End of declarations inserted by SPAG
  !*** subroutine nntc
  !*** numeric solution of the transpose of a sparse nonsymmetric system
  !      of linear equations given lu-factorization (compressed pointer
  !      storage)
  !
  !
  !       input variables..  n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, b
  !       output variables.. z
  !
  !       parameters used internally..
  ! fia   - tmp   - temporary vector which gets result of solving ut y = b
  !       -           size = n.
  !
  !  internal variables..
  !    jmin, jmax - indices of the first and last positions in a row of
  !      u or l  to be used.
  !
        INTEGER R(*) , C(*) , Il(*) , Jl(*) , Ijl(*) , Iu(*) , Ju(*) ,    &
       &        Iju(*)
  !     real l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk,sum
        DOUBLE PRECISION L(*) , D(*) , U(*) , B(*) , Z(*) , Tmp(*) ,      &
       &                 tmpk , sum
  !
  !  ******  set tmp to reordered b  *************************************
        DO k = 1 , N
           Tmp(k) = B(C(k))
        ENDDO
  !  ******  solve  ut y = b  by forward substitution  *******************
        DO k = 1 , N
           jmin = Iu(k)
           jmax = Iu(k+1) - 1
           tmpk = -Tmp(k)
           IF ( jmin.LE.jmax ) THEN
              mu = Iju(k) - jmin
              DO j = jmin , jmax
                 Tmp(Ju(mu+j)) = Tmp(Ju(mu+j)) + tmpk*U(j)
              ENDDO
           ENDIF
        ENDDO
  !  ******  solve  lt x = y  by back substitution  **********************
        k = N
        DO i = 1 , N
           sum = -Tmp(k)
           jmin = Il(k)
           jmax = Il(k+1) - 1
           IF ( jmin.LE.jmax ) THEN
              ml = Ijl(k) - jmin
              DO j = jmin , jmax
                 sum = sum + L(j)*Tmp(Jl(ml+j))
              ENDDO
           ENDIF
           Tmp(k) = -sum*D(k)
           Z(R(k)) = Tmp(k)
           k = k - 1
        ENDDO
        END
  !*==DSTODA.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DSTODA
        SUBROUTINE DSTODA(Neq,Y,Yh,Nyh,Yh1,Ewt,Savf,Acor,Wm,Iwm,F,JAC,    &
       &                  PJAC,SLVS)
        IMPLICIT NONE
  !*--DSTODA4189
  !*** Start of declarations inserted by SPAG
        !INTEGER JAC
  !*** End of declarations inserted by SPAG
        EXTERNAL F , JAC , PJAC , SLVS
        INTEGER Neq , Nyh , Iwm
        DOUBLE PRECISION Y , Yh , Yh1 , Ewt , Savf , Acor , Wm
        DIMENSION Neq(*) , Y(*) , Yh(Nyh,*) , Yh1(*) , Ewt(*) , Savf(*) , &
       &          Acor(*) , Wm(*) , Iwm(*)
        INTEGER IOWnd , IALth , IPUp , LMAx , MEO , NQNyh , NSLp , ICF ,  &
       &        IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , LEWt ,  &
       &        LACor , LSAvf , LWM , LIWm , METh , MITer , MAXord ,      &
       &        MAXcor , MSBp , MXNcf , N , NQ , NST , NFE , NJE , NQU
        INTEGER IOWnd2 , ICOunt , IRFlag , JTYp , MUSed , MXOrdn , MXOrds
        DOUBLE PRECISION CONit , CRAte , EL , ELCo , HOLd , RMAx , TESco ,&
       &                 CCMax , EL0 , H , HMIn , HMXi , HU , RC , TN ,   &
       &                 UROund
        DOUBLE PRECISION ROWnd2 , CM1 , CM2 , PDEst , PDLast , RATio ,    &
       &                 PDNorm
        COMMON /DLS001/ CONit , CRAte , EL(13) , ELCo(13,12) , HOLd ,     &
       &                RMAx , TESco(3,12) , CCMax , EL0 , H , HMIn ,     &
       &                HMXi , HU , RC , TN , UROund , IOWnd(6) , IALth , &
       &                IPUp , LMAx , MEO , NQNyh , NSLp , ICF , IERpj ,  &
       &                IERsl , JCUr , JSTart , KFLag , L , LYH , LEWt ,  &
       &                LACor , LSAvf , LWM , LIWm , METh , MITer ,       &
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        COMMON /DLSA01/ ROWnd2 , CM1(12) , CM2(5) , PDEst , PDLast ,      &
       &                RATio , PDNorm , IOWnd2(3) , ICOunt , IRFlag ,    &
       &                JTYp , MUSed , MXOrdn , MXOrds
        INTEGER i , i1 , iredo , iret , j , jb , m , ncf , newq
        INTEGER lm1 , lm1p1 , lm2 , lm2p1 , nqm1 , nqm2
        DOUBLE PRECISION dcon , ddn , del , delp , dsm , dup , exdn ,     &
       &                 exsm , exup , r , rh , rhdn , rhsm , rhup , told
        DOUBLE PRECISION alpha , dm1 , dm2 , exm1 , exm2 , pdh , pnorm ,  &
       &                 rate , rh1 , rh1it , rh2 , rm , sm1(12)
        SAVE sm1
        DATA sm1/0.5D0 , 0.575D0 , 0.55D0 , 0.45D0 , 0.35D0 , 0.25D0 ,    &
       &     0.20D0 , 0.15D0 , 0.10D0 , 0.075D0 , 0.050D0 , 0.025D0/
  !-----------------------------------------------------------------------
  ! DSTODA performs one step of the integration of an initial value
  ! problem for a system of ordinary differential equations.
  ! Note: DSTODA is independent of the value of the iteration method
  ! indicator MITER, when this is .ne. 0, and hence is independent
  ! of the type of chord method used, or the Jacobian structure.
  ! Communication with DSTODA is done with the following variables:
  !
  ! Y      = an array of length .ge. N used as the Y argument in
  !          all calls to F and JAC.
  ! NEQ    = integer array containing problem size in NEQ(1), and
  !          passed as the NEQ argument in all calls to F and JAC.
  ! YH     = an NYH by LMAX array containing the dependent variables
  !          and their approximate scaled derivatives, where
  !          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
  !          j-th derivative of y(i), scaled by H**j/factorial(j)
  !          (j = 0,1,...,NQ).  On entry for the first step, the first
  !          two columns of YH must be set from the initial values.
  ! NYH    = a constant integer .ge. N, the first dimension of YH.
  ! YH1    = a one-dimensional array occupying the same space as YH.
  ! EWT    = an array of length N containing multiplicative weights
  !          for local error measurements.  Local errors in y(i) are
  !          compared to 1.0/EWT(i) in various error tests.
  ! SAVF   = an array of working storage, of length N.
  ! ACOR   = a work array of length N, used for the accumulated
  !          corrections.  On a successful return, ACOR(i) contains
  !          the estimated one-step local error in y(i).
  ! WM,IWM = real and integer work arrays associated with matrix
  !          operations in chord iteration (MITER .ne. 0).
  ! PJAC   = name of routine to evaluate and preprocess Jacobian matrix
  !          and P = I - H*EL0*Jac, if a chord method is being used.
  !          It also returns an estimate of norm(Jac) in PDNORM.
  ! SLVS   = name of routine to solve linear system in chord iteration.
  ! CCMAX  = maximum relative change in H*EL0 before PJAC is called.
  ! H      = the step size to be attempted on the next step.
  !          H is altered by the error control algorithm during the
  !          problem.  H can be either positive or negative, but its
  !          sign must remain constant throughout the problem.
  ! HMIN   = the minimum absolute value of the step size H to be used.
  ! HMXI   = inverse of the maximum absolute value of H to be used.
  !          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
  !          HMIN and HMXI may be changed at any time, but will not
  !          take effect until the next change of H is considered.
  ! TN     = the independent variable. TN is updated on each step taken.
  ! JSTART = an integer used for input only, with the following
  !          values and meanings:
  !               0  perform the first step.
  !           .gt.0  take a new step continuing from the last.
  !              -1  take the next step with a new value of H,
  !                    N, METH, MITER, and/or matrix parameters.
  !              -2  take the next step with a new value of H,
  !                    but with other inputs unchanged.
  !          On return, JSTART is set to 1 to facilitate continuation.
  ! KFLAG  = a completion code with the following meanings:
  !               0  the step was succesful.
  !              -1  the requested error could not be achieved.
  !              -2  corrector convergence could not be achieved.
  !              -3  fatal error in PJAC or SLVS.
  !          A return with KFLAG = -1 or -2 means either
  !          ABS(H) = HMIN or 10 consecutive failures occurred.
  !          On a return with KFLAG negative, the values of TN and
  !          the YH array are as of the beginning of the last
  !          step, and H is the last step size attempted.
  ! MAXORD = the maximum order of integration method to be allowed.
  ! MAXCOR = the maximum number of corrector iterations allowed.
  ! MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0).
  ! MXNCF  = maximum number of convergence failures allowed.
  ! METH   = current method.
  !          METH = 1 means Adams method (nonstiff)
  !          METH = 2 means BDF method (stiff)
  !          METH may be reset by DSTODA.
  ! MITER  = corrector iteration method.
  !          MITER = 0 means functional iteration.
  !          MITER = JT .gt. 0 means a chord iteration corresponding
  !          to Jacobian type JT.  (The DLSODA/DLSODAR argument JT is
  !          communicated here as JTYP, but is not used in DSTODA
  !          except to load MITER following a method switch.)
  !          MITER may be reset by DSTODA.
  ! N      = the number of first-order differential equations.
  !-----------------------------------------------------------------------
        KFLag = 0
        told = TN
        ncf = 0
        IERpj = 0
        IERsl = 0
        JCUr = 0
        ICF = 0
        delp = 0.0D0
        IF ( JSTart.GT.0 ) GOTO 500
        IF ( JSTart.EQ.-1 ) THEN
  !-----------------------------------------------------------------------
  ! The following block handles preliminaries needed when JSTART = -1.
  ! IPUP is set to MITER to force a matrix update.
  ! If an order increase is about to be considered (IALTH = 1),
  ! IALTH is reset to 2 to postpone consideration one more step.
  ! If the caller has changed METH, DCFODE is called to reset
  ! the coefficients of the method.
  ! If H is to be changed, YH must be rescaled.
  ! If H or METH is being changed, IALTH is reset to L = NQ + 1
  ! to prevent further changes in H for that many steps.
  !-----------------------------------------------------------------------
           IPUp = MITer
           LMAx = MAXord + 1
           IF ( IALth.EQ.1 ) IALth = 2
           IF ( METh.EQ.MUSed ) GOTO 200
           CALL DCFODE(METh,ELCo,TESco)
           IALth = L
           iret = 1
        ELSE
           IF ( JSTart.EQ.-2 ) GOTO 200
  !-----------------------------------------------------------------------
  ! On the first call, the order is set to 1, and other variables are
  ! initialized.  RMAX is the maximum ratio by which H can be increased
  ! in a single step.  It is initially 1.E4 to compensate for the small
  ! initial H, but then is normally equal to 10.  If a failure
  ! occurs (in corrector convergence or error test), RMAX is set at 2
  ! for the next increase.
  ! DCFODE is called to get the needed coefficients for both methods.
  !-----------------------------------------------------------------------
           LMAx = MAXord + 1
           NQ = 1
           L = 2
           IALth = 2
           RMAx = 10000.0D0
           RC = 0.0D0
           EL0 = 1.0D0
           CRAte = 0.7D0
           HOLd = H
           NSLp = 0
           IPUp = MITer
           iret = 3
  ! Initialize switching parameters.  METH = 1 is assumed initially. -----
           ICOunt = 20
           IRFlag = 0
           PDEst = 0.0D0
           PDLast = 0.0D0
           RATio = 5.0D0
           CALL DCFODE(2,ELCo,TESco)
           DO i = 1 , 5
              CM2(i) = TESco(2,i)*ELCo(i+1,i)
           ENDDO
           CALL DCFODE(1,ELCo,TESco)
           DO i = 1 , 12
              CM1(i) = TESco(2,i)*ELCo(i+1,i)
           ENDDO
        ENDIF
  !-----------------------------------------------------------------------
  ! The el vector and related constants are reset
  ! whenever the order NQ is changed, or at the start of the problem.
  !-----------------------------------------------------------------------
   100  DO i = 1 , L
           EL(i) = ELCo(i,NQ)
        ENDDO
        NQNyh = NQ*Nyh
        RC = RC*EL(1)/EL0
        EL0 = EL(1)
        CONit = 0.5D0/(NQ+2)
        GOTO (200,300,500) , iret
  !-----------------------------------------------------------------------
  ! If H is being changed, the H ratio RH is checked against
  ! RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
  ! L = NQ + 1 to prevent a change of H for that many steps, unless
  ! forced by a convergence or error test failure.
  !-----------------------------------------------------------------------
   200  IF ( H.EQ.HOLd ) GOTO 500
        rh = H/HOLd
        H = HOLd
        iredo = 3
        GOTO 400
   300  rh = MAX(rh,HMIn/ABS(H))
   400  rh = MIN(rh,RMAx)
        rh = rh/MAX(1.0D0,ABS(H)*HMXi*rh)
  !-----------------------------------------------------------------------
  ! If METH = 1, also restrict the new step size by the stability region.
  ! If this reduces H, set IRFLAG to 1 so that if there are roundoff
  ! problems later, we can assume that is the cause of the trouble.
  !-----------------------------------------------------------------------
        IF ( METh.NE.2 ) THEN
           IRFlag = 0
           pdh = MAX(ABS(H)*PDLast,0.000001D0)
           IF ( rh*pdh*1.00001D0.GE.sm1(NQ) ) THEN
              rh = sm1(NQ)/pdh
              IRFlag = 1
           ENDIF
        ENDIF
        r = 1.0D0
        DO j = 2 , L
           r = r*rh
           DO i = 1 , N
              Yh(i,j) = Yh(i,j)*r
           ENDDO
        ENDDO
        H = H*rh
        RC = RC*rh
        IALth = L
        IF ( iredo.EQ.0 ) THEN
           RMAx = 10.0D0
           GOTO 1400
        ENDIF
  !-----------------------------------------------------------------------
  ! This section computes the predicted values by effectively
  ! multiplying the YH array by the Pascal triangle matrix.
  ! RC is the ratio of new to old values of the coefficient  H*EL(1).
  ! When RC differs from 1 by more than CCMAX, IPUP is set to MITER
  ! to force PJAC to be called, if a Jacobian is involved.
  ! In any case, PJAC is called at least every MSBP steps.
  !-----------------------------------------------------------------------
   500  IF ( ABS(RC-1.0D0).GT.CCMax ) IPUp = MITer
        IF ( NST.GE.NSLp+MSBp ) IPUp = MITer
        TN = TN + H
        i1 = NQNyh + 1
        DO jb = 1 , NQ
           i1 = i1 - Nyh
  !DIR$ IVDEP
           DO i = i1 , NQNyh
              Yh1(i) = Yh1(i) + Yh1(i+Nyh)
           ENDDO
        ENDDO
        pnorm = DMNORM(N,Yh1,Ewt)
  !-----------------------------------------------------------------------
  ! Up to MAXCOR corrector iterations are taken.  A convergence test is
  ! made on the RMS-norm of each correction, weighted by the error
  ! weight vector EWT.  The sum of the corrections is accumulated in the
  ! vector ACOR(i).  The YH array is not altered in the corrector loop.
  !-----------------------------------------------------------------------
   600  m = 0
        rate = 0.0D0
        del = 0.0D0
        DO i = 1 , N
           Y(i) = Yh(i,1)
        ENDDO
        CALL F(Neq,TN,Y,Savf)
        NFE = NFE + 1
        IF ( IPUp.GT.0 ) THEN
  !-----------------------------------------------------------------------
  ! If indicated, the matrix P = I - H*EL(1)*J is reevaluated and
  ! preprocessed before starting the corrector iteration.  IPUP is set
  ! to 0 as an indicator that this has been done.
  !-----------------------------------------------------------------------
           CALL PJAC(Neq,Y,Yh,Nyh,Ewt,Acor,Savf,Wm,Iwm,F,JAC)
           IPUp = 0
           RC = 1.0D0
           NSLp = NST
           CRAte = 0.7D0
           IF ( IERpj.NE.0 ) GOTO 900
        ENDIF
        DO i = 1 , N
           Acor(i) = 0.0D0
        ENDDO
   700  IF ( MITer.NE.0 ) THEN
  !-----------------------------------------------------------------------
  ! In the case of the chord method, compute the corrector error,
  ! and solve the linear system with that as right-hand side and
  ! P as coefficient matrix.
  !-----------------------------------------------------------------------
           DO i = 1 , N
              Y(i) = H*Savf(i) - (Yh(i,2)+Acor(i))
           ENDDO
           CALL SLVS(Wm,Iwm,Y,Savf)
           IF ( IERsl.LT.0 ) GOTO 900
           IF ( IERsl.GT.0 ) GOTO 800
           del = DMNORM(N,Y,Ewt)
           DO i = 1 , N
              Acor(i) = Acor(i) + Y(i)
              Y(i) = Yh(i,1) + EL(1)*Acor(i)
           ENDDO
        ELSE
  !-----------------------------------------------------------------------
  ! In the case of functional iteration, update Y directly from
  ! the result of the last function evaluation.
  !-----------------------------------------------------------------------
           DO i = 1 , N
              Savf(i) = H*Savf(i) - Yh(i,2)
              Y(i) = Savf(i) - Acor(i)
           ENDDO
           del = DMNORM(N,Y,Ewt)
           DO i = 1 , N
              Y(i) = Yh(i,1) + EL(1)*Savf(i)
              Acor(i) = Savf(i)
           ENDDO
        ENDIF
  !-----------------------------------------------------------------------
  ! Test for convergence.  If M .gt. 0, an estimate of the convergence
  ! rate constant is stored in CRATE, and this is used in the test.
  !
  ! We first check for a change of iterates that is the size of
  ! roundoff error.  If this occurs, the iteration has converged, and a
  ! new rate estimate is not formed.
  ! In all other cases, force at least two iterations to estimate a
  ! local Lipschitz constant estimate for Adams methods.
  ! On convergence, form PDEST = local maximum Lipschitz constant
  ! estimate.  PDLAST is the most recent nonzero estimate.
  !-----------------------------------------------------------------------
        IF ( del.LE.100.0D0*pnorm*UROund ) GOTO 1000
        IF ( m.NE.0 .OR. METh.NE.1 ) THEN
           IF ( m.NE.0 ) THEN
              rm = 1024.0D0
              IF ( del.LE.1024.0D0*delp ) rm = del/delp
              rate = MAX(rate,rm)
              CRAte = MAX(0.2D0*CRAte,rm)
           ENDIF
           dcon = del*MIN(1.0D0,1.5D0*CRAte)/(TESco(2,NQ)*CONit)
           IF ( dcon.LE.1.0D0 ) THEN
              PDEst = MAX(PDEst,rate/ABS(H*EL(1)))
              IF ( PDEst.NE.0.0D0 ) PDLast = PDEst
              GOTO 1000
           ENDIF
        ENDIF
        m = m + 1
        IF ( m.NE.MAXcor ) THEN
           IF ( m.LT.2 .OR. del.LE.2.0D0*delp ) THEN
              delp = del
              CALL F(Neq,TN,Y,Savf)
              NFE = NFE + 1
              GOTO 700
           ENDIF
        ENDIF
  !-----------------------------------------------------------------------
  ! The corrector iteration failed to converge.
  ! If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
  ! the next try.  Otherwise the YH array is retracted to its values
  ! before prediction, and H is reduced, if possible.  If H cannot be
  ! reduced or MXNCF failures have occurred, exit with KFLAG = -2.
  !-----------------------------------------------------------------------
   800  IF ( MITer.NE.0 .AND. JCUr.NE.1 ) THEN
           ICF = 1
           IPUp = MITer
           GOTO 600
        ENDIF
   900  ICF = 2
        ncf = ncf + 1
        RMAx = 2.0D0
        TN = told
        i1 = NQNyh + 1
        DO jb = 1 , NQ
           i1 = i1 - Nyh
  !DIR$ IVDEP
           DO i = i1 , NQNyh
              Yh1(i) = Yh1(i) - Yh1(i+Nyh)
           ENDDO
        ENDDO
        IF ( IERpj.LT.0 .OR. IERsl.LT.0 ) THEN
           KFLag = -3
           GOTO 1500
        ELSEIF ( ABS(H).LE.HMIn*1.00001D0 ) THEN
           KFLag = -2
           GOTO 1500
        ELSEIF ( ncf.EQ.MXNcf ) THEN
           KFLag = -2
           GOTO 1500
        ELSE
           rh = 0.25D0
           IPUp = MITer
           iredo = 1
           GOTO 300
        ENDIF
  !-----------------------------------------------------------------------
  ! The corrector has converged.  JCUR is set to 0
  ! to signal that the Jacobian involved may need updating later.
  ! The local error test is made and control passes to statement 500
  ! if it fails.
  !-----------------------------------------------------------------------
   1000 JCUr = 0
        IF ( m.EQ.0 ) dsm = del/TESco(2,NQ)
        IF ( m.GT.0 ) dsm = DMNORM(N,Acor,Ewt)/TESco(2,NQ)
        IF ( dsm.GT.1.0D0 ) THEN
  !-----------------------------------------------------------------------
  ! The error test failed.  KFLAG keeps track of multiple failures.
  ! Restore TN and the YH array to their previous values, and prepare
  ! to try the step again.  Compute the optimum step size for this or
  ! one lower order.  After 2 or more failures, H is forced to decrease
  ! by a factor of 0.2 or less.
  !-----------------------------------------------------------------------
           KFLag = KFLag - 1
           TN = told
           i1 = NQNyh + 1
           DO jb = 1 , NQ
              i1 = i1 - Nyh
  !DIR$ IVDEP
              DO i = i1 , NQNyh
                 Yh1(i) = Yh1(i) - Yh1(i+Nyh)
              ENDDO
           ENDDO
           RMAx = 2.0D0
           IF ( ABS(H).LE.HMIn*1.00001D0 ) THEN
  !-----------------------------------------------------------------------
  ! All returns are made through this section.  H is saved in HOLD
  ! to allow the caller to change H on the next step.
  !-----------------------------------------------------------------------
              KFLag = -1
              GOTO 1500
           ELSEIF ( KFLag.LE.-3 ) THEN
  !-----------------------------------------------------------------------
  ! Control reaches this section if 3 or more failures have occured.
  ! If 10 failures have occurred, exit with KFLAG = -1.
  ! It is assumed that the derivatives that have accumulated in the
  ! YH array have errors of the wrong order.  Hence the first
  ! derivative is recomputed, and the order is set to 1.  Then
  ! H is reduced by a factor of 10, and the step is retried,
  ! until it succeeds or H reaches HMIN.
  !-----------------------------------------------------------------------
              IF ( KFLag.EQ.-10 ) THEN
                 KFLag = -1
                 GOTO 1500
              ELSE
                 rh = 0.1D0
                 rh = MAX(HMIn/ABS(H),rh)
                 H = H*rh
                 DO i = 1 , N
                    Y(i) = Yh(i,1)
                 ENDDO
                 CALL F(Neq,TN,Y,Savf)
                 NFE = NFE + 1
                 DO i = 1 , N
                    Yh(i,2) = H*Savf(i)
                 ENDDO
                 IPUp = MITer
                 IALth = 5
                 IF ( NQ.EQ.1 ) GOTO 500
                 NQ = 1
                 L = 2
                 iret = 3
                 GOTO 100
              ENDIF
           ELSE
              iredo = 2
              rhup = 0.0D0
           ENDIF
        ELSE
  !-----------------------------------------------------------------------
  ! After a successful step, update the YH array.
  ! Decrease ICOUNT by 1, and if it is -1, consider switching methods.
  ! If a method switch is made, reset various parameters,
  ! rescale the YH array, and exit.  If there is no switch,
  ! consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
  ! If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
  ! use in a possible order increase on the next step.
  ! If a change in H is considered, an increase or decrease in order
  ! by one is considered also.  A change in H is made only if it is by a
  ! factor of at least 1.1.  If not, IALTH is set to 3 to prevent
  ! testing for that many steps.
  !-----------------------------------------------------------------------
           KFLag = 0
           iredo = 0
           NST = NST + 1
           HU = H
           NQU = NQ
           MUSed = METh
           DO j = 1 , L
              DO i = 1 , N
                 Yh(i,j) = Yh(i,j) + EL(j)*Acor(i)
              ENDDO
           ENDDO
           ICOunt = ICOunt - 1
           IF ( ICOunt.LT.0 ) THEN
              IF ( METh.EQ.2 ) THEN
  !-----------------------------------------------------------------------
  ! We are currently using a BDF method.  Consider switching to Adams.
  ! Compute the step size we could have (ideally) used on this step,
  ! with the current (BDF) method, and also that for the Adams.
  ! If NQ .gt. MXORDN, we consider changing to order MXORDN on switching.
  ! Compare the two step sizes to decide whether to switch.
  ! The step size advantage must be at least 5/RATIO = 1 to switch.
  ! If the step size for Adams would be so small as to cause
  ! roundoff pollution, we stay with BDF.
  !-----------------------------------------------------------------------
                 exsm = 1.0D0/L
                 IF ( MXOrdn.GE.NQ ) THEN
                    dm1 = dsm*(CM2(NQ)/CM1(NQ))
                    rh1 = 1.0D0/(1.2D0*dm1**exsm+0.0000012D0)
                    nqm1 = NQ
                    exm1 = exsm
                 ELSE
                    nqm1 = MXOrdn
                    lm1 = MXOrdn + 1
                    exm1 = 1.0D0/lm1
                    lm1p1 = lm1 + 1
                    dm1 = DMNORM(N,Yh(1,lm1p1),Ewt)/CM1(MXOrdn)
                    rh1 = 1.0D0/(1.2D0*dm1**exm1+0.0000012D0)
                 ENDIF
                 rh1it = 2.0D0*rh1
                 pdh = PDNorm*ABS(H)
                 IF ( pdh*rh1.GT.0.00001D0 ) rh1it = sm1(nqm1)/pdh
                 rh1 = MIN(rh1,rh1it)
                 rh2 = 1.0D0/(1.2D0*dsm**exsm+0.0000012D0)
                 IF ( rh1*RATio.GE.5.0D0*rh2 ) THEN
                    alpha = MAX(0.001D0,rh1)
                    dm1 = (alpha**exm1)*dm1
                    IF ( dm1.GT.1000.0D0*UROund*pnorm ) THEN
  ! The switch test passed.  Reset relevant quantities for Adams. --------
                       rh = rh1
                       ICOunt = 20
                       METh = 1
                       MITer = 0
                       PDLast = 0.0D0
                       NQ = nqm1
                       L = NQ + 1
                       GOTO 300
                    ENDIF
                 ENDIF
  !-----------------------------------------------------------------------
  ! We are currently using an Adams method.  Consider switching to BDF.
  ! If the current order is greater than 5, assume the problem is
  ! not stiff, and skip this section.
  ! If the Lipschitz constant and error estimate are not polluted
  ! by roundoff, go to 470 and perform the usual test.
  ! Otherwise, switch to the BDF methods if the last step was
  ! restricted to insure stability (irflag = 1), and stay with Adams
  ! method if not.  When switching to BDF with polluted error estimates,
  ! in the absence of other information, double the step size.
  !
  ! When the estimates are OK, we make the usual test by computing
  ! the step size we could have (ideally) used on this step,
  ! with the current (Adams) method, and also that for the BDF.
  ! If NQ .gt. MXORDS, we consider changing to order MXORDS on switching.
  ! Compare the two step sizes to decide whether to switch.
  ! The step size advantage must be at least RATIO = 5 to switch.
  !-----------------------------------------------------------------------
              ELSEIF ( NQ.LE.5 ) THEN
                 IF ( dsm.GT.100.0D0*pnorm*UROund .AND. PDEst.NE.0.0D0 )  &
       &              THEN
                    exsm = 1.0D0/L
                    rh1 = 1.0D0/(1.2D0*dsm**exsm+0.0000012D0)
                    rh1it = 2.0D0*rh1
                    pdh = PDLast*ABS(H)
                    IF ( pdh*rh1.GT.0.00001D0 ) rh1it = sm1(NQ)/pdh
                    rh1 = MIN(rh1,rh1it)
                    IF ( NQ.LE.MXOrds ) THEN
                       dm2 = dsm*(CM1(NQ)/CM2(NQ))
                       rh2 = 1.0D0/(1.2D0*dm2**exsm+0.0000012D0)
                       nqm2 = NQ
                    ELSE
                       nqm2 = MXOrds
                       lm2 = MXOrds + 1
                       exm2 = 1.0D0/lm2
                       lm2p1 = lm2 + 1
                       dm2 = DMNORM(N,Yh(1,lm2p1),Ewt)/CM2(MXOrds)
                       rh2 = 1.0D0/(1.2D0*dm2**exm2+0.0000012D0)
                    ENDIF
                    IF ( rh2.LT.RATio*rh1 ) GOTO 1050
                 ELSE
                    IF ( IRFlag.EQ.0 ) GOTO 1050
                    rh2 = 2.0D0
                    nqm2 = MIN(NQ,MXOrds)
                 ENDIF
  ! THE SWITCH TEST PASSED.  RESET RELEVANT QUANTITIES FOR BDF. ----------
                 rh = rh2
                 ICOunt = 20
                 METh = 2
                 MITer = JTYp
                 PDLast = 0.0D0
                 NQ = nqm2
                 L = NQ + 1
                 GOTO 300
              ENDIF
           ENDIF
  !
  ! No method switch is being made.  Do the usual step/order selection. --
   1050    IALth = IALth - 1
           IF ( IALth.EQ.0 ) THEN
  !-----------------------------------------------------------------------
  ! Regardless of the success or failure of the step, factors
  ! RHDN, RHSM, and RHUP are computed, by which H could be multiplied
  ! at order NQ - 1, order NQ, or order NQ + 1, respectively.
  ! In the case of failure, RHUP = 0.0 to avoid an order increase.
  ! The largest of these is determined and the new order chosen
  ! accordingly.  If the order is to be increased, we compute one
  ! additional scaled derivative.
  !-----------------------------------------------------------------------
              rhup = 0.0D0
              IF ( L.NE.LMAx ) THEN
                 DO i = 1 , N
                    Savf(i) = Acor(i) - Yh(i,LMAx)
                 ENDDO
                 dup = DMNORM(N,Savf,Ewt)/TESco(3,NQ)
                 exup = 1.0D0/(L+1)
                 rhup = 1.0D0/(1.4D0*dup**exup+0.0000014D0)
              ENDIF
           ELSE
              IF ( IALth.LE.1 ) THEN
                 IF ( L.NE.LMAx ) THEN
                    DO i = 1 , N
                       Yh(i,LMAx) = Acor(i)
                    ENDDO
                 ENDIF
              ENDIF
              GOTO 1400
           ENDIF
        ENDIF
        exsm = 1.0D0/L
        rhsm = 1.0D0/(1.2D0*dsm**exsm+0.0000012D0)
        rhdn = 0.0D0
        IF ( NQ.NE.1 ) THEN
           ddn = DMNORM(N,Yh(1,L),Ewt)/TESco(1,NQ)
           exdn = 1.0D0/NQ
           rhdn = 1.0D0/(1.3D0*ddn**exdn+0.0000013D0)
        ENDIF
  ! If METH = 1, limit RH according to the stability region also. --------
        IF ( METh.NE.2 ) THEN
           pdh = MAX(ABS(H)*PDLast,0.000001D0)
           IF ( L.LT.LMAx ) rhup = MIN(rhup,sm1(L)/pdh)
           rhsm = MIN(rhsm,sm1(NQ)/pdh)
           IF ( NQ.GT.1 ) rhdn = MIN(rhdn,sm1(NQ-1)/pdh)
           PDEst = 0.0D0
        ENDIF
        IF ( rhsm.GE.rhup ) THEN
           IF ( rhsm.GE.rhdn ) THEN
              newq = NQ
              rh = rhsm
              GOTO 1100
           ENDIF
        ELSEIF ( rhup.GT.rhdn ) THEN
           newq = L
           rh = rhup
           IF ( rh.LT.1.1D0 ) THEN
              IALth = 3
              GOTO 1400
           ELSE
              r = EL(L)/L
              DO i = 1 , N
                 Yh(i,newq+1) = Acor(i)*r
              ENDDO
              GOTO 1300
           ENDIF
        ENDIF
        newq = NQ - 1
        rh = rhdn
        IF ( KFLag.LT.0 .AND. rh.GT.1.0D0 ) rh = 1.0D0
  ! If METH = 1 and H is restricted by stability, bypass 10 percent test.
   1100 IF ( METh.NE.2 ) THEN
           IF ( rh*pdh*1.00001D0.GE.sm1(newq) ) GOTO 1200
        ENDIF
        IF ( KFLag.EQ.0 .AND. rh.LT.1.1D0 ) THEN
           IALth = 3
           GOTO 1400
        ENDIF
   1200 IF ( KFLag.LE.-2 ) rh = MIN(rh,0.2D0)
  !-----------------------------------------------------------------------
  ! If there is a change of order, reset NQ, L, and the coefficients.
  ! In any case H is reset according to RH and the YH array is rescaled.
  ! Then exit from 690 if the step was OK, or redo the step otherwise.
  !-----------------------------------------------------------------------
        IF ( newq.EQ.NQ ) GOTO 300
   1300 NQ = newq
        L = NQ + 1
        iret = 2
        GOTO 100
   1400 r = 1.0D0/TESco(2,NQU)
        DO i = 1 , N
           Acor(i) = Acor(i)*r
        ENDDO
   1500 HOLd = H
        JSTart = 1
  !----------------------- End of Subroutine DSTODA ----------------------
        END
  !*==DPRJA.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DPRJA
        SUBROUTINE DPRJA(Neq,Y,Yh,Nyh,Ewt,Ftem,Savf,Wm,Iwm,F,JAC)
        IMPLICIT NONE
  !*--DPRJA4888
        EXTERNAL F , JAC
        INTEGER Neq , Nyh , Iwm
        DOUBLE PRECISION Y , Yh , Ewt , Ftem , Savf , Wm
        DIMENSION Neq(*) , Y(*) , Yh(Nyh,*) , Ewt(*) , Ftem(*) , Savf(*) ,&
       &          Wm(*) , Iwm(*)
        INTEGER IOWnd , IOWns , ICF , IERpj , IERsl , JCUr , JSTart ,     &
       &        KFLag , L , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,     &
       &        METh , MITer , MAXord , MAXcor , MSBp , MXNcf , N , NQ ,  &
       &        NST , NFE , NJE , NQU
        INTEGER IOWnd2 , IOWns2 , JTYp , MUSed , MXOrdn , MXOrds
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        DOUBLE PRECISION ROWnd2 , ROWns2 , PDNorm
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , IOWnd(6) , IOWns(6) , ICF ,    &
       &                IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , &
       &                LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,&
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        COMMON /DLSA01/ ROWnd2 , ROWns2(20) , PDNorm , IOWnd2(3) ,        &
       &                IOWns2(2) , JTYp , MUSed , MXOrdn , MXOrds
        INTEGER i , i1 , i2 , ier , ii , j , j1 , jj , lenp , mba ,       &
       &        mband , meb1 , meband , ml , ml3 , mu , np1
        DOUBLE PRECISION con , fac , hl0 , r , r0 , srur , yi , yj , yjj
  !-----------------------------------------------------------------------
  ! DPRJA is called by DSTODA to compute and process the matrix
  ! P = I - H*EL(1)*J , where J is an approximation to the Jacobian.
  ! Here J is computed by the user-supplied routine JAC if
  ! MITER = 1 or 4 or by finite differencing if MITER = 2 or 5.
  ! J, scaled by -H*EL(1), is stored in WM.  Then the norm of J (the
  ! matrix norm consistent with the weighted max-norm on vectors given
  ! by DMNORM) is computed, and J is overwritten by P.  P is then
  ! subjected to LU decomposition in preparation for later solution
  ! of linear systems with P as coefficient matrix.  This is done
  ! by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
  !
  ! In addition to variables described previously, communication
  ! with DPRJA uses the following:
  ! Y     = array containing predicted values on entry.
  ! FTEM  = work array of length N (ACOR in DSTODA).
  ! SAVF  = array containing f evaluated at predicted y.
  ! WM    = real work space for matrices.  On output it contains the
  !         LU decomposition of P.
  !         Storage of matrix elements starts at WM(3).
  !         WM also contains the following matrix-related data:
  !         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
  ! IWM   = integer work space containing pivot information, starting at
  !         IWM(21).   IWM also contains the band parameters
  !         ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
  ! EL0   = EL(1) (input).
  ! PDNORM= norm of Jacobian matrix. (Output).
  ! IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
  !         P matrix found to be singular.
  ! JCUR  = output flag = 1 to indicate that the Jacobian matrix
  !         (or approximation) is now current.
  ! This routine also uses the Common variables EL0, H, TN, UROUND,
  ! MITER, N, NFE, and NJE.
  !-----------------------------------------------------------------------
        NJE = NJE + 1
        IERpj = 0
        JCUr = 1
        hl0 = H*EL0
        GOTO (100,200,400,500,600) , MITer
  ! If MITER = 1, call JAC and multiply by scalar. -----------------------
   100  lenp = N*N
        DO i = 1 , lenp
           Wm(i+2) = 0.0D0
        ENDDO
        CALL JAC(Neq,TN,Y,0,0,Wm(3),N)
        con = -hl0
        DO i = 1 , lenp
           Wm(i+2) = Wm(i+2)*con
        ENDDO
        GOTO 300
  ! If MITER = 2, make N calls to F to approximate J. --------------------
   200  fac = DMNORM(N,Savf,Ewt)
        r0 = 1000.0D0*ABS(H)*UROund*N*fac
        IF ( r0.EQ.0.0D0 ) r0 = 1.0D0
        srur = Wm(1)
        j1 = 2
        DO j = 1 , N
           yj = Y(j)
           r = MAX(srur*ABS(yj),r0/Ewt(j))
           Y(j) = Y(j) + r
           fac = -hl0/r
           CALL F(Neq,TN,Y,Ftem)
           DO i = 1 , N
              Wm(i+j1) = (Ftem(i)-Savf(i))*fac
           ENDDO
           Y(j) = yj
           j1 = j1 + N
        ENDDO
        NFE = NFE + N
  ! Compute norm of Jacobian. --------------------------------------------
   300  PDNorm = DFNORM(N,Wm(3),Ewt)/ABS(hl0)
  ! Add identity matrix. -------------------------------------------------
        j = 3
        np1 = N + 1
        DO i = 1 , N
           Wm(j) = Wm(j) + 1.0D0
           j = j + np1
        ENDDO
  ! Do LU decomposition on P. --------------------------------------------
        CALL DGEFA(Wm(3),N,N,Iwm(21),ier)
        IF ( ier.NE.0 ) IERpj = 1
        RETURN
  ! Dummy block only, since MITER is never 3 in this routine. ------------
   400  RETURN
  ! If MITER = 4, call JAC and multiply by scalar. -----------------------
   500  ml = Iwm(1)
        mu = Iwm(2)
        ml3 = ml + 3
        mband = ml + mu + 1
        meband = mband + ml
        lenp = meband*N
        DO i = 1 , lenp
           Wm(i+2) = 0.0D0
        ENDDO
        CALL JAC(Neq,TN,Y,ml,mu,Wm(ml3),meband)
        con = -hl0
        DO i = 1 , lenp
           Wm(i+2) = Wm(i+2)*con
        ENDDO
        GOTO 700
  ! If MITER = 5, make MBAND calls to F to approximate J. ----------------
   600  ml = Iwm(1)
        mu = Iwm(2)
        mband = ml + mu + 1
        mba = MIN(mband,N)
        meband = mband + ml
        meb1 = meband - 1
        srur = Wm(1)
        fac = DMNORM(N,Savf,Ewt)
        r0 = 1000.0D0*ABS(H)*UROund*N*fac
        IF ( r0.EQ.0.0D0 ) r0 = 1.0D0
        DO j = 1 , mba
           DO i = j , N , mband
              yi = Y(i)
              r = MAX(srur*ABS(yi),r0/Ewt(i))
              Y(i) = Y(i) + r
           ENDDO
           CALL F(Neq,TN,Y,Ftem)
           DO jj = j , N , mband
              Y(jj) = Yh(jj,1)
              yjj = Y(jj)
              r = MAX(srur*ABS(yjj),r0/Ewt(jj))
              fac = -hl0/r
              i1 = MAX(jj-mu,1)
              i2 = MIN(jj+ml,N)
              ii = jj*meb1 - ml + 2
              DO i = i1 , i2
                 Wm(ii+i) = (Ftem(i)-Savf(i))*fac
              ENDDO
           ENDDO
        ENDDO
        NFE = NFE + mba
  ! Compute norm of Jacobian. --------------------------------------------
   700  PDNorm = DBNORM(N,Wm(ml+3),meband,ml,mu,Ewt)/ABS(hl0)
  ! Add identity matrix. -------------------------------------------------
        ii = mband + 2
        DO i = 1 , N
           Wm(ii) = Wm(ii) + 1.0D0
           ii = ii + meband
        ENDDO
  ! Do LU decomposition of P. --------------------------------------------
        CALL DGBFA(Wm(3),meband,N,ml,mu,Iwm(21),ier)
        IF ( ier.NE.0 ) IERpj = 1
  !----------------------- End of Subroutine DPRJA -----------------------
        END
  !*==DMNORM.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DMNORM
        DOUBLE PRECISION FUNCTION DMNORM(N,V,W)
        IMPLICIT NONE
  !*--DMNORM5063
  !-----------------------------------------------------------------------
  ! This function routine computes the weighted max-norm
  ! of the vector of length N contained in the array V, with weights
  ! contained in the array w of length N:
  !   DMNORM = MAX(i=1,...,N) ABS(V(i))*W(i)
  !-----------------------------------------------------------------------
        INTEGER N , i
        DOUBLE PRECISION V , W , vm
        DIMENSION V(N) , W(N)
        vm = 0.0D0
        DO i = 1 , N
           vm = MAX(vm,ABS(V(i))*W(i))
        ENDDO
        DMNORM = vm
  !----------------------- End of Function DMNORM ------------------------
        END
  !*==DFNORM.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DFNORM
        DOUBLE PRECISION FUNCTION DFNORM(N,A,W)
        IMPLICIT NONE
  !*--DFNORM5084
  !-----------------------------------------------------------------------
  ! This function computes the norm of a full N by N matrix,
  ! stored in the array A, that is consistent with the weighted max-norm
  ! on vectors, with weights stored in the array W:
  !   DFNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
  !-----------------------------------------------------------------------
        INTEGER N , i , j
        DOUBLE PRECISION A , W , an , sum
        DIMENSION A(N,N) , W(N)
        an = 0.0D0
        DO i = 1 , N
           sum = 0.0D0
           DO j = 1 , N
              sum = sum + ABS(A(i,j))/W(j)
           ENDDO
           an = MAX(an,sum*W(i))
        ENDDO
        DFNORM = an
  !----------------------- End of Function DFNORM ------------------------
        END
  !*==DBNORM.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DBNORM
        DOUBLE PRECISION FUNCTION DBNORM(N,A,Nra,Ml,Mu,W)
        IMPLICIT NONE
  !*--DBNORM5109
  !-----------------------------------------------------------------------
  ! This function computes the norm of a banded N by N matrix,
  ! stored in the array A, that is consistent with the weighted max-norm
  ! on vectors, with weights stored in the array W.
  ! ML and MU are the lower and upper half-bandwidths of the matrix.
  ! NRA is the first dimension of the A array, NRA .ge. ML+MU+1.
  ! In terms of the matrix elements a(i,j), the norm is given by:
  !   DBNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
  !-----------------------------------------------------------------------
        INTEGER N , Nra , Ml , Mu
        INTEGER i , i1 , jlo , jhi , j
        DOUBLE PRECISION A , W
        DOUBLE PRECISION an , sum
        DIMENSION A(Nra,N) , W(N)
        an = 0.0D0
        DO i = 1 , N
           sum = 0.0D0
           i1 = i + Mu + 1
           jlo = MAX(i-Ml,1)
           jhi = MIN(i+Mu,N)
           DO j = jlo , jhi
              sum = sum + ABS(A(i1-j,j))/W(j)
           ENDDO
           an = MAX(an,sum*W(i))
        ENDDO
        DBNORM = an
  !----------------------- End of Function DBNORM ------------------------
        END
  !*==DSRCMA.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DSRCMA
        SUBROUTINE DSRCMA(Rsav,Isav,Job)
        IMPLICIT NONE
  !*--DSRCMA5142
  !-----------------------------------------------------------------------
  ! This routine saves or restores (depending on JOB) the contents of
  ! the Common blocks DLS001, DLSA01, which are used
  ! internally by one or more ODEPACK solvers.
  !
  ! RSAV = real array of length 240 or more.
  ! ISAV = integer array of length 46 or more.
  ! JOB  = flag indicating to save or restore the Common blocks:
  !        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
  !        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
  !        A call with JOB = 2 presumes a prior call with JOB = 1.
  !-----------------------------------------------------------------------
        INTEGER Isav , Job
        INTEGER ILS , ILSa
        INTEGER i , lenrls , lenils , lenrla , lenila
        DOUBLE PRECISION Rsav
        DOUBLE PRECISION RLS , RLSa
        DIMENSION Rsav(*) , Isav(*)
        SAVE lenrls , lenils , lenrla , lenila
        COMMON /DLS001/ RLS(218) , ILS(37)
        COMMON /DLSA01/ RLSa(22) , ILSa(9)
        DATA lenrls/218/ , lenils/37/ , lenrla/22/ , lenila/9/
  !
        IF ( Job.EQ.2 ) THEN
  !
           DO i = 1 , lenrls
              RLS(i) = Rsav(i)
           ENDDO
           DO i = 1 , lenrla
              RLSa(i) = Rsav(lenrls+i)
           ENDDO
  !
           DO i = 1 , lenils
              ILS(i) = Isav(i)
           ENDDO
           DO i = 1 , lenila
              ILSa(i) = Isav(lenils+i)
           ENDDO
           GOTO 99999
        ENDIF
        DO i = 1 , lenrls
           Rsav(i) = RLS(i)
        ENDDO
        DO i = 1 , lenrla
           Rsav(lenrls+i) = RLSa(i)
        ENDDO
  !
        DO i = 1 , lenils
           Isav(i) = ILS(i)
        ENDDO
        DO i = 1 , lenila
           Isav(lenils+i) = ILSa(i)
        ENDDO
  !
        RETURN
  !
  !----------------------- End of Subroutine DSRCMA ----------------------
  99999 END
  !*==DRCHEK.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DRCHEK
        SUBROUTINE DRCHEK(Job,G,Neq,Y,Yh,Nyh,G0,G1,Gx,Jroot,Irt)
        IMPLICIT NONE
  !*--DRCHEK5205
        EXTERNAL G
        INTEGER Job , Neq , Nyh , Jroot , Irt
        DOUBLE PRECISION Y , Yh , G0 , G1 , Gx
        DIMENSION Neq(*) , Y(*) , Yh(Nyh,*) , G0(*) , G1(*) , Gx(*) ,     &
       &          Jroot(*)
        INTEGER IOWnd , IOWns , ICF , IERpj , IERsl , JCUr , JSTart ,     &
       &        KFLag , L , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,     &
       &        METh , MITer , MAXord , MAXcor , MSBp , MXNcf , N , NQ ,  &
       &        NST , NFE , NJE , NQU
        INTEGER IOWnd3 , IOWnr3 , IRFnd , ITAskc , NGC , NGE
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        DOUBLE PRECISION ROWnr3 , T0 , TLAst , TOUtc
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , IOWnd(6) , IOWns(6) , ICF ,    &
       &                IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , &
       &                LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,&
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        COMMON /DLSR01/ ROWnr3(2) , T0 , TLAst , TOUtc , IOWnd3(3) ,      &
       &                IOWnr3(2) , IRFnd , ITAskc , NGC , NGE
        INTEGER i , iflag , jflag
        DOUBLE PRECISION hming , t1 , temp1 , temp2 , x
        LOGICAL zroot
  !-----------------------------------------------------------------------
  ! This routine checks for the presence of a root in the vicinity of
  ! the current T, in a manner depending on the input flag JOB.  It calls
  ! Subroutine DROOTS to locate the root as precisely as possible.
  !
  ! In addition to variables described previously, DRCHEK
  ! uses the following for communication:
  ! JOB    = integer flag indicating type of call:
  !          JOB = 1 means the problem is being initialized, and DRCHEK
  !                  is to look for a root at or very near the initial T.
  !          JOB = 2 means a continuation call to the solver was just
  !                  made, and DRCHEK is to check for a root in the
  !                  relevant part of the step last taken.
  !          JOB = 3 means a successful step was just taken, and DRCHEK
  !                  is to look for a root in the interval of the step.
  ! G0     = array of length NG, containing the value of g at T = T0.
  !          G0 is input for JOB .ge. 2, and output in all cases.
  ! G1,GX  = arrays of length NG for work space.
  ! IRT    = completion flag:
  !          IRT = 0  means no root was found.
  !          IRT = -1 means JOB = 1 and a root was found too near to T.
  !          IRT = 1  means a legitimate root was found (JOB = 2 or 3).
  !                   On return, T0 is the root location, and Y is the
  !                   corresponding solution vector.
  ! T0     = value of T at one endpoint of interval of interest.  Only
  !          roots beyond T0 in the direction of integration are sought.
  !          T0 is input if JOB .ge. 2, and output in all cases.
  !          T0 is updated by DRCHEK, whether a root is found or not.
  ! TLAST  = last value of T returned by the solver (input only).
  ! TOUTC  = copy of TOUT (input only).
  ! IRFND  = input flag showing whether the last step taken had a root.
  !          IRFND = 1 if it did, = 0 if not.
  ! ITASKC = copy of ITASK (input only).
  ! NGC    = copy of NG (input only).
  !-----------------------------------------------------------------------
        Irt = 0
        DO i = 1 , NGC
           Jroot(i) = 0
        ENDDO
        hming = (ABS(TN)+ABS(H))*UROund*100.0D0
  !
        GOTO (100,200,300) , Job
  !
  ! Evaluate g at initial T, and check for zero values. ------------------
   100  T0 = TN
        CALL G(Neq,T0,Y,NGC,G0)
        NGE = 1
        zroot = .FALSE.
        DO i = 1 , NGC
           IF ( ABS(G0(i)).LE.0.0D0 ) zroot = .TRUE.
        ENDDO
        IF ( zroot ) THEN
  ! g has a zero at T.  Look at g at T + (small increment). --------------
           temp2 = MAX(hming/ABS(H),0.1D0)
           temp1 = temp2*H
           T0 = T0 + temp1
           DO i = 1 , N
              Y(i) = Y(i) + temp2*Yh(i,2)
           ENDDO
           CALL G(Neq,T0,Y,NGC,G0)
           NGE = NGE + 1
           zroot = .FALSE.
           DO i = 1 , NGC
              IF ( ABS(G0(i)).LE.0.0D0 ) zroot = .TRUE.
           ENDDO
           IF ( zroot ) THEN
  ! g has a zero at T and also close to T.  Take error return. -----------
              Irt = -1
              RETURN
           ENDIF
        ENDIF
  !
        RETURN
  !
  !
   200  IF ( IRFnd.NE.0 ) THEN
  ! If a root was found on the previous step, evaluate G0 = g(T0). -------
           CALL DINTDY(T0,0,Yh,Nyh,Y,iflag)
           CALL G(Neq,T0,Y,NGC,G0)
           NGE = NGE + 1
           zroot = .FALSE.
           DO i = 1 , NGC
              IF ( ABS(G0(i)).LE.0.0D0 ) zroot = .TRUE.
           ENDDO
           IF ( zroot ) THEN
  ! g has a zero at T0.  Look at g at T + (small increment). -------------
              temp1 = SIGN(hming,H)
              T0 = T0 + temp1
              IF ( (T0-TN)*H.LT.0.0D0 ) THEN
                 CALL DINTDY(T0,0,Yh,Nyh,Y,iflag)
              ELSE
                 temp2 = temp1/H
                 DO i = 1 , N
                    Y(i) = Y(i) + temp2*Yh(i,2)
                 ENDDO
              ENDIF
              CALL G(Neq,T0,Y,NGC,G0)
              NGE = NGE + 1
              zroot = .FALSE.
              DO i = 1 , NGC
                 IF ( ABS(G0(i)).LE.0.0D0 ) THEN
                    Jroot(i) = 1
                    zroot = .TRUE.
                 ENDIF
              ENDDO
              IF ( zroot ) THEN
  ! g has a zero at T0 and also close to T0.  Return root. ---------------
                 Irt = 1
                 RETURN
              ENDIF
           ENDIF
        ENDIF
  ! G0 has no zero components.  Proceed to check relevant interval. ------
        IF ( TN.EQ.TLAst ) GOTO 99999
  !
  ! Set T1 to TN or TOUTC, whichever comes first, and get g at T1. -------
   300  IF ( ITAskc.NE.2 .AND. ITAskc.NE.3 .AND. ITAskc.NE.5 ) THEN
           IF ( (TOUtc-TN)*H.LT.0.0D0 ) THEN
              t1 = TOUtc
              IF ( (t1-T0)*H.LE.0.0D0 ) GOTO 99999
              CALL DINTDY(t1,0,Yh,Nyh,Y,iflag)
              GOTO 400
           ENDIF
        ENDIF
        t1 = TN
        DO i = 1 , N
           Y(i) = Yh(i,1)
        ENDDO
   400  CALL G(Neq,t1,Y,NGC,G1)
        NGE = NGE + 1
  ! Call DROOTS to search for root in interval from T0 to T1. ------------
        jflag = 0
   500  CALL DROOTS(NGC,hming,jflag,T0,t1,G0,G1,Gx,x,Jroot)
        IF ( jflag.GT.1 ) THEN
           T0 = x
           CALL DCOPY(NGC,Gx,1,G0,1)
           IF ( jflag.EQ.4 ) GOTO 99999
        ELSE
           CALL DINTDY(x,0,Yh,Nyh,Y,iflag)
           CALL G(Neq,x,Y,NGC,Gx)
           NGE = NGE + 1
           GOTO 500
        ENDIF
  ! Found a root.  Interpolate to X and return. --------------------------
        CALL DINTDY(x,0,Yh,Nyh,Y,iflag)
        Irt = 1
        RETURN
  !
  !----------------------- End of Subroutine DRCHEK ----------------------
  99999 END
  !*==DROOTS.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DROOTS
        SUBROUTINE DROOTS(Ng,Hmin,Jflag,X0,X1,G0,G1,Gx,X,Jroot)
        IMPLICIT NONE
  !*--DROOTS5384
        INTEGER Ng , Jflag , Jroot
        DOUBLE PRECISION Hmin , X0 , X1 , G0 , G1 , Gx , X
        DIMENSION G0(Ng) , G1(Ng) , Gx(Ng) , Jroot(Ng)
        INTEGER IOWnd3 , IMAx , LASt , IDUm3
        DOUBLE PRECISION ALPha , X2 , RDUm3
        COMMON /DLSR01/ ALPha , X2 , RDUm3(3) , IOWnd3(3) , IMAx , LASt , &
       &                IDUm3(4)
  !-----------------------------------------------------------------------
  ! This subroutine finds the leftmost root of a set of arbitrary
  ! functions gi(x) (i = 1,...,NG) in an interval (X0,X1).  Only roots
  ! of odd multiplicity (i.e. changes of sign of the gi) are found.
  ! Here the sign of X1 - X0 is arbitrary, but is constant for a given
  ! problem, and -leftmost- means nearest to X0.
  ! The values of the vector-valued function g(x) = (gi, i=1...NG)
  ! are communicated through the call sequence of DROOTS.
  ! The method used is the Illinois algorithm.
  !
  ! Reference:
  ! Kathie L. Hiebert and Lawrence F. Shampine, Implicitly Defined
  ! Output Points for Solutions of ODEs, Sandia Report SAND80-0180,
  ! February 1980.
  !
  ! Description of parameters.
  !
  ! NG     = number of functions gi, or the number of components of
  !          the vector valued function g(x).  Input only.
  !
  ! HMIN   = resolution parameter in X.  Input only.  When a root is
  !          found, it is located only to within an error of HMIN in X.
  !          Typically, HMIN should be set to something on the order of
  !               100 * UROUND * MAX(ABS(X0),ABS(X1)),
  !          where UROUND is the unit roundoff of the machine.
  !
  ! JFLAG  = integer flag for input and output communication.
  !
  !          On input, set JFLAG = 0 on the first call for the problem,
  !          and leave it unchanged until the problem is completed.
  !          (The problem is completed when JFLAG .ge. 2 on return.)
  !
  !          On output, JFLAG has the following values and meanings:
  !          JFLAG = 1 means DROOTS needs a value of g(x).  Set GX = g(X)
  !                    and call DROOTS again.
  !          JFLAG = 2 means a root has been found.  The root is
  !                    at X, and GX contains g(X).  (Actually, X is the
  !                    rightmost approximation to the root on an interval
  !                    (X0,X1) of size HMIN or less.)
  !          JFLAG = 3 means X = X1 is a root, with one or more of the gi
  !                    being zero at X1 and no sign changes in (X0,X1).
  !                    GX contains g(X) on output.
  !          JFLAG = 4 means no roots (of odd multiplicity) were
  !                    found in (X0,X1) (no sign changes).
  !
  ! X0,X1  = endpoints of the interval where roots are sought.
  !          X1 and X0 are input when JFLAG = 0 (first call), and
  !          must be left unchanged between calls until the problem is
  !          completed.  X0 and X1 must be distinct, but X1 - X0 may be
  !          of either sign.  However, the notion of -left- and -right-
  !          will be used to mean nearer to X0 or X1, respectively.
  !          When JFLAG .ge. 2 on return, X0 and X1 are output, and
  !          are the endpoints of the relevant interval.
  !
  ! G0,G1  = arrays of length NG containing the vectors g(X0) and g(X1),
  !          respectively.  When JFLAG = 0, G0 and G1 are input and
  !          none of the G0(i) should be zero.
  !          When JFLAG .ge. 2 on return, G0 and G1 are output.
  !
  ! GX     = array of length NG containing g(X).  GX is input
  !          when JFLAG = 1, and output when JFLAG .ge. 2.
  !
  ! X      = independent variable value.  Output only.
  !          When JFLAG = 1 on output, X is the point at which g(x)
  !          is to be evaluated and loaded into GX.
  !          When JFLAG = 2 or 3, X is the root.
  !          When JFLAG = 4, X is the right endpoint of the interval, X1.
  !
  ! JROOT  = integer array of length NG.  Output only.
  !          When JFLAG = 2 or 3, JROOT indicates which components
  !          of g(x) have a root at X.  JROOT(i) is 1 if the i-th
  !          component has a root, and JROOT(i) = 0 otherwise.
  !-----------------------------------------------------------------------
        INTEGER i , imxold , nxlast
        DOUBLE PRECISION t2 , tmax , fracint , fracsub , zero , half ,    &
       &                 tenth , five
        LOGICAL zroot , sgnchg , xroot
        SAVE zero , half , tenth , five
        DATA zero/0.0D0/ , half/0.5D0/ , tenth/0.1D0/ , five/5.0D0/
  !
        IF ( Jflag.EQ.1 ) THEN
  ! Check to see in which interval g changes sign. -----------------------
           imxold = IMAx
           IMAx = 0
           tmax = zero
           zroot = .FALSE.
           DO i = 1 , Ng
              IF ( ABS(Gx(i)).LE.zero ) THEN
                 zroot = .TRUE.
  ! Neither G0(i) nor GX(i) can be zero at this point. -------------------
              ELSEIF ( SIGN(1.0D0,G0(i)).NE.SIGN(1.0D0,Gx(i)) ) THEN
                 t2 = ABS(Gx(i)/(Gx(i)-G0(i)))
                 IF ( t2.GT.tmax ) THEN
                    tmax = t2
                    IMAx = i
                 ENDIF
              ENDIF
           ENDDO
           IF ( IMAx.GT.0 ) THEN
              sgnchg = .TRUE.
           ELSE
              sgnchg = .FALSE.
              IMAx = imxold
           ENDIF
           nxlast = LASt
           IF ( sgnchg ) THEN
  ! Sign change between X0 and X2, so replace X1 with X2. ----------------
              X1 = X2
              CALL DCOPY(Ng,Gx,1,G1,1)
              LASt = 1
              xroot = .FALSE.
           ELSEIF ( .NOT.zroot ) THEN
  ! No sign change between X0 and X2.  Replace X0 with X2. ---------------
              CALL DCOPY(Ng,Gx,1,G0,1)
              X0 = X2
              LASt = 0
              xroot = .FALSE.
           ELSE
  ! Zero value at X2 and no sign change in (X0,X2), so X2 is a root. -----
              X1 = X2
              CALL DCOPY(Ng,Gx,1,G1,1)
              xroot = .TRUE.
           ENDIF
           IF ( ABS(X1-X0).LE.Hmin ) xroot = .TRUE.
        ELSE
  ! JFLAG .ne. 1.  Check for change in sign of g or zero at X1. ----------
           IMAx = 0
           tmax = zero
           zroot = .FALSE.
           DO i = 1 , Ng
              IF ( ABS(G1(i)).LE.zero ) THEN
                 zroot = .TRUE.
  ! At this point, G0(i) has been checked and cannot be zero. ------------
              ELSEIF ( SIGN(1.0D0,G0(i)).NE.SIGN(1.0D0,G1(i)) ) THEN
                 t2 = ABS(G1(i)/(G1(i)-G0(i)))
                 IF ( t2.GT.tmax ) THEN
                    tmax = t2
                    IMAx = i
                 ENDIF
              ENDIF
           ENDDO
           IF ( IMAx.GT.0 ) THEN
              sgnchg = .TRUE.
           ELSE
              sgnchg = .FALSE.
           ENDIF
           IF ( .NOT.sgnchg ) THEN
  !
  ! No sign change in the interval.  Check for zero at right endpoint. ---
              IF ( zroot ) GOTO 100
  !
  ! No sign changes in this interval.  Set X = X1, return JFLAG = 4. -----
              CALL DCOPY(Ng,G1,1,Gx,1)
              X = X1
              Jflag = 4
              GOTO 99999
           ELSE
  ! There is a sign change.  Find the first root in the interval. --------
              xroot = .FALSE.
              nxlast = 0
              LASt = 1
           ENDIF
        ENDIF
  !
  ! Repeat until the first root in the interval is found.  Loop point. ---
        IF ( xroot ) THEN
  !
  ! Return with X1 as the root.  Set JROOT.  Set X = X1 and GX = G1. -----
           Jflag = 2
           X = X1
           CALL DCOPY(Ng,G1,1,Gx,1)
           DO i = 1 , Ng
              Jroot(i) = 0
              IF ( ABS(G1(i)).GT.zero ) THEN
                 IF ( SIGN(1.0D0,G0(i)).NE.SIGN(1.0D0,G1(i)) ) Jroot(i)   &
       &              = 1
              ELSE
                 Jroot(i) = 1
              ENDIF
           ENDDO
           RETURN
        ELSE
           IF ( nxlast.NE.LASt ) THEN
              ALPha = 1.0D0
           ELSEIF ( LASt.EQ.0 ) THEN
              ALPha = 2.0D0*ALPha
           ELSE
              ALPha = 0.5D0*ALPha
           ENDIF
           X2 = X1 - (X1-X0)*G1(IMAx)/(G1(IMAx)-ALPha*G0(IMAx))
  ! If X2 is too close to X0 or X1, adjust it inward, by a fractional ----
  ! distance that is between 0.1 and 0.5. --------------------------------
           IF ( ABS(X2-X0)<half*Hmin ) THEN
              fracint = ABS(X1-X0)/Hmin
              fracsub = tenth
              IF ( fracint.LE.five ) fracsub = half/fracint
              X2 = X0 + fracsub*(X1-X0)
           ENDIF
           IF ( ABS(X1-X2)<half*Hmin ) THEN
              fracint = ABS(X1-X0)/Hmin
              fracsub = tenth
              IF ( fracint.LE.five ) fracsub = half/fracint
              X2 = X1 - fracsub*(X1-X0)
           ENDIF
           Jflag = 1
           X = X2
  ! Return to the calling routine to get a value of GX = g(X). -----------
           RETURN
        ENDIF
  !
  ! Zero value at X1 and no sign change in (X0,X1).  Return JFLAG = 3. ---
   100  X = X1
        CALL DCOPY(Ng,G1,1,Gx,1)
        DO i = 1 , Ng
           Jroot(i) = 0
           IF ( ABS(G1(i)).LE.zero ) Jroot(i) = 1
        ENDDO
        Jflag = 3
        RETURN
  !----------------------- End of Subroutine DROOTS ----------------------
  99999 END
  !*==DSRCAR.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DSRCAR
        SUBROUTINE DSRCAR(Rsav,Isav,Job)
        IMPLICIT NONE
  !*--DSRCAR5617
  !-----------------------------------------------------------------------
  ! This routine saves or restores (depending on JOB) the contents of
  ! the Common blocks DLS001, DLSA01, DLSR01, which are used
  ! internally by one or more ODEPACK solvers.
  !
  ! RSAV = real array of length 245 or more.
  ! ISAV = integer array of length 55 or more.
  ! JOB  = flag indicating to save or restore the Common blocks:
  !        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
  !        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
  !        A call with JOB = 2 presumes a prior call with JOB = 1.
  !-----------------------------------------------------------------------
        INTEGER Isav , Job
        INTEGER ILS , ILSa , ILSr
        INTEGER i , ioff , lenrls , lenils , lenrla , lenila , lenrlr ,   &
       &        lenilr
        DOUBLE PRECISION Rsav
        DOUBLE PRECISION RLS , RLSa , RLSr
        DIMENSION Rsav(*) , Isav(*)
        SAVE lenrls , lenils , lenrla , lenila , lenrlr , lenilr
        COMMON /DLS001/ RLS(218) , ILS(37)
        COMMON /DLSA01/ RLSa(22) , ILSa(9)
        COMMON /DLSR01/ RLSr(5) , ILSr(9)
        DATA lenrls/218/ , lenils/37/ , lenrla/22/ , lenila/9/
        DATA lenrlr/5/ , lenilr/9/
  !
        IF ( Job.EQ.2 ) THEN
  !
           DO i = 1 , lenrls
              RLS(i) = Rsav(i)
           ENDDO
           DO i = 1 , lenrla
              RLSa(i) = Rsav(lenrls+i)
           ENDDO
           ioff = lenrls + lenrla
           DO i = 1 , lenrlr
              RLSr(i) = Rsav(ioff+i)
           ENDDO
  !
           DO i = 1 , lenils
              ILS(i) = Isav(i)
           ENDDO
           DO i = 1 , lenila
              ILSa(i) = Isav(lenils+i)
           ENDDO
           ioff = lenils + lenila
           DO i = 1 , lenilr
              ILSr(i) = Isav(ioff+i)
           ENDDO
           GOTO 99999
        ENDIF
        DO i = 1 , lenrls
           Rsav(i) = RLS(i)
        ENDDO
        DO i = 1 , lenrla
           Rsav(lenrls+i) = RLSa(i)
        ENDDO
        ioff = lenrls + lenrla
        DO i = 1 , lenrlr
           Rsav(ioff+i) = RLSr(i)
        ENDDO
  !
        DO i = 1 , lenils
           Isav(i) = ILS(i)
        ENDDO
        DO i = 1 , lenila
           Isav(lenils+i) = ILSa(i)
        ENDDO
        ioff = lenils + lenila
        DO i = 1 , lenilr
           Isav(ioff+i) = ILSr(i)
        ENDDO
  !
        RETURN
  !
  !----------------------- End of Subroutine DSRCAR ----------------------
  99999 END
  !*==DSTODPK.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DSTODPK
        SUBROUTINE DSTODPK(Neq,Y,Yh,Nyh,Yh1,Ewt,Savf,Savx,Acor,Wm,Iwm,F,  &
       &                   JAC,PSOL)
        IMPLICIT NONE
  !*--DSTODPK5700
  !*** Start of declarations inserted by SPAG
        !INTEGER JAC
        !REAL PSOL
  !*** End of declarations inserted by SPAG
        EXTERNAL F , JAC , PSOL
        INTEGER Neq , Nyh , Iwm
        DOUBLE PRECISION Y , Yh , Yh1 , Ewt , Savf , Savx , Acor , Wm
        DIMENSION Neq(*) , Y(*) , Yh(Nyh,*) , Yh1(*) , Ewt(*) , Savf(*) , &
       &          Savx(*) , Acor(*) , Wm(*) , Iwm(*)
        INTEGER IOWnd , IALth , IPUp , LMAx , MEO , NQNyh , NSLp , ICF ,  &
       &        IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , LEWt ,  &
       &        LACor , LSAvf , LWM , LIWm , METh , MITer , MAXord ,      &
       &        MAXcor , MSBp , MXNcf , N , NQ , NST , NFE , NJE , NQU
        INTEGER JPRe , JACflg , LOCwp , LOCiwp , LSAvx , KMP , MAXl ,     &
       &        MNEwt , NNI , NLI , NPS , NCFn , NCFl
        DOUBLE PRECISION CONit , CRAte , EL , ELCo , HOLd , RMAx , TESco ,&
       &                 CCMax , EL0 , H , HMIn , HMXi , HU , RC , TN ,   &
       &                 UROund
        DOUBLE PRECISION DELt , EPCon , SQRtn , RSQrtn
        COMMON /DLS001/ CONit , CRAte , EL(13) , ELCo(13,12) , HOLd ,     &
       &                RMAx , TESco(3,12) , CCMax , EL0 , H , HMIn ,     &
       &                HMXi , HU , RC , TN , UROund , IOWnd(6) , IALth , &
       &                IPUp , LMAx , MEO , NQNyh , NSLp , ICF , IERpj ,  &
       &                IERsl , JCUr , JSTart , KFLag , L , LYH , LEWt ,  &
       &                LACor , LSAvf , LWM , LIWm , METh , MITer ,       &
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        COMMON /DLPK01/ DELt , EPCon , SQRtn , RSQrtn , JPRe , JACflg ,   &
       &                LOCwp , LOCiwp , LSAvx , KMP , MAXl , MNEwt ,     &
       &                NNI , NLI , NPS , NCFn , NCFl
  !-----------------------------------------------------------------------
  ! DSTODPK performs one step of the integration of an initial value
  ! problem for a system of Ordinary Differential Equations.
  !-----------------------------------------------------------------------
  ! The following changes were made to generate Subroutine DSTODPK
  ! from Subroutine DSTODE:
  ! 1. The array SAVX was added to the call sequence.
  ! 2. PJAC and SLVS were replaced by PSOL in the call sequence.
  ! 3. The Common block /DLPK01/ was added for communication.
  ! 4. The test constant EPCON is loaded into Common below statement
  !    numbers 125 and 155, and used below statement 400.
  ! 5. The Newton iteration counter MNEWT is set below 220 and 400.
  ! 6. The call to PJAC was replaced with a call to DPKSET (fixed name),
  !    with a longer call sequence, called depending on JACFLG.
  ! 7. The corrector residual is stored in SAVX (not Y) at 360,
  !    and the solution vector is in SAVX in the 380 loop.
  ! 8. SLVS was renamed DSOLPK and includes NEQ, SAVX, EWT, F, and JAC.
  !    SAVX was added because DSOLPK now needs Y and SAVF undisturbed.
  ! 9. The nonlinear convergence failure count NCFN is set at 430.
  !-----------------------------------------------------------------------
  ! Note: DSTODPK is independent of the value of the iteration method
  ! indicator MITER, when this is .ne. 0, and hence is independent
  ! of the type of chord method used, or the Jacobian structure.
  ! Communication with DSTODPK is done with the following variables:
  !
  ! NEQ    = integer array containing problem size in NEQ(1), and
  !          passed as the NEQ argument in all calls to F and JAC.
  ! Y      = an array of length .ge. N used as the Y argument in
  !          all calls to F and JAC.
  ! YH     = an NYH by LMAX array containing the dependent variables
  !          and their approximate scaled derivatives, where
  !          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
  !          j-th derivative of y(i), scaled by H**j/factorial(j)
  !          (j = 0,1,...,NQ).  On entry for the first step, the first
  !          two columns of YH must be set from the initial values.
  ! NYH    = a constant integer .ge. N, the first dimension of YH.
  ! YH1    = a one-dimensional array occupying the same space as YH.
  ! EWT    = an array of length N containing multiplicative weights
  !          for local error measurements.  Local errors in y(i) are
  !          compared to 1.0/EWT(i) in various error tests.
  ! SAVF   = an array of working storage, of length N.
  !          Also used for input of YH(*,MAXORD+2) when JSTART = -1
  !          and MAXORD .lt. the current order NQ.
  ! SAVX   = an array of working storage, of length N.
  ! ACOR   = a work array of length N, used for the accumulated
  !          corrections.  On a successful return, ACOR(i) contains
  !          the estimated one-step local error in y(i).
  ! WM,IWM = real and integer work arrays associated with matrix
  !          operations in chord iteration (MITER .ne. 0).
  ! CCMAX  = maximum relative change in H*EL0 before DPKSET is called.
  ! H      = the step size to be attempted on the next step.
  !          H is altered by the error control algorithm during the
  !          problem.  H can be either positive or negative, but its
  !          sign must remain constant throughout the problem.
  ! HMIN   = the minimum absolute value of the step size H to be used.
  ! HMXI   = inverse of the maximum absolute value of H to be used.
  !          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
  !          HMIN and HMXI may be changed at any time, but will not
  !          take effect until the next change of H is considered.
  ! TN     = the independent variable. TN is updated on each step taken.
  ! JSTART = an integer used for input only, with the following
  !          values and meanings:
  !               0  perform the first step.
  !           .gt.0  take a new step continuing from the last.
  !              -1  take the next step with a new value of H, MAXORD,
  !                    N, METH, MITER, and/or matrix parameters.
  !              -2  take the next step with a new value of H,
  !                    but with other inputs unchanged.
  !          On return, JSTART is set to 1 to facilitate continuation.
  ! KFLAG  = a completion code with the following meanings:
  !               0  the step was succesful.
  !              -1  the requested error could not be achieved.
  !              -2  corrector convergence could not be achieved.
  !              -3  fatal error in DPKSET or DSOLPK.
  !          A return with KFLAG = -1 or -2 means either
  !          ABS(H) = HMIN or 10 consecutive failures occurred.
  !          On a return with KFLAG negative, the values of TN and
  !          the YH array are as of the beginning of the last
  !          step, and H is the last step size attempted.
  ! MAXORD = the maximum order of integration method to be allowed.
  ! MAXCOR = the maximum number of corrector iterations allowed.
  ! MSBP   = maximum number of steps between DPKSET calls (MITER .gt. 0).
  ! MXNCF  = maximum number of convergence failures allowed.
  ! METH/MITER = the method flags.  See description in driver.
  ! N      = the number of first-order differential equations.
  !-----------------------------------------------------------------------
        INTEGER i , i1 , iredo , iret , j , jb , m , ncf , newq
        DOUBLE PRECISION dcon , ddn , del , delp , dsm , dup , exdn ,     &
       &                 exsm , exup , r , rh , rhdn , rhsm , rhup , told
  !
        KFLag = 0
        told = TN
        ncf = 0
        IERpj = 0
        IERsl = 0
        JCUr = 0
        ICF = 0
        delp = 0.0D0
        IF ( JSTart.GT.0 ) GOTO 500
        IF ( JSTart.EQ.-1 ) THEN
  !-----------------------------------------------------------------------
  ! The following block handles preliminaries needed when JSTART = -1.
  ! IPUP is set to MITER to force a matrix update.
  ! If an order increase is about to be considered (IALTH = 1),
  ! IALTH is reset to 2 to postpone consideration one more step.
  ! If the caller has changed METH, DCFODE is called to reset
  ! the coefficients of the method.
  ! If the caller has changed MAXORD to a value less than the current
  ! order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
  ! If H is to be changed, YH must be rescaled.
  ! If H or METH is being changed, IALTH is reset to L = NQ + 1
  ! to prevent further changes in H for that many steps.
  !-----------------------------------------------------------------------
           IPUp = MITer
           LMAx = MAXord + 1
           IF ( IALth.EQ.1 ) IALth = 2
           IF ( METh.NE.MEO ) THEN
              CALL DCFODE(METh,ELCo,TESco)
              MEO = METh
              IF ( NQ.LE.MAXord ) THEN
                 IALth = L
                 iret = 1
                 GOTO 100
              ENDIF
           ELSEIF ( NQ.LE.MAXord ) THEN
              GOTO 200
           ENDIF
           NQ = MAXord
           L = LMAx
           DO i = 1 , L
              EL(i) = ELCo(i,NQ)
           ENDDO
           NQNyh = NQ*Nyh
           RC = RC*EL(1)/EL0
           EL0 = EL(1)
           CONit = 0.5D0/(NQ+2)
           EPCon = CONit*TESco(2,NQ)
           ddn = DVNORM(N,Savf,Ewt)/TESco(1,L)
           exdn = 1.0D0/L
           rhdn = 1.0D0/(1.3D0*ddn**exdn+0.0000013D0)
           rh = MIN(rhdn,1.0D0)
           iredo = 3
           IF ( H.EQ.HOLd ) GOTO 300
           rh = MIN(rh,ABS(H/HOLd))
           H = HOLd
           GOTO 400
        ELSE
           IF ( JSTart.EQ.-2 ) GOTO 200
  !-----------------------------------------------------------------------
  ! On the first call, the order is set to 1, and other variables are
  ! initialized.  RMAX is the maximum ratio by which H can be increased
  ! in a single step.  It is initially 1.E4 to compensate for the small
  ! initial H, but then is normally equal to 10.  If a failure
  ! occurs (in corrector convergence or error test), RMAX is set at 2
  ! for the next increase.
  !-----------------------------------------------------------------------
           LMAx = MAXord + 1
           NQ = 1
           L = 2
           IALth = 2
           RMAx = 10000.0D0
           RC = 0.0D0
           EL0 = 1.0D0
           CRAte = 0.7D0
           HOLd = H
           MEO = METh
           NSLp = 0
           IPUp = MITer
           iret = 3
  !-----------------------------------------------------------------------
  ! DCFODE is called to get all the integration coefficients for the
  ! current METH.  Then the EL vector and related constants are reset
  ! whenever the order NQ is changed, or at the start of the problem.
  !-----------------------------------------------------------------------
           CALL DCFODE(METh,ELCo,TESco)
        ENDIF
   100  DO i = 1 , L
           EL(i) = ELCo(i,NQ)
        ENDDO
        NQNyh = NQ*Nyh
        RC = RC*EL(1)/EL0
        EL0 = EL(1)
        CONit = 0.5D0/(NQ+2)
        EPCon = CONit*TESco(2,NQ)
        GOTO (200,300,500) , iret
  !-----------------------------------------------------------------------
  ! If H is being changed, the H ratio RH is checked against
  ! RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
  ! L = NQ + 1 to prevent a change of H for that many steps, unless
  ! forced by a convergence or error test failure.
  !-----------------------------------------------------------------------
   200  IF ( H.EQ.HOLd ) GOTO 500
        rh = H/HOLd
        H = HOLd
        iredo = 3
        GOTO 400
   300  rh = MAX(rh,HMIn/ABS(H))
   400  rh = MIN(rh,RMAx)
        rh = rh/MAX(1.0D0,ABS(H)*HMXi*rh)
        r = 1.0D0
        DO j = 2 , L
           r = r*rh
           DO i = 1 , N
              Yh(i,j) = Yh(i,j)*r
           ENDDO
        ENDDO
        H = H*rh
        RC = RC*rh
        IALth = L
        IF ( iredo.EQ.0 ) THEN
           RMAx = 10.0D0
           GOTO 1300
        ENDIF
  !-----------------------------------------------------------------------
  ! This section computes the predicted values by effectively
  ! multiplying the YH array by the Pascal triangle matrix.
  ! The flag IPUP is set according to whether matrix data is involved
  ! (JACFLG .ne. 0) or not (JACFLG = 0), to trigger a call to DPKSET.
  ! IPUP is set to MITER when RC differs from 1 by more than CCMAX,
  ! and at least every MSBP steps, when JACFLG = 1.
  ! RC is the ratio of new to old values of the coefficient  H*EL(1).
  !-----------------------------------------------------------------------
   500  IF ( JACflg.NE.0 ) THEN
           IF ( ABS(RC-1.0D0).GT.CCMax ) IPUp = MITer
           IF ( NST.GE.NSLp+MSBp ) IPUp = MITer
        ELSE
           IPUp = 0
           CRAte = 0.7D0
        ENDIF
        TN = TN + H
        i1 = NQNyh + 1
        DO jb = 1 , NQ
           i1 = i1 - Nyh
  !DIR$ IVDEP
           DO i = i1 , NQNyh
              Yh1(i) = Yh1(i) + Yh1(i+Nyh)
           ENDDO
        ENDDO
  !-----------------------------------------------------------------------
  ! Up to MAXCOR corrector iterations are taken.  A convergence test is
  ! made on the RMS-norm of each correction, weighted by the error
  ! weight vector EWT.  The sum of the corrections is accumulated in the
  ! vector ACOR(i).  The YH array is not altered in the corrector loop.
  !-----------------------------------------------------------------------
   600  m = 0
        MNEwt = 0
        DO i = 1 , N
           Y(i) = Yh(i,1)
        ENDDO
        CALL F(Neq,TN,Y,Savf)
        NFE = NFE + 1
        IF ( IPUp.GT.0 ) THEN
  !-----------------------------------------------------------------------
  ! If indicated, DPKSET is called to update any matrix data needed,
  ! before starting the corrector iteration.
  ! IPUP is set to 0 as an indicator that this has been done.
  !-----------------------------------------------------------------------
           CALL DPKSET(Neq,Y,Yh1,Ewt,Acor,Savf,Wm,Iwm,F,JAC)
           IPUp = 0
           RC = 1.0D0
           NSLp = NST
           CRAte = 0.7D0
           IF ( IERpj.NE.0 ) GOTO 900
        ENDIF
        DO i = 1 , N
           Acor(i) = 0.0D0
        ENDDO
   700  IF ( MITer.NE.0 ) THEN
  !-----------------------------------------------------------------------
  ! In the case of the chord method, compute the corrector error,
  ! and solve the linear system with that as right-hand side and
  ! P as coefficient matrix.
  !-----------------------------------------------------------------------
           DO i = 1 , N
              Savx(i) = H*Savf(i) - (Yh(i,2)+Acor(i))
           ENDDO
           CALL DSOLPK(Neq,Y,Savf,Savx,Ewt,Wm,Iwm,F,PSOL)
           IF ( IERsl.LT.0 ) GOTO 900
           IF ( IERsl.GT.0 ) GOTO 800
           del = DVNORM(N,Savx,Ewt)
           DO i = 1 , N
              Acor(i) = Acor(i) + Savx(i)
              Y(i) = Yh(i,1) + EL(1)*Acor(i)
           ENDDO
        ELSE
  !-----------------------------------------------------------------------
  ! In the case of functional iteration, update Y directly from
  ! the result of the last function evaluation.
  !-----------------------------------------------------------------------
           DO i = 1 , N
              Savf(i) = H*Savf(i) - Yh(i,2)
              Y(i) = Savf(i) - Acor(i)
           ENDDO
           del = DVNORM(N,Y,Ewt)
           DO i = 1 , N
              Y(i) = Yh(i,1) + EL(1)*Savf(i)
              Acor(i) = Savf(i)
           ENDDO
        ENDIF
  !-----------------------------------------------------------------------
  ! Test for convergence.  If M .gt. 0, an estimate of the convergence
  ! rate constant is stored in CRATE, and this is used in the test.
  !-----------------------------------------------------------------------
        IF ( m.NE.0 ) CRAte = MAX(0.2D0*CRAte,del/delp)
        dcon = del*MIN(1.0D0,1.5D0*CRAte)/EPCon
        IF ( dcon.LE.1.0D0 ) THEN
  !-----------------------------------------------------------------------
  ! The corrector has converged.  JCUR is set to 0
  ! to signal that the Jacobian involved may need updating later.
  ! The local error test is made and control passes to statement 500
  ! if it fails.
  !-----------------------------------------------------------------------
           JCUr = 0
           IF ( m.EQ.0 ) dsm = del/TESco(2,NQ)
           IF ( m.GT.0 ) dsm = DVNORM(N,Acor,Ewt)/TESco(2,NQ)
           IF ( dsm.GT.1.0D0 ) THEN
  !-----------------------------------------------------------------------
  ! The error test failed.  KFLAG keeps track of multiple failures.
  ! Restore TN and the YH array to their previous values, and prepare
  ! to try the step again.  Compute the optimum step size for this or
  ! one lower order.  After 2 or more failures, H is forced to decrease
  ! by a factor of 0.2 or less.
  !-----------------------------------------------------------------------
              KFLag = KFLag - 1
              TN = told
              i1 = NQNyh + 1
              DO jb = 1 , NQ
                 i1 = i1 - Nyh
  !DIR$ IVDEP
                 DO i = i1 , NQNyh
                    Yh1(i) = Yh1(i) - Yh1(i+Nyh)
                 ENDDO
              ENDDO
              RMAx = 2.0D0
              IF ( ABS(H).LE.HMIn*1.00001D0 ) THEN
  !-----------------------------------------------------------------------
  ! All returns are made through this section.  H is saved in HOLD
  ! to allow the caller to change H on the next step.
  !-----------------------------------------------------------------------
                 KFLag = -1
                 GOTO 1400
              ELSEIF ( KFLag.LE.-3 ) THEN
  !-----------------------------------------------------------------------
  ! Control reaches this section if 3 or more failures have occured.
  ! If 10 failures have occurred, exit with KFLAG = -1.
  ! It is assumed that the derivatives that have accumulated in the
  ! YH array have errors of the wrong order.  Hence the first
  ! derivative is recomputed, and the order is set to 1.  Then
  ! H is reduced by a factor of 10, and the step is retried,
  ! until it succeeds or H reaches HMIN.
  !-----------------------------------------------------------------------
                 IF ( KFLag.EQ.-10 ) THEN
                    KFLag = -1
                    GOTO 1400
                 ELSE
                    rh = 0.1D0
                    rh = MAX(HMIn/ABS(H),rh)
                    H = H*rh
                    DO i = 1 , N
                       Y(i) = Yh(i,1)
                    ENDDO
                    CALL F(Neq,TN,Y,Savf)
                    NFE = NFE + 1
                    DO i = 1 , N
                       Yh(i,2) = H*Savf(i)
                    ENDDO
                    IPUp = MITer
                    IALth = 5
                    IF ( NQ.EQ.1 ) GOTO 500
                    NQ = 1
                    L = 2
                    iret = 3
                    GOTO 100
                 ENDIF
              ELSE
                 iredo = 2
                 rhup = 0.0D0
                 GOTO 1000
              ENDIF
           ELSE
  !-----------------------------------------------------------------------
  ! After a successful step, update the YH array.
  ! Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
  ! If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
  ! use in a possible order increase on the next step.
  ! If a change in H is considered, an increase or decrease in order
  ! by one is considered also.  A change in H is made only if it is by a
  ! factor of at least 1.1.  If not, IALTH is set to 3 to prevent
  ! testing for that many steps.
  !-----------------------------------------------------------------------
              KFLag = 0
              iredo = 0
              NST = NST + 1
              HU = H
              NQU = NQ
              DO j = 1 , L
                 DO i = 1 , N
                    Yh(i,j) = Yh(i,j) + EL(j)*Acor(i)
                 ENDDO
              ENDDO
              IALth = IALth - 1
              IF ( IALth.EQ.0 ) THEN
  !-----------------------------------------------------------------------
  ! Regardless of the success or failure of the step, factors
  ! RHDN, RHSM, and RHUP are computed, by which H could be multiplied
  ! at order NQ - 1, order NQ, or order NQ + 1, respectively.
  ! In the case of failure, RHUP = 0.0 to avoid an order increase.
  ! the largest of these is determined and the new order chosen
  ! accordingly.  If the order is to be increased, we compute one
  ! additional scaled derivative.
  !-----------------------------------------------------------------------
                 rhup = 0.0D0
                 IF ( L.NE.LMAx ) THEN
                    DO i = 1 , N
                       Savf(i) = Acor(i) - Yh(i,LMAx)
                    ENDDO
                    dup = DVNORM(N,Savf,Ewt)/TESco(3,NQ)
                    exup = 1.0D0/(L+1)
                    rhup = 1.0D0/(1.4D0*dup**exup+0.0000014D0)
                 ENDIF
                 GOTO 1000
              ELSE
                 IF ( IALth.LE.1 ) THEN
                    IF ( L.NE.LMAx ) THEN
                       DO i = 1 , N
                          Yh(i,LMAx) = Acor(i)
                       ENDDO
                    ENDIF
                 ENDIF
                 GOTO 1300
              ENDIF
           ENDIF
        ELSE
           m = m + 1
           IF ( m.NE.MAXcor ) THEN
              IF ( m.LT.2 .OR. del.LE.2.0D0*delp ) THEN
                 MNEwt = m
                 delp = del
                 CALL F(Neq,TN,Y,Savf)
                 NFE = NFE + 1
                 GOTO 700
              ENDIF
           ENDIF
        ENDIF
  !-----------------------------------------------------------------------
  ! The corrector iteration failed to converge.
  ! If MITER .ne. 0 and the Jacobian is out of date, DPKSET is called for
  ! the next try.  Otherwise the YH array is retracted to its values
  ! before prediction, and H is reduced, if possible.  If H cannot be
  ! reduced or MXNCF failures have occurred, exit with KFLAG = -2.
  !-----------------------------------------------------------------------
   800  IF ( MITer.NE.0 .AND. JCUr.NE.1 .AND. JACflg.NE.0 ) THEN
           ICF = 1
           IPUp = MITer
           GOTO 600
        ENDIF
   900  ICF = 2
        ncf = ncf + 1
        NCFn = NCFn + 1
        RMAx = 2.0D0
        TN = told
        i1 = NQNyh + 1
        DO jb = 1 , NQ
           i1 = i1 - Nyh
  !DIR$ IVDEP
           DO i = i1 , NQNyh
              Yh1(i) = Yh1(i) - Yh1(i+Nyh)
           ENDDO
        ENDDO
        IF ( IERpj.LT.0 .OR. IERsl.LT.0 ) THEN
           KFLag = -3
           GOTO 1400
        ELSEIF ( ABS(H).LE.HMIn*1.00001D0 ) THEN
           KFLag = -2
           GOTO 1400
        ELSEIF ( ncf.EQ.MXNcf ) THEN
           KFLag = -2
           GOTO 1400
        ELSE
           rh = 0.5D0
           IPUp = MITer
           iredo = 1
           GOTO 300
        ENDIF
   1000 exsm = 1.0D0/L
        rhsm = 1.0D0/(1.2D0*dsm**exsm+0.0000012D0)
        rhdn = 0.0D0
        IF ( NQ.NE.1 ) THEN
           ddn = DVNORM(N,Yh(1,L),Ewt)/TESco(1,NQ)
           exdn = 1.0D0/NQ
           rhdn = 1.0D0/(1.3D0*ddn**exdn+0.0000013D0)
        ENDIF
        IF ( rhsm.GE.rhup ) THEN
           IF ( rhsm.GE.rhdn ) THEN
              newq = NQ
              rh = rhsm
              GOTO 1100
           ENDIF
        ELSEIF ( rhup.GT.rhdn ) THEN
           newq = L
           rh = rhup
           IF ( rh.LT.1.1D0 ) THEN
              IALth = 3
              GOTO 1300
           ELSE
              r = EL(L)/L
              DO i = 1 , N
                 Yh(i,newq+1) = Acor(i)*r
              ENDDO
              GOTO 1200
           ENDIF
        ENDIF
        newq = NQ - 1
        rh = rhdn
        IF ( KFLag.LT.0 .AND. rh.GT.1.0D0 ) rh = 1.0D0
   1100 IF ( (KFLag.EQ.0) .AND. (rh.LT.1.1D0) ) THEN
           IALth = 3
           GOTO 1300
        ELSE
           IF ( KFLag.LE.-2 ) rh = MIN(rh,0.2D0)
  !-----------------------------------------------------------------------
  ! If there is a change of order, reset NQ, L, and the coefficients.
  ! In any case H is reset according to RH and the YH array is rescaled.
  ! Then exit from 690 if the step was OK, or redo the step otherwise.
  !-----------------------------------------------------------------------
           IF ( newq.EQ.NQ ) GOTO 300
        ENDIF
   1200 NQ = newq
        L = NQ + 1
        iret = 2
        GOTO 100
   1300 r = 1.0D0/TESco(2,NQU)
        DO i = 1 , N
           Acor(i) = Acor(i)*r
        ENDDO
   1400 HOLd = H
        JSTart = 1
  !----------------------- End of Subroutine DSTODPK ---------------------
        END
  !*==DPKSET.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DPKSET
        SUBROUTINE DPKSET(Neq,Y,Ysv,Ewt,Ftem,Savf,Wm,Iwm,F,JAC)
        IMPLICIT NONE
  !*--DPKSET6275
  !*** Start of declarations inserted by SPAG
        !REAL F
  !*** End of declarations inserted by SPAG
        EXTERNAL F , JAC
        INTEGER Neq , Iwm
        DOUBLE PRECISION Y , Ysv , Ewt , Ftem , Savf , Wm
        DIMENSION Neq(*) , Y(*) , Ysv(*) , Ewt(*) , Ftem(*) , Savf(*) ,   &
       &          Wm(*) , Iwm(*)
        INTEGER IOWnd , IOWns , ICF , IERpj , IERsl , JCUr , JSTart ,     &
       &        KFLag , L , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,     &
       &        METh , MITer , MAXord , MAXcor , MSBp , MXNcf , N , NQ ,  &
       &        NST , NFE , NJE , NQU
        INTEGER JPRe , JACflg , LOCwp , LOCiwp , LSAvx , KMP , MAXl ,     &
       &        MNEwt , NNI , NLI , NPS , NCFn , NCFl
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        DOUBLE PRECISION DELt , EPCon , SQRtn , RSQrtn
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , IOWnd(6) , IOWns(6) , ICF ,    &
       &                IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , &
       &                LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,&
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        COMMON /DLPK01/ DELt , EPCon , SQRtn , RSQrtn , JPRe , JACflg ,   &
       &                LOCwp , LOCiwp , LSAvx , KMP , MAXl , MNEwt ,     &
       &                NNI , NLI , NPS , NCFn , NCFl
  !-----------------------------------------------------------------------
  ! DPKSET is called by DSTODPK to interface with the user-supplied
  ! routine JAC, to compute and process relevant parts of
  ! the matrix P = I - H*EL(1)*J , where J is the Jacobian df/dy,
  ! as need for preconditioning matrix operations later.
  !
  ! In addition to variables described previously, communication
  ! with DPKSET uses the following:
  ! Y     = array containing predicted values on entry.
  ! YSV   = array containing predicted y, to be saved (YH1 in DSTODPK).
  ! FTEM  = work array of length N (ACOR in DSTODPK).
  ! SAVF  = array containing f evaluated at predicted y.
  ! WM    = real work space for matrices.
  !         Space for preconditioning data starts at WM(LOCWP).
  ! IWM   = integer work space.
  !         Space for preconditioning data starts at IWM(LOCIWP).
  ! IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
  !         JAC returned an error flag.
  ! JCUR  = output flag = 1 to indicate that the Jacobian matrix
  !         (or approximation) is now current.
  ! This routine also uses Common variables EL0, H, TN, IERPJ, JCUR, NJE.
  !-----------------------------------------------------------------------
        INTEGER ier
        DOUBLE PRECISION hl0
  !
        IERpj = 0
        JCUr = 1
        hl0 = EL0*H
        CALL JAC(F,Neq,TN,Y,Ysv,Ewt,Savf,Ftem,hl0,Wm(LOCwp),Iwm(LOCiwp),  &
       &         ier)
        NJE = NJE + 1
        IF ( ier.EQ.0 ) RETURN
        IERpj = 1
  !----------------------- End of Subroutine DPKSET ----------------------
        END
  !*==DSOLPK.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DSOLPK
        SUBROUTINE DSOLPK(Neq,Y,Savf,X,Ewt,Wm,Iwm,F,PSOL)
        IMPLICIT NONE
  !*--DSOLPK6341
  !*** Start of declarations inserted by SPAG
        !REAL F , PSOL
  !*** End of declarations inserted by SPAG
        EXTERNAL F , PSOL
        INTEGER Neq , Iwm
        DOUBLE PRECISION Y , Savf , X , Ewt , Wm
        DIMENSION Neq(*) , Y(*) , Savf(*) , X(*) , Ewt(*) , Wm(*) , Iwm(*)
        INTEGER IOWnd , IOWns , ICF , IERpj , IERsl , JCUr , JSTart ,     &
       &        KFLag , L , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,     &
       &        METh , MITer , MAXord , MAXcor , MSBp , MXNcf , N , NQ ,  &
       &        NST , NFE , NJE , NQU
        INTEGER JPRe , JACflg , LOCwp , LOCiwp , LSAvx , KMP , MAXl ,     &
       &        MNEwt , NNI , NLI , NPS , NCFn , NCFl
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        DOUBLE PRECISION DELt , EPCon , SQRtn , RSQrtn
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , IOWnd(6) , IOWns(6) , ICF ,    &
       &                IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , &
       &                LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,&
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        COMMON /DLPK01/ DELt , EPCon , SQRtn , RSQrtn , JPRe , JACflg ,   &
       &                LOCwp , LOCiwp , LSAvx , KMP , MAXl , MNEwt ,     &
       &                NNI , NLI , NPS , NCFn , NCFl
  !-----------------------------------------------------------------------
  ! This routine interfaces to one of DSPIOM, DSPIGMR, DPCG, DPCGS, or
  ! DUSOL, for the solution of the linear system arising from a Newton
  ! iteration.  It is called if MITER .ne. 0.
  ! In addition to variables described elsewhere,
  ! communication with DSOLPK uses the following variables:
  ! WM    = real work space containing data for the algorithm
  !         (Krylov basis vectors, Hessenberg matrix, etc.)
  ! IWM   = integer work space containing data for the algorithm
  ! X     = the right-hand side vector on input, and the solution vector
  !         on output, of length N.
  ! IERSL = output flag (in Common):
  !         IERSL =  0 means no trouble occurred.
  !         IERSL =  1 means the iterative method failed to converge.
  !                    If the preconditioner is out of date, the step
  !                    is repeated with a new preconditioner.
  !                    Otherwise, the stepsize is reduced (forcing a
  !                    new evaluation of the preconditioner) and the
  !                    step is repeated.
  !         IERSL = -1 means there was a nonrecoverable error in the
  !                    iterative solver, and an error exit occurs.
  ! This routine also uses the Common variables TN, EL0, H, N, MITER,
  ! DELT, EPCON, SQRTN, RSQRTN, MAXL, KMP, MNEWT, NNI, NLI, NPS, NCFL,
  ! LOCWP, LOCIWP.
  !-----------------------------------------------------------------------
        INTEGER iflag , lb , ldl , lhes , liom , lgmr , lpcg , lp , lq ,  &
       &        lr , lv , lw , lwk , lz , maxlp1 , npsl
        DOUBLE PRECISION delta , hl0
  !
        IERsl = 0
        hl0 = H*EL0
        delta = DELt*EPCon
        GOTO (100,200,300,400,500,500,500,500,500) , MITer
  !-----------------------------------------------------------------------
  ! Use the SPIOM algorithm to solve the linear system P*x = -f.
  !-----------------------------------------------------------------------
   100  lv = 1
        lb = lv + N*MAXl
        lhes = lb + N
        lwk = lhes + MAXl*MAXl
        CALL DCOPY(N,X,1,Wm(lb),1)
        CALL DSCAL(N,RSQrtn,Ewt,1)
        CALL DSPIOM(Neq,TN,Y,Savf,Wm(lb),Ewt,N,MAXl,KMP,delta,hl0,JPRe,   &
       &            MNEwt,F,PSOL,npsl,X,Wm(lv),Wm(lhes),Iwm,liom,Wm(LOCwp)&
       &            ,Iwm(LOCiwp),Wm(lwk),iflag)
        NNI = NNI + 1
        NLI = NLI + liom
        NPS = NPS + npsl
        CALL DSCAL(N,SQRtn,Ewt,1)
        IF ( iflag.NE.0 ) NCFl = NCFl + 1
        IF ( iflag.GE.2 ) IERsl = 1
        IF ( iflag.LT.0 ) IERsl = -1
        RETURN
  !-----------------------------------------------------------------------
  ! Use the SPIGMR algorithm to solve the linear system P*x = -f.
  !-----------------------------------------------------------------------
   200  maxlp1 = MAXl + 1
        lv = 1
        lb = lv + N*MAXl
        lhes = lb + N + 1
        lq = lhes + MAXl*maxlp1
        lwk = lq + 2*MAXl
        ldl = lwk + MIN(1,MAXl-KMP)*N
        CALL DCOPY(N,X,1,Wm(lb),1)
        CALL DSCAL(N,RSQrtn,Ewt,1)
        CALL DSPIGMR(Neq,TN,Y,Savf,Wm(lb),Ewt,N,MAXl,maxlp1,KMP,delta,hl0,&
       &             JPRe,MNEwt,F,PSOL,npsl,X,Wm(lv),Wm(lhes),Wm(lq),lgmr,&
       &             Wm(LOCwp),Iwm(LOCiwp),Wm(lwk),Wm(ldl),iflag)
        NNI = NNI + 1
        NLI = NLI + lgmr
        NPS = NPS + npsl
        CALL DSCAL(N,SQRtn,Ewt,1)
        IF ( iflag.NE.0 ) NCFl = NCFl + 1
        IF ( iflag.GE.2 ) IERsl = 1
        IF ( iflag.LT.0 ) IERsl = -1
        RETURN
  !-----------------------------------------------------------------------
  ! Use DPCG to solve the linear system P*x = -f
  !-----------------------------------------------------------------------
   300  lr = 1
        lp = lr + N
        lw = lp + N
        lz = lw + N
        lwk = lz + N
        CALL DCOPY(N,X,1,Wm(lr),1)
        CALL DPCG(Neq,TN,Y,Savf,Wm(lr),Ewt,N,MAXl,delta,hl0,JPRe,MNEwt,F, &
       &          PSOL,npsl,X,Wm(lp),Wm(lw),Wm(lz),lpcg,Wm(LOCwp),        &
       &          Iwm(LOCiwp),Wm(lwk),iflag)
        NNI = NNI + 1
        NLI = NLI + lpcg
        NPS = NPS + npsl
        IF ( iflag.NE.0 ) NCFl = NCFl + 1
        IF ( iflag.GE.2 ) IERsl = 1
        IF ( iflag.LT.0 ) IERsl = -1
        RETURN
  !-----------------------------------------------------------------------
  ! Use DPCGS to solve the linear system P*x = -f
  !-----------------------------------------------------------------------
   400  lr = 1
        lp = lr + N
        lw = lp + N
        lz = lw + N
        lwk = lz + N
        CALL DCOPY(N,X,1,Wm(lr),1)
        CALL DPCGS(Neq,TN,Y,Savf,Wm(lr),Ewt,N,MAXl,delta,hl0,JPRe,MNEwt,F,&
       &           PSOL,npsl,X,Wm(lp),Wm(lw),Wm(lz),lpcg,Wm(LOCwp),       &
       &           Iwm(LOCiwp),Wm(lwk),iflag)
        NNI = NNI + 1
        NLI = NLI + lpcg
        NPS = NPS + npsl
        IF ( iflag.NE.0 ) NCFl = NCFl + 1
        IF ( iflag.GE.2 ) IERsl = 1
        IF ( iflag.LT.0 ) IERsl = -1
        RETURN
  !-----------------------------------------------------------------------
  ! Use DUSOL, which interfaces to PSOL, to solve the linear system
  ! (no Krylov iteration).
  !-----------------------------------------------------------------------
   500  lb = 1
        lwk = lb + N
        CALL DCOPY(N,X,1,Wm(lb),1)
        CALL DUSOL(Neq,TN,Y,Savf,Wm(lb),Ewt,N,delta,hl0,MNEwt,PSOL,npsl,X,&
       &           Wm(LOCwp),Iwm(LOCiwp),Wm(lwk),iflag)
        NNI = NNI + 1
        NPS = NPS + npsl
        IF ( iflag.NE.0 ) NCFl = NCFl + 1
        IF ( iflag.EQ.3 ) IERsl = 1
        IF ( iflag.LT.0 ) IERsl = -1
  !----------------------- End of Subroutine DSOLPK ----------------------
        END
  !*==DSPIOM.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DSPIOM
        SUBROUTINE DSPIOM(Neq,Tn,Y,Savf,B,Wght,N,Maxl,Kmp,Delta,Hl0,Jpre, &
       &                  Mnewt,F,PSOL,Npsl,X,V,Hes,Ipvt,Liom,Wp,Iwp,Wk,  &
       &                  Iflag)
        IMPLICIT NONE
  !*--DSPIOM6503
  !*** Start of declarations inserted by SPAG
        !REAL F
  !*** End of declarations inserted by SPAG
        EXTERNAL F , PSOL
        INTEGER Neq , N , Maxl , Kmp , Jpre , Mnewt , Npsl , Ipvt , Liom ,&
       &        Iwp , Iflag
        DOUBLE PRECISION Tn , Y , Savf , B , Wght , Delta , Hl0 , X , V , &
       &                 Hes , Wp , Wk
        DIMENSION Neq(*) , Y(*) , Savf(*) , B(*) , Wght(*) , X(*) , V(N,*)&
       &          , Hes(Maxl,Maxl) , Ipvt(*) , Wp(*) , Iwp(*) , Wk(*)
  !-----------------------------------------------------------------------
  ! This routine solves the linear system A * x = b using a scaled
  ! preconditioned version of the Incomplete Orthogonalization Method.
  ! An initial guess of x = 0 is assumed.
  !-----------------------------------------------------------------------
  !
  !      On entry
  !
  !          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
  !
  !           TN = current value of t.
  !
  !            Y = array containing current dependent variable vector.
  !
  !         SAVF = array containing current value of f(t,y).
  !
  !         B    = the right hand side of the system A*x = b.
  !                B is also used as work space when computing the
  !                final approximation.
  !                (B is the same as V(*,MAXL+1) in the call to DSPIOM.)
  !
  !         WGHT = array of length N containing scale factors.
  !                1/WGHT(i) are the diagonal elements of the diagonal
  !                scaling matrix D.
  !
  !         N    = the order of the matrix A, and the lengths
  !                of the vectors Y, SAVF, B, WGHT, and X.
  !
  !         MAXL = the maximum allowable order of the matrix HES.
  !
  !          KMP = the number of previous vectors the new vector VNEW
  !                must be made orthogonal to.  KMP .le. MAXL.
  !
  !        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
  !
  !          HL0 = current value of (step size h) * (coefficient l0).
  !
  !         JPRE = preconditioner type flag.
  !
  !        MNEWT = Newton iteration counter (.ge. 0).
  !
  !           WK = real work array of length N used by DATV and PSOL.
  !
  !           WP = real work array used by preconditioner PSOL.
  !
  !          IWP = integer work array used by preconditioner PSOL.
  !
  !      On return
  !
  !         X    = the final computed approximation to the solution
  !                of the system A*x = b.
  !
  !         V    = the N by (LIOM+1) array containing the LIOM
  !                orthogonal vectors V(*,1) to V(*,LIOM).
  !
  !         HES  = the LU factorization of the LIOM by LIOM upper
  !                Hessenberg matrix whose entries are the
  !                scaled inner products of A*V(*,k) and V(*,i).
  !
  !         IPVT = an integer array containg pivoting information.
  !                It is loaded in DHEFA and used in DHESL.
  !
  !         LIOM = the number of iterations performed, and current
  !                order of the upper Hessenberg matrix HES.
  !
  !         NPSL = the number of calls to PSOL.
  !
  !        IFLAG = integer error flag:
  !                0 means convergence in LIOM iterations, LIOM.le.MAXL.
  !                1 means the convergence test did not pass in MAXL
  !                  iterations, but the residual norm is .lt. 1,
  !                  or .lt. norm(b) if MNEWT = 0, and so X is computed.
  !                2 means the convergence test did not pass in MAXL
  !                  iterations, residual .gt. 1, and X is undefined.
  !                3 means there was a recoverable error in PSOL
  !                  caused by the preconditioner being out of date.
  !               -1 means there was a nonrecoverable error in PSOL.
  !
  !-----------------------------------------------------------------------
        INTEGER i , ier , info , j , k , ll , lm1
        DOUBLE PRECISION bnrm , bnrm0 , prod , rho , snormw , tem
  !
        Iflag = 0
        Liom = 0
        Npsl = 0
  !-----------------------------------------------------------------------
  ! The initial residual is the vector b.  Apply scaling to b, and test
  ! for an immediate return with X = 0 or X = b.
  !-----------------------------------------------------------------------
        DO i = 1 , N
           V(i,1) = B(i)*Wght(i)
        ENDDO
        bnrm0 = DNRM2(N,V,1)
        bnrm = bnrm0
        IF ( bnrm0.GT.Delta ) THEN
  ! Apply inverse of left preconditioner to vector b. --------------------
           ier = 0
           IF ( Jpre.NE.0 .AND. Jpre.NE.2 ) THEN
              CALL PSOL(Neq,Tn,Y,Savf,Wk,Hl0,Wp,Iwp,B,1,ier)
              Npsl = 1
              IF ( ier.NE.0 ) GOTO 300
  ! Calculate norm of scaled vector V(*,1) and normalize it. -------------
              DO i = 1 , N
                 V(i,1) = B(i)*Wght(i)
              ENDDO
              bnrm = DNRM2(N,V,1)
              Delta = Delta*(bnrm/bnrm0)
           ENDIF
           tem = 1.0D0/bnrm
           CALL DSCAL(N,tem,V(1,1),1)
  ! Zero out the HES array. ----------------------------------------------
           DO j = 1 , Maxl
              DO i = 1 , Maxl
                 Hes(i,j) = 0.0D0
              ENDDO
           ENDDO
  !-----------------------------------------------------------------------
  ! Main loop on LL = l to compute the vectors V(*,2) to V(*,MAXL).
  ! The running product PROD is needed for the convergence test.
  !-----------------------------------------------------------------------
           prod = 1.0D0
           DO ll = 1 , Maxl
              Liom = ll
  !-----------------------------------------------------------------------
  ! Call routine DATV to compute VNEW = Abar*v(l), where Abar is
  ! the matrix A with scaling and inverse preconditioner factors applied.
  ! Call routine DORTHOG to orthogonalize the new vector vnew = V(*,l+1).
  ! Call routine DHEFA to update the factors of HES.
  !-----------------------------------------------------------------------
              CALL DATV(Neq,Y,Savf,V(1,ll),Wght,X,F,PSOL,V(1,ll+1),Wk,Wp, &
       &                Iwp,Hl0,Jpre,ier,Npsl)
              IF ( ier.NE.0 ) GOTO 300
              CALL DORTHOG(V(1,ll+1),V,Hes,N,ll,Maxl,Kmp,snormw)
              CALL DHEFA(Hes,Maxl,ll,Ipvt,info,ll)
              lm1 = ll - 1
              IF ( ll.GT.1 .AND. Ipvt(lm1).EQ.lm1 )                       &
       &           prod = prod*Hes(ll,lm1)
              IF ( info.NE.ll ) THEN
  !-----------------------------------------------------------------------
  ! Update RHO, the estimate of the norm of the residual b - A*x(l).
  ! test for convergence.  If passed, compute approximation x(l).
  ! If failed and l .lt. MAXL, then continue iterating.
  !-----------------------------------------------------------------------
                 rho = bnrm*snormw*ABS(prod/Hes(ll,ll))
                 IF ( rho.LE.Delta ) GOTO 200
                 IF ( ll.EQ.Maxl ) GOTO 50
              ELSE
  !-----------------------------------------------------------------------
  ! The last pivot in HES was found to be zero.
  ! If vnew = 0 or l = MAXL, take an error return with IFLAG = 2.
  ! otherwise, continue the iteration without a convergence test.
  !-----------------------------------------------------------------------
                 IF ( snormw.EQ.0.0D0 ) GOTO 100
                 IF ( ll.EQ.Maxl ) GOTO 100
              ENDIF
  ! If l .lt. MAXL, store HES(l+1,l) and normalize the vector v(*,l+1).
              Hes(ll+1,ll) = snormw
              tem = 1.0D0/snormw
              CALL DSCAL(N,tem,V(1,ll+1),1)
           ENDDO
  !-----------------------------------------------------------------------
  ! l has reached MAXL without passing the convergence test:
  ! If RHO is not too large, compute a solution anyway and return with
  ! IFLAG = 1.  Otherwise return with IFLAG = 2.
  !-----------------------------------------------------------------------
   50      IF ( rho.LE.1.0D0 ) THEN
              Iflag = 1
              GOTO 200
           ELSEIF ( rho.LE.bnrm .AND. Mnewt.EQ.0 ) THEN
              Iflag = 1
              GOTO 200
           ENDIF
        ELSEIF ( Mnewt.GT.0 ) THEN
           DO i = 1 , N
              X(i) = 0.0D0
           ENDDO
           RETURN
        ELSE
           CALL DCOPY(N,B,1,X,1)
           RETURN
        ENDIF
   100  Iflag = 2
        RETURN
  !-----------------------------------------------------------------------
  ! Compute the approximation x(l) to the solution.
  ! Since the vector X was used as work space, and the initial guess
  ! of the Newton correction is zero, X must be reset to zero.
  !-----------------------------------------------------------------------
   200  ll = Liom
        DO k = 1 , ll
           B(k) = 0.0D0
        ENDDO
        B(1) = bnrm
        CALL DHESL(Hes,Maxl,ll,Ipvt,B)
        DO k = 1 , N
           X(k) = 0.0D0
        ENDDO
        DO i = 1 , ll
           CALL DAXPY(N,B(i),V(1,i),1,X,1)
        ENDDO
        DO i = 1 , N
           X(i) = X(i)/Wght(i)
        ENDDO
        IF ( Jpre.LE.1 ) RETURN
        CALL PSOL(Neq,Tn,Y,Savf,Wk,Hl0,Wp,Iwp,X,2,ier)
        Npsl = Npsl + 1
        IF ( ier.EQ.0 ) RETURN
  !-----------------------------------------------------------------------
  ! This block handles error returns forced by routine PSOL.
  !-----------------------------------------------------------------------
   300  IF ( ier.LT.0 ) Iflag = -1
        IF ( ier.GT.0 ) Iflag = 3
  !----------------------- End of Subroutine DSPIOM ----------------------
        END
  !*==DATV.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DATV
        SUBROUTINE DATV(Neq,Y,Savf,V,Wght,Ftem,F,PSOL,Z,Vtem,Wp,Iwp,Hl0,  &
       &                Jpre,Ier,Npsl)
        IMPLICIT NONE
  !*--DATV6733
        EXTERNAL F , PSOL
        INTEGER Neq , Iwp , Jpre , Ier , Npsl
        DOUBLE PRECISION Y , Savf , V , Wght , Ftem , Z , Vtem , Wp , Hl0
        DIMENSION Neq(*) , Y(*) , Savf(*) , V(*) , Wght(*) , Ftem(*) ,    &
       &          Z(*) , Vtem(*) , Wp(*) , Iwp(*)
        INTEGER IOWnd , IOWns , ICF , IERpj , IERsl , JCUr , JSTart ,     &
       &        KFLag , L , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,     &
       &        METh , MITer , MAXord , MAXcor , MSBp , MXNcf , N , NQ ,  &
       &        NST , NFE , NJE , NQU
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , IOWnd(6) , IOWns(6) , ICF ,    &
       &                IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , &
       &                LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,&
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
  !-----------------------------------------------------------------------
  ! This routine computes the product
  !
  !   (D-inverse)*(P1-inverse)*(I - hl0*df/dy)*(P2-inverse)*(D*v),
  !
  ! where D is a diagonal scaling matrix, and P1 and P2 are the
  ! left and right preconditioning matrices, respectively.
  ! v is assumed to have WRMS norm equal to 1.
  ! The product is stored in z.  This is computed by a
  ! difference quotient, a call to F, and two calls to PSOL.
  !-----------------------------------------------------------------------
  !
  !      On entry
  !
  !          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
  !
  !            Y = array containing current dependent variable vector.
  !
  !         SAVF = array containing current value of f(t,y).
  !
  !            V = real array of length N (can be the same array as Z).
  !
  !         WGHT = array of length N containing scale factors.
  !                1/WGHT(i) are the diagonal elements of the matrix D.
  !
  !         FTEM = work array of length N.
  !
  !         VTEM = work array of length N used to store the
  !                unscaled version of V.
  !
  !           WP = real work array used by preconditioner PSOL.
  !
  !          IWP = integer work array used by preconditioner PSOL.
  !
  !          HL0 = current value of (step size h) * (coefficient l0).
  !
  !         JPRE = preconditioner type flag.
  !
  !
  !      On return
  !
  !            Z = array of length N containing desired scaled
  !                matrix-vector product.
  !
  !          IER = error flag from PSOL.
  !
  !         NPSL = the number of calls to PSOL.
  !
  ! In addition, this routine uses the Common variables TN, N, NFE.
  !-----------------------------------------------------------------------
        INTEGER i
        DOUBLE PRECISION fac , rnorm , tempn
  !
  ! Set VTEM = D * V.
        DO i = 1 , N
           Vtem(i) = V(i)/Wght(i)
        ENDDO
        Ier = 0
        IF ( Jpre.GE.2 ) THEN
  !
  ! JPRE = 2 or 3.  Apply inverse of right preconditioner to VTEM.
           CALL PSOL(Neq,TN,Y,Savf,Ftem,Hl0,Wp,Iwp,Vtem,2,Ier)
           Npsl = Npsl + 1
           IF ( Ier.NE.0 ) RETURN
  ! Calculate L-2 norm of (D-inverse) * VTEM.
           DO i = 1 , N
              Z(i) = Vtem(i)*Wght(i)
           ENDDO
           tempn = DNRM2(N,Z,1)
           rnorm = 1.0D0/tempn
  ! Save Y in Z and increment Y by VTEM/norm.
           CALL DCOPY(N,Y,1,Z,1)
           DO i = 1 , N
              Y(i) = Z(i) + Vtem(i)*rnorm
           ENDDO
           fac = Hl0*tempn
        ELSE
  !
  ! JPRE = 0 or 1.  Save Y in Z and increment Y by VTEM.
           CALL DCOPY(N,Y,1,Z,1)
           DO i = 1 , N
              Y(i) = Z(i) + Vtem(i)
           ENDDO
           fac = Hl0
        ENDIF
  !
  ! For all JPRE, call F with incremented Y argument, and restore Y.
        CALL F(Neq,TN,Y,Ftem)
        NFE = NFE + 1
        CALL DCOPY(N,Z,1,Y,1)
  ! Set Z = (identity - hl0*Jacobian) * VTEM, using difference quotient.
        DO i = 1 , N
           Z(i) = Ftem(i) - Savf(i)
        ENDDO
        DO i = 1 , N
           Z(i) = Vtem(i) - fac*Z(i)
        ENDDO
  ! Apply inverse of left preconditioner to Z, if nontrivial.
        IF ( Jpre.NE.0 .AND. Jpre.NE.2 ) THEN
           CALL PSOL(Neq,TN,Y,Savf,Ftem,Hl0,Wp,Iwp,Z,1,Ier)
           Npsl = Npsl + 1
           IF ( Ier.NE.0 ) RETURN
        ENDIF
  ! Apply D-inverse to Z and return.
        DO i = 1 , N
           Z(i) = Z(i)*Wght(i)
        ENDDO
  !----------------------- End of Subroutine DATV ------------------------
        END
  !*==DORTHOG.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DORTHOG
        SUBROUTINE DORTHOG(Vnew,V,Hes,N,Ll,Ldhes,Kmp,Snormw)
        IMPLICIT NONE
  !*--DORTHOG6864
        INTEGER N , Ll , Ldhes , Kmp
        DOUBLE PRECISION Vnew , V , Hes , Snormw
        DIMENSION Vnew(*) , V(N,*) , Hes(Ldhes,*)
  !-----------------------------------------------------------------------
  ! This routine orthogonalizes the vector VNEW against the previous
  ! KMP vectors in the V array.  It uses a modified Gram-Schmidt
  ! orthogonalization procedure with conditional reorthogonalization.
  ! This is the version of 28 may 1986.
  !-----------------------------------------------------------------------
  !
  !      On entry
  !
  !         VNEW = the vector of length N containing a scaled product
  !                of the Jacobian and the vector V(*,LL).
  !
  !         V    = the N x l array containing the previous LL
  !                orthogonal vectors v(*,1) to v(*,LL).
  !
  !         HES  = an LL x LL upper Hessenberg matrix containing,
  !                in HES(i,k), k.lt.LL, scaled inner products of
  !                A*V(*,k) and V(*,i).
  !
  !        LDHES = the leading dimension of the HES array.
  !
  !         N    = the order of the matrix A, and the length of VNEW.
  !
  !         LL   = the current order of the matrix HES.
  !
  !          KMP = the number of previous vectors the new vector VNEW
  !                must be made orthogonal to (KMP .le. MAXL).
  !
  !
  !      On return
  !
  !         VNEW = the new vector orthogonal to V(*,i0) to V(*,LL),
  !                where i0 = MAX(1, LL-KMP+1).
  !
  !         HES  = upper Hessenberg matrix with column LL filled in with
  !                scaled inner products of A*V(*,LL) and V(*,i).
  !
  !       SNORMW = L-2 norm of VNEW.
  !
  !-----------------------------------------------------------------------
        INTEGER i , i0
        DOUBLE PRECISION arg , sumdsq , tem , vnrm
  !
  ! Get norm of unaltered VNEW for later use. ----------------------------
        vnrm = DNRM2(N,Vnew,1)
  !-----------------------------------------------------------------------
  ! Do modified Gram-Schmidt on VNEW = A*v(LL).
  ! Scaled inner products give new column of HES.
  ! Projections of earlier vectors are subtracted from VNEW.
  !-----------------------------------------------------------------------
        i0 = MAX(1,Ll-Kmp+1)
        DO i = i0 , Ll
           Hes(i,Ll) = DDOT(N,V(1,i),1,Vnew,1)
           tem = -Hes(i,Ll)
           CALL DAXPY(N,tem,V(1,i),1,Vnew,1)
        ENDDO
  !-----------------------------------------------------------------------
  ! Compute SNORMW = norm of VNEW.
  ! If VNEW is small compared to its input value (in norm), then
  ! reorthogonalize VNEW to V(*,1) through V(*,LL).
  ! Correct if relative correction exceeds 1000*(unit roundoff).
  ! finally, correct SNORMW using the dot products involved.
  !-----------------------------------------------------------------------
        Snormw = DNRM2(N,Vnew,1)
        IF ( vnrm+0.001D0*Snormw.NE.vnrm ) RETURN
        sumdsq = 0.0D0
        DO i = i0 , Ll
           tem = -DDOT(N,V(1,i),1,Vnew,1)
           IF ( Hes(i,Ll)+0.001D0*tem.NE.Hes(i,Ll) ) THEN
              Hes(i,Ll) = Hes(i,Ll) - tem
              CALL DAXPY(N,tem,V(1,i),1,Vnew,1)
              sumdsq = sumdsq + tem**2
           ENDIF
        ENDDO
        IF ( sumdsq.EQ.0.0D0 ) RETURN
        arg = MAX(0.0D0,Snormw**2-sumdsq)
        Snormw = SQRT(arg)
  !
  !----------------------- End of Subroutine DORTHOG ---------------------
        END
  !*==DSPIGMR.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DSPIGMR
        SUBROUTINE DSPIGMR(Neq,Tn,Y,Savf,B,Wght,N,Maxl,Maxlp1,Kmp,Delta,  &
       &                   Hl0,Jpre,Mnewt,F,PSOL,Npsl,X,V,Hes,Q,Lgmr,Wp,  &
       &                   Iwp,Wk,Dl,Iflag)
        IMPLICIT NONE
  !*--DSPIGMR6954
  !*** Start of declarations inserted by SPAG
        !REAL F
  !*** End of declarations inserted by SPAG
        EXTERNAL F , PSOL
        INTEGER Neq , N , Maxl , Maxlp1 , Kmp , Jpre , Mnewt , Npsl ,     &
       &        Lgmr , Iwp , Iflag
        DOUBLE PRECISION Tn , Y , Savf , B , Wght , Delta , Hl0 , X , V , &
       &                 Hes , Q , Wp , Wk , Dl
        DIMENSION Neq(*) , Y(*) , Savf(*) , B(*) , Wght(*) , X(*) , V(N,*)&
       &          , Hes(Maxlp1,*) , Q(*) , Wp(*) , Iwp(*) , Wk(*) , Dl(*)
  !-----------------------------------------------------------------------
  ! This routine solves the linear system A * x = b using a scaled
  ! preconditioned version of the Generalized Minimal Residual method.
  ! An initial guess of x = 0 is assumed.
  !-----------------------------------------------------------------------
  !
  !      On entry
  !
  !          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
  !
  !           TN = current value of t.
  !
  !            Y = array containing current dependent variable vector.
  !
  !         SAVF = array containing current value of f(t,y).
  !
  !            B = the right hand side of the system A*x = b.
  !                B is also used as work space when computing
  !                the final approximation.
  !                (B is the same as V(*,MAXL+1) in the call to DSPIGMR.)
  !
  !         WGHT = the vector of length N containing the nonzero
  !                elements of the diagonal scaling matrix.
  !
  !            N = the order of the matrix A, and the lengths
  !                of the vectors WGHT, B and X.
  !
  !         MAXL = the maximum allowable order of the matrix HES.
  !
  !       MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES.
  !
  !          KMP = the number of previous vectors the new vector VNEW
  !                must be made orthogonal to.  KMP .le. MAXL.
  !
  !        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
  !
  !          HL0 = current value of (step size h) * (coefficient l0).
  !
  !         JPRE = preconditioner type flag.
  !
  !        MNEWT = Newton iteration counter (.ge. 0).
  !
  !           WK = real work array used by routine DATV and PSOL.
  !
  !           DL = real work array used for calculation of the residual
  !                norm RHO when the method is incomplete (KMP .lt. MAXL).
  !                Not needed or referenced in complete case (KMP = MAXL).
  !
  !           WP = real work array used by preconditioner PSOL.
  !
  !          IWP = integer work array used by preconditioner PSOL.
  !
  !      On return
  !
  !         X    = the final computed approximation to the solution
  !                of the system A*x = b.
  !
  !         LGMR = the number of iterations performed and
  !                the current order of the upper Hessenberg
  !                matrix HES.
  !
  !         NPSL = the number of calls to PSOL.
  !
  !         V    = the N by (LGMR+1) array containing the LGMR
  !                orthogonal vectors V(*,1) to V(*,LGMR).
  !
  !         HES  = the upper triangular factor of the QR decomposition
  !                of the (LGMR+1) by lgmr upper Hessenberg matrix whose
  !                entries are the scaled inner-products of A*V(*,i)
  !                and V(*,k).
  !
  !         Q    = real array of length 2*MAXL containing the components
  !                of the Givens rotations used in the QR decomposition
  !                of HES.  It is loaded in DHEQR and used in DHELS.
  !
  !        IFLAG = integer error flag:
  !                0 means convergence in LGMR iterations, LGMR .le. MAXL.
  !                1 means the convergence test did not pass in MAXL
  !                  iterations, but the residual norm is .lt. 1,
  !                  or .lt. norm(b) if MNEWT = 0, and so x is computed.
  !                2 means the convergence test did not pass in MAXL
  !                  iterations, residual .gt. 1, and X is undefined.
  !                3 means there was a recoverable error in PSOL
  !                  caused by the preconditioner being out of date.
  !               -1 means there was a nonrecoverable error in PSOL.
  !
  !-----------------------------------------------------------------------
        INTEGER i , ier , info , ip1 , i2 , j , k , ll , llp1
        DOUBLE PRECISION bnrm , bnrm0 , c , dlnrm , prod , rho , s ,      &
       &                 snormw , tem
  !
        Iflag = 0
        Lgmr = 0
        Npsl = 0
  !-----------------------------------------------------------------------
  ! The initial residual is the vector b.  Apply scaling to b, and test
  ! for an immediate return with X = 0 or X = b.
  !-----------------------------------------------------------------------
        DO i = 1 , N
           V(i,1) = B(i)*Wght(i)
        ENDDO
        bnrm0 = DNRM2(N,V,1)
        bnrm = bnrm0
        IF ( bnrm0.GT.Delta ) THEN
  ! Apply inverse of left preconditioner to vector b. --------------------
           ier = 0
           IF ( Jpre.NE.0 .AND. Jpre.NE.2 ) THEN
              CALL PSOL(Neq,Tn,Y,Savf,Wk,Hl0,Wp,Iwp,B,1,ier)
              Npsl = 1
              IF ( ier.NE.0 ) GOTO 300
  ! Calculate norm of scaled vector V(*,1) and normalize it. -------------
              DO i = 1 , N
                 V(i,1) = B(i)*Wght(i)
              ENDDO
              bnrm = DNRM2(N,V,1)
              Delta = Delta*(bnrm/bnrm0)
           ENDIF
           tem = 1.0D0/bnrm
           CALL DSCAL(N,tem,V(1,1),1)
  ! Zero out the HES array. ----------------------------------------------
           DO j = 1 , Maxl
              DO i = 1 , Maxlp1
                 Hes(i,j) = 0.0D0
              ENDDO
           ENDDO
  !-----------------------------------------------------------------------
  ! Main loop to compute the vectors V(*,2) to V(*,MAXL).
  ! The running product PROD is needed for the convergence test.
  !-----------------------------------------------------------------------
           prod = 1.0D0
           DO ll = 1 , Maxl
              Lgmr = ll
  !-----------------------------------------------------------------------
  ! Call routine DATV to compute VNEW = Abar*v(ll), where Abar is
  ! the matrix A with scaling and inverse preconditioner factors applied.
  ! Call routine DORTHOG to orthogonalize the new vector VNEW = V(*,LL+1).
  ! Call routine DHEQR to update the factors of HES.
  !-----------------------------------------------------------------------
              CALL DATV(Neq,Y,Savf,V(1,ll),Wght,X,F,PSOL,V(1,ll+1),Wk,Wp, &
       &                Iwp,Hl0,Jpre,ier,Npsl)
              IF ( ier.NE.0 ) GOTO 300
              CALL DORTHOG(V(1,ll+1),V,Hes,N,ll,Maxlp1,Kmp,snormw)
              Hes(ll+1,ll) = snormw
              CALL DHEQR(Hes,Maxlp1,ll,Q,info,ll)
              IF ( info.EQ.ll ) GOTO 100
  !-----------------------------------------------------------------------
  ! Update RHO, the estimate of the norm of the residual b - A*xl.
  ! If KMP .lt. MAXL, then the vectors V(*,1),...,V(*,LL+1) are not
  ! necessarily orthogonal for LL .gt. KMP.  The vector DL must then
  ! be computed, and its norm used in the calculation of RHO.
  !-----------------------------------------------------------------------
              prod = prod*Q(2*ll)
              rho = ABS(prod*bnrm)
              IF ( (ll.GT.Kmp) .AND. (Kmp.LT.Maxl) ) THEN
                 IF ( ll.EQ.Kmp+1 ) THEN
                    CALL DCOPY(N,V(1,1),1,Dl,1)
                    DO i = 1 , Kmp
                       ip1 = i + 1
                       i2 = i*2
                       s = Q(i2)
                       c = Q(i2-1)
                       DO k = 1 , N
                          Dl(k) = s*Dl(k) + c*V(k,ip1)
                       ENDDO
                    ENDDO
                 ENDIF
                 s = Q(2*ll)
                 c = Q(2*ll-1)/snormw
                 llp1 = ll + 1
                 DO k = 1 , N
                    Dl(k) = s*Dl(k) + c*V(k,llp1)
                 ENDDO
                 dlnrm = DNRM2(N,Dl,1)
                 rho = rho*dlnrm
              ENDIF
  !-----------------------------------------------------------------------
  ! Test for convergence.  If passed, compute approximation xl.
  ! if failed and LL .lt. MAXL, then continue iterating.
  !-----------------------------------------------------------------------
              IF ( rho.LE.Delta ) GOTO 200
              IF ( ll.EQ.Maxl ) GOTO 50
  !-----------------------------------------------------------------------
  ! Rescale so that the norm of V(1,LL+1) is one.
  !-----------------------------------------------------------------------
              tem = 1.0D0/snormw
              CALL DSCAL(N,tem,V(1,ll+1),1)
           ENDDO
   50      IF ( rho.LE.1.0D0 ) THEN
              Iflag = 1
              GOTO 200
           ELSEIF ( rho.LE.bnrm .AND. Mnewt.EQ.0 ) THEN
              Iflag = 1
              GOTO 200
           ENDIF
        ELSEIF ( Mnewt.GT.0 ) THEN
           DO i = 1 , N
              X(i) = 0.0D0
           ENDDO
           RETURN
        ELSE
           CALL DCOPY(N,B,1,X,1)
           RETURN
        ENDIF
   100  Iflag = 2
        RETURN
  !-----------------------------------------------------------------------
  ! Compute the approximation xl to the solution.
  ! Since the vector X was used as work space, and the initial guess
  ! of the Newton correction is zero, X must be reset to zero.
  !-----------------------------------------------------------------------
   200  ll = Lgmr
        llp1 = ll + 1
        DO k = 1 , llp1
           B(k) = 0.0D0
        ENDDO
        B(1) = bnrm
        CALL DHELS(Hes,Maxlp1,ll,Q,B)
        DO k = 1 , N
           X(k) = 0.0D0
        ENDDO
        DO i = 1 , ll
           CALL DAXPY(N,B(i),V(1,i),1,X,1)
        ENDDO
        DO i = 1 , N
           X(i) = X(i)/Wght(i)
        ENDDO
        IF ( Jpre.LE.1 ) RETURN
        CALL PSOL(Neq,Tn,Y,Savf,Wk,Hl0,Wp,Iwp,X,2,ier)
        Npsl = Npsl + 1
        IF ( ier.EQ.0 ) RETURN
  !-----------------------------------------------------------------------
  ! This block handles error returns forced by routine PSOL.
  !-----------------------------------------------------------------------
   300  IF ( ier.LT.0 ) Iflag = -1
        IF ( ier.GT.0 ) Iflag = 3
  !
  !----------------------- End of Subroutine DSPIGMR ---------------------
        END
  !*==DPCG.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DPCG
        SUBROUTINE DPCG(Neq,Tn,Y,Savf,R,Wght,N,Maxl,Delta,Hl0,Jpre,Mnewt, &
       &                F,PSOL,Npsl,X,P,W,Z,Lpcg,Wp,Iwp,Wk,Iflag)
        IMPLICIT NONE
  !*--DPCG7208
  !*** Start of declarations inserted by SPAG
        !REAL F
  !*** End of declarations inserted by SPAG
        EXTERNAL F , PSOL
        INTEGER Neq , N , Maxl , Jpre , Mnewt , Npsl , Lpcg , Iwp , Iflag
        DOUBLE PRECISION Tn , Y , Savf , R , Wght , Delta , Hl0 , X , P , &
       &                 W , Z , Wp , Wk
        DIMENSION Neq(*) , Y(*) , Savf(*) , R(*) , Wght(*) , X(*) , P(*) ,&
       &          W(*) , Z(*) , Wp(*) , Iwp(*) , Wk(*)
  !-----------------------------------------------------------------------
  ! This routine computes the solution to the system A*x = b using a
  ! preconditioned version of the Conjugate Gradient algorithm.
  ! It is assumed here that the matrix A and the preconditioner
  ! matrix M are symmetric positive definite or nearly so.
  !-----------------------------------------------------------------------
  !
  !      On entry
  !
  !          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
  !
  !           TN = current value of t.
  !
  !            Y = array containing current dependent variable vector.
  !
  !         SAVF = array containing current value of f(t,y).
  !
  !            R = the right hand side of the system A*x = b.
  !
  !         WGHT = array of length N containing scale factors.
  !                1/WGHT(i) are the diagonal elements of the diagonal
  !                scaling matrix D.
  !
  !            N = the order of the matrix A, and the lengths
  !                of the vectors Y, SAVF, R, WGHT, P, W, Z, WK, and X.
  !
  !         MAXL = the maximum allowable number of iterates.
  !
  !        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
  !
  !          HL0 = current value of (step size h) * (coefficient l0).
  !
  !         JPRE = preconditioner type flag.
  !
  !        MNEWT = Newton iteration counter (.ge. 0).
  !
  !           WK = real work array used by routine DATP.
  !
  !           WP = real work array used by preconditioner PSOL.
  !
  !          IWP = integer work array used by preconditioner PSOL.
  !
  !      On return
  !
  !         X    = the final computed approximation to the solution
  !                of the system A*x = b.
  !
  !         LPCG = the number of iterations performed, and current
  !                order of the upper Hessenberg matrix HES.
  !
  !         NPSL = the number of calls to PSOL.
  !
  !        IFLAG = integer error flag:
  !                0 means convergence in LPCG iterations, LPCG .le. MAXL.
  !                1 means the convergence test did not pass in MAXL
  !                  iterations, but the residual norm is .lt. 1,
  !                  or .lt. norm(b) if MNEWT = 0, and so X is computed.
  !                2 means the convergence test did not pass in MAXL
  !                  iterations, residual .gt. 1, and X is undefined.
  !                3 means there was a recoverable error in PSOL
  !                  caused by the preconditioner being out of date.
  !                4 means there was a zero denominator in the algorithm.
  !                  The system matrix or preconditioner matrix is not
  !                  sufficiently close to being symmetric pos. definite.
  !               -1 means there was a nonrecoverable error in PSOL.
  !
  !-----------------------------------------------------------------------
        INTEGER i , ier
        DOUBLE PRECISION alpha , beta , bnrm , ptw , rnrm , ztr , ztr0
  !
        Iflag = 0
        Npsl = 0
        Lpcg = 0
        DO i = 1 , N
           X(i) = 0.0D0
        ENDDO
        bnrm = DVNORM(N,R,Wght)
  ! Test for immediate return with X = 0 or X = b. -----------------------
        IF ( bnrm.GT.Delta ) THEN
  !
           ztr = 0.0D0
        ELSE
           IF ( Mnewt.GT.0 ) RETURN
           CALL DCOPY(N,R,1,X,1)
           RETURN
        ENDIF
  ! Loop point for PCG iterations. ---------------------------------------
   100  Lpcg = Lpcg + 1
        CALL DCOPY(N,R,1,Z,1)
        ier = 0
        IF ( Jpre.NE.0 ) THEN
           CALL PSOL(Neq,Tn,Y,Savf,Wk,Hl0,Wp,Iwp,Z,3,ier)
           Npsl = Npsl + 1
           IF ( ier.NE.0 ) THEN
  !-----------------------------------------------------------------------
  ! This block handles error returns from PSOL.
  !-----------------------------------------------------------------------
              IF ( ier.LT.0 ) Iflag = -1
              IF ( ier.GT.0 ) Iflag = 3
              RETURN
           ENDIF
        ENDIF
        ztr0 = ztr
        ztr = DDOT(N,Z,1,R,1)
        IF ( Lpcg.EQ.1 ) THEN
           CALL DCOPY(N,Z,1,P,1)
        ELSEIF ( ztr0.EQ.0.0D0 ) THEN
  !-----------------------------------------------------------------------
  ! This block handles division by zero errors.
  !-----------------------------------------------------------------------
           Iflag = 4
           GOTO 99999
        ELSE
           beta = ztr/ztr0
           DO i = 1 , N
              P(i) = Z(i) + beta*P(i)
           ENDDO
        ENDIF
  !-----------------------------------------------------------------------
  !  Call DATP to compute A*p and return the answer in W.
  !-----------------------------------------------------------------------
        CALL DATP(Neq,Y,Savf,P,Wght,Hl0,Wk,F,W)
  !
        ptw = DDOT(N,P,1,W,1)
        IF ( ptw.EQ.0.0D0 ) THEN
           Iflag = 4
        ELSE
           alpha = ztr/ptw
           CALL DAXPY(N,alpha,P,1,X,1)
           alpha = -alpha
           CALL DAXPY(N,alpha,W,1,R,1)
           rnrm = DVNORM(N,R,Wght)
           IF ( rnrm.LE.Delta ) RETURN
           IF ( Lpcg.LT.Maxl ) GOTO 100
           Iflag = 2
           IF ( rnrm.LE.1.0D0 ) Iflag = 1
           IF ( rnrm.LE.bnrm .AND. Mnewt.EQ.0 ) Iflag = 1
           RETURN
        ENDIF
  !----------------------- End of Subroutine DPCG ------------------------
  99999 END
  !*==DPCGS.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DPCGS
        SUBROUTINE DPCGS(Neq,Tn,Y,Savf,R,Wght,N,Maxl,Delta,Hl0,Jpre,Mnewt,&
       &                 F,PSOL,Npsl,X,P,W,Z,Lpcg,Wp,Iwp,Wk,Iflag)
        IMPLICIT NONE
  !*--DPCGS7365
  !*** Start of declarations inserted by SPAG
        !REAL F
  !*** End of declarations inserted by SPAG
        EXTERNAL F , PSOL
        INTEGER Neq , N , Maxl , Jpre , Mnewt , Npsl , Lpcg , Iwp , Iflag
        DOUBLE PRECISION Tn , Y , Savf , R , Wght , Delta , Hl0 , X , P , &
       &                 W , Z , Wp , Wk
        DIMENSION Neq(*) , Y(*) , Savf(*) , R(*) , Wght(*) , X(*) , P(*) ,&
       &          W(*) , Z(*) , Wp(*) , Iwp(*) , Wk(*)
  !-----------------------------------------------------------------------
  ! This routine computes the solution to the system A*x = b using a
  ! scaled preconditioned version of the Conjugate Gradient algorithm.
  ! It is assumed here that the scaled matrix D**-1 * A * D and the
  ! scaled preconditioner D**-1 * M * D are close to being
  ! symmetric positive definite.
  !-----------------------------------------------------------------------
  !
  !      On entry
  !
  !          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
  !
  !           TN = current value of t.
  !
  !            Y = array containing current dependent variable vector.
  !
  !         SAVF = array containing current value of f(t,y).
  !
  !            R = the right hand side of the system A*x = b.
  !
  !         WGHT = array of length N containing scale factors.
  !                1/WGHT(i) are the diagonal elements of the diagonal
  !                scaling matrix D.
  !
  !            N = the order of the matrix A, and the lengths
  !                of the vectors Y, SAVF, R, WGHT, P, W, Z, WK, and X.
  !
  !         MAXL = the maximum allowable number of iterates.
  !
  !        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
  !
  !          HL0 = current value of (step size h) * (coefficient l0).
  !
  !         JPRE = preconditioner type flag.
  !
  !        MNEWT = Newton iteration counter (.ge. 0).
  !
  !           WK = real work array used by routine DATP.
  !
  !           WP = real work array used by preconditioner PSOL.
  !
  !          IWP = integer work array used by preconditioner PSOL.
  !
  !      On return
  !
  !         X    = the final computed approximation to the solution
  !                of the system A*x = b.
  !
  !         LPCG = the number of iterations performed, and current
  !                order of the upper Hessenberg matrix HES.
  !
  !         NPSL = the number of calls to PSOL.
  !
  !        IFLAG = integer error flag:
  !                0 means convergence in LPCG iterations, LPCG .le. MAXL.
  !                1 means the convergence test did not pass in MAXL
  !                  iterations, but the residual norm is .lt. 1,
  !                  or .lt. norm(b) if MNEWT = 0, and so X is computed.
  !                2 means the convergence test did not pass in MAXL
  !                  iterations, residual .gt. 1, and X is undefined.
  !                3 means there was a recoverable error in PSOL
  !                  caused by the preconditioner being out of date.
  !                4 means there was a zero denominator in the algorithm.
  !                  the scaled matrix or scaled preconditioner is not
  !                  sufficiently close to being symmetric pos. definite.
  !               -1 means there was a nonrecoverable error in PSOL.
  !
  !-----------------------------------------------------------------------
        INTEGER i , ier
        DOUBLE PRECISION alpha , beta , bnrm , ptw , rnrm , ztr, ztr0
  !
        Iflag = 0
        Npsl = 0
        Lpcg = 0
        DO i = 1 , N
           X(i) = 0.0D0
        ENDDO
        bnrm = DVNORM(N,R,Wght)
  ! Test for immediate return with X = 0 or X = b. -----------------------
        IF ( bnrm.GT.Delta ) THEN
  !
           ztr = 0.0D0
        ELSE
           IF ( Mnewt.GT.0 ) RETURN
           CALL DCOPY(N,R,1,X,1)
           RETURN
        ENDIF
  ! Loop point for PCG iterations. ---------------------------------------
   100  Lpcg = Lpcg + 1
        CALL DCOPY(N,R,1,Z,1)
        ier = 0
        IF ( Jpre.NE.0 ) THEN
           CALL PSOL(Neq,Tn,Y,Savf,Wk,Hl0,Wp,Iwp,Z,3,ier)
           Npsl = Npsl + 1
           IF ( ier.NE.0 ) THEN
  !-----------------------------------------------------------------------
  ! This block handles error returns from PSOL.
  !-----------------------------------------------------------------------
              IF ( ier.LT.0 ) Iflag = -1
              IF ( ier.GT.0 ) Iflag = 3
              RETURN
           ENDIF
        ENDIF
        ztr0 = ztr
        ztr = 0.0D0
        DO i = 1 , N
           ztr = ztr + Z(i)*R(i)*Wght(i)**2
        ENDDO
        IF ( Lpcg.EQ.1 ) THEN
           CALL DCOPY(N,Z,1,P,1)
        ELSEIF ( ztr0.EQ.0.0D0 ) THEN
  !-----------------------------------------------------------------------
  ! This block handles division by zero errors.
  !-----------------------------------------------------------------------
           Iflag = 4
           GOTO 99999
        ELSE
           beta = ztr/ztr0
           DO i = 1 , N
              P(i) = Z(i) + beta*P(i)
           ENDDO
        ENDIF
  !-----------------------------------------------------------------------
  !  Call DATP to compute A*p and return the answer in W.
  !-----------------------------------------------------------------------
        CALL DATP(Neq,Y,Savf,P,Wght,Hl0,Wk,F,W)
  !
        ptw = 0.0D0
        DO i = 1 , N
           ptw = ptw + P(i)*W(i)*Wght(i)**2
        ENDDO
        IF ( ptw.EQ.0.0D0 ) THEN
           Iflag = 4
        ELSE
           alpha = ztr/ptw
           CALL DAXPY(N,alpha,P,1,X,1)
           alpha = -alpha
           CALL DAXPY(N,alpha,W,1,R,1)
           rnrm = DVNORM(N,R,Wght)
           IF ( rnrm.LE.Delta ) RETURN
           IF ( Lpcg.LT.Maxl ) GOTO 100
           Iflag = 2
           IF ( rnrm.LE.1.0D0 ) Iflag = 1
           IF ( rnrm.LE.bnrm .AND. Mnewt.EQ.0 ) Iflag = 1
           RETURN
        ENDIF
  !----------------------- End of Subroutine DPCGS -----------------------
  99999 END
  !*==DATP.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DATP
        SUBROUTINE DATP(Neq,Y,Savf,P,Wght,Hl0,Wk,F,W)
        IMPLICIT NONE
  !*--DATP7528
        EXTERNAL F
        INTEGER Neq
        DOUBLE PRECISION Y , Savf , P , Wght , Hl0 , Wk , W
        DIMENSION Neq(*) , Y(*) , Savf(*) , P(*) , Wght(*) , Wk(*) , W(*)
        INTEGER IOWnd , IOWns , ICF , IERpj , IERsl , JCUr , JSTart ,     &
       &        KFLag , L , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,     &
       &        METh , MITer , MAXord , MAXcor , MSBp , MXNcf , N , NQ ,  &
       &        NST , NFE , NJE , NQU
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , IOWnd(6) , IOWns(6) , ICF ,    &
       &                IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , &
       &                LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,&
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
  !-----------------------------------------------------------------------
  ! This routine computes the product
  !
  !              w = (I - hl0*df/dy)*p
  !
  ! This is computed by a call to F and a difference quotient.
  !-----------------------------------------------------------------------
  !
  !      On entry
  !
  !          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
  !
  !            Y = array containing current dependent variable vector.
  !
  !         SAVF = array containing current value of f(t,y).
  !
  !            P = real array of length N.
  !
  !         WGHT = array of length N containing scale factors.
  !                1/WGHT(i) are the diagonal elements of the matrix D.
  !
  !           WK = work array of length N.
  !
  !      On return
  !
  !
  !            W = array of length N containing desired
  !                matrix-vector product.
  !
  ! In addition, this routine uses the Common variables TN, N, NFE.
  !-----------------------------------------------------------------------
        INTEGER i
        DOUBLE PRECISION fac , pnrm , rpnrm
  !
        pnrm = DVNORM(N,P,Wght)
        rpnrm = 1.0D0/pnrm
        CALL DCOPY(N,Y,1,W,1)
        DO i = 1 , N
           Y(i) = W(i) + P(i)*rpnrm
        ENDDO
        CALL F(Neq,TN,Y,Wk)
        NFE = NFE + 1
        CALL DCOPY(N,W,1,Y,1)
        fac = Hl0*pnrm
        DO i = 1 , N
           W(i) = P(i) - fac*(Wk(i)-Savf(i))
        ENDDO
  !----------------------- End of Subroutine DATP ------------------------
        END
  !*==DUSOL.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DUSOL
        SUBROUTINE DUSOL(Neq,Tn,Y,Savf,B,Wght,N,Delta,Hl0,Mnewt,PSOL,Npsl,&
       &                 X,Wp,Iwp,Wk,Iflag)
        IMPLICIT NONE
  !*--DUSOL7599
        EXTERNAL PSOL
        INTEGER Neq , N , Mnewt , Npsl , Iwp , Iflag
        DOUBLE PRECISION Tn , Y , Savf , B , Wght , Delta , Hl0 , X , Wp ,&
       &                 Wk
        DIMENSION Neq(*) , Y(*) , Savf(*) , B(*) , Wght(*) , X(*) ,       &
       &          Wp(*) , Iwp(*) , Wk(*)
  !-----------------------------------------------------------------------
  ! This routine solves the linear system A * x = b using only a call
  ! to the user-supplied routine PSOL (no Krylov iteration).
  ! If the norm of the right-hand side vector b is smaller than DELTA,
  ! the vector X returned is X = b (if MNEWT = 0) or X = 0 otherwise.
  ! PSOL is called with an LR argument of 0.
  !-----------------------------------------------------------------------
  !
  !      On entry
  !
  !          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
  !
  !           TN = current value of t.
  !
  !            Y = array containing current dependent variable vector.
  !
  !         SAVF = array containing current value of f(t,y).
  !
  !            B = the right hand side of the system A*x = b.
  !
  !         WGHT = the vector of length N containing the nonzero
  !                elements of the diagonal scaling matrix.
  !
  !            N = the order of the matrix A, and the lengths
  !                of the vectors WGHT, B and X.
  !
  !        DELTA = tolerance on residuals b - A*x in weighted RMS-norm.
  !
  !          HL0 = current value of (step size h) * (coefficient l0).
  !
  !        MNEWT = Newton iteration counter (.ge. 0).
  !
  !           WK = real work array used by PSOL.
  !
  !           WP = real work array used by preconditioner PSOL.
  !
  !          IWP = integer work array used by preconditioner PSOL.
  !
  !      On return
  !
  !         X    = the final computed approximation to the solution
  !                of the system A*x = b.
  !
  !         NPSL = the number of calls to PSOL.
  !
  !        IFLAG = integer error flag:
  !                0 means no trouble occurred.
  !                3 means there was a recoverable error in PSOL
  !                  caused by the preconditioner being out of date.
  !               -1 means there was a nonrecoverable error in PSOL.
  !
  !-----------------------------------------------------------------------
        INTEGER i , ier
        DOUBLE PRECISION bnrm
  !
        Iflag = 0
        Npsl = 0
  !-----------------------------------------------------------------------
  ! Test for an immediate return with X = 0 or X = b.
  !-----------------------------------------------------------------------
        bnrm = DVNORM(N,B,Wght)
        IF ( bnrm.GT.Delta ) THEN
  ! Make call to PSOL and copy result from B to X. -----------------------
           ier = 0
           CALL PSOL(Neq,Tn,Y,Savf,Wk,Hl0,Wp,Iwp,B,0,ier)
           Npsl = 1
           IF ( ier.NE.0 ) THEN
  !-----------------------------------------------------------------------
  ! This block handles error returns forced by routine PSOL.
  !-----------------------------------------------------------------------
              IF ( ier.LT.0 ) Iflag = -1
              IF ( ier.GT.0 ) Iflag = 3
              GOTO 99999
           ENDIF
        ELSEIF ( Mnewt.GT.0 ) THEN
           DO i = 1 , N
              X(i) = 0.0D0
           ENDDO
           RETURN
        ELSE
           CALL DCOPY(N,B,1,X,1)
           RETURN
        ENDIF
        CALL DCOPY(N,B,1,X,1)
        RETURN
  !----------------------- End of Subroutine DUSOL -----------------------
  99999 END
  !*==DSRCPK.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DSRCPK
        SUBROUTINE DSRCPK(Rsav,Isav,Job)
        IMPLICIT NONE
  !*--DSRCPK7697
  !-----------------------------------------------------------------------
  ! This routine saves or restores (depending on JOB) the contents of
  ! the Common blocks DLS001, DLPK01, which are used
  ! internally by the DLSODPK solver.
  !
  ! RSAV = real array of length 222 or more.
  ! ISAV = integer array of length 50 or more.
  ! JOB  = flag indicating to save or restore the Common blocks:
  !        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
  !        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
  !        A call with JOB = 2 presumes a prior call with JOB = 1.
  !-----------------------------------------------------------------------
        INTEGER Isav , Job
        INTEGER ILS , ILSp
        INTEGER i , lenilp , lenrlp , lenils , lenrls
        DOUBLE PRECISION Rsav , RLS , RLSp
        DIMENSION Rsav(*) , Isav(*)
        SAVE lenrls , lenils , lenrlp , lenilp
        COMMON /DLS001/ RLS(218) , ILS(37)
        COMMON /DLPK01/ RLSp(4) , ILSp(13)
        DATA lenrls/218/ , lenils/37/ , lenrlp/4/ , lenilp/13/
  !
        IF ( Job.EQ.2 ) THEN
  !
           CALL DCOPY(lenrls,Rsav,1,RLS,1)
           CALL DCOPY(lenrlp,Rsav(lenrls+1),1,RLSp,1)
           DO i = 1 , lenils
              ILS(i) = Isav(i)
           ENDDO
           DO i = 1 , lenilp
              ILSp(i) = Isav(lenils+i)
           ENDDO
           GOTO 99999
        ENDIF
        CALL DCOPY(lenrls,RLS,1,Rsav,1)
        CALL DCOPY(lenrlp,RLSp,1,Rsav(lenrls+1),1)
        DO i = 1 , lenils
           Isav(i) = ILS(i)
        ENDDO
        DO i = 1 , lenilp
           Isav(lenils+i) = ILSp(i)
        ENDDO
        RETURN
  !----------------------- End of Subroutine DSRCPK ----------------------
  99999 END
  !*==DHEFA.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DHEFA
        SUBROUTINE DHEFA(A,Lda,N,Ipvt,Info,Job)
        IMPLICIT NONE
  !*--DHEFA7747
        INTEGER Lda , N , Ipvt(*) , Info , Job
        DOUBLE PRECISION A(Lda,*)
  !-----------------------------------------------------------------------
  !     This routine is a modification of the LINPACK routine DGEFA and
  !     performs an LU decomposition of an upper Hessenberg matrix A.
  !     There are two options available:
  !
  !          (1)  performing a fresh factorization
  !          (2)  updating the LU factors by adding a row and a
  !               column to the matrix A.
  !-----------------------------------------------------------------------
  !     DHEFA factors an upper Hessenberg matrix by elimination.
  !
  !     On entry
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
  !        JOB     INTEGER
  !                JOB = 1    means that a fresh factorization of the
  !                           matrix A is desired.
  !                JOB .ge. 2 means that the current factorization of A
  !                           will be updated by the addition of a row
  !                           and a column.
  !
  !     On return
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
  !                = k  if  U(k,k) .eq. 0.0 .  This is not an error
  !                     condition for this subroutine, but it does
  !                     indicate that DHESL will divide by zero if called.
  !
  !     Modification of LINPACK, by Peter Brown, LLNL.
  !     Written 7/20/83.  This version dated 6/20/01.
  !
  !     BLAS called: DAXPY, IDAMAX
  !-----------------------------------------------------------------------
        INTEGER j , k , km1 , kp1 , l , nm1
        DOUBLE PRECISION t
  !
        IF ( Job.GT.1 ) THEN
  !
  ! The old factorization of A will be updated.  A row and a column
  ! has been added to the matrix A.
  ! N-1 is now the old order of the matrix.
  !
           nm1 = N - 1
  !
  ! Perform row interchanges on the elements of the new column, and
  ! perform elimination operations on the elements using the multipliers.
  !
           IF ( nm1.GT.1 ) THEN
              DO k = 2 , nm1
                 km1 = k - 1
                 l = Ipvt(km1)
                 t = A(l,N)
                 IF ( l.NE.km1 ) THEN
                    A(l,N) = A(km1,N)
                    A(km1,N) = t
                 ENDIF
                 A(k,N) = A(k,N) + A(k,km1)*t
              ENDDO
           ENDIF
  !
  ! Complete update of factorization by decomposing last 2x2 block.
  !
           Info = 0
  !
  !        Find L = pivot index
  !
           l = IDAMAX(2,A(nm1,nm1),1) + nm1 - 1
           Ipvt(nm1) = l
  !
  !        Zero pivot implies this column already triangularized
  !
           IF ( A(l,nm1).EQ.0.0D0 ) THEN
              Info = nm1
           ELSE
  !
  !           Interchange if necessary
  !
              IF ( l.NE.nm1 ) THEN
                 t = A(l,nm1)
                 A(l,nm1) = A(nm1,nm1)
                 A(nm1,nm1) = t
              ENDIF
  !
  !           Compute multipliers
  !
              t = -1.0D0/A(nm1,nm1)
              A(N,nm1) = A(N,nm1)*t
  !
  !           Row elimination with column indexing
  !
              t = A(l,N)
              IF ( l.NE.nm1 ) THEN
                 A(l,N) = A(nm1,N)
                 A(nm1,N) = t
              ENDIF
              A(N,N) = A(N,N) + t*A(N,nm1)
           ENDIF
           Ipvt(N) = N
           IF ( A(N,N).EQ.0.0D0 ) Info = N
        ELSE
  !
  ! A new facorization is desired.  This is essentially the LINPACK
  ! code with the exception that we know there is only one nonzero
  ! element below the main diagonal.
  !
  !     Gaussian elimination with partial pivoting
  !
           Info = 0
           nm1 = N - 1
           IF ( nm1.GE.1 ) THEN
              DO k = 1 , nm1
                 kp1 = k + 1
  !
  !        Find L = pivot index
  !
                 l = IDAMAX(2,A(k,k),1) + k - 1
                 Ipvt(k) = l
  !
  !        Zero pivot implies this column already triangularized
  !
                 IF ( A(l,k).EQ.0.0D0 ) THEN
                    Info = k
                 ELSE
  !
  !           Interchange if necessary
  !
                    IF ( l.NE.k ) THEN
                       t = A(l,k)
                       A(l,k) = A(k,k)
                       A(k,k) = t
                    ENDIF
  !
  !           Compute multipliers
  !
                    t = -1.0D0/A(k,k)
                    A(k+1,k) = A(k+1,k)*t
  !
  !           Row elimination with column indexing
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
           RETURN
        ENDIF
  !----------------------- End of Subroutine DHEFA -----------------------
        END
  !*==DHESL.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DHESL
        SUBROUTINE DHESL(A,Lda,N,Ipvt,B)
        IMPLICIT NONE
  !*--DHESL7928
        INTEGER Lda , N , Ipvt(*)
        DOUBLE PRECISION A(Lda,*) , B(*)
  !-----------------------------------------------------------------------
  ! This is essentially the LINPACK routine DGESL except for changes
  ! due to the fact that A is an upper Hessenberg matrix.
  !-----------------------------------------------------------------------
  !     DHESL solves the real system A * x = b
  !     using the factors computed by DHEFA.
  !
  !     On entry
  !
  !        A       DOUBLE PRECISION(LDA, N)
  !                the output from DHEFA.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  A .
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !        IPVT    INTEGER(N)
  !                the pivot vector from DHEFA.
  !
  !        B       DOUBLE PRECISION(N)
  !                the right hand side vector.
  !
  !     On return
  !
  !        B       the solution vector  x .
  !
  !     Modification of LINPACK, by Peter Brown, LLNL.
  !     Written 7/20/83.  This version dated 6/20/01.
  !
  !     BLAS called: DAXPY
  !-----------------------------------------------------------------------
        INTEGER k , kb , l , nm1
        DOUBLE PRECISION t
  !
        nm1 = N - 1
  !
  !        Solve  A * x = b
  !        First solve  L*y = b
  !
        IF ( nm1.GE.1 ) THEN
           DO k = 1 , nm1
              l = Ipvt(k)
              t = B(l)
              IF ( l.NE.k ) THEN
                 B(l) = B(k)
                 B(k) = t
              ENDIF
              B(k+1) = B(k+1) + t*A(k+1,k)
           ENDDO
        ENDIF
  !
  !        Now solve  U*x = y
  !
        DO kb = 1 , N
           k = N + 1 - kb
           B(k) = B(k)/A(k,k)
           t = -B(k)
           CALL DAXPY(k-1,t,A(1,k),1,B(1),1)
        ENDDO
  !----------------------- End of Subroutine DHESL -----------------------
        END
  !*==DHEQR.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DHEQR
        SUBROUTINE DHEQR(A,Lda,N,Q,Info,Ijob)
        IMPLICIT NONE
  !*--DHEQR7998
        INTEGER Lda , N , Info , Ijob
        DOUBLE PRECISION A(Lda,*) , Q(*)
  !-----------------------------------------------------------------------
  !     This routine performs a QR decomposition of an upper
  !     Hessenberg matrix A.  There are two options available:
  !
  !          (1)  performing a fresh decomposition
  !          (2)  updating the QR factors by adding a row and a
  !               column to the matrix A.
  !-----------------------------------------------------------------------
  !     DHEQR decomposes an upper Hessenberg matrix by using Givens
  !     rotations.
  !
  !     On entry
  !
  !        A       DOUBLE PRECISION(LDA, N)
  !                the matrix to be decomposed.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  A .
  !
  !        N       INTEGER
  !                A is an (N+1) by N Hessenberg matrix.
  !
  !        IJOB    INTEGER
  !                = 1     means that a fresh decomposition of the
  !                        matrix A is desired.
  !                .ge. 2  means that the current decomposition of A
  !                        will be updated by the addition of a row
  !                        and a column.
  !     On return
  !
  !        A       the upper triangular matrix R.
  !                The factorization can be written Q*A = R, where
  !                Q is a product of Givens rotations and R is upper
  !                triangular.
  !
  !        Q       DOUBLE PRECISION(2*N)
  !                the factors c and s of each Givens rotation used
  !                in decomposing A.
  !
  !        INFO    INTEGER
  !                = 0  normal value.
  !                = k  if  A(k,k) .eq. 0.0 .  This is not an error
  !                     condition for this subroutine, but it does
  !                     indicate that DHELS will divide by zero
  !                     if called.
  !
  !     Modification of LINPACK, by Peter Brown, LLNL.
  !     Written 1/13/86.  This version dated 6/20/01.
  !-----------------------------------------------------------------------
        INTEGER i , iq , j , k , km1 , kp1 , nm1
        DOUBLE PRECISION c , s , t , t1 , t2
  !
        IF ( Ijob.GT.1 ) THEN
  !
  ! The old factorization of A will be updated.  A row and a column
  ! has been added to the matrix A.
  ! N by N-1 is now the old size of the matrix.
  !
           nm1 = N - 1
  !
  ! Multiply the new column by the N previous Givens rotations.
  !
           DO k = 1 , nm1
              i = 2*(k-1) + 1
              t1 = A(k,N)
              t2 = A(k+1,N)
              c = Q(i)
              s = Q(i+1)
              A(k,N) = c*t1 - s*t2
              A(k+1,N) = s*t1 + c*t2
           ENDDO
  !
  ! Complete update of decomposition by forming last Givens rotation,
  ! and multiplying it times the column vector (A(N,N), A(N+1,N)).
  !
           Info = 0
           t1 = A(N,N)
           t2 = A(N+1,N)
           IF ( t2.EQ.0.0D0 ) THEN
              c = 1.0D0
              s = 0.0D0
           ELSEIF ( ABS(t2).LT.ABS(t1) ) THEN
              t = t2/t1
              c = 1.0D0/SQRT(1.0D0+t*t)
              s = -c*t
           ELSE
              t = t1/t2
              s = -1.0D0/SQRT(1.0D0+t*t)
              c = -s*t
           ENDIF
           iq = 2*N - 1
           Q(iq) = c
           Q(iq+1) = s
           A(N,N) = c*t1 - s*t2
           IF ( A(N,N).EQ.0.0D0 ) Info = N
        ELSE
  !
  ! A new facorization is desired.
  !
  !     QR decomposition without pivoting
  !
           Info = 0
           DO k = 1 , N
              km1 = k - 1
              kp1 = k + 1
  !
  !           Compute kth column of R.
  !           First, multiply the kth column of A by the previous
  !           k-1 Givens rotations.
  !
              IF ( km1.GE.1 ) THEN
                 DO j = 1 , km1
                    i = 2*(j-1) + 1
                    t1 = A(j,k)
                    t2 = A(j+1,k)
                    c = Q(i)
                    s = Q(i+1)
                    A(j,k) = c*t1 - s*t2
                    A(j+1,k) = s*t1 + c*t2
                 ENDDO
              ENDIF
  !
  !           Compute Givens components c and s
  !
              iq = 2*km1 + 1
              t1 = A(k,k)
              t2 = A(kp1,k)
              IF ( t2.EQ.0.0D0 ) THEN
                 c = 1.0D0
                 s = 0.0D0
              ELSEIF ( ABS(t2).LT.ABS(t1) ) THEN
                 t = t2/t1
                 c = 1.0D0/SQRT(1.0D0+t*t)
                 s = -c*t
              ELSE
                 t = t1/t2
                 s = -1.0D0/SQRT(1.0D0+t*t)
                 c = -s*t
              ENDIF
              Q(iq) = c
              Q(iq+1) = s
              A(k,k) = c*t1 - s*t2
              IF ( A(k,k).EQ.0.0D0 ) Info = k
           ENDDO
           RETURN
        ENDIF
  !----------------------- End of Subroutine DHEQR -----------------------
        END
  !*==DHELS.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DHELS
        SUBROUTINE DHELS(A,Lda,N,Q,B)
        IMPLICIT NONE
  !*--DHELS8153
        INTEGER Lda , N
        DOUBLE PRECISION A(Lda,*) , B(*) , Q(*)
  !-----------------------------------------------------------------------
  ! This is part of the LINPACK routine DGESL with changes
  ! due to the fact that A is an upper Hessenberg matrix.
  !-----------------------------------------------------------------------
  !     DHELS solves the least squares problem
  !
  !           min (b-A*x, b-A*x)
  !
  !     using the factors computed by DHEQR.
  !
  !     On entry
  !
  !        A       DOUBLE PRECISION(LDA, N)
  !                the output from DHEQR which contains the upper
  !                triangular factor R in the QR decomposition of A.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  A .
  !
  !        N       INTEGER
  !                A is originally an (N+1) by N matrix.
  !
  !        Q       DOUBLE PRECISION(2*N)
  !                The coefficients of the N givens rotations
  !                used in the QR factorization of A.
  !
  !        B       DOUBLE PRECISION(N+1)
  !                the right hand side vector.
  !
  !     On return
  !
  !        B       the solution vector  x .
  !
  !     Modification of LINPACK, by Peter Brown, LLNL.
  !     Written 1/13/86.  This version dated 6/20/01.
  !
  !     BLAS called: DAXPY
  !-----------------------------------------------------------------------
        INTEGER iq , k , kb , kp1
        DOUBLE PRECISION c , s , t , t1 , t2
  !
  !        Minimize (b-A*x, b-A*x)
  !        First form Q*b.
  !
        DO k = 1 , N
           kp1 = k + 1
           iq = 2*(k-1) + 1
           c = Q(iq)
           s = Q(iq+1)
           t1 = B(k)
           t2 = B(kp1)
           B(k) = c*t1 - s*t2
           B(kp1) = s*t1 + c*t2
        ENDDO
  !
  !        Now solve  R*x = Q*b.
  !
        DO kb = 1 , N
           k = N + 1 - kb
           B(k) = B(k)/A(k,k)
           t = -B(k)
           CALL DAXPY(k-1,t,A(1,k),1,B(1),1)
        ENDDO
  !----------------------- End of Subroutine DHELS -----------------------
        END
  !*==DLHIN.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DLHIN
        SUBROUTINE DLHIN(Neq,N,T0,Y0,Ydot,F,Tout,Uround,Ewt,Itol,Atol,Y,  &
       &                 Temp,H0,Niter,Ier)
        IMPLICIT NONE
  !*--DLHIN8226
        EXTERNAL F
        DOUBLE PRECISION T0 , Y0 , Ydot , Tout , Uround , Ewt , Atol , Y ,&
       &                 Temp , H0
        INTEGER Neq , N , Itol , Niter , Ier
        DIMENSION Neq(*) , Y0(*) , Ydot(*) , Ewt(*) , Atol(*) , Y(*) ,    &
       &          Temp(*)
  !-----------------------------------------------------------------------
  ! Call sequence input -- NEQ, N, T0, Y0, YDOT, F, TOUT, UROUND,
  !                        EWT, ITOL, ATOL, Y, TEMP
  ! Call sequence output -- H0, NITER, IER
  ! Common block variables accessed -- None
  !
  ! Subroutines called by DLHIN: F, DCOPY
  ! Function routines called by DLHIN: DVNORM
  !-----------------------------------------------------------------------
  ! This routine computes the step size, H0, to be attempted on the
  ! first step, when the user has not supplied a value for this.
  !
  ! First we check that TOUT - T0 differs significantly from zero.  Then
  ! an iteration is done to approximate the initial second derivative
  ! and this is used to define H from WRMS-norm(H**2 * yddot / 2) = 1.
  ! A bias factor of 1/2 is applied to the resulting h.
  ! The sign of H0 is inferred from the initial values of TOUT and T0.
  !
  ! Communication with DLHIN is done with the following variables:
  !
  ! NEQ    = NEQ array of solver, passed to F.
  ! N      = size of ODE system, input.
  ! T0     = initial value of independent variable, input.
  ! Y0     = vector of initial conditions, input.
  ! YDOT   = vector of initial first derivatives, input.
  ! F      = name of subroutine for right-hand side f(t,y), input.
  ! TOUT   = first output value of independent variable
  ! UROUND = machine unit roundoff
  ! EWT, ITOL, ATOL = error weights and tolerance parameters
  !                   as described in the driver routine, input.
  ! Y, TEMP = work arrays of length N.
  ! H0     = step size to be attempted, output.
  ! NITER  = number of iterations (and of f evaluations) to compute H0,
  !          output.
  ! IER    = the error flag, returned with the value
  !          IER = 0  if no trouble occurred, or
  !          IER = -1 if TOUT and t0 are considered too close to proceed.
  !-----------------------------------------------------------------------
  !
  ! Type declarations for local variables --------------------------------
  !
        DOUBLE PRECISION afi , atoli , delyi , half , hg , hlb , hnew ,   &
       &                 hrat , hub , hun , pt1 , t1 , tdist , tround ,   &
       &                 two , yddnrm
        INTEGER i , iter
  !-----------------------------------------------------------------------
  ! The following Fortran-77 declaration is to cause the values of the
  ! listed (local) variables to be saved between calls to this integrator.
  !-----------------------------------------------------------------------
        SAVE half , hun , pt1 , two
        DATA half/0.5D0/ , hun/100.0D0/ , pt1/0.1D0/ , two/2.0D0/
  !
        Niter = 0
        tdist = ABS(Tout-T0)
        tround = Uround*MAX(ABS(T0),ABS(Tout))
        IF ( tdist.LT.two*tround ) THEN
  ! Error return for TOUT - T0 too small. --------------------------------
           Ier = -1
           GOTO 99999
        ELSE
  !
  ! Set a lower bound on H based on the roundoff level in T0 and TOUT. ---
           hlb = hun*tround
  ! Set an upper bound on H based on TOUT-T0 and the initial Y and YDOT. -
           hub = pt1*tdist
           atoli = Atol(1)
           DO i = 1 , N
              IF ( Itol.EQ.2 .OR. Itol.EQ.4 ) atoli = Atol(i)
              delyi = pt1*ABS(Y0(i)) + atoli
              afi = ABS(Ydot(i))
              IF ( afi*hub.GT.delyi ) hub = delyi/afi
           ENDDO
  !
  ! Set initial guess for H as geometric mean of upper and lower bounds. -
           iter = 0
           hg = SQRT(hlb*hub)
  ! If the bounds have crossed, exit with the mean value. ----------------
           IF ( hub.LT.hlb ) THEN
              H0 = hg
              GOTO 200
           ENDIF
  !
  ! Looping point for iteration. -----------------------------------------
  ! Estimate the second derivative as a difference quotient in f. --------
   50      t1 = T0 + hg
           DO i = 1 , N
              Y(i) = Y0(i) + hg*Ydot(i)
           ENDDO
           CALL F(Neq,t1,Y,Temp)
           DO i = 1 , N
              Temp(i) = (Temp(i)-Ydot(i))/hg
           ENDDO
           yddnrm = DVNORM(N,Temp,Ewt)
  ! Get the corresponding new value of H. --------------------------------
           IF ( yddnrm*hub*hub.GT.two ) THEN
              hnew = SQRT(two/yddnrm)
           ELSE
              hnew = SQRT(hg*hub)
           ENDIF
           iter = iter + 1
  !-----------------------------------------------------------------------
  ! Test the stopping conditions.
  ! Stop if the new and previous H values differ by a factor of .lt. 2.
  ! Stop if four iterations have been done.  Also, stop with previous H
  ! if hnew/hg .gt. 2 after first iteration, as this probably means that
  ! the second derivative value is bad because of cancellation error.
  !-----------------------------------------------------------------------
           IF ( iter.LT.4 ) THEN
              hrat = hnew/hg
              IF ( (hrat.LE.half) .OR. (hrat.GE.two) ) THEN
                 IF ( (iter.GE.2) .AND. (hnew.GT.two*hg) ) THEN
                    hnew = hg
                    GOTO 100
                 ENDIF
                 hg = hnew
                 GOTO 50
              ENDIF
           ENDIF
  !
  ! Iteration done.  Apply bounds, bias factor, and sign. ----------------
   100     H0 = hnew*half
           IF ( H0.LT.hlb ) H0 = hlb
           IF ( H0.GT.hub ) H0 = hub
        ENDIF
   200  H0 = SIGN(H0,Tout-T0)
  ! Restore Y array from Y0, then exit. ----------------------------------
        CALL DCOPY(N,Y0,1,Y,1)
        Niter = iter
        Ier = 0
        RETURN
  !----------------------- End of Subroutine DLHIN -----------------------
  99999 END
  !*==DSTOKA.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DSTOKA
        SUBROUTINE DSTOKA(Neq,Y,Yh,Nyh,Yh1,Ewt,Savf,Savx,Acor,Wm,Iwm,F,   &
       &                  JAC,PSOL)
        IMPLICIT NONE
  !*--DSTOKA8370
  !*** Start of declarations inserted by SPAG
        !INTEGER JAC
        !REAL PSOL
  !*** End of declarations inserted by SPAG
        EXTERNAL F , JAC , PSOL
        INTEGER Neq , Nyh , Iwm
        DOUBLE PRECISION Y , Yh , Yh1 , Ewt , Savf , Savx , Acor , Wm
        DIMENSION Neq(*) , Y(*) , Yh(Nyh,*) , Yh1(*) , Ewt(*) , Savf(*) , &
       &          Savx(*) , Acor(*) , Wm(*) , Iwm(*)
        INTEGER IOWnd , IALth , IPUp , LMAx , MEO , NQNyh , NSLp , ICF ,  &
       &        IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , LEWt ,  &
       &        LACor , LSAvf , LWM , LIWm , METh , MITer , MAXord ,      &
       &        MAXcor , MSBp , MXNcf , N , NQ , NST , NFE , NJE , NQU
        INTEGER NEWt , NSFi , NSLj , NJEv
        INTEGER JPRe , JACflg , LOCwp , LOCiwp , LSAvx , KMP , MAXl ,     &
       &        MNEwt , NNI , NLI , NPS , NCFn , NCFl
        DOUBLE PRECISION CONit , CRAte , EL , ELCo , HOLd , RMAx , TESco ,&
       &                 CCMax , EL0 , H , HMIn , HMXi , HU , RC , TN ,   &
       &                 UROund
        DOUBLE PRECISION STIfr
        DOUBLE PRECISION DELt , EPCon , SQRtn , RSQrtn
        COMMON /DLS001/ CONit , CRAte , EL(13) , ELCo(13,12) , HOLd ,     &
       &                RMAx , TESco(3,12) , CCMax , EL0 , H , HMIn ,     &
       &                HMXi , HU , RC , TN , UROund , IOWnd(6) , IALth , &
       &                IPUp , LMAx , MEO , NQNyh , NSLp , ICF , IERpj ,  &
       &                IERsl , JCUr , JSTart , KFLag , L , LYH , LEWt ,  &
       &                LACor , LSAvf , LWM , LIWm , METh , MITer ,       &
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        COMMON /DLS002/ STIfr , NEWt , NSFi , NSLj , NJEv
        COMMON /DLPK01/ DELt , EPCon , SQRtn , RSQrtn , JPRe , JACflg ,   &
       &                LOCwp , LOCiwp , LSAvx , KMP , MAXl , MNEwt ,     &
       &                NNI , NLI , NPS , NCFn , NCFl
  !-----------------------------------------------------------------------
  ! DSTOKA performs one step of the integration of an initial value
  ! problem for a system of Ordinary Differential Equations.
  !
  ! This routine was derived from Subroutine DSTODPK in the DLSODPK
  ! package by the addition of automatic functional/Newton iteration
  ! switching and logic for re-use of Jacobian data.
  !-----------------------------------------------------------------------
  ! Note: DSTOKA is independent of the value of the iteration method
  ! indicator MITER, when this is .ne. 0, and hence is independent
  ! of the type of chord method used, or the Jacobian structure.
  ! Communication with DSTOKA is done with the following variables:
  !
  ! NEQ    = integer array containing problem size in NEQ(1), and
  !          passed as the NEQ argument in all calls to F and JAC.
  ! Y      = an array of length .ge. N used as the Y argument in
  !          all calls to F and JAC.
  ! YH     = an NYH by LMAX array containing the dependent variables
  !          and their approximate scaled derivatives, where
  !          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
  !          j-th derivative of y(i), scaled by H**j/factorial(j)
  !          (j = 0,1,...,NQ).  On entry for the first step, the first
  !          two columns of YH must be set from the initial values.
  ! NYH    = a constant integer .ge. N, the first dimension of YH.
  ! YH1    = a one-dimensional array occupying the same space as YH.
  ! EWT    = an array of length N containing multiplicative weights
  !          for local error measurements.  Local errors in y(i) are
  !          compared to 1.0/EWT(i) in various error tests.
  ! SAVF   = an array of working storage, of length N.
  !          Also used for input of YH(*,MAXORD+2) when JSTART = -1
  !          and MAXORD .lt. the current order NQ.
  ! SAVX   = an array of working storage, of length N.
  ! ACOR   = a work array of length N, used for the accumulated
  !          corrections.  On a successful return, ACOR(i) contains
  !          the estimated one-step local error in y(i).
  ! WM,IWM = real and integer work arrays associated with matrix
  !          operations in chord iteration (MITER .ne. 0).
  ! CCMAX  = maximum relative change in H*EL0 before DSETPK is called.
  ! H      = the step size to be attempted on the next step.
  !          H is altered by the error control algorithm during the
  !          problem.  H can be either positive or negative, but its
  !          sign must remain constant throughout the problem.
  ! HMIN   = the minimum absolute value of the step size H to be used.
  ! HMXI   = inverse of the maximum absolute value of H to be used.
  !          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
  !          HMIN and HMXI may be changed at any time, but will not
  !          take effect until the next change of H is considered.
  ! TN     = the independent variable. TN is updated on each step taken.
  ! JSTART = an integer used for input only, with the following
  !          values and meanings:
  !               0  perform the first step.
  !           .gt.0  take a new step continuing from the last.
  !              -1  take the next step with a new value of H, MAXORD,
  !                    N, METH, MITER, and/or matrix parameters.
  !              -2  take the next step with a new value of H,
  !                    but with other inputs unchanged.
  !          On return, JSTART is set to 1 to facilitate continuation.
  ! KFLAG  = a completion code with the following meanings:
  !               0  the step was succesful.
  !              -1  the requested error could not be achieved.
  !              -2  corrector convergence could not be achieved.
  !              -3  fatal error in DSETPK or DSOLPK.
  !          A return with KFLAG = -1 or -2 means either
  !          ABS(H) = HMIN or 10 consecutive failures occurred.
  !          On a return with KFLAG negative, the values of TN and
  !          the YH array are as of the beginning of the last
  !          step, and H is the last step size attempted.
  ! MAXORD = the maximum order of integration method to be allowed.
  ! MAXCOR = the maximum number of corrector iterations allowed.
  ! MSBP   = maximum number of steps between DSETPK calls (MITER .gt. 0).
  ! MXNCF  = maximum number of convergence failures allowed.
  ! METH/MITER = the method flags.  See description in driver.
  ! N      = the number of first-order differential equations.
  !-----------------------------------------------------------------------
        INTEGER i , i1 , iredo , iret , j , jb , jok , m , ncf , newq ,   &
       &        nslow
        DOUBLE PRECISION dcon , ddn , del , delp , drc , dsm , dup ,      &
       &                 exdn , exsm , exup , dfnorm , r , rh , rhdn ,    &
       &                 rhsm , rhup , roc , stiff , told
  !
        KFLag = 0
        told = TN
        ncf = 0
        IERpj = 0
        IERsl = 0
        JCUr = 0
        ICF = 0
        delp = 0.0D0
        IF ( JSTart.GT.0 ) GOTO 500
        IF ( JSTart.EQ.-1 ) THEN
  !-----------------------------------------------------------------------
  ! The following block handles preliminaries needed when JSTART = -1.
  ! IPUP is set to MITER to force a matrix update.
  ! If an order increase is about to be considered (IALTH = 1),
  ! IALTH is reset to 2 to postpone consideration one more step.
  ! If the caller has changed METH, DCFODE is called to reset
  ! the coefficients of the method.
  ! If the caller has changed MAXORD to a value less than the current
  ! order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
  ! If H is to be changed, YH must be rescaled.
  ! If H or METH is being changed, IALTH is reset to L = NQ + 1
  ! to prevent further changes in H for that many steps.
  !-----------------------------------------------------------------------
           IPUp = MITer
           LMAx = MAXord + 1
           IF ( IALth.EQ.1 ) IALth = 2
           IF ( METh.NE.MEO ) THEN
              CALL DCFODE(METh,ELCo,TESco)
              MEO = METh
              IF ( NQ.LE.MAXord ) THEN
                 IALth = L
                 iret = 1
                 GOTO 100
              ENDIF
           ELSEIF ( NQ.LE.MAXord ) THEN
              GOTO 200
           ENDIF
           NQ = MAXord
           L = LMAx
           DO i = 1 , L
              EL(i) = ELCo(i,NQ)
           ENDDO
           NQNyh = NQ*Nyh
           RC = RC*EL(1)/EL0
           EL0 = EL(1)
           CONit = 0.5D0/(NQ+2)
           EPCon = CONit*TESco(2,NQ)
           ddn = DVNORM(N,Savf,Ewt)/TESco(1,L)
           exdn = 1.0D0/L
           rhdn = 1.0D0/(1.3D0*ddn**exdn+0.0000013D0)
           rh = MIN(rhdn,1.0D0)
           iredo = 3
           IF ( H.EQ.HOLd ) GOTO 300
           rh = MIN(rh,ABS(H/HOLd))
           H = HOLd
           GOTO 400
        ELSE
           IF ( JSTart.EQ.-2 ) GOTO 200
  !-----------------------------------------------------------------------
  ! On the first call, the order is set to 1, and other variables are
  ! initialized.  RMAX is the maximum ratio by which H can be increased
  ! in a single step.  It is initially 1.E4 to compensate for the small
  ! initial H, but then is normally equal to 10.  If a failure
  ! occurs (in corrector convergence or error test), RMAX is set at 2
  ! for the next increase.
  !-----------------------------------------------------------------------
           LMAx = MAXord + 1
           NQ = 1
           L = 2
           IALth = 2
           RMAx = 10000.0D0
           RC = 0.0D0
           EL0 = 1.0D0
           CRAte = 0.7D0
           HOLd = H
           MEO = METh
           NSLp = 0
           NSLj = 0
           IPUp = 0
           iret = 3
           NEWt = 0
           STIfr = 0.0D0
  !-----------------------------------------------------------------------
  ! DCFODE is called to get all the integration coefficients for the
  ! current METH.  Then the EL vector and related constants are reset
  ! whenever the order NQ is changed, or at the start of the problem.
  !-----------------------------------------------------------------------
           CALL DCFODE(METh,ELCo,TESco)
        ENDIF
   100  DO i = 1 , L
           EL(i) = ELCo(i,NQ)
        ENDDO
        NQNyh = NQ*Nyh
        RC = RC*EL(1)/EL0
        EL0 = EL(1)
        CONit = 0.5D0/(NQ+2)
        EPCon = CONit*TESco(2,NQ)
        GOTO (200,300,500) , iret
  !-----------------------------------------------------------------------
  ! If H is being changed, the H ratio RH is checked against
  ! RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
  ! L = NQ + 1 to prevent a change of H for that many steps, unless
  ! forced by a convergence or error test failure.
  !-----------------------------------------------------------------------
   200  IF ( H.EQ.HOLd ) GOTO 500
        rh = H/HOLd
        H = HOLd
        iredo = 3
        GOTO 400
   300  rh = MAX(rh,HMIn/ABS(H))
   400  rh = MIN(rh,RMAx)
        rh = rh/MAX(1.0D0,ABS(H)*HMXi*rh)
        r = 1.0D0
        DO j = 2 , L
           r = r*rh
           DO i = 1 , N
              Yh(i,j) = Yh(i,j)*r
           ENDDO
        ENDDO
        H = H*rh
        RC = RC*rh
        IALth = L
        IF ( iredo.EQ.0 ) THEN
           RMAx = 10.0D0
           GOTO 1300
        ENDIF
  !-----------------------------------------------------------------------
  ! This section computes the predicted values by effectively
  ! multiplying the YH array by the Pascal triangle matrix.
  ! The flag IPUP is set according to whether matrix data is involved
  ! (NEWT .gt. 0 .and. JACFLG .ne. 0) or not, to trigger a call to DSETPK.
  ! IPUP is set to MITER when RC differs from 1 by more than CCMAX,
  ! and at least every MSBP steps, when JACFLG = 1.
  ! RC is the ratio of new to old values of the coefficient  H*EL(1).
  !-----------------------------------------------------------------------
   500  IF ( NEWt.EQ.0 .OR. JACflg.EQ.0 ) THEN
           drc = 0.0D0
           IPUp = 0
           CRAte = 0.7D0
        ELSE
           drc = ABS(RC-1.0D0)
           IF ( drc.GT.CCMax ) IPUp = MITer
           IF ( NST.GE.NSLp+MSBp ) IPUp = MITer
        ENDIF
        TN = TN + H
        i1 = NQNyh + 1
        DO jb = 1 , NQ
           i1 = i1 - Nyh
  !DIR$ IVDEP
           DO i = i1 , NQNyh
              Yh1(i) = Yh1(i) + Yh1(i+Nyh)
           ENDDO
        ENDDO
  !-----------------------------------------------------------------------
  ! Up to MAXCOR corrector iterations are taken.  A convergence test is
  ! made on the RMS-norm of each correction, weighted by the error
  ! weight vector EWT.  The sum of the corrections is accumulated in the
  ! vector ACOR(i).  The YH array is not altered in the corrector loop.
  ! Within the corrector loop, an estimated rate of convergence (ROC)
  ! and a stiffness ratio estimate (STIFF) are kept.  Corresponding
  ! global estimates are kept as CRATE and stifr.
  !-----------------------------------------------------------------------
   600  m = 0
        MNEwt = 0
        stiff = 0.0D0
        roc = 0.05D0
        nslow = 0
        DO i = 1 , N
           Y(i) = Yh(i,1)
        ENDDO
        CALL F(Neq,TN,Y,Savf)
        NFE = NFE + 1
        IF ( NEWt.NE.0 .AND. IPUp.GT.0 ) THEN
  !-----------------------------------------------------------------------
  ! If indicated, DSETPK is called to update any matrix data needed,
  ! before starting the corrector iteration.
  ! JOK is set to indicate if the matrix data need not be recomputed.
  ! IPUP is set to 0 as an indicator that the matrix data is up to date.
  !-----------------------------------------------------------------------
           jok = 1
           IF ( NST.EQ.0 .OR. NST.GT.NSLj+50 ) jok = -1
           IF ( ICF.EQ.1 .AND. drc.LT.0.2D0 ) jok = -1
           IF ( ICF.EQ.2 ) jok = -1
           IF ( jok.EQ.-1 ) THEN
              NSLj = NST
              NJEv = NJEv + 1
           ENDIF
           CALL DSETPK(Neq,Y,Yh1,Ewt,Acor,Savf,jok,Wm,Iwm,F,JAC)
           IPUp = 0
           RC = 1.0D0
           drc = 0.0D0
           NSLp = NST
           CRAte = 0.7D0
           IF ( IERpj.NE.0 ) GOTO 900
        ENDIF
        DO i = 1 , N
           Acor(i) = 0.0D0
        ENDDO
   700  IF ( NEWt.NE.0 ) THEN
  !-----------------------------------------------------------------------
  ! In the case of the chord method, compute the corrector error,
  ! and solve the linear system with that as right-hand side and
  ! P as coefficient matrix.  STIFF is set to the ratio of the norms
  ! of the residual and the correction vector.
  !-----------------------------------------------------------------------
           DO i = 1 , N
              Savx(i) = H*Savf(i) - (Yh(i,2)+Acor(i))
           ENDDO
           dfnorm = DVNORM(N,Savx,Ewt)
           CALL DSOLPK(Neq,Y,Savf,Savx,Ewt,Wm,Iwm,F,PSOL)
           IF ( IERsl.LT.0 ) GOTO 900
           IF ( IERsl.GT.0 ) GOTO 800
           del = DVNORM(N,Savx,Ewt)
           IF ( del.GT.1.0D-8 ) stiff = MAX(stiff,dfnorm/del)
           DO i = 1 , N
              Acor(i) = Acor(i) + Savx(i)
              Y(i) = Yh(i,1) + EL(1)*Acor(i)
           ENDDO
        ELSE
  !-----------------------------------------------------------------------
  ! In the case of functional iteration, update Y directly from
  ! the result of the last function evaluation, and STIFF is set to 1.0.
  !-----------------------------------------------------------------------
           DO i = 1 , N
              Savf(i) = H*Savf(i) - Yh(i,2)
              Y(i) = Savf(i) - Acor(i)
           ENDDO
           del = DVNORM(N,Y,Ewt)
           DO i = 1 , N
              Y(i) = Yh(i,1) + EL(1)*Savf(i)
              Acor(i) = Savf(i)
           ENDDO
           stiff = 1.0D0
        ENDIF
  !-----------------------------------------------------------------------
  ! Test for convergence.  If M .gt. 0, an estimate of the convergence
  ! rate constant is made for the iteration switch, and is also used
  ! in the convergence test.   If the iteration seems to be diverging or
  ! converging at a slow rate (.gt. 0.8 more than once), it is stopped.
  !-----------------------------------------------------------------------
        IF ( m.NE.0 ) THEN
           roc = MAX(0.05D0,del/delp)
           CRAte = MAX(0.2D0*CRAte,roc)
        ENDIF
        dcon = del*MIN(1.0D0,1.5D0*CRAte)/EPCon
        IF ( dcon.LE.1.0D0 ) THEN
  !-----------------------------------------------------------------------
  ! The corrector has converged.  JCUR is set to 0 to signal that the
  ! preconditioner involved may need updating later.
  ! The stiffness ratio STIFR is updated using the latest STIFF value.
  ! The local error test is made and control passes to statement 500
  ! if it fails.
  !-----------------------------------------------------------------------
           JCUr = 0
           IF ( NEWt.GT.0 ) STIfr = 0.5D0*(STIfr+stiff)
           IF ( m.EQ.0 ) dsm = del/TESco(2,NQ)
           IF ( m.GT.0 ) dsm = DVNORM(N,Acor,Ewt)/TESco(2,NQ)
           IF ( dsm.GT.1.0D0 ) THEN
  !-----------------------------------------------------------------------
  ! The error test failed.  KFLAG keeps track of multiple failures.
  ! Restore TN and the YH array to their previous values, and prepare
  ! to try the step again.  Compute the optimum step size for this or
  ! one lower order.  After 2 or more failures, H is forced to decrease
  ! by a factor of 0.2 or less.
  !-----------------------------------------------------------------------
              KFLag = KFLag - 1
              TN = told
              i1 = NQNyh + 1
              DO jb = 1 , NQ
                 i1 = i1 - Nyh
  !DIR$ IVDEP
                 DO i = i1 , NQNyh
                    Yh1(i) = Yh1(i) - Yh1(i+Nyh)
                 ENDDO
              ENDDO
              RMAx = 2.0D0
              IF ( ABS(H).LE.HMIn*1.00001D0 ) THEN
  !-----------------------------------------------------------------------
  ! All returns are made through this section.  H is saved in HOLD
  ! to allow the caller to change H on the next step.
  !-----------------------------------------------------------------------
                 KFLag = -1
                 GOTO 1400
              ELSEIF ( KFLag.LE.-3 ) THEN
  !-----------------------------------------------------------------------
  ! Control reaches this section if 3 or more failures have occured.
  ! If 10 failures have occurred, exit with KFLAG = -1.
  ! It is assumed that the derivatives that have accumulated in the
  ! YH array have errors of the wrong order.  Hence the first
  ! derivative is recomputed, and the order is set to 1.  Then
  ! H is reduced by a factor of 10, and the step is retried,
  ! until it succeeds or H reaches HMIN.
  !-----------------------------------------------------------------------
                 IF ( KFLag.EQ.-10 ) THEN
                    KFLag = -1
                    GOTO 1400
                 ELSE
                    rh = 0.1D0
                    rh = MAX(HMIn/ABS(H),rh)
                    H = H*rh
                    DO i = 1 , N
                       Y(i) = Yh(i,1)
                    ENDDO
                    CALL F(Neq,TN,Y,Savf)
                    NFE = NFE + 1
                    DO i = 1 , N
                       Yh(i,2) = H*Savf(i)
                    ENDDO
                    IPUp = MITer
                    IALth = 5
                    IF ( NQ.EQ.1 ) GOTO 500
                    NQ = 1
                    L = 2
                    iret = 3
                    GOTO 100
                 ENDIF
              ELSE
                 iredo = 2
                 rhup = 0.0D0
                 GOTO 1000
              ENDIF
           ELSE
  !-----------------------------------------------------------------------
  ! After a successful step, update the YH array.
  ! If Newton iteration is being done and STIFR is less than 1.5,
  ! then switch to functional iteration.
  ! Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
  ! If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
  ! use in a possible order increase on the next step.
  ! If a change in H is considered, an increase or decrease in order
  ! by one is considered also.  A change in H is made only if it is by a
  ! factor of at least 1.1.  If not, IALTH is set to 3 to prevent
  ! testing for that many steps.
  !-----------------------------------------------------------------------
              KFLag = 0
              iredo = 0
              NST = NST + 1
              IF ( NEWt.EQ.0 ) NSFi = NSFi + 1
              IF ( NEWt.GT.0 .AND. STIfr.LT.1.5D0 ) NEWt = 0
              HU = H
              NQU = NQ
              DO j = 1 , L
                 DO i = 1 , N
                    Yh(i,j) = Yh(i,j) + EL(j)*Acor(i)
                 ENDDO
              ENDDO
              IALth = IALth - 1
              IF ( IALth.EQ.0 ) THEN
  !-----------------------------------------------------------------------
  ! Regardless of the success or failure of the step, factors
  ! RHDN, RHSM, and RHUP are computed, by which H could be multiplied
  ! at order NQ - 1, order NQ, or order NQ + 1, respectively.
  ! in the case of failure, RHUP = 0.0 to avoid an order increase.
  ! the largest of these is determined and the new order chosen
  ! accordingly.  If the order is to be increased, we compute one
  ! additional scaled derivative.
  !-----------------------------------------------------------------------
                 rhup = 0.0D0
                 IF ( L.NE.LMAx ) THEN
                    DO i = 1 , N
                       Savf(i) = Acor(i) - Yh(i,LMAx)
                    ENDDO
                    dup = DVNORM(N,Savf,Ewt)/TESco(3,NQ)
                    exup = 1.0D0/(L+1)
                    rhup = 1.0D0/(1.4D0*dup**exup+0.0000014D0)
                 ENDIF
                 GOTO 1000
              ELSE
                 IF ( IALth.LE.1 ) THEN
                    IF ( L.NE.LMAx ) THEN
                       DO i = 1 , N
                          Yh(i,LMAx) = Acor(i)
                       ENDDO
                    ENDIF
                 ENDIF
                 GOTO 1300
              ENDIF
           ENDIF
        ELSE
           m = m + 1
           IF ( m.NE.MAXcor ) THEN
              IF ( m.LT.2 .OR. del.LE.2.0D0*delp ) THEN
                 IF ( roc.LE.10.0D0 ) THEN
                    IF ( roc.GT.0.8D0 ) nslow = nslow + 1
                    IF ( nslow.LT.2 ) THEN
                       MNEwt = m
                       delp = del
                       CALL F(Neq,TN,Y,Savf)
                       NFE = NFE + 1
                       GOTO 700
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
  !-----------------------------------------------------------------------
  ! The corrector iteration failed to converge.
  ! If functional iteration is being done (NEWT = 0) and MITER .gt. 0
  ! (and this is not the first step), then switch to Newton
  ! (NEWT = MITER), and retry the step.  (Setting STIFR = 1023 insures
  ! that a switch back will not occur for 10 step attempts.)
  ! If Newton iteration is being done, but using a preconditioner that
  ! is out of date (JACFLG .ne. 0 .and. JCUR = 0), then signal for a
  ! re-evalutation of the preconditioner, and retry the step.
  ! In all other cases, the YH array is retracted to its values
  ! before prediction, and H is reduced, if possible.  If H cannot be
  ! reduced or MXNCF failures have occurred, exit with KFLAG = -2.
  !-----------------------------------------------------------------------
   800  ICF = 1
        IF ( NEWt.EQ.0 ) THEN
           IF ( NST.EQ.0 ) GOTO 900
           IF ( MITer.EQ.0 ) GOTO 900
           NEWt = MITer
           STIfr = 1023.0D0
           IPUp = MITer
           GOTO 600
        ENDIF
        IF ( JCUr.NE.1 .AND. JACflg.NE.0 ) THEN
           IPUp = MITer
           GOTO 600
        ENDIF
   900  ICF = 2
        ncf = ncf + 1
        NCFn = NCFn + 1
        RMAx = 2.0D0
        TN = told
        i1 = NQNyh + 1
        DO jb = 1 , NQ
           i1 = i1 - Nyh
  !DIR$ IVDEP
           DO i = i1 , NQNyh
              Yh1(i) = Yh1(i) - Yh1(i+Nyh)
           ENDDO
        ENDDO
        IF ( IERpj.LT.0 .OR. IERsl.LT.0 ) THEN
           KFLag = -3
           GOTO 1400
        ELSEIF ( ABS(H).LE.HMIn*1.00001D0 ) THEN
           KFLag = -2
           GOTO 1400
        ELSEIF ( ncf.EQ.MXNcf ) THEN
           KFLag = -2
           GOTO 1400
        ELSE
           rh = 0.5D0
           IPUp = MITer
           iredo = 1
           GOTO 300
        ENDIF
   1000 exsm = 1.0D0/L
        rhsm = 1.0D0/(1.2D0*dsm**exsm+0.0000012D0)
        rhdn = 0.0D0
        IF ( NQ.NE.1 ) THEN
           ddn = DVNORM(N,Yh(1,L),Ewt)/TESco(1,NQ)
           exdn = 1.0D0/NQ
           rhdn = 1.0D0/(1.3D0*ddn**exdn+0.0000013D0)
        ENDIF
        IF ( rhsm.GE.rhup ) THEN
           IF ( rhsm.GE.rhdn ) THEN
              newq = NQ
              rh = rhsm
              GOTO 1100
           ENDIF
        ELSEIF ( rhup.GT.rhdn ) THEN
           newq = L
           rh = rhup
           IF ( rh.LT.1.1D0 ) THEN
              IALth = 3
              GOTO 1300
           ELSE
              r = EL(L)/L
              DO i = 1 , N
                 Yh(i,newq+1) = Acor(i)*r
              ENDDO
              GOTO 1200
           ENDIF
        ENDIF
        newq = NQ - 1
        rh = rhdn
        IF ( KFLag.LT.0 .AND. rh.GT.1.0D0 ) rh = 1.0D0
   1100 IF ( (KFLag.EQ.0) .AND. (rh.LT.1.1D0) ) THEN
           IALth = 3
           GOTO 1300
        ELSE
           IF ( KFLag.LE.-2 ) rh = MIN(rh,0.2D0)
  !-----------------------------------------------------------------------
  ! If there is a change of order, reset NQ, L, and the coefficients.
  ! In any case H is reset according to RH and the YH array is rescaled.
  ! Then exit from 690 if the step was OK, or redo the step otherwise.
  !-----------------------------------------------------------------------
           IF ( newq.EQ.NQ ) GOTO 300
        ENDIF
   1200 NQ = newq
        L = NQ + 1
        iret = 2
        GOTO 100
   1300 r = 1.0D0/TESco(2,NQU)
        DO i = 1 , N
           Acor(i) = Acor(i)*r
        ENDDO
   1400 HOLd = H
        JSTart = 1
  !----------------------- End of Subroutine DSTOKA ----------------------
        END
  !*==DSETPK.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DSETPK
        SUBROUTINE DSETPK(Neq,Y,Ysv,Ewt,Ftem,Savf,Jok,Wm,Iwm,F,JAC)
        IMPLICIT NONE
  !*--DSETPK8992
  !*** Start of declarations inserted by SPAG
        !REAL F
  !*** End of declarations inserted by SPAG
        EXTERNAL F , JAC
        INTEGER Neq , Jok , Iwm
        DOUBLE PRECISION Y , Ysv , Ewt , Ftem , Savf , Wm
        DIMENSION Neq(*) , Y(*) , Ysv(*) , Ewt(*) , Ftem(*) , Savf(*) ,   &
       &          Wm(*) , Iwm(*)
        INTEGER IOWnd , IOWns , ICF , IERpj , IERsl , JCUr , JSTart ,     &
       &        KFLag , L , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,     &
       &        METh , MITer , MAXord , MAXcor , MSBp , MXNcf , N , NQ ,  &
       &        NST , NFE , NJE , NQU
        INTEGER JPRe , JACflg , LOCwp , LOCiwp , LSAvx , KMP , MAXl ,     &
       &        MNEwt , NNI , NLI , NPS , NCFn , NCFl
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        DOUBLE PRECISION DELt , EPCon , SQRtn , RSQrtn
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , IOWnd(6) , IOWns(6) , ICF ,    &
       &                IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , &
       &                LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,&
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        COMMON /DLPK01/ DELt , EPCon , SQRtn , RSQrtn , JPRe , JACflg ,   &
       &                LOCwp , LOCiwp , LSAvx , KMP , MAXl , MNEwt ,     &
       &                NNI , NLI , NPS , NCFn , NCFl
  !-----------------------------------------------------------------------
  ! DSETPK is called by DSTOKA to interface with the user-supplied
  ! routine JAC, to compute and process relevant parts of
  ! the matrix P = I - H*EL(1)*J , where J is the Jacobian df/dy,
  ! as need for preconditioning matrix operations later.
  !
  ! In addition to variables described previously, communication
  ! with DSETPK uses the following:
  ! Y     = array containing predicted values on entry.
  ! YSV   = array containing predicted y, to be saved (YH1 in DSTOKA).
  ! FTEM  = work array of length N (ACOR in DSTOKA).
  ! SAVF  = array containing f evaluated at predicted y.
  ! JOK   = input flag showing whether it was judged that Jacobian matrix
  !         data need not be recomputed (JOK = 1) or needs to be
  !         (JOK = -1).
  ! WM    = real work space for matrices.
  !         Space for preconditioning data starts at WM(LOCWP).
  ! IWM   = integer work space.
  !         Space for preconditioning data starts at IWM(LOCIWP).
  ! IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
  !         JAC returned an error flag.
  ! JCUR  = output flag to indicate whether the matrix data involved
  !         is now current (JCUR = 1) or not (JCUR = 0).
  ! This routine also uses Common variables EL0, H, TN, IERPJ, JCUR, NJE.
  !-----------------------------------------------------------------------
        INTEGER ier
        DOUBLE PRECISION hl0
  !
        IERpj = 0
        JCUr = 0
        IF ( Jok.EQ.-1 ) JCUr = 1
        hl0 = EL0*H
        CALL JAC(F,Neq,TN,Y,Ysv,Ewt,Savf,Ftem,hl0,Jok,Wm(LOCwp),          &
       &         Iwm(LOCiwp),ier)
        NJE = NJE + 1
        IF ( ier.EQ.0 ) RETURN
        IERpj = 1
  !----------------------- End of Subroutine DSETPK ----------------------
        END
  !*==DSRCKR.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DSRCKR
        SUBROUTINE DSRCKR(Rsav,Isav,Job)
        IMPLICIT NONE
  !*--DSRCKR9062
  !-----------------------------------------------------------------------
  ! This routine saves or restores (depending on JOB) the contents of
  ! the Common blocks DLS001, DLS002, DLSR01, DLPK01, which
  ! are used internally by the DLSODKR solver.
  !
  ! RSAV = real array of length 228 or more.
  ! ISAV = integer array of length 63 or more.
  ! JOB  = flag indicating to save or restore the Common blocks:
  !        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
  !        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
  !        A call with JOB = 2 presumes a prior call with JOB = 1.
  !-----------------------------------------------------------------------
        INTEGER Isav , Job
        INTEGER ILS , ILS2 , ILSr , ILSp
        INTEGER i , ioff , lenilp , lenrlp , lenils , lenrls , lenilr ,   &
       &        lenrlr
        DOUBLE PRECISION Rsav , RLS , RLS2 , RLSr , RLSp
        DIMENSION Rsav(*) , Isav(*)
        SAVE lenrls , lenils , lenrlp , lenilp , lenrlr , lenilr
        COMMON /DLS001/ RLS(218) , ILS(37)
        COMMON /DLS002/ RLS2 , ILS2(4)
        COMMON /DLSR01/ RLSr(5) , ILSr(9)
        COMMON /DLPK01/ RLSp(4) , ILSp(13)
        DATA lenrls/218/ , lenils/37/ , lenrlp/4/ , lenilp/13/
        DATA lenrlr/5/ , lenilr/9/
  !
        IF ( Job.EQ.2 ) THEN
  !
           CALL DCOPY(lenrls,Rsav,1,RLS,1)
           RLS2 = Rsav(lenrls+1)
           CALL DCOPY(lenrlr,Rsav(lenrls+2),1,RLSr,1)
           CALL DCOPY(lenrlp,Rsav(lenrls+lenrlr+2),1,RLSp,1)
           DO i = 1 , lenils
              ILS(i) = Isav(i)
           ENDDO
           ILS2(1) = Isav(lenils+1)
           ILS2(2) = Isav(lenils+2)
           ILS2(3) = Isav(lenils+3)
           ILS2(4) = Isav(lenils+4)
           ioff = lenils + 2
           DO i = 1 , lenilr
              ILSr(i) = Isav(ioff+i)
           ENDDO
           ioff = ioff + lenilr
           DO i = 1 , lenilp
              ILSp(i) = Isav(ioff+i)
           ENDDO
           GOTO 99999
        ENDIF
        CALL DCOPY(lenrls,RLS,1,Rsav,1)
        Rsav(lenrls+1) = RLS2
        CALL DCOPY(lenrlr,RLSr,1,Rsav(lenrls+2),1)
        CALL DCOPY(lenrlp,RLSp,1,Rsav(lenrls+lenrlr+2),1)
        DO i = 1 , lenils
           Isav(i) = ILS(i)
        ENDDO
        Isav(lenils+1) = ILS2(1)
        Isav(lenils+2) = ILS2(2)
        Isav(lenils+3) = ILS2(3)
        Isav(lenils+4) = ILS2(4)
        ioff = lenils + 2
        DO i = 1 , lenilr
           Isav(ioff+i) = ILSr(i)
        ENDDO
        ioff = ioff + lenilr
        DO i = 1 , lenilp
           Isav(ioff+i) = ILSp(i)
        ENDDO
        RETURN
  !----------------------- End of Subroutine DSRCKR ----------------------
  99999 END
  !*==DAINVG.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DAINVG
        SUBROUTINE DAINVG(RES,ADDA,Neq,T,Y,Ydot,Miter,Ml,Mu,Pw,Ipvt,Ier)
        IMPLICIT NONE
  !*--DAINVG9138
        EXTERNAL RES , ADDA
        INTEGER Neq , Miter , Ml , Mu , Ipvt , Ier
        INTEGER i , lenpw , mlp1 , nrowpw
        DOUBLE PRECISION T , Y , Ydot , Pw
        DIMENSION Y(*) , Ydot(*) , Pw(*) , Ipvt(*)
  !-----------------------------------------------------------------------
  ! This subroutine computes the initial value
  ! of the vector YDOT satisfying
  !     A * YDOT = g(t,y)
  ! when A is nonsingular.  It is called by DLSODI for
  ! initialization only, when ISTATE = 0 .
  ! DAINVG returns an error flag IER:
  !   IER  =  0  means DAINVG was successful.
  !   IER .ge. 2 means RES returned an error flag IRES = IER.
  !   IER .lt. 0 means the a-matrix was found to be singular.
  !-----------------------------------------------------------------------
  !
        IF ( Miter.GE.4 ) THEN
  !
  ! Band matrix case -----------------------------------------------------
  !
           nrowpw = 2*Ml + Mu + 1
           lenpw = Neq*nrowpw
           DO i = 1 , lenpw
              Pw(i) = 0.0D0
           ENDDO
  !
           Ier = 1
           CALL RES(Neq,T,Y,Pw,Ydot,Ier)
           IF ( Ier.GT.1 ) RETURN
  !
           mlp1 = Ml + 1
           CALL ADDA(Neq,T,Y,Ml,Mu,Pw(mlp1),nrowpw)
           CALL DGBFA(Pw,nrowpw,Neq,Ml,Mu,Ipvt,Ier)
           IF ( Ier.EQ.0 ) THEN
              CALL DGBSL(Pw,nrowpw,Neq,Ml,Mu,Ipvt,Ydot,0)
              GOTO 99999
           ENDIF
        ELSE
  !
  ! Full matrix case -----------------------------------------------------
  !
           lenpw = Neq*Neq
           DO i = 1 , lenpw
              Pw(i) = 0.0D0
           ENDDO
  !
           Ier = 1
           CALL RES(Neq,T,Y,Pw,Ydot,Ier)
           IF ( Ier.GT.1 ) RETURN
  !
           CALL ADDA(Neq,T,Y,0,0,Pw,Neq)
           CALL DGEFA(Pw,Neq,Neq,Ipvt,Ier)
           IF ( Ier.EQ.0 ) THEN
              CALL DGESL(Pw,Neq,Neq,Ipvt,Ydot,0)
              RETURN
           ELSE
              Ier = -Ier
              RETURN
           ENDIF
        ENDIF
        Ier = -Ier
        RETURN
  !----------------------- End of Subroutine DAINVG ----------------------
  99999 END
  !*==DSTODI.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DSTODI
        SUBROUTINE DSTODI(Neq,Y,Yh,Nyh,Yh1,Ewt,Savf,Savr,Acor,Wm,Iwm,RES, &
       &                  ADDA,JAC,PJAC,SLVS)
        IMPLICIT NONE
  !*--DSTODI9209
  !*** Start of declarations inserted by SPAG
        !REAL ADDA
        !INTEGER JAC
  !*** End of declarations inserted by SPAG
        EXTERNAL RES , ADDA , JAC , PJAC , SLVS
        INTEGER Neq , Nyh , Iwm
        DOUBLE PRECISION Y , Yh , Yh1 , Ewt , Savf , Savr , Acor , Wm
        DIMENSION Neq(*) , Y(*) , Yh(Nyh,*) , Yh1(*) , Ewt(*) , Savf(*) , &
       &          Savr(*) , Acor(*) , Wm(*) , Iwm(*)
        INTEGER IOWnd , IALth , IPUp , LMAx , MEO , NQNyh , NSLp , ICF ,  &
       &        IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , LEWt ,  &
       &        LACor , LSAvf , LWM , LIWm , METh , MITer , MAXord ,      &
       &        MAXcor , MSBp , MXNcf , N , NQ , NST , NFE , NJE , NQU
        DOUBLE PRECISION CONit , CRAte , EL , ELCo , HOLd , RMAx , TESco ,&
       &                 CCMax , EL0 , H , HMIn , HMXi , HU , RC , TN ,   &
       &                 UROund
        COMMON /DLS001/ CONit , CRAte , EL(13) , ELCo(13,12) , HOLd ,     &
       &                RMAx , TESco(3,12) , CCMax , EL0 , H , HMIn ,     &
       &                HMXi , HU , RC , TN , UROund , IOWnd(6) , IALth , &
       &                IPUp , LMAx , MEO , NQNyh , NSLp , ICF , IERpj ,  &
       &                IERsl , JCUr , JSTart , KFLag , L , LYH , LEWt ,  &
       &                LACor , LSAvf , LWM , LIWm , METh , MITer ,       &
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        INTEGER i , i1 , iredo , ires , iret , j , jb , kgo , m , ncf ,   &
       &        newq
        DOUBLE PRECISION dcon , ddn , del , delp , dsm , dup , eljh ,     &
       &                 el1h , exdn , exsm , exup , r , rh , rhdn ,      &
       &                 rhsm , rhup , told
  !-----------------------------------------------------------------------
  ! DSTODI performs one step of the integration of an initial value
  ! problem for a system of Ordinary Differential Equations.
  ! Note: DSTODI is independent of the value of the iteration method
  ! indicator MITER, and hence is independent
  ! of the type of chord method used, or the Jacobian structure.
  ! Communication with DSTODI is done with the following variables:
  !
  ! NEQ    = integer array containing problem size in NEQ(1), and
  !          passed as the NEQ argument in all calls to RES, ADDA,
  !          and JAC.
  ! Y      = an array of length .ge. N used as the Y argument in
  !          all calls to RES, JAC, and ADDA.
  ! NEQ    = integer array containing problem size in NEQ(1), and
  !          passed as the NEQ argument in all calls tO RES, G, ADDA,
  !          and JAC.
  ! YH     = an NYH by LMAX array containing the dependent variables
  !          and their approximate scaled derivatives, where
  !          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
  !          j-th derivative of y(i), scaled by H**j/factorial(j)
  !          (j = 0,1,...,NQ).  On entry for the first step, the first
  !          two columns of YH must be set from the initial values.
  ! NYH    = a constant integer .ge. N, the first dimension of YH.
  ! YH1    = a one-dimensional array occupying the same space as YH.
  ! EWT    = an array of length N containing multiplicative weights
  !          for local error measurements.  Local errors in y(i) are
  !          compared to 1.0/EWT(i) in various error tests.
  ! SAVF   = an array of working storage, of length N. also used for
  !          input of YH(*,MAXORD+2) when JSTART = -1 and MAXORD is less
  !          than the current order NQ.
  !          Same as YDOTI in the driver.
  ! SAVR   = an array of working storage, of length N.
  ! ACOR   = a work array of length N used for the accumulated
  !          corrections. On a succesful return, ACOR(i) contains
  !          the estimated one-step local error in y(i).
  ! WM,IWM = real and integer work arrays associated with matrix
  !          operations in chord iteration.
  ! PJAC   = name of routine to evaluate and preprocess Jacobian matrix.
  ! SLVS   = name of routine to solve linear system in chord iteration.
  ! CCMAX  = maximum relative change in H*EL0 before PJAC is called.
  ! H      = the step size to be attempted on the next step.
  !          H is altered by the error control algorithm during the
  !          problem.  H can be either positive or negative, but its
  !          sign must remain constant throughout the problem.
  ! HMIN   = the minimum absolute value of the step size H to be used.
  ! HMXI   = inverse of the maximum absolute value of H to be used.
  !          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
  !          HMIN and HMXI may be changed at any time, but will not
  !          take effect until the next change of H is considered.
  ! TN     = the independent variable. TN is updated on each step taken.
  ! JSTART = an integer used for input only, with the following
  !          values and meanings:
  !               0  perform the first step.
  !           .gt.0  take a new step continuing from the last.
  !              -1  take the next step with a new value of H, MAXORD,
  !                    N, METH, MITER, and/or matrix parameters.
  !              -2  take the next step with a new value of H,
  !                    but with other inputs unchanged.
  !          On return, JSTART is set to 1 to facilitate continuation.
  ! KFLAG  = a completion code with the following meanings:
  !               0  the step was succesful.
  !              -1  the requested error could not be achieved.
  !              -2  corrector convergence could not be achieved.
  !              -3  RES ordered immediate return.
  !              -4  error condition from RES could not be avoided.
  !              -5  fatal error in PJAC or SLVS.
  !          A return with KFLAG = -1, -2, or -4 means either
  !          ABS(H) = HMIN or 10 consecutive failures occurred.
  !          On a return with KFLAG negative, the values of TN and
  !          the YH array are as of the beginning of the last
  !          step, and H is the last step size attempted.
  ! MAXORD = the maximum order of integration method to be allowed.
  ! MAXCOR = the maximum number of corrector iterations allowed.
  ! MSBP   = maximum number of steps between PJAC calls.
  ! MXNCF  = maximum number of convergence failures allowed.
  ! METH/MITER = the method flags.  See description in driver.
  ! N      = the number of first-order differential equations.
  !-----------------------------------------------------------------------
        KFLag = 0
        told = TN
        ncf = 0
        IERpj = 0
        IERsl = 0
        JCUr = 0
        ICF = 0
        delp = 0.0D0
        IF ( JSTart.GT.0 ) GOTO 500
        IF ( JSTart.EQ.-1 ) THEN
  !-----------------------------------------------------------------------
  ! The following block handles preliminaries needed when JSTART = -1.
  ! IPUP is set to MITER to force a matrix update.
  ! If an order increase is about to be considered (IALTH = 1),
  ! IALTH is reset to 2 to postpone consideration one more step.
  ! If the caller has changed METH, DCFODE is called to reset
  ! the coefficients of the method.
  ! If the caller has changed MAXORD to a value less than the current
  ! order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
  ! If H is to be changed, YH must be rescaled.
  ! If H or METH is being changed, IALTH is reset to L = NQ + 1
  ! to prevent further changes in H for that many steps.
  !-----------------------------------------------------------------------
           IPUp = MITer
           LMAx = MAXord + 1
           IF ( IALth.EQ.1 ) IALth = 2
           IF ( METh.NE.MEO ) THEN
              CALL DCFODE(METh,ELCo,TESco)
              MEO = METh
              IF ( NQ.LE.MAXord ) THEN
                 IALth = L
                 iret = 1
                 GOTO 100
              ENDIF
           ELSEIF ( NQ.LE.MAXord ) THEN
              GOTO 200
           ENDIF
           NQ = MAXord
           L = LMAx
           DO i = 1 , L
              EL(i) = ELCo(i,NQ)
           ENDDO
           NQNyh = NQ*Nyh
           RC = RC*EL(1)/EL0
           EL0 = EL(1)
           CONit = 0.5D0/(NQ+2)
           ddn = DVNORM(N,Savf,Ewt)/TESco(1,L)
           exdn = 1.0D0/L
           rhdn = 1.0D0/(1.3D0*ddn**exdn+0.0000013D0)
           rh = MIN(rhdn,1.0D0)
           iredo = 3
           IF ( H.EQ.HOLd ) GOTO 300
           rh = MIN(rh,ABS(H/HOLd))
           H = HOLd
           GOTO 400
        ELSE
           IF ( JSTart.EQ.-2 ) GOTO 200
  !-----------------------------------------------------------------------
  ! On the first call, the order is set to 1, and other variables are
  ! initialized.  RMAX is the maximum ratio by which H can be increased
  ! in a single step.  It is initially 1.E4 to compensate for the small
  ! initial H, but then is normally equal to 10.  If a failure
  ! occurs (in corrector convergence or error test), RMAX is set at 2
  ! for the next increase.
  !-----------------------------------------------------------------------
           LMAx = MAXord + 1
           NQ = 1
           L = 2
           IALth = 2
           RMAx = 10000.0D0
           RC = 0.0D0
           EL0 = 1.0D0
           CRAte = 0.7D0
           HOLd = H
           MEO = METh
           NSLp = 0
           IPUp = MITer
           iret = 3
  !-----------------------------------------------------------------------
  ! DCFODE is called to get all the integration coefficients for the
  ! current METH.  Then the EL vector and related constants are reset
  ! whenever the order NQ is changed, or at the start of the problem.
  !-----------------------------------------------------------------------
           CALL DCFODE(METh,ELCo,TESco)
        ENDIF
   100  DO i = 1 , L
           EL(i) = ELCo(i,NQ)
        ENDDO
        NQNyh = NQ*Nyh
        RC = RC*EL(1)/EL0
        EL0 = EL(1)
        CONit = 0.5D0/(NQ+2)
        GOTO (200,300,500) , iret
  !-----------------------------------------------------------------------
  ! If H is being changed, the H ratio RH is checked against
  ! RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
  ! L = NQ + 1 to prevent a change of H for that many steps, unless
  ! forced by a convergence or error test failure.
  !-----------------------------------------------------------------------
   200  IF ( H.EQ.HOLd ) GOTO 500
        rh = H/HOLd
        H = HOLd
        iredo = 3
        GOTO 400
   300  rh = MAX(rh,HMIn/ABS(H))
   400  rh = MIN(rh,RMAx)
        rh = rh/MAX(1.0D0,ABS(H)*HMXi*rh)
        r = 1.0D0
        DO j = 2 , L
           r = r*rh
           DO i = 1 , N
              Yh(i,j) = Yh(i,j)*r
           ENDDO
        ENDDO
        H = H*rh
        RC = RC*rh
        IALth = L
        IF ( iredo.EQ.0 ) THEN
           RMAx = 10.0D0
           GOTO 1500
        ENDIF
  !-----------------------------------------------------------------------
  ! This section computes the predicted values by effectively
  ! multiplying the YH array by the Pascal triangle matrix.
  ! RC is the ratio of new to old values of the coefficient  H*EL(1).
  ! When RC differs from 1 by more than CCMAX, IPUP is set to MITER
  ! to force PJAC to be called.
  ! In any case, PJAC is called at least every MSBP steps.
  !-----------------------------------------------------------------------
   500  IF ( ABS(RC-1.0D0).GT.CCMax ) IPUp = MITer
        IF ( NST.GE.NSLp+MSBp ) IPUp = MITer
        TN = TN + H
        i1 = NQNyh + 1
        DO jb = 1 , NQ
           i1 = i1 - Nyh
  !DIR$ IVDEP
           DO i = i1 , NQNyh
              Yh1(i) = Yh1(i) + Yh1(i+Nyh)
           ENDDO
        ENDDO
  !-----------------------------------------------------------------------
  ! Up to MAXCOR corrector iterations are taken.  A convergence test is
  ! made on the RMS-norm of each correction, weighted by H and the
  ! error weight vector EWT.  The sum of the corrections is accumulated
  ! in ACOR(i).  The YH array is not altered in the corrector loop.
  !-----------------------------------------------------------------------
   600  m = 0
        DO i = 1 , N
           Savf(i) = Yh(i,2)/H
           Y(i) = Yh(i,1)
        ENDDO
        IF ( IPUp.GT.0 ) THEN
  !-----------------------------------------------------------------------
  ! If indicated, the matrix P = A - H*EL(1)*dr/dy is reevaluated and
  ! preprocessed before starting the corrector iteration.  IPUP is set
  ! to 0 as an indicator that this has been done.
  !-----------------------------------------------------------------------
           CALL PJAC(Neq,Y,Yh,Nyh,Ewt,Acor,Savr,Savf,Wm,Iwm,RES,JAC,ADDA)
           IPUp = 0
           RC = 1.0D0
           NSLp = NST
           CRAte = 0.7D0
           IF ( IERpj.EQ.0 ) GOTO 700
           IF ( IERpj.LT.0 ) GOTO 1100
           ires = IERpj
           GOTO (1000,1100,1000) , ires
        ENDIF
  ! Get residual at predicted values, if not already done in PJAC. -------
        ires = 1
        CALL RES(Neq,TN,Y,Savf,Savr,ires)
        NFE = NFE + 1
        kgo = ABS(ires)
        GOTO (700,1100,1000) , kgo
   700  DO i = 1 , N
           Acor(i) = 0.0D0
        ENDDO
  !-----------------------------------------------------------------------
  ! Solve the linear system with the current residual as
  ! right-hand side and P as coefficient matrix.
  !-----------------------------------------------------------------------
   800  CALL SLVS(Wm,Iwm,Savr,Savf)
        IF ( IERsl.LT.0 ) GOTO 1000
        IF ( IERsl.LE.0 ) THEN
           el1h = EL(1)*H
           del = DVNORM(N,Savr,Ewt)*ABS(H)
           DO i = 1 , N
              Acor(i) = Acor(i) + Savr(i)
              Savf(i) = Acor(i) + Yh(i,2)/H
              Y(i) = Yh(i,1) + el1h*Acor(i)
           ENDDO
  !-----------------------------------------------------------------------
  ! Test for convergence.  If M .gt. 0, an estimate of the convergence
  ! rate constant is stored in CRATE, and this is used in the test.
  !-----------------------------------------------------------------------
           IF ( m.NE.0 ) CRAte = MAX(0.2D0*CRAte,del/delp)
           dcon = del*MIN(1.0D0,1.5D0*CRAte)/(TESco(2,NQ)*CONit)
           IF ( dcon.LE.1.0D0 ) THEN
  !-----------------------------------------------------------------------
  ! The corrector has converged.  JCUR is set to 0
  ! to signal that the Jacobian involved may need updating later.
  ! The local error test is made and control passes to statement 500
  ! if it fails.
  !-----------------------------------------------------------------------
              JCUr = 0
              IF ( m.EQ.0 ) dsm = del/TESco(2,NQ)
              IF ( m.GT.0 ) dsm = ABS(H)*DVNORM(N,Acor,Ewt)/TESco(2,NQ)
              IF ( dsm.GT.1.0D0 ) THEN
  !-----------------------------------------------------------------------
  ! The error test failed.  KFLAG keeps track of multiple failures.
  ! restore TN and the YH array to their previous values, and prepare
  ! to try the step again.  Compute the optimum step size for this or
  ! one lower order.  After 2 or more failures, H is forced to decrease
  ! by a factor of 0.1 or less.
  !-----------------------------------------------------------------------
                 KFLag = KFLag - 1
                 TN = told
                 i1 = NQNyh + 1
                 DO jb = 1 , NQ
                    i1 = i1 - Nyh
  !DIR$ IVDEP
                    DO i = i1 , NQNyh
                       Yh1(i) = Yh1(i) - Yh1(i+Nyh)
                    ENDDO
                 ENDDO
                 RMAx = 2.0D0
                 IF ( ABS(H).LE.HMIn*1.00001D0 ) THEN
  !-----------------------------------------------------------------------
  ! All returns are made through this section.  H is saved in HOLD
  ! to allow the caller to change H on the next step.
  !-----------------------------------------------------------------------
                    KFLag = -1
                    GOTO 1600
                 ELSEIF ( KFLag.LE.-7 ) THEN
                    KFLag = -1
                    GOTO 1600
                 ELSE
                    iredo = 2
                    rhup = 0.0D0
                    GOTO 1200
                 ENDIF
              ELSE
  !-----------------------------------------------------------------------
  ! After a successful step, update the YH array.
  ! Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
  ! If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
  ! use in a possible order increase on the next step.
  ! If a change in H is considered, an increase or decrease in order
  ! by one is considered also.  A change in H is made only if it is by a
  ! factor of at least 1.1.  If not, IALTH is set to 3 to prevent
  ! testing for that many steps.
  !-----------------------------------------------------------------------
                 KFLag = 0
                 iredo = 0
                 NST = NST + 1
                 HU = H
                 NQU = NQ
                 DO j = 1 , L
                    eljh = EL(j)*H
                    DO i = 1 , N
                       Yh(i,j) = Yh(i,j) + eljh*Acor(i)
                    ENDDO
                 ENDDO
                 IALth = IALth - 1
                 IF ( IALth.EQ.0 ) THEN
  !-----------------------------------------------------------------------
  ! Regardless of the success or failure of the step, factors
  ! RHDN, RHSM, and RHUP are computed, by which H could be multiplied
  ! at order NQ - 1, order NQ, or order NQ + 1, respectively.
  ! In the case of failure, RHUP = 0.0 to avoid an order increase.
  ! The largest of these is determined and the new order chosen
  ! accordingly.  If the order is to be increased, we compute one
  ! additional scaled derivative.
  !-----------------------------------------------------------------------
                    rhup = 0.0D0
                    IF ( L.NE.LMAx ) THEN
                       DO i = 1 , N
                          Savf(i) = Acor(i) - Yh(i,LMAx)
                       ENDDO
                       dup = ABS(H)*DVNORM(N,Savf,Ewt)/TESco(3,NQ)
                       exup = 1.0D0/(L+1)
                       rhup = 1.0D0/(1.4D0*dup**exup+0.0000014D0)
                    ENDIF
                    GOTO 1200
                 ELSE
                    IF ( IALth.LE.1 ) THEN
                       IF ( L.NE.LMAx ) THEN
                          DO i = 1 , N
                             Yh(i,LMAx) = Acor(i)
                          ENDDO
                       ENDIF
                    ENDIF
                    GOTO 1500
                 ENDIF
              ENDIF
           ELSE
              m = m + 1
              IF ( m.NE.MAXcor ) THEN
                 IF ( m.LT.2 .OR. del.LE.2.0D0*delp ) THEN
                    delp = del
                    ires = 1
                    CALL RES(Neq,TN,Y,Savf,Savr,ires)
                    NFE = NFE + 1
                    kgo = ABS(ires)
                    GOTO (800,1100,900) , kgo
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
  !-----------------------------------------------------------------------
  ! The correctors failed to converge, or RES has returned abnormally.
  ! on a convergence failure, if the Jacobian is out of date, PJAC is
  ! called for the next try.  Otherwise the YH array is retracted to its
  ! values before prediction, and H is reduced, if possible.
  ! take an error exit if IRES = 2, or H cannot be reduced, or MXNCF
  ! failures have occurred, or a fatal error occurred in PJAC or SLVS.
  !-----------------------------------------------------------------------
   900  ICF = 1
        IF ( JCUr.NE.1 ) THEN
           IPUp = MITer
           GOTO 600
        ENDIF
   1000 ICF = 2
        ncf = ncf + 1
        RMAx = 2.0D0
   1100 TN = told
        i1 = NQNyh + 1
        DO jb = 1 , NQ
           i1 = i1 - Nyh
  !DIR$ IVDEP
           DO i = i1 , NQNyh
              Yh1(i) = Yh1(i) - Yh1(i+Nyh)
           ENDDO
        ENDDO
        IF ( ires.EQ.2 ) THEN
           KFLag = -1 - ires
        ELSEIF ( IERpj.LT.0 .OR. IERsl.LT.0 ) THEN
           KFLag = -5
        ELSE
           IF ( ABS(H).GT.HMIn*1.00001D0 ) THEN
              IF ( ncf.NE.MXNcf ) THEN
                 rh = 0.25D0
                 IPUp = MITer
                 iredo = 1
                 GOTO 300
              ENDIF
           ENDIF
           IF ( ires.EQ.3 ) THEN
              KFLag = -1 - ires
           ELSE
              KFLag = -2
           ENDIF
        ENDIF
        GOTO 1600
   1200 exsm = 1.0D0/L
        rhsm = 1.0D0/(1.2D0*dsm**exsm+0.0000012D0)
        rhdn = 0.0D0
        IF ( NQ.NE.1 ) THEN
           ddn = DVNORM(N,Yh(1,L),Ewt)/TESco(1,NQ)
           exdn = 1.0D0/NQ
           rhdn = 1.0D0/(1.3D0*ddn**exdn+0.0000013D0)
        ENDIF
        IF ( rhsm.GE.rhup ) THEN
           IF ( rhsm.GE.rhdn ) THEN
              newq = NQ
              rh = rhsm
              GOTO 1300
           ENDIF
        ELSEIF ( rhup.GT.rhdn ) THEN
           newq = L
           rh = rhup
           IF ( rh.LT.1.1D0 ) THEN
              IALth = 3
              GOTO 1500
           ELSE
              r = H*EL(L)/L
              DO i = 1 , N
                 Yh(i,newq+1) = Acor(i)*r
              ENDDO
              GOTO 1400
           ENDIF
        ENDIF
        newq = NQ - 1
        rh = rhdn
        IF ( KFLag.LT.0 .AND. rh.GT.1.0D0 ) rh = 1.0D0
   1300 IF ( (KFLag.EQ.0) .AND. (rh.LT.1.1D0) ) THEN
           IALth = 3
           GOTO 1500
        ELSE
           IF ( KFLag.LE.-2 ) rh = MIN(rh,0.1D0)
  !-----------------------------------------------------------------------
  ! If there is a change of order, reset NQ, L, and the coefficients.
  ! In any case H is reset according to RH and the YH array is rescaled.
  ! Then exit from 690 if the step was OK, or redo the step otherwise.
  !-----------------------------------------------------------------------
           IF ( newq.EQ.NQ ) GOTO 300
        ENDIF
   1400 NQ = newq
        L = NQ + 1
        iret = 2
        GOTO 100
   1500 r = H/TESco(2,NQU)
        DO i = 1 , N
           Acor(i) = Acor(i)*r
        ENDDO
   1600 HOLd = H
        JSTart = 1
  !----------------------- End of Subroutine DSTODI ----------------------
        END
  !*==DPREPJI.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DPREPJI
        SUBROUTINE DPREPJI(Neq,Y,Yh,Nyh,Ewt,Rtem,Savr,S,Wm,Iwm,RES,JAC,   &
       &                   ADDA)
        IMPLICIT NONE
  !*--DPREPJI9730
        EXTERNAL RES , JAC , ADDA
        INTEGER Neq , Nyh , Iwm
        DOUBLE PRECISION Y , Yh , Ewt , Rtem , Savr , S , Wm
        DIMENSION Neq(*) , Y(*) , Yh(Nyh,*) , Ewt(*) , Rtem(*) , S(*) ,   &
       &          Savr(*) , Wm(*) , Iwm(*)
        INTEGER IOWnd , IOWns , ICF , IERpj , IERsl , JCUr , JSTart ,     &
       &        KFLag , L , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,     &
       &        METh , MITer , MAXord , MAXcor , MSBp , MXNcf , N , NQ ,  &
       &        NST , NFE , NJE , NQU
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , IOWnd(6) , IOWns(6) , ICF ,    &
       &                IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , &
       &                LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,&
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        INTEGER i , i1 , i2 , ier , ii , ires , j , j1 , jj , lenp , mba ,&
       &        mband , meb1 , meband , ml , ml3 , mu
        DOUBLE PRECISION con , fac , hl0 , r , srur , yi , yj , yjj
  !-----------------------------------------------------------------------
  ! DPREPJI is called by DSTODI to compute and process the matrix
  ! P = A - H*EL(1)*J , where J is an approximation to the Jacobian dr/dy,
  ! where r = g(t,y) - A(t,y)*s.  Here J is computed by the user-supplied
  ! routine JAC if MITER = 1 or 4, or by finite differencing if MITER =
  ! 2 or 5.  J is stored in WM, rescaled, and ADDA is called to generate
  ! P. P is then subjected to LU decomposition in preparation
  ! for later solution of linear systems with P as coefficient
  ! matrix.  This is done by DGEFA if MITER = 1 or 2, and by
  ! DGBFA if MITER = 4 or 5.
  !
  ! In addition to variables described previously, communication
  ! with DPREPJI uses the following:
  ! Y     = array containing predicted values on entry.
  ! RTEM  = work array of length N (ACOR in DSTODI).
  ! SAVR  = array used for output only.  On output it contains the
  !         residual evaluated at current values of t and y.
  ! S     = array containing predicted values of dy/dt (SAVF in DSTODI).
  ! WM    = real work space for matrices.  On output it contains the
  !         LU decomposition of P.
  !         Storage of matrix elements starts at WM(3).
  !         WM also contains the following matrix-related data:
  !         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
  ! IWM   = integer work space containing pivot information, starting at
  !         IWM(21).  IWM also contains the band parameters
  !         ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
  ! EL0   = el(1) (input).
  ! IERPJ = output error flag.
  !         = 0 if no trouble occurred,
  !         = 1 if the P matrix was found to be singular,
  !         = IRES (= 2 or 3) if RES returned IRES = 2 or 3.
  ! JCUR  = output flag = 1 to indicate that the Jacobian matrix
  !         (or approximation) is now current.
  ! This routine also uses the Common variables EL0, H, TN, UROUND,
  ! MITER, N, NFE, and NJE.
  !-----------------------------------------------------------------------
        NJE = NJE + 1
        hl0 = H*EL0
        IERpj = 0
        JCUr = 1
        GOTO (100,200,400,500,600) , MITer
  ! If MITER = 1, call RES, then JAC, and multiply by scalar. ------------
   100  ires = 1
        CALL RES(Neq,TN,Y,S,Savr,ires)
        NFE = NFE + 1
        IF ( ires.GT.1 ) GOTO 800
        lenp = N*N
        DO i = 1 , lenp
           Wm(i+2) = 0.0D0
        ENDDO
        CALL JAC(Neq,TN,Y,S,0,0,Wm(3),N)
        con = -hl0
        DO i = 1 , lenp
           Wm(i+2) = Wm(i+2)*con
        ENDDO
        GOTO 300
  ! If MITER = 2, make N + 1 calls to RES to approximate J. --------------
   200  ires = -1
        CALL RES(Neq,TN,Y,S,Savr,ires)
        NFE = NFE + 1
        IF ( ires.GT.1 ) GOTO 800
        srur = Wm(1)
        j1 = 2
        DO j = 1 , N
           yj = Y(j)
           r = MAX(srur*ABS(yj),0.01D0/Ewt(j))
           Y(j) = Y(j) + r
           fac = -hl0/r
           CALL RES(Neq,TN,Y,S,Rtem,ires)
           NFE = NFE + 1
           IF ( ires.GT.1 ) GOTO 800
           DO i = 1 , N
              Wm(i+j1) = (Rtem(i)-Savr(i))*fac
           ENDDO
           Y(j) = yj
           j1 = j1 + N
        ENDDO
        ires = 1
        CALL RES(Neq,TN,Y,S,Savr,ires)
        NFE = NFE + 1
        IF ( ires.GT.1 ) GOTO 800
  ! Add matrix A. --------------------------------------------------------
   300  CALL ADDA(Neq,TN,Y,0,0,Wm(3),N)
  ! Do LU decomposition on P. --------------------------------------------
        CALL DGEFA(Wm(3),N,N,Iwm(21),ier)
        IF ( ier.NE.0 ) IERpj = 1
        RETURN
  ! Dummy section for MITER = 3
   400  RETURN
  ! If MITER = 4, call RES, then JAC, and multiply by scalar. ------------
   500  ires = 1
        CALL RES(Neq,TN,Y,S,Savr,ires)
        NFE = NFE + 1
        IF ( ires.GT.1 ) GOTO 800
        ml = Iwm(1)
        mu = Iwm(2)
        ml3 = ml + 3
        mband = ml + mu + 1
        meband = mband + ml
        lenp = meband*N
        DO i = 1 , lenp
           Wm(i+2) = 0.0D0
        ENDDO
        CALL JAC(Neq,TN,Y,S,ml,mu,Wm(ml3),meband)
        con = -hl0
        DO i = 1 , lenp
           Wm(i+2) = Wm(i+2)*con
        ENDDO
        GOTO 700
  ! If MITER = 5, make ML + MU + 2 calls to RES to approximate J. --------
   600  ires = -1
        CALL RES(Neq,TN,Y,S,Savr,ires)
        NFE = NFE + 1
        IF ( ires.GT.1 ) GOTO 800
        ml = Iwm(1)
        mu = Iwm(2)
        ml3 = ml + 3
        mband = ml + mu + 1
        mba = MIN(mband,N)
        meband = mband + ml
        meb1 = meband - 1
        srur = Wm(1)
        DO j = 1 , mba
           DO i = j , N , mband
              yi = Y(i)
              r = MAX(srur*ABS(yi),0.01D0/Ewt(i))
              Y(i) = Y(i) + r
           ENDDO
           CALL RES(Neq,TN,Y,S,Rtem,ires)
           NFE = NFE + 1
           IF ( ires.GT.1 ) GOTO 800
           DO jj = j , N , mband
              Y(jj) = Yh(jj,1)
              yjj = Y(jj)
              r = MAX(srur*ABS(yjj),0.01D0/Ewt(jj))
              fac = -hl0/r
              i1 = MAX(jj-mu,1)
              i2 = MIN(jj+ml,N)
              ii = jj*meb1 - ml + 2
              DO i = i1 , i2
                 Wm(ii+i) = (Rtem(i)-Savr(i))*fac
              ENDDO
           ENDDO
        ENDDO
        ires = 1
        CALL RES(Neq,TN,Y,S,Savr,ires)
        NFE = NFE + 1
        IF ( ires.GT.1 ) GOTO 800
  ! Add matrix A. --------------------------------------------------------
   700  CALL ADDA(Neq,TN,Y,ml,mu,Wm(ml3),meband)
  ! Do LU decomposition of P. --------------------------------------------
        CALL DGBFA(Wm(3),meband,N,ml,mu,Iwm(21),ier)
        IF ( ier.NE.0 ) IERpj = 1
        RETURN
  ! Error return for IRES = 2 or IRES = 3 return from RES. ---------------
   800  IERpj = ires
  !----------------------- End of Subroutine DPREPJI ---------------------
        END
  !*==DAIGBT.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DAIGBT
        SUBROUTINE DAIGBT(RES,ADDA,Neq,T,Y,Ydot,Mb,Nb,Pw,Ipvt,Ier)
        IMPLICIT NONE
  !*--DAIGBT9913
        EXTERNAL RES , ADDA
        INTEGER Neq , Mb , Nb , Ipvt , Ier
        INTEGER i , lenpw , lblox , lpb , lpc
        DOUBLE PRECISION T , Y , Ydot , Pw
        DIMENSION Y(*) , Ydot(*) , Pw(*) , Ipvt(*) , Neq(*)
  !-----------------------------------------------------------------------
  ! This subroutine computes the initial value
  ! of the vector YDOT satisfying
  !     A * YDOT = g(t,y)
  ! when A is nonsingular.  It is called by DLSOIBT for
  ! initialization only, when ISTATE = 0 .
  ! DAIGBT returns an error flag IER:
  !   IER  =  0  means DAIGBT was successful.
  !   IER .ge. 2 means RES returned an error flag IRES = IER.
  !   IER .lt. 0 means the A matrix was found to have a singular
  !              diagonal block (hence YDOT could not be solved for).
  !-----------------------------------------------------------------------
        lblox = Mb*Mb*Nb
        lpb = 1 + lblox
        lpc = lpb + lblox
        lenpw = 3*lblox
        DO i = 1 , lenpw
           Pw(i) = 0.0D0
        ENDDO
        Ier = 1
        CALL RES(Neq,T,Y,Pw,Ydot,Ier)
        IF ( Ier.GT.1 ) RETURN
        CALL ADDA(Neq,T,Y,Mb,Nb,Pw(1),Pw(lpb),Pw(lpc))
        CALL DDECBT(Mb,Nb,Pw,Pw(lpb),Pw(lpc),Ipvt,Ier)
        IF ( Ier.EQ.0 ) THEN
           CALL DSOLBT(Mb,Nb,Pw,Pw(lpb),Pw(lpc),Ydot,Ipvt)
           GOTO 99999
        ENDIF
        Ier = -Ier
        RETURN
  !----------------------- End of Subroutine DAIGBT ----------------------
  99999 END
  !*==DPJIBT.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DPJIBT
        SUBROUTINE DPJIBT(Neq,Y,Yh,Nyh,Ewt,Rtem,Savr,S,Wm,Iwm,RES,JAC,    &
       &                  ADDA)
        IMPLICIT NONE
  !*--DPJIBT9956
        EXTERNAL RES , JAC , ADDA
        INTEGER Neq , Nyh , Iwm
        DOUBLE PRECISION Y , Yh , Ewt , Rtem , Savr , S , Wm
        DIMENSION Neq(*) , Y(*) , Yh(Nyh,*) , Ewt(*) , Rtem(*) , S(*) ,   &
       &          Savr(*) , Wm(*) , Iwm(*)
        INTEGER IOWnd , IOWns , ICF , IERpj , IERsl , JCUr , JSTart ,     &
       &        KFLag , L , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,     &
       &        METh , MITer , MAXord , MAXcor , MSBp , MXNcf , N , NQ ,  &
       &        NST , NFE , NJE , NQU
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , IOWnd(6) , IOWns(6) , ICF ,    &
       &                IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , &
       &                LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,&
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        INTEGER i , ier , iia , iib , iic , ipa , ipb , ipc , ires , j ,  &
       &        j1 , j2 , k , k1 , lenp , lblox , lpb , lpc , mb , mbsq , &
       &        mwid , nb
        DOUBLE PRECISION con , fac , hl0 , r , srur
  !-----------------------------------------------------------------------
  ! DPJIBT is called by DSTODI to compute and process the matrix
  ! P = A - H*EL(1)*J , where J is an approximation to the Jacobian dr/dy,
  ! and r = g(t,y) - A(t,y)*s.  Here J is computed by the user-supplied
  ! routine JAC if MITER = 1, or by finite differencing if MITER = 2.
  ! J is stored in WM, rescaled, and ADDA is called to generate P.
  ! P is then subjected to LU decomposition by DDECBT in preparation
  ! for later solution of linear systems with P as coefficient matrix.
  !
  ! In addition to variables described previously, communication
  ! with DPJIBT uses the following:
  ! Y     = array containing predicted values on entry.
  ! RTEM  = work array of length N (ACOR in DSTODI).
  ! SAVR  = array used for output only.  On output it contains the
  !         residual evaluated at current values of t and y.
  ! S     = array containing predicted values of dy/dt (SAVF in DSTODI).
  ! WM    = real work space for matrices.  On output it contains the
  !         LU decomposition of P.
  !         Storage of matrix elements starts at WM(3).
  !         WM also contains the following matrix-related data:
  !         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
  ! IWM   = integer work space containing pivot information, starting at
  !         IWM(21).  IWM also contains block structure parameters
  !         MB = IWM(1) and NB = IWM(2).
  ! EL0   = EL(1) (input).
  ! IERPJ = output error flag.
  !         = 0 if no trouble occurred,
  !         = 1 if the P matrix was found to be unfactorable,
  !         = IRES (= 2 or 3) if RES returned IRES = 2 or 3.
  ! JCUR  = output flag = 1 to indicate that the Jacobian matrix
  !         (or approximation) is now current.
  ! This routine also uses the Common variables EL0, H, TN, UROUND,
  ! MITER, N, NFE, and NJE.
  !-----------------------------------------------------------------------
        NJE = NJE + 1
        hl0 = H*EL0
        IERpj = 0
        JCUr = 1
        mb = Iwm(1)
        nb = Iwm(2)
        mbsq = mb*mb
        lblox = mbsq*nb
        lpb = 3 + lblox
        lpc = lpb + lblox
        lenp = 3*lblox
        GOTO (100,200) , MITer
  ! If MITER = 1, call RES, then JAC, and multiply by scalar. ------------
   100  ires = 1
        CALL RES(Neq,TN,Y,S,Savr,ires)
        NFE = NFE + 1
        IF ( ires.GT.1 ) GOTO 400
        DO i = 1 , lenp
           Wm(i+2) = 0.0D0
        ENDDO
        CALL JAC(Neq,TN,Y,S,mb,nb,Wm(3),Wm(lpb),Wm(lpc))
        con = -hl0
        DO i = 1 , lenp
           Wm(i+2) = Wm(i+2)*con
        ENDDO
        GOTO 300
  !
  ! If MITER = 2, make 3*MB + 1 calls to RES to approximate J. -----------
   200  ires = -1
        CALL RES(Neq,TN,Y,S,Savr,ires)
        NFE = NFE + 1
        IF ( ires.GT.1 ) GOTO 400
        mwid = 3*mb
        srur = Wm(1)
        DO i = 1 , lenp
           Wm(2+i) = 0.0D0
        ENDDO
        DO k = 1 , 3
           DO j = 1 , mb
  !         Increment Y(I) for group of column indices, and call RES. ----
              j1 = j + (k-1)*mb
              DO i = j1 , N , mwid
                 r = MAX(srur*ABS(Y(i)),0.01D0/Ewt(i))
                 Y(i) = Y(i) + r
              ENDDO
              CALL RES(Neq,TN,Y,S,Rtem,ires)
              NFE = NFE + 1
              IF ( ires.GT.1 ) GOTO 400
              DO i = 1 , N
                 Rtem(i) = Rtem(i) - Savr(i)
              ENDDO
              k1 = k
              DO i = j1 , N , mwid
  !           Get Jacobian elements in column I (block-column K1). -------
                 Y(i) = Yh(i,1)
                 r = MAX(srur*ABS(Y(i)),0.01D0/Ewt(i))
                 fac = -hl0/r
  !           Compute and load elements PA(*,J,K1). ----------------------
                 iia = i - j
                 ipa = 2 + (j-1)*mb + (k1-1)*mbsq
                 DO j2 = 1 , mb
                    Wm(ipa+j2) = Rtem(iia+j2)*fac
                 ENDDO
                 IF ( k1.GT.1 ) THEN
  !           Compute and load elements PB(*,J,K1-1). --------------------
                    iib = iia - mb
                    ipb = ipa + lblox - mbsq
                    DO j2 = 1 , mb
                       Wm(ipb+j2) = Rtem(iib+j2)*fac
                    ENDDO
                 ENDIF
                 IF ( k1.LT.nb ) THEN
  !           Compute and load elements PC(*,J,K1+1). --------------------
                    iic = iia + mb
                    ipc = ipa + 2*lblox + mbsq
                    DO j2 = 1 , mb
                       Wm(ipc+j2) = Rtem(iic+j2)*fac
                    ENDDO
                 ENDIF
                 IF ( k1.EQ.3 ) THEN
  !           Compute and load elements PC(*,J,1). -----------------------
                    ipc = ipa - 2*mbsq + 2*lblox
                    DO j2 = 1 , mb
                       Wm(ipc+j2) = Rtem(j2)*fac
                    ENDDO
                 ENDIF
                 IF ( k1.EQ.nb-2 ) THEN
  !           Compute and load elements PB(*,J,NB). ----------------------
                    iib = N - mb
                    ipb = ipa + 2*mbsq + lblox
                    DO j2 = 1 , mb
                       Wm(ipb+j2) = Rtem(iib+j2)*fac
                    ENDDO
                 ENDIF
                 k1 = k1 + 3
              ENDDO
           ENDDO
        ENDDO
  ! RES call for first corrector iteration. ------------------------------
        ires = 1
        CALL RES(Neq,TN,Y,S,Savr,ires)
        NFE = NFE + 1
        IF ( ires.GT.1 ) GOTO 400
  ! Add matrix A. --------------------------------------------------------
   300  CALL ADDA(Neq,TN,Y,mb,nb,Wm(3),Wm(lpb),Wm(lpc))
  ! Do LU decomposition on P. --------------------------------------------
        CALL DDECBT(mb,nb,Wm(3),Wm(lpb),Wm(lpc),Iwm(21),ier)
        IF ( ier.NE.0 ) IERpj = 1
        RETURN
  ! Error return for IRES = 2 or IRES = 3 return from RES. ---------------
   400  IERpj = ires
  !----------------------- End of Subroutine DPJIBT ----------------------
        END
  !*==DSLSBT.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DSLSBT
        SUBROUTINE DSLSBT(Wm,Iwm,X,Tem)
        IMPLICIT NONE
  !*--DSLSBT10129
        INTEGER Iwm
        INTEGER lblox , lpb , lpc , mb , nb
        DOUBLE PRECISION Wm , X , Tem
        DIMENSION Wm(*) , Iwm(*) , X(*) , Tem(*)
  !-----------------------------------------------------------------------
  ! This routine acts as an interface between the core integrator
  ! routine and the DSOLBT routine for the solution of the linear system
  ! arising from chord iteration.
  ! Communication with DSLSBT uses the following variables:
  ! WM    = real work space containing the LU decomposition,
  !         starting at WM(3).
  ! IWM   = integer work space containing pivot information, starting at
  !         IWM(21).  IWM also contains block structure parameters
  !         MB = IWM(1) and NB = IWM(2).
  ! X     = the right-hand side vector on input, and the solution vector
  !         on output, of length N.
  ! TEM   = vector of work space of length N, not used in this version.
  !-----------------------------------------------------------------------
        mb = Iwm(1)
        nb = Iwm(2)
        lblox = mb*mb*nb
        lpb = 3 + lblox
        lpc = lpb + lblox
        CALL DSOLBT(mb,nb,Wm(3),Wm(lpb),Wm(lpc),X,Iwm(21))
  !----------------------- End of Subroutine DSLSBT ----------------------
        END
  !*==DDECBT.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DDECBT
        SUBROUTINE DDECBT(M,N,A,B,C,Ip,Ier)
        IMPLICIT NONE
  !*--DDECBT10160
        INTEGER M , N , Ip(M,N) , Ier
        DOUBLE PRECISION A(M,M,N) , B(M,M,N) , C(M,M,N)
  !-----------------------------------------------------------------------
  ! Block-tridiagonal matrix decomposition routine.
  ! Written by A. C. Hindmarsh.
  ! Latest revision:  November 10, 1983 (ACH)
  ! Reference:  UCID-30150
  !             Solution of Block-Tridiagonal Systems of Linear
  !             Algebraic Equations
  !             A.C. Hindmarsh
  !             February 1977
  ! The input matrix contains three blocks of elements in each block-row,
  ! including blocks in the (1,3) and (N,N-2) block positions.
  ! DDECBT uses block Gauss elimination and Subroutines DGEFA and DGESL
  ! for solution of blocks.  Partial pivoting is done within
  ! block-rows only.
  !
  ! Note: this version uses LINPACK routines DGEFA/DGESL instead of
  ! of dec/sol for solution of blocks, and it uses the BLAS routine DDOT
  ! for dot product calculations.
  !
  ! Input:
  !     M = order of each block.
  !     N = number of blocks in each direction of the matrix.
  !         N must be 4 or more.  The complete matrix has order M*N.
  !     A = M by M by N array containing diagonal blocks.
  !         A(i,j,k) contains the (i,j) element of the k-th block.
  !     B = M by M by N array containing the super-diagonal blocks
  !         (in B(*,*,k) for k = 1,...,N-1) and the block in the (N,N-2)
  !         block position (in B(*,*,N)).
  !     C = M by M by N array containing the subdiagonal blocks
  !         (in C(*,*,k) for k = 2,3,...,N) and the block in the
  !         (1,3) block position (in C(*,*,1)).
  !    IP = integer array of length M*N for working storage.
  ! Output:
  ! A,B,C = M by M by N arrays containing the block-LU decomposition
  !         of the input matrix.
  !    IP = M by N array of pivot information.  IP(*,k) contains
  !         information for the k-th digonal block.
  !   IER = 0  if no trouble occurred, or
  !       = -1 if the input value of M or N was illegal, or
  !       = k  if a singular matrix was found in the k-th diagonal block.
  ! Use DSOLBT to solve the associated linear system.
  !
  ! External routines required: DGEFA and DGESL (from LINPACK) and
  ! DDOT (from the BLAS, or Basic Linear Algebra package).
  !-----------------------------------------------------------------------
        INTEGER nm1 , nm2 , km1 , i , j , k
        DOUBLE PRECISION dp
        IF ( M.LT.1 .OR. N.LT.4 ) THEN
           Ier = -1
           GOTO 99999
        ELSE
           nm1 = N - 1
           nm2 = N - 2
  ! Process the first block-row. -----------------------------------------
           CALL DGEFA(A,M,M,Ip,Ier)
           k = 1
           IF ( Ier.EQ.0 ) THEN
              DO j = 1 , M
                 CALL DGESL(A,M,M,Ip,B(1,j,1),0)
                 CALL DGESL(A,M,M,Ip,C(1,j,1),0)
              ENDDO
  ! Adjust B(*,*,2). -----------------------------------------------------
              DO j = 1 , M
                 DO i = 1 , M
                    dp = DDOT(M,C(i,1,2),M,C(1,j,1),1)
                    B(i,j,2) = B(i,j,2) - dp
                 ENDDO
              ENDDO
  ! Main loop.  Process block-rows 2 to N-1. -----------------------------
              DO k = 2 , nm1
                 km1 = k - 1
                 DO j = 1 , M
                    DO i = 1 , M
                       dp = DDOT(M,C(i,1,k),M,B(1,j,km1),1)
                       A(i,j,k) = A(i,j,k) - dp
                    ENDDO
                 ENDDO
                 CALL DGEFA(A(1,1,k),M,M,Ip(1,k),Ier)
                 IF ( Ier.NE.0 ) GOTO 100
                 DO j = 1 , M
                    CALL DGESL(A(1,1,k),M,M,Ip(1,k),B(1,j,k),0)
                 ENDDO
              ENDDO
  ! Process last block-row and return. -----------------------------------
              DO j = 1 , M
                 DO i = 1 , M
                    dp = DDOT(M,B(i,1,N),M,B(1,j,nm2),1)
                    C(i,j,N) = C(i,j,N) - dp
                 ENDDO
              ENDDO
              DO j = 1 , M
                 DO i = 1 , M
                    dp = DDOT(M,C(i,1,N),M,B(1,j,nm1),1)
                    A(i,j,N) = A(i,j,N) - dp
                 ENDDO
              ENDDO
              CALL DGEFA(A(1,1,N),M,M,Ip(1,N),Ier)
              k = N
              IF ( Ier.EQ.0 ) RETURN
           ENDIF
        ENDIF
  ! Error returns. -------------------------------------------------------
   100  Ier = k
        RETURN
  !----------------------- End of Subroutine DDECBT ----------------------
  99999 END
  !*==DSOLBT.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DSOLBT
        SUBROUTINE DSOLBT(M,N,A,B,C,Y,Ip)
        IMPLICIT NONE
  !*--DSOLBT10273
        INTEGER M , N , Ip(M,N)
        DOUBLE PRECISION A(M,M,N) , B(M,M,N) , C(M,M,N) , Y(M,N)
  !-----------------------------------------------------------------------
  ! Solution of block-tridiagonal linear system.
  ! Coefficient matrix must have been previously processed by DDECBT.
  ! M, N, A,B,C, and IP  must not have been changed since call to DDECBT.
  ! Written by A. C. Hindmarsh.
  ! Input:
  !     M = order of each block.
  !     N = number of blocks in each direction of matrix.
  ! A,B,C = M by M by N arrays containing block LU decomposition
  !         of coefficient matrix from DDECBT.
  !    IP = M by N integer array of pivot information from DDECBT.
  !     Y = array of length M*N containg the right-hand side vector
  !         (treated as an M by N array here).
  ! Output:
  !     Y = solution vector, of length M*N.
  !
  ! External routines required: DGESL (LINPACK) and DDOT (BLAS).
  !-----------------------------------------------------------------------
  !
        INTEGER nm1 , nm2 , i , k , kb , km1 , kp1
        DOUBLE PRECISION dp
        nm1 = N - 1
        nm2 = N - 2
  ! Forward solution sweep. ----------------------------------------------
        CALL DGESL(A,M,M,Ip,Y,0)
        DO k = 2 , nm1
           km1 = k - 1
           DO i = 1 , M
              dp = DDOT(M,C(i,1,k),M,Y(1,km1),1)
              Y(i,k) = Y(i,k) - dp
           ENDDO
           CALL DGESL(A(1,1,k),M,M,Ip(1,k),Y(1,k),0)
        ENDDO
        DO i = 1 , M
           dp = DDOT(M,C(i,1,N),M,Y(1,nm1),1)                             &
       &        + DDOT(M,B(i,1,N),M,Y(1,nm2),1)
           Y(i,N) = Y(i,N) - dp
        ENDDO
        CALL DGESL(A(1,1,N),M,M,Ip(1,N),Y(1,N),0)
  ! Backward solution sweep. ---------------------------------------------
        DO kb = 1 , nm1
           k = N - kb
           kp1 = k + 1
           DO i = 1 , M
              dp = DDOT(M,B(i,1,k),M,Y(1,kp1),1)
              Y(i,k) = Y(i,k) - dp
           ENDDO
        ENDDO
        DO i = 1 , M
           dp = DDOT(M,C(i,1,1),M,Y(1,3),1)
           Y(i,1) = Y(i,1) - dp
        ENDDO
  !----------------------- End of Subroutine DSOLBT ----------------------
        END
  !*==DIPREPI.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DIPREPI
        SUBROUTINE DIPREPI(Neq,Y,S,Rwork,Ia,Ja,Ic,Jc,Ipflag,RES,JAC,ADDA)
        IMPLICIT NONE
  !*--DIPREPI10334
  !*** Start of declarations inserted by SPAG
        !REAL ADDA , RES
        !INTEGER JAC
  !*** End of declarations inserted by SPAG
        EXTERNAL RES , JAC , ADDA
        INTEGER Neq , Ia , Ja , Ic , Jc , Ipflag
        DOUBLE PRECISION Y , S , Rwork
        DIMENSION Neq(*) , Y(*) , S(*) , Rwork(*) , Ia(*) , Ja(*) ,       &
       &          Ic(*) , Jc(*)
        INTEGER IOWnd , IOWns , ICF , IERpj , IERsl , JCUr , JSTart ,     &
       &        KFLag , L , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,     &
       &        METh , MITer , MAXord , MAXcor , MSBp , MXNcf , N , NQ ,  &
       &        NST , NFE , NJE , NQU
        INTEGER IPLost , IESp , ISTatc , IYS , IBA , IBIan , IBJan ,      &
       &        IBJgp , IPIan , IPJan , IPJgp , IPIgp , IPR , IPC , IPIc ,&
       &        IPIsp , IPRsp , IPA , LENyh , LENyhm , LENwk , LREq ,     &
       &        LRAt , LREst , LWMin , MOSs , MSBj , NSLj , NGP , NLU ,   &
       &        NNZ , NSP , NZL , NZU
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        DOUBLE PRECISION RLSs
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , IOWnd(6) , IOWns(6) , ICF ,    &
       &                IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , &
       &                LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,&
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        COMMON /DLSS01/ RLSs(6) , IPLost , IESp , ISTatc , IYS , IBA ,    &
       &                IBIan , IBJan , IBJgp , IPIan , IPJan , IPJgp ,   &
       &                IPIgp , IPR , IPC , IPIc , IPIsp , IPRsp , IPA ,  &
       &                LENyh , LENyhm , LENwk , LREq , LRAt , LREst ,    &
       &                LWMin , MOSs , MSBj , NSLj , NGP , NLU , NNZ ,    &
       &                NSP , NZL , NZU
        INTEGER i , imax , lewtn , lyhd , lyhn
        INTEGER, DIMENSION(LENwk) :: Rwork_i
  !-----------------------------------------------------------------------
  ! This routine serves as an interface between the driver and
  ! Subroutine DPREPI.  Tasks performed here are:
  !  * call DPREPI,
  !  * reset the required WM segment length LENWK,
  !  * move YH back to its final location (following WM in RWORK),
  !  * reset pointers for YH, SAVR, EWT, and ACOR, and
  !  * move EWT to its new position if ISTATE = 0 or 1.
  ! IPFLAG is an output error indication flag.  IPFLAG = 0 if there was
  ! no trouble, and IPFLAG is the value of the DPREPI error flag IPPER
  ! if there was trouble in Subroutine DPREPI.
  !-----------------------------------------------------------------------
        Ipflag = 0
  ! Call DPREPI to do matrix preprocessing operations. -------------------
        Rwork_i = Rwork(LWM)
        CALL DPREPI(Neq,Y,S,Rwork(LYH),Rwork(LSAvf),Rwork(LEWt),          &
       &            Rwork(LACor),Ia,Ja,Ic,Jc,Rwork(LWM),Rwork_i,Ipflag,&
       &            RES,JAC,ADDA)
        LENwk = MAX(LREq,LWMin)
        IF ( Ipflag.LT.0 ) RETURN
  ! If DPREPI was successful, move YH to end of required space for WM. ---
        lyhn = LWM + LENwk
        IF ( lyhn.GT.LYH ) RETURN
        lyhd = LYH - lyhn
        IF ( lyhd.NE.0 ) THEN
           imax = lyhn - 1 + LENyhm
           DO i = lyhn , imax
              Rwork(i) = Rwork(i+lyhd)
           ENDDO
           LYH = lyhn
        ENDIF
  ! Reset pointers for SAVR, EWT, and ACOR. ------------------------------
        LSAvf = LYH + LENyh
        lewtn = LSAvf + N
        LACor = lewtn + N
        IF ( ISTatc.NE.3 ) THEN
  ! If ISTATE = 1, move EWT (left) to its new position. ------------------
           IF ( lewtn.GT.LEWt ) RETURN
           DO i = 1 , N
              Rwork(i+lewtn-1) = Rwork(i+LEWt-1)
           ENDDO
        ENDIF
        LEWt = lewtn
  !----------------------- End of Subroutine DIPREPI ---------------------
        END
  !*==DPREPI.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DPREPI
        SUBROUTINE DPREPI(Neq,Y,S,Yh,Savr,Ewt,Rtem,Ia,Ja,Ic,Jc,Wk,Iwk,    &
       &                  Ipper,RES,JAC,ADDA)
        IMPLICIT NONE
  !*--DPREPI10418
        EXTERNAL RES , JAC , ADDA
        INTEGER Neq , Ia , Ja , Ic , Jc , Iwk , Ipper
        DOUBLE PRECISION Y , S , Yh , Savr , Ewt , Rtem , Wk
        DIMENSION Neq(*) , Y(*) , S(*) , Yh(*) , Savr(*) , Ewt(*) ,       &
       &          Rtem(*) , Ia(*) , Ja(*) , Ic(*) , Jc(*) , Wk(*) , Iwk(*)
        INTEGER IOWnd , IOWns , ICF , IERpj , IERsl , JCUr , JSTart ,     &
       &        KFLag , L , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,     &
       &        METh , MITer , MAXord , MAXcor , MSBp , MXNcf , N , NQ ,  &
       &        NST , NFE , NJE , NQU
        INTEGER IPLost , IESp , ISTatc , IYS , IBA , IBIan , IBJan ,      &
       &        IBJgp , IPIan , IPJan , IPJgp , IPIgp , IPR , IPC , IPIc ,&
       &        IPIsp , IPRsp , IPA , LENyh , LENyhm , LENwk , LREq ,     &
       &        LRAt , LREst , LWMin , MOSs , MSBj , NSLj , NGP , NLU ,   &
       &        NNZ , NSP , NZL , NZU
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        DOUBLE PRECISION RLSs
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , IOWnd(6) , IOWns(6) , ICF ,    &
       &                IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , &
       &                LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,&
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        COMMON /DLSS01/ RLSs(6) , IPLost , IESp , ISTatc , IYS , IBA ,    &
       &                IBIan , IBJan , IBJgp , IPIan , IPJan , IPJgp ,   &
       &                IPIgp , IPR , IPC , IPIc , IPIsp , IPRsp , IPA ,  &
       &                LENyh , LENyhm , LENwk , LREq , LRAt , LREst ,    &
       &                LWMin , MOSs , MSBj , NSLj , NGP , NLU , NNZ ,    &
       &                NSP , NZL , NZU
        INTEGER i , ibr , ier , ipil , ipiu , iptt1 , iptt2 , j , k ,     &
       &        knew , kamax , kamin , kcmax , kcmin , ldif , lenigp ,    &
       &        lenwk1 , liwk , ljfo , maxg , np1 , nzsut
        DOUBLE PRECISION erwt , fac , yj
  !-----------------------------------------------------------------------
  ! This routine performs preprocessing related to the sparse linear
  ! systems that must be solved.
  ! The operations that are performed here are:
  !  * compute sparseness structure of the iteration matrix
  !      P = A - con*J  according to MOSS,
  !  * compute grouping of column indices (MITER = 2),
  !  * compute a new ordering of rows and columns of the matrix,
  !  * reorder JA corresponding to the new ordering,
  !  * perform a symbolic LU factorization of the matrix, and
  !  * set pointers for segments of the IWK/WK array.
  ! In addition to variables described previously, DPREPI uses the
  ! following for communication:
  ! YH     = the history array.  Only the first column, containing the
  !          current Y vector, is used.  Used only if MOSS .ne. 0.
  ! S      = array of length NEQ, identical to YDOTI in the driver, used
  !          only if MOSS .ne. 0.
  ! SAVR   = a work array of length NEQ, used only if MOSS .ne. 0.
  ! EWT    = array of length NEQ containing (inverted) error weights.
  !          Used only if MOSS = 2 or 4 or if ISTATE = MOSS = 1.
  ! RTEM   = a work array of length NEQ, identical to ACOR in the driver,
  !          used only if MOSS = 2 or 4.
  ! WK     = a real work array of length LENWK, identical to WM in
  !          the driver.
  ! IWK    = integer work array, assumed to occupy the same space as WK.
  ! LENWK  = the length of the work arrays WK and IWK.
  ! ISTATC = a copy of the driver input argument ISTATE (= 1 on the
  !          first call, = 3 on a continuation call).
  ! IYS    = flag value from ODRV or CDRV.
  ! IPPER  = output error flag , with the following values and meanings:
  !        =   0  no error.
  !        =  -1  insufficient storage for internal structure pointers.
  !        =  -2  insufficient storage for JGROUP.
  !        =  -3  insufficient storage for ODRV.
  !        =  -4  other error flag from ODRV (should never occur).
  !        =  -5  insufficient storage for CDRV.
  !        =  -6  other error flag from CDRV.
  !        =  -7  if the RES routine returned error flag IRES = IER = 2.
  !        =  -8  if the RES routine returned error flag IRES = IER = 3.
  !-----------------------------------------------------------------------
        IBIan = LRAt*2
        IPIan = IBIan + 1
        np1 = N + 1
        IPJan = IPIan + np1
        IBJan = IPJan - 1
        lenwk1 = LENwk - N
        liwk = LENwk*LRAt
        IF ( MOSs.EQ.0 ) liwk = liwk - N
        IF ( MOSs.EQ.1 .OR. MOSs.EQ.2 ) liwk = lenwk1*LRAt
        IF ( IPJan+N-1.GT.liwk ) GOTO 600
        IF ( MOSs.NE.0 ) THEN
  !
           IF ( ISTatc.NE.3 ) THEN
  ! ISTATE = 1 and MOSS .ne. 0.  Perturb Y for structure determination.
  ! Initialize S with random nonzero elements for structure determination.
              DO i = 1 , N
                 erwt = 1.0D0/Ewt(i)
                 fac = 1.0D0 + 1.0D0/(i+1.0D0)
                 Y(i) = Y(i) + fac*SIGN(erwt,Y(i))
                 S(i) = 1.0D0 + fac*erwt
              ENDDO
              GOTO (100,200,300,400) , MOSs
           ENDIF
  !
  ! ISTATE = 3 and MOSS .ne. 0. Load Y from YH(*,1) and S from YH(*,2). --
           DO i = 1 , N
              Y(i) = Yh(i)
              S(i) = Yh(N+i)
           ENDDO
           GOTO (100,200,300,400) , MOSs
        ENDIF
  !
  ! MOSS = 0. Process user's IA,JA and IC,JC. ----------------------------
        knew = IPJan
        kamin = Ia(1)
        kcmin = Ic(1)
        Iwk(IPIan) = 1
        DO j = 1 , N
           DO i = 1 , N
              Iwk(liwk+i) = 0
           ENDDO
           kamax = Ia(j+1) - 1
           IF ( kamin.LE.kamax ) THEN
              DO k = kamin , kamax
                 i = Ja(k)
                 Iwk(liwk+i) = 1
                 IF ( knew.GT.liwk ) GOTO 600
                 Iwk(knew) = i
                 knew = knew + 1
              ENDDO
           ENDIF
           kamin = kamax + 1
           kcmax = Ic(j+1) - 1
           IF ( kcmin.LE.kcmax ) THEN
              DO k = kcmin , kcmax
                 i = Jc(k)
                 IF ( Iwk(liwk+i).EQ.0 ) THEN
                    IF ( knew.GT.liwk ) GOTO 600
                    Iwk(knew) = i
                    knew = knew + 1
                 ENDIF
              ENDDO
           ENDIF
           Iwk(IPIan+j) = knew + 1 - IPJan
           kcmin = kcmax + 1
        ENDDO
        GOTO 500
  !
  ! MOSS = 1. Compute structure from user-supplied Jacobian routine JAC. -
  ! A dummy call to RES allows user to create temporaries for use in JAC.
   100  ier = 1
        CALL RES(Neq,TN,Y,S,Savr,ier)
        IF ( ier.GT.1 ) GOTO 1000
        DO i = 1 , N
           Savr(i) = 0.0D0
           Wk(lenwk1+i) = 0.0D0
        ENDDO
        k = IPJan
        Iwk(IPIan) = 1
        DO j = 1 , N
           CALL ADDA(Neq,TN,Y,j,Iwk(IPIan),Iwk(IPJan),Wk(lenwk1+1))
           CALL JAC(Neq,TN,Y,S,j,Iwk(IPIan),Iwk(IPJan),Savr)
           DO i = 1 , N
              ljfo = lenwk1 + i
              IF ( Wk(ljfo).EQ.0.0D0 ) THEN
                 IF ( Savr(i).EQ.0.0D0 ) GOTO 150
                 Savr(i) = 0.0D0
              ELSE
                 Wk(ljfo) = 0.0D0
                 Savr(i) = 0.0D0
              ENDIF
              IF ( k.GT.liwk ) GOTO 600
              Iwk(k) = i
              k = k + 1
   150     ENDDO
           Iwk(IPIan+j) = k + 1 - IPJan
        ENDDO
        GOTO 500
  !
  ! MOSS = 2. Compute structure from results of N + 1 calls to RES. ------
   200  DO i = 1 , N
           Wk(lenwk1+i) = 0.0D0
        ENDDO
        k = IPJan
        Iwk(IPIan) = 1
        ier = -1
        IF ( MITer.EQ.1 ) ier = 1
        CALL RES(Neq,TN,Y,S,Savr,ier)
        IF ( ier.GT.1 ) GOTO 1000
        DO j = 1 , N
           CALL ADDA(Neq,TN,Y,j,Iwk(IPIan),Iwk(IPJan),Wk(lenwk1+1))
           yj = Y(j)
           erwt = 1.0D0/Ewt(j)
           Y(j) = yj + SIGN(erwt,yj)
           CALL RES(Neq,TN,Y,S,Rtem,ier)
           IF ( ier.GT.1 ) RETURN
           Y(j) = yj
           DO i = 1 , N
              ljfo = lenwk1 + i
              IF ( Wk(ljfo).NE.0.0D0 ) THEN
                 Wk(ljfo) = 0.0D0
              ELSEIF ( Rtem(i).EQ.Savr(i) ) THEN
                 GOTO 250
              ENDIF
              IF ( k.GT.liwk ) GOTO 600
              Iwk(k) = i
              k = k + 1
   250     ENDDO
           Iwk(IPIan+j) = k + 1 - IPJan
        ENDDO
        GOTO 500
  !
  ! MOSS = 3. Compute structure from the user's IA/JA and JAC routine. ---
  ! A dummy call to RES allows user to create temporaries for use in JAC.
   300  ier = 1
        CALL RES(Neq,TN,Y,S,Savr,ier)
        IF ( ier.GT.1 ) GOTO 1000
        DO i = 1 , N
           Savr(i) = 0.0D0
        ENDDO
        knew = IPJan
        kamin = Ia(1)
        Iwk(IPIan) = 1
        DO j = 1 , N
           CALL JAC(Neq,TN,Y,S,j,Iwk(IPIan),Iwk(IPJan),Savr)
           kamax = Ia(j+1) - 1
           IF ( kamin.LE.kamax ) THEN
              DO k = kamin , kamax
                 i = Ja(k)
                 Savr(i) = 0.0D0
                 IF ( knew.GT.liwk ) GOTO 600
                 Iwk(knew) = i
                 knew = knew + 1
              ENDDO
           ENDIF
           kamin = kamax + 1
           DO i = 1 , N
              IF ( Savr(i).NE.0.0D0 ) THEN
                 Savr(i) = 0.0D0
                 IF ( knew.GT.liwk ) GOTO 600
                 Iwk(knew) = i
                 knew = knew + 1
              ENDIF
           ENDDO
           Iwk(IPIan+j) = knew + 1 - IPJan
        ENDDO
        GOTO 500
  !
  ! MOSS = 4. Compute structure from user's IA/JA and N + 1 RES calls. ---
   400  knew = IPJan
        kamin = Ia(1)
        Iwk(IPIan) = 1
        ier = -1
        IF ( MITer.EQ.1 ) ier = 1
        CALL RES(Neq,TN,Y,S,Savr,ier)
        IF ( ier.GT.1 ) GOTO 1000
        DO j = 1 , N
           yj = Y(j)
           erwt = 1.0D0/Ewt(j)
           Y(j) = yj + SIGN(erwt,yj)
           CALL RES(Neq,TN,Y,S,Rtem,ier)
           IF ( ier.GT.1 ) RETURN
           Y(j) = yj
           kamax = Ia(j+1) - 1
           IF ( kamin.LE.kamax ) THEN
              DO k = kamin , kamax
                 i = Ja(k)
                 Rtem(i) = Savr(i)
                 IF ( knew.GT.liwk ) GOTO 600
                 Iwk(knew) = i
                 knew = knew + 1
              ENDDO
           ENDIF
           kamin = kamax + 1
           DO i = 1 , N
              IF ( Rtem(i).NE.Savr(i) ) THEN
                 IF ( knew.GT.liwk ) GOTO 600
                 Iwk(knew) = i
                 knew = knew + 1
              ENDIF
           ENDDO
           Iwk(IPIan+j) = knew + 1 - IPJan
        ENDDO
  !
   500  IF ( MOSs.NE.0 .AND. ISTatc.NE.3 ) THEN
  ! If ISTATE = 0 or 1 and MOSS .ne. 0, restore Y from YH. ---------------
           DO i = 1 , N
              Y(i) = Yh(i)
           ENDDO
        ENDIF
        NNZ = Iwk(IPIan+N) - 1
        Ipper = 0
        NGP = 0
        lenigp = 0
        IPIgp = IPJan + NNZ
        IF ( MITer.EQ.2 ) THEN
  !
  ! Compute grouping of column indices (MITER = 2). ----------------------
  !
           maxg = np1
           IPJgp = IPJan + NNZ
           IBJgp = IPJgp - 1
           IPIgp = IPJgp + N
           iptt1 = IPIgp + np1
           iptt2 = iptt1 + N
           LREq = iptt2 + N - 1
           IF ( LREq.GT.liwk ) GOTO 700
           CALL JGROUP(N,Iwk(IPIan),Iwk(IPJan),maxg,NGP,Iwk(IPIgp),       &
       &               Iwk(IPJgp),Iwk(iptt1),Iwk(iptt2),ier)
           IF ( ier.NE.0 ) GOTO 700
           lenigp = NGP + 1
        ENDIF
  !
  ! Compute new ordering of rows/columns of Jacobian. --------------------
        IPR = IPIgp + lenigp
        IPC = IPR
        IPIc = IPC + N
        IPIsp = IPIc + N
        IPRsp = (IPIsp-2)/LRAt + 2
        IESp = LENwk + 1 - IPRsp
        IF ( IESp.LT.0 ) GOTO 800
        ibr = IPR - 1
        DO i = 1 , N
           Iwk(ibr+i) = i
        ENDDO
        NSP = liwk + 1 - IPIsp
        CALL ODRV(N,Iwk(IPIan),Iwk(IPJan),Wk,Iwk(IPR),Iwk(IPIc),NSP,      &
       &          Iwk(IPIsp),1,IYS)
        IF ( IYS.EQ.11*N+1 ) THEN
  !
           Ipper = -4
           RETURN
        ELSE
           IF ( IYS.NE.0 ) GOTO 800
  !
  ! Reorder JAN and do symbolic LU factorization of matrix. --------------
           IPA = LENwk + 1 - NNZ
           NSP = IPA - IPRsp
           LREq = MAX(12*N/LRAt,6*N/LRAt+2*N+NNZ) + 3
           LREq = LREq + IPRsp - 1 + NNZ
           IF ( LREq.GT.LENwk ) GOTO 900
           IBA = IPA - 1
           DO i = 1 , NNZ
              Wk(IBA+i) = 0.0D0
           ENDDO
           IPIsp = LRAt*(IPRsp-1) + 1
           CALL CDRV(N,Iwk(IPR),Iwk(IPC),Iwk(IPIc),Iwk(IPIan),Iwk(IPJan), &
       &             Wk(IPA),Wk(IPA),Wk(IPA),NSP,Iwk(IPIsp),Wk(IPRsp),    &
       &             IESp,5,IYS)
           LREq = LENwk - IESp
           IF ( IYS.EQ.10*N+1 ) GOTO 900
           IF ( IYS.NE.0 ) THEN
  !
              Ipper = -6
              LREq = LENwk
              RETURN
           ELSE
              ipil = IPIsp
              ipiu = ipil + 2*N + 1
              NZU = Iwk(ipil+N) - Iwk(ipil)
              NZL = Iwk(ipiu+N) - Iwk(ipiu)
              IF ( LRAt.LE.1 ) THEN
                 CALL ADJLR(N,Iwk(IPIsp),ldif)
                 LREq = LREq + ldif
              ENDIF
              IF ( LRAt.EQ.2 .AND. NNZ.EQ.N ) LREq = LREq + 1
              NSP = NSP + LREq - LENwk
              IPA = LREq + 1 - NNZ
              IBA = IPA - 1
              Ipper = 0
              RETURN
           ENDIF
        ENDIF
  !
   600  Ipper = -1
        LREq = 2 + (2*N+1)/LRAt
        LREq = MAX(LENwk+1,LREq)
        RETURN
  !
   700  Ipper = -2
        LREq = (LREq-1)/LRAt + 1
        RETURN
  !
   800  Ipper = -3
        CALL CNTNZU(N,Iwk(IPIan),Iwk(IPJan),nzsut)
        LREq = LENwk - IESp + (3*N+4*nzsut-1)/LRAt + 1
        RETURN
  !
   900  Ipper = -5
        RETURN
  !
   1000 Ipper = -ier - 5
        LREq = 2 + (2*N+1)/LRAt
  !----------------------- End of Subroutine DPREPI ----------------------
        END
  !*==DAINVGS.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DAINVGS
        SUBROUTINE DAINVGS(Neq,T,Y,Wk,Iwk,Tem,Ydot,Ier,RES,ADDA)
        IMPLICIT NONE
  !*--DAINVGS10811
        EXTERNAL RES , ADDA
        INTEGER Neq , Iwk , Ier
        INTEGER IPLost , IESp , ISTatc , IYS , IBA , IBIan , IBJan ,      &
       &        IBJgp , IPIan , IPJan , IPJgp , IPIgp , IPR , IPC , IPIc ,&
       &        IPIsp , IPRsp , IPA , LENyh , LENyhm , LENwk , LREq ,     &
       &        LRAt , LREst , LWMin , MOSs , MSBj , NSLj , NGP , NLU ,   &
       &        NNZ , NSP , NZL , NZU
        INTEGER i , imul , j , k , kmin , kmax
        DOUBLE PRECISION T , Y , Wk , Tem , Ydot
        DOUBLE PRECISION RLSs
        DIMENSION Y(*) , Wk(*) , Iwk(*) , Tem(*) , Ydot(*)
        COMMON /DLSS01/ RLSs(6) , IPLost , IESp , ISTatc , IYS , IBA ,    &
       &                IBIan , IBJan , IBJgp , IPIan , IPJan , IPJgp ,   &
       &                IPIgp , IPR , IPC , IPIc , IPIsp , IPRsp , IPA ,  &
       &                LENyh , LENyhm , LENwk , LREq , LRAt , LREst ,    &
       &                LWMin , MOSs , MSBj , NSLj , NGP , NLU , NNZ ,    &
       &                NSP , NZL , NZU
  !-----------------------------------------------------------------------
  ! This subroutine computes the initial value of the vector YDOT
  ! satisfying
  !     A * YDOT = g(t,y)
  ! when A is nonsingular.  It is called by DLSODIS for initialization
  ! only, when ISTATE = 0.  The matrix A is subjected to LU
  ! decomposition in CDRV.  Then the system A*YDOT = g(t,y) is solved
  ! in CDRV.
  ! In addition to variables described previously, communication
  ! with DAINVGS uses the following:
  ! Y     = array of initial values.
  ! WK    = real work space for matrices.  On output it contains A and
  !         its LU decomposition.  The LU decomposition is not entirely
  !         sparse unless the structure of the matrix A is identical to
  !         the structure of the Jacobian matrix dr/dy.
  !         Storage of matrix elements starts at WK(3).
  !         WK(1) = SQRT(UROUND), not used here.
  ! IWK   = integer work space for matrix-related data, assumed to
  !         be equivalenced to WK.  In addition, WK(IPRSP) and WK(IPISP)
  !         are assumed to have identical locations.
  ! TEM   = vector of work space of length N (ACOR in DSTODI).
  ! YDOT  = output vector containing the initial dy/dt. YDOT(i) contains
  !         dy(i)/dt when the matrix A is non-singular.
  ! IER   = output error flag with the following values and meanings:
  !       = 0  if DAINVGS was successful.
  !       = 1  if the A-matrix was found to be singular.
  !       = 2  if RES returned an error flag IRES = IER = 2.
  !       = 3  if RES returned an error flag IRES = IER = 3.
  !       = 4  if insufficient storage for CDRV (should not occur here).
  !       = 5  if other error found in CDRV (should not occur here).
  !-----------------------------------------------------------------------
  !
        DO i = 1 , NNZ
           Wk(IBA+i) = 0.0D0
        ENDDO
  !
        Ier = 1
        CALL RES(Neq,T,Y,Wk(IPA),Ydot,Ier)
        IF ( Ier.GT.1 ) RETURN
  !
        kmin = Iwk(IPIan)
        DO j = 1 , Neq
           kmax = Iwk(IPIan+j) - 1
           DO k = kmin , kmax
              i = Iwk(IBJan+k)
              Tem(i) = 0.0D0
           ENDDO
           CALL ADDA(Neq,T,Y,j,Iwk(IPIan),Iwk(IPJan),Tem)
           DO k = kmin , kmax
              i = Iwk(IBJan+k)
              Wk(IBA+k) = Tem(i)
           ENDDO
           kmin = kmax + 1
        ENDDO
        NLU = NLU + 1
        Ier = 0
        DO i = 1 , Neq
           Tem(i) = 0.0D0
        ENDDO
  !
  ! Numerical factorization of matrix A. ---------------------------------
        CALL CDRV(Neq,Iwk(IPR),Iwk(IPC),Iwk(IPIc),Iwk(IPIan),Iwk(IPJan),  &
       &          Wk(IPA),Tem,Tem,NSP,Iwk(IPIsp),Wk(IPRsp),IESp,2,IYS)
        IF ( IYS.EQ.0 ) THEN
  !
  ! Solution of the linear system. ---------------------------------------
           CALL CDRV(Neq,Iwk(IPR),Iwk(IPC),Iwk(IPIc),Iwk(IPIan),Iwk(IPJan)&
       &             ,Wk(IPA),Ydot,Ydot,NSP,Iwk(IPIsp),Wk(IPRsp),IESp,4,  &
       &             IYS)
           IF ( IYS.NE.0 ) Ier = 5
           GOTO 99999
        ENDIF
        imul = (IYS-1)/Neq
        Ier = 5
        IF ( imul.EQ.8 ) Ier = 1
        IF ( imul.EQ.10 ) Ier = 4
        RETURN
  !----------------------- End of Subroutine DAINVGS ---------------------
  99999 END
  !*==DPRJIS.spg  processed by SPAG 6.72Dc at 17:55 on 29 Mar 2017
  !DECK DPRJIS
        SUBROUTINE DPRJIS(Neq,Y,Yh,Nyh,Ewt,Rtem,Savr,S,Wk,Iwk,RES,JAC,    &
       &                  ADDA)
        IMPLICIT NONE
  !*--DPRJIS10913
        EXTERNAL RES , JAC , ADDA
        INTEGER Neq , Nyh , Iwk
        DOUBLE PRECISION Y , Yh , Ewt , Rtem , Savr , S , Wk
        DIMENSION Neq(*) , Y(*) , Yh(Nyh,*) , Ewt(*) , Rtem(*) , S(*) ,   &
       &          Savr(*) , Wk(*) , Iwk(*)
        INTEGER IOWnd , IOWns , ICF , IERpj , IERsl , JCUr , JSTart ,     &
       &        KFLag , L , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,     &
       &        METh , MITer , MAXord , MAXcor , MSBp , MXNcf , N , NQ ,  &
       &        NST , NFE , NJE , NQU
        INTEGER IPLost , IESp , ISTatc , IYS , IBA , IBIan , IBJan ,      &
       &        IBJgp , IPIan , IPJan , IPJgp , IPIgp , IPR , IPC , IPIc ,&
       &        IPIsp , IPRsp , IPA , LENyh , LENyhm , LENwk , LREq ,     &
       &        LRAt , LREst , LWMin , MOSs , MSBj , NSLj , NGP , NLU ,   &
       &        NNZ , NSP , NZL , NZU
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        DOUBLE PRECISION RLSs
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , IOWnd(6) , IOWns(6) , ICF ,    &
       &                IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH , &
       &                LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,&
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
        COMMON /DLSS01/ RLSs(6) , IPLost , IESp , ISTatc , IYS , IBA ,    &
       &                IBIan , IBJan , IBJgp , IPIan , IPJan , IPJgp ,   &
       &                IPIgp , IPR , IPC , IPIc , IPIsp , IPRsp , IPA ,  &
       &                LENyh , LENyhm , LENwk , LREq , LRAt , LREst ,    &
       &                LWMin , MOSs , MSBj , NSLj , NGP , NLU , NNZ ,    &
       &                NSP , NZL , NZU
        INTEGER i , imul , ires , j , jj , jmax , jmin , k , kmax , kmin ,&
       &        ng
        DOUBLE PRECISION con , fac , hl0 , r , srur
  !-----------------------------------------------------------------------
  ! DPRJIS is called to compute and process the matrix
  ! P = A - H*EL(1)*J, where J is an approximation to the Jacobian dr/dy,
  ! where r = g(t,y) - A(t,y)*s.  J is computed by columns, either by
  ! the user-supplied routine JAC if MITER = 1, or by finite differencing
  ! if MITER = 2.  J is stored in WK, rescaled, and ADDA is called to
  ! generate P.  The matrix P is subjected to LU decomposition in CDRV.
  ! P and its LU decomposition are stored separately in WK.
  !
  ! In addition to variables described previously, communication
  ! with DPRJIS uses the following:
  ! Y     = array containing predicted values on entry.
  ! RTEM  = work array of length N (ACOR in DSTODI).
  ! SAVR  = array containing r evaluated at predicted y. On output it
  !         contains the residual evaluated at current values of t and y.
  ! S     = array containing predicted values of dy/dt (SAVF in DSTODI).
  ! WK    = real work space for matrices.  On output it contains P and
  !         its sparse LU decomposition.  Storage of matrix elements
  !         starts at WK(3).
  !         WK also contains the following matrix-related data.
  !         WK(1) = SQRT(UROUND), used in numerical Jacobian increments.
  ! IWK   = integer work space for matrix-related data, assumed to be
  !         equivalenced to WK.  In addition,  WK(IPRSP) and IWK(IPISP)
  !         are assumed to have identical locations.
  ! EL0   = EL(1) (input).
  ! IERPJ = output error flag (in COMMON).
  !         =  0 if no error.
  !         =  1 if zero pivot found in CDRV.
  !         = IRES (= 2 or 3) if RES returned IRES = 2 or 3.
  !         = -1 if insufficient storage for CDRV (should not occur).
  !         = -2 if other error found in CDRV (should not occur here).
  ! JCUR  = output flag = 1 to indicate that the Jacobian matrix
  !         (or approximation) is now current.
  ! This routine also uses other variables in Common.
  !-----------------------------------------------------------------------
        hl0 = H*EL0
        con = -hl0
        JCUr = 1
        NJE = NJE + 1
        GOTO (100,200) , MITer
  !
  ! If MITER = 1, call RES, then call JAC and ADDA for each column. ------
   100  ires = 1
        CALL RES(Neq,TN,Y,S,Savr,ires)
        NFE = NFE + 1
        IF ( ires.GT.1 ) GOTO 400
        kmin = Iwk(IPIan)
        DO j = 1 , N
           kmax = Iwk(IPIan+j) - 1
           DO i = 1 , N
              Rtem(i) = 0.0D0
           ENDDO
           CALL JAC(Neq,TN,Y,S,j,Iwk(IPIan),Iwk(IPJan),Rtem)
           DO i = 1 , N
              Rtem(i) = Rtem(i)*con
           ENDDO
           CALL ADDA(Neq,TN,Y,j,Iwk(IPIan),Iwk(IPJan),Rtem)
           DO k = kmin , kmax
              i = Iwk(IBJan+k)
              Wk(IBA+k) = Rtem(i)
           ENDDO
           kmin = kmax + 1
        ENDDO
        GOTO 300
  !
  ! If MITER = 2, make NGP + 1 calls to RES to approximate J and P. ------
   200  ires = -1
        CALL RES(Neq,TN,Y,S,Savr,ires)
        NFE = NFE + 1
        IF ( ires.GT.1 ) GOTO 400
        srur = Wk(1)
        jmin = Iwk(IPIgp)
        DO ng = 1 , NGP
           jmax = Iwk(IPIgp+ng) - 1
           DO j = jmin , jmax
              jj = Iwk(IBJgp+j)
              r = MAX(srur*ABS(Y(jj)),0.01D0/Ewt(jj))
              Y(jj) = Y(jj) + r
           ENDDO
           CALL RES(Neq,TN,Y,S,Rtem,ires)
           NFE = NFE + 1
           IF ( ires.GT.1 ) GOTO 400
           DO j = jmin , jmax
              jj = Iwk(IBJgp+j)
              Y(jj) = Yh(jj,1)
              r = MAX(srur*ABS(Y(jj)),0.01D0/Ewt(jj))
              fac = -hl0/r
              kmin = Iwk(IBIan+jj)
              kmax = Iwk(IBIan+jj+1) - 1
              DO k = kmin , kmax
                 i = Iwk(IBJan+k)
                 Rtem(i) = (Rtem(i)-Savr(i))*fac
              ENDDO
              CALL ADDA(Neq,TN,Y,jj,Iwk(IPIan),Iwk(IPJan),Rtem)
              DO k = kmin , kmax
                 i = Iwk(IBJan+k)
                 Wk(IBA+k) = Rtem(i)
              ENDDO
           ENDDO
           jmin = jmax + 1
        ENDDO
        ires = 1
        CALL RES(Neq,TN,Y,S,Savr,ires)
        NFE = NFE + 1
        IF ( ires.GT.1 ) GOTO 400
  !
  ! Do numerical factorization of P matrix. ------------------------------
   300  NLU = NLU + 1
        IERpj = 0
        DO i = 1 , N
           Rtem(i) = 0.0D0
        ENDDO
        CALL CDRV(N,Iwk(IPR),Iwk(IPC),Iwk(IPIc),Iwk(IPIan),Iwk(IPJan),    &
       &          Wk(IPA),Rtem,Rtem,NSP,Iwk(IPIsp),Wk(IPRsp),IESp,2,IYS)
        IF ( IYS.EQ.0 ) RETURN
        imul = (IYS-1)/N
        IERpj = -2
        IF ( imul.EQ.8 ) IERpj = 1
        IF ( imul.EQ.10 ) IERpj = -1
        RETURN
  ! Error return for IRES = 2 or IRES = 3 return from RES. ---------------
   400  IERpj = ires
  !----------------------- End of Subroutine DPRJIS ----------------------
        END

end module ode_lsoda_aux1
