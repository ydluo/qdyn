module fftsg

! This module is a f90 wrapper for Takuya Ooura's FFT routines 
! General Purpose FFT (Fast Fourier/Cosine/Sine Transform) Package
! Copyright(C) 1996-2001 Takuya OOURA
! email: ooura@kurims.kyoto-u.ac.jp
! download: http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html

  private 

 ! working arrays for Ooura's fft
  type OouraFFT_type
    private
    integer :: n=-1, nw=-1, n1=-1, n2=-1
    integer, allocatable, dimension(:) :: iw
    double precision, allocatable, dimension(:) :: rw, tw
  end type OouraFFT_type

  public :: cdft,rdft,ddct,ddst,dfct,dfst,my_rdft,my_rdft2,conj2d,OouraFFT_type

contains

! Fast Fourier/Cosine/Sine Transform
!     dimension   :one
!     data length :power of 2
!     decimation  :frequency
!     radix       :split-radix
!     data        :inplace
!     table       :use
! subroutines
!     cdft: Complex Discrete Fourier Transform
!     rdft: Real Discrete Fourier Transform
!     ddct: Discrete Cosine Transform
!     ddst: Discrete Sine Transform
!     dfct: Cosine Transform of RDFT (Real Symmetric DFT)
!     dfst: Sine Transform of RDFT (Real Anti-symmetric DFT)
!
!
! -------- Complex DFT (Discrete Fourier Transform) --------
!     [definition]
!         <case1>
!             X(k) = sum_j=0^n-1 x(j)*exp(2*pi*i*j*k/n), 0<=k<n
!         <case2>
!             X(k) = sum_j=0^n-1 x(j)*exp(-2*pi*i*j*k/n), 0<=k<n
!         (notes: sum_j=0^n-1 is a summation from j=0 to n-1)
!     [usage]
!         <case1>
!             ip(0) = 0  ! first time only
!             call cdft(2*n, 1, a, ip, w)
!         <case2>
!             ip(0) = 0  ! first time only
!             call cdft(2*n, -1, a, ip, w)
!     [parameters]
!         2*n          :data length (integer)
!                       n >= 1, n = power of 2
!         a(0:2*n-1)   :input/output data (real*8)
!                       input data
!                           a(2*j) = Re(x(j)), 
!                           a(2*j+1) = Im(x(j)), 0<=j<n
!                       output data
!                           a(2*k) = Re(X(k)), 
!                           a(2*k+1) = Im(X(k)), 0<=k<n
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n)
!                       strictly, 
!                       length of ip >= 
!                           2+2**(int(log(n+0.5)/log(2.0))/2).
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n/2-1)   :cos/sin table (real*8)
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             call cdft(2*n, -1, a, ip, w)
!         is 
!             call cdft(2*n, 1, a, ip, w)
!             do j = 0, 2 * n - 1
!                 a(j) = a(j) / n
!             end do
!         .
!
!
! -------- Real DFT / Inverse of Real DFT --------
!     [definition]
!         <case1> RDFT
!             R(k) = sum_j=0^n-1 a(j)*cos(2*pi*j*k/n), 0<=k<=n/2
!             I(k) = sum_j=0^n-1 a(j)*sin(2*pi*j*k/n), 0<k<n/2
!         <case2> IRDFT (excluding scale)
!             a(k) = (R(0) + R(n/2)*cos(pi*k))/2 + 
!                    sum_j=1^n/2-1 R(j)*cos(2*pi*j*k/n) + 
!                    sum_j=1^n/2-1 I(j)*sin(2*pi*j*k/n), 0<=k<n
!     [usage]
!         <case1>
!             ip(0) = 0  ! first time only
!             call rdft(n, 1, a, ip, w)
!         <case2>
!             ip(0) = 0  ! first time only
!             call rdft(n, -1, a, ip, w)
!     [parameters]
!         n            :data length (integer)
!                       n >= 2, n = power of 2
!         a(0:n-1)     :input/output data (real*8)
!                       <case1>
!                           output data
!                               a(2*k) = R(k), 0<=k<n/2
!                               a(2*k+1) = I(k), 0<k<n/2
!                               a(1) = R(n/2)
!                       <case2>
!                           input data
!                               a(2*j) = R(j), 0<=j<n/2
!                               a(2*j+1) = I(j), 0<j<n/2
!                               a(1) = R(n/2)
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n/2)
!                       strictly, 
!                       length of ip >= 
!                           2+2**(int(log(n/2+0.5)/log(2.0))/2).
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n/2-1)   :cos/sin table (real*8)
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             call rdft(n, 1, a, ip, w)
!         is 
!             call rdft(n, -1, a, ip, w)
!             do j = 0, n - 1
!                 a(j) = a(j) * 2 / n
!             end do
!         .
!
!
! -------- DCT (Discrete Cosine Transform) / Inverse of DCT --------
!     [definition]
!         <case1> IDCT (excluding scale)
!             C(k) = sum_j=0^n-1 a(j)*cos(pi*j*(k+1/2)/n), 0<=k<n
!         <case2> DCT
!             C(k) = sum_j=0^n-1 a(j)*cos(pi*(j+1/2)*k/n), 0<=k<n
!     [usage]
!         <case1>
!             ip(0) = 0  ! first time only
!             call ddct(n, 1, a, ip, w)
!         <case2>
!             ip(0) = 0  ! first time only
!             call ddct(n, -1, a, ip, w)
!     [parameters]
!         n            :data length (integer)
!                       n >= 2, n = power of 2
!         a(0:n-1)     :input/output data (real*8)
!                       output data
!                           a(k) = C(k), 0<=k<n
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n/2)
!                       strictly, 
!                       length of ip >= 
!                           2+2**(int(log(n/2+0.5)/log(2.0))/2).
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n*5/4-1) :cos/sin table (real*8)
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             call ddct(n, -1, a, ip, w)
!         is 
!             a(0) = a(0) / 2
!             call ddct(n, 1, a, ip, w)
!             do j = 0, n - 1
!                 a(j) = a(j) * 2 / n
!             end do
!         .
!
!
! -------- DST (Discrete Sine Transform) / Inverse of DST --------
!     [definition]
!         <case1> IDST (excluding scale)
!             S(k) = sum_j=1^n A(j)*sin(pi*j*(k+1/2)/n), 0<=k<n
!         <case2> DST
!             S(k) = sum_j=0^n-1 a(j)*sin(pi*(j+1/2)*k/n), 0<k<=n
!     [usage]
!         <case1>
!             ip(0) = 0  ! first time only
!             call ddst(n, 1, a, ip, w)
!         <case2>
!             ip(0) = 0  ! first time only
!             call ddst(n, -1, a, ip, w)
!     [parameters]
!         n            :data length (integer)
!                       n >= 2, n = power of 2
!         a(0:n-1)     :input/output data (real*8)
!                       <case1>
!                           input data
!                               a(j) = A(j), 0<j<n
!                               a(0) = A(n)
!                           output data
!                               a(k) = S(k), 0<=k<n
!                       <case2>
!                           output data
!                               a(k) = S(k), 0<k<n
!                               a(0) = S(n)
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n/2)
!                       strictly, 
!                       length of ip >= 
!                           2+2**(int(log(n/2+0.5)/log(2.0))/2).
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n*5/4-1) :cos/sin table (real*8)
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             call ddst(n, -1, a, ip, w)
!         is 
!             a(0) = a(0) / 2
!             call ddst(n, 1, a, ip, w)
!             do j = 0, n - 1
!                 a(j) = a(j) * 2 / n
!             end do
!         .
!
!
! -------- Cosine Transform of RDFT (Real Symmetric DFT) --------
!     [definition]
!         C(k) = sum_j=0^n a(j)*cos(pi*j*k/n), 0<=k<=n
!     [usage]
!         ip(0) = 0  ! first time only
!         call dfct(n, a, t, ip, w)
!     [parameters]
!         n            :data length - 1 (integer)
!                       n >= 2, n = power of 2
!         a(0:n)       :input/output data (real*8)
!                       output data
!                           a(k) = C(k), 0<=k<=n
!         t(0:n/2)     :work area (real*8)
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n/4)
!                       strictly, 
!                       length of ip >= 
!                           2+2**(int(log(n/4+0.5)/log(2.0))/2).
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n*5/8-1) :cos/sin table (real*8)
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             a(0) = a(0) / 2
!             a(n) = a(n) / 2
!             call dfct(n, a, t, ip, w)
!         is 
!             a(0) = a(0) / 2
!             a(n) = a(n) / 2
!             call dfct(n, a, t, ip, w)
!             do j = 0, n
!                 a(j) = a(j) * 2 / n
!             end do
!         .
!
!
! -------- Sine Transform of RDFT (Real Anti-symmetric DFT) --------
!     [definition]
!         S(k) = sum_j=1^n-1 a(j)*sin(pi*j*k/n), 0<k<n
!     [usage]
!         ip(0) = 0  ! first time only
!         call dfst(n, a, t, ip, w)
!     [parameters]
!         n            :data length + 1 (integer)
!                       n >= 2, n = power of 2
!         a(0:n-1)     :input/output data (real*8)
!                       output data
!                           a(k) = S(k), 0<k<n
!                       (a(0) is used for work area)
!         t(0:n/2-1)   :work area (real*8)
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n/4)
!                       strictly, 
!                       length of ip >= 
!                           2+2**(int(log(n/4+0.5)/log(2.0))/2).
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n*5/8-1) :cos/sin table (real*8)
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             call dfst(n, a, t, ip, w)
!         is 
!             call dfst(n, a, t, ip, w)
!             do j = 1, n - 1
!                 a(j) = a(j) * 2 / n
!             end do
!         .
!
!
! Appendix :
!     The cos/sin table is recalculated when the larger table required.
!     w() and ip() are compatible with all routines.
!
!

! One-dimensional real FFT
subroutine my_rdft(isgn, a, m_fft)

  integer, intent(in) :: isgn
  real*8, intent(inout) :: a(:)
  type (OouraFFT_type), intent(inout) :: m_fft

  integer :: n

  n = size(a)
  if (n < 2**nextpow2(n)) stop 'FATAL ERROR in my_rdft: length must be a power of 2'

  ! Initialize or re-initialize working arrays if needed
  ! Handle 2 possible states:
  !   1. m_fft has never been intialized yet. This is indicated by m_fft%n<0 
  !      (by default m_fft%n= -1, see the definition of OouraFFT_type)
  !   2. m_fft has been initialized before but for a different length n.
  if (m_fft%n /= n) then
    if (m_fft%n >0) deallocate(m_fft%iw, m_fft%rw)
    m_fft%n = n 
    m_fft%nw = 2+ceiling(dsqrt(dble(m_fft%n/2.d0)))
    allocate (m_fft%iw(0:m_fft%nw-1))
    allocate (m_fft%rw(0:m_fft%n/2-1))
    m_fft%iw = 0
    m_fft%rw = 0d0
  endif
  
  ! compute fft
  call rdft(n,isgn,a,m_fft%iw,m_fft%rw)

end subroutine my_rdft

! Two-dimensional real FFT
subroutine my_rdft2(isgn, a, m_fft)

  integer, intent(in) :: isgn
  real*8, intent(inout) :: a(:,:)
  type(OouraFFT_type), intent(inout) :: m_fft

  integer, dimension(2) :: shp
  integer :: n1, n2, n

  ! Initialize working arrays with same scheme as before
  ! See comments in Ooura rdft2 code for allocation sizes
  shp = shape(a)
  n1 = shp(1)
  n2 = shp(2)
  if (n1 < 2**nextpow2(n1)) stop 'FATAL ERROR in my_rdft2: 1st size must be a power of 2'
  if (n2 < 2**nextpow2(n2)) stop 'FATAL ERROR in my_rdft2: 2nd size must be a power of 2'
  if (m_fft%n1 /= n1 .or. m_fft%n2 /= n2) then
    if (m_fft%n1 > 0) deallocate(m_fft%iw, m_fft%rw, m_fft%tw)
    m_fft%n1 = n1
    m_fft%n2 = n2
    n = max(n1/2, n2)
    m_fft%nw = 2 + ceiling(dsqrt(dble(n)))
    allocate(m_fft%tw(8*n2))
    allocate(m_fft%iw(m_fft%nw))
    n = max(n1/4, n2/2) + n1/4
    allocate(m_fft%rw(n))
  end if

  ! Compute FFT
  call rdft2d(n1, n1, n2, isgn, a, m_fft%tw, m_fft%iw, m_fft%rw)

end subroutine my_rdft2

! same as matlab's nextpow2
integer function nextpow2(n)
  
  integer, intent(in) :: n
  
  integer :: n2
  
  n = abs(n)
  nextpow2 = 1
  n2 = 2
  do while( n > n2 )
    nextpow2 = nextpow2 + 1
    n2 = 2*n2
  end do 
  
end function nextpow2
    
! Complex conjugate of 2D real FFT array
! Output array MUST already have the correct shape
subroutine conj2d(a, b)

  real*8, intent(in) :: a(:,:)
  real*8, intent(inout) :: b(:,:)

  integer :: n2
  integer, dimension(2) :: shp

  ! Determine size of input array and allocate output
  shp = shape(a)
  n2 = shp(2)

  ! Assign values to output
  b = a
  b(1,n2/2+2:) = -1.0d0 * b(1,n2/2+2:)
  b(2,2:n2/2)  = -1.0d0 * b(2,2:n2/2)
  b(4::2,:)    = -1.0d0 * b(4::2,:)

end subroutine conj2d
  

!_________________________________________________________

      subroutine cdft(n, isgn, a, ip, w)
      integer n, isgn, ip(0 : *), nw
      real*8 a(0 : n - 1), w(0 : *)
      nw = ip(0)
      if (n .gt. 4 * nw) then
          nw = n / 4
          call makewt(nw, ip, w)
      end if
      if (isgn .ge. 0) then
          call cftfsub(n, a, ip, nw, w)
      else
          call cftbsub(n, a, ip, nw, w)
      end if
      end subroutine cdft
!
      subroutine rdft(n, isgn, a, ip, w)
      integer n, isgn, ip(0 : *), nw, nc
      real*8 a(0 : n - 1), w(0 : *), xi
      nw = ip(0)
      if (n .gt. 4 * nw) then
          nw = n / 4
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. 4 * nc) then
          nc = n / 4
          call makect(nc, ip, w(nw))
      end if
      if (isgn .ge. 0) then
          if (n .gt. 4) then
              call cftfsub(n, a, ip, nw, w)
              call rftfsub(n, a, nc, w(nw))
          else if (n .eq. 4) then
              call cftfsub(n, a, ip, nw, w)
          end if
          xi = a(0) - a(1)
          a(0) = a(0) + a(1)
          a(1) = xi
      else
          a(1) = 0.5d0 * (a(0) - a(1))
          a(0) = a(0) - a(1)
          if (n .gt. 4) then
              call rftbsub(n, a, nc, w(nw))
              call cftbsub(n, a, ip, nw, w)
          else if (n .eq. 4) then
              call cftbsub(n, a, ip, nw, w)
          end if
      end if
      end subroutine rdft
!
      subroutine ddct(n, isgn, a, ip, w)
      integer n, isgn, ip(0 : *), j, nw, nc
      real*8 a(0 : n - 1), w(0 : *), xr
      nw = ip(0)
      if (n .gt. 4 * nw) then
          nw = n / 4
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. nc) then
          nc = n
          call makect(nc, ip, w(nw))
      end if
      if (isgn .lt. 0) then
          xr = a(n - 1)
          do j = n - 2, 2, -2
              a(j + 1) = a(j) - a(j - 1)
              a(j) = a(j) + a(j - 1)
          end do
          a(1) = a(0) - xr
          a(0) = a(0) + xr
          if (n .gt. 4) then
              call rftbsub(n, a, nc, w(nw))
              call cftbsub(n, a, ip, nw, w)
          else if (n .eq. 4) then
              call cftbsub(n, a, ip, nw, w)
          end if
      end if
      call dctsub(n, a, nc, w(nw))
      if (isgn .ge. 0) then
          if (n .gt. 4) then
              call cftfsub(n, a, ip, nw, w)
              call rftfsub(n, a, nc, w(nw))
          else if (n .eq. 4) then
              call cftfsub(n, a, ip, nw, w)
          end if
          xr = a(0) - a(1)
          a(0) = a(0) + a(1)
          do j = 2, n - 2, 2
              a(j - 1) = a(j) - a(j + 1)
              a(j) = a(j) + a(j + 1)
          end do
          a(n - 1) = xr
      end if
      end subroutine ddct
!
      subroutine ddst(n, isgn, a, ip, w)
      integer n, isgn, ip(0 : *), j, nw, nc
      real*8 a(0 : n - 1), w(0 : *), xr
      nw = ip(0)
      if (n .gt. 4 * nw) then
          nw = n / 4
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. nc) then
          nc = n
          call makect(nc, ip, w(nw))
      end if
      if (isgn .lt. 0) then
          xr = a(n - 1)
          do j = n - 2, 2, -2
              a(j + 1) = -a(j) - a(j - 1)
              a(j) = a(j) - a(j - 1)
          end do
          a(1) = a(0) + xr
          a(0) = a(0) - xr
          if (n .gt. 4) then
              call rftbsub(n, a, nc, w(nw))
              call cftbsub(n, a, ip, nw, w)
          else if (n .eq. 4) then
              call cftbsub(n, a, ip, nw, w)
          end if
      end if
      call dstsub(n, a, nc, w(nw))
      if (isgn .ge. 0) then
          if (n .gt. 4) then
              call cftfsub(n, a, ip, nw, w)
              call rftfsub(n, a, nc, w(nw))
          else if (n .eq. 4) then
              call cftfsub(n, a, ip, nw, w)
          end if
          xr = a(0) - a(1)
          a(0) = a(0) + a(1)
          do j = 2, n - 2, 2
              a(j - 1) = -a(j) - a(j + 1)
              a(j) = a(j) - a(j + 1)
          end do
          a(n - 1) = -xr
      end if
      end subroutine ddst
!
      subroutine dfct(n, a, t, ip, w)
      integer n, ip(0 : *), j, k, l, m, mh, nw, nc
      real*8 a(0 : n), t(0 : n / 2), w(0 : *), xr, xi, yr, yi
      nw = ip(0)
      if (n .gt. 8 * nw) then
          nw = n / 8
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. 2 * nc) then
          nc = n / 2
          call makect(nc, ip, w(nw))
      end if
      m = n / 2
      yi = a(m)
      xi = a(0) + a(n)
      a(0) = a(0) - a(n)
      t(0) = xi - yi
      t(m) = xi + yi
      if (n .gt. 2) then
          mh = m / 2
          do j = 1, mh - 1
              k = m - j
              xr = a(j) - a(n - j)
              xi = a(j) + a(n - j)
              yr = a(k) - a(n - k)
              yi = a(k) + a(n - k)
              a(j) = xr
              a(k) = yr
              t(j) = xi - yi
              t(k) = xi + yi
          end do
          t(mh) = a(mh) + a(n - mh)
          a(mh) = a(mh) - a(n - mh)
          call dctsub(m, a, nc, w(nw))
          if (m .gt. 4) then
              call cftfsub(m, a, ip, nw, w)
              call rftfsub(m, a, nc, w(nw))
          else if (m .eq. 4) then
              call cftfsub(m, a, ip, nw, w)
          end if
          a(n - 1) = a(0) - a(1)
          a(1) = a(0) + a(1)
          do j = m - 2, 2, -2
              a(2 * j + 1) = a(j) + a(j + 1)
              a(2 * j - 1) = a(j) - a(j + 1)
          end do
          l = 2
          m = mh
          do while (m .ge. 2)
              call dctsub(m, t, nc, w(nw))
              if (m .gt. 4) then
                  call cftfsub(m, t, ip, nw, w)
                  call rftfsub(m, t, nc, w(nw))
              else if (m .eq. 4) then
                  call cftfsub(m, t, ip, nw, w)
              end if
              a(n - l) = t(0) - t(1)
              a(l) = t(0) + t(1)
              k = 0
              do j = 2, m - 2, 2
                  k = k + 4 * l
                  a(k - l) = t(j) - t(j + 1)
                  a(k + l) = t(j) + t(j + 1)
              end do
              l = 2 * l
              mh = m / 2
              do j = 0, mh - 1
                  k = m - j
                  t(j) = t(m + k) - t(m + j)
                  t(k) = t(m + k) + t(m + j)
              end do
              t(mh) = t(m + mh)
              m = mh
          end do
          a(l) = t(0)
          a(n) = t(2) - t(1)
          a(0) = t(2) + t(1)
      else
          a(1) = a(0)
          a(2) = t(0)
          a(0) = t(1)
      end if
      end subroutine dfct
!
      subroutine dfst(n, a, t, ip, w)
      integer n, ip(0 : *), j, k, l, m, mh, nw, nc
      real*8 a(0 : n - 1), t(0 : n / 2 - 1), w(0 : *), xr, xi, yr, yi
      nw = ip(0)
      if (n .gt. 8 * nw) then
          nw = n / 8
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. 2 * nc) then
          nc = n / 2
          call makect(nc, ip, w(nw))
      end if
      if (n .gt. 2) then
          m = n / 2
          mh = m / 2
          do j = 1, mh - 1
              k = m - j
              xr = a(j) + a(n - j)
              xi = a(j) - a(n - j)
              yr = a(k) + a(n - k)
              yi = a(k) - a(n - k)
              a(j) = xr
              a(k) = yr
              t(j) = xi + yi
              t(k) = xi - yi
          end do
          t(0) = a(mh) - a(n - mh)
          a(mh) = a(mh) + a(n - mh)
          a(0) = a(m)
          call dstsub(m, a, nc, w(nw))
          if (m .gt. 4) then
              call cftfsub(m, a, ip, nw, w)
              call rftfsub(m, a, nc, w(nw))
          else if (m .eq. 4) then
              call cftfsub(m, a, ip, nw, w)
          end if
          a(n - 1) = a(1) - a(0)
          a(1) = a(0) + a(1)
          do j = m - 2, 2, -2
              a(2 * j + 1) = a(j) - a(j + 1)
              a(2 * j - 1) = -a(j) - a(j + 1)
          end do
          l = 2
          m = mh
          do while (m .ge. 2)
              call dstsub(m, t, nc, w(nw))
              if (m .gt. 4) then
                  call cftfsub(m, t, ip, nw, w)
                  call rftfsub(m, t, nc, w(nw))
              else if (m .eq. 4) then
                  call cftfsub(m, t, ip, nw, w)
              end if
              a(n - l) = t(1) - t(0)
              a(l) = t(0) + t(1)
              k = 0
              do j = 2, m - 2, 2
                  k = k + 4 * l
                  a(k - l) = -t(j) - t(j + 1)
                  a(k + l) = t(j) - t(j + 1)
              end do
              l = 2 * l
              mh = m / 2
              do j = 1, mh - 1
                  k = m - j
                  t(j) = t(m + k) + t(m + j)
                  t(k) = t(m + k) - t(m + j)
              end do
              t(0) = t(m + mh)
              m = mh
          end do
          a(l) = t(0)
      end if
      a(0) = 0
      end subroutine dfst

! -------- Two-dimensional routines -------

      subroutine cdft2d(n1max, n1, n2, isgn, a, t, ip, w)
      integer n1max, n1, n2, isgn, ip(0 : *), n, j
      real*8 a(0 : n1max - 1, 0 : n2 - 1), t(0 : 8 * n2 - 1), w(0 : *)
      n = max(n1, 2 * n2)
      if (n .gt. 4 * ip(0)) then
          call makewt(n / 4, ip, w)
      end if
      do j = 0, n2 - 1
          call cdft(n1, isgn, a(0, j), ip, w)
      end do
      call cdft2d_sub(n1max, n1, n2, isgn, a, t, ip, w)
      end subroutine cdft2d
!
      subroutine rdft2d(n1max, n1, n2, isgn, a, t, ip, w)
      integer n1max, n1, n2, isgn, ip(0 : *), n, nw, nc, j
      real*8 a(0 : n1max - 1, 0 : n2 - 1), t(0 : 8 * n2 - 1), w(0 : *)
      n = max(n1, 2 * n2)
      nw = ip(0)
      if (n .gt. 4 * nw) then
          nw = n / 4
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n1 .gt. 4 * nc) then
          nc = n1 / 4
          call makect(nc, ip, w(nw))
      end if
      if (isgn .lt. 0) then
          call rdft2d_sub(n1max, n1, n2, isgn, a)
          call cdft2d_sub(n1max, n1, n2, isgn, a, t, ip, w)
      end if
      do j = 0, n2 - 1
          call rdft(n1, isgn, a(0, j), ip, w)
      end do
      if (isgn .ge. 0) then
          call cdft2d_sub(n1max, n1, n2, isgn, a, t, ip, w)
          call rdft2d_sub(n1max, n1, n2, isgn, a)
      end if
      end subroutine rdft2d
!
      subroutine rdft2dsort(n1max, n1, n2, isgn, a)
      integer n1max, n1, n2, isgn, n2h, j
      real*8 a(0 : n1max - 1, 0 : n2 - 1), x, y
      n2h = n2 / 2
      if (isgn .lt. 0) then
          do j = n2h + 1, n2 - 1
              a(0, j) = a(n1 + 1, j)
              a(1, j) = a(n1, j)
          end do
          a(1, 0) = a(n1, 0)
          a(1, n2h) = a(n1, n2h)
      else
          do j = n2h + 1, n2 - 1
              y = a(0, j)
              x = a(1, j)
              a(n1, j) = x
              a(n1 + 1, j) = y
              a(n1, n2 - j) = x
              a(n1 + 1, n2 - j) = -y
              a(0, j) = a(0, n2 - j)
              a(1, j) = -a(1, n2 - j)
          end do
          a(n1, 0) = a(1, 0)
          a(n1 + 1, 0) = 0
          a(1, 0) = 0
          a(n1, n2h) = a(1, n2h)
          a(n1 + 1, n2h) = 0
          a(1, n2h) = 0
      end if
      end subroutine rdft2dsort
!
      subroutine ddct2d(n1max, n1, n2, isgn, a, t, ip, w)
      integer n1max, n1, n2, isgn, ip(0 : *), n, nw, nc, j
      real*8 a(0 : n1max - 1, 0 : n2 - 1), t(0 : 4 * n2 - 1), w(0 : *)
      n = max(n1, n2)
      nw = ip(0)
      if (n .gt. 4 * nw) then
          nw = n / 4
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. nc) then
          nc = n
          call makect(nc, ip, w(nw))
      end if
      do j = 0, n2 - 1
          call ddct(n1, isgn, a(0, j), ip, w)
      end do
      call ddxt2d_sub(n1max, n1, n2, 0, isgn, a, t, ip, w)
      end subroutine ddct2d
!
      subroutine ddst2d(n1max, n1, n2, isgn, a, t, ip, w)
      integer n1max, n1, n2, isgn, ip(0 : *), n, nw, nc, j
      real*8 a(0 : n1max - 1, 0 : n2 - 1), t(0 : 4 * n2 - 1), w(0 : *)
      n = max(n1, n2)
      nw = ip(0)
      if (n .gt. 4 * nw) then
          nw = n / 4
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. nc) then
          nc = n
          call makect(nc, ip, w(nw))
      end if
      do j = 0, n2 - 1
          call ddst(n1, isgn, a(0, j), ip, w)
      end do
      call ddxt2d_sub(n1max, n1, n2, 1, isgn, a, t, ip, w)
      end subroutine ddst2d
!
! -------- child routines --------
!
      subroutine cdft2d_sub(n1max, n1, n2, isgn, a, t, ip, w)
      integer n1max, n1, n2, isgn, ip(0 : *), i, j
      real*8 a(0 : n1max - 1, 0 : n2 - 1), t(0 : 8 * n2 - 1), w(0 : *)
      if (n1 .gt. 4) then
          do i = 0, n1 - 8, 8
              do j = 0, n2 - 1
                  t(2 * j) = a(i, j)
                  t(2 * j + 1) = a(i + 1, j)
                  t(2 * n2 + 2 * j) = a(i + 2, j)
                  t(2 * n2 + 2 * j + 1) = a(i + 3, j)
                  t(4 * n2 + 2 * j) = a(i + 4, j)
                  t(4 * n2 + 2 * j + 1) = a(i + 5, j)
                  t(6 * n2 + 2 * j) = a(i + 6, j)
                  t(6 * n2 + 2 * j + 1) = a(i + 7, j)
              end do
              call cdft(2 * n2, isgn, t, ip, w)
              call cdft(2 * n2, isgn, t(2 * n2), ip, w)
              call cdft(2 * n2, isgn, t(4 * n2), ip, w)
              call cdft(2 * n2, isgn, t(6 * n2), ip, w)
              do j = 0, n2 - 1
                  a(i, j) = t(2 * j)
                  a(i + 1, j) = t(2 * j + 1)
                  a(i + 2, j) = t(2 * n2 + 2 * j)
                  a(i + 3, j) = t(2 * n2 + 2 * j + 1)
                  a(i + 4, j) = t(4 * n2 + 2 * j)
                  a(i + 5, j) = t(4 * n2 + 2 * j + 1)
                  a(i + 6, j) = t(6 * n2 + 2 * j)
                  a(i + 7, j) = t(6 * n2 + 2 * j + 1)
              end do
          end do
      else if (n1 .eq. 4) then
          do j = 0, n2 - 1
              t(2 * j) = a(0, j)
              t(2 * j + 1) = a(1, j)
              t(2 * n2 + 2 * j) = a(2, j)
              t(2 * n2 + 2 * j + 1) = a(3, j)
          end do
          call cdft(2 * n2, isgn, t, ip, w)
          call cdft(2 * n2, isgn, t(2 * n2), ip, w)
          do j = 0, n2 - 1
              a(0, j) = t(2 * j)
              a(1, j) = t(2 * j + 1)
              a(2, j) = t(2 * n2 + 2 * j)
              a(3, j) = t(2 * n2 + 2 * j + 1)
          end do
      else if (n1 .eq. 2) then
          do j = 0, n2 - 1
              t(2 * j) = a(0, j)
              t(2 * j + 1) = a(1, j)
          end do
          call cdft(2 * n2, isgn, t, ip, w)
          do j = 0, n2 - 1
              a(0, j) = t(2 * j)
              a(1, j) = t(2 * j + 1)
          end do
      end if
      end subroutine cdft2d_sub
!
      subroutine rdft2d_sub(n1max, n1, n2, isgn, a)
      integer n1max, n1, n2, isgn, n2h, i, j
      real*8 a(0 : n1max - 1, 0 : n2 - 1), xi
      n2h = n2 / 2
      if (isgn .lt. 0) then
          do i = 1, n2h - 1
              j = n2 - i
              xi = a(0, i) - a(0, j)
              a(0, i) = a(0, i) + a(0, j)
              a(0, j) = xi
              xi = a(1, j) - a(1, i)
              a(1, i) = a(1, i) + a(1, j)
              a(1, j) = xi
          end do
      else
          do i = 1, n2h - 1
              j = n2 - i
              a(0, j) = 0.5d0 * (a(0, i) - a(0, j))
              a(0, i) = a(0, i) - a(0, j)
              a(1, j) = 0.5d0 * (a(1, i) + a(1, j))
              a(1, i) = a(1, i) - a(1, j)
          end do
      end if
      end subroutine rdft2d_sub
!
      subroutine ddxt2d_sub(n1max, n1, n2, ics, isgn, a, t, ip, w)
      integer n1max, n1, n2, ics, isgn, ip(0 : *), i, j
      real*8 a(0 : n1max - 1, 0 : n2 - 1), t(0 : 4 * n2 - 1), w(0 : *)
      if (n1 .gt. 2) then
          do i = 0, n1 - 4, 4
              do j = 0, n2 - 1
                  t(j) = a(i, j)
                  t(n2 + j) = a(i + 1, j)
                  t(2 * n2 + j) = a(i + 2, j)
                  t(3 * n2 + j) = a(i + 3, j)
              end do
              if (ics .eq. 0) then
                  call ddct(n2, isgn, t, ip, w)
                  call ddct(n2, isgn, t(n2), ip, w)
                  call ddct(n2, isgn, t(2 * n2), ip, w)
                  call ddct(n2, isgn, t(3 * n2), ip, w)
              else
                  call ddst(n2, isgn, t, ip, w)
                  call ddst(n2, isgn, t(n2), ip, w)
                  call ddst(n2, isgn, t(2 * n2), ip, w)
                  call ddst(n2, isgn, t(3 * n2), ip, w)
              end if
              do j = 0, n2 - 1
                  a(i, j) = t(j)
                  a(i + 1, j) = t(n2 + j)
                  a(i + 2, j) = t(2 * n2 + j)
                  a(i + 3, j) = t(3 * n2 + j)
              end do
          end do
      else if (n1 .eq. 2) then
          do j = 0, n2 - 1
              t(j) = a(0, j)
              t(n2 + j) = a(1, j)
          end do
          if (ics .eq. 0) then
              call ddct(n2, isgn, t, ip, w)
              call ddct(n2, isgn, t(n2), ip, w)
          else
              call ddst(n2, isgn, t, ip, w)
              call ddst(n2, isgn, t(n2), ip, w)
          end if
          do j = 0, n2 - 1
              a(0, j) = t(j)
              a(1, j) = t(n2 + j)
          end do
      end if
      end subroutine ddxt2d_sub

!
! -------- initializing routines --------
!
      subroutine makewt(nw, ip, w)
      integer nw, ip(0 : *), j, nwh, nw0, nw1
      real*8 w(0 : nw - 1), delta, wn4r, wk1r, wk1i, wk3r, wk3i
      ip(0) = nw
      ip(1) = 1
      if (nw .gt. 2) then
          nwh = nw / 2
          delta = atan(1.0d0) / nwh
          wn4r = cos(delta * nwh)
          w(0) = 1
          w(1) = wn4r
          if (nwh .eq. 4) then
              w(2) = cos(delta * 2)
              w(3) = sin(delta * 2)
          else if (nwh .gt. 4) then
              call makeipt(nw, ip)
              w(2) = 0.5d0 / cos(delta * 2)
              w(3) = 0.5d0 / cos(delta * 6)
              do j = 4, nwh - 4, 4
                  w(j) = cos(delta * j)
                  w(j + 1) = sin(delta * j)
                  w(j + 2) = cos(3 * delta * j)
                  w(j + 3) = -sin(3 * delta * j)
              end do
          end if
          nw0 = 0
          do while (nwh .gt. 2)
              nw1 = nw0 + nwh
              nwh = nwh / 2
              w(nw1) = 1
              w(nw1 + 1) = wn4r
              if (nwh .eq. 4) then
                  wk1r = w(nw0 + 4)
                  wk1i = w(nw0 + 5)
                  w(nw1 + 2) = wk1r
                  w(nw1 + 3) = wk1i
              else if (nwh .gt. 4) then
                  wk1r = w(nw0 + 4)
                  wk3r = w(nw0 + 6)
                  w(nw1 + 2) = 0.5d0 / wk1r
                  w(nw1 + 3) = 0.5d0 / wk3r
                  do j = 4, nwh - 4, 4
                      wk1r = w(nw0 + 2 * j)
                      wk1i = w(nw0 + 2 * j + 1)
                      wk3r = w(nw0 + 2 * j + 2)
                      wk3i = w(nw0 + 2 * j + 3)
                      w(nw1 + j) = wk1r
                      w(nw1 + j + 1) = wk1i
                      w(nw1 + j + 2) = wk3r
                      w(nw1 + j + 3) = wk3i
                  end do
              end if
              nw0 = nw1
          end do
      end if
      end subroutine makewt
!
      subroutine makeipt(nw, ip)
      integer nw, ip(0 : *), j, l, m, m2, p, q
      ip(2) = 0
      ip(3) = 16
      m = 2
      l = nw
      do while (l .gt. 32)
          m2 = 2 * m
          q = 8 * m2
          do j = m, m2 - 1
              p = 4 * ip(j)
              ip(m + j) = p
              ip(m2 + j) = p + q
          end do
          m = m2
          l = l / 4
      end do
      end subroutine makeipt
!
      subroutine makect(nc, ip, c)
      integer nc, ip(0 : *), j, nch
      real*8 c(0 : nc - 1), delta
      ip(1) = nc
      if (nc .gt. 1) then
          nch = nc / 2
          delta = atan(1.0d0) / nch
          c(0) = cos(delta * nch)
          c(nch) = 0.5d0 * c(0)
          do j = 1, nch - 1
              c(j) = 0.5d0 * cos(delta * j)
              c(nc - j) = 0.5d0 * sin(delta * j)
          end do
      end if
      end subroutine makect
!
! -------- child routines --------
!
      subroutine cftfsub(n, a, ip, nw, w)
      integer n, ip(0 : *), nw
      real*8 a(0 : n - 1), w(0 : nw - 1)
      if (n .gt. 8) then
          if (n .gt. 32) then
              call cftf1st(n, a, w(nw - n / 4))
              if (n .gt. 512) then
                  call cftrec4(n, a, nw, w)
              else if (n .gt. 128) then
                  call cftleaf(n, 1, a, nw, w)
              else
                  call cftfx41(n, a, nw, w)
              end if
              call bitrv2(n, ip, a)
          else if (n .eq. 32) then
              call cftf161(a, w(nw - 8))
              call bitrv216(a)
          else
              call cftf081(a, w)
              call bitrv208(a)
          end if
      else if (n .eq. 8) then
          call cftf040(a)
      else if (n .eq. 4) then
          call cftx020(a)
      end if
      end subroutine cftfsub
!
      subroutine cftbsub(n, a, ip, nw, w)
      integer n, ip(0 : *), nw
      real*8 a(0 : n - 1), w(0 : nw - 1)
      if (n .gt. 8) then
          if (n .gt. 32) then
              call cftb1st(n, a, w(nw - n / 4))
              if (n .gt. 512) then
                  call cftrec4(n, a, nw, w)
              else if (n .gt. 128) then
                  call cftleaf(n, 1, a, nw, w)
              else
                  call cftfx41(n, a, nw, w)
              end if
              call bitrv2conj(n, ip, a)
          else if (n .eq. 32) then
              call cftf161(a, w(nw - 8))
              call bitrv216neg(a)
          else
              call cftf081(a, w)
              call bitrv208neg(a)
          end if
      else if (n .eq. 8) then
          call cftb040(a)
      else if (n .eq. 4) then
          call cftx020(a)
      end if
      end subroutine cftbsub
!
      subroutine bitrv2(n, ip, a)
      integer n, ip(0 : *), j, j1, k, k1, l, m, nh, nm
      real*8 a(0 : n - 1), xr, xi, yr, yi
      m = 1
      l = n / 4
      do while (l .gt. 8)
          m = m * 2
          l = l / 4
      end do
      nh = n / 2
      nm = 4 * m
      if (l .eq. 8) then
          do k = 0, m - 1
              do j = 0, k - 1
                  j1 = 4 * j + 2 * ip(m + k)
                  k1 = 4 * k + 2 * ip(m + j)
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nm
                  k1 = k1 + 2 * nm
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nm
                  k1 = k1 - nm
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nm
                  k1 = k1 + 2 * nm
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nh
                  k1 = k1 + 2
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nm
                  k1 = k1 - 2 * nm
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nm
                  k1 = k1 + nm
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nm
                  k1 = k1 - 2 * nm
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + 2
                  k1 = k1 + nh
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nm
                  k1 = k1 + 2 * nm
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nm
                  k1 = k1 - nm
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nm
                  k1 = k1 + 2 * nm
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nh
                  k1 = k1 - 2
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nm
                  k1 = k1 - 2 * nm
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nm
                  k1 = k1 + nm
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nm
                  k1 = k1 - 2 * nm
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
              k1 = 4 * k + 2 * ip(m + k)
              j1 = k1 + 2
              k1 = k1 + nh
              xr = a(j1)
              xi = a(j1 + 1)
              yr = a(k1)
              yi = a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
              j1 = j1 + nm
              k1 = k1 + 2 * nm
              xr = a(j1)
              xi = a(j1 + 1)
              yr = a(k1)
              yi = a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
              j1 = j1 + nm
              k1 = k1 - nm
              xr = a(j1)
              xi = a(j1 + 1)
              yr = a(k1)
              yi = a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
              j1 = j1 - 2
              k1 = k1 - nh
              xr = a(j1)
              xi = a(j1 + 1)
              yr = a(k1)
              yi = a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
              j1 = j1 + nh + 2
              k1 = k1 + nh + 2
              xr = a(j1)
              xi = a(j1 + 1)
              yr = a(k1)
              yi = a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
              j1 = j1 - nh + nm
              k1 = k1 + 2 * nm - 2
              xr = a(j1)
              xi = a(j1 + 1)
              yr = a(k1)
              yi = a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
          end do
      else
          do k = 0, m - 1
              do j = 0, k - 1
                  j1 = 4 * j + ip(m + k)
                  k1 = 4 * k + ip(m + j)
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nm
                  k1 = k1 + nm
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nh
                  k1 = k1 + 2
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nm
                  k1 = k1 - nm
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + 2
                  k1 = k1 + nh
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nm
                  k1 = k1 + nm
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nh
                  k1 = k1 - 2
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nm
                  k1 = k1 - nm
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
              k1 = 4 * k + ip(m + k)
              j1 = k1 + 2
              k1 = k1 + nh
              xr = a(j1)
              xi = a(j1 + 1)
              yr = a(k1)
              yi = a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
              j1 = j1 + nm
              k1 = k1 + nm
              xr = a(j1)
              xi = a(j1 + 1)
              yr = a(k1)
              yi = a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
          end do
      end if
      end subroutine bitrv2
!
      subroutine bitrv2conj(n, ip, a)
      integer n, ip(0 : *), j, j1, k, k1, l, m, nh, nm
      real*8 a(0 : n - 1), xr, xi, yr, yi
      m = 1
      l = n / 4
      do while (l .gt. 8)
          m = m * 2
          l = l / 4
      end do
      nh = n / 2
      nm = 4 * m
      if (l .eq. 8) then
          do k = 0, m - 1
              do j = 0, k - 1
                  j1 = 4 * j + 2 * ip(m + k)
                  k1 = 4 * k + 2 * ip(m + j)
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nm
                  k1 = k1 + 2 * nm
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nm
                  k1 = k1 - nm
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nm
                  k1 = k1 + 2 * nm
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nh
                  k1 = k1 + 2
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nm
                  k1 = k1 - 2 * nm
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nm
                  k1 = k1 + nm
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nm
                  k1 = k1 - 2 * nm
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + 2
                  k1 = k1 + nh
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nm
                  k1 = k1 + 2 * nm
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nm
                  k1 = k1 - nm
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nm
                  k1 = k1 + 2 * nm
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nh
                  k1 = k1 - 2
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nm
                  k1 = k1 - 2 * nm
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nm
                  k1 = k1 + nm
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nm
                  k1 = k1 - 2 * nm
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
              k1 = 4 * k + 2 * ip(m + k)
              j1 = k1 + 2
              k1 = k1 + nh
              a(j1 - 1) = -a(j1 - 1)
              xr = a(j1)
              xi = -a(j1 + 1)
              yr = a(k1)
              yi = -a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
              a(k1 + 3) = -a(k1 + 3)
              j1 = j1 + nm
              k1 = k1 + 2 * nm
              xr = a(j1)
              xi = -a(j1 + 1)
              yr = a(k1)
              yi = -a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
              j1 = j1 + nm
              k1 = k1 - nm
              xr = a(j1)
              xi = -a(j1 + 1)
              yr = a(k1)
              yi = -a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
              j1 = j1 - 2
              k1 = k1 - nh
              xr = a(j1)
              xi = -a(j1 + 1)
              yr = a(k1)
              yi = -a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
              j1 = j1 + nh + 2
              k1 = k1 + nh + 2
              xr = a(j1)
              xi = -a(j1 + 1)
              yr = a(k1)
              yi = -a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
              j1 = j1 - nh + nm
              k1 = k1 + 2 * nm - 2
              a(j1 - 1) = -a(j1 - 1)
              xr = a(j1)
              xi = -a(j1 + 1)
              yr = a(k1)
              yi = -a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
              a(k1 + 3) = -a(k1 + 3)
          end do
      else
          do k = 0, m - 1
              do j = 0, k - 1
                  j1 = 4 * j + ip(m + k)
                  k1 = 4 * k + ip(m + j)
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nm
                  k1 = k1 + nm
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nh
                  k1 = k1 + 2
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nm
                  k1 = k1 - nm
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + 2
                  k1 = k1 + nh
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + nm
                  k1 = k1 + nm
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nh
                  k1 = k1 - 2
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 - nm
                  k1 = k1 - nm
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
              k1 = 4 * k + ip(m + k)
              j1 = k1 + 2
              k1 = k1 + nh
              a(j1 - 1) = -a(j1 - 1)
              xr = a(j1)
              xi = -a(j1 + 1)
              yr = a(k1)
              yi = -a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
              a(k1 + 3) = -a(k1 + 3)
              j1 = j1 + nm
              k1 = k1 + nm
              a(j1 - 1) = -a(j1 - 1)
              xr = a(j1)
              xi = -a(j1 + 1)
              yr = a(k1)
              yi = -a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
              a(k1 + 3) = -a(k1 + 3)
          end do
      end if
      end subroutine bitrv2conj
!
      subroutine bitrv216(a)
      real*8 a(0 : 31), x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i
      real*8 x5r, x5i, x7r, x7i, x8r, x8i, x10r, x10i
      real*8 x11r, x11i, x12r, x12i, x13r, x13i, x14r, x14i
      x1r = a(2)
      x1i = a(3)
      x2r = a(4)
      x2i = a(5)
      x3r = a(6)
      x3i = a(7)
      x4r = a(8)
      x4i = a(9)
      x5r = a(10)
      x5i = a(11)
      x7r = a(14)
      x7i = a(15)
      x8r = a(16)
      x8i = a(17)
      x10r = a(20)
      x10i = a(21)
      x11r = a(22)
      x11i = a(23)
      x12r = a(24)
      x12i = a(25)
      x13r = a(26)
      x13i = a(27)
      x14r = a(28)
      x14i = a(29)
      a(2) = x8r
      a(3) = x8i
      a(4) = x4r
      a(5) = x4i
      a(6) = x12r
      a(7) = x12i
      a(8) = x2r
      a(9) = x2i
      a(10) = x10r
      a(11) = x10i
      a(14) = x14r
      a(15) = x14i
      a(16) = x1r
      a(17) = x1i
      a(20) = x5r
      a(21) = x5i
      a(22) = x13r
      a(23) = x13i
      a(24) = x3r
      a(25) = x3i
      a(26) = x11r
      a(27) = x11i
      a(28) = x7r
      a(29) = x7i
      end subroutine bitrv216
!
      subroutine bitrv216neg(a)
      real*8 a(0 : 31), x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i
      real*8 x5r, x5i, x6r, x6i, x7r, x7i, x8r, x8i
      real*8 x9r, x9i, x10r, x10i, x11r, x11i, x12r, x12i
      real*8 x13r, x13i, x14r, x14i, x15r, x15i
      x1r = a(2)
      x1i = a(3)
      x2r = a(4)
      x2i = a(5)
      x3r = a(6)
      x3i = a(7)
      x4r = a(8)
      x4i = a(9)
      x5r = a(10)
      x5i = a(11)
      x6r = a(12)
      x6i = a(13)
      x7r = a(14)
      x7i = a(15)
      x8r = a(16)
      x8i = a(17)
      x9r = a(18)
      x9i = a(19)
      x10r = a(20)
      x10i = a(21)
      x11r = a(22)
      x11i = a(23)
      x12r = a(24)
      x12i = a(25)
      x13r = a(26)
      x13i = a(27)
      x14r = a(28)
      x14i = a(29)
      x15r = a(30)
      x15i = a(31)
      a(2) = x15r
      a(3) = x15i
      a(4) = x7r
      a(5) = x7i
      a(6) = x11r
      a(7) = x11i
      a(8) = x3r
      a(9) = x3i
      a(10) = x13r
      a(11) = x13i
      a(12) = x5r
      a(13) = x5i
      a(14) = x9r
      a(15) = x9i
      a(16) = x1r
      a(17) = x1i
      a(18) = x14r
      a(19) = x14i
      a(20) = x6r
      a(21) = x6i
      a(22) = x10r
      a(23) = x10i
      a(24) = x2r
      a(25) = x2i
      a(26) = x12r
      a(27) = x12i
      a(28) = x4r
      a(29) = x4i
      a(30) = x8r
      a(31) = x8i
      end subroutine bitrv216neg
!
      subroutine bitrv208(a)
      real*8 a(0 : 15), x1r, x1i, x3r, x3i, x4r, x4i, x6r, x6i
      x1r = a(2)
      x1i = a(3)
      x3r = a(6)
      x3i = a(7)
      x4r = a(8)
      x4i = a(9)
      x6r = a(12)
      x6i = a(13)
      a(2) = x4r
      a(3) = x4i
      a(6) = x6r
      a(7) = x6i
      a(8) = x1r
      a(9) = x1i
      a(12) = x3r
      a(13) = x3i
      end subroutine bitrv208
!
      subroutine bitrv208neg(a)
      real*8 a(0 : 15), x1r, x1i, x2r, x2i, x3r, x3i, x4r, x4i
      real*8 x5r, x5i, x6r, x6i, x7r, x7i
      x1r = a(2)
      x1i = a(3)
      x2r = a(4)
      x2i = a(5)
      x3r = a(6)
      x3i = a(7)
      x4r = a(8)
      x4i = a(9)
      x5r = a(10)
      x5i = a(11)
      x6r = a(12)
      x6i = a(13)
      x7r = a(14)
      x7i = a(15)
      a(2) = x7r
      a(3) = x7i
      a(4) = x3r
      a(5) = x3i
      a(6) = x5r
      a(7) = x5i
      a(8) = x1r
      a(9) = x1i
      a(10) = x6r
      a(11) = x6i
      a(12) = x2r
      a(13) = x2i
      a(14) = x4r
      a(15) = x4i
      end subroutine bitrv208neg
!
      subroutine cftf1st(n, a, w)
      integer n, j, j0, j1, j2, j3, k, m, mh
      real*8 a(0 : n - 1), w(0 : *)
      real*8 wn4r, csc1, csc3, wk1r, wk1i, wk3r, wk3i
      real*8 wd1r, wd1i, wd3r, wd3i
      real*8 x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
      real*8 y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i
      mh = n / 8
      m = 2 * mh
      j1 = m
      j2 = j1 + m
      j3 = j2 + m
      x0r = a(0) + a(j2)
      x0i = a(1) + a(j2 + 1)
      x1r = a(0) - a(j2)
      x1i = a(1) - a(j2 + 1)
      x2r = a(j1) + a(j3)
      x2i = a(j1 + 1) + a(j3 + 1)
      x3r = a(j1) - a(j3)
      x3i = a(j1 + 1) - a(j3 + 1)
      a(0) = x0r + x2r
      a(1) = x0i + x2i
      a(j1) = x0r - x2r
      a(j1 + 1) = x0i - x2i
      a(j2) = x1r - x3i
      a(j2 + 1) = x1i + x3r
      a(j3) = x1r + x3i
      a(j3 + 1) = x1i - x3r
      wn4r = w(1)
      csc1 = w(2)
      csc3 = w(3)
      wd1r = 1
      wd1i = 0
      wd3r = 1
      wd3i = 0
      k = 0
      do j = 2, mh - 6, 4
          k = k + 4
          wk1r = csc1 * (wd1r + w(k))
          wk1i = csc1 * (wd1i + w(k + 1))
          wk3r = csc3 * (wd3r + w(k + 2))
          wk3i = csc3 * (wd3i + w(k + 3))
          wd1r = w(k)
          wd1i = w(k + 1)
          wd3r = w(k + 2)
          wd3i = w(k + 3)
          j1 = j + m
          j2 = j1 + m
          j3 = j2 + m
          x0r = a(j) + a(j2)
          x0i = a(j + 1) + a(j2 + 1)
          x1r = a(j) - a(j2)
          x1i = a(j + 1) - a(j2 + 1)
          y0r = a(j + 2) + a(j2 + 2)
          y0i = a(j + 3) + a(j2 + 3)
          y1r = a(j + 2) - a(j2 + 2)
          y1i = a(j + 3) - a(j2 + 3)
          x2r = a(j1) + a(j3)
          x2i = a(j1 + 1) + a(j3 + 1)
          x3r = a(j1) - a(j3)
          x3i = a(j1 + 1) - a(j3 + 1)
          y2r = a(j1 + 2) + a(j3 + 2)
          y2i = a(j1 + 3) + a(j3 + 3)
          y3r = a(j1 + 2) - a(j3 + 2)
          y3i = a(j1 + 3) - a(j3 + 3)
          a(j) = x0r + x2r
          a(j + 1) = x0i + x2i
          a(j + 2) = y0r + y2r
          a(j + 3) = y0i + y2i
          a(j1) = x0r - x2r
          a(j1 + 1) = x0i - x2i
          a(j1 + 2) = y0r - y2r
          a(j1 + 3) = y0i - y2i
          x0r = x1r - x3i
          x0i = x1i + x3r
          a(j2) = wk1r * x0r - wk1i * x0i
          a(j2 + 1) = wk1r * x0i + wk1i * x0r
          x0r = y1r - y3i
          x0i = y1i + y3r
          a(j2 + 2) = wd1r * x0r - wd1i * x0i
          a(j2 + 3) = wd1r * x0i + wd1i * x0r
          x0r = x1r + x3i
          x0i = x1i - x3r
          a(j3) = wk3r * x0r + wk3i * x0i
          a(j3 + 1) = wk3r * x0i - wk3i * x0r
          x0r = y1r + y3i
          x0i = y1i - y3r
          a(j3 + 2) = wd3r * x0r + wd3i * x0i
          a(j3 + 3) = wd3r * x0i - wd3i * x0r
          j0 = m - j
          j1 = j0 + m
          j2 = j1 + m
          j3 = j2 + m
          x0r = a(j0) + a(j2)
          x0i = a(j0 + 1) + a(j2 + 1)
          x1r = a(j0) - a(j2)
          x1i = a(j0 + 1) - a(j2 + 1)
          y0r = a(j0 - 2) + a(j2 - 2)
          y0i = a(j0 - 1) + a(j2 - 1)
          y1r = a(j0 - 2) - a(j2 - 2)
          y1i = a(j0 - 1) - a(j2 - 1)
          x2r = a(j1) + a(j3)
          x2i = a(j1 + 1) + a(j3 + 1)
          x3r = a(j1) - a(j3)
          x3i = a(j1 + 1) - a(j3 + 1)
          y2r = a(j1 - 2) + a(j3 - 2)
          y2i = a(j1 - 1) + a(j3 - 1)
          y3r = a(j1 - 2) - a(j3 - 2)
          y3i = a(j1 - 1) - a(j3 - 1)
          a(j0) = x0r + x2r
          a(j0 + 1) = x0i + x2i
          a(j0 - 2) = y0r + y2r
          a(j0 - 1) = y0i + y2i
          a(j1) = x0r - x2r
          a(j1 + 1) = x0i - x2i
          a(j1 - 2) = y0r - y2r
          a(j1 - 1) = y0i - y2i
          x0r = x1r - x3i
          x0i = x1i + x3r
          a(j2) = wk1i * x0r - wk1r * x0i
          a(j2 + 1) = wk1i * x0i + wk1r * x0r
          x0r = y1r - y3i
          x0i = y1i + y3r
          a(j2 - 2) = wd1i * x0r - wd1r * x0i
          a(j2 - 1) = wd1i * x0i + wd1r * x0r
          x0r = x1r + x3i
          x0i = x1i - x3r
          a(j3) = wk3i * x0r + wk3r * x0i
          a(j3 + 1) = wk3i * x0i - wk3r * x0r
          x0r = y1r + y3i
          x0i = y1i - y3r
          a(j3 - 2) = wd3i * x0r + wd3r * x0i
          a(j3 - 1) = wd3i * x0i - wd3r * x0r
      end do
      wk1r = csc1 * (wd1r + wn4r)
      wk1i = csc1 * (wd1i + wn4r)
      wk3r = csc3 * (wd3r - wn4r)
      wk3i = csc3 * (wd3i - wn4r)
      j0 = mh
      j1 = j0 + m
      j2 = j1 + m
      j3 = j2 + m
      x0r = a(j0 - 2) + a(j2 - 2)
      x0i = a(j0 - 1) + a(j2 - 1)
      x1r = a(j0 - 2) - a(j2 - 2)
      x1i = a(j0 - 1) - a(j2 - 1)
      x2r = a(j1 - 2) + a(j3 - 2)
      x2i = a(j1 - 1) + a(j3 - 1)
      x3r = a(j1 - 2) - a(j3 - 2)
      x3i = a(j1 - 1) - a(j3 - 1)
      a(j0 - 2) = x0r + x2r
      a(j0 - 1) = x0i + x2i
      a(j1 - 2) = x0r - x2r
      a(j1 - 1) = x0i - x2i
      x0r = x1r - x3i
      x0i = x1i + x3r
      a(j2 - 2) = wk1r * x0r - wk1i * x0i
      a(j2 - 1) = wk1r * x0i + wk1i * x0r
      x0r = x1r + x3i
      x0i = x1i - x3r
      a(j3 - 2) = wk3r * x0r + wk3i * x0i
      a(j3 - 1) = wk3r * x0i - wk3i * x0r
      x0r = a(j0) + a(j2)
      x0i = a(j0 + 1) + a(j2 + 1)
      x1r = a(j0) - a(j2)
      x1i = a(j0 + 1) - a(j2 + 1)
      x2r = a(j1) + a(j3)
      x2i = a(j1 + 1) + a(j3 + 1)
      x3r = a(j1) - a(j3)
      x3i = a(j1 + 1) - a(j3 + 1)
      a(j0) = x0r + x2r
      a(j0 + 1) = x0i + x2i
      a(j1) = x0r - x2r
      a(j1 + 1) = x0i - x2i
      x0r = x1r - x3i
      x0i = x1i + x3r
      a(j2) = wn4r * (x0r - x0i)
      a(j2 + 1) = wn4r * (x0i + x0r)
      x0r = x1r + x3i
      x0i = x1i - x3r
      a(j3) = -wn4r * (x0r + x0i)
      a(j3 + 1) = -wn4r * (x0i - x0r)
      x0r = a(j0 + 2) + a(j2 + 2)
      x0i = a(j0 + 3) + a(j2 + 3)
      x1r = a(j0 + 2) - a(j2 + 2)
      x1i = a(j0 + 3) - a(j2 + 3)
      x2r = a(j1 + 2) + a(j3 + 2)
      x2i = a(j1 + 3) + a(j3 + 3)
      x3r = a(j1 + 2) - a(j3 + 2)
      x3i = a(j1 + 3) - a(j3 + 3)
      a(j0 + 2) = x0r + x2r
      a(j0 + 3) = x0i + x2i
      a(j1 + 2) = x0r - x2r
      a(j1 + 3) = x0i - x2i
      x0r = x1r - x3i
      x0i = x1i + x3r
      a(j2 + 2) = wk1i * x0r - wk1r * x0i
      a(j2 + 3) = wk1i * x0i + wk1r * x0r
      x0r = x1r + x3i
      x0i = x1i - x3r
      a(j3 + 2) = wk3i * x0r + wk3r * x0i
      a(j3 + 3) = wk3i * x0i - wk3r * x0r
      end subroutine cftf1st
!
      subroutine cftb1st(n, a, w)
      integer n, j, j0, j1, j2, j3, k, m, mh
      real*8 a(0 : n - 1), w(0 : *)
      real*8 wn4r, csc1, csc3, wk1r, wk1i, wk3r, wk3i
      real*8 wd1r, wd1i, wd3r, wd3i
      real*8 x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
      real*8 y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i
      mh = n / 8
      m = 2 * mh
      j1 = m
      j2 = j1 + m
      j3 = j2 + m
      x0r = a(0) + a(j2)
      x0i = -a(1) - a(j2 + 1)
      x1r = a(0) - a(j2)
      x1i = -a(1) + a(j2 + 1)
      x2r = a(j1) + a(j3)
      x2i = a(j1 + 1) + a(j3 + 1)
      x3r = a(j1) - a(j3)
      x3i = a(j1 + 1) - a(j3 + 1)
      a(0) = x0r + x2r
      a(1) = x0i - x2i
      a(j1) = x0r - x2r
      a(j1 + 1) = x0i + x2i
      a(j2) = x1r + x3i
      a(j2 + 1) = x1i + x3r
      a(j3) = x1r - x3i
      a(j3 + 1) = x1i - x3r
      wn4r = w(1)
      csc1 = w(2)
      csc3 = w(3)
      wd1r = 1
      wd1i = 0
      wd3r = 1
      wd3i = 0
      k = 0
      do j = 2, mh - 6, 4
          k = k + 4
          wk1r = csc1 * (wd1r + w(k))
          wk1i = csc1 * (wd1i + w(k + 1))
          wk3r = csc3 * (wd3r + w(k + 2))
          wk3i = csc3 * (wd3i + w(k + 3))
          wd1r = w(k)
          wd1i = w(k + 1)
          wd3r = w(k + 2)
          wd3i = w(k + 3)
          j1 = j + m
          j2 = j1 + m
          j3 = j2 + m
          x0r = a(j) + a(j2)
          x0i = -a(j + 1) - a(j2 + 1)
          x1r = a(j) - a(j2)
          x1i = -a(j + 1) + a(j2 + 1)
          y0r = a(j + 2) + a(j2 + 2)
          y0i = -a(j + 3) - a(j2 + 3)
          y1r = a(j + 2) - a(j2 + 2)
          y1i = -a(j + 3) + a(j2 + 3)
          x2r = a(j1) + a(j3)
          x2i = a(j1 + 1) + a(j3 + 1)
          x3r = a(j1) - a(j3)
          x3i = a(j1 + 1) - a(j3 + 1)
          y2r = a(j1 + 2) + a(j3 + 2)
          y2i = a(j1 + 3) + a(j3 + 3)
          y3r = a(j1 + 2) - a(j3 + 2)
          y3i = a(j1 + 3) - a(j3 + 3)
          a(j) = x0r + x2r
          a(j + 1) = x0i - x2i
          a(j + 2) = y0r + y2r
          a(j + 3) = y0i - y2i
          a(j1) = x0r - x2r
          a(j1 + 1) = x0i + x2i
          a(j1 + 2) = y0r - y2r
          a(j1 + 3) = y0i + y2i
          x0r = x1r + x3i
          x0i = x1i + x3r
          a(j2) = wk1r * x0r - wk1i * x0i
          a(j2 + 1) = wk1r * x0i + wk1i * x0r
          x0r = y1r + y3i
          x0i = y1i + y3r
          a(j2 + 2) = wd1r * x0r - wd1i * x0i
          a(j2 + 3) = wd1r * x0i + wd1i * x0r
          x0r = x1r - x3i
          x0i = x1i - x3r
          a(j3) = wk3r * x0r + wk3i * x0i
          a(j3 + 1) = wk3r * x0i - wk3i * x0r
          x0r = y1r - y3i
          x0i = y1i - y3r
          a(j3 + 2) = wd3r * x0r + wd3i * x0i
          a(j3 + 3) = wd3r * x0i - wd3i * x0r
          j0 = m - j
          j1 = j0 + m
          j2 = j1 + m
          j3 = j2 + m
          x0r = a(j0) + a(j2)
          x0i = -a(j0 + 1) - a(j2 + 1)
          x1r = a(j0) - a(j2)
          x1i = -a(j0 + 1) + a(j2 + 1)
          y0r = a(j0 - 2) + a(j2 - 2)
          y0i = -a(j0 - 1) - a(j2 - 1)
          y1r = a(j0 - 2) - a(j2 - 2)
          y1i = -a(j0 - 1) + a(j2 - 1)
          x2r = a(j1) + a(j3)
          x2i = a(j1 + 1) + a(j3 + 1)
          x3r = a(j1) - a(j3)
          x3i = a(j1 + 1) - a(j3 + 1)
          y2r = a(j1 - 2) + a(j3 - 2)
          y2i = a(j1 - 1) + a(j3 - 1)
          y3r = a(j1 - 2) - a(j3 - 2)
          y3i = a(j1 - 1) - a(j3 - 1)
          a(j0) = x0r + x2r
          a(j0 + 1) = x0i - x2i
          a(j0 - 2) = y0r + y2r
          a(j0 - 1) = y0i - y2i
          a(j1) = x0r - x2r
          a(j1 + 1) = x0i + x2i
          a(j1 - 2) = y0r - y2r
          a(j1 - 1) = y0i + y2i
          x0r = x1r + x3i
          x0i = x1i + x3r
          a(j2) = wk1i * x0r - wk1r * x0i
          a(j2 + 1) = wk1i * x0i + wk1r * x0r
          x0r = y1r + y3i
          x0i = y1i + y3r
          a(j2 - 2) = wd1i * x0r - wd1r * x0i
          a(j2 - 1) = wd1i * x0i + wd1r * x0r
          x0r = x1r - x3i
          x0i = x1i - x3r
          a(j3) = wk3i * x0r + wk3r * x0i
          a(j3 + 1) = wk3i * x0i - wk3r * x0r
          x0r = y1r - y3i
          x0i = y1i - y3r
          a(j3 - 2) = wd3i * x0r + wd3r * x0i
          a(j3 - 1) = wd3i * x0i - wd3r * x0r
      end do
      wk1r = csc1 * (wd1r + wn4r)
      wk1i = csc1 * (wd1i + wn4r)
      wk3r = csc3 * (wd3r - wn4r)
      wk3i = csc3 * (wd3i - wn4r)
      j0 = mh
      j1 = j0 + m
      j2 = j1 + m
      j3 = j2 + m
      x0r = a(j0 - 2) + a(j2 - 2)
      x0i = -a(j0 - 1) - a(j2 - 1)
      x1r = a(j0 - 2) - a(j2 - 2)
      x1i = -a(j0 - 1) + a(j2 - 1)
      x2r = a(j1 - 2) + a(j3 - 2)
      x2i = a(j1 - 1) + a(j3 - 1)
      x3r = a(j1 - 2) - a(j3 - 2)
      x3i = a(j1 - 1) - a(j3 - 1)
      a(j0 - 2) = x0r + x2r
      a(j0 - 1) = x0i - x2i
      a(j1 - 2) = x0r - x2r
      a(j1 - 1) = x0i + x2i
      x0r = x1r + x3i
      x0i = x1i + x3r
      a(j2 - 2) = wk1r * x0r - wk1i * x0i
      a(j2 - 1) = wk1r * x0i + wk1i * x0r
      x0r = x1r - x3i
      x0i = x1i - x3r
      a(j3 - 2) = wk3r * x0r + wk3i * x0i
      a(j3 - 1) = wk3r * x0i - wk3i * x0r
      x0r = a(j0) + a(j2)
      x0i = -a(j0 + 1) - a(j2 + 1)
      x1r = a(j0) - a(j2)
      x1i = -a(j0 + 1) + a(j2 + 1)
      x2r = a(j1) + a(j3)
      x2i = a(j1 + 1) + a(j3 + 1)
      x3r = a(j1) - a(j3)
      x3i = a(j1 + 1) - a(j3 + 1)
      a(j0) = x0r + x2r
      a(j0 + 1) = x0i - x2i
      a(j1) = x0r - x2r
      a(j1 + 1) = x0i + x2i
      x0r = x1r + x3i
      x0i = x1i + x3r
      a(j2) = wn4r * (x0r - x0i)
      a(j2 + 1) = wn4r * (x0i + x0r)
      x0r = x1r - x3i
      x0i = x1i - x3r
      a(j3) = -wn4r * (x0r + x0i)
      a(j3 + 1) = -wn4r * (x0i - x0r)
      x0r = a(j0 + 2) + a(j2 + 2)
      x0i = -a(j0 + 3) - a(j2 + 3)
      x1r = a(j0 + 2) - a(j2 + 2)
      x1i = -a(j0 + 3) + a(j2 + 3)
      x2r = a(j1 + 2) + a(j3 + 2)
      x2i = a(j1 + 3) + a(j3 + 3)
      x3r = a(j1 + 2) - a(j3 + 2)
      x3i = a(j1 + 3) - a(j3 + 3)
      a(j0 + 2) = x0r + x2r
      a(j0 + 3) = x0i - x2i
      a(j1 + 2) = x0r - x2r
      a(j1 + 3) = x0i + x2i
      x0r = x1r + x3i
      x0i = x1i + x3r
      a(j2 + 2) = wk1i * x0r - wk1r * x0i
      a(j2 + 3) = wk1i * x0i + wk1r * x0r
      x0r = x1r - x3i
      x0i = x1i - x3r
      a(j3 + 2) = wk3i * x0r + wk3r * x0i
      a(j3 + 3) = wk3i * x0i - wk3r * x0r
      end subroutine cftb1st
!
      subroutine cftrec4(n, a, nw, w)
      integer n, nw, isplt, j, k, m
      real*8 a(0 : n - 1), w(0 : nw - 1)
      m = n
      do while (m .gt. 512)
          m = m / 4
          call cftmdl1(m, a(n - m), w(nw - m / 2))
      end do
      call cftleaf(m, 1, a(n - m), nw, w)
      k = 0
      do j = n - m, m, -m
          k = k + 1
          isplt = cfttree(m, j, k, a, nw, w)
          call cftleaf(m, isplt, a(j - m), nw, w)
      end do
      end subroutine cftrec4
!
      integer function cfttree(n, j, k, a, nw, w)
      integer n, j, k, nw, i, isplt, m
      real*8 a(0 : j - 1), w(0 : nw - 1)
      if (mod(k, 4) .ne. 0) then
          isplt = mod(k, 2)
          if (isplt .ne. 0) then
              call cftmdl1(n, a(j - n), w(nw - n / 2))
          else
              call cftmdl2(n, a(j - n), w(nw - n))
          end if
      else
          m = n
          i = k
          do while (mod(i, 4) .eq. 0)
              m = m * 4
              i = i / 4
          end do
          isplt = mod(i, 2)
          if (isplt .ne. 0) then
              do while (m .gt. 128)
                  call cftmdl1(m, a(j - m), w(nw - m / 2))
                  m = m / 4
              end do
          else
              do while (m .gt. 128)
                  call cftmdl2(m, a(j - m), w(nw - m))
                  m = m / 4
              end do
          end if
      end if
      cfttree = isplt
      end function cfttree
!
      subroutine cftleaf(n, isplt, a, nw, w)
      integer n, isplt, nw
      real*8 a(0 : n - 1), w(0 : nw - 1)
      if (n .eq. 512) then
          call cftmdl1(128, a, w(nw - 64))
          call cftf161(a, w(nw - 8))
          call cftf162(a(32), w(nw - 32))
          call cftf161(a(64), w(nw - 8))
          call cftf161(a(96), w(nw - 8))
          call cftmdl2(128, a(128), w(nw - 128))
          call cftf161(a(128), w(nw - 8))
          call cftf162(a(160), w(nw - 32))
          call cftf161(a(192), w(nw - 8))
          call cftf162(a(224), w(nw - 32))
          call cftmdl1(128, a(256), w(nw - 64))
          call cftf161(a(256), w(nw - 8))
          call cftf162(a(288), w(nw - 32))
          call cftf161(a(320), w(nw - 8))
          call cftf161(a(352), w(nw - 8))
          if (isplt .ne. 0) then
              call cftmdl1(128, a(384), w(nw - 64))
              call cftf161(a(480), w(nw - 8))
          else
              call cftmdl2(128, a(384), w(nw - 128))
              call cftf162(a(480), w(nw - 32))
          end if
          call cftf161(a(384), w(nw - 8))
          call cftf162(a(416), w(nw - 32))
          call cftf161(a(448), w(nw - 8))
      else
          call cftmdl1(64, a, w(nw - 32))
          call cftf081(a, w(nw - 8))
          call cftf082(a(16), w(nw - 8))
          call cftf081(a(32), w(nw - 8))
          call cftf081(a(48), w(nw - 8))
          call cftmdl2(64, a(64), w(nw - 64))
          call cftf081(a(64), w(nw - 8))
          call cftf082(a(80), w(nw - 8))
          call cftf081(a(96), w(nw - 8))
          call cftf082(a(112), w(nw - 8))
          call cftmdl1(64, a(128), w(nw - 32))
          call cftf081(a(128), w(nw - 8))
          call cftf082(a(144), w(nw - 8))
          call cftf081(a(160), w(nw - 8))
          call cftf081(a(176), w(nw - 8))
          if (isplt .ne. 0) then
              call cftmdl1(64, a(192), w(nw - 32))
              call cftf081(a(240), w(nw - 8))
          else
              call cftmdl2(64, a(192), w(nw - 64))
              call cftf082(a(240), w(nw - 8))
          end if
          call cftf081(a(192), w(nw - 8))
          call cftf082(a(208), w(nw - 8))
          call cftf081(a(224), w(nw - 8))
      end if
      end subroutine cftleaf
!
      subroutine cftmdl1(n, a, w)
      integer n, j, j0, j1, j2, j3, k, m, mh
      real*8 a(0 : n - 1), w(0 : *)
      real*8 wn4r, wk1r, wk1i, wk3r, wk3i
      real*8 x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
      mh = n / 8
      m = 2 * mh
      j1 = m
      j2 = j1 + m
      j3 = j2 + m
      x0r = a(0) + a(j2)
      x0i = a(1) + a(j2 + 1)
      x1r = a(0) - a(j2)
      x1i = a(1) - a(j2 + 1)
      x2r = a(j1) + a(j3)
      x2i = a(j1 + 1) + a(j3 + 1)
      x3r = a(j1) - a(j3)
      x3i = a(j1 + 1) - a(j3 + 1)
      a(0) = x0r + x2r
      a(1) = x0i + x2i
      a(j1) = x0r - x2r
      a(j1 + 1) = x0i - x2i
      a(j2) = x1r - x3i
      a(j2 + 1) = x1i + x3r
      a(j3) = x1r + x3i
      a(j3 + 1) = x1i - x3r
      wn4r = w(1)
      k = 0
      do j = 2, mh - 2, 2
          k = k + 4
          wk1r = w(k)
          wk1i = w(k + 1)
          wk3r = w(k + 2)
          wk3i = w(k + 3)
          j1 = j + m
          j2 = j1 + m
          j3 = j2 + m
          x0r = a(j) + a(j2)
          x0i = a(j + 1) + a(j2 + 1)
          x1r = a(j) - a(j2)
          x1i = a(j + 1) - a(j2 + 1)
          x2r = a(j1) + a(j3)
          x2i = a(j1 + 1) + a(j3 + 1)
          x3r = a(j1) - a(j3)
          x3i = a(j1 + 1) - a(j3 + 1)
          a(j) = x0r + x2r
          a(j + 1) = x0i + x2i
          a(j1) = x0r - x2r
          a(j1 + 1) = x0i - x2i
          x0r = x1r - x3i
          x0i = x1i + x3r
          a(j2) = wk1r * x0r - wk1i * x0i
          a(j2 + 1) = wk1r * x0i + wk1i * x0r
          x0r = x1r + x3i
          x0i = x1i - x3r
          a(j3) = wk3r * x0r + wk3i * x0i
          a(j3 + 1) = wk3r * x0i - wk3i * x0r
          j0 = m - j
          j1 = j0 + m
          j2 = j1 + m
          j3 = j2 + m
          x0r = a(j0) + a(j2)
          x0i = a(j0 + 1) + a(j2 + 1)
          x1r = a(j0) - a(j2)
          x1i = a(j0 + 1) - a(j2 + 1)
          x2r = a(j1) + a(j3)
          x2i = a(j1 + 1) + a(j3 + 1)
          x3r = a(j1) - a(j3)
          x3i = a(j1 + 1) - a(j3 + 1)
          a(j0) = x0r + x2r
          a(j0 + 1) = x0i + x2i
          a(j1) = x0r - x2r
          a(j1 + 1) = x0i - x2i
          x0r = x1r - x3i
          x0i = x1i + x3r
          a(j2) = wk1i * x0r - wk1r * x0i
          a(j2 + 1) = wk1i * x0i + wk1r * x0r
          x0r = x1r + x3i
          x0i = x1i - x3r
          a(j3) = wk3i * x0r + wk3r * x0i
          a(j3 + 1) = wk3i * x0i - wk3r * x0r
      end do
      j0 = mh
      j1 = j0 + m
      j2 = j1 + m
      j3 = j2 + m
      x0r = a(j0) + a(j2)
      x0i = a(j0 + 1) + a(j2 + 1)
      x1r = a(j0) - a(j2)
      x1i = a(j0 + 1) - a(j2 + 1)
      x2r = a(j1) + a(j3)
      x2i = a(j1 + 1) + a(j3 + 1)
      x3r = a(j1) - a(j3)
      x3i = a(j1 + 1) - a(j3 + 1)
      a(j0) = x0r + x2r
      a(j0 + 1) = x0i + x2i
      a(j1) = x0r - x2r
      a(j1 + 1) = x0i - x2i
      x0r = x1r - x3i
      x0i = x1i + x3r
      a(j2) = wn4r * (x0r - x0i)
      a(j2 + 1) = wn4r * (x0i + x0r)
      x0r = x1r + x3i
      x0i = x1i - x3r
      a(j3) = -wn4r * (x0r + x0i)
      a(j3 + 1) = -wn4r * (x0i - x0r)
      end subroutine cftmdl1
!
      subroutine cftmdl2(n, a, w)
      integer n, j, j0, j1, j2, j3, k, kr, m, mh
      real*8 a(0 : n - 1), w(0 : *)
      real*8 wn4r, wk1r, wk1i, wk3r, wk3i, wd1r, wd1i, wd3r, wd3i
      real*8 x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
      real*8 y0r, y0i, y2r, y2i
      mh = n / 8
      m = 2 * mh
      wn4r = w(1)
      j1 = m
      j2 = j1 + m
      j3 = j2 + m
      x0r = a(0) - a(j2 + 1)
      x0i = a(1) + a(j2)
      x1r = a(0) + a(j2 + 1)
      x1i = a(1) - a(j2)
      x2r = a(j1) - a(j3 + 1)
      x2i = a(j1 + 1) + a(j3)
      x3r = a(j1) + a(j3 + 1)
      x3i = a(j1 + 1) - a(j3)
      y0r = wn4r * (x2r - x2i)
      y0i = wn4r * (x2i + x2r)
      a(0) = x0r + y0r
      a(1) = x0i + y0i
      a(j1) = x0r - y0r
      a(j1 + 1) = x0i - y0i
      y0r = wn4r * (x3r - x3i)
      y0i = wn4r * (x3i + x3r)
      a(j2) = x1r - y0i
      a(j2 + 1) = x1i + y0r
      a(j3) = x1r + y0i
      a(j3 + 1) = x1i - y0r
      k = 0
      kr = 2 * m
      do j = 2, mh - 2, 2
          k = k + 4
          wk1r = w(k)
          wk1i = w(k + 1)
          wk3r = w(k + 2)
          wk3i = w(k + 3)
          kr = kr - 4
          wd1i = w(kr)
          wd1r = w(kr + 1)
          wd3i = w(kr + 2)
          wd3r = w(kr + 3)
          j1 = j + m
          j2 = j1 + m
          j3 = j2 + m
          x0r = a(j) - a(j2 + 1)
          x0i = a(j + 1) + a(j2)
          x1r = a(j) + a(j2 + 1)
          x1i = a(j + 1) - a(j2)
          x2r = a(j1) - a(j3 + 1)
          x2i = a(j1 + 1) + a(j3)
          x3r = a(j1) + a(j3 + 1)
          x3i = a(j1 + 1) - a(j3)
          y0r = wk1r * x0r - wk1i * x0i
          y0i = wk1r * x0i + wk1i * x0r
          y2r = wd1r * x2r - wd1i * x2i
          y2i = wd1r * x2i + wd1i * x2r
          a(j) = y0r + y2r
          a(j + 1) = y0i + y2i
          a(j1) = y0r - y2r
          a(j1 + 1) = y0i - y2i
          y0r = wk3r * x1r + wk3i * x1i
          y0i = wk3r * x1i - wk3i * x1r
          y2r = wd3r * x3r + wd3i * x3i
          y2i = wd3r * x3i - wd3i * x3r
          a(j2) = y0r + y2r
          a(j2 + 1) = y0i + y2i
          a(j3) = y0r - y2r
          a(j3 + 1) = y0i - y2i
          j0 = m - j
          j1 = j0 + m
          j2 = j1 + m
          j3 = j2 + m
          x0r = a(j0) - a(j2 + 1)
          x0i = a(j0 + 1) + a(j2)
          x1r = a(j0) + a(j2 + 1)
          x1i = a(j0 + 1) - a(j2)
          x2r = a(j1) - a(j3 + 1)
          x2i = a(j1 + 1) + a(j3)
          x3r = a(j1) + a(j3 + 1)
          x3i = a(j1 + 1) - a(j3)
          y0r = wd1i * x0r - wd1r * x0i
          y0i = wd1i * x0i + wd1r * x0r
          y2r = wk1i * x2r - wk1r * x2i
          y2i = wk1i * x2i + wk1r * x2r
          a(j0) = y0r + y2r
          a(j0 + 1) = y0i + y2i
          a(j1) = y0r - y2r
          a(j1 + 1) = y0i - y2i
          y0r = wd3i * x1r + wd3r * x1i
          y0i = wd3i * x1i - wd3r * x1r
          y2r = wk3i * x3r + wk3r * x3i
          y2i = wk3i * x3i - wk3r * x3r
          a(j2) = y0r + y2r
          a(j2 + 1) = y0i + y2i
          a(j3) = y0r - y2r
          a(j3 + 1) = y0i - y2i
      end do
      wk1r = w(m)
      wk1i = w(m + 1)
      j0 = mh
      j1 = j0 + m
      j2 = j1 + m
      j3 = j2 + m
      x0r = a(j0) - a(j2 + 1)
      x0i = a(j0 + 1) + a(j2)
      x1r = a(j0) + a(j2 + 1)
      x1i = a(j0 + 1) - a(j2)
      x2r = a(j1) - a(j3 + 1)
      x2i = a(j1 + 1) + a(j3)
      x3r = a(j1) + a(j3 + 1)
      x3i = a(j1 + 1) - a(j3)
      y0r = wk1r * x0r - wk1i * x0i
      y0i = wk1r * x0i + wk1i * x0r
      y2r = wk1i * x2r - wk1r * x2i
      y2i = wk1i * x2i + wk1r * x2r
      a(j0) = y0r + y2r
      a(j0 + 1) = y0i + y2i
      a(j1) = y0r - y2r
      a(j1 + 1) = y0i - y2i
      y0r = wk1i * x1r - wk1r * x1i
      y0i = wk1i * x1i + wk1r * x1r
      y2r = wk1r * x3r - wk1i * x3i
      y2i = wk1r * x3i + wk1i * x3r
      a(j2) = y0r - y2r
      a(j2 + 1) = y0i - y2i
      a(j3) = y0r + y2r
      a(j3 + 1) = y0i + y2i
      end subroutine cftmdl2
!
      subroutine cftfx41(n, a, nw, w)
      integer n, nw
      real*8 a(0 : n - 1), w(0 : nw - 1)
      if (n .eq. 128) then
          call cftf161(a, w(nw - 8))
          call cftf162(a(32), w(nw - 32))
          call cftf161(a(64), w(nw - 8))
          call cftf161(a(96), w(nw - 8))
      else
          call cftf081(a, w(nw - 8))
          call cftf082(a(16), w(nw - 8))
          call cftf081(a(32), w(nw - 8))
          call cftf081(a(48), w(nw - 8))
      end if
      end subroutine cftfx41
!
      subroutine cftf161(a, w)
      real*8 a(0 : 31), w(0 : *), wn4r, wk1r, wk1i
      real*8 x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
      real*8 y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i
      real*8 y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i
      real*8 y8r, y8i, y9r, y9i, y10r, y10i, y11r, y11i
      real*8 y12r, y12i, y13r, y13i, y14r, y14i, y15r, y15i
      wn4r = w(1)
      wk1r = w(2)
      wk1i = w(3)
      x0r = a(0) + a(16)
      x0i = a(1) + a(17)
      x1r = a(0) - a(16)
      x1i = a(1) - a(17)
      x2r = a(8) + a(24)
      x2i = a(9) + a(25)
      x3r = a(8) - a(24)
      x3i = a(9) - a(25)
      y0r = x0r + x2r
      y0i = x0i + x2i
      y4r = x0r - x2r
      y4i = x0i - x2i
      y8r = x1r - x3i
      y8i = x1i + x3r
      y12r = x1r + x3i
      y12i = x1i - x3r
      x0r = a(2) + a(18)
      x0i = a(3) + a(19)
      x1r = a(2) - a(18)
      x1i = a(3) - a(19)
      x2r = a(10) + a(26)
      x2i = a(11) + a(27)
      x3r = a(10) - a(26)
      x3i = a(11) - a(27)
      y1r = x0r + x2r
      y1i = x0i + x2i
      y5r = x0r - x2r
      y5i = x0i - x2i
      x0r = x1r - x3i
      x0i = x1i + x3r
      y9r = wk1r * x0r - wk1i * x0i
      y9i = wk1r * x0i + wk1i * x0r
      x0r = x1r + x3i
      x0i = x1i - x3r
      y13r = wk1i * x0r - wk1r * x0i
      y13i = wk1i * x0i + wk1r * x0r
      x0r = a(4) + a(20)
      x0i = a(5) + a(21)
      x1r = a(4) - a(20)
      x1i = a(5) - a(21)
      x2r = a(12) + a(28)
      x2i = a(13) + a(29)
      x3r = a(12) - a(28)
      x3i = a(13) - a(29)
      y2r = x0r + x2r
      y2i = x0i + x2i
      y6r = x0r - x2r
      y6i = x0i - x2i
      x0r = x1r - x3i
      x0i = x1i + x3r
      y10r = wn4r * (x0r - x0i)
      y10i = wn4r * (x0i + x0r)
      x0r = x1r + x3i
      x0i = x1i - x3r
      y14r = wn4r * (x0r + x0i)
      y14i = wn4r * (x0i - x0r)
      x0r = a(6) + a(22)
      x0i = a(7) + a(23)
      x1r = a(6) - a(22)
      x1i = a(7) - a(23)
      x2r = a(14) + a(30)
      x2i = a(15) + a(31)
      x3r = a(14) - a(30)
      x3i = a(15) - a(31)
      y3r = x0r + x2r
      y3i = x0i + x2i
      y7r = x0r - x2r
      y7i = x0i - x2i
      x0r = x1r - x3i
      x0i = x1i + x3r
      y11r = wk1i * x0r - wk1r * x0i
      y11i = wk1i * x0i + wk1r * x0r
      x0r = x1r + x3i
      x0i = x1i - x3r
      y15r = wk1r * x0r - wk1i * x0i
      y15i = wk1r * x0i + wk1i * x0r
      x0r = y12r - y14r
      x0i = y12i - y14i
      x1r = y12r + y14r
      x1i = y12i + y14i
      x2r = y13r - y15r
      x2i = y13i - y15i
      x3r = y13r + y15r
      x3i = y13i + y15i
      a(24) = x0r + x2r
      a(25) = x0i + x2i
      a(26) = x0r - x2r
      a(27) = x0i - x2i
      a(28) = x1r - x3i
      a(29) = x1i + x3r
      a(30) = x1r + x3i
      a(31) = x1i - x3r
      x0r = y8r + y10r
      x0i = y8i + y10i
      x1r = y8r - y10r
      x1i = y8i - y10i
      x2r = y9r + y11r
      x2i = y9i + y11i
      x3r = y9r - y11r
      x3i = y9i - y11i
      a(16) = x0r + x2r
      a(17) = x0i + x2i
      a(18) = x0r - x2r
      a(19) = x0i - x2i
      a(20) = x1r - x3i
      a(21) = x1i + x3r
      a(22) = x1r + x3i
      a(23) = x1i - x3r
      x0r = y5r - y7i
      x0i = y5i + y7r
      x2r = wn4r * (x0r - x0i)
      x2i = wn4r * (x0i + x0r)
      x0r = y5r + y7i
      x0i = y5i - y7r
      x3r = wn4r * (x0r - x0i)
      x3i = wn4r * (x0i + x0r)
      x0r = y4r - y6i
      x0i = y4i + y6r
      x1r = y4r + y6i
      x1i = y4i - y6r
      a(8) = x0r + x2r
      a(9) = x0i + x2i
      a(10) = x0r - x2r
      a(11) = x0i - x2i
      a(12) = x1r - x3i
      a(13) = x1i + x3r
      a(14) = x1r + x3i
      a(15) = x1i - x3r
      x0r = y0r + y2r
      x0i = y0i + y2i
      x1r = y0r - y2r
      x1i = y0i - y2i
      x2r = y1r + y3r
      x2i = y1i + y3i
      x3r = y1r - y3r
      x3i = y1i - y3i
      a(0) = x0r + x2r
      a(1) = x0i + x2i
      a(2) = x0r - x2r
      a(3) = x0i - x2i
      a(4) = x1r - x3i
      a(5) = x1i + x3r
      a(6) = x1r + x3i
      a(7) = x1i - x3r
      end subroutine cftf161
!
      subroutine cftf162(a, w)
      real*8 a(0 : 31), w(0 : *)
      real*8 wn4r, wk1r, wk1i, wk2r, wk2i, wk3r, wk3i
      real*8 x0r, x0i, x1r, x1i, x2r, x2i
      real*8 y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i
      real*8 y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i
      real*8 y8r, y8i, y9r, y9i, y10r, y10i, y11r, y11i
      real*8 y12r, y12i, y13r, y13i, y14r, y14i, y15r, y15i
      wn4r = w(1)
      wk1r = w(4)
      wk1i = w(5)
      wk3r = w(6)
      wk3i = -w(7)
      wk2r = w(8)
      wk2i = w(9)
      x1r = a(0) - a(17)
      x1i = a(1) + a(16)
      x0r = a(8) - a(25)
      x0i = a(9) + a(24)
      x2r = wn4r * (x0r - x0i)
      x2i = wn4r * (x0i + x0r)
      y0r = x1r + x2r
      y0i = x1i + x2i
      y4r = x1r - x2r
      y4i = x1i - x2i
      x1r = a(0) + a(17)
      x1i = a(1) - a(16)
      x0r = a(8) + a(25)
      x0i = a(9) - a(24)
      x2r = wn4r * (x0r - x0i)
      x2i = wn4r * (x0i + x0r)
      y8r = x1r - x2i
      y8i = x1i + x2r
      y12r = x1r + x2i
      y12i = x1i - x2r
      x0r = a(2) - a(19)
      x0i = a(3) + a(18)
      x1r = wk1r * x0r - wk1i * x0i
      x1i = wk1r * x0i + wk1i * x0r
      x0r = a(10) - a(27)
      x0i = a(11) + a(26)
      x2r = wk3i * x0r - wk3r * x0i
      x2i = wk3i * x0i + wk3r * x0r
      y1r = x1r + x2r
      y1i = x1i + x2i
      y5r = x1r - x2r
      y5i = x1i - x2i
      x0r = a(2) + a(19)
      x0i = a(3) - a(18)
      x1r = wk3r * x0r - wk3i * x0i
      x1i = wk3r * x0i + wk3i * x0r
      x0r = a(10) + a(27)
      x0i = a(11) - a(26)
      x2r = wk1r * x0r + wk1i * x0i
      x2i = wk1r * x0i - wk1i * x0r
      y9r = x1r - x2r
      y9i = x1i - x2i
      y13r = x1r + x2r
      y13i = x1i + x2i
      x0r = a(4) - a(21)
      x0i = a(5) + a(20)
      x1r = wk2r * x0r - wk2i * x0i
      x1i = wk2r * x0i + wk2i * x0r
      x0r = a(12) - a(29)
      x0i = a(13) + a(28)
      x2r = wk2i * x0r - wk2r * x0i
      x2i = wk2i * x0i + wk2r * x0r
      y2r = x1r + x2r
      y2i = x1i + x2i
      y6r = x1r - x2r
      y6i = x1i - x2i
      x0r = a(4) + a(21)
      x0i = a(5) - a(20)
      x1r = wk2i * x0r - wk2r * x0i
      x1i = wk2i * x0i + wk2r * x0r
      x0r = a(12) + a(29)
      x0i = a(13) - a(28)
      x2r = wk2r * x0r - wk2i * x0i
      x2i = wk2r * x0i + wk2i * x0r
      y10r = x1r - x2r
      y10i = x1i - x2i
      y14r = x1r + x2r
      y14i = x1i + x2i
      x0r = a(6) - a(23)
      x0i = a(7) + a(22)
      x1r = wk3r * x0r - wk3i * x0i
      x1i = wk3r * x0i + wk3i * x0r
      x0r = a(14) - a(31)
      x0i = a(15) + a(30)
      x2r = wk1i * x0r - wk1r * x0i
      x2i = wk1i * x0i + wk1r * x0r
      y3r = x1r + x2r
      y3i = x1i + x2i
      y7r = x1r - x2r
      y7i = x1i - x2i
      x0r = a(6) + a(23)
      x0i = a(7) - a(22)
      x1r = wk1i * x0r + wk1r * x0i
      x1i = wk1i * x0i - wk1r * x0r
      x0r = a(14) + a(31)
      x0i = a(15) - a(30)
      x2r = wk3i * x0r - wk3r * x0i
      x2i = wk3i * x0i + wk3r * x0r
      y11r = x1r + x2r
      y11i = x1i + x2i
      y15r = x1r - x2r
      y15i = x1i - x2i
      x1r = y0r + y2r
      x1i = y0i + y2i
      x2r = y1r + y3r
      x2i = y1i + y3i
      a(0) = x1r + x2r
      a(1) = x1i + x2i
      a(2) = x1r - x2r
      a(3) = x1i - x2i
      x1r = y0r - y2r
      x1i = y0i - y2i
      x2r = y1r - y3r
      x2i = y1i - y3i
      a(4) = x1r - x2i
      a(5) = x1i + x2r
      a(6) = x1r + x2i
      a(7) = x1i - x2r
      x1r = y4r - y6i
      x1i = y4i + y6r
      x0r = y5r - y7i
      x0i = y5i + y7r
      x2r = wn4r * (x0r - x0i)
      x2i = wn4r * (x0i + x0r)
      a(8) = x1r + x2r
      a(9) = x1i + x2i
      a(10) = x1r - x2r
      a(11) = x1i - x2i
      x1r = y4r + y6i
      x1i = y4i - y6r
      x0r = y5r + y7i
      x0i = y5i - y7r
      x2r = wn4r * (x0r - x0i)
      x2i = wn4r * (x0i + x0r)
      a(12) = x1r - x2i
      a(13) = x1i + x2r
      a(14) = x1r + x2i
      a(15) = x1i - x2r
      x1r = y8r + y10r
      x1i = y8i + y10i
      x2r = y9r - y11r
      x2i = y9i - y11i
      a(16) = x1r + x2r
      a(17) = x1i + x2i
      a(18) = x1r - x2r
      a(19) = x1i - x2i
      x1r = y8r - y10r
      x1i = y8i - y10i
      x2r = y9r + y11r
      x2i = y9i + y11i
      a(20) = x1r - x2i
      a(21) = x1i + x2r
      a(22) = x1r + x2i
      a(23) = x1i - x2r
      x1r = y12r - y14i
      x1i = y12i + y14r
      x0r = y13r + y15i
      x0i = y13i - y15r
      x2r = wn4r * (x0r - x0i)
      x2i = wn4r * (x0i + x0r)
      a(24) = x1r + x2r
      a(25) = x1i + x2i
      a(26) = x1r - x2r
      a(27) = x1i - x2i
      x1r = y12r + y14i
      x1i = y12i - y14r
      x0r = y13r - y15i
      x0i = y13i + y15r
      x2r = wn4r * (x0r - x0i)
      x2i = wn4r * (x0i + x0r)
      a(28) = x1r - x2i
      a(29) = x1i + x2r
      a(30) = x1r + x2i
      a(31) = x1i - x2r
      end subroutine cftf162
!
      subroutine cftf081(a, w)
      real*8 a(0 : 15), w(0 : *)
      real*8 wn4r, x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
      real*8 y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i
      real*8 y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i
      wn4r = w(1)
      x0r = a(0) + a(8)
      x0i = a(1) + a(9)
      x1r = a(0) - a(8)
      x1i = a(1) - a(9)
      x2r = a(4) + a(12)
      x2i = a(5) + a(13)
      x3r = a(4) - a(12)
      x3i = a(5) - a(13)
      y0r = x0r + x2r
      y0i = x0i + x2i
      y2r = x0r - x2r
      y2i = x0i - x2i
      y1r = x1r - x3i
      y1i = x1i + x3r
      y3r = x1r + x3i
      y3i = x1i - x3r
      x0r = a(2) + a(10)
      x0i = a(3) + a(11)
      x1r = a(2) - a(10)
      x1i = a(3) - a(11)
      x2r = a(6) + a(14)
      x2i = a(7) + a(15)
      x3r = a(6) - a(14)
      x3i = a(7) - a(15)
      y4r = x0r + x2r
      y4i = x0i + x2i
      y6r = x0r - x2r
      y6i = x0i - x2i
      x0r = x1r - x3i
      x0i = x1i + x3r
      x2r = x1r + x3i
      x2i = x1i - x3r
      y5r = wn4r * (x0r - x0i)
      y5i = wn4r * (x0r + x0i)
      y7r = wn4r * (x2r - x2i)
      y7i = wn4r * (x2r + x2i)
      a(8) = y1r + y5r
      a(9) = y1i + y5i
      a(10) = y1r - y5r
      a(11) = y1i - y5i
      a(12) = y3r - y7i
      a(13) = y3i + y7r
      a(14) = y3r + y7i
      a(15) = y3i - y7r
      a(0) = y0r + y4r
      a(1) = y0i + y4i
      a(2) = y0r - y4r
      a(3) = y0i - y4i
      a(4) = y2r - y6i
      a(5) = y2i + y6r
      a(6) = y2r + y6i
      a(7) = y2i - y6r
      end subroutine cftf081
!
      subroutine cftf082(a, w)
      real*8 a(0 : 15), w(0 : *)
      real*8 wn4r, wk1r, wk1i, x0r, x0i, x1r, x1i
      real*8 y0r, y0i, y1r, y1i, y2r, y2i, y3r, y3i
      real*8 y4r, y4i, y5r, y5i, y6r, y6i, y7r, y7i
      wn4r = w(1)
      wk1r = w(2)
      wk1i = w(3)
      y0r = a(0) - a(9)
      y0i = a(1) + a(8)
      y1r = a(0) + a(9)
      y1i = a(1) - a(8)
      x0r = a(4) - a(13)
      x0i = a(5) + a(12)
      y2r = wn4r * (x0r - x0i)
      y2i = wn4r * (x0i + x0r)
      x0r = a(4) + a(13)
      x0i = a(5) - a(12)
      y3r = wn4r * (x0r - x0i)
      y3i = wn4r * (x0i + x0r)
      x0r = a(2) - a(11)
      x0i = a(3) + a(10)
      y4r = wk1r * x0r - wk1i * x0i
      y4i = wk1r * x0i + wk1i * x0r
      x0r = a(2) + a(11)
      x0i = a(3) - a(10)
      y5r = wk1i * x0r - wk1r * x0i
      y5i = wk1i * x0i + wk1r * x0r
      x0r = a(6) - a(15)
      x0i = a(7) + a(14)
      y6r = wk1i * x0r - wk1r * x0i
      y6i = wk1i * x0i + wk1r * x0r
      x0r = a(6) + a(15)
      x0i = a(7) - a(14)
      y7r = wk1r * x0r - wk1i * x0i
      y7i = wk1r * x0i + wk1i * x0r
      x0r = y0r + y2r
      x0i = y0i + y2i
      x1r = y4r + y6r
      x1i = y4i + y6i
      a(0) = x0r + x1r
      a(1) = x0i + x1i
      a(2) = x0r - x1r
      a(3) = x0i - x1i
      x0r = y0r - y2r
      x0i = y0i - y2i
      x1r = y4r - y6r
      x1i = y4i - y6i
      a(4) = x0r - x1i
      a(5) = x0i + x1r
      a(6) = x0r + x1i
      a(7) = x0i - x1r
      x0r = y1r - y3i
      x0i = y1i + y3r
      x1r = y5r - y7r
      x1i = y5i - y7i
      a(8) = x0r + x1r
      a(9) = x0i + x1i
      a(10) = x0r - x1r
      a(11) = x0i - x1i
      x0r = y1r + y3i
      x0i = y1i - y3r
      x1r = y5r + y7r
      x1i = y5i + y7i
      a(12) = x0r - x1i
      a(13) = x0i + x1r
      a(14) = x0r + x1i
      a(15) = x0i - x1r
      end subroutine cftf082
!
      subroutine cftf040(a)
      real*8 a(0 : 7), x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
      x0r = a(0) + a(4)
      x0i = a(1) + a(5)
      x1r = a(0) - a(4)
      x1i = a(1) - a(5)
      x2r = a(2) + a(6)
      x2i = a(3) + a(7)
      x3r = a(2) - a(6)
      x3i = a(3) - a(7)
      a(0) = x0r + x2r
      a(1) = x0i + x2i
      a(2) = x1r - x3i
      a(3) = x1i + x3r
      a(4) = x0r - x2r
      a(5) = x0i - x2i
      a(6) = x1r + x3i
      a(7) = x1i - x3r
      end subroutine cftf040
!
      subroutine cftb040(a)
      real*8 a(0 : 7), x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
      x0r = a(0) + a(4)
      x0i = a(1) + a(5)
      x1r = a(0) - a(4)
      x1i = a(1) - a(5)
      x2r = a(2) + a(6)
      x2i = a(3) + a(7)
      x3r = a(2) - a(6)
      x3i = a(3) - a(7)
      a(0) = x0r + x2r
      a(1) = x0i + x2i
      a(2) = x1r + x3i
      a(3) = x1i - x3r
      a(4) = x0r - x2r
      a(5) = x0i - x2i
      a(6) = x1r - x3i
      a(7) = x1i + x3r
      end subroutine cftb040
!
      subroutine cftx020(a)
      real*8 a(0 : 3), x0r, x0i
      x0r = a(0) - a(2)
      x0i = a(1) - a(3)
      a(0) = a(0) + a(2)
      a(1) = a(1) + a(3)
      a(2) = x0r
      a(3) = x0i
      end subroutine cftx020
!
      subroutine rftfsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks, m
      real*8 a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr, xi, yr, yi
      m = n / 2
      ks = 2 * nc / m
      kk = 0
      do j = 2, m - 2, 2
          k = n - j
          kk = kk + ks
          wkr = 0.5d0 - c(nc - kk)
          wki = c(kk)
          xr = a(j) - a(k)
          xi = a(j + 1) + a(k + 1)
          yr = wkr * xr - wki * xi
          yi = wkr * xi + wki * xr
          a(j) = a(j) - yr
          a(j + 1) = a(j + 1) - yi
          a(k) = a(k) + yr
          a(k + 1) = a(k + 1) - yi
      end do
      end subroutine rftfsub
!
      subroutine rftbsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks, m
      real*8 a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr, xi, yr, yi
      m = n / 2
      ks = 2 * nc / m
      kk = 0
      do j = 2, m - 2, 2
          k = n - j
          kk = kk + ks
          wkr = 0.5d0 - c(nc - kk)
          wki = c(kk)
          xr = a(j) - a(k)
          xi = a(j + 1) + a(k + 1)
          yr = wkr * xr + wki * xi
          yi = wkr * xi - wki * xr
          a(j) = a(j) - yr
          a(j + 1) = a(j + 1) - yi
          a(k) = a(k) + yr
          a(k + 1) = a(k + 1) - yi
      end do
      end subroutine rftbsub
!
      subroutine dctsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks, m
      real*8 a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr
      m = n / 2
      ks = nc / n
      kk = 0
      do j = 1, m - 1
          k = n - j
          kk = kk + ks
          wkr = c(kk) - c(nc - kk)
          wki = c(kk) + c(nc - kk)
          xr = wki * a(j) - wkr * a(k)
          a(j) = wkr * a(j) + wki * a(k)
          a(k) = xr
      end do
      a(m) = c(0) * a(m)
      end subroutine dctsub
!
      subroutine dstsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks, m
      real*8 a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr
      m = n / 2
      ks = nc / n
      kk = 0
      do j = 1, m - 1
          k = n - j
          kk = kk + ks
          wkr = c(kk) - c(nc - kk)
          wki = c(kk) + c(nc - kk)
          xr = wki * a(k) - wkr * a(j)
          a(k) = wkr * a(k) + wki * a(j)
          a(j) = xr
      end do
      a(m) = c(0) * a(m)
      end subroutine dstsub
!
end module fftsg 

