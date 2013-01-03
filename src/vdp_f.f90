!! Copyright (c) 2012 Michel Crucifix <michel.crucifix@uclouvain.be>
!!
!! Permission is hereby granted, free of charge, to any person obtaining
!! a copy of this software and associated documentation files (the
!! "Software"), to deal in the Software without restriction, including
!! without limitation the rights to use, copy, modify, merge, publish,
!! distribute, sublicense, and/or sell copies of the Software, and to
!! permit persons to whom the Software is furnished to do so, subject
!!
!! the following conditions:
!!
!! The above copyright notice and this permission notice shall be
!! incluudedin all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND INFRINGEMENT
!! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
!!
!! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
!! TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
!! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!!
!! ------------------------------------------------------------------
!! Fortran 95 norm
!! ------------------------------------------------------------------

! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!  VANDERPOL WITH BETA PARAMETER
!  REDUCED TO 2D
!  RUNGE-KUTTA INTEGRATION
!  TIME SCALE : 1 TIME UNIT = 10 KA
!  DETERMINISTIC
!  USES OBLIQUITY AND PRECESSION SEPARATELY
!  AS FORCINGS (GAMMAP AND GAMMAO)
!  AUTHOR : MICHEL CRUCIFIX
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

subroutine vdp_f(state,n,ndim,npar,par,forcing, dforcingdt, deltat,dy,ds,dl)
!     par   (n,6)  : parameters
!--     alpha    = par(:,1)
!--     beta     = par(:,2)
!--     gammapre = par(:,3)
!--     gammaobl = par(:,4)
!--     omega    = par(:,5)
!--     asym     = par(:,6)
!
  implicit none
  integer, intent(in) ::  n
  integer, intent(in) :: ndim, npar !! expect ndim=2, npar=6
  double precision, intent(in) :: par(n,npar),state(n,ndim) 
  double precision, intent(in) :: deltat
  double precision, intent(in) :: forcing(3),dforcingdt(3)
  double precision, intent(inout) :: ds(n,ndim),dl(n,ndim)
  double precision, intent(out)   :: dy(n,ndim)

  integer :: i,k,m 
  double precision da(n,ndim)
  double precision p(n), pp(n), dt(n)


!     variables for runge-kutta integration
  double precision kk1(n,ndim), kk2(n,ndim), kk3(n,ndim)
  double precision kk4(n,ndim)
  double precision u(n,ndim), dadu(n,ndim,ndim)
  double precision dadt(n,ndim), du(n,ndim)
  double precision da1dt(n)


  double precision alpha(n), beta(n), gammapre(n) 
  double precision gammaobl(n), omega(n), asym(n),f(n)
      

  u  = state(:,:)
  du = ds(:,:)

  alpha    = par(:,1)
  beta     = par(:,2)
  gammapre = par(:,3)
  gammaobl = par(:,4)
  omega    = par(:,5)
  asym     = par(:,6)

  dt = deltat / omega

  do m=1,n
     if (asym(m).gt. 0.) then
       da1dt(m) = -exp(asym(m)*(gammapre(m)*forcing(1) + gammaobl(m)*forcing(3)))
        f(m)     = (da1dt(m) + 1.)/asym(m)
      else
       da1dt(m) =  - 1.
       f(m)     =  - ( gammapre(m)*forcing(1) + gammaobl(m)*forcing(3))
      endif
  enddo

  call p_(n,u(:,2),p)
  call pp_(n,u(:,2),pp)

  !  first-order time derivatives
  kk1(:,1)   =   -  (u(:,2) + beta) + f
  kk1(:,2)   =    alpha * ( u(:,1) - p)


  !  differential elements for calculating lyapunov

  da (:,1) = - du(:,2) 
  da (:,2) = alpha * (du(:,1) - pp*du(:,2))

  !  jacobian for runge-kutta

  dadu(:,1,1)  =  0
  dadu(:,1,2)  =  -1
  dadu(:,2,1)  =  alpha
  dadu(:,2,2)  =  (- alpha * pp)

  !  time derivatives for runge-kutta
  dadt(:,1)  =  da1dt * (gammapre * dforcingdt(1) + gammaobl * dforcingdt(3)) 
  dadt(:,2)  = 0. 

    call rk4(kk1, dt, dadu, dadt, n, ndim, dy)
  ! the system itself is integreted with runge-kutta
  ! its tangent is integrated with euler.
  ! this could be fixed in a nicer version

  do k=1,ndim
  !    dy(:,k) = dt(:) * (  kk1(:,k) ) 
    dl(:,k) = 1. * dt(:) * da(:,k)
  enddo  

  ! local subroutines (might beter be inlined)
  contains
    subroutine p_(n,x,p) ! phi(x)
     integer n
     double precision p(n),x(n)
      p=x*x*x/3.-x
    end subroutine

     subroutine pp_(n,x,pp) ! dphidx(x)
     integer n
     double precision pp(n),x(n)
      pp=x*x-1
    end subroutine

      
end
