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
!  Dufing van der pol oscillation
!  Ongoingy study Crucifix - Ditlevsen - Mitsui
!  RUNGE-KUTTA INTEGRATION
!  TIME SCALE : 1 TIME UNIT = 10 KA
!  DETERMINISTIC
!  USES OBLIQUITY AND PRECESSION SEPARATELY
!  AS FORCINGS (GAMMAP AND GAMMAO)
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

subroutine duffing_vdp_f(state,n,ndim,npar,par,forcing, dforcingdt, deltat,dy,ds,dl)
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
  double precision dt(n)


!     variables for runge-kutta integration
  double precision kk1(n,ndim), kk2(n,ndim), kk3(n,ndim)
  double precision kk4(n,ndim)
  double precision u(n,ndim), dadu(n,ndim,ndim)
  double precision dadt(n,ndim), du(n,ndim)

  double precision, dimension(n) :: kappa, mu, omega0, beta, omega,&
                                    gammapre, gammaobl, f
      

  u  = state(:,:)
  du = ds(:,:)

  kappa    = par(:,1)
  mu       = par(:,2)
  omega0   = par(:,3)
  beta     = par(:,4)
  gammapre = par(:,5)
  gammaobl = par(:,6)
  omega    = par(:,7)

  dt = deltat / omega

  f =  gammapre * forcing(1) + gammaobl * forcing(3)

  !  first-order time derivatives
  kk1(:,1)   =   u(:,2) + kappa * u(:,1) - mu / 3. * u(:,1)*u(:,1)*u(:,1) 
  kk1(:,2)   =   omega0 * u(:,1) - beta * u(:,1)*u(:,1)*u(:,1) + f


  !  jacobian for runge-kutta

  dadu(:,1,1)  =  kappa - mu * u(:,1)*u(:,1)
  dadu(:,1,2)  =  1.
  dadu(:,2,1)  =  omega0 - 3 * u(:,1)*u(:,1) * beta
  dadu(:,2,2)  =  0.

  !  differential elements for calculating lyapunov

  da(:,1) = dadu(:,1,1) * du(:,1) + dadu(:,1,2) * du(:,2)
  da(:,2) = dadu(:,2,1) * du(:,1) + dadu(:,2,2) * du(:,2)

  !  time derivatives for runge-kutta
  dadt(:,2)  =  (gammapre * dforcingdt(1) + gammaobl * dforcingdt(3))  !!! suppressed da1dt on 08.03.2018
  dadt(:,1)  = 0. 

  ! runge-kutta fourth order scheme

  do k=1,2
     kk2 (:,k) = kk1(:,k) + 0.5*dadt(:,k) * dt
     do i=1,2
      kk2(:,k) = kk2(:,k) + 0.5 * dadu(:,k,i) * kk1(:,i) * dt
     enddo
  enddo

  do k=1,2
     kk3 (:,k) = kk1(:,k) + 0.5*dadt(:,k) * dt
     do i=1,2
      kk3(:,k) = kk3(:,k) + 0.5*dadu(:,k,i) * kk2(:,i) * dt
     enddo
  enddo

  do k=1,2
     kk4 (:,k) = kk1(:,k) + dadt(:,k) * dt
    do i=1,2
      kk4(:,k) = kk4(:,k) +  dadu(:,k,i) * kk3(:,i) * dt
    enddo
  enddo


  ! the system itself is integreted with runge-kutta
  ! its tangent is integrated with euler.
  ! this could be fixed in a nicer version

  do k=1,2
    dy(:,k) = 1./6. * dt(:) * (  kk1(:,k) + 2*kk2(:,k) + 2*kk3(:,k) +  kk4(:,k))
    dl(:,k) = 1.    * dt(:) * da(:,k)
  enddo  

end
      
