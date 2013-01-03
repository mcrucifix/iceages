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
!  SALTZMAN91-like model
!  RUNGE-KUTTA INTEGRATION
!  TIME SCALE : 1 TIME UNIT = 10 KA
!  DETERMINISTIC
!  USES OBLIQUITY AND PRECESSION SEPARATELY
!  AS FORCINGS (GAMMAP AND GAMMAO)
!  Written by this author based on :
!  B. Saltzman and K. A. Maasch, A first-order global model of late
!  Cenozoic climate. II Further analysis based on a simplification of the
!  CO_2 dynamics,  Climate Dynamics, 5, 201-210  1991
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

subroutine s91_f(state,n,ndim,npar,par,f,dfdt,deltat,dy,ds,dl)
  implicit none
  integer k,i,m
  integer, intent(in)  :: n,ndim,npar
  double precision par(n,npar),state(n,ndim), ds(n,ndim),dl(n,ndim)
  double precision deltat
  double precision da(ndim)
  double precision f(3), dfdt(3)
  double precision dy(n,ndim)
  double precision dt

  double precision gp,go
  parameter (gp = 0.7639)
  parameter (go = 0.4756)

  ! variables for runge-kutta integration
  double precision kk1(ndim), kk2(ndim), kk3(ndim), kk4(ndim)
  double precision u(ndim), dadu(ndim,ndim), dadt(ndim), du(ndim)


  double precision gamma, omega
  
  !! parameters propre to saltzman and maasch
  double precision p, q, r, v,  s,uu

  parameter (p=1.0)
  parameter (q=2.5)
  parameter (r=1.3)
  parameter (v=0.2)
  parameter (s=0.6)
  parameter (uu=0.6)

  do m=1,n

   u  = state(m,:)
   du = ds(m,:)


   gamma = par(m,1)
   omega = par(m,2)

   dt = deltat / omega

    ! first-order time derivatives
    kk1(1) =  - u(1) - u(2) - v*u(3) &
             - uu * gamma * ( gp * f(1) - go * f(3))
    kk1(2) =  -p * u(3) + r * u(2) - s * u(2)*u(2) - u(2)*u(2)*u(2) 
    kk1(3) =  -q * (u(1) + u(3)) 

   ! differential elements for calculating lyapunov

    da (1) = - du(1) - du(2) - v * du(3)
    da (2) = - p * du(3) + r * du(2) + 2 * s * u(2) * du(2)  &
            - 3 * du(2) * u(2) * u(2)
    da (3) = -q * (du(1) + du(3))

    ! jacobian for runge-kutta
    dadu(1,1)  =  -1
    dadu(1,2)  =  -1
    dadu(1,3)  =  -v
    dadu(2,1)  =  0
    dadu(2,2)  = r - 2 * s * u(2) - 3 * u(2)*u(2)
    dadu(2,3)  = -p 
    dadu(3,1)  = -q
    dadu(3,2)  = 0
    dadu(3,3)  = -q

   ! time derivatives for runge-kutta
    dadt(1)  = - gamma * (gp * dfdt(1) + go * dfdt(3))
    dadt(2)  = 0. 
    dadt(3)  = 0. 

    call rk4(kk1, dt, dadu, dadt, n, ndim, dy)

  !  the system itself is integreted with runge-kutta
  !  its tangent is integrated with euler.
  !  this could be fixed in a nicer version

    do k=1,ndim
     dl(m,k) = 1. * dt * da(k)
    enddo  

  enddo

eND
