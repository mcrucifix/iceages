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
!  saltzman90 model
!  runge-kutta integration
!  time scale : 1 time unit = 10 ka
!  deterministic
!  uses obliquity and precession separately
!  as forcings (gammap and gammao)
!  written by this author based on 
!  B. Saltzman and K. A. Maasch, A first-order global model of late
!  Cenozoic climate,  Transactions of the Royal Society of Edinburgh Earth
!  Sciences, 81, 315-325  1990
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! #define force_euler

subroutine s90_f(state,n,ndim,npar,par,f,dfdt, deltat,dy,ds,dl)
implicit none
integer k,i,j,l,m
integer, intent (in) :: n, ndim, npar ! expect ndim=3, npar=2
double precision par(n,npar),state(n,ndim), ds(n,ndim),dl(n,ndim)
double precision x,y
double precision deltat,sdeltat
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

double precision dxt, dyt

double precision alpha, beta, gamma, omega
double precision  dx

!! parameters propre to saltzman and maasch
double precision p, q, r, v, w, s, uu

parameter (p=1.0)
parameter (q=2.5)
parameter (r=0.9) 
parameter (v=0.2)
! parameter (w=0.5)  ! w now input as a parameter
parameter (uu=0.6)
 parameter (s=1.0) 

do m=1,n

 u  = state(m,:)
 du = ds(m,:)


 gamma = par(m,1)
 omega = par(m,2)
 w     = par(m,3)

 dt = deltat / omega

  !       first-order time derivatives
  kk1(1) =  - u(1) - u(2) - v*u(3) - uu * gamma * ( gp* f(1) + go*f(3) )
  kk1(2) =  -p * u(3) + r * u(2) + s * u(3)*u(3) - w*u(2)*u(3) - u(2)*u(3)*u(3)
  kk1(3) =  -q * (u(1) + u(3))

  ! differential elements for calculating lyapunov

  da (1) = - du(1) - du(2) - v * du(3)
  da (2) = - p * du(3) + r * du(3) + 2* s * u(3) * du(3) - & 
             w * (u(2) * du(3) + du(2) * u(3)) - du(2) * u(3) * u(3) - &
             2 * u(2) * u(3) * du(3)
  da (3) = -q * (du(1) + du(3))

  ! jacobian for runge-kutta
  dadu(1,1)  =  -1
  dadu(1,2)  =  -1
  dadu(1,3)  =  -v
  dadu(2,1)  =  0
  dadu(2,2)  = r - w * u(3) - u(3)*u(3)
  dadu(2,3)  = -p + 2 * s * u(3) - w * u(2) - 2 * u(2) * u(3)             
  dadu(3,1)  = -q
  dadu(3,2)  = 0
  dadu(3,3)  = -q

  ! time derivatives for runge-kutta
  dadt(1)  = - uu* gamma * (gp * dfdt(1) + go * dfdt(3)) 
  dadt(2)  = 0. 
  dadt(3)  = 0. 

  ! runge-kutta fourth order scheme

  call rk4(kk1, dt, dadu, dadt, n, ndim, dy)

!  the system itself is integreted with runge-kutta
!  its tangent is integrated with euler.
!  this could be fixed in a nicer version

  do k=1,ndim
   dl(m,k) = 1. * dt * da(k)
  enddo  
enddo

end
      
