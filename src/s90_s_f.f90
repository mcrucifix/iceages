! Copyright (c) 2012 Michel Crucifix <michel.crucifix@uclouvain.be>
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

subroutine s90_s (n,ndim,npar,state,par,dt, sdt, f,dfdt, u1, u2, dy)
implicit none
integer k,i,j,l,m
integer, intent (in) :: n, ndim, npar ! expect ndim=3, npar=3
double precision, intent(in), dimension(n,ndim), target ::  par,state
double precision, intent(in), dimension (n) ::  dt,sdt
double precision da(n,ndim)
double precision, intent(in), dimension(3) :: f, dfdt
double precision, intent(in), dimension(n,ndim) :: u1, u2
double precision dy(n,ndim)

double precision gp,go
parameter (gp = 0.7639)
parameter (go = 0.4756)

! internal, will be passed to stochastic integrator
double precision, dimension(n,ndim) :: dw
double precision, dimension(n,ndim) :: a, b, dadt
double precision, dimension(n,ndim,ndim) :: dadx

double precision, pointer, dimension(:) ::   gamma, w,  omega, sigma

!! parameters propre to saltzman and maasch
double precision p, q, r, v,  s, uu

parameter (p=1.0)
parameter (q=2.5)
parameter (r=0.9) 
parameter (v=0.2)
! parameter (w=0.5)  ! w now input as a parameter
parameter (uu=0.6)
parameter (s=1.0) 


 ! redefine parameter and state names for readability
  gamma    => par(:,1)
  omega    => par(:,2)
  w        => par(:,3)
  sigma    => par(:,4)

  
  ! model (drift terms)
  a(:,1) =  - state(:,1) - state(:,2) - v*state(:,3) - uu * gamma * ( gp* f(1) + go*f(3) )
  a(:,2) =  -p * state(:,3) + r * state(:,2) + s * state(:,3)*state(:,3) &
            - w*state(:,2)*state(:,3) - state(:,2)*state(:,3)*state(:,3)
  a(:,3) =  -q * (state(:,1) + state(:,3))

  ! 
  dadx(:,1,1)  =  -1
  dadx(:,1,2)  =  -1
  dadx(:,1,3)  =  -v
  dadx(:,2,1)  =  0
  dadx(:,2,2)  = r - w * state(:,3) - state(:,3)*state(:,3)
  dadx(:,2,3)  = -p + 2 * s * state(:,3) - w * state(:,2) - 2 * state(:,2) * state(:,3)             
  dadx(:,3,1)  = -q
  dadx(:,3,2)  = 0
  dadx(:,3,3)  = -q

  ! time derivatives for runge-kutta
  dadt(:,1)  = - uu* gamma * (gp * dfdt(1) + go * dfdt(3)) 
  dadt(:,2)  = 0. 
  dadt(:,3)  = 0. 

! diffusion terms 

   b(:,1)   = 0
   b(:,2)   = 0
   b(:,3)   = sigma


!  call semi-implicit stochastic integrator
  dw = u1
  call si_stoch(n,ndim,npar,state,par,a, dadx, dadt, b,dt, sdt, dw, dy)


end
      
