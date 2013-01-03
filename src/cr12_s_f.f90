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

! --------------------------------
! STOCHASTIC VERSION OF CR12_F
! See documentation in that file
! --------------------------------

subroutine cr12_s (n, ndim, npar, state, par, dt,sdt, f, dfdt,u1, u2, dy)

  ! interface
  implicit none
  integer, intent(in) ::  n, ndim, npar ! (expect : ndim=2, npar = 6)
  double precision, intent(in), dimension(n,ndim), target ::  par,state
  double precision, intent(in), dimension (n) ::  dt,sdt, u1, u2
  double precision, intent(in), dimension(3) :: f, dfdt
  double precision, intent(out) ::  dy(n,ndim)
  ! end interface

  ! internal, will be passed to stochastic integrator
  double precision, dimension(n) :: dw
  double precision, dimension(n,ndim) :: a, b, dadt
  double precision, dimension(n,ndim,ndim) :: dadx



  ! aliases for parameters and state
  double precision, pointer, dimension(:) :: alpha,  beta0, beta1, beta2, &
                                            delta,  gammapre, gammaobl,  &  
                                            sigmax, sigmay, omega, x, y 

  dadx=0.

  b = 0.

 ! redefine parameter and state names for readability
  alpha     => par(:,1)
  beta0     => par(:,2)
  beta1     => par(:,3)
  beta2     => par(:,4)
  delta     => par(:,5)
  gammapre  => par(:,6)
  gammaobl  => par(:,7)
  omega     => par(:,8)
  sigmax    => par(:,9)
  sigmay    => par(:,10)

  x    => state(:,1)
  y    => state(:,2)

!  the model : drift terms
   a(:,1)  =  beta0 + beta1 * x - beta2 * (x*x - 1)*x - delta * y 
   a(:,2)  =    alpha * delta * ( y * ( 1 - y*y / 3.)  + x )

!  forcing

   a(:,1) = a(:,1) + f(1) * gammapre + f(3) * gammaobl

! diffusion terms 

   b(:,1)   = sigmax
   b(:,2)   = sigmay

!  first order derivatives of drift with respect to state
   dadx(:,1,1)  =  beta1 + beta2 * (3*x*x - 1)
   dadx(:,1,2)  =  - delta
   dadx(:,2,1)  =  alpha * delta
   dadx(:,2,2)  =  - alpha * delta * (y*y - 1)

!                     ... of drift with respect to time
   dadt(:,1)  = gammapre * dfdt(1) + gammaobl * dfdt(3)

!  call semi-implicit stochastic integrator
  dw = u1
  call si_stoch(n,ndim,npar,state,par,a, dadx, dadt, b,dt,sdt,dw, dy)

end subroutine cr12_s


