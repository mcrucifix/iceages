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
! modif 20.03.13 : u1, u2 and dw go n x ndim
! stocastic version of vdp_f.
! see references in that file


subroutine daruka_s (n, ndim, npar, state, par, dt,sdt, f, dfdt, u1, u2, dy)

  implicit none
  ! interface
  integer, intent(in) ::  n, ndim, npar ! (expect : ndim=2, npar = 6)
  double precision, intent(in), dimension(n,ndim), target ::  par,state
  double precision, intent(in), dimension (n) ::  dt,sdt
  double precision, intent(in), dimension(3) :: f, dfdt
  double precision, intent(in), dimension(n,ndim) :: u1, u2
  double precision, intent(out) ::  dy(n,ndim)
  ! end interface

  ! internal, will be passed to stochastic integrator
  double precision, dimension(n,ndim) :: dw
  double precision, dimension(n,ndim) :: a, b, dadt
  double precision, dimension(n,ndim,ndim) :: dadx

  ! aliases for parameters and state
  double precision, pointer, dimension(:) ::  alpha,kappa, lambda, gammapre, gammaobl, omega, sigmax, sigmay
  double precision, dimension(n) :: Vprime, Vprimeprime, forcing

  dadx=0.
  dadt=0.

  b = 0.

 ! redefine parameter and state names for readability
  alpha    => par(:,1)
  kappa    => par(:,2)
  lambda   => par(:,3)
  gammapre => par(:,4)
  gammaobl => par(:,5)
  omega    => par(:,6)
  sigmax    => par(:,7)
  sigmay    => par(:,8)

  forcing =  gammapre * f(1) + gammaobl * f(3)
  Vprime = alpha - state(:,2) + state(:,2)**3
  Vprimeprime = - 1 + 3*state(:,2)**2

   a(:,1)   =  ( -Vprime - kappa * state(:,1) - state(:,2)*forcing  )  
   a(:,2)   =  lambda * state(:,1) 



! diffusion terms 

   b(:,1)   = sigmax
   b(:,2)   = sigmay

!  first order derivatives of drift with respect to x
   dadx(:,1,1)  =  - kappa  
   dadx(:,1,2)  =  - ( Vprimeprime - forcing ) 
   dadx(:,2,1)  =  0
   dadx(:,2,2)  =  lambda 

   dadt(:,1)  =  -state(:,2) * (gammapre * dfdt(1) + gammaobl * dfdt(3))  

!  call semi-implicit stochastic integrator

  dw = u1
  call si_stoch(n,ndim,npar,state,par,a, dadx, dadt, b,dt,sdt,dw, dy)
end subroutine daruka_s


