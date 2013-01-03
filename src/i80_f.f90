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
!  Written by this author based on:
! Imbrie and Imbrie, Science, 1980.
!  The different model parameters use units
!  such that one time unit = 10 ka

! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


subroutine i80_f(state,n,ndim,npar,par,f, dfdt, deltat,dy,ds,dl)
  implicit none
  integer k,j,l,m
  integer, intent(in) :: n,ndim, npar
  double precision, intent(in) ::  par(n,npar),state(n,ndim)
  double precision, intent(in) ::  deltat
  double precision, intent(in) ::  f(3), dfdt(3)
  double precision, intent(inout) :: ds(n,ndim), dl(n,ndim)
  double precision, intent(out) ::  dy(n,ndim)

  double precision x,y
  double precision a(ndim), da(ndim)
  
  double precision, parameter :: gp = 0.7639
  double precision, parameter :: go = 0.4756

  double precision alpha, beta, gamma, omega, i0, i1, v0
  
  double precision,  parameter ::  b  = 0.6,  tm = 1.7

  do m=1,n
    y=state(m,1)

    gamma  = par(m,1)
    omega  = par(m,2)

    x = - (gp*f(1) + go*f(3))

    dy(m,1) = y + (1 + sign(b, x-y))/tm  * (x - y)*deltat/omega
    if (dy(m,1) .lt. 0) dy(m,1) = 0.
  enddo

  dy = dy - state

eND
