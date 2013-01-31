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
!! ATTENTION 
!  SEE README.TXT
!  The different model parameters use units
!  such that one time unit = 10 ka
!  Written by this author based on : 
! E. Tziperman et al., Consequences of pacing the {P}leistocene 100 kyr ice
! ages by nonlinear phase locking to {M}ilankovitch forcing,  Paleoceanography,
! 21, PA4206  2006
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


subroutine t06_f(state,n,ndim,npar,par,f, dfdt, deltat,dy,ds,dl)

  implicit none
  integer, intent(in) ::  n,ndim,npar ! expect :ndim=2, npar=2
  double precision, intent(in) ::  f(3), dfdt(3)
  double precision, intent(in) ::  deltat
  double precision, intent(in) ::  par(n,npar),state(n,ndim) 
  double precision, intent(inout) :: ds(n,ndim),dl(n,ndim)
  double precision, intent (out) ::  dy(n,2)

  double precision :: x,y
  double precision a(ndim), da(ndim)

  integer k,j,l,m

  double precision, parameter :: gp = 0.7639, go = 0.4756

  double precision alpha, beta, gamma, omega, i0, i1, v0
  double precision  dx, asi 
  
  double precision, parameter :: cfactor = 86400. * 360. * 10000. *1.e-9
  double precision :: p0,s,sm
 ! double precision, parameter :: p0 = 0.26*1e-9*86400*360*10000
  double precision, parameter :: kk  = 0.7 / 4.
 ! double precision, parameter :: s  = 0.23*1e-9*86400*360*10000
 ! double precision, parameter :: sm = 0.03*1e-9*86400*360*10000


 ! to do : include comment about the units of p0, s,  and sm 
 ! and the role of th scaling factor. 

  do m=1,n

    x=state(m,1)
    y=state(m,2)

    gamma  = par(m,1)
    omega  = par(m,2)
    p0     = par(m,3) * cfactor
    s      = par(m,4) * cfactor
    sm     = par(m,5) * cfactor

    if (y .ne. 0 .and. y.ne. 1) y = 1.

    if (y .eq. 0. .and. (x.gt.45)) then 
      y = 1.
    elseif (y .eq. 1. .and. (x.lt.3)) then
      y = 0.
    endif

    asi = 0.

    if (y .eq. 0.) asi = 0.
    if (y .eq. 1.) asi = 0.46

    dx = (p0 - kk* x)*(1-asi) - (s + sm * gamma* (gp * f(1) + go * f(3)))

     !! attention : this is an iteration

     dy(m,1) = x + dx * deltat / omega
     if (dy(m,1) .lt. 0) dy(m,1) = 0.
     dy(m,2) = y
  enddo

  !! i repeat : attention : this is now an iteration
  dy = dy - state

end
