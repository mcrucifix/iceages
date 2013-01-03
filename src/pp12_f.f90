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

!  rewritten by this author based on : 
! 
!  PARRENIN AND PAILLARD 2012
!  Parrenin, F. and Paillard, D.: Terminations VI and VIII 
!  (∼ 530 and ∼ 720 kyr BP)
!  tell us the importance of obliquity and precession in the triggering of
!  deglaciations, Clim. Past, 8, 2031-2037, doi:10.5194/cp-8-2031-2012, 2012.  
!  Terminations VI and VII (-530 and -720 kyrBP)
!  The time unit is 10 ka
! 
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function  f ( x )
  double precision f,x
  double precision, parameter  :: a    = 0.68034
  if (x > 0) then 
     f = x
  else
     f = x + sqrt ( 4 * a * a + x * x ) - 2 * a
  endif
  f =  f 
end function

function dfdx ( x )
  double precision dfdx, x
  double precision, parameter  :: a    = 0.68034
  if (x > 0) then 
     dfdx = 1
  else
     dfdx = 1 + x / sqrt ( 4 * a * a + x * x ) 
  endif
  dfdx =  dfdx 
end function


subroutine pp12_f(state,n,ndim,npar,par,forcing, dforcingdt, deltat,dy,ds,dl)

  implicit none
  integer k,j,l,m
  integer, intent(in)  :: n, ndim, npar
  double precision, intent(in) :: par(n,npar),state(n,ndim) 
  double precision, intent(in) :: deltat
  double precision, intent(in) :: forcing(3), dforcingdt(3)
  double precision, intent(inout) :: ds(n,ndim),dl(n,ndim)
  double precision, intent(out)   :: dy(n,2)


  ! constants to multiply precession and obliquity
  ! to obtain normalised values

  double precision, parameter :: gp = 0.8536269 
  double precision, parameter :: go = 1.094381 

  !     parameters (time unit = 10 ka ! ) 
  double precision, parameter ::  aesi = 14.561
  double precision, parameter ::  aeco = 3.87
  double precision, parameter ::  ao   = 11.3662
  double precision, parameter ::  ag   = 9.7825
  double precision, parameter ::  ad   = -7.469
  double precision, parameter ::  tdm1 = 1.199  ! 1. / td
  double precision, parameter ::  kesi = 14.6348
  double precision, parameter ::  keco = 2.28061
  double precision, parameter ::  ko   = 18.5162
  double precision, parameter ::  v0   = 122.918
  double precision, parameter ::  v1   = 3.10301

  !     local variables

  double precision f

  double precision localstate(n,ndim)
  double precision v(n), s (n)

  double precision cfobl, dfdx
  double precision, dimension(n) ::  gamma, omega, ggp
  double precision, dimension (n) :: dt, kk1, kk2, kk3, kk4, dkk1dt, dkk1dv
  double precision, dimension (n) ::  threshold, ifpre, ifcpre, ifobl
  double precision, dimension(n) ::  difpredt, difcpredt, difobldt

  localstate = state

  v = state(:,1)
  s = state(:,2)
   
  ! needed to reproduce results exactly (small shift of obliquity
  ! because the oblquity means were not calculated over the same periods. 
  cfobl = forcing(3) - 0.0298181

  gamma  = par(:,1)
  omega  = par(:,2)

  ggp = gamma * gp

  !     non-linear transformation of the precession forcing

  ifpre =  ggp * ( f(forcing(1))  - 0.14308 ) / 1.22139
  ifcpre = ggp * ( f(forcing(2)) - 0.14308 ) / 1.22139
  ifobl  = gamma * go * ( f(cfobl) - 0.1298 ) / 1.0786

  difpredt = ggp * (dforcingdt(1) * dfdx(forcing(1)) ) / 1.22139
  difcpredt = ggp * (dforcingdt(2) * dfdx(forcing(2)) ) / 1.22139
  difobldt = gamma * go * dforcingdt(3) * dfdx(forcing(3))

  !     s can only be 0 (glaciation) or 1 (deglaciation)

  threshold = kesi * ggp * forcing(1) + keco * ggp * forcing(2) & 
              + ko * go * gamma * cfobl

  
  where (s > .5 .and. (threshold  < v1) )  s = 0.     ! "enter glaciation"
  where (s < .5 .and. (threshold + v >  v0) )  s = 1. ! "enter deglaciation"

  kk1 = - aesi * ifpre - aeco * ifcpre - ao * ifobl
  dkk1dt = -aesi * difpredt - aeco * difcpredt - ao * difobldt

  where (s >  .5)  ! "deglaciation mode"
    kk1 = kk1 + ad - v * tdm1
    dkk1dv = -tdm1
  elsewhere
    kk1 = kk1 + ag ! "glaciation mode"
    dkk1dv = 0
  endwhere

  dt = deltat / omega

  kk2 = kk1 + 0.5 * dt * (dkk1dt + dkk1dv * kk1)
  kk3 = kk1 + 0.5 * dt * (dkk1dt + dkk1dv * kk2)
  kk4 = kk1 + 0.5 * dt * (dkk1dt + dkk1dv * kk3)

  localstate(:,1) = v + 1./6. * dt * (kk1 + 2 * (kk2  + kk3 ) + kk4)

  where (localstate (:,1) < 0) localstate(:,1) = 0

  localstate(:,2) = s

  dy = localstate - state

eND
