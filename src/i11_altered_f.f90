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

!  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!! ATTENTION 
!  SEE README.TXT
!  John Z. Imbrie, Annabel Imbrie-Moore, and Lorraine E. Lisiecki, 
!  A phase-space model for Pleistocene ice volume,  
!  Earth and Planetary Science Letters, 307, 94--102  2011
!  Note : the stability function was fitted with a 3-degree polynomial
!         for robustess, ease of coding and computational efficiency
!  The different model parameters use units
!  such that one time unit = 10 ka
!  Author of code : Michel Crucifix 2012
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

subroutine i11_altered_f(state,n,ndim,npar,par,forcing,dforcingdt, deltat,dy,ds,dl)

  implicit none
  integer, intent(in) :: n
  integer, intent(in) :: ndim, npar
  double precision, intent(in) :: par(n,npar),state(n,ndim) 
  double precision, intent(in) :: deltat
  double precision, intent(in) :: forcing(3),dforcingdt(3)
  double precision, intent(inout) :: ds(n,ndim),dl(n,ndim)
  double precision, intent(out)   :: dy(n,2)

  integer k,j,l,m

  !     constants 
  double precision gp,go, pi180

  ! scale parameters so that precession and obl are normalised
  parameter (gp = 0.8536269)
  parameter (go = 1.094381 )
  parameter (pi180 = 0.01745329)

  !     local variables
  double precision localstate(n,ndim)
  double precision y(n), yp (n)


  double precision gamma(n), omega(n)
  double precision h1obl(n), a(n), lag(n)
  double precision h2si(n), h3co(n), f(n), d(n), trans(n)

  localstate = state

  y =state(:,1)
  yp=state(:,2)

  gamma  = par(:,1)
  omega  = par(:,2)

  h1obl = (0.05 - y * 0.005) * go * forcing(3)

  where (y .lt. 0)
      trans = 0.135 + 0.07 * y
  elsewhere
      trans = 0.135
  endwhere
  a     = 0.07 
  lag   = 10 * pi180

  h2si  = a * gp * forcing(1) * sin(lag) 
  h3co  = a * gp * forcing(2) * cos(lag) 
  f     = gamma*(h1obl + h2si + h3co)

  where ( yp .gt. trans)
     d = -1. + 3. * (f + 0.28) / 0.51
     dy(:,1) = yp
     dy(:,2) = 0.5 * (0.25 * 0.25 * (d - y) - yp*yp / (d-y)  )

     localstate(:,1) = state(:,1) + deltat / omega * dy(:,1) * 10.
     localstate(:,2) = state(:,2) + deltat / omega * dy(:,2) * 10.

  elsewhere 

     dy(:,1) =  &
      -0.036375 + y * ( 0.014416 + y *  & 
      (0.001121 - y * 0.008264)) + 0.5*f 
      localstate(:,1) =  state(:,1) + dy(:,1) * deltat / omega * 10.
      localstate(:,2) =  dy(:,1)
  endwhere

  where (y .gt. d .and. yp .gt. trans)
      localstate(:,2) =  0.
  endwhere

  dy = localstate - state

eND
