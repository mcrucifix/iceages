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

subroutine vcv18_f(state,n,ndim,npar,par,f,dfdt,deltat,dy,ds,dl)
 
!     TODO: mettre la reference Verbitsky et al. Climate of the past 2018
!     ------------------------------------------
!
!     NPAR = 2 
!     NDIM = 3
!     ------------------------------------------

  implicit none

  integer, intent(in) :: n
  integer, intent(in) :: ndim, npar
  double precision, intent(in) :: par(n,npar),state(n,ndim) 
  double precision, intent(in) :: deltat
  double precision, intent (inout) :: ds(n,3), dl(n,3)
  double precision, intent (in)    :: f(3),dfdt(3)
  double precision, intent(out)    :: dy(n,3)

  ! variables for runge-kutta integration
  double precision kk1(n,ndim)
  double precision u(n,ndim), dadu(n,ndim,ndim), dadt(n,ndim), du(n,ndim)
  ! -- local variables
    integer m
    double precision gamma(n), omega(n)
    double precision I(n), T(n), O(n)
    double precision h(n), ff(n), dt(n)
    double precision da(3,n), i65(n), di65dt(n)
         
  ! -- parameters 
    ! insolation
    ! todo: corriger: 
    double precision, parameter :: gp = 0.7060
    double precision, parameter :: gc = -0.4066
    double precision, parameter ::  go = 0.4386

    double precision, parameter ::  z = 1.
    double precision, parameter ::  aa = 0.065
    double precision, parameter ::  k = 0.005
    double precision, parameter ::  cc = 0.042
    double precision, parameter ::  alpha = 2.
    double precision, parameter ::  beta = 2.
    double precision, parameter ::  g1 = 0.
    double precision, parameter ::  g2 = 0.21
    double precision, parameter ::  g3 = 0.30
    double precision, parameter ::  S0 = 12
!    double precision, parameter ::  gamma = 0.11

  I=state(:,1)
  T=state(:,2)
  O=state(:,3)

  gamma  = par(:,1)
  omega  = par(:,2)

  dt = deltat / omega * 10.   ! units are 10ka (hence the /10 on second line)

  i65 = gamma * ( gp * f(1) + gc * f(2) + go * f(3) )
  di65dt =  0* gamma * ( gp * dfdt(1) - gc * dfdt (2) + go * dfdt(3) ) / 10.


  !       first-order time derivatives
  kk1(:,1) =  0.8*(z)*I**0.75*(aa-i65-k*O-cc*T)
  where (I + kk1(:,1)*dt < 0.1) 
    kk1(:,1) = 0
    I(:) = 0.1
  endwhere
  kk1(:,2) =  (1/z)*I**(-0.25)*(aa-i65-k*O)*(alpha*O+beta*(I-S0)-T)
  kk1(:,3) =  g1 - g2*(I-S0)-g3*O
  

  ! -----------------------------------------
  ! test for non-negative ice area and valume
  ! -----------------------------------------



  ! dadu and dadt are not specified
  dadu  = 0
  dadt  = 0 

  dadu(:,1,1)  =  0.75*kk1(:,1)*I**(-1)
  dadu(:,1,2)  =  -cc*0.8*(z)*I**0.75
  dadu(:,1,3)  =  -k*0.8*(z)*I**0.75
  dadu(:,2,1)  =  -0.25*kk1(:,2)*I**(-1)
  dadu(:,2,2)  = (1/z)*I**(-0.25)*(aa-i65-k*O)*(-1)
  dadu(:,2,3)  = (1/z)*I**(-0.25)*(-k*(alpha*O+beta*(I-S0)-T) + alpha*(aa-i65-k*O))
  dadu(:,3,1)  = -g2
  dadu(:,3,2)  = 0
  dadu(:,3,3)  = -g3

  dadt(:,1) = 0.8*(z)*I**0.75*( - di65dt ) 
  dadt(:,2) = (1/z)*I**(-0.25)*( - di65dt )*(alpha*O+beta*(I-S0)-T)
  dadt(:,3) = 0 



  ! runge-kutta fourth order scheme

  call rk4(kk1, dt, dadu, dadt, n, ndim, dy)

end subroutine
    
