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

!! modif history
!!  18 march 2013 : add rmean optional parameter 
!!                  useful for proposals in particle filter algorithms

!!  20 march 2013 : go for multidimensional u1 and u2 (n x ndim)
!!                  and rmean is accordingly a ndim-vector
!!                  isum facility temporarily discarded

!!  22 march 2013 : output generated wiener processes

subroutine propagate_s (n, state, par, scaletime, t0, t1, deltat, ix, nap, nao,  &
                        amppre, omepre, angpre,  &
                        ampobl, omeobl, angobl,  & 
                        ndim, npar, model_s, isum, rmean, dw)

!     MAIN ROUTINE CALLED FROM R 
!     IN PARAMETERS 
!     -------------
!     T0 , T1      : INITIAL AND FINAL TIME
!     NAP, NAO     : NUMBER OF TERMS RETAINED IN COEF. EXPANSION
!                     PRECESSION AND OBLIQUITY
!     SEED    : RANDOM NUMBER SEED. MUST NOT BE 0. 
!     NDIM    : Number of dimensions
!     NPAR    : Number of parameters

!     ISUM    : contructs every random number as a sum
!               of isum random numbers, then divided
!               by sqrt(isum)
!               (see example use for 
!               generation of Brownian tree)
!     INOUT
!     -----
!     state (N,2)  : initial conditions
  ! interface
  implicit none
  integer, intent(in) :: n, ix, nap, nao, ndim, npar, isum
  double precision, intent(in) :: t0,t1, deltat
  double precision, intent(in) :: par(n,ndim), scaletime(n)
  double precision, intent(in), dimension (nap)  ::  amppre, omepre, angpre
  double precision, intent(in), dimension (nao)  ::  ampobl, omeobl, angobl
  double precision, intent(in), dimension (ndim), optional :: rmean
  double precision, intent(inout) :: state(n,ndim)
  double precision, intent(inout) :: dw(n,ndim)

  external model_s

  interface
    subroutine model_s (n, ndim, npar, state, par, dt,sdt, f, dfdt, u1, u2, dy)

    ! interface
      implicit none
      integer, intent(in) ::  n, ndim, npar ! (expect : ndim=2, npar = 6)
      double precision, intent(in), dimension(n,ndim), target ::  par,state
      double precision, intent(in), dimension (n) ::  dt, sdt
      double precision, intent(in), dimension(3) :: f, dfdt
      double precision, intent(in), dimension(n,ndim) :: u1, u2
      double precision, intent(out) ::  dy(n,ndim)
    end subroutine model_s
  end interface

  integer i,j,k, imax

  ! internal (time step)
  double precision :: t, deltat_adjusted
  double precision, dimension (n) :: dt(n), sdt(n)
  double precision f(3), dfdt(3),  dy(n,ndim)


  ! variables for generation of random numbers
  ! note that u2 is presently a dummy (not in use)
  integer, parameter :: nwork=1e4
  double precision dwork(nwork)
  double precision, dimension(n*isum, ndim) :: u1, u2


  ! ------- end variable declaration -----------  ! 


  ! adjusts imax to have an int. numb. of timesteps

  imax = max(int(((t1-t0)/deltat)+0.5),1)
  deltat_adjusted  = (t1-t0)/imax     

  ! main time loop

  dt = deltat / scaletime
!  dt = deltat 
  sdt = sqrt(dt)

! dw = sum of random numbers generated
!        will be used to compute prior / proposal probabilities
  
  do i=1,imax

    t = t0+(i-1)*deltat_adjusted
    call astro(nap, nao, t,amppre, omepre, angpre, ampobl, omeobl, angobl, f, dfdt )
    call rannw(0.d0, 1.d0, ix, u1, n*ndim, dwork, nwork)

   
    ! ------------------------------------------
    ! NOTE : ISUM FACILITY NEEDS TO BE RE-CHECKED
    ! ------------------------------------------

    ! if the user wants to do a brownian tree 
    if (isum .gt. 1) then
     do k=1, ndim
      do j=1,n
       u1(j,k) = sum(u1( ( (j-1)*isum+1 ) : j*isum , k) ) 
      enddo
     enddo 
     u1 = u1 / sqrt(real(isum))
    endif 


    if(present(rmean)) then
     do k=1,ndim
      u1(1:n, k)= u1(1:n,k) + rmean(k)
     enddo
    endif



    do k=1,ndim
      dw(:,k) = dw(:,k) + u1(1:n,k)*sdt
    enddo 

    call model_s (n, ndim, npar, state, par, dt,sdt, f, dfdt, u1, u2, dy)
    state = state + dy

  enddo ! i
 end 
 


