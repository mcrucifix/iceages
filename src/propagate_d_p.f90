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

subroutine propagate_d (n, state, par, t0, t1,icalclyap, ds, &
                       lyap, nap, nao, &
                       amppre, omepre, angpre,  &
                       ampobl, omeobl, angobl,  & 
                       ndim, npar, ddtlin_f, deltat)  

!     MAIN ROUTINE CALLED FROM R 
!     IN PARAMETERS 
!     -------------
!     T0 , T1      : INITIAL AND FINAL TIME
!     NAP, NAO     : NUMBER OF TERMS RETAINED IN COEF. EXPANSION
!                     PRECESSION AND OBLIQUITY
!     ICALCLYAP    : SET 1 IF YOU WANT THE GREATEST LYAP. EXP. TO BE CALCUATED. 
!     NPAR         : NUMBER OF INPUT PARAME
!     NDIM         : STATE VECTOR DIMENSIO
!     INOUT
!     -----
!     STATE (N,2)  : INITIAL CONDITIONS
!     LYAP         : GREATEST LYAPUNOV COEFFICIENT. SET EQ 1 
!     DS           : DIRECTION OF THE MOST UNSTABLE DIRECTION
!                    (IF NOT KNOWN : MAY BE OBTAINED BY SPIN-UP)


  integer i, imax 
  integer,intent (in) ::  n, npar, ndim, nap, nao
  integer,intent (in) ::  icalclyap
  double precision, intent(in) :: amppre(nap), omepre(nap), angpre(nap)
  double precision, intent(in) :: ampobl(nao), omeobl(nao), angobl(nao)
  double precision, intent(in) ::  t0,t1, par(n,npar), deltat
  double precision, intent(inout) :: state(n,ndim),ds(n,ndim), lyap(n)
  double precision :: dy(n,ndim)
  double precision :: norm(n), dl(n,ndim)
  double precision :: t,deltat_
  double precision :: forcing(3),dforcingdt(3)

  interface
    subroutine ddtlin_f(state,n,ndim,npar,par,forcing, dforcingdt, deltat,dy,ds,dl)
      integer, intent(in) :: n
      integer :: ndim, npar
      double precision, intent(in) :: deltat, forcing(3), dforcingdt(3)
      double precision, intent(in) :: state(n,ndim), par(n,npar),ds(n,ndim)
      double precision, intent(out) :: dy(n,ndim)
      double precision, intent(out) :: dl(n,ndim)
    end subroutine ddtlin_f
  end interface

  !! adjusts imax to have an int. numb. of timesteps
  imax = max(int(((t1-t0)/deltat)+0.5),1)
  deltat_  = (t1-t0)/imax     


!!  call read_precession(nap,amppre,omepre, angpre)
!!  call read_obliquity(nao,ampobl,omeobl, angobl)


  do i=1,(imax)
    t = t0+(i-1)*deltat
    call astro(nap, nao, t,amppre, omepre, angpre, ampobl, omeobl, &
               angobl, forcing, dforcingdt)

    call ddtlin_f(state,n,ndim,npar,par,forcing, dforcingdt, deltat_,dy,ds,dl)

    state = state +  dy
    if (icalclyap.eq.1) then 
       ds = ds +  dl
       if (mod(i,10).eq.0) then
          norm = ds(:,1)*ds(:,1)+ds(:,2)*ds(:,2)
          lyap=lyap+dlog(norm)
          forall (j = 1:n) ds(j,:) = ds(j,:) / norm(j)
       endif
    end if
  enddo

  if (icalclyap.eq.1) then 
     norm = ds(:,1)*ds(:,1)+ds(:,2)*ds(:,2)
     lyap=lyap + dlog(norm)
     forall (j = 1:n) ds(j,:) = ds(j,:) / norm(j)
     lyap = lyap / 2.d0
  endif

end 


