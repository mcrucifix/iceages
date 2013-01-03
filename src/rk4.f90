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

subroutine rk4 (kk1, dt, dadu, dadt, n, ndim, dy )

! runge-kutta fourth order scheme
! works 'n' integration in parellel, each with its specific time step
! kk1 : first order derivatives
! dadu : jacobian ; dadt : time derivatives of kk1
! ndim : ndimensions of the system
! dy   : increment (output)

! author : M. Crucifix (2012)


  integer, intent(in) :: n,ndim
  double precision, intent(in) :: dt(n)
  double precision, intent(in) :: kk1(n,ndim)
  double precision, intent(in) :: dadt(n,ndim), dadu(n,ndim,ndim)
  double precision, intent(out) :: dy(n,ndim)
  double precision             :: kk2(n,ndim), kk3(n,ndim), kk4(n,ndim)

  do k=1,ndim
     kk2 (:,k) = kk1(:,k) + 0.5*dadt(:,k) * dt
     do i=1,ndim
      kk2(:,k) = kk2(:,k) + 0.5 * dadu(:,k,i) * kk1(:,i) * dt
     enddo
  enddo

  do k=1,ndim
     kk3 (:,k) = kk1(:,k) + 0.5*dadt(:,k) * dt
     do i=1,ndim
      kk3(:,k) = kk3(:,k) + 0.5*dadu(:,k,i) * kk2(:,i) * dt
     enddo
  enddo

  do k=1,ndim
     kk4 (:,k) = kk1(:,k) + dadt(:,k) * dt
    do i=1,ndim
      kk4(:,k) = kk4(:,k) +  dadu(:,k,i) * kk3(:,i) * dt
    enddo
  enddo

  do k=1,ndim
    dy(:,k) = 1./6. * dt(:) * (  kk1(:,k) + 2*kk2(:,k) + 2*kk3(:,k) +  kk4(:,k))
  enddo

end subroutine


