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

! -----------------------------------
! implicit order 1.5 scheme
! reference : Kloeden ... check page number
! author : Crucifix (2011)
! modif : 21.03.13 : dw go n x nidm
! ---------------------------------
!       construct the matrix which is on the left hand side of the
!       x[t+1] vector, on the left hand side of the equation


  subroutine si_stoch(n,ndim, npar, state, par, a, dadx, dadt, b, dt, sdt, dw, dy  )

  ! interface
  implicit none
  integer, intent(in) ::  n, ndim, npar ! 
  double precision, intent(in), dimension(n,ndim) ::  state,dadt, a, b
  double precision, intent(in), dimension(n,npar) ::  par
  double precision, intent(in), dimension(n,ndim,ndim) ::  dadx
  double precision, intent(in), dimension (n) :: dt,sdt
  double precision, intent(in), dimension (n,ndim) :: dw 
  double precision, intent(out), dimension(n,ndim) :: dy

  ! auxiliary variables
  integer i,j,m, info
  integer, dimension(n) :: ipiv
  double precision, dimension(ndim) :: rhstmp
  double precision, dimension(n,ndim) :: rhs
  double precision, dimension(n,ndim,ndim) :: amat
  double precision, dimension(ndim,ndim) :: amattmp

    amat = 0.
    forall (m=1:n)
      forall (i=1:ndim) amat(m,i,i)=1.
        amat(m,:,:) = amat(m,:,:) - 0.5*dadx(m,:,:) * dt(m)  
    end forall


!    now constructs the right hand side vector
    forall (i=1:ndim)  
     rhs (:,i) =  state(:,i) + (a(:,i) + 0.5 * (dadt(:,i)) * dt) * dt
     forall (j=1:ndim)
       rhs (:,i) =  rhs(:,i) - 0.5*(dadx(:,i,j)*state(:,j)) * dt 
     end forall 
     rhs (:,i) =  rhs(:,i) + b(:,i) * dw(:,i) *sdt
    end forall
!     now the equation to be solved is  amat * x = rhs
!     use lapack

   do m=1,n
      amattmp = amat(m,:,:)
      rhstmp  = rhs(m,:)
      call dgesv( ndim,  1,  amattmp,  ndim,  ipiv,  rhstmp,  ndim,  info)
      rhs(m,:)=rhstmp(:)
   enddo

    dy = rhs - state
    end
 
