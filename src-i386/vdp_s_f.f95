! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! STOCHASTIC PROPAGATION OF NOISE VAN DER POL
! OSCILLATOR, WITH DIFFERENT NUMERICAL
! SOLUTIONS TO BE SELECTED AT COMPILATION
!

! with IMPLICIT you may ask or not the LAPACK
! routine, in which case you need to compile
! with -lblas -llapack
! Probably does not improve performance in the
! 2-d case but could be useful for higher-dim
! problems 
! 
! Most of the code was crafted for constant additive
! noise (B CONSTANT), but there are some commented lines to be
! activated for more general settings. 

! ANY stochastic differential equations is a particular case
! which explains while, unlike ODEs, there is no general
! solver

! MUST BE LINKED WITH RANUT AND RANNW, available in the directory
! for generation of normal-distributed random numbers

! see ./compile.sh for a sample of compile script
! I haven't made a proper Makefile

! NOTE : THE TIME UNIT FOR THIS PROGRAMME IS
! 1 TIME UNIT = 10 ka (see scale time)

! original code Michel Crucifix, 1 july 2011 
! adapted : 26 septembre 2011
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

!     PAR   (NPAR,5)  : PARAMETERS
!--     ALPHA    = PAR(:,1)
!--     BETA     = PAR(:,2)
!--     GAMMAPRE = PAR(:,3)
!--     GAMMAOBL = PAR(:,4)
!--     OMEGA    = PAR(:,5)
!--     ASYM     = PAR(:,6)
!
subroutine vdp_s (n, state,par,f,dfdt, deltat,u1, u2, dy, ndim, npar)
  implicit none
  integer k,m, i
  integer, intent(in) ::  ndim, npar ! (expect : ndim=2, npar = 6)
  integer n

!     variables specitif to implicit scheme
  double precision  amat(n,ndim,ndim), rhs(n,ndim)
  integer           ipiv(ndim), info

  double precision par(n,npar),state(n,ndim)
  double precision deltat
  double precision f(3),dfdt(3)
  double precision a(n,ndim)
  double precision dadx(n,ndim,ndim), dadt(n,ndim)
  double precision amattmp(ndim,ndim), rhstmp(ndim)

  double precision b(n,ndim), dbdt(n,ndim)
  double precision dy(n,ndim)


!     variables for time-steps
  double precision dt(n), sdt(n), dt2(n), dt32(n), sq3
 
!     dummay variables for generation of random numbers
  integer nwork
  double precision u1(n), u2(n), dw(n), dw2(n), dz(n)

!     model parameters
  double precision alpha(n), beta(n), omega(n), sigma(n)
  double precision gammapre(n), gammaobl(n)

  sq3  = sqrt(1./3.)

!     random numbers
  dw  = u1

  dadx=0.
  dbdt=0.

  b = 0.

  alpha    = par(:,1)
  beta     = par(:,2)
  gammapre = par(:,3)
  gammaobl = par(:,4)
  omega    = par(:,5)
  sigma    = par(:,6)

  dt = deltat / omega
  dt2 = dt*dt
  sdt = sqrt(dt)
  dt32 = dt**(1.5)

!  the model : drift terms
   a(:,1)  = -  (state(:,2) + beta) - gammapre*f(1) - gammaobl*f(3)
   a(:,2)  =    alpha * ( state(:,1) - &
     (  state(:,2)*state(:,2)*state(:,2) / 3.  - state(:,2)) )

! diffusion terms 

   b(:,2)   = sigma

!  first order derivatives of drift with respect to x
   dadx(:,1,2)  =  -1
   dadx(:,2,1)  =  alpha
   dadx(:,2,2)  =  - alpha * (state(:,2) * state(:,2) - 1)

!                     ... of drift with respect to time
   dadt(:,1)  = gammapre * dfdt(1) + gammaobl * dfdt(3)


!   ... if you are order15 you also need
!   the second-order derivatives of the drift
!   except in the implicit or pred-corr schemes
!   here, given that we only have stoch. noise in the
!   second eq. you only need the second-order derivatives
!   with respect to x(2)

!   diffusion first order derivatives
!   default to zero. no need to update 

!   dbdt(1) = 0
!   dbdt(2) = 0


! -----------------------------------
! implicit scheme
! ---------------------------------
!       construct the matrix which is on the left hand side of the
!       x[t+1] vector, on the left hand side of the equation

    amat = 0.

    forall (m=1:n)
      forall (i=1:ndim) amat(m,i,i)=1.
        amat(m,:,:) = amat(m,:,:) - 0.5*dadx(m,:,:) * dt(m)
    end forall


!    now constructs the right hand side vector
      
    forall (i=1:ndim)  
     rhs (:,i) =  state(:,i) + (a(:,i) - 0.5*(dadx(:,i,1)*state(:,1) &
                     + dadx(:,i,2)*state(:,2) -dadt(:,i) * dt)) * dt &
                     + b(:,i) * dw(:) *sdt
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
    
