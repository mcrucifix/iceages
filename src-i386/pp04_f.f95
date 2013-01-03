subroutine pp04_f(state,n,ndim,npar,par,f,dfdt,deltat,dy,ds,dl)
 
!     ------------------------------------------
!     D. PAILLARD AND F. PARRENIN, THE {A}NTARCTIC ICE SHEET AND THE
!     TRIGGERING OF DEGLACIATIONS,  EARTH PLANET. SC. LETT., 227, 263-271
!     (2004)
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

  ! -- local variables
    integer m
    double precision i65, i60
    double precision gamma(n), omega(n)
    double precision a(n), v(n), c(n)
    double precision h(n), ff(n), dt(n)
    double precision da(3,n)
         
  ! -- parameters 
    ! insolation
    double precision, parameter  :: gp = 0.7636
    double precision, parameter ::  go = 0.4752
    ! relaxation times
    double precision, parameter :: tauv = 1.5
    double precision, parameter :: tauc = .5
    double precision, parameter :: taua = 1.2
    ! other constants
    double precision, parameter :: xx   = 1.3
    double precision, parameter :: yy   = 0.5
    double precision, parameter :: zz   = 0.8

    double precision, parameter :: alpha = 0.15
    double precision, parameter :: beta  = 0.5
    double precision, parameter :: ggamma = 0.5
    double precision, parameter :: ddelta = 0.4
    ! original value = 0.3, but needed to be adjusted
    ! to approximately reproduce their results
    ! cf. pers. communication by Paillard
    ! acknowledging error in manuscript
    double precision, parameter :: aa = 0.40
    double precision, parameter :: bb = 0.7
    double precision, parameter :: cc = 0.01
    double precision, parameter :: d = 0.27

  !  for test only:
  !  double precision, parameter :: d = 0.25

  v=state(:,1)
  a=state(:,2)
  c=state(:,3)

  gamma  = par(:,1)
  omega  = par(:,2)

  dt = deltat / omega 

  i65 = gp * f(1) + go * f(3)
  i60 = -0.421820 * f(1) + 0.716977 * f(2) + 0.247518 * f(3)

  ! -----------------------------------------
  ! ocean circulation switch
  ! -----------------------------------------

  ff = aa*v - bb*a + d - gamma * cc * i60
  h = ggamma
  where (ff > 0) h = 0.

  ! -----------------------------------------
  ! relaxation equations
  ! -----------------------------------------
  ! simple euler forward integration
  ! no calculation of lyapunov exponent
  !--------------------------------


  dy(:,1) = (- xx * c - yy * gamma * i65 + zz - v)/tauv * dt
  dy(:,2) = (v - a ) / taua * dt
  dy(:,3) = (alpha * gamma * i65 - beta*v + h + ddelta - c) / tauc * dt


  ! -----------------------------------------
  ! test for non-negative ice area and valume
  ! -----------------------------------------

  where ((dy(:,1) ) < -v)  dy(:,1) = -v
  where ((dy(:,2) ) < -a)  dy(:,2) = -a
  where ((dy(:,3) ) < -c)  dy(:,3) = -c

end subroutine
    
