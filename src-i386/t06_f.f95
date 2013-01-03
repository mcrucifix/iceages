! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!! ATTENTION 
!  SEE README.TXT
!  The different model parameters use units
!  such that one time unit = 10 ka

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
  
  double precision, parameter :: p0 = 0.26*1e-9*86400*360*10000
  double precision, parameter :: kk  = 0.7 / 4.
  double precision, parameter :: s  = 0.23*1e-9*86400*360*10000
  double precision, parameter :: sm = 0.03*1e-9*86400*360*10000

  do m=1,n

    x=state(m,1)
    y=state(m,2)

    gamma  = par(m,1)
    omega  = par(m,2)

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
