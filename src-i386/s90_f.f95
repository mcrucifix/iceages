! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!  saltzman90-like model
!  reduced to 2d
!  runge-kutta integration
!  time scale : 1 time unit = 10 ka
!  deterministic
!  uses obliquity and precession separately
!  as forcings (gammap and gammao)
!  author : michel crucifix
!!      force euler by passing -wl,-dforce_euler at compilation
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! #define force_euler

subroutine s90_f(state,n,ndim,npar,par,f,dfdt, deltat,dy,ds,dl)
implicit none
integer k,i,j,l,m
integer, intent (in) :: n, ndim, npar ! expect ndim=3, npar=2
double precision par(n,npar),state(n,ndim), ds(n,ndim),dl(n,ndim)
double precision x,y
double precision deltat,sdeltat
double precision da(ndim)
double precision f(3), dfdt(3)
double precision dy(n,ndim)
double precision dt

double precision gp,go
parameter (gp = 0.7639)
parameter (go = 0.4756)

! variables for runge-kutta integration
double precision kk1(ndim), kk2(ndim), kk3(ndim), kk4(ndim)
double precision u(ndim), dadu(ndim,ndim), dadt(ndim), du(ndim)

double precision dxt, dyt

double precision alpha, beta, gamma, omega
double precision  dx

!! parameters propre to saltzman and maasch
double precision p, q, r, v, w, s, uu

parameter (p=1.0)
parameter (q=2.5)
parameter (r=0.9)
parameter (v=0.2)
parameter (w=0.5)
parameter (uu=0.6)
parameter (s=1.0)

do m=1,n

 u  = state(m,:)
 du = ds(m,:)


 gamma = par(m,1)
 omega = par(m,2)

 dt = deltat / omega

  !       first-order time derivatives
  kk1(1) =  - u(1) - u(2) - v*u(3) - uu * gamma * ( gp* f(1) + go*f(3) )
  kk1(2) =  -p * u(3) + r * u(2) + s * u(3)*u(3) - w*u(2)*u(3) - u(2)*u(3)*u(3)
  kk1(3) =  -q * (u(1) + u(3))

  ! differential elements for calculating lyapunov

  da (1) = - du(1) - du(2) - v * du(3)
  da (2) = - p * du(3) + r * du(3) + 2* s * u(3) * du(3) - & 
             w * (u(2) * du(3) + du(2) * u(3)) - du(2) * u(3) * u(3) - &
             2 * u(2) * u(3) * du(3)
  da (3) = -q * (du(1) + du(3))

  ! jacobian for runge-kutta
  dadu(1,1)  =  -1
  dadu(1,2)  =  -1
  dadu(1,3)  =  -v
  dadu(2,1)  =  0
  dadu(2,2)  = r - w * u(3) - u(3)*u(3)
  dadu(2,3)  = -p + 2 * s * u(3) - w * u(2) - 2 * u(2) * u(3)             
  dadu(3,1)  = -q
  dadu(3,2)  = 0
  dadu(3,3)  = -q

  ! time derivatives for runge-kutta
  dadt(1)  = - uu* gamma * (gp * dfdt(1) + go * dfdt(3)) 
  dadt(2)  = 0. 
  dadt(3)  = 0. 

  ! runge-kutta fourth order scheme

  do k=1,ndim
    kk2 (k) = kk1(k) + 0.5*dadt(k) * dt
    do i=1,ndim
     kk2(k) = kk2(k) + 0.5 * dadu(k,i) * kk1(i) * dt
    enddo
   enddo

  do k=1,ndim
    kk3 (k) = kk1(k) + 0.5*dadt(k) * dt
    do i=1,ndim
     kk3(k) = kk3(k) + 0.5*dadu(k,i) * kk2(i) * dt
    enddo
   enddo

  do k=1,ndim
    kk4 (k) = kk1(k) + dadt(k) * dt
    do i=1,ndim
     kk4(k) = kk4(k) +  dadu(k,i) * kk3(i) * dt
    enddo
  enddo


!  the system itself is integreted with runge-kutta
!  its tangent is integrated with euler.
!  this could be fixed in a nicer version

  do k=1,ndim
!   dy(m,k) = 1./6. * dt * ( kk1(k) + 2*kk2(k) + 2*kk3(k) +  kk4(k))
   dy(m,k) = 1./6. * dt * ( 6 * kk1(k) )
   dl(m,k) = 1. * dt * da(k)
  enddo  
enddo

end
      
