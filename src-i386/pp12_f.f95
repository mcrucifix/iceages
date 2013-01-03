!  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!  PARRENIN AND PAILLARD 2012
!  Terminations VI and VII (-530 and -720 kyrBP)
!  tell us the importance of obliquity and precession in
!  the triggering of deglaciations
!  Climate of the Past Discussions
!  doi: 10.5194/cpd-8-3143-2012
!  The time unit is 10 ka
! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function  f ( x )
  double precision f,x
  double precision, parameter  :: a    = 0.68034
  if (x > 0) then 
     f = x
  else
     f = x + sqrt ( 4 * a * a + x * x ) - 2 * a
  endif
  f =  f 
end function

subroutine pp12_f(state,n,ndim,npar,par,forcing, dforcingdt, deltat,dy,ds,dl)

  implicit none
  integer k,j,l,m
  integer, intent(in)  :: n, ndim, npar
  double precision, intent(in) :: par(n,npar),state(n,ndim) 
  double precision, intent(in) :: deltat
  double precision, intent(in) :: forcing(3), dforcingdt(3)
  double precision, intent(inout) :: ds(n,ndim),dl(n,ndim)
  double precision, intent(out)   :: dy(n,2)


  ! constants to multiply precession and obliquity
  ! to obtain normalised values

  double precision, parameter :: gp = 0.8536269 
  double precision, parameter :: go = 1.094381 

  !     parameters (time unit = 10 ka ! ) 
  double precision, parameter ::  aesi = 14.561
  double precision, parameter ::  aeco = 3.87
  double precision, parameter ::  ao   = 11.3662
  double precision, parameter ::  ag   = 9.7825
  double precision, parameter ::  ad   = -7.469
  double precision, parameter ::  tdm1 = 1.199  ! 1. / td
  double precision, parameter ::  kesi = 14.6348
  double precision, parameter ::  keco = 2.28061
  double precision, parameter ::  ko   = 18.5162
  double precision, parameter ::  v0   = 122.918
  double precision, parameter ::  v1   = 3.10301

  !     local variables

  double precision f

  double precision localstate(n,ndim)
  double precision v(n), s (n)

  double precision cfobl
  double precision gamma(n), omega(n), dvdt(n), ggp(n)
  double precision threshold (n), ifpre(n), ifcpre(n), ifobl(n)

  localstate = state

  v = state(:,1)
  s = state(:,2)
   
  !     non-linear transformation of the precession forcing

  cfobl = forcing(3) - 0.0298181

  gamma  = par(:,1)
  omega  = par(:,2)

  ggp = gamma * gp

  ifpre =  ggp * ( f(forcing(1))  - 0.14308 ) / 1.22139
  ifcpre = ggp * ( f(forcing(2)) - 0.14308 ) / 1.22139
  ifobl  = gamma * go * ( f(cfobl) - 0.1298 ) / 1.0786

  !     s can only be 0 (glacial) or 1 (interglacial)

  threshold = kesi * ggp * forcing(1) + keco * ggp * forcing(2) & 
              + ko * go * gamma * cfobl

  where (s < .5 .and. (threshold + v >  v0) )  s = 1.
  where (s > .5 .and. (threshold  < v1) )  s = 0.

  dvdt = - aesi * ifpre - aeco * ifcpre - ao * ifobl
  where (s >  .5)
    dvdt = dvdt + ad - v * tdm1
  elsewhere
    dvdt = dvdt + ag
  endwhere

  localstate(:,1) = v + deltat / omega * dvdt

  where (localstate (:,1) < 0) localstate(:,1) = 0

  localstate(:,2) = s

  dy = localstate - state

eND
