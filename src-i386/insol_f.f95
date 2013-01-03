! Fortran subroutine for reading precession and obliquity
! Precession and obliquity are computed after Berger (1978)
! Precession is reconstructed as a d'Alambert series. This
! is not quite the proceduce advised by Berger, but is the most
! efficient for the purpose of investigating the dynamical
! response of a simple model

! NAP : Number of terms in the d'Alember decomp. of precession : 
!       chose at least 20
! NAO : Same but for obliquity : at least 20

! USage : execute first 'READ_PRECESSION' and 'READ_OBLIQUITY'
! with the desired NAO and NAP. Theses routines output
! AMPlitudes, OMEgas (angular velocity) and phase (ANG)
! as arrays of dimension NAP and NAO, respectively

! SCALE TIME defaults to 10^4 : this means that time units
! are 10,000 years

! Then calles ASTRO, which returns a normalised version
! of PRECESSION, COPRECESSION and OBLIQUITY and their
! derivatives. 

! Compile in FORTRAN90, e.g. : gfortan -c -o insol.o insol.f 

! supplied : obliquity and precession must be in a 
!            a Data directory.


! With this procedure a good approximation of 
! Summer Solstice Insolation at 65 degrees North
! may be obtained as follows (normalised)

! 0.7639 * FORCING(1) + 0.4756 * FORCING(5)

! Author : M. Crucifix, Summer 2011. 


subroutine read_precession(nap,amppre,omepre, angpre, filename, ic)

  implicit none
  integer, intent(in) :: nap,ic
  character(len=ic), intent(in) :: filename
  integer i    ! counter
  double precision, intent(out) :: amppre(nap),omepre(nap), angpre(nap)
  double precision tmp1,tmp2,tmp3 ! tmp
  double precision, parameter :: scaletime = 1.d4
  double precision, parameter :: scaleprec = 1./0.01864561223205 
  double precision, parameter :: twopi = 2*3.1415926535897  
  double precision, parameter :: twopi360  = twopi / 360.
 
  open (10, file=filename, status='old')
  do i=1,nap
    read (10,*) tmp1,tmp2,tmp3
    amppre(i)=tmp1*scaleprec
    omepre(i)=tmp2*(twopi360/(60.*60.))*scaletime
    angpre(i)=tmp3*twopi360
  enddo
  close(10)
  end

subroutine read_obliquity(nao,ampobl,omeobl, angobl, filename, ic)

  implicit none
  integer, intent(in) ::  nao,ic
  character(len=ic), intent(in) :: filename
  double precision, intent(out) ::  ampobl(nao),omeobl(nao), angobl(nao)

  integer i    ! counter
  double precision :: tmp1,tmp2,tmp3 ! tmp
  double precision, parameter :: scaletime = 1.d4
  double precision, parameter ::  scaleobl  = 1./2462.2214466 
  double precision, parameter :: twopi = 2*3.1415926535897  
  double precision, parameter :: twopi360  = twopi / 360. 

  open (10, file=filename, status='old')
  do i=1,nao
    read (10,*) tmp1,tmp2,tmp3
    ampobl(i)=tmp1*scaleobl
    omeobl(i)=tmp2*(twopi360/(60.*60.))*scaletime
    angobl(i)=tmp3*twopi360
  enddo
  close(10)

end

subroutine astro(nap, nao, t,amppre, omepre, angpre, ampobl, & 
                 omeobl, angobl, forcing, dforcingdt)

   ! t : time in 10 ka  (i. e.: t=234 -> + 2340 ka ap) 
   ! nap, nao, ampre,  etc: as defined above 
   ! output: 
   ! forcing : precession, coprecession, obliquity


  implicit none
  integer, intent(in) ::  nap, nao
  double precision,intent(in)  :: t
  double precision, intent(in) :: amppre(nap), omepre(nap), angpre(nap)  
  double precision, intent(in) :: ampobl(nao), omeobl(nao), angobl(nao)  
  double precision, intent(out) :: forcing(3), dforcingdt(3)

   ! local variables
  double precision somep(nap),comep(nap)             ! tmp
  double precision someo(nao),comeo(nao)             ! tmp

  ! precession
  somep = dsin(omepre*t + angpre)
  comep = dcos(omepre*t + angpre)

  forcing(1)  = sum(amppre*somep)
  dforcingdt(1) = sum(omepre*(amppre*comep))

  forcing(2)  = sum(amppre*comep)
  dforcingdt(2)  = -sum(omepre*(amppre*somep))

  ! obliquity
  someo = dsin(omeobl*t + angobl)
  comeo = dcos(omeobl*t + angobl)

  forcing(3)  = sum(ampobl*comeo)
  dforcingdt(3) = - sum(omeobl*(ampobl*someo))

end  subroutine astro


