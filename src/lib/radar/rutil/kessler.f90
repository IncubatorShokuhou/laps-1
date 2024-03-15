
subroutine kessler_mr2z(nx,ny,nz &
                       ,rho &                                ! kg/m^3
                       ,rainmr,icemr,snowmr,graupelmr  &
                       ,refl)

! subroutine to compute estimated radar reflectivity (z) from
! the precipitation mixing ratios.  the estimation
! is done using formulae from kessler (1969) and 
! rogers and yau (1989).  

! adapted from usaf weather agency routine.  
! brent shaw, noaa forecast system lab, dec 2000

! put into subroutine by steve albers
  
implicit none

integer :: nx,ny,nz,i,j,k
real, parameter :: svnfrth=7.0/4.0
real, dimension(nx,ny,nz) :: rho,rainmr,icemr,snowmr,graupelmr,refl

refl=0.0

do j=1,ny
do i=1,nx
   do k=1,nz

!     compute the basic reflectivity 
      refl(i,j,k) =17300.0 * &
                  (rho(i,j,k) * 1000.0 * &
                   max(0.0,rainmr(i,j,k)))**svnfrth

!     add the ice component
      refl(i,j,k)=refl(i,j,k) + &
                  38000.0*(rho(i,j,k) * 1000.0 * &
                  max(0.0,icemr(i,j,k)+snowmr(i,j,k)+graupelmr(i,j,k)))**2.2

   enddo

enddo
enddo

return
end

