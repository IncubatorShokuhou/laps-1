module grid_utils
  
! this module contains utilities related to grid manipulation
! (e.g., creating data for the b and c grids on an arakawa-c
!  type stagger when only the a grid is available)


contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine arakawa_c_n2t(datain, nx, ny, nz, dataout)
    
    ! staggers a 3d array of data from the non-staggered points
    ! to the mass grid of an arakawa c stagger.

    implicit none
    integer, intent(in)                :: nx
    integer, intent(in)                :: ny
    integer, intent(in)                :: nz
    real, intent(in)                   :: datain(nx,ny,nz)
    real, intent(out)                  :: dataout(nx,ny,nz)
   
    integer                            :: i,j,k
    print *, 'staggering to t grid (arakawa c)'
    do k = 1, nz
      do j = 1, ny-1
        do i = 1, nx-1
          dataout(i,j,k) = 0.25 * ( datain(i,j,k)    + &
                                    datain(i+1,j,k)  + &
                                    datain(i+1,j+1,k)+ &
                                    datain(i,j+1,k) )
        enddo
        ! fill unused rightmost column
        dataout(nx,j,k) = dataout(nx-1,j,k)
      enddo
      ! fill unused uppermost row
      dataout(:,ny,k) = dataout(:,ny-1,k)
    enddo
    return
  end subroutine arakawa_c_n2t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine arakawa_c_n2u(datain, nx, ny, nz, dataout)

    ! staggers a 3d array of data from the non-staggered points
    ! to the u grid of an arakawa c stagger.

    implicit none
    integer, intent(in)                :: nx
    integer, intent(in)                :: ny
    integer, intent(in)                :: nz
    real, intent(in)                   :: datain(nx,ny,nz)
    real, intent(out)                  :: dataout(nx,ny,nz)

    integer                            :: i,j,k
    print *, 'staggering to u grid (arakawa c)'
    do k = 1, nz
      do j = 1, ny-1
        do i = 1, nx
          dataout(i,j,k) = 0.50 * ( datain(i,j,k)    + &
                                    datain(i,j+1,k) )
        enddo
      enddo
      ! fill unused uppermost row
      dataout(:,ny,k) = dataout(:,ny-1,k)
    enddo
    return
  end subroutine arakawa_c_n2u                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine arakawa_c_n2v(datain, nx, ny, nz, dataout)

    ! staggers a 3d array of data from the non-staggered points
    ! to the v grid of an arakawa c stagger.

    implicit none
    integer, intent(in)                :: nx
    integer, intent(in)                :: ny
    integer, intent(in)                :: nz
    real, intent(in)                   :: datain(nx,ny,nz)
    real, intent(out)                  :: dataout(nx,ny,nz)

    integer                            :: i,j,k
    print *, 'staggering to v grid (arakawa c)'
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx-1
          dataout(i,j,k) = 0.50 * ( datain(i,j,k)    + &
                                    datain(i+1,j,k) )
        enddo
        ! fill unused right column
        dataout(nx,j,k) = dataout(nx-1,j,k)
      enddo
    enddo
    return
  end subroutine arakawa_c_n2v            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine arakawa_c_t2n(datain,nx_t,ny_t,nz,dataout)

    ! destaggers a wrf arakawa c from the staggered thermodynamic
    ! points to the non-staggered points.   note that this routine
    ! returns an array that is one element larger in each direction
    ! than the input array.
    !
    !  example:  input 3x3 "t" points, return 4x4 "n" points
    !
    !           n   n   n   n
    !             t   t   t
    !           n   n   n   n 
    !             t   t   t
    !           n   n   n   n
    !             t   t   t 
    !           n   n   n   n

    implicit none

    integer, intent(in)       :: nx_t
    integer, intent(in)       :: ny_t    
    integer, intent(in)       :: nz
    real, intent(in)          :: datain(nx_t,ny_t,nz)       
    real, intent(out)         :: dataout(nx_t+1,ny_t+1,nz)   

    integer                   :: i,j,k

    vertical_loop:  do k=1,nz

      ! first, compute all of the interior points

      do j = 2, ny_t
        do i = 1, nx_t

           dataout(i,j,k) = 0.25*( datain(i-1,j-1,k) + datain(i-1,j,k) + &
                                  datain(i,j,k) + datain(i,j-1,k) )

        enddo
      enddo

      ! now, extrapolate upper and lower rows, except corner points

      do i = 2, nx_t
 
        dataout(i,1,k) = 2.0* dataout(i,2,k)-dataout(i,3,k)
        dataout(i,ny_t+1,k) = 2.0*dataout(i,ny_t,k)-dataout(i,ny_t-1,k)

      enddo

      ! extrapolate left and right columns, except corner points

      do j = 2, ny_t
    
        dataout(1,j,k) = 2.0*dataout(2,j,k)-dataout(3,j,k)
        dataout(nx_t+1,j,k) = 2.0*dataout(nx_t,j,k)-dataout(nx_t-1,j,k)

      enddo

      ! compute corner point values by solving for 4 point average

      dataout(1,1,k) = 4.0 * datain(1,1,k) - &
                             dataout(1,2,k) - &
                             dataout(2,2,k) - &
                             dataout(2,1,k)
      dataout(1,ny_t+1,k) = 4.0 * datain(1,ny_t,k) - &
                                  dataout(1,ny_t,k) - &
                                  dataout(2,ny_t+1,k) - &
                                  dataout(2,ny_t,k)
      dataout(nx_t+1,ny_t+1,k) = 4.0 * datain(nx_t,ny_t,k) - &
                                       dataout(nx_t,ny_t,k) - &
                                       dataout(nx_t,ny_t+1,k) - &
                                       dataout(nx_t+1,ny_t,k)
      dataout(nx_t+1,1,k) = 4.0 * datain(nx_t,1,k) - &
                                dataout(nx_t,1,k) - &
                                dataout(nx_t,2,k) - &
                                dataout(nx_t+1,2,k)

    enddo vertical_loop       
    return
  end subroutine arakawa_c_t2n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine arakawa_c_u2n(datain,nx_u,ny_u,nz,dataout)

    ! destaggers from the wrf "u" grid to the non-staggered grid.  the
    ! return array will be 1 element larger in the y direction than
    ! the input array.  example:
    !
    !   input 4x3 "u" grid, return 4x4 "n" grid:
    !
    !            n   n   n   n
    !            u   u   u   u
    !            n   n   n   n
    !            u   u   u   u
    !            n   n   n   n
    !            u   u   u   u
    !            n   n   n   n

    implicit none

    integer,  intent(in)   :: nx_u
    integer,  intent(in)   :: ny_u
    integer,  intent(in)   :: nz
    real,     intent(in)   :: datain(nx_u,ny_u,nz)
    real,     intent(out)  :: dataout(nx_u,ny_u+1,nz)

    integer                :: i,k
    integer                :: nx,ny

    nx = nx_u
    ny = ny_u+1

    do k = 1, nz
   
      ! linear interpolation along each column, except top/bottom rows
      do i = 1, nx
   
        ! average of points above and below to fill interior rows 
        dataout(i,2:ny_u,k) = 0.5*(datain(i,1:ny_u-1,k)+datain(i,2:ny_u,k))
 
        ! fill bottom row
        dataout(i,1,k) = 2.0 * dataout(i,2,k) - dataout(i,3,k)

        ! fill top row
        dataout(i,ny,k) = 2.0 * dataout(i,ny-1,k) - dataout(i,ny-2,k)

      enddo
    enddo 
    return
  end subroutine arakawa_c_u2n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine arakawa_c_v2n(datain,nx_v,ny_v,nz,dataout)

    ! destaggers from the wrf "v" grid to the non-staggered grid.  the
    ! return array will be 1 element larger in the x direction than
    ! the input array.  example:
    !
    !   input 3x4 "v" grid, return 4x4 "n" grid:
    !
    !            n v n v n v n
    !           
    !            n v n v n v n
    !       
    !            n v n v n v n
    !        
    !            n v n v n v n

    implicit none
    integer,  intent(in)   :: nx_v
    integer,  intent(in)   :: ny_v
    integer,  intent(in)   :: nz
    real,     intent(in)   :: datain(nx_v,ny_v,nz)
    real,     intent(out)  :: dataout(nx_v+1,ny_v,nz)

    integer                :: j,k
    integer                :: nx,ny

    nx = nx_v + 1
    ny = ny_v

    do k = 1, nz

      ! linear interpolation along each row, except left/right columns
      do j = 1, ny

        ! average of points above and below to fill interior rows 
        dataout(2:nx_v,j,k) = 0.5*(datain(1:nx_v-1,j,k)+datain(2:nx_v,j,k))

        ! fill left column
        dataout(1,j,k) = 2.0 * dataout(2,j,k) - dataout(3,j,k)

        ! fill top row
        dataout(nx,j,k) = 2.0 * dataout(nx-1,j,k) - dataout(nx-2,j,k)

      enddo
    enddo
    return
  end subroutine arakawa_c_v2n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine arakawa_c_v2t(datain,nx_v,ny_v,nz,dataout)

    ! destaggers from the wrf "v" grid to the mass grid.  the
    ! return array will be 1 element smaller in the y direction than
    ! the input array.  example:
    !
    !   input 3x4 "v" grid, return 3x3 "t" grid:
    !
    !              v   v   v  
    !              t   t   t
    !              v   v   v  
    !              t   t   t
    !              v   v   v  
    !              t   t   t
    !              v   v   v  

    implicit none
    integer,  intent(in)   :: nx_v
    integer,  intent(in)   :: ny_v
    integer,  intent(in)   :: nz
    real,     intent(in)   :: datain(nx_v,ny_v,nz)
    real,     intent(out)  :: dataout(nx_v,ny_v-1,nz)

    integer                :: j,k
    integer                :: nx,ny

    nx = nx_v 
    ny = ny_v - 1

    do k = 1, nz

      ! linear interpolation along each column
      do j = 1, ny

        dataout(:,j,k) = 0.5*(datain(:,j,k) + datain(:,j+1,k))

      enddo
    enddo
    return
  end subroutine arakawa_c_v2t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine arakawa_c_u2t(datain,nx_u,ny_u,nz,dataout)

    ! destaggers from the wrf "u" grid to the mass grid.  the
    ! return array will be 1 element smaller in the x direction than
    ! the input array.  example:
    !
    !   input 4x3 "u" grid, return 3x3 "t" grid:
    !
    !                           
    !            u t u t u t u
    !                             
    !            u t u t u t u
    !                             
    !            u t u t u t u
    !                           
    implicit none
    integer,  intent(in)   :: nx_u
    integer,  intent(in)   :: ny_u
    integer,  intent(in)   :: nz
    real,     intent(in)   :: datain(nx_u,ny_u,nz)
    real,     intent(out)  :: dataout(nx_u-1,ny_u,nz)

    integer                :: i,k
    integer                :: nx,ny

    nx = nx_u - 1
    ny = ny_u

    do k = 1, nz

      ! linear interpolation along each row     
      do i = 1, nx

        dataout(i,:,k) = 0.5*(datain(i,:,k) + datain(i+1,:,k))

      enddo
    enddo
    return
  end subroutine arakawa_c_u2t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine wlevs2hlevs(datain,nx,ny,nz_in,dataout)

    ! vertically destaggers an array from the full w levels to
    ! the half levels.  output array is one less in the z dimension
    ! than the input array

    implicit none
    integer, intent(in)   :: nx,ny,nz_in
    real, intent(in)      :: datain(nx,ny,nz_in)
    real, intent(out)     :: dataout(nx,ny,nz_in-1)

    integer :: k

    do k = 1,nz_in-1

      dataout(:,:,k) = 0.5*(datain(:,:,k)+datain(:,:,k+1))
  
    enddo
    return
  end subroutine wlevs2hlevs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module grid_utils
