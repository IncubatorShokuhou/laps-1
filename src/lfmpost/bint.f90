 function bint(ri, rj, data, nx, ny)

    ! performs bilinear interpolation to point i,j from data array
    ! with dimensions nx,ny

    implicit none

    real :: ri, rj
    integer :: nx, ny
    real :: data(nx, ny)
    real :: bint
    real :: t, u
    integer :: i1, i2, j1, j2

    i1 = int(ri)
    i2 = int(ri + 1.)
    j1 = int(rj)
    j2 = int(rj + 1.)

    if (ri .ne. float(i1)) then
       t = (ri - float(i1))/float(i2 - i1)
    else
       t = 0.
       i2 = i1
    end if

    if (rj .ne. float(j1)) then
       u = (rj - float(j1))/float(j2 - j1)
    else
       u = 0.
       j2 = j1
    end if
    bint = (1.-t)*(1.-u)*data(i1, j1) + &
           t*(1.-u)*data(i2, j1) + &
           t*u*data(i2, j2) + &
           (1.-t)*u*data(i1, j2)
    return
 end function bint
