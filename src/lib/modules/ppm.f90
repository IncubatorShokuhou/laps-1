! 
! write portable pixmap image file (ppm, pbm, pgm) 
! with fortran 
!
!   charles o'neill oct 23, 2009
!    charles.oneill@gmail.com
!    www.caselab.okstate.edu
!
! copyright (c) 2009 charles o'neill
!
! permission is hereby granted, free of charge, to any person
! obtaining a copy of this software and associated documentation
! files (the "software"), to deal in the software without
! restriction, including without limitation the rights to use,
! copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the software, and to permit persons to whom the
! software is furnished to do so, subject to the following
! conditions:
!
! the above copyright notice and this permission notice shall be
! included in all copies or substantial portions of the software.
!
! the software is provided "as is", without warranty of any kind,
! express or implied, including but not limited to the warranties
! of merchantability, fitness for a particular purpose and
! noninfringement. in no event shall the authors or copyright
! holders be liable for any claim, damages or other liability,
! whether in an action of contract, tort or otherwise, arising
! from, out of or in connection with the software or the use or
! other dealings in the software.
!

module ppm
  implicit none

contains

  !--------------------------------------------------------------
  ! portable pixmap type 1 (black and white)
  subroutine writeppm1matrix(m,text)
    integer :: m(:,:)
    character(len=*) :: text
    integer :: cols,rows
    integer :: i,j
 
    ! open file   
    open(unit=100, file=trim(text)//".pbm", status='unknown')
    
    ! write header and ppm file type
    write(100,'( a )') "p1"
    write(100,'( a )') "# ppm type 1 file (generated with fortran)"
 
    ! write image size
    cols = size(m,2)
    rows = size(m,1)
    write(100,'( i6, 1x, i6 )') cols, rows
    
    ! write image
    do i=1,rows
      do j=1,cols
        write(100,'( i1 )', advance='no') m(i,j)
      enddo
      write(100,*) ! endline
    enddo
  end subroutine

  !--------------------------------------------------------------
  ! portable pixmap type 2 (grayscale)
  subroutine writeppm2matrix(m,text)
    integer :: m(:,:)
    character(len=*) :: text
    integer :: cols,rows
    integer :: i,j
    integer :: maxvalue

    ! open file   
    open(unit=100, file=trim(text)//".pgm", status='unknown')
    
    ! write header and ppm file type
    write(100,'( a )') "p2"
    write(100,'( a )') "# ppm type 2 file (generated with fortran)"
 
    ! write image size
    cols = size(m,2)
    rows = size(m,1)
    write(100,'( i6, 1x, i6 )') cols, rows
    
    ! write maximum value
    maxvalue = maxval(maxval(m,dim=1),dim=1)
    write(100,'( i6 )') maxvalue
    
    ! write image
    do i=1,rows
      do j=1,cols
        write(100,'( i5,1x )', advance='no') m(i,j)
      enddo
      write(100,*) ! endline
    enddo
  end subroutine

  !--------------------------------------------------------------
  ! portable pixmap type 3 (rgb color)
  subroutine writeppm3matrix(r,g,b,text)
    integer :: r(:,:),g(:,:),b(:,:)
    character(len=*) :: text
    integer :: cols,rows
    integer :: i,j
    integer :: maxvalue,maxr,maxg,maxb
 
    ! open file   
    open(unit=100, file=trim(text)//".ppm", status='unknown')
    
    ! write header and ppm file type
    write(100,'( a )') "p3"
    write(100,'( a )') "# ppm type 2 file (generated with fortran)"
 
    ! write image size
    cols = size(r,2)
    rows = size(r,1)
    write(100,'( i6, 1x, i6 )') cols, rows
    
    ! write maximum value
    maxr = maxval(maxval(r,dim=1),dim=1)
    maxg = maxval(maxval(g,dim=1),dim=1)
    maxb = maxval(maxval(b,dim=1),dim=1)
    maxvalue = max(maxr,maxg,maxb)
    write(6,*)' image maxvalue is ',maxr,maxg,maxb,maxvalue

!   set brightness scaling in the image at a floor value
!   image will always have a max >= 85 counts
!   if(maxvalue .lt. 64)then
!       maxvalue = 4 * maxvalue
    if(maxvalue .lt. 85)then
        maxvalue = 3 * maxvalue
    else
        maxvalue = 255
    endif

    write(6,*)' image scaled to   ',maxvalue
    write(100,'( i6 )') maxvalue
    
    ! write image
    do i=1,rows
      do j=1,cols
        write(100,'( 3(i5,1x) )') r(i,j),g(i,j),b(i,j)
      enddo
    enddo

    close(100)
  end subroutine

  !--------------------------------------------------------------
  ! portable pixmap type 3 (rgb color)
  subroutine writeppm3_16bit(r,g,b,text)
    integer :: r(:,:),g(:,:),b(:,:)
    character(len=*) :: text
    integer :: cols,rows
    integer :: i,j
    integer :: maxvalue
 
    ! open file   
    open(unit=100, file=trim(text)//".ppm", status='unknown')
    
    ! write header and ppm file type
    write(100,'( a )') "p3"
    write(100,'( a )') "# ppm type 2 file (generated with fortran)"
 
    ! write image size
    cols = size(r,2)
    rows = size(r,1)
    write(100,'( i6, 1x, i6 )') cols, rows
    
    ! write maximum value
    maxvalue = max( maxval(maxval(r,dim=1),dim=1)&
                   ,maxval(maxval(g,dim=1),dim=1)&
                   ,maxval(maxval(b,dim=1),dim=1))
    write(6,*)' image maxvalue is ',maxvalue
    maxvalue = 65535
    write(6,*)' image scaled to   ',maxvalue
    write(100,'( i6 )') maxvalue
    
    ! write image
    do i=1,rows
      do j=1,cols
        write(100,'( 3(i5,1x) )') r(i,j),g(i,j),b(i,j)
      enddo
    enddo
  end subroutine

  !--------------------------------------------------------------
  ! test module
  subroutine testppm
    integer,parameter :: n = 100
    integer :: a(n,n)
    integer :: r(n,n)
    integer :: g(n,n)
    integer :: b(n,n)
    real :: aa(n,n)
    integer :: i,j

    ! show the pixmap format with a simple case
    open(unit=100, file="test.ppm", status='unknown')
    write(100,'( a )') "p1"
    write(100,'( a )') "# this is an example bitmap"
    write(100,'( a )') "18 10"
    write(100,'( a )') "0 0 0 0 0 0   0 0 0 0 0 0   0 0 0 0 0 0"
    write(100,'( a )') "0 0 1 1 0 0   0 0 1 1 1 0   0 1 0 0 1 0"
    write(100,'( a )') "0 1 0 0 1 0   0 1 0 0 1 0   0 1 0 0 1 0"
    write(100,'( a )') "0 1 0 0 1 0   0 1 0 0 0 0   0 1 0 0 1 0"
    write(100,'( a )') "0 1 0 0 1 0   0 0 1 0 0 0   0 1 0 0 1 0"
    write(100,'( a )') "0 1 0 0 1 0   0 0 0 1 0 0   0 1 0 0 1 0"
    write(100,'( a )') "0 1 0 0 1 0   0 0 0 0 1 0   0 1 0 0 1 0"
    write(100,'( a )') "0 1 0 0 1 0   0 1 0 0 1 0   0 1 0 0 1 0"
    write(100,'( a )') "0 0 1 1 0 0   0 1 1 1 0 0   0 0 1 1 0 0"
    write(100,'( a )') "0 0 0 0 0 0   0 0 0 0 0 0   0 0 0 0 0 0"
    close(100)
        
    ! get a matrix of random numbers
    call random_seed()
    call random_number(aa)
    
    ! setup and write type 1
    a = aa + 0.5
    call writeppm1matrix(a,"subtest")
    write(*,*) "writeppm1matrix"
  
    ! setup and write type 2
    a = aa * 256
    call writeppm2matrix(a,"subtest2")
    write(*,*) "writeppm2matrix"

    ! setup and write type 3 
    do i=1,n
      do j=1,n
        r(i,j) = abs( mod(i+j,n/5) )
        g(i,j) = abs( mod(i-j,n/5) )
        b(i,j) = abs( mod(i*j,n/5) )
      enddo
    enddo
    call writeppm3matrix(r,g,b,"subtest3")
    write(*,*) "writeppm3matrix"
  end subroutine

  subroutine read_ppm(u,img,ncol,iwidth,iheight)

    integer ncol,iwidth,iheight
    integer ncol2,iwidth2,iheight2

    integer img(ncol,iwidth,iheight)
    integer u
    character(2) :: sign

    read(u, '(a2)') sign
    read(u, *) iwidth2,iheight2
    read(u, *) ncol2

    write(6,*) sign
    write(6,*) iwidth2,iheight2
    write(6,*) ncol2

    if (ncol2 .ne. 255 .and. ncol2 .ne. 65535)then
        write(6,*)' error in read_ppm - ncol2 appears incorrect: ',ncol2
        return
    endif

    read(u,*)img

    write(6,*)' image has been read in subroutine read_ppm'

  end subroutine

end module
