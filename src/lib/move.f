cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis
c
c
        subroutine move(a,b,imax,jmax)
c
c.....  routine to move array 'a' into array 'b'.
c
        real a(imax,jmax), b(imax,jmax)
c
        do j=1,jmax
        do i=1,imax
          b(i,j) = a(i,j)
        enddo !i
        enddo !j
c
        return
        end
c
c-----------------------------------------------------
c
        subroutine move_i(a,b,imax,jmax)
c
c.....  routine to move array 'a' into array 'b'.
c
        integer a(imax,jmax), b(imax,jmax)
c
        do j=1,jmax
        do i=1,imax
          b(i,j) = a(i,j)
        enddo !i
        enddo !j
c
        return
        end

c
c-----------------------------------------------------
c
        subroutine move_3d(a,b,imax,jmax,kmax)
c
c.....  routine to move (copy) array 'a' into array 'b'.
c
        real a(imax,jmax,kmax), b(imax,jmax,kmax)
c
        do k=1,kmax
        do j=1,jmax
        do i=1,imax
          b(i,j,k) = a(i,j,k)
        enddo !i
        enddo !j
        enddo !k
c
        return
        end
c
c
        subroutine move_2dto3d(a,b,index,imax,jmax,kmax)
c
c.....  routine to move (copy) the 2d array 'a' into one level
c.....  of the 3d array 'b'.  the level is defined by 'index'.
c
c       original:  p. stamus  noaa/fsl  15 apr 1997
c
        real a(imax,jmax), b(imax,jmax,kmax)
c
        do j=1,jmax
        do i=1,imax
          b(i,j,index) = a(i,j)
        enddo !i
        enddo !j
c
        return
        end
c
c
        subroutine move_3dto2d(a,index,b,imax,jmax,kmax)
c
c.....  routine to move (copy) one level of the 3d array 'a' into 
c.....  the 2d array 'b'.  the level is defined by 'index'.
c
c       original:  p. stamus  noaa/fsl  15 apr 1997
c
        real a(imax,jmax,kmax), b(imax,jmax)
c
        do j=1,jmax
        do i=1,imax
          b(i,j) = a(i,j,index)
        enddo !i
        enddo !j
c
        return
        end
