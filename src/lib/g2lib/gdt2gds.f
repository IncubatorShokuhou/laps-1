        subroutine gdt2gds(igds,igdstmpl,idefnum,ideflist,kgds,
     &                     igrid,iret)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    gdt2gds
c   prgmmr: gilbert        org: w/np11    date: 2003-06-17
c
c abstract: this routine converts grid information from a grib2
c   grid description section as well as its
c   grid definition template to grib1 gds info.  in addition,
c   a check is made to determine if the grid is an ncep
c   predefined grid.
c
c program history log:
c 2003-06-17  gilbert
c 2004-04-27  gilbert - added support for gaussian grids.
c
c usage:    call gdt2gds(igds,igdstmpl,idefnum,ideflist,kgds,igrid,iret)
c   input argument list:
c     igds()   - contains information read from the appropriate grib grid
c                definition section 3 for the field being returned.
c                must be dimensioned >= 5.
c                igds(1)=source of grid definition (see code table 3.0)
c                igds(2)=number of grid points in the defined grid.
c                igds(3)=number of octets needed for each
c                            additional grid points definition.
c                            used to define number of
c                            points in each row ( or column ) for
c                            non-regular grids.
c                            = 0, if using regular grid.
c                igds(4)=interpretation of list for optional points
c                            definition.  (code table 3.11)
c                igds(5)=grid definition template number (code table 3.1)
c     igdstmpl() - grid definition template values for gdt 3.igds(5)
c     idefnum    - the number of entries in array ideflist.  
c                  i.e. number of rows ( or columns )
c                  for which optional grid points are defined.
c     ideflist() - optional integer array containing
c                  the number of grid points contained in each row (or column).
c
c   output argument list:      (including work arrays)
c     kgds()   - grib1 gds as described in w3fi63 format.
c     igrid    - ncep predifined grib1 grid number
c                set to 255, if not ncep grid
c     iret     - error return value:
c                  0  = successful
c                  1  = unrecognized grib2 gdt number 3.igds(5)
c
c remarks: list caveats, other helpful hints or information
c
c attributes:
c   language: indicate extensions, compiler options
c   machine:  ibm sp
c
c$$$
! 
        integer,intent(in) :: idefnum
        integer,intent(in) :: igds(*),igdstmpl(*),ideflist(*)
        integer,intent(out) :: kgds(*),igrid,iret

        integer :: kgds72(200),kgds71(200),idum(200),jdum(200)

        iret=0
        if (igds(5).eq.0) then       !  lat/lon grid
           kgds(1)=0
           kgds(2)=igdstmpl(8)            ! ni
           kgds(3)=igdstmpl(9)            ! nj
           kgds(4)=igdstmpl(12)/1000      ! lat of 1st grid point
           kgds(5)=igdstmpl(13)/1000      ! long of 1st grid point
           kgds(6)=0                      ! resolution and component flags
           if (igdstmpl(1)==2 ) kgds(6)=64
           if ( btest(igdstmpl(14),4).or.btest(igdstmpl(14),5) ) 
     &         kgds(6)=kgds(6)+128
           if ( btest(igdstmpl(14),3) ) kgds(6)=kgds(6)+8
           kgds(7)=igdstmpl(15)/1000      ! lat of last grid point
           kgds(8)=igdstmpl(16)/1000      ! long of last grid point
           kgds(9)=igdstmpl(17)/1000      ! di
           kgds(10)=igdstmpl(18)/1000     ! dj
           kgds(11)=igdstmpl(19)          ! scanning mode
           kgds(12)=0
           kgds(13)=0
           kgds(14)=0
           kgds(15)=0
           kgds(16)=0
           kgds(17)=0
           kgds(18)=0
           kgds(19)=0
           kgds(20)=255
           kgds(21)=0
           kgds(22)=0
           !
           !  process irreg grid stuff, if necessary
           !
           if ( idefnum.ne.0 ) then
              if ( igdstmpl(8).eq.-1 ) then
                 kgds(2)=65535
                 kgds(9)=65535
              endif
              if ( igdstmpl(9).eq.-1 ) then
                 kgds(3)=65535
                 kgds(10)=65535
              endif
              kgds(19)=0
              kgds(20)=33
              if ( kgds(1).eq.1.or.kgds(1).eq.3 ) kgds(20)=43
              kgds(21)=igds(2)                   ! num of grid points
              do j=1,idefnum
                 kgds(21+j)=ideflist(j)
              enddo
           endif
        elseif (igds(5).eq.10) then       !  mercator grid
           kgds(1)=1                 ! grid definition template number
           kgds(2)=igdstmpl(8)            ! ni
           kgds(3)=igdstmpl(9)            ! nj
           kgds(4)=igdstmpl(10)/1000      ! lat of 1st grid point
           kgds(5)=igdstmpl(11)/1000      ! long of 1st grid point
           kgds(6)=0                      ! resolution and component flags
           if (igdstmpl(1)==2 ) kgds(6)=64
           if ( btest(igdstmpl(12),4).or.btest(igdstmpl(12),5) ) 
     &         kgds(6)=kgds(6)+128
           if ( btest(igdstmpl(12),3) ) kgds(6)=kgds(6)+8
           kgds(7)=igdstmpl(14)/1000      ! lat of last grid point
           kgds(8)=igdstmpl(15)/1000      ! long of last grid point
           kgds(9)=igdstmpl(13)/1000      ! lat intersects earth
           kgds(10)=0
           kgds(11)=igdstmpl(16)          ! scanning mode
           kgds(12)=igdstmpl(18)/1000     ! di
           kgds(13)=igdstmpl(19)/1000     ! dj
           kgds(14)=0
           kgds(15)=0
           kgds(16)=0
           kgds(17)=0
           kgds(18)=0
           kgds(19)=0
           kgds(20)=255
           kgds(21)=0
           kgds(22)=0
        elseif (igds(5).eq.30) then       ! lambert conformal grid
           kgds(1)=3
           kgds(2)=igdstmpl(8)            ! nx
           kgds(3)=igdstmpl(9)            ! ny
           kgds(4)=igdstmpl(10)/1000      ! lat of 1st grid point
           kgds(5)=igdstmpl(11)/1000      ! long of 1st grid point
           kgds(6)=0                      ! resolution and component flags
           if (igdstmpl(1)==2 ) kgds(6)=64
           if ( btest(igdstmpl(12),4).or.btest(igdstmpl(12),5) ) 
     &         kgds(6)=kgds(6)+128
           if ( btest(igdstmpl(12),3) ) kgds(6)=kgds(6)+8
           kgds(7)=igdstmpl(14)/1000      ! lon of orientation
           kgds(8)=igdstmpl(15)/1000      ! dx
           kgds(9)=igdstmpl(16)/1000      ! dy
           kgds(10)=igdstmpl(17)          ! projection center flag
           kgds(11)=igdstmpl(18)          ! scanning mode
           kgds(12)=igdstmpl(19)/1000     ! lat in 1
           kgds(13)=igdstmpl(20)/1000     ! lat in 2
           kgds(14)=igdstmpl(21)/1000     ! lat of s. pole of projection
           kgds(15)=igdstmpl(22)/1000     ! lon of s. pole of projection
           kgds(16)=0
           kgds(17)=0
           kgds(18)=0
           kgds(19)=0
           kgds(20)=255
           kgds(21)=0
           kgds(22)=0
        elseif (igds(5).eq.40) then       !  gaussian lat/lon grid
           kgds(1)=4
           kgds(2)=igdstmpl(8)            ! ni
           kgds(3)=igdstmpl(9)            ! nj
           kgds(4)=igdstmpl(12)/1000      ! lat of 1st grid point
           kgds(5)=igdstmpl(13)/1000      ! long of 1st grid point
           kgds(6)=0                      ! resolution and component flags
           if (igdstmpl(1)==2 ) kgds(6)=64
           if ( btest(igdstmpl(14),4).or.btest(igdstmpl(14),5) ) 
     &         kgds(6)=kgds(6)+128
           if ( btest(igdstmpl(14),3) ) kgds(6)=kgds(6)+8
           kgds(7)=igdstmpl(15)/1000      ! lat of last grid point
           kgds(8)=igdstmpl(16)/1000      ! long of last grid point
           kgds(9)=igdstmpl(17)/1000      ! di
           kgds(10)=igdstmpl(18)          ! n - number of parallels
           kgds(11)=igdstmpl(19)          ! scanning mode
           kgds(12)=0
           kgds(13)=0
           kgds(14)=0
           kgds(15)=0
           kgds(16)=0
           kgds(17)=0
           kgds(18)=0
           kgds(19)=0
           kgds(20)=255
           kgds(21)=0
           kgds(22)=0
        elseif (igds(5).eq.20) then       ! polar stereographic grid
           kgds(1)=5
           kgds(2)=igdstmpl(8)            ! nx
           kgds(3)=igdstmpl(9)            ! ny
           kgds(4)=igdstmpl(10)/1000      ! lat of 1st grid point
           kgds(5)=igdstmpl(11)/1000      ! long of 1st grid point
           kgds(6)=0                      ! resolution and component flags
           if (igdstmpl(1)==2 ) kgds(6)=64
           if ( btest(igdstmpl(12),4).or.btest(igdstmpl(12),5) ) 
     &         kgds(6)=kgds(6)+128
           if ( btest(igdstmpl(12),3) ) kgds(6)=kgds(6)+8
           kgds(7)=igdstmpl(14)/1000      ! lon of orientation
           kgds(8)=igdstmpl(15)/1000      ! dx
           kgds(9)=igdstmpl(16)/1000      ! dy
           kgds(10)=igdstmpl(17)          ! projection center flag
           kgds(11)=igdstmpl(18)          ! scanning mode
           kgds(12)=0
           kgds(13)=0
           kgds(14)=0
           kgds(15)=0
           kgds(16)=0
           kgds(17)=0
           kgds(18)=0
           kgds(19)=0
           kgds(20)=255
           kgds(21)=0
           kgds(22)=0
        else
           print *,'gdt2gds: unrecognized grib2 gdt = 3.',igds(5)
           iret=1
           kgds(1:22)=0
           return
        endif
!
!   can we determine ncep grid number ?
!
        igrid=255
        do j=254,1,-1
        !do j=225,225
           kgds71=0
           kgds72=0
           call w3fi71(j,kgds71,ierr)
           if ( ierr.ne.0 ) cycle
           ! convert w to e for longitudes
           if ( kgds71(3).eq.0 ) then    ! lat/lon
              if ( kgds71(7).lt.0 ) kgds71(7)=360000+kgds71(7)
              if ( kgds71(10).lt.0 ) kgds71(10)=360000+kgds71(10)
           elseif ( kgds71(3).eq.1 ) then    ! mercator
              if ( kgds71(7).lt.0 ) kgds71(7)=360000+kgds71(7)
              if ( kgds71(10).lt.0 ) kgds71(10)=360000+kgds71(10)
           elseif ( kgds71(3).eq.3 ) then     ! lambert conformal
              if ( kgds71(7).lt.0 ) kgds71(7)=360000+kgds71(7)
              if ( kgds71(9).lt.0 ) kgds71(9)=360000+kgds71(9)
              if ( kgds71(18).lt.0 ) kgds71(18)=360000+kgds71(18)
           elseif ( kgds71(3).eq.4 ) then     ! guassian lat/lon
              if ( kgds71(7).lt.0 ) kgds71(7)=360000+kgds71(7)
              if ( kgds71(10).lt.0 ) kgds71(10)=360000+kgds71(10)
           elseif ( kgds71(3).eq.5 ) then     ! polar stereographic
              if ( kgds71(7).lt.0 ) kgds71(7)=360000+kgds71(7)
              if ( kgds71(9).lt.0 ) kgds71(9)=360000+kgds71(9)
           endif
           call r63w72(idum,kgds,jdum,kgds72)
           if ( kgds72(3).eq.3 ) kgds72(14)=0    ! lambert conformal fix
           if ( kgds72(3).eq.1 ) kgds72(15:18)=0    ! mercator fix
           if ( kgds72(3).eq.5 ) kgds72(14:18)=0    ! polar str fix
           !print *,'sagt71:',(kgds71(k),k=1,30)
           !print *,'sagt72:',(kgds72(k),k=1,30)
           if ( all(kgds71.eq.kgds72) ) then
              igrid=j
              return
           endif
        enddo

        return
        end
