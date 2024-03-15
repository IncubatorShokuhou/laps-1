      subroutine get_raob_data_af(i4time_sys,ilaps_cycle_time,nx_l,ny_l
     1     ,i4time_raob_earliest,i4time_raob_latest,a9time_file
     1     ,filename       
     1     ,istatus)

!     ken dritz     28-jul-1997       added nx_l, ny_l to dummy argument list.
!     ken dritz     28-jul-1997       added call to get_r_missing_data.
!     ken dritz     28-jul-1997       changed laps_domain_file to 'nest7grid'.
!     ken dritz     28-jul-1997       removed include of lapsparms.for.
!     ken dritz     28-jul-1997       removed comment about "non-automatic
!                                     declarations" (above arrays dimensioned
!                                     by nx_l, ny_l); they are now automatic.

c     fortran template for file= 9614912000300o                          
      parameter (nvars=39) !number of variables
      parameter (nrec=   200)   !change this to generalize
c     variable ids run sequentially from 1 to nvars= 39
      integer rcode
      integer recdim
c     ****variables for this netcdf file****
c

      character*170 filename
      real lat_a(nx_l,ny_l)
      real lon_a(nx_l,ny_l)
      real topo_a(nx_l,ny_l)

      character*5 c5_staid
      character*11 a11_raob_reltime,a11_pibal_reltime
      character*9 a9time_raob, a9time_file
      character*8 c8_obstype

      call get_r_missing_data(r_missing_data,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'error getting r_missing_data'
          return
      endif

      c8_obstype = 'raob'

      nrecs = 0

      open(21,file=filename,status='old')

!     read initial header and ignore it
      do i = 1,6
          read(21,*,err=997,end=998)
      enddo ! i

 400  read(21,501,err=997,end=998)iwmostanum,c5_staid,n_good_levels       
     1       ,stalat,stalon,a11_raob_reltime,a11_pibal_reltime
     1       ,iheight,pressure,temp,dewpoint
     1       ,iwind_ind,idir,ispd
 501  format(i7,1x,a5,6x,i8
     1         ,4x,f6.2,3x,f7.2,2x,a11,2x,a11
     1         /3x,i9,3x,f9.0,2x,f10.0,2x,f10.0
     1         ,2x,i9,2x,i9,2x,i9)

      staelev = 0.

      write(6,*)' reltimes ',a11_raob_reltime,a11_pibal_reltime
 
      a9time_raob = a9time_file

      write(6,511,err=997)
     1             iwmostanum,n_good_levels,stalat
     1            ,stalon,staelev,c5_staid,a9time_raob,c8_obstype
      write(11,511,err=997)
     1             iwmostanum,n_good_levels,stalat
     1            ,stalon,staelev,c5_staid,a9time_raob,c8_obstype

  511 format(i12,i12,f11.4,f15.4,f15.0,1x,a5,3x,a9,1x,a8)


      if(.true.)then ! write out stuff for first level
              rheight = iheight
              if(abs(rheight) .gt. 1e6)then
                  rheight = r_missing_data
              endif

              if(abs(temp) .gt. 1e6)then
                  temp = r_missing_data
              else
                  temp = temp - 273.15
              endif

              if(abs(dewpoint) .gt. 1e6)then
                  dewpoint = r_missing_data
              elseif(temp .eq. r_missing_data)then
                  dewpoint = r_missing_data
              else
                  dewpoint = temp - dewpoint
              endif

              dir = idir
              spd = ispd

              if(abs(dir) .gt. 1e6 .or. abs(spd) .gt. 1e6)then
                  dir = r_missing_data
                  spd = r_missing_data
              endif

              write(6,*) rheight,pressure
     1              ,temp
     1              ,dewpoint
     1              ,iwind_ind,dir,spd,ilvl

              write(11,*) rheight,pressure
     1              ,temp
     1              ,dewpoint
     1              ,dir,spd

      endif


      if(n_good_levels .gt. 1)then
          do ilvl = 2,n_good_levels
              read(21,501,err=997,end=998)iwmostanum,c5_staid
     1       ,n_good_levels       
     1       ,stalat,stalon,a11_raob_reltime,a11_pibal_reltime
     1       ,iheight,pressure,temp,dewpoint
     1       ,iwind_ind,idir,ispd

              rheight = iheight
              if(abs(rheight) .gt. 1e6)then
                  rheight = r_missing_data
              endif

              if(abs(temp) .gt. 1e6)then
                  temp = r_missing_data
              else
                  temp = temp - 273.15
              endif

              if(abs(dewpoint) .gt. 1e6)then
                  dewpoint = r_missing_data
              elseif(temp .eq. r_missing_data)then
                  dewpoint = r_missing_data
              else
                  dewpoint = temp - dewpoint
              endif

              dir = idir
              spd = ispd

              if(abs(dir) .gt. 1e6 .or. abs(spd) .gt. 1e6)then
                  dir = r_missing_data
                  spd = r_missing_data
              endif

              write(6,*) rheight,pressure
     1              ,temp
     1              ,dewpoint
     1              ,iwind_ind,dir,spd,ilvl

              write(11,*) rheight,pressure
     1              ,temp
     1              ,dewpoint
     1              ,dir,spd

          enddo ! ilvl
      endif ! > 1 lvl

      nrecs = nrecs + 1

      write(6,*)' looping back to search for another sounding'

      go to 400

 997  write(6,*)' error reading raob file'

      go to 999

 998  write(6,*)' end of raob file'

 999  continue

      return
      end



