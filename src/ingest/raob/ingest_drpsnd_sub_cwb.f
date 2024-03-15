      subroutine get_drpsnd_data_cwb ( i4time_sys, ilaps_cycle_time,
     ~         nx_l, ny_l, i4time_drpsnd_earliest, i4time_drpsnd_latest,
     ~         a9_time, filename,lun_out, istatus )

      integer   loopnum, levelnum  
      parameter ( loopnum=100, levelnum=90 )

      character(*)   filename
      character(3)   reportflag
      character(5)   taskname(loopnum)
      character(2)   yy, mo, dd, hh, mn, flag
      character(9)   a9time(loopnum), a9timedummy, a10_to_a9, a9_time
      character(10)  time
  
      real latitude_out(loopnum,levelnum)
      real longitude_out(loopnum,levelnum)
      character*9 a9time_out(loopnum,levelnum)
      character c8_obstype(loopnum)*8
      character c5_staid(loopnum)*5

      real  lat_a(nx_l,ny_l), lon_a(nx_l,ny_l), topo_a(nx_l,ny_l)
c wen modify 
c      real  latitudedummy, longitudedummy
       real  lat1, lon1
       integer  latitudedummy, longitudedummy
c
      real  elevation(loopnum), latitude(loopnum), longitude(loopnum)
      real  pressure(loopnum,levelnum), height(loopnum,levelnum)
      real  temperature(loopnum,levelnum)
      real  tempdewdiff(loopnum,levelnum), dewpoint(loopnum,levelnum)
      real  winddir(loopnum,levelnum), windspeed(loopnum,levelnum)
      real  prshtdbtl(loopnum), prshtdbth(loopnum)
      real  prstmdbtl(loopnum), prstmdbth(loopnum)
      real  staelev(loopnum)

      integer wmoid(loopnum), layernum(loopnum)
      integer heightqua(loopnum,levelnum), dewpointqua(loopnum,levelnum)       
      integer temperaturequa(loopnum,levelnum),windqua(loopnum,levelnum)
      integer recnum, innum, jumpnum, logicrecnum
      integer d(12)

      data  d / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

      call get_r_missing_data(r_missing_data,istatus)
      if ( istatus /= 1 ) then
         write (6,*) ' error getting r_missing_data'
         return
      endif

      recnum=    0
      innum=     0        ! innum : the record number within time window
      istatus=   0
      wmoid=     0
      elevation= 0.

      open ( 1, file=filename, status='old', err=1000 )

      istatus= 1

      do i= 1,loopnum
         read (1,5,end=99,err=29) reportflag, 
     ~                            latitudedummy, longitudedummy,
     ~                            iy, m1, id, ih, m2, logicrecnum
c wen add 5        format ( a3, i5, f4.0, 2f5.2, 2x, 5i2, i3 )
5        format ( a3,9x,i5,i5,2x, 5i2, i3 )
ccccc               write(*,*)'!!!!1111!',reportflag,iy,m1,id,ih,m2,logicrecnum

         if ( reportflag /= '*15' )  then
            write (6,*) 
     ~            ' error reading drpsnd data of identification -reject'
	    write (6,*) reportflag, latitudedummy, longitudedummy

            jumpnum= logicrecnum -1
            do 11 k= 1,jumpnum
11             read (1,*)
            go to 51
         endif

c               ------ creat a9time in yydddhhmm format ------
         if ( m1 == 2  .and.  mod(iy,4) == 0 )  d(m1)= d(m1) +1
	  
         if ( m2 /= -9 )  then   ! minus 30 mins to obtain the time in the air
            m2= m2 -30

            if ( m2 .lt. 0 )  then
               m2= m2 +60
               ih= ih -1

               if ( ih .lt. 0 )  then
                  ih= 23
                  id= id -1

                  if ( id .lt. 1 )  then
                     m1= m1 -1
                         
                     if ( m1 .lt. 1 )  then
                        m1= 12
                        iy= iy -1
                     endif

                     id= d(m1)
                  endif
               endif
            endif

         else         ! 00:-9 23:-9 -> 00:00 as the time in the air
            m2= 0     ! 12:-9 11:-9 -> 12:00 
	    if ( ih == 11  .or.  ih == 23 )  ih= ih +1

            if ( ih >= 24 )  then
               ih= 0
               id= id +1

               if ( id > d(m1) )  then
                  id= 1
                  m1= m1 +1
                      
                  if ( m1 > 12 )  then
                     m1= 1
                     iy= iy +1
                  endif
               endif
            endif

         endif
            
         call i2a ( iy, yy )
         call i2a ( m1, mo )
         call i2a ( id, dd )
         call i2a ( ih, hh )
         call i2a ( m2, mn )

         time= yy//mo//dd//hh//mn
         a9timedummy= a10_to_a9(time,istatus)
         call cv_asc_i4time( a9timedummy, i4time_drpsnd )

c          ----------    test if drpsnd is within time window    ----------
         if ( i4time_drpsnd /= 0 ) then    
            if ( i4time_drpsnd >= i4time_drpsnd_earliest .and.
     ~           i4time_drpsnd <= i4time_drpsnd_latest )  then
	       write (6,*) reportflag, latitudedummy, longitudedummy,
     ~                     yy, mo, dd, hh, mn, logicrecnum,
     ~                     ' inside time window'
	       innum= innum +1
cc  wen modify
c	       latitude(innum)= latitudedummy
c	       longitude(innum)= longitudedummy
	       lat1= latitudedummy/100.
	       lon1= longitudedummy/100.

	       latitude(innum)= lat1
	       longitude(innum)= lon1
	       a9time(innum)= a9timedummy

               layernum(innum)= logicrecnum -5
               do j= 1,layernum(innum)
                  read (1,15,err=19,end=99) pressure(innum,j),
     ~              height(innum,j), heightqua(innum,j),      
     ~              temperature(innum,j), temperaturequa(innum,j),        
     ~              tempdewdiff(innum,j), dewpointqua(innum,j),
     ~              winddir(innum,j),windspeed(innum,j),windqua(innum,j)      
15                format ( 2x, f5.1, f5.0, i2, 2(f4.1,i2), 2f3.0, i2 )
	          go to 20

19                write (6,*)' error reading variables of drpsnd data'
                  do k= 1,j
                     write (6,*) pressure(innum,k),
     ~                    height(innum,k), heightqua(innum,k),       
     ~                    temperature(innum,k), temperaturequa(innum,k),
     ~                    tempdewdiff(innum,k), dewpointqua(innum,k),
     ~                    winddir(innum,k), windspeed(innum,k),
     ~                    windqua(innum,k)
                  enddo
20             enddo

               read (1,*)
               read (1,21) taskname(innum), 
     ~                     prshtdbtl(innum), prshtdbth(innum)
               read (1,22) prstmdbtl(innum), prstmdbth(innum)
               read (1,*)
21             format ( 5x, a5, 10x, 2f5.1 )
22             format ( 2f5.1 )
    	       goto 50

            else
               write (6,*) reportflag, latitudedummy, longitudedummy,
     ~                     yy, mo, dd, hh, mn, logicrecnum,
     ~                     ' outside time window -reject'
    	       goto 40

            endif
         endif

29       write (6,*) ' error reading drpsnd codes of stations -reject'
	 write (6,*) reportflag, latitudedummy, longitudedummy,
     ~               iy, m1, id, ih, m2, logicrecnum
	 do k= 1,levelnum
            read (1,'(a2)') flag
	    if ( flag == '25' )  go to 50
	 enddo

40       jumpnum= logicrecnum -1
         do 41 k= 1,jumpnum
41          read (1,*) 

50       recnum= recnum +1
51    enddo

c      ----------     examing data quality and changing units     ---------    
c      when elevation is missing, return -999. without change for the sake of 
c      format f15.0 in snd files
99    do 100 i= 1,innum
      do 100 j= 1,layernum(i)

         if ( pressure(i,j) == -999. )  pressure(i,j)= r_missing_data
         if ( heightqua(i,j) /= 1 )  height(i,j)= r_missing_data
         if ( temperaturequa(i,j) /= 1 ) temperature(i,j)=r_missing_data

         tmpflag= 0
         if     ( prshtdbtl(i) >= pressure(i,j)  .and.
     ~            prshtdbth(i) <= pressure(i,j) ) then
            height(i,j)= r_missing_data
         elseif ( prstmdbtl(i) >= pressure(i,j)  .and.
     ~            prstmdbth(i) <= pressure(i,j) ) then
            temperature(i,j)= r_missing_data
            tmpflag= 1.
         endif

         if ( temperaturequa(i,j) == 1 .and. dewpointqua(i,j) == 1 .and.
     ~        tmpflag == 0 ) then
               dewpoint(i,j)= temperature(i,j) -tempdewdiff(i,j)
            else
               dewpoint(i,j)= r_missing_data
         endif

c wen modi          if ( windqua(i,j) /= 1 )  then
         if ( windqua(i,j) .eq. 11 ) go to 100
         if ( windqua(i,j) .eq. 21 ) go to 100
         if ( windqua(i,j) .eq. 31 ) go to 100
             write(*,*)' wen test windqua',windqua(i,j),winddir(i,j)

            winddir(i,j)= r_missing_data
            windspeed(i,j)= r_missing_data
c wen modi         endif
100   continue

!      do 900 i= 1,innum
!	 write(*,895) wmoid(i), layernum(i), latitude(i), longitude(i),
!     ~                 elevation(i), taskname(i), a9time(i), 'drpsnd'
!895      format (i12, i12, f11.4, f15.4, f15.0, 1x, a5, 3x, a9, 1x, a8)
!
!         do 900 j= 1,layernum(i)
!            write (*,*) height(i,j), pressure(i,j), temperature(i,j),
!     ~                   dewpoint(i,j), winddir(i,j), windspeed(i,j)
!900   continue
       do 900 i= 1,innum
          latitude_out(i,:) = latitude(i)
          longitude_out(i,:) = longitude(i)
          a9time_out(i,:) = a9time(i)
          c8_obstype(i) = 'dropsnd '            ! note revised spelling
          c5_staid(i) = '     '
          staelev(i) = -999.                    ! dummy values
900   continue

!     call write_snd routine
      call write_snd(      lun_out                            ! i
     1                    ,loopnum,levelnum,innum             ! i
     1                    ,wmoid                              ! i
     1                    ,latitude_out,longitude_out,staelev ! i
     1                    ,c5_staid,a9time_out,c8_obstype     ! i
     1                    ,layernum                           ! i
     1                    ,height                             ! i
     1                    ,pressure                           ! i
     1                    ,temperature                        ! i
     1                    ,dewpoint                           ! i
     1                    ,winddir                            ! i
     1                    ,windspeed                          ! i
     1                    ,istatus)                           ! o



      write (6,*) ' found', innum, 
     ~            'stations available within time window in',
     ~            recnum, 'drpsnd stations'

1000  return
      end



