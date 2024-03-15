      subroutine get_cloud_drift_cwb_satob 
     ~           (i4time_sys, i4_window, nx_l, ny_l, filename, istatus)     

      parameter ( loopnum=35 )

      character*(*)  filename
      character*3    reportflag
      character*2    yy, mo, dd, hh, mn
      character*9    a9timeobs(loopnum), a10_to_a9
      character*10   time

      real  lat_a(nx_l,ny_l), lon_a(nx_l,ny_l), topo_a(nx_l,ny_l)
      real  latitude(loopnum), longitude(loopnum)
      real  pressure(loopnum), winddir(loopnum), windspeed(loopnum)

      integer pressurequa(loopnum), windqua(loopnum), recnum, innum

      call get_r_missing_data(r_missing_data,istatus)
      if ( istatus .ne. 1 )  then
         write (6,*) 'error getting r_missing_data'
         return
      endif

      recnum= 0
      innum= 0              !  innum : the record number within time window
      istatus= 0

      open ( 1, file=filename, status='old', err=1000 )

      istatus= 1

      do i= 1,loopnum
         read (1,40,end=99,err=8) reportflag, yy, mo, dd, hh, mn

         if ( reportflag .ne. '*61' )  then
            write (6,*) 'error reading satob code heading ', reportflag
	    read (1,*)
            go to 10
         endif

c               ------ creat a9timeobs in yydddhhmm format ------
         if ( yy(1:1) .eq. ' ' )  yy= '0'//yy(2:2)
         if ( mo(1:1) .eq. ' ' )  mo= '0'//mo(2:2)
         if ( dd(1:1) .eq. ' ' )  dd= '0'//dd(2:2)
         if ( hh(1:1) .eq. ' ' )  hh= '0'//hh(2:2)
         if ( mn(1:1) .eq. ' ' )  mn= '0'//mn(2:2)
 
         time= yy//mo//dd//hh//mn
         a9timeobs(i)= a10_to_a9(time,istatus)
	 if ( istatus .ne. 1 )  then
	    write (6,*) 'bad observation time - reject ', time
	 endif  

         call cv_asc_i4time( a9timeobs(i), i4time_obs )
         i4_residue= abs( i4time_obs -i4time_sys )

c          ----------    test if raob is within time window    ----------
	 if ( i4_residue .le. i4_window )  then
            write (6,*) 'inside time window ', time, i4_window
            innum= innum +1

            read (1,50,end=99,err=9) latitude(innum), longitude(innum),
     ~                  pressure(innum), pressurequa(innum),      
     ~                  winddir(innum), windspeed(innum), windqua(innum)
         else
            write (6,*) 'outside time window -reject ', time, i4_window
	    read (1,*)
         endif
         go to 10

8        write (6,*) 'error reading actual time of satob code ',
     ~               reportflag, yy, mo, dd, hh, mn
	 read (1,*)
	 go to 10

9        write (6,*) 'error reading variables of satob code '
         write (6,*) latitude(innum), longitude(innum),
     ~               pressure(innum), pressurequa(innum),
     ~               winddir(innum), windspeed(innum), windqua(innum)

10       recnum= recnum +1       
      enddo

40    format ( a3, 21x, 5a2 )
50    format ( 2f4.1, 2x, f3.0, i1, 4x, 2f3.0, i1 )

c      ----------       examine data quality and change units       ---------
99    do 100 i= 1,innum
         if ( pressurequa(i) .eq. 1 )  then
	    pressure(i)= pressure(i) *100.          !  unit: mb --> hpa
          else
	    pressure(i)= r_missing_data
         endif

         if ( windqua(i) .ne. 1 )  then
            winddir(i)= r_missing_data
            windspeed(i)= r_missing_data
         endif
100   enddo

      do 900 i= 1,innum
         call open_ext(11,i4time_sys,'cdw',istatus)
900      write (11,21) latitude(i), longitude(i), pressure(i), 
     ~                 winddir(i), windspeed(i), a9timeobs(i)
 21      format(f8.3,f10.3,f8.0,f6.0,f6.1,2x,a9)

      write (6,*) 'found', innum,'satob data within time window in',      
     ~            recnum, ' satob codes'
       
1000  return
      end
