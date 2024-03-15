      subroutine get_sounding_info_cdf(csat_id,
     &                         c_sounding_path,
     &                         i4time_latest,
     &                         c_filename_sat,
     &                         ndimx,ndimy,ndimch,
     &                         istatus)
c
c subroutine determines if goes8/9 satellite sounder data are available
c in fsl's /public data base. these are raw gvar satellite sounder
c files, 16 bit, 8 km resolution.
c it is then determined which sat files have not yet been processed
c for laps and returns these files for further mapping/navigation.
c it turns out that only the goes9 eastern sector covers the roc-laps
c domain.
c 
c any port of this software to another location, or perhaps a larger
c size than roc-laps will require adjustments (software) to determine
c which sector properly maps into the domain.
c
c there is no preliminary installation step atm to make this sector
c determination.
c
c j. smart   2/96     original version
c j. smart   10/97    modify for newlaps structure;
c                     modify for dynamic memory;
c                     modify to accept both goes8 and goes9.
c j. smart   02/02    add allocate statement for raw sounder data array
c                     allowing dynamic adjustment depending on file for
c                     given satellites.
      implicit none

      integer     max_files
      parameter  (max_files=300)
      integer     i_snddata_age
      parameter  (i_snddata_age=20000)
      integer     nch

      character*255 c_filename_sat
      character*255 c_filename_lsr
      character*255 c_filename(max_files)
      character*150 c_lsr_dir
      character*255 pathname
      character*255 c_sounding_path
      character*9   c_filetime_sat
      character*9   c_filetime_lsr
      character*9   a9_time_sys
      character*9   c_fname
      character*6   csat_id

      integer     i,j,k,n,jj
      integer     istatus
      integer     mstatus
      integer     gfn_status
      integer     nfiles_sat
      integer     numoffiles
      integer     lend
      integer     laps_cycle_time
      integer     ndimx,ndimy,ndimch
      integer     imci4
      integer     itstatus
      integer     init_timer
      integer     ishow_timer

      integer     i4time
      integer     i4time_min
      integer     i4time_proc
      integer     i4time_latest
      integer     i4time_nearest_lsr
      integer     i4time_now_gg
      integer     i4time_sys

c
      character*1   imc(4)
c
c ----------------------------- start --------------------------------------
c
      istatus = 1
c
c find files for goes-8 sounding data sector over colorado. this is
c the _sdr file at ~ 45 or 46 past the hour..
c
      call get_systime(i4time_sys,a9_time_sys,istatus)
      i4time_proc=i4time_now_gg()

c this is designed to allow archive data runs!
      call get_laps_cycle_time(laps_cycle_time,istatus)
      if(i4time_proc-i4time_sys .gt. 2*laps_cycle_time)then
         print*,'set i4time to contents of systime.dat'
         i4time_proc=i4time_sys
      endif


      call get_directory('lsr',c_lsr_dir,lend)
      c_lsr_dir=c_lsr_dir(1:lend)//'/'//csat_id//'/'

      n=index(c_sounding_path,' ')-1
      c_filename_sat=c_sounding_path(1:n)//'*_sdr'
      call get_file_names(c_filename_sat,numoffiles,c_filename
     1        ,max_files,istatus)
      i4time_min=99999
      do i = 1,numoffiles
         j=index(c_filename(i),' ')-5
         jj=j-8
         call cv_asc_i4time(c_filename(i)(jj:j),i4time)
         if(abs(i4time-i4time_proc).lt.i4time_min)then
            i4time_min=abs(i4time-i4time_proc)
            i4time_latest = i4time
         endif
      enddo

c     do i = 1,numoffiles
c        j=index(c_filename(i),' ')-5
c        jj=j-8
c        if(csat_id.eq.'goes08')then
c           if(c_filename(i)(j-1:j-1).eq.'4')then
c              call cv_asc_i4time(c_filename(i)(jj:j),i4time)
c              if(i4time.gt.i4time_latest)then
c                 i4time_latest = i4time
c              endif
c           endif
c        elseif(csat_id.eq.'goes09')then
c           if(c_filename(i)(j-1:j-1).eq.'0')then
c              call cv_asc_i4time(c_filename(i)(jj:j),i4time)
c              if(i4time.gt.i4time_latest)then
c                 i4time_latest = i4time
c              endif
c           endif
c        else
c           write(6,*)'unknown satellite type - returning to main'
c           goto 1000
c        endif
c     enddo

      if(abs(i4time-i4time_latest).gt.i_snddata_age)then
c
c put in wait for data here
c
         write(6,*)'data is old'
         write(6,*)'returning to main without new sounder data'
         istatus = -1
         goto 1000
      else
         print*,'i4time age difference = ',i4time_proc-i4time_latest
      endif
c
c we have "current" data. determine if it has already been processed.
c get latest lsr file (ie., in lapsprd/lsr)
c
      call make_fnam_lp(i4time_latest,c_filetime_sat
     +,mstatus)
      c_filename_sat=c_sounding_path(1:n)//c_filetime_sat//'_sdr'

      call get_file_time(c_lsr_dir,i4time_proc,
     +i4time_nearest_lsr)

      if(i4time_latest-i4time_nearest_lsr.eq.0)then
         write(6,*)'lsr has already processed this sounder data.'
         write(6,*)'return to main without new data.'
         istatus = -1
         goto 1000
      endif

      call make_fnam_lp (i4time_nearest_lsr, c_filetime_lsr, mstatus)
c
c get the times for the current sounding data
c
      n=index(c_filename_lsr,' ')
      write(6,*)'filetime lsr: ',c_filetime_lsr
      n=index(c_filename_sat,' ')
      write(6,*)'filetime sat: ',c_filetime_sat
c
c----------------------------------------------
c !for sounder data the image motion compensation is off (ie.,=1)
c
      imci4 = 1

      call get_line_elem_sounder_cdf(c_filename_sat
     &,ndimx,ndimy,ndimch,istatus)
      if(istatus.eq.0)then
         print*,'error getting sounder file line/elem dims'
         print*,'terminating'
         stop
      endif

      if(istatus.ne.1)then
         write(6,*)'error occurred reading sounding db'
         write(6,*)'getlsrdata.f'
         return
      endif

      write(6,*)'attributes for sounding data obtained'

      goto 1000

996   write(6,*)'error'

1000  return
      end
