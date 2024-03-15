      subroutine update_gvarimg_parms(cd6,
     &                        cstype,l_cell_afwa,chtype,
     &                        cdir_path,
     &                        ewcycles,
     &                        ewincs,
     &                        nscycles,
     &                        nsincs,
     &                        framestarttime,
     &                        imc,x_res,y_res,
     &                        nw_vis_pix,
     &                        nw_vis_line,
     &                        orbitattitude,
     &                        satsublat,
     &                        satsublon,
     &                        x_step,y_step,
     &                        nx,ny,
     &                        istatus)
c
c program reads the netcdf gvarimage satellite data file header
c or the afwa file header and returns the relevant information for updating
c the satellite namelist (/static/satellite.nl).
c
      implicit none

      integer   max_files
      parameter   (max_files = 20000)

c     integer   i
c     integer   max_sat
c     integer   max_channels

      integer   nw_vis_pix
      integer   nw_vis_line
      integer   nw_vis_pix_gwc
      integer   nw_vis_line_gwc
      integer   se_vis_pix
      integer   se_vis_line
      integer   image_depth
      integer   image_width
      integer   strbdy1,strbdy2
      integer   stpbdy1,stpbdy2
      integer   strtpix,strtline
      integer   stoppix,stopline
      integer   bepixfc,bescnfc,fsci
      integer   decimat
      integer   x_step,y_step

      real*8      framestarttime 
      real*8      getftime
      real*8      satsublat,satsublon
      real        golonsbp
      real        golatsbp
      real        goalpha
      character*1 imc(4)
      character*4 c_imc
      integer   ewcycles                 
      integer   ewincs                   
      integer   nscycles               
      integer   nsincs                 
      integer   x_res,y_res
      integer   imci4
c     integer   nch

      logical   l_cell_afwa

      real*8      orbitattitude(336)

      integer   i4time_now_gg
      integer   i4time_cur
      integer   i4time_nearest
      integer   i_obstime
      integer   ld
      integer   numoffiles
      integer   nf,nl,nn
c     integer   nsat
c     integer   n,ns,nc
      integer   n
      integer   istatus
      integer   gstatus
      integer   lstatus
      integer   nx,ny
c     integer   lend
c     integer   nstypes
c     integer   indx
c     integer   idum

c     character*150 dir
      character*150 c_filespec
      character*255 filename_cdf
      character*200 cdir_path
      character*255 fname_sat
      character*255 c_filenames(max_files)
c     character*200 cdum
      character     cfname*100
      character     c_afwa_fname*100
c     character*10  cmode
      character*9   c_fname_cur
      character*9   c_fname
      character*3   chtype
      character*3   cstype
      character*6   cd6
      character*2   image_type
c     character*1   ct
c
c ========================================================
c
      write(6,*)' subroutine update_gvarimg_parms...'

      istatus = -1  !bad status return

      i4time_cur = i4time_now_gg()
      call make_fnam_lp(i4time_cur,c_fname_cur,lstatus)
      if(lstatus.ne.1)goto 997
c
c find latest gvarimage data.  cstype = 'gvr' is fsl gvar data
c
      if(cstype.eq.'gvr')then

         nn=index(chtype,' ')-1
         if(nn.le.0)nn=3

         n=index(cdir_path,' ')-1
         fname_sat=cdir_path(1:n)//'*_'//chtype(1:nn)
         call get_file_time(fname_sat,i4time_cur,i4time_nearest)
c
c just in case there is no current data, skip this channel.
c
         if(i4time_nearest.gt.0)then

            call make_fnam_lp(i4time_nearest,c_fname,lstatus)
            if(lstatus.ne.1)goto 997
            fname_sat=cdir_path(1:n)//c_fname//'*_'//chtype(1:nn)

            call get_file_names(fname_sat,numoffiles,c_filenames
     1,max_files,lstatus)
            if(lstatus.ne.1)goto 998
            if(numoffiles.gt.0)then
               if(numoffiles.eq.1)then
                  filename_cdf=c_filenames(numoffiles)
               else
                  write(6,*)'warning: numoffiles > 1: = ',numoffiles
                  write(6,*)'fname_sat = ',trim(fname_sat)
                  filename_cdf=c_filenames(1)
               endif
            else
               nf=index(filename_cdf,' ')
               write(6,*)'no files found ',filename_cdf(1:nf)
               goto 899
            endif

         else

            nf=index(cdir_path,' ')
            write(6,*)'no data in ',cdir_path(1:nf)
            istatus = 1  !status flag indicating gvar header not read
            goto 900

         endif
c
c read header of current 11u gvarimage
c
         write(6,*)'reading gvar cdf header = ',chtype(1:nn)

         call rd_gvarimg_cdf_header(filename_cdf,
     &                        nw_vis_pix,
     &                        nw_vis_line,
     &                        se_vis_pix,
     &                        se_vis_line,
     &                        ewcycles,
     &                        ewincs,
     &                        nscycles,
     &                        nsincs,
     &                        x_res,y_res,
     &                        framestarttime,
     &                        imc,
     &                        orbitattitude,
     &                        nx,ny,
     &                        x_step,y_step,
     &                        lstatus)

         if(lstatus .ne. 1)then
            write(6,*)'error reading gvar cdf header ',chtype(1:nn)
            goto 996
         else
            write(6,*)'header info obtained = ',chtype(1:nn)
         endif
c
c compute satellite sub latitude and sub longitude.
c
c        if(imc(4).eq.' ')imc(4)='0'
c        read(imc(4),'(i1)')imci4
cc
c the bottom line:
c
         imci4 = 0

         write(6,*)'compute satsublat/satsublon - call sat_sublatlon...'       
         call sat_sublatlon(ewcycles,ewincs,nscycles,nsincs,
     &framestarttime,imci4,orbitattitude,satsublat,satsublon,lstatus)
         if(lstatus.ne.1)then
            goto 999
         endif
c
c ====================== gwc switch =======================
c =========== gwc switch removed on 12-2007: jrs ==========
c
      endif
      istatus = 0

      goto 900
 
101   write(6,*)'error getting ri/rj luts'
      goto 900

996   write(6,*)'error in rd_gvarimg_cdf_header: terminating'
      istatus = 1
      goto 900

997   write(6,*)'error in make_fname_lp'
      goto 900

998   write(6,*)'error in get_file_names'
      goto 900

999   write(6,*)'error in satsublatlon'
      istatus = 1
      goto 900

899   write(6,*)'parmfile not updated'

900   return
      end

c ===============================================================
c
      subroutine get_wfo_nav_parms(cpath,chtype,centerlat,
     &centerlon,rlat,rlon,rlatnxny,rlonnxny,rlatdxdy,rlondxdy,
     &dx,dy,nx,ny,istatus)

      include 'netcdf.inc'

      character*(*) cpath
      character     cfname_sat*200
      character     fname9_to_wfo_fname13*13
      character     cfname9*9
      character     cfname13*13
      character     chtype*3
      integer       istatus
      integer       i4time_cur
      integer       nx,ny
      real          centerlat
      real          centerlon
      real          rlat
      real          rlon
      real          rlatdxdy
      real          rlondxdy
      real          rlatnxny
      real          rlonnxny
      real          dx,dy

      istatus = -1  !bad status return

      i4time_cur = i4time_now_gg()
c
c find latest wfo satellite data.
c
      nn=index(chtype,' ')-1
      if(nn.le.0)nn=3
      n=index(cpath,' ')-1
      cfname_sat=cpath(1:n)//'*'

      call get_file_time(cfname_sat,i4time_cur,i4time_nearest)

      if(i4time_nearest.gt.0)then
         call make_fnam_lp(i4time_nearest,cfname9,lstatus)
         cfname13=fname9_to_wfo_fname13(cfname9)
         if(lstatus.ne.1)goto 997
         cfname_sat=cpath(1:n)//cfname13
      else
         write(6,*)'no data in ',cpath(1:n)
         goto 900
      endif
c
c read header of current 11u gvarimage
c
      lenf=index(cfname_sat,' ')-1
      print*,'opening ',cfname_sat(1:lenf)
      nf_status = nf_open(cfname_sat,nf_nowrite,nf_fid)

      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ',cfname_sat(1:lenf)
        istatus=-1
        goto 900
      endif

      write(6,*)'calling get_attribute_wfo ',chtype(1:nn)
      call get_attribute_wfo(nf_fid,centerlat,centerlon,rlat,
     &rlon,rlatnxny,rlonnxny,rlatdxdy,rlondxdy,dx,dy,nx,ny,lstatus)

      if(lstatus .lt. 0)then
         write(6,*)'no attributes returned: get_attribute_wfo ',chtype
     &(1:nn)
         goto 900
      endif

      istatus = 0
      goto 1000

900   write(6,*)'returning without new attributes'
      goto 1000

997   write(6,*)'error returned from make_fnam_lp'

1000  return
      end
