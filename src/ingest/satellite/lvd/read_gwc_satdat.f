      subroutine read_afgwc_satdat(c_filename,
     &                       isat,jtype,l_cell_afwa,
     &                       chtype,i_delta_t,
     &                       i4time_current,
     &                       nlines,nelems,
     &                       image_data,
     &                       i4time_data,
     &                       istatus)
c
      implicit none

      integer     cell_depth
      integer     cell_width
      parameter  (cell_depth=64,
     &            cell_width=256)

      integer     nlines,nelems

      integer     n
      integer     i4time_current
      integer     i4time_data
      integer     i4time_diff
      integer     istatus
      integer     i_delta_t
      integer     isat,jtype
      integer     nefi,nlfi
      integer     lun

      integer decimat
      integer strpix     !start pixel
      integer strscnl    !start scanline 
      integer stppix     !stop pixel
      integer stpscnl    !stop scanline
      integer reqobstm   !requested observation time
      integer bepixfc
      integer bescnfc
      integer fsci
      integer strbdy1    !requested start boundary 1
      integer strbdy2    !requested start boundary 2
      integer stpbdy1    !requested stop boundary 1 
      integer stpbdy2    !requested stop boundary 2
      real   golonsbp      !goes longitude subpoint
      real   golatsbp      !goes latitude subpoint
      real   goalpha       !goes alpha

      real   image_data  (nelems,nlines)

      integer width      !tracks in width of image
      integer depth      !tracks in depth of image
      character*2 imgtype  !image type

      character     c_filetime*9
      character     cfname*9
      character     chtype*3
      character     c_filename*255

      logical       lopen,lext
      logical       l_cell_afwa

      istatus=0
 
      n=index(c_filename,' ')-1

      inquire(file=c_filename,exist=lext,opened=lopen,number=lun)
      if(.not.lext)then
         print*,'file does not exist: ',c_filename(1:n)
         goto 1000
      endif
      if(lopen)then
         print*,'file is already open: ',c_filename(1:n)
         goto 1000
      endif
      call read_gwc_header(c_filename(1:n),l_cell_afwa,strpix,strscnl,
     +   stppix,stpscnl,reqobstm,imgtype,golatsbp,golonsbp,width,depth, 
     +   goalpha,strbdy1,strbdy2,stpbdy1,stpbdy2,bepixfc,bescnfc,
     +   fsci,decimat,istatus)
      if(istatus.eq.0)then
         write(6,*)'got gwc header info'
      else
         write(6,*)'error reading gwc header '
         goto 1000
      endif

c     write(c_filetime(1:7),111)reqobstm
c111   format(i7)
c     do i=1,7
c        if(c_filetime(i:i).eq.' ')then
c           c_filetime(i:i)='0'
c        endif
c     enddo
c     call make_fnam_lp(i4time_current,cfname,istatus)
c     c_yr = cfname(1:2)
c     c_jday = c_filetime(1:3)
c     c_hhmm = c_filetime(4:7)
c     c_filetime=c_yr//c_jday//c_hhmm
c     call cv_asc_i4time(c_filetime,i4time_data)

c
c put check in here for latest file time in lvd subdirectory.
c
      i4time_data=reqobstm

      i4time_diff = 0  !i4time_current-i4time_data    <---- don't forget to un-comment.
      if(i4time_diff.gt.i_delta_t)then
         write(6,*)'data is too old!'
         write(6,*)'data    time: ',c_filetime
         call make_fnam_lp(i4time_current,c_filetime,istatus)
         write(6,*)'current time: ',c_filetime
         goto 902
      elseif(i4time_diff.lt.0)then
         write(6,*)'data time is more recent than current time!?'
         write(6,*)'this does not make sense!'
      else
         write(6,*)'found afwa data'
      endif
c
c read afgwc binary data file
c
      if(l_cell_afwa)then
         nlfi=depth*cell_depth
         nefi=width*cell_width
      else
         nlfi=stpscnl-strscnl+1
         nefi=stppix-strpix+1
      endif

      call process_sdhs_gvar_sub(c_filename,chtype,isat,jtype,
     &nelems,nlines,image_data,nlfi,nefi,depth,width,istatus)

      if(istatus.eq.1)then
         write(6,*)'got afwa satellite data'
      else
         write(6,*)'error reading gwc sat data - process_sdhs_sub'
      endif

      goto 1000

902   write(6,*)'returning without new data'

1000  return
      end
