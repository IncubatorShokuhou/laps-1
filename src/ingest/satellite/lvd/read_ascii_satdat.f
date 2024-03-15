      subroutine read_ascii_satdat(c_filename,
     &                       i4time_current,
     &                       c_filetime,i_delta_t,
     &                       nlines,nelems,
     &                       i,j,
     &                       rlat,rlon,
     &                       img_line,img_elem,
     &                       image_data,
     &                       i4time_data,
     &                       grid_spacing_km,
     &                       istatus)
c
      implicit none
      real cosd 
      integer     nlines,nelems

      integer     i,j
      integer     itb
      integer     img_lin
      integer     img_lin_prev
      integer     img_ele
      integer     img_line(nelems,nlines)
      integer     img_elem(nelems,nlines)
      integer     image_data(nelems,nlines)
      integer     i4time_current
      integer     i4time_data
      integer     i4time_diff
      integer     istatus
      integer     i_delta_t

      logical       eof

      real        xlat,xlon
      real        rlat(nelems,nlines)
      real        rlon(nelems,nlines)

      real        grid_spacing_deg
      real        grid_spacing_km

      integer     n_vars_req
      character*100 c_values_req
      character*40  c_vars_req
      character     c_date*5
      character     c_time*4
      character     c_filetime*9
      character     c_filename*255

      istatus=1
c
c open ascii file
c
      open(19,file=c_filename,form='formatted',status='old',
     &err=900)
c
c first line is date/time
c
      read(19,102,err=901)c_date,c_time
102   format(1x,a5,1x,a4)
      c_filetime=c_date//c_time
      write(6,*)c_filetime
c
c check that time in file is "current"
c
      call cv_asc_i4time(c_filetime,i4time_data)

      i4time_diff = i4time_current-i4time_data
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
         write(6,*)'found current data'
      endif
c
c read second line to set up the line and element counting
c
      read(19,*,err=901)xlat,xlon,img_lin_prev,img_ele,itb
      eof=.false.
      i=1
      j=1
      rlat(i,j)=xlat
      rlon(i,j)=xlon
      img_line(i,j)=img_lin_prev
      img_elem(i,j)=img_ele

      do while (.not.eof)

         read(19,*,end=29)xlat,xlon,img_lin,img_ele,itb
         if(img_lin.ne.img_lin_prev)then
            j=j+1
            i=0
            img_lin_prev=img_lin
         endif
 
         i=i+1
         rlat(i,j)=xlat
         rlon(i,j)=xlon
         img_line(i,j)=img_lin
         img_elem(i,j)=img_ele
         image_data(i,j)=itb

         goto 200

29       eof=.true.

200   enddo
100   format(f7.4,1x,f8.4,1x,i5,1x,i5,2x,i4)
c
c define 0.5 grid window for remapping
c
      grid_spacing_deg = sqrt( 
     1    (  rlat(1,2) - rlat(1,1)                   )**2
     1  + ( (rlon(1,2) - rlon(1,1))*cosd(rlat(1,1))  )**2 )
      grid_spacing_km = grid_spacing_deg*111.1               !1 deg lat = 111.1 km
c
      write(6,*)'grid spacing degrees ',grid_spacing_deg
      write(6,*)'grid spacing km :',grid_spacing_km
      write(6,*)

      goto 995

900   write(6,*)'error opening laps.asc data file'
      istatus=-1
      goto 1000

901   write(6,*)'error - initial read laps.asc'
      istatus=-1
      goto 1000

902   write(6,*)'old data. done in rd_atl_satdat '
      istatus=-1
      goto 1000

995   write(6,*)'finished. i/j totals: ',i,j
      write(6,*)
      write(6,*)'lat(1,1)/lon(1,1) ',rlat(1,1),rlon(1,1)
      write(6,*)'lat(i,1)/lon(i,1) ',rlat(i,1),rlon(i,1)
      write(6,*)'lat(1,j)/lon(1,j) ',rlat(1,j),rlon(1,j)
      write(6,*)'lat(i,j)/lon(i,j) ',rlat(i,j),rlon(i,j)

1000  return
      end
