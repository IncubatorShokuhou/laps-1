      subroutine read_sounder_db_cdf(filename,
     &                               imax,jmax,nch,
     &                               sounding,
     &                               wavelength,
     &                               scalingbias,
     &                               scalinggain,
     &                               northwest_sdr_pixel,
     &                               northwest_sdr_line,
     &                               southeast_sdr_pixel,
     &                               southeast_sdr_line,
     &                               eastwestcycles,
     &                               eastwestincs,
     &                               northsouthcycles,
     &                               northsouthincs,
     &                               framestarttime,
     &                               linetimebegin,linetimeend,
     &                               imc,ires_x,ires_y,
     &                               orbitattitude,
     &                               istatus)
c
c
      implicit none
      include 'netcdf.inc'
      integer irec_max
      parameter (irec_max=5000)

      integer  nvars
      parameter (nvars=60) !number of variables
c     variable ids run sequentially from 1 to nvars= 60
      integer rcode
c     ****variables for this netcdf file****
c
c the following declarations are for testing using the original netcdf read
c routine. ie., read the entire array all at once
c

      integer   imax,jmax,nch
      integer   sounding(imax,jmax,nch)
      integer   i,j,k,n
      integer   dim_id_y

      real*8      wavelength      ( nch )
      real      scalingbias     (jmax,nch)
      real      scalinggain     (jmax,nch)
      real      scaling_rec     (irec_max)

      integer   northwest_sdr_pixel
      integer   northwest_sdr_line
      integer   southeast_sdr_pixel
      integer   southeast_sdr_line
      integer   varid
      character*1 imc                      (   4)
      real      imcenabletime
      integer   eastwestcycles                 
      integer   eastwestincs                   
      integer   northsouthcycles               
      integer   northsouthincs                 
      integer   ires_x,ires_y
      real*8      framestarttime                 
      real*8      orbitattitude            (336)
      real*8      linetimebegin(jmax,nch)
      real*8      linetimeend(jmax,nch)
      real*8      ltrec(irec_max)

      integer istatus
      integer ndsize_sb(nch)
      integer ndsize_ltb(nch)
      integer ndsize_lte(nch)
      integer ndsize, ncid

      integer nvdim
      integer ntp,nvs
      integer lenstr
      integer start(10)
      integer count(10)
      integer vdims(10) !allow up to 10 dimensions
      character*31 dummy
      character*255 filename
c
      istatus = 1

      rcode=nf_open(filename,nf_nowrite,ncid)
      if(rcode.ne.0)then
         n=index(filename,' ')
         write(6,*)'error openning netcdf file'
         write(6,*)'filename: ',filename(1:n-1)
         istatus = -1
         return
      endif
c
c code to get dimension size and read individual element of sounding array
c get dimensions for sounding array (x,y,lambda) [lambda is # of wavelengths]
c this code has now been subroutine-ized; rdimg_line_elem_sub.f.
c
c
      call rddata_line_elem(ncid,imax,jmax,nch,sounding,istatus)

      if(istatus .ne. 1)then
         write(6,*)'error reading sounding - rddata_img_line_elem'
         return
      endif
c -------------------------------------
c
c    statements to fill wavelength                     
c
      rcode=nf_inq_varid(ncid,'wavelength',varid)
      if(rcode.ne.0) return
      call ncvinq(ncid, varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do  20 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
  20  continue
      rcode=nf_get_vara_double(ncid,varid,start,count,wavelength)
c
c    statements to fill northwest_sdr_pixel            
c
      rcode=nf_inq_varid(ncid,'northwest_sdr_pixel',varid)
      if(rcode.ne.0) return
      call ncvinq(ncid, varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do  50 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
  50  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,northwest_sdr_pixel)
c
c    statements to fill northwest_sdr_line             
c
      rcode=nf_inq_varid(ncid,'northwest_sdr_line',varid)
      if(rcode.ne.0) return
      call ncvinq(ncid, varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do  60 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
  60  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,northwest_sdr_line)
c
c    statements to fill framestarttime                 
c
      rcode=nf_inq_varid(ncid,'framestarttime',varid)
      if(rcode.ne.0) return
      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 70 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 70   continue
      rcode=nf_get_vara_double(ncid,varid,start,count,framestarttime)
c
c    statements to fill imc                            
c
      rcode=nf_inq_varid(ncid,'imc',varid)
      if(rcode.ne.0) return
      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 180 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 180  continue
      rcode= nf_get_vara_text(ncid,varid,start,count,imc)
c
c    statements to fill imcenabletime                  
c
      rcode=nf_inq_varid(ncid,'imcenabletime',varid)
      if(rcode.ne.0) return
      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 260 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 260  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,imcenabletime)
c
c    statements to fill eastwestcycles                 
c
      rcode=nf_inq_varid(ncid,'eastwestcycles',varid)
      if(rcode.ne.0) return
      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 300 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 300  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,eastwestcycles)
c
c    statements to fill eastwestincs                   
c
      rcode=nf_inq_varid(ncid,'eastwestincs',varid)
      if(rcode.ne.0) return
      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 310 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 310  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,eastwestincs)
c
c    statements to fill northsouthcycles               
c
      rcode=nf_inq_varid(ncid,'northsouthcycles',varid)
      if(rcode.ne.0) return
      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 320 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 320  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,northsouthcycles)
c
c    statements to fill northsouthincs                 
c
      rcode=nf_inq_varid(ncid,'northsouthincs',varid)
      if(rcode.ne.0) return
      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 330 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 330  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,northsouthincs)
c
c    statements to fill orbitattitude                  
c
      rcode=nf_inq_varid(ncid,'orbitattitude',varid)
      if(rcode.ne.0) return
      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 340 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 340  continue
      rcode=nf_get_vara_double(ncid,varid,start,count,orbitattitude)
c
c new: 11-21-96. jrsmart. retrieve scalinggain and scalingbias. these needed
c to convert counts to useable brightness temps and radiances.
c
      dim_id_y = ncdid(ncid, 'y', rcode)

      rcode=nf_inq_varid(ncid,'scalingbias',varid)
      if(rcode.ne.0)then
         write(6,*)'error getting scalingbias '
         istatus = -1
      endif

      if(rcode.ne.0)then
         write(6,*)'error getting scalingbias id code - returning'
         istatus = -1
         return
      endif

      write(6,*)'reading scalingbias'
      do k = 1,nch  !# of channels in sounding database (= 19).

c        write(6,*)'reading scalingbias for channel ',k

         call ncdinq(ncid,dim_id_y,dummy,ndsize_sb(k),rcode)
         if(rcode.ne.0)then
            write(6,*)'error getting sb dimension - ndsize_sb'
         endif

         if(ndsize_sb(k).gt.jmax.or.ndsize_sb(k).le.0)then
            write(6,*)'sb dimension size > nlines_max or <= 0'
            write(6,*)'ndsize_db(k) = ', ndsize_sb(k)
            istatus = -1
            write(6,*)'returning to main, istatus = ',istatus
            return
         endif

         count(1)=ndsize_sb(k)
         start(2)=k
c
c read record
c
      rcode=nf_get_vara_real(ncid,varid,start,count,scaling_rec)
         if(rcode.ne.0)then 
            write(6,*)'error reading scaling database'
            istatus = -1
            write(6,*)'returning to main, istatus = ',istatus
         endif
c
c load values for all lines into the 3d output array
c
         do i = 1,ndsize_sb(k)
            scalingbias(i,k) = scaling_rec(i)
         enddo

      enddo
c
c  statements to fill scalinggain
c
      rcode=nf_inq_varid(ncid,'scalinggain',varid)
      if(rcode.ne.0)then
         write(6,*)'error getting scalinggain '
         istatus = -1
      endif

      write(6,*)'reading scalinggain'
      do k = 1,nch  !# of channels in sounding database (= 19).

c        write(6,*)'reading scalinggain for channel ',k

         call ncdinq(ncid,dim_id_y,dummy,ndsize_sb(k),rcode)
         if(rcode.ne.0)then
            write(6,*)'error getting sg dimension - ndsize_sb'
         endif

         if(ndsize_sb(k).gt.jmax.or.ndsize_sb(k).le.0)then
            write(6,*)'sb dimension size > jmax or <= 0'
            write(6,*)'nsize_sb(k) = ',ndsize_sb(k)
            istatus = -1
            write(6,*)'returning to main, istatus = ',istatus
            return
         endif

         count(1)=ndsize_sb(k)
         start(2)=k
c
c read record
c
      rcode=nf_get_vara_real(ncid,varid,start,count,scaling_rec)
         if(rcode.ne.0)then
            write(6,*)'error reading scaling database'
            istatus = -1
            write(6,*)'returning to main, istatus = ',istatus
         endif
c
c load values for all lines into the 3d output array
c
         do i = 1,ndsize_sb(k)
            scalinggain(i,k) = scaling_rec(i)
         enddo

      enddo
c
c 11-14-97 (j.smart). acquire linetimebegin and linetimeend.
c
      rcode=nf_inq_varid(ncid,'linetimebegin',varid)
      if(rcode.ne.0)then
         write(6,*)'error getting linetimebegin '
         istatus = -1
      endif
      write(6,*)'reading linetimebegin'
      do k = 1,nch  !# of channels in sounding database (= 19).
         call ncdinq(ncid,dim_id_y,dummy,ndsize_ltb(k),rcode)
         if(rcode.ne.0)then
            write(6,*)'error getting ltb dimension - ndsize_ltb'
         endif
         if(ndsize_ltb(k).gt.jmax.or.ndsize_ltb(k).le.0)then
            write(6,*)'ltb dimension size > jmax or <= 0'
            write(6,*)'nsize_ltb(k) = ',ndsize_ltb(k)
            istatus = -1
            write(6,*)'returning to main, istatus = ',istatus
            return
         endif
         count(1)=ndsize_ltb(k)
         start(2)=k
c
c read record, use scaling_rec for the i/o.
c
      rcode=nf_get_vara_double(ncid,varid,start,count,ltrec)
         if(rcode.ne.0)then
            write(6,*)'error reading scaling database'
            istatus = -1
            write(6,*)'returning to main, istatus = ',istatus
         endif
c
c load values for all lines into the 3d output array
c
         do i = 1,ndsize_ltb(k)
            linetimebegin(i,k) = ltrec(i)
         enddo
      enddo
c -----------
c linetimeend
c
      rcode=nf_inq_varid(ncid,'linetimeend',varid)
      if(rcode.ne.0)then
         write(6,*)'error getting linetimebegin '
         istatus = -1
      endif
      write(6,*)'reading linetimeend'
      do k = 1,nch  !# of channels in sounding database (= 19).
         call ncdinq(ncid,dim_id_y,dummy,ndsize_lte(k),rcode)
         if(rcode.ne.0)then
            write(6,*)'error getting lte dimension - ndsize_lte'
         endif
         if(ndsize_lte(k).gt.jmax.or.ndsize_lte(k).le.0)then
            write(6,*)'lte dimension size > jmax or <= 0'
            write(6,*)'nsize_lte(k) = ',ndsize_lte(k)
            istatus = -1
            write(6,*)'returning to main, istatus = ',istatus
            return
         endif
         count(1)=ndsize_lte(k)
         start(2)=k
c
c read record, use scaling_rec for the i/o.
c
      rcode=nf_get_vara_double(ncid,varid,start,count,ltrec)
         if(rcode.ne.0)then
            write(6,*)'error reading scaling database'
            istatus = -1
            write(6,*)'returning to main, istatus = ',istatus
         endif
c
c load values for all lines into the 3d output array
c
         do i = 1,ndsize_lte(k)
            linetimeend(i,k) = ltrec(i)
         enddo
      enddo
c
c    statements to fill x_resolution and y_resolution
c
      rcode=nf_inq_varid(ncid,'x_resolution',varid)
      if(rcode.ne.0) return
      call ncvinq(ncid, varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do  53 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
  53  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,ires_x)
c
      rcode=nf_inq_varid(ncid,'y_resolution',varid)
      if(rcode.ne.0) return
      call ncvinq(ncid, varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do  54 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
  54  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,ires_y)
c
c    statements to fill southeast_sdr_pixel
c
      rcode=nf_inq_varid(ncid,'southeast_sdr_pixel',varid)
      if(rcode.ne.0) return
      call ncvinq(ncid, varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do  80 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
  80  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,southeast_sdr_pixel)
c
c    statements to fill southeast_sdr_line
c
      rcode=nf_inq_varid(ncid,'southeast_sdr_line',varid)
      if(rcode.ne.0) return
      call ncvinq(ncid, varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do  90 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
  90  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,southeast_sdr_line)

      return
      end
c
c get x/y dimensions for satellite data
c
      subroutine get_line_elem_sounder_cdf(filename,nelems,nlines
     &,nchannels,istatus)

      implicit none

      integer nlines
      integer nelems
      integer nchannels
      integer istatus
      integer ln
      integer ncid
      integer dim_id_x
      integer dim_id_y
      integer rcode

      character    filename*(*)
      character*31 dummy

      include 'netcdf.inc'
c
c open file for reading
c
      istatus = 0

      call s_len(filename,ln)
      print*,'open file for x/y/ch dimension read'
      print*,'filename: ',filename(1:ln)
      rcode=nf_open(filename,nf_nowrite,ncid)
      if(rcode.ne.0)then
         print*,'error openning netcdf file'
         print*,'filename: ',filename(1:ln)
         return
      endif
c
c this is the number of lines
c
      dim_id_y = ncdid(ncid, 'y', rcode)
      if(rcode.ne.0)then
         print*,'error getting y id code - returning'
         return
      endif
      call ncdinq(ncid, dim_id_y,dummy,nlines,rcode)
      if(rcode.ne.0)then
         print*,'error getting y dimension - nlines'
         return
      endif
c
c get x dimension id
c
      dim_id_x = ncdid(ncid, 'x', rcode)
      if(rcode.ne.0)then
         print*,'error getting x id code - returning'
         return
      endif

      call ncdinq(ncid,dim_id_x,dummy,nelems,rcode)
      if(rcode.ne.0)then
         print*,'error getting x dimension - nelems'
         return
      endif
c
c get wavelength dimension id (this is the number of channels)
c
      dim_id_x = ncdid(ncid, 'wavelength', rcode)
      if(rcode.ne.0)then
         print*,'error getting wavelength id code - returning'
         return
      endif

      call ncdinq(ncid,dim_id_x,dummy,nchannels,rcode)
      if(rcode.ne.0)then
         print*,'error getting wavelength dimension - nchannels'
         return
      endif

      rcode = nf_close(ncid)
      if(rcode.ne.nf_noerr) then
        print *, nf_strerror(rcode)
        print *,'nf_close'
        return
      endif

      istatus = 1

      return
      end
