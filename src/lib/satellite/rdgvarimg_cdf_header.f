      subroutine rd_gvarimg_cdf_header(filename,
     &                              northwest_vis_pixel,
     &                              northwest_vis_line,
     &                              southeast_vis_pixel,
     &                              southeast_vis_line,
     &                              eastwestcycles,
     &                              eastwestincs,
     &                              northsouthcycles,
     &                              northsouthincs,
     &                              x_resolution,
     &                              y_resolution,
     &                              framestarttime,
     &                              imc,
     &                              orbitattitude,
     &                              nx4,ny4,
     &                              x_step,y_step,
     &                              istatus)
c
      parameter (nvars=67) !number of variables
      include 'netcdf.inc'
c
      integer rcode
      integer varid

c     ****variables for this netcdf file****
c
      character*1 wavelength                     (   1,   3)
      integer   northwest_vis_pixel            (   1)
      integer   northwest_vis_line             (   1)
      integer   southeast_vis_pixel            (   1)
      integer   southeast_vis_line             (   1)
      integer   x_step,y_step
      real*8      framestarttime                 
      character*1 elem_dim                       (   1,   8)
      character*1 line_dim                       (   1,   8)
      character*1 imc                            (   4)
      integer   eastwestcycles                 
      integer   eastwestincs                   
      integer   northsouthcycles               
      integer   northsouthincs                 
      real*8      orbitattitude                  (   1, 336)
c     character*1 x_dim                          (   1,   8)
c     character*1 y_dim                          (   1,   8)
      integer   nx                             (   1)
      integer   ny                             (   1)
      real        lap                            (   1)
      real        lop                            (   1)
      integer   x_resolution
      integer   y_resolution
      real        dx                             (   1)
      real        dy                             (   1)
      character   filename*255
c
c*************************************
c
      integer start(10)
      integer count(10)
      integer nx4,ny4
      integer vdims(10) !allow up to 10 dimensions
      character*31 dummy
c
c ======  start =======
c
      rcode=nf_open(filename,nf_nowrite,ncid)

c
c    statements to fill wavelength                     
c
       rcode=nf_inq_varid(ncid,'wavelength',varid)
      if(rcode.ne.0)return
      call ncvinq(ncid, varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do  90 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
  90  continue
      rcode= nf_get_vara_text(ncid, varid,start,count,wavelength)
c
c    statements to fill northwest_vis_pixel            
c
       rcode=nf_inq_varid(ncid,'northwest_vis_pixel',varid)
      if(rcode.ne.0)return
      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 100 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 100  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,northwest_vis_pixel)
c
c    statements to fill northwest_vis_line             
c
       rcode=nf_inq_varid(ncid,'northwest_vis_line',varid)
      if(rcode.ne.0)return
      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 110 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 110  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,northwest_vis_line)
c
c    statements to fill southeast_vis_pixel            
c
       rcode=nf_inq_varid(ncid,'southeast_vis_pixel',varid)
      if(rcode.ne.0)return
      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 120 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 120  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,southeast_vis_pixel)
c
c    statements to fill southeast_vis_line             
c
       rcode=nf_inq_varid(ncid,'southeast_vis_line',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 130 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 130  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,southeast_vis_line)
c
c    statements to fill framestarttime                 
c
       rcode=nf_inq_varid(ncid,'framestarttime',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 150 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 150  continue
      rcode=nf_get_vara_double(ncid,varid,start,count,framestarttime)
c
c    statements to fill elem_dim                       
c
       rcode=nf_inq_varid(ncid,'elem_dim',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 230 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 230  continue
      rcode= nf_get_vara_text(ncid,varid,start,count,elem_dim)
c
c    statements to fill line_dim                       
c
       rcode=nf_inq_varid(ncid,'line_dim',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 240 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 240  continue
      rcode= nf_get_vara_text(ncid,varid,start,count,line_dim)
c
c    statements to fill imc                            
c
       rcode=nf_inq_varid(ncid,'imc',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 250 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 250  continue
      rcode= nf_get_vara_text(ncid,varid,start,count,imc)
c
c    statements to fill eastwestcycles                 
c
       rcode=nf_inq_varid(ncid,'eastwestcycles',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 390 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 390  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,eastwestcycles)
c
c    statements to fill eastwestincs                   
c
       rcode=nf_inq_varid(ncid,'eastwestincs',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 400 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 400  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,eastwestincs)
c
c    statements to fill northsouthcycles               
c
       rcode=nf_inq_varid(ncid,'northsouthcycles',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 410 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 410  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,northsouthcycles)
c
c    statements to fill northsouthincs                 
c
       rcode=nf_inq_varid(ncid,'northsouthincs',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 420 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 420  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,northsouthincs)
c
c    statements to fill orbitattitude                  
c
       rcode=nf_inq_varid(ncid,'orbitattitude',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 430 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 430  continue
      rcode=nf_get_vara_double(ncid,varid,start,count,orbitattitude)
c
c    statements to fill nx                             
c
       rcode=nf_inq_varid(ncid,'nx',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 510 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 510  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,nx)
      nx4 = nx(1)
c
c    statements to fill ny                             
c
       rcode=nf_inq_varid(ncid,'ny',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 520 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 520  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,ny)
      ny4=ny(1)
c
c    statements to fill lap                            
c
       rcode=nf_inq_varid(ncid,'lap',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 530 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 530  continue
      rcode=nf_get_vara_real(ncid,varid,start,count,lap)
c
c    statements to fill lop                            
c
       rcode=nf_inq_varid(ncid,'lop',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 540 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 540  continue
      rcode=nf_get_vara_real(ncid,varid,start,count,lop)
c
c    statements to fill x_resolution                     
c
       rcode=nf_inq_varid(ncid,'x_resolution',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 550 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 550  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,x_resolution)
c
c    statements to fill y_resolution
c
       rcode=nf_inq_varid(ncid,'y_resolution',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 551 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 551  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,y_resolution)
c
c    statements to fill dx                             
c
       rcode=nf_inq_varid(ncid,'dx',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 580 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 580  continue
      rcode=nf_get_vara_real(ncid,varid,start,count,dx)
c
c    statements to fill dy                             
c
       rcode=nf_inq_varid(ncid,'dy',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 590 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 590  continue
      rcode=nf_get_vara_real(ncid,varid,start,count,dy)
c
c    statements to fill x_step
c
      rcode=nf_inq_varid(ncid,'x_step',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 680 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 680  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,x_step)
c
c    statements to fill y_step 
c
      rcode=nf_inq_varid(ncid,'y_step',varid)
      if(rcode.ne.0)return

      call ncvinq(ncid,varid,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do 690 j=1,nvdim
      call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
      lenstr=lenstr*ndsize
      start(j)=1
      count(j)=ndsize
 690  continue
      rcode=nf_get_vara_int(ncid,varid,start,count,y_step)

      return
      end
