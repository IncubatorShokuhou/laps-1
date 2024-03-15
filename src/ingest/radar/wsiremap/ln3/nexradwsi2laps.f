cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis 
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps 
cdis 
cdis    this software and its documentation are in the public domain and 
cdis    are furnished "as is."  the united states government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  they assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  if significant modifications or enhancements 
cdis    are made to this software, the fsl software policy manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
       subroutine nexradwsi_to_laps(filename,
     &                    msng_radar,
     &                    nlines,nelements,nlev,
     &                    imax,jmax,
     &                    lat,lon,
     &                    validtime,
     &                    remapped_prod,
     &                    istatus)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c routine reads netcdf wsi-nexrad data using subroutine read_wsi_cdf.
c nowrad data is then remapped to laps domain given lat/lon of domain.
c
       integer   nlev

       real  lat(imax,jmax)
       real  lon(imax,jmax)
       real  ri(imax,jmax),rj(imax,jmax)
       real  remapped_prod(imax,jmax)
       integer sum(imax,jmax)
       integer count(imax,jmax)
       integer i_max_value(imax,jmax)
       integer validtime
       integer istatus
       integer jstatus
       integer lstatus
       integer data_levels(nlev)
       integer il(nlev)
       integer ul
       integer ndvalue
       integer num_levels
       integer msng_radar

       integer ickint,itotwait,iageth

       character filename*200
       character c_data_type*2
       character level_prefix(nlev)*50
       character lprefix(nlev)*2
       character cpathwsi3d*200

       real dgtord
       real rlat1, rlon1
       real rlat2, rlon2
       real rdlat, rdlon
       real starti,startj
       real endi,  endj

       real r_missing_data

       integer image(nelements,nlines)
       integer imdata

       real  rlatc,rlonc
       real  nw(2),se(2)

       common /cegrid/nx,ny,nz,nw,se,rlatc,rlonc
c
c ***************************************************************************
c **************************  start routine *********************************
c
      istatus=1

      n=index(filename,' ')
      c_data_type=filename(n-2:n-1)

      rdtodg=180.0/3.1415926
      dgtord=0.017453292
c
      call get_ref_base(ref_base,iostatus)
      call get_ref_base_useable(ref_base_useable,iostatus)
      call get_r_missing_data(r_missing_data,iostatus)
c
c read nimbus netcdf file
c
      write(6,*)'reading: ',filename(1:n)

      call rd_wsi_3dradar_cdf(filename,nlines,nelements,
     &rdlat,rdlon,rlat1,rlon1,rlat2,rlon2,centerlon,
     &validtime,nlev,data_levels,num_levels,level_prefix,
     &image,jstatus)

      if(jstatus.eq.-1)then
         write(6,*)'error reading wsi-nexrad data'
         goto 19
      end if
      write(6,*)'found 3d wsi data. valid time: ',validtime
      print*
c
c load common block for cylindrical equidistant grid
c
      rdlat=rdlat*rdtodg
      rdlon=rdlon*rdtodg
      rlatc=rlat1-(rdlat*nlines*0.5)
      nx=nelements
      ny=nlines
      nz=1
c     rlatc=25.0
      rlonc=centerlon*rdtodg
      nw(1)=rlat1
      nw(2)=rlon1
      se(1)=rlat2
      se(2)=rlon2
c
c compute i/j starting/ending points for this domain
c
      call latlon_2_ceij(imax*jmax,lat,lon,ri,rj)

      starti= 99999.
      startj= 99999.
      endi  =-99999.
      endj  =-99999. 
      do j=1,jmax
         starti=min(ri(1,j),starti)
         endi  =max(ri(imax,j),endi)
      enddo
      do i=1,imax
         startj=min(rj(i,jmax),startj)
         endj  =max(rj(i,1),endj)
      enddo
      istart=max(1,nint(starti+0.5))
      jstart=max(1,nint(startj+0.5))
      iend  =min(nelements,nint(endi+0.5))
      jend  =min(nlines,nint(endj+0.5))

      write(6,*)' istart/iend ',istart,iend
      write(6,*)' jstart/jend ',jstart,jend
      write(6,*)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      print*,'useable levels = ',num_levels
      write(6,*)(data_levels(i),i=1,num_levels)
c
c  rlon1 and rlat1 already in degrees
c
c  remap nowrad data onto laps grid.  we now have the 
c  section of the nowrad grid in units of dbz which 
c  is in and near the laps domain.  this area is defined by
c  jstart, jend, istart, iend.
c
       do j=1,jmax
       do i=1,imax
          count(i,j)=0
          sum(i,j)=0
          i_max_value(i,j)=0
       enddo
       enddo

       iobcnd=0
       iobcm=0
       iobcn=0
       iobcu=0
       itot =0
       igood=0
       incnt=0
       do j=jstart,jend
          wsi_lat=rlat1-rdlat*(j)
          do i=istart,iend
             wsi_lon=rlon1+rdlon*(i)
c
c correct for signed integer
c
             if(image(i,j).lt.0)image(i,j)=256+image(i,j)

             call latlon_to_rlapsgrid(wsi_lat,wsi_lon,
     &                                lat,lon,imax,jmax,rii,rjj,
     &                                jstatus)

             if(jstatus.eq.1)then

               ii=nint(rii)
               jj=nint(rjj)

               itot=itot+1
               imdata=image(i,j)+1
               if(imdata.le.num_levels.and.imdata.ne.msng_radar+1)then

                   igood=igood+1
                   sum(ii,jj)=sum(ii,jj)+data_levels(imdata)
                   count(ii,jj)=count(ii,jj)+1
                   i_max_value(ii,jj)=max(i_max_value(ii,jj),
     &                                    data_levels(imdata))
               elseif(imdata.eq.msng_radar+1)then
                 iobcm=iobcm+1
               elseif(imdata.gt.num_levels.and.
     &                imdata.lt.msng_radar+1)then
                 iobcn=iobcn+1
                 incnt=imdata+incnt
               else
                 iobcu=iobcu+1 
               end if

             end if
          end do       !all i within window
       end do       !all j within window.
c
c check for echo tops or vil product
c
       if(c_data_type.ne.'et'.and.c_data_type.ne.'vi')then
          do j = 1,jmax
          do i = 1,imax
             if(count(i,j) .gt. 0)then
                if(sum(i,j) .gt. 0)then
                   remapped_prod(i,j)=float(sum(i,j))/
     &                                float(count(i,j))

c               remapped_prod(i,j)=float(i_max_value(i,j))

                else
                   remapped_prod(i,j)=ref_base
                endif
             else
                remapped_prod(i,j)=r_missing_data
             end if
          end do
          end do
       else
          do j = 1,jmax
          do i = 1,imax
             if(count(i,j) .gt. 0)then
                if(i_max_value(i,j).gt.0)then
                   remapped_prod(i,j)=float(i_max_value(i,j))
                else
                   remapped_prod(i,j)=r_missing_data
                endif
              else
                remapped_prod(i,j)=r_missing_data
              endif
          end do
          end do
       endif
       write(6,*)
       write(6,*)'       wsi data attributes:'
       write(6,*)' ------------------------------------'
c      write(6,*)'  pixels no data (nd) : ',iobcnd
c      write(6,*)'    (%): ',float(iobcnd)/float(itot)
       write(6,*)'  pixels = missing: ',iobcm
       write(6,*)'    (%): ',float(iobcm)/float(itot)
       write(6,*)'  pixels > nlevs : ',iobcn
       write(6,*)'    (%): ',float(iobcn)/float(itot)
       if(iobcn.gt.0)then
          write(6,*)'      avg pix value: ',
     &              float(incnt)/float(iobcn) 
       endif
       write(6,*)'  pixels unaccounted for: ',iobcu
       write(6,*)'    (%): ',float(iobcu)/float(itot)
       write(6,*)'  pixels useable : ',igood
       write(6,*)'    (%): ',float(igood)/float(itot)
       write(6,*)' ------------------------------------'

       goto 16
19     write(6,*)'error in nowrad_2_laps, terminating'
       goto 16
14     write(6,*)'nowrad data not found for given time'

16     write(6,*)'finished in nowrad_2_laps'
       return
       end
