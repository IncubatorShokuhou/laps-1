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
       subroutine satdat2laps_vis(
     &                  r_grid_ratio,
     &                  r_llij_lut_ri,
     &                  r_llij_lut_rj,
     &                  imax,jmax,
     &                  line_dim,elem_dim, ! image_vis array dimensions
     &                  image_vis,         ! satellite grid
     &                  sv,                ! model grid
     &                  istatus)
c
c.....       this is the chandran version returing visible values
c.....       modified by steve albers to average a window the size of a laps
c.....       grid box, giving it better resolution and faster run times.
c           changes:  p. stamus  11-13-92  install in operational sfc code.
c
c                     j. smart   3-1-94    install for using ispan feed for
c                                          satellite data.
c                                          allow analysis technique to vary
c                                          as function of input/output ratio
       implicit none

       integer max_elem,max_line
       parameter (max_elem = 100)
       parameter (max_line = 100)

       integer imax, jmax
       integer line_dim, elem_dim

       real sv(imax,jmax)                  ! model grid
       real t_array(max_line*max_elem)
       real image_vis(elem_dim,line_dim)   ! satellite grid
       real r_llij_lut_ri(imax,jmax)
       real r_llij_lut_rj(imax,jmax)
       real elem_mn,elem_mx
       real line_mn,line_mx

c      integer i_s(imax*jmax)
c      integer j_s(imax*jmax)

       integer npix
       integer maxpix
       integer i,j,ii,jj
       integer istart
       integer iend
       integer jstart
       integer jend
       integer istatus
       integer qcstatus
       integer insufdata
       integer icnt

       real r_missing_data
       real temp
       real r_grid_ratio
       real result
       logical lforce_switch, l_rhombus /.false./
c
c -----------------------------begin--------------------------------
c
       call get_r_missing_data(r_missing_data, istatus)
       call zero(sv,imax,jmax)
       istatus = 0
       qcstatus =0
c
c       grid_spacing_deg = sqrt( 
c    1       (  xlat(1,2) - xlat(1,1)                   )**2
c    1     + ( (xlon(1,2) - xlon(1,1))*cosd(xlat(1,1))  )**2 
c    1                         )   

c.....       define half of the window dimensions

c      wdw_lat =  grid_spacing_deg / 2.
c      wdw_lon = (grid_spacing_deg / 2.) / cosd(xlat(1,1))
c      write(6,*)' get vis: wdw_lat, wdw_lon = ',wdw_lat,wdw_lon

       insufdata=0
       lforce_switch=.false.
       icnt=0
       do j=1,line_dim
       do i=1,elem_dim
          if(image_vis(i,j).eq.r_missing_data)then
             icnt=icnt+1
          endif
       enddo
       enddo
       if(icnt.gt.(.1*elem_dim*line_dim))then
          lforce_switch=.true.
          print*,'more than 10% of data missing: '
     &,float(icnt)/float(imax*jmax)
          print*,'force grid point averaging in satir2laps'
       endif

       if(r_grid_ratio .lt. 0.5.or.lforce_switch)then

c      ------------------------------
c in this block the average pixel value is used for remapping the visible
c satellite to the output laps grid.
c
          write(6,*)'image ratio .lt. 0.5 ',r_grid_ratio
          write(6,*)'use pixel avg for vis count'
          do 10 j=1,jmax
          do 10 i=1,imax
            if(sv(i,j).ne.0.) go to 10
c
c compute the line and element for a window surrounding laps grid point,
c using a simple assumption that the satellite grid has the same orientation
c as the model with a 1:1 aspect ratio.
c
c this might also be done by projecting the model grid box onto the satellite
c grid as a rhombus, and determining whether each satellite point is included in
c the rhombus.
c
            if(r_llij_lut_ri(i,j).ne.r_missing_data.and.
     &         r_llij_lut_rj(i,j).ne.r_missing_data)then

              elem_mx = r_llij_lut_ri(i,j) + ((1./r_grid_ratio) * 0.5)
              elem_mn = r_llij_lut_ri(i,j) - ((1./r_grid_ratio) * 0.5)
              line_mx = r_llij_lut_rj(i,j) + ((1./r_grid_ratio) * 0.5)
              line_mn = r_llij_lut_rj(i,j) - ((1./r_grid_ratio) * 0.5)
              istart = nint(elem_mn+0.5)
              iend   = int(elem_mx)
              jstart = nint(line_mn+0.5)
              jend   = int(line_mx)

              if(istart .le. 0        .or. jstart .le. 0 .or.
     &           iend   .gt. elem_dim .or. jend   .gt. line_dim)then
c               write(*,*)'outside visible lat/lon sector'
c               write(*,1020)i,j
c1020	        format(1x,'laps grid (i,j) = ',i3,1x,i3)

                insufdata=insufdata+1
                sv(i,j)=r_missing_data

c               i_s(insufdata)=i
c               j_s(insufdata)=j

   	      else
 
c **** find the average visible pixels around grid point
c
                temp = 0.0
                npix = 0
                maxpix = 0

                if(.not. l_rhombus)then
                  do jj=jstart,jend
                  do ii=istart,iend
                    if(image_vis(ii,jj) .ne. r_missing_data)then
                      npix = npix + 1
                      t_array(npix) = image_vis(ii,jj)
                    endif
                  enddo ! ii
                  enddo ! jj
                endif ! l_rhombus

                if(npix .ge. maxpix)maxpix=npix
c
c added 11-6-95... a quality control test on the group of pixels used to
c                  derive the laps average pixel value.
c
                if(npix.ge.2)then

                  do ii=1,npix
                     temp = temp + t_array(ii)
                  enddo

                  sv(i,j) = temp/ npix

                elseif(npix.eq.1)then

                  sv(i,j) = t_array(npix)

                else

                  sv(i,j) = r_missing_data

                endif

cd              if(i .eq. i/10*10 .and. j .eq. j/10*10)then
cd                write(6,5555)i,j,wm,wc,npix,nwarm,sc(i,j)
cd5555            format(1x,2i4,2f10.2,2i5,f10.2)
cd              endif

              end if

            else
              sv(i,j)=r_missing_data
            endif

   10     continue ! i,j

          write(6,*)'max num vis pixels for avg: ',maxpix
          write(6,*)'number of vis satellite pixels modified'
          write(6,*)'           by statistical qc: ',qcstatus

       else
c      ----------------------------
c this block bilinearly interpolates four surrounding grid points to the
c output grid point.
c
          write(6,*)'grid ratio .ge. 0.5 ',r_grid_ratio
          write(6,*)'use bilinear interp for vis count'
          do 20 j=1,jmax
          do 20 i=1,imax

            if(sv(i,j).ne.0.) go to 20
c
c line/elem are floating point i/j positions in ispan grid for input lat/lon
c 
c bilinear_interp_extrap checks on boundary conditions and
c uses r_missing_data if out of bounds.
c

            if(r_llij_lut_ri(i,j).ne.r_missing_data.and.
     &         r_llij_lut_rj(i,j).ne.r_missing_data)then

               call bilinear_laps(
     &           r_llij_lut_ri(i,j),
     &           r_llij_lut_rj(i,j),
     &           elem_dim,line_dim,image_vis,
     &           result,istatus)

               if(result .lt. r_missing_data .and.
     &            result .gt. 0.0)then

                  sv(i,j) = result

               else

                  sv(i,j) = r_missing_data

               endif

            else
               sv(i,j)=r_missing_data
            endif

   20     continue ! i,j

        end if

c       write(6,1234) ib,i4vtime,ict
1234       format(1x,'band ',i4,' count for i4time ',i10,' is ',i8)

        if(insufdata.gt.0)then
           print*,'found ',insufdata,' points that are too'
           print*,'close to data edge to compute average'

c          do i=1,insufdata,10
c             write(6,55)(i_s(j),j_s(j),j=i,i+9)
c55            format(1x,'i/j ',10(i3,',',i3,1x))
c          enddo

        endif

        istatus = 1
c
        return
        end
