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
      subroutine process_ir_satellite(i4time,
     &           ni,nj,lat,lon,
     &           n_ir_lines,n_ir_elem,
     &           r_grid_ratio,
     &           image_ir,
     &           r_llij_lut_ri,
     &           r_llij_lut_rj,
     &           c_type,
     &           ta8,tb8,tc8,
     &           istatus)
c
c**************************************************************************
c
c       routine to collect satellite data for the laps analyses.
c
c       changes:
c        p.a. stamus       12-01-92       original (from 'get_vas_bt' routine).
c                     01-11-93       write output in laps standard format.
c                     02-01-93       add snooze call for vis data.
c                     09-16-93       add bands 3,4,5,12 data.
c       j.r. smart    03-01-94	     implement on the sun.  adapt to ispan
c                             grids. this required removing all references to
c                             goes mdals/ground station satellite receiving.
c          "          10-26-94       modified for goes-8 data. no need to compute
c                                    radiance and btemp as this conversion is included
c                                    in icnt_lut.
c          "          11-28-95       extract ir processing from original sub (process_satellite)
c                                    to isolate the ir processing in one subroutine.
c          "           3-8-96        removed icnt_lut.
c
c       notes:
c       this program gets satellite data from ispan satellite database, does
c       some processing, then writes the data on the laps grids to the
c       lvd file.  lvd is written in standard laps netcdf form.
c       the ispan data will only initially contain  band-8 (11.2)
c       and visible data.  more bands should be available when goes-i becomes
c       operational.
c
c
c       variables:
c       ta8              ra       o       band 8 brightness temps (averaged)
c       tb8              ra       o       band 8 brightness temps (warm pixel)
c       tc8              ra       o       band 8 brightness temps (filtered)
c
c
c       note: for details on the filtering and averaging methods, see the
c       vasdat2 routine or talk to s. albers.
c
c****************************************************************************
c
       implicit none
c
       integer ni,nj
c
c..... grids to put the satellite data on.
c
       real ta8(ni,nj)
       real tb8(ni,nj)
       real tc8(ni,nj)
c 
c..... laps lat/lon files.
c
       real lat(ni,nj)
       real lon(ni,nj)
c
       integer n_ir_lines,n_ir_elem

       real r_llij_lut_ri(ni,nj)
       real r_llij_lut_rj(ni,nj)
       real r_grid_ratio
       real image_ir(n_ir_elem,n_ir_lines)
c
       integer istatus_a, istatus_f
       integer istatus_r
       integer istatus_w, istatus
       integer i,j,ik
       integer nn

       real badlow,badhigh
       real r_missing_data
c
       integer i4time,imax,jmax

       character c_type*3
c
c using the original lvd output required 14 fields. we will stay with this
c so to keep the lvd output and reduce (eliminate) impact on other processes.
c
       istatus = -1

c      bad = 1.e6 - 2.
c
c these values are bases upon the gvar cnt-to-btemp lookup tables.
c they are good estimates for satellites other than goes.
c
       if(c_type.eq.'wv '.or.c_type.eq.'iwv'.or.c_type.eq.'wvp')then
          badlow=148.486
          badhigh=291.182
       elseif(c_type.eq.'4u'.or.c_type.eq.'i39')then
c per eric gregow 22 aug 2008 -- get the ir channels ingested without any exclusion
c of points (that were previously set to missing values)
c         badlow=205.908
          badlow=185.
          badhigh=341.786
       elseif(c_type.eq.'11u'.or.c_type.eq.'i11')then
          badlow=112.105
          badhigh=341.25  !in line with goes08 cnt-to-btemp lut max value.
       else          !must be the 12u data
          badlow=110.611
          badhigh=336.347
       endif

       imax = ni
       jmax = nj

       call get_r_missing_data(r_missing_data,istatus)
       if(istatus.ne.1)then
          write(6,*)'error getting r_missing_data'
          goto 999
       endif

       do j = 1,jmax
       do i = 1,imax
          ta8(i,j) = r_missing_data
          tb8(i,j) = r_missing_data
          tc8(i,j) = r_missing_data
       enddo
       enddo
c
c.....  call the satellite data processing subroutine for each band required.
c
       write(6,900)c_type 
900    format('proc channel type ',a3,' satellite data.')
c --------------------------------------------------------------------------
       call satdat2laps_ir(imax,jmax,
     &                     r_grid_ratio,
     &                     r_missing_data,
     &                     image_ir,
     &                     r_llij_lut_ri,
     &                     r_llij_lut_rj,
     &                     n_ir_lines,n_ir_elem,
     &                     tb8,tc8,ta8,
     &                     istatus_r)

       if(istatus_r .ne. 1) then
          write(6,920)istatus_r, c_type
920       format(' +++ warning. bad status',i3          
     &          ,' from satdat2laps_ir for band ',a2)
       endif
c
c----------------------------------------------------------------------------
c.....       do a quick check on the data. 
c
       ik=0
       do j=1,jmax
       do i=1,imax
          if(ta8(i,j).le.badlow .or. 
     &       ta8(i,j).gt.badhigh)then
             ik=ik+1
             if(ik.le.25)print*,'ta8(',i,',',j,')= ',ta8(i,j)
             ta8(i,j) = r_missing_data
          endif
          if(tb8(i,j).le.badlow .or.
     &       tb8(i,j).gt.badhigh) tb8(i,j) = r_missing_data
          if(tc8(i,j).le.badlow .or.
     &       tc8(i,j).gt.badhigh) tc8(i,j) = r_missing_data
       enddo !i
       enddo !j
       call check(ta8,r_missing_data,istatus_a,imax,jmax)
       call check(tb8,r_missing_data,istatus_w,imax,jmax)
       call check(tc8,r_missing_data,istatus_f,imax,jmax)

       istatus = 1

       if(istatus_a .lt. 0)then
          write(6,910) istatus_a
          if(istatus_a .eq. -(ni*nj))then ! entire array is missing
              write(6,*)' entire ta8 array is missing'
              istatus = 0
          endif
       else
          write(6,*)'ta8 checked out ok'
       end if
       if(istatus_w .lt. 0)then
          write(6,911) istatus_w
          if(istatus_w .eq. -(ni*nj))then ! entire array is missing
              write(6,*)' entire tb8 array is missing'
              istatus = 0
          endif
       else
          write(6,*)'tb8 checked out ok'
       end if
       if(istatus_f .lt. 0)then
          write(6,912) istatus_f
          if(istatus_f .eq. -(ni*nj))then ! entire array is missing
              write(6,*)' entire tc8 array is missing'
              istatus = 0
          endif
       else
          write(6,*)'tc8 checked out ok'
       end if
910    format(' warning! ta8 check istatus_a = ',i8)
911    format(' warning! tb8 check istatus_w = ',i8)
912    format(' warning! tc8 check istatus_f = ',i8)

999   return
      end
