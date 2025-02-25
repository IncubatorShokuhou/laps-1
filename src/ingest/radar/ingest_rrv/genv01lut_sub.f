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
      subroutine genv01lut_sub(nx_l,ny_l,max_radars,radlat,radlon,
     &                      c_radar_name,radar_la1,radar_lo1,
     &                 nxv01,nyv01,dxv01,dyv01,ri,rj,istatus)
c
      implicit none

      integer*4 nx_l,ny_l
      integer*4 max_radars
      integer*4 nxv01,nyv01

      real*4    lat(nx_l,ny_l)
      real*4    lon(nx_l,ny_l)
      real*4    grid_spacing
      real*4    data(nx_l,ny_l,2)
      real*4    ri(nx_l,ny_l,max_radars)
      real*4    rj(nx_l,ny_l,max_radars)
      real*4    radar_la1(max_radars)
      real*4    radar_lo1(max_radars)
      real*4    dxv01,dyv01
      real*4    radlat(max_radars)
      real*4    radlon(max_radars)

      integer*4 i,j
      integer*4 istatus
      integer*4 lend
c
c dimensions for lat/lon
c
      character*125 comment_ll(2)
      character*10 units_ll(2)
      character*3 var_ll(2)
      character*150 dir_static
      character c_radar_name(max_radars)*4
c
c ======================== start ==============================
c
c acquire laps latitude and longitude arrays.
c -------------------------------------------------------------
      call get_directory('static',dir_static,lend)
      var_ll(1) = 'lat'
      var_ll(2) = 'lon'

      write(6,*)'get laps lat/lon grid'
      call rd_laps_static(dir_static,'nest7grid',nx_l,ny_l, 2,
     &     var_ll, units_ll, comment_ll, data, grid_spacing,
     &     istatus)

      if(istatus.eq.1)then
         write(6,*)'laps lat/lon grid obtained'
         write(6,*)
         do j=1,ny_l
            do i=1,nx_l
               lat(i,j)=data(i,j,1)
               lon(i,j)=data(i,j,2)
            end do
         end do
      else
         write(6,*)'error - unable to get lat/lon data'
         goto 900 
      end if
c-----------------------------------------------------------------
c
      do i = 1,max_radars

         call gen_llij_lut_v01(c_radar_name(i),radlat(i),radlon(i),
     &nx_l,ny_l,lat,lon,radar_la1(i),radar_lo1(i),nxv01,nyv01,
     &dxv01,dyv01,ri(1,1,i),rj(1,1,i),istatus)

         if(istatus.eq.1)then
            write(6,*)'lut generated ',c_radar_name(i)
            write(6,*)'********************************'
            write(6,*)
         else
            write(6,*)'lut not generated ',c_radar_name(i)
            write(6,*)
         endif

      enddo


900   return 
      end
