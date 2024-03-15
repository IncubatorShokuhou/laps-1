
      program genlvdlut

      implicit none

      integer  nx_l,ny_l
      integer  istatus

      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'
c
c ======================== start ==============================
c
      call get_grid_dim_xy(nx_l,ny_l,istatus)
      if(istatus.ne.1)then
         write(6,*)'error getting nx_l/ny_l'
         goto 1000
      else
         write(6,*)'laps nx_l and ny_l obtained'
      endif
c
c---------------------------------------------------------------
c
      call config_satellite_lvd(istatus)
      if(istatus.ne.1)then
         write(*,*)'error configuring satellite-master common'
         goto 1000
      endif
c
c------------------- main program -----------------------
c
      call genlvdlut_sub(nx_l,ny_l,istatus)
      if(istatus.lt.0)then
         write(*,*)'error generating lut'
         goto 1000
      endif

c ---------------------------------------------------------
      call rewrite_satellite_lvd_nl(istatus)
      if(istatus.ne.1)then
         write(*,*)'error returned from rewrite_satellite_master_nl'
         write(*,*)'check satellite_master.nl namelist in data/static'
         goto 1000
      endif

      print*, 'finished in genlvdlut'
1000  stop
      end
