      subroutine wrt_laps_static (dir, laps_dom_file, imax, jmax, 
     1                            n_grids, dx, dy, lov, latin1, 
     1                            latin2, origin, var, comment, 
     1                            data, model, grid_spacing,map_proj,
     1                            la1,lo1,la2,lo2,
     1                            center_lat, center_lon,lli,llj,uri,
     1                            urj,parent_id,ratio_2_parent,
     1                            status)

      implicit none

c this routine calls functions in static_routines.c.  to create the new
c static file, a cdl file named <laps_dom_file>.cdl (ie. nest7grid.cdl)
c must be located in the laps cdl directory and is written out to the
c directory specified in "dir".

c passed in variables
      character      dir*(*),laps_dom_file*(*)
      integer        imax, jmax, n_grids
      real         dx, dy, lov, latin1, latin2
      character*(*)  origin
      character*(*)  var(n_grids)
      character*(*)  comment(n_grids)
      real         data(imax,jmax,n_grids)
      character*(*)  model
      real         grid_spacing
      character*(*)  map_proj
      real           la1,lo1,la2,lo2
      real           center_lat,center_lon
      integer        lli,llj,uri,urj,parent_id,ratio_2_parent
      integer        status

!mp
      integer n
!mp

c local variables

      integer        i4time_now_gg, dom_len, map_len
      integer        origin_len, unixtime, nx_lp, ny_lp
      integer        var_len, com_len, cdl_dir_len , asc_len
      character*150  file_name, cdl_dir
      character*24   asctime
      character*30   map_projection
      
      integer        i4time,
     1               error(2),
     1               flag

      integer        f_len

      common         /prt/flag

      include       'grid_fname.cmn'

      error(1)=1
      error(2)=0

c  begin subroutine

      print*,'wrt_laps_static'
      print*,'make_static_fname'

	write(6,*) 'dir= ', dir
	write(6,*) 'laps_dom_file= ', laps_dom_file
	write(6,*) 'file_name= ', file_name
        write(6,*)' n_grids= ',n_grids
	write(6,*) 'f_len= ', f_len
      call make_static_fname(dir,laps_dom_file,file_name,f_len,status)
	write(6,*) 'after make_static_fname'
	write(6,*) 'dir= ', dir
	write(6,*) 'laps_dom_file= ', laps_dom_file
	write(6,*) 'file_name= ', file_name
	write(6,*) 'f_len= ', f_len
      if (status .eq. error(2)) goto 990

      print*,'call i4time_now_gg'
      i4time = i4time_now_gg()
      print*,'i4time = ',i4time
      call cv_i4tim_asc_lp(i4time,asctime,status)
      if(status .ne. 1)then
         print*,'error returned: cv_i4tim_asc_lp: ',i4time
         return
      endif

      unixtime = i4time - 315619200

c get nx_l and ny_l from namelist (or common if namelist read already).
      if(nest.eq.1)then
         write(6,*) 'call get_grid_dim_xy: ' 
         call get_grid_dim_xy(nx_lp, ny_lp, status)
	 write(6,*) 'return get_grid_dim_xy, status: ', status
         if (status .ne. 1) goto 930
      else
         nx_lp=imax
         ny_lp=jmax
      endif

      print*,'dimensions used in wrt_cdf_static: ',nx_lp,ny_lp
        
      asc_len = len(asctime)
      call s_len(laps_dom_file, dom_len)
      var_len = len(var(1))
      com_len = len(comment(1))
      origin_len = len(origin)

      map_len = len(map_projection)
      map_projection = ' '
      if (map_proj .eq. 'plrstr') 
     1  map_projection = 'polar stereographic'
      if (map_proj .eq. 'lambrt' ) then
        if (latin1 .eq. latin2) then
          map_projection = 'tangential lambert conformal'
        else
          map_projection = 'secant lambert conformal'
        endif
      endif
      if (map_proj .eq. 'merctr') map_projection = 'mercator'
      if (map_proj .eq. 'rotlat') map_projection = 'rotated lat-lon'

	write(6,*) 'map_projection: ', map_projection
	
	write(6,*) 'lov = ', lov

      if (lov .lt. 0.0) lov = 360.0 + lov
     
      call get_directory('cdl',cdl_dir, cdl_dir_len)

	write(6,*) 'call write_cdf_static using: '
	write(6,*) 'laps_dom_file= ', laps_dom_file
      call write_cdf_static(file_name,f_len,asctime,asc_len,
     1                      cdl_dir,cdl_dir_len,var,var_len,
     1                      comment,
     1                      com_len,laps_dom_file,dom_len,imax,
     1                      jmax,n_grids,nx_lp,ny_lp,data,model,
     1                      grid_spacing,dx,dy,lov,latin1,latin2,
     1                      la1,lo1,la2,lo2,
     1                      center_lat,center_lon,lli,llj,uri,urj,
     1                      parent_id,ratio_2_parent,
     1                      origin,origin_len,map_projection,map_len,
     1                      unixtime, status)

      if (status .eq. -2) goto 940
      if (status .eq. -3) goto 950
      if (status .eq. -4) goto 960
      if (status .eq. -5) goto 970
      if (status .eq. -6) goto 980
c
c ****  return normally.
c
      status = error(1)
999   return
c
c ****  error trapping.
c
930   if (flag .ne. 1)
     1write(6,*) 
     1'error getting info from ',laps_dom_file(1:dom_len),
     1'.parms..write aborted.'
      status = error(2)
      goto 999
c
940   if (flag .ne. 1)
     1write(6,*) 'error opening file to be written to...write aborted.'
      status = error(2)
      goto 999
c
950   if (flag .ne. 1)
     1   write(6,*) 'error in imax,jmax, or n_grids...write aborted.'
      status = error(2)
      goto 999
c
960   if (flag .ne. 1)
     1   write(6,*) 'error writing data to file...write aborted.'
      status = error(2)
      goto 999
c
970   if (flag .ne. 1)
     1write(6,*) 'error writing file header info...write aborted.'
      status = error(2)
      goto 999
c
980   if (flag .ne. 1)then
       if(trim(grid_fnam_common) == 'nest7grid')then
            write(6,*) 'x and y values in ',laps_dom_file(1:dom_len),
     1'.cdl do not match imax and jmax in ',laps_dom_file(1:dom_len),
     1'.parms ... write aborted.'
       else
            write(6,*) 'x and y values in ',laps_dom_file(1:dom_len),
     1'.cdl do not match imax and jmax computed from '
     1,laps_dom_file(1:5),'.nl ... write aborted.'
            write(6,*)'nx_lp,ny_lp | imax,jmax ',nx_lp,ny_lp,imax,jmax
       endif

       status = error(2)
       goto 999
      endif
c
990   if (flag .ne. 1) then
        write(6,*) 'length of dir+file-name is greater than 150 char.'
        write(6,*) 'static file cannot be accessed.'
      endif
      status = error(2)
      goto 999
c
      end

