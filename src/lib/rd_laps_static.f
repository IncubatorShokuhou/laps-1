      subroutine rd_laps_static (dir,laps_dom_file,imax,jmax,n_grids,
     1                           var,units,comment,data,grid_spacing,
     1                           status)

      implicit none

c this routine calls functions in static_routines.c.  to create the new
c static file, a cdl file named <laps_dom_file>.cdl (ie. nest7grid.cdl)
c must be located in the laps cdl directory and is written out to the
c directory specified in "dir".

c passed in variables
      character      dir*(*), laps_dom_file*(*)
      integer        imax, jmax, n_grids
      character*(*)  var(n_grids)
      character*(*)  units(n_grids)
      character*(*)  comment(n_grids)
      real         data(imax,jmax,n_grids), grid_spacing
      integer        status

c local variables

      character*150   file_name
      integer        var_len, com_len, unit_len,
     1               error(3),
     1               no_laps_diag,
     1               flag

      integer        f_len

      logical        l_some_missing

      common         /prt/flag

      error(1)=1
      error(2)=0
      error(3)=-2
      l_some_missing = .false.

c  begin subroutine

      call make_static_fname(dir,laps_dom_file,file_name,f_len,status)
      if (status .eq. error(2)) goto 980
      
      var_len = len(var(1))
      com_len = len(comment(1))
      unit_len = len(units(1))

      call read_cdf_static(file_name,f_len,var,var_len,comment,com_len,
     1                     units,unit_len,imax,jmax,n_grids,data,
     1                     grid_spacing,no_laps_diag,status)

      if (status .ge. 0) then       !return from read with no errors

c       do i = 1, n_grids
c         comment(i) = c_comment(i)
c         units(i) = c_units(i)
c       enddo

        if (status .gt. 0) l_some_missing = .true.

        if(l_some_missing) then
          status=error(3) != -2
        else
          status=error(1) !=  1
        endif
      endif

      if (status .eq. -1) goto 940 !error opening file
      if (status .eq. -3) goto 950 !error in imax,jmax or n_grids
      if (status .eq. -4) goto 960 !error reading file data
      if (status .eq. -5) goto 970 !error reading header
c
c ****  return normally.
c
      status = error(1)
999   return
c
c ****  error trapping.
c
940   if (flag .ne. 1)
     1   write(6,*) 'error opening file to be read.'
      status = error(2)
      goto 999
c
950   if (flag .ne. 1)
     1write(6,*) 'error in imax,jmax, or n_grids...read aborted.'
      status = error(2)
      goto 999
c
960   if (flag .ne. 1)
     1   write(6,*) 'error reading file data.'
      status = error(2)
      goto 999
c
970   if (flag .ne. 1)
     1   write(6,*) 'error reading header info.'
      status = error(2)
      goto 999
c
980   if (flag .ne. 1) then
        write(6,*) 'length of dir+file-name is greater than 150 char.'
        write(6,*) 'static file cannot be accessed.'
      endif
      status = error(2)
      goto 999
c

      end

c########################################################################
      subroutine make_static_fname(dir,laps_dom_file,file_name,
     1                             f_len,status)
c
c**********************************************************************
c
c      subroutine make_static_fname
c
c      author:    linda wharton 4/93
c
c      inputed dir and laps_dom_file are converted to ascii values and
c      b_filename is created.
c
c**********************************************************************

      implicit  none

      integer end_dir, end_dom,
     1          error(2),
     1          status

      integer f_len

      character       dir*(*)       !directory to read data from
      character       laps_dom_file*(*)
      character*30    dn_laps_dom   !downcase of laps_dom_file

      character*(*)    file_name

      error(1)=1
      error(2)=0

      call downcase(laps_dom_file,dn_laps_dom)

c ******  find end of dir
c
      call s_len(dir,end_dir)
c
c ******  find end of laps_dom_file
c
      call s_len(laps_dom_file,end_dom)
c
c ******  find end of file_name
c
      f_len = len(file_name)
c
c ****  make file_name
c
      if (end_dir+end_dom+7 .gt. f_len) then
        status = error(2)
        goto 999
      else
        file_name = dir(1:end_dir)//'static.'//dn_laps_dom(1:end_dom)
        f_len = end_dir+7+end_dom
      endif
c
c ****  return normally.
c
      status = error(1)
999   return
      end

c##########################################################################
