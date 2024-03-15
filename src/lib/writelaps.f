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
      subroutine write_laps(i4_reftime,i4_valtime,dir,ext,
     1   imax,jmax,kmax,kdim,var,lvl,lvl_coord,units,
     1   comment,data,istatus)

c**********************************************************************
c
c      this file contains the following fortran subroutines:
c            write_laps
c
c      the write_laps_data subroutine reads the following fortran
c      subroutines from the readlapsdata.f file:
c            cvt_fname_v3
c
c      the write_laps_data subroutine reads the following c subroutines
c      from the rwl_v3.c file:
c            write_cdf_v3
c
c**********************************************************************
c
c      subroutine write_laps
c
c      author:    john snook
c      modified:  to write netcdf data files  1/93 linda wharton
c                 to remove byte arrays       4/94 linda wharton
c                 to accept netcdf ver. 3 data files  9/97 linda wharton
c                 to read levels from pressure.nl    2/01 linda wharton
c
c      writes data in arrays data and comment to the netcdf file name
c      specified by i4time, dir and ext.  the data in var, lvl, lvl_coord,
c      imax, jmax, kmax, kdim and units are stored into the netcdf file
c      when it is created.  istatus is returned.
c
c**********************************************************************
c
        implicit  none
c
      include       'grid_fname.cmn'

      integer      i4_reftime,           !input i4time of run
     1               i4_valtime,           !input i4time data is valid
     1               imax,jmax,kmax,       !input # cols, # rows, # fields
     1               kdim,                 !input k dimension of data array
     1               lvl(kdim),            !input level of each field 
     1               istatus               !output

      real         data(imax,jmax,kdim)  !input raw data to be written
      character*(*)  dir                   !input directory to be written to
      character*(*)  ext                   !input file name ext
      character*(*)  var(kdim)             !input 3 letter id of each field
      character*(*)  lvl_coord(kdim)       !input vertical coordinate of fields
      character*(*)  units(kdim)           !input units of each field
      character*(*)  comment(kdim)         !input comments for each field
c
      integer      flag,                 !print flag (1 = off)
     1               i_reftime,            !unix time of data
     1               i_valtime,            !unix time of data
     1               error(2),
     1               i,j,n7g_nx, n7g_ny,
     1               lgfc,
     1               ldf_len,leng,
     1               fn_length,
     1               var_len,
     1               comm_len,
     1               ext_len,
     1               asc_len,
     1               lvl_coord_len,
     1               units_len,
     1               cdl_path_len,
     1               stat_len,
     1               n_levels,
     1               max_levels,
     1               called_from,          !0=fortran, 1=c
     1               append                !0=no, 1=yes
c
      parameter (max_levels = 100)
      real         pr(max_levels),       !pressures read from get_pres_1d
     1               cdl_levels(max_levels),
     1               bott, top
c
      character*5    fcst_hh_mm
      character*9    gtime
      character*150  file_name
      character*150  cdl_path
      character*150  static_path
      character*9    laps_dom_file
      character*24   asctime
      character*40   v_g
c
      common         /prt/flag
c
c-------------------------------------------------------------------------------
c
      error(1)=1
      error(2)=0
c
c
c ****  call get_config to setup to read nest7grid.parms
c
      call get_config(istatus)

c ****  setup laps domain name
c
      call s_len(grid_fnam_common,lgfc)
      if (lgfc .gt. 0) then
        laps_dom_file = grid_fnam_common(1:lgfc)
      else
        laps_dom_file = grid_fnam_common
        write(6,*) 'domain name not retrieved by get_config'
        write(6,*) 'navigation info will not be written to',
     1' output file'
      endif

      call get_laps_dimensions(n_levels,istatus)
      if (istatus .ne. 1) then
         write (6,*) 'error getting vertical domain dimension'
         return
      endif

      call get_grid_dim_xy(n7g_nx, n7g_ny, istatus)
      if (istatus .ne. 1) then
          write(6,*) 'return get_grid_dim_xy, status: ', istatus
          return
      endif

      call get_vertical_grid(v_g,istatus)
      if (istatus .ne. 1) then
         write (6,*) 'error getting vertical grid'
         return
      endif

      call upcase(v_g, v_g)
  
      call s_len(v_g,leng)
      if (v_g(1:leng) .eq. 'pressure') then     
        call get_pres_1d(i4_valtime,n_levels,pr,istatus)
        do j = 1,n_levels
          pr(j)=pr(j)/100.
        enddo

      elseif (v_g(1:leng) .eq. 'sigma_p') then
        call get_sigma_1d(n_levels,pr,istatus)
        do j = 1,n_levels
          pr(j)=pr(j)*1000.
        enddo

      elseif (v_g(1:leng) .eq. 'sigma_ht') then
        call get_ht_1d(n_levels,pr,istatus) ! sigma heights fill the pr array
        if(istatus .ne. 1)then
            write(6,*)' error returned from get_ht_1d'
            goto 920
        endif

      else
        goto 920

      endif

c **** special case where write_laps is called with fua or fsf extension
      if (ext .eq. 'fua') then
        bott = pr(1)
        top = pr(n_levels)
        if (bott .gt. top) then
          do i = n_levels, 1, -1
            cdl_levels(i) = pr(i)
          enddo
        else
          do i = 1, n_levels
            cdl_levels(i) = pr(i)
          enddo
        endif
      endif

      if (ext .eq. 'fsf') then
        n_levels = 1
        cdl_levels(1) = 0
      endif

c ****  various checks on input data.
c
      if (kmax .gt. kdim) then
        if (flag .ne. 1)
     1write (6,*) 'illegal k dimension in data array...write aborted.'
        istatus=error(2)
        return
      endif
c
      if (imax .ne. n7g_nx) then
        if (flag .ne. 1) then
          if (laps_dom_file(1:lgfc) .eq. 'nest7grid') then
            write (6,*)
     1'imax passed in does not match nest7grid.parms...write aborted.'
          else
            write (6,*) 'imax passed in does not match '
     1,laps_dom_file,'.nl...write aborted.'
          endif
        endif
        istatus=error(2)
        return
      endif
c
      if (jmax .ne. n7g_ny) then
        if (flag .ne. 1) then
          if (laps_dom_file(1:lgfc) .eq. 'nest7grid') then
            write (6,*)
     1'jmax passed in does not match nest7grid.parms...write aborted.'
          else
            write (6,*) 'jmax passed in does not match '
     1,laps_dom_file,'.nl...write aborted.'
          endif
        endif
        istatus=error(2)
        return
      endif
c
c ****  get cdl_path
c
      call get_directory('cdl',cdl_path, cdl_path_len)
c
c ****  get static_path
c
      call get_directory('static',static_path, stat_len)
c
c ****  specify file name
c
      call make_fnam_lp(i4_reftime,gtime,istatus)
      if (istatus .ne. 1) then
        write (6,*)
     1'error converting i4time to file name...write aborted.'
        istatus=error(2)
        return
      endif

      call make_fcst_time(i4_valtime,i4_reftime,
     1                    fcst_hh_mm,istatus)

c
c ****  create ascii time variables.
c
      call cv_i4tim_asc_lp(i4_valtime,asctime,istatus)

      call s_len(ext, ext_len)

      call cvt_fname_v3(dir,gtime,fcst_hh_mm,ext,ext_len,
     1                  file_name,fn_length,istatus)
      if (istatus .eq. error(2)) goto 930


      called_from = 0    !called from fortran
      append = 0         ! only one analysis time allowed per file

      call s_len(laps_dom_file,ldf_len)
      var_len = len(var(1))
      comm_len = len(comment(1))
      lvl_coord_len = len(lvl_coord(1))
      units_len = len(units(1))
      asc_len = len(asctime)
      i_reftime = i4_reftime - 315619200
      i_valtime = i4_valtime - 315619200
c
c **** write out netcdf file
c
      print*,'writing file: ',file_name

!     if (v_g(1:leng) .eq. 'sigma_ht') then
!         write(6,*)' cdl_levels = ',cdl_levels
!     endif

      call write_cdf_v3 (file_name,ext,var,comment,asctime,cdl_path, 
     1                   static_path,laps_dom_file,ldf_len,fn_length,
     1                   ext_len,var_len, 
     1                   comm_len, asc_len, cdl_path_len, stat_len,
     1                   i_reftime, i_valtime,imax, jmax, kmax, kdim, 
     1                   lvl, data, pr, n_levels, cdl_levels,
     1                   called_from, append, istatus)
c
      if (istatus .gt. 0) goto 980
      if (istatus .eq. -2) goto 940
      if (istatus .eq. -3) goto 950
      if (istatus .eq. -4) goto 960
      if (istatus .eq. -5) goto 970
      if (istatus .eq. -6) goto 990
c
c ****  return normally.
c
        istatus=error(1)
999     return
c
c ****  error trapping.
c
920     if (flag .ne. 1) then
          write(6,*) 'write_laps aborted!'
          write(6,*) ' laps will currently only work on a pressure'
     1,' vertical grid'
          if (laps_dom_file(1:lgfc) .eq. 'nest7grid') then
            write(6,*) ' make sure vertical_grid is set to pressure'
     1,' in nest7grid.parms'
          else
            write(6,*) ' make sure vertical_grid is set to pressure'
     1,' in ',laps_dom_file,'.nl'
          endif
        endif
        istatus=error(2)
        goto 999

930     if (flag .ne. 1)
     1    write (6,*) 'file_name variable too short...write aborted.'
        istatus=error(2)
        goto 999
c
940     if (flag .ne. 1)
     1    write (6,*) 'error opening file to be written to...write abort
     1ed.'
        istatus=error(2)
        goto 999
c
950     if (flag .ne. 1)
     1    write (6,*) 'error in imax,jmax,or n_levels...write aborted
     1.'
        istatus=error(2)
        goto 999
c
960     if (flag .ne. 1)
     1    write (6,*) 'error writing data to file...write aborted.'
        istatus=error(2)
        goto 999
c
970     if (flag .ne. 1)
     1    write (6,*) 'error writing header info into file...write abort
     1ed.'
        istatus=error(2)
        goto 999
c
980     if (flag .ne. 1)
     1    write (6,*) 'some grids not written....could not convert laps
     1variables.', istatus
        istatus=error(2)
        goto 999
c
990     if (flag .ne. 1)
     1    write (6,*) 'file already exists for analysis time...write
     1aborted.'
        istatus=error(2)
        goto 999
c
        end

c##########################################################################
