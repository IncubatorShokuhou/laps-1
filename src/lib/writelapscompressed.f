
      subroutine write_laps_compressed(i4time,dir,ext,imax,jmax,
     1   kmax,kdim,var,lvl,lvl_coord,units,comment,data,
     1   istatus)

c**********************************************************************
c
!       implicit  none
c
      integer      i4time,               !input i4time of data
     1               i4_valtime,
     1               imax,jmax,kmax,       !input # cols, # rows, # fields
     1               kdim,                 !input k dimension of data array
     1               lvl(kdim),            !input level of each field 
     1               istatus               !output

      real         data(imax,jmax,kdim)    !input raw data to be written

      integer, allocatable, dimension(:) :: array1 !local compressed array
      real,    allocatable, dimension(:) :: array2 !local compressed array

      character*(*)  dir                     !input directory to be written to
      character*(*)  ext                     !input file name ext
      character*(*)  var(kdim)               !input 3 letter id of each field
      character*(*)  lvl_coord(kdim)         !input vertical coordinate of fields
      character*(*)  units(kdim)             !input units of each field
      character*(*)  comment(kdim)           !input comments for each field
c
      integer      flag,                 !print flag (1 = off)
     1               i_reftime,            !unix time of data
     1               i_valtime,            !unix time of data
     1               error(2),
     1               i,j,n7g_nx, n7g_ny,
     1               lgfc,
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
     1               max_levels,	   !maximum vertical levels
     1               called_from,          !0=fortran, 1=c
     1               append                !0=no, 1=yes
c
      parameter (max_levels=100)
      real         pr(max_levels),       !pressures read from get_pres_1d
     1               cdl_levels(max_levels)
c
      logical l_check_encoding
c
      character*5    fcst_hh_mm
      character*9    gtime
      character*150  file_name
      character*150  cdl_path
      character*150  static_path
      character*24   asctime
      character*20   v_g
c
      common         /prt/flag
c
c
c-------------------------------------------------------------------------------
c
      error(1)=1
      error(2)=0

c
c ****  various checks on input data.
c
      if (kmax .gt. kdim) then
        if (flag .ne. 1)
     1write (6,*) 'illegal k dimension in data array...write aborted.'
        istatus=error(2)
        return
      endif
c
c ****  specify file name
c
      call make_fnam_lp(i4time,gtime,istatus)
      if (istatus .ne. 1) then
        write (6,*)
     1'error converting i4time to file name...write aborted.'
        istatus=error(2)
        return
      endif
c
c **** get actual reftime from gtime...
c
      i_reftime = i4time - 315619200
      i_valtime = i_reftime

c
c ****  create ascii time variables.
c
      i4_valtime = i_valtime +  315619200
      call cv_i4tim_asc_lp(i4_valtime,asctime,istatus)

      call s_len(ext, ext_len)

c fcst_hh_mm: hard wired as a place holder - will be used in filename  only
c   if write_laps_compressed is called on lga, lgb, fua, fsf, ram, rsf
c to fix this, call write_laps instead
      fcst_hh_mm = '0000'

      call cvt_fname_v3(dir,gtime,fcst_hh_mm,ext,ext_len,
     1                  file_name,fn_length,istatus)
      if (istatus .eq. error(2)) goto 930

      called_from = 0    !called from fortran
      append = 0         ! only one analysis time allowed per file

      var_len = len(var(1))
      comm_len = len(comment(1))
      lvl_coord_len = len(lvl_coord(1))
      units_len = len(units(1))
      asc_len = len(asctime)

!     allocate arrays then do run-length encoding
!     n_cmprs_max = imax*jmax*kdim
      n_cmprs_max = 2000000  
 800  write(6,*)' allocate arrays with size of ',n_cmprs_max

      allocate( array1(n_cmprs_max), stat=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' error: could not allocate array1'
          write(6,*)' try reducing n_cmprs_max from ',n_cmprs_max
          goto 950
      endif

      allocate( array2(n_cmprs_max), stat=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' error: could not allocate array2'
          write(6,*)' try reducing n_cmprs_max from ',n_cmprs_max
          goto 950
      endif

      ngrids = imax*jmax*kdim

      call runlength_encode(ngrids,n_cmprs_max,data           ! i
     1                     ,n_cmprs,array1,array2,istatus)    ! o
      if(istatus .eq. -1)then ! try increasing array allocations
          deallocate(array1)
          deallocate(array2)
          n_cmprs_max = n_cmprs_max * 2
          goto 800
      elseif(istatus .ne. 1)then
          goto 950
      endif

      l_check_encoding = .false.
      if(l_check_encoding)then ! just for debugging purposes
          call runlength_decode(ngrids,n_cmprs,array1,array2        ! i
     1                         ,data                                ! o
     1                         ,istatus)                            ! o
          if(istatus .ne. 1)then
              write(6,*)
     1            ' 1st decoding test of compressed data unsuccessful'
              goto 980
          else
              write(6,*)
     1            ' 1st decoding test of compressed data successful'
          endif
      endif
c
c **** write out compressed file
c
      lun = 65
      call open_lapsprd_file(lun,i4time,ext,istatus)
      if(istatus .ne. 1)goto 940

      write(lun,*)kdim

      do k = 1,kdim
          write(lun,1)comment(k)
1         format(a)
      enddo ! k

      write(lun,*)n_cmprs

      icheck_sum = 0

      do i = 1,n_cmprs
          write(lun,*)array1(i),array2(i)
          icheck_sum = icheck_sum + array1(i)
      enddo ! i

      close(lun)

      deallocate(array1)
      deallocate(array2)

!     second "internal" checksum test
      if(icheck_sum .ne. ngrids)then
          write(6,*)' 2nd checksum test discrepancy: '
     1             ,icheck_sum, ngrids
          go to 980
      endif
c
c ****  return normally.
c
        istatus=error(1)
999     return
c
c ****  error trapping.
c
920     if (flag .ne. 1) then
          write(6,*) ' write_laps_compressed aborted!'
          write(6,*) ' laps will currently only work on a pressure'
     1              ,' vertical grid'
          write(6,*) ' make sure vertical_grid is set to pressure'
     1              ,' in nest7grid.parms'
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
     1    write (6,*) 'error in/near runlength_encode...write aborted'       
        istatus=error(2)
        goto 999
c
960     if (flag .ne. 1)
     1    write (6,*) 'error writing data to file...write aborted.'
        istatus=error(2)
        goto 999
c
970     if (flag .ne. 1)
     1    write (6,*) 
     1 'error writing header info into file...write aborted.'
        istatus=error(2)
        goto 999
c
980     if (flag .ne. 1)
     1    write (6,*) 
     1   'checksum error with runlength encoding'
        istatus=error(2)
        goto 999
c
990     if (flag .ne. 1)
     1    write (6,*) 
     1 'file already exists for analysis time...write aborted.'
        istatus=error(2)
        goto 999
c
        end


        subroutine runlength_encode(ngrids,n_cmprs_max,data           ! i
     1                             ,n_cmprs,array1,array2,istatus)    ! o

        integer array1(n_cmprs_max)
        real array2(n_cmprs_max)
        real data(ngrids)

!       setup for first point
        n_cmprs = 0
        i_count_same = 1

        do i = 2,ngrids-1

            if(data(i) .eq. data(i-1))then    
                i_count_same = i_count_same + 1
            else
                n_cmprs = n_cmprs + 1
                if(n_cmprs .le. n_cmprs_max)then
                    array1(n_cmprs) = i_count_same
                    array2(n_cmprs) = data(i-1)
                    i_count_same = 1
                else
                    write(6,*)' error, increase n_cmprs_max',n_cmprs_max
                    istatus = -1
                    return
                endif
            endif

        enddo ! i

!       take care of the last point
        i = ngrids

        if(data(i) .eq. data(i-1))then    
            i_count_same = i_count_same + 1
        else
            i_count_same = 1
        endif

        n_cmprs = n_cmprs + 1
        if(n_cmprs .le. n_cmprs_max)then
            array1(n_cmprs) = i_count_same
            array2(n_cmprs) = data(i)
        else
            write(6,*)' error, increase n_cmprs_max',n_cmprs_max
            istatus = -1
            return
        endif

        write(6,*)' end of runlength_encode, number of pts = '
     1           ,n_cmprs,ngrids       
        write(6,*)' compression ratio = ',float(n_cmprs)/float(ngrids)

        istatus = 1
        return
        end
