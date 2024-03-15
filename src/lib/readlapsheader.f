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
      subroutine read_laps_header(i4time,dir,ext,imax,jmax,kmax,
     1                         laps_dom_file,asctime,version,
     1                         model,origin,var,lvl,num_variables,
     1                         var_avail,laps_var_avail,num_levels,
     1                         lvl_avail,lvl_coord,units,
     1                         comment,l_packed_data,istatus)
c
c**********************************************************************
c
c      this file contains the following fortran subroutines:
c            readlapsheader
c
c      the readlapsheader subroutine reads the following fortran
c       subroutines from the readlapsdata.for file:
c            cvt_fname_data
c
c      the readlapsheader subroutine reads the following c subroutines
c      from the readwritelaps.c file:
c            make_c_fname
c            read_cdf_header
c            cdf_retrieve_header
c            itoa
c
c**********************************************************************
c
c      subroutine read_laps_header
c
c      author:    steve albers
c      modified:  to accept netcdf data files          1/93 linda wharton
c                 to accept extended format filenames  2/96  linda wharton
c                 will not make correct filetimes with
c                 extension lga
c
c      creates filename using cvt_fname_data.
c      reads header variables imax, jmax, kmax, laps_dom_file, asctime,
c      version, model, origin, and num_variables.
c
c**********************************************************************
        implicit        none

        integer         fn_length
c
        integer       i4time,         ! input i4time of data
     1          flag,           ! print flag (1 = off)
     1          no_laps_diag,   !if = 0, print diagnostic output
     1          rec,
     1          rcdl,
     1          imax,jmax,kmax, ! output dimensions/#fields from header
     1          lvl(200),      ! output levels from header
     1          i,k,
     1          error(2),
     1                num_variables,
     1                num_levels,
     1          istatus
c
        real          msg_flag
c
        character*150   dir             ! input directory to read data from
        character*31    ext             ! input file name ext (up to 31 chars)
        character*31    ext_i           !input input file name ext
        character*3     var(200)       ! output variables
        character*4     lvl_coord(200) ! output comments
        character*10    units(200)     ! output units
        character*125   comment(200)   ! output comments
        character*9     gtime
        character*5     fcst_hh_mm
        character*91    file_name
        character*4     mark
        character*4     cimax,cjmax,ckmax
        character*24    asctime
        character*18    asct
        character*4     version
        character*5     vern
        character*131   model           !meteorological model in file
        character*132   modl
        character*131   origin          !location where file was created
        character*132   orign
        character*11    laps_dom_file   !name of domain file e.g. nest7grid
        character*12    ldf
        character*1     hmark
        character*4     clvl
        character*4     cstart_rec
        character*3     laps_var_avail(200)
        character*4     lvar_a(200)
        character*19    var_avail(200)
        character*20    var_a(200)
c
        logical         l_packed_data
c
        integer       lvl_avail(200)
c
        data            msg_flag/9.e30/
c
        common          /prt/flag
        common          /laps_diag/no_laps_diag

c
c-------------------------------------------------------------------------------
c
        error(1)=1
        error(2)=0
c
c ****  create file name.
c
        call make_fnam_lp(i4time,gtime,istatus)
        if (istatus .ne. 1) then
                istatus=error(2)
                return
        endif

c fcst_hh_mm: hard wired as a place holder - will be used in filename 
c   only if read_laps_header is called on lga, lgb, fua, fsf, ram, rsf
        fcst_hh_mm = '0000'

        call upcase(ext,ext_i)

        call cvt_fname_data(dir,gtime,fcst_hh_mm,ext_i,file_name,
     1                      fn_length,istatus)

        call read_cdf_header(file_name, fn_length,imax, jmax, kmax,
     1                 num_variables,ldf,asct,vern,modl,orign,
     1                 var_a,lvar_a,num_levels, lvl_avail,
     1                 no_laps_diag,istatus)

        if (istatus .eq. 1) then
          laps_dom_file = ldf
          asctime = asct
          version = vern
          model = modl
          origin = orign
          do i = 1, num_variables
            laps_var_avail(i) = lvar_a(i)
            var_avail(i) = var_a(i)
          enddo

         istatus=error(1)
         call setup_var_lvl(ext_i,num_levels,lvl_avail,
     1                      num_variables,laps_var_avail,
     1                      var,lvl,kmax,istatus)
         if (istatus .ne. 1) then
            write (6,*) 'error in setup_var_lvl'
            istatus=error(2)
         endif
         l_packed_data = .false.
         goto 999
        else
          if (istatus .eq. -2) goto 905

!         if you get to here, istatus = -1 so file not netcdf

c         if (ext(1:3) .eq. 'lc3' .or. ext(1:3) .eq. 'sc3' .or.
c    1    ext(1:3) .eq. 'lps' .or. ext(1:3) .eq. 'sps' .or.
c    1    ext(1:3) .eq. 'lrp' .or. ext(1:3) .eq. 'srp' .or.
c    1    ext(1:3) .eq. 'lty' .or. ext(1:3) .eq. 'sty' .or.
c    1    ext(1:3) .eq. 'lh3' .or. ext(1:3) .eq. 'sh3' .or.
c    1    ext(1:3) .eq. 'lmd' .or. ext(1:3) .eq. 'smd'            )then
            l_packed_data = .false.
c         else
c           l_packed_data = .false.
c         endif
          if(l_packed_data)then
             goto 905
          endif
c
c ****  open file and read first header record.
c
        open(1,file=file_name,status='old',err=940)
        read(1,890,err=895) mark,cimax,cjmax,ckmax,asctime,version
890     format(a4,3a4,a17,a4)
        close(1)
        if (mark(1:1) .ne. '*') then
          goto 895
        endif
c
c ****  decode charater variables.
c
        read(cimax,900)  imax
        read(cjmax,900)  jmax
        read(ckmax,900)  kmax
900     format(i4)
c
        if(kmax .gt. 200)goto 960
c
c ****  open file for direct access.
c
        rcdl=imax+1
        if (rcdl .lt. 40) rcdl=40
        open(1,file=file_name,status='old',
     1     access='direct',
     1     recl=rcdl,err=950)
c
c ****  read header files for each field.
c
        rec=1
        do k=1,kmax
                rec=rec+1
                read(1,rec=rec,err=905) hmark,var(k),clvl,lvl_coord(k),u
     1nits(k),
     1              cstart_rec,comment(k)
                if (hmark .ne. '*') then
                  goto 905
                endif
                read(clvl,900)  lvl(k)
        enddo
906     format(i3)

910     format(1x,a3,i4,' field not found.')

        endif  !if status .eq. 1
c
c ****  return normally.
c
        l_packed_data = .false.
        istatus=error(1)
998     close(1,err=999)
999     return
c
c ****  error trapping.
c
895     if (flag .ne. 1)
     1    write (6,*) 'error reading main header...read aborted.'
        istatus=error(2)
        goto 998
c
905     if (flag .ne. 1)
     1    write (6,*) 'error reading field header...read aborted.'
        istatus=error(2)
        goto 998
c
940     if (flag .ne. 1)
     1    write (6,*) 'error opening file for given i4time...read aborte
     1d.'
        istatus=error(2)
        goto 998
c
950     if (flag .ne. 1)
     1    write (6,*) 'error opening file as direct access...read aborte
     1d.'
        istatus=error(2)
        goto 998
c
960     if (flag .ne. 1)
     1    write (6,*) 'kmax greater than 200...arrays to small'
        istatus=error(2)
        goto 998
c
        end

c########################################################################
      subroutine setup_var_lvl(ext_in,num_levels,lvl_avail,
     1                         num_variables,laps_var_avail,
     1                         var,lvl,kdim,istatus)

c**********************************************************************
c
c      subroutine setup_var_lvl
c
c      author:    linda wharton
c
c      takes data read in by readlapsheader and fills var and lvl
c      so readlapsdata may be called.  this subroutines is used only
c      if the file is in netcdf format.
c
c**********************************************************************
c
      implicit  none

      integer lvl(*),
     1          i, j,
     1          kdim,
     1          num_levels,
     1          num_variables,
     1          error(3),
     1          pos,
     1          istatus                     !output

      integer  lvl_avail(*)

      character*3       laps_var_avail(*)
      character*3       var(*)
      character*31      ext_in              !input input file name ext

c
c-------------------------------------------------------------------------------
c
      error(1)=1
      error(2)=0
      error(3)=-2

      istatus = error(1)

      if (ext_in .eq. 'lmr') then
         var(1) = 'r00'
         var(2) = 'r06'
         var(3) = 'r12'
         lvl(1) = 0
         lvl(2) = 0
         lvl(3) = 0
         pos = 4
      else
         if (ext_in .eq. 'lf1') then
            var(1) = 'h00'
            var(2) = 'h06'
            var(3) = 'h12'
            lvl(1) = 0
            lvl(2) = 0
            lvl(3) = 0
            pos = 4
          else
           if (ext_in .eq. 'lhe') then
              var(1) = 'lhe'
              lvl(1) = 0
              pos = 2
           else
              pos = 1

              do i = 1, num_variables
                 if (laps_var_avail(i)(1:1) .ne. char(0)) then
                    do j = 1, num_levels
                       var(pos) = laps_var_avail(i)
                       lvl(pos) = lvl_avail(j)
                       pos = pos + 1
                    enddo
                 endif
              enddo
            endif
         endif
      endif

      if ((pos - 1) .ne. kdim) istatus = error(2)

      return
      end

c########################################################################
