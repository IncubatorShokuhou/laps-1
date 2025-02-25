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

        subroutine get_directory(ext_in,directory_out,len_dir)

!       1993         s. albers

!       this routine is a multipurpose routine that accepts an extension.
!       if the extension is a domain name, the location of the static file
!       directory is returned.

!       added wrfsi as a possible domain name (namelist) in addition to nes7grid.
!       thus, made the data root a bit more generic. retain functionality for
!       calls independent of get_laps_config; ie.,  still gets the env or command
!       line input
!       2000         j. smart

!       if the extension is a product file, the domain is tested, the extension
!       is modified for the domain if necessary, and the directory containing
!       the products is returned. this routine is designed for a unix
!       environment.

        character*(*) ext_in         ! input
        character*80 ext             ! local
        character*(*) directory_out  ! output
        integer len_dir            ! output

        include 'grid_fname.cmn'

        character*201 directory

        call s_len(grid_fnam_common,len_dir)

cc        len_dir = index(grid_fnam_common,'/',.true.)
        if(len_dir.eq.0) then
           call get_config(istatus)
           if(istatus.ne.1)then
              print*,'error: get_config failed'
              return
           endif
        endif

!       test if the extension is the domain name. if so, then return the
!       directory containing the static files.

        call s_len(grid_fnam_common,len_grid_fnam)
        call s_len(ext_in,len_ext_in)
        call s_len(generic_data_root,len_root)

        if(len_root .gt. 200) then
          write(6,*)'get_directory error: the env var for the data root'
     1             ,' is too long; shorten to <200'
          stop
        endif 
c        print *, generic_data_root(1:len_root),len_root

        if(generic_data_root(len_root:len_root).ne.'/') then
           len_root=len_root+1
           generic_data_root(len_root:len_root)='/'
        endif
c
c both laps and wrfsi have static files in 'static' at the moment.
        if(len_ext_in .eq. len_grid_fnam)then
          if(len_ext_in .ne. 0)then
            if(ext_in(1:len_ext_in).eq.grid_fnam_common(1:len_grid_fnam)
     1)then
              if(ext_in(1:len_ext_in) .eq. 'nest7grid'
     1 .or.      ext_in(1:len_ext_in) .eq. 'wrfsi')then
                directory = generic_data_root(1:len_root)//'static/'
                goto 999
              endif
            endif ! ext_in .eq. grid_fnam_common
          endif ! len_ext .ne. 0
        endif ! lens are = ?

        call downcase(ext_in,ext)

!       test if the extension is for the thermo directory.
!       if so, then return the directory containing the static thermo files.
!       added direct test for static dir
        if(ext(1:3) .eq. 'dat' .or. ext(1:6) .eq. 'static')then
           directory = generic_data_root(1:len_root)//'static/'
           goto 999
        endif

        if(ext(1:3) .eq. 'etc' .or. ext(1:4).eq.'time')then
c           directory = generic_data_root(1:len_root)//'etc/'
           directory = generic_data_root(1:len_root)//'time/'
           goto 999
        endif

        if(ext(1:3) .eq. 'cdl')then
           directory = generic_data_root(1:len_root)//'cdl/'
           goto 999
        endif

        if(ext(1:3) .eq. 'log')then
           directory = generic_data_root(1:len_root)//'log/'
           goto 999
        endif

        if(ext(1:4) .eq. 'root' .or. ext(1:8) .eq. 'dataroot')then
           directory = generic_data_root(1:len_root)//'/'
           goto 999
        endif

!       in this section, we assume the extension points to a product file
        call s_len(ext,len_ext)


        directory = generic_data_root(1:len_root)
     +       //'lapsprd/'//ext(1:len_ext)//'/'

!       get length of directory

 999    call s_len(directory,len_dir)

!       test length of input directory
        if(len_dir .gt. len(directory_out))then
            write(6,*)' get_directory error: input directory is too'     
            write(6,*)' short, lengthen from ',len(directory_out)
     1               ,' to ',len_dir     
            stop

        elseif(len_dir .eq. 201)then
            write(6,*)' get_directory error: laps_data_root is too long'       
     1               ,' total directory length >= ',len_dir
            stop ! consider increasing max string lengths in this routine

        else ! len_dir fits all constraints
            directory_out = directory(1:len_dir)

        endif ! len_dir
        return
        end
