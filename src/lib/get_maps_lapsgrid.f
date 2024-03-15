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
        subroutine get_modelfg_3d(i4time,var_2d,imax,jmax,kmax
     1                          ,field_3d_laps,istatus)

!       this routine reads the model background from rams or lga/f and
!       interpolates in time to the requested i4time. the time interpolation
!       occurs for lga/f, but the i4time must be a whole hour for a successful
!       read of the .ram files to occur.

!          ~1990  steve albers - original version
!       sep 1994      "        - read in .ram file first, then default over
!                                to lga/f.
!       jun 1995      "        - slight cleanup and added comments
!       jul 1995      "        - now handles lga/f files that come as frequently
!                                as the laps cycle. this is backwards compatable
!                                as it still can do the time interpolation with
!                                three hourly forecast files.
!                     "        - new subroutine 'get_fcst_filename', modularizes
!                                the handling of forecast file naming convention.
!       jan 1996               - now will handle 9 or 13 character versions
!                                of lga filenames
!       feb 1997               - added entry get_modelfg_3d.
!       oct 1998 linda wharton - removed variables never used: a9_filename,
!                                ext_f, field_3d_maps, l_fill, field_3d_maps_1,
!                                field_3d_maps_2
!
!       feb 1999 john smart    - add subroutine get_modelfg_3d_sub which allows
!                                specific model extension.
!
!       jun 2000     "         - add subroutine get_best_fcst so that this code
!                                can be use elsewhere (eg., lga, dprep, laps_sfc, osse)
!          "         "         - incorporate namelist parameter fdda_model_source to
!                                determine the subdirectory in which fdda model bckgnd
!                                exist.  remove 'ram' as fdda extension.
!       jan 2001     "         - obtain 2d fields from fsf or lgb when kmax = 1.
!                                new jacket routine get_modelfg_2d.

        include 'bgdata.inc'

        real field_3d_laps(imax,jmax,kmax)       ! output array

        logical      lgab

        character*3  var_2d
        character*3  c_bkgd_ext(2)
        character*9  a9_time
        character*31 ext_a(maxbgmodels)
        character*31 subdir(maxbgmodels)
        character*9 fdda_model_source(maxbgmodels)
        character*5 bgmodelnames(maxbgmodels)
 
        integer     nbgm

!       ****************** model selection section *********************

!       restrictions:

!       1) no time interpolation is performed
!       2) fdda/lga file must be available valid at i4time/i4time_needed

        call bgmodel_name(maxbgmodels,nbgm,bgmodelnames,istatus)
!       do i=1,nbgm
!        print*,'bg model derived from from cmodel = ',bgmodelnames(i)
!       enddo

        if(kmax.eq.1)then
           c_bkgd_ext(1)='fsf'
           c_bkgd_ext(2)='lgb'
        else
           c_bkgd_ext(1)='fua'
           c_bkgd_ext(2)='lga'
        endif

        i4time_needed = i4time

        call get_fdda_model_source(fdda_model_source,n_fdda_models
     +,istatus)
        if(istatus .ne. 1)then
           print*,'error getting fdda_model_source'
           return
        endif

        lgab=.false.
        if(n_fdda_models .eq. 0)then
           n_fdda_models = 1
c          subdir(n_fdda_models)=bgmodelnames(1:3) #replaces line below (needs mod) when lga subdirs is activated
           subdir(n_fdda_models)=c_bkgd_ext(2)
           ext_a(n_fdda_models) =c_bkgd_ext(2)
           lgab=.true.

c if lga is the first in the list then we'll force this to be the
c only possible background returned from get_modelfg_2/3d..
c removed functionality 12-06-01 (jrs)
c       elseif(fdda_model_source(1).eq.'lga')then
c          n_fdda_models = 1
c          subdir(n_fdda_models)=c_bkgd_ext(2)
c          ext_a(n_fdda_models) =c_bkgd_ext(2)
c          lgab=.true.

        else
           do i=1,n_fdda_models
              subdir(i)=fdda_model_source(i)
              if(subdir(i).eq.'lga')then        !.or.subdir(i).eq.'lgb')then
                 ext_a(i) =c_bkgd_ext(2)
                 lgab=.true.
              else
                 ext_a(i) =c_bkgd_ext(1)
              endif
           enddo
        endif
c this part adds lga/b to the list since it isn't already  part of it.
        if(.not.lgab)then
           if(n_fdda_models.lt.maxbgmodels)then
              n_fdda_models = n_fdda_models + 1
c             subdir(n_fdda_models)=bgmodelnames(1:3) #replaces line below (needs mod) when lga subdirs is activated
              subdir(n_fdda_models)=c_bkgd_ext(2)
              ext_a(n_fdda_models) =c_bkgd_ext(2)
           else
              print*,'*** warning *** '
              print*,'cannot add lga/b to model background list in'
     .,'get_modelfg_3d, n_fdda_models > maxbgmodels'
           endif
        endif
           
        do isource = 1,n_fdda_models

           call make_fnam_lp(i4time_needed,a9_time,istatus)

           write(6,*)
           write(6,*)' searching for model background valid at: '
     1                        ,a9_time,' ',ext_a(isource)(1:6),var_2d

           call get_modelfg_3d_sub(i4time_needed,var_2d,subdir(isource),
     1                      ext_a(isource),imax,jmax,kmax,field_3d_laps,
     1                      istatus)

           if(istatus.eq.1)then
              print*,'file obtained in get_modelfg_3d_sub - return'
              return
           endif

        enddo

        return
        end

!***************** start new section ******************************************
! start subroutine

        subroutine get_modelfg_3d_sub(i4time_needed,var_2d,subdir,ext_a
     1                         ,imax,jmax,kmax,field_3d_laps,istatus)
!
!
        real field_3d_laps(imax,jmax,kmax)       ! output array

        character*(*) var_2d
        character*(*) subdir
        character*9 a9_time
        character*14 a_filename

        character*125 comment_2d
        character*10 units_2d
        character*31 ext,ext_a
        character*150  directory
        character*255 c_filespec

        integer max_files
        parameter (max_files = 20000)
        character c_fnames(max_files)*180


        call get_directory(ext_a,directory,lend)
c here -> reverse the directory order from "model/fua" to "fua/model"
        if(ext_a.ne.'lga'.and.ext_a.ne.'lgb')then
           call s_len(subdir,ld)
           if(ld .le. 0)then
               write(6,*)' subdir has zero length in get_modelfg_3d_sub'
               istatus = 0
               return
           endif
           directory = directory(1:lend)//subdir(1:ld)//'/'
           lend=lend+ld+1
        endif
        c_filespec=directory(1:lend)

        c_filespec=c_filespec(1:lend)//'*.'//ext_a(1:3)

!       obtain list of analysis/forecast filenames
        call get_file_names(c_filespec,
     1                      i_nbr_files_ret,
     1                      c_fnames,
     1                      max_files,
     1                      istatus )

!       determine which file having the proper valid time has the
!       most recent initialization time.

        call get_best_fcst(max_files,i4time_needed,i_nbr_files_ret
     1,c_fnames,i_best_file)

        if(i_best_file .gt. 0)then ! file for this ext exists with proper
           i = i_best_file
           call get_directory_length(c_fnames(i),lend)
           call get_time_length(c_fnames(i),lenf)                               ! valid time.
           ext = ext_a
           a_filename = c_fnames(i)(lend+1:lenf)

           write(6,*)' found file for: ',c_fnames(i)(lend+1:lenf)
     1                                       ,' ',ext(1:6),var_2d

           call get_fcst_times(a_filename,i4_initial,i4_valid
     1                        ,i4_fn)

           if(lenf - lend .eq. 9)then

c             call get_laps_3d(i4_fn,imax,jmax,kmax,ext,var_2d
c    1                ,units_2d,comment_2d,field_3d_laps,istatus)
c
c
c new
c
              call get_3d_dir_time(directory,i4_fn
     1                      ,ext,var_2d,units_2d,comment_2d
     1                      ,imax,jmax,kmax,field_3d_laps,istatus)

           elseif(lenf - lend .eq. 13 .or. lenf - lend .eq. 14)then       

               if(kmax.gt.1)then
c
c new: changed variable name "ext" to "directory".
                  call get_lapsdata_3d(i4_initial,i4_valid,imax
     1                 ,jmax,kmax,directory,var_2d
     1                 ,units_2d,comment_2d,field_3d_laps,istatus)

                else

                  call get_lapsdata_2d(i4_initial,i4_valid,directory
     1              ,var_2d,units_2d,comment_2d,imax,jmax
     1              ,field_3d_laps(1,1,1),istatus)

c                 call get_2dgrid_dname(directory
c    1         ,i4_fn,0,i4time_nearest,ext,var_2d,units_2d
c    1         ,comment_2d,imax,jmax,field_3d_laps(1,1,1),0,istatus)


                endif

           else
               write(6,*)' error, illegal length of bckgd filename'
     1                 ,lend,lenf
               istatus = 0

           endif

           if(istatus .ne. 1)then
              write(6,*)'get_modelfg_3d_sub: warning - could not read'
     1                 ,' model file'
           elseif(kmax.gt.1)then ! istatus = 1
               call qc_field_3d(var_2d,field_3d_laps
     1                         ,imax,jmax,kmax,istatus)            
           else
               call qc_field_2d(var_2d,field_3d_laps(1,1,1)
     1                         ,imax,jmax,istatus)
           endif

       else ! i_best_file = 0
           write(6,*)' no file with proper valid time'
           istatus = 0

       endif

       if(istatus .eq. 1)then
          call s_len(comment_2d,lenc)
          lenc = max(lenc,1)
          write(6,*)' successfully obtained: '
     1                        ,c_fnames(i)(lend+1:lenf),' '
     1                        ,ext(1:6),var_2d,comment_2d(1:lenc)
          write(6,*)' exiting get_modelfg_3d_sub'
          write(6,*)
          return
       else
          write(6,*)' no good ',ext_a(1:6),' files available.'          
       endif

       write(6,*)
       write(6,*)' no good files: exiting get_modelfg_3d_sub'

       istatus = 0
       return
       end


        subroutine get_fcst_times(a_filename,i4_initial,i4_valid,
     1                            i4_filename)

!       a_filename    (character)                                         i
!       i4_initial                                                        o
!       i4_valid                                                          o
!       i4_filename   (this is what you might pass into read_laps_data)   o

!       this routine deals with either c9 or c13 forecast file naming
!       convention

        character*(*) a_filename
        character*2 c2_hr, c2_mn
        character*3 c3_hr

        call s_len(a_filename,i_len)

        if(i_len .eq. 9)then
            call cv_asc_i4time(a_filename,i4_filename)
            i4_initial = (i4_filename/3600) * 3600
            i4_fcst = (i4_filename - i4_initial) * 60
            i4_valid = i4_initial + i4_fcst

        elseif(i_len .eq. 13)then
            call cv_asc_i4time(a_filename(1:9),i4_filename)
            i4_initial = i4_filename
            c2_hr = a_filename(10:11)
            c2_mn = a_filename(12:13)
            read(c2_hr,1)i4_fcst_hr
            read(c2_mn,1)i4_fcst_mn
 1          format(i2)
            i4_fcst =  i4_fcst_hr * 3600 + i4_fcst_mn * 60
            i4_valid = i4_initial + i4_fcst

        elseif(i_len .eq. 14)then
            call cv_asc_i4time(a_filename(1:9),i4_filename)
            i4_initial = i4_filename
            c3_hr = a_filename(10:12)
            c2_mn = a_filename(13:14)
            read(c3_hr,3)i4_fcst_hr
            read(c2_mn,1)i4_fcst_mn
 3          format(i3)
            i4_fcst =  i4_fcst_hr * 3600 + i4_fcst_mn * 60
            i4_valid = i4_initial + i4_fcst

        else
            write(6,*)' get_fcst_time: illegal value of i_len',i_len

        endif

        return
        end

c --------------------------------------------------------------------------
        subroutine get_best_fcst(maxfiles,i4time_needed
     1,i_nbr_files,c_fnames,i_best_file)

c
c determine the best file that matches the i4time_needed input
c j. smart 6-20-00:  pulled this section of software out of
c                    get_modelfg_3d_sub for use elsewhere in laps.
c
        implicit  none

        integer i4time_needed
        integer i_best_file,i
        integer i_nbr_files
        integer i4_fcst_time_min
        integer i4_valid,i4_fcst_time,i4_fn
        integer i4_initial
        integer lend,lenf
        integer maxfiles

        character*(*)  c_fnames(maxfiles)
        
        i_best_file = 0
        i4_fcst_time_min = 9999999

        do i=1,i_nbr_files
            call get_directory_length(c_fnames(i),lend)
            call get_time_length(c_fnames(i),lenf)
            if(lenf.le.0)then
               call get_fname_length(c_fnames(i),lenf)
               lenf=lend+lenf
            endif
            if(lenf.eq.lend)then
               print*,'no filenames in c_fnames array: get_best_fcst'
               return
            endif

            call get_fcst_times(c_fnames(i)(lend+1:lenf)
     1                 ,i4_initial,i4_valid,i4_fn)
            if(i4_valid .eq. i4time_needed)then
               i4_fcst_time = i4_valid - i4_initial

               if(i4_fcst_time .lt. i4_fcst_time_min)then

                  i4_fcst_time = i4_fcst_time_min
                  i_best_file = i

               endif ! smallest forecast time?
            endif ! correct valid time
        enddo ! i

        return
        end
c
c-------------------------------------------------------------------------
        subroutine get_modelfg_2d(i4time,var_2d,imax,jmax
     1                          ,field_2d_laps,istatus)

c
c routine requires kmax = 1
c
        include 'bgdata.inc'

        integer imax,jmax

        real field_2d_laps(imax,jmax)       ! output array

        character*3  var_2d
        character*9  a9_time
        character*9  fdda_model_source(maxbgmodels)
        character*31 ext_a(maxbgmodels)
        character*31 subdir(maxbgmodels)

        call get_modelfg_3d(i4time,var_2d,imax,jmax,1
     1                          ,field_2d_laps,istatus)

        if(istatus.ne.1)then
           print*,'warning: no 2d field from get_modelfg_2d'
           return
        endif

        return
        end
