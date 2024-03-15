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
       subroutine read_laps_data(i4time,dir,ext,iimax,jjmax,kkmax,kdim,
     1                    var_req,lvl_req,lvl_coord_req,
     1                    units_req,comment_req,data,istatus)
c**********************************************************************
c
c      this file contains the following fortran subroutines:
c            read_laps_data
c            cvt_fname_v3
c            make_fcst_time
c
c      the read_laps_data subroutine calls the following c function
c      from the rwl_v3.c file:
c            read_cdf_v3
c
c**********************************************************************
c
c      subroutine read_laps_data
c
c      author:    linda wharton
c      modified:  to accept netcdf ver. 2 data files  1/93 linda wharton
c                 to accept netcdf ver. 3 data files...capablility to 
c                   read binary laps files separated out into 
c                   read_old_laps.f  9/97 linda wharton
c
c      reads data requested by arrays var_req and lvl_req for the
c      i4time, dir and ext specified.  returns lvl_coord-req,
c      units_req, comment_req, data and istatus
c
c**********************************************************************
c
      implicit  none
c
      integer       i4time,              !input i4time of data
     1                iimax,jjmax,kkmax,   !input # cols, # rows, # fields
     1                kdim,                !input k dimension of data array
     1                lvl_req(kdim),       !input requested levels
     1                istatus              !output
      character*(*)   dir                  !input directory to read data from
      character*(*)   ext                  !input file name ext 
      character*(*)   var_req(kdim)        !input 3 letter id of requested fields
      character*(*)   lvl_coord_req(kdim)  !output vertical coordinate of fields
      character*(*)   units_req(kdim)      !output units of requested fields
      character*(*)   comment_req(kdim)    !output comments for requested fields
      real        data(iimax,jjmax,kdim) !output data
c
      integer fn_length,
     1          i_reftime,              !unix time of data
     1          i_valtime,              !unix time of data
     1          flag,                   !print flag (1 = off)
     1          error(3),
     1          called_from,
     1          var_len,
     1          comm_len,
     1          ext_len,
     1          lvl_coord_len,
     1          units_len
  
c
      character*5       fcst_hh_mm
      character*9       gtime
      character*150     file_name
c
      common            /prt/flag
c
c-------------------------------------------------------------------------------
c
      error(1)=1
      error(2)=0
      error(3)=-2
c
c ****  create file name.
c
      call make_fnam_lp(i4time,gtime,istatus)
      if (istatus .ne. 1) then
        write (6,*) 
     1'error converting i4time to file name...read aborted.'
        istatus=error(2)
        return
      endif

      call s_len(ext, ext_len)

c hard wired as a place holder - will be used in filename only if read_laps_data
c   is called on lga, lgb, fua, fsf, ram, rsf
c to fix this, call read_laps instead
      fcst_hh_mm = '0000'

      call cvt_fname_v3(dir,gtime,fcst_hh_mm,ext,ext_len,
     1                  file_name,fn_length,istatus)
      if (istatus .eq. error(2)) goto 930
c
c **** get actual reftime from gtime...
c
      i_reftime = i4time - 315619200
      i_valtime = i_reftime

      called_from = 0   !called from fortran

      var_len = len(var_req(1))
      comm_len = len(comment_req(1))
      lvl_coord_len = len(lvl_coord_req(1))
      units_len = len(units_req(1))

      call read_cdf_v3 (file_name, ext, var_req, comment_req, 
     1                  lvl_coord_req, units_req, var_len, comm_len, 
     1                  fn_length, ext_len, lvl_coord_len, units_len, 
     1                  i_reftime, i_valtime, iimax, jjmax, kkmax,  
     1                  kdim, lvl_req, data, called_from, istatus)


      if (istatus .ge. 0) then   !return from read with no errors
                                 !convert byte data to characters

         if (istatus .gt. 0) then
            istatus=error(3)
         else
            istatus=error(1)
         endif
         goto 999
      endif

      if (istatus .eq. -5)  goto 992  !file not there
      if (istatus .eq. -4)  goto 990  !error in version 
      if (istatus .eq. -3)  goto 980  !error in dimensioning arrays
      if (istatus .eq. -2)  goto 970  !error retrieving data
      if (istatus .eq. -1)  goto 950  !error opening file as netcdf
c
c ****  return normally.
c
        istatus=error(1)
999     return
c
c ****  error trapping.
c
930     if (flag .ne. 1)
     1    write (6,*) 'file_name variable too short...read aborted.'
        istatus=error(2)
        goto 999
c
950     if (flag .ne. 1)
     1    write (6,*) 'error opening netcdf file...read aborted.'
        istatus=error(2)
        goto 999
c
970     if (flag .ne. 1)
     1    write (6,*) 'error retrieving data...read aborted.'
        istatus=error(2)
        goto 999
c
980     if (flag .ne. 1) then
          write (6,*) 'error in array dimensioning...read aborted.'
          write (6,*) 'iimax/jjmax/kkmax/kdim = ',iimax,jjmax,kkmax,kdim
        endif
        istatus=error(2)
        goto 999
c
990     if (flag .ne. 1)
     1    write (6,*) 'error in version, not a valid laps file... '
     1                ,'read aborted.'
        istatus=error(2)
        goto 999
c
992     continue
        istatus=error(2)
        goto 999
c
        end

c########################################################################
      subroutine cvt_fname_v3(dir,gtime,fhh,ext,ext_len,file_name,
     1                        fn_length,istatus)
c**********************************************************************
c
c      subroutine cvt_fname_v3
c
c      author:    linda wharton 1/93
c
c      inputed dir, gtime and ext are converted to ascii values in a byte
c      array.  b_fname is created in c function make_c_fname.  file_name
c      is also created.  b_ext is returned for use in other functions.
c
c      modified:  linda wharton 4/94
c      changed to remove references to byte data type, file_name is
c      created from dir, gtime and ext.  no ascii conversions done.
c
c**********************************************************************

      implicit  none

      character*(*)     dir          !directory to read data from
      character*(*)     gtime
      character*(*)     ext          !file name ext 
      character*(*)     file_name
      character*5       fhh

      integer         fn_length,
     1                  ext_len,fhh_len,
     1                  istatus

      integer         end_dir, end_ext, error(2)

c#ifdefined nodynamic
c      character*31  ext_dn
c#else
      character*(10) ext_dn
c#endif

      call downcase(ext,ext_dn)

      error(1)=1
      error(2)=0

c **** convert string data so it can be used by c programs
c ******  find end of dir 
c
      call s_len(dir,end_dir)
c
c ******  find end of ext_dn
c
      call s_len(ext_dn,end_ext)
c
c ******  find end of file_name
c
      fn_length = len(file_name)
c
c ******  find end of fhh
c
      call s_len(fhh,fhh_len)
c
c ****  make fortran file_name
c

      if (end_dir+end_ext+fhh_len+10 .gt. fn_length) then
        write (6,*) 'length of dir+file-name exceeds file_name length.'
        istatus = error(2)
        goto 999
      else
        if (ext_dn(1:2) .eq. 'lg' .or.
     +     ext_dn(1:3).eq.'fua' .or.
     +     ext_dn(1:4).eq.'cont' .or.
     +     ext_dn(1:3).eq.'fsf' .or.
     +     ext_dn(1:3).eq.'ram' .or.
     +     ext_dn(1:3).eq.'rsf') then
           file_name=dir(1:end_dir)//gtime//fhh(1:fhh_len)//'.'//
     +ext_dn(1:end_ext)
        else
          file_name = dir(1:end_dir)//gtime//'.'//ext_dn(1:end_ext)
        endif

        call s_len(file_name,fn_length)
      endif

c
c ****  return normally.
c
        istatus=error(1)
999     return
        end
c########################################################################
      subroutine make_fcst_time(valtime,reftime,fcst_hh_mm,istatus)

      implicit none

      integer       valtime, reftime, istatus
      integer       fcst_hr, fcst_min, fcst_min_sec, fcst_sec
      integer       error(3), extended, interim_fh
      character*5   fcst_hh_mm
      character*1   h1, h2, h3, m1, m2

      error(1)=1
      error(2)=0

      fcst_sec = valtime - reftime

c see if fcst_sec > 356400 seconds (99 hours)
      if (fcst_sec .gt. 356400) then
        extended = 1  ! fcst_hh_mm format is hhhmm if hh > 99
      else
        extended = 0  ! fcst_hh_mm format is hhmm if hh <= 99
      endif

      if (fcst_sec .eq. 0) then
        fcst_hh_mm = '0000'
        goto 998
      endif

c ****  fcst_sec > 0 .... create fcst_hh_mm

      fcst_min_sec = mod(fcst_sec,3600)
      fcst_hr = (fcst_sec - fcst_min_sec) / 3600
      fcst_min = fcst_min_sec / 60

c     fcst_hr can be between 0 and 999
c     fcst_min can be between 0 and 59

      if ((fcst_hr .lt. 0) .or. (fcst_hr .gt. 999)) then
       write(6,*) ' forecast hour cannot exceed 999: ',fcst_hr
        goto 997
      else
        if ((fcst_min .lt. 0) .or. (fcst_min .gt. 59)) then
          write(6,*) ' forecast minute cannot exceed 59: ',fcst_min
          goto 997
        endif
      endif

      if (extended .eq. 1) then  !fcst_hr > 99, so 3 digit hr is valid
        interim_fh = fcst_hr - mod(fcst_hr,100) 
        h1 = char((interim_fh/100) + 48)
        interim_fh = fcst_hr - interim_fh 
        h2 = char(((interim_fh - mod(interim_fh,10)) / 10) + 48)
        h3 = char(mod(interim_fh,10) + 48)
      else ! fcst_hr <= 99, so 2 digit hr is valid
        if (fcst_hr .lt. 10) then
          h1 = '0'
          h2 = char(fcst_hr + 48)
        else
          h1 = char(((fcst_hr - mod(fcst_hr,10)) / 10) + 48)
          h2 = char(mod(fcst_hr,10) + 48)
        endif
      endif

      if (fcst_min .lt. 10) then
        m1 = '0'
        m2 = char(fcst_min + 48)
      else
        m1 = char(((fcst_min - mod(fcst_min,10)) / 10) + 48)
        m2 = char(mod(fcst_min,10) + 48)
      endif

      if (extended .eq. 1) then
        fcst_hh_mm = h1//h2//h3//m1//m2
      else
        fcst_hh_mm = h1//h2//m1//m2
      endif
      goto 998
c
c ****  error return.
c

997   istatus=error(2)
      goto 999

c
c ****  return normally.
c

998   istatus=error(1)

999   return
      end
c########################################################################
