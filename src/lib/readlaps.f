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
       subroutine read_laps(reftime,valtime,dir,ext,iimax,jjmax,
     1                      kkmax,kdim,var_req,lvl_req,lvl_coord_req,
     1                      units_req,comment_req,data,istatus)
c**********************************************************************
c
c      this file contains the following fortran subroutines:
c            read_laps
c
c      the read_laps routine reads the following fortran
c      routines from the readlapsdata.f file:
c            make_fcst_time
c            cvt_fname_v3
c
c      the read_laps routine calls the following c function
c      from the rwl_v3.c file:
c            read_cdf_v3
c
c**********************************************************************
c
c      author:    linda wharton
c      modified:  to accept netcdf ver. 2 data files  1/93 linda wharton
c                 to accept valtime & reftime  2/96 linda wharton
c                 to accept netcdf ver. 3 data files  9/97 linda wharton 
c
cdoc   reads data requested by arrays var_req and lvl_req for the
cdoc   i4time, dir and ext specified.  returns lvl_coord-req,
cdoc   units_req, comment_req, data and istatus
cdoc   reftime is the time of the model run.
cdoc   valtime is the valid time of the data.
cdoc   for analysis, valtime = reftime.
c
c**********************************************************************
c
      implicit  none
c
      integer       reftime,             !input i4time of model run
     1                valtime,             !input i4time data is valid
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
      call make_fnam_lp(reftime,gtime,istatus)
      if (istatus .ne. 1) then
        write (6,*) 
     1'error converting reftime to file name...read aborted.'
        istatus=error(2)
        return
      endif

      call make_fcst_time(valtime,reftime,fcst_hh_mm,istatus)

      call s_len(ext, ext_len)

      call cvt_fname_v3(dir,gtime,fcst_hh_mm,ext,ext_len,
     1                  file_name,fn_length,istatus)
      if (istatus .eq. error(2)) goto 930


      i_reftime = reftime - 315619200
      i_valtime = valtime - 315619200
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
     1    write (6,*) 'error opening netcdf file...read aborted...'
     1               ,file_name
        istatus=error(2)
        goto 999
c
970     if (flag .ne. 1)
     1    write (6,*) 'error retrieving data...read aborted.'
        istatus=error(2)
        goto 999
c
980     if (flag .ne. 1)
     1    write (6,*) 'error in array dimensioning...read aborted.'
        istatus=error(2)
        goto 999
c
990     if (flag .ne. 1)
     1    write (6,*) 'error in version, not a valid laps file...
     1read aborted.'
        istatus=error(2)
        goto 999
c
992     continue
        istatus=error(2)
        goto 999
c
        end

c########################################################################

