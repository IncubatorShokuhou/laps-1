cdis
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
c       prof_cdf_close.for              michael barth           11-aug-1993
c
c       this subroutine is used to close a netcdf file of profiler data.
c
c input:
c
c       cdfid           integer:  the netcdf id of the file to be closed.  it
c                       must previously have been opened by prof_cdf_open.
c
c       status          integer:  returned status:  0 is good, -1 and +n are
c                       errors returned by netcdf.  the definitions of these can
c                       be found in netcdf.inc.  in addition, these error
c                       returns are from prof_cdf_close:
c
c                       -5      netcdf file not open.
c
c
c modifications:
c       28-nov-1995/lab put prof_cdf_common in it's own include file.
c       07-aug-2000/lsw modified to netcdf version 3.4

c
        subroutine prof_cdf_close(cdfid,status)
c
        implicit none
c
        integer cdfid,status
c
c       internal prof_cdf information:
c
c       open_list       logical flags indicating which of the four file slots
c                       are in use.
c
c       cdfid_list      netcdf id's for the close files.
c
c       stanam_list     station names for the records in the close files.
c
c       wmoid_list      wmo id's for the records in the close files.
c
c       nrec_list       number of records (stations) in each close file.
c
c       error_flag      error handling flag.
c
        integer max_profilers
        parameter (max_profilers = 200)
c
        include 'prof_cdf_common.inc'
        include 'netcdf.inc'
c
        integer i
c
        do i = 1, 4
c
c       find the user's cdfid in our internal list.
c
           if(cdfid_list(i).eq.cdfid.and.open_list(i))then
c
                status = nf_close(cdfid_list(i))
c
                open_list(i) = .false.                  ! indicate internally.
                go to 900
           endif
        enddo
c
        status = -5                                     ! file not found.
c
c       use the non-default error handling for prof_cdf errors.
c
        if(error_flag.ne.0)then
                write(*,*)'prof_cdf_close:  netcdf file not open.'
                if(error_flag.eq.2)stop
        endif
c
900     return
        end

c   prof_cdf_open.for           michael barth           11-aug-1993
c
c       this subroutine is used to open a netcdf file holding profiler data.
c       it will also read in the station names and wmo id's for later use when 
c       read routines attempt to locate data.
c
c       a user of the prof_cdf routines is allowed to have up to 4 open
c       netcdf files at a time.
c
c input:
c
c       file            character string giving the name (including any path
c                       information) of the netcdf file to be opened.
c
c output:
c
c       cdfid           integer:  the netcdf id for the open file.
c
c       status          integer:  returned status:  0 is good, -1 and +n are
c                       errors returned by netcdf.  the definitions of these can
c                       be found in netcdf.inc.  in addition, these error
c                       returns are from prof_cdf_open:
c
c                       -2      attempt to open more than 4 files.
c
c modifications:
c
c       28-nov-1995/lab put prof_cdf_common in it's own include file.
c
c       07-mar-1996/mfb added the use of recdim_list to determine if a
c                       one-dimensional variable's dimension is the record
c                       dimension or not -- the hyperslab count
c                       depends on this distinction. ***** note that the
c                       recdim_list calculation is currently set to
c                       compensate for a bug in the old vms netcdf library
c                       in use.  see below where recdim_list is calculated
c                       when porting this routine to unix (or in changing
c                       netcdf library versions).
c
c       21-mar-1996/mfb finally went out and got the latest version of
c                       netcdf (2.4.1) and installed on my cluster.  this
c                       eliminated the need to compensate for the netcdf
c                       bug described in the 3/7/96 modification note.
c                       also removed declarations that duplicate those in new
c                       netcdf.inc.
c
c       01-may-1996/mfb went back to the old version of the netcdf library
c                       due to performance problems with the new one.  thus
c                       the bug in ncinq mentioned in the 3/7/96 mod note
c                       again has to be accommodated.  see the code below
c                       for more explanation.
c
c       26-feb-1997/mfb fixed a bug that was exposed if an existing file
c                       is opened, but it has no records in it.  this should
c                       be legal, and now that i've tried it and found the
c                       bug, it is legal...
c
c       28-feb-1997/lab removed a forgotten write statement
c
c       17-mar-1997/mfb removed the recdim fix for vms's old netcdf library.
c                       this should fix problems currently encountered with
c                       single dimensional variables.
c
c       16-apr-1997/mfb removed duplicate definitions of netcdf functions
c                       that are contained in netcdf.inc.
c
c       27-may-1997/mfb renamed "ncopts" common and variable so it doesn't
c                       conflict with symbols used in c netcdf interface.
c
c       07-aug-2000/lsw modified to netcdf version 3.4
c
        subroutine prof_cdf_open(file,cdfid,status)
c
        implicit none
c
        character*(*) file
        integer cdfid,status
c
c       internal prof_cdf information:
c
c       open_list       logical flags indicating which of the four file slots
c                       are in use.
c
c       cdfid_list      netcdf id's for the open files.
c
c       staname_list    station names for the records in the open files.
c
c       wmoid_list      wmo id's for the records in the open files.
c
c       nrec_list       number of records (stations) in each open file.
c
c       recdim_list     record dimension id's for each file.
c
c       error_flag      error handling flag.
c
        integer max_profilers
        parameter (max_profilers = 200)
c
        include 'prof_cdf_common.inc'
        include 'netcdf.inc'
        character*128 dimname                   ! must match netcdf.inc's
c                                               ! maxncnam parameter.
        integer ncopts_val
        common/ncopts_cmn/ncopts_val            ! netcdf error handling flag.
c
        integer i,j,varid,start(2),count(2)
        integer ndims,nvars,ngatts,recdim

        logical l_exist

        write(6,*)' subroutine prof_cdf_open'

        inquire(file=file,exist=l_exist)
        if(.not. l_exist)then
            write(6,*)' warning in prof_cdf_open, file does not exist: '
     1               ,file 
            status=-7
            return
        endif ! l_exist
c
c       unless non-default error handling has been set by a call to
c       prof_cdf_set_error set the default error handling.
c
c
        do i = 1, 4
c
c       find an open slot.
c
           if(.not.open_list(i))then
c
c       open the file.
c
              status = nf_open(file,nf_nowrite,cdfid)
              if(status.eq.nf_noerr)then
c
c       clear out the data structures for this slot.
c
                 do j = 1, max_profilers
                    staname_list(j,i) = '      '
                    wmoid_list(j,i) = 0
                 enddo
                 nrec_list(i) = 0
c
c       get the number of records (stations) in the file.
c
                 status = nf_inq_dimid(cdfid, 'recnum',varid)
                 if (status.eq.nf_noerr)then
                   status = nf_inq_dimlen(cdfid,varid,nrec_list(i))

                 
                   if(nrec_list(i) .gt. max_profilers)then
                     write(6,*)' too many profilers in prof_cdf_open '       
     1                        ,nrec_list(i),max_profilers
                     status=-7
                     return
                   endif
                 endif

                 if(status.eq.nf_noerr.and.nrec_list(i).gt.0)then
c
c       get the station names.  (have to do this in a loop as i can't get
c       ncvgtc to work on one call for the whole enchilada.)
c
                    status = nf_inq_varid(cdfid,'staname',varid)

                    start(1) = 1
                    count(1) = 6 ! 6 characters
                    count(2) = 1 ! 1 record at a time
                    do j = 1, nrec_list(i)
                       start(2) = j
                       status = nf_get_vara_text(cdfid,varid,start,
     $                                           count,
     $                                           staname_list(j,i))
                    enddo
                    if(status.eq.0)then
c
c       get the wmo id's.
c
                       start(1) = 1 ! start at the beginning
                       start(2) = 1
                       count(1) = nrec_list(i) ! and get all records.
                       status = nf_inq_varid(cdfid,'wmostanum',varid)
                       status = nf_get_vara_int(cdfid,varid,start,count,
     $                                          wmoid_list(1,i))
                       if(status.ne.nf_noerr)go to 900
                    endif
                 endif
c
c       get the dimension id of the record dimension into the list where it
c       can be used by prof_cdf_read.  note the +1 --> in the old netcdf
c       bug-riddled library ncinq returns 0-(ndims-1) but the fortran interface
c       requires 1-ndim.  note that in the current version of the netcdf
c       library (4/96) this bug has been fixed, but other problems exist.  thus
c       we're sticking with the old version until porting this to unix.
c
c       3/17/97: have ported to unix, and have removed the +1.  if this rears
c       it's head again the old code was: recdim_list(i) = recdim + 1
c
                status = nf_inq(cdfid,ndims,nvars,ngatts,recdim)
                if(status.eq.nf_noerr)recdim_list(i) = recdim
c
              endif
              go to 900                 ! done processing with the open slot.
           endif
        enddo
c
        status = -2                     ! no open slots.
c
900     if(status.eq.0)then
c
           open_list(i) = .true.
           cdfid_list(i) = cdfid
c
        else
c
           cdfid = 0
c
           if(error_flag.gt.0.and.status.eq.-2)then
c
              write(*,*)'prof_cdf_open:  attempt to open more than 4'
     $                  //' files.'
c
           endif
c
 
           if(error_flag.eq.2)stop
 
c
        endif
c
        return
        end
c       prof_cdf_read.for               michael barth           11-aug-1993
c
c       this subroutine is used to read all of one variable for one station
c       from a netcdf file holding profiler data.
c
c       the user must have previously opened the file via prof_cdf_open.
c
c input:
c
c       cdfid           integer:  the netcdf id for the file.
c
c       staname         character*6:  the name of the station (e.g., 'pltc2 ' --
c                       note the trailing blank).  this is the primary way the
c                       user can request a station.  if this is all blanks,
c                       the second method, wmoid, will be used.
c
c                       ***** note that the station name must be all capital
c                       letters (and a number).
c
c       wmoid           integer:  the wmo block and station number (e.g.,
c                       74533).  this is the alternative method of selecting
c                       the station.
c
c                       ***** note that if staname is not all blank, wmoid
c                       will be ignored.
c
c       varname         character*(*):  a character string giving the name of
c                       the desired variable (use the names in the netcdl
c                       file describing the dataset in use).  this is
c                       case-sensitive so the exact string must be used.
c
c       dtype           integer:  0 - the variable is numeric data.
c                                   1 - the variable is character data.
c
c output:
c
c       array           array of data returned to the caller.  the data type
c                       and dimensions are not used in this routine, or
c                       validated.  the caller must use the appropriate type
c                       (see the netcdl file) and know the expected amount
c                       of data.
c
c       status          integer:  returned status:  0 is good, -1 and +n
c                       are errors returned by netcdf.  the definitions of
c                       these can be found in netcdf.inc.  in addition, these
c                       error returns are from prof_cdf_read:
c
c                       -3      station not found.
c                       -4      invalid dtype.
c                       -5      netcdf file not open.
c                       -9      too many dimensions in variable (16 maximum).
c
c modifications:
c
c       15-sep-1994/mfb fixed a really stupid oversight:  the routine didn't
c                       work with 2-dimensional arrays for a single record
c                       (e.g., moments - rec_num, level, beam).  part of this
c                       fix was to remove the narray argument into this
c                       routine.  that was the size of the expected array.
c                       as the fix included the need to determine that on a
c                       dimension-by-dimension basis, narray became redundant.
c
c       11-sep-1995/mfb change the dimensions of start & count to use maximum
c                       number of allowable netcdf dimensions.  also changed
c                       how character data are read to fix a bug revealed by
c                       the new blp cdl.
c
c       28-nov-1995/lab put prof_cdf_common in it's own include file.
c
c       29-feb-1996/mfb changed declaration of array from byte to integer to
c                       make the code more portable.
c
c       05-mar-1996/mfb limited the check of staname to the first five
c                       characters.  this was done to make the code compatible
c                       with all known profiler-based data formats in fsl.
c                       (some have a blank in the 6th character, some have a
c                       zero.)
c
c       07-mar-1996/mfb added the use of recdim_list to determine if a
c                       one-dimensional variable's dimension is the record
c                       dimension or not -- the hyperslab count
c                       depends on this distinction.
c
c       21-mar-1996/mfb removed declarations that duplicate those in new
c                       netcdf.inc.
c
c       07-aug-2000/lsw modified to netcdf version 3.4
c
        subroutine prof_cdf_read(cdfid,staname,wmoid,varname,dtype,
     $                           array,status)
c
        implicit none
c
        integer cdfid,wmoid,dtype,status
        character*6 staname
        character*(*) varname
        integer array(*)
c
c       internal prof_cdf information:
c
c       open_list       logical flags indicating which of the four file slots
c                       are in use.
c
c       cdfid_list      netcdf id's for the read files.
c
c       staname_list    station names for the records in the read files.
c
c       wmoid_list      wmo id's for the records in the read files.
c
c       nrec_list       number of records (stations) in each read file.
c
c       recdim_list     record dimension id's for each file.
c
c       error_flag      error handling flag.
c
        integer max_profilers
        parameter (max_profilers = 200)
 
        include 'prof_cdf_common.inc'
        include 'netcdf.inc'
c
        character*21 error_text(3)
        data error_text/'station not found.   ',
     $                  'invalid dtype.       ',
     $                  'netcdf file not open.'/
        character*(maxncnam) name
        integer vartyp,nvdims,vdims(maxvdims),nvatts,dimsiz
        integer i,j,k,varid,start(maxvdims),count(maxvdims)
        integer count_total
        character*5 c5_1,c5_2
c
        do i = 1, 4
c
c       find the user's cdfid in our internal list.
c
           if(cdfid_list(i).eq.cdfid.and.open_list(i))then
c
c       locate the station.
c
              do j = 1, nrec_list(i)
c
                 c5_1 = staname_list(j,i)(1:5)
                 c5_2 = staname(1:5)
c
                 if((staname.ne.'      '.and.c5_1.eq.c5_2).or.
     $              (staname.eq.'      '.and.
     $               wmoid_list(j,i).eq.wmoid))then
c
c       get the variable id.
c
                    status = nf_inq_varid(cdfid,varname,varid)
                    if(status.ne.nf_noerr)go to 900
c
c       find out how many dimensions are in the variable.
c
                    status = nf_inq_var(cdfid,varid,name,vartyp,
     $                                  nvdims,vdims,nvatts)
                    if(status.ne.nf_noerr)go to 900
c
c       use this information to, one dimension at a time, fill in the
c       dimensions of the hyperslab to be obtained.
c
                    count_total = 1     ! must calculate product of count
c                                       ! vector for character data type.
                    do k = 1, nvdims
                       status = nf_inq_dim(cdfid,vdims(k),name,dimsiz)
                       if(status.ne.nf_noerr)go to 900
                       if(k.lt.nvdims.or.
     $                    (nvdims.eq.1.and.vdims(k).ne.
     $                     recdim_list(i)))then
                                start(k) = 1
                                count(k) = dimsiz
                                count_total = count_total * dimsiz
                       else
                                start(k) = j    ! rec_num dimension
                                count(k) = 1
                       endif
                    enddo
c
c       call the appropriate routine to get the variable.
c
                    if(dtype.eq.0)then  !reading an integer
c
                      status = nf_get_vara_int(cdfid,varid,start,count,
     $                                         array)
c
                    else if(dtype.eq.1)then   !reading text
c
                      status = nf_get_vara_text(cdfid,varid,start,count,
     $                                          array)
c
                    else if(dtype.eq.2)then   !reading real/float
c
                      status = nf_get_vara_real(cdfid,varid,start,count,
     $                                          array)
                    else
c
                             status = -4                ! invalid dtype.
c
                    endif
                    if (status.ne.nf_noerr) print *, nf_strerror(status)
c
                    go to 900                   ! end of found station.
c
                 endif
              enddo
c
              status = -3                       ! station not found.
              go to 900
c
           endif                                ! end of found open file.
c
        enddo
c
        status = -5                             ! file not found.
c
900     if(status.lt.-1.and.error_flag.ne.0)then
c
c       use the non-default error handling for prof_cdf errors.
c
           
                write(*,*)'prof_cdf_read:  '//
     $                    error_text(abs(status)-2)
                if(error_flag.eq.2)stop
        endif
c
        return
        end
c       prof_cdf_set_error.for          michael barth           12-aug-1993
c
c       this subroutine is used to set error handling for prof_cdf routines.
c       the caller doesn't have to use this feature, if no call is made the
c       default error handling will be performed.
c
c input:
c
c       error_code      0       return status codes -- this is the default.
c                       1       return status codes and write an error message
c                               to standard output (sys$output on vms).
c                       2       write an error message to standard output and
c                               exit the program.
c
c output:
c
c       status          integer:  returned status:  0 is good, -1 and +n are
c                       errors returned by netcdf.  the definitions of these can
c                       be found in netcdf.inc.  in addition, these error
c                       returns are from prof_cdf_set_error:
c
c                       -6      invalid error_code.
c
c modifications:
c
c       28-nov-1995/lab put prof_cdf_common in it's own include file.
c
c       27-may-1997/mfb renamed "ncopts" common and variable so it doesn't
c                       conflict with symbols used in c netcdf interface.
c
c       07-aug-2000/lsw modified to netcdf version 3.4 ncpopt is no longer
c                       supported, and as such is no longer called.  
c
        subroutine prof_cdf_set_error(error_code,status)
c
        implicit none
c
        integer error_code,status
c
c       internal prof_cdf information:
c
c       open_list       logical flags indicating which of the four file slots
c                       are in use.
c
c       cdfid_list      netcdf id's for the open files.
c
c       stanam_list     station names for the records in the open files.
c
c       wmoid_list      wmo id's for the records in the open files.
c
c       nrec_list       number of records (stations) in each open file.
c
c       error_flag      error handling flag.
c
        integer max_profilers
        parameter (max_profilers = 200)
c
        include 'prof_cdf_common.inc'
        include 'netcdf.inc'
c
        integer ncopts_val
        common/ncopts_cmn/ncopts_val            ! netcdf error handling flag.
c
c       validate the error_code and then save it in the prof_cdf common.  then
c       set the netcdf version (ncopts_val) accordingly.
c
        status = 0                              ! assume success.
c
        if(error_code.eq.0)then
c
                error_flag = 0
 
c
        else if(error_code.eq.1)then
c
                error_flag = 1
c
        else if(error_code.eq.2)then
c
                error_flag = 2
        else
                status = -6                     ! invalid error_code.
        endif
c
        return
        end
c       prof_cdf_get_stations.for       michael barth           19-oct-1995
c
c       this subroutine is used to get a list of stations from a netcdf file
c       holding profiler data.
c
c       the user must have previously opened the file via prof_cdf_open.
c
c input:
c
c       cdfid           integer:  the netcdf id for the file.
c
c output:
c
c       n_stations      number of stations currently in the file.
c
c       staname         character*6 array:  the names of the stations (e.g.,
c                       'pltc2 ').
c
c       wmoid           integer array :  the wmo block and station numbers
c                       (e.g., 74533).
c
c       status          integer:  returned status:  0 is good.  no netcdf
c                       calls are performed, so no errors are returned by
c                       netcdf.  in addition, these error returns are from
c                       prof_cdf_get_stations:
c
c                       -5      netcdf file not open.
c
c modfications:
c       28-nov-1995/lab put prof_cdf_common in it's own include file.
c
c       07-aug-2000/lsw modified to netcdf version 3.4
c
        subroutine prof_cdf_get_stations(cdfid,n_stations,staname,
     $                                   wmoid,status)
c
        implicit none
c
        integer cdfid,n_stations,wmoid(*),status
        character*6 staname(*)
c
c       internal prof_cdf information:
c
c       open_list       logical flags indicating which of the four file slots
c                       are in use.
c
c       cdfid_list      netcdf id's for the read files.
c
c       stanam_list     station names for the records in the read files.
c
c       wmoid_list      wmo id's for the records in the read files.
c
c       nrec_list       number of records (stations) in each read file.
c
c       error_flag      error handling flag.
c
        integer max_profilers
        parameter (max_profilers = 200)
 
        include 'prof_cdf_common.inc'
 
c
        character*21 error_text
        data error_text/'netcdf file not open.'/
        integer i,j
c
        do i = 1, 4
c
c       find the user's cdfid in our internal list.
c
           if(cdfid_list(i).eq.cdfid.and.open_list(i))then
c
              status = 0
              n_stations = nrec_list(i)
              if(n_stations.ne.0)then
c
                 do j = 1, nrec_list(i)
                    staname(j) = staname_list(j,i)
                    wmoid(j) = wmoid_list(j,i)
                 enddo
c
              endif
c
              go to 900
c
           endif                                ! end of found open file.
c
        enddo
c
        status = -5                             ! file not found.
c
900     if(status.ne.0.and.error_flag.ne.0)then
c
c       use the non-default error handling for prof_cdf errors.
c
                write(*,*)'prof_cdf_get_stations:  '//error_text
                if(error_flag.eq.2)stop
        endif
c
        return
        end


        subroutine prof_i4_avg_wdw(i4_avg_wdw_sec, cdfid, istatus)
!  modified to pass in cdfid, and actually read file for data lw 8-27-98

        integer cdfid
	character*20 timestr
        character*13 attname
	integer    lenstr, sp_loc, units_loc, to_seconds, itime
        integer      tlen, ttype      ! attribute type and length

        include 'netcdf.inc'

c       netcdf file is opened, and accessed via cdfid
c	read "avgtimeperiod" global attribute into timestr 

        lenstr = len(timestr)  
        attname = 'avgtimeperiod'

c determine type of attribute 'avgtimeperiod'
        istatus = nf_inq_atttype(cdfid,nf_global,'avgtimeperiod',ttype)
        if (istatus .ne. nf_noerr) then  !error retrieving info about avgtimeperiod
          istatus = 0
          return
        endif

        if (ttype .ne. ncchar) then   !read avgtimeperiod as an integer

          istatus = nf_get_att_int(cdfid,nf_global,'avgtimeperiod',
     $                             itime)
	  if (istatus .ne. nf_noerr) then   !error retrieving avgtimeperiod from file
            istatus = 0
            return
          endif
	  
c  units assumed to be minutes
	  to_seconds = 60

        else  !looking for a character string

c	  format of avgtimestring should be a number followed by a space and
c           then a lower case string of units, ie. "60 minutes" or "6 minutes"
c	    verify that the units are minutes, seconds or hour, and convert
c	    the number string into an integer.  return value is via i4_avg_wdw_sec
c  	    and is in seconds 

          istatus = nf_get_att_text(cdfid,nf_global,'avgtimeperiod',
     $                              timestr)
	  if (istatus .ne. nf_noerr) then   !error retrieving avgtimeperiod from file
            istatus = 0
            return
          endif

	  sp_loc = index(timestr,' ')  !everything to the left should be number
	
c	  convert numerical string to number
101	  format(i1)
102	  format(i2)
103	  format(i3)
104	  format(i4)
          itime = -1 ! check at end to see if set

	  if (sp_loc .eq. 2) then !one digit number
 	    read(timestr(1:1),101)itime
          endif

	  if (sp_loc .eq. 3) then !two digit number
 	    read(timestr(1:2),102)itime
          endif

	  if (sp_loc .eq. 4) then !three digit number
 	    read(timestr(1:3),103)itime
          endif

	  if (sp_loc .eq. 5) then !four digit number
 	    read(timestr(1:4),104)itime
          endif

	  if (itime .eq. -1) then   ! couldn't convert timestr to itime, return error
            write(6,*)
     1       ' error: couldnt convert avgtimeperiod string to integer: '  
     1       , timestr     
	    istatus = 0
            return
	  endif

c	  determine units in timestr
          to_seconds = 0 !will check after looking for units to see if set
   	  units_loc = index(timestr,'minute')
          if (units_loc .gt. 0) to_seconds = 60
 
	  units_loc = index(timestr,'second')
          if (units_loc .gt. 0) to_seconds = 1
 
	  units_loc = index(timestr,'hour')
          if (units_loc .gt. 0) to_seconds = 3600

        endif

	if (to_seconds .eq. 0) then   ! couldn't identify units, return error
          write(6,*)
     1      ' error: couldnt decode avgtimeperiod from file ', timestr       
	  istatus = 0
          return
        else
          i4_avg_wdw_sec = itime * to_seconds
          istatus = 1
	endif
        
        return
        end
