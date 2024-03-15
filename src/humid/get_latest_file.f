cdis   
cdis    open source license/disclaimer, forecast systems laboratory
cdis    noaa/oar/fsl, 325 broadway boulder, co 80305
cdis    
cdis    this software is distributed under the open source definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    in particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - all modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - if significant modifications or enhancements are made to this
cdis    software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    this software and its documentation are in the public domain
cdis    and are furnished "as is."  the authors, the united states
cdis    government, its instrumentalities, officers, employees, and
cdis    agents make no warranty, express or implied, as to the usefulness
cdis    of the software and documentation for any purpose.  they assume
cdis    no responsibility (1) for the use of the software and
cdis    documentation; or (2) to provide technical support to users.
cdis   
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
        subroutine get_latest_file (path,i4time,filename,istatus)


c       $log: get_latest_file.for,v $
c revision 1.1  1996/08/30  20:45:49  birk
c initial revision
c

c this routine returns the latest filename in a given path sans the .ext.
c it is useful to determine i4times from existing files on the system.

c author dan birkenheuer

c date: 8/14/96

        implicit none

        character*256 path
        character*9 filename,f_save
        integer i4time
        integer istatus
        
        integer numoffiles,maxfiles
        parameter (maxfiles = 3000)
        character*256 c_filenames(maxfiles),dum
        integer i4time_test,i,mintime
        
        
        
        
        call get_file_names (path, numoffiles,
     1       c_filenames, maxfiles, istatus)
        
        if (istatus.eq.1 .and. numoffiles.gt.0) then
           
           mintime = 1000000000
           
           do i = 1,numoffiles
              
              dum = c_filenames(i)
              
              filename = dum (index(dum,' ')-13:index(dum,' ')-4)
              
              call i4time_fname_lp (filename, i4time_test, istatus)
              
              mintime = min(abs(i4time_test-i4time),mintime)
              
              if(mintime.eq.abs(i4time_test-i4time) ) f_save = filename
              
           enddo
           
           filename = f_save

           write (6,*) 'file found is ',mintime,'(sec) old'
           
           if (mintime.lt.3600) return
           
        endif

        write (6,*) 'file found is ',mintime,'(sec) old'
        
        istatus = 0
        
        return
        end

