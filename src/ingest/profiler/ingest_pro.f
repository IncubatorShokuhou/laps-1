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

!       1997 jul      ken dritz        added call to get_grid_dim_xy.
!       1997 jul      ken dritz        pass nx_l, ny_l to ingest_pro.

        character*9 a9_time
        character*8 c8_project,c8_blp_format,c8_blp_format_in
        character*3 ext

        call getenv('laps_a9time',a9_time)
        call s_len(a9_time,ilen)

        lun_out = 1

        call get_c8_project(c8_project,istatus)
        if(istatus .ne. 1)goto999

        call get_c8_blpfmt(c8_blp_format_in,istatus)
        if(istatus .ne. 1)goto999

        if(ilen .eq. 9)then
            write(6,*)' systime (from env) = ',a9_time
            call i4time_fname_lp(a9_time,i4time,istatus)
        else
            call get_systime(i4time,a9_time,istatus)
            if(istatus .ne. 1)go to 999
            write(6,*)' systime = ',a9_time
        endif

        call get_grid_dim_xy(nx_l,ny_l,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'error getting horizontal domain dimensions'
           go to 999
        endif

!       if(i4time .eq. (i4time / 3600) * 3600)then
            write(6,*)
            write(6,*)' running wpdn (nimbus/wfo) profiler ingest'       
            call ingest_pro(i4time,nx_l,ny_l,lun_out,j_status)
            write(6,*)' return from wpdn (nimbus/wfo) profiler ingest'              

!       else
!           write(6,*)' not on the hour, no wpdn profiler ingest run'       

!       endif

        if(c8_blp_format_in .eq. 'default')then
            c8_blp_format = c8_project
        else
            c8_blp_format = c8_blp_format_in
        endif            

        if(c8_blp_format .eq. 'rsa')then
            write(6,*)
            write(6,*)' running rsa/ldad local wind profile ingest '
            call ingest_rsapro(i4time,nx_l,ny_l,lun_out,j_status)
            write(6,*)' return from rsa/ldad local wind profile ingest'
        elseif(c8_blp_format .eq. 'wfo' .or. 
     1         c8_blp_format .eq. 'madis')then
            write(6,*)
            write(6,*)' running madis (wfo) multi-agency profile ingest'       
            ext = 'pro'
            call ingest_madis_map(i4time,nx_l,ny_l,ext,lun_out
     1                           ,j_status)
            write(6,*)' return from madis (wfo) map ingest'
        else
            write(6,*)
            write(6,*)' running blp (nimbus) local profiler ingest'
            call ingest_blppro(i4time,nx_l,ny_l,lun_out,j_status)
            write(6,*)' return from blp (nimbus) local profiler ingest'
        endif

        write(6,*)
        write(6,*)' running vad (nimbus) ingest'
        call ingest_vad(istatus)
        write(6,*)' return from vad (nimbus) ingest'

 999    continue

        close(lun_out)

        end


       subroutine get_pro_parms(c8_blp_format,istatus)

!      this subroutine and namelist isn't being used at the present time.

       character*8 c8_blp_format

       namelist /pro_nl/ c8_blp_format
 
       character*150 static_dir,filename
 
       call get_directory('static',static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/ingest_pro.nl'
 
       open(51,file=filename,status='old',err=900)
       read(51,pro_nl,err=901)
       close(51)

       print*,'success reading pro_nl in ',filename
       write(*,pro_nl)

       istatus = 1
       return

  900  print*,'error opening file ',filename
       istatus = 0
       return

  901  print*,'error reading pro_nl in ',filename
       write(*,pro_nl)
       istatus = 0
       return

       end
