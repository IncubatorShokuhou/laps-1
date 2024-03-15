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

!       1997 jul 31 k. dritz  - added call to get_grid_dim_xy to get the
!                               values of nx_l, ny_l.
!       1997 jul 31 k. dritz  - now pass nx_l, ny_l as arguments to laps_cloud.
!       1997 jul 31 k. dritz  - added call to get_meso_sao_pirep to get the
!                               value of n_pirep, which is passed to
!                               laps_cloud.
!       1997 jul 31 k. dritz  - added call to get_maxstns, and pass the value
!                               of maxstns to laps_cloud.
!       1997 jul 31 k. dritz  - compute max_cld_snd as maxstns + n_pirep and
!                               pass to laps_cloud.

        use mem_namelist, only: read_namelist_laps

        use mem_namelist, only: nx_l, ny_l, nk_laps, maxstns, n_pirep       

        integer j_status(20),iprod_number(20)
        character*150 static_dir,filename
        character*9 a9time

!       read global parameters into module memory structure
        call get_directory('static',static_dir,len_dir)
        filename = static_dir(1:len_dir)//'/nest7grid.parms'
        call read_namelist_laps('lapsparms',filename)

!       read cloud parameters into module memory structure
        filename = static_dir(1:len_dir)//'/cloud.nl'
        call read_namelist_laps('cloud_anal',filename)

        call get_systime(i4time,a9time,istatus)
        if(istatus .ne. 1)go to 999

        write(6,*)' systime = ',a9time

        isplit = 1

        max_cld_snd = maxstns + n_pirep
          
        call laps_cloud_sub(i4time,
     1                  nx_l,ny_l,
     1                  nk_laps,
     1                  n_pirep,
     1                  maxstns,
     1                  max_cld_snd,
     1                  i_diag,
     1                  n_prods,
     1                  iprod_number,
     1                  isplit,
     1                  j_status)

999     continue

        end

 



