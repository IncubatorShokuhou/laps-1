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

        subroutine put_wind_2d(i4time,directory,ext,var,units,
     1                  comment,wind_2d,imax,jmax,istatus)

        character*150 directory
        character*31 ext

        character*125 comment,comment_2d(2)
        character*10 units,units_2d(2)
        character*3 var,var_2d(2)
        integer lvl,lvl_2d(2)
        character*4 lvl_coord,lvl_coord_2d(2)

        real wind_2d(imax,jmax,2)

        write(6,11)directory,ext(1:5)
11      format(' writing 2d wind ',a50,1x,a5,1x,a3)

        lvl = 0
        lvl_coord = 'msl'

        var_2d(1) = 'u'
        var_2d(2) = 'v'

        do k = 1,2
            comment_2d(k) = comment
            lvl_2d(k) = lvl
            lvl_coord_2d(k) = lvl_coord
            units_2d(k) = 'm/s'
        enddo

        call write_laps_data(i4time,directory,ext,imax,jmax,
     1  2,2,var_2d,lvl_2d,lvl_coord_2d,units_2d,
     1                     comment_2d,wind_2d,istatus)

        return
        end
