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


        function aski4t()

        character*11 asc9_tim
        integer aski4t,lent

        write(6,*)'enter time, <yydddhhmm>, or <return> for latest:'

        read(5,1)asc9_tim
1       format(a)

        call s_len(asc9_tim,lent)

        if(lent .eq. 0)then
            i4time = i4time_now_gg()
        elseif(lent .eq. 11)then
            call cv_asc_i4time(asc9_tim(1:9),i4time)
            read(asc9_tim,2)isec
2           format(9x,i2.2)
            i4time = i4time + isec
            write(6,*)' 11 character input, i4time = ',i4time
        else ! lent = 9
            call cv_asc_i4time(asc9_tim(1:9),i4time)
        endif

        aski4t = i4time

        return
        end
