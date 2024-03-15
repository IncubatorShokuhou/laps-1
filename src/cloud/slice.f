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

        subroutine slice (a,imax,jmax,kmax,b,ib,jb,hz,ew,ns,kk,jj,ii)
c this subroutne allows a slice to be taken through a 3-d array
c control is by setting integers hz,ew,ns, to 1 or 0.  kk,jj,ii
c then determine which slice is taken. result is put into b.
        integer hz,ew,ns
        dimension a(imax,jmax,kmax),b(ib,jb)

        if (hz.eq.1) then
             do j=1,jmax
             do i=1,imax
             b(i,j)=a(i,j,kk)
             enddo
             enddo
        endif

        if(ew.eq.1) then
                do k=1,kmax
                do i=1,imax
                b(i,k)=a(i,jj,k)
                enddo
                enddo
        endif

        if(ns.eq.1) then
                do k=1,kmax
                do j=1,jmax
                b(j,k)=a(ii,j,k)
                enddo
                enddo
        endif

        return
        end
