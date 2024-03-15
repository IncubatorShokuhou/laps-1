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
cdis cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis



        subroutine slv_laplc (data,mask, nx, ny)

c       $log: slv_laplc.for,v $
c revision 1.1  1996/08/30  20:57:55  birk
c initial revision
c

        implicit none

        integer nx,ny
        integer mask(nx,ny)
        real data(nx,ny),error
        real maxerror
        real typical_data

        integer i,j,k



        do k = 1,6000

        maxerror =0.0


        do j = 2,ny -1
           do i = 2,nx -1

              if(mask(i,j) .eq. 0) then

                 error = 0.25 * ( data(i+1,j) + data(i-1,j) +
     1                data(i,j+1) + data (i,j-1) ) - data (i,j)

                 data (i,j) = error + data(i,j)

                 maxerror = max(maxerror,abs(error))

              else 
                 typical_data = data(i,j)


              endif

           enddo
        enddo

c       print*, maxerror

        if(typical_data > 1.e-9) then
           if (maxerror/typical_data.le. 1.e-3) go to 22
        else
           if (maxerror.le. 1.e-3) then
             write(6,*) 'typical_data = 0.0, divide avoided'
             go to 22
           endif
        endif


        enddo

        write(6,*) 'diriclet terminated on iterations ', k

22      write (6,*) 'max error solving dirichlet problem ', 
     1       maxerror,'/',typical_data,k

c     fill boarders (normally zero) with nearest neighbor values

        do j = 1,ny
           do i = 1,nx
              if(i .eq. 1 ) then !boarder
                 data(i,j) = data(i+1, j)
              elseif(i .eq. nx) then
                 data(i,j) = data(i-1,j)
              endif
           enddo
        enddo
     
        do j = 1,ny
           do i = 1,nx
              if(j .eq. 1 ) then !boarder
                 data(i,j) = data(i, j+1)
              elseif(j .eq. ny) then
                 data(i,j) = data(i,j-1)
              endif
           enddo
        enddo
                   

        return
        end

