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
        subroutine glst (i4time,t,ii,jj,istatus)

c       $log: glst.for,v $
c revision 1.3  1996/03/13  22:43:39  birk
c changed read_laps_data to read_laps to allow a sub-hourly cycle
c
c revision 1.2  1995/09/13  21:35:22  birk
c added disclaimer to files
c
c revision 1.1  1994/04/25  15:13:22  birk
c initial revision
c

c       gets laps surface temperature

c       input i4time

        implicit none

c        include 'lapsparms.for'
c        include 'parmtrs.inc'

c parameter variables

      integer i4time
      integer ii,jj
      integer istatus
      real t(ii,jj)

c variables dynamically dependent on lapsparms.inc

      real data(ii,jj,1)


c internal variables

      integer i,j

      character*256 dir
      character*31 ext
      character*3 var (1)
      integer lvl (1)
      character*4 lvl_coord  (1)
      character*10 units (1)
      character*125 comment (1)
      integer kmax, len

      real numel,sum
      data ext /'lsx'/
c code--------------------------------------
      call get_directory(ext,dir,len)

        numel = 0.0
        sum = 0.0

        istatus = 0

        lvl (1) = 0


        var(1) = 't'
        kmax = 1

        call read_laps (i4time,i4time,dir,ext,ii,jj,
     1  kmax,1,var,lvl,
     1  lvl_coord,units,comment,data,istatus)

        if (istatus.ne.1) then

        write(6,*) 'error in routine glst'
        return
        endif

        write (6,*) istatus
        write (6,*) lvl_coord,units,comment


        do j = 1,jj
        do i = 1,ii

        numel = numel + 1.
        sum = sum + data(i,j,1)
        t(i,j) = data(i,j,1)


        enddo
        enddo

        sum = sum /numel
        write (6,*) '  average value of field is ',sum


        istatus = 1

        return

        end
