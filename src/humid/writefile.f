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
        subroutine writefile (i4time,commentline,lvl,data,
     1                         ii,jj,kk,istatus)

c       $log: writefile.for,v $
c revision 1.4  1996/03/05  21:12:06  birk
c modified i/o interface to enable subhourly cycle
c
c revision 1.3  1995/09/13  21:36:17  birk
c added disclaimer to files
c
c revision 1.2  1994/10/28  21:36:21  birk
c removed unused variables
c placed data in parameter list
c
c revision 1.1  1994/04/25  15:05:17  birk
c initial revision
c

c     module writes laps data base file of specific humidity at the
c     various levels

c     author d. birkenheuer
c     changed to pressure coords  14 aug. 89

c     24 april 1989
c     updated for unix 19 oct 1992  db

        implicit none


c     include 'lapsparms.for'
c     include 'parmtrs.inc'



c     parameter variables

        integer ii,jj,kk
        real data(ii,jj,kk)
        integer i4time 
        character*125 commentline
        integer lvl (kk)
        integer istatus

c     variables relying on dynamic initialization

        character  var(kk)*3,
     1       lvl_coord(kk)*4,
     1       units(kk)*10,
     1       comment(kk)*125

c        data var/kk*'sh '/
c        data lvl_coord/kk*'hpa'/
c        data units/kk*'          '/


c     internal variables

        integer kmax
        character
     1       dir*250,
     1       extlt1*31,ext*50,rhext*50,extpw*50,ext3*50
        integer k

        data extpw/'lh1'/
        data ext3/'lh2'/
        data extlt1/'lt1'/
        data ext /'lq3'/
        data rhext /'lh3'/
        integer len


        do k = 1,kk
           var(k) = 'sh '
           units(k) = '          '
           lvl_coord(k) = 'hpa'
        enddo


c     call get_directory(extpw,dirpw,len)
c     call get_directory(ext3,dir3,len)
c     call get_directory(extlt1,dirlt1,len)
        call get_directory(ext,dir,len)
c     call get_directory(rhext,rhdir,len)


        kmax = kk

        do k  = 1,kk
           comment(k) = commentline
        enddo

        call write_laps (i4time,i4time,
     1       dir,
     1       ext,
     1       ii,
     1       jj,
     1       kk,
     1       kk,
     1       var,
     1       lvl,
     1       lvl_coord,
     1       units,
     1       comment,
     1       data,
     1       istatus)

        print*, ext, istatus, dir(1:len), len

        if(istatus.ne.1) then
        istatus = 134316524
        return
        endif

        istatus = 1

        print*, 'write laps data module complete',istatus


        return

        end

