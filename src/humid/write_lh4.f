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
      subroutine write_lh4 (i4time,tpw,bias_one,ii,jj,istatus)


      implicit none

c     parameter variables

      integer ii,jj
      real tpw (ii,jj)      
      real bias_one
      integer istatus
      
c     internal variables
      
      integer
     1     i4time,
     1     kmax,
     1     lvl(1)
      real mdf
      
      integer i,j, len
      
      character
     1     tpdir*250,
     1     tpext*31,
     1     var(1)*3,
     1     lvl_coord(1)*4,
     1     units(1)*10,
     1     comment(1)*125
      
      data var/1*'tpw'/
      data lvl_coord/1*'  '/
      data units/1*'m  '/
      data tpext /'lh4'/
c     set missing data flag
      call get_r_missing_data(mdf, istatus)
      if (istatus.ne.1) then
         write(6,*) 'fatal error in assigning missing data flag'
         write(6,*) 'assigning default value of 1e+37'
         mdf = 1.e37
         write(6,*) 'continuing using default'
      endif
      
      call get_directory(tpext,tpdir,len)
      
      istatus = 0               ! bad
      
      if(tpw(1,1).eq.mdf) then  ! field is not valid
         write(6,*) 'tpw field not valid'
         return
      endif
      
      write(comment(1),1)  bias_one
 1    format('vsn 5.1 tpw via sh integration, bias correction = ',
     1     e20.10)
      
      kmax = 1
      
c     convert from cm to meters for archive
      
      do i  = 1,ii
         do j  = 1,jj
            
            tpw(i,j) = tpw(i,j) * 1.e-2
c     place mod for missing data flag
            if (tpw(i,j) .lt. 0.0) tpw (i,j) = mdf
            
         enddo
      enddo
      
c     initialize 

      lvl(1) = 0 ! assigned for ibm problems (7/21/2000) db
      
      call write_laps (i4time,i4time,
     1     tpdir,
     1     tpext,
     1     ii,
     1     jj,
     1     1,
     1     1,
     1     var,
     1     lvl,
     1     lvl_coord,
     1     units,
     1     comment,
     1     tpw,
     1     istatus)
      
      if(istatus.ne.1) then
         istatus = 134316524
         return
      endif
      
      istatus = 1               ! success
      
      print*, 'write laps tpw complete',istatus
      
      return
      end
