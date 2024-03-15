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
      subroutine draw_county_map(sw,ne,jproj,polat,polon,rrot,jdot
     1                          ,icol_sta,icol_cou,ni,nj,namelist_parms)       

      include 'lapsplot.inc'
c
      real sw(2),ne(2),pl1(2),pl2(2),pl3(2),pl4(2),
     +       polat,polon,rrot

c
!abdel      
      double precision grsp,dpolat,dpolon,dsw(2),dne(2)
      integer irgl
      
      data grsp,irgl / 1.d0 , 0  /
      parameter (ncra=10000,ngps=10,lrwk=2*ncra)

      dimension xcra(ncra),ycra(ncra)
      dimension rwrk(lrwk)     
      equivalence (rwrk(1),xcra(1)),(rwrk(ncra+1),ycra(1))
! abdel     

      integer jproj,jjlts,jgrid,jus,jdot,ier
c
      common/supmp9/ds,di,dsrdi
      common /zoom/       zoom 
!     di = 50.
!     polat=90.

!     rrot=0.
      pl1(1)=sw(1)
      pl2(1)=sw(2)
      pl3(1)=ne(1)
      pl4(1)=ne(2)
      jjlts=-2
      jgrid=0
 

!     call get_lapsplot_parms(namelist_parms,istatus)       

!     1 means use local version of supmap, 2 means use the ncarglib version
!     3 means to try newer (ezmap) routines
      mode_supmap = namelist_parms%mode_supmap

      call get_grid_spacing(grid_spacing_m,istatus)
      if(istatus .ne. 1)then
          write(6,*)' no grid spacing, stop in draw_county_map'
          stop
      else
          write(6,*)
          write(6,*)' subroutine draw_county_map...',mode_supmap,jproj
      endif

      domsize = (float(nj)-1.) * grid_spacing_m / zoom

!     plot counties
      if(jdot .eq. 1)then
          call gsln(3) ! dotted
      else
          call gsln(1) ! solid
      endif

      jgrid=namelist_parms%latlon_int        ! draw lat/lon lines?

      if(domsize .le. 1500e3 .and. 
     1   namelist_parms%continent_line_width .gt. 0)then
          write(6,*)' plotting counties ',domsize,mode_supmap,jgrid
          call setusv_dum(2hin,icol_cou)

          if(mode_supmap .eq. 1)then
              jus=-4
              call supmap_local(jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,
     +                          jjlts,jgrid,jus,jdot,ier)
          elseif(mode_supmap .eq. 2)then
              iout = 0
              call supmap      (jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,
     +                          jjlts,jgrid*1000,iout,jdot,ier)
          elseif(mode_supmap .eq. 3)then
              write(6,*)' calling mplndr, etc. for counties...'
!             call mapdrw()
              call mapint
!             call maplot
              call mplndr ('earth..2',5)
              if(jgrid .gt. 0)then ! draw lat/lon lines
                  call mpsetr('gr',float(jgrid))
                  call mapgrd()
              endif
c abdel	      
          elseif(mode_supmap .eq. 4)then
              write(6,*)' calling sub submap=4 for europe...'
              dpolat=polat
              dpolon=polon
              dsw(1)=sw(1)
              dsw(2)=sw(2)
              dne(1)=ne(1)
              dne(2)=ne(2)
              call mdproj ('st',dpolat,dpolon,0.d0)
              call mdpset ('co',dsw(1),dsw(2),dne(1),dne(2))
              call mapstd ('gr',grsp)
              call mdrgol (irgl,rwrk,lrwk)    
          
          elseif(mode_supmap .eq. 5)then
              write(6,*)' accessing rangs database...'
              write(6,*)' not yet supported - stop' 
              stop
          endif
          if(ier .ne. 0)write(6,*)' ier = ',ier

          call sflush

          jgrid=0                        ! do not draw subsequent lat/lon lines

      else
          write(6,*)' omitting counties ',domsize
     1             ,namelist_parms%continent_line_width

      endif

      call gsln(1)
      call setusv_dum('in',namelist_parms%icol_state)

      call gslwsc(namelist_parms%continent_line_width)

      if(mode_supmap .eq. 1)then
          write(6,*)' plotting states from counties ',mode_supmap,jgrid
          jus=-8
          call supmap_local(jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,
     +                      jjlts,jgrid,jus,jdot,ier)
      elseif(mode_supmap .eq. 2)then
          write(6,*)' plotting states from counties ',mode_supmap,jgrid
          iout = 0
          call supmap      (jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,
     +                      jjlts,jgrid*1000,iout,jdot,ier)
      elseif(mode_supmap .eq. 3)then
          call mapint

          if(namelist_parms%state_line_width .gt. 0. .or.
     1       namelist_parms%country_line_width .gt. 0.     )then
              write(6,*)' calling mapdrw, etc. for countries/states...'
     1                 ,namelist_parms%state_line_width
     1                 ,namelist_parms%icol_country
     1                 ,icol_sta
              call gslwsc(namelist_parms%state_line_width)
              call setusv_dum('in',namelist_parms%icol_state) 
              call mplndr ('earth..2',4) ! states & countries
          else
              write(6,*)' skip plot of countries/states'
              if(namelist_parms%continent_line_width .gt. 0.)then
                  write(6,*)' calling mapdrw just for continents...'
                  call gslwsc(namelist_parms%continent_line_width)
                  call setusv_dum('in',namelist_parms%icol_continent) 
                  call mplndr ('earth..2',2) ! continents
              endif
          endif

          if( (namelist_parms%icol_country .ne. 
     1         namelist_parms%icol_state) 
     1                     .or.
     1        (namelist_parms%country_line_width .eq. 0. 
     1                     .and.
     1         namelist_parms%state_line_width .gt. 0.)       
     1                                                 )then
              write(6,*)
     1              ' replotting countries in separate color/width: '
     1                 ,0
     1                 ,namelist_parms%state_line_width
              call gslwsc(namelist_parms%state_line_width)
              if(namelist_parms%country_line_width .eq. 0.)then
                  call setusv_dum('in',0) 
              else
                  call setusv_dum('in',namelist_parms%icol_country) 
              endif
              call mplndr ('earth..2',3) ! countries
          endif

          if(jgrid .gt. 0)then ! draw lat/lon lines
              call setusv_dum(2hin,icol_cou)
              call mpsetr('gr',float(jgrid))
              call mapgrd()
          endif
c abdel	  
      elseif(mode_supmap .eq. 4)then
            write(6,*)' calling sub submap=4 for europe...'
            dpolat=polat
            dpolon=polon
            dsw(1)=sw(1)
            dsw(2)=sw(2)
            dne(1)=ne(1)
            dne(2)=ne(2)
            call mdproj ('st',dpolat,dpolon,0.d0)
            call mdpset ('co',dsw(1),dsw(2),dne(1),dne(2))
            call mapstd ('gr',grsp)
            call mdrgol (irgl,rwrk,lrwk)	  
	  
      else
          write(6,*)' accessing rangs database...'
          write(6,*)' not yet supported - stop' 
          stop

      endif

      if(ier .ne. 0)write(6,*)' ier = ',ier

      call sflush

      call gsln(1)
      call setusv_dum(2hin,icol_sta)

      jgrid=0                                ! do not draw lat/lon lines
       
      if(mode_supmap .eq. 1)then
          write(6,*)' plotting continents ',mode_supmap,jgrid
          jus=-1
          call supmap_local(jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,
     +                      jjlts,jgrid,jus,jdot,ier)
      elseif(mode_supmap .eq. 2)then
          write(6,*)' plotting continents ',mode_supmap,jgrid
          iout = 2
          call supmap      (jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,
     +                      jjlts,jgrid*1000,iout,jdot,ier)
      endif
      if(ier .ne. 0)write(6,*)' ier = ',ier

      call gslwsc(1.0)

      call sflush

      write(6,*)

      return
      end

cabdel       
              subroutine mdrgdi (dinm)
c
c this is a user-replaceable routine that returns the name of the
c directory in which the rangs/gshhs data files have been placed.
c
       character*(*) dinm

c fitxer bo
c return the name of the directory where the rangs/gshhs data reside.
c
       dinm='/usr/local/ncarg/lib/ncarg/database/rangs_gshhs'
c
c done.
c
       return

       end
