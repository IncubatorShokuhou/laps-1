



c
c********************** map projection routines ****************************
c
c
c
        subroutine xytops(x,y,pla,plo,erad)
c
c     this convert x,y-polar stereographic coordinates to
c     lat/lon values.
c     longitude:   0 - 180  ; greenwich and east, 0 - -180 greenwich and west
c     latitude : -90 -  90  ; positive for northern hemisphere
c     it is assumed that the x-axis point towards the east why the
c     longitude is rotated relative to the 'standard pol.ste.' location
c            
      pi180=3.14159265/180.
      vdist = erad*2.0
c
c     calculate distance from (0,0) and define (0,0) as 90,0 (90,-90 in
c     the rotated system)
c
      dist=sqrt(x**2+y**2)
      if(dist.eq.0) then
         pla= 90.0
         plo=-90.0
      else
c
c     calculate the latitude by means of atan
c
         pla=atan(dist/vdist)/pi180
         pla=90.0-2.0*pla
c
c     calculate the longitude taking the directions into account
c
         if(x.eq.0.0) then
            if(y.gt.0.0) then
               plo= 90.0
            else
               plo=-90.0
            end if
         else
            if(x.gt.0.0) then
               plo=atan(y/x)/pi180
            else
               plo=atan(y/x)/pi180+180.0
            end if
         end if
      end if
c
c     rotate the longitude
c
      plo=amod(plo+450.0,360.0)
      return
      end
c
c     *******************************************************************
c
      subroutine pstoge(pla,plo,glat,glon,rlat,wlon1)
c
c     convert polar stereographic coordinates to geographical lat/lon
c     ditto with the pol.ste. pole at rlat,wlon1
c     
c     longitude:   0 ; 360 positive east (on input)
c               -180 ; 180 positive east (on output)
c     latitude : -90 ;  90 posive on northern hemisphere
c     it is assumed that the polar stereographic coordinates have been
c     rotated to the standard format with 0 degrees longitude along wlon1
c
c     tsp 21 june 89
c
c     set flag for n/s hemisphere
c
      pi180=3.14159265/180.
c      
      if(rlat.ge.0.0) then
         hsign= 1.0
      else
         hsign=-1.0
      end if
c
c     test for a n/s pole case
c
      if(rlat.eq.90.0) then
	 glat=pla
         glon=mod(plo+wlon1,360.0)
         go to 2000
      end if
      if(rlat.eq.-90.0) then
         glat=-pla
         glon=mod(plo+wlon1,360.0)
         go to 2000
      end if
c
c     test for longitude on 'greenwich or date line'
c
      if(plo.eq.0) then
         glat=rlat-90.0+pla
         if(glat.lt.-90.0) then
            glat=-180.0-glat
            glon=mod(wlon1+180.0,360.0)
         else
            glon=wlon1
         end if
         go to 2000
      end if      
      if(plo.eq.180.0) then
         glat=rlat+90.0-pla
         if(glat.gt.90.0) then
            glat=180.0-glat
            glon=mod(wlon1+180.0,360.0)
         else
            glon=wlon1
         end if
         go to 2000         
      end if
c
c     determine longitude distance relative to wlon1 so it belongs to
c     the absolute interval 0 - 180
c
      argu1=plo
      if(plo.gt.180.0) argu1 = plo-360.0
c
c     get the latitude, the help circle bb and the longitude by first
c     calculating the argument and legalize it - then take the inverse fct.
c
      if(hsign.gt.0.0) then
         arg2a = sin(pla*pi180)*sin(hsign*rlat*pi180)+
     +        cos(pla*pi180)*cos(rlat*pi180)*cos((180.0-argu1)*pi180)
      else
         arg2a = sin(pla*pi180)*sin(hsign*rlat*pi180)+
     +        cos(pla*pi180)*cos(rlat*pi180)*cos(argu1*pi180)
      end if
      arg2a = min(arg2a, 1.0)
      arg2a = max(arg2a,-1.0)
      glat  = hsign*asin(arg2a)
c
      if(hsign.gt.0.0) then
         arg2a = cos(rlat*pi180)*sin(pla*pi180)+
     +        sin(rlat*pi180)*cos(pla*pi180)*cos(argu1*pi180)
      else
         arg2a = cos(rlat*pi180)*sin(pla*pi180)+
     +       sin(-rlat*pi180)*cos(pla*pi180)*cos((180.0-argu1)*pi180)
      end if
      arg2a = min(arg2a, 1.0)
      arg2a = max(arg2a,-1.0)      
      bb    = acos(arg2a)
c
!     arg2a = cos(glat)*cos(bb)/(1.0-sin(glat)**2)
      arg2a = cos(bb)/cos(glat)              ! 1997 steve albers
      arg2a = min(arg2a, 1.0)
      arg2a = max(arg2a,-1.0)      
      glon  = acos(arg2a)
c     
c     convert the radians to degrees 
c
        glat = glat/pi180
        glon = glon/pi180
c
c       the solution is symmetric so the direction must be if'ed
c
        if(argu1.lt.0.0) then
           glon = 360.0-glon
        end if
        glon=amod(glon+wlon1,360.0)
c
 2000 continue
c
c     the resultant longitude must be in the interval from -180, 180
c      
      if(glon.gt.180.0) glon=glon-360.0
      return
      end




c     ****************************************************************
c
      subroutine getops(pla_r4,plo_r4,glat_r4,glon_r4,rlat_r4,wlon1_r4)       
c                        out     out    in      in      in       in
c
c     convert geographical lat/lon coordinates to polar stereographic
c     ditto with the pol.ste. pole at rlat,wlon1
c     
c     longitude:-180 ; 180 positive east (on input)
c              :   0 ; 360 positive east (on output)
c     latitude : -90 ;  90 posive on northern hemisphere
c     the result is rotated 270 degrees relative to 'standard pol.ste.'
c     wlon1 is defined in the same way as the input
c     approach so as to get the x-axis to point towards the east, and the
c     y-axis towards the north along 0 degrees (at np south along 180)
c
c     tsp 20/06-89
c
c     constants
c
      implicit real*8 (a-h,o-z)

      real pla_r4,plo_r4,glat_r4,glon_r4,rlat_r4,wlon1_r4

      pla = pla_r4
      plo = plo_r4
      glat = glat_r4
      glon = glon_r4
      rlat = rlat_r4
      wlon1 = wlon1_r4

      pi180 = 3.1415926535897932/180.0
c
c     set flag for n/s hemisphere and convert longitude to <0 ; 360> interval
c
      if(rlat.ge.0.0) then
         hsign= 1.0
      else
         hsign=-1.0
      end if
      glor=glon
      if(glor.lt.0.0) glor=360.0+glor
      rwlon1=wlon1
      if(rwlon1.lt.0.0) rwlon1=360.0+wlon1
c
c     test for a n/s pole case
c
      if(rlat.eq.90.0) then
         pla=glat
         plo=dmod(glor+270.0d0,360.0d0)
         go to 2000
      end if
      if(rlat.eq.-90.0) then
         pla=-glat
         plo=dmod(glor+270.0d0,360.0d0)
         go to 2000
      end if
c
c     test for ge coordinates at ps pole (steve albers - 1998)
c
      if(glat .eq. rlat .and. glon .eq. wlon1)then
         pla = 90.0
         plo = 0.
         goto2000
      endif
c
c     test for longitude on 'greenwich or date line'
c
      if(glor.eq.rwlon1) then
         if(glat.gt.rlat) then
            pla=90.0-glat+rlat
            plo=90.0
         else
            pla=90.0-rlat+glat
            plo=270.0
         end if
         go to 2000
      end if      
      if(dmod(glor+180.0d0,360.0d0).eq.rwlon1) then
         pla=rlat-90.0+glat
         if(pla.lt.-90.0) then
            pla=-180.0-pla
            plo=270.0
         else
            plo= 90.0
         end if
         go to 2000         
      end if
c
c     determine longitude distance relative to rwlon1 so it belongs to
c     the absolute interval 0 - 180
c
      argu1 = glor-rwlon1
      if(argu1.gt. 180.0) argu1 = argu1-360.0
      if(argu1.lt.-180.0) argu1 = argu1+360.0
c
c     1. get the help circle bb and angle alpha (legalize arguments)
c
      arg2a = cos(glat*pi180)*cos(argu1*pi180)
      arg2a = max(arg2a,-1.0d0)
      arg2a = min(arg2a, 1.0d0)         
      bb    = acos(arg2a)
c
      arg2a = sin(hsign*glat*pi180)/sin(bb)
      arg2a = max(arg2a,-1.0d0)
      arg2a = min(arg2a, 1.0d0)
      alpha = asin(arg2a)
c
c     2. get pla and plo (still legalizing arguments)
c
      arg2a = cos(rlat*pi180)*cos(bb)+
     +        sin(hsign*rlat*pi180)*sin(hsign*glat*pi180)
      arg2a = max(arg2a,-1.0d0)
      arg2a = min(arg2a, 1.0d0)         
      pla   = asin(arg2a)
c
      arg2a = sin(bb)*cos(alpha)/cos(pla)
      arg2a = max(arg2a,-1.0d0)
      arg2a = min(arg2a, 1.0d0)
      plo   = asin(arg2a)
c
c    test for passage of the 90 degree longitude (duallity in plo)
c         get pla for which plo=90 when glat is the latitude
c
      arg2a = sin(hsign*glat*pi180)/sin(hsign*rlat*pi180)
      arg2a = max(arg2a,-1.0d0)
      arg2a = min(arg2a, 1.0d0)         
      pla90 = asin(arg2a)
c
c         get help arc bb and angle alpha
c
      arg2a = cos(rlat*pi180)*sin(pla90)
      arg2a = max(arg2a,-1.0d0)
      arg2a = min(arg2a, 1.0d0)
      bb    = acos(arg2a)

      arg2a = sin(hsign*glat*pi180)/sin(bb)
      arg2a = max(arg2a,-1.0d0)
      arg2a = min(arg2a, 1.0d0)        
      alpha = asin(arg2a)
c
c         get glolim - it is nesc. to test for the existence of solution
c
      argu2  = cos(glat*pi180)*cos(bb)/
     +            (1.-sin(hsign*glat*pi180)*sin(bb)*sin(alpha))
      if(abs(argu2).gt.1.0) then
      glolim = 999.0
      else
        glolim = acos(argu2)/pi180
      end if
c
c     modify (if nesc.) the plo solution
c
      if((abs(argu1).gt.glolim.and.glat.le.rlat).or.
     +   glat.gt.rlat) then
            plo = pi180*180.0 - plo
      end if
c
c     the solution is symmetric so the direction must be if'ed
c
      if(argu1.lt.0.0) then
         plo = -plo
      end if
c
c     convert the radians to degrees
c
      pla = pla/pi180        
      plo = plo/pi180
c
c     to obtain a rotated value (ie so x-axis in pol.ste. points east)
c     add 270 to longitude
c
      plo=dmod(plo+270.0d0,360.0d0)
c
 2000 continue      

      pla_r4 = pla
      plo_r4 = plo
!     glat_r4 = glat
!     glon_r4 = glon
      rlat_r4 = rlat
      wlon1_r4 = wlon1

      return
      end                                  
c
c     ******************************************************************
c
      subroutine pstoxy(x,y,pla,plo,erad)
c
c     this program convert polar stereographic coordinates to x,y ditto
c     longitude:   0 - 360  ; positive to the east
c     latitude : -90 -  90  ; positive for northern hemisphere
c     it is assumed that the x-axis point towards the east and
c     corresponds to longitude = 0
c
c     tsp 20/06-89
c
c     constants and functions
c            
      fac(pla) = erad*2.0/(1.0+sin(pla*pi180))*cos(pla*pi180)
      xc(pla,plo) = fac(pla)*cos(plo*pi180)
      yc(pla,plo) = fac(pla)*sin(plo*pi180)      
      pi180=3.14159265/180.0
c
c     calculate the coordinates
c
      x = xc(pla,plo)
      y = yc(pla,plo)
c
      return
      end

      subroutine latlon_to_xy_old(glat,glon,rlat,wlon1,erad,x,y)

!     combines getops + pstoxy

!                  o   o   i    i    i     i
      call getops(pla,plo,glat,glon,rlat,wlon1)

      plo = plo - wlon1

!                 o o  i   i    i
      call pstoxy(x,y,pla,plo,erad)

      return
      end

      subroutine latlon_to_xy(glat,glon,erad,x,y)

!     combines getops + pstoxy

      call latlon_to_uv(glat,glon,u,v,istatus)

      call uv_to_xy(u,v,erad,x,y)

      return
      end

      subroutine xy_to_latlon_old(x,y,erad,rlat,wlon1,glat,glon)

!     combines xytops + pstoge

!                 i i  o   o    i
      call xytops(x,y,pla,plo,erad)

!                  i   i    o    o    i    i
      call pstoge(pla,plo,glat,glon,rlat,wlon1)

      return
      end


      subroutine xy_to_latlon(x,y,erad,glat,glon)

!     combines xytops + pstoge

      call xy_to_uv(x,y,erad,u,v)

      call uv_to_latlon(u,v,glat,glon,istatus)

      return
      end


      subroutine xy_to_uv(x,y,erad,u,v)

      use mem_namelist, only: c6_maproj,slat1=>standard_latitude
     1                                 ,slat2=>standard_latitude2
     1                                 ,grid_spacing_m

      include 'trigd.inc'

      real n 

      if(c6_maproj .eq. 'plrstr')then     ! haltiner & williams 1-21
          call get_ps_parms(slat1,slat2,grid_spacing_m,phi0
     1                                 ,grid_spacing_proj_m)
          factor = 1. + sind(phi0)

          u = x / (factor * erad)
          v = y / (factor * erad)

      elseif(c6_maproj .eq. 'lambrt')then 
          call lambert_parms(slat1,slat2,n,s,rconst)

          u = x / (erad * rconst)
          v = y / (erad * rconst)

      elseif(c6_maproj .eq. 'merctr')then ! haltiner & williams 1-8-2
          u = x / erad
          v = y / erad

      elseif(c6_maproj .eq. 'latlon')then 
          u = x / erad
          v = y / erad

      else
          write(6,*)'xy_to_uv - error: invalid map projection '
     1             ,c6_maproj             
          stop

      endif

      return
      end


      subroutine uv_to_xy(u,v,erad,x,y)

      use mem_namelist, only: c6_maproj,slat1=>standard_latitude
     1                                 ,slat2=>standard_latitude2
     1                                 ,grid_spacing_m
      include 'trigd.inc'

      real n 

      if(c6_maproj .eq. 'plrstr')then     ! haltiner & williams 1-21
          call get_ps_parms(slat1,slat2,grid_spacing_m,phi0
     1                                 ,grid_spacing_proj_m)
          factor = 1. + sind(phi0)

          x = u * factor * erad
          y = v * factor * erad

      elseif(c6_maproj .eq. 'lambrt')then 
          call lambert_parms(slat1,slat2,n,s,rconst)

          x = u * (erad * rconst)
          y = v * (erad * rconst)

      elseif(c6_maproj .eq. 'merctr')then ! haltiner & williams 1-8-2
          x = u * erad
          y = v * erad

      elseif(c6_maproj .eq. 'latlon')then 
          x = u * erad
          y = v * erad

      else
          write(6,*)'uv_to_xy - error: invalid map projection '
     1             ,c6_maproj       
          stop

      endif

      return
      end
