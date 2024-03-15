c
c 5-11-99:  j. smart  afwa code modified for inclusion into laps satellite
c                     library to compute meteosat line/pixel values for
c                     given laps lat/lon.
c   4.4.1.2.2.11.7.1 fc01_conv_pts_el_to_met_1 
ca
c   fc01_conv_pts_el_to_met_1 converts latitude/longitude in radians to
c   pixel scanline number using algorithm described in annex "e"
c   to the meteosat-5 calibration report (april 1991). prepared by the
c   meteosat exploitation project at european space operations center.
ce
c---------------------  arguments passed  ---------------------------
cca
c
c      name                    i/o        description
c
c      latitude                 i   input latitude
c
c      longitude                i   input longitude
c
c      x_coord                  o   output x coordinate
c
c      y_coord                  o   output y coordinate
c
c      cstatus                  o   error code c
ccae
c-------------------  software unit history  ------------------------
ch
c
c       06dec95 : created, scr95124, sfl;
c       26feb97 : renamed rs to rs1 to avoid conflict with ingest
c                  include file, sdhsup3, kla;
c
c      25mar97 : sdhsu   placed key in comments for
c                        sdd sec 4 generation; rfb
che

c--------------------------------------------------------------------
c
c     procedure fc01_conv_pts_el_to_met_1 is
c
cp    begin
ce    end fc01_conv_pts_el_to_met_1 ;
c
c--------------------------------------------------------------------

      subroutine fc01_conv_pts_el_to_met_1(golonsbp,
     &                          stop_pix,fsci,
     &                          start_line,decimat,
     &                          ct,
     &                          latitude,
     &                          longitude,
     &                          x_coord,
     &                          y_coord,
     &                          cstatus)

      implicit none

c     include 'data_type_code_tables.inc'
c     include 'message_data_variables.inc'
c     include 'is1_ingest_parameters_5.inc'

      character ct*(*)

      real   latitude,
     &       longitude,
     &       x_coord,
     &       y_coord

      integer   cstatus

      real   h
      parameter (h = 35785.845)
      real   re
      parameter (re = 6378.155)
      real   a
      parameter (a = 1.0/297.0)
      real   rp
      parameter (rp = re/(1.0+a))
      real   pi
      parameter (pi = 3.141592653)
      real   cdr
      parameter (cdr = pi/180.0)
      real   crd 
      parameter (crd = 180.0/pi)
      integer   lpsi2
      parameter (lpsi2 = 1)
      real   deltax
      parameter (deltax = 18.0/2500.0)
      real   deltay 
      parameter (deltay = 18.0/2500.0)

      real   xfi,
     &       xla,
     &       rom,
     &       y,
     &       r1,
     &       r2,
     &       rs1,
     &       xr,
     &       yr,
     &       xt,
     &       yt,
     &       zt,
     &       px,
     &       py,
     &       teta
      real   radjust
      real      golonsbp
      integer   start_line
      integer   stop_pix
      integer   fsci
      integer   decimat
      
c     include 'data_areas_imd_hdr_offsets.inc'
c     include 'data_areas_imd_hdr_data_names.inc'

c     common /block_13/ image_list_char

c-----------------------------------------------------------
c--   begin
c-----------------------------------------------------------
c#
      cstatus = 0
      radjust=1.0
      if(ct.eq.'vis')radjust=2.0
      xfi = latitude * cdr
      xla = (longitude - golonsbp) * cdr
      rom = (re * rp)/sqrt(rp**2 * cos(xfi)**2 + re**2 * sin(xfi)**2)
      y = sqrt(h**2 + rom**2 - 2 * h * rom * cos(xfi) * cos(xla))
      r1 = y**2 + rom**2
      r2 = h**2

      if (r1 .gt. r2) then
         x_coord = 9999.
         y_coord = 9999.
         cstatus = 1
      else
         rs1 = re + h
         teta = atan((rp / re) * tan(xfi))
         xt = re * cos(teta) * cos(xla)
         yt = re * cos(teta) * sin(xla)
         zt = rp * sin(teta)

         px = atan(yt / (xt-rs1))

         py = atan(zt * (-1.0 / (xt-rs1))
     &                * cos(px))

         px = px * crd
         py = py * crd

         xr = px / (deltax * lpsi2)
         yr = py / (deltay * lpsi2)

         if(xr .ge. 0.0) xr = int(px / (deltax * lpsi2)) + 0.5
         if(xr .lt. 0.0) xr = int(px / (deltax * lpsi2)) - 0.5
         if(yr .ge. 0.0) yr = int(py / (deltay * lpsi2)) + 0.5
         if(yr .lt. 0.0) yr = int(py / (deltay * lpsi2)) - 0.5

         x_coord = (xr + 1250.0)*radjust
         y_coord = (yr + 1251.0)*radjust
c
c compute file relative x/y coordinate positions.
c
         x_coord = x_coord/decimat
         x_coord = stop_pix-x_coord

         y_coord = y_coord - fsci
         y_coord = y_coord/decimat + 0.5
         y_coord = y_coord - start_line

c        y_coord = ((y_coord-(fsci/radjust))/decimat)+0.5-start_line

      end if

99999 return
      end
