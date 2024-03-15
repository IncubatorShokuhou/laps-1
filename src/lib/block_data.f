      block data

c     routine supmap.f

        common/supmp1/pi,tovpi,dtr,rtd,eps,ov90,con1,con2,part
        common/supmp2/npts,maxlat,minlat,maxlon,minlon,pts(200)
        common/supmp3/polong,cone,rlat,rlon,jgr,ilf,sgn
        common/supmp4/ifst,igo,igold,icross,iout,uold,vold
        common/supmp5/phioc,sino,coso,sinr,cosr,iproj
        common/supmp6/umin,umax,vmin,vmax,ueps,veps
        common/supmp7/phio,phia,igrid,idot,ilts
        common/supmp8/u,v,u1,v1,u2,v2
        common/supmpa/iier
        common/mapcol/mpcol1,mpcol2,mpcol3,mpcol4
        common/mapdas/ldash1,ldash2,ldash3,ldash4
c        data                    !default line intensities and dash patterns
c     1          mpcol1,ldash1   /255,'1777'o/,  !map lines
c     2          mpcol2,ldash2   /128,'1756'o/,  !grid lines
c     3          mpcol3,ldash3   /192,'1777'o/,  !limb lines
c     4          mpcol4,ldash4   /255,'1777'o/   !perimeter
        data                    !default line intensities and dash patterns
     1          mpcol1,ldash1   /255,1023/,  !map lines
     2          mpcol2,ldash2   /128,1006/,  !grid lines
     3          mpcol3,ldash3   /192,1023/,  !limb lines
     4          mpcol4,ldash4   /255,1023/   !perimeter


c       common/supmp1/pi,tovpi,dtr,rtd,eps,ov90,con1,con2,part
c       common/supmp4/ifst,igo,igold,icross,iout,uold,vold
        common/supmp9/ds,di,dsrdi

        data   con1 / 1.00001/
        data   con2 / 179.99999/
        data     di / 16./
        data    dtr / 1.7453292519943e-2/
        data    eps / 1.e-6/
        data   ov90 / 1.11111111111111e-2/
        data     pi / 3.1415926535898/
        data    rtd / 57.295779513082/
        data  tovpi / 0.63661977236758/
        data   uold / 0.0 /
        data   vold / 0.0 /
        data part   / 1.0 /          !size of picture (90% of screen)
  
c     routine xsect.f
        logical l_convert
        common/lapsplot_omega/l_convert
        data l_convert /.true./

c     routine config_satellite_lvd (in file lapsgrid.f).
c       include 'satellite_dims_lvd.inc'
c       include 'satellite_common_lvd.inc'
c       include 'sat_data_static_lvd.inc'
c       data iflag_lvd_common /0/ 

      end

