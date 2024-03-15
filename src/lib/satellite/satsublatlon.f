      subroutine sat_sublatlon(ewc4,ewi4,nsc4,nsi4,f_time,imc,orbat,
     &satsublat,satsublon,istatus)
c
c
c
      implicit none

      include 'instco.inc'

      real*8        orbat(336)
      real          time_50,time50
      real*8        t50_8
      real*8        t
      real*8        f_time
      real*8        satsublat,satsublon
      real          pi,radtodeg

      integer     ewc4,ewi4
      integer     nsc4,nsi4
      integer     instr
      integer     time_spec(2)
      integer     imc
      integer     istatus

      istatus = -1
      instr=1          !1=imager, 2=sounder

      pi=3.141592653589793
      radtodeg=180.0/pi

      call bcd_to_int(orbat(12),time_spec)
      time_50 = time50(time_spec)
      t50_8=time_50
      t = f_time /60. + 7305. * 24. * 60.

      call setcon(instr,nsc4,nsi4,ewc4,ewi4)
      call lmodel(t,t50_8,orbat,imc,satsublat,satsublon)

      write(6,*)'  sat subpoint lat (deg) ',satsublat*radtodeg
      write(6,*)'  sat subpoint lon (deg) ',satsublon*radtodeg

      istatus=1

      return
      end
