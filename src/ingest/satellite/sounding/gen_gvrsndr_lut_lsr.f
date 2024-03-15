       subroutine gen_gvrsndr_lut_lsr(c_filename_sat,nlines,nelems,
     &wavelength,resolution_x,resolution_y,rsndr_res_km,start_pix,
     &start_line,end_pix,end_line,ewcycles,ewincs,nscycles,nsincs,
     &f_time,orbat,nch,nxl,nyl,lat,lon,rpix,rline,pct_req_lsr,
     &istatus)
c
c
c
      implicit none
 
      include       'instco.inc'

      integer     nlines,nelems
      integer     nxl,nyl
      integer     nch
c
      real*8        wavelength(nch)
      real        scalingbias(nlines,nch)
      real        scalinggain(nlines,nch)
      real        lat(nxl,nyl)
      real        lon(nxl,nyl)
      real        rline(nxl,nyl)
      real        rpix(nxl,nyl)
      real        pi
      real        rd2dg
      real        rsndr_res_m
      real        rsndr_res_km
      real        r_missing_data
      real        rnx,rny,rnp
      real        pct_covered
      real        pct_req_lsr
      real        time_50,time50

      real*8        r8lat(2),r8lon(2)
      real*8        r8sl,r8sp,r8el,r8ep
      real*8        elev,scan
      real*8        rl
      real*8        rp
      real*8        orbat(336)
      real*8        t50_8
      real*8        t
      real*8        f_time
      real*8        satsublat,satsublon
      real*8        evln,scpx
      real*8        rl_div
      real*8        rp_div
c
      integer     ewcycles,ewincs
      integer     nscycles,nsincs
      integer     start_line
      integer     start_pix
      integer     end_line
      integer     end_pix
      integer     resolution_x
      integer     resolution_y

      integer     i1,j1
      integer     cstatus
      integer     istatus
      integer     mstatus
      integer     nstatus
      integer     wstatus
      integer     ierr
      integer     instr
      integer     i,j,k,n
      integer     n2,nn
      integer     time_spec(2)
      integer     imc
      integer     npoints_out
      integer     i4time_nearest_sat
      integer     isat
      integer     nrl,nrp

      character     filename_parm*255
      character     c_filename_sat*255
      character*255 c_sounding_path
      character     table_path*200
      character*9   c_ftimesat
      character     ctype*4
      character     c_imc*4
c ---------------------------------------------------------
c read static information for gvarimage navigation
c
c     if(istatus.ne.1)then
c        print*,'error returned get_r_missing_data'
c        return
c     endif

      istatus = -1

      rsndr_res_m=(resolution_y+resolution_x)/2.0
      rl_div = resolution_y/8.
      rp_div = resolution_x/8.
      if(resolution_x.eq.0)then
         rp_div=1.
      endif
      if(resolution_y.eq.0)then
         rl_div=1.
      endif

c     read(c_imc(2:4),111)imc
c111   format(i3)
c     if(imc.ne.0)imc=1

      imc=1
      instr=2          !1=imager, 2=sounder
      pi=3.141592653589793
      rd2dg=180.0/pi

      call bcd_to_int(orbat(12),time_spec)
      time_50 = time50(time_spec)
      t50_8=time_50
      t = f_time /60. + 7305. * 24. * 60.

      call setcon(instr,nscycles,nsincs,ewcycles,ewincs)
      call lmodel(t,t50_8,orbat,imc,satsublat,satsublon)

      write(6,*)'sat subpoint lat (deg) ',satsublat*rd2dg
      write(6,*)'sat subpoint lon (deg) ',satsublon*rd2dg
      print*,'-------------------------------------------'
      print*
      print*,'this satellite sector lat/lon corner points'
      print*,'-------------------------------------------'
      r8sl=float(start_line)
      r8el=float(end_line)
      r8sp=float(start_pix)
      r8ep=float(end_pix)
      elev = evln(instr,r8sl)
      scan = scpx(instr,r8sp)
      call lpoint(elev,scan,r8lat(1),r8lon(1),ierr)
      elev = evln(instr,r8sl)
      scan = scpx(instr,r8ep)
      call lpoint(elev,scan,r8lat(2),r8lon(2),ierr)

      r8lat(1)=r8lat(1)*rd2dg
      r8lon(1)=r8lon(1)*rd2dg
      r8lat(2)=r8lat(2)*rd2dg
      r8lon(2)=r8lon(2)*rd2dg
      print*,'      nw                  ne'
      write(6,101)r8lat(1),r8lon(1),r8lat(2),r8lon(2)
101   format(5x,2(f8.2),2(f9.2))

      elev = evln(instr,r8el)
      scan = scpx(instr,r8sp)
      call lpoint(elev,scan,r8lat(1),r8lon(1),ierr)
      elev = evln(instr,r8el)
      scan = scpx(instr,r8ep)
      call lpoint(elev,scan,r8lat(2),r8lon(2),ierr)
      r8lat(1)=r8lat(1)*rd2dg
      r8lon(1)=r8lon(1)*rd2dg
      r8lat(2)=r8lat(2)*rd2dg
      r8lon(2)=r8lon(2)*rd2dg
      print*,'      sw                  se'
      write(6,101)r8lat(1),r8lon(1),r8lat(2),r8lon(2)

      cstatus = 0
      npoints_out = 0
      do j=1,nyl
      do k=1,nxl
 
         r8lat(1)=lat(k,j)*pi/180.d0
         r8lon(1)=lon(k,j)*pi/180.d0

         call gpoint(r8lat(1),r8lon(1),elev,scan,ierr)

         if(ierr.ne.0)then
c           write(6,*)'error computing elev/scan in gpoint from'
c           write(6,*)'lat/lon ', lat(k,j),lon(k,j)
            rline(k,j)=-2.
            rpix(k,j)=-2.
            cstatus=cstatus-1
         else

            call evsc2l(instr,elev,scan,rl,rp)

            nrl=nint(rl)
            nrp=nint(rp)
            if( (nrl.gt.0.0.and.nrp.gt.0.0).and.
     &          (nrl.ge.start_line.and.nrp.ge.start_pix).and.
     &          (nrl.le.end_line.and.nrp.le.end_pix) )then

               rline(k,j)=(rl-start_line+rl_div)*rl_div
               rpix(k,j)= (rp-start_pix+rp_div )*rp_div
               i1=k
               j1=j
            else

               npoints_out = npoints_out+1
               rline(k,j)=r_missing_data
               rpix(k,j) =r_missing_data

            endif

         endif

      enddo
      enddo

      if(cstatus.lt.0)then

        write(6,*)'warning! some rl/rp values not computed '
        write(6,*)'for lut ',ctype,' status = ',cstatus

      endif

      if(npoints_out .gt. 0)then

         rnx=float(nxl)
         rny=float(nyl)
         rnp=float(npoints_out)
         write(6,*)'warning! some rl/rp values out of domain'
         pct_covered=((rnx*rny)-rnp)/(rnx*rny)
         if(pct_covered.lt.pct_req_lsr)then
            print*,'exceeded domain cover namelist parameter',
     +' pct_req_lsr: % req/% covered:',pct_req_lsr,'/',pct_covered
            return
         else
            write(6,*)'percent of domain covered = ',pct_covered
         endif
      else
         write(6,*)'100% percent of domain covered '
      endif
c
      if(.false.)then
         do j = 1,nyl,20
          do k = 1,nxl,20
           write(6,*)'i,j,rpix,rline: ',k,j,rpix(k,j),rline(k,j)
          enddo
         enddo
      endif
c
c compute image resolution in meters
c
      if(rpix(nxl/2,nyl/2).ne.r_missing_data.and.
     &   rline(nxl/2,nyl/2).ne.r_missing_data)then
         i1=nxl/2
         j1=nyl/2
      endif

      if(i1.gt.0 .and. i1.lt.nelems .and.
     +   j1.gt.0 .and. j1.le.nlines)then
         call compute_sat_res_m(rp_div,rl_div,
     &rpix(i1,j1),rline(i1,j1),start_pix,start_line,instr,
     &rsndr_res_m,istatus)
      else
         rsndr_res_m=(resolution_y+resolution_x)/2.0 
      endif
      rsndr_res_km=rsndr_res_m/1000.

      istatus = 1
      goto 1000

1000  return
      end
