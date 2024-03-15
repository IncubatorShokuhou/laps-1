c*****************************************************************
c
c the .swp disc file is read by first reading the header and then
c reading and unpacking the dbz data one row at a time.
c subroutine readhead may be used to open and read the radar
c header.  subroutines passxcmp and unpakdbz2d may be used to read
c and unpack the data respectively.  the readhead, passxcmp and
c unpakdbz2d subroutines are given below(variables beginning w/
c ijklmn are integer*4 and the rest are real unless explicitly
c declared).
c
c the size of the maps in bins is denoted by the variables names
c "imax" and "jmax".  the spacing of the data is denoted by the
c variable names "sx" and "sy".
c
c *****************************************************************

      subroutine readhead(luout, name, lunit, ierr, keyword, fltname,
     +  stmname, radarid, trackname, createtime, ramname, imax, jmax,
     +  kmax, nsweeps, nmosmflag, iunfoldflag, intattflag, ieditflag,
     +  iextra2, iextra3, stime, etime, olat, olon, sx, sy, sz, xdis,
     +  ydis, z0, rot, radaralt, calibco1, calibco2, azmcor, elcor,
     +  thresh, dbzrs, pcor, dcor, rcor, starthorelev, htbb, dbb)
c     opens radar data file and reads radar header data values.
c     paul a. leighton, usdc/noaa/aoml/hrd, 4 jun 1991
      character      keyword*4, fltname*8, stmname*12, radarid*4,
     +               trackname*32, createtime*32, ramname*28,
     +               name*140, jfile*140

c     close previously opened files
      close(lunit)

c     open the file and read the header variables.
      open(lunit, err = 991, file = name, iostat = ierr,
     +  status = 'old', form = 'unformatted') 
   10 read(lunit, iostat = ierr, err = 993) keyword, fltname, stmname,
     +  radarid, trackname, createtime, ramname, imax, jmax, kmax,
     +  nsweeps, nmosmflag, iunfoldflag, intattflag, ieditflag,
     +  iextra2, iextra3, stime, etime, olat, olon, sx, sy, sz,
     +  xdis, ydis, z0, rot, radaralt, calibco1, calibco2,
     +  azmcor, elcor, thresh, dbzrs, pcor, dcor, rcor,
     +  starthorelev, htbb, dbb
   
      return

c     error messages: 
  991 write(luout, '(''error '', i4, '' on file '', a40,
     +  ''trying /hrd/dat/'')') ierr, name 
      jfile = '/hrd/dat/' // name 
      open(lunit, err = 995, file = jfile, iostat = ierr,
     +  status = 'old', form = 'unformatted') 
      goto 10

  993 write(luout, '(//, "***********************", //,
     + "read error on header", //, "**********************")')

  995 write(luout, '(//, "**********************", //, "error ", i4,
     +  " on file ", a40, "not on /hrd/dat/", //,
     +  "**********************")') ierr, jfile

      return
      end

c *****************************************************************

      subroutine passxcmp(lunit, imax, jmax, max_x, max_y, zarray)
c     reads and unpacks composite data values stored two per word in kpac and
c     returns them in zarray for later??? use.  assumes header is previously
c     read.  written by paul a. leighton, usdc/noaa/aoml/hrd, 12 aug 92.
c     luout - printed output lu number for error message.
c     lunit - disc lu number of composite file opened previously in readhead.
c     imax - row dimension of z array,
c          - maximum number of z values per row in kpac
c     jmax - column dimension of z array,
c          - column dimension of stored composite file being read into kpac.
c     kpac - packed array of dbz composite values 
c     zarray - buffer of unpacked dbz values to be returned
      integer*2    kpac(32767)   ! dbz composite values stored 2 per word
      integer*2    zarray(0:max_x-1,0:max_y-1) !row of dbz values to be 
                                               !returned
      data luout/6/
c check array boundaries:
      max=(imax*jmax)/2
      if (max .gt. 32767) then
         write(luout,'(" the boundaries for imax or jmax are to big.")')
         stop
      endif
c read packed dbz values from disc:
      do jx=1, jmax
         j=jmax+1-jx
         n1=((j-1)*imax+2)/2
         n2=(imax+(j-1)*imax+1)/2
         read(lunit,iostat=ierr,err=991) (kpac(i),i=n1,n2)
      end do
c unpack data:
      do j = 0, jmax-1
         do i = 0, imax-1
            call unpakdbz2d(i+1, j+1, kpac, imax, jmax, z) ! unpacks kpac.
            zarray(i,j) = nint(z)  ! store dbz value in output array.
         enddo
      end do
      goto 9999  ! end subroutine

c     error messages: 
  991 write(luout, '(''error on read:'', i4)') ierr
      write(luout, '("aborting program!!!!")')
      stop

 9999 return  ! to calling subroutine
      end

c *****************************************************************

      subroutine unpakdbz2d(i, j, kpac, imax, jmax, z)

      integer*2 kpac(32767)
      integer   i4_kpac

      n = (i + (j-1)*imax + 1) / 2 
      i4_kpac = kpac(n)
      if (mod(i, 2).eq.0)  then
          l = iand(i4_kpac, 255)
      else
          l = ishft(i4_kpac, -8)
      endif
      z = (l - 64) / 2.
      if (l.eq.0)  z = -999.0 
      return
      end
      
c *****************************************************************
