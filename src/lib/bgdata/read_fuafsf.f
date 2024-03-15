      subroutine read_fuafsf_cdf(fullname, x, y, z, 
     +     ht, level, om, sh, t3, u3, v3,
     +     usfc, vsfc, tsfc, dsfc, psfc, mslp, zsfc, r01, rto,
     +     lmr, llr, s8a, swi, tpw,
     +     istatus)
c
      implicit none
      include 'netcdf.inc'
      character*(*) fullname
      character*200 cfname_int
      character*80 c8_proj
      integer x, y, z, nf_fid, nf_vid, nf_status
      integer lname,lname_in,lend,l,lenc8
      integer istatus

      real ht( x,  y,  z, 1)
     +,    om( x,  y,  z, 1)
     +,    sh( x,  y,  z, 1)
     +,    t3( x,  y,  z, 1)
     +,    u3( x,  y,  z, 1)
     +,    v3( x,  y,  z, 1)
     +,    w3( x,  y,  z, 1)

      real level(z)

      real tsfc( x,  y)
     +,    psfc( x,  y)
     +,    dsfc( x,  y)
     +,    mslp( x,  y)
     +,    usfc( x,  y)
     +,    vsfc( x,  y)
     +,    zsfc( x,  y)
     +,    r01( x,  y)
     +,    rto( x,  y)
     +,    lmr( x,  y)
     +,    llr( x,  y)
     +,    s8a( x,  y)
     +,    swi( x,  y)
     +,    tpw( x,  y)

      call get_c8_project(c8_proj,istatus)

      istatus = 1 ! good status

      call downcase(c8_proj,c8_proj)
      call s_len(c8_proj,lenc8)

      call s_len(fullname,lname_in)

      if(x*y .eq. 0)then
          write(6,*)' error in read_fuafsf_cdf input x/y ',x,y
          istatus = 0
          return
      endif

      if(z .le. 1)then
          write(6,*)' skip fua read, lname_in,x,y = ',lname_in,x,y
          cfname_int=fullname(1:lname_in)//'.fsf'
          goto100
      endif
c
c  open netcdf file for reading
c
      cfname_int=fullname(1:lname_in)//'.fua'
      call s_len(cfname_int,lname)
      print*,'open/read ',cfname_int(1:lname)
      nf_status = nf_open(cfname_int,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ',cfname_int(1:lname)
        istatus = 0; return
      endif
c
c     variable        netcdf long name
c      ht           "laps fcst height"
c
      nf_status = nf_inq_varid(nf_fid,'ht',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ht'
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,ht)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ht'
        istatus = -1            
      endif

c     if(c8_proj(1:lenc8).eq.'airdrop')
      call swap_array_k(ht,x,y,z)
c
c     variable        netcdf long name
c      level        "level of data"
c
      nf_status = nf_inq_varid(nf_fid,'level',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var level'
        istatus = -1             
      else ! found id
        nf_status = nf_get_var_real(nf_fid,nf_vid,level)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var level'
          istatus = -1           
        else ! found data
c         if(c8_proj(1:lenc8).eq.'airdrop')
          call swap_array_k(level,1,1,z)
        endif
      endif
c
c     variable        netcdf long name
c      om           "laps fcst omega wind component"
c
      nf_status = nf_inq_varid(nf_fid,'om',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var om'
        istatus = -1            
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,om)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var om'
        istatus = -1          
      endif

c     if(c8_proj(1:lenc8).eq.'airdrop')
      call swap_array_k(om,x,y,z)
c
c     variable        netcdf long name
c      sh           "laps fcst specific humidity"
c
      nf_status = nf_inq_varid(nf_fid,'sh',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var sh'
        istatus = -1          
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,sh)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var sh'
        istatus = -1                
      endif

c     if(c8_proj(1:lenc8).eq.'airdrop')
      call swap_array_k(sh,x,y,z)
c
c     variable        netcdf long name
c      t3           "laps fcst temperature"
c
      nf_status = nf_inq_varid(nf_fid,'t3',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var t3'
        istatus = -1            
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,t3)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var t3'
        istatus = -1           
      endif

c     if(c8_proj(1:lenc8).eq.'airdrop')
      call swap_array_k(t3,x,y,z)
c
c     variable        netcdf long name
c      u3           "laps fcst u wind component"
c
      nf_status = nf_inq_varid(nf_fid,'u3',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var u3'
        istatus = -1          
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,u3)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var u3'
        istatus = -1           
      endif

c     if(c8_proj(1:lenc8).eq.'airdrop')
      call swap_array_k(u3,x,y,z)
c
c     variable        netcdf long name
c      v3           "laps fcst v wind component"
c
      nf_status = nf_inq_varid(nf_fid,'v3',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var v3'
        istatus = -1              
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,v3)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var v3'
        istatus = -1           
      endif

c     if(c8_proj(1:lenc8).eq.'airdrop')
      call swap_array_k(v3,x,y,z)

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
        istatus = -1           
      endif
c
c now fsf file:
c
      call get_directory_length(cfname_int,lend)
      search_dir: do l=lend,1,-1
         if(cfname_int(l:l).eq.'f')then
            if(cfname_int(l:l+2).eq.'fua')then
               exit search_dir
            endif
         endif
      enddo search_dir

      if(l.le.1)then
         print*,'didnt determine location of fua in'
         print*,'string for fsf. return with no data'
         print*,'read_fuafsf.f: abort'
         istatus = 0; return
      endif

      cfname_int(l:l+2)='fsf'

 100  continue

      cfname_int=cfname_int(1:lname_in)//'.fsf'

      call s_len(cfname_int,lname)
      print*,'open/read ',cfname_int(1:lname)
      nf_status = nf_open(cfname_int,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ',cfname_int(1:lname)
        istatus = 0; return
      endif

      write(6,*)' reading usf'
c
c     variable        netcdf long name
c      usf          "laps fcst sfc u wind component"
c
      nf_status = nf_inq_varid(nf_fid,'usf',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var usf'
        istatus = -1          
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,usfc)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var usf'
        istatus = -1           
      endif
c
c     variable        netcdf long name
c      vsf          "laps fcst sfc v wind component"
c
      nf_status = nf_inq_varid(nf_fid,'vsf',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vsf'
        istatus = -1          
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,vsfc)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vsf'
        istatus = -1          
      endif
c
c     variable        netcdf long name
c      vsf          "laps fcst sfc temperature"
c
      nf_status = nf_inq_varid(nf_fid,'tsf',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tsf'
        istatus = -1          
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,tsfc)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tsf'
        istatus = -1           
      endif
c
c     variable        netcdf long name
c      vsf          "laps fcst sfc dew point temperature"
c
      nf_status = nf_inq_varid(nf_fid,'dsf',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dsf'
        istatus = -1           
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,dsfc)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dsfc'
        istatus = -1            
      endif
c
c     variable        netcdf long name
c      vsf          "laps fcst mean sea level pressure"
c
      nf_status = nf_inq_varid(nf_fid,'slp',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var slp'
        istatus = -1            
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,mslp)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var mslp'
        istatus = -1           
      endif
c
c     variable        netcdf long name
c      vsf          "laps fcst sfc pressure"
c
      nf_status = nf_inq_varid(nf_fid,'psf',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var psf'
        istatus = -1          
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,psfc)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var psfc'
        istatus = -1           
      endif

      write(6,*)' reading zsfc'
cc
c     variable        netcdf long name
c      ter          "laps sfc terrain"
c
      nf_status = nf_inq_varid(nf_fid,'ter',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ter'
        istatus = -1           
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,zsfc)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var zsfc'
        istatus = -1            
      endif

      write(6,*)' reading r01'
c
c     variable        netcdf long name
c      r01          "laps incremental precip"
c
      nf_status = nf_inq_varid(nf_fid,'r01',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var r01'
        istatus = -1          
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,r01)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var r01'
        istatus = -1          
      endif

      write(6,*)' reading rto'
c
c     variable        netcdf long name
c      rto          "laps storm-total precip"
c
      nf_status = nf_inq_varid(nf_fid,'rto',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rto'
        istatus = -1          
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,rto)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rto'
        istatus = -1          
      endif

      write(6,*)' reading lmr'
c
c     variable        netcdf long name
c      lmr          "laps column max reflectivity"
c
      nf_status = nf_inq_varid(nf_fid,'lmr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lmr'
        istatus = -1          
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,lmr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lmr'
        istatus = -1          
      endif

      write(6,*)' reading swi'
c
c     variable        netcdf long name
c      swi          "laps ghi"
c
      nf_status = nf_inq_varid(nf_fid,'swi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var swi'
        istatus = -1          
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,swi)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var swi'
        istatus = -1          
      endif

      write(6,*)' reading s8a'
c
c     variable        netcdf long name
c      s8a          "laps ghi"
c
      nf_status = nf_inq_varid(nf_fid,'s8a',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var s8a'
        istatus = -1          
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,s8a)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var s8a'
        istatus = -1          
      endif

      write(6,*)' reading tpw'
c
c     variable        netcdf long name
c      tpw          "laps tpw"
c
      nf_status = nf_inq_varid(nf_fid,'tpw',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var tpw'
        istatus = -1          
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,tpw)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var tpw'
          istatus = -1          
        endif
      endif

      write(6,*)' closing file ',nf_fid,z

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
        return
      endif

      write(6,*)' returning from read_fuafsf_cdf...'

      return
      end
c
c =================================================
c
      subroutine swap_array_k(a1,nx,ny,nz)
c
      implicit none
      integer nx,ny,nz,k
      real  a1(nx,ny,nz)
      real, allocatable:: a2(:,:,:)

      allocate (a2(nx,ny,nz))

      do k=1,nz
         a2(:,:,k)=a1(:,:,nz-k+1)
      enddo
      a1=a2
      deallocate (a2)
      return
      end

