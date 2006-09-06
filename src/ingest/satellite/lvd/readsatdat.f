      subroutine readsatdat(csat_id,
     &                      csat_type,
     &                      c_dir_path,
     &                      c_fname_data,
     &                      c_type,
     &                      ntf,ntm,
     &                      maxchannels,
     &                      max_files,
     &                      maximages,
     &                      nlineswv,nelemwv,
     &                      nlinesir,nelemir,
     &                      nlinesvis,nelemvis,
     &                      image_67,
     &                      image_ir,image_12,
     &                      image_39,image_vis,
     &                      istatus)
c
c
c
      implicit none

      include 'netcdf.inc'

      integer i,j,n
      integer ntf
      integer nlinesir,nelemir
      integer nlinesvis,nelemvis
      integer nlineswv,nelemwv
      integer max_files,maximages,maxchannels
      integer ntm(max_files)
      integer istatus
      integer istatus_ir
      integer istatus_vis
      integer istatus_wv
      integer wstatus
      integer istat
      integer ispec

      real*4    image_ir(nelemir,nlinesir,maximages)
      real*4    image_12(nelemir,nlinesir,maximages)
      real*4    image_39(nelemir,nlinesir,maximages)
      real*4    image_67(nelemwv,nlineswv,maximages)
      real*4    image_vis(nelemvis,nlinesvis,maximages)
      real*4    wv_image(nelemwv,nlineswv)
      real*4    ir_image(nelemir,nlinesir)
      real*4    vis_image(nelemvis,nlinesvis)

      character c_fname_data(max_files)*9
      character c_type(maxchannels,max_files)*3
      character c_dir_path(maxchannels)*200
      character c_filename*200
      character c_wfo_fname*13
      character fname9_to_wfo_fname13*13
      character csat_id*6
      character csat_type*3
c
      INTEGER   Nx
      INTEGER   Ny 
      INTEGER   record
      INTEGER   ncid
      INTEGER   rcode
      INTEGER   ivalidTime
      doubleprecision validtime
      REAL*4    dummy
      REAL*4    la1_vis,la1_ir,la1_wv
      REAL*4    lo1_vis,lo1_ir,lo1_wv
      REAL*4    dx_vis,dx_ir,dx_wv
      REAL*4    dy_vis,dy_ir,dy_wv
      REAL*4    latin_vis,latin_ir,latin_wv
      REAL*4    lov_vis,lov_ir,lov_wv

      istatus=1

      do i=1,ntf
         do j=1,ntm(i)

            call lvd_file_specifier(c_type(j,i),ispec,istat)

            if(csat_type.eq.'wfo')then
               n=index(c_dir_path(ispec),' ')-1
               c_wfo_fname = fname9_to_wfo_fname13(c_fname_data(i))
               c_filename=c_dir_path(ispec)(1:n)//c_wfo_fname

               n=index(c_filename,' ')-1
               print*,'opening ',c_filename(1:n)
               rcode = NF_OPEN(c_filename,NF_NOWRITE,ncid)

               if(rcode.ne.NF_NOERR) then
                  print *, NF_STRERROR(rcode)
                  print *,'NF_OPEN ',c_filename(1:n)
                  istatus=-1
                  return
               endif

               print*,'calling get_attribute_wfo ',c_type(i,j)
               call get_attribute_wfo(ncid,dummy,dummy,dummy,
     &dummy,dummy,dummy,dummy,dummy,dummy,dummy,nx,ny,wstatus)

               if(wstatus .lt. 0)then
                  print*,'No attributes: get_attribute_wfo ',c_type(j,i)
                  return
               endif
               record=1
               n=index(c_filename,' ')
               write(6,*)'Reading: ',c_filename(1:n)
            else
               n=index(c_dir_path(1),' ')-1
               c_filename=c_dir_path(1)(1:n)//c_fname_data(i)//
     &'_'//c_type(j,i)
               n=index(c_filename,' ')
               write(6,*)'Reading: ',c_filename(1:n)

               rcode=NF_OPEN(c_filename,NF_NOWRITE,NCID)
               if(rcode.ne.nf_noerr) return

               call get_cdf_dims(ncid,record,nx,ny,istatus)
               record = 1
               if(istatus.eq.1)then
                  print*,'Error reading cdf dimensions'
                  return
               endif
               istatus =1

            endif

            if(ispec.ne.1.and.ispec.ne.3)then    !check for visible and water vapor

               call readcdf(csat_id,
     &                    csat_type,
     &                    c_type(j,i),
     &                    nx,ny,
     &                    record,
     &                    nelemir,nlinesir,
     &                    ir_image,
     &                    la1_ir,lo1_ir,
     &                    Dx_ir,Dy_ir,
     &                    Latin_ir,Lov_ir,
     &                    ivalidTime,
     &                    ncid,
     &                    istatus_ir)
               if(istatus_ir.eq.1)then
                  write(6,*)'Successful'
               else
                  Write(6,*)'NOT Successful!!'
                  istatus=-1
                  goto 125
               endif


cisido          
		write(6,*)'la1_ir',la1_ir,'lo1_ir',lo1_ir,
     &                     'Dx_ir',Dx_ir,'Dy_ir',Dy_ir,
     &                     'Latin_ir',Latin_ir,'Lov_ir',Lov_ir 
cisid

c
               if(ispec .eq. 2)then 
                  call move(ir_image,image_39(1,1,i),nelemir,nlinesir)
               elseif(ispec .eq. 4)then
                  call move(ir_image,image_ir(1,1,i),nelemir,nlinesir)
               elseif(ispec .eq. 5)then
                  call move(ir_image,image_12(1,1,i),nelemir,nlinesir)
               endif
c
            elseif(ispec.eq.1)then

               call readcdf(csat_id,
     &                    csat_type,
     &                    c_type(j,i),
     &                    nx,ny,
     &                    record,
     &                    nelemvis,nlinesvis,
     &                    vis_image,
     &                    la1_vis,lo1_vis,
     &                    Dx_vis,Dy_vis,
     &                    Latin_vis,Lov_vis,
     &                    ivalidTime,
     &                    ncid,
     &                    istatus_vis)
               if(istatus_vis.eq.1)then
                  write(6,*)'Successful'
                  call move(vis_image,image_vis(1,1,i),nelemvis,
     &                      nlinesvis)
               else
                  Write(6,*)'NOT Successful!!'
                  istatus=-1
               endif

cisido
                write(6,*)'la1_vis',la1_vis,'lo1_vis',lo1_vis,
     &                     'Dx_vis',Dx_vis,'Dy_vis',Dy_vis,
     &                     'Latin_vis',Latin_vis,'Lov_vis',Lov_vis
cisid





c
c load vis attributes
c
            elseif(ispec.eq.3)then

               call readcdf(csat_id,
     &                    csat_type,
     &                    c_type(j,i),
     &                    nx,ny,
     &                    record,
     &                    nelemwv,nlineswv,
     &                    wv_image,
     &                    la1_wv,lo1_wv,
     &                    Dx_wv,Dy_wv,
     &                    Latin_wv,Lov_wv,
     &                    ivalidTime,
     &                    ncid,
     &                    istatus_wv)
               if(istatus_wv.eq.1)then
                  write(6,*)'Successful'
                  call move(wv_image,image_67(1,1,i),nelemwv,
     &                      nlineswv)
               else
                  Write(6,*)'NOT Successful!!'
                  istatus=-1
               endif

cisido
                write(6,*)'la1_wv',la1_wv,'lo1_wv',lo1_wv,
     &                     'Dx_wv',Dx_wv,'Dy_wv',Dy_wv,
     &                     'Latin_wv',Latin_wv,'Lov_wv',Lov_wv
cisid





c
c load vis attributes
c
            endif
125      enddo
      enddo
c
      return
      end
