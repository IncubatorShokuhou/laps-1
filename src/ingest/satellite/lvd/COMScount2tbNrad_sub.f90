
subroutine comscount2tbnrad_sub(path_to_raw_sat,max_files &
                               ,n_lines_ir,n_pixels_ir &
                               ,n_vis_lines,n_vis_elem &
                               ,n_wv_lines,n_wv_elem &
                               ,r_missing_data &
                               ,image_lat_ir,image_lon_ir &
                               ,ir1_tb_out,vis_rad_out &
                               ,ir2_tb_out,swir_tb_out,wv_tb_out &
                               ,i4time_data,istatus)


!****************************************************************************
!*this fortran code converts from coms count value (binary)                 *
!*                                     to temperature and radiance (ascii). *
!*                                                                          *
!*input                                                                     *
!*nx and ny : array size of coms data                                       *
!* fd :  1375 * 1429                                                        * 
!* enh : 1934 * 1544                                                        *
!*input units number 10 to 11 : temperature and radiance conversion tables  *
!*input units number 20 to 24 : coms count data                             *
!*each unit is assigned chanel in ascending order(ir1, ir2, wv, swir, vis). *
!*example code uses jan 1, 2016, 02:45 utc data.                            *
!*input units number 40 (optional) : latitude and longitude conversion table*
!*                                                                          *
!*output                                                                    *
!*chanel_tb.dat : temperature of each chanel obtained from count value      *
!*chanel_rad.dat : radiance of each chanel obtained from count value        *
!*their unit numbers are 30 to 39.                                          *
!*                                                                          *
!*optional commands                                                         *
!*there are several lines for additional analysis.                          *
!*line number 42, 43 and 154 to 162 are commented.                          *
!*if users want to know latitude and longitude of each cell,                *
!*                                     relevant to same row of latticepoint.*
!*each column of latticepoint denotes column number and row number of cell. *
!*and each column of latlon represents latitude and longitude of cell,      *
!*                                     relevant to same row of latticepoint.*
!*                                                                          *
!*author : gyeong-gyeun ha                                                  *
!*affiliation : national meteorological satellite center (of kma)           *
!*date: nov.2016                                                            *
!****************************************************************************

integer :: i, j
integer, parameter :: nx=1934, ny=1544, nbyte=2, numofch=5, sizeofctab=1024
integer*2, dimension(nx,ny) :: ir1, ir2, wv, swir, vis
real, dimension(sizeofctab,numofch) :: tb, rad
real, dimension(nx,ny) :: ir1_tb, ir2_tb, wv_tb, swir_tb, vis_tb, ir1_rad, ir2_rad, wv_rad, swir_rad, vis_rad
integer, dimension(nx*ny,2) :: latticepoint
real, dimension(nx*ny,2) :: latlon
real, dimension(nx,ny) :: image_lat_ir,image_lon_ir
character*12 a12time_coms
character*13 a13time,fname9_to_wfo_fname13
integer cvt_wfo_fname13_i4time

! inputs
character*200 path_to_raw_sat

! outputs
integer i4time_data(max_files)
real, dimension(n_pixels_ir,n_lines_ir) :: ir1_tb_out
real, dimension(n_pixels_ir,n_lines_ir) :: ir2_tb_out
real, dimension(n_pixels_ir,n_lines_ir) :: wv_tb_out
real, dimension(n_pixels_ir,n_lines_ir) :: swir_tb_out
real, dimension(n_pixels_ir,n_lines_ir) :: vis_rad_out

write(6,*)' subroutine comscounttbnrad ',nx,ny,n_pixels_ir,n_lines_ir

if(.true.)then
  call get_systime(i4time,a9_time,istatus)
  i4time_data(1) = i4time
  a13time = fname9_to_wfo_fname13(a9_time)
  a12time_coms = a13time(1:8)//a13time(10:13)
else
  a12time_coms = '201701010000'
  a13time = a12time_coms(1:8)//'_'//a12time_coms(9:12)
  i4time_data(1)=cvt_wfo_fname13_i4time(a13time)
endif

write(6,*)' read binary data for ',a12time_coms,' ',a13time,i4time_data(1)

! read binary data

open(20,file=trim(path_to_raw_sat)//'/coms_le1b_ir1_ch1_cn_'//a12time_coms//'.bin',& 
        form='unformatted', recl=nx*ny*nbyte, access='direct',&
        status='old', convert='big_endian')
open(21,file=trim(path_to_raw_sat)//'/coms_le1b_ir2_ch2_cn_'//a12time_coms//'.bin',& 
        form='unformatted', recl=nx*ny*nbyte, access='direct',&
        status='old', convert='big_endian')
open(22,file=trim(path_to_raw_sat)//'/coms_le1b_wv_ch3_cn_'//a12time_coms//'.bin',& 
        form='unformatted', recl=nx*ny*nbyte, access='direct',&
        status='old', convert='big_endian')
open(23,file=trim(path_to_raw_sat)//'/coms_le1b_swir_ch4_cn_'//a12time_coms//'.bin',& 
        form='unformatted', recl=nx*ny*nbyte, access='direct',&
        status='old', convert='big_endian')
open(24,file=trim(path_to_raw_sat)//'/coms_le1b_vis_ch5_cn_'//a12time_coms//'.bin',& 
        form='unformatted', recl=nx*ny*nbyte, access='direct',&
        status='old', convert='big_endian')

read(20, rec=1) ((ir1(i,j), i=1,nx), j=1,ny)
read(21, rec=1) ((ir2(i,j), i=1,nx), j=1,ny)
read(22, rec=1) ((wv(i,j), i=1,nx), j=1,ny)
read(23, rec=1) ((swir(i,j), i=1,nx), j=1,ny)
read(24, rec=1) ((vis(i,j), i=1,nx), j=1,ny)
close(20)
close(21)
close(22)
close(23)
close(24)


write(6,*)' read temperature and radiance table'

! read temperature and radiance table

open(10,file=trim(path_to_raw_sat)//'/coms_mi_conversion_table_tb.dat')
open(11,file=trim(path_to_raw_sat)//'/coms_mi_conversion_table_radiance.dat')

do i=1, sizeofctab
	read(10,*) tb(i,1), tb(i,2), tb(i,3), tb(i,4), tb(i,5)
	read(11,*) rad(i,1), rad(i,2), rad(i,3), rad(i,4), rad(i,5)
enddo 

close(10)
close(11)

write(6,*)' convert count value to tb and albedo'

! convert count value to tb and albedo

do i=1,nx
	do j=1,ny
		ir1_tb(i,j)=tb(ir1(i,j)+1,1)      ! add 1 to assign array index start from 1 in fortran 
		ir1_rad(i,j)=rad(ir1(i,j)+1,1)

		ir2_tb(i,j)=tb(ir2(i,j)+1,2)
		ir2_rad(i,j)=rad(ir2(i,j)+1,2)

		wv_tb(i,j)=tb(wv(i,j)+1,3)
		wv_rad(i,j)=rad(wv(i,j)+1,3)

		swir_tb(i,j)=tb(swir(i,j)+1,4)
		swir_rad(i,j)=rad(swir(i,j)+1,4)

		vis_tb(i,j)=tb(vis(i,j)+1,5)
		vis_rad(i,j)=rad(vis(i,j)+1,5)
	enddo
enddo

iwrite = 0 

if(iwrite .eq. 1)then

! save output in ascii (example) 

    open(30,file='./ir1_tb.dat')
    open(31,file='./ir1_rad.dat')
    open(32,file='./ir2_tb.dat')
    open(33,file='./ir2_rad.dat')
    open(34,file='./wv_tb.dat')
    open(35,file='./wv_rad.dat')
    open(36,file='./swir_tb.dat')
    open(37,file='./swir_rad.dat')
    open(38,file='./vis_tb.dat')
    open(39,file='./vis_rad.dat')

    do i=1,nx
	write(30, '(1543(f8.4,x),f8.4)') ir1_tb(i,:)    ! to deal with fd or la, put 'ny-1' value in 1543.
	write(31, '(1543(es11.5,x),es11.5)') ir1_rad(i,:)

	write(32, '(1543(f8.4,x),f8.4)') ir2_tb(i,:)
	write(33, '(1543(es11.5,x),es11.5)') ir2_rad(i,:)

	write(34, '(1543(f8.4,x),f8.4)') wv_tb(i,:)
	write(35, '(1543(es11.5,x),es11.5)') wv_rad(i,:)

	write(36, '(1543(f8.4,x),f8.4)') swir_tb(i,:)
	write(37, '(1543(es11.5,x),es11.5)') swir_rad(i,:)

	write(38, '(1543(f8.4,x),f8.4)') vis_tb(i,:)
	write(39, '(1543(es11.5,x),es11.5)') vis_rad(i,:)
    enddo

    close(30)
    close(31)
    close(32)
    close(33)
    close(34)
    close(35)
    close(36)
    close(37)
    close(38)
    close(39)

else
    ir1_tb_out(:,:) = ir1_tb(:,:) 
    ir2_tb_out(:,:) = ir2_tb(:,:) 
!   wv_tb_out(:,:) = wv_tb(:,:) 
    swir_tb_out(:,:) = swir_tb(:,:) 
    vis_rad_out(:,:) = vis_tb(:,:) / 100.

endif

! read latitude and longitude (optional)
! to deal with fd or la, put proper latitude-longitude table in open command.
! example code uses latitude-longitude table for enh mode.  

write(6,*)' initialize latlon data'
image_lat_ir = r_missing_data
image_lon_ir = r_missing_data

write(6,*)' open latlon data ',trim(path_to_raw_sat)//'/cn_latlon.txt'
open(40,file=trim(path_to_raw_sat)//'/cn_latlon.txt',status='old',err=998)

write(6,*)' read latlon data'
do i=1, nx*ny 
  read(40,*)lattice1, lattice2, rlat, rlon
! write(6,*)lattice1,lattice2,rlat,rlon
  ii = lattice1+1
  jj = lattice2+1
! write(6,*)lattice1,lattice2,ii,jj,rlat,rlon
  image_lat_ir(ii,jj) = rlat
  image_lon_ir(ii,jj) = rlon
enddo 
close(40)

write(6,*)' center ir_tb_out value ',ir1_tb_out(nx/2,ny/2)
write(6,*)' corner ir_tb_out value ',ir1_tb_out(1,1)
write(6,*)' center vis_rad_out value ',vis_rad_out(nx/2,ny/2)
write(6,*)' corner vis_rad_out value ',vis_rad_out(1,1)
write(6,*)' center lat/lon ',image_lat_ir(nx/2,ny/2),image_lon_ir(nx/2,ny/2)
write(6,*)' corner lat/lon ',image_lat_ir(1,1),image_lon_ir(1,1)
write(6,*)' range lat      ',minval(image_lat_ir),maxval(image_lat_ir)
write(6,*)' range lon      ',minval(image_lon_ir),maxval(image_lon_ir)

istatus = 1
return ! normal return

998 write(6,*)' error reading coms lat/lon cn_latlon.txt data'
istatus = 0
return ! error return

end subroutine comscount2tbnrad_sub
