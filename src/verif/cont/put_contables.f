
        subroutine put_contables(i4_initial,i4_valid,nthr
     1                          ,cont_4d,ni,nj,nk,cont_dir)  

        character*150 cont_dir,filename
        character*31 ext

        character*125 comment_3d(nk*nthr),comment_2d
        character*10 units_3d(nk*nthr),units_2d
        character*3 var_3d(nk*nthr),var_2d
        integer lvl_3d(nk*nthr)
        character*4 lvl_coord_3d(nk*nthr)
        logical ltest_vertical_grid

        real o(ni,nj,nk)
        real f(ni,nj,nk)
        real cont_3d(ni,nj,nk)

        integer contable(0:1,0:1)

!       write contingency tables as an ascii free format file
!       write(6,*)' write out 3-d contingency table (ascii) in '
!    1             ,cont_dir

!       call s_len(cont_dir,lend)
!       filename = cont_dir(1:lend)//'cont3d'
!       open(12,file=filename,status = 'unknown')
!       write(12,*)cont_4d
!       close(12)

!       write netcdf contingency table with a call to 'write_laps'
        if(.true.)then
            ext = 'cont'

            units_2d = ' '
            comment_2d = ' laps radar verification contingency table'

            do ithr = 1,nthr
              idbz = 10 + ithr*10
              write(var_2d,1)idbz
 1            format('t',i2.2)

              write(6,*)' 3d contingency table var = ',var_2d

              do k = 1,nk

                kk = k + (ithr-1)*nk ! calcualate 1d subscript
                units_3d(kk)   = units_2d
                comment_3d(kk) = comment_2d
                if(ltest_vertical_grid('height'))then  
                    lvl_3d(kk) = zcoord_of_level(k)/10
                    lvl_coord_3d(kk) = 'msl'
                elseif(ltest_vertical_grid('pressure'))then 
                    lvl_3d(kk) = nint(zcoord_of_level(k))/100
                    lvl_coord_3d(kk) = 'hpa'
                else
                    write(6,*)' error, vertical grid not supported,'
     1                     ,' this routine supports pressure or height'     
                    istatus = 0
                    return  
                endif

                var_3d(kk) = var_2d

              enddo ! k
            enddo ! ithr

            nf = nk*nthr

            write(6,*)' calling write_laps where ext is ',ext

            call write_laps(i4_valid,i4_valid,cont_dir,ext,
!           call write_laps(i4_initial,i4_valid,cont_dir,ext,
     .                      ni,nj,nf,nf,var_3d,
     .                      lvl_3d,lvl_coord_3d,units_3d,comment_3d,
     .                      cont_4d,istatus)
            if(istatus .ne. 1)then
                print *
     1           ,'error writing interpolated data to laps database.'      
                return
            endif
        endif

        return
        end
