
        subroutine put_vrc(i4time,comment_2d 
     1                         ,rlat_radar,rlon_radar,rheight_radar
     1                         ,lat,lon,topo
     1                         ,dbz_2d,imax,jmax,i_vrc        
     1                         ,vrc_outdir,r_missing_data,istatus)

!       stuff from 'put_laps_2d' except how do we handle radar subdir?

        character*150 vrc_full_path
        character*31 ext

        character*125 comment_2d
        character*125 comments_2d(2)
        character*10 units(2)
        character*3 vars(2)

        integer lvl_2d(2)

        real fields_2d(imax,jmax,2)
        real lat(imax,jmax)
        real lon(imax,jmax)
        real topo(imax,jmax)
        real dist(imax,jmax)
        real dbz_2d(imax,jmax)

        character*9 a9time
        character*8 radar_subdir
        character*3 c3_radar_subdir
        character*7 vrc_outdir
        character*4 lvl_coord_2d(2)

        call make_fnam_lp(i4time,a9time,istatus)
        if(istatus .ne. 1)return

        write(6,*)' subroutine put_vrc for ',a9time,' ',vrc_outdir

        call get_ref_base(ref_base,istatus)
        if(istatus .ne. 1)return

!       calculate closest radar array
        write(6,*)' calculating closest radar array (dist to vrc radar)'       
        do i = 1,imax
        do j = 1,jmax
            call latlon_to_radar(lat(i,j),lon(i,j),topo(i,j)
     1                          ,azimuth,dist(i,j),elev
     1                          ,rlat_radar,rlon_radar,rheight_radar)       
        enddo ! j
        enddo ! i

!       call vrc_clutter_thresh(      fields_2d(1,1,1)                   ! i/o
!    1                               ,dist                               ! i
!    1                               ,imax,jmax,ref_base,r_missing_data) ! i

!       set up 2 fields for output
        do i = 1,imax
        do j = 1,jmax
            fields_2d(i,j,1) = dbz_2d(i,j)
            if(fields_2d(i,j,1) .ne. r_missing_data)then
                fields_2d(i,j,2) = dist(i,j)
            else
                fields_2d(i,j,2) = r_missing_data
            endif
        enddo ! j
        enddo ! i

        ext = 'vrc'

        vars(1) = 'ref'
        vars(2) = 'dis'

        units(1) = 'dbz'
        units(2) = 'm'

        lvl_2d(1) = 0
        lvl_2d(2) = 0
        
        lvl_coord_2d(1) = 'msl'
        lvl_coord_2d(2) = 'msl'

        comments_2d(1) = comment_2d
        comments_2d(2) = comment_2d

        call get_vrc_full_path(i_vrc,vrc_outdir
     1                        ,vrc_full_path,lenv,istatus)
        

        write(6,11)vrc_full_path(1:lenv),ext(1:5),vars
11      format(' writing 2d ',a,1x,a5,2(1x,a3))

        call write_laps_data(i4time,vrc_full_path,ext,imax,jmax,
     1                       2,2,vars,lvl_2d,lvl_coord_2d,units,
     1                       comments_2d,fields_2d,istatus)
        if(istatus .eq. 1)then
            write(6,*)' vrc successfully written'
        else
            write(6,*)' vrc not successfully written', istatus
        endif

        return
        end


        subroutine get_vrc_full_path(i_vrc,vrc_outdir
     1                              ,vrc_full_path,lenv,istatus)

        character*(*) vrc_full_path
        character*(*) vrc_outdir
        character*150 directory1
        character*3 c3_radar_subdir

        write(c3_radar_subdir,801)i_vrc
 801    format(i3.3)

        call s_len(vrc_outdir,lenp)
        if(lenp .gt. 0)then
            write(6,*)'vrc_outdir = ',vrc_outdir(1:lenp)
        else
            write(6,*)'get_vrc_full_path: vrc_outdir has zero length'
            istatus = 0
            return
        endif

        if(vrc_outdir .eq. 'rdr')then
            write(6,*)' radar_subdir = ',c3_radar_subdir

            call get_directory('rdr',directory1,len_dir1)

            vrc_full_path = directory1(1:len_dir1)//c3_radar_subdir
     1                                            //'/vrc/'  
            call s_len(vrc_full_path,lenv)

        else ! 'lapsprd'
            if(i_vrc .gt. 1)then
                write(6,*)' warning - lapsprd output when i_vrc > 1'
                write(6,*)' i_vrc = ',i_vrc
                istatus = 0
                return
            endif
            call get_directory('vrc',vrc_full_path,lenv)

        endif            

        istatus = 1
        return
        end
