
      subroutine write_orb_att(path,sat_name_in,num_atts,orb_att)
c
c**********************************************************************
c
        implicit  none
c
        character*(*)  path
        character*(*)  sat_name_in
        character*15   sat_name_down

        integer          num_atts
        real*8         orb_att(num_atts)
c
        integer        sat_name_len,
     1                 end_path,
     1                 noerror, error,
     1                 istatus,
     1                 cdl_len,
     1                 filename_len
c
        character*150  filename
        character*150  cdlfile
c
c ****  begin
c
        noerror=1   !success
        error=0     !error

c ****  downcase sat_name_in for making file name
        call downcase(sat_name_in, sat_name_down)
        call s_len(sat_name_down,sat_name_len)

c ****  make filename: path/satname_orbatt.dat

        call s_len(path,end_path)
        if (path(end_path:end_path) .eq. '/') then
        else
          path = path//'/'
          end_path = end_path + 1
        endif
        filename = 
     1path(1:end_path)//sat_name_down(1:sat_name_len)//'_orbatt.dat'
        call s_len(filename,filename_len) 
c
c **** make full cdl file name        
c
        call get_directory('cdl',cdlfile, cdl_len)
        if (cdlfile(cdl_len:cdl_len) .eq. '/') then
        else
          cdlfile = cdlfile(1:cdl_len)//'/'
          cdl_len = cdl_len + 1
        endif
        cdlfile = cdlfile(1:cdl_len)//'orb.cdl'
        cdl_len = cdl_len + 7
c
c **** call c program write_att_c to write file
c
        print*,'calling write_att_c'
        call write_att_c(filename, filename_len, cdlfile, cdl_len, 
     1                   sat_name_in, sat_name_len, num_atts, orb_att, 
     1                   istatus)

        if (istatus .eq. noerror)  goto 999  !success
        if (istatus .eq. -1) goto 950 
        if (istatus .eq. -2) goto 960
        if (istatus .eq. -3) goto 970
        if (istatus .eq. 2) goto 980
c
c ****  return normally.
c
        istatus=noerror
999     return
c
c ****  error trapping.
c
950     write (6,*) 'error opening netcdf file...write aborted.'
        istatus=error
        goto 999
c
960     write (6,*) 'error creating netcdf file...write aborted.'
        istatus=error
        goto 999
c
970     write (6,*) 'error writing netcdf file...write aborted.'
        istatus=error
        goto 999
c
980     write (6,*) 'error in array dimensioning...write aborted.'
        istatus=error
        goto 999
c
992     continue
        istatus=error
        goto 999
c
        end

