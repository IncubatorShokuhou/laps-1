        function delta_t(t) ! tdt - ut1  (days)

        implicit real*8 (a-z)

        integer*4 n_entries
        parameter (n_entries = 34)

        real ar(2,n_entries)
        save ar

        integer*4 init
        save init
        data init/0/

        integer*4 iby,i,len_dir
        character*150 static_dir,filename

        include '../../include/astparms.for'

        if(init .eq. 0)then
            call get_directory('static',static_dir,len_dir)
            filename = static_dir(1:len_dir)//'/delta_t.parms'          
            open(1,file=trim(filename),status='old',form='formatted')
            read(1,*)ar
            init = 1
        endif

        by = 1950.d0 + (t - t1950)/ 365.2421988d0

        if(by .lt. ar(1,1) .or. by .gt. ar(1,n_entries))then
            tcent = (by-1900.d0)/100.d0
            delta_t_min =
     1       0.41d0 + 1.2053d0 * tcent + 0.4992 * tcent**2
            delta_t = delta_t_min / 1440d0
c           write(6,*)' general formula'
c           write(6,*)by,tcent,delta_t_min,delta_t*86400.
            return

        else
            do i = n_entries-1,1,-1
                if(by .gt. ar(1,i))then
                    frac = (by - ar(1,i)) / (ar(1,i+1) - ar(1,i))
                    delta_t_sec = ar(2,i)   *  (1.-frac)
     1                  + ar(2,i+1) *      frac
                    delta_t = delta_t_sec / 86400d0
c                   write(6,*)' interpolation'
c                   write(6,*)by,i,frac,delta_t*86400.
                    return
                endif
            enddo
        endif

        write(6,*)' error in deltat ',t
        stop

        end
