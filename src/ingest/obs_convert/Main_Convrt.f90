program obs_convert

!==========================================================
!doc  this program reads observation data through laps ingest
!doc  and converts the data into requested format, e.g., bufr
!doc
!doc  history:
!doc	creation:	yuanfu xie	may 2007
!==========================================================

  use laps_params

  implicit none

  ! local variables:
  character :: suffix*14,wtfile*201,rdfile*201
  integer   :: length(2),ibfmsg(10000)
  
  ! laps configuration:
  call laps_config

  ! open file to write:
  if (format_request .eq. 'bufr') then
    suffix = 'bufr'
    call get_directory(suffix,wtfile,length(1))
    ! open bufr table:
    call get_directory('static',rdfile,length(2))
    rdfile(length(2)+1:length(2)+22) = 'prepobs_prep.bufrtable'
    open(unit=bfrtbl_channel,file=rdfile(1:length(2)+22),status='old')

    write(wtfile(length(1)+1:length(1)+14),1) system_asctime,'.bufr'
1   format(a9,a5)
    length = length(1)+14

    open(unit=output_channel,file=wtfile(1:length(1)),form='unformatted')
    call openbf(output_channel,'out',bfrtbl_channel)
  else if (format_request .eq. 'wrf') then
    suffix = 'wrf'
    call get_directory(suffix,wtfile,length(1))
    write(wtfile(length(1)+1:length(1)+13),2) system_asctime,'.wrf'
2   format(a9,a4)
    length = length(1)+13

    open(unit=output_channel,file=wtfile(1:length(1)),form='formatted')
  else
    write(6,*) 'main: unknown data format'
    stop
  endif

  ! laps observation data ingest:
  call laps_ingest

  ! close file:
  if (format_request .eq. 'bufr') then
    close(bfrtbl_channel)
    ! force flush the last bufr message:
    call writsa(-output_channel,ibfmsg,length(1))
  endif
  close(output_channel)

end program obs_convert
