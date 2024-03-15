module mem_grid

real   , allocatable, dimension(:,:) :: lat, lon, topo, ldf
integer, allocatable, dimension(:,:) :: istartiend
integer, allocatable, dimension(:)   :: recvcounts,displs
integer                              :: npes,rank

end module
