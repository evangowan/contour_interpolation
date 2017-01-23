module direction_grid_mod
	implicit none
	integer, parameter :: dir_grid_unit=70

	double precision, dimension(:,:), allocatable :: direction_grid

	character (len=256), parameter :: dir_grid_file = "direction_grid.bin"

	logical, dimension(:,:), allocatable :: direction_mask

contains

! initialize memory
subroutine read_direction_grid_init(number_x_grid,number_y_grid)
	implicit none

	integer, intent(in) :: number_x_grid,number_y_grid

	integer :: istat

	allocate (direction_mask(number_x_grid,number_y_grid), direction_grid(number_x_grid,number_y_grid), stat=istat)
	if (istat /= 0) THEN
		WRITE(6,*) "allocation error in read_direction_grid_init"
		stop
	ENDIF

end subroutine read_direction_grid_init


! clear memory
subroutine read_direction_grid_clear()
	implicit none

	integer :: istat
	deallocate( direction_mask, direction_grid, stat=istat)
	if (istat /= 0) THEN
		WRITE(6,*) "deallocation error in read_direction_grid_clear", istat
		stop
	ENDIF
end subroutine read_direction_grid_clear


! write out file
subroutine write_direction_grid(number_x_grid,number_y_grid)
	implicit none

	integer, intent(in) :: number_x_grid,number_y_grid
	integer :: istat, rec_num, x_counter, y_counter


	open(unit=dir_grid_unit, file=dir_grid_file, access="direct", form="unformatted", recl=8, status="replace", iostat=istat)
	if (istat /= 0) THEN
		WRITE(6,*) "error opening direction grid file for writing: ", trim(adjustl(dir_grid_file)), istat
		stop
	ENDIF

	rec_num = 0
	do x_counter = 1, number_x_grid
		do y_counter = 1, number_y_grid
			rec_num = rec_num + 1
			if( .not. direction_mask(x_counter, y_counter)) THEN

				write(dir_grid_unit, rec=rec_num) 0.
			else
				write(dir_grid_unit,rec=rec_num) direction_grid(x_counter,y_counter)
			endif
		end do
	end do

	close(unit=dir_grid_unit)


end subroutine write_direction_grid


! read in file

subroutine read_direction_grid(number_x_grid,number_y_grid)
	implicit none

	integer, intent(in) :: number_x_grid,number_y_grid
	integer :: istat, rec_num, x_counter, y_counter


	open(unit=dir_grid_unit, file=dir_grid_file, access="direct", form="unformatted", recl=8, status="old", iostat=istat)
	if (istat /= 0) THEN
		WRITE(6,*) "error opening direction grid file for reading: ", trim(adjustl(dir_grid_file)), istat
		stop
	ENDIF

	rec_num = 0
	do x_counter = 1, number_x_grid
		do y_counter = 1, number_y_grid
			rec_num = rec_num + 1
			read(dir_grid_unit,rec=rec_num) direction_grid(x_counter,y_counter)

		end do
	end do

	close(unit=dir_grid_unit)


end subroutine read_direction_grid

end module direction_grid_mod
