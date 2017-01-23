module boundary_mask_mod
	implicit none

	integer, parameter :: grid_unit=50

	integer(kind=1), dimension(:,:,:), allocatable ::  mask

	character (len=256), dimension(2), parameter :: mask_files = (/ "boundary_mask_1.bin", "boundary_mask_2.bin" /)





contains

! allocate memory
subroutine allocate_mask(number_x_grid,number_y_grid)
	implicit none
	integer, intent(in) :: number_x_grid,number_y_grid
	integer :: istat

	allocate (mask(2,number_x_grid,number_y_grid), stat=istat)
	if (istat /= 0) THEN
		WRITE(6,*) "allocation error for boundary mask", istat
		stop
	ENDIF

end subroutine allocate_mask


! clear memory
subroutine boundary_mask_clear()
	implicit none
	integer :: istat

	deallocate (mask, stat=istat)
	if (istat /= 0) THEN
		WRITE(6,*) "deallocation error for boundary mask", istat
		stop
	ENDIF

end subroutine boundary_mask_clear


subroutine write_mask(number_x_grid,number_y_grid)
	implicit none

	integer, intent(in) :: number_x_grid,number_y_grid
	integer :: istat, rec_num, x_counter, y_counter, counter

	do counter = 1, 2, 1

		open(unit=grid_unit, file=mask_files(counter), access="direct", form="unformatted", recl=1, status="replace",&
		 iostat=istat)
		if (istat /= 0) THEN
			WRITE(6,*) "error opening mask file: ", trim(adjustl(mask_files(counter))), istat
			stop
		ENDIF

		rec_num = 0
		do x_counter = 1, number_x_grid
			do y_counter = 1, number_y_grid
				rec_num = rec_num + 1
				write(grid_unit,rec=rec_num) mask(counter,x_counter,y_counter)

			end do
		end do

		close(unit=grid_unit)


	end do

end subroutine write_mask

! read mask from file
subroutine read_mask(number_x_grid,number_y_grid)
	implicit none

	integer, intent(in) :: number_x_grid,number_y_grid
	integer :: istat, rec_num, counter, x_counter, y_counter

	do counter = 1, 2, 1

		open(unit=grid_unit, file=mask_files(counter), access="direct", form="unformatted", recl=1, status="old",&
		 iostat=istat)
		if (istat /= 0) THEN
			WRITE(6,*) "error opening mask file: ", trim(adjustl(mask_files(counter))), istat
			stop
		ENDIF

		rec_num = 0
		do x_counter = 1, number_x_grid
			do y_counter = 1, number_y_grid
				rec_num = rec_num + 1
				read(grid_unit,rec=rec_num) mask(counter,x_counter,y_counter)

			end do
		end do

		close(unit=grid_unit)


	end do

end subroutine read_mask


end module boundary_mask_mod
