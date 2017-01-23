program flowlines2

! read in the boundaries, boundary masks, direction field and interpolate between the boundaries

	use read_polygons
	use read_minmax
	use boundary_mask_mod
	use direction_grid_mod

	implicit none

	integer :: commandline_count, polygon_counter, points_counter, x_grid_index, y_grid_index, x_dummy_index, y_dummy_index
	integer :: check_peak, x_increment, y_increment, index_next, extra_counter, extra_total
	integer, parameter :: gmt_unit = 90

	double precision :: grid_spacing, current_x, current_y, current_direction, next_x, next_y, grid_x, grid_y, temp_x, temp_y
	
	double precision :: crossover_x, crossover_y, angle
	double precision, parameter :: r_increment = 0.25, fining_increment=5


	double precision, parameter :: pi = 3.141592653589793, to_round_down = 0.00001

	character(len=256) :: arg_in
	character (len=1), parameter :: divider_character = ">"
	character(len=256), parameter :: gmt_file = "flowlines.txt"


	logical :: found_point, end_line, round_down_x, round_down_y

	commandline_count = command_argument_count()

	if(commandline_count /= 1) THEN
		write(6,*) "correct usage of this program:"
		write(6,*) "flowlines grid_spacing"
		stop
	endif

	call GET_COMMAND_ARGUMENT(1,arg_in) 
	read(arg_in,*) grid_spacing

	! read polygons
	call read_polygons_init()
	! read grid parameters
	call run_read_minmax(grid_spacing)

	! read in the boundary mask
	call allocate_mask(number_x_grid,number_y_grid)
	call read_mask(number_x_grid,number_y_grid)

	! read direction grid
	call read_direction_grid_init(number_x_grid,number_y_grid)
	call read_direction_grid(number_x_grid,number_y_grid)


	! write GMT file for plotting/debugging

	open(unit=gmt_unit, file=gmt_file, access="sequential", form="formatted", status="replace")


	! at each initial boundary point, go out until it reaches the secondary boundary
	do polygon_counter = 1, number_polygons(1)
		do points_counter = 1, polygon_points(1, polygon_counter)
	
			if(points_counter < polygon_points(1, polygon_counter) ) then
				index_next = points_counter + 1
			else
				index_next = 1
			endif
			! add some points in between if there is a larg space
			angle = atan2((y_coordinates(1,polygon_counter,index_next) - y_coordinates(1,polygon_counter,points_counter)), &
				  (x_coordinates(1,polygon_counter,index_next) - x_coordinates(1,polygon_counter,points_counter)))

			extra_total = int(sqrt((y_coordinates(1,polygon_counter,index_next) - &
					            y_coordinates(1,polygon_counter,points_counter))**2 + &
					  (x_coordinates(1,polygon_counter,index_next) - x_coordinates(1,polygon_counter,points_counter))**2)&
					  /fining_increment) + 1

			do extra_counter = 1, extra_total, 1
				current_x = x_coordinates(1,polygon_counter,points_counter) + &
				     dble(extra_counter-1)*fining_increment * cos(angle)
				current_y = y_coordinates(1,polygon_counter,points_counter) + &
				     dble(extra_counter-1)*fining_increment * sin(angle)


				write(gmt_unit,'(A1)') divider_character

				write(gmt_unit,*) current_x, current_y
			!	write(666,'(A1)') divider_character



				call find_grid_index(current_x, current_y, grid_spacing, x_grid_index, y_grid_index)

				flowline_loop: do
					! find grid points


					call find_grid_location(x_grid_index, y_grid_index, grid_spacing, grid_x, grid_y)




					current_direction = direction_grid(x_grid_index, y_grid_index)

					next_x = current_x + r_increment * cos(current_direction)
					next_y = current_y + r_increment * sin(current_direction)


					call find_grid_index(next_x, next_y, grid_spacing, x_dummy_index, y_dummy_index) ! includes inside grid check

					temp_x = next_x
					temp_y = next_y

					end_line = .false.

					if(mask(2,x_grid_index, y_grid_index) == 1 .or. mask(2,x_dummy_index, y_dummy_index) == 1) THEN ! check to see if the line crosses over the boundary

						call crossover_point(current_x, current_y, next_x, next_y, end_line, grid_spacing)
				
					endif

			!		write(666,*) current_x, current_y, temp_x, temp_y, next_x,next_y, mask(2,x_grid_index, y_grid_index) &
			!			 == 1, mask(2,x_dummy_index, y_dummy_index) == 1

					write(gmt_unit,*) next_x, next_y


					if(end_line) then
						exit flowline_loop
					endif

					x_grid_index = x_dummy_index
					y_grid_index = y_dummy_index

					current_x = next_x
					current_y = next_y
				end do flowline_loop

				if(points_counter == 8) THEn
				!	stop
				endif
			end do
		end do

	end do

	close(unit=gmt_unit)

	! clear memory
	call read_polygons_clear()
	call boundary_mask_clear()
	call read_direction_grid_clear()

contains

subroutine crossover_point(current_x, current_y, next_x, next_y, end_line, grid_spacing)

	! find out if there is an overlap between the 

	use read_polygons

	implicit none

	double precision, intent(in) :: current_x, current_y, grid_spacing
	double precision, intent(inout) :: next_x, next_y
	logical, intent(inout) :: end_line



	integer :: polygon_counter, points_counter, next_index, x_index1, y_index1, x_index2, y_index2

	double precision :: slope1, slope2, intercept1, intercept2, temp_x, temp_y1, temp_y2
	double precision :: seg_x1, seg_x2, seg_y1, seg_y2, min_cell_x, min_cell_y, max_cell_x, max_cell_y

	double precision, parameter :: small = 1e-8

	logical :: vertical_line1, vertical_line2

	if(current_x /= next_x) THEN

		slope1 = (current_y - next_y) / (current_x - next_x) 
		intercept1 = current_y - current_x * slope1
		vertical_line1 = .false.
	else

		vertical_line2 = .true.

	endif

	min_cell_x = min(current_x, next_x)
	min_cell_y = min(current_y, next_y)
	max_cell_x = max(current_x, next_x)
	max_cell_y = max(current_y, next_y)

	polygon_loop: do polygon_counter = 1, number_polygons(2), 1
		points_loop: do points_counter = 1, polygon_points(2, polygon_counter), 1

			next_index = points_counter + 1
			if(next_index > polygon_points(2, polygon_counter)) then
				next_index = 1
			endif

			seg_x1 = x_coordinates(2,polygon_counter,points_counter)
			seg_y1 = y_coordinates(2,polygon_counter,points_counter)

			seg_x2 = x_coordinates(2,polygon_counter,next_index)
			seg_y2 = y_coordinates(2,polygon_counter,next_index)


			if(seg_x2 /= seg_x1) THEN
				slope2 = (seg_y2 - seg_y1) / (seg_x2 - seg_x1)
				intercept2 = seg_y2 - slope2 * seg_x2
				vertical_line2 = .false.
			else
				vertical_line2 = .true.
			endif


			if( vertical_line1 .and. .not. vertical_line2 ) THEN

				temp_x = current_x
				temp_y = slope2 * temp_x + intercept2

			elseif(vertical_line2 .and. .not. vertical_line1) THEN

				temp_x = seg_x1
				temp_y = slope1 * temp_x + intercept1

			elseif (vertical_line1 .and. vertical_line2) THEN

				temp_x = current_x


				if(min(seg_y1,seg_y2) < min_cell_y) then
					temp_y = (max(seg_y1,seg_y2) + max_cell_y) / 2.
				else
					temp_y = (min(seg_y1,seg_y2) + min_cell_y) / 2.
				endif

			else
				temp_x = (intercept2 - intercept1) / (slope1 - slope2)
				

				temp_y = slope1 * temp_x + intercept1
			endif

			if( temp_x >= min_cell_x .and. temp_x <=max_cell_x .and. temp_y >= min_cell_y .and. temp_y <= max_cell_y) THEN ! crossover
				next_x = temp_x
				next_y = temp_y
				end_line = .true.
				write(546,*) next_x, next_y
				exit polygon_loop
			endif


		end do points_loop
	end do polygon_loop


end subroutine crossover_point

end program flowlines2
