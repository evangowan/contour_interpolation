program flowlines2

! read in the boundaries, boundary masks, direction field and interpolate between the boundaries

	use read_polygons
	use read_minmax
	use boundary_mask_mod
	use direction_grid_mod

	implicit none

	integer :: commandline_count, polygon_counter, points_counter, x_grid_index, y_grid_index, x_dummy_index, y_dummy_index
	integer :: check_peak, x_increment, y_increment, index_next, extra_counter, extra_total, flowline_point_count, flow_counter
	integer :: max_line, max_point
	integer, parameter :: gmt_unit = 90, max_flowline_points = 100000

	double precision :: grid_spacing, current_x, current_y, current_direction, next_x, next_y, grid_x, grid_y, temp_x, temp_y
	
	double precision :: crossover_x, crossover_y, angle, closeness_threshold
	double precision, parameter :: r_increment = 0.25, threshold_factor = 100.

	double precision, dimension(max_flowline_points) :: x_flowline_store, y_flowline_store


	double precision, parameter :: pi = 3.141592653589793, to_round_down = 0.00001

	character(len=256) :: arg_in
	character (len=1), parameter :: divider_character = ">"
	character(len=256), parameter :: gmt_file = "flowlines.txt"


	logical :: found_point, end_line, round_down_x, round_down_y, looping, hit_saddle

	closeness_threshold = r_increment / threshold_factor


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

	max_line = 0
	max_point = 0
	! at each initial boundary point, go out until it reaches the secondary boundary
	do polygon_counter = 1, number_polygons(1)
		do points_counter = 1, polygon_points(1, polygon_counter)
	
			if(points_counter < polygon_points(1, polygon_counter) ) then
				index_next = points_counter + 1
			else
				index_next = 1
			endif




			flowline_point_count = 1

			x_flowline_store(flowline_point_count) = x_coordinates(1,polygon_counter,points_counter) 
			y_flowline_store(flowline_point_count) = y_coordinates(1,polygon_counter,points_counter)

			call find_grid_index(x_flowline_store(flowline_point_count), y_flowline_store(flowline_point_count),&
			   grid_spacing, x_grid_index, y_grid_index)



			flowline_loop: do
				! find grid points
				flowline_point_count = flowline_point_count + 1



				call find_grid_location(x_grid_index, y_grid_index, grid_spacing, grid_x, grid_y)

				current_direction = direction_grid(x_grid_index, y_grid_index)

				x_flowline_store(flowline_point_count) = x_flowline_store(flowline_point_count-1) +&
				  r_increment * cos(current_direction)
				y_flowline_store(flowline_point_count) = y_flowline_store(flowline_point_count-1) +&
				  r_increment * sin(current_direction)


				if(isnan(x_flowline_store(flowline_point_count))) THEN
					write(6,*) current_direction
					stop
				end if

				call find_grid_index(x_flowline_store(flowline_point_count), y_flowline_store(flowline_point_count),&
				  grid_spacing, x_dummy_index, y_dummy_index) ! includes inside grid check


				! check to see if the line crosses over the boundary

				end_line = .false.
				looping = .false.
				hit_saddle = .false.

				if(mask(2,x_grid_index, y_grid_index) == 1 .or. mask(2,x_dummy_index, y_dummy_index) == 1) THEN
					call cross_polygon(x_flowline_store(flowline_point_count-1), &
					  y_flowline_store(flowline_point_count-1), x_flowline_store(flowline_point_count), &
					  y_flowline_store(flowline_point_count), end_line, grid_spacing)

				endif

				if(end_line) then
					write(549,*) "> ", points_counter, x_flowline_store(flowline_point_count), &
					  y_flowline_store(flowline_point_count)
					exit flowline_loop
				endif

				! check to see if the line crosses over itself



				if(flowline_point_count > 4) THEN

					call cross_itself(x_flowline_store(1:flowline_point_count), &
					  y_flowline_store(1:flowline_point_count), flowline_point_count, looping)

				endif

				if(looping) then

					write(549,*) ">>", points_counter, x_flowline_store(flowline_point_count), &
					  y_flowline_store(flowline_point_count)
					exit flowline_loop
				endif


				if(flowline_point_count > 3) THEN ! need three points to check this

					if(abs(x_flowline_store(flowline_point_count) - x_flowline_store(flowline_point_count-2)) &
						< closeness_threshold .and. &
					   abs(y_flowline_store(flowline_point_count) - y_flowline_store(flowline_point_count-2)) &
						< closeness_threshold) THEN
						hit_saddle = .true.
					endif
				endif

				if(hit_saddle) then

					write(549,*) ">>>", points_counter, x_flowline_store(flowline_point_count), &
					  y_flowline_store(flowline_point_count), x_flowline_store(flowline_point_count-2), &
					  y_flowline_store(flowline_point_count-2)
					exit flowline_loop
				endif

				if(flowline_point_count == max_flowline_points) THEN
					write(6,*) "Warning: exceeded number of points, maybe went out of bounds?"

					do flow_counter = 1, flowline_point_count, 1

						write(887,*) x_flowline_store(flow_counter), y_flowline_store(flow_counter)

					end do

				endif

				x_grid_index = x_dummy_index
				y_grid_index = y_dummy_index

			end do flowline_loop

			if(.not. looping .and. .not. hit_saddle) THEN
				if(flowline_point_count > max_point) THEN
					max_point = flowline_point_count
					max_line = points_counter
				endif

				write(gmt_unit,'(A1,I7,I7)') divider_character, points_counter, flowline_point_count

				do flow_counter = 1, flowline_point_count, 1

					write(gmt_unit,*) x_flowline_store(flow_counter), y_flowline_store(flow_counter)

				end do
			endif


!			end do
		end do

	end do

	write(6,*) max_line, max_point

	close(unit=gmt_unit)

	! clear memory
	call read_polygons_clear()
	call boundary_mask_clear()
	call read_direction_grid_clear()

contains

subroutine cross_polygon(current_x, current_y, next_x, next_y, end_line, grid_spacing)

	! find out if there is an overlap between the current line segment and the inner polygons

	use read_polygons

	implicit none

	double precision, intent(in) :: current_x, current_y, grid_spacing
	double precision, intent(inout) :: next_x, next_y
	logical, intent(inout) :: end_line



	integer :: polygon_counter, points_counter, next_index


	double precision :: seg_x1, seg_x2, seg_y1, seg_y2, temp_x, temp_y

	logical :: is_crossover

	end_line = .false.


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


			call crossover_point(current_x, current_y, next_x, next_y, seg_x1, seg_y1, seg_x2, seg_y2, &
			  temp_x, temp_y, is_crossover)

			if( is_crossover ) THEN ! crossover
				next_x = temp_x
				next_y = temp_y
				end_line = .true.
				write(546,*) next_x, next_y
				exit polygon_loop
			endif


		end do points_loop
	end do polygon_loop


end subroutine cross_polygon



subroutine cross_itself(x_flowline_store, y_flowline_store, flowline_point_count, end_line)
	! checks to see if the line crosses over itself, will adjust flowline_point_count if it does

	implicit none
	integer, intent(inout) :: flowline_point_count
	double precision, dimension(flowline_point_count), intent(inout) :: x_flowline_store, y_flowline_store
	logical, intent(out) :: end_line

	integer :: point_counter, end_index

	double precision :: a_x1, a_y1, a_x2, a_y2, b_x1, b_y1, b_x2, b_y2, crossover_x, crossover_y
	logical :: is_crossover

	! for ease of reading

	a_x1 = x_flowline_store(flowline_point_count-1)
	a_y1 = y_flowline_store(flowline_point_count-1)
	a_x2 = x_flowline_store(flowline_point_count)
	a_y2 = y_flowline_store(flowline_point_count)
	is_crossover = .false.

	point_loop: do point_counter = 1, flowline_point_count-3, 1 ! shouldn't need to check the last line, unless there is a freak case where the line goes completely the opposite way

		b_x1 = x_flowline_store(point_counter)
		b_y1 = y_flowline_store(point_counter)
		b_x2 = x_flowline_store(point_counter+1)
		b_y2 = y_flowline_store(point_counter+1)
		
		call crossover_point(a_x1, a_y1, a_x2, a_y2, b_x1, b_y1, b_x2, b_y2, crossover_x, crossover_y, is_crossover)

		if(is_crossover) THEN
			write(548,*) ">", point_counter, flowline_point_count
			write(548,*) a_x1, a_y1,b_x1, b_y1,  crossover_x, crossover_y
			write(548,*) a_x2, a_y2,b_x2, b_y2
			end_index = point_counter
			exit point_loop

		end if

	end do point_loop

	if(is_crossover) THEN

	!	write(6,*) end_index, flowline_point_count
		flowline_point_count = end_index
		x_flowline_store(end_index) = crossover_x
		y_flowline_store(end_index) = crossover_y
		end_line = .true.
		write(547,*) crossover_x, crossover_y
	else
		end_line = .false.
	end if

end subroutine cross_itself


subroutine crossover_point(a_x1, a_y1, a_x2, a_y2, b_x1, b_y1, b_x2, b_y2, crossover_x, crossover_y, is_crossover)

	! checks if the two given line segments overlap. If they do, this subroutine returns the crossover point

	implicit none

	double precision, intent(in) :: a_x1, a_y1, a_x2, a_y2, b_x1, b_y1, b_x2, b_y2
	double precision, intent(out) :: crossover_x, crossover_y
	logical, intent(out) :: is_crossover

	double precision :: slope1, slope2, intercept1, intercept2, a_min_cell_x, a_min_cell_y, a_max_cell_x, a_max_cell_y
	double precision :: b_min_cell_x, b_min_cell_y, b_max_cell_x, b_max_cell_y
	double precision :: temp_x, temp_y1, temp_y2
	double precision, parameter :: threshold = 1e6

	logical :: vertical_line1, vertical_line2

	if(a_x1 /= a_x2) THEN
		
		
		vertical_line1 = .false.
	else
		vertical_line1 = .true.
	endif


	a_min_cell_x = min(a_x1, a_x2)
	a_min_cell_y = min(a_y1, a_y2)
	a_max_cell_x = max(a_x1, a_x2)
	a_max_cell_y = max(a_y1, a_y2)

	b_min_cell_x = min(b_x1, b_x2)
	b_min_cell_y = min(b_y1, b_y2)
	b_max_cell_x = max(b_x1, b_x2)
	b_max_cell_y = max(b_y1, b_y2)


	is_crossover = .false.



	slope1 = (a_y2 - a_y1) / (a_x2 - a_x1) 
	intercept1 = a_y2 - a_x2 * slope1

	slope2 = (b_y2 - b_y1) / (b_x2 - b_x1) 
	intercept2 = b_y2 - b_x2 * slope2

	temp_x = (intercept2 - intercept1) / (slope1 - slope2)
	temp_y = slope1 * temp_x + intercept1

	if(a_x2 == a_x1) then
		temp_x = a_x2
		temp_y = slope2 * temp_x + intercept2
	elseif(b_x2 == b_x1) THEN
		temp_x = b_x2
		temp_y = slope1 * temp_x + intercept1
	endif


	is_crossover = .false.
	if(temp_x > a_min_cell_x) THEN
	  if(temp_x < a_max_cell_x) THEN
	     if(temp_x > b_min_cell_x) THEN
		 if(temp_x < b_max_cell_x) THEN
	   	   if(temp_y > a_min_cell_y) THEN
		     if(temp_y < a_max_cell_y) THEN
	   		 if(temp_y > b_min_cell_y) THEN
			   if(temp_y < b_max_cell_y) THEN

				is_crossover = .true.
			   endif
			 endif
		     endif
		   endif
		 endif
	     endif
	   endif
	endif





	crossover_x = temp_x
	crossover_y = temp_y
end subroutine crossover_point

end program flowlines2
