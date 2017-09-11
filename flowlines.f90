program flowlines2

! read in the boundaries, boundary masks, direction field and interpolate between the boundaries

	use read_polygons
	use read_minmax
	use boundary_mask_mod
	use direction_grid_mod
	use crossover

	implicit none

	integer :: commandline_count, polygon_counter, points_counter, x_grid_index, y_grid_index, x_dummy_index, y_dummy_index
	integer :: check_peak, x_increment, y_increment, index_next, extra_counter, extra_total, flowline_point_count, flow_counter
	integer :: max_line, max_point
	integer, parameter :: gmt_unit = 90, max_flowline_points = 100000, discard_unit=200

	double precision :: grid_spacing, current_x, current_y, current_direction, next_x, next_y, grid_x, grid_y, temp_x, temp_y
	
	double precision :: crossover_x, crossover_y, angle, closeness_threshold, r_increment
	double precision, parameter :: r_increment_minimum = 0.02, threshold_factor = 100.

	double precision, dimension(max_flowline_points) :: x_flowline_store, y_flowline_store


	double precision, parameter :: pi = 3.141592653589793, to_round_down = 0.00001

	integer :: last_x_index, last_y_index
	integer, allocatable, dimension(:,:) :: visited_grid

	character(len=256) :: arg_in
	character (len=1), parameter :: divider_character = ">"
	character(len=256), parameter :: gmt_file = "flowlines.txt", discard_file="discarded_flowlines.txt"


	logical :: found_point, end_line, round_down_x, round_down_y, looping, hit_saddle, oscillating, return_status, outside




	commandline_count = command_argument_count()

	if(commandline_count /= 1) THEN
		write(6,*) "correct usage of this program:"
		write(6,*) "flowlines grid_spacing"
		stop
	endif

	call GET_COMMAND_ARGUMENT(1,arg_in) 
	read(arg_in,*) grid_spacing

	r_increment = grid_spacing * r_increment_minimum

	closeness_threshold = r_increment / threshold_factor

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

	allocate(visited_grid(number_x_grid,number_y_grid))

	! write GMT file for plotting/debugging

	open(unit=gmt_unit, file=gmt_file, access="sequential", form="formatted", status="replace")

	max_line = 0 
	max_point = 0

	open(unit=discard_unit, file=discard_file, access="sequential", form="formatted", status="replace")

	! at each initial boundary point, go out until it reaches the secondary boundary
	do polygon_counter = 1, number_polygons(1)
		do points_counter = 1, polygon_points(1, polygon_counter)
			write(6,*) "polygon: ", polygon_counter, "/", number_polygons(1), " point:", points_counter, "/", &
				polygon_points(1, polygon_counter)
			if(points_counter < polygon_points(1, polygon_counter) ) then
				index_next = points_counter + 1
			else
				index_next = 1
			endif




			flowline_point_count = 1

			x_flowline_store(flowline_point_count) = x_coordinates(1,polygon_counter,points_counter) 
			y_flowline_store(flowline_point_count) = y_coordinates(1,polygon_counter,points_counter)

			call find_grid_index(x_flowline_store(flowline_point_count), y_flowline_store(flowline_point_count),&
			   grid_spacing, x_grid_index, y_grid_index, return_status)

			visited_grid = 0
			last_x_index = x_grid_index
			last_y_index = y_grid_index

			oscillating = .false.
			outside = .false.
			flowline_loop: do
				! find grid points
				flowline_point_count = flowline_point_count + 1



				call find_grid_location(x_grid_index, y_grid_index, grid_spacing, grid_x, grid_y)


				current_direction = direction_grid(x_grid_index, y_grid_index)

				x_flowline_store(flowline_point_count) = x_flowline_store(flowline_point_count-1) +&
				  r_increment * cos(current_direction)
				y_flowline_store(flowline_point_count) = y_flowline_store(flowline_point_count-1) +&
				  r_increment * sin(current_direction)

				if(x_grid_index /= last_x_index .or. y_grid_index /= last_y_index) THEN
					visited_grid(last_x_index,last_y_index) = visited_grid(last_x_index,last_y_index) + 1
				endif
				
				last_x_index = x_grid_index
				last_y_index = y_grid_index


				if(visited_grid(x_grid_index,y_grid_index) > 2) THEN
					write(6,*) "likely oscillating"
					oscillating = .true.
					exit flowline_loop
				endif
				

				if(isnan(x_flowline_store(flowline_point_count))) THEN
					write(6,*) current_direction
					stop
				end if

				call find_grid_index(x_flowline_store(flowline_point_count), y_flowline_store(flowline_point_count),&
				  grid_spacing, x_dummy_index, y_dummy_index, return_status) ! includes inside grid check


				if(.not. return_status) THEN

					outside =.true.
					exit flowline_loop

				endif

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
!					write(549,*) "> ", points_counter, x_flowline_store(flowline_point_count), &
!					  y_flowline_store(flowline_point_count)
					exit flowline_loop
				endif

				! check to see if the line crosses over itself



				if(flowline_point_count > 4) THEN

					call cross_itself(x_flowline_store(1:flowline_point_count), &
					  y_flowline_store(1:flowline_point_count), flowline_point_count, looping)

				endif

				if(looping) then

!					write(549,*) ">>", points_counter, x_flowline_store(flowline_point_count), &
!					  y_flowline_store(flowline_point_count)
					write(6,*) "detected looping"
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
					write(6,*) "detected a saddle point"
!					write(549,*) ">>>", points_counter, x_flowline_store(flowline_point_count), &
!					  y_flowline_store(flowline_point_count), x_flowline_store(flowline_point_count-2), &
!					  y_flowline_store(flowline_point_count-2)
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

!			if( .not. hit_saddle) THEN
				if(flowline_point_count > max_point) THEN
					max_point = flowline_point_count
					max_line = points_counter
				endif
				if (oscillating .or. outside .or. looping .or. hit_saddle) THEN
				!	write(gmt_unit,'(A1,I7,I7)') divider_character, points_counter, 0
					write(discard_unit,'(A1,I7,I7)') divider_character, points_counter, flowline_point_count
				else
					write(gmt_unit,'(A1,I7,I7)') divider_character, points_counter, flowline_point_count
				endif

				do flow_counter = 1, flowline_point_count, 1
					if(oscillating .or. outside .or. looping .or. hit_saddle) THEN
						write(discard_unit,*) x_flowline_store(flow_counter), y_flowline_store(flow_counter)
					else
						write(gmt_unit,*) x_flowline_store(flow_counter), y_flowline_store(flow_counter)
					endif

				end do
!			endif


!			end do
		end do

	end do

	write(6,*) max_line, max_point

	close(unit=gmt_unit)

	deallocate(visited_grid)

	! clear memory
	call read_polygons_clear()
	call boundary_mask_clear()
	call read_direction_grid_clear()


end program flowlines2
