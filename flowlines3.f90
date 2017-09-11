program flowlines3

! read in the boundaries, boundary masks, direction field and interpolate between the boundaries

	use read_polygons
	use read_minmax
	use boundary_mask_mod
	use direction_grid_mod


	implicit none

	integer :: commandline_count, polygon_counter, points_counter
	integer ::  index_next, extra_counter, extra_total, flowline_point_count
	integer, parameter :: gmt_unit = 90, max_flowline_points = 100000, discard_unit=200



	double precision :: grid_spacing

	
	double precision ::  r_increment
	double precision, parameter :: r_increment_minimum = 0.1

	double precision, dimension(max_flowline_points) :: x_flowline_store, y_flowline_store, distance_store

	double precision, parameter :: pi = 3.141592653589793, to_round_down = 0.00001

	character(len=256) :: arg_in

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

			call flowline_loop(x_flowline_store,y_flowline_store,distance_store,grid_spacing, r_increment,&
	                   oscillating,hit_saddle,outside,flowline_point_count)

			call write_flowline(x_flowline_store,y_flowline_store,distance_store,flowline_point_count,gmt_unit,&
				  discard_unit,hit_saddle, oscillating, outside)
			


		end do

	end do



	close(unit=gmt_unit)



	! clear memory
	call read_polygons_clear()
	call boundary_mask_clear()
	call read_direction_grid_clear()


contains

double precision function bicubic(x,y,corner_x,corner_y,corner_values)

	! returns the bilinearly interpolated value of a given x,y value within the given square

	implicit none

	double precision, intent(in) :: x, y
	double precision, dimension(2), intent(in) :: corner_x,corner_y
	double precision, dimension(2,2), intent(in) :: corner_values

      

	!  bilinear interpolation alorithm
	bicubic = (corner_values(1,1)*(corner_x(2)-x)*(corner_y(2)-y) + corner_values(2,1)*(x-corner_x(1))*(corner_y(2)-y) &
		    + corner_values(1,2)*(corner_x(2)-x)*(y-corner_y(1)) + corner_values(2,2)*(x-corner_x(1))*(y-corner_y(1))) &
			/ ((corner_x(2)-corner_x(1))*(corner_y(2)-corner_y(1)))

end function bicubic


subroutine interpolate_between_points(x1, y1, z1, x2, y2, z2, x_out, y_out, z_expected, return_status)

	implicit none
	! direction is from point1 to point2

	double precision, intent(in) :: x1, y1, z1, x2, y2, z2, z_expected
	double precision, intent(out) :: x_out, y_out
	logical, intent(out) :: return_status

	double precision :: angle, intermediate_fraction, dx, dy, distance


	return_status = .true.

	intermediate_fraction = (z_expected - z1) / (z2 - z1)
	
	if(intermediate_fraction > 1.0 .or. intermediate_fraction < 0.0) THEN
		return_status = .false.
	endif


	if(x1 == x2) THEN ! straight up-down line
		x_out = x1
		y_out = y1 + (y2-y1) * intermediate_fraction

	else if (y1 == y2) THEN ! straight line left-right

		x_out = x1 + (x2-x1) * intermediate_fraction
		y_out = y1
	else

		dx = x2-x1
		dy = y2-y1
		distance = z2 - z1
		angle = atan2(dy,dx)
		x_out = x1 + distance * intermediate_fraction * cos(angle)
		y_out = y1 + distance * intermediate_fraction * sin(angle)	
	endif	

!	write(6,*) ">"
!	write(6,*) x1, y1, z1
!	write(6,*) x_out, y_out, z_expected
!	write(6,*) x2, y2, z2


end subroutine interpolate_between_points



subroutine flowline_loop(x_flowline_store,y_flowline_store,distance_store,grid_spacing, r_increment,&
	                   oscillating,hit_saddle,outside,flowline_point_count)

	use direction_grid_mod
	use boundary_mask_mod
	use crossover

	implicit none
	integer, parameter :: max_flowline_points = 100000
	double precision, intent(in) :: grid_spacing, r_increment
	integer, intent(out) :: flowline_point_count
	double precision, dimension(max_flowline_points), intent(inout) :: x_flowline_store, y_flowline_store
	double precision, dimension(max_flowline_points), intent(out) :: distance_store
	logical, intent(out) :: hit_saddle, oscillating, outside

	double precision :: grid_x(2), grid_y(2), distance, dx, dy, total_distance
	integer :: x_grid_index, y_grid_index, x_grid_point, y_grid_point
	logical :: return_status
	double precision, dimension(2,2) :: corner_values, corner_values_x, corner_values_y

	double precision, parameter :: increment_minimum = 10e-5

	total_distance = 0.
	flowline_point_count = 1

	oscillating = .false.
	hit_saddle = .false.
	outside = .false.

	loop: do
		! find grid points

		call find_grid_corner(x_flowline_store(flowline_point_count), y_flowline_store(flowline_point_count),&
		   grid_spacing, x_grid_index, y_grid_index, return_status)



		call find_grid_location(x_grid_index, y_grid_index, grid_spacing, grid_x(1), grid_y(1))

		grid_x(2) = grid_x(1) + grid_spacing
		grid_y(2) = grid_y(1) + grid_spacing


		corner_values(1,1) = direction_grid(x_grid_index, y_grid_index)
		corner_values(2,1) = direction_grid(x_grid_index+1, y_grid_index)
		corner_values(2,2) = direction_grid(x_grid_index+1, y_grid_index+1)
		corner_values(1,2) = direction_grid(x_grid_index, y_grid_index+1)

		corner_values_x = cos(corner_values)
		corner_values_y = sin(corner_values)
		

		dx = bicubic(x_flowline_store(flowline_point_count),y_flowline_store(flowline_point_count),&
			grid_x,grid_y,corner_values_x)
		dy = bicubic(x_flowline_store(flowline_point_count),y_flowline_store(flowline_point_count),&
			grid_x,grid_y,corner_values_y)

		distance = sqrt(dx**2 + dy**2)



		if(distance < increment_minimum) THEN
			! likely hit a discontinuity in the field
			write(6,*) "hit discontinuity"
			hit_saddle = .true.
			exit loop
		endif

		flowline_point_count = flowline_point_count + 1

		if(flowline_point_count > max_flowline_points) THEN
			write(6,*) "possibly oscillating"
			oscillating = .true.
			exit loop
		endif

		x_flowline_store(flowline_point_count) = x_flowline_store(flowline_point_count-1) +&
		  r_increment * dx
		y_flowline_store(flowline_point_count) = y_flowline_store(flowline_point_count-1) +&
		  r_increment * dy




		if(.not. return_status) THEN

			outside =.true.
			exit loop

		endif

		! check to see if the line crosses over the boundary

		end_line = .false.


		call find_grid_index(x_flowline_store(flowline_point_count), y_flowline_store(flowline_point_count),&
		   grid_spacing, x_grid_point, y_grid_point, return_status)

		if(mask(2,x_grid_point, y_grid_point) == 1 ) THEN
			call cross_polygon(x_flowline_store(flowline_point_count-1), &
			  y_flowline_store(flowline_point_count-1), x_flowline_store(flowline_point_count), &
			  y_flowline_store(flowline_point_count), end_line, grid_spacing)

		endif


		total_distance = total_distance + sqrt((x_flowline_store(flowline_point_count)-&
		  x_flowline_store(flowline_point_count-1))**2 + (y_flowline_store(flowline_point_count)&
		  -y_flowline_store(flowline_point_count-1))**2)

		distance_store(flowline_point_count) = total_distance


		if(end_line) then
			exit loop
		endif


	end do loop



end subroutine flowline_loop



subroutine write_flowline(x_flowline_store,y_flowline_store,distance_store,flowline_point_count,gmt_unit,&
				  discard_unit,hit_saddle, oscillating, outside)

	implicit none


	integer, parameter :: max_flowline_points = 100000
	double precision, dimension(max_flowline_points), intent(in) :: x_flowline_store, y_flowline_store,distance_store
	integer, intent(in) :: flowline_point_count,gmt_unit,discard_unit
	logical, intent(in) :: hit_saddle, oscillating, outside

	double precision, dimension(max_flowline_points) :: distance_store_normalized

	integer :: coarse_counter, flow_counter
	integer, parameter :: coarse_factor = 21 ! every 0.05 units

	double precision :: total_distance
	double precision :: x1, y1, z1, x2, y2, z2, x, y
	double precision, dimension(coarse_factor) :: coarse_x, coarse_y
	character (len=1), parameter :: divider_character = ">"

	logical :: 	return_status

	if (oscillating .or. outside  .or. hit_saddle) THEN

		write(discard_unit,'(A1,I7,I7)') divider_character
	else
		write(gmt_unit,'(A1,I7,I7)') divider_character
	endif

	! reduce the amount of points

	total_distance = distance_store(flowline_point_count)

	if(flowline_point_count > 1) THEN
		coarse_counter = 1
		coarse_x(coarse_counter) = x_flowline_store(1)
		coarse_y(coarse_counter) = y_flowline_store(1)
	
		distance_store_normalized = distance_store / total_distance * dble(coarse_factor-1) ! normalizes the distance, then puts it into fractions that can be used for finding the points

		flow_counter = 2
		do while (flow_counter < flowline_point_count)
	
			if(int(distance_store_normalized(flow_counter)) >= coarse_counter) THEN ! add point


				x1 = x_flowline_store(flow_counter-1)
				y1 = y_flowline_store(flow_counter-1)
				z1 = distance_store_normalized(flow_counter-1)

				x2 = x_flowline_store(flow_counter)
				y2 = y_flowline_store(flow_counter)
				z2 = distance_store_normalized(flow_counter)


				call interpolate_between_points(x1, y1, z1, x2, y2, z2, x, y, dble(coarse_counter), return_status)

				if(.not. return_status) THEN
					write(6,*) "possible error in the code"
					write(6,*) x1, y1, z1
					write(6,*) x2, y2, z2
					write(6,*) coarse_counter
				endif

				coarse_counter = coarse_counter + 1
				coarse_x(coarse_counter) = x
				coarse_y(coarse_counter) = y

			else
	!			write(gmt_unit,*) x_flowline_store(flow_counter),y_flowline_store(flow_counter)

				flow_counter = flow_counter + 1

			endif

			

		end do

		coarse_counter = coarse_factor
		coarse_x(coarse_counter) = x_flowline_store(flowline_point_count)
		coarse_y(coarse_counter) = y_flowline_store(flowline_point_count)

		do coarse_counter = 1, coarse_factor, 1

			if(oscillating .or. outside .or. hit_saddle) THEN
				write(discard_unit,*) coarse_x(coarse_counter), coarse_y(coarse_counter)
			else

				write(gmt_unit,*) coarse_x(coarse_counter), coarse_y(coarse_counter)
			endif

		end do

	else
	 	write(discard_unit,*) x_flowline_store(1), y_flowline_store(1)
	endif


end subroutine write_flowline


end program flowlines3
