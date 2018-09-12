program flowlines3

! Calculates the path from the outer boundary to the inner boundary through the vector field defined in interpolate_direction
! currently very slow

! read in the boundaries, boundary masks, direction field and interpolate between the boundaries

	use read_polygons
	use read_minmax
	use boundary_mask_mod
	use direction_grid_mod


	implicit none

	integer :: commandline_count, polygon_counter, points_counter, dummy_counter
	integer ::  index_next, extra_counter, extra_total, flowline_point_count
	integer, parameter :: gmt_unit = 90, max_flowline_points = 100000, discard_unit=200



	double precision :: grid_spacing

	
	double precision ::  r_increment
	double precision, parameter :: r_increment_minimum = 0.1

	double precision, dimension(max_flowline_points) :: x_flowline_store, y_flowline_store, distance_store



	character(len=256) :: arg_in

	character(len=256), parameter :: gmt_file = "flowlines.txt", discard_file="discarded_flowlines.txt"


	logical :: found_point, end_line, round_down_x, round_down_y, looping, hit_saddle, oscillating, return_status, outside
	logical :: reverse
	integer :: polygon_set, polygon_compare


	commandline_count = command_argument_count()

	if(commandline_count /= 2) THEN
		write(6,*) "not enough args: ", commandline_count
		write(6,*) "correct usage of this program:"
		write(6,*) "flowlines3 grid_spacing reverse(T/F)"
		stop
	endif

	call GET_COMMAND_ARGUMENT(1,arg_in) 
	read(arg_in,*) grid_spacing

	call GET_COMMAND_ARGUMENT(2,arg_in) 

	if(trim(adjustl(arg_in)) == "T") THEN
		reverse = .true.
		polygon_set = 2
		polygon_compare = 1
	elseif(trim(adjustl(arg_in)) == "F") THEN
		reverse = .false.
		polygon_set = 1
		polygon_compare = 2
	else
		write(6,*) "problem in argument 2: ", trim(adjustl(arg_in))
		write(6,*) "correct usage of this program:"
		write(6,*) "flowlines grid_spacing reverse(T/F)"
		stop
	endif

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
	do polygon_counter = 1, number_polygons(polygon_set)

		do points_counter = 1, polygon_points(polygon_set, polygon_counter)
			write(6,*) "polygon: ", polygon_counter, "/", number_polygons(polygon_set), " point:", points_counter, "/", &
				polygon_points(polygon_set, polygon_counter), reverse
			if(points_counter < polygon_points(polygon_set, polygon_counter) ) then
				index_next = points_counter + 1
			else
				index_next = 1
			endif




			flowline_point_count = 1

			x_flowline_store(flowline_point_count) = x_coordinates(polygon_set,polygon_counter,points_counter) 
			y_flowline_store(flowline_point_count) = y_coordinates(polygon_set,polygon_counter,points_counter)

!			call flowline_loop(x_flowline_store,y_flowline_store,distance_store,grid_spacing, r_increment,&
!	                   oscillating,hit_saddle,outside,flowline_point_count,reverse, polygon_compare)
			call flowline_loop_runga_kutta(x_flowline_store,y_flowline_store,distance_store,grid_spacing, r_increment,&
	                   oscillating,hit_saddle,outside,flowline_point_count,reverse, polygon_compare)


			if(flowline_point_count > 2) THEN

				call write_flowline(x_flowline_store,y_flowline_store,distance_store,flowline_point_count,gmt_unit,&
				  discard_unit,hit_saddle, oscillating, outside,reverse)
			else
				write(6,*) "did not proceed far enough"
			endif
			


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

!	write(667,*) "intermediate_fraction:", intermediate_fraction

	if(x1 == x2) THEN ! straight up-down line
		x_out = x1
		y_out = y1 + (y2-y1) * intermediate_fraction

	else if (y1 == y2) THEN ! straight line left-right

		x_out = x1 + (x2-x1) * intermediate_fraction
		y_out = y1
	else

		dx = x2-x1
		dy = y2-y1
		distance = sqrt(dx**2+dy**2)
		angle = atan2(dy,dx)


		x_out = x1 + distance * intermediate_fraction * cos(angle)
		y_out = y1 + distance * intermediate_fraction * sin(angle)	

!		write(667,*) "dx:", dx
!		write(667,*) "dy:", dy
!		write(667,*) "distance:", distance
!		write(667,*) "angle:", angle
!		write(667,*) "x_distance", distance * intermediate_fraction * cos(angle)
!		write(667,*) "y_distance", distance * intermediate_fraction * sin(angle)	
	endif	

!	write(6,*) ">"
!	write(6,*) x1, y1, z1
!	write(6,*) x_out, y_out, z_expected
!	write(6,*) x2, y2, z2


end subroutine interpolate_between_points


subroutine flowline_loop_runga_kutta(x_flowline_store,y_flowline_store,distance_store,grid_spacing, r_increment,&
	                   oscillating,hit_saddle,outside,flowline_point_count,reverse, polygon_compare)

! found on these online notes
! http://graphics.cs.ucdavis.edu/~joy/ecs277/other-notes/Numerical-Methods-for-Particle-Tracing-in-Vector-Fields.pdf

	use direction_grid_mod
	use boundary_mask_mod
	use crossover

	implicit none
	integer, parameter :: max_flowline_points = 100000
	double precision, intent(in) :: grid_spacing, r_increment
	integer, intent(out) :: flowline_point_count
	double precision, dimension(max_flowline_points), intent(inout) :: x_flowline_store, y_flowline_store
	double precision, dimension(max_flowline_points), intent(out) :: distance_store
	logical, intent(in) :: reverse
	logical, intent(out) :: hit_saddle, oscillating, outside
	integer, intent(in) :: polygon_compare

	integer :: oscillation_check

	double precision, parameter :: pi = 3.141592653589793

	double precision :: grid_x(2), grid_y(2), distance, dx, dy, total_distance, dummy_x, dummy_y, angle
	integer :: x_grid_index, y_grid_index, x_grid_point, y_grid_point, counter1, counter2
	integer :: dummy_x_grid_index, dummy_y_grid_index
	logical :: return_status, dummy_return_status
	double precision, dimension(2,2) :: corner_values, corner_values_x, corner_values_y

	double precision :: k1_x, k1_y, k2_x, k2_y, k3_x, k3_y, k4_x, k4_y
	double precision :: temp_x, temp_y

	double precision :: increment_minimum



	integer :: x_index_check, y_index_check

	integer :: previous_mask

	increment_minimum = r_increment * 1.0e-5

	total_distance = 0.
	flowline_point_count = 1

	oscillating = .false.
	hit_saddle = .false.
	outside = .false.

	x_index_check = 0
	y_index_check = 0

	if(reverse) THEN
		write(6,*) "detected reverse"
	else
		write(6,*) "normal direction"
	endif

	previous_mask = 0

	loop: do
		! find grid points


		! first step of the runga kutta algorithm

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
		if(reverse) then
			dx = -dx
			dy = -dy
		endif


		angle = atan2(dy,dx)
		
		k1_x = r_increment * cos(angle)
		k1_y = r_increment * sin(angle)

		! second step of the runga kutta algorithm

		temp_x = x_flowline_store(flowline_point_count) + k1_x / 2.0
		temp_y = y_flowline_store(flowline_point_count) + k1_y / 2.0

		call find_grid_corner(temp_x, temp_y,&
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
		if(reverse) then
			dx = -dx
			dy = -dy
		endif

		angle = atan2(dy,dx)

		k2_x = r_increment * cos(angle)
		k2_y = r_increment * sin(angle)

		! third step of the Runga Kutta

		temp_x = x_flowline_store(flowline_point_count) + k2_x / 2.0
		temp_y = y_flowline_store(flowline_point_count) + k2_y / 2.0


		call find_grid_corner(temp_x, temp_y,&
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
		if(reverse) then
			dx = -dx
			dy = -dy
		endif

		angle = atan2(dy,dx)

		k3_x = r_increment * cos(angle)
		k3_y = r_increment * sin(angle)
		

		! fourth step of the Runga Kutta

		temp_x = (x_flowline_store(flowline_point_count) + k3_x)
		temp_y = (y_flowline_store(flowline_point_count) + k3_y) 


		call find_grid_corner(temp_x, temp_y,&
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
		if(reverse) then
			dx = -dx
			dy = -dy
		endif

		angle = atan2(dy,dx)
		k4_x = r_increment * cos(angle)
		k4_y = r_increment * sin(angle)


		flowline_point_count = flowline_point_count + 1

		if(flowline_point_count > max_flowline_points) THEN
			write(6,*) "possibly oscillating"
			oscillating = .true.
			exit loop
		endif

		x_flowline_store(flowline_point_count) = x_flowline_store(flowline_point_count-1) + (k1_x + 2.0*k2_x + 2.0*k3_x + k4_x) / 6.0
		y_flowline_store(flowline_point_count) = y_flowline_store(flowline_point_count-1) + (k1_y + 2.0*k2_y + 2.0*k3_y + k4_y) / 6.0


		! check to see if the line crosses over the boundary

		! TODO This is not working properly

		end_line = .false.


		call find_grid_index(x_flowline_store(flowline_point_count), y_flowline_store(flowline_point_count),&
		   grid_spacing, x_grid_point, y_grid_point, return_status)

		if(mask(polygon_compare,x_grid_point, y_grid_point) == 1 .or. previous_mask == 1 ) THEN
			call cross_polygon(x_flowline_store(flowline_point_count-1), &
			  y_flowline_store(flowline_point_count-1), x_flowline_store(flowline_point_count), &
			  y_flowline_store(flowline_point_count), end_line, grid_spacing,polygon_compare)

		endif
 

		total_distance = total_distance + sqrt((x_flowline_store(flowline_point_count)-&
		  x_flowline_store(flowline_point_count-1))**2 + (y_flowline_store(flowline_point_count)&
		  -y_flowline_store(flowline_point_count-1))**2)



		distance_store(flowline_point_count) = total_distance


		previous_mask = mask(polygon_compare,x_grid_point, y_grid_point)


		if(.not. return_status .or. (mask(1,x_grid_point, y_grid_point) == 0 .and. mask(2,x_grid_point, y_grid_point) == 0)) THEN

			outside =.true.
			exit loop

		endif


		if(mask(1,x_grid_point, y_grid_point) == 2 .and. mask(2,x_grid_point, y_grid_point) == 2) THEN
			! inside both
			outside =.true.
			exit loop

		endif


		if(end_line) then
			exit loop
		endif

		! checks to see if the point crosses over itself, indicating it hit a peak

		call cross_itself(x_flowline_store(1:flowline_point_count), y_flowline_store(1:flowline_point_count), &
			 flowline_point_count, end_line)

		if(end_line) then
			write(6,*) "crossed over"
			exit loop
		endif


	end do loop

	write(6,*) "number of points:", flowline_point_count

end subroutine flowline_loop_runga_kutta


subroutine flowline_loop(x_flowline_store,y_flowline_store,distance_store,grid_spacing, r_increment,&
	                   oscillating,hit_saddle,outside,flowline_point_count,reverse, polygon_compare)

	use direction_grid_mod
	use boundary_mask_mod
	use crossover

	implicit none
	integer, parameter :: max_flowline_points = 100000
	double precision, intent(in) :: grid_spacing, r_increment
	integer, intent(out) :: flowline_point_count
	double precision, dimension(max_flowline_points), intent(inout) :: x_flowline_store, y_flowline_store
	double precision, dimension(max_flowline_points), intent(out) :: distance_store
	logical, intent(in) :: reverse
	logical, intent(out) :: hit_saddle, oscillating, outside
	integer, intent(in) :: polygon_compare

	double precision :: grid_x(2), grid_y(2), distance, dx, dy, total_distance, dummy_x, dummy_y
	integer :: x_grid_index, y_grid_index, x_grid_point, y_grid_point, counter1, counter2
	integer :: dummy_x_grid_index, dummy_y_grid_index
	logical :: return_status, dummy_return_status
	double precision, dimension(2,2) :: corner_values, corner_values_x, corner_values_y

	double precision :: increment_minimum

	double precision, parameter :: pi = 3.141592653589793


	integer :: x_index_check, y_index_check

	integer :: previous_mask

	increment_minimum = r_increment * 1.0e-5

	total_distance = 0.
	flowline_point_count = 1

	oscillating = .false.
	hit_saddle = .false.
	outside = .false.

	x_index_check = 0
	y_index_check = 0

	if(reverse) THEN
		write(6,*) "detected reverse"
	else
		write(6,*) "normal direction"
	endif

	previous_mask = 0

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

		if(reverse) then
			dx = -dx
			dy = -dy
		endif


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


		! check to see if the line crosses over the boundary

		end_line = .false.


		call find_grid_index(x_flowline_store(flowline_point_count), y_flowline_store(flowline_point_count),&
		   grid_spacing, x_grid_point, y_grid_point, return_status)

		if(mask(polygon_compare,x_grid_point, y_grid_point) == 1 .or. previous_mask == 1 ) THEN
			call cross_polygon(x_flowline_store(flowline_point_count-1), &
			  y_flowline_store(flowline_point_count-1), x_flowline_store(flowline_point_count), &
			  y_flowline_store(flowline_point_count), end_line, grid_spacing,polygon_compare)

		endif
 

		total_distance = total_distance + sqrt((x_flowline_store(flowline_point_count)-&
		  x_flowline_store(flowline_point_count-1))**2 + (y_flowline_store(flowline_point_count)&
		  -y_flowline_store(flowline_point_count-1))**2)



		distance_store(flowline_point_count) = total_distance


		previous_mask = mask(polygon_compare,x_grid_point, y_grid_point)


		if(.not. return_status .or. (mask(1,x_grid_point, y_grid_point) == 0 .and. mask(2,x_grid_point, y_grid_point) == 0)) THEN

			outside =.true.
			exit loop

		endif


		if(mask(1,x_grid_point, y_grid_point) == 2 .and. mask(2,x_grid_point, y_grid_point) == 2) THEN
			! inside both
			outside =.true.
			exit loop

		endif


		if(end_line) then
			exit loop
		endif




	end do loop

!	stop

end subroutine flowline_loop



subroutine write_flowline(x_flowline_store,y_flowline_store,distance_store,flowline_point_count,gmt_unit,&
				  discard_unit,hit_saddle, oscillating, outside, reverse)

	implicit none


	integer, parameter :: max_flowline_points = 100000
	double precision, dimension(max_flowline_points), intent(in) :: x_flowline_store, y_flowline_store,distance_store
	integer, intent(in) :: flowline_point_count,gmt_unit,discard_unit
	logical, intent(in) :: hit_saddle, oscillating, outside
	logical, intent(in) :: reverse
	double precision, dimension(max_flowline_points) :: distance_store_normalized

	integer :: coarse_counter, flow_counter
	integer, parameter :: coarse_factor = 21 ! every 0.05 units

	double precision :: total_distance
	double precision :: x1, y1, z1, x2, y2, z2, x, y
	double precision, dimension(coarse_factor) :: coarse_x, coarse_y
	character (len=1), parameter :: divider_character = ">"

	logical :: 	return_status, kill

!	if (oscillating .or. outside  .or. hit_saddle) THEN
	if (oscillating .or. outside ) THEN
		write(discard_unit,'(A1,1X,L1,1X,L1,1X,L1)') divider_character, oscillating, outside, hit_saddle
	else
		write(gmt_unit,'(A1)') divider_character
	endif
	kill = .false.
	! reduce the amount of points


	total_distance = distance_store(flowline_point_count)

	if(flowline_point_count > 1) THEN
		coarse_counter = 1
		coarse_x(coarse_counter) = x_flowline_store(1)
		coarse_y(coarse_counter) = y_flowline_store(1)
	
		distance_store_normalized = distance_store / total_distance * dble(coarse_factor-1) ! normalizes the distance, then puts it into fractions that can be used for finding the points

		if(int(distance_store_normalized(flowline_point_count)) > 20) THEN
			write(6,*) "normalization failed"
			write(6,*) distance_store_normalized(flowline_point_count)
			write(6,*) distance_store(flowline_point_count)
			write(6,*) total_distance
			write(6,*) dble(coarse_factor-1)
			stop
		endif




		flow_counter = 2
		do while (flow_counter <= flowline_point_count)

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


				flow_counter = flow_counter + 1

			endif

			

		end do

		coarse_counter = coarse_factor
		coarse_x(coarse_counter) = x_flowline_store(flowline_point_count)
		coarse_y(coarse_counter) = y_flowline_store(flowline_point_count)




		if(reverse) THEN
			do coarse_counter = coarse_factor, 1, -1
				if(coarse_y(coarse_counter) < 40.) THEN
					kill = .true.
				endif

!				if(oscillating .or. outside .or. hit_saddle) THEN
				if(oscillating .or. outside) THEN
					write(discard_unit,*) coarse_x(coarse_counter), coarse_y(coarse_counter)
				else

					write(gmt_unit,*) coarse_x(coarse_counter), coarse_y(coarse_counter), &
						dble(coarse_factor-coarse_counter) / dble(coarse_factor-1)
				endif

			end do

		else
			do coarse_counter = 1, coarse_factor, 1
				if(coarse_y(coarse_counter) < 40.) THEN
					kill = .true.
				endif

!				if(oscillating .or. outside .or. hit_saddle) THEN
				if(oscillating .or. outside) THEN
					write(discard_unit,*) coarse_x(coarse_counter), coarse_y(coarse_counter)
				else

					write(gmt_unit,*) coarse_x(coarse_counter), coarse_y(coarse_counter), &
					  dble(coarse_counter-1) / dble(coarse_factor-1)
				endif

			end do
		endif

	else
	 	write(discard_unit,*) x_flowline_store(1), y_flowline_store(1)
	endif



	if(kill) THEN
		write(6,*) "killing"
	!	stop
	endif

	if(oscillating) THEN
		write(6,*) "oscillating"
	!	stop
	endif
	if(outside) THEN
		write(6,*) "outside"
	!	stop
	endif

end subroutine write_flowline


end program flowlines3
