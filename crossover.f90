module crossover


	use read_polygons
! make sure the main program includes read_polygons.o and read in those polygons

contains 

subroutine cross_polygon(current_x, current_y, next_x, next_y, end_line, grid_spacing,compare_in)



	! find out if there is an overlap between the current line segment and the inner polygons


	implicit none

	double precision, intent(in) :: current_x, current_y, grid_spacing
	double precision, intent(inout) :: next_x, next_y
	logical, intent(inout) :: end_line
   	integer,optional :: compare_in

	integer :: polygon_set


	integer :: polygon_counter, points_counter, next_index


	double precision :: seg_x1, seg_x2, seg_y1, seg_y2, temp_x, temp_y

	logical :: is_crossover



	if(present(compare_in))then
	  polygon_set=compare_in
	else
	  polygon_set=2
	endif

	end_line = .false.


	polygon_loop: do polygon_counter = 1, number_polygons(polygon_set), 1
		points_loop: do points_counter = 1, polygon_points(polygon_set, polygon_counter), 1

			next_index = points_counter + 1
			if(next_index > polygon_points(polygon_set, polygon_counter)) then
				next_index = 1
			endif

			seg_x1 = x_coordinates(polygon_set,polygon_counter,points_counter)
			seg_y1 = y_coordinates(polygon_set,polygon_counter,points_counter)

			seg_x2 = x_coordinates(polygon_set,polygon_counter,next_index)
			seg_y2 = y_coordinates(polygon_set,polygon_counter,next_index)


			call crossover_point(current_x, current_y, next_x, next_y, seg_x1, seg_y1, seg_x2, seg_y2, &
			  temp_x, temp_y, is_crossover)

			if( is_crossover ) THEN ! crossover
				next_x = temp_x
				next_y = temp_y
				end_line = .true.
			!	write(546,*) next_x, next_y
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
	double precision :: temp_x, temp_y


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





	if(a_x2 == a_x1) then
		temp_x = a_x2
		temp_y = slope2 * temp_x + intercept2
	elseif(b_x2 == b_x1) THEN
		temp_x = b_x2
		temp_y = slope1 * temp_x + intercept1
	else 

		slope1 = (a_y2 - a_y1) / (a_x2 - a_x1) 
		intercept1 = a_y2 - a_x2 * slope1

		slope2 = (b_y2 - b_y1) / (b_x2 - b_x1) 
		intercept2 = b_y2 - b_x2 * slope2

		temp_x = (intercept2 - intercept1) / (slope1 - slope2)
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


end module crossover
