module polygon_storage

	double precision, allocatable, dimension(:,:), save :: x_store, y_store
	double precision, allocatable, dimension(:,:) :: distance_start, distance_end
	double precision, allocatable, dimension(:) :: x_temp, y_temp
	logical, allocatable, dimension(:), save :: line_mask
	integer, allocatable, dimension(:), save :: points_store
	integer, save :: number_lines, max_points, total_points

	character(len=80), parameter :: out_file = "out_polygon.gmt"
	integer, parameter :: out_unit = 30
	

contains


subroutine read_lines(file_name)


	implicit none


	character(len=256), intent(in) :: file_name

	integer, parameter :: file_unit = 10

	integer :: istat, points_count, line_count

	character(len=1) :: carrot, hash, test_first

	carrot = ">"
	hash = "#"

	open(unit=file_unit, file=file_name, status="old", form="formatted", access="sequential")

	number_lines = 0
	max_points = 0
	points_count = 0
	read_init: do

		read(file_unit,*,iostat=istat) test_first
		if(istat /= 0) THEN
			exit read_init
		endif

		if (test_first == carrot) THEN
			number_lines = number_lines + 1
			if(points_count > max_points) THEN
				max_points = points_count
			endif

			points_count = 0

		else if (test_first /= hash) THEN

			points_count = points_count + 1

		endif




	end do read_init


	if(points_count > max_points) THEN
		max_points = points_count
	endif


	allocate(x_store(number_lines,max_points), y_store(number_lines,max_points), points_store(number_lines), &
		  line_mask(number_lines), distance_start(number_lines,2), distance_end(number_lines,2), stat=istat)
	if(istat /=0) THEN
		write(6,*) "problem allocating arrays in read_lines"
		stop
	endif


	rewind(unit=file_unit)

	line_mask = .true.

	write(6,*) number_lines, max_points


	x_store = 0
	y_store = 0
	points_store = 0

	line_count = 0
	total_points = 0

	read_second: do

		read(file_unit,*,iostat=istat) test_first
		if(istat /= 0) THEN
			exit read_second
		endif

		if (test_first == carrot) THEN

			if(line_count > 0) THEN
				points_store(line_count) = points_count
			endif

			line_count = line_count + 1
			points_count = 0



		else if (test_first /= hash) THEN

			backspace(unit=file_unit)

			points_count = points_count + 1
			total_points = total_points + 1

			read(file_unit,*) x_store(line_count,points_count), y_store(line_count,points_count) 


		endif


	end do read_second


	points_store(line_count) = points_count

	total_points = total_points * 3 ! just so I don't run into problems later with rounding errors
	allocate(x_temp(total_points), y_temp(total_points), stat=istat)
	if(istat /=0) THEN
		write(6,*) "problem allocating arrays in read_lines"
		stop
	endif


	close(unit = file_unit)

end subroutine read_lines


subroutine read_lines_clear()
! clears the memory from this subroutine

	implicit none

	integer :: istat

	deallocate(x_store, y_store, points_store,line_mask,x_temp,y_temp,distance_start,distance_end, stat=istat)
	if(istat /=0) THEN
		write(6,*) "problem deallocating arrays in read_line"
		stop
	endif

end subroutine read_lines_clear



subroutine check_out_open()

	logical :: open_status

	inquire(unit=out_unit, opened=open_status) 
	if ( .not. open_status ) then
		open(unit=out_unit,file=out_file, status="replace", form="formatted", access="sequential")
		write(out_unit,'(A20)') "# @VGMT1.0 @GPOLYGON"
	end if


end subroutine check_out_open


subroutine close_out()

	close(unit=out_unit)

end subroutine close_out

subroutine write_polygon(polygon_x, polygon_y,polygon_point_count)


	implicit none

	integer, intent(in) :: polygon_point_count
	double precision, dimension(polygon_point_count), intent(in) :: polygon_x, polygon_y

	integer :: counter, end_count


	call check_out_open()


	if(polygon_x(1) == polygon_x(polygon_point_count) .and. polygon_y(1) == polygon_y(polygon_point_count)) THEN
		end_count = polygon_point_count - 1
	else
		end_count = polygon_point_count
	end if


	if (end_count > 2) THEN

		write(out_unit,'(A1)') ">"

		write(out_unit,*) polygon_x(1), polygon_y(1)

		do counter = 2, end_count, 1

			if(polygon_x(counter) /= polygon_x(counter-1) .and. polygon_y(counter) /= polygon_y(counter-1)) THEN
				write(out_unit,*) polygon_x(counter), polygon_y(counter)
			endif

		end do

		write(out_unit,*) polygon_x(1), polygon_y(1)
	
	end if

end subroutine write_polygon


subroutine remove_existing_polygons()

	implicit none
	integer :: counter


	do counter = 1, number_lines, 1

		if(points_store(counter) < 2) then ! ignore

			line_mask(counter) = .false.
		
		else if(close_number(x_store(counter,1),x_store(counter,points_store(counter))) .and. &
			  close_number(y_store(counter,1),y_store(counter,points_store(counter))) ) THEN

			! for debugging, commented out
			call write_polygon(x_store(counter,1:points_store(counter)),y_store(counter,1:points_store(counter)),&
			  points_store(counter))
			line_mask(counter) = .false.
		end if


	end do

end subroutine remove_existing_polygons

logical function close_number(a, b)

	implicit none

	double precision, intent(in) :: a, b

	double precision, parameter :: fraction_threshold = 10e-6

	if(abs(a-b) / a < fraction_threshold) THEN ! assume they are very close
		close_number = .true.
	else
		close_number = .false.
	endif


end function close_number


subroutine create_polygons()

	integer :: temp_points, starting_line, start_index, end_index, counter, min_distance_index
	integer, dimension(4) :: min_index
	double precision, dimension(5) :: min_distance

	logical :: counted_first, incomplete

	x_temp = 0
	y_temp = 0
	temp_points = 0


	! find the first line segement that will serve as the starting point

	starting_line = 1

	find_first: do

		if(line_mask(starting_line) .and. points_store(starting_line) > 2) THEN
			exit find_first
		else
			starting_line = starting_line + 1

		endif

		if(starting_line == number_lines) THEN
			line_mask = .false.
			return
		endif

	end do find_first

 	! hopefully I don't run into issues doing this, might only be a problem if there is only one polygon and encompasses all points
	start_index = total_points / 2 - points_store(starting_line) / 2

	end_index = start_index - 1

	! transfer first line's points to temp array


	do counter = 1, points_store(starting_line)
		end_index = end_index + 1
		x_temp(end_index) = x_store(starting_line,counter)
		y_temp(end_index) = y_store(starting_line,counter)

	end do

	line_mask(starting_line) = .false.

	counted_first = .false.

	incomplete = .true.

	do while(incomplete)

		incomplete = .false.
		do counter = 1, number_lines, 1

			if(line_mask(counter)) THEN

				distance_start(counter,1) = distance(x_temp(start_index),y_temp(start_index),x_store(counter,1),y_store(counter,1))
				distance_start(counter,2) = distance(x_temp(start_index),y_temp(start_index),&
				  x_store(counter,points_store(counter)),y_store(counter,points_store(counter)))

				distance_end(counter,1) = distance(x_temp(end_index),y_temp(end_index),x_store(counter,1),y_store(counter,1))
				distance_end(counter,2) = distance(x_temp(end_index),y_temp(end_index),&
				  x_store(counter,points_store(counter)),y_store(counter,points_store(counter)))


				incomplete = .true. ! if there are available line segments, otherwise the polygon will be complete
			endif

		end do

		if(incomplete) THEN
			min_index(1) = minloc(distance_start(:,1),1,line_mask)	
			min_index(2) = minloc(distance_start(:,2),1,line_mask)
			min_index(3) = minloc(distance_end(:,1),1,line_mask)	
			min_index(4) = minloc(distance_end(:,2),1,line_mask)

			min_distance(1) = distance_start(min_index(1),1)
			min_distance(2) = distance_start(min_index(2),2)
			min_distance(3) = distance_end(min_index(3),1)
			min_distance(4) = distance_end(min_index(4),2)
			min_distance(5) = distance(x_temp(start_index),y_temp(start_index),x_temp(end_index),y_temp(end_index))

			if( end_index - start_index == 1  ) THEN ! should at least link two lines
				min_distance(5) = 1e10
				counted_first = .true.
			endif


			min_distance_index = minloc(min_distance,1)

			if(min_distance_index == 1) THEN ! connects the start point with the start point of the min line

				do counter = 1, points_store(min_index(1)), 1

					start_index = start_index - 1

					x_temp(start_index) = x_store(min_index(1),counter)
					y_temp(start_index) = y_store(min_index(1),counter)
					line_mask(min_index(1)) = .false.


				end do
			else if(min_distance_index == 2) THEN ! connect the start point with the end point of the min line
				do counter =  points_store(min_index(2)), 1, -1

					start_index = start_index - 1

					x_temp(start_index) = x_store(min_index(2),counter)
					y_temp(start_index) = y_store(min_index(2),counter)
					line_mask(min_index(2)) = .false.


				end do
			else if(min_distance_index == 3) THEN ! connects the end point with the start point of the min line

				do counter = 1, points_store(min_index(3)), 1

					end_index = end_index + 1

					x_temp(end_index) = x_store(min_index(3),counter)
					y_temp(end_index) = y_store(min_index(3),counter)
					line_mask(min_index(3)) = .false.


				end do
			else if(min_distance_index == 4) THEN ! connects the end point with the end point of the min line

				do counter = points_store(min_index(4)), 1, -1

					end_index = end_index + 1

					x_temp(end_index) = x_store(min_index(4),counter)
					y_temp(end_index) = y_store(min_index(4),counter)
					line_mask(min_index(4)) = .false.


				end do

			else if(min_distance_index == 5) THEN ! complete polygon

				write(6,*) "completed polygon"
				incomplete = .false.


			else
				write(6,*) "something messed up"

			endif


		endif
 
	end do
	
	end_index = end_index + 1
	x_temp(end_index) = x_temp(start_index)
	y_temp(end_index) = y_temp(start_index)

	write(6,*) x_temp(start_index), y_temp(start_index)
	write(6,*) x_temp(end_index), y_temp(end_index)

	call write_polygon(x_temp(start_index:end_index),y_temp(start_index:end_index),end_index-start_index+1)


end subroutine create_polygons


double precision function distance(x1,y1,x2,y2)

	double precision, intent(in) :: x1, y1, x2, y2

	distance = sqrt((x2-x1)**2+(y2-y1)**2)

end function distance

end module polygon_storage


program merge_polygons

	use polygon_storage
	implicit none

	integer :: commandline_count

	character(len=256) :: file_name

	character(len=256) :: arg_in


	commandline_count = command_argument_count()

	if(commandline_count /= 1) THEN

		write(6,*) "correct usage of this program:"
		write(6,*) "merge_polygons filename"
		stop
	endif

	call GET_COMMAND_ARGUMENT(1,arg_in) 
	read(arg_in,*) file_name


	call read_lines(file_name)

	call remove_existing_polygons()


	do while (any(line_mask))

		call create_polygons()
	end do

	call read_lines_clear()


end program merge_polygons
