program contour_creation

	use read_polygons
	use spline 

	implicit none
	double precision :: distance_factor, target_distance, current_distance, distance_add, angle, x, y, part_distance
	double precision :: contour_density, distance, slope, t, interval

	double precision, dimension(:), allocatable :: total_distance
	double precision, dimension(:,:), allocatable :: x_array, y_array, x_final, y_final

	double precision, dimension(spline_points) :: x_points, y_points

	integer :: commandline_count, number_lines, istat, max_points_contour, points_counter, line, line_step, step_counter, next_index
	integer :: half_size, fill_counter, temp_index, add_points, add_counter
	integer, parameter :: gmt_unit = 90, out_unit=100, spline_unit=110
	integer, dimension(:), allocatable :: points_array

	character(len=1) :: test_divider
	character(len=256) :: arg_in

	character (len=1), parameter :: divider_character = ">"
	character(len=256), parameter :: gmt_file = "flowlines.txt", out_file="final_contour.txt", spline_file="spline_contour.txt"


	commandline_count = command_argument_count()

	if(commandline_count /= 1) THEN
		write(6,*) "correct usage of this program:"
		write(6,*) "contour_creation distance_factor"
		stop
	endif

	! the density of points along the flowline. Note it must be between 0 and 1
	call GET_COMMAND_ARGUMENT(1,arg_in) 
	read(arg_in,*) distance_factor

	if (distance_factor >= 1 .or. distance_factor <= 0) THEN
		write(6,*) "distance_factor must be between 0 and 1"
		stop
	endif

	call read_params() ! just need to get the fining_increment

	line_step = int(1./distance_factor)

	open(unit=gmt_unit, file=gmt_file, access="sequential", form="formatted", status="old")

	! find the number of lines and points
	number_lines = 0
	max_points_contour = 0
	read_file1: do

		read(gmt_unit,*,iostat=istat) test_divider
		if(istat /= 0) THEN
			exit read_file1
		endif

		if(test_divider == divider_character) THEN
			number_lines = number_lines + 1
			points_counter = 0
		else
			points_counter = points_counter + 1
			if(points_counter > max_points_contour) THEN
				max_points_contour = points_counter
			endif
		endif

	end do read_file1
	
	rewind(unit=gmt_unit)

	! allocate memory

	allocate(x_array(number_lines,max_points_contour), y_array(number_lines,max_points_contour), points_array(number_lines), &
		total_distance(number_lines), x_final(number_lines,line_step-1), y_final(number_lines,line_step-1), stat=istat)
	if (istat /= 0) THEN
		WRITE(6,*) "allocation error for contour_creation", istat
		stop
	ENDIF

	! read each flowline

	line = 0
	read_file2: do

		read(gmt_unit,*,iostat=istat) test_divider
		if(istat /= 0) THEN
			exit read_file2
		endif

		if(test_divider == divider_character) THEN
			line = line + 1
			points_counter = 0

		else
			backspace(unit=gmt_unit)
			points_counter = points_counter + 1

			read(gmt_unit,*) x_array(line,points_counter), y_array(line,points_counter)
			points_array(line) = points_counter

		endif

	end do read_file2

	close(gmt_unit)

	! determine the total distance for each flowline



	total_distance = 0

	do line = 1, number_lines
	
		if(points_array(line) >= 2) THEN
			do points_counter = 2, points_array(line)

				total_distance(line) = total_distance(line) + sqrt((x_array(line,points_counter)- &
				    x_array(line,points_counter-1))**2 + (y_array(line,points_counter)-y_array(line,points_counter-1))**2)
	
			end do

		endif

	end do

	! find the point on the flowline at the target distance

	open(unit=out_unit, file=out_file, access="sequential", form="formatted", status="replace")

	do line = 1, number_lines

		write(out_unit,*) x_array(line,1), y_array(line,1), 0.

		do step_counter = 1, line_step-1, 1

			interval = (distance_factor * dble(step_counter))

			target_distance = total_distance(line) * interval


			if(target_distance > 0.) THEN
				current_distance = 0.

				go_the_distance: do points_counter = 2, points_array(line)
					distance_add = sqrt((x_array(line,points_counter)-x_array(line,points_counter-1))**2 &
						         + (y_array(line,points_counter)-y_array(line,points_counter-1))**2)

					if(current_distance + distance_add > target_distance) THEN

						part_distance = target_distance - current_distance

						angle = atan2((y_array(line,points_counter) - y_array(line,points_counter-1)), &
							  (x_array(line,points_counter) - x_array(line,points_counter-1)))
						x = cos(angle) * part_distance + x_array(line,points_counter-1)
						y = sin(angle) * part_distance + y_array(line,points_counter-1)


						exit go_the_distance


					else

						current_distance = current_distance + distance_add
					endif

				end do go_the_distance

			else

				x = x_array(line,1)
				y = y_array(line,1)

			end if

			x_final(line,step_counter) = x
			y_final(line,step_counter) = y

			write(out_unit,*) x, y, interval
		end do

		write(out_unit,*) x_array(line,points_array(line)), y_array(line,points_array(line)), 1.

	end do

	close(unit=out_unit)


	! loop through the results, and produce contours using splines

	open(unit=spline_unit, file=spline_file, access="sequential", form="formatted", status="replace")

	half_size = int(spline_points / 2)

	if (number_lines >= spline_points) THEN ! doesn't work if there are fewer lines than the number of spline points, probably should change later to be allocatable

		do step_counter = 1, line_step - 1, 1

			interval = (distance_factor * dble(step_counter))
			write(spline_unit,'(A1,F5.3)') divider_character, interval
			do line = 1, number_lines

				if(line < number_lines) THEN
					next_index = 1  
				else
					next_index = line + 1
				endif


				write(spline_unit,*) x_final(line,step_counter), y_final(line,step_counter), interval

				distance = sqrt((x_final(line,step_counter) - x_final(next_index,step_counter))**2 +&
					          (y_final(line,step_counter) - y_final(next_index,step_counter))**2)
				

				if(distance > fining_increment) THEN ! add some points

					! determine how many points should be added

					add_points = int((distance*2.)/fining_increment)

					write(531,*) ">", x_final(line,step_counter), y_final(line,step_counter), distance, add_points

					! populate temperary arrays
				
					do fill_counter = 1, spline_points

						temp_index = line + (fill_counter - half_size)

						if(temp_index < 1) then
							temp_index = number_lines +  (fill_counter - half_size)

						endif

						x_points(fill_counter) = x_final(temp_index,step_counter)
						y_points(fill_counter) = y_final(temp_index,step_counter)

						if(fill_counter == half_size) THEN
							write(531,*) fill_counter, x_points(fill_counter), y_points(fill_counter), "#"
						else
							write(531,*) fill_counter, x_points(fill_counter), y_points(fill_counter)
						endif

					end do



					do add_counter = 1, add_points

						t = dble(add_counter) / dble(add_points + 1)

						call cubic_spline(x_points,y_points,t,x,y,slope) ! slope isn't used here, though
						write(spline_unit,*) x, y, interval, " >" 

					end do

				endif

			end do
		end do
	end if


	deallocate(x_array, y_array, points_array, total_distance, x_final, y_final, stat=istat)
	if (istat /= 0) THEN
		WRITE(6,*) "deallocation error for contour_creation", istat
		stop
	ENDIF

end program contour_creation
