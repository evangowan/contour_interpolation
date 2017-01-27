module read_polygons

! reads the GMT formatted multisegment files that contain the starting and ending polygons.

! required: file called "params.txt" that contains the number of polygons for each step and the file names of a GMT-style multi-segment
! files (where the individual polygons are separated by a "greater than" character > ) corresponding to those steps. See below:
!
! 1 2
! step1.txt
! step2.txt

	implicit none

	integer :: max_polygons
	integer, parameter ::  max_points = 100000, step_unit = 20, param_unit=10
	integer, dimension(2) :: number_polygons
	integer, dimension(:,:), allocatable :: polygon_points

	double precision, dimension(:,:,:), allocatable :: x_coordinates, y_coordinates
	double precision, save :: fining_increment
	character (len=256), dimension(2) :: step_file


contains

subroutine read_params()

	open(unit=param_unit, file="params.txt", access="sequential", form="formatted", status="old")

	read(param_unit,*) number_polygons(1), number_polygons(2)
	read(param_unit,*) step_file(1)
	read(param_unit,*) step_file(2)
	read(param_unit,*) fining_increment

	close(unit=param_unit)

end subroutine read_params

subroutine read_polygons_init()

	implicit none

	integer :: istat, counter, polygon_counter, add_points, add_counter


	double precision :: x, y, distance, angle

	character (len=1) :: divider ! make sure that there are no spaces before the ">" character or this will mess up
	character (len=1), parameter :: divider_character = ">", ignore_character = "#"



	max_polygons = max(number_polygons(1),number_polygons(2))
	allocate(x_coordinates(2,max_polygons,max_points), y_coordinates(2,max_polygons,max_points), &
		   polygon_points(2,max_polygons), stat=istat)
	if(istat /=0) THEN
		write(6,*) "problem allocating arrays in read_polygons"
		stop
	endif

	call read_params()

	polygon_points = 0
	read_files: do counter = 1, 2, 1

		open(unit=step_unit, file=step_file(counter), access="sequential", form="formatted", status="old")
		polygon_counter = 0


		read_polygons: do 
			
			read(step_unit,*, iostat=istat) divider
			if(istat /=0) THEN
				exit read_polygons
			endif

			if(divider == divider_character) THEN
				polygon_counter = polygon_counter + 1
				cycle read_polygons
			elseif (divider == ignore_character) THEN ! anything starting with '#' should be ignored
				cycle read_polygons
			else
				backspace(unit=step_unit)

				polygon_points(counter,polygon_counter) = polygon_points(counter,polygon_counter) + 1
				call check_array(polygon_points(counter,polygon_counter))

				read(step_unit,*) x, y
				if (polygon_points(counter,polygon_counter) > 1) THEN ! check if points need to be added

					distance = sqrt(&
					  (x_coordinates(counter,polygon_counter,polygon_points(counter,polygon_counter)-1)-x)**2 + &
					  (y_coordinates(counter,polygon_counter,polygon_points(counter,polygon_counter)-1)-y)**2)

					if(distance > fining_increment) THEN ! add points
					! currently this just does linear interpolation, probably better in the future to do a spline

						add_points = int(distance / fining_increment)

						angle = atan2(&
						 y - y_coordinates(counter,polygon_counter,polygon_points(counter,polygon_counter)-1), &
						 x - x_coordinates(counter,polygon_counter,polygon_points(counter,polygon_counter)-1))

						do add_counter = 1, add_points

							

							x_coordinates(counter,polygon_counter,polygon_points(counter,polygon_counter)) = &
						x_coordinates(counter,polygon_counter,polygon_points(counter,polygon_counter)-add_counter) + &
							  cos(angle) * dble(add_counter) / dble(add_points+1) * distance
							y_coordinates(counter,polygon_counter,polygon_points(counter,polygon_counter)) = &
						y_coordinates(counter,polygon_counter,polygon_points(counter,polygon_counter)-add_counter) + &
							  sin(angle) * dble(add_counter) / dble(add_points+1) * distance


							polygon_points(counter,polygon_counter) = polygon_points(counter,polygon_counter) + 1
							call check_array(polygon_points(counter,polygon_counter))

						end do
					endif


				end if

				x_coordinates(counter,polygon_counter,polygon_points(counter,polygon_counter)) = x
				y_coordinates(counter,polygon_counter,polygon_points(counter,polygon_counter)) = y
			endif

		end do read_polygons


		close(unit=step_unit)
	end do read_files

! check to see if the start point is the same as the end point. The point_in_polygon routine does not work if it is, so this loop removes it

	do counter = 1, 2, 1

		do polygon_counter = 1, number_polygons(counter), 1

			if (x_coordinates(counter,polygon_counter,1) == &
			    x_coordinates(counter,polygon_counter,polygon_points(counter,polygon_counter)) .and. &
			    y_coordinates(counter,polygon_counter,1) == &
			    y_coordinates(counter,polygon_counter,polygon_points(counter,polygon_counter))) THEN

				polygon_points(counter,polygon_counter) = polygon_points(counter,polygon_counter) - 1

			endif

		end do

	end do

end subroutine read_polygons_init

subroutine read_polygons_clear()
! clears the memory from this subroutine

	implicit none

	integer :: istat

	deallocate(x_coordinates, y_coordinates, polygon_points, stat=istat)
	if(istat /=0) THEN
		write(6,*) "problem deallocating arrays in read_polygons"
		stop
	endif

end subroutine read_polygons_clear

subroutine check_array(point_count)

	integer, intent(in) :: point_count

	if(point_count > max_points) THEN
		write(6,*) "Number of points in polygon exceeds internal memory"
		write(6,*) "If you want to proceed, you must increase max_points"
		write(6,*) "and recompile"
		close(unit=step_unit)
		stop
	end if

end subroutine check_array

end module read_polygons
