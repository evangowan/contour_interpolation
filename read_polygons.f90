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
	integer, parameter ::  max_points = 100000
	integer, dimension(2) :: number_polygons
	integer, dimension(:,:), allocatable :: polygon_points

	double precision, dimension(:,:,:), allocatable :: x_coordinates, y_coordinates

	character (len=256), dimension(2) :: step_file


contains

subroutine read_polygons_init()

	implicit none

	integer :: istat, counter, polygon_counter
	integer, parameter :: param_unit=10, step_unit = 20

	character (len=1) :: divider ! make sure that there are no spaces before the ">" character or this will mess up
	character (len=1), parameter :: divider_character = ">", ignore_character = "#"

	open(unit=param_unit, file="params.txt", access="sequential", form="formatted", status="old")

	read(param_unit,*) number_polygons(1), number_polygons(2)
	read(param_unit,*) step_file(1)
	read(param_unit,*) step_file(2)

	close(unit=param_unit)

	max_polygons = max(number_polygons(1),number_polygons(2))
	allocate(x_coordinates(2,max_polygons,max_points), y_coordinates(2,max_polygons,max_points), &
		   polygon_points(2,max_polygons), stat=istat)
	if(istat /=0) THEN
		write(6,*) "problem allocating arrays in read_polygons"
		stop
	endif



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
				if(polygon_points(counter,polygon_counter) > max_points) THEN
					write(6,*) "Number of points in polygon exceeds internal memory"
					write(6,*) "If you want to proceed, you must increase max_points"
					write(6,*) "and recompile"
					close(unit=step_unit)
					stop
				end if

				read(step_unit,*) x_coordinates(counter,polygon_counter,polygon_points(counter,polygon_counter)), &
							y_coordinates(counter,polygon_counter,polygon_points(counter,polygon_counter))
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

end module read_polygons
