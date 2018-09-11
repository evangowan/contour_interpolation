!	copyright Evan. J. Gowan, 2016
!
!	This file is part of countour_interpolation
!
!	ICESHEET is free software: you can redistribute it and/or modify
!	it under the terms of the GNU General Public License as published by
!	the Free Software Foundation, version 3 of the License
!
!	ICESHEET is distributed in the hope that it will be useful,
!	but WITHOUT ANY WARRANTY; without even the implied warranty of
!	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!	GNU General Public License for more details.
!
!	You should have received a copy of the GNU General Public License
!	along with countour_interpolation.  If not, see <http://www.gnu.org/licenses/>.

!***********************************************************************************

program find_direction

! This program takes one or more input polygons from two steps (could be time steps for ice sheet margins, or successive elevation contours)
! and finds the direction that they should be going. Notably, the input file will have two sets of polygons, one representing the initial step
! and the second is the second step (i.e. for an ice sheet boundary, time slice #1 and time slice #2). In order for this to work, polygons
! in a single set cannot overlap in any way, though they can overlap between the sets.

! The algorithm works in this way:
! 1) reads in all the polygons for each set
! 2) does a "point in polygon" routine for each polygon in each step. This loop will be nxm steps, where n and m are the number of polygons
!    set one and set two respectively. It has to check both ways
! 3) The direction of the vector is found by setting the vector to be towards the middle of the polygon if the set one polygon points are outside
!    the set two polygons, and towards the outside of the polygon if the set two polygon points are outside of the set one polygons



! See the run script for a quick way to generate this file from existing multi-segment files

	use pip
	use direction_mod
	use read_polygons

	implicit none

	integer :: counter, istat, polygon_counter, points_counter, opposite_step, polygon_counter2
	integer :: previous_point, next_point, direction, record_counter
	integer, parameter :: dir_unit = 30, minmax_unit = 40, dir_text_unit = 50



	double precision :: slope, mid_x, mid_y, current_x, current_y, intercept, temp_x, temp_y, angle,  last_x, last_y
	double precision :: next_x, next_y
	double precision :: min_x, min_y, max_x, max_y, x_component, y_component, direction1, direction2
	double precision, parameter :: arbitrary_increment = 0.01




	character (len=256), parameter :: dir_file = "direction_file.bin", minmax_file = "minmax.txt"
	character (len=256), parameter :: dir_text_file = "direction_file.txt"

	logical :: true_inside, is_inside
	logical, dimension(:,:,:), allocatable :: inside 

	! time values

      integer,dimension(8) :: values


! 1)  read in the polygons

	call read_polygons_init()


	allocate(inside(2,max_polygons,max_points), stat=istat)
	if(istat /=0) THEN
		write(6,*) "problem allocating arrays in find_direction"
		stop
	endif


	call date_and_time(VALUES=values)
	write(6,*) "finished reading in polygon: ", values(5), values(6), values(7)
	
	call date_and_time(VALUES=values)
	write(6,'(A30,I2,A1,I2,A1,I2,A1)') "finished reading in polygon: ", values(5), ":", values(6), ":", values(7)


	
! 2) find out if the points of the polygon are inside another polygon

	inside = .false.

	find_inside: do counter = 1, 2, 1

		if(counter == 1) then
			opposite_step = 2
		else
			opposite_step = 1
		endif

		polygon_loop: do polygon_counter = 1, number_polygons(counter)

			points_loop: do points_counter = 1, polygon_points(counter,polygon_counter), 1

				polygon_loop2: do polygon_counter2 = 1, number_polygons(opposite_step)


					if( x_coordinates(counter,polygon_counter,points_counter) >= &
							poly_x_min(opposite_step,polygon_counter2) .and. &
					    x_coordinates(counter,polygon_counter,points_counter) <= &
							poly_x_max(opposite_step,polygon_counter2) .and. &
					    y_coordinates(counter,polygon_counter,points_counter) <= &
							poly_y_max(opposite_step,polygon_counter2) .and. &
					    y_coordinates(counter,polygon_counter,points_counter) >= &
							poly_y_min(opposite_step,polygon_counter2)  ) THEN

						! checks to see if this point is inside of another polygon, if so, it is flagged
						inside(counter,polygon_counter,points_counter) = point_in_polygon( &
						  x_coordinates(opposite_step,polygon_counter2,1:polygon_points(opposite_step,polygon_counter2)), &
						  y_coordinates(opposite_step,polygon_counter2,1:polygon_points(opposite_step,polygon_counter2)), &
						  x_coordinates(counter,polygon_counter,points_counter), &
						  y_coordinates(counter,polygon_counter,points_counter), &
						  polygon_points(opposite_step,polygon_counter2))
					else
						inside(counter,polygon_counter,points_counter) = .false.
					endif



					if(inside(counter,polygon_counter,points_counter)) THEN
						cycle points_loop
					endif


				end do polygon_loop2
			end do points_loop
		end do polygon_loop
	end do find_inside


	call date_and_time(VALUES=values)
	write(6,'(A27,I2,A1,I2,A1,I2,A1)') "finished point in polygon: ", values(5), ":", values(6), ":", values(7)
	

! 3) find the direction vectors

	! might as well find the range of the resulting polygons while I am at it
	min_x = x_coordinates(1,1,1)
	min_y = y_coordinates(1,1,1)
	max_x = x_coordinates(1,1,1)
	max_y = y_coordinates(1,1,1)

	! save into a binary file that can be read into GMT
	open(unit=dir_unit, file=dir_file, access="direct", form="unformatted", recl=8, status="replace")
	open(unit=dir_text_unit, file=dir_text_file, access="sequential", form="formatted", status="replace")

	record_counter = 0

	find_direction_loop: do counter = 1, 2, 1

		if (counter == 1) THEN
			true_inside = .false.
		else
			true_inside = .true.
		endif

		do polygon_counter = 1, number_polygons(counter), 1

			do points_counter = 1, polygon_points(counter,polygon_counter), 1

				previous_point = points_counter - 1
				next_point = points_counter + 1

				if(points_counter == 1) THEN
					previous_point = polygon_points(counter,polygon_counter)
				endif

				if(points_counter == polygon_points(counter,polygon_counter) + 1 ) THEN
					next_point = 1
				endif

				! this is taken from ICESHEET, read_icefile.f90

				current_x = x_coordinates(counter,polygon_counter,points_counter)
				current_y = y_coordinates(counter,polygon_counter,points_counter)



				if(current_x < min_x) THEN
					min_x = current_x
				endif
				if(current_y < min_y) THEN
					min_y = current_y
				endif
				if(current_x > max_x) THEN
					max_x = current_x
				endif
				if(current_y > max_y) THEN
					max_y = current_y
				endif
			


				last_x = x_coordinates(counter,polygon_counter,previous_point)
				last_y = y_coordinates(counter,polygon_counter,previous_point)

				next_x = x_coordinates(counter,polygon_counter,next_point)
				next_y = y_coordinates(counter,polygon_counter,next_point)


				direction1 = atan2(last_y-current_y,last_x-current_x)
				direction2 = atan2(next_y-current_y, next_x-current_x)

				angle = average_direction(direction1, direction2)

				temp_x = current_x + cos(angle) * arbitrary_increment
				temp_y = current_y + sin(angle) * arbitrary_increment



				! note that this will fail miserably if the points are closer together than the arbitrary increment
				! I could add a check, but that would increase the overhead a lot, and really the user
				! should split the polygon in two if the polygon closes together that closely
				
				is_inside = point_in_polygon( &
					 x_coordinates(counter,polygon_counter,1:polygon_points(counter,polygon_counter)), &
					 y_coordinates(counter,polygon_counter,1:polygon_points(counter,polygon_counter)),  &
					 temp_x, temp_y, polygon_points(counter,polygon_counter))




				if (.not. is_inside) THEN

					angle = check_angle(angle + pi)

				endif



				if(inside(counter,polygon_counter,points_counter) .neqv. true_inside) THEN ! the vector should be pointing in the opposite direction

					angle = check_angle(angle + pi)
!
				endif


				! create a unit vector

	
				x_component =  cos(angle)
				y_component =  sin(angle)



				record_counter=record_counter+1
				write(dir_unit,rec=record_counter) current_x
				record_counter=record_counter+1
				write(dir_unit,rec=record_counter) current_y
				record_counter=record_counter+1
				write(dir_unit,rec=record_counter) atan2(y_component,x_component)


				write(dir_text_unit,*) current_x, current_y, atan2(y_component,x_component) * 180. / pi, 0.5 ! 0.5 is the length of vector for plotting in GMT

			end do


		end do


	end do find_direction_loop


	close (unit=dir_unit)
	close (unit=dir_text_unit)

	open(unit=minmax_unit, file=minmax_file, access="sequential", form="formatted", status="replace")


	write(minmax_unit,*) min_x, max_x
	write(minmax_unit,*) min_y, max_y

	close(unit=minmax_unit)

	call read_polygons_clear()

	deallocate(inside, stat=istat)
	if(istat /=0) THEN
		write(6,*) "problem deallocating arrays in find_direction"
		stop
	endif

	call date_and_time(VALUES=values)
	write(6,'(A27,I2,A1,I2,A1,I2,A1)') "finished program: ", values(5), ":", values(6), ":", values(7)



end program find_direction
