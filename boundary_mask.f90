program boundary_mask

! reads in the boundary polygons and creates two masks 

! must run find_direction before running this

	use read_polygons
	use read_minmax
	use boundary_mask_mod
	use pip

	integer :: counter, polygon_counter, istat, points_counter, next_point, current_index, commandline_count
	integer :: min_x_index, min_y_index, max_x_index, max_y_index, x_counter, y_counter
	integer, parameter ::   gmt_unit=80, gmt_unit2=90




	double precision :: grid_spacing, current_x, current_y, next_x, next_y, slope, intercept
	double precision :: min_x_segment, min_y_segment, max_x_segment, max_y_segment, cell_x, cell_y
	double precision :: min_x, max_x, min_y, max_y

	character(len=256) :: arg_in

	character (len=256), dimension(2) :: gmt_files, gmt_files2

	logical :: inside, was_inside

	commandline_count = command_argument_count()

	if(commandline_count /= 1) THEN
		write(6,*) "correct usage of this program:"
		write(6,*) "boundary_mask grid_spacing"
		stop
	endif

	call GET_COMMAND_ARGUMENT(1,arg_in) 
	read(arg_in,*) grid_spacing


	! for plotting, testing
	gmt_files(1) = "boundary_mask_1.txt"
	gmt_files(2) = "boundary_mask_2.txt"

	gmt_files2(1) = "inside_mask_1.txt"
	gmt_files2(2) = "inside_mask_2.txt"

	! read in the polygons
	call read_polygons_init()
	call run_read_minmax(grid_spacing)

	! allocate memory
	call allocate_mask(number_x_grid,number_y_grid)



	mask = 0 ! mask is true when it has the value of 1

	do counter = 1, 2, 1

		do polygon_counter = 1, number_polygons(counter), 1


			do points_counter = 1, polygon_points(counter,polygon_counter), 1


				

				if(points_counter == polygon_points(counter,polygon_counter) + 1 ) THEN
					next_point = 1
				else

					next_point = points_counter + 1
				endif


				current_x = x_coordinates(counter,polygon_counter,points_counter)
				current_y = y_coordinates(counter,polygon_counter,points_counter)

				next_x = x_coordinates(counter,polygon_counter,next_point)
				next_y = y_coordinates(counter,polygon_counter,next_point)

				! flare out a bit to make sure each cell is checked
				min_x_segment = min(current_x,next_x) - grid_spacing/2.
				max_x_segment = max(current_x,next_x) + grid_spacing/2.
				min_y_segment = min(current_y,next_y) - grid_spacing/2.
				max_y_segment = max(current_y,next_y) + grid_spacing/2.


				min_x_index = max(1, nint((min_x_segment-min_x_grid) / grid_spacing)+1)
				min_y_index = max(1, nint((min_y_segment-min_y_grid) / grid_spacing)+1)
				max_x_index = min(number_x_grid, nint((max_x_segment-min_x_grid) / grid_spacing)+1)
				max_y_index = min(number_y_grid, nint((max_y_segment-min_y_grid) / grid_spacing)+1)


				if(min_x_segment == max_x_segment) THEN ! just go between the y points

					do y_counter = min_y_index, max_y_index

						current_index = nint((min_x_segment-min_x_grid) / grid_spacing)+1
						mask(counter,current_index, y_counter) = 1


					end do
				else if (min_y_segment == max_y_segment) THEn ! just go between x points

					do x_counter = min_y_index, max_y_index

						current_index = nint((max_y_segment-min_y_grid) / grid_spacing)+1
						mask(counter,x_counter, current_index) = 1

					end do
				else ! go along the line in both the x and y directions

					do x_counter = min_x_index, max_x_index
						do y_counter = min_y_index, max_y_index

							cell_x = dble(x_counter-1) * grid_spacing + min_x_grid
							cell_y = dble(y_counter-1) * grid_spacing + min_y_grid
							inside = inside_cell(current_x, next_x, current_y, next_y, cell_x, cell_y, grid_spacing)


							if(inside) then
								mask(counter,x_counter,y_counter) = 1
							endif
						end do
					end do

				endif
				

			end do
		end do

	end do

	! add a mask for the interior of the polygons

	do counter = 1, 2, 1

		do polygon_counter = 1, number_polygons(counter), 1
			min_x = minval( x_coordinates(counter,polygon_counter,1:polygon_points(counter,polygon_counter)))
			max_x = maxval( x_coordinates(counter,polygon_counter,1:polygon_points(counter,polygon_counter)))
			min_y = minval( y_coordinates(counter,polygon_counter,1:polygon_points(counter,polygon_counter)))
			max_y = maxval( y_coordinates(counter,polygon_counter,1:polygon_points(counter,polygon_counter)))

			do x_counter = 1, number_x_grid

				do y_counter = 1, number_y_grid

					if(mask(counter,x_counter,y_counter) == 0) THEN 

						cell_x = dble(x_counter-1) * grid_spacing + min_x_grid
						cell_y = dble(y_counter-1) * grid_spacing + min_y_grid

						if(cell_x > min_x .and. cell_x < max_x .and. cell_y > min_y .and. cell_y < max_y) THEN ! should check

							inside = point_in_polygon(&
							  x_coordinates(counter,polygon_counter,1:polygon_points(counter,polygon_counter)),&
							  y_coordinates(counter,polygon_counter,1:polygon_points(counter,polygon_counter)),&
							  cell_x, cell_y, polygon_points(counter,polygon_counter))

							if(inside) THEN
								mask(counter,x_counter,y_counter) = 2
							endif

						endif



					endif


				end do
			end do
		end do
	end do

	! write out the files

	call write_mask(number_x_grid,number_y_grid)

	! for plotting/debugging, write out gmt file
	do counter = 1, 2, 1


		open(unit=gmt_unit, file=gmt_files(counter), access="sequential", form="formatted",  status="replace")
		open(unit=gmt_unit2, file=gmt_files2(counter), access="sequential", form="formatted",  status="replace")

		do x_counter = 1, number_x_grid
			do y_counter = 1, number_y_grid


				if(mask(counter,x_counter, y_counter) == 1) THEN
					write(gmt_unit, *) dble(x_counter-1) * grid_spacing + min_x_grid, &
							 	 dble(y_counter-1) * grid_spacing + min_y_grid
				endif
				if(mask(counter,x_counter, y_counter) == 2) THEN
					write(gmt_unit2, *) dble(x_counter-1) * grid_spacing + min_x_grid, &
							 	 dble(y_counter-1) * grid_spacing + min_y_grid
				endif

			end do
		end do
		close(unit=gmt_unit)
		close(unit=gmt_unit2)

	end do

	! clear memory
	call read_polygons_clear()
	call boundary_mask_clear()



contains

	logical function inside_cell(x1, x2, y1, y2, cell_x, cell_y, grid_spacing)
		! I didn't think this would be necessary to be this rigorous and check each cell, but it was missing
		! some points by simply going back and forth along the line at grid spacing intervals

		implicit none

		double precision, intent(in) :: x1, x2, y1, y2, cell_x, cell_y, grid_spacing
		
		double precision :: min_cell_x, min_cell_y, max_cell_x, max_cell_y, slope, intercept, crossover_point

		inside_cell = .false.

		min_cell_x = cell_x - grid_spacing / 2.
		min_cell_y = cell_y - grid_spacing / 2.
		max_cell_x = cell_x + grid_spacing / 2.
		max_cell_y = cell_y + grid_spacing / 2.


		if(x1 < min_cell_x .and. x2 < min_cell_x) THEN
			return
		endif

		if(y1 < min_cell_y .and. y2 < min_cell_y) THEN
			return
		endif

		if(x1 > max_cell_x .and. x2 > max_cell_x) THEN
			return
		endif

		if(y1 > max_cell_y .and. y2 > max_cell_y) THEN
			return
		endif
		
		slope = (y2-y1) / (x2-x1)
		intercept = y1 - slope * x1

		! min_cell_x
		crossover_point = min_cell_x * slope + intercept

		if (crossover_point >= min_cell_y .and. crossover_point <= max_cell_y) THEN
			inside_cell = .true.
			return
		endif

		! max_cell_x
		crossover_point = max_cell_x * slope + intercept

		if (crossover_point >= min_cell_y .and. crossover_point <= max_cell_y) THEN
			inside_cell = .true.
			return
		endif


	

		slope = (x2-x1) / (y2-y1)
		intercept = x1 - slope * y1



		! min_cell_x
		crossover_point = min_cell_y * slope + intercept

		if (crossover_point >= min_cell_x .and. crossover_point <= max_cell_x) THEN
			inside_cell = .true.
			return
		endif

		! max_cell_x
		crossover_point = max_cell_y * slope + intercept

		if (crossover_point >= min_cell_x .and. crossover_point <= max_cell_x) THEN
			inside_cell = .true.
			return
		endif


	end function inside_cell

end program boundary_mask
