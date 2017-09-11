program eliminate_outside


	use pip
	use read_polygons

	implicit none

	integer :: counter, score, poly_counter, istat
	integer, parameter :: contour_unit=25, out_unit=35

	double precision :: x, y, z

	character (len=1) :: divider ! make sure that there are no spaces before the ">" character or this will mess up
	character(len=80) :: all_line
	character (len=1), parameter :: divider_character = ">", ignore_character = "#"
	character(len=256), parameter :: contour_file="contours.txt", out_contour_file="contours_elim.txt"

	logical :: inside

! 1)  read in the polygons

	call read_polygons_init()


	open(unit=contour_unit, file = contour_file, access="sequential", form="formatted", status="old")
	open(unit=out_unit, file =out_contour_file, access="sequential", form="formatted", status="replace")

	write(out_unit,'(A)') ">"

	read_file: do

		read(contour_unit,*, iostat=istat) divider
		if(istat /=0) THEN
			exit read_file
		endif

		if(divider == divider_character) THEN
			backspace(unit=contour_unit)
			read(contour_unit,'(A)') all_line
			write(out_unit,'(A)') adjustl(trim(all_line))
			cycle read_file
		elseif (divider == ignore_character) THEN ! anything starting with '#' should be ignored
			cycle read_file
		else
			backspace(unit=contour_unit)
	
			read(contour_unit,*) x, y, z
		endif

		score = 0
		find_inside: do counter = 1, 2, 1

			do poly_counter = 1, number_polygons(counter), 1

				inside = point_in_polygon( &
					  x_coordinates(counter,poly_counter,1:polygon_points(counter,poly_counter)), &
					  y_coordinates(counter,poly_counter,1:polygon_points(counter,poly_counter)), &
					  x, y, polygon_points(counter,poly_counter))

				if(inside) then
					score = score + 1
					cycle find_inside
				endif

			end do


		end do find_inside

		if(score == 1) THEN
			write(out_unit,*) x, y, z
		else
			write(out_unit,'(A)') adjustl(trim(all_line))
		endif

	end do read_file

	call read_polygons_clear()
	close(unit=contour_unit)
	close(unit=out_unit)

end program eliminate_outside
