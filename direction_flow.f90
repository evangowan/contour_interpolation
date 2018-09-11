program direction_flow

	IMPLICIT NONE

	integer :: number_of_records, counter, rec_num, commandline_count, line_counter, total_lines, point_counter, istat
	integer, parameter :: dir_unit = 30, supp_unit = 40
	double precision :: x, y, dir, x_diff, y_diff

	integer, parameter :: max_lines = 10000, max_points = 10000
	integer, dimension(max_lines) :: line_points

	double precision, dimension(max_lines,max_points) :: x_coordinates, y_coordinates

	character(len=256) :: arg_in, supplement_gmt_file

	character (len=256), parameter :: dir_file = "direction_file.bin"


	character (len=1) :: divider ! make sure that there are no spaces before the ">" character or this will mess up
	character (len=1), parameter :: divider_character = ">", ignore_character = "#"

	commandline_count = command_argument_count()

	if(commandline_count /= 1) THEN
		write(6,*) "not enough args: ", commandline_count
		write(6,*) "correct usage of this program:"
		write(6,*) "direction_flow [file] "
		stop
	endif

	call GET_COMMAND_ARGUMENT(1,arg_in) 
	read(arg_in,*) supplement_gmt_file


	! open the gmt formatted file

	open(unit=supp_unit, file=supplement_gmt_file, access="sequential", form="formatted", status="old")


	open(unit=dir_unit, file=dir_file, access="direct", form="unformatted", recl=8, status="old")

	inquire(dir_unit, size=number_of_records)
	number_of_records = number_of_records / 8

!	write(6,*) number_of_records

	! read in the supplementary direction file, assuming that the lines point in the direction of retreat/advance (depending on application)
	line_counter = 0
	supp_loop: do

		read(supp_unit,*, iostat=istat) divider
		if(istat /=0) THEN
			exit supp_loop
		endif

		if(divider == divider_character) THEN
			line_counter = line_counter + 1
			line_points(line_counter) = 0
			cycle supp_loop

		elseif (divider == ignore_character) THEN ! anything starting with '#' should be ignored
				
			cycle supp_loop
		else
			backspace(unit=supp_unit)

			line_points(line_counter) = line_points(line_counter) + 1

			read(supp_unit,*) x, y

			x_coordinates(line_counter,line_points(line_counter)) = x
			y_coordinates(line_counter,line_points(line_counter)) = y

		endif


	end do supp_loop

	close(supp_unit)

	total_lines = line_counter

	! add the supplementary directions to the dir_file

	do line_counter = 1, total_lines, 1
!		write(6,*) line_counter,  line_points(line_counter)
		do point_counter = 1, line_points(line_counter) - 1, 1

			x = (x_coordinates(line_counter,point_counter) + &
			     x_coordinates(line_counter,point_counter+1) ) / 2.0
			y = (y_coordinates(line_counter,point_counter) + &
			     y_coordinates(line_counter,point_counter+1) ) / 2.0

			x_diff = x_coordinates(line_counter,point_counter+1) - &
				   x_coordinates(line_counter,point_counter)
			y_diff = y_coordinates(line_counter,point_counter+1) - &
				   y_coordinates(line_counter,point_counter)

			dir = atan2( y_diff, x_diff)


			number_of_records=number_of_records+1
			write(dir_unit,rec=number_of_records) x
			number_of_records=number_of_records+1
			write(dir_unit,rec=number_of_records) y
			number_of_records=number_of_records+1
			write(dir_unit,rec=number_of_records) dir


		end do
	end do

!	do counter = 1, number_of_records / 3
!		rec_num = (counter-1) * 3 + 1
!		read(dir_unit,rec=rec_num) x
!		read(dir_unit,rec=rec_num+1) y
!		read(dir_unit,rec=rec_num+2) dir

!		write(666,*) x, y, dir
!	end do


	close(unit=dir_unit)

end program direction_flow
