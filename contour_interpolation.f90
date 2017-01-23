program contour_interpolation


	double precision, dimension(:), allocatable :: contour1_x, contour2_x, contour1_y, contour2_y

	integer :: countour1_points, contour2_points

	integer, parameter :: contour_file_unit = 10, param_unit=20

	character(len=256) :: countour1_file, contour2_file

! read in parameters

	open(unit=param_unit, file="params.txt", access="sequential", form="formatted", status="old")

	read(param_unit,*) contour1_file
	read(param_unit,*) contour2_file


! first step, read in the contour lines



end program contour_interpolation
