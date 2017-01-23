module read_minmax

	implicit none

	integer :: number_x_grid, number_y_grid

	double precision :: min_x_grid, max_x_grid, min_y_grid, max_y_grid


! used only after running find_direction

contains

subroutine run_read_minmax(grid_spacing)

	implicit none

	double precision, intent(in) :: grid_spacing

	integer :: istat
	integer, parameter :: minmax_unit = 40

	double precision :: min_x, max_x, min_y, max_y

	character (len=256), parameter :: minmax_file = "minmax.txt"


	! extreme points in the contours of interest, file created when find_direction is run
	open(unit=minmax_unit, file=minmax_file, access="sequential", form="formatted", status="old", iostat=istat)
	if(istat /=0) THEN
		write(6,*) "error opening ", trim(adjustl(minmax_file))
		write(6,*) "did you run 'find_direction' first?"
		stop
	endif

	read(minmax_unit,*) min_x, max_x
	read(minmax_unit,*) min_y, max_y

	close(unit=minmax_unit)

	close (unit=11)


	! grid creation step
	! buffer by one grid unit
	min_x_grid = dble(floor(min_x/grid_spacing)) * grid_spacing - grid_spacing
	min_y_grid = dble(floor(min_y/grid_spacing)) * grid_spacing - grid_spacing
	max_x_grid = dble(ceiling(max_x/grid_spacing)) * grid_spacing + grid_spacing
	max_y_grid = dble(ceiling(max_y/grid_spacing)) * grid_spacing + grid_spacing


	number_x_grid = nint((max_x_grid-min_x_grid)/grid_spacing)
	number_y_grid = nint((max_y_grid-min_y_grid)/grid_spacing)
end subroutine run_read_minmax


subroutine find_grid_index(x, y, grid_spacing, x_grid_index, y_grid_index)
	implicit none

	double precision, intent(in) :: x, y, grid_spacing
	integer, intent(out) :: x_grid_index, y_grid_index

	x_grid_index = nint((x-min_x_grid) / grid_spacing)+1
	y_grid_index = nint((y-min_y_grid) / grid_spacing)+1

	if(x_grid_index < 1 .or. x_grid_index > number_x_grid) then
		write(6,*) "warning: x-coordinate out of bounds"
		stop
	endif

	if(y_grid_index < 1 .or. y_grid_index > number_y_grid) then
		write(6,*) "warning: y-coordinate out of bounds"
		stop
	endif

end subroutine find_grid_index


subroutine find_grid_location(x_grid_index, y_grid_index, grid_spacing, x, y)
! not done

	double precision, intent(in) :: grid_spacing
	integer, intent(in) :: x_grid_index, y_grid_index
	double precision, intent(out) :: x, y

	x = min_x_grid + dble(x_grid_index-1) * grid_spacing
	y = min_y_grid + dble(y_grid_index-1) * grid_spacing

end subroutine find_grid_location

end module read_minmax
