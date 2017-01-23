!	copyright Evan. J. Gowan, 2013, 2016
!
!	This file is part of countour_interpolation
!
!	countour_interpolation is free software: you can redistribute it and/or modify
!	it under the terms of the GNU General Public License as published by
!	the Free Software Foundation, version 3 of the License
!
!	countour_interpolation is distributed in the hope that it will be useful,
!	but WITHOUT ANY WARRANTY; without even the implied warranty of
!	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!	GNU General Public License for more details.
!
!	You should have received a copy of the GNU General Public License
!	along with countour_interpolation.  If not, see <http://www.gnu.org/licenses/>.

!***********************************************************************************


module pip

contains

logical function point_in_polygon(x_boundary, y_boundary, x, y, number_points)

! this function determines whether a give point (x,y) is within a polygon defined by x_boundary and y_boundary
	
	implicit none

	integer, intent(in) :: number_points
	double precision, intent(in) :: x, y
	double precision, dimension(number_points), intent(in) :: x_boundary, y_boundary

	integer :: current_point, next_point, last_point, crossover_counter
	logical :: found_first, found_last, inside

	found_first = .false.
	found_last = .false.
	inside = .false.

	current_point = 1
	search_boundary: do

		next_point = current_point + 1
	
		if (next_point > number_points) THEN
			next_point = 1
			found_last = .true.
		endif

! even-odd rule algorithm to determine if the point is inside or outside

		if (min(y_boundary(current_point), y_boundary(next_point)) < y .and.&
		    max(y_boundary(current_point), y_boundary(next_point)) >= y) THEN

			if (x_boundary(current_point) + (y - y_boundary(current_point)) /&
			    (y_boundary(next_point) - y_boundary(current_point)) * &
			    (x_boundary(next_point) - x_boundary(current_point)) < x) THEN

				inside = .not.(inside)

			endif

		endif

		current_point = current_point + 1

		if (found_last) THEN
			exit search_boundary
		endif
		
	
	end do search_boundary

	point_in_polygon = inside

	return
end function point_in_polygon




end module pip
