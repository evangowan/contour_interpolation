! code written by Evan Gowan (evangowan@gmail.com)
module spline
! cubic_spline storage
	integer, parameter :: spline_points = 20
	double precision, dimension(spline_points), save :: x_store, y_store
	double precision, dimension(2,spline_points), save :: parametric_parameters
contains

subroutine cubic_spline(x_points,y_points,t,x,y,slope)
! this function returns the local slope. Input includes the "spline_points" points needed to do the 
! interpolation (the interval of interest is between the n/2 and n/2+1 points).
! t is the distance between those points in the parametric equations, and is an input as well, must be in the range [0,1]
! (i.e., t=0 is the point n/2, and t=1 is the point n/2+1). x, y, and slope are the output.
! note that the values of x_points and y_points are stored in this module so that repeat calculations of
! the spline parameters are unnecessary for the same interval.
! see this website: 
! Weisstein, Eric W. "Cubic Spline." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/CubicSpline.html 
	double precision, dimension(spline_points), intent(in) :: x_points, y_points
	double precision, intent(in) :: t
	double precision, intent(out) :: x, y, slope
	double precision, dimension(2,spline_points) :: D_vect
	double precision, dimension(2,spline_points) :: point_equations
	double precision, dimension(2,spline_points) :: matrix_values
	double precision, dimension(2,spline_points) :: points
	double precision, dimension(2,spline_points) :: right_solved
	double precision, dimension(4) :: results
	integer :: counter, counter2
	
	
	if(sum(x_points - x_store) /= 0 .or. sum(y_points - y_store) /= 0) THEN ! calculate the parameters
		! store the points in a 2xspline_points array so that the code doesn't have to be repeated for both x and y
		points(1,:) = x_points
		points(2,:) = y_points
		do counter = 1, 2, 1
		! first calculate the point equations (right hand side of equation 18 on the website)

			point_equations(counter,1) = 3. * (points(counter,2) - points(counter,1))
			point_equations(counter,spline_points) = 3. * (points(counter,spline_points) - points(counter,spline_points-1))
			do counter2 = 2, spline_points-1, 1
				point_equations(counter,counter2) = 3. * (points(counter,counter2+1) - points(counter,counter2-1))
			end do
		! next step is to calculate the values of the diagonal matrix found in a RREF version of the matrix on the left side in equation 18
		! after solving the equations by hand, it was a pretty simple pattern.
			matrix_values(counter,1) = 2.
			matrix_values(counter,2) = 4. * matrix_values(counter,1) - 1
			do counter2 = 3, spline_points-1, 1
				matrix_values(counter,counter2) = 4. * matrix_values(counter,counter2-1) - matrix_values(counter,counter2-2)
			end do
			matrix_values(counter,spline_points) = 2. * matrix_values(counter,counter2-1) - matrix_values(counter,counter2-2)
		! after converting the above matrix to RREF, the right hand side of equation 18 also changes
			right_solved(counter,1) = point_equations(counter,1)
			do counter2 = 1, spline_points, 1
				right_solved(counter,counter2) = matrix_values(counter,counter2-1) * point_equations(counter,counter2) - &
					right_solved(counter,counter2-1)
			end do

		! calculate the "D" coefficients used to solve the parametric equations (the left hand side of equation 18 on the website)
			D_vect(counter,spline_points) = right_solved(counter,spline_points) / matrix_values(counter,spline_points)
			do counter2 = spline_points-1, 2, -1
				D_vect(counter,counter2) = (right_solved(counter,counter2) - &
					matrix_values(counter,counter2-1)*D_vect(counter,counter2+1)) / matrix_values(counter,counter2)
			end do
			D_vect(counter,1) = (right_solved(counter,1) - D_vect(counter,2)) / matrix_values(counter,1)
			do counter2 = 1, spline_points
			end do

		! now determine the parametric equation parameters, which will be stored for future use (these are equations 6-9 on the website)

			parametric_parameters(counter,1) = points(counter,spline_points/2)
			parametric_parameters(counter,2) = D_vect(counter,spline_points/2)
			parametric_parameters(counter,3) = 3. * (points(counter,spline_points/2+1)-points(counter,spline_points/2)) - 2. &
				* D_vect(counter,spline_points/2) - D_vect(counter,spline_points/2+1) 
			parametric_parameters(counter,4) = 2.*(points(counter,spline_points/2) - points(counter,spline_points/2+1)) + &
				D_vect(counter,spline_points/2) + D_vect(counter,spline_points/2+1)
		end do
	endif
	do counter = 1, 2
		! solve the parametric equations (equation 1 on the website)
		results(counter) = parametric_parameters(counter,1) + parametric_parameters(counter,2) * t + &
			parametric_parameters(counter,3) * t**2 + parametric_parameters(counter,4) * t**3
		! this is the derivative of the parametric equation, which can be used to find the local slope
		results(counter+2) = parametric_parameters(counter,2) + &
			2.*parametric_parameters(counter,3) * t + 3.*parametric_parameters(counter,4) * t**2
	end do
	x = results(1)
	y = results(2)
	slope = results(4) / results(3)


end subroutine cubic_spline

end module spline
