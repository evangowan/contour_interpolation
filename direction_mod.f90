module direction_mod

! taken from global_parameters.f90 in ICESHEET

	implicit none

	double precision, parameter :: pi = 3.141592653589793


contains


double precision function average_direction(direction1, direction2)


	implicit none

	double precision, intent(in) :: direction1, direction2


	if((direction1 > pi/2.0 .and. direction2 <  -pi/2.0) ) THEN
		average_direction = check_angle( (direction1 + (2.0*pi+direction2))/2.0)

	elseif (direction1 < -pi/2.0 .and. direction2 > pi/2.0) THEN

		average_direction = check_angle( (direction2 + (2.0*pi+direction1))/2.0)

	elseif (direction1 > pi/2.0 .and. (direction2 >-pi/2.0 .and. direction2 < 0.0 )) THEN
		
		if (direction1-direction2 <= pi) THEN

			average_direction = (direction1 + direction2) / 2.0

		else

			average_direction = check_angle( (direction1 + (2.0*pi+direction2))/2.0)


		endif

	elseif (direction2 > pi/2.0 .and. (direction1 >-pi/2.0 .and. direction1 < 0.0 )) THEN

		if (direction2-direction1 <= pi) THEN

			average_direction = (direction1 + direction2) / 2.0

		else

			average_direction = check_angle( (direction2 + (2.0*pi+direction1))/2.0)

		endif

	else 
		average_direction = (direction1 + direction2) / 2.0


	endif

	

end function average_direction


double precision function check_angle(angle)

	implicit none
	double precision, intent(in) :: angle
!	double precision ::  check_angle

	if (angle > pi) THEN
		check_angle = angle - 2.0 * pi
	elseif (angle < -pi) THEN
		check_angle = angle + 2.0 * pi
	else
		check_angle = angle
	end if


end function check_angle


double precision function general_direction(direction1, direction2, rotate_factor)


	implicit none

	double precision, intent(in) :: direction1, direction2, rotate_factor


	if((direction1 >= pi/2.0 .and. direction2 <=  -pi/2.0) ) THEN
		general_direction = check_angle( ((2*pi+direction2)-direction1)*rotate_factor + direction1)

	elseif (direction1 <= -pi/2.0 .and. direction2 >= pi/2.0) THEN

		general_direction = check_angle( direction1 - ((2.0*pi+direction1)-direction2)*rotate_factor )

	elseif (direction1 >= pi/2.0 .and. (direction2 >=-pi/2.0 .and. direction2 < 0.0 )) THEN
		
		if (direction1-direction2 <= pi) THEN

			general_direction = direction1 - (direction1 - direction2) *rotate_factor

		else

			general_direction = check_angle( direction1+(2.0*pi -(direction1-direction2))*rotate_factor)


		endif

	elseif (direction2 >= pi/2.0 .and. (direction1 >=-pi/2.0 .and. direction1 <= 0.0 )) THEN

		if (direction2-direction1 <= pi) THEN

			general_direction = direction1 + (direction2 - direction1)*rotate_factor

		else

			general_direction = check_angle(direction1 - (2.0*pi -(direction2 - direction1))*rotate_factor)

		endif

	else 
		general_direction = check_angle(direction1 - (direction1-direction2)*rotate_factor)

	endif

	

end function general_direction

end module direction_mod
