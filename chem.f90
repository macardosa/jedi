module chem
	implicit none
	integer, parameter :: dp=kind(0.d0)
	real(dp) :: PI = 4*ATAN(1.d0)
	
	private
	public altitude, projection_onto_plane, normal_to_plane, discriminant, cross_product, distance, compute_length, & 
	compute_angle, compute_dihedral, compute_improper
	
	contains
		function altitude(A, B, C)
			! computes altitude of triangle abc with respect to base bc
			real(dp), dimension(:), intent(in) :: A, B, C
			real(dp), dimension(3) :: altitude, v, M
			real(dp) :: ba, altitude_norm, bm

			ba = NORM2(B-A)
			v = C - B
			altitude_norm = NORM2(cross_product( (A-C), v))/NORM2(v)

			bm = SQRT( ba**2 - altitude_norm**2 )
			M = B + (bm/NORM2(v))*v

			altitude = A - M
		end function

		function projection_onto_plane(v, n)
			real(dp), dimension(:), intent(in) :: v, n
			real(dp), dimension(3) :: pvector
			real(dp) :: projection_onto_plane

			pvector = v - ( DOT_PRODUCT(v,n)/NORM2(n)**2 )*n
			projection_onto_plane = NORM2(pvector)**2
		end function

		function normal_to_plane(a, b, c)
			real(dp), dimension(:), intent(in) :: a, b, c
			real(dp), dimension(3) :: normal_to_plane, u, v

			u = a - b
			v = c - b

			normal_to_plane = cross_product(u,v)
			normal_to_plane = normal_to_plane/NORM2(normal_to_plane)

		end function

		function discriminant(a, b, c, d)
			real(dp) :: a, b, c, d
			real(dp) :: discriminant

			discriminant = 18*a*b*c*d - 4*(b**3)*d + (b**2)*(c**2) - 4*a*(c**3) - 27*(a**2)*(d**2)
		end function

		function cross_product(u, v)
			real(dp), dimension(3) :: u, v, cross_product

			cross_product(1) = u(2)*v(3) - u(3)*v(2)
			cross_product(2) = -u(1)*v(3) + u(3)*v(1)
			cross_product(3) = u(1)*v(2) - u(2)*v(1)

		end function

		function distance(b, c)
			real(dp) :: distance
			real(dp), dimension(3), intent(in) :: b, c
			integer :: i

			distance = 0
			do i=1,3
				distance = distance + (c(i)-b(i))**2
			enddo
			!distance = sqrt(distance)
		end function
		
		function compute_length(u, v)
			real(dp), dimension(:) :: u, v
			real(dp) :: compute_length

			compute_length = NORM2( v - u )
		end function

		function compute_angle(u, v, w)
			real(dp), dimension(:) :: u, v, w
			real(dp) :: compute_angle

			compute_angle = ACOS( DOT_PRODUCT((u - v), (w - v))/(NORM2(u-v)*NORM2(w-v)) )
			
		end function

		function compute_dihedral(u, v, w, x)
			real(dp), dimension(:) :: u, v, w, x
			real(dp), dimension(3) :: n1, n2, m
			real(dp) :: compute_dihedral, s, t

			n1 = cross_product( (u-v), (w-v) )
			n1 = n1/NORM2(n1)
			n2 = cross_product( (v-w), (x-w) )
			n2 = n2/NORM2(n2)

			m = (w-v)/NORM2(w-v)
			m = cross_product(n1,m)

			s = DOT_PRODUCT(n1,n2)
			t = DOT_PRODUCT(m,n2)

			compute_dihedral = -ATAN2(t,s)
			
		end function

		function compute_improper(u,v,w,x)
			real(dp), dimension(:) :: u, v, w, x
			real(dp), dimension(3) :: n1, n2, m
			real(dp) :: compute_improper, s, t

			n1 = cross_product( (v-u), (w-u) )
			n1 = n1/NORM2(n1)

			n2 = cross_product( (w-v), (x-v) )
			n2 = n2/NORM2(n2)

			compute_improper = ACOS( DOT_PRODUCT(n1,n2) )
		end function
		
end module chem
