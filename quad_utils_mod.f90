module quad_utils_mod
  implicit none
  contains

  subroutine quad_bilinear_interp(lon_in, lat_in, x_corners_in, y_corners_in, p, expected_obs, d, do_rotate)
    implicit none
    real(8), intent(in) :: lon_in, lat_in, x_corners_in(4), y_corners_in(4), p(4)
    real(8), intent(out) :: expected_obs
    real(8), intent(out) :: d
    logical, intent(in) :: do_rotate
    real(8) :: m(3, 3), v(3), r(3), a, b(2), c(2)
    real(8) :: x_corners(4), lon, y_corners(4), lat, angle
    
    integer :: i

    lon = lon_in
    x_corners = x_corners_in
    lat = lat_in
    y_corners = y_corners_in

    if (do_rotate) then
        ! Rotate the quadrilateral
        do i = 2, 4
            x_corners(i) = x_corners(i) - x_corners(1)
            y_corners(i) = y_corners(i) - y_corners(1)
        end do
        lon = lon - x_corners(1)
        lat = lat - y_corners(1)
        x_corners(1) = 0.d0
        y_corners(1) = 0.d0

        b(1) = x_corners(2)
        b(2) = y_corners(2)

        ! Avoid degenerate cases where the grid is rotated exactly +/- 90 degrees
        if (abs(x_corners(2)) > 0.d001) then
            c(1) = x_corners(2)
            c(2) = 0.d0
        else
            c(1) = 0.d0
            c(2) = y_corners(2)
        end if

        angle = angle2(b, c)

        if (abs(angle) > 0.d001) then
            do i = 2, 4
                b(1) = x_corners(i)
                b(2) = y_corners(i)
                b = rotate2(b, angle)
                x_corners(i) = b(1)
                y_corners(i) = b(2)
            end do
            b(1) = lon
            b(2) = lat
            b = rotate2(b, angle)
            lon = b(1)
            lat = b(2)
        end if
    end if

    ! Fit a surface and interpolate; solve for 3x3 matrix
    do i = 1, 3
        m(i, 1) = x_corners(i) - x_corners(i + 1)
        m(i, 2) = y_corners(i) - y_corners(i + 1)
        m(i, 3) = x_corners(i) * y_corners(i) - x_corners(i + 1) * y_corners(i + 1)
        v(i) = p(i) - p(i + 1)
    end do

    d = deter3(m)
    ! Solve the matrix for b, c, and d
    call mat3x3(m, v, r)

    ! r contains b, c, and d; solve for a
    a = p(4) - r(1) * x_corners(4) - r(2) * y_corners(4) - r(3) * x_corners(4) * y_corners(4)

    ! Now do the interpolation
    expected_obs = a + r(1) * lon + r(2) * lat + r(3) * lon * lat
  end subroutine quad_bilinear_interp

  subroutine mat3x3(m, v, r)
    real(8), intent(in) :: m(3, 3), v(3)
    real(8), intent(out) :: r(3)
    real(8) :: m_sub(3, 3), numer, denom
    integer :: i

    denom = deter3(m)

    do i = 1, 3
      m_sub = m
      m_sub(:, i) = v
      numer = deter3(m_sub)
      r(i) = numer / denom
    end do

  end subroutine mat3x3

  function deter3(m)
    real(8), intent(in) :: m(3, 3)
    real(8) :: deter3

    deter3 = m(1, 1) * m(2, 2) * m(3, 3) + m(1, 2) * m(2, 3) * m(3, 1) + &
             m(1, 3) * m(2, 1) * m(3, 2) - m(3, 1) * m(2, 2) * m(1, 3) - &
             m(1, 1) * m(2, 3) * m(3, 2) - m(3, 3) * m(2, 1) * m(1, 2)

  end function deter3

  ! rotate vector a counterclockwise by angle theta (in radians)
  function rotate2(a, theta)
    real(8), intent(in) :: a(2)
    real(8), intent(in) :: theta
    real(8)             :: rotate2(2)
  
  real(8) :: r(2,2)
  
  r(1,1) = cos(theta)
  r(1,2) = sin(theta)
  r(2,1) = sin(-theta)
  r(2,2) = cos(theta)
  
  rotate2(1) = r(1,1)*a(1) + r(1,2)*a(2)
  rotate2(2) = r(2,1)*a(1) + r(2,2)*a(2)
  
  end function rotate2

  ! compute the angle between two 2-vectors
  function angle2(a, b)
    real(8), intent(in) :: a(2), b(2)
    real(8)             :: angle2
  
  angle2 = acos(dot2(a,b) / (mag2(a) * mag2(b)))
  
  end function angle2

  ! compute the magnitude of a 2-vector
    function mag2(a)
      real(8), intent(in) :: a(2)
      real(8)             :: mag2
    
    mag2 = sqrt(a(1)*a(1) + a(2)*a(2))
    
    end function mag2

   !Computes dot product of two 2-vectors

   function dot2(a, b)
    real(8), intent(in) :: a(2), b(2)
    real(8)             :: dot2
   
   dot2 = a(1)*b(1) + a(2)*b(2)
   
   end function dot2
   

end module quad_utils_mod
