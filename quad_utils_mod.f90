module quad_utils_mod
  implicit none
  contains

  subroutine quad_bilinear_interp(lon_in, lat_in, x_corners_in, y_corners_in, p, expected_obs)
    implicit none
    real(8), intent(in) :: lon_in, lat_in, x_corners_in(4), y_corners_in(4), p(4)
    real(8), intent(out) :: expected_obs
    real(8) :: m(3, 3), v(3), r(3), a, b(2), c(2), d
    real(8) :: x_corners(4), lon, y_corners(4), lat
    integer :: i

    lon = lon_in
    x_corners = x_corners_in
    lat = lat_in
    y_corners = y_corners_in

    ! See if the side wraps around in longitude.
    if (maxval(x_corners) - minval(x_corners) > 180.0_8) then
      if (lon < 180.0_8) lon = lon + 360.0_8
      do i = 1, 4
        if (x_corners(i) < 180.0_8) x_corners(i) = x_corners(i) + 360.0_8
      end do
    end if

    ! Fit a surface and interpolate; solve for 3x3 matrix
    do i = 1, 3
      m(i, 1) = x_corners(i) - x_corners(i + 1)
      m(i, 2) = y_corners(i) - y_corners(i + 1)
      m(i, 3) = x_corners(i) * y_corners(i) - x_corners(i + 1) * y_corners(i + 1)
      v(i) = p(i) - p(i + 1)
    end do

    ! Solve the matrix for b, c and d
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

end module quad_utils_mod
