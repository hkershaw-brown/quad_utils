module quad_barycentric_mod
  implicit none
  contains

  ! Subroutine to compute barycentric interpolation for a quadrilateral
  subroutine quad_barycentric_interp(lon, lat, x_corners, y_corners, values, interpolated_value)
    implicit none
    real(8), intent(in) :: lon, lat                  ! Interpolation point (longitude, latitude)
    real(8), intent(in) :: x_corners(4), y_corners(4) ! Quadrilateral corners
    real(8), intent(in) :: values(4)                ! Values at the corners
    real(8), intent(out) :: interpolated_value      ! Interpolated value at (lon, lat)
    real(8) :: lambda(4)                            ! Barycentric coordinates
    integer :: i

    ! Compute barycentric coordinates
    call compute_barycentric_coords(lon, lat, x_corners, y_corners, lambda)

    ! Interpolate the value using the barycentric coordinates
    interpolated_value = 0.0d0
    do i = 1, 4
      interpolated_value = interpolated_value + lambda(i) * values(i)
    end do
  end subroutine quad_barycentric_interp

  ! Subroutine to compute barycentric coordinates for a quadrilateral
  subroutine compute_barycentric_coords(lon, lat, x_corners, y_corners, lambda)
    implicit none
    real(8), intent(in) :: lon, lat                  ! Interpolation point (longitude, latitude)
    real(8), intent(in) :: x_corners(4), y_corners(4) ! Quadrilateral corners
    real(8), intent(out) :: lambda(4)               ! Barycentric coordinates
    real(8) :: area_total, area_sub(4)              ! Total area and sub-areas
    integer :: i

    ! Compute the total area of the quadrilateral using the shoelace formula
    area_total = compute_quad_area(x_corners, y_corners)

    ! Compute the sub-areas formed by the interpolation point and each triangle
    do i = 1, 4
      area_sub(i) = compute_triangle_area(lon, lat, x_corners(mod(i, 4) + 1), y_corners(mod(i, 4) + 1), &
                                          x_corners(mod(i + 1, 4) + 1), y_corners(mod(i + 1, 4) + 1))
    end do

    ! Compute barycentric coordinates as the ratio of sub-areas to the total area
    do i = 1, 4
      lambda(i) = area_sub(i) / area_total
    end do
  end subroutine compute_barycentric_coords

  ! Function to compute the area of a quadrilateral using the shoelace formula
  function compute_quad_area(x_corners, y_corners) result(area)
    implicit none
    real(8), intent(in) :: x_corners(4), y_corners(4) ! Quadrilateral corners
    real(8) :: area                                   ! Total area of the quadrilateral
    integer :: i

    area = 0.0d0
    do i = 1, 4
      area = area + (x_corners(i) * y_corners(mod(i, 4) + 1) - y_corners(i) * x_corners(mod(i, 4) + 1))
    end do
    area = 0.5d0 * abs(area)
  end function compute_quad_area

  ! Function to compute the area of a triangle given its three vertices
  function compute_triangle_area(x1, y1, x2, y2, x3, y3) result(area)
    implicit none
    real(8), intent(in) :: x1, y1, x2, y2, x3, y3 ! Triangle vertices
    real(8) :: area                               ! Area of the triangle

    area = 0.5d0 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
  end function compute_triangle_area

end module quad_barycentric_mod