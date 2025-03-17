module quad_jacobian_mod
  implicit none
  contains

  ! Subroutine to compute the Jacobian determinant for a quadrilateral
  subroutine compute_jacobian_det(x_corners, y_corners, xi, eta, jacobian_det)
    implicit none
    real(8), intent(in) :: x_corners(4), y_corners(4) ! Quadrilateral corners
    real(8), intent(in) :: xi, eta                   ! Local coordinates (-1 <= xi, eta <= 1)
    real(8), intent(out) :: jacobian_det             ! Jacobian determinant
    real(8) :: dN_dxi(4), dN_deta(4)                ! Partial derivatives of shape functions
    real(8) :: dx_dxi, dx_deta, dy_dxi, dy_deta     ! Partial derivatives of x and y

    ! Compute partial derivatives of shape functions
    dN_dxi(1) = -0.25d0 * (1.0d0 - eta)
    dN_dxi(2) =  0.25d0 * (1.0d0 - eta)
    dN_dxi(3) =  0.25d0 * (1.0d0 + eta)
    dN_dxi(4) = -0.25d0 * (1.0d0 + eta)

    dN_deta(1) = -0.25d0 * (1.0d0 - xi)
    dN_deta(2) = -0.25d0 * (1.0d0 + xi)
    dN_deta(3) =  0.25d0 * (1.0d0 + xi)
    dN_deta(4) =  0.25d0 * (1.0d0 - xi)

    ! Compute partial derivatives of x and y
    dx_dxi = sum(dN_dxi(:) * x_corners(:))
    dx_deta = sum(dN_deta(:) * x_corners(:))
    dy_dxi = sum(dN_dxi(:) * y_corners(:))
    dy_deta = sum(dN_deta(:) * y_corners(:))

    ! Compute the Jacobian determinant
    jacobian_det = dx_dxi * dy_deta - dx_deta * dy_dxi
  end subroutine compute_jacobian_det

end module quad_jacobian_mod