program test_jacobian
  use quad_jacobian_mod
  implicit none
  real(8) :: x_corners(4), y_corners(4)
  real(8) :: xi, eta, jacobian_det

  ! Define the quadrilateral corners
  x_corners = [253.978338312172, 254.698325218660, 256.152459376863, 255.464594303337]
  y_corners = [74.4775053981321, 74.2803767364551, 74.6534788353970, 74.8496146897753]

  ! Define the local coordinates
  xi = 0.0d0
  eta = 0.0d0

  ! Compute the Jacobian determinant
  call compute_jacobian_det(x_corners, y_corners, xi, eta, jacobian_det)

  ! Print the result
  print *, "Jacobian determinant at (xi, eta) = (", xi, ",", eta, "): ", jacobian_det
end program test_jacobian