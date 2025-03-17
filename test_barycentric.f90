program test_barycentric
  use quad_barycentric_mod
  implicit none
  real(8) :: x_corners(4), y_corners(4), values(4)
  real(8) :: lon, lat, interpolated_value

  ! Define the quadrilateral corners
  x_corners = [253.978338312172, 254.698325218660, 256.152459376863, 255.464594303337]
  y_corners = [74.4775053981321, 74.2803767364551, 74.6534788353970, 74.8496146897753]

  ! Define the values at the corners
  values = [0.04222093679422526, 0.02939444349332336, 0.02361927128232367, 0.02457275269829678]

  ! Define the interpolation point
  lon = 255.070312500000
  lat = 74.4735717773438

  ! Perform barycentric interpolation
  call quad_barycentric_interp(lon, lat, x_corners, y_corners, values, interpolated_value)

  ! Print the result
  print *, "Interpolated value at (", lon, ",", lat, "): ", interpolated_value
end program test_barycentric