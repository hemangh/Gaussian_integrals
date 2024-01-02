module conversions
    use precision, only: idp
    implicit none
    
contains

function cart2sph (xyz_coord) result(rthetaphi)
          
    real(idp)::xyz_coord (3), rthetaphi (3)

    real(idp), parameter :: epsilon = 1.d-15

    rthetaphi(1) = sqrt( dot_product (xyz_coord, xyz_coord))
    if (abs(xyz_coord(3)) <= epsilon .and. abs(rthetaphi(1)) <= epsilon) then
      rthetaphi(2) = 0.d0
    else
      rthetaphi(2) = acos (xyz_coord(3)/rthetaphi(1))
    endif
    
    rthetaphi(3) = atan2 (xyz_coord(2),xyz_coord(1))

end function  

function sph2cart (rthetaphi) result(xyz_coord)

real(idp)::xyz_coord (3), rthetaphi (3)

xyz_coord(1) = sin(rthetaphi(2)) * cos (rthetaphi(3)) * rthetaphi(1)
xyz_coord(2) = sin(rthetaphi(2)) * sin (rthetaphi(3)) * rthetaphi(1)
xyz_coord(3) = cos(rthetaphi(2)) * rthetaphi(1)

end function
    
end module conversions