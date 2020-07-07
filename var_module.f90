module var_module

    implicit none
    save

    integer :: ncomp, ots, option_dsbc
    real(kind=4) :: cfl, qus, f, yds, rhs1, rhs2, c11, c12, c21, c22, us, thes
    real(kind=4) :: vv, dtini

end module var_module
