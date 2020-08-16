module var_module_R01

    implicit none
    save

    integer :: ntim, ncomp, ots, option_dsbc, saveFrequency, noLatFlow, noQSKtable, igate, boundaryFileMaxEntry
	integer :: ppp, qqq
    real(kind=4) :: cfl, qus, f, yds, rhs1, rhs2, c11, c12, c21, c22, us, thes
    real(kind=4) :: vv, dtini, skk, thetas, thesinv, yy, qq, timesDepth

end module var_module_R01
