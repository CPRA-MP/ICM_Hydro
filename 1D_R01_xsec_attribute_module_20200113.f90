module xsec_attribute_module_R01

    implicit none
    save

    real(kind=4), allocatable :: xsec_tab(:, :, :)

    contains

    ! Allocate storage for all of the arrays in this module based on the number
    ! of time steps and spatial points
    subroutine setup_xsec_attribute_module(elements, num_points)
        implicit none

        ! Input
        integer, intent(in) :: elements, num_points

        allocate(xsec_tab(10, elements, num_points))
		
		xsec_tab=0.

        ! Headers:
        ! Elev(m)    Area(m2)     Peri(m)      Radi(m)   Conv(m3/s)    topWidth(m)    newI1(m3)    dPdA(1/m)     dbdxp(-)    dbdxc(-)      I2p     I2c



    end subroutine
end module xsec_attribute_module_R01
