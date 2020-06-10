module arrays_module

    implicit none
    save

    real(kind=4), allocatable :: area(:), bo(:) !, y(:, :), q(:, :)
    real(kind=4), allocatable :: areap(:), qp(:), z(:), dqp(:)
    real(kind=4), allocatable :: av11(:), av12(:), av21(:), av22(:)
    real(kind=4), allocatable :: dqc(:), dap(:), dac(:), ci1(:), ci2(:)
    real(kind=4), allocatable :: aso(:), f1(:), f2(:), depth(:)
    real(kind=4), allocatable :: g11inv(:), g12inv(:), g21inv(:), g22inv(:)
    real(kind=4), allocatable :: b11(:), b12(:), b21(:), b22(:)
    real(kind=4), allocatable :: eps2(:), eps4(:), d1(:), d2(:), u(:), c(:)
    real(kind=4), allocatable :: sk(:), co(:), gso(:), dbdx(:)
    real(kind=4), allocatable :: dt_1d(:), dx(:), froud(:), courant(:)

    real(kind=4), allocatable :: USBoundary(:,:), DSBoundary(:,:), Q_Sk_Table(:,:)
! change for unsteady flow
    real(kind=4), allocatable :: pere(:),dpda(:)

    real(kind=4), allocatable :: oldQ(:), newQ(:), oldArea(:), newArea(:), oldY(:), newY(:)

    integer, allocatable :: ityp(:), latFlowLocations(:), dataInEachLatFlow(:), latFlowType(:), latFlowXsecs(:)

    real(kind=4), allocatable :: lateralFlowTable(:,:,:), lateralFlow(:)

contains

    ! Allocate storage for all of the arrays in this module based on the number
    ! of time steps and spatial points
    subroutine setup_arrays(num_time, num_points, maxTableEntry1, maxTableEntry2, totalLatFlow)

        implicit none

        ! Input
        integer, intent(in) :: num_time, num_points, maxTableEntry1, maxTableEntry2, totalLatFlow

        allocate(area(num_points))

! change for unsteady flow

        allocate(bo(num_points))

        allocate(pere(num_points))
        allocate(dpda(num_points))


        allocate(areap(num_points))
        allocate(qp(num_points))
        allocate(z(num_points))
        allocate(dqp(num_points))
        allocate(av11(num_points))
        allocate(av12(num_points))
        allocate(av21(num_points))
        allocate(av22(num_points))
        allocate(dqc(num_points))
        allocate(dap(num_points))
        allocate(dac(num_points))
        allocate(ci1(num_points))
        allocate(ci2(num_points))
        allocate(aso(num_points))
        allocate(depth(num_points))
        allocate(f1(num_points))
        allocate(f2(num_points))
        allocate(g11inv(num_points))
        allocate(g12inv(num_points))
        allocate(g21inv(num_points))
        allocate(g22inv(num_points))
        allocate(b11(num_points))
        allocate(b12(num_points))
        allocate(b21(num_points))
        allocate(b22(num_points))
        allocate(eps2(num_points))
        allocate(eps4(num_points))
        allocate(d1(num_points))
        allocate(d2(num_points))
        allocate(u(num_points))
        allocate(c(num_points))
        allocate(sk(num_points))
        allocate(co(num_points))
        allocate(gso(num_points))
        allocate(dbdx(num_points))
        allocate(dt_1d(num_points))
        allocate(ityp(num_points))
        allocate(dx(num_points-1))

        allocate(froud(num_points))

        allocate(Q_Sk_Table(2, maxTableEntry1))
        allocate(USBoundary(2, maxTableEntry2))
        allocate(DSBoundary(2, maxTableEntry2))

        allocate(courant(num_points-1))

        allocate(oldQ(num_points))
        allocate(newQ(num_points))
        allocate(oldArea(num_points))
        allocate(newArea(num_points))
        allocate(oldY(num_points))
        allocate(newY(num_points))

        allocate(lateralFlowTable(2, maxTableEntry2, totalLatFlow))
        allocate(dataInEachLatFlow(totalLatFlow))
        allocate(lateralFlow(num_points))

    end subroutine setup_arrays

end module arrays_module
