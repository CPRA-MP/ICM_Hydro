module arrays_module_R02

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
    real(kind=4), allocatable :: dt(:), dx(:), froud(:), courant(:)

    real(kind=4), allocatable :: USBoundary(:,:), DSBoundary(:,:)
! change for unsteady flow
    real(kind=4), allocatable :: pere(:),dpda(:), hy(:)

    real(kind=4), allocatable :: oldQ(:), newQ(:), oldArea(:), newArea(:), oldY(:), newY(:)

    integer, allocatable :: ityp(:), latFlowLocations(:), dataInEachLatFlow(:), latFlowType(:), latFlowXsecs(:)

    real(kind=4), allocatable :: lateralFlowTable(:,:,:), lateralFlow(:)

    integer, allocatable :: Q_SK_tableEntry(:)
    real(kind=4), allocatable :: eachQSKtableNodeRange(:,:), Q_SK_Table(:,:,:)

contains

    ! Allocate storage for all of the arrays in this module based on the number
    ! of time steps and spatial points
    subroutine setup_arrays(num_time, num_points, maxTableEntry1, maxTableEntry2, totalLatFlow, totalQSKtable)

        implicit none

        ! Input
        integer, intent(in) :: num_time, num_points, maxTableEntry1, maxTableEntry2, totalLatFlow, totalQSKtable

        allocate(area(num_points))
		area=0.
! change for unsteady flow

        allocate(bo(num_points))
		bo=0.

        allocate(pere(num_points))
		pere=0.
        allocate(dpda(num_points))
		dpda=0.
        allocate(hy(num_points))
		hy=0.


        allocate(areap(num_points))
		areap=0.
        allocate(qp(num_points))
		qp=0.
        allocate(z(num_points))
		z=0.
        allocate(dqp(num_points))
		dqp=0.
        allocate(av11(num_points))
		av11=0.
        allocate(av12(num_points))
		av12=0.
        allocate(av21(num_points))
		av21=0.
        allocate(av22(num_points))
		av22=0.
        allocate(dqc(num_points))
		dqc=0.
        allocate(dap(num_points))
		dap=0.
        allocate(dac(num_points))
		dac=0.
        allocate(ci1(num_points))
		ci1=0.
        allocate(ci2(num_points))
		ci2=0.
        allocate(aso(num_points))
		aso=0.
        allocate(depth(num_points))
		depth=0.
        allocate(f1(num_points))
		f1=0.
        allocate(f2(num_points))
		f2=0.
        allocate(g11inv(num_points))
		g11inv=0.
        allocate(g12inv(num_points))
		g12inv=0.
        allocate(g21inv(num_points))
		g21inv=0.
        allocate(g22inv(num_points))
		g22inv=0.
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
        allocate(dt(num_points))
        allocate(ityp(num_points))
        allocate(dx(num_points-1))
		b11=0.
        b12=0.
        b21=0.
        b22=0.
        eps2=0.
        eps4=0.
        d1=0.
        d2=0.
        u=0.
        c=0.
        sk=0.
        co=0.
        gso=0.
        dbdx=0.
        dt=0.
        ityp=0
        dx=0.

        allocate(froud(num_points))
		froud=0.

        allocate(Q_SK_Table(2, maxTableEntry1, totalQSKtable))
        allocate(Q_SK_tableEntry(totalQSKtable))

		allocate(USBoundary(2, maxTableEntry2))
		USBoundary=0.
        allocate(DSBoundary(2, maxTableEntry2))
		DSBoundary=0.

        allocate(courant(num_points-1))
		courant=0.

        allocate(oldQ(num_points))
		oldQ=0.
        allocate(newQ(num_points))
		newQ=0.
        allocate(oldArea(num_points))
		oldArea=0.
        allocate(newArea(num_points))
		newArea=0.
        allocate(oldY(num_points))
		oldY=0.
        allocate(newY(num_points))
		newY=0.

        allocate(lateralFlowTable(2, maxTableEntry2, totalLatFlow))
		lateralFlowTable=0.
        allocate(dataInEachLatFlow(totalLatFlow))
		dataInEachLatFlow=0
        allocate(lateralFlow(num_points))
		lateralFlow=0.

    end subroutine setup_arrays

end module arrays_module_R02
