module arrays_section_module_R07

    implicit none
    save

    real(kind=4), allocatable :: elevTable(:),areaTable(:)
    real(kind=4), allocatable :: pereTable(:),rediTable(:)
    real(kind=4), allocatable :: convTable(:),topwTable(:)
    real(kind=4), allocatable :: nwi1Table(:),dPdATable(:)
    real(kind=4), allocatable :: ncompElevTable(:), ncompAreaTable(:)
! only in version 2 20191511
    real(kind=4), allocatable :: I2Tablep(:),I2Tablec(:)
    real(kind=4), allocatable :: upstreamI2Tablec(:), downstreamI2Tablep(:)
    real(kind=4), allocatable :: currentSquareDepth(:), downstreamSquareDepth(:), upstreamSquareDepth(:)

    integer :: maxTableLength, nel

    character*4 :: file_num
    character(len=128) :: xSection_path

contains

    subroutine setup_arrays_section

        implicit none

        allocate(elevTable(nel))
		elevTable=0.
        allocate(areaTable(nel))
		areaTable=0.
        allocate(pereTable(nel))
		pereTable=0.
        allocate(rediTable(nel))
		rediTable=0.
        allocate(convTable(nel))
		convTable=0.
        allocate(topwTable(nel))
		topwTable=0.
        allocate(nwi1Table(nel))
		nwi1Table=0.
        allocate(dPdATable(nel))
		dPdATable=0.

		allocate(ncompElevTable(nel))
		ncompElevTable=0.
        allocate(ncompAreaTable(nel))
		ncompAreaTable=0.

! only in version 2 20151115
        allocate(I2Tablep(nel))
		I2Tablep=0.
        allocate(I2Tablec(nel))
		I2Tablec=0.
        allocate(upstreamI2Tablec(nel))
		upstreamI2Tablec=0.
        allocate(downstreamI2Tablep(nel))
		downstreamI2Tablep=0.
        allocate(currentSquareDepth(nel))
		currentSquareDepth=0.
        allocate(downstreamSquareDepth(nel))
		downstreamSquareDepth=0.
        allocate(upstreamSquareDepth(nel))
		upstreamSquareDepth=0.

    end subroutine setup_arrays_section
end module
