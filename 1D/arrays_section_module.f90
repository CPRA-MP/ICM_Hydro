module arrays_section_module

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

    integer :: maxTableLength, nel, Q_sk_tableEntry

    character*4 :: file_num
    character(len=128) :: xSection_path

contains

    subroutine setup_arrays_section

        implicit none

        allocate(elevTable(nel))
        allocate(areaTable(nel))
        allocate(pereTable(nel))
        allocate(rediTable(nel))
        allocate(convTable(nel))
        allocate(topwTable(nel))
        allocate(nwi1Table(nel))
        allocate(dPdATable(nel))

		allocate(ncompElevTable(nel))
        allocate(ncompAreaTable(nel))

! only in version 2 20151115
        allocate(I2Tablep(nel))
        allocate(I2Tablec(nel))
        allocate(upstreamI2Tablec(nel))
        allocate(downstreamI2Tablep(nel))
        allocate(currentSquareDepth(nel))
        allocate(downstreamSquareDepth(nel))
        allocate(upstreamSquareDepth(nel))

    end subroutine setup_arrays_section
end module
