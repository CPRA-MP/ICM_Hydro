module common_array_R

    implicit none
    save

    real(kind=4), allocatable :: y_R(:,:), area_R(:,:), q_R(:,:), hy_R(:,:), wl_lat_R(:,:), Q_lat_from_ICM_R(:,:)
	real(kind=4), allocatable :: Q_terminal_R(:), WL_terminal_from_ICM_R(:), Q_upstream_from_ICM_R(:)
	integer, allocatable :: ncomp_R(:), ndt_R(:), ioutf_R(:), nlat_R(:), n_R(:), n_SAL_R(:), n_TMP_R(:), n_FINE_R(:), n_SAND_R(:)
	character(len=128), allocatable :: Input_file_R(:), Input_SAL_R(:), Input_TMP_R(:), Input_FINE_R(:), Input_SAND_R(:)

! SAL TMP, FINE, SAND
	integer, allocatable :: ndt_SAL_R(:), ndt_TMP_R(:), ndt_FINE_R(:), ndt_SAND_R(:), ndt_total(:)
	integer, allocatable :: Nctr_SAL_R(:), Nctr_TMP_R(:), Nctr_FINE_R(:), Nctr_SAND_R(:)

	real(kind=4), allocatable :: sal_R(:,:), SAL_lat_from_ICM_R(:,:)
	real(kind=4), allocatable :: SAL_upstream_from_ICM_R(:), SAL_terminal_from_ICM_R(:)
	real(kind=4), allocatable :: tmp_R(:,:), TMP_lat_from_ICM_R(:,:)
	real(kind=4), allocatable :: TMP_upstream_from_ICM_R(:), TMP_terminal_from_ICM_R(:)
	real(kind=4), allocatable :: fine_R(:,:), FINE_lat_from_ICM_R(:,:)
	real(kind=4), allocatable :: FINE_upstream_from_ICM_R(:), FINE_terminal_from_ICM_R(:)
	real(kind=4), allocatable :: sand_R(:,:), SAND_lat_from_ICM_R(:,:)
	real(kind=4), allocatable :: SAND_upstream_from_ICM_R(:), SAND_terminal_from_ICM_R(:)
	

contains

    subroutine setup_common_var(num_reg)

        implicit none

        ! Input
        integer, intent(in) :: num_reg

		allocate(ncomp_R(num_reg))
		ncomp_R=0
		allocate(ndt_R(num_reg))
		ndt_R=1
		allocate(ioutf_R(num_reg))
		ioutf_R=0	
		allocate(nlat_R(num_reg))
		nlat_R=0			
		allocate(n_R(num_reg))
		n_R=0			
		allocate(n_SAL_R(num_reg))
		n_SAL_R=0			
		allocate(n_TMP_R(num_reg))
		n_TMP_R=0			
		allocate(n_FINE_R(num_reg))
		n_FINE_R=0			
		allocate(n_SAND_R(num_reg))
		n_SAND_R=0			
		allocate(Input_file_R(num_reg))
		allocate(Input_SAL_R(num_reg))
		allocate(Input_TMP_R(num_reg))
		allocate(Input_FINE_R(num_reg))
		allocate(Input_SAND_R(num_reg))

!SAL, TMP, FINE, SAND
		allocate(ndt_SAL_R(num_reg))
		ndt_SAL_R=1
		allocate(ndt_TMP_R(num_reg))
		ndt_TMP_R=1
		allocate(ndt_FINE_R(num_reg))
		ndt_FINE_R=1
		allocate(ndt_SAND_R(num_reg))
		ndt_SAND_R=1
		allocate(ndt_total(5*num_reg))
		ndt_total=1
		
		allocate(Nctr_SAL_R(num_reg))
		Nctr_SAL_R=0
		allocate(Nctr_TMP_R(num_reg))
		Nctr_TMP_R=0
		allocate(Nctr_FINE_R(num_reg))
		Nctr_FINE_R=0
		allocate(Nctr_SAND_R(num_reg))
		Nctr_SAND_R=0

		
	end subroutine setup_common_var
	
	

    subroutine setup_common_arrays_R(num_reg, num_points, num_lat)

        implicit none

        ! Input
        integer, intent(in) :: num_reg, num_points, num_lat

! R00
		allocate(y_R(num_reg, num_points))
		y_R=0.
		allocate(area_R(num_reg, num_points))
		area_R=0.
		allocate(q_R(num_reg, num_points))
		q_R=0.
		allocate(hy_R(num_reg, num_points))
		hy_R=0.
		
		allocate(wl_lat_R(num_reg, num_lat))
		wl_lat_R=-99999.
		allocate(Q_lat_from_ICM_R(num_reg, num_lat))
		Q_lat_from_ICM_R=-99999.
		allocate(WL_terminal_from_ICM_R(num_reg))
		WL_terminal_from_ICM_R=-99999.
		allocate(Q_upstream_from_ICM_R(num_reg))
		Q_upstream_from_ICM_R = -99999.
		allocate(Q_terminal_R(num_reg))
		Q_terminal_R = -99999.

		allocate(sal_R(num_reg, num_points))
		sal_R=0.
		allocate(tmp_R(num_reg, num_points))
		tmp_R=0.
		allocate(fine_R(num_reg, num_points))
		fine_R=0.
		allocate(sand_R(num_reg, num_points))
		sand_R=0.

		allocate(SAL_lat_from_ICM_R(num_reg, num_lat))
		SAL_lat_from_ICM_R=-99999.
		allocate(TMP_lat_from_ICM_R(num_reg, num_lat))
		TMP_lat_from_ICM_R=-99999.
		allocate(FINE_lat_from_ICM_R(num_reg, num_lat))
		FINE_lat_from_ICM_R=-99999.
		allocate(SAND_lat_from_ICM_R(num_reg, num_lat))
		SAND_lat_from_ICM_R=-99999.

		allocate(SAL_upstream_from_ICM_R(num_reg))
		SAL_upstream_from_ICM_R = -99999.
		allocate(TMP_upstream_from_ICM_R(num_reg))
		TMP_upstream_from_ICM_R = -99999.
		allocate(FINE_upstream_from_ICM_R(num_reg))
		FINE_upstream_from_ICM_R = -99999.
		allocate(SAND_upstream_from_ICM_R(num_reg))
		SAND_upstream_from_ICM_R = -99999.

		allocate(SAL_terminal_from_ICM_R(num_reg))
		SAL_terminal_from_ICM_R = -99999.
		allocate(TMP_terminal_from_ICM_R(num_reg))
		TMP_terminal_from_ICM_R = -99999.
		allocate(FINE_terminal_from_ICM_R(num_reg))
		FINE_terminal_from_ICM_R = -99999.
		allocate(SAND_terminal_from_ICM_R(num_reg))
		SAND_terminal_from_ICM_R = -99999.

		
	end subroutine setup_common_arrays_R
	
end module common_array_R

