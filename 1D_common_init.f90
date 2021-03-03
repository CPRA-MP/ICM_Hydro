subroutine common_init(nr)

    use common_array_R
    use params
    implicit none
    integer, intent(out) :: nr
    integer :: i, j, NGCD_array, NGCD, max_nlat

! read dimensions for each region
	open(999,file='../region_input.txt')
!	read(999,*)nr, days_all, ndt_ICM
	read(999,*)nr
    days_all = simdays

	call setup_common_var(nr)

	do i=1, nr
		read(999,*)ncomp_R(i), ndt_R(i), nlat_R(i)
		!if(nlat_R(i).le.0)nlat_R(i)=1
		
		ndt_SAL_R(i)=ndt_R(i)
		ndt_TMP_R(i)=ndt_R(i)
		ndt_FINE_R(i)=ndt_R(i)
		ndt_SAND_R(i)=ndt_R(i)

		
		read(999,*)Input_file_R(i)
        write(*,*) i,Input_file_R(i)
        
		read(999,*)Nctr_SAL_R(i), Nctr_TMP_R(i), Nctr_FINE_R(i), Nctr_SAND_R(i)
		if(Nctr_SAL_R(i).eq.1)then
			read(999,*)ndt_SAL_R(i)
			read(999,*)Input_SAL_R(i)
		endif
		if(Nctr_TMP_R(i).eq.1)then
			read(999,*)ndt_TMP_R(i)
			read(999,*)Input_TMP_R(i)
		endif
		if(Nctr_FINE_R(i).eq.1)then
			read(999,*)ndt_FINE_R(i)
			read(999,*)Input_FINE_R(i)
		endif
		if(Nctr_SAND_R(i).eq.1)then
			read(999,*)ndt_SAND_R(i)
			read(999,*)Input_SAND_R(i)
		endif
		
		ndt_total((i-1)*5+1)=ndt_R(i)
		ndt_total((i-1)*5+2)=ndt_SAL_R(i)
		ndt_total((i-1)*5+3)=ndt_TMP_R(i)
		ndt_total((i-1)*5+4)=ndt_FINE_R(i)
		ndt_total((i-1)*5+5)=ndt_SAND_R(i)
	enddo

! change with R
	nlat_R_ori=nlat_R
	max_nlat=maxval(nlat_R)
	if(max_nlat.le.0)max_nlat=1
	call setup_common_arrays_R(nr, maxval(ncomp_R), max_nlat)
!!!!!

!	ndt_all=NGCD_array(nr,ndt_R)
! ndt_all and ntim_all are set in main.f for 2D modeling only - they will be updated here if there are 1D reaches included in the model run
	ndt_all=NGCD_array(5*nr,ndt_total)
    ndt_all=NGCD(ndt_all, ndt_ICM)
	ntim_all=int(days_all*86400/ndt_all)
	do i=1, nr
		ioutf_R(i)=600+8*(i-1)
	enddo

	
! Read lateral info
	read(999,*)
	
	do i=1, nr

    if(nlat_R(i)>0)then
        write(*,*) 'region', i,nlat_R(i)
		read(999,*) (latFlowLoc_R(i,j), j=1, nlat_R(i))
		do j=1,nlat_R(i)
			if ((latFlowLoc_R(i,j)-1)*(latFlowLoc_R(i,j)-ncomp_R(i)) .eq. 0) then
				print*, 'Region ',i, ' ERROR: Lateral flow cannot be applied at the boundaries'
				stop
			endif
		enddo
		
		read(999,*) (latFlowType_R(i,j), j=1, nlat_R(i))
		do j=1,nlat_R(i)
			if (latFlowType_R(i,j) .eq. 1) then
				print*, 'Region ',i,' Lateral flow at node = ', latFlowLoc_R(i,j), ', is a time series'
			elseif (latFlowType_R(i,j) .eq. 2) then
				print*, 'Region ',i,' Lateral flow at node = ', latFlowLoc_R(i,j), ', is a function of upstream flow'
			elseif (latFlowType_R(i,j) .eq. 3) then
                print*, 'Region ',i,' Lateral flow at node = ', latFlowLoc_R(i,j), ', is a coupling point with ICM model'
			else
				print*, 'Region ',i,' Wrong lateral flow type is provided. Type ', latFlowType_R(i,j), 'is not a valid type'
				stop
			endif
		enddo

		read(999,*) (latFlowXsec_R(i,j), j=1, nlat_R(i))	
	else
		nlat_R(i)=1
	endif
	

	enddo

	close(999)


end subroutine common_init
