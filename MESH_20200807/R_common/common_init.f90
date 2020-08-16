subroutine common_init(nr)

	use common_array_R
    implicit none
	integer, intent(out) :: nr
	integer :: i, NGCD_array, NGCD

! read dimensions for each region
	open(999,file='..\region_input.txt')
	read(999,*)nr, days_all, ndt_ICM

	call setup_common_var(nr)

	do i=1, nr
		read(999,*)ncomp_R(i), ndt_R(i), nlat_R(i)
		if(nlat_R(i).le.0)nlat_R(i)=1
		
		ndt_SAL_R(i)=ndt_R(i)
		ndt_TMP_R(i)=ndt_R(i)
		ndt_FINE_R(i)=ndt_R(i)
		ndt_SAND_R(i)=ndt_R(i)

		
		read(999,*)Input_file_R(i)
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
	close(999)

! change with R
	call setup_common_arrays_R(nr, maxval(ncomp_R), maxval(nlat_R))
!!!!!

!	ndt_all=NGCD_array(nr,ndt_R)
	ndt_all=NGCD_array(5*nr,ndt_total)
    ndt_all=NGCD(ndt_all, ndt_ICM)
	ntim_all=int(days_all*86400/ndt_all)

	do i=1, nr
		ioutf_R(i)=600+8*(i-1)
	enddo

end subroutine common_init
