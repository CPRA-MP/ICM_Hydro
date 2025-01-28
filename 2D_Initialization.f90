!> @file
!> @brief This is the subroutine to initialize 1D & 2D models.
!> @details This is the subroutine to initialize 1D & 2D models 
!> (assigning initial conditions for 2D link & compartments)

      subroutine Initialization

      use params
      use common_array_R  ! 1D-ICM coupling

      implicit none
      integer :: i,j,jj,jjj,sedclass,iir,ichem,me       !iterators
     
!>> Initialize 1D model subroutines
      if (n1D > 0) then
          call common_init(n_region)              ! n1d is read from RuncontrolR.data and n_region is read in from region_input.txt during configuration of 1D model
          if (n1D /= n_region) then
              write(*,*) 'Number of 1D regions defined in ICM (RuncontrolR.dat) and MESH (region_input.txt) are not consistent.'
              write(1,*) 'Number of 1D regions defined in ICM (RuncontrolR.dat) and MESH (region_input.txt) are not consistent.'
              stop !pause
          endif
          if (nlc /= sum(nlat_R(1:n_region))) then
              write(*,*) 'Number of 1D lateral connections defined in ICM (RuncontrolR.dat) and MESH (region_input.txt) are not consistent.'
              write(1,*) 'Number of 1D lateral connections defined in ICM (RuncontrolR.dat) and MESH (region_input.txt) are not consistent.'
              stop !pause
          endif	  
          do iir=1, n_region
              write(*,*) 'Initializing 1D model for region: ',iir
              write(1,*) 'Initializing 1D model for region: ',iir
              call init_R(iir, ncomp_R(iir),ioutf_R(iir), Input_file_R(iir), days_all, ndt_R(iir), nlat_R_ori(iir), nlat_R(iir), latFlowLoc_R(iir,1:nlat_R(iir)), latFlowType_R(iir,1:nlat_R(iir)), latFlowXsec_R(iir,1:nlat_R(iir)), y_R(iir,1:ncomp_R(iir)), q_R(iir,1:ncomp_R(iir)), area_R(iir,1:ncomp_R(iir)), hy_R(iir,1:ncomp_R(iir)), wl_lat_R(iir,1:nlat_R(iir)), Q_terminal_R(iir))
              if (Nctr_SAL_R(iir) .eq. 1)call init_SAL_R(iir, ncomp_R(iir), ioutf_R(iir)+4, Input_SAL_R(iir), days_all, ndt_SAL_R(iir), nlat_R_ori(iir), nlat_R(iir), latFlowLoc_R(iir,1:nlat_R(iir)), latFlowType_R(iir,1:nlat_R(iir)), latFlowXsec_R(iir,1:nlat_R(iir)))
              if (Nctr_TMP_R(iir) .eq. 1)call init_TMP_R(iir, ncomp_R(iir), ioutf_R(iir)+5, Input_TMP_R(iir),days_all, ndt_TMP_R(iir), nlat_R_ori(iir), nlat_R(iir), latFlowLoc_R(iir,1:nlat_R(iir)), latFlowType_R(iir,1:nlat_R(iir)), latFlowXsec_R(iir,1:nlat_R(iir)))
              if (Nctr_FINE_R(iir) .eq. 1)call init_FINE_R(iir, ncomp_R(iir), ioutf_R(iir)+6, Input_FINE_R(iir), days_all, ndt_FINE_R(iir), nlat_R_ori(iir), nlat_R(iir), latFlowLoc_R(iir,1:nlat_R(iir)), latFlowType_R(iir,1:nlat_R(iir)), latFlowXsec_R(iir,1:nlat_R(iir)))
              if (Nctr_SAND_R(iir) .eq. 1)call init_SAND_R(iir, ncomp_R(iir), ioutf_R(iir)+7, Input_SAND_R(iir),days_all, ndt_SAND_R(iir), nlat_R_ori(iir), nlat_R(iir), latFlowLoc_R(iir,1:nlat_R(iir)), latFlowType_R(iir,1:nlat_R(iir)), latFlowXsec_R(iir,1:nlat_R(iir)))
          enddo
          ndt_all_ICM = ndt_all
          ntim_all_ICM = ntim_all          
      else
          write(*,*) 'No 1D regions defined in ICM (RuncontrolR.dat, n1D=0) - only 2D compartments will be modeled.'
          write(1,*) 'No 1D regions defined in ICM (RuncontrolR.dat, n1D=0) - only 2D compartments will be modeled.'
      endif

!>> Set initial conditions for 2D links
      do i=1,M
          fa(i) = fa_def*fa_mult(i) !Set array of initial upwind factor to default value
          Q(i,1)=0.0	!YW!
          Q(i,2)=Q(i,1)

          ! MP2023 zw added 04/06/2020
          EAOL(i)=0.0
          SlkAve(i) = 0
          SL(i,1) = 0  ! SL is for links
          SL(i,2)=0
          TL(i,1)=0
          asedout(i)=0.0
      enddo

      if(nlinksw > 0) then               !YW!
          do jj = 1,nlinksw
              FLO(jj) = 0.0
          enddo
      endif

!>> Initialize 2D compartment related arrays
      do j=1,N
          As(j,2)= As(j,1)				    ! Surface water area of cells (m2)
          Es(j,2)= Es(j,1)				    ! Stage in storage cells (m)
          ds(j,1)= Es(j,1)-Bed(j)			! Depth in storage cells (m)
          !Eh(j,1) = BedM(j) + 0.1           ! Initial marsh depth (override hotstart file read in above)
          Eh(j,2)=Eh(j,1)					! Stage in Marsh storage (m)	!JAM Oct 2010
          BCnosurge(j,1) = 0.0              ! Initialize no surge BC to zero for all compartments - only BC nodes will be updated - rest of array will be 0.0
          BCnosurge(j,2) = BCnosurge(j,1)   ! boundary conditions stage(m) before surge is added
          BCsurge(j,1) = 0.0                ! Initialize surge BC to zero for all compartments - only BC nodes will be updated - rest of array will be 0.0 -YW
          BCsurge(j,2) = BCsurge(j,1)
          ESMX(j,2)=ES(j,1)
          ESMN(j,2)=ES(j,1)
          dailyHW(j)=0.0
          dailyLW(j)=0.0
          ESAV(j,1) = ES(j,1)*dt/(3600.*24.)
          EHAV(j,1) = EH(j,1)*dt/(3600.*24.)
          floodf(j)=0.0
          
          Qmarsh(j,1) = 0.0				    ! Flow into/out of marsh area
          Qmarsh(j,2) = Qmarsh(j,1)		    
          Qmarshmax(j) = 0.0
          QmarshAve(j) = Qmarsh(j,1)*dt/(3600.*24.)

          S(j,2) = S(j,1)
          SALAV(j) = S(j,1)*dt/(3600.*24.)
          sal_ave(j)=0.0
          Tempw(j,2) = Tempw(j,1)
          TempwAve(j) = Tempw(j,1)*dt/(3600.*24.)

          do sedclass=1,4
               CSS(j,2,sedclass) = CSS(j,1,sedclass)
               CSSh(j,1,sedclass) = CSS(j,1,sedclass)
               CSSh(j,2,sedclass) = CSS(j,1,sedclass)
          enddo
          accsed(j)=0.0
          cumul_retreat(j) = 0.0
          Sacc(j,1)=0.0
          Sacch_int(j,1)=0.0
          Sacch_edge(j,1)=0.0
          Sandacc(j,1) = 0.0
          Siltacc(j,1) = 0.0
          Clayacc(j,1) = 0.0
          CSSvRs(j,1)= 0.0
          Sacc(j,2)=Sacc(j,1)
          Sacch_int(j,2)=Sacch_int(j,1)
          Sacch_edge(j,2)=Sacch_edge(j,1)
          Sandacc(j,2) = Sandacc(j,1)
          Siltacc(j,2) = Siltacc(j,1)
          Clayacc(j,2) = Clayacc(j,1)
          CSSvRs(j,2)= CSSvRs(j,1)
          TSSave(j) = ( CSS(j,1,1) + CSS(j,1,2) + CSS(j,1,3) + CSS(j,1,4) )*dt/(3600.*24.)

          do ichem = 1,14
              Chem(j,ichem,2)=Chem(j,ichem,1)
              ChemAve(j,ichem) = Chem(j,ichem,1)*dt/(3600.*24.)
              ! initially set GrowAlgae array equal to zero
              do me=1,14
                  GrowAlgae(j,ichem,me) = 0.
              enddo
          enddo
          denit(j,1)=0
          denit(j,2)=0

      enddo

!>> initialize BCnosurge and BCsurge with initial tide and surge       -YW
      do jj=1,tidegages
          BCnosurge(transposed_tide(jj,1),1) = TideData(1,jj)
          do jjj=1,Mds
              if (KBC(jjj)==transposed_tide(jj,1)) then
                  BCsurge(transposed_tide(jj,1),1)= Surge(1,jjj)
              endif
          enddo
      enddo
              
!>> Take first timestep of imported wind data and save into windx and windy arrays.
!>> These arrays will be overwritten at a delta t that matches the wind data timestep (this update occurs immediately prior to calling hydrod)
      do jjj = 1,N
!          windx(jjj) = max(0.01,windx_data(1,jwind(jjj)))
!          windy(jjj) = max(0.01,windy_data(1,jwind(jjj)))
          windx(jjj) = windx_data(1,jwind(jjj))
          windy(jjj) = windy_data(1,jwind(jjj))
      enddo

      return
      end
