!> @file
!> @brief This is the subroutine to open input/output files.
!> @details This is the subroutine to open the input/output hydro files

      Subroutine OpenFiles

      use params

      implicit none

!>> Open input text files.

      open (unit=32, file= 'Cells.csv', status = 'unknown')
      open (unit=323, file ='Fetch.csv', status = 'unknown')
      open (unit=33, file= 'Links.csv', status = 'unknown')         ! JAM Oct 2010
      open (unit=34, file= 'LinksClosedHours.dat', status = 'unknown') !-EDW !value of 1 means link is closed for the hour, zero means it is open
      open (unit=35,file='LockControlObservedData.csv',status='unknown')
!      open (unit=36, file= 'SWR.dat', status = 'unknown')
      open (unit=74, file= 'MissRToC.csv', status = 'unknown')      ! JAM Oct 2010
      open (unit=39, file= 'TribQ.csv', status = 'unknown')   !Tributary flow (m3/s)
      open (unit=40, file= 'PET.csv', status = 'unknown')
      open (unit=42, file= 'Precip.csv', status='unknown')
      open (unit=45, file= 'Meteorology.csv', status = 'unknown')
      open (unit=43, file= 'WindVectorsX.csv',status= 'unknown')
      open (unit=46, file= 'WindVectorsY.csv',status= 'unknown')
      open (unit=47, file= 'TideData.csv',status='unknown')
      open (unit=48, file= 'TideTranspose.csv',status='unknown')
      open (unit=49, file= 'TideWeight.csv',status='unknown')
      open (unit=77, file= 'QMult.csv', form= 'formatted',status = 'unknown')
      open (unit=55, file= 'TribS.csv', status = 'unknown') ! Tributary sand concentration (mg/L)
      open (unit=555,file= 'TribF.csv',status='unknown')      ! Tributary fines concentration (mg/L)
      open (unit=56, file= 'SBC.dat', status = 'unknown')
      open (unit=110, file= 'surge.csv', form = 'formatted')
      open (unit=125, file= 'KBC.dat', status = 'unknown')              !node numbers of open boundary
      open (unit=101, file= 'BCToC2.dat', form = 'formatted')

!==== conditional input files
      if (Ndiv>0) then
          open (unit=86, file= 'DivQm.csv', status ='unknown')  ! diversion flow multiplier on Miss Riv flow
          open (unit=88, file= 'QMult_div.csv', form= 'formatted',status ='unknown')
          open (unit=89, file= 'DivSW.csv', status = 'unknown')
      endif
!      open (unit=117, file= 'AsedOW.csv',form ='formatted',status ='unknown')      ! Sediment Accretion  !Status='unknown' added by Joao Pereira 5/17/2011
      if (nlinksw>0) then
          open (unit=126, file= 'links_to_write.csv',status='unknown')      !input csv file with the link ID numbers of links to write flowrate to output file
      endif
      if (nstghr>0) then
          open (unit=127, file='hourly_stage_to_write.csv',status='unknown')  !input csv file with the ID numbers of compartments to write hourly stage to output file
      endif
      if (nlinklimiter>0) then
          open (unit=500, file= 'links_to_apply.csv',status='unknown')
      endif

!     If WQ modeling is excluded (iWQ=0), WQ inputs are disabled      
      if (iWQ>0) then
          open (unit=44, file= 'AnthL.csv', status = 'unknown')         ! Farm and Urban WW Loads (kg/d)
          open (unit=80, file= 'NO2NO3Data.csv', status = 'unknown')
          open (unit=81, file= 'NH4Data.csv', status = 'unknown')
          open (unit=82, file= 'OrgNData.csv', status = 'unknown')
          open (unit=83, file= 'PhosphorusData.csv', status ='unknown')
          open (unit=84, file= 'AtmChemData.csv', status ='unknown')
          open (unit=85, file= 'Decay.csv', status ='unknown')
          open (unit=87, file= 'DivWQ.csv', status ='unknown')
!          open (unit=118, file= 'UplandNP.dat', form ='formatted')
      endif
      
!     If 1D2D coupling enabled, 1D2Dcoupling input files
      if (n1D>0) then
          if (ntc>0) then
              open (unit=402, file= '1D2Dcoupling_tc.csv', status = 'unknown')
          endif
          if (nlc>0) then
              open (unit=403, file= '1D2Dcoupling_lc.csv', status = 'unknown')
          endif
          if (nuc>0) then
              open (unit=404, file= '1D2Dcoupling_uc.csv', status = 'unknown')
          endif
      endif   
      
!     If idt_schem = 2 user specified varying time step input file in daily interval
      if(idt_schem == 2) then
          open (unit=900, file= 'Dt_Varying_user.csv', status = 'unknown')
      endif

!>> Open output text files (in append mode, if needed).
      open (unit=75,file='SAL.out',form ='formatted',position='append')         ! Salinity.out
      open (unit=96,file='TSS.out',form='formatted',position='append')
      open (unit=100,file='TMP.out',form='formatted',position='append')        
      open (unit=103,file='SedAcc.out',form='formatted',position='append')      ! Last row will be used to compute open water Acc. 
      open (unit=104,file='SedAcc_MarshInt.out',form='formatted',position='append')     ! Last row will be used to compute interior marsh Acc. 
      open (unit=1045,file='SedAcc_MarshEdge.out',form='formatted',position='append')       ! Last row will be used to compute marsh edge Acc. 
      open (unit=105,file='fflood.out',form='formatted',position='append')
      open (unit=111,file='STG.out',form='formatted',position='append')         ! ESAVE.OUT
      open (unit=112,file='TRG.out',form='formatted',position='append')         ! Range.out
      open (unit=124,file='FLOm.out',form='formatted',position='append')
      open(unit=210,file='STGhr.out',form='formatted',position='append')                !output file for hourly water level in Boundary Condition cells
      open(unit=211,file='FLO.out',form='formatted',position='append')      !output file for flowrate   
      open(unit=212,file='STGm.out',form='formatted',position='append')
      open (unit=93,file='O2Sat.out',form='formatted',position='append')

!     If WQ modeling is excluded (iWQ=0), WQ outputs are disabled     
      if (iWQ>0) then
          open (unit=70,file='DIN.out',form ='formatted',position='append')
          open (unit=71,file='OrgN.out',form='formatted',position='append')
          open (unit=72,file='TPH.out',form='formatted',position='append')          ! TP.out
          open (unit=73,file='TOC.out',form='formatted',position='append')
          open (unit=91,file='NO3.out',form='formatted',position='append')          ! NO2NO3.out
          open (unit=92,file='NH4.out',form = 'formatted',position='append')        
          open (unit=94,file='ALG.out',form='formatted',position='append')
          open (unit=95,file='DO.out',form='formatted',position='append')
          open (unit=97,file='DET.out',form='formatted',position='append')          ! DeadAlgae.out
          open (unit=113,file='DON.out',form='formatted',position='append')
          open (unit=119,file='SPH.out',form='formatted',position='append')         ! SRP.out
          open (unit=121,file='NRM.out',form='formatted',position='append')         ! NRAcc.out -> Denitrification
          open (unit=123,file='TKN.out',form='formatted',position='append')
      endif
! read in information for grid cells used to pass data to other ICM routines !-EDW
      open (unit=200, file='grid_lookup_500m.csv', form='formatted')              ! compartment and link lookup table for 500-m grid cells
      open (unit=201, file='grid_interp_dist_500m.csv',form='formatted')          ! distance from each 500-m grid cell centroid to the compartment and link centroids
      open (unit=202, file='grid_data_500m.csv', form='formatted')                ! mean elevation for 500 m grid cells     
      open (unit=203, file='grid_IDs_Veg_matrix.csv', form='formatted')           ! 500m grid cell names formatted in the matrix used by Vegetation ICM routine

      open (unit=204, file='grid_500m_out.csv', form='formatted')                 ! output file for 500 m grid cells - in list form
      open (unit=205, file='compartment_out.csv',form='formatted')                ! output file for hydro compartments - summary values for ICM in list form

      open(unit=206,file='sal_monthly_ave_500m.out',form='formatted')
      open(unit=207,file='tmp_monthly_ave_500m.out',form='formatted')
      open(unit=208,file='tkn_monthly_ave_500m.out',form='formatted')
      open(unit=209,file='TSS_monthly_ave_500m.out',form='formatted')

! output files for use in the Vegetation ICM routine !-EDW
! these are written in append mode. ICM checks when first run as to whethere these files exist.
! If Ecohydro is run outside of the ICM these files may be erroneously appended to if they contain data and the model is re-run.
      open (unit=301, file=TRIM(ADJUSTL(VegMeanSalFile)),form='formatted', position='append')         ! mean salinity formatted for input into Vegetation ICM routine     
      open (unit=302, file=TRIM(ADJUSTL(VegSummerSalFile)),form='formatted', position='append')       ! mean summertime salinity for input into Vegetation ICM routine     
      open (unit=303, file=TRIM(ADJUSTL(VegSummerDepthFile)),form='formatted', position='append')     ! mean summertime water depth for input into Vegetation ICM routine   
      open (unit=304, file=TRIM(ADJUSTL(VegWaveAmpFile)),form='formatted', position='append')         ! variance in daily water depth for input into Vegetation ICM routine    
      open (unit=305, file=TRIM(ADJUSTL(VegSummerTempFile)),form='formatted', position='append')      ! mean summertime water temperature for input into Vegetation ICM routine
      open (unit=306, file=TRIM(ADJUSTL(VegTreeEstCondFile)),form='formatted', position='append')     ! tree establishment conditions for input into Vegetation ICM routine    
      open (unit=307, file=TRIM(ADJUSTL(VegBIHeightFile)),form='formatted', position='append')      ! tree establishment conditions for input into Vegetation ICM routine     
      open (unit=308, file=TRIM(ADJUSTL(VegPerLandFile)),form='formatted', position='append')      ! tree establishment conditions for input into Vegetation ICM routine       
! hotstart files
      open(unit=400,file='hotstart_in.dat',form='formatted')
      open(unit=401,file='hotstart_out.dat',form='formatted')


      return
      end
