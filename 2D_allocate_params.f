      subroutine allocate_params

      use params

      implicit none
      integer :: windsteps,tidesteps,maxconnectuse

      cells=N
      links=M
      windsteps = simdays*24/dtwind
      tidesteps = simdays*24/dttide+1    !YW! Tide transpose is assumed been handled in the input files. +1 is to include the final row
      maxconnectuse = max(maxconnect,100) ! upper limit on memory allocation for link connectivity matrices, icc and sicc
      numChem=14  !zw added 04/07/2020 to replace fixed variable dimensions related to chemicals

      WRITE(1,*) '----------------------------------------------------'
      WRITE(1,*) ' GRID AND RUN INFORMATION'
      WRITE(1,*) '----------------------------------------------------'
      write(1,*)
      write(1,1) ' Simulation time (days)          : ', simdays
      write(1,1) ' Hydrologic compartments         : ', cells
      write(1,1) ' Hydraulic links                 : ',links
      write(1,1) ' 500m ICM grid cells (with data) : ',n_500m_cells
      write(1,1) ' LA-Veg Model rows               : ',veg_matrix_rows
      write(1,1) ' LA-Veg Model columns            : ',veg_matrix_cols
      write(1,*)
      write(1,1) ' Links allowed per compartment   : ',maxconnect
      write(1,*)
      write(1,2) ' Water area required (sq km)     : ',minwater/1000000.
      write(1,*)
      write(1,2) ' Maximum change in stage allowed : ',maxdz


      WRITE(*,*) '----------------------------------------------------'
      WRITE(*,*) ' GRID AND RUN INFORMATION'
      WRITE(*,*) '----------------------------------------------------'
      write(*,*)
      write(*,1) ' Simulation time (days)          : ', simdays
      write(*,1) ' Hydrologic compartments         : ', cells
      write(*,1) ' Hydraulic links                 : ',links
      write(*,1) ' 500m ICM grid cells (with data) : ',n_500m_cells
      write(*,1) ' LA-Veg Model rows               : ',veg_matrix_rows
      write(*,1) ' LA-Veg Model columns            : ',veg_matrix_cols
      write(*,*)
      write(*,1) ' Links allowed per compartment   : ',maxconnect
      write(*,*)
      write(*,2) ' Water area required (sq km)     : ',minwater/1000000.
      write(*,*)
      write(*,2) ' Maximum change in stage allowed : ',maxdz
1     format(A,I0)
2     format(A,F0.2)

      ! lookup tables (dimension is number of grid cells) !-EDW
      allocate(grid_lookup_500m(n_500m_cells,22))
      allocate(grid_interp_dist_500m(n_500m_cells,22))
      allocate(bed_elev_500m(n_500m_cells))
      allocate(land_elev_500m(n_500m_cells))
      allocate(per_land_500m(n_500m_cells))
      allocate(per_water_500m(n_500m_cells))
      allocate(height_500m(n_500m_cells))

! output summary arrays for compartments !-EDW
      allocate(sal_2wk_ave_max(cells))
      allocate(sal_daily(simdays,cells))
      allocate(sal_ave(cells))
      allocate(sal_summer(125,cells))
      allocate(sal_ave_summer(cells))
      allocate(stage_daily(simdays,cells))
      allocate(stage_ave(cells))
      allocate(stage_max(cells))
      allocate(octapr_stage(cells))
      allocate(sepmar_stage(cells))
      allocate(stage_summer(125,cells))
      allocate(stage_ave_summer(cells))
      allocate(difmean(125,cells))
      allocate(stage_var_summer(cells))
      allocate(stage_stdv_summer(cells))
      allocate(trg_summer(125,cells))
      allocate(trg_ave_summer(cells))
      allocate(stage_wlv_summer(cells))
      allocate(tmp_daily(simdays,cells))
      allocate(tmp_ave(cells))
      allocate(tmp_summer(125,cells))
      allocate(tmp_ave_summer(cells))
      allocate(SedOW(cells))
      allocate(SedMarshInt(cells))
      allocate(SedMarshEdge(cells))
      allocate(tidal_range_daily(simdays,cells))
      allocate(dailyLW(cells))
      allocate(dailyHW(cells))
      allocate(TRsum(cells))
      allocate(TRave(cells))
      allocate(tidalprism_ave(cells))
      allocate(tkn_daily(simdays,cells))
      allocate(tss_daily(simdays,cells))
      allocate(tss_ave(cells))
      allocate(sal_month_ave(12,cells))
      allocate(stg_month_ave(12,cells))
      allocate(tmp_month_ave(12,cells))
      allocate(tkn_month_ave(12,cells))
      allocate(tss_month_ave(12,cells))
      allocate(tss_var_annual(cells))
      allocate(difmean_tss(simdays,cells))



!output summary arrays for links  !-EDW
      allocate(sal_2wk_ave_max_links(links))
      allocate(sal_daily_links(simdays,links))
      allocate(sal_ave_links(links))
      allocate(sal_summer_links(125,links))
      allocate(sal_ave_summer_links(links))
      allocate(tmp_daily_links(simdays,links))
      allocate(tmp_ave_links(links))
      allocate(tmp_summer_links(125,links))
      allocate(tmp_ave_summer_links(links))
      allocate(tkn_daily_links(simdays,links))
      allocate(tss_daily_links(simdays,links))
      allocate(sal_month_ave_links(12,links))
      allocate(tmp_month_ave_links(12,links))
      allocate(tkn_month_ave_links(12,links))
      allocate(tss_month_ave_links(12,links))



      allocate(SlkAve(links))
      allocate(TlkAve(links))

!output summary arrays for 500 m grid !-EDW
      allocate(pct_sand_bed_500m(n_500m_cells))
      allocate(salinity_500m(n_500m_cells))
      allocate(salinity_IDW_500m(n_500m_cells))
      allocate(salinity_summer_500m(n_500m_cells))
      allocate(salinity_summer_IDW_500m(n_500m_cells))
      allocate(sal_thresh_500m(n_500m_cells))
      allocate(sal_thresh_IDW_500m(n_500m_cells))
      allocate(stage_500m(n_500m_cells))
      allocate(stage_summer_500m(n_500m_cells))
      allocate(stage_wlv_summer_500m(n_500m_cells))
      allocate(depth_500m(n_500m_cells))
      allocate(depth_summer_500m(n_500m_cells))
      allocate(tmp_500m(n_500m_cells))
      allocate(tmp_IDW_500m(n_500m_cells))
      allocate(tmp_summer_500m(n_500m_cells))
      allocate(tmp_summer_IDW_500m(n_500m_cells))
      allocate(tree_est(n_500m_cells))

      allocate(sal_IDW_500m_month(12,n_500m_cells))
      allocate(tmp_IDW_500m_month(12,n_500m_cells))
      allocate(tkn_500m_month(12,n_500m_cells))
      allocate(tss_500m_month(12,n_500m_cells))



!output summary arrays for 500 m grid in Veg model matrix
      allocate(veg_grid_IDs(veg_matrix_cols,veg_matrix_rows))
      allocate(depth_summer_forVeg(veg_matrix_cols,veg_matrix_rows))
      allocate(salinity_forVeg(veg_matrix_cols,veg_matrix_rows))
      allocate(salinity_summer_forVeg(veg_matrix_cols,veg_matrix_rows))
      allocate(stage_var_forVeg(veg_matrix_cols,veg_matrix_rows))
      allocate(tmp_summer_forVeg(veg_matrix_cols,veg_matrix_rows))
      allocate(tree_est_forVeg(veg_matrix_cols,veg_matrix_rows))
      allocate(ht_abv_water_forVeg(veg_matrix_cols,veg_matrix_rows))
      allocate(per_land_forVeg(veg_matrix_cols,veg_matrix_rows))

      ! sediment resupsension/erodible bed terms
      allocate(depo_on_off(cells))
      allocate(erBedAvail(cells,4))
      allocate(erBedDepth(cells))
      allocate(erBedBD(cells))
      allocate(CSSresusOff(4))

      ! variables of (:) dimensions

      ! arrays of various length
      allocate(stghrwrite(nstghr))
      allocate(EShrly(nstghr))
      allocate(linkswrite(nlinksw))
      allocate(FLO(nlinksw))
      allocate(Dgr(4))    !sediment class parameters
      allocate(D90(4))    !sediment class parameters
      allocate(rhoSed(4)) !sediment class parameters
      allocate(SGsed(4))  !sediment class parameters

      allocate(dCSS(4))
      allocate(dCSSh(4))
      allocate(SedAccumRate(4))
      allocate(SedAccumRate_h(4))
      allocate(dSacc(4))
      allocate(dSacch(4))
      allocate(QSmarsh(4))     !used to be an array the size of compartments - now only sized for the 4 sediment classes
      allocate(QStrib(4))
      allocate(QSdiv(4))
      allocate(QSsum(4))
      allocate(QSsumh(4))
      allocate(netflux(4))
      allocate(daccm(4))
      allocate(BCSedRatio(4))
      allocate(MEESedRatio(4))
      allocate(MEESedRate(4))
      allocate(month_DOY(12))
      allocate(KBC(mds))  !allocate(KBC(20))  !zw change to mds instead of 20  04/07/2020
      allocate(Sal(1,1))
      allocate(SWR(50))

      !zw change to numChem instead of 15 04/07/2020
      allocate(cChemface(numChem))  !allocate(cChemface(15))
      allocate(DChem(numChem))  !allocate(DChem(15))
      allocate(QChemSUM(numChem)) !allocate(QChemSUM(15))
      allocate(QChemSUManth(numChem)) !allocate(QChemSUManth(15))
      allocate(QChemSUMtrib(numChem)) !allocate(QChemSUMtrib(15))
      allocate(QChemSUMdiv(numChem)) !allocate(QChemSUMdiv(15))
      allocate(QChemSUMflows(numChem)) !allocate(QChemSUMflows(15))
      allocate(QChemSUMatm(numChem)) !allocate(QChemSUMatm(15))
      allocate(Ws(5))

      ! arrays of length equal to number of days simulated
      allocate(ParSandD(simdays))
      allocate(lnO2Sat(simdays))
      allocate(O2Sat(simdays))
      allocate(ta(simdays))
      allocate(T7(simdays))
      allocate(ta_k(simdays))
      allocate(TempMR(simdays))
      allocate(tw(simdays))
      allocate(wd(simdays))
!      allocate(wspd(simdays))

      ! arrays equal to number of gages
      allocate(windx_data(windsteps,windgages))
      allocate(windy_data(windsteps,windgages))
      allocate(PET(simdays,etgages))
      allocate(ETA(etgages))
      allocate(rain(simdays,raingages))
      allocate(TideData(tidesteps,tidegages))
      allocate(transposed_tide(tidegages,2))
      allocate(weighted_tide(Mds-tidegages,5))
      allocate(surge(tidesteps,Mds)) !-EDW


      ! arrays of length equal to number of cells
      allocate(nlink2cell(cells))
      allocate(ka(cells))
      allocate(cf(cells))
      allocate(sedn(cells))
      allocate(sedcalib(cells))
      allocate(alphaSed(cells))
      allocate(accsed(cells))
!      allocate(acss(cells))
      allocate(acssh(cells))
      allocate(Ahf(cells))
      allocate(Ahydro(cells))
      allocate(Apctupland(cells))
      allocate(Apctmarsh(cells))
      allocate(Apctwater(cells))
      allocate(Atotal(cells))
      allocate(AM1(cells))
      allocate(AM2(cells))
      allocate(AnthL(cells))
      allocate(ASandA(cells))
!      allocate(bcss(cells))
      allocate(Bed(cells))
      allocate(BedM(cells))
      allocate(BedMOrig(cells))
      allocate(BedMSD(cells))
      allocate(BedMAdj(cells))
      allocate(ChemDONo(cells))
      allocate(ChemDOPo(cells))
      allocate(ChLAo(cells))
      allocate(CSSo(cells))
      allocate(dAdz(cells))
      allocate(CSSos(cells))
      allocate(Esho(cells))
      allocate(Eso(cells))
!      allocate(flood(cells))
      allocate(floodf(cells))
      allocate(hLength(cells))
!      allocate(han(cells))
      allocate(HDTo(cells))
!      allocate(hkm(cells))
      allocate(hWidth(cells))
      allocate(ar_int(cells))
      allocate(ar_ed(cells))
      allocate(jbar(cells))
      allocate(jdiv(cells))
      allocate(Jrain(cells))
      allocate(Jwind(cells))
      allocate(Jet(cells))
      allocate(jtrib(Ntrib))  !ZW 12/12/2023 jtrib is the tributary IDs
      allocate(Percent(cells))
      allocate(phz(cells))
      allocate(por(cells))
      allocate(SBC(cells))
      allocate(SOD(cells))
      allocate(SourceBM(cells))
      allocate(SRPo(cells))
      allocate(SSource(cells))      !this was original SSource(150) in comdeck.h - effected calculation results
!      allocate(tauc(cells))
      allocate(TM1(cells))
      allocate(TM2(cells))
      allocate(TSS(cells))
      allocate(TSSAve(cells))
!      allocate(Vss(cells))
!      allocate(Vsh(cells))
      allocate(windx(cells))
      allocate(windy(cells))
      allocate(Qsum_in(cells))
      allocate(Qsum_out(cells))
!      allocate(Qsum_out_links(cells))
      allocate(Qsum_abs(cells))
      allocate(KKa(cells))
      allocate(KKdepth(cells))
      allocate(MEE(cells))
      allocate(cumul_retreat(cells))
      allocate(Sacch_int(cells,2))
      allocate(Sacch_edge(cells,2))
      allocate(flag_offbc(cells))  !zw offshore bc cells flag 04/07/2020
      allocate(adaption_coeff(cells))   !non-equilibrium adaption coefficient for sand resuspension/deposition source term


      ! arrays of length equal to number of links
      allocate(USx(links))
      allocate(USy(links))
      allocate(DSx(links))
      allocate(DSy(links))
      allocate(linkt(links))
      allocate(Latr1(links))
      allocate(Latr2(links))
      allocate(Latr3(links))
      allocate(Latr4(links))
      allocate(Latr5(links))
      allocate(Latr6(links))
      allocate(Latr7(links))
      allocate(Latr8(links))
      allocate(Latr9(links))
      allocate(Latr10(links))
      allocate(Latr11(links))      
      allocate(fa_mult(links))
      allocate(hourclosed(links,24))

      allocate(fa(links))
      allocate(fb(links))
      allocate(ACCSEDj(links))
!      allocate(Achan(links))
!   allocate(an(links))
      allocate(asedout(links))  !this was originally asedout(cells) in comdeck.h - should be allocated for links, not cells
      allocate(BCage(cells))
      allocate(BCDA(cells))
      allocate(BCDO(cells))
      allocate(BCDIN(cells))
      allocate(BCLA(cells))
      allocate(BCNH4(cells))
      allocate(BCNO3(cells))
      allocate(BCON(cells))
      allocate(BCTOC(cells))
      allocate(BCTP(cells))
      allocate(BCTSS(cells))
!      allocate(Cinvert(links))
      allocate(CssL(links))     !this was originally CssL(30) - should be allocated for links
!      allocate(Deptho(links))
      allocate(EAOL(links))
      allocate(Exy(links))
!      allocate(itype(links))
      allocate(jds(links))
      allocate(jus(links))
!      allocate(Ken(links))
!      allocate(Km(links))
!      allocate(Kx(links))
!      allocate(Length(links))
!      allocate(NR(links))
!      allocate(Pmsh(links))
      allocate(Qsed(links))     !this was originally Qsed(30) - should be allocated for links
      allocate(Qsedsm(links))
!      allocate(Resist(links))    !no longer an array - now a local calculation
!      allocate(Width(links))

! variables of (:,:) dimensions

      ! arrays of various length
!      allocate(BL(150,7))
      allocate(cssFines(Ntrib))
      allocate(cssT(Ntrib,simdays,4))
      allocate(cSSTdiv(Ndiv,simdays,4))
      allocate(SWRsand(Ndiv))
      allocate(SWRfines(Ndiv))
      allocate(DivMult(Ndiv))
      allocate(daymo(12,2))
      allocate(decay(numChem,numChem))  !allocate(decay(25,25)) !zw changed to numChem instead of 25 04/07/2020
      allocate(Qdiv(Ndiv,simdays))
      allocate(QssT(Ntrib,simdays))
      allocate(QssTdiv(Ndiv,simdays))
      allocate(Qtrib(Ntrib,simdays))
      allocate(Strib(Ntrib,simdays))

      ! arrays of length equal to number of cells
      allocate(kinvisc(cells))
      allocate(group_vel(cells,2))
      allocate(Hs(cells,2))
      allocate(Uorb(cells,2))
      allocate(wave_energy(cells,2))
      allocate(wave_frequency(cells,2))
      allocate(wavelength(cells,2))
      allocate(wave_period(cells,2))
      allocate(velset(cells,4))   !sediment class parameters

      allocate(Aas(cells,3))
      allocate(As(cells,3))
      allocate(ASandD(cells,simdays))
      allocate(ASandT(cells,simdays))
      allocate(barp(cells,simdays))
      allocate(ChemDON(cells,3))
      allocate(ChemDOP(cells,3))
      allocate(ChemPOP(cells,3))
      allocate(ChLa(cells,3))
      allocate(CSS(cells,3,4)) !third array dimension is for 4 different sediment classes
      allocate(CSSf(cells,3))
      allocate(CSSh(cells,3,4)) !third array dimension is for 4 different sediment classes
      allocate(denit(cells,2))
      allocate(Eh(cells,3))
      allocate(ds(cells,3))
      allocate(Es(cells,3))
      allocate(EHAV(cells,3))
      allocate(BCnosurge(cells,2))
      allocate(BCsurge(cells,2))   !YW! added for tide calculation
      allocate(ESAV(cells,3))
      allocate(ESMN(cells,3))
      allocate(ESMX(cells,3))
      allocate(Fetch(cells,16)) !allocate(Fetch(cells,10))!zw changed to 16 instead of 10 04/07/2020
      allocate(HDT(cells,3))
      allocate(hmarsho(cells,3))
      allocate(icc(cells,maxconnectuse))
      allocate(Qmarsh(cells,2))
      allocate(Qmarshmax(cells))
      allocate(QmarshAve(cells))
      allocate(QMult(cells,Ntrib))
      allocate(qmultdiv(cells,Ndiv))
      allocate(EsRange(cells,2))
      allocate(S(cells,3))
      allocate(SALAV(cells))  !zw added 3/22/2015 for daily average salinity from each time step
!      allocate(SA(cells,12))   !not used in program
      allocate(Sacc(cells,2))
      allocate(Sacch(cells,2))
      allocate(Sandacc(cells,2))
      allocate(Siltacc(cells,2))
      allocate(Clayacc(cells,2))
      allocate(pct_sand_bed(cells))
      allocate(Sh(cells,3))
      allocate(Shg(cells,2))
      allocate(sicc(cells,maxconnectuse))
      allocate(SRP(cells,3))
      allocate(STEMP(cells,4))
!      allocate(tau(cells,3))
      allocate(Tempair(cells,simdays))
      allocate(Tempe(cells,simdays))
      allocate(resuspension(cells,4))     !sediment array - 4 columns for 4 sediment classes
      allocate(deposition(cells,4))     !sediment array - 4 columns for 4 sediment classes
      allocate(CSSvRs(cells,2))
      allocate(QRain(cells))            !ZW 1/31/2024 save total runoff volume (m3/s) at each time step 


      ! !EDW Tempw() used to be used to read in the boundary condition water temps, and then also used as the cell values of temperature at each timestep
      ! !EDW Essentially the first 7 columns were used to read in the BC, then columns 1 and 2 were then re-used as the temperature array for compartments
      ! !EDW Created a new array so that they are now separate.

      allocate(TempwBC(mds,simdays))     ! !BUG! this was originally (cells,3) - but error was thrown on infile - should be a timeseries
      allocate(Tempw(cells,3))     ! !BUG! this was originally (cells,3) - but error was thrown on infile - should be a timeseries
      allocate(TempwAve(cells))   !-EDW array to save average temperature values for the day


      allocate(Tmp(cells,2))
      allocate(Tmpe(cells,2))     ! !-EDW I don't think this is used
      allocate(Tmpm(cells,2))     ! !-EDW I don't think this is used
      !allocate(TMtrib(30,simdays))!-EDW   array size of 30 is too small when running in Debug mode
      allocate(TMtrib(cells,simdays))!-EDW

      allocate(ave_vel_sum(cells))
      allocate(ave_vel_cnt(cells))
      allocate(ave_vel(cells))
      allocate(max_vel(cells))
      allocate(min_vel(cells))
      

      ! arrays of length equal to number of links
      allocate(Depth(links,3))
      allocate(Q(links,3))
!      allocate(QSal(links,3))
!      allocate(Qss(links,3))
      allocate(SL(links,2))
      allocate(Tcoef(links,25))
      allocate(TL(links,3))
      allocate(link_vel(links))
!      allocate(uplandNP(links,14))

! variables of (:,:,:) dimensions

! arrays of various length realted to number of chemicals
!zw changed to numChem instead of numbers such as 15 or 25, 04/07/2020
      allocate(cCHEM(Ntrib,numChem,simdays))  !allocate(cCHEM(Ntrib,15,simdays))
      allocate(cChemdiv(1,numChem,simdays))  !allocate(cChemdiv(1,15,simdays))
      allocate(dccc1(1,numChem,2))  !allocate(dccc1(1,25,2))
      allocate(dccc2(1,numChem,2))   !allocate(dccc2(1,25,2))

      allocate(QAtm(1,numChem,simdays))  !allocate(QAtm(1,15,simdays))
      allocate(QChemdiv(Ndiv,numChem,simdays))   !allocate(QChemdiv(Ndiv,15,simdays))

! arrays of length equal to number of cells
      allocate(Chem(cells,numChem,2))                !allocate(Chem(cells,15,2))
      allocate(ChemAve(cells,numChem))               !allocate(ChemAve(cells,15))
      allocate(GrowAlgae(cells,numChem,numChem))     !allocate(GrowAlgae(cells,25,25))
      allocate(GrowChlA(cells,numChem,numChem))      !allocate(GrowChlA(cells,25,25))
      allocate(QCHEM(Ntrib,numChem,simdays))         !allocate(QCHEM(Ntrib,15,simdays))

 6666 format(A,I3,A,I3)

      allocate(linkslimiter(nlinklimiter))  !YW! store link numbers to apply flow limiter
      allocate(flag_apply(links))           !zW! store flag to apply flow limiter in links

!1D-ICM coupling variables
      allocate(tcr2D(ntc))  ! terminal connection ICM receiving compartment list
      allocate(tcf2D(ntc))  ! terminal connection ICM connecting compartment list
      allocate(tcl2D(ntc))  ! terminal connection ICM link list      
      allocate(tcn1D(ntc))  ! terminal connection 1D node list
      allocate(tcr1D(ntc))  ! terminal connection 1D region list      
      
      allocate(lcr2D(nlc))  ! lateral connection ICM receiving compartment list
      allocate(lcf2D(nlc))  ! lateral connection ICM connecting compartment list
      allocate(lcl2D(nlc))  ! lateral connection ICM link list
      allocate(lcn1D(nlc))  ! lateral connection 1D node list
      allocate(lcr1D(nlc))  ! lateral connection 1D region list      
      
      allocate(ucr2D(nuc))  ! upstream connection ICM upstream compartment list
      allocate(ucf2D(nuc))  ! upstream connection ICM connecting compartment list
      allocate(ucl2D(nuc))  ! upstream connection ICM link list
      allocate(ucn1D(nuc))  ! upstream connection 1D node list  
      allocate(ucr1D(nuc))  ! upstream connection 1D regio list        
      
      allocate(tcH(ntc)) ! terminal connection stage (applied to the connecting compartment from ICM)
      allocate(tcQ(ntc)) ! terminal connection discharge (applied to the connecting link from 1D)
      allocate(lcH(nlc)) ! lateral connection stage (applied to the connecting compartment from 1D)
      allocate(lcQ(nlc)) ! lateral connection discharge (calculated)
      allocate(ucH(nuc)) ! upstream connection stage (applied to the connecting compartment from 1D)
      allocate(ucQ(nuc)) ! upstream connection discharge (calculated)    
      
      allocate(NTs_Ratio(n_region)) ! time step ratio between 1D river and ICM

!     If idt_schm = 2 user specified varying time step input file in daily interval
      if(idt_schm == 2) then
          allocate(dt_var_user(simdays))  !variable time step in daily interval
      endif

      return
      end
