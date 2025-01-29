!> @file
!> @brief This is the subroutine to write output for other ICM routines.
!> @details This is the subroutine to post-process hydro model outputs to be read into
!> other ICM routines.

      subroutine ICM_PP

      use params

      implicit none
      integer :: j,k,kk,kj,kl      !iterators
     
!>> Call 'ICM_Summaries' subroutine which calculates the average output values used by the other ICM routines (Wetland Morph, Vegetation, and Ecosytems)
      call ICM_Summaries

!>> Call 'ICM_SalinityThreshold_WM' subroutine which calculates the maximum two-week salinity average for each hydro compartment (used in  Marsh Collapse calculations in Wetland Morph)
      call ICM_SalinityThreshold_WM

!>> Call 'ICM_MapToGrid' subroutine which maps (without interpolation) the hydro compartment results to the gridded structrure used by the other ICM routines.
      call ICM_MapToGrid(1,500)               ! Map annual mean salinity to grid
      call ICM_MapToGrid(2,500)               ! Map summer salinity to grid
      call ICM_MapToGrid(3,500)               ! Interpolate annual mean temperature to grid
      call ICM_MapToGrid(4,500)               ! Interpolate summer mean temperature to grid
      call ICM_MapToGrid(5,500)               ! Map annual mean water stage to grid
      call ICM_MapToGrid(6,500)               ! Map summer water stage to grid
      call ICM_MapToGrid(7,500)               ! Map summer water stage variance to grid
      call ICM_MapToGrid(8,500)               ! Map summer maximum 2-wk mean salinity to grid
      call ICM_MapToGrid(9,500)               ! Map percent sand in bed to grid

      do kk=301,312
          call ICM_MapMonthlytoGrid(kk,500)       ! Map monthly TKN values to grid
      enddo

      do kk=401,412
          call ICM_MapMonthlytoGrid(kk,500)       ! Map monthly TSS values to grid
      enddo

!>> Calculate water depth for 500 m-grid cells.
      do j=1,n_500m_cells
          depth_500m(j) = stage_500m(j)-bed_elev_500m(j)
          depth_summer_500m(j) = stage_summer_500m(j)-bed_elev_500m(j)
          height_500m(j) = land_elev_500m(j)-stage_500m(j)
      enddo

!>> Call 'ICM_InterpolateToGrid' subroutine which interpolates hydro compartment and link results to the gridded structure used by the other ICM routines.
      call ICM_InterpolateToGrid(1,500)       ! Interpolate annual mean salinity to grid
      call ICM_InterpolateToGrid(2,500)       ! Interpolate summer mean salinity to grid
      call ICM_InterpolateToGrid(3,500)       ! Interpolate annual mean temperature to grid
      call ICM_InterpolateToGrid(4,500)       ! Interpolate summer mean temperature to grid
      call ICM_InterpolateToGrid(5,500)       ! Interpolate summer maximum 2-wk mean salinity to grid

      do kk=101,112
          call ICM_InterpolateMonthlyToGrid(kk,500)       ! Interpolate monthly mean salinity values to grid
      enddo

      do kk=201,212
          call ICM_InterpolateMonthlyToGrid(kk,500)       ! Interpolate monthly mean temperatures to grid
      enddo

!>> Call 'ICM_TreeConditions_Veg' subroutine which determines if tree establishment conditions are met for each grid cell of the Veg model.
      call ICM_TreeConditions_Veg

!>> Call 'ICM_Formatting' subroutine which formats the gridded output into versions directly digestable by other ICM routines.
      call ICM_Formatting

!>> Calculate average tidal prism for each compartment.
      do j = 1,N
          TRsum(j) = 0.0
          do kl=1,simdays
              TRsum(j) = TRsum(j) + tidal_range_daily(kl,j)
          enddo
          TRave(j) = TRsum(j)/simdays
          tidalprism_ave(j) = TRave(j)*As(j,1)
      enddo

!>> Write compartment summary output to file in list form - one row for each compartment
!>> THESE HEADERS ARE USED BY OTHER ICM ROUTINES - DO NOT CHANGE WITHOUT UPDATING ICM.PY, HSI.PY, & WM.PY
      write(205,1117)'ICM_ID',      &
             'max_annual_stage',        &
             'ave_annual_stage',        &
             'ave_stage_summer',        &
             'var_stage_summer',        &
             'ave_annual_salinity',     &
             'ave_salinity_summer',     &
             'sal_2wk_max',     &
             'ave_tmp',     &
             'ave_tmp_summer',      &
             'openwater_sed_accum',     &
             'marsh_int_sed_accum',     &
             'marsh_edge_sed_accum',        &
             'tidal_prism_ave',     &
             'ave_sepmar_stage',        &
             'ave_octapr_stage',        &
             'marsh_edge_erosion_rate',     &
             'ave_annual_tss',      &
             'stdev_annual_tss',        &
             'totalland_m2'

      do kj=1,N
          write(205,1118) kj,       &
             max(-stagemax,min(stagemax,stage_max(kj))),        &
             max(-stagemax,min(stagemax,stage_ave(kj))),        &
             max(-stagemax,min(stagemax,stage_ave_summer(kj))),     &
             !max(-rangemax,min(rangemax,stage_var_summer(kj))),        &
             max(-rangemax,min(rangemax,stage_wlv_summer(kj))),     &
             min(salmax,sal_ave(kj)),       &
             min(salmax,sal_ave_summer(kj)),        &
             min(salmax,sal_2wk_ave_max(kj)),       &
             min(tmpmax,tmp_ave(kj)),       &
             min(tmpmax,tmp_ave_summer(kj)),        &
             max(-sedmaxow,min(sedmaxow,Sacc(kj,2))),       &
             max(-sedmaxmi,min(sedmaxmi,Sacch_int(kj,2))),      &
             max(-sedmaxme,min(sedmaxme,Sacch_edge(kj,2))),     &
             tidalprism_ave(kj),        &
             max(-stagemax,min(stagemax,sepmar_stage(kj))),     &
             max(-stagemax,min(stagemax,octapr_stage(kj))),     &
             MEE(kj),       &
             tss_ave(kj),       &
             tss_var_annual(kj)**0.5,       &  !stdev = sqrt(variance)
             max(0.0,(Atotal(kj)-As(kj,1)))
      enddo
      close(205)
!>> Write gridded output to file in list form - one row for each grid cell.
!>> Write header rows in grid output files
!>> THESE HEADERS ARE USED BY OTHER ICM ROUTINES - DO NOT CHANGE WITHOUT UPDATING ICM.PY, HSI.PY, & WM

      write(204,1115)'GRID_ID',     &
        'compartment_ave_salinity_ppt',     &
         'IDW_ave_salinity_ppt',        &
         'compartment_ave_summer_salinity_ppt',     &
         'IDW_ave_summer_salinity_ppt',     &
         'compartment_max_2wk_summer_salinity_ppt',     &
         'IDW_max_2wk_summer_salinity_ppt',     &
         'bed_pct_sand',        &
         'compartment_ave_temp',        &
         'IDW_ave_temp',        &
         'compartment_ave_summer_temp',     &
         'IDW_ave_summer_temp',     &
         'stage_ave',       &
         'stage_summer_ave',        &
         'WLV_stage_summer',        &
         'ave_depth_summer',        &
         'ave_depth'
      write(206,1119)'GRID_ID',     &
         'sal_ave_jan',     &
         'sal_ave_feb',     &
         'sal_ave_mar',     &
        'sal_ave_apr',      &
         'sal_ave_may',     &
         'sal_ave_jun',     &
         'sal_ave_jul',     &
        'sal_ave_aug',      &
         'sal_ave_sep',     &
         'sal_ave_oct',     &
         'sal_ave_nov',     &
        'sal_ave_dec'
      write(207,1119)'GRID_ID',     &
         'tmp_ave_jan',     &
         'tmp_ave_feb',     &
        'tmp_ave_mar',      &
         'tmp_ave_apr',     &
         'tmp_ave_may',     &
         'tmp_ave_jun',     &
        'tmp_ave_jul',      &
         'tmp_ave_aug',     &
         'tmp_ave_sep',     &
         'tmp_ave_oct',     &
        'tmp_ave_nov',      &
         'tmp_ave_dec'
      write(208,1119)'GRID_ID',     &
         'tkn_ave_jan',     &
         'tkn_ave_feb',     &
        'tkn_ave_mar',      &
         'tkn_ave_apr',     &
         'tkn_ave_may',     &
         'tkn_ave_jun',     &
        'tkn_ave_jul',      &
         'tkn_ave_aug',     &
         'tkn_ave_sep',     &
         'tkn_ave_oct',     &
        'tkn_ave_nov',      &
         'tkn_ave_dec'
      write(209,1119)'GRID_ID',     &
         'TSS_ave_jan',     &
         'TSS_ave_feb',     &
         'TSS_ave_mar',     &
        'TSS_ave_apr',      &
         'TSS_ave_may',     &
         'TSS_ave_jun',     &
         'TSS_ave_jul',     &
        'TSS_ave_aug',      &
         'TSS_ave_sep',     &
         'TSS_ave_oct',     &
         'TSS_ave_nov',     &
        'TSS_ave_dec'

    ! write various summary results in 500m grid output file
      do k=1,n_500m_cells
          write(204,1116) grid_lookup_500m(k,1),     &
           min(salmax,salinity_500m(k)),        &
           min(salmax,salinity_IDW_500m(k)),        &
           min(salmax,salinity_summer_500m(k)),     &
           min(salmax,salinity_summer_IDW_500m(k)),     &
           min(salmax,sal_thresh_500m(k)),      &
           min(salmax,sal_thresh_IDW_500m(k)),      &
           pct_sand_bed_500m(k),        &
           min(tmpmax,tmp_500m(k)),     &
           min(tmpmax,tmp_IDW_500m(k)),     &
           min(tmpmax,tmp_summer_500m(k)),      &
           min(tmpmax,tmp_summer_IDW_500m(k)),      &
           max(-stagemax,min(stagemax,stage_500m(k))),      &
           max(-stagemax,min(stagemax,stage_summer_500m(k))),       &
           max(-rangemax,min(rangemax,stage_wlv_summer_500m(k))),       &
           max(-depthmax,min(depthmax,depth_summer_500m(k))),       &
           max(-depthmax,min(depthmax,depth_500m(k)))

          write(206,1120) k,     &
           min(salmax,sal_IDW_500m_month(1,k)),     &
           min(salmax,sal_IDW_500m_month(2,k)),     &
           min(salmax,sal_IDW_500m_month(3,k)),     &
           min(salmax,sal_IDW_500m_month(4,k)),     &
           min(salmax,sal_IDW_500m_month(5,k)),     &
           min(salmax,sal_IDW_500m_month(6,k)),     &
           min(salmax,sal_IDW_500m_month(7,k)),     &
           min(salmax,sal_IDW_500m_month(8,k)),     &
           min(salmax,sal_IDW_500m_month(9,k)),     &
           min(salmax,sal_IDW_500m_month(10,k)),        &
           min(salmax,sal_IDW_500m_month(11,k)),        &
           min(salmax,sal_IDW_500m_month(12,k))

          write(207,1120) k,     &
           min(tmpmax,tmp_IDW_500m_month(1,k)),     &
           min(tmpmax,tmp_IDW_500m_month(2,k)),     &
           min(tmpmax,tmp_IDW_500m_month(3,k)),     &
           min(tmpmax,tmp_IDW_500m_month(4,k)),     &
           min(tmpmax,tmp_IDW_500m_month(5,k)),     &
           min(tmpmax,tmp_IDW_500m_month(6,k)),     &
           min(tmpmax,tmp_IDW_500m_month(7,k)),     &
           min(tmpmax,tmp_IDW_500m_month(8,k)),     &
           min(tmpmax,tmp_IDW_500m_month(9,k)),     &
           min(tmpmax,tmp_IDW_500m_month(10,k)),        &
           min(tmpmax,tmp_IDW_500m_month(11,k)),        &
           min(tmpmax,tmp_IDW_500m_month(12,k))

          write(208,1120) k,     &
           min(tknmax,tkn_500m_month(1,k)),     &
           min(tknmax,tkn_500m_month(2,k)),     &
           min(tknmax,tkn_500m_month(3,k)),     &
           min(tknmax,tkn_500m_month(4,k)),     &
           min(tknmax,tkn_500m_month(5,k)),     &
           min(tknmax,tkn_500m_month(6,k)),     &
           min(tknmax,tkn_500m_month(7,k)),     &
           min(tknmax,tkn_500m_month(8,k)),     &
           min(tknmax,tkn_500m_month(9,k)),     &
           min(tknmax,tkn_500m_month(10,k)),        &
           min(tknmax,tkn_500m_month(11,k)),        &
           min(tknmax,tkn_500m_month(12,k))

          write(209,1120) k,TSS_500m_month(1,k),     &
           min(tssmax,TSS_500m_month(2,k)),     &
           min(tssmax,TSS_500m_month(3,k)),     &
           min(tssmax,TSS_500m_month(4,k)),     &
           min(tssmax,TSS_500m_month(5,k)),     &
           min(tssmax,TSS_500m_month(6,k)),     &
           min(tssmax,TSS_500m_month(7,k)),     &
           min(tssmax,TSS_500m_month(8,k)),     &
           min(tssmax,TSS_500m_month(9,k)),     &
           min(tssmax,TSS_500m_month(10,k)),        &
           min(tssmax,TSS_500m_month(11,k)),        &
           min(tssmax,TSS_500m_month(12,k))
      enddo


1115  format(A,16(',',A))
1116  format(I0,16(',',F0.4))
1117  format(A,19(',',A))
1118  format(I0,19(',',F0.4))
1119  format(A,12(',',A))
1120  format(I0,12(',',F0.4))

      return
      end
