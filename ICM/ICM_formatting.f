!> @file
!> @brief This subroutine formats the summarized data for use in the other ICM routines.
!> @details The data that was summarized at each grid cell for use in the ICM is formatted here
!> so that it can be directly passed to the non-Fortran based portions of the ICM.
!> Specifically, this subroutine prepares text files that are formatted as required for
!> direct ingestion into the Vegetation model.
      
!> @author Eric White - The Water Institute of the Gulf      

!> @param[in]     depth_summer_500m(grid_cell)      
!> @param[in]     salinity_IDW_500m(grid_cell)                            Interpolated salinity values for each 500-m grid cell
!> @param[in]     salinity_summer_IDW_500m(grid_cell)                     Interpolated summertime salinity values for each 500-m grid cell
!> @param[in]     stage_var_summer_500m(grid_cell)                        Variance in summertime stage mapped to each 500-m grid cell
!> @param[in]     tmp_summer_IDW_500m(grid_cell)                          Interpolated summertime water temperature for each 500-m grid cell
!> @param[in]     veg_matrix_cols                                         Number of columns in the matrix structure of the Vegetation I/O files
!> @param[in]     veg_matrix_rows                                         Number of rows in the matrix structure of the Vegetation I/O files      

!> @param[out]    depth_summer_forVeg(j,k)                                Mean summertime depth formatted for Vegetation ICM routine
!> @param[out]    salinity_forVeg(veg_matrix_cols,veg_matrix_rows)        Mean salinity formatted for Vegetation ICM routine
!> @param[out]    salinity_summer_forVeg(veg_matrix_cols,veg_matrix_rows) Mean summer salinity formatted for Vegetation ICM routine
!> @param[out]    stage_var_forVeg(j,k)                                   Variance in summertime stage formatted for Vegetation ICM routine
!> @param[out]    tmp_summer_forVeg(j,k)                                  Mean summertime water temperature formatted for Vegetation ICM routine

!> @param         grid_cell                                               temporary storage of 500-m grid cell ID number used to lookup respective value
!> @param         j                                                       DO loop counter
!> @param         k                                                       DO loop counter
!> @param         veg_grid_IDs(veg_matrix_cols,veg_matrix_rows)           matrix of 500-m grid cell ID numbers
      
      subroutine ICM_Formatting
			
	use params

      implicit none
      integer :: j,k,grid_cell
      
      write(1,*)' Preparing input files for Vegetation model.'
      write(*,*)' Preparing input files for Vegetation model.'

!>@par General Structure of Subroutine Logic:      

!>> Read in matrix of Grid IDs for 500x500 m cells used by Vegetation ICM routine.
!>> Structure of matrix matches the grid matrix in the Vegetation I/O files.
      read(203,*) veg_grid_IDs

    
! Put data for each grid cell into matrix format required for Vegetation ICM routine
!>> Loop over columns and rows of grid cell IDs matrix.
      do j=1,veg_matrix_cols
          do k=1,veg_matrix_rows
!>> -- Lookup grid cell ID from input matrix.          
              grid_cell = veg_grid_IDs(j,k)
!>> -- Find grid cell's respective data value from the Fortran arrays generated in 'ICM_summaries' & 'ICM_InterpolateToGrid' subroutines.              
              ! Flag 'no data' values in input grid matrix as 'no data' in output array (e.g. = -9999)
              if (grid_cell == -9999) then
                  salinity_forVeg(j,k) = -9999.
                  salinity_summer_forVeg(j,k) = -9999.
                  tmp_summer_forVeg(j,k) = -9999.
                  stage_var_forVeg(j,k) = -9999.
                  depth_summer_forVeg(j,k) = -9999.
                  tree_est_forVeg(j,k) = -9999.
                  ht_abv_water_forVeg(j,k) = -9999.
                  per_land_forVeg(j,k) = -9999.
              else    
!>> -- Low-pass filter on mean annual salinity for Veg model input files
                  salinity_forVeg(j,k) = min(salmax,
     &                         salinity_IDW_500m(grid_cell))
!>> -- Low-pass filter on mean summer salinity for Veg model input files
                  salinity_summer_forVeg(j,k) = min(salmax,
     &                            salinity_summer_IDW_500m(grid_cell))
!>> -- Low-pass filter on mean summer temperature for Veg model input files
                  tmp_summer_forVeg(j,k)= min(tmpmax,
     &                            tmp_summer_IDW_500m(grid_cell))
!>> -- High-pass & Low-pass filter on summer water level standard deviation for Veg model input files (ICM_summaries calculates WSEL variance, convert to st dev here)
                  stage_var_forVeg(j,k)=max(-rangemax,min(rangemax,
     &                           stage_var_summer_500m(grid_cell)**0.5))
!>> -- High-pass & Low-pass filter on summer water depth for Veg model input files 
                  depth_summer_forVeg(j,k)= max(-depthmax,min(depthmax,
     &                            depth_summer_500m(grid_cell)))
!>> -- High-pass & Low-pass filter on height above mean water for Veg model input files 
                  ht_abv_water_forVeg(j,k) = max(-stagemax,min(stagemax,
     &                            height_500m(grid_cell)))
                  tree_est_forVeg(j,k) = tree_est(grid_cell)
!>> -- High-pass & low-pass filter on percent land to avoid any near zero/hundred errors, but keep -9999 NoData values
                  if (per_land_500m(grid_cell) == -9999.) then
                      per_land_forVeg(j,k) = per_land_500m(grid_cell)
                  else
                      per_land_forVeg(j,k) = max(0.0,min(100.0,
     &                                    per_land_500m(grid_cell)))
                  endif
              endif
          enddo
      enddo
!>> End loop of grid cell IDs matrix.

!>> Write these new, matrix-formatted arrays to file to be read directly into Vegetation model.
!>> THESE HEADERS ARE USED BY OTHER ICM ROUTINES - DO NOT CHANGE WITHOUT UPDATING OTHER ICM ROUTINES   
      write(301,1117)'# Year =',year
      write(301,1117)'NROWS',veg_matrix_rows
      write(301,1117)'NCOLS',veg_matrix_cols
      write(301,1117)'XLLCORNER',veg_xllcorner
      write(301,1117)'YLLCORNER',veg_yllcorner
      write(301,1117)'CELLSIZE 500'
      write(301,1117)'NODATA_VALUE -9999.00'
      write(301,1118) salinity_forVeg
      
      write(302,1117)'# Year =',year
      write(302,1117)'NROWS',veg_matrix_rows
      write(302,1117)'NCOLS',veg_matrix_cols
      write(302,1117)'XLLCORNER',veg_xllcorner
      write(302,1117)'YLLCORNER',veg_yllcorner
      write(302,1117)'CELLSIZE 500'
      write(302,1117)'NODATA_VALUE -9999.00'
      write(302,1118) salinity_summer_forVeg       
          
      write(303,1117)'# Year =',year
      write(303,1117)'NROWS',veg_matrix_rows
      write(303,1117)'NCOLS',veg_matrix_cols
      write(303,1117)'XLLCORNER',veg_xllcorner
      write(303,1117)'YLLCORNER',veg_yllcorner
      write(303,1117)'CELLSIZE 500'
      write(303,1117)'NODATA_VALUE -9999.00'
      write(303,1118) depth_summer_forVeg  
      
      write(304,1117)'# Year =',year
      write(304,1117)'NROWS',veg_matrix_rows
      write(304,1117)'NCOLS',veg_matrix_cols
      write(304,1117)'XLLCORNER',veg_xllcorner
      write(304,1117)'YLLCORNER',veg_yllcorner
      write(304,1117)'CELLSIZE 500'
      write(304,1117)'NODATA_VALUE -9999.00'
      write(304,1118) stage_var_forVeg
      
      write(305,1117)'# Year =',year
      write(305,1117)'NROWS',veg_matrix_rows
      write(305,1117)'NCOLS',veg_matrix_cols
      write(305,1117)'XLLCORNER',veg_xllcorner
      write(305,1117)'YLLCORNER',veg_yllcorner
      write(305,1117)'CELLSIZE 500'
      write(305,1117)'NODATA_VALUE -9999.00'
      write(305,1118) tmp_summer_forVeg  

      write(306,1117)'# Year =',year
      write(306,1117)'NROWS',veg_matrix_rows
      write(306,1117)'NCOLS',veg_matrix_cols
      write(306,1117)'XLLCORNER',veg_xllcorner
      write(306,1117)'YLLCORNER',veg_yllcorner
      write(306,1117)'CELLSIZE 500'
      write(306,1117)'NODATA_VALUE -9999.00'
      write(306,1119) tree_est_forVeg

      write(307,1117)'# Year =',year
      write(307,1117)'NROWS',veg_matrix_rows
      write(307,1117)'NCOLS',veg_matrix_cols
      write(307,1117)'XLLCORNER',veg_xllcorner
      write(307,1117)'YLLCORNER',veg_yllcorner
      write(307,1117)'CELLSIZE 500'
      write(307,1117)'NODATA_VALUE -9999.00'
      write(307,1118) ht_abv_water_forVeg
      
      write(308,1117)'# Year =',year
      write(308,1117)'NROWS',veg_matrix_rows
      write(308,1117)'NCOLS',veg_matrix_cols
      write(308,1117)'XLLCORNER',veg_xllcorner
      write(308,1117)'YLLCORNER',veg_yllcorner
      write(308,1117)'CELLSIZE 500'
      write(308,1117)'NODATA_VALUE -9999.00'
      write(308,1118) per_land_forVeg
      
 1117 format(A,1x,I0)
 1118 format(<veg_matrix_cols-1>(F0.4,1x),F0.4)
 1119 format(<veg_matrix_cols-1>(I0,1x),I0)
      
      
      return
      end 
			
			