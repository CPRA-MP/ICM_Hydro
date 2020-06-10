!> @file
!> @brief This subroutine interpolates model output data to the pre-defined grid structure.
!> @details This subroutine uses an inverse distance weighting (IDW) method to interpolate model output 
!> as summarized for use in the ICM (e.g. mean salinity) to the grid used by other ICM routines.
!> The model output used for interpolation comes from each hydro compartment and all primary hydraulic links in each compartment.
!> Look-up tables that link each grid cell to the appropriate hydro compartments/links and the distances between all respective centroids are read into this subroutine.
!> This subroutine uses dynamically allocated arrays that are deallocated at the end of the subroutine.
!> This subroutine is called for each dataset that requires interpolation; the 'output_flag' and 'gridsize' parameters are passed into this subroutine.
!> These two input parameters control the array size needed for interpolation and the correct assignment from the temporary interpolation array into the permanent storage array that is eventually used to write output to a file.

!> @author  Eric White - The Water Institute of the Gulf
      
!> @param[in]     N                                       number of hydro compartments in model, N
!> @param[in]     M                                       number of links in model, M
!> @param[in]     grid_lookup_500m(n_500m_cells,16)       lookup table for 500 m grid - matches compartment and link to grid cell
!> @param[in]     gridsize                                size of grid cell output (30m or 500m)
!> @param[in]     n_500m_cells                            number of 500 m grid cells
!> @param[in]     output_flag                             flag for identifying data value to be interpolated
!> @parma[in]     pct_sand_bed                            percentage of open water bed sediments that are sand
!> @param[in]     sal_ave(N)                              mean salinity for compartments
!> @param[in]     sal_ave_summmer(N)                      mean summertime salinity for compartments
!> @param[in]     stage_ave(N)                            mean stage for compartments
!> @param[in]     stage_ave_summer(N)                     mean summertime stage for compartments
!> @param[in]     stage_var_summer(N)                     variance in summertime stage for compartments
!> @param[in]     tmp_ave(N)                              mean water temperature for compartments
!> @param[in]     tmp_ave_summer(N)                       mean summertime water temperature for compartments
      

!> @param[out]    pct_sand_bed_500m(n_500m_cells)         global array for compartment-to-grid overlay percent sand in bed sediments  for 500m grid
!> @param[out]    salinity_500m(n_500m_cells)             global array for compartment-to-grid overlay salinities for 500m grid
!> @param[out]    salinity_summer_500m(n_500m_cells)      global array for compartment-to-grid overlay salinities for 500m grid
!> @param[out]    stage_500m(n_500m_cells)                global array for compartment-to-grid overlay stage for 500m grid
!> @param[out]    stage_summer_500m(n_500m_cells)         global array for compartment-to-grid overlay summertime stage  for 500m grid
!> @param[out]    stage_var_summer_500m(n_500m_cells)     global array for compartment-to-grid overlay summertime stage variance for 500m grid
!> @param[out]    tmp_500m(n_500m_cells)                  global array for compartment-to-grid overlay water temp for 500m grid
!> @param[out]    tmp_summer_500m(n_500m_cells)           global array for compartment-to-grid overlay summertime water temp  for 500m grid

!> @param         k                                       DO loop counter
!> @param         n_cells                                 local storage for number of grid cells
!> @param         map_comp_input(N)                       local storage for incoming compartment data to be used
!> @param         grid_no_interp(n_cells)                 local storage for compartment-to-grid overlay (no interpolation)
!> @param         grid_lookup(n_cells,16)                 local storage for grid-to-compartment-to-link lookup table
!> @param         grid_interp_dist(n_cells,16)            local storage for distance lookup table
      
      subroutine ICM_MapToGrid(output_flag,gridsize)
      
      use params      

      implicit none
      integer :: k,output_flag,gridsize,n_cells
      real :: dist
      real(dp),dimension(:),allocatable :: map_comp_input
      real(dp),dimension(:),allocatable :: grid_no_interp
      integer, dimension(:,:), allocatable :: grid_lookup

    
!>@par General Structure of Subroutine Logic:

!>> Allocate temporary IDW arrays to be of length equal to number of grid cells - these are deallocated at end of this subroutine      
      if(gridsize == 500) then
          
          n_cells = n_500m_cells
          
          ! these array sizes change depending on grid cell size
          allocate(grid_lookup(n_cells,22))
          allocate(grid_no_interp(n_cells))
          
          ! these array sizes do not change - however they are allocated here to ensure
          ! that the interpolation scheme starts with a fresh array
          allocate(map_comp_input(N))             

          
!>> Error checks - print error message if interpolation arrays were not read into subroutine correctly
          if(size(grid_lookup)==size(grid_lookup_500m)) then
              grid_lookup = grid_lookup_500m
          else
              write(1,*)'***************ERROR**********************'
              write(1,*)'Unequal array dimensions for grid lookups.'
              write(1,*)'***************ERROR**********************'
              
              write(*,*)'***************ERROR**********************'
              write(*,*)'Unequal array dimensions for grid lookups.'
              write(*,*)'***************ERROR**********************'
          endif
      
      else
          write(1,*)'*********ERROR**********************'
          write(1,*)'Mapping grid size is not 500 m!'
          write(1,*)'Only 500 m grid cells are currently supported.'
      
          write(*,*)'*********ERROR**********************'
          write(*,*)'Mapping grid size is not 500 m!'
          write(*,*)'Only 500 m grid cells are currently supported.'
      
      endif
 
      
!>> Assign model output data to be interpolated to the temporary IDW arrays based on the 'output_flag' value passed to this subroutine.
!>> (e.g. output_flag 1 = mean salinity data, output_flag 2 = mean summer salinity, etc.)
      if(output_flag == 1) then
          write(1,*)' Mapping mean annual salinity output to grid.'
          write(*,*)' Mapping mean annual salinity output to grid.'
          map_comp_input = sal_ave

      elseif(output_flag == 2) then
          write(1,*)' Mapping mean summer salinity output to grid.'
          write(*,*)' Mapping mean summer salinity output to grid.'
          map_comp_input = sal_ave_summer
          
      elseif(output_flag == 3) then
          write(1,*)' Mapping mean annual temperature output to grid.'
          write(*,*)' Mapping mean annual temperature output to grid.'
          map_comp_input = tmp_ave
          
      elseif(output_flag == 4) then
          write(1,*)' Mapping mean summer temperature output to grid.'
          write(*,*)' Mapping mean summer temperature output to grid.'
          map_comp_input = tmp_ave_summer
          
      elseif(output_flag == 5) then
          write(1,*)' Mapping mean stage output to grid.'
          write(*,*)' Mapping mean stage output to grid.'
          map_comp_input = stage_ave

      elseif(output_flag == 6) then
          write(1,*)' Mapping mean summer stage output to grid.'
          write(*,*)' Mapping mean summer stage output to grid.'
          map_comp_input = stage_ave_summer

      elseif(output_flag == 7) then
          write(1,*)' Mapping mean summer stage output to grid.'
          write(*,*)' Mapping mean summer stage output to grid.'
          map_comp_input = stage_var_summer
      
      elseif(output_flag == 8) then
          write(1,*)' Mapping max 2-wk summer salinity output to grid.'
          write(*,*)' Mapping max 2-wk summer salinity output to grid.'
          map_comp_input = sal_2wk_ave_max
                    
      elseif(output_flag == 9) then
          write(1,*)' Mapping percent sand in OW bed output to grid.'
          write(*,*)' Mapping percent sand in OW bed output to grid.'
          map_comp_input = pct_sand_bed
      endif
      

                                          
!>> First DO is looped over the number of grid cells that will have an interpolated data value calculated
      do k=1,n_cells
          
!>> -- Map compartment output values to grid cells - no interpolation (simple, non-weighted overlay)
          grid_no_interp(k) = map_comp_input(grid_lookup(k,2))
       
      enddo
!>> End first DO loop
      
!>> Set output arrays equal to the temporary interpolation array
      if (output_flag == 1) then
          salinity_500m = grid_no_interp
      elseif (output_flag == 2) then
          salinity_summer_500m = grid_no_interp
      elseif (output_flag == 3) then
          tmp_500m = grid_no_interp
      elseif (output_flag == 4) then
          tmp_summer_500m = grid_no_interp
      elseif (output_flag == 5) then
          stage_500m = grid_no_interp
      elseif (output_flag == 6) then
          stage_summer_500m = grid_no_interp
      elseif (output_flag == 7) then
          stage_var_summer_500m = grid_no_interp
      elseif (output_flag == 8) then
          sal_thresh_500m = grid_no_interp
      elseif (output_flag == 9) then
         pct_sand_bed_500m = grid_no_interp
      endif
      
!>> Deallocate temporary interpolation arrays      
      deallocate(grid_lookup)
      deallocate(grid_no_interp)
      deallocate(map_comp_input)

      return
      end