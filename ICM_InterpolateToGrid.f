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
!> @param[in]     grid_interp_dist_500m(n_500m_cells,16)  lookup table for distances between grid cell and compartments and links
!> @param[in]     grid_lookup_500m(n_500m_cells,16)       lookup table for 500 m grid - matches compartment and link to grid cell
!> @param[in]     gridsize                                size of grid cell output (30m or 500m)
!> @param[in]     n_500m_cells                            number of 500 m grid cells
!> @param[in]     output_flag                             flag for identifying data value to be interpolated
!> @param[in]     sal_ave(N)                              mean salinity for compartments
!> @param[in]     sal_ave_links(M)                        mean salinity for links
!> @param[in]     sal_ave_summmer(N)                      mean summertime salinity for compartments
!> @param[in]     sal_ave_summer_links(M)                 mean summertime salinity for links
!> @param[in]     tmp_ave(N)                              mean water temperature for compartments
!> @param[in]     tmp_ave_links(M)                        mean water temperature for links
!> @param[in]     tmp_ave_summmer(N)                      mean summertime temperature for compartments
!> @param[in]     tmp_ave_summer_links(M)                 mean summertime temperature for links

     
!> @param[out]    salinity_IDW_500m(n_500m_cells)         global array for interpolated compartment and link salinities for 500m grid
!> @param[out]    salinity_summer_IDW_500m(n_500m_cells)  global array for interpolated compartment and link summertime salinities for 500m grid
!> @param[out]    tmp_IDW_500m(n_500m_cells)              global array for interpolated compartment and link water temp for 500m grid
!> @param[out]    tmp_summer_IDW_500m(n_500m_cells)       global array for interpolated compartment and link summertime water temp  for 500m grid
     
!> @param         dist                                    distance correction factor to avoid DIV/0 errors
!> @param         k                                       DO loop counter
!> @param         kk                                      DO loop counter
!> @param         n_cells                                 local storage for number of grid cells
!> @param         IDW_comp_input(N)                       local storage for incoming compartment data to be used
!> @param         IDW_link_input(M)                       local storage for incoming link data to be used 
!> @param         IDW_denom(n_cells)                      local storage for IDW denominator
!> @param         IDW_exp                                 exponent to be used in IDW equation
!> @param         IDW_numer(n_cells)                      local storage for IDW numerator
!> @param         IDW_output(n_cells)                     local storage for IDW output values
!> @param         grid_lookup(n_cells,16)                 local storage for grid-to-compartment-to-link lookup table
!> @param         grid_interp_dist(n_cells,16)            local storage for distance lookup table
      
      subroutine ICM_InterpolateToGrid(output_flag,gridsize)
      
      use params      

      implicit none
      integer :: k,kk,output_flag,gridsize,n_cells
      real :: dist
      real(dp),dimension(:),allocatable :: IDW_comp_input
      real(dp),dimension(:),allocatable :: IDW_link_input
      real(dp),dimension(:),allocatable :: IDW_denom
      real(dp),dimension(:),allocatable :: IDW_numer
      real(dp),dimension(:),allocatable :: IDW_output
      integer, dimension(:,:), allocatable :: grid_lookup
      real(sp), dimension(:,:), allocatable :: grid_interp_dist
    
!>@par General Structure of Subroutine Logic:

!>> Allocate temporary IDW arrays to be of length equal to number of grid cells - these are deallocated at end of this subroutine      
      if(gridsize == 500) then
          
          n_cells = n_500m_cells
          
          ! these array sizes change depending on interpolated grid cell size
          allocate(grid_interp_dist(n_cells,22))
          allocate(grid_lookup(n_cells,22))
          allocate(IDW_denom(n_cells))
          allocate(IDW_numer(n_cells))
          allocate(IDW_output(n_cells))
          
          ! these array sizes do not change - however they are allocated here to ensure
          ! that the interpolation scheme starts with a fresh array
          allocate(IDW_comp_input(N))             
          allocate(IDW_link_input(M))
          
!>> Error checks - print error message if interpolation arrays were not read into subroutine correctly
          if(size(grid_lookup)==size(grid_lookup_500m)) then
              grid_lookup = grid_lookup_500m
          else
              write(1,*)'***************ERROR**********************'
              write(1,*)' Unequal array dimensions for grid lookups.'
              write(1,*)'***************ERROR**********************'
              
              write(*,*)'***************ERROR**********************'
              write(*,*)' Unequal array dimensions for grid lookups.'
              write(*,*)'***************ERROR**********************'
          endif
          
          if(size(grid_interp_dist)==size(grid_interp_dist_500m)) then
              grid_interp_dist = grid_interp_dist_500m
          else
              write(1,*)'***************ERROR**********************'
              write(1,*)' Unequal array dimensions for grid lookups.'
              write(1,*)'***************ERROR**********************'
          
              write(*,*)'***************ERROR**********************'
              write(*,*)' Unequal array dimensions for grid lookups.'
              write(*,*)'***************ERROR**********************'
          endif
      
      else
          write(1,*)'*********ERROR**********************'
          write(1,*)' Interpolation grid size is not 500 m!'
          write(1,*)' Only 500 m grid cells are currently supported.'
          
          write(*,*)'*********ERROR**********************'
          write(*,*)' Interpolation grid size is not 500 m!'
          write(*,*)' Only 500 m grid cells are currently supported.'
      
      endif
 
      
!>> Assign model output data to be interpolated to the temporary IDW arrays based on the 'output_flag' value passed to this subroutine.
!>> (e.g. output_flag 1 = mean salinity data, output_flag 2 = mean summer salinity, etc.)
      if(output_flag == 1) then
          write(1,*)' Interpolating mean annual salinity to grid.'
          write(*,*)' Interpolating mean annual salinity to grid.'
          IDW_comp_input = sal_ave
          IDW_link_input = sal_ave_links
      elseif(output_flag == 2) then
          write(1,*)' Interpolating mean summer salinity to grid.'
          write(*,*)' Interpolating mean summer salinity to grid.'
          IDW_comp_input = sal_ave_summer
          IDW_link_input = sal_ave_summer_links
      elseif(output_flag == 3) then
          write(1,*)' Interpolating mean annual temperature to grid.'
          write(*,*)' Interpolating mean annual temperature to grid.'
          IDW_comp_input = tmp_ave
          IDW_link_input = tmp_ave_links
      elseif(output_flag == 4) then
          write(1,*)' Interpolating mean summer temperature to grid.'
          write(*,*)' Interpolating mean summer temperature to grid.'
          IDW_comp_input = tmp_ave_summer
          IDW_link_input = tmp_ave_summer_links 
      elseif(output_flag == 5) then
          write(1,*)' Interpolating 2-wk max summer salinity to grid.'
          write(*,*)' Interpolating 2-wk max summer salinity to grid.'
          IDW_comp_input = sal_2wk_ave_max
          IDW_link_input = sal_2wk_ave_max_links 
      endif
      IDW_exp = 1.                                            ! exponent to use in inverse distance weighted interpolation 
!>> First DO is looped over the number of grid cells that will have an interpolated data value calculated
      do k=1,n_cells
          
!>> -- Start Inverse Distance Weighting interpolation
          
          IDW_denom(k) = 0.
          IDW_numer(k) = 0.                                   ! reset to zero before starting DO loop
      
!>> -- Second DO is looped over the 15 potential data points at each grid cell (1 compartment, 14 potential links) - only links that are active and not overland marsh or ridge links are used for interpolation.
          do kk=2,16
              if (kk==2) then                                         ! lookup compartment value - 
                  dist = max(0.0001,grid_interp_dist(k,kk))           ! if distance is zero, set to near-zero value to prevent and div-by-zero NaNs
!>> - Calculate IDW numerator - function of distances and model output values
                  IDW_numer(k) = IDW_numer(k) + 
     &                IDW_comp_input(grid_lookup(k,kk))/(dist**IDW_exp)
!>> - Calculate IDW denominator - only a function of distances
                  IDW_denom(k) = IDW_denom(k) + 1./(dist**IDW_exp)
              else                                                    ! lookup link salinity values
                  if (grid_interp_dist(k,kk) >= 0) then               !skips grid cells that have -9999 distances
                      dist = max(0.0001,grid_interp_dist(k,kk))
                  endif
                  if (grid_lookup(k,kk) >= 0) then                    !skips grid cells that have -9999 distances
                      if (linkt(grid_lookup(k,kk))  > 0) then         !skips links that are inactive due to negative link type flag   
                          if (linkt(grid_lookup(k,kk)) /= 8) then     !skips links that are overland marsh links (surge only)
                              if (linkt(grid_lookup(k,kk)) /= 9) then  !skips linkst hat are overland ridge/levee links (surge only)
                                  IDW_numer(k) = IDW_numer(k) + 
     &                                IDW_link_input(grid_lookup(k,kk))/
     &                                (dist**IDW_exp)
                                  IDW_denom(k) = IDW_denom(k) + 
     &                                 1./(dist**IDW_exp)
                              endif
                          endif
                      endif
                  endif
              endif
          enddo
!>> -- End second DO loop
      
!>> -- IDW value calculated from 'numer' and 'denom' for grid cell, k. 'numer' and 'denom' will be reset to zero for next grid cell, k+1.
         IDW_output(k) = IDW_numer(k)/IDW_denom(k)
      
      enddo
!>> End first DO loop
      
!>> Set output arrays equal to the temporary interpolation array
      if (output_flag == 1) then
          salinity_IDW_500m = IDW_output
      elseif (output_flag == 2) then
          salinity_summer_IDW_500m = IDW_output
      elseif (output_flag == 3) then
          tmp_IDW_500m = IDW_output
      elseif (output_flag == 4) then
          tmp_summer_IDW_500m = IDW_output
      elseif (output_flag == 5) then
          sal_thresh_IDW_500m = IDW_output
      endif
      
!>> Deallocate temporary interpolation arrays      
      deallocate(grid_interp_dist)
      deallocate(grid_lookup)
      deallocate(IDW_denom)
      deallocate(IDW_numer)
      deallocate(IDW_output)
      deallocate(IDW_comp_input)
      deallocate(IDW_link_input)
      
      return
      end