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
!> @param[in]     sal_ave(N)                              mean salinity for compartments
!> @param[in]     sal_ave_summmer(N)                      mean summertime salinity for compartments
!> @param[in]     stage_ave(N)                            mean stage for compartments
!> @param[in]     stage_ave_summer(N)                     mean summertime stage for compartments
!> @param[in]     stage_var_summer(N)                     variance in summertime stage for compartments
!> @param[in]     tmp_ave(N)                              mean water temperature for compartments
!> @param[in]     tmp_ave_summer(N)                       mean summertime water temperature for compartments
      
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
      
      subroutine ICM_MapMonthlyToGrid(output_flag,gridsize)
      
      use params      

      implicit none
      integer :: k,kk,mmm,output_flag,gridsize,n_cells
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
              write(1,*)' Unequal array dimensions for grid lookups.'
              write(1,*)'***************ERROR**********************'
          
              write(*,*)'***************ERROR**********************'
              write(*,*)' Unequal array dimensions for grid lookups.'
              write(*,*)'***************ERROR**********************'
          endif
      
      else
          write(1,*)'*********ERROR**********************'
          write(1,*)' Mapping grid size is not 500 m!'
          write(1,*)' Only 500 m grid cells are currently supported.'
                
          write(*,*)'*********ERROR**********************'
          write(*,*)' Mapping grid size is not 500 m!'
          write(*,*)' Only 500 m grid cells are currently supported.'
      
      endif
 
      
!>> Assign model output data to be interpolated to the temporary IDW arrays based on the 'output_flag' value passed to this subroutine.
!>> (e.g. output_flag 1 = mean salinity data, output_flag 2 = mean summer salinity, etc.)

!!! salinity and temperature are mapped to grid using IDW interpolation subroutine      
!!!   do mmm=101,112
!!!         if(output_flag == mmm) then
!!!               write(*,5555)' Mapping monthly mean salinity',
!!!  &                      ' values to grid. Month: ',mmm-100
!!!	   	    do kk=1,N
!!!	   		    map_comp_input(kk) = sal_month_ave(mmm-100,kk)
!!!	   	    enddo
!!!       endif
!!!	enddo
!!!
!!!	do mmm=201,212
!!!          if(output_flag == mmm) then
!!!               write(*,5555)' Mapping monthly mean temperature',
!!!  &                      ' values to grid. Month: ',mmm-200
!!!	        do kk=1,N
!!!	     		map_comp_input(kk) = tmp_month_ave(mmm-200,kk)
!!!	    	enddo
!!!          endif
!!!      enddo
!!! salinity and temperature are mapped to grid using IDW interpolation subroutine
      
   	do mmm=301,312
	    if(output_flag == mmm) then
               write(1,5555)' Mapping monthly mean TKN',
     &                      ' values to grid. Month: ',mmm-300
               write(*,5555)' Mapping monthly mean TKN',
     &                      ' values to grid. Month: ',mmm-300
              do kk=1,N
	     	    map_comp_input(kk) = tkn_month_ave(mmm-300,kk)
              enddo
          endif
      enddo

      do mmm=401,412
	   	if(output_flag == mmm) then
              write(1,5555)' Mapping monthly mean TSS',
     &                      ' values to grid. Month: ',mmm-400
              write(*,5555)' Mapping monthly mean TSS',
     &                      ' values to grid. Month: ',mmm-400
              do kk=1,N
	     		map_comp_input(kk) = tss_month_ave(mmm-400,kk)
	    	enddo
          endif
      enddo

                                          
!>> First DO is looped over the number of grid cells that will have an interpolated data value calculated
      do k=1,n_cells
          
!>> -- Map compartment output values to grid cells - no interpolation (simple, non-weighted overlay)
          grid_no_interp(k) = map_comp_input(grid_lookup(k,2))
       
      
!>> End first DO loop
      
!>> Set output arrays equal to the temporary interpolation array
!!! salinity and temperature are mapped to grid using IDW interpolation subroutine
!!!      do mmm=101,112
!!!    	    if (output_flag == mmm) then
!!!               salinity_500m_month(mmm-100,k) = grid_no_interp(k)
!!!           endif
!!!      endo	
!!!      
!!!       do mmm=201,212
!!!     	    if (output_flag == mmm) then
!!!               tmp_500m_month(mmm-200,k) = grid_no_interp(k)
!!!           endif
!!!       enddo
!!! salinity and temperature are mapped to grid using IDW interpolation subroutine      

          do mmm=301,312
      	    if (output_flag == mmm) then
                  tkn_500m_month(mmm-300,k) = grid_no_interp(k)
              endif
          enddo
      
          do mmm=401,412
      	    if (output_flag == mmm) then
                  TSS_500m_month(mmm-400,k) = grid_no_interp(k)
       	    endif
          enddo
      
      enddo
      
!>> Deallocate temporary interpolation arrays      
      deallocate(grid_lookup)
      deallocate(grid_no_interp)
      deallocate(map_comp_input)

5555  format(A,A,I3)

      return
      end