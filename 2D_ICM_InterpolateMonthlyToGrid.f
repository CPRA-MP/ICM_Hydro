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
!> @param         IDW_comp_input12(N)                     local storage for incoming compartment monthly array of data to be used
!> @param         IDW_link_input12(M)                     local storage for incoming link monthly array data to be used 
!> @param         IDW_denom(n_cells)                      local storage for IDW denominator
!> @param         IDW_exp                                 exponent to be used in IDW equation
!> @param         IDW_numer(n_cells)                      local storage for IDW numerator
!> @param         IDW_output(n_cells)                     local storage for IDW output values
!> @param         grid_lookup(n_cells,16)                 local storage for grid-to-compartment-to-link lookup table
!> @param         grid_interp_dist(n_cells,16)            local storage for distance lookup table
      
      subroutine ICM_InterpolateMonthlyToGrid(output_flag,gridsize)
      
      use params      

      implicit none
      integer :: jj,k,kk,mmm,output_flag,gridsize,n_cells
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
          
          ! these array dimensions do not change
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
  
!>> filter through output_flag values 101-112, 201-212, 301-312, 401-412 
!>> choose corresponding month dimension of monthly output arrays       
      do mmm=101,112
          if(output_flag == mmm) then
              write(1,5555)' Interpolating monthly mean salinity',
     &                       ' values to grid. Month: ',mmm-100
              
              write(*,5555)' Interpolating monthly mean salinity',
     &                       ' values to grid. Month: ',mmm-100
              do kk=1,N
                  IDW_comp_input(kk) = sal_month_ave(mmm-100,kk)
              enddo
              do jj=1,M
                  IDW_link_input(jj) = sal_month_ave_links(mmm-100,jj)
              enddo
          endif
      enddo
          
      do mmm=201,212
          if(output_flag == mmm) then
              write(1,5555)' Interpolating monthly mean',
     &                    ' temperature values to grid. Month: ',mmm-200
                  
              write(*,5555)' Interpolating monthly mean',
     &                    ' temperature values to grid. Month: ',mmm-200
     
              do kk=1,N
                  IDW_comp_input(kk) = tmp_month_ave(mmm-200,kk)
              enddo
              do jj=1,M
                  IDW_link_input(jj) = tmp_month_ave_links(mmm-200,jj)
              enddo
          endif
      enddo

!!! Algae and TSS not interpolated to grid - they are mapped in the MapMonthlyToGrid subroutine      
!!!      do mmm=301,312
!!!         if(output_flag == mmm) then
!!!               write(*,*)'Interpolating monthly mean Chlorophyll A',
!!!  &                       'values to grid. Month: ',mmm-300
!!!             do kk=1,N
!!!                 IDW_comp_input(kk) = alg_month_ave(mmm-300,kk)
!!!             enddo
!!!             do jj=1,M
!!!                 IDW_link_input(jj) = alg_month_ave_links(mmm-300,jj)
!!!             enddo
!!!         endif
!!!     enddo
!!!
!!!      do mmm=401,412
!!!         if(output_flag == mmm) then
!!!                  write(*,*)'Interpolating monthly mean TSS',
!!!  &                       'values to grid. Month: ',mmm-400
!!!             do kk=1,N
!!!                 IDW_comp_input(kk) = tss_month_ave(mmm-400,kk)
!!!             enddo
!!!             do jj=1,M
!!!                 IDW_link_input(jj) = tss_month_ave_links(mmm-400,jj)
!!!             enddo
!!!         endif
!!!     enddo
!!! Algae and TSS not interpolated to grid - they are mapped in the MapMonthlyToGrid subroutine      
          
          
      IDW_exp = 1.                                            ! exponent to use in inverse distance weighted interpolation 
!>> First DO is looped over the number of grid cells that will have an interpolated data value calculated
      do k=1,n_cells
          
!>> -- Start Inverse Distance Weighting interpolation
          
          IDW_denom(k) = 0.
          IDW_numer(k) = 0.                                   ! reset to zero before starting DO loop
      
!>> -- Second DO is looped over the 15 potential data points at each grid cell (1 compartment, 14 potential links)

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
                              if (linkt(grid_lookup(k,kk)) /= 9) then  !skips links that are overland ridge/levee links (surge only)
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

!>> Set output arrays equal to the temporary interpolation array      
          do mmm=101,112
              if (output_flag == mmm) then
                  sal_IDW_500m_month(mmm-100,k) = IDW_output(k)
              endif
          enddo 
      
          do mmm=201,212
              if (output_flag == mmm) then
                  tmp_IDW_500m_month(mmm-200,k) = IDW_output(k)
              endif
          enddo
!!! Algae and TSS not interpolated to grid - they are mapped in the MapMonthlyToGrid subroutine      
!!!          do mmm=301,312
!!!              if (output_flag == mmm) then
!!!                  Alg_IDW_500m_month(mmm-300,k) = IDW_output(k)
!!!              endif
!!!          enddo
!!!      
!!!          do mmm=401,412
!!!              if (output_flag == mmm) then
!!!                  TSS_IDW_500m_month(mmm-400,k) = IDW_output(k)
!!!              endif
!!!          enddo
!!! Algae and TSS not interpolated to grid - they are mapped in the MapMonthlyToGrid subroutine      
         
      enddo
!>> End first DO loop
      




!>> Deallocate temporary interpolation arrays      
      deallocate(grid_interp_dist)
      deallocate(grid_lookup)
      deallocate(IDW_denom)
      deallocate(IDW_numer)
      deallocate(IDW_output)
      deallocate(IDW_comp_input)
      deallocate(IDW_link_input)

5555  format(A,A,I3)      
      
      return
      end