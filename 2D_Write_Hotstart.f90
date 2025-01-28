!> @file
!> @brief This is the subroutine to write hotstart values.
!> @details This is the subroutine to save the last timestep compartment values to a hotstart output file
!> which will be used for the next year's initial conditions

      subroutine Write_Hotstart

      use params

      implicit none
      integer :: j       !iterators
      
!>> Add error terms to last timestep value before writing hotstart file
      write(1,*) 'Adding error adjustment to last timestep values.'
      write(*,*) 'Adding error adjustment to last timestep values.'

      do j=1,N
          Es(j,2) = Es(j,2) + stage_error

          if(S(j,2) < 1.0) then
              S(j,2) = max(0.0,S(j,2) + sal_0_1_error)
          elseif(S(j,2) < 5.0) then
              S(j,2) = max(0.0,S(j,2) + sal_1_5_error)
          elseif(S(j,2) < 20.0) then
              S(j,2) = max(0.0,S(j,2) + sal_5_20_error)
          else
              S(j,2) = max(0.0,S(j,2) + sal_20_35_error)
          endif

          Css(j,2,1) = max(0.0,Css(j,2,1)*(1.0+tss_error))
          Css(j,2,2) = max(0.0,Css(j,2,2)*(1.0+tss_error))
          Css(j,2,3) = max(0.0,Css(j,2,3)*(1.0+tss_error))
          Css(j,2,4) = max(0.0,Css(j,2,4)*(1.0+tss_error))

      enddo

!>> Save output values of last simulation timestep in a hotstart file
      write(1,*)'Saving values from last timestep to a hotstart file.'
      write(*,*)'Saving values from last timestep to a hotstart file.'

      write(401,11140)'Compartment',	&
                      'Stage',		    &
                      'Salinity',		&
                      'CSS_sand',		&
                      'CSS_silt',		&
                      'CSS_clay',		&
                      'CSS_clayfloc',	&
                      'WaterTemp',		&
                      'NO3_NO2',		&
                      'NH4',		    &
                      'DIN',		    &
                      'OrgN',		    &
                      'TIP',		    &
                      'TOC',		    &
                      'DO',		        &
                      'LiveAlg',	    &
                      'DeadAlg',	    &
                      'DON',		    &
                      'DOP',		    &
                      'DIP',		    &
                      'ChlA',		    &
                      'TKN',            &
                      'Marsh_stage'

      do j=1,N
	      write(401,11142) j,		                    &
                          Es(j,2),		                &
                          S(j,2),		                &
                          Css(j,2,1),		            &
                          Css(j,2,2),		            &
                          Css(j,2,3),		            &
                          Css(j,2,4),		            &
                          Tempw(j,2),		            &
                          min(Chem(j,1,2),1000.0),		&
                          min(Chem(j,2,2),1000.0),		&
                          min(Chem(j,3,2),1000.0),		&
                          min(Chem(j,4,2),1000.0),		&
                          min(Chem(j,5,2),1000.0),		&
                          min(Chem(j,6,2),1000.0),		&
                          min(Chem(j,7,2),1000.0),		&
                          min(Chem(j,8,2),1000.0),		&
                          min(Chem(j,9,2),1000.0),		&
                          min(Chem(j,10,2),1000.0),		&
                          min(Chem(j,11,2),1000.0),		&
                          min(Chem(j,12,2),1000.0),		&
                          min(Chem(j,13,2),1000.0),		&
                          min(Chem(j,14,2),1000.0),     & ! chem unit = mg/L
                          Eh(j,2)
      enddo
      close(401)

11140 format(A,22(' ',A))
11142 format(I8,22(F20.4))

      return
      end
