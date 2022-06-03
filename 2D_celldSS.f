!       Subroutine CellDSS(Dz,j,CSSTRIBj,dref,Tres)

! QSSUM, QSSumh, kthr and kday now global parameters - no longer needed to be passed into subroutine      
  	Subroutine CelldSS(j,kday,kthr,CSSTRIBj,dref,Tres)
cJAM      Tributary and resuspension/deposition contributions to SS

!> @param     QSsum(k)        sediment flux in Open Water from all links - negative value is flux INTO open water (g/s)
!> @param     QSsumh(k)       sediment flux in Marsh from all marsh links - negative value is flux INTO marsh (g/s)
!> @param     MEEsedRate(k)   sediment flux from erosion of marsh face - negative value is flux INTO open water (g/s)
!> @param     QSmarsh(k)      sediment flux via exchange flow between marsh and open water (Kadlec-Knight flow) - negative value is flux FROM marsh INTO open water (g/s)
!> @param     QStrib(k)       sediment flux due to tributary flows - negative value is flux INTO compartment from tributaries (g/s)
!> @param     QSdiv(k)        sediment flux due to diversion flows - negative value is flux INTO compartment from diversions (g/s)
!> @param     SedAccumRate(k) sediment flux due to depostion/resuspension - positive value is flux FROM water column ONTO bed (g/s)
!> @param     sed_avail_h     sediment suspended in marsh availabe for deposition (g)
!> @param     depo_avail_h    sediment flux available to deposit in marsh during timestep based on available sediment (g/s)
!> @param     depo_settling_h sediment flux available to desposit in marsh based on settling velocity (g/s)
      
	use params      
      
      
      real :: ddy_1,ddy_2,ddh_1,ddh_2
      real :: sed_avail_h,depo_avail_h,depo_settling_h
      real :: marshface,MEErho,e,insta_retreat,MEEvol,MEEom,MEheight
      real :: dSacch_int,dSacch_edge
      real :: dry_depth
    
      !>> Define depth, in meters, for dry cells that will turn off TSS change calculations 
      !      this is used in other celldXXX subroutines but each subroutine may have a separate dry depth value assigned - double check for consistency
      dry_depth = 0.05
      
!>> Initialize local sediment flux variables to zero          
      do k=1,4
          QSsum(k)=0.0
          QSsumh(k)=0.0
          QSmarsh(k) = 0.0
          QStrib(k) = 0.0
          QSdiv(k) = 0.0
          MEESedRate(k) = 0.0
      enddo

     
!>> Sediment from marsh edge erosion is total mass - must split into sediment classes
      MEESedRatio(1) = 1./3. 
      MEESedRatio(2) = 1./3.
      MEESedRatio(3) = 1./6.
      MEESedRatio(4) = 1./6.
      MEErho = 2650.          ! bulk density of marsh sediments (kg/m3)
      MEEom = 0.6             ! organic matter of marsh sediments (ratio)
   

!>> Calculate density of water (corrected for salinity concentration)
      rhow=1000.*(1+S(j,1)/1000.)

!>> Calculate depth in open water - depth is not allowed to be less than 0.0001 in celldQ, so this value will always be non-zero
      ddy_1 = Es(j,1) - Bed(j)
      ddy_2 = Es(j,2) - Bed(j)
!>> Calculate depth in marsh - marsh depth is not allowed to be less than marsh elevation, but it can be zero, so add a minimum depth to avoid div_by_zero errors
      ddh_1 = Eh(j,1) - BedM(j) + 0.0001
      ddh_2 = Eh(j,2) - BedM(j) + 0.0001
      
!>> Calculate height of marsh edge erosion face from equation 1 in Marsh Edge Erosion methodology report.
!>> Marsh face height comes from Wilson & Allison (2008) 

      !MEE is in m/yr, NTs is number of timesteps per year
      
      MEheight = 1.2
      insta_retreat = MEE(j)/float(NTs)
      MEEvol = insta_retreat*MEheight*hwidth(j)
      
! Loop over sediment size classes      
      do k=1,4
!>> Calculate sediment flux from eroded marsh face - multiply by -1 to keep sign convention of negative sediment flux is INTO open water
! multiply by 1000 to convert from kg/sec to g/sec      
          MEEsedRate(k) = -1000.0*MEESedRatio(k)
     &            *MEEvol*MEErho*(1-MEEom)/dt
            
!>> Calculate sediment flux from marsh-to-open water exchange flows
!>> Negative flux is FROM marsh TO open water
          if(QMarsh(j,2) >= 0.0) then
!>> if flow is into marsh and there are suspended solid in open water determine maximum possible sediment flux (g/sec) based on available sediment
              if (Css(j,1,k) > CSSmin) then
                  QSmarsh_avail = CSS(j,1,k)*As(j,1)*ddy_1/dt
                  QSmarsh(k) = min(Qmarsh(j,2)*Css(j,1,k),QSmarsh_avail)	!g/s !going into marsh
              else
                  QSmarsh(k) = 0.0
              endif
          else
!>> if flow is out of marsh and there are suspended solids in marsh determine maximum possible sediment flux (g/sec) based on available sediment
              if (Cssh(j,1,k) > CSSmin) then
                  QSmarsh_avail = -CSSh(j,1,k)*Ahf(j)*ddh_1/dt
	            QSmarsh(k) = max(Qmarsh(j,2)*Cssh(j,1,k),QSmarsh_avail)	! QSmarsh is negative here (since Qmarsh is negative), therefore the max() is on negative values and will take the value with the smaller magnitude

              else
                  QSmarsh(k) = 0.0
              endif
          endif
    
!>> Calculate sediment flux from tributaries
!>> Negative flux is FROM tributary TO open water
!>> positive tributary flow is INTO compartment - multiply by -1 to keep sign convention of negative sediment flux is INTO open water
          do jn =1,ntrib
              jnt = jtrib(jn)
              tribflow = -Qmult(j,jn)*QTrib(jnt,kday)
              if (tribflow <= 0.0) then
                  QStrib(k) = QStrib(k)
     &                  + tribflow*cssT(jnt,kday,k)
              else
!>> if tributary flow is leaving compartment and there are suspended solids in compartment, determine maximum possible sediment flux (g/sec) based on available sediment
                  if (Css(j,1,k) > CSSmin) then
                      QStrib_avail = CSS(j,1,k)*As(j,1)*ddy_1/dt
                      QStrib(k) = QStrib(k) 
     &                      + min(tribflow*CSS(j,1,k),QStrib_avail)  ! QStrib is positive here
                  else
                      QStrib(k) = QStrib(k) + 0.0
                  endif
              endif
          enddo
!>> Calculate sediment flux from diversions 
!>> Negative flux is FROM diversion TO open water
!>> positive diversion flow is INTO compartment - multiply by -1 to keep sign convention of negative sediment flux is INTO open water
          do it = 1,ndiv
              jit = jdiv(it)
              divflow = -Qmultdiv(j,it)*Qdiv(it,kday)
              if (divflow <= 0.0) then
                  QSdiv(k) = Qsdiv(k)
     &                 + divflow*cssTdiv(jit,kday,k)
              else
!>> if tributary flow is leaving compartment and there are suspended solids in compartment, determine maximum possible sediment flux (g/sec) based on available sediment
                  if(CSS(j,1,k) > CSSmin) then
                      QSdiv_avail = CSS(j,1,k)*As(j,1)*ddy_1/dt
                      QSdiv(k) = QSdiv(k)
     &                     + min(divflow*css(j,1,k),QSdiv_avail)
                  else
                      QSdiv(k) = QSdiv(k) + 0.0
                  endif
              endif
          enddo
      enddo
! End loop over sediment classes

!>> Calculate sediment flux and velocities from link network
!>> note that negative QSSum values for sand (flux OUT of compartment)are included in the SedAccumRate for sand via the van Rijn equilibrium CSS calculation
      ave_vel_sum(j) = 0.0
      ave_vel_cnt(j) = 0.0
      ave_vel(j) = 0.0
      max_vel(j) = 0.0
      min_vel(j) = 0.0

      do k = 1,nlink2cell(j)
          if(icc(j,k) /= 0) then
              if(icc(j,k) < 0)then
                  jnb=jus(abs(icc(j,k)))
              else
                  jnb=jds(abs(icc(j,k)))
              endif  
          endif
              iab=abs(icc(j,k))
!>> if link value in connectivity matrix is non-zero - call subroutine that will accumulate sediment flux from all links
          if (iab /= 0) then
!>> calculate velocity magnitude of all flows into/out of compartment
              if (linkt(iab) > 0) then
                  if (linkt(iab) /= 8) then
                      if (linkt(iab) /= 9) then
                          ave_vel_sum(j) = ave_vel_sum(j) + abs(link_vel(iab))
                          ave_vel_cnt(j) = ave_vel_cnt(j) + 1.0
                          max_vel(j) = max( max_vel(j),abs(link_vel(iab)) )
                          min_vel(j) = min( min_vel(j),abs(link_vel(iab)) )
                      endif
                  endif
              endif  
!>> call link suspended solids computations for non-sand particles (sand link flow is calculated in van Rijn subroutine)
              do sedclass=1,4
                  call TSSOLIDS(mm,iab,jnb,j,k,dz,dzh,dref,sedclass)
                  !if (j == 115) then
                  !    write(*,*) sedclass,iab,jnb,j,QSsum(sedclass)
                  !end if
              enddo
          endif
      enddo
      ave_vel(j) = ave_vel_sum(j)/ave_vel_cnt(j)

!>> Call sediment resuspension and deposition subroutine - this calculates the net sediment accumulation rate (g/sec)      
      call vanRijnSediment(j,kday)

!>> set delta of each accumulation type to zero - will be used to sum over sediment classes, k
!>> marsh accumulation deltas are running sums (since deposition zone varies based on available edge-vs-interior area)
      dSacch_edge = 0.0
      dSacch_int = 0.0
      
!      if (j == 113) then
!          do k=1,4 
!              !write(*,'(4I,1I,F,F,F,F)')  j,k,SedAccumRate(k),insta_retreat,MEE(j),MEESedRate(k)
!              write(*,'(I,2(1x,F),A,6(1x,F))') k,CSS(j,1,k),CSS(j,1,k)*As(j,1)*ddy_1/dt, ' dep:',
!     &         SedAccumRate(k), MEESedRate(k), QSmarsh(k), QSsum(k), QStrib(k), QSdiv(k)
!          end do
!          pause
!      end if

      
!>> open water accumulation delta are kept by sediment class to account for erodible bed of each sediment class
      do k=1,4
          dSacc(k) = 0.0
!>> loop over sediment classes to 
!>> Calculate change in open water CSS concentration for each sediment class,  based on net accumulation rates
          CSS(j,2,k) = ((CSS(j,1,k)*As(j,1)*ddy_1/dt)
     &                - SedAccumRate(k)   
     &                - MEESedRate(k)
     &                - QSmarsh(k)
     &                - QSsum(k)
     &                - QStrib(k)
     &                - QSdiv(k))
     &                *dt/(As(j,1)*ddy_2)

          
!>> Update CSS for each sediment class - filter based on available sediment
          CSS(j,2,k)=min(max(CSSmin,CSS(j,2,k)),CSSmax)
                    
!>> Calculate total mass accumulated on open water bed
!>> If available sediment in erodible bed has already been lost, SedAccumRate was set to zero and only use positive SedAccumRate will be used (deposition)         
          dSacc(k) = SedAccumRate(k)*dt/As(j,1)
                    
!>> Calculate marsh settling and change to CSS in marsh
          if (Ahf(j) > 0.0) then
!>> Calculate sediment suspended in marsh during during timestep (g)
! Negative sediment flux via links (QSsumh) is entering marsh
              sed_avail_h = CSSh(j,1,k)*Ahf(j)*ddh_1
!>> Calculate deposition rate if all sediment in marsh were to settle out (g/sec)
              depo_avail_h = max(sed_avail_h,0.0)/dt
!>> Calculate deposition rate based on current settling velocity (g/sec)
              depo_settling_h = CSSh(j,1,k)*velset(j,k)*Ahf(j)
!>> Determine actual amount of sediment accumulated on marsh based on available sediment              
              SedAccumRate_h(k) = min(depo_settling_h,depo_avail_h)

!>> Calculate change in marsh CSS concentration in each sediment class
              CSSh(j,2,k) = (CSSh(j,1,k)*Ahf(j)*ddh_1/dt
     &            + QSmarsh(k)
     &            - QSsumh(k)
     &            - SedAccumRate_h(k))
     &            *dt/(ddh_2*Ahf(j))
              
              CSSh(j,2,k) = min(max(CSSmin,CSSh(j,2,k)),CSSmax)
              
!>> Calculate mass deposited in marsh zones (in g/m2)
!>> Sand particles deposit in marsh edge zone (first 30 meters) 
              if (k == 1) then
                  if(ar_ed(j) /= 0.0) then
                      !dSacch(k)=SedAccumRate_h(k)*(dt/ar_ed)
                      dSacch_edge = dSacch_edge
     &                        + SedAccumRate_h(k)*(dt/ar_ed(j))
                  endif
!>> Silt particles will deposit on edge if marsh depth is less than some threshold depth
              elseif (k == 2) then
                  if(ar_ed(j) /= 0.0) then
                      if (ddh_2 < 0.3) then
                          dSacch_edge = dSacch_edge
     &                        + SedAccumRate_h(k)*(dt/ar_ed(j))
!>> if marsh depth is greater than threshold, deposit silt on interior (unless all marsh is edge)                      
                      else
                          if(ar_int(j) /= 0.0) then
                              dSachh_int = dSacch_int
     &                             + SedAccumRate_h(k)*(dt/ar_int(j))
                          else
                              dSacch_edge = dSacch_edge
     &                        + SedAccumRate_h(k)*(dt/ar_ed(j))
                          endif
                      endif
                  endif
                  
!>> Calculate amount deposited in interior of marsh (in g/m2)
              else
                  if(ar_int(j) /= 0.0) then
                      dSacch_int = dSacch_int
     &                            + SedAccumRate_h(k)*(dt/ar_int(j))
!>> If there all marsh area is edge, deposit fine material on edge instead of interior
                  elseif (ar_ed(j) /= 0.0) then
                      dSacch_edge = dSacch_edge
     &                    + SedAccumRate_h(k)*(dt/ar_ed(j))
                  endif
              endif              
!>> If compartment has no marsh area, do not calculate values          
          else
              dSacch_edge = 0.0
              dSacch_int = 0.0
              CSSh(j,2,k) = CSSh(j,1,k)
          endif

!>> Check if depths are at or below dry-depth threshold
          if (ddy_2 <= dry_depth) then
              CSS(j,2,k) = 0.0
              !CSS(j,2,k) = CSS(j,1,k)
          endif
          
          if (ddh_2 <= dry_depth) then
              CSSh(j,2,k) = 0.0
              !CSSh(j,2,k) = CSSh(j,1,k)
          endif
          
!>> Check if compartment has an offshore TSS concentration assigned in Cells.csv
!>> If set to -9999 then the CSS calculated above will be used, otherwise the value will be set to 
          if (CSSos(j) >= 0)  then
              Css(j,2,k) = CSSos(j)*BCSedRatio(k)
          endif
          
      enddo
      
     
!>> Calculate total sediment accumulation in each zone (cumulative up to current timestep)
      Sandacc(j,2) = Sandacc(j,1) + dSacc(1)
      Siltacc(j,2) = Siltacc(j,1) + dSacc(2)
      Clayacc(j,2) = Clayacc(j,1) + dSacc(3) + dSacc(4)
      Sacc(j,2) = Sacc(j,1) + dSacc(1) + dSacc(2) + dSacc(3) + dSacc(4)
      Sacch_edge(j,2) = Sacch_edge(j,1) + dSacch_edge !dSacch(1) + dSacch(2)
      Sacch_int(j,2)= Sacch_int(j,1) + dSacch_int !dSacch(3) + dSacch(4)
      Sacch(j,2) = Sacch_int(j,2) + Sacch_edge(j,2)

! Miscellaneous calculations
!      ASandA(j)=ASandT(j,kday)*dt+ASandA(j)					!kg/s*s*1000 = kg*1000=g June 2011  JAM Dec 12, 2010
!      ASandA(j)=ASandD(j,kday)*dt+ASandA(j)						!JAM Dec 12 2010 *** June 21, 2011
!      sbm=-SourceBM(j)*Ahf(j)*1000./(365.25*24*3600)				!Biomass SS source in marsh !JAM Oct 2010

     
!     if (daystep == lastdaystep) then
!         if(j == 638)then
!             write(*,*) 'CSS:'
!             write(*,*) CSS(j,2,1),CSS(j,2,2),
!    &                      CSS(j,2,3),CSS(j,2,4)
!             write(*,*) 'CSSh:'
!             write(*,*) CSSh(j,2,1),CSSh(j,2,2),
!    &                            CSSh(j,2,3),CSSh(j,2,4)
!             write(*,*) 'Sacc:',sacc(j,2)
!             write(*,*) 'SedAccumrate',SedAccumRate(1),
!    &                SedAccumRate(2),SedAccumRate(3),SedAccumRate(4)
!             write(*,*) 'sacch_edge:',sacch_edge(j,2)
!             write(*,*) SedAccumRate_h(1),SedAccumRate_h(2)
!             write(*,*) 'sacch_int:',sacch_int(j,2)
!             write(*,*) SedAccumRate_h(3),SedAccumRate_h(4)
!             write(*,*) 'Qsmarsh:',QSmarsh(1),QSmarsh(2),
!    &                            QSMarsh(3),QSMarsh(4)
!             write(*,*) 'Qmarsh:',Qmarsh(j,2)
!             write(*,*) 'depth marsh:',ddh_2
!             write(*,*) 'area marsh:', Ahf(j)
!         endif
!     endif

     
      return 
	end
