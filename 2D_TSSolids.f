!cc**********************Start Subroutine for TSS Solids*****************************************
    
!      Subroutine TSSOLIDS(mm,iab,jnb,j,k,dz,dzh,dref,sedclass)		!face densities from node densities
      Subroutine TSSOLIDS(iab,jnb,j,k,sedclass)		!face densities from node densities

      use params
      
      implicit none
      integer :: iab,jnb,j,k,sedclass
      real :: Cssface,diffus,Qlink,d1,d2
 
!>> Check if link is marsh overland flow link, if so use CSS in marsh and add sediment to marsh cumulative sediment flux term
      if (linkt(iab) == 8) then
!>> Only calculate marsh sediment exchange if depth in marsh is greater than 0.0 (or set to 0.3 m for assumed surge conditions)
          if (sedclass == 1) then
!>> Overland marsh links do not transport sand          
!              QSsumh(sedclass) = 0.0 ! do not transport sand across marsh links
              QSsumh(sedclass) = QSsumh(sedclass) + 0.0 ! do not transport sand across marsh links
          else          
!              if (Eh(j,1)-Bedm(j) > 0.0) then !0.3) then
!>> calculate average concentration at interface between marshes (directionality of flow is needed to apply upwind factor correctly)
              if (Q(iab,2) >= 0.0) then
                  Cssface = CSSh(jus(abs(icc(j,k))),1,sedclass)
              else    
                  Cssface = CSSh(jds(abs(icc(j,k))),1,sedclass)
              endif
              diffus = EAOL(iab) 
              Qlink=Q(iab,2)
              
              if (Ahf(j)>0) then			  
!>> add sediment flux from this marsh link to the cumulative sediment flux for marsh in this comparment (QSsumh)
                  QSSumh(sedclass) = QSsumh(sedclass)
     &              + sicc(j,k)*Qlink*Cssface
     &              + fe*diffus*(CSSh(j,1,sedclass)-CSSh(jnb,1,sedclass))			!Diffusion
              else
                  QSSum(sedclass) = QSSum(sedclass) 
     &              + sicc(j,k)*Qlink*Cssface				!face exchange
     &              + fe*diffus*(CSSh(j,1,sedclass)-CSSh(jnb,1,sedclass))			!Diffusion
              endif
          endif
!>> If link is not marsh overland flow link, use CSS in open water and add sediment flux to open water cumulative sediment flux term
      else
!>> calculate concentration at interface between compartments (directionality of flow is needed)
          if (Q(iab,2) >= 0.0) then
              Cssface = CSS(jus(abs(icc(j,k))),1,sedclass)
          else
              Cssface = CSS(jds(abs(icc(j,k))),1,sedclass)
          endif
          diffus = EAOL(iab)              
          Qlink = Q(iab,2)

! diffusion term reinforcement, although EAOL has been dealed with in hydrod.f
! no diffusion associated with weir links unless submerged         
          if(linkt(iab)==2) then
              d1 = Es(j,1)-Latr1(iab) 
              d2 = Es(jnb,1)-Latr1(iab) 
              if((d1<0) .or. (d2<0)) then
                  diffus = 0.0
              endif
! pump link has no diffusion
          elseif (linkt(iab) == 7) then
              diffus = 0.0
! ridge link distribute water to OW & marsh proportionally based on area 
! no diffusion associated with ridge links unless submerged         
          elseif(linkt(iab)==9) then
              Qlink = Q(iab,2)*As(j,1)/(As(j,1)+Ahf(j))
              Qlink2 = Q(iab,2)*Ahf(j)/(As(j,1)+Ahf(j))
              d1 = Es(j,1)-Latr1(iab) 
              d2 = Es(jnb,1)-Latr1(iab) 
              if((d1<0) .or. (d2<0)) then
                  diffus = 0.0
              endif
!>> add sediment flux from ridge link to the cumulative sediment flux for marsh in this comparment (QSsumh)
              if(sedclass/=1) then
                  QSSumh(sedclass) = QSsumh(sedclass)
     &            + sicc(j,k)*Qlink2*Cssface
     &            + fe*diffus*(CSS(j,1,sedclass)-CSS(jnb,1,sedclass))			!Diffusion
              endif
          endif
          
          !! update cssface for Atchafalya diversion links to have a SWR of 0.5 of sand - MP project # 03b.DI.04
          !if (iab == 3859) then
          !    Cssface_new = 0.5*Cssface
          !    Cssface = Cssface_new
          !endif
          !
          
!>> Calculate incoming flux for sand
!>> van Rijn suspended sand flux (in g/sec) entering compartment from current link inflows and their respective van Rijn equilibrium flow CSS from previous timestep
!>> sand flux out of compartment is set by equilibrium rate in vanRijn subroutine and this will only collect incoming sand fluxes
!          if (sedclass == 1) then
!              QSsum(sedclass) =  QSsum(sedclass) 
!     &            + sicc(j,k)*Q(abs(icc(j,k)),2)*Cssface
              
!          else
!>> add sediment flux from this link to the cumulative sediment flux for the compartment (QSsum)			    
          QSSum(sedclass) = QSSum(sedclass) 
     &            + sicc(j,k)*Qlink*Cssface				!face exchange
     &            + fe*diffus*(CSS(j,1,sedclass)-CSS(jnb,1,sedclass))			!Diffusion
              
!          endif
      endif

      
	
      return
      end

!c***********************End Subroutine for TSS SOLIDS*******************************************
