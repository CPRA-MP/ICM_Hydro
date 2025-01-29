!cc**********************Start Subroutine for Change in Cell Chemistry*******JAM Oct 2010********

!      Subroutine chemical(mm,iab,jnb,j,k,ichem)    ! face values from node values **** Collect Advective diffusion terms
      Subroutine chemical(iab,jnb,j,k,ichem)    ! face values from node values **** Collect Advective diffusion terms

      use params      

      implicit none
      integer :: iab,jnb,j,k,ichem
      real :: diffus,Qlink,d1,d2
      
      if(abs(Q(iab,2)) > 0) then
          if(Q(iab,2) >= 0.0) then
              Cchemface(ichem) = Chem(jus(iab),ichem,1)
!              Cchemface(ichem)= ((fa(iab)*
!     &            Chem(jus(abs(icc(j,k))),ichem,1)
!     &         + fb(iab)*Chem(jds(abs(icc(j,k))),ichem,1)))
          else
              Cchemface(ichem) = Chem(jds(iab),ichem,1)
!           Cchemface(ichem) =
!     &        ((fa(iab)*Chem(jds(abs(icc(j,k))),ichem,1)
!     &        + fb(iab)*Chem(jus(abs(icc(j,k))),ichem,1)))                 
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
! exclude marsh link flow contribution to openwater temperature
!ZW 1/13/24          elseif ((linkt(iab) == 8) .and. (Ahf(j) > 0)) then
!              CTMPface = 0.0
!              diffus = 0.0
! ridge link distribute water to OW & marsh proportionally based on area 
! no diffusion associated with ridge links unless submerged         
          elseif(linkt(iab)==9) then
!ZW 1/13/24              Qlink = Q(iab,2)*As(j,1)/(As(j,1)+Ahf(j))
              d1 = Es(j,1)-Latr1(iab) 
              d2 = Es(jnb,1)-Latr1(iab) 
              if((d1<0) .or. (d2<0)) then
                  diffus = 0.0
              endif
          endif
          
          QChemSUMflows(ichem)=QChemSUMflows(ichem)    
     &          +sicc(j,k)*Qlink*Cchemface(ichem)
     &          +fe*diffus*(Chem(j,ichem,1)-Chem(jnb,ichem,1))  !zw 4/28/2015 delete /2.
          
      endif
    
      return
      end

!c***********************End Subroutine for Chemical Conditions**********************************