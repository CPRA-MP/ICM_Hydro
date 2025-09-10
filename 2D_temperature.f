!cc**********************Start Subroutine for face temperature****JAM Oct 2010*******************
    
!      Subroutine Temperature(mm,iab,jnb,j,k,QTMPsum)		!face densities from node densities
      Subroutine Temperature(iab,jnb,j,k,QTMPsum)		!face densities from node densities

!cc      TempMR(100,366*20) Miss River Temperature 
!cc      Tempair(100,366*20) Air Temperature
!cc      Tempw(100,366*20) Cell water Temperature
!cc      Tempe(100,366*20) Cell equilibrium Temperature
!cc      TL(200,3) Link water Temperature

      use params
	  
      implicit none
      integer :: iab,jnb,j,k
      real :: QTMPsum,CTMPface,diffus,Qlink,d1,d2
      real :: fa_DS

      if(abs(Q(iab,2)) == 0.0) then
          TL(iab,1)=0
      else
!          if (linkt(iab) == 8) then
!              if (Eh(j,1) - Bedm(j) > 0.3) then
!                  if(Q(iab,1) >= 0.0) then    
!                      CTMPface=Tempw(jus(abs(icc(j,k))),1)
!                  else
!                      CTMPface=Tempw(jds(abs(icc(j,k))),1)
!                  endif
!                  diffus = EAOL(iab)
!              else
!                  CTMPface = 0.0
!                  diffus = 0.0
!              endif
          diffus = EAOL(iab)
          Qlink = Q(iab,2)

!==ZW 2/1/2024 add for Blended Differencing scheme
          if (iAdvTrans==1) then
              if ((linkt(iab) == 1) .or. (linkt(iab) == 3) 
     &           .or. (linkt(iab) == 6) .or. (linkt(iab) == 11)  
     &           .or. (linkt(iab) == 12).or. (linkt(iab) == 8)) then
 !                if(fa(iab)<1)then
 !                    fa_DS=2.0-r_BD-fa(iab)
 !                else
 !                    fa_DS=fa(iab)
 !                endif
                 if(Q(iab,2) > 0.0) then			 
                     fa_US=1.0-r_BD*(1-fx_marsh)  !this is fa for flow from US to DS (Q>0)
				 else
                     fa_DS=1.0-r_BD*fx_marsh  !this is fa for flow from DS to US (Q<0)
				 endif
              else
                 fa_US=fa(iab)
                 fa_DS=fa(iab)
              endif
          else
              fa_US=fa(iab)
              fa_DS=fa(iab)
          endif
!===

          if(Q(iab,2) >= 0.0) then
!              CTMPface=Tempw(jus(iab),1)
              CTMPface= fa_US*Tempw(jus(iab),1)				!cell face values
     &                  +(1-fa_US)*Tempw(jds(iab),1)
          else
!              CTMPface=Tempw(jds(iab),1)
              CTMPface= fa_DS*Tempw(jds(iab),1)
     &                  +(1.0-fa_DS)*Tempw(jus(iab),1)
          endif    

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

          QTMPSum=QTMPSum + sicc(j,k)*Qlink*CTMPface
     &            +fe*diffus*(Tempw(j,1)-Tempw(jnb,1))
                  
          TL(iab,1)=CTMPface
              
      endif

      return
      end

!c***********************End Subroutine face temperature****JAM Oct 2010*************************