!cc**********************Start Subroutine for face temperature****JAM Oct 2010*******************
    
!      Subroutine Temperature(mm,iab,jnb,j,k,QTMPsum)		!face densities from node densities
      Subroutine Temperature(iab,jnb,j,k,QTMPsum)		!face densities from node densities

!cc      TempMR(100,366*20) Miss River Temperature 
!cc      Tempair(100,366*20) Air Temperature
!cc      Tempw(100,366*20) Cell water Temperature
!cc      Tempe(100,366*20) Cell equilibrium Temperature
!cc      TL(200,3) Link water Temperature

      use params


      if(iab /= 0) then
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
          if(Q(iab,2) >= 0.0) then
              CTMPface=Tempw(jus(iab),1)
          else
              CTMPface=Tempw(jds(iab),1)
          endif    
          diffus = EAOL(iab)
          Qlink = Q(iab,2)

! exclude marsh link flow contribution to openwater temperature
          if ((linkt(iab) == 8) .and. (Ahf(j) > 0)) then
              CTMPface = 0.0
              diffus = 0.0
          endif

! ridge link distribute water to OW & marsh proportionally based on area 
! no diffusion associated with ridge links unless submerged         
          if(linkt(iab)==9) then
              Qlink = Q(iab,2)*As(j,1)/(As(j,1)+Ahf(j))
              d1 = Es(j,1)-Latr1(iab) 
              d2 = Es(jnb,1)-Latr1(iab) 
              if((d1<0) .or. (d2<0)) then
                  diffus = 0.0
              endif
          endif

! pump link has no diffusion
          if (linkt(iab) == 7) then
              diffus = 0.0
          endif

          QTMPSum=QTMPSum + sicc(j,k)*Qlink*CTMPface
     &            +fe*diffus*(Tempw(j,1)-Tempw(jnb,1))
                  
          TL(iab,1)=CTMPface
              
      endif

      return
      end

!c***********************End Subroutine face temperature****JAM Oct 2010*************************