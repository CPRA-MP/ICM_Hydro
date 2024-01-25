!c***********************Start Subroutine for face salinity**************************************
    
!      Subroutine salinity(mm,iab,jnb,j,k,Qsalsum)		!face densities from node densities
      Subroutine salinity(iab,jnb,j,k,Qsalsum)		!face densities from node densities

      use params
	  
      implicit none
      integer :: iab,jnb,j,k
      real :: Qsalsum,Csalface,diffus,Qlink,d1,d2,cfacemax

!>@par General Structure of Subroutine Logic:
!>>
!      if(iab /= 0.0)then
!          if(Q(iab,1) >= 0.0) then
!              Csalface = S(jus(abs(icc(j,k))),1)
!          else
!              Csalface = S(jds(abs(icc(j,k))),1)
!          endif
!          
!          diffus = EAOL(iab)
!          
!          QSalSum=QSalSum + sicc(j,k)*(Q(abs(icc(j,k)),1))*Csalface
!     &                +fe*diffus*(S(j,1)-S(jnb,1))
!    
!          SL(iab,2)=Csalface
!      endif
      
      if(abs(Q(iab,2)) == 0.0) then
          SL(iab,2)=0
          QSalSum=QSalSum+0.0
      else
!	    if (linkt(iab) == 8) then
!              if(Ahf(j) == 0) then
!!              if (Eh(j,1) - Bedm(j) > 0.3) then
!                  if(Q(iab,2) >= 0.0) then
!                      Csalface=((fa(iab)*S(jus(abs(icc(j,k))),1)+
!     &					fb(iab)*S(jds(abs(icc(j,k))),1)))
!                  else
!                      Csalface=((fa(iab)*S(jds(abs(icc(j,k))),1)+
!     &					fb(iab)*S(jus(abs(icc(j,k))),1)))
!                  endif
!                  diffus = EAOL(iab)
!              else
!                  Csalface = 0.0
!                  diffus = 0.0
!              endif  
!          else    
          diffus = EAOL(iab)              
          Qlink = Q(iab,2)
          if(Q(iab,2) > 0.0) then
!              Csalface= ((fa(iab)*S(jus(iab),1)				!cell face values
!     &                  +fb(iab)*S(jds(iab),1)))
              Csalface= (0.75*S(jus(iab),1)				!cell face values
     &                  +0.25*S(jds(iab),1))
          else
!              Csalface= ((fa(iab)*S(jds(iab),1)
!     &                  +fb(iab)*S(jus(iab),1)))
              Csalface= (0.75*S(jds(iab),1)
     &                  +0.25*S(jus(iab),1))
          endif
   		!endif
   
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

!check EAOL
          if (isNan(diffus)) then
              write(1,*)'EAOL missing in link: ',iab
              diffus=0.0
          endif       

          QSalSum=QSalSum + sicc(j,k)*Qlink*Csalface
     &                +fe*diffus*(S(j,1)-S(jnb,1))

          SL(iab,2)=Csalface

      endif
          
      return
	  end


!c***********************End Subrout