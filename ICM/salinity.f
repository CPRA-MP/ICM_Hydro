c***********************Start Subroutine for face salinity**************************************
    
      Subroutine salinity(mm,iab,jnb,j,k,Qsalsum)		!face densities from node densities

      use params
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
      
      
      if(iab /= 0.0)then											!ne = 0 there no link available go to 22
	    if (linkt(iab) == 8) then
              if(Ahf(j) == 0) then
!              if (Eh(j,1) - Bedm(j) > 0.3) then
                  if(Q(iab,1) >= 0.0) then
                      Csalface=((fa(iab)*S(jus(abs(icc(j,k))),1)+
     &					fb(iab)*S(jds(abs(icc(j,k))),1)))
                  else
                      Csalface=((fa(iab)*S(jds(abs(icc(j,k))),1)+
     &					fb(iab)*S(jus(abs(icc(j,k))),1)))
                  endif
                  diffus = EAOL(iab)
              else
                  Csalface = 0.0
                  diffus = 0.0
              endif  
          else    
              if(Q(iab,1) >= 0.0) then
				Csalface= ((fa(iab)*S(jus(abs(icc(j,k))),1)				!cell face values
     &				+fb(iab)*S(jds(abs(icc(j,k))),1)))
              else
				Csalface= ((fa(iab)*S(jds(abs(icc(j,k))),1)+
     &					fb(iab)*S(jus(abs(icc(j,k))),1)))
			endif
              diffus = EAOL(iab)              
   		endif
   
          QSalSum=QSalSum + sicc(j,k)*(Q(abs(icc(j,k)),1))*Csalface
     &                +fe*diffus*(S(j,1)-S(jnb,1))
     
          SL(iab,2)=Csalface
      endif
      
          
      return
	end


c***********************End Subrout