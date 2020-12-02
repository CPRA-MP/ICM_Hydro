cc**********************Start Subroutine for face temperature****JAM Oct 2010*******************
    
      Subroutine Temperature(mm,iab,jnb,j,k,QTMPsum)		!face densities from node densities

cc      TempMR(100,366*20) Miss River Temperature 
cc      Tempair(100,366*20) Air Temperature
cc      Tempw(100,366*20) Cell water Temperature
cc      Tempe(100,366*20) Cell equilibrium Temperature
cc      TL(200,3) Link water Temperature

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
!          else    
              if(Q(iab,2) >= 0.0) then
                  CTMPface=Tempw(jus(abs(icc(j,k))),1)
              else
                  CTMPface=Tempw(jds(abs(icc(j,k))),1)
              endif    
              diffus = EAOL(iab)
!          endif

          QTMPSum=QTMPSum + sicc(j,k)*(Q(abs(icc(j,k)),2))*CTMPface
     &			    +fe*diffus*(Tempw(j,1)-Tempw(jnb,1))
                  
          TL(iab,1)=CTMPface
              
      endif
	return
	end

c***********************End Subroutine face temperature****JAM Oct 2010*************************