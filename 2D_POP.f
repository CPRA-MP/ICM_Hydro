!      Subroutine POP(DChemSum,ichem,mex,j,k)
      
! kday now global parameter - no longer needed to be passed into subroutine   
      Subroutine POP(DChemSum,ichem,mex,j,k,kday)			!POP= chemical PArtition SRP = ParP*TP
	!JAM Oct 2010 Chem #14  April 16 2011

	  use params
      

	  FNoC=1./5.681
	  FPoC=1/41.
	  fnodin=Chem(j,1,1)/(Chem(j,3,1)+10E-08)
	  DChemSUM=DChemSUM+decay(ichem,11)*Chem(j,9,1)*FPoC	! Dissolution = anaerobic decomposition of organic material  
	  return
	  end

!c***********************End Subroutine for chemical PArtition SRP = ParP*TP --> POP*************