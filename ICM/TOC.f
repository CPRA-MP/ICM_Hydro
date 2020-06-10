!      Subroutine TOC(DChemSum,ichem,mex,j,k)
      
! kday now global parameter - no longer needed to be passed into subroutine   
      Subroutine TOC(DChemSum,ichem,mex,j,k,kday)	! TOC as C = C in live algae + C in dead algae
	!JAM Oct 2010 Chem #6
	
      use params
      
      DChemSUM=DChemSUM+decay(ichem,8)*GrowAlgae(j,8,8)*Chem(j,8,1)

	return
	end

c***********************End Subroutine for chemical Total Organic Carbon (TOC)******************