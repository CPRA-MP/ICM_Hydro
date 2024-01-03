!      Subroutine OrgN(DChemSum,ichem,mex,j,k)
      
! kday now global parameter - no longer needed to be passed into subroutine   
    Subroutine OrgN(DChemSum,ichem,mex,j,k,kday) !Organic Nitrogen as N 
	!JAM Oct 2010 Chem #4
	
      ! Subroutine no longer used - organic Nitrogen calcualted from Algae and Detritus concentrations in celldChem !-EDW
      
      
    use params
      
 
	FNoC=1./5.681											! Nitrogen to carbon ratio in algae
	DChemSUM=DChemSUM+decay(ichem,ichem)*Chem(j,ichem,1)	! Loss to mineralization
	DChemSUM=DChemSUM+decay(ichem,2)*Chem(j,ichem,1)		! Loss due to anaerobic decomposition ON --> NH4
    DChemSUM=DChemSUM+decay(ichem,8)*GrowAlgae(j,8,8)*		! Addition by uptake by live algae
     &	Chem(j,8,1)*FNoC
	
	return
	end

!c***********************End Subroutine for chemical OrgN****************************************