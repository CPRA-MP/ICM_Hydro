!cc**********************Start Subroutine for Change in Cell Chemistry*******JAM Oct 2010********

!      Subroutine chemical(mm,iab,jnb,j,k,ichem)	! face values from node values **** Collect Advective diffusion terms
    Subroutine chemical(iab,jnb,j,k,ichem)	! face values from node values **** Collect Advective diffusion terms

	use params      

	if(iab /= 0) then
        if(Q(iab,2) >= 0.0) then
		    Cchemface(ichem) = Chem(jus(iab),ichem,1)
!              Cchemface(ichem)= ((fa(iab)*
!     &            Chem(jus(abs(icc(j,k))),ichem,1)
!     &			+ fb(iab)*Chem(jds(abs(icc(j,k))),ichem,1)))
        else
		    Cchemface(ichem) = Chem(jds(iab),ichem,1)
!		    Cchemface(ichem) =
!     &        ((fa(iab)*Chem(jds(abs(icc(j,k))),ichem,1)
!     &        + fb(iab)*Chem(jus(abs(icc(j,k))),ichem,1)))					
        endif
          
        diffus = EAOL(iab)
                  
        QChemSUMflows(ichem)=QChemSUMflows(ichem) + sicc(j,k)*   
     &			(Q(iab,2))*Cchemface(ichem)
     &			+fe*diffus*(Chem(j,ichem,1)-Chem(jnb,ichem,1))  !zw 4/28/2015 delete /2.
          
    endif
	
	return
	end

!c***********************End Subroutine for Chemical Conditions**********************************