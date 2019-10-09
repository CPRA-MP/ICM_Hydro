!      Subroutine DIP(DChemSum,ichem,mex,j,k)

      ! kday now global parameter - no longer needed to be passed into subroutine 
      Subroutine DIP(DChemSum,ichem,mex,j,k,kday)			!DIP=PArtition SRP = ParP*TP	

      use params
      
      !JAM Oct 2010 Chem #12

	FNoC=1./5.681
	fnodin=Chem(j,1,1)/(Chem(j,3,1)+10E-08)
      DChemSUM=DChemSUM+decay(ichem,11)*Chem(j,11,1)		!dissolution = anaerobic decomposition of organic material   

	return
	end
 
c***********************End Subroutine for chemical PArtition SRP = ParP*TP --> DIP*************