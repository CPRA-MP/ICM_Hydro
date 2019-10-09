!      Subroutine DOx(DChemSum,ichem,mex,j,k)
      ! kday now global parameter - no longer needed to be passed into subroutine   
	Subroutine DOx(DChemSum,ichem,mex,j,k,kday)		! DO as O
	!JAM Oct 2010 Chem #7
      
      use params
      

	FOoC=2.*16/14.								! JAM guess
	FOoN=3*16/12.								! JAM guess
	REO=.1										! Re-aeration rate per day
	DDO=O2Sat(kday)-Chem(j,ichem,1)/consd
	IF(DDO.LT.0.0)ddo=0.0
	windspd = sqrt(windx(j)**2+windy(j)**2)
      DChemSUM=DChemSUM+decay(ichem,ichem)*windspd*windspd*DDO	! Re-aeration
!      DChemSUM=DChemSUM+decay(ichem,ichem)*wspd(kday)*wspd(kday)*DDO	! Re-aeration
	DChemSUM=DChemSUM+decay(ichem,2)*Chem(j,2,1)*FOoN				! Nitrification 
      DChemSUM=DChemSUM+decay(ichem,8)*Chem(j,8,1)					! Net photosynthesis and respiration  
      DChemSUM=DChemSUM+decay(ichem,9)*Chem(j,9,1)*FOoC				! Loss by oxydation of TOC as Dead Algae C

	return
	end
c______________________________________________________________________________________
cc! NO3	NH4	DIN	ON	 TP	TOC	DO	  LA-C	DA-C	DON	DOP	SRP	ChLa POP
cc   1	2	3	4	 5	 6	7	  8	     9	   10	11	12	13	  14 ***********
c*********************************************************************************
cc-0.0815	0.1	0	0	 0	 0	0    -0.5	 0	   0	 0	 0	 0	  0		!1 NO2+NO3
cc  0 -0.21	0	0    0	 0	0    -0.5	 0	   0.01	 0	 0	 0	  0		!2 NH4
c    0	0	0	0	 0	 0	0	  0	     0	   0	 0	 0	 0	  0		!3 DIN
cc 0-0.00001 0 -0.000001 0 0	0     1	-0.00001   0	 0	 0	 0	  0		!4 ON
cc 0	    0	0	0 -0.002 0	0     -1	  0    0     0	 0	 0	  0		!5 TP
c   0	   0    0	0	0	0 -0.005 0.265 -0.055  0	 0	 0	 0	  0		!8 LiveAlgae
c   0	   0    0	0	0	0 -0.005-0.055  0.055 -0.01-0.005 0	 0	  0		!9 DeadAlgae
cc  0 -0.01	0	0	0	0	0    0.0	0.005	0	 0	 0	 0	  0		!10 DON
cc  0	   0	0	0	0	0	0   0.01	  0	   0  0.01 -0.01 0	  0		!11 DOP
cc  0	   0	0	0	0	0	0	  0	      0	   0  -0.01	0.01 0	  0		!12 SRP
cc  0	   0	0	0	0	0	0	0.06667	  0	   0	0	0  0.01	  0		!13 Chla
cc  0	   0	0	0	0	0	0	  1	      0	-0.01	 0	 0	 0	  0		!14 POP
c**********************************************************************************
cc! NO3	NH4	DIN	ON	 TP	TOC	DO	  LA-C	DA-C	DON	DOP	SRP	ChLa POP
c______________________________________________________________________________________

c***********************End Subroutine for chemical DO as O ************************************