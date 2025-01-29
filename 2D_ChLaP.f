!      Subroutine ChLaP(DChemSum,ichem,mex,j,k)
      
! kday now global parameter - no longer needed to be passed into subroutine   
      Subroutine ChLaP(DChemSum,ichem,mex,j,k,kday)     ! PArtition Chla= ParCla*LivA   
    !JAM Oct 2010 Chem #13

      use params      
      

      FNoC=1./5.681
      fdead=1.0
      fnodin=Chem(j,1,1)/(Chem(j,3,1)+10E-08)
      fPL=1.0                                               ! Modified CHLA JAM March 27, 2011
      fChLa=1./20.*(1.+1./(s(j,1)+1.))
      GM7=1+Chem(j,8,1)/(0.00251+Chem(j,8,1))
      if(j.eq.14)then                                       ! Changed from 2012 SMP Cell 6 to Lower B&B Cell 14 Oct 16 2013 JKS
          fLP=2. 
      endif
      ffc=fPL*fChla
      DChemSUM=DChemSUM+GrowAlgae(j,ichem,ichem)*Chem(j,ichem,1)*ffc        ! increase due to uptake
     &  +decay(ichem,9)*Chem(j,ichem,1)*Tcorn/1.5*gm7*fdead             ! *Chem(j,9,1)/        ! net loss ~ death rate & elimination rate  = - DA

      return
      end
!c______________________________________________________________________________________
!cc! NO3    NH4 DIN ON   TP TOC DO    LA-C  DA-C    DON DOP SRP ChLa POP
!cc   1 2   3   4    5   6  7     8      9     10   11  12  13    14 ***********
!c*********************************************************************************
!cc-0.0815  0.1 0   0    0   0  0    -0.5    0     0     0   0   0    0     !1 NO2+NO3
!cc  0 -0.21    0   0    0   0  0    -0.5    0     0.01  0   0   0    0     !2 NH4
!c    0 0   0   0    0   0  0     0      0     0     0   0   0    0     !3 DIN
!cc 0-0.00001 0 -0.000001 0 0   0     1 -0.00001   0     0   0   0    0     !4 ON
!cc 0       0   0   0 -0.002 0  0     -1      0    0     0   0   0    0     !5 TP
!c   0     0    0   0   0   0 -0.005 0.265 -0.055  0     0   0   0    0     !8 LiveAlgae
!c   0     0    0   0   0   0 -0.005-0.055  0.055 -0.01-0.005 0  0    0     !9 DeadAlgae
!cc  0 -0.01    0   0   0   0   0    0.0    0.005   0    0   0   0    0     !10 DON
!cc  0     0    0   0   0   0   0   0.01      0    0  0.01 -0.01 0    0     !11 DOP
!cc  0     0    0   0   0   0   0     0       0    0  -0.01 0.01 0    0     !12 SRP
!cc  0     0    0   0   0   0   0   0.06667   0    0    0   0  0.01   0     !13 Chla
!cc  0     0    0   0   0   0   0     1       0 -0.01    0   0   0    0     !14 POP
!c**********************************************************************************
!cc! NO3    NH4 DIN ON   TP TOC DO    LA-C  DA-C    DON DOP SRP ChLa POP
!c______________________________________________________________________________________

!c***********************End Subroutine for chemical PArtition DeadA*StchP --> ChLaP*************
