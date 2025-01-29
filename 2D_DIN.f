!      Subroutine DIN(DChemSum,ichem,mex,j,k)
! kday now global parameter - no longer needed to be passed into subroutine   
      Subroutine DIN(DChemSum,ichem,mex,j,k,kday)       ! Nitrate + nitrite as N        
    !JAM Chem # 3  JAM March 12 2011

      use params      
      
      rca=0.004
      fixN=0.00000001
      fixNO3=rca*Chem(j,8,1)*fixN
      Tcorn=1.085**(Tempw(j,2)-20.)                         ! Temperature Correction  Nitrification
      Tcord=1.05**(Tempw(j,2)-20.)                          ! Temperature Correction deNitrification
      ALAG=0.2
      PP=1-0.2*cos(2*PI*(Float(kday)/365.25-ALAG)) 
      fnodin=Chem(j,1,1)/(Chem(j,3,1)+10E-08)                   ! Assumes proportionate uptake of NO3 and NH4
      FNoC=1./5.681
      cxx=1/3600/24/1000000.
      fcN=1.
      FSEASON2=1.0+0.2*cos(2*PI*(Float(kday)/365.25+0.25))  ! JAM Jan 09 11
      frr=0.82
      fmm=1
!   dz=ES(j,2)-ES(j,1)
      Chem(j,3,2)=Chem(j,1,2)+Chem(j,2,2)
      QChemsum(3)=QChemsum(1)+QChemsum(2)                       ! JAM March 12 2011 
    
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

!c***********************End Subroutine for DIN**************************************************