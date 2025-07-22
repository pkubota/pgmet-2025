!
!  $Author: pkubota $
!  $Date: 2008/09/23 17:51:54 $
!  $Revision: 1.9 $
!
MODULE Constants

  IMPLICIT NONE
  PRIVATE       

  ! Selecting Kinds
  INTEGER,PUBLIC, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER,PUBLIC, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER,PUBLIC, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER,PUBLIC, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER,PUBLIC, PARAMETER :: r16 = SELECTED_REAL_KIND(15)! Kind for 128-bits Real Numbers
  ! Selecting Unit

  INTEGER,PUBLIC, PARAMETER :: nfprt=0
  REAL (KIND=r8), PARAMETER :: Avogrado=6022.52_r8  ! Avogrado's Constant (10**-20/kmol)
  REAL (KIND=r8), PARAMETER :: Boltzmann =1.38054_r8! Boltzmann's Constant (J/K)
  REAL (KIND=r8), PARAMETER :: REstar=Avogrado*Boltzmann    ! Universal Gas Constant (10**23 J/K/kmol)
  REAL (KIND=r8), PARAMETER :: MWWater=18.016_r8   ! Molecular Weight of Water (kg/kmol)
  REAL (KIND=r8), PARAMETER :: MWN2 =28.0134_r8     ! Molecular Weight of Nitrogen (kg/kmol)
  REAL (KIND=r8), PARAMETER :: PN2 =0.7809_r8      ! Atmosphere Nitrogen Percetage
  REAL (KIND=r8), PARAMETER :: MWO2=31.9988_r8      ! Molecular Weight of Oxigen (kg/kmol)
  REAL (KIND=r8), PARAMETER :: PO2=0.2095_r8       ! Atmosphere Oxigen Percetage
  REAL (KIND=r8), PARAMETER :: MWAr=39.9480_r8      ! Molecular Weight of Argon (kg/kmol)
  REAL (KIND=r8), PARAMETER :: PAr=0.0093_r8       ! Atmosphere Argon Percetage
  REAL (KIND=r8), PARAMETER :: MWCO2=44.0103_r8     ! Molecular Weight of Carbon Gas (kg/kmol)
  REAL (KIND=r8), PARAMETER :: PCO2 =0.0003_r8     ! Atmosphere Carbon Gas Percetage
  REAL (KIND=r8), PARAMETER :: MWDAir=PN2*MWN2+PO2*MWO2+PAr*MWAr+PCO2*MWCO2! Mean Molecular Weight of Dry Air (kg/kmol)
  REAL (KIND=r8),PUBLIC, PARAMETER :: Rd=REstar/MWDAir     ! Dry Air Gas Constant (m2/s2/K)
  REAL (KIND=r8),PUBLIC, PARAMETER :: Cp=3.5_r8*Rd     ! Dry Air Gas Specific Heat Capcity at Pressure Constant (m2/s2/K)
  REAL (KIND=r8),PUBLIC, PARAMETER :: Rv=REstar/MWWater

  REAL(kind=r8),PUBLIC, PARAMETER :: pi     =4.0*atan(1.0_r8)
  REAL(kind=r8),PUBLIC, PARAMETER :: Deg2Rad=pi/180.0_r8
  REAL(kind=r8),PUBLIC, PARAMETER :: r_earth=6.37e6_r8 !m
  REAL(KIND=r8),PUBLIC, PARAMETER :: omega = 7.27e-5_r8
  REAL(KIND=r8),PUBLIC, PARAMETER :: kappa = Rd/Cp
  REAL(KIND=r8),PUBLIC, PARAMETER ::  Eps=Rd/Rv
  REAL(KIND=r8),PUBLIC, PARAMETER ::  Eps1=1.0_r8-Eps

  REAL (KIND=r8),PUBLIC, PARAMETER :: CTv=Eps1/Eps  ! Constant Used to Convert Tv into T, or vice-versa.

  PUBLIC :: InitClassModuleConstants

CONTAINS

 SUBROUTINE InitClassModuleConstants () 
    IMPLICIT NONE


 END SUBROUTINE InitClassModuleConstants

END MODULE Constants
