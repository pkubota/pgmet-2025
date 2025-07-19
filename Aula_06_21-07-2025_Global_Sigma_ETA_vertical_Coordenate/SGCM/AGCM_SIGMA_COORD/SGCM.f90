MODULE MODEL_SGCM
!   Simplified Global Circulation Model
!
!     by Paulo Kubota  January 30, 2019


IMPLICIT NONE
PRIVATE
INTEGER,PUBLIC, PARAMETER :: r4=4
INTEGER,PUBLIC :: iDimLon=144
INTEGER,PUBLIC :: jDimLat=72
INTEGER,PUBLIC :: kDimLev=11
INTEGER,PUBLIC :: IX=480*100
INTEGER,PUBLIC :: IO=480
REAL   ,PUBLIC :: DeltaT=0.180E3
NAMELIST /RESOLUTION/ iDimLon,jDimLat,kDimLev,IX,IO,DeltaT

PUBLIC :: Init_SGCM
PUBLIC :: PARAM
PUBLIC :: BASIC,POLES,DIAGX,DIAGH,DIAGR,TPGEV,TENDZ,PGRAD,INTEL,Diffusion,Damping,FORCE,WRDAT
PUBLIC :: BSCPRN,DiagVertSigVelocity,ADVEC,RSGN,ZERO,INITC,OFLOPN,Convergence,CORIO,ADBTC,INTEE,DIAGS,DIAGZ
CONTAINS

SUBROUTINE Init_SGCM()
  IMPLICIT NONE
  OPEN(1,FILE='NameList',STATUS='OLD',ACTION='READ')
  READ(1,NML=RESOLUTION)
  PRINT*,iDimLon,jDimLat,kDimLev
END SUBROUTINE Init_SGCM


! ########################################
!   SUBROUTINES
! ########################################
! ========================================

 SUBROUTINE CONST (PI     ,&!REAL(KIND=r4), INTENT(OUT  ) :: PI 
                   RE     ,&!REAL(KIND=r4), INTENT(OUT  ) :: RE 
                   OM     ,&!REAL(KIND=r4), INTENT(OUT  ) :: OM 
                   GC     ,&!REAL(KIND=r4), INTENT(OUT  ) :: GC 
                   RC     ,&!REAL(KIND=r4), INTENT(OUT  ) :: RC 
                   T00    ,&!REAL(KIND=r4), INTENT(OUT  ) :: T00
                   P00     )!REAL(KIND=r4), INTENT(OUT  ) :: P00
   IMPLICIT NONE
   REAL(KIND=r4), INTENT(OUT  ) :: PI 
   REAL(KIND=r4), INTENT(OUT  ) :: RE 
   REAL(KIND=r4), INTENT(OUT  ) :: OM 
   REAL(KIND=r4), INTENT(OUT  ) :: GC 
   REAL(KIND=r4), INTENT(OUT  ) :: RC 
   REAL(KIND=r4), INTENT(OUT  ) :: T00
   REAL(KIND=r4), INTENT(OUT  ) :: P00

      PI  = 0.31415926535E1_r4   !geometric constant PI 
      RE  = 0.637E7_r4           !radius of the earth [meters]
      OM  = 0.727E-4_r4          ! 
      GC  = 0.9807E1_r4          ! gravity constat
      RC  = 0.2870E3_r4          ! gas constat 
      T00 = 0.298E3_r4
      P00 = 0.101325E6_r4

   RETURN
 END SUBROUTINE CONST 

! ========================================

  SUBROUTINE SIGMA (kDimLev     ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                    DZ          ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                    ilev        ,&!INTEGER      , INTENT(IN   ) :: ilev
                    si            )!REAL(KIND=r4), INTENT(OUT  ) :: si
    INTEGER      , INTENT(IN   ) :: kDimLev
    REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
    INTEGER      , INTENT(IN   ) :: ilev
    REAL(KIND=r4), INTENT(OUT  ) :: si
    INTEGER :: kk
      IF (ilev.GE.kDimLev+1) THEN
        si = 1.0E0
        RETURN
      ENDIF
      IF (ilev.LE.0) THEN
        si = 1.0E0
      DO kk=kDimLev,1,-1
        si = si - DZ(kk)
      END DO
        RETURN
      ENDIF
        si = 1.0E0
      DO kk=kDimLev,ilev+1,-1
        si = si - DZ(kk)
      END DO
        si = si - 0.5E0*DZ(ilev)
      RETURN
  END SUBROUTINE SIGMA


! ========================================
!   BASIC FIELDS
! ========================================

 SUBROUTINE PARAM (&
                   jDimLat      ,& !INTEGER      , INTENT(IN   ) :: kDimLev
                   kDimLev      ,& !INTEGER      , INTENT(IN   ) :: iDimLon
                   DX      ,& !REAL(KIND=r4), INTENT(OUT  ) :: DX(jDimLat)
                   DY      ,& !REAL(KIND=r4), INTENT(OUT  ) :: DY(jDimLat)
                   DZ      ,& !REAL(KIND=r4), INTENT(OUT  ) :: DZ(kDimLev)
                   PLEV0   ,& !REAL(KIND=r4), INTENT(OUT  ) :: PLEV0
                   PLEV    ,& !REAL(KIND=r4), INTENT(OUT  ) :: PLEV(kDimLev)
                   FC   )     !REAL(KIND=r4), INTENT(OUT  ) :: FC(jDimLat)!componente horizontal da força de Coriolis
  IMPLICIT NONE
  INTEGER      , INTENT(IN   ) :: jDimLat
  INTEGER      , INTENT(IN   ) :: kDimLev
  REAL(KIND=r4), INTENT(OUT  ) :: DX(jDimLat)
  REAL(KIND=r4), INTENT(OUT  ) :: DY(jDimLat)
  REAL(KIND=r4), INTENT(OUT  ) :: DZ(kDimLev)
  REAL(KIND=r4), INTENT(OUT  ) :: PLEV0
  REAL(KIND=r4), INTENT(OUT  ) :: PLEV(kDimLev)
  REAL(KIND=r4), INTENT(OUT  ) :: FC(jDimLat)!componente horizontal da força de Coriolis
  REAL(KIND=r4) :: PI  !geometric constant PI 
  REAL(KIND=r4) :: RE  !radius of the earth [meters]
  REAL(KIND=r4) :: OM  !angular velocity of the earth
  REAL(KIND=r4) :: GC
  REAL(KIND=r4) :: RC
  REAL(KIND=r4) :: T00
  REAL(KIND=r4) :: P00
  INTEGER :: M
  INTEGER :: N
  REAL(KIND=r4) :: Z1
  REAL(KIND=r4) :: Z2
  REAL(KIND=r4) :: R
  REAL(KIND=r4) :: PS
  REAL(KIND=r4) :: Z
  
  CALL CONST (PI,RE,OM,GC,RC,T00,P00)

  DO  M=1,jDimLat
      DX(M) = PI*RE/REAL(jDimLat) * SIN(PI*(REAL(M)-0.5E0)/REAL(jDimLat))
      DY(M) = PI*RE/REAL(jDimLat)
      FC(M) = - 0.2E1 * OM * COS(PI*(REAL(M)-0.5E0)/REAL(jDimLat))
  END DO
  Z1 = 0.100E0
  Z2 = 1.000E0
  DO  N=1,kDimLev
      R = (Z2/Z1)**(1.0E0/REAL(kDimLev))
      DZ(N) = Z1*(R-1.0E0)*R**(N-1)
  END DO
  PS = 0.101325E4
  Z = Z1
  PLEV0 = PS * Z1
  DO  N=1,kDimLev
      Z = Z + 0.5E0*DZ(N)
      PLEV(N) = PS * Z
      Z = Z + 0.5E0*DZ(N)
  END DO
  RETURN
  END  SUBROUTINE PARAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  SUBROUTINE BASIC (iDimLon     ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat     ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                    kDimLev     ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                    DZ     ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                    PLEV0  ,&!REAL(KIND=r4), INTENT(IN   ) :: PLEV0
                    PTOP   ,&!REAL(KIND=r4), INTENT(OUT  ) :: PTOP
                    U      ,&!REAL(KIND=r4), INTENT(OUT  ) :: U(iDimLon,jDimLat,kDimLev)
                    V      ,&!REAL(KIND=r4), INTENT(OUT  ) :: V(iDimLon,jDimLat,kDimLev)
                    T      ,&!REAL(KIND=r4), INTENT(OUT  ) :: T(iDimLon,jDimLat,kDimLev)
                    P      ,&!REAL(KIND=r4), INTENT(OUT  ) :: P(iDimLon,jDimLat)
                    H      ) !REAL(KIND=r4), INTENT(OUT  ) :: H(iDimLon,jDimLat)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: PLEV0
   REAL(KIND=r4), INTENT(OUT  ) :: PTOP
   REAL(KIND=r4), INTENT(OUT  ) :: U(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: V(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: T(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: P(iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(OUT  ) :: H(iDimLon,jDimLat)

   REAL(KIND=r4) ::  PDAT  (144,73,999)
   REAL(KIND=r4) ::  PLEVS (144,73,999)
   REAL(KIND=r4) ::  UDAT  (144,73,17)
   REAL(KIND=r4) ::  UU    (144,73,999)
   REAL(KIND=r4) ::  VDAT  (144,73,17)
   REAL(KIND=r4) ::  VV    (144,73,999)
   REAL(KIND=r4) ::  TDAT  (144,73,17)
   REAL(KIND=r4) ::  TT    (144,73,999)
   REAL(KIND=r4) ::  HDAT  (144,73,999)
   !REAL(KIND=r4) ::  HHDAT (144,73,999)

   REAL(KIND=r4), PARAMETER :: PSTAN(17)=(/1000.,925.,850.,700.,600.,500.,400.,300.,250.,200., &
                                               150.,100.,70.,50.,30.,20.,10./)
   CHARACTER INFL*80

   INTEGER :: M
   INTEGER :: N
   INTEGER :: L

   WRITE(6,*) 'Basic field data files are automatically designated:'
    101 FORMAT (1X,A18,1X,A9)
!C      INFL = 'TOPO.dat'
!C      WRITE(6,101) 'Topography       =',INFL
!C      CALL RDDAT (144,73, 1,INFL,HDAT)
        DO  M=1,73
           DO  L=1,144
               HDAT(L,M,1) = 0.0E0
           END DO
        END DO
!C      DO 13 M=1,73
!C      MA = MAX(M-1, 1)
!C      MB = MIN(M+1,73)
!C      DO 14 L=1,144
!C      LA = MOD(L-1+144-1,144)+1
!C      LB = MOD(L+1+144-1,144)+1
!C        HHDAT(L,M,1) =   0.5E0*HDAT(L,M,1)
!C     +               + 0.125E0*( HDAT(LA,MA,1)+HDAT(LB,MA,1)
!C     +                          +HDAT(LA,MB,1)+HDAT(LB,MB,1))
!C   14 CONTINUE
!C   13 CONTINUE
!C      DO 15 M=1,73
!C      DO 16 L=1,144
!C        HDAT(L,M,1) = HHDAT(L,M,1)
!C   16 CONTINUE
!C   15 CONTINUE
        INFL = 'U0.dat'
        WRITE(6,101) 'Zonal wind       =',INFL
        CALL RDDAT (144   ,&  !INTEGER         , INTENT(IN   ) :: iDimLon
                    73    ,&  !INTEGER         , INTENT(IN   ) :: jDimLat
                    17    ,&  !INTEGER         , INTENT(IN   ) :: NNX
                    INFL  ,&  !CHARACTER(LEN=*), INTENT(IN   ) :: INFL
                    UDAT   )  !REAL(KIND=r4)   , INTENT(OUT  ) :: G(iDimLon,jDimLat,NNX)
        INFL = 'V0.dat'
        WRITE(6,101) 'Meridional wind  =',INFL
        CALL RDDAT (144   ,&  !INTEGER         , INTENT(IN   ) :: iDimLon
                    73    ,&  !INTEGER         , INTENT(IN   ) :: jDimLat
                    17    ,&  !INTEGER         , INTENT(IN   ) :: NNX
                    INFL  ,&  !CHARACTER(LEN=*), INTENT(IN   ) :: INFL
                    VDAT   )  !REAL(KIND=r4)   , INTENT(OUT  ) :: G(iDimLon,jDimLat,NNX)
        INFL = 'T0.dat'
        WRITE(6,101) 'Temperature      =',INFL
        CALL RDDAT (144  ,&  !INTEGER         , INTENT(IN   ) :: iDimLon
                    73   ,&  !INTEGER         , INTENT(IN   ) :: jDimLat
                    17   ,&  !INTEGER         , INTENT(IN   ) :: NNX
                    INFL ,&  !CHARACTER(LEN=*), INTENT(IN   ) :: INFL
                    TDAT  )  !REAL(KIND=r4)   , INTENT(OUT  ) :: G(iDimLon,jDimLat,NNX)
        INFL = 'P0.dat'
        WRITE(6,101) 'Surface pressure =',INFL
        CALL RDDAT (144  ,&  !INTEGER         , INTENT(IN   ) :: iDimLon
                     73  ,&  !INTEGER         , INTENT(IN   ) :: jDimLat
                      1  ,&  !INTEGER         , INTENT(IN   ) :: NNX
                   INFL  ,&  !CHARACTER(LEN=*), INTENT(IN   ) :: INFL
                   PDAT(1:144,1:73,1)   )  !REAL(KIND=r4)   , INTENT(OUT  ) :: G(iDimLon,jDimLat,NNX)
        WRITE(6,*) '===================================='
        CALL THETA (144    , &  !INTEGER     , INTENT(IN   ) :: iDimLon
                     73    , &  !INTEGER     , INTENT(IN   ) :: jDimLat
                     17    , &  !INTEGER     , INTENT(IN   ) :: NX
                    PSTAN  , &  !REAL(KIND=r4), INTENT(IN   ) :: PSTAN(NX)
                    TDAT    )   !REAL(KIND=r4), INTENT(INOUT) :: T(iDimLon,jDimLat,NX)
        CALL PSMOD (144     ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                    73      ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                    17      ,&!INTEGER      , INTENT(IN   ) :: NX
                    PSTAN   ,&!REAL(KIND=r4), INTENT(IN   ) :: PSTAN(NX)
                    PDAT(1:144,1:73,1)    ,&!REAL(KIND=r4), INTENT(INOUT) :: P(iDimLon,jDimLat)
                    TDAT    ,&!REAL(KIND=r4), INTENT(IN   ) :: T(iDimLon,jDimLat,NX)
                    HDAT(1:144,1:73,1)     )!REAL(KIND=r4), INTENT(IN   ) :: H(iDimLon,jDimLat)
        CALL PCOOR (144     , & !INTEGER      , INTENT(IN   ) :: iDimLon
                    73      , & !INTEGER      , INTENT(IN   ) :: jDimLat
                    kDimLev      , & !INTEGER      , INTENT(IN   ) :: kDimLev
                    DZ      , & !REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                    PLEV0   , & !REAL(KIND=r4), INTENT(IN   ) :: PTOP
                    PDAT(1:144,1:73,1)    , & !REAL(KIND=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
                    PLEVS     ) !REAL(KIND=r4), INTENT(OUT  ) :: PLEVS(iDimLon,jDimLat,999)
        CALL LAYER (144     , & !INTEGER      , INTENT(IN   ) :: iDimLon
                    73      , & !INTEGER      , INTENT(IN   ) :: jDimLat
                    17      , & !INTEGER      , INTENT(IN   ) :: NNX
                    kDimLev      , & !INTEGER      , INTENT(IN   ) :: kDimLev
                    PSTAN   , & !REAL(KIND=r4), INTENT(IN   ) :: PSTAN(NNX)
                    PLEVS   , & !REAL(KIND=r4), INTENT(IN   ) :: PLEVS(iDimLon,jDimLat,999)
                    UDAT    , & !REAL(KIND=r4), INTENT(IN   ) :: GDAT(iDimLon,jDimLat,NNX)
                    UU        ) !REAL(KIND=r4), INTENT(OUT  ) :: G(iDimLon,jDimLat,999)
        CALL LAYER (144     , & !INTEGER      , INTENT(IN   ) :: iDimLon
                    73      , & !INTEGER      , INTENT(IN   ) :: jDimLat
                    17      , & !INTEGER      , INTENT(IN   ) :: NNX
                    kDimLev      , & !INTEGER      , INTENT(IN   ) :: kDimLev
                    PSTAN   , & !REAL(KIND=r4), INTENT(IN   ) :: PSTAN(NNX)
                    PLEVS   , & !REAL(KIND=r4), INTENT(IN   ) :: PLEVS(iDimLon,jDimLat,999)
                    VDAT    , & !REAL(KIND=r4), INTENT(IN   ) :: GDAT(iDimLon,jDimLat,NNX)
                    VV        ) !REAL(KIND=r4), INTENT(OUT  ) :: G(iDimLon,jDimLat,999)
        CALL LAYER (144     , & !INTEGER      , INTENT(IN   ) :: iDimLon
                    73      , & !INTEGER      , INTENT(IN   ) :: jDimLat
                    17      , & !INTEGER      , INTENT(IN   ) :: NNX
                    kDimLev      , & !INTEGER      , INTENT(IN   ) :: kDimLev
                    PSTAN   , & !REAL(KIND=r4), INTENT(IN   ) :: PSTAN(NNX)
                    PLEVS   , & !REAL(KIND=r4), INTENT(IN   ) :: PLEVS(iDimLon,jDimLat,999)
                    TDAT    , & !REAL(KIND=r4), INTENT(IN   ) :: GDAT(iDimLon,jDimLat,NNX)
                    TT        ) !REAL(KIND=r4), INTENT(OUT  ) :: G(iDimLon,jDimLat,999)
        CALL RDCNVU (iDimLon     , & !INTEGER      , INTENT(IN   ) ::  iDimLon
                     jDimLat     , & !INTEGER      , INTENT(IN   ) ::  jDimLat
                     kDimLev     , & !INTEGER      , INTENT(IN   ) ::  kDimLev
                     UU     , & !REAL(KINd=r4), INTENT(IN   ) :: G(144,73,999)
                     U      )   !REAL(KINd=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat,kDimLev)
        CALL RDCNVV (iDimLon     , & !INTEGER     , INTENT(IN   ) :: iDimLon
                     jDimLat     , & !INTEGER     , INTENT(IN   ) :: jDimLat
                     kDimLev     , & !INTEGER     , INTENT(IN   ) :: kDimLev
                     VV     , & !REAL(KIND=r4), INTENT(IN   ) :: G(144,73,999)
                     V       )  !REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat,kDimLev)
        CALL RDCNVG (iDimLon     , & !INTEGER      , INTENT(IN   ) :: iDimLon
                     jDimLat     , & !INTEGER      , INTENT(IN   ) :: jDimLat
                     kDimLev     , & !INTEGER      , INTENT(IN   ) :: kDimLev
                     TT     , & !REAL(KIND=r4), INTENT(IN   ) :: G(144,73,999)
                     T       )  !REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat,kDimLev)
        CALL RDCNVG (iDimLon      ,& !INTEGER      , INTENT(IN   ) :: iDimLon
                     jDimLat      ,& !INTEGER      , INTENT(IN   ) :: jDimLat
                     1       ,& !INTEGER      , INTENT(IN   ) :: kDimLev
                     PDAT(1:144,1:73,1)    ,& !REAL(KIND=r4), INTENT(IN   ) :: G(144,73,999)
                     P        ) !REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat,kDimLev)
        CALL RDCNVG (iDimLon       ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                     jDimLat       ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                     1        ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                     HDAT(1:144,1:73,1)     ,&!REAL(KIND=r4), INTENT(IN   ) :: G(144,73,999)
                     H         )!REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat,kDimLev)

        PTOP = 1.0E2*PLEV0

        DO M=1,jDimLat
           DO L=1,iDimLon
              P(L,M) = 1.0E2*P(L,M)
           END DO
        END DO
        DO N=1,kDimLev
           DO L=1,iDimLon
              V(L,1,N) = 0.0E0
           END DO
        END DO
   RETURN
  END  SUBROUTINE BASIC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  SUBROUTINE RDDAT (iDimLon       ,& !INTEGER         , INTENT(IN   ) :: iDimLon
                    jDimLat       ,& !INTEGER         , INTENT(IN   ) :: jDimLat
                    NNX      ,& !INTEGER         , INTENT(IN   ) :: NNX
                    INFL     ,& !CHARACTER(LEN=*), INTENT(IN   ) :: INFL
                    G        )  !REAL(KIND=r4)   , INTENT(OUT  ) :: G(iDimLon,jDimLat,NNX)
    IMPLICIT NONE
    INTEGER         , INTENT(IN   ) :: iDimLon
    INTEGER         , INTENT(IN   ) :: jDimLat
    INTEGER         , INTENT(IN   ) :: NNX
    CHARACTER(LEN=*), INTENT(IN   ) :: INFL
    REAL(KIND=r4)   , INTENT(OUT) :: G(iDimLon,jDimLat,NNX)
    INTEGER :: lrec
    INQUIRE(IOLENGTH=lrec) G
    PRINT*,lrec
    OPEN(10,FILE=TRIM(INFL),STATUS='OLD', &
          FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec)
    READ(10,REC=1) G
    PRINT*,'MAXVAL(G)=',MAXVAL(G),'MINVAL(G)=',MINVAL(G)
    CLOSE(10)
    RETURN
  END  SUBROUTINE RDDAT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


  SUBROUTINE THETA (iDimLon     ,& !INTEGER     , INTENT(IN   ) :: iDimLon
                    jDimLat     ,& !INTEGER     , INTENT(IN   ) :: jDimLat
                    kDimLev     ,& !INTEGER     , INTENT(IN   ) :: kDimLev
                    PSTAN  ,& !REAL(KIND=r4), INTENT(IN   ) :: PSTAN(kDimLev)
                    T       ) !REAL(KIND=r4), INTENT(INOUT) :: T(iDimLon,jDimLat,kDimLev)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: PSTAN(kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: T(iDimLon,jDimLat,kDimLev)
   INTEGER :: k
   INTEGER :: j
   INTEGER :: i
   DO k=1,kDimLev
      DO j=1,jDimLat
         DO i=1,iDimLon
            T(i,j,k) = T(i,j,k) * (1.0E3/PSTAN(k))**(0.2E1/0.7E1)
         END DO
      END DO
    END DO
   RETURN
  END SUBROUTINE THETA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  SUBROUTINE PSMOD (iDimLon    ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat    ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                    kDimLev    ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                    PSTAN ,&!REAL(KIND=r4), INTENT(IN   ) :: PSTAN(kDimLev)
                    P     ,&!REAL(KIND=r4), INTENT(INOUT) :: P(iDimLon,jDimLat)
                    T     ,&!REAL(KIND=r4), INTENT(IN   ) :: T(iDimLon,jDimLat,kDimLev)
                    H      )!REAL(KIND=r4), INTENT(IN   ) :: H(iDimLon,jDimLat)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: PSTAN(kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: P(iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: T(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: H(iDimLon,jDimLat)

   INTEGER :: M
   INTEGER :: L
   INTEGER :: N

   REAL(KIND=r4) :: PR
   REAL(KIND=r4) :: DPR
   REAL(KIND=r4) :: TM
   REAL(KIND=r4) :: Z
   REAL(KIND=r4) :: Z1
   REAL(KIND=r4) :: PI  !geometric constant PI 
   REAL(KIND=r4) :: RE  !radius of the earth [meters]
   REAL(KIND=r4) :: OM  !angular velocity of the earth
   REAL(KIND=r4) :: GC  ! gravity constat
   REAL(KIND=r4) :: RC  ! gas constat 
   REAL(KIND=r4) :: T00
   REAL(KIND=r4) :: P00 

      CALL CONST (PI  ,&!REAL(KIND=r4), INTENT(OUT  ) :: PI 
                  RE  ,&!REAL(KIND=r4), INTENT(OUT  ) :: RE 
                  OM  ,&!REAL(KIND=r4), INTENT(OUT  ) :: OM 
                  GC  ,&!REAL(KIND=r4), INTENT(OUT  ) :: GC 
                  RC  ,&!REAL(KIND=r4), INTENT(OUT  ) :: RC 
                  T00 ,&!REAL(KIND=r4), INTENT(OUT  ) :: T00
                  P00  )!REAL(KIND=r4), INTENT(OUT  ) :: P00
   DO M=1,jDimLat
      DO L=1,iDimLon

         PR = 1.0E2*P(L,M)
         DPR = 1.0E2*(P(L,M)-PSTAN(1))
         TM = T(L,M,1) * (PR/1.0E5)**(0.2E1/0.7E1)
         Z = RC*TM/GC * DPR / PR
         IF (Z  .GE.  H(L,M)) THEN
            P(L,M) = ((Z-H(L,M))*P(L,M)+H(L,M)*PSTAN(1))/ Z
         ELSE
            DO N=1,kDimLev-1
               PR = 1.0E2*0.5E0*(PSTAN(N)+PSTAN(N+1))
               DPR = 1.0E2*(PSTAN(N)-PSTAN(N+1))
               TM = T(L,M,N) * (PR/1.0E5)**(0.2E1/0.7E1)
               Z1 = Z
               Z = Z + RC*TM/GC * DPR / PR
               IF (Z.GE.H(L,M)) THEN
                  P(L,M) = ((Z-H(L,M))*PSTAN(N)+(H(L,M)-Z1)*PSTAN(N+1))/(Z-Z1)
                  GO TO 29
               ENDIF
               P(L,M) = PSTAN(1)
            END DO
         ENDIF
         29 CONTINUE
      END DO
   END DO
   RETURN
  END SUBROUTINE PSMOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE PCOOR (iDimLon     , & !INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat     , & !INTEGER      , INTENT(IN   ) :: jDimLat
                    kDimLev     , & !INTEGER      , INTENT(IN   ) :: kDimLev
                    DZ     , & !REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                    PTOP   , & !REAL(KIND=r4), INTENT(IN   ) :: PTOP
                    P      , & !REAL(KIND=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
                    PLEVS    ) !REAL(KIND=r4), INTENT(OUT  ) :: PLEVS(iDimLon,jDimLat,999)

   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: PTOP
   REAL(KIND=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(OUT  ) :: PLEVS(iDimLon,jDimLat,999)
   REAL(KIND=r4) :: Z0
   REAL(KIND=r4) :: Z
   INTEGER :: N
   INTEGER :: NN
   INTEGER :: M
   INTEGER :: L
   
   Z0 = 0.0E0
   DO N=1,kDimLev
      Z0=Z0+DZ(N)
   END DO
   DO N=1,kDimLev 
        Z=0.0E0
      DO NN=1,N-1
         Z=Z+DZ(NN)/Z0
      END DO
      Z=Z+0.5E0*DZ(N)/Z0
      DO M=1,jDimLat
         DO L=1,iDimLon
            PLEVS(L,M,N) = (1.0E0-Z)*PTOP + Z*P(L,M)
         END DO
      END DO
   END DO
   RETURN
  END SUBROUTINE PCOOR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE LAYER (iDimLon       ,&!  INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat       ,&!  INTEGER      , INTENT(IN   ) :: jDimLat
                    NNX      ,&!  INTEGER      , INTENT(IN   ) :: NNX
                    kDimLev       ,&!  INTEGER      , INTENT(IN   ) :: kDimLev
                    PSTAN    ,&!  REAL(KIND=r4), INTENT(IN   ) :: PSTAN(NNX)
                    PLEVS    ,&!  REAL(KIND=r4), INTENT(IN   ) :: PLEVS(iDimLon,jDimLat,999)
                    GDAT     ,&!  REAL(KIND=r4), INTENT(IN   ) :: GDAT(iDimLon,jDimLat,NNX)
                    G         )!  REAL(KIND=r4), INTENT(OUT  ) :: G(iDimLon,jDimLat,999)
    INTEGER      , INTENT(IN   ) :: iDimLon
    INTEGER      , INTENT(IN   ) :: jDimLat
    INTEGER      , INTENT(IN   ) :: NNX
    INTEGER      , INTENT(IN   ) :: kDimLev
    REAL(KIND=r4), INTENT(IN   ) :: PSTAN(NNX)
    REAL(KIND=r4), INTENT(IN   ) :: PLEVS(iDimLon,jDimLat,999)
    REAL(KIND=r4), INTENT(IN   ) :: GDAT(iDimLon,jDimLat,NNX)
    REAL(KIND=r4), INTENT(OUT  ) :: G(iDimLon,jDimLat,999)

    REAL(KIND=r4) :: C1
    REAL(KIND=r4) :: C2
    REAL(KIND=r4) :: P
    REAL(KIND=r4) :: P1
    REAL(KIND=r4) :: P2
    
    
    INTEGER :: M
    INTEGER :: L
    INTEGER :: N
    INTEGER :: NN
    INTEGER :: NN1
    INTEGER :: NN2
    
    DO M=1,jDimLat
       DO L=1,iDimLon
          DO N=1,kDimLev
             P = PLEVS(L,M,N)
             DO NN=NNX-1,1,-1
                P1 = PSTAN (NN+1)
                P2 = PSTAN (NN)
                IF ((P.LE.P2).OR.(NN.EQ.1)) THEN
                   NN1 = NN+1
                   NN2 = NN
                   C1 = (P2-P)/(P2-P1)
                   C2 = (P-P1)/(P2-P1)
                   GO TO 39
                ENDIF
             END DO
            
            39 CONTINUE

            G(L,M,N) = C1*GDAT(L,M,NN1) + C2*GDAT(L,M,NN2)

         END DO
      END DO
   END DO
   RETURN
  END SUBROUTINE LAYER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE RDCNVU (iDimLon      ,&!INTEGER      , INTENT(IN   ) ::  iDimLon
                     jDimLat      ,&!INTEGER      , INTENT(IN   ) ::  jDimLat
                     kDimLev      ,&!INTEGER      , INTENT(IN   ) ::  kDimLev
                     G            ,&!REAL(KINd=r4), INTENT(IN   ) :: G(144,73,999)
                     GG            )!REAL(KINd=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat,kDimLev)
   INTEGER      , INTENT(IN   ) ::  iDimLon
   INTEGER      , INTENT(IN   ) ::  jDimLat
   INTEGER      , INTENT(IN   ) ::  kDimLev
   REAL(KINd=r4), INTENT(IN   ) :: G(144,73,999)
   REAL(KINd=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat,kDimLev)
   INTEGER :: k
   INTEGER :: j
   INTEGER :: i
   INTEGER :: yb
   INTEGER :: yf
   INTEGER :: nlat
   nlat=73
   DO k=1,kDimLev
      DO j=1,jDimLat
         yb = nlat+0-j
         yf = nlat+1-j
         DO i=1,iDimLon
            GG(i,j,k) = 0.25E0*(G(i,yb,k)+G(i,yb,k)+G(i,yf,k)+G(i,yf,k))
         END DO
      END DO
   END DO
   RETURN
  END SUBROUTINE RDCNVU


  SUBROUTINE RDCNVV (iDimLon    ,& !INTEGER     , INTENT(IN   ) :: iDimLon
                     jDimLat    ,& !INTEGER     , INTENT(IN   ) :: jDimLat
                     kDimLev    ,& !INTEGER     , INTENT(IN   ) :: kDimLev
                     G     ,& !REAL(KIND=r4), INTENT(IN   ) :: G(144,73,999)
                     GG     ) !REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat,kDimLev)
    INTEGER      , INTENT(IN   ) :: iDimLon
    INTEGER      , INTENT(IN   ) :: jDimLat
    INTEGER      , INTENT(IN   ) :: kDimLev
    REAL(KIND=r4), INTENT(IN   ) :: G(144,73,999)
    REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat,kDimLev)
    INTEGER   :: N
    INTEGER   :: M
    INTEGER   :: L
    INTEGER   :: M1
    INTEGER   :: M2
    INTEGER   :: L1
    INTEGER   :: L2
    DO N=1,kDimLev
       DO M=1,jDimLat
          M1 = 74-M
          M2 = 74-M
          DO L=1,iDimLon
             L1 = L
             L2 = L+1
             IF (L2 .GT. 144) L2 = L2-144
             GG(L,M,N) = 0.25E0*(G(L1,M1,N)+G(L2,M1,N)+G(L1,M2,N)+G(L2,M2,N))
          END DO
       END DO
    END DO
    RETURN
  END SUBROUTINE RDCNVV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE RDCNVG (iDimLon    ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                     jDimLat    ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                     kDimLev    ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                     G     ,&! REAL(KIND=r4), INTENT(IN   ) :: G(144,73,999)
                     GG     )! REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat,kDimLev)
    INTEGER      , INTENT(IN   ) :: iDimLon
    INTEGER      , INTENT(IN   ) :: jDimLat
    INTEGER      , INTENT(IN   ) :: kDimLev
    REAL(KIND=r4), INTENT(IN   ) :: G(144,73,999)
    REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat,kDimLev)
    INTEGER :: k
    INTEGER :: j
    INTEGER :: i
    INTEGER :: j1
    INTEGER :: j2
    INTEGER :: i1
    INTEGER :: i2
    DO k=1,kDimLev
       DO j=1,jDimLat
          j1 = 73-j
          j2 = 74-j
          DO i=1,iDimLon
             i1 = i
             i2 = i+1
             IF (i2.GT.144) i2 = i2-144
             GG(i,j,k) = 0.25E0*(G(i1,j1,k)+G(i2,j1,k)+G(i1,j2,k)+G(i2,j2,k))
          END DO
       END DO
    END DO
      RETURN
  END SUBROUTINE RDCNVG

! ========================================
!   BOUNDARY CONDITIONS
! ========================================

  SUBROUTINE POLES (iDimLon   ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat   ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                    kDimLev   ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                    U    ,&! REAL(KIND=r4), INTENT(INOUT) :: U(iDimLon,jDimLat,kDimLev)
                    V    ,&! REAL(KIND=r4), INTENT(INOUT) :: V(iDimLon,jDimLat,kDimLev)
                    T    ,&! REAL(KIND=r4), INTENT(INOUT) :: T(iDimLon,jDimLat,kDimLev)
                    P    ,&! REAL(KIND=r4), INTENT(INOUT) :: P(iDimLon,jDimLat)
                    H     )! REAL(KIND=r4), INTENT(INOUT) :: H(iDimLon,jDimLat)
      INTEGER      , INTENT(IN   ) :: iDimLon
      INTEGER      , INTENT(IN   ) :: jDimLat
      INTEGER      , INTENT(IN   ) :: kDimLev
      REAL(KIND=r4), INTENT(INOUT) :: U(iDimLon,jDimLat,kDimLev)
      REAL(KIND=r4), INTENT(INOUT) :: V(iDimLon,jDimLat,kDimLev)
      REAL(KIND=r4), INTENT(INOUT) :: T(iDimLon,jDimLat,kDimLev)
      REAL(KIND=r4), INTENT(INOUT) :: P(iDimLon,jDimLat)
      REAL(KIND=r4), INTENT(INOUT) :: H(iDimLon,jDimLat)
      CALL POLESG (iDimLon ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                   kDimLev ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                   U   )!REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
      CALL POLESV (iDimLon ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                   kDimLev ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                   V   )!REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
      CALL POLESG (iDimLon  ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat  ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                   kDimLev  ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                   T    )!REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
      CALL POLESG (iDimLon    ,&!INTEGER    , INTENT(IN   ) :: iDimLon
                   jDimLat    ,&!INTEGER    , INTENT(IN   ) :: jDimLat
                   1     ,&!INTEGER    , INTENT(IN   ) :: kDimLev
                   P      )!REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
      CALL POLESG (iDimLon    ,&!INTEGER    , INTENT(IN   ) :: iDimLon
                   jDimLat    ,&!INTEGER    , INTENT(IN   ) :: jDimLat
                   1     ,&!INTEGER    , INTENT(IN   ) :: kDimLev
                   H      )!REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
      RETURN
  END SUBROUTINE POLES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE POLESG (iDimLon   ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                     jDimLat   ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                     kDimLev   ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                     G     )!REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
      INTEGER      , INTENT(IN   ) :: iDimLon
      INTEGER      , INTENT(IN   ) :: jDimLat
      INTEGER      , INTENT(IN   ) :: kDimLev
      REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
      INTEGER :: M
      INTEGER :: M1
      INTEGER :: M2
      INTEGER :: M3
      INTEGER :: M4
      INTEGER :: M5

      M1 = 12
      M2 = 6
      M3 = 3
      M4 = 2
      M5 = 1

   !  north polo

      !  7  ,12  
   DO M=M2+1,M1
      CALL POLESA (iDimLon    ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat    ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                   kDimLev    ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                   M          ,&!INTEGER      , INTENT(IN   ) :: yc
                   2          ,&!INTEGER      , INTENT(IN   ) :: LR --> (iDimLon-1)/LR+1
                   G           )!REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
   END DO
      !  4  ,6  
   DO M=M3+1,M2
      CALL POLESA (iDimLon    ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat    ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                   kDimLev    ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                   M          ,&!INTEGER      , INTENT(IN   ) :: yc
                   4          ,&!INTEGER      , INTENT(IN   ) :: LR --> (iDimLon-1)/LR+1
                   G           )!REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
   END DO
      !  3  ,3  
   DO M=M4+1,M3
      CALL POLESA (iDimLon    ,&!INTEGER     , INTENT(IN   ) :: iDimLon
                   jDimLat    ,&!INTEGER     , INTENT(IN   ) :: jDimLat
                   kDimLev    ,&!INTEGER     , INTENT(IN   ) :: kDimLev
                   M          ,&!INTEGER     , INTENT(IN   ) :: yc
                   8          ,&!INTEGER     , INTENT(IN   ) :: LR --> (iDimLon-1)/LR+1
                   G           )!REAL(KIND=r4),INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
   END DO
      !  2  ,2  
   DO M=M5+1,M4
      CALL POLESA (iDimLon    ,&!INTEGER     , INTENT(IN   ) :: iDimLon
                   jDimLat    ,&!INTEGER     , INTENT(IN   ) :: jDimLat
                   kDimLev    ,&!INTEGER     , INTENT(IN   ) :: kDimLev
                   M          ,&!INTEGER     , INTENT(IN   ) :: yc
                   16         ,&!INTEGER     , INTENT(IN   ) :: LR --> (iDimLon-1)/LR+1
                   G           )!REAL(KIND=r4),INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
   END DO
      ! 1,1
   DO M=1,1
      CALL POLESA (iDimLon    ,&!INTEGER    , INTENT(IN   ) :: iDimLon
                   jDimLat    ,&!INTEGER    , INTENT(IN   ) :: jDimLat
                   kDimLev    ,&!INTEGER    , INTENT(IN   ) :: kDimLev
                   M          ,&!INTEGER    , INTENT(IN   ) :: yc
                   iDimLon    ,&!INTEGER    , INTENT(IN   ) :: LR --> (iDimLon-1)/LR+1
                   G           )!REAL(KIND=r4), INTENT(INOUT) ::  G(iDimLon,jDimLat,kDimLev)
   END DO
   !  north polo

   !    jDimLat-11,jDimLat-6
   DO M=jDimLat-M1+1,jDimLat-M2
      CALL POLESA (iDimLon    ,&!INTEGER     , INTENT(IN   ) :: iDimLon
                   jDimLat    ,&!INTEGER     , INTENT(IN   ) :: jDimLat
                   kDimLev    ,&!INTEGER     , INTENT(IN   ) :: kDimLev
                   M          ,&!INTEGER     , INTENT(IN   ) :: yc
                   2          ,&!INTEGER     , INTENT(IN   ) :: LR --> (iDimLon-1)/LR+1
                   G           )!REAL(KIND=r4), INTENT(INOUT) ::  G(iDimLon,jDimLat,kDimLev)
   END DO
   !    jDimLat-5,jDimLat-3
   DO M=jDimLat-M2+1,jDimLat-M3
      CALL POLESA (iDimLon    ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat    ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                   kDimLev    ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                   M          ,&!INTEGER      , INTENT(IN   ) :: yc
                   4          ,&!INTEGER      , INTENT(IN   ) ::  LR --> (iDimLon-1)/LR+1
                   G           )!REAL(KIND=r4), INTENT(INOUT) ::  G(iDimLon,jDimLat,kDimLev)
   END DO
   !    jDimLat-2,jDimLat-2
   DO M=jDimLat-M3+1,jDimLat-M4
      CALL POLESA (iDimLon    ,&!INTEGER     , INTENT(IN   ) :: iDimLon
                   jDimLat    ,&!INTEGER     , INTENT(IN   ) :: jDimLat
                   kDimLev    ,&!INTEGER     , INTENT(IN   ) :: kDimLev
                   M          ,&!INTEGER     , INTENT(IN   ) :: yc
                   8          ,&!INTEGER     , INTENT(IN   ) ::  LR --> (iDimLon-1)/LR+1
                   G           )!REAL(KIND=r4), INTENT(INOUT) ::  G(iDimLon,jDimLat,kDimLev)
   END DO
   !    jDimLat-1,jDimLat-1
   DO M=jDimLat-M4+1,jDimLat-M5
      CALL POLESA (iDimLon    ,&!INTEGER    , INTENT(IN   ) :: iDimLon
                   jDimLat    ,&!INTEGER    , INTENT(IN   ) :: jDimLat
                   kDimLev    ,&!INTEGER    , INTENT(IN   ) :: kDimLev
                   M          ,&!INTEGER    , INTENT(IN   ) :: yc
                   16         ,&!INTEGER    , INTENT(IN   ) ::  LR --> (iDimLon-1)/LR+1
                   G           )!REAL(KIND=r4), INTENT(INOUT) ::  G(iDimLon,jDimLat,kDimLev)
   END DO
   !    jDimLat,jDimLat
   DO M=jDimLat,jDimLat
      CALL POLESA (iDimLon   ,&!INTEGER    , INTENT(IN   ) :: iDimLon
                   jDimLat   ,&!INTEGER    , INTENT(IN   ) :: jDimLat
                   kDimLev   ,&!INTEGER    , INTENT(IN   ) :: kDimLev
                   M         ,&!INTEGER    , INTENT(IN   ) :: yc
                   iDimLon   ,&!INTEGER    , INTENT(IN   ) ::  LR --> (iDimLon-1)/LR+1
                   G          )!REAL(KIND=r4), INTENT(INOUT) ::  G(iDimLon,jDimLat,kDimLev)
   END DO
      RETURN
  END SUBROUTINE POLESG 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE POLESV (iDimLon,jDimLat,kDimLev,G)
      INTEGER      , INTENT(IN   ) :: iDimLon
      INTEGER      , INTENT(IN   ) :: jDimLat
      INTEGER      , INTENT(IN   ) :: kDimLev
      REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
      INTEGER :: M
      INTEGER :: M1
      INTEGER :: M2
      INTEGER :: M3
      INTEGER :: M4
      INTEGER :: M5

      M1 = 12
      M2 = 6
      M3 = 3
      M4 = 2
      M5 = 1
      DO M=M2+1,M1
        CALL POLESA (iDimLon    ,&!INTEGER     , INTENT(IN   ) :: iDimLon
                     jDimLat    ,&!INTEGER     , INTENT(IN   ) :: jDimLat
                     kDimLev    ,&!INTEGER     , INTENT(IN   ) :: kDimLev
                     M     ,&!INTEGER     , INTENT(IN   ) :: M
                     2     ,&!INTEGER     , INTENT(IN   ) :: LR
                     G      )!REAL(KIND=r4), INTENT(INOUT) ::  G(iDimLon,jDimLat,kDimLev)
      END DO
      DO M=M3+1,M2
        CALL POLESA (iDimLon    ,&!INTEGER     , INTENT(IN   ) :: iDimLon
                     jDimLat    ,&!INTEGER     , INTENT(IN   ) :: jDimLat
                     kDimLev    ,&!INTEGER     , INTENT(IN   ) :: kDimLev
                     M     ,&!INTEGER     , INTENT(IN   ) :: M
                     4     ,&!INTEGER     , INTENT(IN   ) :: LR
                     G      )!REAL(KIND=r4), INTENT(INOUT) ::  G(iDimLon,jDimLat,kDimLev)
      END DO
      DO M=M4+1,M3
        CALL POLESA (iDimLon    ,&!INTEGER     , INTENT(IN   ) :: iDimLon
                     jDimLat    ,&!INTEGER     , INTENT(IN   ) :: jDimLat
                     kDimLev    ,&!INTEGER     , INTENT(IN   ) :: kDimLev
                     M     ,&!INTEGER     , INTENT(IN   ) :: M
                     8     ,&!INTEGER     , INTENT(IN   ) :: LR
                     G      )!REAL(KIND=r4), INTENT(INOUT) ::  G(iDimLon,jDimLat,kDimLev)
      END DO
      DO M=M5+1,M4
        CALL POLESA (iDimLon    ,&!INTEGER     , INTENT(IN   ) :: iDimLon
                     jDimLat    ,&!INTEGER     , INTENT(IN   ) :: jDimLat
                     kDimLev    ,&!INTEGER     , INTENT(IN   ) :: kDimLev
                     M     ,&!INTEGER     , INTENT(IN   ) :: M
                     16    ,&!INTEGER     , INTENT(IN   ) :: LR
                     G      )!REAL(KIND=r4), INTENT(INOUT) ::  G(iDimLon,jDimLat,kDimLev)
      END DO
      DO M=jDimLat-M1+1,jDimLat-M2
        CALL POLESA (iDimLon    ,&!INTEGER     , INTENT(IN   ) :: iDimLon
                     jDimLat    ,&!INTEGER     , INTENT(IN   ) :: jDimLat
                     kDimLev    ,&!INTEGER     , INTENT(IN   ) :: kDimLev
                     M+1   ,&!INTEGER     , INTENT(IN   ) :: M
                     2     ,&!INTEGER     , INTENT(IN   ) :: LR
                     G      )!REAL(KIND=r4), INTENT(INOUT) ::  G(iDimLon,jDimLat,kDimLev)
      END DO
      DO M=jDimLat-M2+1,jDimLat-M3
        CALL POLESA (iDimLon     ,&!INTEGER      , INTENT(IN    ) :: iDimLon
                     jDimLat     ,&!INTEGER      , INTENT(IN    ) :: jDimLat
                     kDimLev     ,&!INTEGER      , INTENT(IN    ) :: kDimLev
                     M+1    ,&!INTEGER      , INTENT(IN    ) :: M
                     4      ,&!INTEGER      , INTENT(IN    ) :: LR
                     G       )!REAL(KIND=r4), INTENT(INOUT) ::  G(iDimLon,jDimLat,kDimLev)
      END DO
      DO M=jDimLat-M3+1,jDimLat-M4
        CALL POLESA (iDimLon     ,&!INTEGER      , INTENT(IN    ) :: iDimLon
                     jDimLat     ,&!INTEGER      , INTENT(IN    ) :: jDimLat
                     kDimLev     ,&!INTEGER      , INTENT(IN    ) :: kDimLev
                     M+1    ,&!INTEGER      , INTENT(IN    ) :: M
                     8      ,&!INTEGER      , INTENT(IN    ) :: LR
                     G       )!REAL(KIND=r4), INTENT(INOUT) ::  G(iDimLon,jDimLat,kDimLev)
      END DO
      DO M=jDimLat-M4+1,jDimLat-M5
        CALL POLESA (iDimLon   ,&!INTEGER      , INTENT(IN    ) :: iDimLon
                     jDimLat   ,&!INTEGER      , INTENT(IN    ) :: jDimLat
                     kDimLev   ,&!INTEGER      , INTENT(IN    ) :: kDimLev
                     M+1  ,&!INTEGER      , INTENT(IN    ) :: M
                     16   ,&!INTEGER      , INTENT(IN    ) :: LR
                     G       )!REAL(KIND=r4), INTENT(INOUT) ::  G(iDimLon,jDimLat,kDimLev)
      END DO
      RETURN
  END SUBROUTINE POLESV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE POLESA (iDimLon        ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                     jDimLat        ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                     kDimLev        ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                     yc             ,&!INTEGER      , INTENT(IN   ) :: yc
                     LR             ,&!INTEGER      , INTENT(IN   ) :: LR
                     G               )!REAL(KIND=r4), INTENT(INOUT) ::  G(iDimLon,jDimLat,kDimLev)
    INTEGER      , INTENT(IN   ) :: iDimLon
    INTEGER      , INTENT(IN   ) :: jDimLat
    INTEGER      , INTENT(IN   ) :: kDimLev
    INTEGER      , INTENT(IN   ) :: yc
    INTEGER      , INTENT(IN   ) :: LR
    REAL(KIND=r4), INTENT(INOUT) ::  G(iDimLon,jDimLat,kDimLev)

    INTEGER       :: ii1
    INTEGER       :: ii2
    REAL(KIND=r4) :: SUM
    INTEGER       :: i,xc,k
    
    DO k=1,kDimLev
       DO i=1,(iDimLon-1)/LR+1

          ii1 = 1+LR*(i-1)
          ii2 = MIN0(LR+LR*(i-1),iDimLon)
          SUM = 0.0E0

          DO xc=ii1,ii2
             SUM = SUM + G(xc,yc,k)
          END DO

          DO xc=ii1,ii2
             G(xc,yc,k) = SUM / REAL(ii2-ii1+1)
          END DO

       END DO

    END DO
    RETURN
  END  SUBROUTINE POLESA 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE DIAGX (iDimLon     ,&!  INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat     ,&!  INTEGER      , INTENT(IN   ) :: jDimLat
                    kDimLev     ,&!  INTEGER      , INTENT(IN   ) :: kDimLev
                    DZ     ,&!  REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                    PTOP   ,&!  REAL(KIND=r4), INTENT(IN   ) :: PTOP
                    P      ,&!  REAL(KIND=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
                    EX      )!  REAL(KIND=r4), INTENT(OUT  ) :: EX(iDimLon,jDimLat,kDimLev)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: PTOP
   REAL(KIND=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(OUT  ) :: EX(iDimLon,jDimLat,kDimLev)

   REAL(KIND=r4) :: PI 
   REAL(KIND=r4) :: RE 
   REAL(KIND=r4) :: OM 
   REAL(KIND=r4) :: GC 
   REAL(KIND=r4) :: RC 
   REAL(KIND=r4) :: T00
   REAL(KIND=r4) :: P00
   REAL(KIND=r4) :: Z
   REAL(KIND=r4) :: Z1
   REAL(KIND=r4) :: Z0
   REAL(KIND=r4) :: PR
   REAL(KIND=r4) :: Z1Z0I
   INTEGER       :: M
   INTEGER       :: N
   INTEGER       :: L
   CALL CONST (PI   ,&!REAL(KIND=r4), INTENT(OUT  ) :: PI 
               RE   ,&!REAL(KIND=r4), INTENT(OUT  ) :: RE 
               OM   ,&!REAL(KIND=r4), INTENT(OUT  ) :: OM 
               GC   ,&!REAL(KIND=r4), INTENT(OUT  ) :: GC 
               RC   ,&!REAL(KIND=r4), INTENT(OUT  ) :: RC 
               T00  ,&!REAL(KIND=r4), INTENT(OUT  ) :: T00
               P00   )!REAL(KIND=r4), INTENT(OUT  ) :: P00
!  P : Surface pressure.
!  EX: Exner function (ratio of temperature to potential
!      temperature). 
   CALL SIGMA (kDimLev    ,&!INTEGER      , INTENT(IN   ) :: kDimLev
               DZ    ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
               kDimLev+1  ,&!INTEGER      , INTENT(IN   ) :: N
               Z1     )!REAL(KIND=r4), INTENT(OUT  ) :: Z

   CALL SIGMA (kDimLev    ,&!INTEGER      , INTENT(IN   ) :: kDimLev
               DZ    ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
               0     ,&!INTEGER      , INTENT(IN   ) :: N
               Z0     )!REAL(KIND=r4), INTENT(OUT  ) :: Z
   Z1Z0I = 1.0E0/(Z1-Z0)
   DO M=1,jDimLat
      DO L=1,iDimLon
         DO N=kDimLev,1,-1
            CALL SIGMA (kDimLev    ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                        DZ    ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                        N     ,&!INTEGER      , INTENT(IN   ) :: N
                        Z      )!REAL(KIND=r4), INTENT(OUT  ) :: Z

            PR = (P(L,M)*(Z-Z0)+PTOP*(Z1-Z))*Z1Z0I
            EX(L,M,N) = (PR/1.0E5)**(0.2E1/0.7E1)
         END DO
      END DO
   END DO
   RETURN
  END SUBROUTINE DIAGX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DIAGH (iDimLon     ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat     ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                    kDimLev     ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                    DZ     ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ (kDimLev)
                    TPG    ,&!REAL(KIND=r4), INTENT(IN   ) :: TPG(iDimLon,jDimLat)
                    PTOP   ,&!REAL(KIND=r4), INTENT(IN   ) :: PTOP
                    P      ,&!REAL(KIND=r4), INTENT(IN   ) :: P  (iDimLon,jDimLat)
                    EX     ,&!REAL(KIND=r4), INTENT(IN   ) :: EX (iDimLon,jDimLat,kDimLev)
                    T      ,&!REAL(KIND=r4), INTENT(IN   ) :: T  (iDimLon,jDimLat,kDimLev)
                    H       )!REAL(KIND=r4), INTENT(OUT  ) :: H  (iDimLon,jDimLat,kDimLev)

   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: DZ (kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: TPG(iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: PTOP
   REAL(KIND=r4), INTENT(IN   ) :: P  (iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: EX (iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: T  (iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: H  (iDimLon,jDimLat,kDimLev)

   INTEGER :: M
   INTEGER :: N
   INTEGER :: L

   REAL(KIND=r4) :: PI 
   REAL(KIND=r4) :: RE 
   REAL(KIND=r4) :: OM 
   REAL(KIND=r4) :: GC 
   REAL(KIND=r4) :: RC 
   REAL(KIND=r4) :: T00
   REAL(KIND=r4) :: P00

   REAL(KIND=r4) :: PR
   REAL(KIND=r4) :: DPR
   REAL(KIND=r4) :: TM
   REAL(KIND=r4) :: DH
   REAL(KIND=r4) :: Z1Z0I
   REAL(KIND=r4) :: RCGC
   REAL(KIND=r4) :: Z1
   REAL(KIND=r4) :: Z0
   REAL(KIND=r4) :: Z
   REAL(KIND=r4) :: HH
!  EX: Exner function (ratio of temperature to potential temperature). 
!  T : Temperature.
!  H : Geopotential height of the sigma-surface.

   CALL CONST (PI   ,&!REAL(KIND=r4), INTENT(OUT  ) :: PI 
               RE   ,&!REAL(KIND=r4), INTENT(OUT  ) :: RE 
               OM   ,&!REAL(KIND=r4), INTENT(OUT  ) :: OM 
               GC   ,&!REAL(KIND=r4), INTENT(OUT  ) :: GC 
               RC   ,&!REAL(KIND=r4), INTENT(OUT  ) :: RC 
               T00  ,&!REAL(KIND=r4), INTENT(OUT  ) :: T00
               P00   )!REAL(KIND=r4), INTENT(OUT  ) :: P00


   RCGC = RC / GC
   CALL SIGMA (kDimLev      ,&!INTEGER      , INTENT(IN   ) :: kDimLev
               DZ      ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
               kDimLev+1    ,&!INTEGER      , INTENT(IN   ) :: N
               Z1       )!REAL(KIND=r4), INTENT(OUT  ) :: Z
   CALL SIGMA (kDimLev      ,&!INTEGER      , INTENT(IN   ) :: kDimLev
               DZ      ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
               0       ,&!INTEGER      , INTENT(IN   ) :: N
               Z0       )!REAL(KIND=r4), INTENT(OUT  ) :: Z
   Z1Z0I = 1.0E0/(Z1-Z0)
   DO M=1,jDimLat
      DO L=1,iDimLon
         HH = TPG(L,M)
         DO N=kDimLev,1,-1
            CALL SIGMA (kDimLev   ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                        DZ   ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                        N    ,&!INTEGER      , INTENT(IN   ) :: N
                        Z     )!REAL(KIND=r4), INTENT(OUT  ) :: Z
!           ZA = Z-0.5E0*DZ(N)
!           ZB = Z+0.5E0*DZ(N)
            PR       = (P(L,M)*(Z-Z0)+PTOP*(Z1-Z))*Z1Z0I
            DPR      = (P(L,M)-PTOP)*DZ(N)*Z1Z0I
            TM       = T(L,M,N) * EX(L,M,N)
            DH       = RCGC * TM * DPR / PR
            H(L,M,N) = HH + 0.5*DH
            HH       = HH + DH
!           PRA = (P(L,M)*(ZA-Z0)+PTOP*(Z1-ZA))*Z1Z0I
!           PRB = (P(L,M)*(ZB-Z0)+PTOP*(Z1-ZB))*Z1Z0I
!           DH = RCGC * TM * (-ALOG(PR/P(L,M))+ALOG(PRB/P(L,M)))
!           H(L,M,N) = HH + DH
!           DH = RCGC * TM * (-ALOG(PRA/P(L,M))+ALOG(PRB/P(L,M)))
!           HH = HH + DH
         END DO
      END DO
   END DO
      RETURN
  END SUBROUTINE DIAGH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DIAGR (iDimLon    ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat    ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                    kDimLev    ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                    DZ    ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ  (kDimLev)
                    PTOP  ,&!REAL(KIND=r4), INTENT(IN   ) :: PTOP
                    P     ,&!REAL(KIND=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
                    T     ,&!REAL(KIND=r4), INTENT(IN   ) :: T(iDimLon,jDimLat,kDimLev)
                    Scale_height    )!REAL(KIND=r4), INTENT(OUT  ) :: Scale_height (iDimLon,jDimLat,kDimLev)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: DZ  (kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: PTOP
   REAL(KIND=r4), INTENT(IN   ) :: P   (iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: T   (iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: Scale_height (iDimLon,jDimLat,kDimLev)
!  P  : Surface pressure.
!  T  : Temperature.
!  Scale_height: Scale height.
   REAL(KIND=r4) :: PI 
   REAL(KIND=r4) :: RE 
   REAL(KIND=r4) :: OM 
   REAL(KIND=r4) :: GC 
   REAL(KIND=r4) :: RC 
   REAL(KIND=r4) :: T00
   REAL(KIND=r4) :: P00
   REAL(KIND=r4) :: Z
   REAL(KIND=r4) :: Z1
   REAL(KIND=r4) :: Z0
   REAL(KIND=r4) :: Z1Z0I
   REAL(KIND=r4) :: PR
   REAL(KIND=r4) :: TM

   INTEGER :: N,M,L
   
      CALL CONST (PI      ,&!REAL(KIND=r4), INTENT(OUT  ) :: PI 
                  RE      ,&!REAL(KIND=r4), INTENT(OUT  ) :: RE 
                  OM      ,&!REAL(KIND=r4), INTENT(OUT  ) :: OM 
                  GC      ,&!REAL(KIND=r4), INTENT(OUT  ) :: GC 
                  RC      ,&!REAL(KIND=r4), INTENT(OUT  ) :: RC 
                  T00     ,&!REAL(KIND=r4), INTENT(OUT  ) :: T00
                  P00      )!REAL(KIND=r4), INTENT(OUT  ) :: P00

      CALL SIGMA (kDimLev      ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DZ      ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                  kDimLev+1    ,&!INTEGER      , INTENT(IN   ) :: N
                  Z1       )!REAL(KIND=r4), INTENT(OUT  ) :: Z

      CALL SIGMA (kDimLev      ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DZ      ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                  0       ,&!INTEGER      , INTENT(IN   ) :: N
                  Z0       )!REAL(KIND=r4), INTENT(OUT  ) :: Z
    Z1Z0I = 1.0E0/(Z1-Z0)
    DO N=1,kDimLev
       DO M=1,jDimLat
          DO L=1,iDimLon

             CALL SIGMA (kDimLev    ,&!INTEGER       , INTENT(IN   ) :: kDimLev
                         DZ    ,&!REAL(KIND=r4) , INTENT(IN   ) :: DZ(kDimLev)
                         N     ,&!INTEGER       , INTENT(IN   ) :: N
                         Z      )!REAL(KIND=r4) , INTENT(OUT  ) :: Z

             PR = (P(L,M)*(Z-Z0)+PTOP*(Z1-Z))*Z1Z0I
             TM = T(L,M,N)*(PR/1.0E5)**(0.2E1/0.7E1)
            Scale_height(L,M,N) = RC*TM/GC
          END DO
       END DO
    END DO
   RETURN
  END SUBROUTINE DIAGR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE BSCPRN (iDimLon      ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                     jDimLat      ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                     kDimLev      ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                     DZ      ,&! REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                     PLEV    ,&! REAL(KIND=r4), INTENT(IN   ) :: PLEV(kDimLev)
                     H       ,&! REAL(KIND=r4), INTENT(IN   ) :: H(iDimLon,jDimLat,kDimLev)
                     U       ,&! REAL(KIND=r4), INTENT(IN   ) :: U(iDimLon,jDimLat,kDimLev)
                     T        )! REAL(KIND=r4), INTENT(IN   ) :: T(iDimLon,jDimLat,kDimLev)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: PLEV(kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: H(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: U(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: T(iDimLon,jDimLat,kDimLev)

   REAL(KIND=r4) :: H0
   REAL(KIND=r4) :: U0
   REAL(KIND=r4) :: T1
   REAL(KIND=r4) :: T2

   INTEGER :: N,L
   WRITE(6,*) 'Basic field data:'
   WRITE(6,'(1X,A5,1X,A5,1X,A5,1X,A6,1X,A6,2(1X,A6))')  &
   'Lyr #','d sgm','Pres.','Z(45N)','U(45N)','th(Eq)','th(NP)'
   DO N=1,kDimLev
      H0 = 0.0E0
      U0 = 0.0E0
      T1 = 0.0E0
      T2 = 0.0E0
      DO L=1,iDimLon
         H0 = H0 + 0.5E0*(H(L,18,N)+H(L,19,N))
         U0 = U0 + 0.5E0*(U(L,18,N)+U(L,19,N))
         T1 = T1 + T(L, 1,N)
         T2 = T2 + T(L,jDimLat,N)
      END DO
      H0 = H0 / REAL(iDimLon)
      U0 = U0 / REAL(iDimLon)
      T1 = T1 / REAL(iDimLon)
      T2 = T2 / REAL(iDimLon)
      WRITE(6,'(1X,I5,1X,F5.3,1X,F5.0,1X,F6.0,1X,F6.1,2(1X,F6.1))') &
        N,DZ(N),PLEV(N),H0,U0,T1,T2
   END DO
   WRITE(6,*) '===================================='
   RETURN
  END SUBROUTINE BSCPRN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ========================================
!   DIAGNOSES OF DEPENDENT VARIABLES
! ========================================

  SUBROUTINE DiagVertSigVelocity (iDimLon      ,&!INTEGER, INTENT(IN   ) :: iDimLon
                    jDimLat      ,&!INTEGER, INTENT(IN   ) :: jDimLat
                    kDimLev      ,&!INTEGER, INTENT(IN   ) :: kDimLev
                    DX      ,&!REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
                    DY      ,&!REAL(KIND=r4), INTENT(IN   ) :: DY  (jDimLat)
                    DZ      ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ  (kDimLev)
                    PTOP    ,&!REAL(KIND=r4), INTENT(IN   ) :: PTOP
                    P      ,& !REAL(KIND=r4), INTENT(IN   ) :: P   (iDimLon,jDimLat)
                    U      ,& !REAL(KIND=r4), INTENT(IN   ) :: U   (iDimLon,jDimLat,kDimLev)
                    V      ,& !REAL(KIND=r4), INTENT(IN   ) :: V   (iDimLon,jDimLat,kDimLev)
                    W      )  !REAL(KIND=r4), INTENT(OUT  ) :: W   (iDimLon,jDimLat,kDimLev)
   INTEGER, INTENT(IN   ) :: iDimLon
   INTEGER, INTENT(IN   ) :: jDimLat
   INTEGER, INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: DY  (jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: DZ  (kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: PTOP
   REAL(KIND=r4), INTENT(IN   ) :: P   (iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: U   (iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: V   (iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: W   (iDimLon,jDimLat,kDimLev)

   REAL(KIND=r4) :: Z
   REAL(KIND=r4) :: Z0
   REAL(KIND=r4) :: Z1
   REAL(KIND=r4) :: Z1Z0I
   REAL(KIND=r4) :: PRI
   REAL(KIND=r4) :: W0
   REAL(KIND=r4) :: W1
   REAL(KIND=r4) :: DX1
   REAL(KIND=r4) :: DX2
   REAL(KIND=r4) :: DXDY
   REAL(KIND=r4) :: DXDYI
   REAL(KIND=r4) :: FX1
   REAL(KIND=r4) :: FX2
   REAL(KIND=r4) :: FY1
   REAL(KIND=r4) :: FY2
   
   INTEGER :: MA
   INTEGER :: MB
   INTEGER :: LA 
   INTEGER :: LB 
   INTEGER :: M1
   INTEGER :: M2   
   INTEGER :: N, M, L,NN
!  P0: Surface pressure of the basic field.
!  U0: Zonal wind of the basic field.
!  V0: Meridional wind of the basic field.
!  P : Surface pressure of the anomaly field.
!  U : Zonal wind of the anomaly field.
!  V : Meridional wind of the anomaly field.
!  W : Diagnosed vertical sigma-velocity of the anomaly field.

   DO N=1,kDimLev
      DO M=1,jDimLat
         DO L=1,iDimLon
            W(L,M,N) = 0.0E0
         END DO
      END DO
   END DO

      CALL SIGMA (kDimLev      ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DZ      ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                  kDimLev+1    ,&!INTEGER      , INTENT(IN   ) :: N
                  Z1       )!REAL(KIND=r4), INTENT(OUT  ) :: Z

      CALL SIGMA (kDimLev    ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DZ    ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                  0     ,&!INTEGER      , INTENT(IN   ) :: N
                  Z0     )!REAL(KIND=r4), INTENT(OUT  ) :: Z

      Z1Z0I = 1.0E0/(Z1-Z0)

   DO M=1,jDimLat
      MA = MOD(M-1+jDimLat-1,jDimLat)+1
      MB = MOD(M+1+jDimLat-1,jDimLat)+1
      M1 = MAX0(M-1, 1)
      M2 = MIN0(M+1,jDimLat)
      DO L=1,iDimLon

         LA = MOD(L-1+iDimLon-1,iDimLon)+1
         LB = MOD(L+1+iDimLon-1,iDimLon)+1
         PRI = 1.0E0/P(L,M)
         W0 = 0.0E0

         DO N=1,kDimLev
            DX1   = 0.5E0*(DX(M)+DX(M1))
            DX2   = 0.5E0*(DX(M)+DX(M2))
            DXDY  = DX(M) * DY(M)
            DXDYI = 1.0E0/DXDY

            FX1   = DY(M) * U(L ,M ,N) * (0.5E0 * (P(L ,M )+P(LA,M ))-PTOP)
            FX2   = DY(M) * U(LB,M ,N) * (0.5E0 * (P(LB,M )+P(L ,M ))-PTOP)

            FY1   = DX1   * V(L ,M ,N) * (0.5E0 * (P(L ,M )+P(L ,MA))-PTOP)
            FY2   = DX2   * V(L ,MB,N) * (0.5E0 * (P(L ,MB)+P(L ,M ))-PTOP)

            W0    = W0 + (FX2 - FX1 + FY2 - FY1 ) * DXDYI * DZ(N) * Z1Z0I
         END DO

         DO N=1,kDimLev-1

            CALL SIGMA (kDimLev      ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                        DZ           ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                        N            ,&!INTEGER      , INTENT(IN   ) :: N
                        Z             )!REAL(KIND=r4), INTENT(OUT  ) :: Z

            Z  = Z + 0.5E0*DZ(N)
            W1 = 0.0E0

            DO NN=1,N

               FX1 = DY(M) * U(L ,M ,NN) * (0.5E0 * (P(L ,M ) + P(LA,M )) - PTOP)
               FX2 = DY(M) * U(LB,M ,NN) * (0.5E0 * (P(LB,M ) + P(L ,M )) - PTOP)

               FY1 = DX1   * V(L ,M ,NN) * (0.5E0 * (P(L ,M ) + P(L ,MA)) - PTOP)
               FY2 = DX2   * V(L ,MB,NN) * (0.5E0 * (P(L ,MB) + P(L ,M )) - PTOP)

               W1 = W1 - (FX2 - FX1 + FY2 - FY1) * DXDYI * DZ(NN) * Z1Z0I

            END DO

            W(L,M,N+1)= W(L,M,N+1) + ( W1 + (Z - Z0) * Z1Z0I * W0 ) * PRI

         END DO
      END DO
   END DO
   RETURN
  END SUBROUTINE DiagVertSigVelocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ========================================

  SUBROUTINE ADVEC (iDimLon       ,&! INTEGER      , INTENT(IN   ) :: iDimLon 
                    jDimLat       ,&! INTEGER      , INTENT(IN   ) :: jDimLat 
                    kDimLev       ,&! INTEGER      , INTENT(IN   ) :: kDimLev 
                    DX       ,&! REAL(KIND=r4), INTENT(IN   ) ::  DX(jDimLat)
                    DY       ,&! REAL(KIND=r4), INTENT(IN   ) ::  DY(jDimLat)
                    DZ       ,&! REAL(KIND=r4), INTENT(IN   ) ::  DZ(kDimLev)
                    U        ,&! REAL(KIND=r4), INTENT(IN   ) ::  U (iDimLon,jDimLat,kDimLev)
                    V        ,&! REAL(KIND=r4), INTENT(IN   ) ::  V (iDimLon,jDimLat,kDimLev)
                    W        ,&! REAL(KIND=r4), INTENT(IN   ) ::  W (iDimLon,jDimLat,kDimLev)
                    T        ,&! REAL(KIND=r4), INTENT(IN   ) ::  T (iDimLon,jDimLat,kDimLev)
                    DU       ,&! REAL(KIND=r4), INTENT(INOUT) ::  DU(iDimLon,jDimLat,kDimLev)
                    DV       ,&! REAL(KIND=r4), INTENT(INOUT) ::  DV(iDimLon,jDimLat,kDimLev)
                    DT        )! REAL(KIND=r4), INTENT(INOUT) ::  DT(iDimLon,jDimLat,kDimLev)
    INTEGER      , INTENT(IN   ) :: iDimLon 
    INTEGER      , INTENT(IN   ) :: jDimLat 
    INTEGER      , INTENT(IN   ) :: kDimLev 
    REAL(KIND=r4), INTENT(IN   ) ::  DX(jDimLat)
    REAL(KIND=r4), INTENT(IN   ) ::  DY(jDimLat)
    REAL(KIND=r4), INTENT(IN   ) ::  DZ(kDimLev)
    REAL(KIND=r4), INTENT(IN   ) ::  U (iDimLon,jDimLat,kDimLev)
    REAL(KIND=r4), INTENT(IN   ) ::  V (iDimLon,jDimLat,kDimLev)
    REAL(KIND=r4), INTENT(IN   ) ::  W (iDimLon,jDimLat,kDimLev)
    REAL(KIND=r4), INTENT(IN   ) ::  T (iDimLon,jDimLat,kDimLev)
    REAL(KIND=r4), INTENT(INOUT) ::  DU(iDimLon,jDimLat,kDimLev)
    REAL(KIND=r4), INTENT(INOUT) ::  DV(iDimLon,jDimLat,kDimLev)
    REAL(KIND=r4), INTENT(INOUT) ::  DT(iDimLon,jDimLat,kDimLev)
    CALL ADVECU (iDimLon       ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                 jDimLat       ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                 kDimLev       ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                 DX       ,&! REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
                 DY       ,&! REAL(KIND=r4), INTENT(IN   ) :: DY(jDimLat)
                 DZ       ,&! REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                 U        ,&! REAL(KIND=r4), INTENT(IN   ) :: U(iDimLon,jDimLat,kDimLev)
                 V        ,&! REAL(KIND=r4), INTENT(IN   ) :: V(iDimLon,jDimLat,kDimLev)
                 W        ,&! REAL(KIND=r4), INTENT(IN   ) :: W(iDimLon,jDimLat,kDimLev)
                 U        ,&! REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                 DU        )! REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)
    CALL ADVECV (iDimLon       ,&!  INTEGER      , INTENT(IN   ) :: iDimLon
                 jDimLat       ,&!  INTEGER      , INTENT(IN   ) :: jDimLat
                 kDimLev       ,&!  INTEGER      , INTENT(IN   ) :: kDimLev
                 DX       ,&!  REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
                 DY       ,&!  REAL(KIND=r4), INTENT(IN   ) :: DY(jDimLat)
                 DZ       ,&!  REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                 U        ,&!  REAL(KIND=r4), INTENT(IN   ) :: U(iDimLon,jDimLat,kDimLev)
                 V        ,&!  REAL(KIND=r4), INTENT(IN   ) :: V(iDimLon,jDimLat,kDimLev)
                 W        ,&!  REAL(KIND=r4), INTENT(IN   ) :: W(iDimLon,jDimLat,kDimLev)
                 V        ,&!  REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                 DV        )!  REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)
    CALL ADVECG (iDimLon       ,& ! INTEGER      , INTENT(IN   ) :: iDimLon 
                 jDimLat       ,& ! INTEGER      , INTENT(IN   ) :: jDimLat
                 kDimLev       ,& ! INTEGER      , INTENT(IN   ) :: kDimLev
                 DX       ,& ! REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
                 DY       ,& ! REAL(KIND=r4), INTENT(IN   ) :: DY(jDimLat)
                 DZ       ,& ! REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                 U        ,& ! REAL(KIND=r4), INTENT(IN   ) :: U(iDimLon,jDimLat,kDimLev)
                 V        ,& ! REAL(KIND=r4), INTENT(IN   ) :: V(iDimLon,jDimLat,kDimLev)
                 W        ,& ! REAL(KIND=r4), INTENT(IN   ) :: W(iDimLon,jDimLat,kDimLev)
                 T        ,& ! REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                 DT        ) ! REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)
    RETURN
  END SUBROUTINE ADVEC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ADVECU (iDimLon       ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                     jDimLat       ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                     kDimLev       ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                     DX       ,&! REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
                     DY       ,&! REAL(KIND=r4), INTENT(IN   ) :: DY(jDimLat)
                     DZ       ,&! REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                     U        ,&! REAL(KIND=r4), INTENT(IN   ) :: U(iDimLon,jDimLat,kDimLev)
                     V        ,&! REAL(KIND=r4), INTENT(IN   ) :: V(iDimLon,jDimLat,kDimLev)
                     W        ,&! REAL(KIND=r4), INTENT(IN   ) :: W(iDimLon,jDimLat,kDimLev)
                     G        ,&! REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                     DG        )! REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: DY(jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: U(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: V(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: W(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)
   INTEGER :: NA
   INTEGER :: NB
   INTEGER :: MA
   INTEGER :: MB
   INTEGER :: LA
   INTEGER :: LB

   REAL(KIND=r4) :: DX1
   REAL(KIND=r4) :: DX2
   REAL(KIND=r4) :: SX1
   REAL(KIND=r4) :: SX2
   REAL(KIND=r4) :: SY1
   REAL(KIND=r4) :: SY2
   REAL(KIND=r4) :: SZ1
   REAL(KIND=r4) :: SZ2
   REAL(KIND=r4) :: DXDYDZ
   REAL(KIND=r4) :: DXDYDZI
   REAL(KIND=r4) :: U1 
   REAL(KIND=r4) :: U2 
   REAL(KIND=r4) :: V1 
   REAL(KIND=r4) :: V2 
   REAL(KIND=r4) :: W1 
   REAL(KIND=r4) :: W2 
   REAL(KIND=r4) :: FX1
   REAL(KIND=r4) :: FX2
   REAL(KIND=r4) :: FY1
   REAL(KIND=r4) :: FY2
   REAL(KIND=r4) :: FZ1
   REAL(KIND=r4) :: FZ2
   
   INTEGER :: N,M,L 
   DO N=1,kDimLev
      NA = MOD(N-1+kDimLev-1,kDimLev)+1
      NB = MOD(N+1+kDimLev-1,kDimLev)+1
      DO M=1,jDimLat
         MA  = MOD(M-1+jDimLat-1,jDimLat)+1
         MB  = MOD(M+1+jDimLat-1,jDimLat)+1
         DX1 = 0.5E0*(DX(M)+DX(MA))
         DX2 = 0.5E0*(DX(M)+DX(MB))
         SX1 = DY(M) * DZ(N)
         SX2 = SX1
         SY1 = DX1 * DZ(N)
         SY2 = DX2 * DZ(N)
         SZ1 = DX(M) * DY(M)
         SZ2 = SZ1
         DXDYDZ = SZ1 * DZ(N)
         DXDYDZI = 1.0E0/DXDYDZ
         DO L=1,iDimLon
            LA = MOD(L-1+iDimLon-1,iDimLon)+1
            LB = MOD(L+1+iDimLon-1,iDimLon)+1
!  -- Technical modification --
!      U1 = 0.5E0*(U(LA,M ,N )+U(L ,M ,N ))
!      U2 = 0.5E0*(U(L ,M ,N )+U(LB,M ,N ))
!      V1 = 0.5E0*(V(LA,M ,N )+V(L ,M ,N ))
!      V2 = 0.5E0*(V(LA,MB,N )+V(L ,MB,N ))
!      W1 = 0.5E0*(W(LA,M ,N )+W(L ,M ,N ))
!      W2 = 0.5E0*(W(LA,M ,NB)+W(L ,M ,NB))
!  --
            U1 = (U(LA,M ,N )+U(L ,M ,N ))
            U2 = (U(L ,M ,N )+U(LB,M ,N ))
            V1 = (V(LA,M ,N )+V(L ,M ,N ))
            V2 = (V(LA,MB,N )+V(L ,MB,N ))
            W1 = (W(LA,M ,N )+W(L ,M ,N ))
            W2 = (W(LA,M ,NB)+W(L ,M ,NB))
!  --
            FX1 = SX1 * U1 * (G(L ,M ,N )-G(LA,M ,N ))
            FX2 = SX2 * U2 * (G(LB,M ,N )-G(L ,M ,N ))
            FY1 = SY1 * V1 * (G(L ,M ,N )-G(L ,MA,N ))
            FY2 = SY2 * V2 * (G(L ,MB,N )-G(L ,M ,N ))
            FZ1 = SZ1 * W1 * (G(L ,M ,N )-G(L ,M ,NA))
            FZ2 = SZ2 * W2 * (G(L ,M ,NB)-G(L ,M ,N ))
            DG(L,M,N) = DG(L,M,N)                        &
!  -- Technical modification --
!     +    - 0.5E0 * (FX2+FX1+FY2+FY1+FZ2+FZ1) * DXDYDZI
!  --
          - (FX2+FX1+FY2+FY1+FZ2+FZ1) * DXDYDZI
!  --
         END DO
      END DO
   END DO
   RETURN
  END SUBROUTINE ADVECU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE ADVECV (iDimLon      ,&!  INTEGER      , INTENT(IN   ) :: iDimLon
                     jDimLat      ,&!  INTEGER      , INTENT(IN   ) :: jDimLat
                     kDimLev      ,&!  INTEGER      , INTENT(IN   ) :: kDimLev
                     DX      ,&!  REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
                     DY      ,&!  REAL(KIND=r4), INTENT(IN   ) :: DY(jDimLat)
                     DZ      ,&!  REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                     U       ,&!  REAL(KIND=r4), INTENT(IN   ) :: U(iDimLon,jDimLat,kDimLev)
                     V       ,&!  REAL(KIND=r4), INTENT(IN   ) :: V(iDimLon,jDimLat,kDimLev)
                     W       ,&!  REAL(KIND=r4), INTENT(IN   ) :: W(iDimLon,jDimLat,kDimLev)
                     G       ,&!  REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                     DG       )!  REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: DY(jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: U(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: V(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: W(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)


   REAL(KIND=r4) :: DX0
   REAL(KIND=r4) :: DY0
   REAL(KIND=r4) :: SX1
   REAL(KIND=r4) :: SX2
   REAL(KIND=r4) :: SY1
   REAL(KIND=r4) :: SY2
   REAL(KIND=r4) :: SZ1
   REAL(KIND=r4) :: SZ2
   REAL(KIND=r4) :: DXDYDZ 
   REAL(KIND=r4) :: DXDYDZI
   REAL(KIND=r4) :: U1 
   REAL(KIND=r4) :: U2 
   REAL(KIND=r4) :: V1 
   REAL(KIND=r4) :: V2 
   REAL(KIND=r4) :: W1 
   REAL(KIND=r4) :: W2 
   REAL(KIND=r4) :: FX1
   REAL(KIND=r4) :: FX2
   REAL(KIND=r4) :: FY1
   REAL(KIND=r4) :: FY2
   REAL(KIND=r4) :: FZ1
   REAL(KIND=r4) :: FZ2

   INTEGER :: kb
   INTEGER :: kc
   INTEGER :: kf
   INTEGER :: yb
   INTEGER :: yc
   INTEGER :: yf
   INTEGER :: xb
   INTEGER :: xc
   INTEGER :: xf
   INTEGER :: k,j,i


   DO k=1,kDimLev
      kb = MOD(k-1+kDimLev-1,kDimLev)+1
      kc = k
      kf = MOD(k+1+kDimLev-1,kDimLev)+1
      DO j=2,jDimLat
         yb = MOD(j-1+jDimLat-1,jDimLat)+1
         yc = j
         yf = MOD(j+1+jDimLat-1,jDimLat)+1

         DX0 = 0.5E0*(DX(yc)+DX(yb))
         DY0 = 0.5E0*(DY(yc)+DY(yb))
         !
         !
         !
         ! SX1  = 0.5E0*(DY(yc)+DY(yb))* DZ(kc)
         !
         SX1 = ( 0.5E0*(DY(yc)+DY(yb))) * DZ(kc)
         SX2 = ( 0.5E0*(DY(yc)+DY(yb))) * DZ(kc)

         SY1 = DX(yb) * DZ(kc)
         SY2 = DX(yc) * DZ(kc)
         !
         !SZ1 =  0.5E0*(DX(yc)+DX(yb)) * 0.5E0*(DY(yc)+DY(yb))

         !SZ1 =  0.25E0*((DX(yc)+DX(yb))*(DY(yc)+DY(yb)))
         !
         SZ1 = DX0 * ( 0.5E0*(DY(yc)+DY(yb)))
         SZ2 = SZ1
         !
         !DXDYDZ =  0.25E0*((DX(yc)+DX(yb))*(DY(yc)+DY(yb)))* DZ(kc)
         !
         DXDYDZ  = SZ1 * DZ(kc)

         DXDYDZI = 1.0E0/DXDYDZ

         DO i=1,iDimLon
            xb = MOD(i-1+iDimLon-1,iDimLon)+1
            xc = i
            xf = MOD(i+1+iDimLon-1,iDimLon)+1
!            PRINT*,'iDimLon=',iDimLon ,'xc=',xc ,'xb=',xb ,'xf=',xf 
!  -- Technical modification --
!      U1 = 0.5E0*(U(xc ,yb,kc )+U(xc ,yc ,kc ))
!      U2 = 0.5E0*(U(xf,yb,kc )+U(xf,yc ,kc ))
!      V1 = 0.5E0*(V(xc ,yb,kc )+V(xc ,yc ,kc ))
!      V2 = 0.5E0*(V(xc ,yc ,kc )+V(xc ,yf,kc ))
!      W1 = 0.5E0*(W(xc ,yb,kc )+W(xc ,yc ,kc ))
!      W2 = 0.5E0*(W(xc ,yb,kf)+W(xc ,yc ,kf))
!  --
            !                                     
            ! d(v)         d(v)         d(v)         d(v)   
            ! -----  + u * ------ + v * ------- w * ------- =0
            ! dt           dx            dy           dz    
            !

            U1 = (U(xc ,yb ,kc )+U(xc ,yc ,kc ))
            U2 = (U(xf ,yb ,kc )+U(xf ,yc ,kc ))

            V1 = (V(xc ,yb ,kc )+V(xc ,yc ,kc ))
            V2 = (V(xc ,yc ,kc )+V(xc ,yf ,kc ))

            W1 = (W(xc ,yb ,kc )+W(xc ,yc ,kc ))
            W2 = (W(xc ,yb ,kf )+W(xc ,yc ,kf ))
!  --
           !          1
           !  SX1=   --- * (DY(yc) + DY(yb)) * DZ(kc)
           !          2
           !
            FX1 = SX1 * U1 * (G(xc ,yc ,kc )-G(xb ,yc ,kc ))
           !
           !          1
           !  SX2=   --- * (DY(yc) + DY(yb)) * DZ(kc)
           !          2
           !
            FX2 = SX2 * U2 * (G(xf ,yc ,kc )-G(xc ,yc ,kc ))

            FY1 = SY1 * V1 * (G(xc ,yc ,kc )-G(xc ,yb ,kc ))
            FY2 = SY2 * V2 * (G(xc ,yf ,kc )-G(xc ,yc ,kc ))

            FZ1 = SZ1 * W1 * (G(xc ,yc ,kc )-G(xc ,yc ,kb))
            FZ2 = SZ2 * W2 * (G(xc ,yc ,kf )-G(xc ,yc ,kc ))
            !
            !          1
            !  SX1=   --- * (DY(yc) + DY(yb)) * DZ(kc)
            !          2
            !
            !              1
            !  DXDYDZ =  -----*((DX(yc)+DX(yb))*(DY(yc)+DY(yb)))* DZ(kc)
            !              4
            !
            !  SX1            2 * (DY(yc) + DY(yb)) * DZ(kc)                     2
            !------  =  ------------------------------------------- =      ---------------
            !DXDYDZI      ((DX(yc)+DX(yb))*(DY(yc)+DY(yb)))* DZ(kc)         (DX(yc)+DX(yb))
            !

            !
            !          1
            !  SX2=   --- * (DY(yc) + DY(yb)) * DZ(kc)
            !          2
            !
            !              1
            !  DXDYDZ =  -----*((DX(yc)+DX(yb))*(DY(yc)+DY(yb)))* DZ(kc)
            !              4
            !
            !  SX2            2 * (DY(yc) + DY(yb)) * DZ(kc)                     2
            !------  =  ------------------------------------------- =      ---------------
            !DXDYDZI      ((DX(yc)+DX(yb))*(DY(yc)+DY(yb)))* DZ(kc)         (DX(yc)+DX(yb))
            !
            DG(xc,yc,kc) = DG(xc,yc,kc) - (FX2+FX1+FY2+FY1+FZ2+FZ1) * (1.0E0/DXDYDZ)

!  -- Technical modification --
!         DG(xc,yc,kc) = DG(xc,yc,kc) - 0.5E0 * (FX2+FX1+FY2+FY1+FZ2+FZ1) * DXDYDZI
!  --
!  --
         END DO
      END DO
   END DO

   RETURN
  END SUBROUTINE ADVECV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ADVECG (iDimLon      ,&! INTEGER      , INTENT(IN   ) :: iDimLon 
                     jDimLat      ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                     kDimLev      ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                     DX      ,&! REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
                     DY      ,&! REAL(KIND=r4), INTENT(IN   ) :: DY(jDimLat)
                     DZ      ,&! REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                     U       ,&! REAL(KIND=r4), INTENT(IN   ) :: U(iDimLon,jDimLat,kDimLev)
                     V       ,&! REAL(KIND=r4), INTENT(IN   ) :: V(iDimLon,jDimLat,kDimLev)
                     W       ,&! REAL(KIND=r4), INTENT(IN   ) :: W(iDimLon,jDimLat,kDimLev)
                     G       ,&! REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                     DG       )! REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)
    INTEGER      , INTENT(IN   ) :: iDimLon 
    INTEGER      , INTENT(IN   ) :: jDimLat
    INTEGER      , INTENT(IN   ) :: kDimLev
    REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
    REAL(KIND=r4), INTENT(IN   ) :: DY(jDimLat)
    REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
    REAL(KIND=r4), INTENT(IN   ) :: U(iDimLon,jDimLat,kDimLev)
    REAL(KIND=r4), INTENT(IN   ) :: V(iDimLon,jDimLat,kDimLev)
    REAL(KIND=r4), INTENT(IN   ) :: W(iDimLon,jDimLat,kDimLev)
    REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
    REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)
    REAL(KIND=r4) :: FX1
    REAL(KIND=r4) :: FX2
    REAL(KIND=r4) :: FY1
    REAL(KIND=r4) :: FY2
    REAL(KIND=r4) :: FZ1
    REAL(KIND=r4) :: FZ2
    INTEGER :: kb
    INTEGER :: kc
    INTEGER :: kf
    INTEGER :: yb
    INTEGER :: yc
    INTEGER :: yf
    INTEGER :: xb
    INTEGER :: xc
    INTEGER :: xf

    INTEGER :: k,j,i

    DO k=1,kDimLev
       kb = MOD(k-1+kDimLev-1,kDimLev)+1
       kc = k
       kf = MOD(k+1+kDimLev-1,kDimLev)+1
       DO j=1,jDimLat
          yb = MOD(j-1+jDimLat-1,jDimLat)+1
          yc = j
          yf = MOD(j+1+jDimLat-1,jDimLat)+1
          DO i=1,iDimLon
             xb = MOD(i-1+iDimLon-1,iDimLon)+1
             xc = i
             xf = MOD(i+1+iDimLon-1,iDimLon)+1


             FX1 = (DY(yc) * DZ(kc))                * U(xc ,yc ,kc ) * (G(xc ,yc ,kc )-G(xb ,yc ,kc ))
             FX2 = (DY(yc) * DZ(kc))                * U(xf ,yc ,kc ) * (G(xf ,yc ,kc )-G(xc ,yc ,kc ))

             FY1 = (0.5E0*(DX(yc)+DX(yb))) * DZ(kc) * V(xc ,yc ,kc ) * (G(xc ,yc ,kc )-G(xc ,yb ,kc ))
             FY2 = (0.5E0*(DX(yc)+DX(yf))) * DZ(kc) * V(xc ,yf ,kc ) * (G(xc ,yf ,kc )-G(xc ,yc ,kc ))

             FZ1 = (DX(yc) * DY(yc))                * W(xc ,yc ,kc ) * (G(xc ,yc ,kc )-G(xc ,yc ,kb ))
             FZ2 = (DX(yc) * DY(yc))                * W(xc ,yc ,kf ) * (G(xc ,yc ,kf )-G(xc ,yc ,kc ))

             DG(xc,yc,kc) = DG(xc,yc,kc) - 0.5E0 * (FX2+FX1+FY2+FY1+FZ2+FZ1) * (1.0/ ( DX(yc) * DY(yc)) * DZ(kc))
          END DO
       END DO
    END DO
    RETURN
  END  SUBROUTINE ADVECG 

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE RSGN (iDimLon        ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat        ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                   kDimLev        ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                   G          )!REAL(KIND=r4), INTENT(INOUT) ::   G(iDimLon,jDimLat,kDimLev)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)

   INTEGER :: i,j,k

   DO k=1,kDimLev
      DO j=1,jDimLat
         DO i=1,iDimLon
            G(i,j,k) = - G(i,j,k)
         END DO
      END DO
   END DO

   RETURN
  END SUBROUTINE RSGN
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ZERO (iDimLon     ,&!  INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat     ,&!  INTEGER      , INTENT(IN   ) :: jDimLat
                   kDimLev     ,&!  INTEGER      , INTENT(IN   ) :: kDimLev
                   G       )!  REAL(KIND=r4), INTENT(OUT  ) :: G(iDimLon,jDimLat,kDimLev)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(OUT  ) :: G(iDimLon,jDimLat,kDimLev)
   INTEGER :: i,j,k
   DO k=1,kDimLev
      DO j=1,jDimLat
         DO i=1,iDimLon
            G(i,j,k) = 0.0E0
         END DO
      END DO
   END DO

   RETURN
  END SUBROUTINE ZERO 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE INITC (iDimLon    ,&! INTEGER     , INTENT(IN    ) :: iDimLon
                    jDimLat    ,&! INTEGER     , INTENT(IN    ) :: jDimLat
                    kDimLev    ,&! INTEGER     , INTENT(IN    ) :: kDimLev
                    PLEV  ,&! REAL(KIND=r4), INTENT(IN   ) :: PLEV(kDimLev)
                    TPG   ,&! REAL(KIND=r4), INTENT(IN   ) :: TPG(iDimLon,jDimLat)
                    U     ,&! REAL(KIND=r4), INTENT(OUT  ) :: U (iDimLon,jDimLat,kDimLev)
                    V     ,&! REAL(KIND=r4), INTENT(OUT  ) :: V(iDimLon,jDimLat,kDimLev)
                    T     ,&! REAL(KIND=r4), INTENT(OUT  ) :: T(iDimLon,jDimLat,kDimLev)
                    P     ,&! REAL(KIND=r4), INTENT(OUT  ) :: P(iDimLon,jDimLat)
                    UP    ,&! REAL(KIND=r4), INTENT(OUT  ) :: UP(iDimLon,jDimLat,kDimLev)
                    VP    ,&! REAL(KIND=r4), INTENT(OUT  ) :: VP(iDimLon,jDimLat,kDimLev)
                    TP    ,&! REAL(KIND=r4), INTENT(OUT  ) :: TP(iDimLon,jDimLat,kDimLev)
                    PP     )! REAL(KIND=r4), INTENT(OUT  ) :: PP(iDimLon,jDimLat)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: PLEV(kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: TPG(iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(OUT  ) :: U (iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: V(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: T(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: P(iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(OUT  ) :: UP(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: VP(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: TP(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: PP(iDimLon,jDimLat)

   REAL(KIND=r4) :: PI 
   REAL(KIND=r4) :: RE 
   REAL(KIND=r4) :: OM 
   REAL(KIND=r4) :: GC 
   REAL(KIND=r4) :: RC 
   REAL(KIND=r4) :: T00
   REAL(KIND=r4) :: P00

   REAL(KIND=r4) :: RTG1
   REAL(KIND=r4) :: PLEV1
   
   INTEGER k,j,i
   
   CALL CONST (PI   ,&!REAL(KIND=r4), INTENT(OUT  ) :: PI 
               RE   ,&!REAL(KIND=r4), INTENT(OUT  ) :: RE 
               OM   ,&!REAL(KIND=r4), INTENT(OUT  ) :: OM 
               GC   ,&!REAL(KIND=r4), INTENT(OUT  ) :: GC 
               RC   ,&!REAL(KIND=r4), INTENT(OUT  ) :: RC 
               T00  ,&!REAL(KIND=r4), INTENT(OUT  ) :: T00
               P00   )!REAL(KIND=r4), INTENT(OUT  ) :: P00

   RTG1 = RC*T00/GC

   DO k=1,kDimLev
      DO j=1,jDimLat
         DO i=1,iDimLon
            P(i,j)    = P00*EXP(-TPG(i,j)/RTG1)
            PLEV1     = (1.0E2*PLEV(k)) * (P(i,j)/P00)
            T(i,j,k)  = T00/((PLEV1/1.0E5)**(0.2E1/0.7E1))
            U(i,j,k)  = 0.0E0
            V(i,j,k)  = 0.0E0
            PP(i,j)   = P(i,j)
            TP(i,j,k) = T(i,j,k)
            UP(i,j,k) = U(i,j,k)
            VP(i,j,k) = V(i,j,k)
         END DO
      END DO
   END DO
   RETURN
  END SUBROUTINE INITC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE TPGEV (iDimLon    ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat    ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                    kDimLev    ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                    Scale_height   ,&! REAL(KIND=r4), INTENT(IN   ) :: Scale_height(iDimLon,jDimLat,kDimLev)
                    TPG   ,&! REAL(KIND=r4), INTENT(IN   ) :: TPG(iDimLon,jDimLat)
                    TIM   ,&! REAL(KIND=r4), INTENT(IN   ) :: TIM
                    TPG1  ,&! REAL(KIND=r4), INTENT(INOUT) :: TPG1(iDimLon,jDimLat)
                    P     ,&! REAL(KIND=r4), INTENT(INOUT) :: P(iDimLon,jDimLat)
                    PP     )! REAL(KIND=r4), INTENT(INOUT) :: PP(iDimLon,jDimLat)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: Scale_height(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: TPG(iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: TIM
   REAL(KIND=r4), INTENT(INOUT) :: TPG1(iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(INOUT) :: P(iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(INOUT) :: PP(iDimLon,jDimLat)

   REAL(KIND=r4)  :: C
   REAL(KIND=r4)  :: RATIO
   REAL(KIND=r4)  :: TPG1P
   INTEGER        :: M,L

   C = AMIN1(1.0E0,TIM/0.864000E6)
   DO M=1,jDimLat
      DO L=1,iDimLon
         TPG1P = TPG1(L,M)
         TPG1(L,M) = C * TPG(L,M)
         RATIO = 1.0E0 - (TPG1(L,M)-TPG1P)/Scale_height(L,M,kDimLev)
         P (L,M) = P (L,M) * RATIO
         PP(L,M) = PP(L,M) * RATIO
      END DO
   END DO
   RETURN
  END SUBROUTINE TPGEV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ========================================
!   DYNAMICAL PROCESSES
! ========================================

  SUBROUTINE Convergence (iDimLon     ,&!INTEGER     , INTENT(IN ) ::  iDimLon
                    jDimLat     ,&!INTEGER     , INTENT(IN ) ::  jDimLat
                    kDimLev     ,&!INTEGER     , INTENT(IN ) ::  kDimLev
                    DX     ,&!REAL(KIND=r4), INTENT(IN  ) ::  DX(jDimLat)
                    DY     ,&!REAL(KIND=r4), INTENT(IN  ) ::  DY(jDimLat)
                    DZ     ,&!REAL(KIND=r4), INTENT(IN  ) ::  DZ(kDimLev)
                    PTOP   ,&!REAL(KIND=r4), INTENT(IN  ) ::  PTOP
                    P      ,&!REAL(KIND=r4), INTENT(IN  ) ::  P(iDimLon,jDimLat)
                    U      ,&!REAL(KIND=r4), INTENT(IN  ) ::  U(iDimLon,jDimLat,kDimLev)
                    V      ,&!REAL(KIND=r4), INTENT(IN  ) ::  V(iDimLon,jDimLat,kDimLev)
                    DP      )!REAL(KIND=r4), INTENT(INOUT) ::  DP(iDimLon,jDimLat)
   INTEGER      , INTENT(IN   ) ::  iDimLon
   INTEGER      , INTENT(IN   ) ::  jDimLat
   INTEGER      , INTENT(IN   ) ::  kDimLev
   REAL(KIND=r4), INTENT(IN   ) ::  DX(jDimLat)
   REAL(KIND=r4), INTENT(IN   ) ::  DY(jDimLat)
   REAL(KIND=r4), INTENT(IN   ) ::  DZ(kDimLev)
   REAL(KIND=r4), INTENT(IN   ) ::  PTOP
   REAL(KIND=r4), INTENT(IN   ) ::  P(iDimLon,jDimLat)        !  P : Surface pressure of the anomaly field.
   REAL(KIND=r4), INTENT(IN   ) ::  U(iDimLon,jDimLat,kDimLev)!  U : Zonal wind of the anomaly field.
   REAL(KIND=r4), INTENT(IN   ) ::  V(iDimLon,jDimLat,kDimLev)!  V : Meridional wind of the anomaly field.
   REAL(KIND=r4), INTENT(INOUT) ::  DP(iDimLon,jDimLat)       !  DP: Time-derivative of surface pressure.

   REAL(KIND=r4) :: Z1
   REAL(KIND=r4) :: Z0
   REAL(KIND=r4) :: Z1Z0I
   REAL(KIND=r4) :: DX1_bar 
   REAL(KIND=r4) :: DX2_bar 
   REAL(KIND=r4) :: FX1
   REAL(KIND=r4) :: FX2
   REAL(KIND=r4) :: FY1
   REAL(KIND=r4) :: FY2   
   INTEGER        :: i,j,k
   INTEGER        :: xb,xc,xf   
   INTEGER        :: yb,yc,yf   
!  P : Surface pressure of the anomaly field.
!  U : Zonal wind of the anomaly field.
!  V : Meridional wind of the anomaly field.
!  DP: Time-derivative of surface pressure.
      CALL SIGMA (kDimLev          ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DZ          ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                  kDimLev+1        ,&!INTEGER      , INTENT(IN   ) :: N
                  Z1           )!REAL(KIND=r4), INTENT(OUT  ) :: Z
      CALL SIGMA (kDimLev          ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DZ          ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                  0           ,&!INTEGER      , INTENT(IN   ) :: N
                  Z0           )!REAL(KIND=r4), INTENT(OUT  ) :: Z

      Z1Z0I = 1.0E0/(Z1-Z0)

      DO j=1,jDimLat
         yb = MOD(j-1+jDimLat-1,jDimLat)+1
         yc = j
         yf = MOD(j+1+jDimLat-1,jDimLat)+1

         DX1_bar = 0.5E0*(DX(yc)+DX(yb))
         DX2_bar = 0.5E0*(DX(yc)+DX(yf))


         DO i=1,iDimLon

            xb = MOD(i-1+iDimLon-1,iDimLon)+1
            xc = i
            xf = MOD(i+1+iDimLon-1,iDimLon)+1

            DO k=1,kDimLev
               !
               !                 _
               !  FX1  = DY * U* Ps
               !
               !
               FX1 = DY(yc)     * U(xc ,yc ,k) * (0.5E0*(P(xc ,yc ) + P(xb ,yc )) - PTOP)
               FX2 = DY(yc)     * U(xf ,yc ,k) * (0.5E0*(P(xf ,yc ) + P(xc ,yc )) - PTOP)

               FY1 = DX1_bar    * V(xc ,yc ,k) * (0.5E0*(P(xc ,yc ) + P(xc ,yb )) - PTOP)
               FY2 = DX2_bar    * V(xc ,yf ,k) * (0.5E0*(P(xc ,yf ) + P(xc ,yc )) - PTOP)
               !
               !             1
               !        -------
               !        \    --                 --    --        --
               !         \  |  d ps*u      d ps*v |  |            |
               ! DP = -   \ | ------- +   ------- |* | sig2 - sig1|
               !          / | dx          dy      |  |            |
               !         /  |                     |   --        --
               !        /    --                 --
               !        ------
               !           sig_top
               
               DP(xc,yc) = DP(xc,yc) - (FX2 - FX1 + FY2 - FY1 ) * ( 1.0E0/(DX(yc)*DY(yc))) * DZ(k) * Z1Z0I

            END DO

         END DO

      END DO

   RETURN
  END SUBROUTINE Convergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ========================================
!   TEMPORAL EVOLUTION
!========================================

  SUBROUTINE TENDZ (iDimLon    ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat    ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                    kDimLev    ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                    DU    ,&!REAL(KIND=r4), INTENT(OUT  ) :: DU(iDimLon,jDimLat,kDimLev)
                    DV    ,&!REAL(KIND=r4), INTENT(OUT  ) :: DV(iDimLon,jDimLat,kDimLev)
                    DT    ,&!REAL(KIND=r4), INTENT(OUT  ) :: DT(iDimLon,jDimLat,kDimLev)
                    DP     )!REAL(KIND=r4), INTENT(OUT  ) :: DP(iDimLon,jDimLat)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(OUT  ) :: DU(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: DV(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: DT(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: DP(iDimLon,jDimLat)
   INTEGER        :: j,k,i

   DO k=1,kDimLev
      DO j=1,jDimLat
         DO i=1,iDimLon
            DU(i,j,k) = 0.0E0
            DV(i,j,k) = 0.0E0
            DT(i,j,k) = 0.0E0
         END DO
      END DO
   END DO
   DO j=1,jDimLat
      DO i=1,iDimLon
         DP(i,j) = 0.0E0
      END DO
   END DO
      RETURN
  END SUBROUTINE TENDZ 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ========================================

  SUBROUTINE PGRAD (iDimLon       ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat       ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                    kDimLev       ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                    DX       ,&!REAL(KINd=r4), INTENT(IN   ) :: DX(jDimLat)
                    DY       ,&!REAL(KINd=r4), INTENT(IN   ) :: DY(jDimLat)
                    P        ,&!REAL(KINd=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
                    H        ,&!REAL(KINd=r4), INTENT(IN   ) :: H(iDimLon,jDimLat,kDimLev)
                    Scale_height      ,&!REAL(KINd=r4), INTENT(IN   ) :: Scale_height(iDimLon,jDimLat,kDimLev)
                    DU       ,&!REAL(KINd=r4), INTENT(INOUT) :: DU(iDimLon,jDimLat,kDimLev)
                    DV        )!REAL(KINd=r4), INTENT(INOUT) :: DV(iDimLon,jDimLat,kDimLev)

   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KINd=r4), INTENT(IN   ) :: DX(jDimLat)
   REAL(KINd=r4), INTENT(IN   ) :: DY(jDimLat)
   REAL(KINd=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
   REAL(KINd=r4), INTENT(IN   ) :: H(iDimLon,jDimLat,kDimLev)
   REAL(KINd=r4), INTENT(IN   ) :: Scale_height(iDimLon,jDimLat,kDimLev)
   REAL(KINd=r4), INTENT(INOUT) :: DU(iDimLon,jDimLat,kDimLev)
   REAL(KINd=r4), INTENT(INOUT) :: DV(iDimLon,jDimLat,kDimLev)

   REAL(KIND=r4) :: PI 
   REAL(KIND=r4) :: RE 
   REAL(KIND=r4) :: OM 
   REAL(KIND=r4) :: GC 
   REAL(KIND=r4) :: RC 
   REAL(KIND=r4) :: T00
   REAL(KIND=r4) :: P00
   INTEGER       :: N
   INTEGER       :: j
   INTEGER       :: i
   INTEGER       :: xb
   INTEGER       :: xc
   INTEGER       :: yb
   INTEGER       :: yc
   REAL(KIND=r4) :: PP
   CALL CONST (PI   ,& !REAL(KIND=r4), INTENT(OUT  ) :: PI 
               RE   ,& !REAL(KIND=r4), INTENT(OUT  ) :: RE 
               OM   ,& !REAL(KIND=r4), INTENT(OUT  ) :: OM 
               GC   ,& !REAL(KIND=r4), INTENT(OUT  ) :: GC 
               RC   ,& !REAL(KIND=r4), INTENT(OUT  ) :: RC 
               T00  ,& !REAL(KIND=r4), INTENT(OUT  ) :: T00
               P00  )  !REAL(KIND=r4), INTENT(OUT  ) :: P00

!      CALL SIGMA (kDimLev,DZ,kDimLev+1,Z1)
!      CALL SIGMA (kDimLev,DZ,   0,Z0)
!      Z1Z0I = 1.0E0/(Z1-Z0)
   DO N=1,kDimLev
      !      CALL SIGMA (kDimLev,DZ,N,Z)
      !      ZR = (Z-Z0)*Z1Z0I
      !      C1 = (Z1-Z)*Z1Z0I
      !      C2 = (Z-Z0)*Z1Z0I
      DO j=1,jDimLat
         yc = j
         DO i=1,iDimLon
            xc = i
            xb = MOD(i-1+iDimLon-1,iDimLon)+1
            PP = 0.5E0*(P(xb,yc)+P(xc,yc))
            !
            !             d z        R * T       d Ps
            ! DU = - g * ------- - -------- * ---------
            !             dx          Ps          dx
            !
            !
            DU(xc,yc,N) = DU(xc,yc,N) - GC * (  (H(xc,yc,N)-H(xb,yc,N))                                &
                                           + Scale_height(xc,yc,N)*(P(xc,yc)-P(xb,yc))/PP) * (1.0E0/DX(yc))
                                                 
            !        PR  = C1*PTOP+C2*P(xc ,yc)
            !        PRA = C1*PTOP+C2*P(xb,yc)
            !        DU(xc,yc,N) = DU(xc,yc,N)
            !     +    - GC * (  (H(xc,yc,N)-H(xb,yc,N))
            !     +            + Scale_height(xc,yc,N)*ALOG(PR/PRA))
            !     +         * DXI
         END DO
      END DO
   END DO

   DO N=1,kDimLev
      !      CALL SIGMA (kDimLev,DZ,N,Z)
      !      ZR = (Z-Z0)*Z1Z0I
      !      C1 = (Z1-Z)*Z1Z0I
      !      C2 = (Z-Z0)*Z1Z0I
      DO j=2,jDimLat
         yb = MOD(j-1+jDimLat-1,jDimLat)+1
         yc = j
         DO i=1,iDimLon
            xc = i
            PP = 0.5E0*(P(xc,yb)+P(xc,yc))
            !
            !             d z        R * T       d Ps
            ! DV = - g * ------- - -------- * ---------
            !             dy          Ps          dy
            !
            !
            DV(xc,yc,N) = DV(xc,yc,N) - GC * (  (H(xc,yc,N)-H(xc,yb,N))       &
                                      + Scale_height(xc,yc,N)*(P(xc,yc)-P(xc,yb))/PP) * (1.0E0/(0.5E0*(DY(yb)+DY(yc))))
            !        PR  = C1*PTOP+C2*P(xc,yc )
            !        PRA = C1*PTOP+C2*P(xc,yb)
            !        DV(xc,yc,N) = DV(xc,yc,N)
            !     +    - GC * (  (H(xc,yc,N)-H(xc,yb,N))
            !     +            + Scale_height(xc,yc,N)*ALOG(PR/PRA))
            !     +         * DYI
         END DO
      END DO
   END DO
   RETURN
  END SUBROUTINE PGRAD

  SUBROUTINE CORIO (iDimLon     ,&! INTEGER         , INTENT(IN   ) :: iDimLon
                    jDimLat     ,&! INTEGER         , INTENT(IN   ) :: jDimLat
                    kDimLev     ,&! INTEGER         , INTENT(IN   ) :: kDimLev
                    CG     ,&! CHARACTER(LEN=1), INTENT(IN   ) :: CG
                    FC     ,&! REAL(KIND=r4)   , INTENT(IN   ) :: FC(jDimLat)
                    U      ,&! REAL(KIND=r4)   , INTENT(IN   ) :: U(iDimLon,jDimLat,kDimLev)
                    V      ,&! REAL(KIND=r4)   , INTENT(IN   ) :: V(iDimLon,jDimLat,kDimLev)
                    DU     ,&! REAL(KIND=r4)   , INTENT(INOUT) :: DU(iDimLon,jDimLat,kDimLev)
                    DV      )! REAL(KIND=r4)   , INTENT(INOUT) :: DV(iDimLon,jDimLat,kDimLev)
   INTEGER         , INTENT(IN   ) :: iDimLon
   INTEGER         , INTENT(IN   ) :: jDimLat
   INTEGER         , INTENT(IN   ) :: kDimLev
   CHARACTER(LEN=1), INTENT(IN   ) :: CG
   REAL(KIND=r4)   , INTENT(IN   ) :: FC(jDimLat)
   REAL(KIND=r4)   , INTENT(IN   ) :: U(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4)   , INTENT(IN   ) :: V(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4)   , INTENT(INOUT) :: DU(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4)   , INTENT(INOUT) :: DV(iDimLon,jDimLat,kDimLev)

   REAL(KIND=r4) :: PI 
   REAL(KIND=r4) :: RE 
   REAL(KIND=r4) :: OM 
   REAL(KIND=r4) :: GC 
   REAL(KIND=r4) :: RC 
   REAL(KIND=r4) :: T00
   REAL(KIND=r4) :: P00
   REAL(KIND=r4) :: REI
   REAL(KIND=r4) :: TH
   REAL(KIND=r4) :: TANTH
   REAL(KIND=r4) :: TANRE
   REAL(KIND=r4) :: FF
   REAL(KIND=r4) :: FC1
   REAL(KIND=r4) :: VV
   REAL(KIND=r4) :: UU
   INTEGER       :: yb
   INTEGER       :: yc
   INTEGER       :: yf
   INTEGER       :: xb
   INTEGER       :: xc
   INTEGER       :: xf
   INTEGER       :: j
   INTEGER       :: k
   INTEGER       :: i

   CALL CONST (PI    ,& !REAL(KIND=r4), INTENT(OUT  ) :: PI 
               RE    ,& !REAL(KIND=r4), INTENT(OUT  ) :: RE 
               OM    ,& !REAL(KIND=r4), INTENT(OUT  ) :: OM 
               GC    ,& !REAL(KIND=r4), INTENT(OUT  ) :: GC 
               RC    ,& !REAL(KIND=r4), INTENT(OUT  ) :: RC 
               T00   ,& !REAL(KIND=r4), INTENT(OUT  ) :: T00
               P00    ) !REAL(KIND=r4), INTENT(OUT  ) :: P00

   REI = 1.0E0/RE

   DO j=1,jDimLat
      yc  = j
      yf = MOD(j+1+jDimLat-1,jDimLat)+1
      IF (CG.EQ.'H') TH = 0.5E0*PI*(REAL(yc)-0.5E0)/REAL(jDimLat)
      IF (CG.EQ.'G') TH = PI*((REAL(yc)-0.5E0)/REAL(jDimLat)-0.5E0)
      TANTH = TAN(TH)
      TANRE = TANTH * REI
      DO k=1,kDimLev
         DO i=1,iDimLon
            xb = MOD(i-1+iDimLon-1,iDimLon)+1
            xc = i

            FF = FC(yc) + U(xc,yc,k) * TANRE

            VV = 0.25E0*(V(xb,yc,k)+V(xc,yc,k)+V(xb,yf,k)+V(xc,yf,k))
            !
            !              u*v
            ! DU =  f*v + ------ * tan(fhi)
            !              ra
            !
            DU(xc,yc,k) = DU(xc,yc,k) + FF * VV

         END DO
      END DO
   END DO
 
   DO j=2,jDimLat
      yb = MOD(j-1+jDimLat-1,jDimLat)+1
      yc  = j

      IF (CG.EQ.'H') TH = 0.5E0*PI*REAL(yc-1)/REAL(jDimLat)
      IF (CG.EQ.'G') TH = PI*(REAL(yc-1)/REAL(jDimLat)-0.5E0)
      TANTH = TAN(TH)
      TANRE = TANTH * REI
      FC1 = 0.5E0*(FC(yb)+FC(yc))
      DO k=1,kDimLev
         DO i=1,iDimLon

            xc = i
            xf = MOD(xc+1+iDimLon-1,iDimLon)+1
            !
            !       U(xc,yb)+U(xf,yb)+U(xc,yc)+U(xf,yc)
            ! UU = -------------------------------------
            !                       4
            !
            UU = 0.25E0 * ( U(xc,yb,k)+U(xf,yb,k)+U(xc,yc,k)+U(xf,yc,k) )
            !
            !                u
            ! FF = + f   + ------ * tan(fhi)
            !                ra
            !
            FF = FC1 + UU * TANRE
            !
            !                u*u
            ! DV = - f*u - ------ * tan(fhi)
            !                ra
            !
            DV(xc,yc,k) = DV(xc,yc,k) - FF * UU

         END DO
      END DO
   END DO
   RETURN
  END SUBROUTINE CORIO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE INTEL (iDimLon        ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat        ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                    kDimLev        ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                    DeltaT       ,&!REAL(KIND=r4), INTENT(IN   ) :: DeltaT
                    DU        ,&!REAL(KIND=r4), INTENT(IN   ) :: DU(iDimLon,jDimLat,kDimLev)
                    DV        ,&!REAL(KIND=r4), INTENT(IN   ) :: DV(iDimLon,jDimLat,kDimLev)
                    DT        ,&!REAL(KIND=r4), INTENT(IN   ) :: DT(iDimLon,jDimLat,kDimLev)
                    DP        ,&!REAL(KIND=r4), INTENT(IN   ) :: DP(iDimLon,jDimLat)
                    UP        ,&!REAL(KIND=r4), INTENT(INOUT) :: UP(iDimLon,jDimLat,kDimLev)
                    VP        ,&!REAL(KIND=r4), INTENT(INOUT) :: VP(iDimLon,jDimLat,kDimLev)
                    TP        ,&!REAL(KIND=r4), INTENT(INOUT) :: TP(iDimLon,jDimLat,kDimLev)
                    PP        ,&!REAL(KIND=r4), INTENT(INOUT) :: PP(iDimLon,jDimLat)
                    U         ,&!REAL(KIND=r4), INTENT(INOUT) :: U (iDimLon,jDimLat,kDimLev)
                    V         ,&!REAL(KIND=r4), INTENT(INOUT) :: V (iDimLon,jDimLat,kDimLev)
                    T         ,&!REAL(KIND=r4), INTENT(INOUT) :: T (iDimLon,jDimLat,kDimLev)
                    P          )!REAL(KIND=r4), INTENT(INOUT) :: P (iDimLon,jDimLat),
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: DeltaT
   REAL(KIND=r4), INTENT(IN   ) :: DU(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: DV(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: DT(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: DP(iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(INOUT) :: UP(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: VP(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: TP(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: PP(iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(INOUT) :: U (iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: V (iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: T (iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: P (iDimLon,jDimLat)

      CALL INTEL1 (iDimLon      ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat      ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                   kDimLev      ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                   DeltaT     ,&! REAL(KIND=r4), INTENT(IN   ) :: DeltaT
                   DU      ,&! REAL(KIND=r4), INTENT(IN   ) :: DG(iDimLon,jDimLat,kDimLev)
                   UP      ,&! REAL(KIND=r4), INTENT(INOUT) :: GP(iDimLon,jDimLat,kDimLev)
                   U        )! REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
      CALL INTEL1 (iDimLon      ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat      ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                   kDimLev      ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                   DeltaT     ,&! REAL(KIND=r4), INTENT(IN   ) :: DeltaT
                   DV      ,&! REAL(KIND=r4), INTENT(IN   ) :: DG(iDimLon,jDimLat,kDimLev)
                   VP      ,&! REAL(KIND=r4), INTENT(INOUT) :: GP(iDimLon,jDimLat,kDimLev)
                   V        )! REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
      CALL INTEL1 (iDimLon      ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat      ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                   kDimLev      ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                   DeltaT     ,&! REAL(KIND=r4), INTENT(IN   ) :: DeltaT
                   DT      ,&! REAL(KIND=r4), INTENT(IN   ) :: DG(iDimLon,jDimLat,kDimLev)
                   TP      ,&! REAL(KIND=r4), INTENT(INOUT) :: GP(iDimLon,jDimLat,kDimLev)
                   T        )! REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
      CALL INTEL1 (iDimLon      ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat      ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                   1       ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                   DeltaT     ,&! REAL(KIND=r4), INTENT(IN   ) :: DeltaT
                   DP      ,&! REAL(KIND=r4), INTENT(IN   ) :: DG(iDimLon,jDimLat,kDimLev)
                   PP      ,&! REAL(KIND=r4), INTENT(INOUT) :: GP(iDimLon,jDimLat,kDimLev)
                   P        )! REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
      RETURN
  END SUBROUTINE INTEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE INTEL1 (iDimLon     ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                     jDimLat     ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                     kDimLev     ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                     DeltaT    ,&! REAL(KIND=r4), INTENT(IN   ) :: DeltaT
                     DG     ,&! REAL(KIND=r4), INTENT(IN   ) :: DG(iDimLon,jDimLat,kDimLev)
                     GP     ,&! REAL(KIND=r4), INTENT(INOUT) :: GP(iDimLon,jDimLat,kDimLev)
                     G       )! REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: DeltaT
   REAL(KIND=r4), INTENT(IN   ) :: DG(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: GP(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4) ::  DDeltaT 
   REAL(KIND=r4) ::  EF   
   REAL(KIND=r4) ::  EF1  
   REAL(KIND=r4) ::  EF2  
   REAL(KIND=r4) :: GP1
   INTEGER       :: N,M,L
   DDeltaT = 0.2E1*DeltaT
   EF   = 0.5E-1
   EF1  = (1.0E0-0.2E1*EF)/(1.0E0-EF)
   EF2  = EF/(1.0E0-EF)
   DO N=1,kDimLev
      DO M=1,jDimLat
         DO L=1,iDimLon
            GP1       = GP(L,M,N)
            GP(L,M,N) = EF1*G(L,M,N) + EF2*GP(L,M,N)
            G (L,M,N) = GP1 + DDeltaT * DG(L,M,N)
         END DO
      END DO
   END DO
   RETURN
  END SUBROUTINE INTEL1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ========================================
!   PHYSICAL PROCESSES
! ========================================

  SUBROUTINE Diffusion (iDimLon    ,&! INTEGER     , INTENT(IN   ) :: iDimLon
                    jDimLat    ,&! INTEGER     , INTENT(IN   ) :: jDimLat
                    kDimLev    ,&! INTEGER     , INTENT(IN   ) :: kDimLev
                    DX    ,&! REAL(KIND=r4), INTENT(IN  ) :: DX(jDimLat)
                    DY    ,&! REAL(KIND=r4), INTENT(IN  ) :: DY(jDimLat)
                    DZ    ,&! REAL(KIND=r4), INTENT(IN  ) :: DZ(kDimLev)
                    Scale_height   ,&! REAL(KIND=r4), INTENT(IN  ) :: Scale_height(iDimLon,jDimLat,kDimLev
                    U     ,&! REAL(KIND=r4), INTENT(IN  ) :: U (iDimLon,jDimLat,kDimLev)
                    V     ,&! REAL(KIND=r4), INTENT(IN  ) :: V (iDimLon,jDimLat,kDimLev)
                    T     ,&! REAL(KIND=r4), INTENT(IN  ) :: T (iDimLon,jDimLat,kDimLev)
                    DU    ,&! REAL(KIND=r4), INTENT(INOUT) :: DU(iDimLon,jDimLat,kDimLev)
                    DV    ,&! REAL(KIND=r4), INTENT(INOUT) :: DV(iDimLon,jDimLat,kDimLev)
                    DT     )! REAL(KIND=r4), INTENT(INOUT) :: DT(iDimLon,jDimLat,kDimLev)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: DY(jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: Scale_height(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: U (iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: V (iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: T (iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: DU(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: DV(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: DT(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4) :: DFH
   REAL(KIND=r4) :: DFV

   DFH = 1.0E6
   DFV = 1.0E0
   CALL DiffusionU (iDimLon      ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                jDimLat      ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                kDimLev      ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                DX      ,&! REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
                DY      ,&! REAL(KIND=r4), INTENT(IN   ) :: DY(jDimLat)
                DZ      ,&! REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                Scale_height     ,&! REAL(KIND=r4), INTENT(IN   ) :: Scale_height(iDimLon,jDimLat,kDimLev)
                DFH     ,&! REAL(KIND=r4), INTENT(IN   ) :: DFH
                DFV     ,&! REAL(KIND=r4), INTENT(IN   ) :: DFV
                U       ,&! REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                DU       )! REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)

   CALL DiffusionV (iDimLon      ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                jDimLat      ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                kDimLev      ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                DX      ,&!REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
                DY      ,&!REAL(KIND=r4), INTENT(IN   ) :: DY(jDimLat)
                DZ      ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                Scale_height     ,&!REAL(KIND=r4), INTENT(IN   ) :: Scale_height(iDimLon,jDimLat,kDimLev)
                DFH     ,&!REAL(KIND=r4), INTENT(IN   ) :: DFH
                DFV     ,&!REAL(KIND=r4), INTENT(IN   ) :: DFV
                V       ,&!REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                DV       )!REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)

   CALL DiffusionG (iDimLon      ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                jDimLat      ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                kDimLev      ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                DX      ,&! REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
                DY      ,&! REAL(KIND=r4), INTENT(IN   ) :: DY(jDimLat)
                DZ      ,&! REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                Scale_height     ,&! REAL(KIND=r4), INTENT(IN   ) :: Scale_height(iDimLon,jDimLat,kDimLev)
                DFH     ,&! REAL(KIND=r4), INTENT(IN   ) :: DFH
                DFV     ,&! REAL(KIND=r4), INTENT(IN   ) :: DFV
                T       ,&! REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                DT       )! REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)
   RETURN
  END SUBROUTINE Diffusion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE DiffusionU (iDimLon       ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                     jDimLat       ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                     kDimLev       ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                     DX       ,&! REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
                     DY       ,&! REAL(KIND=r4), INTENT(IN   ) :: DY(jDimLat)
                     DZ       ,&! REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                     Scale_height      ,&! REAL(KIND=r4), INTENT(IN   ) :: Scale_height(iDimLon,jDimLat,kDimLev)
                     DFH      ,&! REAL(KIND=r4), INTENT(IN   ) :: DFH
                     DFV      ,&! REAL(KIND=r4), INTENT(IN   ) :: DFV
                     G        ,&! REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                     DG        )! REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)

   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: DY(jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: Scale_height(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: DFH
   REAL(KIND=r4), INTENT(IN   ) :: DFV
   REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)

   REAL(KIND=r4) :: PI 
   REAL(KIND=r4) :: RE 
   REAL(KIND=r4) :: OM 
   REAL(KIND=r4) :: GC 
   REAL(KIND=r4) :: RC 
   REAL(KIND=r4) :: T00
   REAL(KIND=r4) :: P00

   REAL(KIND=r4) :: DZ1 
   REAL(KIND=r4) :: DZ2 
   REAL(KIND=r4) :: DZ1I
   REAL(KIND=r4) :: DZ2I
   REAL(KIND=r4) :: DFV0
   REAL(KIND=r4) :: DFV1
   REAL(KIND=r4) :: DFV2
   REAL(KIND=r4) :: DZI
   
   REAL(KIND=r4) :: Z
   REAL(KIND=r4) :: Z1
   REAL(KIND=r4) :: Z2
   INTEGER :: i,j,k
   INTEGER :: yb
   INTEGER :: yc
   INTEGER :: yf
   INTEGER :: xb
   INTEGER :: xc
   INTEGER :: xf
   INTEGER :: kb
   INTEGER :: kc
   INTEGER :: kf

   CALL CONST (PI   ,& !REAL(KIND=r4), INTENT(OUT  ) :: PI 
               RE   ,& !REAL(KIND=r4), INTENT(OUT  ) :: RE 
               OM   ,& !REAL(KIND=r4), INTENT(OUT  ) :: OM 
               GC   ,& !REAL(KIND=r4), INTENT(OUT  ) :: GC 
               RC   ,& !REAL(KIND=r4), INTENT(OUT  ) :: RC 
               T00  ,& !REAL(KIND=r4), INTENT(OUT  ) :: T00
               P00   ) !REAL(KIND=r4), INTENT(OUT  ) :: P00

   DO k=1,kDimLev
      kb=k-1
      kc=k
      kf=k+1
      DO j=1,jDimLat
         yb   = MOD(j-1+jDimLat-1,jDimLat)+1
         yc   = j
         yf   = MOD(j+1+jDimLat-1,jDimLat)+1
         yb   = MAX(j-1, 1)
         yc   = j
         yf   = MIN(j+1,jDimLat)

         DO i=1,iDimLon
            xb = MOD(i-1+iDimLon-1,iDimLon)+1
            xc = i
            xf = MOD(i+1+iDimLon-1,iDimLon)+1
            !
            !             -                         -
            !            |  d d u             d d u  |
            ! DG = nih * |------------ + ------------|
            !            |  dx dx            dy dy   |
            !             -                         -
            !
            DG(xc,yc,kc) = DG(xc,yc,kc)                                                  &
                       + DFH * (G(xb,yc,kc)+G(xf,yc,kc)-G(xc,yc,kc)-G(xc,yc,kc)) * (1.0E0/(DX(yc)**2))  &
                       + DFH * (G(xc,yb,kc)+G(xc,yf,kc)-G(xc,yc,kc)-G(xc,yc,kc)) * (1.0E0/(DY(yc)**2))

         END DO

      END DO
   END DO

   DO k=1,kDimLev
      kb = MOD(k-1+kDimLev-1,kDimLev)+1
      kc = k
      kf = MOD(k+1+kDimLev-1,kDimLev)+1
      kb = MAX(k-1, 1)
      kc = k
      kf = MIN(k+1,kDimLev)

      CALL SIGMA (kDimLev       ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DZ       ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                  kc        ,&!INTEGER      , INTENT(IN   ) :: kc
                  Z         )!REAL(KIND=r4), INTENT(OUT  ) :: Z

      CALL SIGMA (kDimLev       ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DZ       ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                  kb       ,&!INTEGER      , INTENT(IN   ) :: N
                  Z1        )!REAL(KIND=r4), INTENT(OUT  ) :: Z

      CALL SIGMA (kDimLev       ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DZ       ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                  kf       ,&!INTEGER      , INTENT(IN   ) :: N
                  Z2        )!REAL(KIND=r4), INTENT(OUT  ) :: Z

      DZ1  = 0.5E0*(DZ(kb)+DZ(kc))
      DZ2  = 0.5E0*(DZ(kc)+DZ(kf))
      DZ1I = 1.0E0/DZ1
      DZ2I = 1.0E0/DZ2

      DO j=1,jDimLat
         yc   = j
         DO i=1,iDimLon
            xc = i
            !
            !                        R * T
            !Scale_height         = ----------
            !                          g
            !
            !
            !                                g
            !DFV0                 =DFV * ----------
            !                              R * T
            !
            DFV0 = DFV / Scale_height(xc,yc,kc)
            !
            !                              g * sig1
            !DFV1                 =DFV * -------------
            !                               R * T
            !
            DFV1 = DFV0 * (0.5E0*(Z+Z1))
            !
            !
            !                              g * sig2
            !DFV2                 =DFV * -------------
            !                               R * T
            !
            DFV2 = DFV0 * (0.5E0*(Z+Z2))
            !
            !            g         1
            !DZI =  ---------- * ------
            !          R * T       dz
            !
            DZI  = 1.0E0 / (DZ(kc)*Scale_height(xc,yc,kc))
            !
            !                                -                                  -
            !        g * sig       1        |        g*sig             d u       |
            ! DG = ----------- * -------- * | niv * -------- * ----------------- |
            !        R * T         dsig     |         R * T            dsig      |
            !                                -                                  -
            !
  
            DG(xc,yc,kc) = DG(xc,yc,kc)    +  DZI*(  DFV2*(G(xc,yc,kf)-G(xc,yc,kc))*DZ2I       &
                                           -         DFV1*(G(xc,yc,kb)-G(xc,yc,kc))*DZ1I)  * Z
            IF (kc.EQ.kDimLev) THEN
               DG(xc,yc,kc) = DG(xc,yc,kc) -  DZI*(  DFV2*G(xc,yc,kc)*(0.2E1*DZ2I)) * (0.5E0*(Z+Z2))
            ENDIF

         END DO
      END DO
   END DO

   RETURN
  END SUBROUTINE DiffusionU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DiffusionV (iDimLon      ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                     jDimLat      ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                     kDimLev      ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                     DX      ,&!REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
                     DY      ,&!REAL(KIND=r4), INTENT(IN   ) :: DY(jDimLat)
                     DZ      ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                     Scale_height     ,&!REAL(KIND=r4), INTENT(IN   ) :: Scale_height(iDimLon,jDimLat,kDimLev)
                     DFH     ,&!REAL(KIND=r4), INTENT(IN   ) :: DFH
                     DFV     ,&!REAL(KIND=r4), INTENT(IN   ) :: DFV
                     G       ,&!REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                     DG       )!REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: DY(jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: Scale_height(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: DFH
   REAL(KIND=r4), INTENT(IN   ) :: DFV
   REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)

   REAL(KIND=r4) :: PI 
   REAL(KIND=r4) :: RE 
   REAL(KIND=r4) :: OM 
   REAL(KIND=r4) :: GC 
   REAL(KIND=r4) :: RC 
   REAL(KIND=r4) :: T00
   REAL(KIND=r4) :: P00
   REAL(KIND=r4) :: DX2I
   REAL(KIND=r4) :: DY2I

   REAL(KIND=r4) :: DZ1 
   REAL(KIND=r4) :: DZ2 
   REAL(KIND=r4) :: DZ1I
   REAL(KIND=r4) :: DZ2I
   REAL(KIND=r4) :: DX0 
   REAL(KIND=r4) :: DY0 
   REAL(KIND=r4) :: DFV0
   REAL(KIND=r4) :: DFV1
   REAL(KIND=r4) :: DFV2
   REAL(KIND=r4) :: DZI

   REAL(KIND=r4) :: Z
   REAL(KIND=r4) :: Z1
   REAL(KIND=r4) :: Z2
   INTEGER :: N,M,L
   INTEGER :: MA
   INTEGER :: MB
   INTEGER :: NA
   INTEGER :: NB
   INTEGER :: N1
   INTEGER :: N2
   INTEGER :: LA
   INTEGER :: LB

   CALL CONST (PI   ,& !REAL(KIND=r4), INTENT(OUT  ) :: PI 
               RE   ,& !REAL(KIND=r4), INTENT(OUT  ) :: RE 
               OM   ,& !REAL(KIND=r4), INTENT(OUT  ) :: OM 
               GC   ,& !REAL(KIND=r4), INTENT(OUT  ) :: GC 
               RC   ,& !REAL(KIND=r4), INTENT(OUT  ) :: RC 
               T00  ,& !REAL(KIND=r4), INTENT(OUT  ) :: T00
               P00   ) !REAL(KIND=r4), INTENT(OUT  ) :: P00

   DO N=1,kDimLev
      DO M=2,jDimLat
         MA   = MOD(M-1+jDimLat-1,jDimLat)+1
         MB   = MOD(M+1+jDimLat-1,jDimLat)+1
         DX0  = 0.5E0*(DX(MA)+DX(M))
         DY0  = 0.5E0*(DY(MA)+DY(M))
         DX2I = 1.0E0/(DX0**2)
         DY2I = 1.0E0/(DY0**2)
         DO L=1,iDimLon
            LA = MOD(L-1+iDimLon-1,iDimLon)+1
            LB = MOD(L+1+iDimLon-1,iDimLon)+1
            DG(L,M,N) = DG(L,M,N)                     &
                      + DFH * (G(LA,M,N)+G(LB,M,N)-G(L,M,N)-G(L,M,N)) * DX2I      &
                      + DFH * (G(L,MA,N)+G(L,MB,N)-G(L,M,N)-G(L,M,N)) * DY2I
         END DO
      END DO
   END DO

   DO N=1,kDimLev
      NA = MOD(N-1+kDimLev-1,kDimLev)+1
      NB = MOD(N+1+kDimLev-1,kDimLev)+1
      N1 = MAX(N-1, 1)
      N2 = MIN(N+1,kDimLev)

      CALL SIGMA (kDimLev       ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DZ       ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                  N        ,&!INTEGER      , INTENT(IN   ) :: N
                  Z         )!REAL(KIND=r4), INTENT(OUT  ) :: Z

      CALL SIGMA (kDimLev       ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DZ       ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                  N1       ,&!INTEGER      , INTENT(IN   ) :: N
                  Z1        )!REAL(KIND=r4), INTENT(OUT  ) :: Z

      CALL SIGMA (kDimLev       ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DZ       ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                  N2       ,&!INTEGER      , INTENT(IN   ) :: N
                  Z2        )!REAL(KIND=r4), INTENT(OUT  ) :: Z

      DZ1  = 0.5E0*(DZ(N1)+DZ(N ))
      DZ2  = 0.5E0*(DZ(N )+DZ(N2))
      DZ1I = 1.0E0/DZ1
      DZ2I = 1.0E0/DZ2
      DO M=2,jDimLat
         DO L=1,iDimLon
            DFV0 = DFV / Scale_height(L,M,N)
            DFV1 = DFV0 * (0.5E0*(Z+Z1))
            DFV2 = DFV0 * (0.5E0*(Z+Z2))
            DZI  = 1.0E0 / (DZ(N)*Scale_height(L,M,N))
            DG(L,M,N) = DG(L,M,N)                              &
                      + (  DFV2*(G(L,M,N2)-G(L,M,N))*DZ2I      &
                      -    DFV1*(G(L,M,N1)-G(L,M,N))*DZ1I)        &
                      * DZI * Z
            IF (N.EQ.kDimLev) THEN
               DG(L,M,N) = DG(L,M,N) - DFV2*G(L,M,N)*(0.2E1*DZ2I) * DZI * (0.5E0*(Z+Z2))
            ENDIF
         END DO
      END DO
   END DO
   RETURN
  END SUBROUTINE DiffusionV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DiffusionG (iDimLon     ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                     jDimLat     ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                     kDimLev     ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                     DX     ,&! REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
                     DY     ,&! REAL(KIND=r4), INTENT(IN   ) :: DY(jDimLat)
                     DZ     ,&! REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                     Scale_height    ,&! REAL(KIND=r4), INTENT(IN   ) :: Scale_height(iDimLon,jDimLat,kDimLev)
                     DFH    ,&! REAL(KIND=r4), INTENT(IN   ) :: DFH
                     DFV    ,&! REAL(KIND=r4), INTENT(IN   ) :: DFV
                     G      ,&! REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                     DG      )! REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)

   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: DY(jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: Scale_height(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: DFH
   REAL(KIND=r4), INTENT(IN   ) :: DFV
   REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)

   REAL(KIND=r4) :: PI 
   REAL(KIND=r4) :: RE 
   REAL(KIND=r4) :: OM 
   REAL(KIND=r4) :: GC 
   REAL(KIND=r4) :: RC 
   REAL(KIND=r4) :: T00
   REAL(KIND=r4) :: P00
   REAL(KIND=r4) :: DX2I
   REAL(KIND=r4) :: DY2I
   REAL(KIND=r4) :: Z
   REAL(KIND=r4) :: Z1
   REAL(KIND=r4) :: Z2
   REAL(KIND=r4) :: DZ1 
   REAL(KIND=r4) :: DZ2 
   REAL(KIND=r4) :: DZ1I
   REAL(KIND=r4) :: DZ2I

   REAL(KIND=r4) :: DFV0
   REAL(KIND=r4) :: DFV1
   REAL(KIND=r4) :: DFV2
   INTEGER :: xb,xf,yb,yf
   INTEGER :: k,j,i
!   INTEGER :: MA
!   INTEGER :: MB
!   INTEGER :: M1
!   INTEGER :: M2
   INTEGER :: NA
   INTEGER :: NB
   INTEGER :: N1
   INTEGER :: N2
!   INTEGER :: LA
!   INTEGER :: LB
   
   CALL CONST (PI    ,& !REAL(KIND=r4), INTENT(OUT  ) :: PI 
               RE    ,& !REAL(KIND=r4), INTENT(OUT  ) :: RE 
               OM    ,& !REAL(KIND=r4), INTENT(OUT  ) :: OM 
               GC    ,& !REAL(KIND=r4), INTENT(OUT  ) :: GC 
               RC    ,& !REAL(KIND=r4), INTENT(OUT  ) :: RC 
               T00   ,& !REAL(KIND=r4), INTENT(OUT  ) :: T00
               P00    ) !REAL(KIND=r4), INTENT(OUT  ) :: P00

   DO k=1,kDimLev
      DO j=1,jDimLat
         !yb   = MOD(j-1+jDimLat-1,jDimLat)+1  ! j=1==> jDimLat-1
         !yf   = MOD(j+1+jDimLat-1,jDimLat)+1  ! j=1==> jDimLat+1
         yb   = MAX(j-1, 1)
         yf   = MIN(j+1,jDimLat)
         DX2I = 1.0E0/(DX(j)**2)
         DY2I = 1.0E0/(DY(j)**2)
         DO i=1,iDimLon
            xb = MOD(i-1+iDimLon-1,iDimLon)+1
            xf = MOD(i+1+iDimLon-1,iDimLon)+1
            !
            !             -                         -
            !            |  d d theta      d dthata  |
            ! DG = nih * |------------ + ------------|
            !            |  dx dx            dy dy   |
            !             -                         -
            DG(i,j,k) = DG(i,j,k)                                           &
                      + DFH * (G(xf,j,k) - 2.0*G(i,j,k) + G(xb,j,k)) * DX2I &
                      + DFH * (G(i,yf,k) - 2.0*G(i,j,k) + G(i,yb,k)) * DY2I
         END DO
      END DO
   END DO
   
   DO k=1,kDimLev
      NA = MOD(k-1+kDimLev-1,kDimLev)+1
      NB = MOD(k+1+kDimLev-1,kDimLev)+1
      N1 = MAX(k-1, 1)
      N2 = MIN(k+1,kDimLev)

      CALL SIGMA (kDimLev    ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DZ    ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                  k     ,&!INTEGER      , INTENT(IN   ) :: N
                  Z      )!REAL(KIND=r4), INTENT(OUT  ) :: Z

      CALL SIGMA (kDimLev    ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DZ    ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                  N1    ,&!INTEGER      , INTENT(IN   ) :: N
                  Z1     )!REAL(KIND=r4), INTENT(OUT  ) :: Z

      CALL SIGMA (kDimLev    ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DZ    ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                  N2    ,&!INTEGER      , INTENT(IN   ) :: N
                  Z2     )!REAL(KIND=r4), INTENT(OUT  ) :: Z

      DZ1  = 0.5E0*(DZ(N1)+DZ(k ))
      DZ2  = 0.5E0*(DZ(k )+DZ(N2))
      DZ1I = 1.0E0/DZ1
      DZ2I = 1.0E0/DZ2

      DO j=1,jDimLat
         DO i=1,iDimLon
            !
            !                        R * T
	    !Scale_height         = ----------
            !                          g
            !
            !
            !                                g
	    !DFV0                 =DFV * ----------
            !                              R * T
            !

            DFV0 = DFV / Scale_height(i,j,k)
            !
            !                              g * sig1
	    !DFV1                 =DFV * -------------
            !                               R * T
            !

            DFV1 = DFV0 * (0.5E0*(Z+Z1))
            !
            !                              g * sig2
	    !DFV1                 =DFV * -------------
            !                               R * T
            !
            DFV2 = DFV0 * (0.5E0*(Z+Z2))
            !
	    !                                -                                  -
	    !        g * sig       1        |        g*sig       Theta2 -Theta1  |
	    ! DG = ----------- * -------- * | niv * -------- * ----------------- |
	    !        R * T         Dsig     |         R * T            Dsig      |
	    !                                -                                  -
	    !
            DG(i,j,k) = DG(i,j,k) + (  DFV2 * (G(i,j,N2)-G(i,j,k)) * DZ2I         &
                                     - DFV1 * (G(i,j,N1)-G(i,j,k)) * DZ1I )       &
                                      / (DZ(k)*Scale_height(i,j,k)) * Z
         END DO
      END DO
   END DO
   RETURN
  END SUBROUTINE DiffusionG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ========================================

  SUBROUTINE ADBTC (iDimLon     ,&!INTEGER     , INTENT(IN   ) :: iDimLon
                    jDimLat     ,&!INTEGER     , INTENT(IN   ) :: jDimLat
                    kDimLev     ,&!INTEGER     , INTENT(IN   ) :: kDimLev
                    UF     ,&!REAL(KIND=r4), INTENT(IN  ) :: UF(iDimLon,jDimLat,kDimLev)
                    VF     ,&!REAL(KIND=r4), INTENT(IN  ) :: VF(iDimLon,jDimLat,kDimLev)
                    TF     ,&!REAL(KIND=r4), INTENT(IN  ) :: TF(iDimLon,jDimLat,kDimLev)
                    DU     ,&!REAL(KIND=r4), INTENT(INOUT) :: DU(iDimLon,jDimLat,kDimLev)
                    DV     ,&!REAL(KIND=r4), INTENT(INOUT) :: DV(iDimLon,jDimLat,kDimLev)
                    DT      )!REAL(KIND=r4), INTENT(INOUT) :: DT(iDimLon,jDimLat,kDimLev)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: UF(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: VF(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: TF(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: DU(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: DV(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: DT(iDimLon,jDimLat,kDimLev)
   INTEGER   :: i,j,k
   DO k=1,kDimLev
      DO j=1,jDimLat
         DO i=1,iDimLon
!           DU(i,j,k) = DU(i,j,k) + UF(i,j,k)
!           DV(i,j,k) = DV(i,j,k) + VF(i,j,k)
            DT(i,j,k) = DT(i,j,k) + TF(i,j,k)
         END DO
      END DO
   END DO
   RETURN
  END SUBROUTINE ADBTC
! ========================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE Damping (iDimLon      ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat      ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                    kDimLev      ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                    U0      ,&!REAL(KIND=r4), INTENT(IN   ) :: U0(iDimLon,jDimLat,kDimLev)
                    V0      ,&!REAL(KIND=r4), INTENT(IN   ) :: V0(iDimLon,jDimLat,kDimLev)
                    T0      ,&!REAL(KIND=r4), INTENT(IN   ) :: T0(iDimLon,jDimLat,kDimLev)
                    U       ,&!REAL(KIND=r4), INTENT(IN   ) :: U (iDimLon,jDimLat,kDimLev)
                    V       ,&!REAL(KIND=r4), INTENT(IN   ) :: V (iDimLon,jDimLat,kDimLev)
                    T       ,&!REAL(KIND=r4), INTENT(IN   ) :: T (iDimLon,jDimLat,kDimLev)
                    DU      ,&!REAL(KIND=r4), INTENT(INOUT) :: DU(iDimLon,jDimLat,kDimLev)
                    DV      ,&!REAL(KIND=r4), INTENT(INOUT) :: DV(iDimLon,jDimLat,kDimLev)
                    DT       )!REAL(KIND=r4), INTENT(INOUT) :: DT(iDimLon,jDimLat,kDimLev)

   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: U0(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: V0(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: T0(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: U (iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: V (iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: T (iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: DU(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: DV(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: DT(iDimLon,jDimLat,kDimLev)

   REAL(KIND=r4) :: TAU
   REAL(KIND=r4) :: TAU1
   REAL(KIND=r4) :: TAU2
   REAL(KIND=r4) :: TAUI
   REAL(KIND=r4) :: TAU1I
   REAL(KIND=r4) :: TAU2I

   INTEGER :: N,M,L
   
   TAU  = 0.86400E5*2.0E1
   TAU2 = 0.86400E5*0.5E1
   TAU1 = 0.86400E5*0.2E1
   !      TAU  = 0.86400E5*1.0E1
   !      TAU2 = 0.86400E5*0.3E1
   !      TAU1 = 0.86400E5*0.1E1
   TAUI  = 1.0E0/TAU
   TAU2I = 1.0E0/TAU2
   TAU1I = 1.0E0/TAU1

   DO N=1,kDimLev
      DO M=1,jDimLat
         DO L=1,iDimLon
            IF (N.EQ.kDimLev) THEN
               DU(L,M,N) = DU(L,M,N) - TAU1I * (U(L,M,N)-U0(L,M,N))
               DV(L,M,N) = DV(L,M,N) - TAU1I * (V(L,M,N)-V0(L,M,N))
               DT(L,M,N) = DT(L,M,N) - TAU1I * (T(L,M,N)-T0(L,M,N))
            ELSE
               IF (N.EQ.kDimLev-1) THEN
                  DU(L,M,N) = DU(L,M,N) - TAU2I * (U(L,M,N)-U0(L,M,N))
                  DV(L,M,N) = DV(L,M,N) - TAU2I * (V(L,M,N)-V0(L,M,N))
                  DT(L,M,N) = DT(L,M,N) - TAU2I * (T(L,M,N)-T0(L,M,N))
               ELSE
                  DU(L,M,N) = DU(L,M,N) - TAUI * (U(L,M,N)-U0(L,M,N))
                  DV(L,M,N) = DV(L,M,N) - TAUI * (V(L,M,N)-V0(L,M,N))
                  DT(L,M,N) = DT(L,M,N) - TAUI * (T(L,M,N)-T0(L,M,N))
               ENDIF
            ENDIF
         END DO
      END DO
   END DO
   RETURN
  END SUBROUTINE Damping

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ========================================
!   FORCINGS
! ========================================

  SUBROUTINE FORCE (iDimLon      ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat      ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                    kDimLev      ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                    DZ      ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                    B0      ,&!REAL(KIND=r4), INTENT(IN   ) :: B0
                    A0      ,&!REAL(KIND=r4), INTENT(IN   ) :: A0
                    R0      ,&!REAL(KIND=r4), INTENT(IN   ) :: R0
                    Z1      ,&!REAL(KIND=r4), INTENT(IN   ) :: Z1
                    Z2      ,&!REAL(KIND=r4), INTENT(IN   ) :: Z2
                    C       ,&!REAL(KIND=r4), INTENT(IN   ) :: C 
                    DG       )!REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)

   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: B0
   REAL(KIND=r4), INTENT(IN   ) :: A0
   REAL(KIND=r4), INTENT(IN   ) :: R0
   REAL(KIND=r4), INTENT(IN   ) :: Z1
   REAL(KIND=r4), INTENT(IN   ) :: Z2
   REAL(KIND=r4), INTENT(IN   ) :: C 
   REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)
      
   REAL(KIND=r4) :: PI 
   REAL(KIND=r4) :: RE 
   REAL(KIND=r4) :: OM 
   REAL(KIND=r4) :: GC 
   REAL(KIND=r4) :: RC 
   REAL(KIND=r4) :: T00
   REAL(KIND=r4) :: P00
   REAL(KIND=r4) :: X0
   REAL(KIND=r4) :: Y0
   REAL(KIND=r4) :: Z0
   REAL(KIND=r4) :: A
   REAL(KIND=r4) :: B
   REAL(KIND=r4) :: X
   REAL(KIND=r4) :: Y
   REAL(KIND=r4) :: Z
   REAL(KIND=r4) :: P
   REAL(KIND=r4) :: R
   REAL(KIND=r4) :: GH
   REAL(KIND=r4) :: GV
   INTEGER :: M,N,L
      
   CALL CONST (PI       ,& !REAL(KIND=r4), INTENT(OUT  ) :: PI 
               RE       ,& !REAL(KIND=r4), INTENT(OUT  ) :: RE 
               OM       ,& !REAL(KIND=r4), INTENT(OUT  ) :: OM 
               GC       ,& !REAL(KIND=r4), INTENT(OUT  ) :: GC 
               RC       ,& !REAL(KIND=r4), INTENT(OUT  ) :: RC 
               T00      ,& !REAL(KIND=r4), INTENT(OUT  ) :: T00
               P00       ) !REAL(KIND=r4), INTENT(OUT  ) :: P00

   CALL VECT  (PI/REAL(180)*A0   ,&!REAL(KIND=r4), INTENT(IN   ) :: A
               PI/REAL(180)*B0   ,&!REAL(KIND=r4), INTENT(IN   ) :: B
               X0                ,&!REAL(KIND=r4), INTENT(OUT  ) :: X
               Y0                ,&!REAL(KIND=r4), INTENT(OUT  ) :: Y
               Z0                 )!REAL(KIND=r4), INTENT(OUT  ) :: Z
   DO M=1,jDimLat
      B = REAL(180)*(REAL(M)-0.5E0)/REAL(jDimLat)-REAL(90)
      DO L=1,iDimLon

         A = REAL(360)*(REAL(L)-0.5E0)/REAL(iDimLon)

         CALL VECT (PI/REAL(180)*A   ,& !REAL(KIND=r4), INTENT(IN   ) :: A
                    PI/REAL(180)*B   ,& !REAL(KIND=r4), INTENT(IN   ) :: B
                    X                ,& !REAL(KIND=r4), INTENT(OUT  ) :: X
                    Y                ,& !REAL(KIND=r4), INTENT(OUT  ) :: Y
                    Z                 ) !REAL(KIND=r4), INTENT(OUT  ) :: Z

         P = X0*X + Y0*Y + Z0*Z
         R = REAL(180)/PI * ACOS(P)
         IF (R.LT.R0) THEN
            GH = COS(0.5E0*PI*R/R0)**2
            DO N=1,kDimLev

               CALL SIGMA (kDimLev     ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                           DZ     ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                           N      ,&!INTEGER      , INTENT(IN   ) :: N
                           Z       )!REAL(KIND=r4), INTENT(OUT  ) :: Z

               IF ((Z.GE.Z1).AND.(Z.LE.Z2)) THEN
                  GV = 0.5E0*PI*( SIN(PI*(Z-Z1)/(Z2-Z1))         &
                      +0.5E0*SIN(0.2E1*PI*(Z-Z1)/(Z2-Z1)))

                  DG (L,M,N) = DG (L,M,N) + C * GH * GV

               ENDIF
            END DO
         ENDIF
      END DO
   END DO
   RETURN
  END SUBROUTINE FORCE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE VECT (A        ,&!REAL(KIND=r4), INTENT(IN   ) :: A
                   B        ,&!REAL(KIND=r4), INTENT(IN   ) :: B
                   X        ,&!REAL(KIND=r4), INTENT(OUT  ) :: X
                   Y        ,&!REAL(KIND=r4), INTENT(OUT  ) :: Y
                   Z        ) !REAL(KIND=r4), INTENT(OUT  ) :: Z

   REAL(KIND=r4), INTENT(IN   ) :: A
   REAL(KIND=r4), INTENT(IN   ) :: B
   REAL(KIND=r4), INTENT(OUT  ) :: X
   REAL(KIND=r4), INTENT(OUT  ) :: Y
   REAL(KIND=r4), INTENT(OUT  ) :: Z

      X = COS(B) * COS(A)
      Y = COS(B) * SIN(A)
      Z = SIN(B)

   RETURN
  END SUBROUTINE VECT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE INTEE (iDimLon      ,&!INTEGER    , INTENT(IN   ) :: iDimLon
                    jDimLat      ,&!INTEGER    , INTENT(IN   ) :: jDimLat
                    kDimLev      ,&!INTEGER    , INTENT(IN   ) :: kDimLev
                    DeltaT     ,&!REAL(KIND=r4), INTENT(IN   ) :: DeltaT
                    DU      ,&!REAL(KIND=r4), INTENT(IN   ) :: DU(iDimLon,jDimLat,kDimLev)
                    DV      ,&!REAL(KIND=r4), INTENT(IN   ) :: DV(iDimLon,jDimLat,kDimLev)
                    DT      ,&!REAL(KIND=r4), INTENT(IN   ) :: DT(iDimLon,jDimLat,kDimLev)
                    DP      ,&!REAL(KIND=r4), INTENT(IN   ) :: DP(iDimLon,jDimLat)
                    UP      ,&!REAL(KIND=r4), INTENT(OUT  ) :: UP(iDimLon,jDimLat,kDimLev)
                    VP      ,&!REAL(KIND=r4), INTENT(OUT  ) :: VP(iDimLon,jDimLat,kDimLev)
                    TP      ,&!REAL(KIND=r4), INTENT(OUT  ) :: TP(iDimLon,jDimLat,kDimLev)
                    PP      ,&!REAL(KIND=r4), INTENT(OUT  ) :: PP(iDimLon,jDimLat)
                    U       ,&!REAL(KIND=r4), INTENT(INOUT) :: U (iDimLon,jDimLat,kDimLev)
                    V       ,&!REAL(KIND=r4), INTENT(INOUT) :: V (iDimLon,jDimLat,kDimLev)
                    T       ,&!REAL(KIND=r4), INTENT(INOUT) :: T (iDimLon,jDimLat,kDimLev)
                    P        )!REAL(KIND=r4), INTENT(INOUT) :: P (iDimLon,jDimLat)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: DeltaT
   REAL(KIND=r4), INTENT(IN   ) :: DU(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: DV(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: DT(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: DP(iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(OUT  ) :: UP(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: VP(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: TP(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: PP(iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(INOUT) :: U (iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: V (iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: T (iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: P (iDimLon,jDimLat)

      CALL INTEE1 (iDimLon     ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat     ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                   kDimLev     ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                   DeltaT    ,&!REAL(KIND=r4), INTENT(IN   ) :: DeltaT
                   DU     ,&!REAL(KIND=r4), INTENT(IN   ) :: DG(iDimLon,jDimLat,kDimLev)
                   UP     ,&!REAL(KIND=r4), INTENT(OUT  ) :: GP(iDimLon,jDimLat,kDimLev)
                   U       )!REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)

      CALL INTEE1 (iDimLon     ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat     ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                   kDimLev     ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                   DeltaT    ,&!REAL(KIND=r4), INTENT(IN   ) :: DeltaT
                   DV     ,&!REAL(KIND=r4), INTENT(IN   ) :: DG(iDimLon,jDimLat,kDimLev)
                   VP     ,&!REAL(KIND=r4), INTENT(OUT  ) :: GP(iDimLon,jDimLat,kDimLev)
                   V       )!REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)

      CALL INTEE1 (iDimLon     ,&!INTEGER      , INTENT(IN    ) :: iDimLon
                   jDimLat     ,&!INTEGER      , INTENT(IN    ) :: jDimLat
                   kDimLev     ,&!INTEGER      , INTENT(IN    ) :: kDimLev
                   DeltaT    ,&!REAL(KIND=r4), INTENT(IN    ) :: DeltaT
                   DT     ,&!REAL(KIND=r4), INTENT(IN    ) :: DG(iDimLon,jDimLat,kDimLev)
                   TP     ,&!REAL(KIND=r4), INTENT(OUT  ) :: GP(iDimLon,jDimLat,kDimLev)
                   T       )!REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)

      CALL INTEE1 (iDimLon     ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat     ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                   1      ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                   DeltaT    ,&!REAL(KIND=r4), INTENT(IN   ) :: DeltaT
                   DP     ,&!REAL(KIND=r4), INTENT(IN   ) :: DG(iDimLon,jDimLat,kDimLev)
                   PP     ,&!REAL(KIND=r4), INTENT(OUT  ) :: GP(iDimLon,jDimLat,kDimLev)
                   P       )!REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)

  RETURN
  END SUBROUTINE INTEE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE INTEE1 (iDimLon     ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                     jDimLat     ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                     kDimLev     ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                     DeltaT    ,&!REAL(KIND=r4), INTENT(IN   ) :: DeltaT
                     DG     ,&!REAL(KIND=r4), INTENT(IN   ) :: DG(iDimLon,jDimLat,kDimLev)
                     GP     ,&!REAL(KIND=r4), INTENT(OUT  ) :: GP(iDimLon,jDimLat,kDimLev)
                     G       )!REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: DeltaT
   REAL(KIND=r4), INTENT(IN   ) :: DG(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: GP(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(INOUT) :: G(iDimLon,jDimLat,kDimLev)
   
   REAL(KIND=r4) :: DDeltaT
   REAL(KIND=r4) :: EF 
   REAL(KIND=r4) :: EF1
   REAL(KIND=r4) :: EF2
   INTEGER       :: N,M,L

   DDeltaT = 0.2E1*DeltaT
   EF   = 0.5E-1
   EF1  = 1.0E0-EF
   EF2  = EF
   DO N=1,kDimLev
      DO M=1,jDimLat
         DO L=1,iDimLon
            G (L,M,N) = G(L,M,N) + DDeltaT * DG(L,M,N)
            GP(L,M,N) = EF1*GP(L,M,N) + EF2*G(L,M,N)
         END DO
      END DO
   END DO
   RETURN
  END SUBROUTINE INTEE1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE DIAGS (iDimLon    ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat    ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                    TPG   ,&! REAL(KIND=r4), INTENT(IN   ) :: TPG(iDimLon,jDimLat)
                    P     ,&! REAL(KIND=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
                    SLP    )! REAL(KIND=r4), INTENT(OUT  ) :: SlP(iDimLon,jDimLat)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   REAL(KIND=r4), INTENT(IN   ) :: TPG(iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(OUT  ) :: SlP(iDimLon,jDimLat)

   REAL(KIND=r4) :: GC  
   REAL(KIND=r4) :: RC  
   REAL(KIND=r4) :: T00 
   REAL(KIND=r4) :: P00 
   REAL(KIND=r4) :: RTG1
   INTEGER :: M,L
      GC   = 0.9807E1
      RC   = 0.2870E3
      T00  = 0.298E3
      P00  = 0.101325E6
      RTG1 = RC*T00/GC
   DO M=1,jDimLat
      DO L=1,iDimLon
         SLP(L,M) = P(L,M)/EXP(-TPG(L,M)/RTG1)
      END DO
   END DO
      RETURN
  END SUBROUTINE DIAGS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE DIAGZ (iDimLon        ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat        ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                    kDimLev        ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                    DZ        ,&! REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                    PTOP      ,&! REAL(KIND=r4), INTENT(IN   ) :: PTOP
                    P         ,&! REAL(KIND=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
                    H         ,&! REAL(KIND=r4), INTENT(IN   ) :: H(iDimLon,jDimLat,kDimLev)
                    Z         ) ! REAL(KIND=r4), INTENT(OUT  ) :: Z(iDimLon,jDimLat,kDimLev)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: PTOP
   REAL(KIND=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(IN   ) :: H(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: Z(iDimLon,jDimLat,kDimLev)
!  P : Surface pressure of the anomaly field.
!  H : Geopotential height of the sigma-surface.
!  Z : Diagnosed geopotential height of the pressure surface.
   REAL(KIND=r4) :: PI 
   REAL(KIND=r4) :: RE 
   REAL(KIND=r4) :: OM 
   REAL(KIND=r4) :: GC 
   REAL(KIND=r4) :: RC 
   REAL(KIND=r4) :: T00
   REAL(KIND=r4) :: P00
   REAL(KIND=r4) :: Z1
   REAL(KIND=r4) :: Z0
   REAL(KIND=r4) :: Z1Z0I
   REAL(KIND=r4) :: PR0
   REAL(KIND=r4) :: ZZ
   REAL(KIND=r4) :: ZZ1
   REAL(KIND=r4) :: ZZ2
   REAL(KIND=r4) :: PR1
   REAL(KIND=r4) :: PR2
   REAL(KIND=r4) :: C1
   REAL(KIND=r4) :: C2
   INTEGER :: M,N,L,NN
   
   CALL CONST (PI    ,&!REAL(KIND=r4), INTENT(OUT  ) :: PI 
               RE    ,&!REAL(KIND=r4), INTENT(OUT  ) :: RE 
               OM    ,&!REAL(KIND=r4), INTENT(OUT  ) :: OM 
               GC    ,&!REAL(KIND=r4), INTENT(OUT  ) :: GC 
               RC    ,&!REAL(KIND=r4), INTENT(OUT  ) :: RC 
               T00   ,&!REAL(KIND=r4), INTENT(OUT  ) :: T00
               P00    )!REAL(KIND=r4), INTENT(OUT  ) :: P00

   CALL SIGMA (kDimLev    ,&!INTEGER      , INTENT(IN   ) :: kDimLev
               DZ    ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
               kDimLev+1  ,&!INTEGER      , INTENT(IN   ) :: N
               Z1     )!REAL(KIND=r4), INTENT(OUT  ) :: Z
   CALL SIGMA (kDimLev    ,&!INTEGER      , INTENT(IN   ) :: kDimLev
               DZ    ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
               0     ,&!INTEGER      , INTENT(IN   ) :: N
               Z0     )!REAL(KIND=r4), INTENT(OUT  ) :: Z

   Z1Z0I = 1.0E0/(Z1-Z0)

   DO M=1,jDimLat
      DO L=1,iDimLon
         DO N=1,kDimLev
            CALL SIGMA (kDimLev   ,&!!INTEGER      , INTENT(IN   ) :: kDimLev
                        DZ   ,&!!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                        N    ,&!!INTEGER      , INTENT(IN   ) :: N
                        ZZ    )!!REAL(KIND=r4), INTENT(OUT  ) :: Z

            PR0 = (P00*(ZZ-Z0)+PTOP*(Z1-ZZ))*Z1Z0I

            CALL SIGMA (kDimLev   ,&!!INTEGER      , INTENT(IN   ) :: kDimLev
                        DZ   ,&!!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                        1    ,&!!INTEGER      , INTENT(IN   ) :: N
                        ZZ2   )!!REAL(KIND=r4), INTENT(OUT  ) :: Z

            DO NN=1,kDimLev-1
               ZZ1 = ZZ2

               CALL SIGMA (kDimLev    ,&!!INTEGER      , INTENT(IN   ) :: kDimLev
                           DZ    ,&!!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                           NN+1  ,&!!INTEGER      , INTENT(IN   ) :: N
                           ZZ2    ) !REAL(KIND=r4), INTENT(OUT  ) :: Z

               PR1 = (P(L,M)*(ZZ1-Z0)+PTOP*(Z1-ZZ1))*Z1Z0I
               PR2 = (P(L,M)*(ZZ2-Z0)+PTOP*(Z1-ZZ2))*Z1Z0I
               IF ((PR2.GE.PR0).OR.(NN.EQ.kDimLev-1)) THEN
                  !        C1 = ALOG(PR2/PR0)
                  C1 = (PR2-PR0)/PR0
                  !        C2 = ALOG(PR0/PR1)
                  C2 = (PR0-PR1)/PR0
                  Z(L,M,N) = (C1*H(L,M,NN)+C2*H(L,M,NN+1))/ (C1+C2)
                  EXIT
               ENDIF
            END DO
         END DO
      END DO
   END DO
   RETURN
  END SUBROUTINE DIAGZ


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ========================================
!   OUTPUT
! ========================================

  SUBROUTINE OFLOPN (iDimLon,jDimLat,kDimLev)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev

   CHARACTER(LEN=80) :: OUTFL
      WRITE(6,*) 'Output files are automatically designated:'
  101 FORMAT (1X,A21,1X,A9)
      OUTFL = 'U1.dat'
      WRITE(6,101) 'Zonal wind          =',TRIM(OUTFL)
      OPEN(11,FILE=TRIM(OUTFL),STATUS='UNKNOWN',                              &
           FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*iDimLon*(jDimLat+1)*kDimLev)
      OUTFL = 'V1.dat'
      WRITE(6,101) 'Meridional wind     =',TRIM(OUTFL)
      OPEN(12,FILE=TRIM(OUTFL),STATUS='UNKNOWN',                              &
           FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*iDimLon*(jDimLat+1)*kDimLev)
      OUTFL = 'T1.dat'
      WRITE(6,101) 'Potential temp.     =',TRIM(OUTFL)
      OPEN(13,FILE=TRIM(OUTFL),STATUS='UNKNOWN',                              &
           FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*iDimLon*(jDimLat+1)*kDimLev)
      OUTFL = 'Z1.dat'
      WRITE(6,101) 'Geopotential height =',TRIM(OUTFL)
      OPEN(14,FILE=TRIM(OUTFL),STATUS='UNKNOWN',                              &
           FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*iDimLon*(jDimLat+1)*kDimLev)
      OUTFL = 'P1.dat'
      WRITE(6,101) 'Surface pressure    =',TRIM(OUTFL)
      OPEN(15,FILE=TRIM(OUTFL),STATUS='UNKNOWN',                              &
           FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*iDimLon*(jDimLat+1))
      WRITE(6,*) '===================================='
      RETURN
  END SUBROUTINE OFLOPN
! ========================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE WRDAT (iDimLon      ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat      ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                    kDimLev      ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                    I       ,&!INTEGER      , INTENT(IN   ) :: I
                    U       ,&!REAL(KIND=r4), INTENT(IN   ) :: U(iDimLon,jDimLat,kDimLev)
                    V       ,&!REAL(KIND=r4), INTENT(IN   ) :: V(iDimLon,jDimLat,kDimLev)
                    T       ,&!REAL(KIND=r4), INTENT(IN   ) :: T(iDimLon,jDimLat,kDimLev)
                    Z       ,&!REAL(KIND=r4), INTENT(IN   ) :: Z(iDimLon,jDimLat,kDimLev)
                    P       ,&!REAL(KIND=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
                    GG      ,&!REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat+1,kDimLev)
                    PP       )!REAL(KIND=r4), INTENT(OUT  ) :: PP(iDimLon,jDimLat+1)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   INTEGER      , INTENT(IN   ) :: I
   REAL(KIND=r4), INTENT(IN   ) :: U(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: V(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: T(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: Z(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
   REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat+1,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: PP(iDimLon,jDimLat+1)
   INTEGER :: M,L
      CALL WRCNVU (iDimLon      ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat      ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                   kDimLev      ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                   U       ,&!REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                   GG       )!REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat+1,kDimLev)
      WRITE(11,REC=I) GG
      CALL WRCNVV (iDimLon      ,& ! INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat      ,& ! INTEGER      , INTENT(IN   ) :: jDimLat
                   kDimLev      ,& ! INTEGER      , INTENT(IN   ) :: kDimLev
                   V       ,& ! REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                   GG       ) ! REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat+1,kDimLev)
      WRITE(12,REC=I) GG
      CALL WRCNVG (iDimLon      ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat      ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                   kDimLev      ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                   T       ,&! REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                   GG       )! REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat+1,kDimLev)
      WRITE(13,REC=I) GG
      CALL WRCNVG (iDimLon      ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat      ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                   kDimLev      ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                   Z       ,&! REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                   GG       )! REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat+1,kDimLev)
      WRITE(14,REC=I) GG
      CALL WRCNVG (iDimLon      ,& ! INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat      ,& ! INTEGER      , INTENT(IN   ) :: jDimLat
                   1       ,& ! INTEGER      , INTENT(IN   ) :: kDimLev
                   P       ,& ! REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                   PP       ) ! REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat+1,kDimLev)
      DO M=1,jDimLat+1
         DO L=1,iDimLon
            PP(L,M) = 1.0E-2*PP(L,M)
         END DO
      END DO
      WRITE(15,REC=I) PP
      RETURN
  END SUBROUTINE WRDAT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE WRCNVG (iDimLon    ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                     jDimLat    ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                     kDimLev    ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                     G     ,&! REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                     GG     )! REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat+1,kDimLev)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat+1,kDimLev)

   INTEGER :: N,M,L
   INTEGER :: M1
   INTEGER :: M2
   INTEGER :: MM
   INTEGER :: L1
   INTEGER :: L2
   DO N=1,kDimLev
      DO M=1,jDimLat+1
         M1 = MAX0(M-1, 1)
         M2 = MIN0(M  ,jDimLat)
         MM = (jDimLat+1)+1-M
         DO L=1,iDimLon
            L1 = MOD(iDimLon+L-2,iDimLon)+1
            L2 = MOD(iDimLon+L-1,iDimLon)+1
            GG(L,MM,N) = 0.25E0*( G(L1,M1,N)+G(L2,M1,N) + G(L1,M2,N)+G(L2,M2,N))
         END DO
      END DO
   END DO
   RETURN
  END SUBROUTINE WRCNVG
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  SUBROUTINE WRCNVV (iDimLon   ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                     jDimLat   ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                     kDimLev   ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                     G    ,&! REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                     GG    )! REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat+1,kDimLev)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat+1,kDimLev)

   INTEGER :: N,M,L
   INTEGER :: M1
   INTEGER :: M2
   INTEGER :: MM
   INTEGER :: L1
   INTEGER :: L2
   
   DO N=1,kDimLev
      DO L=1,iDimLon
         GG(L,   1,N) = 0.0E0
         GG(L,jDimLat+1,N) = 0.0E0
      END DO
      
      DO M=2,jDimLat
         M1 = M-1
         M2 = M
         MM = (jDimLat+1)+1-M
         DO L=1,iDimLon
            L1 = MOD(iDimLon+L-2,iDimLon)+1
            L2 = MOD(iDimLon+L-1,iDimLon)+1
            GG(L,MM,N) = 0.25E0*( G(L1,M1,N)+G(L2,M1,N)+G(L1,M2,N)+G(L2,M2,N))
         END DO
      END DO
   END DO

      RETURN
  END SUBROUTINE WRCNVV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE WRCNVU (iDimLon        ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                     jDimLat        ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                     kDimLev        ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                     G         ,&!REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
                     GG         )!REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat+1,kDimLev)
   INTEGER      , INTENT(IN   ) :: iDimLon
   INTEGER      , INTENT(IN   ) :: jDimLat
   INTEGER      , INTENT(IN   ) :: kDimLev
   REAL(KIND=r4), INTENT(IN   ) :: G(iDimLon,jDimLat,kDimLev)
   REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat+1,kDimLev)

   INTEGER :: N,M,L
   INTEGER :: M1
   INTEGER :: M2
   INTEGER :: MM
   INTEGER :: LL

   DO N=1,kDimLev
      DO M=1,jDimLat+1
         M1 = MAX0(M-1, 1)
         M2 = MIN0(M  ,jDimLat)
         MM = (jDimLat+1)+1-M
         DO L=1,iDimLon
            LL = MOD(iDimLon+L-1,iDimLon)+1
            GG(L,MM,N) = 0.5E0*(G(LL,M1,N)+G(LL,M2,N))
         END DO
      END DO
   END DO
   RETURN
  END SUBROUTINE WRCNVU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MODEL_SGCM

PROGRAM Main
 USE MODEL_SGCM, ONLY: Init_SGCM,PARAM,r4,jDimLat,kDimLev,iDimLon,&
                        BASIC,POLES,DIAGX,DIAGH,DIAGR,TPGEV,TENDZ,PGRAD,INTEL,Diffusion,INTEE,DIAGS,WRDAT,&
                        BSCPRN,DiagVertSigVelocity,ADVEC,RSGN,ZERO,INITC,OFLOPN,Convergence,CORIO,ADBTC,Damping,FORCE,DIAGZ
 IMPLICIT NONE

 !INTEGER, PARAMETER :: iDimLon=144
 !INTEGER, PARAMETER :: jDimLat=72
 !INTEGER, PARAMETER :: kDimLev=11

 INTEGER, PARAMETER :: IX=480*100
 INTEGER, PARAMETER :: IO=480
 REAL(KIND=r4),PARAMETER :: DeltaT=0.180E3
 REAL(KIND=r4) :: PTOP
 REAL(KIND=r4), ALLOCATABLE :: DX(:)
 REAL(KIND=r4), ALLOCATABLE :: DY(:)
 REAL(KIND=r4), ALLOCATABLE :: FC(:)
 REAL(KIND=r4), ALLOCATABLE :: DZ(:)
 REAL(KIND=r4), ALLOCATABLE :: PLEV(:)
      
 REAL(KIND=r4), ALLOCATABLE :: TPG(:,:)
 REAL(KIND=r4), ALLOCATABLE :: TPG1(:,:)
     
 REAL(KIND=r4), ALLOCATABLE :: U0(:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: V0(:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: T0(:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: P0(:,:)
     
 REAL(KIND=r4), ALLOCATABLE :: U (:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: V (:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: T (:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: P (:,:)
     
 REAL(KIND=r4), ALLOCATABLE :: UP(:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: VP(:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: TP(:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: PP(:,:)
     
 REAL(KIND=r4), ALLOCATABLE :: DU(:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: DV(:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: DT(:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: DP(:,:)
     
 REAL(KIND=r4), ALLOCATABLE :: UF(:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: VF(:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: TF(:,:,:)
     
 REAL(KIND=r4), ALLOCATABLE :: H0(:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: H (:,:,:)
     
 REAL(KIND=r4), ALLOCATABLE :: W0(:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: W (:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: Z (:,:,:)
     
 REAL(KIND=r4), ALLOCATABLE :: EX(:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: Scale_height(:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: SLP(:,:)
 REAL(KIND=r4), ALLOCATABLE :: DM1(:,:,:)
 REAL(KIND=r4), ALLOCATABLE :: DM2(:,:)
 
 CALL  Init()

 
CONTAINS

 SUBROUTINE Init()
  IMPLICIT NONE
  INTEGER :: I
  REAL    :: PLEV0
  REAL    :: TIM

  CALL  Init_SGCM

  ALLOCATE(  DX(jDimLat) )
  ALLOCATE(  DY(jDimLat) )
  ALLOCATE(  FC(jDimLat) )
  ALLOCATE(  DZ(kDimLev) );DZ=0.0_r4
  ALLOCATE(  PLEV(kDimLev) )
  
  ALLOCATE(  TPG(iDimLon,jDimLat)   )
  ALLOCATE(  TPG1(iDimLon,jDimLat)  )
 
  ALLOCATE(  U0(iDimLon,jDimLat,kDimLev) );U0=0.0_r4
  ALLOCATE(  V0(iDimLon,jDimLat,kDimLev) );V0=0.0_r4
  ALLOCATE(  T0(iDimLon,jDimLat,kDimLev) );T0=0.0_r4
  ALLOCATE(  P0(iDimLon,jDimLat)    );P0=0.0_r4
 
  ALLOCATE(  U (iDimLon,jDimLat,kDimLev) )
  ALLOCATE(  V (iDimLon,jDimLat,kDimLev) )
  ALLOCATE(  T (iDimLon,jDimLat,kDimLev) )
  ALLOCATE(  P (iDimLon,jDimLat) )
 
  ALLOCATE(  UP(iDimLon,jDimLat,kDimLev) )
  ALLOCATE(  VP(iDimLon,jDimLat,kDimLev) )
  ALLOCATE(  TP(iDimLon,jDimLat,kDimLev) )
  ALLOCATE(  PP(iDimLon,jDimLat) )
 
  ALLOCATE(  DU(iDimLon,jDimLat,kDimLev) )
  ALLOCATE(  DV(iDimLon,jDimLat,kDimLev) )
  ALLOCATE(  DT(iDimLon,jDimLat,kDimLev) )
  ALLOCATE(  DP(iDimLon,jDimLat) )
 
  ALLOCATE(  UF(iDimLon,jDimLat,kDimLev) )
  ALLOCATE(  VF(iDimLon,jDimLat,kDimLev) )
  ALLOCATE(  TF(iDimLon,jDimLat,kDimLev) )
 
  ALLOCATE(  H0(iDimLon,jDimLat,kDimLev) )
  ALLOCATE(  H (iDimLon,jDimLat,kDimLev) )
 
  ALLOCATE(  W0(iDimLon,jDimLat,kDimLev) )
  ALLOCATE(  W (iDimLon,jDimLat,kDimLev) )
  ALLOCATE(  Z (iDimLon,jDimLat,kDimLev) )
 
  ALLOCATE(  EX(iDimLon,jDimLat,kDimLev) )
  ALLOCATE(  Scale_height(iDimLon,jDimLat,kDimLev)   )
  ALLOCATE(  SLP(iDimLon,jDimLat)      )
  ALLOCATE(  DM1(iDimLon,jDimLat+1,kDimLev) )
  ALLOCATE(  DM2(iDimLon,jDimLat+1)    )

! ########################################
!   1. PREPARATION
! ########################################
  CALL PARAM (&
             jDimLat      ,& !INTEGER      , INTENT(IN   ) :: kDimLev
             kDimLev      ,& !INTEGER      , INTENT(IN   ) :: iDimLon
             DX      ,& !REAL(KIND=r4), INTENT(OUT  ) :: DX(jDimLat)
             DY      ,& !REAL(KIND=r4), INTENT(OUT  ) :: DY(jDimLat)
             DZ      ,& !REAL(KIND=r4), INTENT(OUT  ) :: DZ(kDimLev)
             PLEV0   ,& !REAL(KIND=r4), INTENT(OUT  ) :: PLEV0
             PLEV    ,& !REAL(KIND=r4), INTENT(OUT  ) :: PLEV(kDimLev)
             FC   )     !REAL(KIND=r4), INTENT(OUT  ) :: FC(jDimLat)!componente horizontal da força de Coriolis

  CALL BASIC (iDimLon    ,&!INTEGER      , INTENT(IN   ) :: iDimLon
              jDimLat    ,&!INTEGER      , INTENT(IN   ) :: jDimLat
              kDimLev    ,&!INTEGER      , INTENT(IN   ) :: kDimLev
              DZ    ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
              PLEV0 ,&!REAL(KIND=r4), INTENT(IN   ) :: PLEV0
              PTOP  ,&!REAL(KIND=r4), INTENT(OUT  ) :: PTOP
              U0    ,&!REAL(KIND=r4), INTENT(OUT  ) :: U(iDimLon,jDimLat,kDimLev)
              V0    ,&!REAL(KIND=r4), INTENT(OUT  ) :: V(iDimLon,jDimLat,kDimLev)
              T0    ,&!REAL(KIND=r4), INTENT(OUT  ) :: T(iDimLon,jDimLat,kDimLev)
              P0    ,&!REAL(KIND=r4), INTENT(OUT  ) :: P(iDimLon,jDimLat)
              TPG    )!REAL(KIND=r4), INTENT(OUT  ) :: H(iDimLon,jDimLat)
 
      CALL POLES (iDimLon    ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                  jDimLat    ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                  kDimLev    ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                  U0    ,&! REAL(KIND=r4), INTENT(INOUT) :: U(iDimLon,jDimLat,kDimLev)
                  V0    ,&! REAL(KIND=r4), INTENT(INOUT) :: V(iDimLon,jDimLat,kDimLev)
                  T0    ,&! REAL(KIND=r4), INTENT(INOUT) :: T(iDimLon,jDimLat,kDimLev)
                  P0    ,&! REAL(KIND=r4), INTENT(INOUT) :: P(iDimLon,jDimLat)
                  TPG    )! REAL(KIND=r4), INTENT(INOUT) :: H(iDimLon,jDimLat)

      CALL DIAGX (iDimLon   ,&!  INTEGER      , INTENT(IN   ) :: iDimLon
                  jDimLat   ,&!  INTEGER      , INTENT(IN   ) :: jDimLat
                  kDimLev   ,&!  INTEGER      , INTENT(IN   ) :: kDimLev
                  DZ   ,&!  REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                  PTOP ,&!  REAL(KIND=r4), INTENT(IN   ) :: PTOP
                  P0   ,&!  REAL(KIND=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
                  EX    )!  REAL(KIND=r4), INTENT(OUT  ) :: EX(iDimLon,jDimLat,kDimLev)

      CALL DIAGH (iDimLon   ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                  jDimLat   ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                  kDimLev   ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DZ   ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ (kDimLev)
                  TPG  ,&!REAL(KIND=r4), INTENT(IN   ) :: TPG(iDimLon,jDimLat)
                  PTOP ,&!REAL(KIND=r4), INTENT(IN   ) :: PTOP
                  P0   ,&!REAL(KIND=r4), INTENT(IN   ) :: P  (iDimLon,jDimLat)
                  EX   ,&!REAL(KIND=r4), INTENT(IN   ) :: EX (iDimLon,jDimLat,kDimLev)
                  T0   ,&!REAL(KIND=r4), INTENT(IN   ) :: T  (iDimLon,jDimLat,kDimLev)
                  H0    )!REAL(KIND=r4), INTENT(OUT  ) :: H  (iDimLon,jDimLat,kDimLev)
       
      CALL DIAGR (iDimLon   ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                  jDimLat   ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                  kDimLev   ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DZ   ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ  (kDimLev)
                  PTOP ,&!REAL(KIND=r4), INTENT(IN   ) :: PTOP
                  P0   ,&!REAL(KIND=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
                  T0   ,&!REAL(KIND=r4), INTENT(IN   ) :: T(iDimLon,jDimLat,kDimLev)
                  Scale_height   )!REAL(KIND=r4), INTENT(OUT  ) :: Scale_height (iDimLon,jDimLat,kDimLev)

      CALL BSCPRN (iDimLon   ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                   jDimLat   ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                   kDimLev   ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                   DZ   ,&! REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                   PLEV ,&! REAL(KIND=r4), INTENT(IN   ) :: PLEV(kDimLev)
                   H0   ,&! REAL(KIND=r4), INTENT(IN   ) :: H(iDimLon,jDimLat,kDimLev)
                   U0   ,&! REAL(KIND=r4), INTENT(IN   ) :: U(iDimLon,jDimLat,kDimLev)
                   T0    )! REAL(KIND=r4), INTENT(IN   ) :: T(iDimLon,jDimLat,kDimLev)


! #### Gravity wave drag and adiabatic heating ####

      CALL DiagVertSigVelocity (iDimLon   ,& !INTEGER, INTENT(IN   ) :: iDimLon
                  jDimLat   ,& !INTEGER, INTENT(IN   ) :: jDimLat
                  kDimLev   ,& !INTEGER, INTENT(IN   ) :: kDimLev
                  DX   ,& !REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
                  DY   ,& !REAL(KIND=r4), INTENT(IN   ) :: DY  (jDimLat)
                  DZ   ,& !REAL(KIND=r4), INTENT(IN   ) :: DZ  (kDimLev)
                  PTOP ,& !REAL(KIND=r4), INTENT(IN   ) :: PTOP
                  P0   ,& !REAL(KIND=r4), INTENT(IN   ) :: P   (iDimLon,jDimLat)
                  U0   ,& !REAL(KIND=r4), INTENT(IN   ) :: U   (iDimLon,jDimLat,kDimLev)
                  V0   ,& !REAL(KIND=r4), INTENT(IN   ) :: V   (iDimLon,jDimLat,kDimLev)
                  W     ) !REAL(KIND=r4), INTENT(OUT  ) :: W   (iDimLon,jDimLat,kDimLev)

      CALL ADVEC (iDimLon   ,& ! INTEGER      , INTENT(IN   ) :: iDimLon 
                  jDimLat   ,& ! INTEGER      , INTENT(IN   ) :: jDimLat 
                  kDimLev   ,& ! INTEGER      , INTENT(IN   ) :: kDimLev 
                  DX   ,& ! REAL(KIND=r4), INTENT(IN   ) ::  DX(jDimLat)
                  DY   ,& ! REAL(KIND=r4), INTENT(IN   ) ::  DY(jDimLat)
                  DZ   ,& ! REAL(KIND=r4), INTENT(IN   ) ::  DZ(kDimLev)
                  U0   ,& ! REAL(KIND=r4), INTENT(IN   ) ::  U (iDimLon,jDimLat,kDimLev)
                  V0   ,& ! REAL(KIND=r4), INTENT(IN   ) ::  V (iDimLon,jDimLat,kDimLev)
                  W    ,& ! REAL(KIND=r4), INTENT(IN   ) ::  W (iDimLon,jDimLat,kDimLev)
                  T    ,& ! REAL(KIND=r4), INTENT(IN   ) ::  T (iDimLon,jDimLat,kDimLev)
                  UF   ,& ! REAL(KIND=r4), INTENT(INOUT) ::  DU(iDimLon,jDimLat,kDimLev)
                  VF   ,& ! REAL(KIND=r4), INTENT(INOUT) ::  DV(iDimLon,jDimLat,kDimLev)
                  TF    ) ! REAL(KIND=r4), INTENT(INOUT) ::  DT(iDimLon,jDimLat,kDimLev)

      CALL RSGN  (iDimLon   ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                  jDimLat   ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                  kDimLev   ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  UF    )!REAL(KIND=r4), INTENT(INOUT) ::   G(iDimLon,jDimLat,kDimLev)

      CALL RSGN  (iDimLon   ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                  jDimLat   ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                  kDimLev   ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  VF    )!REAL(KIND=r4), INTENT(INOUT) ::   G(iDimLon,jDimLat,kDimLev)

      CALL RSGN  (iDimLon   ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                  jDimLat   ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                  kDimLev   ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  TF    )!REAL(KIND=r4), INTENT(INOUT) ::   G(iDimLon,jDimLat,kDimLev)

! #### Isothermal static atmosphere ####

      CALL ZERO  (iDimLon    ,&!  INTEGER      , INTENT(IN   ) :: iDimLon
                  jDimLat    ,&!  INTEGER      , INTENT(IN   ) :: jDimLat
                  1     ,&!  INTEGER      , INTENT(IN   ) :: kDimLev
                  TPG1   )!  REAL(KIND=r4), INTENT(OUT  ) :: G(iDimLon,jDimLat,kDimLev)

      CALL INITC (iDimLon    ,&! INTEGER     , INTENT(IN    ) :: iDimLon
                  jDimLat    ,&! INTEGER     , INTENT(IN    ) :: jDimLat
                  kDimLev    ,&! INTEGER     , INTENT(IN    ) :: kDimLev
                  PLEV  ,&! REAL(KIND=r4), INTENT(IN   ) :: PLEV(kDimLev)
                  TPG1  ,&! REAL(KIND=r4), INTENT(IN   ) :: TPG(iDimLon,jDimLat)
                  U     ,&! REAL(KIND=r4), INTENT(OUT  ) :: U (iDimLon,jDimLat,kDimLev)
                  V     ,&! REAL(KIND=r4), INTENT(OUT  ) :: V(iDimLon,jDimLat,kDimLev)
                  T     ,&! REAL(KIND=r4), INTENT(OUT  ) :: T(iDimLon,jDimLat,kDimLev)
                  P     ,&! REAL(KIND=r4), INTENT(OUT  ) :: P(iDimLon,jDimLat)
                  UP    ,&! REAL(KIND=r4), INTENT(OUT  ) :: UP(iDimLon,jDimLat,kDimLev)
                  VP    ,&! REAL(KIND=r4), INTENT(OUT  ) :: VP(iDimLon,jDimLat,kDimLev)
                  TP    ,&! REAL(KIND=r4), INTENT(OUT  ) :: TP(iDimLon,jDimLat,kDimLev)
                  PP     )! REAL(KIND=r4), INTENT(OUT  ) :: PP(iDimLon,jDimLat)

!   Output files are opened.

      CALL OFLOPN (iDimLon    ,&!INTEGER, INTENT(IN   ) :: iDimLon
                   jDimLat    ,&!INTEGER, INTENT(IN   ) :: jDimLat
                   kDimLev     )!INTEGER, INTENT(IN   ) :: kDimLev

! #### Do loop 11 starts here ####

!   In roop 11, the time derivateves of U, V, T and P are calculated,
!   and time-integrations are performed for them.

      WRITE(6,*) 'Time-integration:'
      WRITE(6,'(1X,A22,1X,F9.3,1X,A4)')                  &
                                'Time step            =',DeltaT,'sec.'
      WRITE(6,'(1X,A22,1X,I6)') 'Number of time steps =',IX
      WRITE(6,'(1X,A22,1X,I6)') 'Interval of output   =',IO
      WRITE(6,'(1X,A22,1X,I6)') 'Number of outputs    =',IX/IO
      WRITE(6,*) '===================================='
      WRITE(6,*) 'Time-integration started...'

      DO 11 I=1,IX

      TIM = DeltaT * REAL(I-1)
!      IF (MOD(I,40).EQ.0) WRITE(6,'(A13,1X,I6)') 'Time step # =',I


      CALL TPGEV (iDimLon        ,& ! INTEGER      , INTENT(IN   ) :: iDimLon
                  jDimLat        ,& ! INTEGER      , INTENT(IN   ) :: jDimLat
                  kDimLev        ,& ! INTEGER      , INTENT(IN   ) :: kDimLev
                  Scale_height       ,& ! REAL(KIND=r4), INTENT(IN   ) :: Scale_height(iDimLon,jDimLat,kDimLev)
                  TPG       ,& ! REAL(KIND=r4), INTENT(IN   ) :: TPG(iDimLon,jDimLat)
                  TIM       ,& ! REAL(KIND=r4), INTENT(IN   ) :: TIM
                  TPG1      ,& ! REAL(KIND=r4), INTENT(INOUT) :: TPG1(iDimLon,jDimLat)
                  P         ,& ! REAL(KIND=r4), INTENT(INOUT) :: P(iDimLon,jDimLat)
                  PP         ) ! REAL(KIND=r4), INTENT(INOUT) :: PP(iDimLon,jDimLat)

! ########################################
!   2. CALCULATION OF DEPENDENT VARIABLES
! ########################################

      CALL DiagVertSigVelocity (iDimLon    ,&!INTEGER, INTENT(IN   ) :: iDimLon
                  jDimLat    ,&!INTEGER, INTENT(IN   ) :: jDimLat
                  kDimLev    ,&!INTEGER, INTENT(IN   ) :: kDimLev
                  DX    ,&!REAL(KIND=r4), INTENT(IN   ) :: DX(jDimLat)
                  DY    ,&!REAL(KIND=r4), INTENT(IN   ) :: DY  (jDimLat)
                  DZ    ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ  (kDimLev)
                  PTOP  ,&!REAL(KIND=r4), INTENT(IN   ) :: PTOP
                  P     ,&!REAL(KIND=r4), INTENT(IN   ) :: P   (iDimLon,jDimLat)
                  U     ,&!REAL(KIND=r4), INTENT(IN   ) :: U   (iDimLon,jDimLat,kDimLev)
                  V     ,&!REAL(KIND=r4), INTENT(IN   ) :: V   (iDimLon,jDimLat,kDimLev)
                  W      )!REAL(KIND=r4), INTENT(OUT  ) :: W   (iDimLon,jDimLat,kDimLev)

      CALL DIAGH (iDimLon     ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                  jDimLat     ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                  kDimLev     ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DZ     ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ (kDimLev)
                  TPG1   ,&!REAL(KIND=r4), INTENT(IN   ) :: TPG(iDimLon,jDimLat)
                  PTOP   ,&!REAL(KIND=r4), INTENT(IN   ) :: PTOP
                  P      ,&!REAL(KIND=r4), INTENT(IN   ) :: P  (iDimLon,jDimLat)
                  EX     ,&!REAL(KIND=r4), INTENT(IN   ) :: EX (iDimLon,jDimLat,kDimLev)
                  T      ,&!REAL(KIND=r4), INTENT(IN   ) :: T  (iDimLon,jDimLat,kDimLev)
                  H       )!REAL(KIND=r4), INTENT(OUT  ) :: H  (iDimLon,jDimLat,kDimLev)

! ########################################
!   3. DYNAMICAL PROCESSES
! ########################################

      CALL TENDZ (iDimLon      ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                  jDimLat      ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                  kDimLev      ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DU      ,&!REAL(KIND=r4), INTENT(OUT  ) :: DU(iDimLon,jDimLat,kDimLev)
                  DV      ,&!REAL(KIND=r4), INTENT(OUT  ) :: DV(iDimLon,jDimLat,kDimLev)
                  DT      ,&!REAL(KIND=r4), INTENT(OUT  ) :: DT(iDimLon,jDimLat,kDimLev)
                  DP       )!REAL(KIND=r4), INTENT(OUT  ) :: DP(iDimLon,jDimLat)

      CALL Convergence (iDimLon      ,&!INTEGER     , INTENT(IN ) ::  iDimLon
                  jDimLat      ,&!INTEGER     , INTENT(IN ) ::  jDimLat
                  kDimLev      ,&!INTEGER     , INTENT(IN ) ::  kDimLev
                  DX      ,&!REAL(KIND=r4), INTENT(IN  ) ::  DX(jDimLat)
                  DY      ,&!REAL(KIND=r4), INTENT(IN  ) ::  DY(jDimLat)
                  DZ      ,&!REAL(KIND=r4), INTENT(IN  ) ::  DZ(kDimLev)
                  PTOP    ,&!REAL(KIND=r4), INTENT(IN  ) ::  PTOP
                  P       ,&!REAL(KIND=r4), INTENT(IN  ) ::  P(iDimLon,jDimLat)
                  U       ,&!REAL(KIND=r4), INTENT(IN  ) ::  U(iDimLon,jDimLat,kDimLev)
                  V       ,&!REAL(KIND=r4), INTENT(IN  ) ::  V(iDimLon,jDimLat,kDimLev)
                  DP       )!REAL(KIND=r4), INTENT(INOUT) ::  DP(iDimLon,jDimLat)

      CALL ADVEC (iDimLon      ,& ! INTEGER      , INTENT(IN   ) :: iDimLon 
                  jDimLat      ,& ! INTEGER      , INTENT(IN   ) :: jDimLat 
                  kDimLev      ,& ! INTEGER      , INTENT(IN   ) :: kDimLev 
                  DX      ,& ! REAL(KIND=r4), INTENT(IN   ) ::  DX(jDimLat)
                  DY      ,& ! REAL(KIND=r4), INTENT(IN   ) ::  DY(jDimLat)
                  DZ      ,& ! REAL(KIND=r4), INTENT(IN   ) ::  DZ(kDimLev)
                  U       ,& ! REAL(KIND=r4), INTENT(IN   ) ::  U (iDimLon,jDimLat,kDimLev)
                  V       ,& ! REAL(KIND=r4), INTENT(IN   ) ::  V (iDimLon,jDimLat,kDimLev)
                  W       ,& ! REAL(KIND=r4), INTENT(IN   ) ::  W (iDimLon,jDimLat,kDimLev)
                  T       ,& ! REAL(KIND=r4), INTENT(IN   ) ::  T (iDimLon,jDimLat,kDimLev)
                  DU      ,& ! REAL(KIND=r4), INTENT(INOUT) ::  DU(iDimLon,jDimLat,kDimLev)
                  DV      ,& ! REAL(KIND=r4), INTENT(INOUT) ::  DV(iDimLon,jDimLat,kDimLev)
                  DT       ) ! REAL(KIND=r4), INTENT(INOUT) ::  DT(iDimLon,jDimLat,kDimLev)

      CALL PGRAD (iDimLon      ,& !INTEGER      , INTENT(IN   ) :: iDimLon
                  jDimLat      ,& !INTEGER      , INTENT(IN   ) :: jDimLat
                  kDimLev      ,& !INTEGER      , INTENT(IN   ) :: kDimLev
                  DX      ,& !REAL(KINd=r4), INTENT(IN   ) :: DX(jDimLat)
                  DY      ,& !REAL(KINd=r4), INTENT(IN   ) :: DY(jDimLat)
                  P       ,& !REAL(KINd=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
                  H       ,& !REAL(KINd=r4), INTENT(IN   ) :: H(iDimLon,jDimLat,kDimLev)
                  Scale_height     ,& !REAL(KINd=r4), INTENT(IN   ) :: Scale_height(iDimLon,jDimLat,kDimLev)
                  DU      ,& !REAL(KINd=r4), INTENT(INOUT) :: DU(iDimLon,jDimLat,kDimLev)
                  DV       ) !REAL(KINd=r4), INTENT(INOUT) :: DV(iDimLon,jDimLat,kDimLev)

      CALL CORIO (iDimLon      ,&! INTEGER         , INTENT(IN   ) :: iDimLon
                  jDimLat      ,&! INTEGER         , INTENT(IN   ) :: jDimLat
                  kDimLev      ,&! INTEGER         , INTENT(IN   ) :: kDimLev
                  'G'     ,&! CHARACTER(LEN=1), INTENT(IN   ) :: CG
                  FC      ,&! REAL(KIND=r4)   , INTENT(IN   ) :: FC(jDimLat)
                  U       ,&! REAL(KIND=r4)   , INTENT(IN   ) :: U(iDimLon,jDimLat,kDimLev)
                  V       ,&! REAL(KIND=r4)   , INTENT(IN   ) :: V(iDimLon,jDimLat,kDimLev)
                  DU      ,&! REAL(KIND=r4)   , INTENT(INOUT) :: DU(iDimLon,jDimLat,kDimLev)
                  DV       )! REAL(KIND=r4)   , INTENT(INOUT) :: DV(iDimLon,jDimLat,kDimLev)

      CALL INTEL (iDimLon      ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                  jDimLat      ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                  kDimLev      ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DeltaT     ,&!REAL(KIND=r4), INTENT(IN   ) :: DeltaT
                  DU      ,&!REAL(KIND=r4), INTENT(IN   ) :: DU(iDimLon,jDimLat,kDimLev)
                  DV      ,&!REAL(KIND=r4), INTENT(IN   ) :: DV(iDimLon,jDimLat,kDimLev)
                  DT      ,&!REAL(KIND=r4), INTENT(IN   ) :: DT(iDimLon,jDimLat,kDimLev)
                  DP      ,&!REAL(KIND=r4), INTENT(IN   ) :: DP(iDimLon,jDimLat)
                  UP      ,&!REAL(KIND=r4), INTENT(INOUT) :: UP(iDimLon,jDimLat,kDimLev)
                  VP      ,&!REAL(KIND=r4), INTENT(INOUT) :: VP(iDimLon,jDimLat,kDimLev)
                  TP      ,&!REAL(KIND=r4), INTENT(INOUT) :: TP(iDimLon,jDimLat,kDimLev)
                  PP      ,&!REAL(KIND=r4), INTENT(INOUT) :: PP(iDimLon,jDimLat),
                  U       ,&!REAL(KIND=r4), INTENT(INOUT) :: U (iDimLon,jDimLat,kDimLev)
                  V       ,&!REAL(KIND=r4), INTENT(INOUT) :: V (iDimLon,jDimLat,kDimLev)
                  T       ,&!REAL(KIND=r4), INTENT(INOUT) :: T (iDimLon,jDimLat,kDimLev)
                  P        )!REAL(KIND=r4), INTENT(INOUT) :: P (iDimLon,jDimLat),

      CALL POLES (iDimLon      ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                  jDimLat      ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                  kDimLev      ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                  U       ,&! REAL(KIND=r4), INTENT(INOUT) :: U(iDimLon,jDimLat,kDimLev)
                  V       ,&! REAL(KIND=r4), INTENT(INOUT) :: V(iDimLon,jDimLat,kDimLev)
                  T       ,&! REAL(KIND=r4), INTENT(INOUT) :: T(iDimLon,jDimLat,kDimLev)
                  P       ,&! REAL(KIND=r4), INTENT(INOUT) :: P(iDimLon,jDimLat)
                  TPG1     )! REAL(KIND=r4), INTENT(INOUT) :: H(iDimLon,jDimLat)

! ########################################
!   4. PHYSICAL PROCESSES
! ########################################

      CALL TENDZ (iDimLon      ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                  jDimLat      ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                  kDimLev      ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  DU      ,&!REAL(KIND=r4), INTENT(OUT  ) :: DU(iDimLon,jDimLat,kDimLev)
                  DV      ,&!REAL(KIND=r4), INTENT(OUT  ) :: DV(iDimLon,jDimLat,kDimLev)
                  DT      ,&!REAL(KIND=r4), INTENT(OUT  ) :: DT(iDimLon,jDimLat,kDimLev)
                  DP       )!REAL(KIND=r4), INTENT(OUT  ) :: DP(iDimLon,jDimLat)

      CALL Diffusion (iDimLon     ,& ! INTEGER      , INTENT(IN  ) :: iDimLon
                  jDimLat     ,& ! INTEGER      , INTENT(IN  ) :: jDimLat
                  kDimLev     ,& ! INTEGER      , INTENT(IN  ) :: kDimLev
                  DX     ,& ! REAL(KIND=r4), INTENT(IN  ) :: DX(jDimLat)
                  DY     ,& ! REAL(KIND=r4), INTENT(IN  ) :: DY(jDimLat)
                  DZ     ,& ! REAL(KIND=r4), INTENT(IN  ) :: DZ(kDimLev)
                  Scale_height    ,& ! REAL(KIND=r4), INTENT(IN  ) :: Scale_height(iDimLon,jDimLat,kDimLev
                  U      ,& ! REAL(KIND=r4), INTENT(IN  ) :: U (iDimLon,jDimLat,kDimLev)
                  V      ,& ! REAL(KIND=r4), INTENT(IN  ) :: V (iDimLon,jDimLat,kDimLev)
                  T      ,& ! REAL(KIND=r4), INTENT(IN  ) :: T (iDimLon,jDimLat,kDimLev)
                  DU     ,& ! REAL(KIND=r4), INTENT(INOUT) :: DU(iDimLon,jDimLat,kDimLev)
                  DV     ,& ! REAL(KIND=r4), INTENT(INOUT) :: DV(iDimLon,jDimLat,kDimLev)
                  DT      ) ! REAL(KIND=r4), INTENT(INOUT) :: DT(iDimLon,jDimLat,kDimLev)

      CALL ADBTC (iDimLon     ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                  jDimLat     ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                  kDimLev     ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  UF     ,&!REAL(KIND=r4), INTENT(IN  ) :: UF(iDimLon,jDimLat,kDimLev)
                  VF     ,&!REAL(KIND=r4), INTENT(IN  ) :: VF(iDimLon,jDimLat,kDimLev)
                  TF     ,&!REAL(KIND=r4), INTENT(IN  ) :: TF(iDimLon,jDimLat,kDimLev)
                  DU     ,&!REAL(KIND=r4), INTENT(INOUT) :: DU(iDimLon,jDimLat,kDimLev)
                  DV     ,&!REAL(KIND=r4), INTENT(INOUT) :: DV(iDimLon,jDimLat,kDimLev)
                  DT      )!REAL(KIND=r4), INTENT(INOUT) :: DT(iDimLon,jDimLat,kDimLev)

      CALL Damping (iDimLon     ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                  jDimLat     ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                  kDimLev     ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                  U0     ,&!REAL(KIND=r4), INTENT(IN   ) :: U0(iDimLon,jDimLat,kDimLev)
                  V0     ,&!REAL(KIND=r4), INTENT(IN   ) :: V0(iDimLon,jDimLat,kDimLev)
                  T0     ,&!REAL(KIND=r4), INTENT(IN   ) :: T0(iDimLon,jDimLat,kDimLev)
                  U      ,&!REAL(KIND=r4), INTENT(IN   ) :: U (iDimLon,jDimLat,kDimLev)
                  V      ,&!REAL(KIND=r4), INTENT(IN   ) :: V (iDimLon,jDimLat,kDimLev)
                  T      ,&!REAL(KIND=r4), INTENT(IN   ) :: T (iDimLon,jDimLat,kDimLev)
                  DU     ,&!REAL(KIND=r4), INTENT(INOUT) :: DU(iDimLon,jDimLat,kDimLev)
                  DV     ,&!REAL(KIND=r4), INTENT(INOUT) :: DV(iDimLon,jDimLat,kDimLev)
                  DT      )!REAL(KIND=r4), INTENT(INOUT) :: DT(iDimLon,jDimLat,kDimLev)


! #### Subroutine for response problem ####
!   For response problem, remove the comment of the following line.
!      CALL FORCE (iDimLon          ,&!INTEGER      , INTENT(IN   ) :: iDimLon
!                  jDimLat          ,&!INTEGER      , INTENT(IN   ) :: jDimLat
!                  kDimLev          ,&!INTEGER      , INTENT(IN   ) :: kDimLev
!                  DZ          ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
!                  20.         ,&!REAL(KIND=r4), INTENT(IN   ) :: B0
!                  130.        ,&!REAL(KIND=r4), INTENT(IN   ) :: A0
!                  10.         ,&!REAL(KIND=r4), INTENT(IN   ) :: R0
!                  0.15        ,&!REAL(KIND=r4), INTENT(IN   ) :: Z1
!                  0.75        ,&!REAL(KIND=r4), INTENT(IN   ) :: Z2
!                  2.*1.0E-5   ,&!REAL(KIND=r4), INTENT(IN   ) :: C 
!                  DT           )!REAL(KIND=r4), INTENT(INOUT) :: DG(iDimLon,jDimLat,kDimLev)

      CALL INTEE (iDimLon     ,&!INTEGER    , INTENT(IN   ) :: iDimLon
                  jDimLat     ,&!INTEGER    , INTENT(IN   ) :: jDimLat
                  kDimLev     ,&!INTEGER    , INTENT(IN   ) :: kDimLev
                  DeltaT    ,&!REAL(KIND=r4), INTENT(IN   ) :: DeltaT
                  DU     ,&!REAL(KIND=r4), INTENT(IN   ) :: DU(iDimLon,jDimLat,kDimLev)
                  DV     ,&!REAL(KIND=r4), INTENT(IN   ) :: DV(iDimLon,jDimLat,kDimLev)
                  DT     ,&!REAL(KIND=r4), INTENT(IN   ) :: DT(iDimLon,jDimLat,kDimLev)
                  DP     ,&!REAL(KIND=r4), INTENT(IN   ) :: DP(iDimLon,jDimLat)
                  UP     ,&!REAL(KIND=r4), INTENT(OUT  ) :: UP(iDimLon,jDimLat,kDimLev)
                  VP     ,&!REAL(KIND=r4), INTENT(OUT  ) :: VP(iDimLon,jDimLat,kDimLev)
                  TP     ,&!REAL(KIND=r4), INTENT(OUT  ) :: TP(iDimLon,jDimLat,kDimLev)
                  PP     ,&!REAL(KIND=r4), INTENT(OUT  ) :: PP(iDimLon,jDimLat)
                  U      ,&!REAL(KIND=r4), INTENT(INOUT) :: U (iDimLon,jDimLat,kDimLev)
                  V      ,&!REAL(KIND=r4), INTENT(INOUT) :: V (iDimLon,jDimLat,kDimLev)
                  T      ,&!REAL(KIND=r4), INTENT(INOUT) :: T (iDimLon,jDimLat,kDimLev)
                  P       )!REAL(KIND=r4), INTENT(INOUT) :: P (iDimLon,jDimLat)

      CALL POLES (iDimLon   ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                  jDimLat   ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                  kDimLev   ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                  U    ,&! REAL(KIND=r4), INTENT(INOUT) :: U(iDimLon,jDimLat,kDimLev)
                  V    ,&! REAL(KIND=r4), INTENT(INOUT) :: V(iDimLon,jDimLat,kDimLev)
                  T    ,&! REAL(KIND=r4), INTENT(INOUT) :: T(iDimLon,jDimLat,kDimLev)
                  P    ,&! REAL(KIND=r4), INTENT(INOUT) :: P(iDimLon,jDimLat)
                  TPG1  )! REAL(KIND=r4), INTENT(INOUT) :: H(iDimLon,jDimLat)

! ########################################
!   5. OUTPUT
! ########################################
      IF (MOD(I,IO).EQ.0) THEN
        CALL DIAGH (iDimLon     ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat     ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                    kDimLev     ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                    DZ     ,&!REAL(KIND=r4), INTENT(IN   ) :: DZ (kDimLev)
                    TPG1   ,&!REAL(KIND=r4), INTENT(IN   ) :: TPG(iDimLon,jDimLat)
                    PTOP   ,&!REAL(KIND=r4), INTENT(IN   ) :: PTOP
                    P      ,&!REAL(KIND=r4), INTENT(IN   ) :: P  (iDimLon,jDimLat)
                    EX     ,&!REAL(KIND=r4), INTENT(IN   ) :: EX (iDimLon,jDimLat,kDimLev)
                    T      ,&!REAL(KIND=r4), INTENT(IN   ) :: T  (iDimLon,jDimLat,kDimLev)
                    H       )!REAL(KIND=r4), INTENT(OUT  ) :: H  (iDimLon,jDimLat,kDimLev)
        CALL DIAGS (iDimLon    ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat    ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                    TPG1  ,&! REAL(KIND=r4), INTENT(IN   ) :: TPG(iDimLon,jDimLat)
                    P     ,&! REAL(KIND=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
                    SLP    )! REAL(KIND=r4), INTENT(OUT  ) :: SlP(iDimLon,jDimLat)
        CALL DIAGZ (iDimLon    ,&! INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat    ,&! INTEGER      , INTENT(IN   ) :: jDimLat
                    kDimLev    ,&! INTEGER      , INTENT(IN   ) :: kDimLev
                    DZ    ,&! REAL(KIND=r4), INTENT(IN   ) :: DZ(kDimLev)
                    PTOP  ,&! REAL(KIND=r4), INTENT(IN   ) :: PTOP
                    P     ,&! REAL(KIND=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
                    H     ,&! REAL(KIND=r4), INTENT(IN   ) :: H(iDimLon,jDimLat,kDimLev)
                    Z      )! REAL(KIND=r4), INTENT(OUT  ) :: Z(iDimLon,jDimLat,kDimLev)
        WRITE(6,'(1X,I4,1X,A1,1X,I4)') I/IO,'/',IX/IO
        CALL WRDAT (iDimLon      ,&!INTEGER      , INTENT(IN   ) :: iDimLon
                    jDimLat      ,&!INTEGER      , INTENT(IN   ) :: jDimLat
                    kDimLev      ,&!INTEGER      , INTENT(IN   ) :: kDimLev
                    I/IO    ,&!INTEGER      , INTENT(IN   ) :: I
                    U       ,&!REAL(KIND=r4), INTENT(IN   ) :: U(iDimLon,jDimLat,kDimLev)
                    V       ,&!REAL(KIND=r4), INTENT(IN   ) :: V(iDimLon,jDimLat,kDimLev)
                    T       ,&!REAL(KIND=r4), INTENT(IN   ) :: T(iDimLon,jDimLat,kDimLev)
                    Z       ,&!REAL(KIND=r4), INTENT(IN   ) :: Z(iDimLon,jDimLat,kDimLev)
                    SLP     ,&!REAL(KIND=r4), INTENT(IN   ) :: P(iDimLon,jDimLat)
                    DM1     ,&!REAL(KIND=r4), INTENT(OUT  ) :: GG(iDimLon,jDimLat+1,kDimLev)
                    DM2      )!REAL(KIND=r4), INTENT(OUT  ) :: PP(iDimLon,jDimLat+1)

PRINT*,'pkubota','TIM=',TIM,'I/IO=',I/IO,'I=',I,'IO=',IO

      ENDIF
PRINT*,'TIM=',TIM,'I/IO=',I/IO,'I=',I,'IO=',IO


   11 CONTINUE

! #### Do loop 11 ends here ####

      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
      CLOSE(15)

      STOP

 END  SUBROUTINE Init
 SUBROUTINE Run()
    IMPLICIT NONE
    !CALL  Init_SGCM
 END  SUBROUTINE Run
 SUBROUTINE Finalize()
    IMPLICIT NONE
    !CALL  Init_SGCM
 END  SUBROUTINE Finalize
END PROGRAM Main

