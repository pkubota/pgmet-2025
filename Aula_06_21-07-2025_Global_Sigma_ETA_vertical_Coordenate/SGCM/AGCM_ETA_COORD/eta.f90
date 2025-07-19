MODULE MODULE_ETA_MODEL
  IMPLICIT NONE
  PRIVATE
  ! Selecting Kinds
  INTEGER, PUBLIC, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PUBLIC, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PUBLIC, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PUBLIC, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PUBLIC, PARAMETER :: r16 = SELECTED_REAL_KIND(15)! Kind for 128-bits Real Numbers


  PUBLIC :: eta,TROCA,FILTER,SWRITE,LEAP,EULER,TVENTO,CAMP,TEMPE,DIAG2,DIAG1,CONDI,FORCE,PREFE,INICIA

CONTAINS


  SUBROUTINE eta()
    !************************************************************************
    !                   PROGRAMA PRINCIPAL modelp2.f  V.2
    !
    ! ********************************************************* 

    !              MODELO TRIDIMENSIONAL BAROCLINICO DE ECUACOES 
    !          PRIMITIVAS  EM COORDENAS ESFERICAS NA HORIZONTAL E COORDE-
    !          NADAS ETA NA VERTICAL. ESTA VERSAO NAO TEM VAPOR DA AGUA
    !************************************************************************
    !                  DEFINICAO DOS PARAMETROS
    !                  ------------------------
    !                JM.......NUMR. PONTOS EM LAT
    !                IM.......NUMR. PONTOS EM LONG.
    !                KM.......NUMR. DE NIVEIS PRINCIPAIS(U,V,T)
    !                KMI......NUMR. DE NIVEIS INTERMEDIOS(PHI,EP(ETA PONTO)
    !                NT.......NUMERO DE ITERACOES NA INTEGRACAO NO TEMPO 
    !                NTMAX....MAXIMA VALOR DE NT
    !                DT.......DELTA DE TEMPO            
    !---------------------------------------------- ------------------
    !                  VARIAVEIS PROGNOSTICAS
    !                  ----------------------
    !              PIA,PIF,PIC....... P*( FUNCAO DE P. DE SUP)
    !              UA,U,UC........... VENTO ZONAL
    !              VA,VF,VC.......... VENTO MERIDIONAL
    !              TA,TF,TC.......... TEMPERATURA
    !              A,F,C............. TEMPO NT-1,NT,NT+1 RESPECT.
    !-------------------------------------------------------------
    !********************************************************************
    !             PROGRAMA PRINCIPAL version 2.0
    !**********************************************************************  
    !                 P  A   R   A   M   E   T   R   O   S
    !**********************************************************************                      

    !       DADOS PARA INICIAR
    !       ================== 
    !       FLAT1=LATITUDE INICIAL DO DOMINIO
    !       FLON1=LONGITUD INICIAL DO DOMINIO
    !       KPAD=1 CON TEMPERATURA PADRAO INICIAL PARA REPOSO
    !       KPAD=2 TEMPERATURA MEDIA ZONAL MERIDIONAL INICIAL PARA REPOSO


    REAL(KIND=r8), PARAMETER :: FLAT1=-55.
    REAL(KIND=r8), PARAMETER :: FLON1=-177.5
    INTEGER      , PARAMETER :: KPAD=2
    !----------------------------------------------------------------------
    !      DADOS PARA TOPOGRAFIA
    !      =====================
    !      KTOP=1(SI TOPOGRAFIA)
    !      KTOP=2(NO TOPOGRAFIA)
    !      L=POSICION CENTRAL DOS ANDES LONGITUD
    !      LY=POSICION CENTRAL DOS ANDES LATITUD
    !      KSIG=1(EN SIGMA)
    !      KSIG=2(EN ETA)


    INTEGER, PARAMETER :: KSIG=2
    INTEGER, PARAMETER :: LLPK=32+11
    INTEGER, PARAMETER :: LY=16
    !--------------------------------------------------------------------
    !        DADOS DE FRONTEIRA
    !        ==================
    !            FRONTERA N/S                     FRONTEIRA E/W
    !            KFRON =1(NEWMAN                  KCI=1(CICLICA)
    !            KFRON =2(ESPONJA                 KCI=2(RADIACIONAL)
    !            KFRON =3(RADiAC)       
    !   nota en esta versao radiacional nao se usa porque e global

    INTEGER, PARAMETER :: KFRON=1
    INTEGER, PARAMETER :: KCI=1
    INTEGER, PARAMETER :: NT1=1
    REAL(KIND=r8), PARAMETER :: CDF=0.
    !-----------------------------------------------------------------------
    !      DADOS FORCANTE
    !      ==============
    !      KFOR=1(SI)
    !      KFOR=2(NO)
    !      NT2=NT INICIAL PARA FORCANTE* 
    !      TMAX= T MAX.INTEGRADA ESP-TEMPO
    ! 
    !      KQ=1   6HORAS  INV. CICLICA
    !      KQ=2   ESTACIONA COM INV. 
    !      KQ=3   12HORAS SENO UN CICLO  
    !      KQ=4  
    !      FLO1........LONGITUD ONDE Q E MAX. para Amazonia
    !      FLA1........LATITUD ONDE Q E MAX. 
    !      LX,LYY......AMPLITUD EN X E Y DE Q
    !      AIN.........ANGULO DE INCLINACAO RESPECTO A X

    !ccc        PARAMETER(FLO1=-58.,FLA1=-10.,TMAX=8.5)
    REAL(KIND=r8), PARAMETER :: FLO1=-58.
    REAL(KIND=r8), PARAMETER :: FLA1=-10.
    REAL(KIND=r8), PARAMETER :: TMAX=10.
    REAL(KIND=r8), PARAMETER :: LX=10.
    REAL(KIND=r8), PARAMETER :: LYY=10.
    REAL(KIND=r8), PARAMETER :: AIN=0.
    !-------------quit-----------------------------------------------------  
    !      DADOS INICIAIS E SAIDA
    !      ===========================
    !      KINI=1(INICIAL MEDIAi ZONAL)
    !      KINI=2(INICIAL REPOSO)
    !      MUL=NT DE SAIDA
    !      NT1=NUMERO DE VECES FORW-BACK ANTES DE AVANZAR NO TEMPO
    !      CF COEF.DE FILTR0(MAX=0.5, 0.2 PARA MEDIA ZONAL O.1 PARA REPOSO)
    !234567~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !     abajo puede ver-se MUL
    !               MUL---------horas de saida
    !                144        12h          
    !                288        24h
    !                576        48h
    !                864        72h
    !                  .         .
    !                 72         6h
    !                 36         3h


    !***********************************************************
    !                     NO TEM MEDIA ZONAL      
    !***********************************************************

    !                      6HORAS INVErTIDA UM CICLO POR 48 HORAS
    !                              ----------------------------
    INTEGER, PARAMETER :: KINI=2  !      KINI=1(INICIAL MEDIAi ZONAL)
    INTEGER, PARAMETER :: KTOP=2
    INTEGER, PARAMETER :: KFOR=1
    INTEGER, PARAMETER :: NT2=1
    INTEGER, PARAMETER :: KQ=1
    INTEGER, PARAMETER :: MUL=36
    INTEGER, PARAMETER :: NTMAX=500
    !2       PARAMETER(KINI=2,KTOP=1,KFOR=1,NT2=1,KQ=4,MUL=72,NTMAX=2*288)

    !2345678                  12HORAS POR 72 HORAS
    !                                  --------------------
    !exp3          PARAMETER(KINI=2,KTOP=1,KFOR=1,NT2=1,KQ=3,MUL=72,NTMAX=5*288)
    !       PARAMETER(KINI=2,KTOP=1,KFOR=1,NT2=1,KQ=3,MUL=72,NTMAX=3*288)

    !                          CICLOS POR 6 DIAS cubica 6h
    !234567           -----------------------------
    !5ff     PARAMETER(KINI=2,KTOP=2,KFOR=1,NT2=1,KQ=1,MUL=72,NTMAX=288*5)
    !6ff      PARAMETER(KINI=2,KTOP=1,KFOR=1,NT2=1,KQ=1,MUL=72,NTMAX=288*6)

    !                            ESTACIONARIA
    !                            -------------
    !7    PARAMETER(KINI=2,KTOP=2,KFOR=1,NT2=1,KQ=2,MUL=288,NTMAX=288*15)
    !8     PARAMETER(KINI=2,KTOP=1,KFOR=1,NT2=1,KQ=2,MUL=288,NTMAX=288*6)

    !
    !********************************************************************        
    !                       TEM MEDIA ZONAL  
    !********************************************************************

    !234567                    CICLICO POR 5 DIAS(144=12HORAS)
    !                    ------------------------------
    !9ff        PARAMETER(KINI=1,KTOP=2,KFOR=1,NT2=1,KQ=1,MUL=72,NTMAX=288*6)
    !10ff      PARAMETER(KINI=1,KTOP=1,KFOR=1,NT2=1,KQ=1,MUL=72,NTMAX=288*6)


    !                  ESTACIONARIO con media zonal 
    !2345678           ----------------------------           

    !11ff      PARAMETER(KINI=1,KTOP=2,KFOR=1,NT2=1,KQ=2,MUL=288,NTMAX=288*8)
    !12ff  PARAMETER(KINI=1,KTOP=1,KFOR=1,NT2=1,KQ=2,MUL=288,NTMAX=288*8)


    !                    MEDIA ZONAL SOLO 
    !                    ------------------------
    !         PARAMETER(KINI=1,KTOP=1,KFOR=1,NT2=1,KQ=1,MUL=288,NTMAX=288*20)
    !c         PARAMETER(KINI=1,KTOP=2,KFOR=2,NT2=1,KQ=0,MUL=18,NTMAX=288*5)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    INTEGER      , PARAMETER :: IM=146
    INTEGER      , PARAMETER :: JM=44
    INTEGER      , PARAMETER :: KM=7
    INTEGER      , PARAMETER :: JM1=JM-1
    INTEGER      , PARAMETER :: IM1=IM
    INTEGER      , PARAMETER :: JM2=JM-2
    INTEGER      , PARAMETER :: IM2=IM-2
    INTEGER      , PARAMETER :: KIM=KM+1
    REAL(KIND=r8), PARAMETER :: PT=20.
    REAL(KIND=r8), PARAMETER :: RT=6.365E+6
    REAL(KIND=r8), PARAMETER :: PIS=3.141596
    REAL(KIND=r8), PARAMETER :: RS=287.
    REAL(KIND=r8), PARAMETER :: RCP=0.286
    REAL(KIND=r8), PARAMETER :: CP=1004.
    REAL(KIND=r8), PARAMETER :: DLAT=2.5
    REAL(KIND=r8), PARAMETER :: DLON=2.5
    REAL(KIND=r8), PARAMETER :: RAD=PIS/180.
    REAL(KIND=r8), PARAMETER :: N=1./RT
    REAL(KIND=r8), PARAMETER :: DX=DLON*RAD
    REAL(KIND=r8), PARAMETER :: DY=DLAT*RAD 
    REAL(KIND=r8), PARAMETER :: OMEGA=7.292E-5
    REAL(KIND=r8), PARAMETER :: GRA=9.8
    REAL(KIND=r8), PARAMETER :: QMXY=PIS/4.
    REAL(KIND=r8), PARAMETER :: QMZ=2./PIS


    REAL(KIND=r8), PARAMETER :: DT=300.
    REAL(KIND=r8) :: ALFAC(IM,JM,KM)
    REAL(KIND=r8) :: ALFA(IM,JM,KM)
    REAL(KIND=r8) :: PIF(IM,JM)
    REAL(KIND=r8) :: UF(IM1,JM1,KM)
    REAL(KIND=r8) :: VF(IM1,JM1,KM)
    REAL(KIND=r8) :: TF(IM,JM,KM)
    REAL(KIND=r8) :: PIA(IM,JM)
    REAL(KIND=r8) :: UA(IM1,JM1,KM)
    REAL(KIND=r8) :: VA(IM1,JM1,KM)
    REAL(KIND=r8) :: TA(IM,JM,KM)
    REAL(KIND=r8) :: PIZ(IM,JM)
    REAL(KIND=r8) :: UZ(IM1,JM1,KM)
    REAL(KIND=r8) :: VZ(IM1,JM1,KM)
    REAL(KIND=r8) :: TZ(IM,JM,KM)

    REAL(KIND=r8) :: PIC(IM,JM)
    REAL(KIND=r8) :: UC(IM1,JM1,KM)
    REAL(KIND=r8) :: VC(IM1,JM1,KM)
    REAL(KIND=r8) :: TC(IM,JM,KM)
    REAL(KIND=r8) :: GC(IM1,JM1,KM)
    REAL(KIND=r8) :: GA(IM1,JM1,KM)
    REAL(KIND=r8) :: GF(IM1,JM1,KM)
    REAL(KIND=r8) :: TPI(IM,JM)
    REAL(KIND=r8) :: TU(IM1,JM1,KM)
    REAL(KIND=r8) :: TV(IM1,JM1,KM)
    REAL(KIND=r8) :: TT(IM,JM,KM)
    REAL(KIND=r8) :: NS(IM,JM)
    REAL(KIND=r8) :: M(JM)
    REAL(KIND=r8) :: MI(JM)
    REAL(KIND=r8) :: FLAT(JM)
    REAL(KIND=r8) :: FLON(IM)
    INTEGER       :: KSE(IM,JM)
    INTEGER       :: KSU(IM1,JM1)
    INTEGER       :: KSV(IM1,JM1)
    REAL(KIND=r8) :: EK(KM)
    REAL(KIND=r8) :: EKI(KIM)
    REAL(KIND=r8) :: PHIS(IM,JM)
    REAL(KIND=r8) :: PRF(IM,JM)
    REAL(KIND=r8) :: PS(IM,JM)
    REAL(KIND=r8) :: EP(IM,JM,KIM)
    REAL(KIND=r8) :: WP(IM,JM,KM)
    REAL(KIND=r8) :: PHI(IM,JM,KIM)
    REAL(KIND=r8) :: PHIC(IM,JM,KIM)
    REAL(KIND=r8) :: GPX(IM1,JM1,KM)
    REAL(KIND=r8) :: GPY(IM1,JM1,KM)
    REAL(KIND=r8) :: WW(JM)
    REAL(KIND=r8) :: Q(IM,JM,KM)
    REAL(KIND=r8) :: CF
    !Auxliar Variable 
    REAL(KIND=r8) :: DIV(IM,JM,KM)
    REAL(KIND=r8) :: EP1(IM,JM,KIM)
    REAL(KIND=r8) :: GPXC(IM1,JM1,KM)
    INTEGER       :: irec1
    INTEGER       :: KFOR1
    INTEGER       :: KINI1
    INTEGER       :: KONT
    INTEGER       :: KTOP1
    INTEGER       :: NT
    REAL(KIND=r8) :: PHI1(IM,JM,KIM)
    REAL(KIND=r8) :: PHIC1(IM,JM,KIM)
    ! REAL(KIND=r8), INTENT(IN   ) :: GRA
    !2345678
    !----------------------------------------------------------
    !      Define as variaveis auxuliares lat., lon.,etas
    !----------------------------------------------------------
    UF=0.0
    VF=0.0
    TF=0.0
    WP=0.0
    EP=0.0
    EP1=0.0
    TPI=0.0
    PIF=0.0
    !CALL INICIA(FLAT1,FLON1,KPAD,KSIG) 
    CALL INICIA(&
         IM        ,&!INTEGER      , INTENT(IN   ) :: IM
         JM        ,&!INTEGER      , INTENT(IN   ) :: JM
         KM        ,&!INTEGER      , INTENT(IN   ) :: KM
         KIM       ,&!INTEGER      , INTENT(IN   ) :: KIM 
         IM1       ,&!INTEGER      , INTENT(IN   ) :: IM1
         JM1       ,&!INTEGER      , INTENT(IN   ) :: JM1
         PT        ,&!REAL(KIND=r8), INTENT(IN   ) :: PT
         RT        ,&!REAL(KIND=r8), INTENT(IN   ) :: RT      !RT=6.365E+6
         RAD       ,&!REAL(KIND=r8), INTENT(IN   ) :: RAD    !RAD=PIS/180.
         DLON      ,&!REAL(KIND=r8), INTENT(IN   ) :: DLON   !DLON=2.5
         DLAT      ,&!REAL(KIND=r8), INTENT(IN   ) :: DLAT   !DLAT=2.5
         FLAT1     ,&!REAL(KIND=r8), INTENT(IN   ) :: FLAT1
         FLON1     ,&!REAL(KIND=r8), INTENT(IN   ) :: FLON1
         KPAD      ,&!INTEGER      , INTENT(IN   ) :: KPAD
         KSIG      ,&!INTEGER      , INTENT(IN   ) :: KSIG
         FLAT      ,&!REAL(KIND=r8), INTENT(OUT  ) :: FLAT(JM)
         FLON      ,&!REAL(KIND=r8), INTENT(OUT  ) :: FLON(IM)
         M         ,&!REAL(KIND=r8), INTENT(OUT  ) :: M (JM)
         MI        ,&!REAL(KIND=r8), INTENT(OUT  ) :: MI(JM)
         EK        ,&!REAL(KIND=r8), INTENT(OUT  ) :: EK(KM)
         EKI       ,&!REAL(KIND=r8), INTENT(OUT  ) :: EKI(KIM)
         PRF       ,&!REAL(KIND=r8), INTENT(OUT  ) :: PRF(IM,JM)
         NS        ,&!REAL(KIND=r8), INTENT(OUT  ) :: NS(IM,JM)
         KSE       ,&!INTEGER      , INTENT(OUT  ) :: KSE(IM,JM)
         KSU       ,&!INTEGER      , INTENT(OUT  ) :: KSU(IM1,JM1)
         KSV       ,&!INTEGER      , INTENT(OUT  ) :: KSV(IM1,JM1)
         PHIS      ,&!REAL(KIND=r8), INTENT(OUT  ) :: PHIS(IM,JM)
         PS        ,&!REAL(KIND=r8), INTENT(OUT  ) :: PS(IM,JM)
         PIF       ,&!REAL(KIND=r8), INTENT(OUT  ) :: PIF(IM,JM)
         TF        ,&!REAL(KIND=r8), INTENT(OUT  ) :: TF(IM,JM,KM)
         WW         &!REAL(KIND=r8), INTENT(OUT  ) :: WW(JM)
         ) 
   PRINT*,'INICIA',MAXVAL(TF),MINVAL(TF)
    !----------------------------------------------------------
    !     Define as presoes de referenca da orografia para
    !     para definir os valores de etas com orografia, asimi-
    !     calcula as alturas geopotenciais sobre a orografia 
    !---------------------------------------------------------

    KTOP1=KTOP
    IF(KTOP1.EQ.1)THEN
       !CALL PREFE(KSIG,L,LY) 
       CALL PREFE(&
            KSIG      ,&!INTEGER      , INTENT(IN   ) :: KSIG
            LLPK      ,&!INTEGER      , INTENT(IN   ) :: LLPK
            LY        ,&!INTEGER      , INTENT(IN   ) :: LY
            IM        ,&!INTEGER      , INTENT(IN   ) :: IM
            JM        ,&!INTEGER      , INTENT(IN   ) :: JM
            KM        ,&!INTEGER      , INTENT(IN   ) :: KM
            KIM       ,&!INTEGER      , INTENT(IN   ) :: KIM
            IM1       ,&!INTEGER      , INTENT(IN   ) :: IM1
            JM1       ,&!INTEGER      , INTENT(IN   ) :: JM1
            PT        ,&!REAL(KIND=r8), INTENT(IN   ) :: PT
            RS        ,&!REAL(KIND=r8), INTENT(IN   ) :: RS      !RS=287.
            EK        ,&!REAL(KIND=r8), INTENT(IN   ) :: EK (KIM)
            EKI       ,&!REAL(KIND=r8), INTENT(IN   ) :: EKI(KIM)
            PRF       ,&!REAL(KIND=r8), INTENT(OUT  ) :: PRF(IM,JM)
            NS        ,&!REAL(KIND=r8), INTENT(OUT  ) :: NS(IM,JM)
            PS        ,&!REAL(KIND=r8), INTENT(OUT  ) :: PS(IM,JM)
            PIF       ,&!REAL(KIND=r8), INTENT(OUT  ) :: PIF(IM,JM)
            KSU       ,&!INTEGER      , INTENT(OUT  ) :: KSU(IM1,JM1)
            KSV       ,&!INTEGER      , INTENT(OUT  ) :: KSV(IM1,JM1)
            KSE       ,&!INTEGER      , INTENT(OUT  ) :: KSE(IM,JM)
            TF        ,&!REAL(KIND=r8), INTENT(INOUT) :: TF(IM,JM,KM)
            PHIS       &!REAL(KIND=r8), INTENT(INOUT) :: PHIS(IM,JM)
            )            !PRESIONES DE REFERENCIA


    END IF
   PRINT*,'PREFE',MAXVAL(TF),MINVAL(TF)

    !-----------------------------------------------------------
    !     Define a distribucao horizontal da forcante termica
    !-----------------------------------------------------------


    KFOR1=KFOR                 
    IF(KFOR1.EQ.1)THEN
       !CALL FORCE(L,LY,FLO1,FLA1,LX,LYY,AIN,ktop)
       CALL FORCE(&
            IM       ,&!INTEGER      , INTENT(IN   ) :: IM
            JM       ,&!INTEGER      , INTENT(IN   ) :: JM
            KM       ,&!INTEGER      , INTENT(IN   ) :: KM
            LLPK     ,&!INTEGER      , INTENT(IN   ) :: LLPK
            LX       ,&!REAL(KIND=r8), INTENT(IN   ) :: LX
            LY       ,&!INTEGER      , INTENT(IN   ) :: LY
            FLO1     ,&!REAL(KIND=r8), INTENT(IN   ) :: FLO1
            FLA1     ,&!REAL(KIND=r8), INTENT(IN   ) :: FLA1
            LYY      ,&!REAL(KIND=r8), INTENT(IN   ) :: LYY
            AIN      ,&!REAL(KIND=r8), INTENT(IN   ) :: AIN
            PIS      ,&!REAL(KIND=r8), INTENT(IN   ) :: PIS !PIS=3.141596
            RAD      ,&!REAL(KIND=r8), INTENT(IN   ) :: RAD !RAD=PIS/180.
            FLON     ,&!REAL(KIND=r8), INTENT(IN   ) :: FLON(IM)
            FLAT     ,&!REAL(KIND=r8), INTENT(IN   ) :: FLAT(JM)
            EK       ,&!REAL(KIND=r8), INTENT(IN   ) :: EK(KM)
            KSE      ,&!INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
            Q        ,&!REAL(KIND=r8), INTENT(OUT  ) :: Q(IM,JM,KM)
            ktop      &!INTEGER      , INTENT(IN   ) :: ktop
            )

    END IF
   PRINT*,'FORCE',KFOR,MAXVAL(Q),MINVAL(Q)

    !----------------------------------------------------------
    !     Define as condicoes iniciais das variaveis prognost.
    !----------------------------------------------------------

    KINI1=KINI
    IF(KINI1.EQ.1)THEN
       !CALL CONDI(KTOP)
       CALL CONDI(&
            IM      ,&! INTEGER      , INTENT(IN   ) :: IM
            JM      ,&! INTEGER      , INTENT(IN   ) :: JM
            KM      ,&! INTEGER      , INTENT(IN   ) :: KM
            KIM     ,&! INTEGER      , INTENT(IN   ) :: KIM
            IM1     ,&! INTEGER      , INTENT(IN   ) :: IM1
            JM1     ,&! INTEGER      , INTENT(IN   ) :: JM1
            DY      ,&! REAL(KIND=r8), INTENT(IN   ) :: DY
            N       ,&! REAL(KIND=r8), INTENT(IN   ) :: N    !N=1./RT
            OMEGA   ,&! REAL(KIND=r8), INTENT(IN   ) :: OMEGA!OMEGA=7.29E-5
            PT      ,&! REAL(KIND=r8), INTENT(IN   ) :: PT
            RAD     ,&! REAL(KIND=r8), INTENT(IN   ) :: RAD   !RAD=PIS/180
            RS      ,&! REAL(KIND=r8), INTENT(IN   ) :: RS    !    RS=287.
            EK      ,&! REAL(KIND=r8), INTENT(IN   ) :: EK(KM)
            EKI     ,&! REAL(KIND=r8), INTENT(IN   ) :: EKI(KM)
            FLAT    ,&! REAL(KIND=r8), INTENT(IN   ) :: FLAT(JM)
            NS      ,&! REAL(KIND=r8), INTENT(IN   ) :: NS(IM,JM)
            PRF     ,&! REAL(KIND=r8), INTENT(IN   ) :: PRF(IM,JM)
            PHI     ,&! REAL(KIND=r8), INTENT(INOUT) :: PHI(IM,JM,KIM)
            PS      ,&! REAL(KIND=r8), INTENT(INOUT) :: PS(IM,JM)
            PIF     ,&! REAL(KIND=r8), INTENT(INOUT) :: PIF(IM,JM)
            PHI1    ,&! REAL(KIND=r8), INTENT(INOUT) :: PHI1(IM,JM,KIM)
            TF      ,&! REAL(KIND=r8), INTENT(INOUT) :: TF(IM,JM,KM)
            UF      ,&! REAL(KIND=r8), INTENT(INOUT) :: UF(IM1,JM1,KM)
            KSE      &! INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
            ) 

       !         open(8,file='inicia.dat',status='old')

       !         DO 400 K=1,KM
       ! c        DO 400 J=1,JM
       !         DO 400 I=1,IM
       !         read(8,98)PS(I,J),TF(I,J,K),UF(I,J,K)
       !98        FORMAT(1X,F12.4,F12.4,F12.4)
       !400       CONTINUE
       !              close(8)
       WRITE(*,*)'pasei1'
    ENDIF
   PRINT*,'CONDI',KINI,MAXVAL(TF),MINVAL(TF)

    !----------------------------------------------------------
    !        inicializa irec1 para gravar os dados
    !        --------------------------------------  

    !----------------------------------------------------------
    !               inicio da integracao
    !               --------------------  
    irec1=0

    DO NT=1,NTMAX

       PRINT*,'NT=',NT,'NTMAX=',NTMAX
       !c                      KONT=144   !com media zonal  y topografia
       KONT=36    !com media zonal  

       !                     KONT=1     ! no tem media zonal  y  no tienen topo

       IF(NT.LT.KONT)THEN

          CF=0.1           

       ELSE

          !                     CF=0.05
          CF=0.10

       ENDIF

       !********************************************************
       !     Calcula velocidade vertical eta e omega,asimismo PIC
       !---------------------------------------------------------

       !CALL  DIAG1(NT,KFRON,KCI,DT,NT1)
       CALL DIAG1(&
            NT              ,&!INTEGER      , INTENT(IN   ) :: NT
            KFRON           ,&!INTEGER      , INTENT(IN   ) :: KFRON 
            KCI             ,&!INTEGER      , INTENT(IN   ) :: KCI   
            DT              ,&!REAL(KIND=r8), INTENT(IN   ) :: DT    
            NT1             ,&!INTEGER      , INTENT(IN   ) :: NT1   
            IM              ,&!INTEGER      , INTENT(IN   ) :: IM
            JM              ,&!INTEGER      , INTENT(IN   ) :: JM
            KM              ,&!INTEGER      , INTENT(IN   ) :: KM
            KIM             ,&!INTEGER      , INTENT(IN   ) :: KIM
            IM1             ,&!INTEGER      , INTENT(IN   ) :: IM1
            JM1             ,&!INTEGER      , INTENT(IN   ) :: JM1
            DX              ,&!REAL(KIND=r8), INTENT(IN   ) :: DX
            DY              ,&!REAL(KIND=r8), INTENT(IN   ) :: DY
            PT              ,&!REAL(KIND=r8), INTENT(IN   ) :: PT
            RT              ,&!REAL(KIND=r8), INTENT(IN   ) :: RT  !RT=6.365E+6,
            N               ,&!REAL(KIND=r8), INTENT(IN   ) :: N   !N=1./RT
            WW              ,&!REAL(KIND=r8), INTENT(IN   ) :: WW(JM)
            M               ,&!REAL(KIND=r8), INTENT(IN   ) :: M(JM)
            MI              ,&!REAL(KIND=r8), INTENT(IN   ) :: MI(JM)
            EK              ,&!REAL(KIND=r8), INTENT(IN   ) :: EK(KIM)
            EKI             ,&!REAL(KIND=r8), INTENT(IN   ) :: EKI(KIM)
            PIA             ,&!REAL(KIND=r8), INTENT(IN   ) :: PIA(IM,JM)
            PIZ             ,&!REAL(KIND=r8), INTENT(IN   ) :: PIZ(IM,JM)
            PS              ,&!REAL(KIND=r8), INTENT(OUT  ) :: PS(IM,JM)
            DIV             ,&!REAL(KIND=r8), INTENT(OUT  ) :: DIV(IM,JM,KM)
            WP              ,&!REAL(KIND=r8), INTENT(OUT  ) :: WP(IM,JM,KM)
            NS              ,&!REAL(KIND=r8), INTENT(IN   ) :: NS(IM,JM)
            PIF             ,&!REAL(KIND=r8), INTENT(IN   ) :: PIF(IM,JM)
            UF              ,&!REAL(KIND=r8), INTENT(IN   ) :: UF(IM1,JM1,KM)
            VF              ,&!REAL(KIND=r8), INTENT(IN   ) :: VF(IM1,JM1,KM)
            TPI             ,&!REAL(KIND=r8), INTENT(INOUT) :: TPI(IM,JM)
            EP              ,&!REAL(KIND=r8), INTENT(INOUT) :: EP (IM,JM,KIM)
            EP1             ,&!REAL(KIND=r8), INTENT(INOUT) :: EP1(IM,JM,KIM)
            PIC             ,&!REAL(KIND=r8), INTENT(INOUT) :: PIC(IM,JM)
            KSE              &!INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
            )

   PRINT*,'DIAG1',NT,MAXVAL(EP),MINVAL(EP)

       !---------------------------------------------------------
       !     Calcula  alturas geopotenciais, gradientes de presao
       !     en x e y, asimismo calcula volumen especifico(alfa)
       !----------------------------------------------------------

       !CALL DIAG2(NT,KSIG) 
       CALL DIAG2(&
            NT      ,&! INTEGER      , INTENT(IN   ) :: NT
            IM      ,&! INTEGER      , INTENT(IN   ) :: IM
            JM      ,&! INTEGER      , INTENT(IN   ) :: JM
            KM      ,&! INTEGER      , INTENT(IN   ) :: KM
            KIM     ,&! INTEGER      , INTENT(IN   ) :: KIM
            IM1     ,&! INTEGER      , INTENT(IN   ) :: IM1
            JM1     ,&! INTEGER      , INTENT(IN   ) :: JM1
            DX      ,&! REAL(KIND=r8), INTENT(IN   ) :: DX
            DY      ,&! REAL(KIND=r8), INTENT(IN   ) :: DY
            PT      ,&! REAL(KIND=r8), INTENT(IN   ) :: PT
            RS      ,&! REAL(KIND=r8), INTENT(IN   ) :: RS        !RS=287.
            RT      ,&! REAL(KIND=r8), INTENT(IN   ) :: RT        !RT=6.365E+6
            NS      ,&! REAL(KIND=r8), INTENT(IN   ) :: NS(IM,JM)
            PS      ,&! REAL(KIND=r8), INTENT(IN   ) :: PS(IM,JM)
            KSE     ,&! INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
            EK      ,&! REAL(KIND=r8), INTENT(IN   ) :: EK(KM)
            EKI     ,&! REAL(KIND=r8), INTENT(IN   ) :: EKI(KIM)
            TF      ,&! REAL(KIND=r8), INTENT(INOUT) :: TF(IM,JM,KM)
            PIF     ,&! REAL(KIND=r8), INTENT(IN   ) :: PIF(IM,JM)
            PHIS    ,&! REAL(KIND=r8), INTENT(IN   ) :: PHIS(IM,JM)
            ALFA    ,&! REAL(KIND=r8), INTENT(OUT  ) :: ALFA(IM,JM,KM)
            PHI     ,&! REAL(KIND=r8), INTENT(OUT  ) :: PHI(IM,JM,KIM)
            GPX     ,&! REAL(KIND=r8), INTENT(OUT  ) :: GPX(IM1,JM1,KM)
            GPY     ,&! REAL(KIND=r8), INTENT(OUT  ) :: GPY(IM1,JM1,KM)
            GF       &! REAL(KIND=r8), INTENT(OUT  ) :: GF(IM1,JM1,KM)
            )

   PRINT*,'DIAG2',NT,MAXVAL(GF),MINVAL(GF)

       !---------------------------------------------------------
       !     Calcula a tendencia da temperatura e TC
       !---------------------------------------------------------
       !2345678

       !CALL TEMPE(NT,KINI,KFRON,KTOP,KCI,DT,NT2,TMAX,KQ,NT1,CDF, L,LY)
       CALL TEMPE( &
            NT      ,&!INTEGER      , INTENT(IN ) :: NT
            IM      ,&!INTEGER      , INTENT(IN ) :: IM
            JM      ,&!INTEGER      , INTENT(IN ) :: JM
            KM      ,&!INTEGER      , INTENT(IN ) :: KM
            KIM     ,&!INTEGER      , INTENT(IN ) :: KIM
            IM1     ,&!INTEGER      , INTENT(IN ) :: IM1
            JM1     ,&!INTEGER      , INTENT(IN ) :: JM1
            KCI     ,&!INTEGER      , INTENT(IN ) :: KCI
            KFRON   ,&!INTEGER      , INTENT(IN ) :: KFRON
            KTOP    ,&!INTEGER      , INTENT(IN ) :: KTOP
            KQ      ,&!INTEGER      , INTENT(IN ) :: KQ
            DT      ,&!REAL(KIND=r8), INTENT(IN ) :: DT
            DX      ,&!REAL(KIND=r8), INTENT(IN ) :: DX
            DY      ,&!REAL(KIND=r8), INTENT(IN ) :: DY
            CP      ,&!REAL(KIND=r8), INTENT(IN ) :: CP
            TMAX    ,&!REAL(KIND=r8), INTENT(IN ) :: TMAX
            QMXY    ,&!REAL(KIND=r8), INTENT(IN ) :: QMXY  !QMXY=PIS/4
            QMZ     ,&!REAL(KIND=r8), INTENT(IN ) :: QMZ   !QMZ=2./PIS
            N       ,&!REAL(KIND=r8), INTENT(IN ) :: N     !N=1./RT
            RAD     ,&!REAL(KIND=r8), INTENT(IN ) :: RAD   !RAD=PIS/180.
            RT      ,&!REAL(KIND=r8), INTENT(IN ) :: RT    !RT=6.365E+6
            NT1     ,&!INTEGER      , INTENT(IN ) :: NT1
            NT2     ,&!INTEGER      , INTENT(IN ) :: NT2
            WW      ,&!REAL(KIND=r8), INTENT(IN ) :: WW(JM)
            FLAT    ,&!REAL(KIND=r8), INTENT(IN ) :: FLAT(JM)
            M       ,&!REAL(KIND=r8), INTENT(IN ) :: M(JM)
            MI      ,&!REAL(KIND=r8), INTENT(IN ) :: MI(JM)
            EKI     ,&!REAL(KIND=r8), INTENT(IN ) :: EKI(KIM)
            KSE     ,&!INTEGER      , INTENT(IN ) :: KSE(IM,JM)
            UF      ,&!REAL(KIND=r8), INTENT(IN ) :: UF(IM1,JM1,KM)
            VF      ,&!REAL(KIND=r8), INTENT(IN ) :: VF(IM1,JM1,KM)
            TF      ,&!REAL(KIND=r8), INTENT(IN ) :: TF(IM,JM,KM)
            Q       ,&!REAL(KIND=r8), INTENT(IN ) :: Q(IM,JM,KM)
            PIF     ,&!REAL(KIND=r8), INTENT(IN ) :: PIF(IM,JM)
            EP      ,&!REAL(KIND=r8), INTENT(IN ) :: EP(IM,JM,KIM)
            WP      ,&!REAL(KIND=r8), INTENT(IN ) :: WP(IM,JM,KM)
            ALFA    ,&!REAL(KIND=r8), INTENT(IN ) :: ALFA(IM,JM,KM) 
            PIC     ,&!REAL(KIND=r8), INTENT(IN ) :: PIC(IM,JM)
            PIA     ,&!REAL(KIND=r8), INTENT(IN ) :: PIA(IM,JM)
            TA      ,&!REAL(KIND=r8), INTENT(IN ) :: TA(IM,JM,KM)
            TZ      ,&!REAL(KIND=r8), INTENT(IN ) :: TZ(IM,JM,KM)
            TC      ,&!REAL(KIND=r8), INTENT(INOUT) :: TC(IM,JM,KM)
            TT       &!REAL(KIND=r8), INTENT(INOUT) :: TT(IM,JM,KM)  
            )

       !---------------------------------------------------------
       !     Calcula as medias ponderadas do gradiente de presao 
       !      e geopotencial seguiendo o esquema de Campana  
       !--------------------------------------------------------  


       !CALL  CAMP(NT) 
       CALL CAMP(&
            NT        ,&!INTEGER     , INTENT(IN   ) :: NT
            IM        ,&!INTEGER     , INTENT(IN   ) :: IM
            JM        ,&!INTEGER     , INTENT(IN   ) :: JM
            KM        ,&!INTEGER     , INTENT(IN   ) :: KM
            KIM       ,&!INTEGER     , INTENT(IN   ) :: KIM
            IM1       ,&!INTEGER     , INTENT(IN   ) :: IM1
            JM1       ,&!INTEGER     , INTENT(IN   ) :: JM1
            PT        ,&!REAL(KIND=r8), INTENT(IN   ) :: PT
            RS        ,&!REAL(KIND=r8), INTENT(IN   ) :: RS
            RT        ,&!REAL(KIND=r8), INTENT(IN   ) :: RT
            DX        ,&!REAL(KIND=r8), INTENT(IN   ) :: DX
            PS        ,&!REAL(KIND=r8), INTENT(OUT  ) :: PS(IM,JM)
            NS        ,&!REAL(KIND=r8), INTENT(IN   ) :: NS(IM,JM)
            PIC       ,&!REAL(KIND=r8), INTENT(IN   ) :: PIC(IM,JM)
            GPXC      ,&!REAL(KIND=r8), INTENT(OUT  ) :: GPXC(IM1,JM1,KM)
            EK        ,&!REAL(KIND=r8), INTENT(IN   ) :: EK(KM)
            EKI       ,&!REAL(KIND=r8), INTENT(IN   ) :: EKI(KIM)
            PHIS      ,&!REAL(KIND=r8), INTENT(IN   ) :: PHIS(IM,JM)
            PHIC1     ,&!REAL(KIND=r8), INTENT(INOUT) :: PHIC1(IM,JM,KIM)
            PHIC      ,&!REAL(KIND=r8), INTENT(INOUT) :: PHIC(IM,JM,KIM)
            TC        ,&!REAL(KIND=r8), INTENT(IN   ) :: TC(IM,JM,KM)
            ALFAC     ,&!REAL(KIND=r8), INTENT(INOUT) :: ALFAC(IM,JM,KM)
            GA        ,&!REAL(KIND=r8), INTENT(IN   ) :: GA(IM1,JM1,KM)
            GC        ,&!REAL(KIND=r8), INTENT(INOUT) :: GC(IM1,JM1,KM)
            GF        ,&!REAL(KIND=r8), INTENT(INOUT) :: GF(IM1,JM1,KM)
            KSE        &!INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
            )
       !--------------------------------------------------------
       !     Calcula as tendencia de vento zonal e meridional
       !--------------------------------------------------------


       !CALL TVENTO(CDF)
       CALL TVENTO(&
            IM       ,&! INTEGER      , INTENT(IN   ) :: IM
            JM       ,&! INTEGER      , INTENT(IN   ) :: JM
            KM       ,&! INTEGER      , INTENT(IN   ) :: KM
            KIM      ,&! INTEGER      , INTENT(IN   ) :: KIM
            IM1      ,&! INTEGER      , INTENT(IN   ) :: IM1
            JM1      ,&! INTEGER      , INTENT(IN   ) :: JM1
            DX       ,&! REAL(KIND=r8), INTENT(IN   ) :: DX
            DY       ,&! REAL(KIND=r8), INTENT(IN   ) :: DY
            GRA      ,&! REAL(KIND=r8), INTENT(IN   ) :: GRA
            OMEGA    ,&! REAL(KIND=r8), INTENT(IN   ) :: OMEGA
            RAD      ,&! REAL(KIND=r8), INTENT(IN   ) :: RAD
            RT       ,&! REAL(KIND=r8), INTENT(IN   ) :: RT
            N        ,&! REAL(KIND=r8), INTENT(IN   ) :: N
            KSE      ,&! INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
            UF       ,&! REAL(KIND=r8), INTENT(IN   ) :: UF(IM1,JM1,KM)
            VF       ,&! REAL(KIND=r8), INTENT(IN   ) :: VF(IM1,JM1,KM)
            PIF      ,&! REAL(KIND=r8), INTENT(IN   ) :: PIF(IM,JM)
            FLAT     ,&! REAL(KIND=r8), INTENT(IN   ) :: FLAT(JM)
            M        ,&! REAL(KIND=r8), INTENT(IN   ) :: M(JM)
            MI       ,&! REAL(KIND=r8), INTENT(IN   ) :: MI(JM)
            EK       ,&! REAL(KIND=r8), INTENT(IN   ) :: EK(KM)
            EKI      ,&! REAL(KIND=r8), INTENT(IN   ) :: EKI(KIM)
            EP       ,&! REAL(KIND=r8), INTENT(IN   ) :: EP(IM,JM,KIM)
            PHIS     ,&! REAL(KIND=r8), INTENT(IN   ) :: PHIS(IM,JM)
            ALFA     ,&! REAL(KIND=r8), INTENT(IN   ) :: ALFA(IM,JM,KM)
            GPY      ,&! REAL(KIND=r8), INTENT(IN   ) :: GPY(IM1,JM1,KM) 
            PHI      ,&! REAL(KIND=r8), INTENT(IN   ) :: PHI(IM,JM,KIM)
            TU       ,&! REAL(KIND=r8), INTENT(OUT  ) :: TU(IM1,JM1,KM)
            TV       ,&! REAL(KIND=r8), INTENT(OUT  ) :: TV(IM1,JM1,KM)
            GF        &! REAL(KIND=r8), INTENT(IN   ) :: GF(IM1,JM1,KM)
            )

       !-------------------------------------------------------
       !     Para o primer step usa-se esquema de euler para
       !     integr.no tempo do vento e para os outros steps 
       !     usa-se leapfrog
       !------------------------------------------------------

       IF(NT.LE.NT1)THEN
          ! CALL EULER(DT,NT)
          CALL EULER(&
               IM       , &!INTEGER      , INTENT(IN   ) :: IM
               JM       , &!INTEGER      , INTENT(IN   ) :: JM
               KM       , &!INTEGER      , INTENT(IN   ) :: KM
               IM1      , &!INTEGER      , INTENT(IN   ) :: IM1
               JM1      , &!INTEGER      , INTENT(IN   ) :: JM1
               NT       , &!INTEGER      , INTENT(IN   ) :: NT
               DT       , &!REAL(KIND=r8), INTENT(IN   ) :: DT
               KSU      , &!INTEGER      , INTENT(IN   ) :: KSU(IM1,JM1)
               KSV      , &!INTEGER      , INTENT(IN   ) :: KSV(IM1,JM1)
               UC       , &!REAL(KIND=r8), INTENT(INOUT) :: UC(IM1,JM1,KM)
               VC       , &!REAL(KIND=r8), INTENT(INOUT) :: VC(IM1,JM1,KM)
               UF       , &!REAL(KIND=r8), INTENT(IN   ) :: UF(IM1,JM1,KM)
               VF       , &!REAL(KIND=r8), INTENT(IN   ) :: VF(IM1,JM1,KM)
               TU       , &!REAL(KIND=r8), INTENT(IN   ) :: TU(IM1,JM1,KM)
               TV       , &!REAL(KIND=r8), INTENT(IN   ) :: TV(IM1,JM1,KM)
               PIC      , &!REAL(KIND=r8), INTENT(IN   ) :: PIC(IM,JM)
               PIF        )!REAL(KIND=r8), INTENT(IN   ) :: PIF(IM,JM)
       ELSE
          !CALL LEAP(KFRON,KCI,DT)
          CALL LEAP( &
               IM          ,&!INTEGER      , INTENT(IN   ) :: IM
               JM          ,&!INTEGER      , INTENT(IN   ) :: JM
               KM          ,&!INTEGER      , INTENT(IN   ) :: KM
               IM1         ,&!INTEGER      , INTENT(IN   ) :: IM1
               JM1         ,&!INTEGER      , INTENT(IN   ) :: JM1
               KFRON       ,&!INTEGER      , INTENT(IN   ) :: KFRON
               KCI         ,&!INTEGER      , INTENT(IN   ) :: KCI
               DT          ,&!REAL(KIND=r8), INTENT(IN   ) :: DT
               KSU         ,&!INTEGER      , INTENT(IN   ) :: KSU(IM1,JM1)
               KSV         ,&!INTEGER      , INTENT(IN   ) :: KSV(IM1,JM1)
               PIC         ,&!REAL(KIND=r8), INTENT(IN   ) :: PIC(IM,JM)
               PIA         ,&!REAL(KIND=r8), INTENT(IN   ) :: PIA(IM,JM)
               UA          ,&!REAL(KIND=r8), INTENT(IN   ) :: UA(IM1,JM1,KM)
               VA          ,&!REAL(KIND=r8), INTENT(IN   ) :: VA(IM1,JM1,KM)
               UC          ,&!REAL(KIND=r8), INTENT(INOUT) :: UC(IM1,JM1,KM)
               VC          ,&!REAL(KIND=r8), INTENT(INOUT) :: VC(IM1,JM1,KM)
               UF          ,&!REAL(KIND=r8), INTENT(IN   ) :: UF(IM1,JM1,KM)
               VF          ,&!REAL(KIND=r8), INTENT(IN   ) :: VF(IM1,JM1,KM)
               UZ          ,&!REAL(KIND=r8), INTENT(IN   ) :: UZ(IM1,JM1,KM)
               VZ          ,&!REAL(KIND=r8), INTENT(IN   ) :: VZ(IM1,JM1,KM)
               TU          ,&!REAL(KIND=r8), INTENT(IN   ) :: TU(IM1,JM1,KM)
               TV          ,&!REAL(KIND=r8), INTENT(IN   ) :: TV(IM1,JM1,KM)
               WW           &!REAL(KIND=r8), INTENT(IN   ) :: WW(JM)
               )

       ENDIF

       !-------------------------------------------------------
       !     saida dos resultados para cada mult de NT
       !-------------------------------------------------------

       IF(MOD(NT,MUL).EQ.0.)THEN 
          !CALL SWRITE(IREC1,L,LY,KTOP)
          CALL SWRITE(&
               IREC1    ,&!INTEGER      , INTENT(INOUT) :: IREC1
               IM       ,&!INTEGER      , INTENT(IN   ) :: IM
               JM       ,&!INTEGER      , INTENT(IN   ) :: JM
               KM       ,&!INTEGER      , INTENT(IN   ) :: KM
               IM1      ,&!INTEGER      , INTENT(IN   ) :: IM1
               JM1      ,&!INTEGER      , INTENT(IN   ) :: JM1
               IM2      ,&!INTEGER      , INTENT(IN   ) :: IM2
               JM2      ,&!INTEGER      , INTENT(IN   ) :: JM2
               KIM      ,&!INTEGER      , INTENT(IN   ) :: KIM
               PT       ,&!REAL(KIND=r8), INTENT(IN   ) :: PT
               KSE      ,&!INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
               KSU      ,&!INTEGER      , INTENT(IN   ) :: KSU(IM1,JM1)
               KSV      ,&!INTEGER      , INTENT(IN   ) :: KSV(IM1,JM1)
               PIC      ,&!REAL(KIND=r8), INTENT(IN   ) :: PIC(IM,JM)
               EKI      ,&!REAL(KIND=r8), INTENT(IN   ) :: EKI(KIM)
               PHIC     ,&!REAL(KIND=r8), INTENT(IN   ) :: PHIC(IM,JM,KIM)
               WP       ,&!REAL(KIND=r8), INTENT(IN   ) :: WP(IM,JM,KM)
               UC       ,&!REAL(KIND=r8), INTENT(IN   ) :: UC(IM1,JM1,KM)
               VC       )!REAL(KIND=r8), INTENT(IN   ) :: VC(IM1,JM1,KM)
          WRITE(*,*)'pasie swrite',NT/MUL
       ENDIF

       !------------------------------------------------------
       !              para gravar las condiciones iniciaes
       !-----------------------------------------------------
       !               if(nt.eq.10)then
       !         open(9,file='inicia.dat',status='new')
       !         
       !         DO 400 K=1,KM
       !         DO 400 J=1,JM
       !         DO 400 I=1,IM
       !         write(9,99)PS(I,J),TF(I,J,K),UF(I,J,K)
       !99        FORMAT(1X,F12.4,F12.4,F12.4)
       !400       CONTINUE
       !           endif                    
       ! 
       !-------------------------------------------------------
       !    Se faz filtro no espacio e no tempo
       !------------------------------------------------------

       IF(NT.NE.1)THEN
          !CALL FILTER(KTOP,KFRON,NT,CF,L,LY)
          CALL FILTER(&
               KTOP      ,&! INTEGER, INTENT(IN   ) ::  KTOP
               CF        ,&! INTEGER, INTENT(IN   ) ::  CF  
               IM        ,&! INTEGER, INTENT(IN   ) ::  IM
               JM        ,&! INTEGER, INTENT(IN   ) ::  JM
               KM        ,&! INTEGER, INTENT(IN   ) ::  KM
               IM1       ,&! INTEGER, INTENT(IN   ) ::  IM1
               JM1       ,&! INTEGER, INTENT(IN   ) ::  JM1
               KSE       ,&! INTEGER      , INTENT(IN   ) :: KSE(IM,JM)  !COMMON/BLOCOS/KSE(IM,JM)
               KSU       ,&! INTEGER      , INTENT(IN   ) :: KSU(IM1,JM1)
               KSV       ,&! INTEGER      , INTENT(IN   ) :: KSV(IM1,JM1)
               TA        ,&! REAL(KIND=r8), INTENT(IN   ) :: TA(IM,JM,KM)
               TC        ,&! REAL(KIND=r8), INTENT(INOUT) :: TC(IM,JM,KM)
               TF        ,&! REAL(KIND=r8), INTENT(INOUT) :: TF(IM,JM,KM)
               UA        ,&! REAL(KIND=r8), INTENT(IN   ) :: UA(IM1,JM1,KM)
               VA        ,&! REAL(KIND=r8), INTENT(IN   ) :: VA(IM1,JM1,KM)
               UC        ,&! REAL(KIND=r8), INTENT(INOUT) :: UC(IM1,JM1,KM)
               VC        ,&! REAL(KIND=r8), INTENT(INOUT) :: VC(IM1,JM1,KM)
               UF        ,&! REAL(KIND=r8), INTENT(INOUT) :: UF(IM1,JM1,KM)
               VF        ,&! REAL(KIND=r8), INTENT(INOUT) :: VF(IM1,JM1,KM)
               PIF       ,&! REAL(KIND=r8), INTENT(INOUT) :: PIF(IM,JM)
               PIC       ,&! REAL(KIND=r8), INTENT(INOUT) :: PIC(IM,JM)
               PIA        )! REAL(KIND=r8), INTENT(IN   ) :: PIA(IM,JM)
       END IF

       !-------------------------------------------------------
       !	      Mudanza de variaveis    
       !-------------------------------------------------------

       !CALL TROCA
       CALL TROCA(&
            IM         ,&!!& INTEGER      , INTENT(IN   ) :: IM
            JM         ,&!!& INTEGER      , INTENT(IN   ) :: JM
            KM         ,&!!& INTEGER      , INTENT(IN   ) :: KM
            IM1        ,&!!& INTEGER      , INTENT(IN   ) :: IM1
            JM1        ,&!!& INTEGER      , INTENT(IN   ) :: JM1
            PIF        ,&!!& REAL(KIND=r8), INTENT(INOUT) :: PIF (IM,JM)        !       COMMON/F/
            TA         ,&!!& REAL(KIND=r8), INTENT(INOUT) :: TA  (IM,JM,KM)
            TC         ,&!!& REAL(KIND=r8), INTENT(IN   ) :: TC  (IM,JM,KM)
            UC         ,&!!& REAL(KIND=r8), INTENT(IN   ) :: UC  (IM1,JM1,KM)
            VC         ,&!!& REAL(KIND=r8), INTENT(IN   ) :: VC  (IM1,JM1,KM)
            GC         ,&!!& REAL(KIND=r8), INTENT(IN   ) :: GC  (IM1,JM1,KM)
            UF         ,&!!& REAL(KIND=r8), INTENT(INOUT) :: UF  (IM1,JM1,KM) ! COMMON/F
            VF         ,&!!& REAL(KIND=r8), INTENT(INOUT) :: VF  (IM1,JM1,KM) ! COMMON/F
            GF         ,&!!& REAL(KIND=r8), INTENT(INOUT) :: GF  (IM1,JM1,KM)
            TF         ,&!!& REAL(KIND=r8), INTENT(INOUT) :: TF  (IM ,JM ,KM) ! COMMON/F/
            UZ         ,&!!& REAL(KIND=r8), INTENT(OUT  ) :: UZ  (IM1,JM1,KM)
            VZ         ,&!!& REAL(KIND=r8), INTENT(OUT  ) :: VZ  (IM1,JM1,KM)
            PIA        ,&!!& REAL(KIND=r8), INTENT(INOUT) :: PIA (IM,JM)        !       COMMON/A/
            UA         ,&!!& REAL(KIND=r8), INTENT(INOUT) :: UA  (IM1,JM1,KM) ! COMMON/A/
            VA         ,&!!& REAL(KIND=r8), INTENT(INOUT) :: VA  (IM1,JM1,KM) ! COMMON/A/
            GA         ,&!!& REAL(KIND=r8), INTENT(OUT  ) :: GA  (IM1,JM1,KM)
            PIZ        ,&!!& REAL(KIND=r8), INTENT(OUT  ) :: PIZ (IM,JM)        !       COMMON/Z/
            TZ         ,&!!& REAL(KIND=r8), INTENT(OUT  ) :: TZ  (IM,JM,KM) ! COMMON/Z/
            PIC         )!  REAL(KIND=r8), INTENT(INOUT) :: PIC (IM,JM)        !       COMMON/C/

       !------------------------------------------------------

       !-----------------------------------------------------
    END DO ! END DO NT=1,NTMAX
  END SUBROUTINE eta
  !*******************************************************************
  SUBROUTINE INICIA(&
       IM        ,&!INTEGER      , INTENT(IN   ) :: IM
       JM        ,&!INTEGER      , INTENT(IN   ) :: JM
       KM        ,&!INTEGER      , INTENT(IN   ) :: KM
       KIM       ,&!INTEGER      , INTENT(IN   ) :: KIM 
       IM1       ,&!INTEGER      , INTENT(IN   ) :: IM1
       JM1       ,&!INTEGER      , INTENT(IN   ) :: JM1
       PT        ,&!REAL(KIND=r8), INTENT(IN   ) :: PT
       RT        ,&!REAL(KIND=r8), INTENT(IN   ) :: RT      !RT=6.365E+6
       RAD       ,&!REAL(KIND=r8), INTENT(IN   ) :: RAD    !RAD=PIS/180.
       DLON      ,&!REAL(KIND=r8), INTENT(IN   ) :: DLON   !DLON=2.5
       DLAT      ,&!REAL(KIND=r8), INTENT(IN   ) :: DLAT   !DLAT=2.5
       FLAT1     ,&!REAL(KIND=r8), INTENT(IN   ) :: FLAT1
       FLON1     ,&!REAL(KIND=r8), INTENT(IN   ) :: FLON1
       KPAD      ,&!INTEGER      , INTENT(IN   ) :: KPAD
       KSIG      ,&!INTEGER      , INTENT(IN   ) :: KSIG
       FLAT      ,&!REAL(KIND=r8), INTENT(OUT  ) :: FLAT(JM)
       FLON      ,&!REAL(KIND=r8), INTENT(OUT  ) :: FLON(IM)
       M         ,&!REAL(KIND=r8), INTENT(OUT  ) :: M (JM)
       MI        ,&!REAL(KIND=r8), INTENT(OUT  ) :: MI(JM)
       EK        ,&!REAL(KIND=r8), INTENT(OUT  ) :: EK(KM)
       EKI       ,&!REAL(KIND=r8), INTENT(OUT  ) :: EKI(KIM)
       PRF       ,&!REAL(KIND=r8), INTENT(OUT  ) :: PRF(IM,JM)
       NS        ,&!REAL(KIND=r8), INTENT(OUT  ) :: NS(IM,JM)
       KSE       ,&!INTEGER      , INTENT(OUT  ) :: KSE(IM,JM)
       KSU       ,&!INTEGER      , INTENT(OUT  ) :: KSU(IM1,JM1)
       KSV       ,&!INTEGER      , INTENT(OUT  ) :: KSV(IM1,JM1)
       PHIS      ,&!REAL(KIND=r8), INTENT(OUT  ) :: PHIS(IM,JM)
       PS        ,&!REAL(KIND=r8), INTENT(OUT  ) :: PS(IM,JM)
       PIF       ,&!REAL(KIND=r8), INTENT(OUT  ) :: PIF(IM,JM)
       TF        ,&!REAL(KIND=r8), INTENT(OUT  ) :: TF(IM,JM,KM)
       WW         &!REAL(KIND=r8), INTENT(OUT  ) :: WW(JM)
       ) 
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: IM
    INTEGER      , INTENT(IN   ) :: JM
    INTEGER      , INTENT(IN   ) :: KM
    INTEGER      , INTENT(IN   ) :: KIM 
    INTEGER      , INTENT(IN   ) :: IM1
    INTEGER      , INTENT(IN   ) :: JM1
    REAL(KIND=r8), INTENT(IN   ) :: PT
    REAL(KIND=r8), INTENT(IN   ) :: RT     !RT=6.365E+6
    REAL(KIND=r8), INTENT(IN   ) :: RAD    !RAD=PIS/180.
    REAL(KIND=r8), INTENT(IN   ) :: DLON   !DLON=2.5
    REAL(KIND=r8), INTENT(IN   ) :: DLAT   !DLAT=2.5
    REAL(KIND=r8), INTENT(IN   ) :: FLAT1
    REAL(KIND=r8), INTENT(IN   ) :: FLON1
    INTEGER      , INTENT(IN   ) :: KPAD
    INTEGER      , INTENT(IN   ) :: KSIG
    REAL(KIND=r8), INTENT(OUT  ) :: FLAT(JM)
    REAL(KIND=r8), INTENT(OUT  ) :: FLON(IM)
    REAL(KIND=r8), INTENT(OUT  ) :: M (JM)
    REAL(KIND=r8), INTENT(OUT  ) :: MI(JM)
    REAL(KIND=r8), INTENT(OUT  ) :: EK(KM)
    REAL(KIND=r8), INTENT(OUT  ) :: EKI(KIM)
    REAL(KIND=r8), INTENT(OUT  ) :: PRF(IM,JM)
    REAL(KIND=r8), INTENT(OUT  ) :: NS(IM,JM)
    INTEGER      , INTENT(OUT  ) :: KSE(IM,JM)
    INTEGER      , INTENT(OUT  ) :: KSU(IM1,JM1)
    INTEGER      , INTENT(OUT  ) :: KSV(IM1,JM1)
    REAL(KIND=r8), INTENT(OUT  ) :: PHIS(IM,JM)
    REAL(KIND=r8), INTENT(OUT  ) :: PS(IM,JM)
    REAL(KIND=r8), INTENT(OUT  ) :: PIF(IM,JM)
    REAL(KIND=r8), INTENT(OUT  ) :: TF(IM,JM,KM)
    REAL(KIND=r8), INTENT(OUT  ) :: WW(JM)
    REAL(KIND=r8) :: AI
    REAL(KIND=r8) :: AJ
    REAL(KIND=r8) :: ANG
    REAL(KIND=r8) :: DD
    INTEGER       :: KPAD1
    INTEGER       :: KSIG1
    INTEGER       :: i,j

    !-----------------------------------------------------------
    !             FLAT....... LATITUDE
    !             FLON....... LONGITUDE             
    !-----------------------------------------------------------
    !2345678
    !   REAL  NS,N,M,MI
    !   DOUBLE PRECISION EP,WP,PHI,TPI,TU,TV,TT,PS                            
    !   DOUBLE PRECISION PIF,UF,VF,TF,PIA,UA,VA,TA,PIC,UC,VC,TC
    !   DOUBLE PRECISION GPX,GPY,GA,GF,GC
    !   DOUBLE PRECISION EKI,EK,PHIC       
    ! 
    !  PARAMETER(IM=146,JM=44,KM=7,JM1=JM-1,IM1=IM,KIM=KM+1,PT=20.)
    !  PARAMETER(PIS=3.141596,RAD=PIS/180.,RT=6.365E+6)         
    !  PARAMETER (DLON=2.5,DLAT=2.5) 
    !  COMMON/F/PIF(IM,JM),UF(IM1,JM1,KM),VF(IM1,JM1,KM),TF(IM,JM,KM)
    !  COMMON/A/PIA(IM,JM),UA(IM1,JM1,KM),VA(IM1,JM1,KM),TA(IM,JM,KM)
    !  COMMON/Z/PIZ(IM,JM),UZ(IM1,JM1,KM),VZ(IM1,JM1,KM),TZ(IM,JM,KM)

    !  COMMON/C/PIC(IM,JM),UC(IM1,JM1,KM),VC(IM1,JM1,KM),TC(IM,JM,KM)
    !  COMMON/G/GC(IM1,JM1,KM),GA(IM1,JM1,KM),GF(IM1,JM1,KM)
    !  COMMON/TEN/TPI(IM,JM),TU(IM1,JM1,KM),TV(IM1,JM1,KM),TT(IM,JM,KM)    
    !  COMMON/COORD/NS(IM,JM),M(JM),MI(JM),FLAT(JM),FLON(IM)
    !  COMMON/BLOCOS/KSE(IM,JM),KSU(IM1,JM1),KSV(IM1,JM1)
    !  COMMON/ETAS/EK(KM),EKI(KIM)
    !  COMMON/TOPOG/PHIS(IM,JM),PRF(IM,JM),PS(IM,JM)
    !  COMMON/DIAG/EP(IM,JM,KIM),WP(IM,JM,KM)
    !  COMMON/PHISS/PHI(IM,JM,KIM),PHIC(IM,JM,KIM)
    !  COMMON/GRAD/GPX(IM1,JM1,KM),GPY(IM1,JM1,KM)                   
    !  COMMON/AUX/WW(JM),Q(IM,JM,KM),ALFA(IM,JM,KM)
    FLAT=0.0
    FLON=0.0
    M =0.0
    MI=0.0
    EK=0.0
    EKI=0.0
    PRF=0.0
    NS=0.0
    KSE=0.0
    KSU=0.0
    KSV=0.0
    PHIS=0.0
    PS=0.0
    PIF=0.0
    TF=0.0
    WW=0.0
    AI=0.0
    AJ=0.0
    ANG=0.0
    DD=0.0
    !2345678  
    !----------------------------------------------------------          
    !           DETERMINACAO DAS LATITUDES E LONGITUDES
    !----------------------------------------------------------
    DO J=1,JM
       AJ=J
       FLAT(J)=FLAT1+(AJ-1.)*DLAT
       ANG=FLAT(J)*RAD
       MI(J)=RT*COS(ANG)
       M(J)=1./MI(J)  ! 1/(RT*COS(ANG))
    END DO
    DO I=1,IM
       AI=I
       FLON(I)=FLON1+(AI-1.)*DLON
    END DO
    !-----------------------------------------------------------
    !            DETERMINACAO DE ETAS
    !------------------------------------------------------------
    DD=1000.-PT
    !               NIVELEIS PRINCIPAIS
    !               -------------------
    EK(1)=(100.-PT)/DD
    EK(2)=(200.-PT)/DD
    EK(3)=(300.-PT)/DD
    EK(4)=(500.-PT)/DD
    EK(5)=(700.-PT)/DD
    EK(6)=(850.-PT)/DD
    EK(7)=(960.-PT)/DD
    !             NIVEIS INTERMEDIOS(K+1/2)
    !             ------------------------
    EKI(1)=0.
    EKI(2)=(150.-PT)/DD
    EKI(3)=(250.-PT)/DD
    EKI(4)=(400.-PT)/DD
    EKI(5)=(620.-PT)/DD
    EKI(6)=(780.-PT)/DD
    EKI(7)=(920.-PT)/DD
    EKI(8)=1.                   

    !------------------------------------------------------
    !     CALCULO DE NS,KSE SIN TOPO.,NS VARIAVEI QUE E
    !     FUNCAO DA PRESAO DE REFERENCIA DA TOPOG.KSE E UMA
    !     VARIAVEI ENTERA QUE VAI INDICAR ATE ONDE DEVE SER 
    !     INTEGRADO NA VERTICAL. PHIS E ALTURA GEOPOTENCIAL
    !     EN ESTA SUBROUTINA INGRESA SEM TOPO. POR TANTO,
    !     NS=1,KSE=KMI,PHIS=0.//PS E PRERSAO DA SUPERICIE           
    !----------------------------------------------------- 
    KSIG1=KSIG

    DO I=1,IM
       DO J=1,JM
          PRF(I,J)=1000.   !SIN TOPO

          IF(KSIG1.EQ.1)THEN
             NS(I,J)=1.
          ELSE
             NS(I,J)=(PRF(I,J)-PT)/DD
          ENDIF

          KSE(I,J)=INT(NS(I,J)*7.1)+1
          PHIS(I,J)=0.
          PS(I,J)=PRF(I,J)
          PIF(I,J)=(PS(I,J)-PT)/NS(I,J)
       END DO
    END DO
    !-------------------------------------------------------
    !    KSU E KSV SAO IDENTIFICADORES DO MAXIMO NIVEL ATE ONDE
    !    DEVE SER INTEGRADO NA VERTICAL, NESTA SUBROUTINA INGR.
    !    SEM TOPO. POR TANTO KSU=KSV=KSE.       
    !-----------------------------------------------------
    DO J=1,JM1
       DO I=1,IM1
          KSU(I,J)=KSE(I,J)
          KSV(I,J)=KSE(I,J)                   
       END DO
    END DO
    !---------------------------------------------------
    !            TEMPERATURA PADRAO
    !   ESTAS TEMPERATURA SERAO USADOS PARA DEFINIR AS
    !   ALTURAS GEOPOTENCIAIS SOBRE A TOPOGRAFIA; ADICIONAL-
    !   PUEDE SER USADO COMO CONDICAO INICIAL EN REPOSO
    !---------------------------------------------------

    KPAD1=KPAD
    IF(KPAD1.EQ.1)THEN
       DO I=1,IM
          DO J=1,JM
             TF(I,J,1)=273.3-55.
             TF(I,J,2)=273.3-55.
             TF(I,J,3)=273.3-44.5         
             TF(I,J,4)=273.3-21.2
             TF(I,J,5)=273.3-4.6
             TF(I,J,6)=273.3+5.5
             TF(I,J,7)=273.3+12.1
          END DO
       END DO
       !-----------------------------------------------------
       ! temperatura media en tropico para geopo da topo phis 
       ! con condiciones iniciales diferencte a padrao/modifi
       !--------------------------------------------------
    ELSE

       DO I=1,IM
          DO J=1,JM
             TF(I,J,1)=273.3-68.06
             TF(I,J,2)=273.3-53.56
             TF(I,J,3)=273.3-32.20         
             TF(I,J,4)=273.3-6.26
             TF(I,J,5)=273.3+9.08
             TF(I,J,6)=273.3+16.88
             TF(I,J,7)=273.3+22.22
          END DO
       END DO
    ENDIF

    !------------------------------------------------------- 
    !               PESOS PARA A FORNTEIRA NORTE E SUL SPONJA
    !-------------------------------------------------------
    WW(1)=0.
    WW(2)=0.4
    WW(3)=0.7
    WW(4)=0.9
    WW(5)=1.
    !---------------------------
    WW(JM)=0.
    WW(JM-1)=0.4
    WW(JM-2)=0.7
    WW(JM-3)=0.9
    WW(JM-4)=1.
    !-----------------------------        
    RETURN
  END SUBROUTINE INICIA

  !**************************************************************
  SUBROUTINE PREFE(&
       KSIG      ,&!INTEGER      , INTENT(IN   ) :: KSIG
       LLPK      ,&!INTEGER      , INTENT(IN   ) :: LLPK
       LY        ,&!INTEGER      , INTENT(IN   ) :: LY
       IM        ,&!INTEGER      , INTENT(IN   ) :: IM
       JM        ,&!INTEGER      , INTENT(IN   ) :: JM
       KM        ,&!INTEGER      , INTENT(IN   ) :: KM
       KIM       ,&!INTEGER      , INTENT(IN   ) :: KIM
       IM1       ,&!INTEGER      , INTENT(IN   ) :: IM1
       JM1       ,&!INTEGER      , INTENT(IN   ) :: JM1
       PT        ,&!REAL(KIND=r8), INTENT(IN   ) :: PT
       RS        ,&!REAL(KIND=r8), INTENT(IN   ) :: RS      !RS=287.
       EK        ,&!REAL(KIND=r8), INTENT(IN   ) :: EK (KIM)
       EKI       ,&!REAL(KIND=r8), INTENT(IN   ) :: EKI(KIM)
       PRF       ,&!REAL(KIND=r8), INTENT(OUT  ) :: PRF(IM,JM)
       NS        ,&!REAL(KIND=r8), INTENT(OUT  ) :: NS(IM,JM)
       PS        ,&!REAL(KIND=r8), INTENT(OUT  ) :: PS(IM,JM)
       PIF       ,&!REAL(KIND=r8), INTENT(OUT  ) :: PIF(IM,JM)
       KSU       ,&!INTEGER      , INTENT(OUT  ) :: KSU(IM1,JM1)
       KSV       ,&!INTEGER      , INTENT(OUT  ) :: KSV(IM1,JM1)
       KSE       ,&!INTEGER      , INTENT(OUT  ) :: KSE(IM,JM)
       TF        ,&!REAL(KIND=r8), INTENT(INOUT) :: TF(IM,JM,KM)
       PHIS       &!REAL(KIND=r8), INTENT(INOUT) :: PHIS(IM,JM)
       )       !PRESIONES DE REFERENCIA
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: KSIG
    INTEGER      , INTENT(IN   ) :: LLPK
    INTEGER      , INTENT(IN   ) :: LY
    INTEGER      , INTENT(IN   ) :: IM
    INTEGER      , INTENT(IN   ) :: JM
    INTEGER      , INTENT(IN   ) :: KM
    INTEGER      , INTENT(IN   ) :: KIM
    INTEGER      , INTENT(IN   ) :: IM1
    INTEGER      , INTENT(IN   ) :: JM1
    REAL(KIND=r8), INTENT(IN   ) :: PT
    REAL(KIND=r8), INTENT(IN   ) :: RS     !RS=287.
    REAL(KIND=r8), INTENT(IN   ) :: EK (KIM)
    REAL(KIND=r8), INTENT(IN   ) :: EKI(KIM)
    REAL(KIND=r8), INTENT(OUT  ) :: PRF(IM,JM)
    REAL(KIND=r8), INTENT(OUT  ) :: NS(IM,JM)
    REAL(KIND=r8), INTENT(OUT  ) :: PS(IM,JM)
    REAL(KIND=r8), INTENT(OUT  ) :: PIF(IM,JM)
    INTEGER      , INTENT(OUT  ) :: KSU(IM1,JM1)
    INTEGER      , INTENT(OUT  ) :: KSV(IM1,JM1)
    INTEGER      , INTENT(OUT  ) :: KSE(IM,JM)
    REAL(KIND=r8), INTENT(INOUT) :: TF(IM,JM,KM)
    REAL(KIND=r8), INTENT(INOUT) :: PHIS(IM,JM)


    REAL(KIND=r8) :: TS(IM,JM)
    REAL(KIND=r8) :: DENO1
    REAL(KIND=r8) :: FF3
    REAL(KIND=r8) :: FF4
    REAL(KIND=r8) :: PK1
    REAL(KIND=r8) :: PK2
    REAL(KIND=r8) :: PK21
    REAL(KIND=r8) :: PK3

    INTEGER :: KK1
    INTEGER :: KK2
    INTEGER :: kmedio
    INTEGER :: kalto
    INTEGER :: KS
    INTEGER :: KSIG1
    INTEGER :: i,j,k
    !----------------------------------------------------------------
    !            L............O VALOR DE I (LON.) ONDE ESTA CENTRADO
    !                         A CORDILHEIRA DOS ANDES SOBRE CHILE
    !            LY...........O VALOR DE J (LAT.) ONDE ESTA 17.5 SUL
    !                         EXACTAMENTE SOBRE BOLIVIA 
    !            RS.......... CONSTANTE DE AR SECO PARA CALCULAR ALT
    !                         GEOPOTENCIAL SOBRE A TOPOGRAFIA
    !            PRF......... PRESAO DE REFERENCA DA TOPOGRAFIA
    !----------------------------------------------------------------              
    !2345678
    ! REAL  NS,N,M,MI
    ! DOUBLE PRECISION EP,WP,PHI,TPI,TU,TV,TT,PS                            
    ! DOUBLE PRECISION PIF,UF,VF,TF,PIA,UA,VA,TA,PIC,UC,VC,TC
    ! DOUBLE PRECISION GPX,GPY,GA,GF,GC
    ! DOUBLE PRECISION EKI,EK,PHIC
    !
    ! PARAMETER(IM=146,JM=44,KM=7,JM1=JM-1,IM1=IM,KIM=KM+1,PT=20.)                              
    ! PARAMETER(RS=287.,JM2=JM-2,IM2=IM-2)
    !COMMON/F/PIF(IM,JM),UF(IM1,JM1,KM),VF(IM1,JM1,KM),TF(IM,JM,KM)
    !COMMON/A/PIA(IM,JM),UA(IM1,JM1,KM),VA(IM1,JM1,KM),TA(IM,JM,KM)
    !COMMON/Z/PIZ(IM,JM),UZ(IM1,JM1,KM),VZ(IM1,JM1,KM),TZ(IM,JM,KM)

    !COMMON/C/PIC(IM,JM),UC(IM1,JM1,KM),VC(IM1,JM1,KM),TC(IM,JM,KM)
    !COMMON/G/GC(IM1,JM1,KM),GA(IM1,JM1,KM),GF(IM1,JM1,KM)
    !COMMON/TEN/TPI(IM,JM),TU(IM1,JM1,KM),TV(IM1,JM1,KM),TT(IM,JM,KM)    
    !COMMON/COORD/NS(IM,JM),M(JM),MI(JM),FLAT(JM),FLON(IM)
    !COMMON/BLOCOS/KSE(IM,JM),KSU(IM1,JM1),KSV(IM1,JM1)
    !COMMON/ETAS/EK(KM),EKI(KIM)
    !COMMON/TOPOG/PHIS(IM,JM),PRF(IM,JM),PS(IM,JM)
    !COMMON/DIAG/EP(IM,JM,KIM),WP(IM,JM,KM)
    !COMMON/PHISS/PHI(IM,JM,KIM),PHIC(IM,JM,KIM)
    !COMMON/GRAD/GPX(IM1,JM1,KM),GPY(IM1,JM1,KM)                   
    !COMMON/AUX/WW(JM),Q(IM,JM,KM),ALFA(IM,JM,KM)

    !DIMENSION PRF1(IM2,JM2),PHIS1(IM2,JM2),TS(IM,JM)
    !--------------------------------------------------------
    !       DEFICAO DAS PRESOES DE REFERENCIA DA TOPOGRAFIA   
    !----------------------------------------------------------
    !       KSE 6666666666666666666666666666666666
    !---------------------------------------------------------                                          
    DO J=LY-10,LY+1
       DO I=LLPK-2,LLPK+2
          PRF(I,J)=920.
       END DO
    END DO
    !--------------------------------------------------------
    DO J=LY-2,LY
       PRF(LLPK+3,J)=920.
    END DO
    !-------------------------------------------------------
    DO J=LY+5,LY+9
       PRF(LLPK-6,J)=920.
    END DO
    !------------------------------------------------------
    DO J=LY+4,LY+9
       PRF(LLPK-5,J)=920.
    END DO
    !------------------------------------------------------
    DO J=LY+3,LY+10 
       PRF(LLPK-4,J)=920.
    END DO
    !-----------------------------------------------------
    DO J=LY+2,LY+11
       PRF(LLPK-3,J)=920.
    END DO
    !------------------------------------------------------
    DO J=LY+2,LY+12
       PRF(LLPK-2,J)=920.
    END DO
    !-----------------------------------------------------
    DO J=LY+9,LY+12      
       PRF(LLPK-1,J)=920.
    END DO
    DO J=LY+10,LY+12   
       PRF(LLPK,J)=920.
    END DO
    !----------------------------------------------------
    DO J=LY+2,LY+4
       PRF(LLPK-1,J)=920.
    END DO
    !----------------------------------------------------
    DO J=LY+2,LY+3
       PRF(LLPK,J)=920.
    END DO
    !---------------------------------------------------- 
    PRF(LLPK+1,LY+2)=920.
    !---------------------------------------------------
    !     KSE 555555555555555555555555555555555555555555
    !---------------------------------------------------

    kmedio=1
    IF(kmedio.EQ.1)THEN
       DO J=LY-8,LY+1
          DO I=LLPK-1,LLPK+1
             PRF(I,J)=780.
          END DO
       END DO
       DO J=LY-2,LY
          PRF(LLPK+2,J)=780.
       END DO
       !2345678-----------------------------------------

       PRF(LLPK,LY+2)=780.

       DO J=LY+2,LY+3
          PRF(LLPK-1,J)=780.
       END DO

       DO J=LY+2,LY+4
          PRF(LLPK-2,J)=780
       END DO


       DO J=LY+3,LY+10
          PRF(LLPK-3,J)=780.
       END DO

       DO J=LY+4,LY+9
          PRF(LLPK-4,J)=780.
       END DO

       DO J=LY+5,LY+8
          PRF(LLPK-5,J)=780.
       END DO
       !-------------------------------------------
       DO J=LY+9,LY+11
          PRF(LLPK-2,J)=780.
       END DO
       DO J=LY+10,LY+11
          PRF(LLPK-1,J)=780.
       END DO
       !---------------------------------------------
    ENDIF
    !--------------------------------------------------
    !    KSE 444444444444444444444444444444444444444444
    !--------------------------------------------------
    kalto=1
    IF(kalto.EQ.1)THEN
       DO J=LY-6,LY+1
          PRF(LLPK,J)=620.
       END DO
       DO J=LY-2,LY
          PRF(LLPK+1,J)=620.
       END DO
       !---------------------------------------------=
       PRF(LLPK-1,LY+2)=620.
       PRF(LLPK-2,LY+3)=620.
       PRF(LLPK-3,LY+4)=620.
       PRF(LLPK-3,LY+9)=620.
       PRF(LLPK-2,LY+10)=620.
       DO J=LY+5,LY+8
          PRF(LLPK-4,J)=620.
       END DO
    ENDIF
    !--------------------------------------------
    DO J=LY+5,LY+9
       PRF(LLPK-6,J)=1000.
    END DO
    DO J=LY+4,LY+9
       PRF(LLPK-5,J)=920.
    END DO
    DO J=LY+4,LY+8
       PRF(LLPK-4,J)=780.
       PRF(LLPK-3,J)=620.
       PRF(LLPK-2,J)=780.
       PRF(LLPK-1,J)=920.
    END DO
    PRF(LLPK-3,LY+8)=780
    PRF(LLPK-2,LY+8)=620
    PRF(LLPK-2,LY+9)=780
    PRF(LLPK-1,LY+8)=780
    PRF(LLPK,LY+8)=920

    PRF(LLPK-4,LY+9)=920.
    PRF(LLPK-4,LY+10)=1000.

    PRF(LLPK-3,LY+9)=780.
    PRF(LLPK-3,LY+10)=920.
    PRF(LLPK-3,LY+11)=1000.

    PRF(LLPK-2,LY+9)=780.
    PRF(LLPK-2,LY+10)=920.
    PRF(LLPK-2,LY+11)=1000.
    PRF(LLPK-2,LY+12)=1000.

    PRF(LLPK-1,LY+9)=780.
    PRF(LLPK-1,LY+10)=780.
    PRF(LLPK-1,LY+11)=920.
    PRF(LLPK-1,LY+12)=1000.

    PRF(LLPK,LY+9)=920.
    PRF(LLPK,LY+10)=780.
    PRF(LLPK,LY+11)=920.
    PRF(LLPK,LY+12)=1000.  
    !==================================
    PRF(LLPK+1,LY+10)=920.
    PRF(LLPK+1,LY+11)=920.
    PRF(LLPK+1,LY+12)=920.         
    !--------------------------------------

    PRF(LLPK-2,LY)=1000.
    PRF(LLPK-1,LY)=920.
    PRF(LLPK,LY)=780.

    PRF(LLPK-2,LY-1)=1000.
    PRF(LLPK-1,LY-1)=920.
    PRF(LLPK,LY-1)=780.


    PRF(LLPK-2,LY-2)=1000.
    PRF(LLPK-1,LY-2)=920.
    PRF(LLPK,LY-2)=780.
    !------------------------------------------------
    !           hasta aqui perfecto
    !-------------------------------------------------


    PRF(LLPK-5,LY+8)=1000.
    PRF(LLPK-5,LY+9)=1000.
    PRF(LLPK-4,LY+8)=920.

    PRF(LLPK-4,LY+9)=1000.
    PRF(LLPK-3,LY+10)=1000.
    PRF(LLPK-3,LY+9)=920.


    !   
    !--------------------------------------------------------
    !---------------------------------------------------------
    !      CALCULO DE NS,KSE E PRESAO DA SUPERFICIE SOBRE A 
    !      TOPOGRAFIA. KSE PARA OS PONTOS ONDE EXISTE TOPOGRAFIA
    !      SERA DIFRERENTE DO CALCULADO NA SUB. INICIA. 
    ! 
    !      PARA TRABALHAR EM COORD. SIGMA SE DEVE MANTER NS=1, MAS
    !      OS DADOS INICIAS DEVEN SER AJUSTADOS A NIVEIS SIGMA
    !
    !      PARA A PRESAO DA SUPERFICIE SOBRE A TOPOGRAFIA PODE SER
    !      USADO AS MESMAS PRESOES DE REFERENCA(COMO HACEMOS EN ESTE)
    !      MAS CUALQUIER MUDANZA PODE-SE FAZER NA SUBROUTINA CONDI
    !------------------------------------------------------------             
    DENO1=1000.-PT
    KSIG1=KSIG

    DO I=1,IM
       DO J=1,JM 

          IF(KSIG1.EQ.1)THEN              
             NS(I,J)=1.

          ELSE            

             NS(I,J)=(PRF(I,J)-PT)/DENO1 

          ENDIF

          KSE(I,J)=INT(NS(I,J)*7.1)+1
          PS(I,J)=PRF(I,J)                         
          PIF(I,J)=(PS(I,J)-PT)/NS(I,J)

          IF(PRF(I,J).EQ.1000.)TS(I,J)=19.3+273.3
          IF(PRF(I,J).EQ.920.)TS(I,J)=9.8+273.3
          IF(PRF(I,J).EQ.780.)TS(I,J)=1.0+273.3
          IF(PRF(I,J).EQ.620.)TS(I,J)=-10.7+273.3

       END DO
    END DO
    !-------------------------------------------------------------
    !
    OPEN(21,FILE='KSE.dat',STATUS='UNKNOWN')
    WRITE(21,99)((KSE(I,J),I=1,IM),J=1,JM)
99  FORMAT(1X,I2)


    !--------------------------------------------------------                                                                 
    !      CALCULO DA ALTURA GEOPOTENCIAL DA TOPOGRAFIA
    !--------------------------------------------------------
    DO J=1,JM
       DO I=1,IM
          KS=KSE(I,J)
          IF(KS.LT.KIM)THEN
             DO K=KS,KIM-1                  
                FF3=PIF(I,J)*(EKI(K+1)-EKI(K))
                FF4=PIF(I,J)*(EKI(K+1)+EKI(K))+2.*PT
                PHIS(I,J)=PHIS(I,J)+2.*RS*TF(I,J,K)*FF3/FF4  
             END DO
          ENDIF
       END DO
    END DO
    !-------------------------------------------------------------


    !          grafico de pHIS
    !---------------------------------------------------------------
    !             DO 393 I=2,IM-1
    !             DO 393 J=2,JM-1
    !          PHIS1(I-1,J-1)=PHIS(I,J)                       
    !393       CONTINUE
    !23456789
    !       OPEN(22,FILE='topo.dat',STATUS='UNKNOWN',
    !     * form='unformatted',access='direct',recl=IM2*JM2)                               
    !       WRITE(22,rec=1)PHIS1

    !---------------------------------------------------------
    !         CALCULO DOS IDENTIFICADORES DO VENTO, DE MANEIRA
    !         QUE DENTRO DA TOPOGRAFIA NAO E RECONHECIDO E ALEM
    !         DE ISSO FAZ QUE AS PERPENDICULARE AOS BLOCOS SEJAM
    !         NULOS
    !-----------------------------------------------------------
    DO J=1,JM1
       DO I=1,IM1
          KSU(I,J)=KSE(I,J)
          KSV(I,J)=KSE(I,J)
       END DO
    END DO
    !----------------------------------------------------------
    DO I=1,IM1
       DO J=1,JM1
          KK1=KSE(I+1,J)-KSE(I,J)
          IF(KK1.LE.0)KSU(I,J)=KSE(I+1,J)
          KK2=KSE(I,J+1)-KSE(I,J)
          IF(KK2.LE.0)KSV(I,J)=KSE(I,J+1)
       END DO
    END DO
    !--------------------------------------------------------
    !*********************************************************** 
    !
    !           TEMperatura INICIAL PRA   SIGMA          
    !
    !***********************************************************
    KSIG1=KSIG
    IF(KSIG1.EQ.1)THEN
       ! 
       DO I=1,IM
          DO J=1,JM
             DO K=1,KM
                PK1=EK(K)*(1000.-PT)+PT
                PK2=EK(K)*(PS(I,J)-PT)+PT
                PK21=PK2/PK1
                TF(I,J,K)=TF(I,J,K)*PK21**0.286
                PK3=EK(7)*(PS(I,J)-PT)+PT

                TF(I,J,7)=TS(I,J)*(PK3/PS(I,J))**0.286 

             END DO
          END DO
       END DO
       !---------------------------------------------------
    ENDIF

    !------------------------------------------------------------             
    RETURN
  END SUBROUTINE PREFE

  !*******************************************************************
  SUBROUTINE FORCE(&
       IM       ,&!INTEGER      , INTENT(IN   ) :: IM
       JM       ,&!INTEGER      , INTENT(IN   ) :: JM
       KM       ,&!INTEGER      , INTENT(IN   ) :: KM
       LLPK     ,&!INTEGER      , INTENT(IN   ) :: LLPK
       LX       ,&!REAL(KIND=r8), INTENT(IN   ) :: LX
       LY       ,&!INTEGER      , INTENT(IN   ) :: LY
       FLO1     ,&!REAL(KIND=r8), INTENT(IN   ) :: FLO1
       FLA1     ,&!REAL(KIND=r8), INTENT(IN   ) :: FLA1
       LYY      ,&!REAL(KIND=r8), INTENT(IN   ) :: LYY
       AIN      ,&!REAL(KIND=r8), INTENT(IN   ) :: AIN
       PIS      ,&!REAL(KIND=r8), INTENT(IN   ) :: PIS !PIS=3.141596
       RAD      ,&!REAL(KIND=r8), INTENT(IN   ) :: RAD !RAD=PIS/180.
       FLON     ,&!REAL(KIND=r8), INTENT(IN   ) :: FLON(IM)
       FLAT     ,&!REAL(KIND=r8), INTENT(IN   ) :: FLAT(JM)
       EK       ,&!REAL(KIND=r8), INTENT(IN   ) :: EK(KM)
       KSE      ,&!INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
       Q        ,&!REAL(KIND=r8), INTENT(OUT  ) :: Q(IM,JM,KM)
       ktop      &!INTEGER      , INTENT(IN   ) :: ktop
       )
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: IM
    INTEGER      , INTENT(IN   ) :: JM
    INTEGER      , INTENT(IN   ) :: KM
    INTEGER      , INTENT(IN   ) :: LLPK
    REAL(KIND=r8), INTENT(IN   ) :: LX
    INTEGER      , INTENT(IN   ) :: LY
    REAL(KIND=r8), INTENT(IN   ) :: FLO1
    REAL(KIND=r8), INTENT(IN   ) :: FLA1
    REAL(KIND=r8), INTENT(IN   ) :: LYY
    REAL(KIND=r8), INTENT(IN   ) :: AIN
    REAL(KIND=r8), INTENT(IN   ) :: PIS !PIS=3.141596
    REAL(KIND=r8), INTENT(IN   ) :: RAD !RAD=PIS/180.
    REAL(KIND=r8), INTENT(IN   ) :: FLON(IM)
    REAL(KIND=r8), INTENT(IN   ) :: FLAT(JM)
    REAL(KIND=r8), INTENT(IN   ) :: EK(KM)
    INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
    REAL(KIND=r8), INTENT(OUT  ) :: Q(IM,JM,KM)
    INTEGER      , INTENT(IN   ) :: ktop


    REAL(KIND=r8) :: AIN2
    REAL(KIND=r8) :: ANGI
    REAL(KIND=r8) :: B
    REAL(KIND=r8) :: B1
    REAL(KIND=r8) :: B2
    REAL(KIND=r8) :: B3
    REAL(KIND=r8) :: B4
    REAL(KIND=r8) :: FF
    REAL(KIND=r8) :: FFX
    REAL(KIND=r8) :: FFY
    REAL(KIND=r8) :: FLA
    REAL(KIND=r8) :: FLA1A
    REAL(KIND=r8) :: FLATC
    REAL(KIND=r8) :: FLO
    REAL(KIND=r8) :: FLO1A
    REAL(KIND=r8) :: FLONC
    REAL(KIND=r8) :: LX1
    REAL(KIND=r8) :: LY1
    REAL(KIND=r8) :: LYY1
    REAL(KIND=r8) :: QXY
    REAL(KIND=r8) :: QZ
    INTEGER       :: KSE1(IM,JM) 
    INTEGER       :: kaf
    INTEGER       :: kam
    INTEGER       :: kin
    INTEGER       :: knino
    INTEGER       :: ksacz
    INTEGER       :: KS
    INTEGER       :: i,j,k
    !-----------------------------------------------------------------
    !             Q.............DISTRIBUCAO HORIZONTAL DA FORCANTE TERM.

    !-------------------------------------------------------------------        
    !2345678
    ! REAL  NS,N,M,MI,LX,LYY
    ! DOUBLE PRECISION EP,WP,PHI,TPI,TU,TV,TT,PS                           
    ! DOUBLE PRECISION PIF,UF,VF,TF,PIA,UA,VA,TA,PIC,UC,VC,TC
    ! DOUBLE PRECISION GPX,GPY,GA,GF,GC
    ! DOUBLE PRECISION EKI,EK,PHIC
    ! 
    !PARAMETER(IM=146,JM=44,KM=7,IM1=IM,JM1=JM-1,KIM=KM+1,PT=20.)      
    !PARAMETER(PIS=3.141596,RAD=PIS/180.,IM2=IM-2,JM2=JM-2)

    !COMMON/F/PIF(IM,JM),UF(IM1,JM1,KM),VF(IM1,JM1,KM),TF(IM,JM,KM)
    !COMMON/A/PIA(IM,JM),UA(IM1,JM1,KM),VA(IM1,JM1,KM),TA(IM,JM,KM)
    !COMMON/Z/PIZ(IM,JM),UZ(IM1,JM1,KM),VZ(IM1,JM1,KM),TZ(IM,JM,KM)

    !COMMON/C/PIC(IM,JM),UC(IM1,JM1,KM),VC(IM1,JM1,KM),TC(IM,JM,KM)
    !COMMON/G/GC(IM1,JM1,KM),GA(IM1,JM1,KM),GF(IM1,JM1,KM)
    !COMMON/TEN/TPI(IM,JM),TU(IM1,JM1,KM),TV(IM1,JM1,KM),TT(IM,JM,KM)    
    !COMMON/COORD/NS(IM,JM),M(JM),MI(JM),FLAT(JM),FLON(IM)
    !COMMON/BLOCOS/KSE(IM,JM),KSU(IM1,JM1),KSV(IM1,JM1)
    !COMMON/ETAS/EK(KM),EKI(KIM)
    !COMMON/TOPOG/PHIS(IM,JM),PRF(IM,JM),PS(IM,JM)
    !COMMON/DIAG/EP(IM,JM,KIM),WP(IM,JM,KM)
    !COMMON/PHISS/PHI(IM,JM,KIM),PHIC(IM,JM,KIM)
    !COMMON/GRAD/GPX(IM1,JM1,KM),GPY(IM1,JM1,KM)                   
    !COMMON/AUX/WW(JM),Q(IM,JM,KM),ALFA(IM,JM,KM)

    !DIMENSION KSE1(IM,JM)                             
    !2345678
    !---------------------------------------------------------------
    !   PARA NO TENER DENTRO DA TOPOGRAFIA QQ CON OU SIN TOPO
    !------------------------------------------------------=------
    IF(KTOP.EQ.1)THEN
       OPEN(19,FILE='KSE.dat',STATUS='OLD')
       READ(19,99)((KSE1(I,J),I=1,IM),J=1,JM)
99     FORMAT(1X,I2)

       CLOSE(19)
    ENDIF
    Q=0.0
    !******************************************************

    knino=1
    IF(knino.EQ.0)THEN
       !---------------------------------------------------
       !                        EL N I N O  3
       !-------------------------------------------------
       !      VARIACAO HORIzONTAL CIRCULAR 
       !------------------------------------------------------------

       FLO1A=-130
       !ccccc              FLA1A=-7.5
       FLA1A=-30
       LX1=25.
       LYY1=5.

       DO I=1,LLPK-5                          
          !cccccccc          DO 100 J=LY-5,JM 
          DO J=1,JM-10   

             FLO=FLON(I)-FLO1A
             FLA=FLAT(J)-FLA1A
             ANGI=AIN*RAD

             FFX=FLO*COS(ANGI)-FLA*SIN(ANGI)
             FFY=FLO*SIN(ANGI)+FLA*COS(ANGI)

             IF(ABS(FLO).LE.30.AND.ABS(FLA).LE.30)THEN             
                FF=(FFX/LX1)**2+(FFY/LYY1)**2
                QXY=EXP(-FF)   ! GAUSSIANA
             ELSE
                QXY=0.
             ENDIF
             !------------------------------------------------------------
             !          VARIACION VERTICAL SINOSOIDAL              
             !-----------------------------------------------------------
             !     SI LA FORCANTE DEBE ESTAR  SIN RESTRICION KSE1=8) 
             IF(KTOP.EQ.1)THEN
                KS=KSE1(I,J)                   
             ELSE
                ks=KSE(i,j)
             ENDIF
             !c               B=2.    !max en 400mb
             B=0.

             DO K=1,KS-1
                B1=PIS**2
                B2=B1*(1.+EXP(-B))  
                B3=EXP(-B*EK(K))        
                B4=(B**2+B1)*B3/B2

                QZ=B4*SIN(PIS*EK(K)) 
                !cccccccc    Q(I,J,K)=QXY*QZ
                Q(I,J,K)=QXY*QZ*2.0
             END DO
          END DO
       END DO

       !----------------------------------------------

    ENDIF
PRINT*,' A M A Z O N A S'
    kam=0
    IF(kam.EQ.0)THEN
       !----------------------------------------------------
       !              A M A Z O N A S
       !----------------------------------------
       !      VARIACAO HORIzONTAL CIRCULAR
       !---------------------------------------
       DO J=LY-2,LY+10
          DO I=LLPK-2,LLPK+20

             FLO=FLON(I)-FLO1
             FLA=FLAT(J)-FLA1
             ANGI=AIN*RAD
             FFX=FLO*COS(ANGI)-FLA*SIN(ANGI)
             FFY=FLO*SIN(ANGI)+FLA*COS(ANGI)

             IF(ABS(FLO).LE.30.AND.ABS(FLA).LE.30)THEN             
                FF=(FFX/LX)**2+(FFY/LYY)**2
                QXY=EXP(-FF)   ! GAUSSIANA
             ELSE
                QXY=0.
             ENDIF
             !------------------------------------------------------------
             !          VARIACION VERTICAL SINOSOIDAL              
             !-----------------------------------------------------------
             !     SI LA FORCANTE DEBE ESTAR  SIN RESTRICION KSE1=8)           
             IF(KTOP.EQ.1)THEN
                KS=KSE1(I,J)                   
             ELSE
                ks=KSE(i,J)
             ENDIF
             !c               B=2.    !max en 400mb
             B=0.

             DO K=1,KS-1
                B1=PIS**2
                B2=B1*(1.+EXP(-B))  
                B3=EXP(-B*EK(K))        
                B4=(B**2+B1)*B3/B2

                QZ=B4*SIN(PIS*EK(K))   
                Q(I,J,K)=QXY*QZ
             END DO
          END DO
       END DO

       !-----------------------------------------------------------------
       DO J=LY-5,LY+2
          DO I=LLPK-7,LLPK-1
             DO K=1,KM
                Q(I,J,K)=0.
             END DO
          END DO
       END DO

       !--------------------------------------------------------
       DO J=LY+8,LY+10
          DO I=LLPK-10,LLPK-1
             DO K=1,KM
                Q(I,J,K)=0.
             END DO
          END DO
       END DO

       !---------------------------------------------------------
    ENDIF
    !------------------------------------------------
    !
    !
    !              Z    A     C     Z   
    !              
    !
    !-------------------------------------------------
PRINT*,' Z    A     C   Z'
    ksacz=1
    IF(ksacz.EQ.1)THEN
       FLATC=-24.375
       FLONC=-45                     !   -45.0
       AIN2=45.0
       lx1=10.
       ly1=7.5

       DO j=1,ly+5
          DO i=LLPK+2,LLPK+35

             FLO=FLON(i)-FLONC
             FLA=FLAT(j)-FLATC
             ANGI=AIN2*RAD
             FFX=FLO*COS(ANGI)-FLA*SIN(ANGI)
             FFY=FLO*SIN(ANGI)+FLA*COS(ANGI)

             IF(ABS(FLO).LE.20.AND.ABS(FLA).LE.20)THEN
                FF=(FFX/LX1)**2+(FFY/LY1)**2
                QXY=EXP(-FF)   ! GAUSSIANA
             ELSE
                QXY=0.
             ENDIF

             B=0.

             DO K=1,KS-1
                B1=PIS**2
                B2=B1*(1.+EXP(-B))  
                B3=EXP(-B*EK(K))        
                B4=(B**2+B1)*B3/B2

                QZ=B4*SIN(PIS*EK(K))   
                Q(I,J,K)=QXY*QZ
             END DO
          END DO
       END DO

       !-----------------------------------------------------------------
    ENDIF

    !                        A F R I C A
    !-----------------------------------------
    !      VARIACAO HORIzONTAL CIRCULAR 
    !------------------------------------------------------------

    kaf=0
    IF(kaf.EQ.0)THEN
       FLO1A=26.25
       FLA1A=-11.25

       DO I=LLPK+21,IM-30
          DO J=LY-5,JM 

             FLO=FLON(I)-FLO1A
             FLA=FLAT(J)-FLA1A
             ANGI=AIN*RAD

             FFX=FLO*COS(ANGI)-FLA*SIN(ANGI)
             FFY=FLO*SIN(ANGI)+FLA*COS(ANGI)

             IF(ABS(FLO).LE.30.AND.ABS(FLA).LE.30)THEN             
                FF=(FFX/LX)**2+(FFY/LYY)**2
                QXY=EXP(-FF)   ! GAUSSIANA
             ELSE
                QXY=0.
             ENDIF
             !------------------------------------------------------------
             !          VARIACION VERTICAL SINOSOIDAL              
             !-----------------------------------------------------------
             !     SI LA FORCANTE DEBE ESTAR  SIN RESTRICION KSE1=8)           
             IF(KTOP.EQ.1)THEN
                KS=KSE1(I,J)                   
             ELSE
                ks=KSE(i,j)
             ENDIF
             !c               B=2.    !max en 400mb
             B=0.

             DO  K=1,KS-1
                B1=PIS**2
                B2=B1*(1.+EXP(-B))  
                B3=EXP(-B*EK(K))        
                B4=(B**2+B1)*B3/B2

                QZ=B4*SIN(PIS*EK(K))   
                Q(I,J,K)=QXY*QZ*0.8
             END DO
          END DO
       END DO

       !-----------------------------------------
    ENDIF
    !----------------------------------------------------


    !---------------------------------------------------
    !                         I N D  O N E S I  A
    !-------------------------------------------------
    kin=0
    IF(kin.EQ.0)THEN

       !---------------------------------------------
       !      VARIACAO HORIzONTAL CIRCULAR 
       !---------------------------------------------------
       FLO1A=135
       !ccc              FLA1A=-7.5
       FLA1A=-30

       LX1=25.
       !c                  LYY1=5.
       LYY1=10.

       DO I=IM-40,Im-5
          DO J=1,JM-10 

             FLO=FLON(I)-FLO1A
             FLA=FLAT(J)-FLA1A
             ANGI=AIN*RAD

             FFX=FLO*COS(ANGI)-FLA*SIN(ANGI)
             FFY=FLO*SIN(ANGI)+FLA*COS(ANGI)

             IF(ABS(FLO).LE.30.AND.ABS(FLA).LE.30)THEN
                FF=(FFX/LX1)**2+(FFY/LYY1)**2
                QXY=EXP(-FF)   ! GAUSSIANA
             ELSE
                QXY=0.
             ENDIF
             !------------------------------------------------------------
             !          VARIACION VERTICAL SINOSOIDAL              
             !-----------------------------------------------------------
             !     SI LA FORCANTE DEBE ESTAR  SIN RESTRICION KSE1=8) 
             IF(KTOP.EQ.1)THEN
                KS=KSE1(I,J)                   
             ELSE
                ks=KSE(i,j)
             ENDIF
             !c               B=2.    !max en 400mb
             B=0.

             DO  K=1,KS-1
                B1=PIS**2
                B2=B1*(1.+EXP(-B))  
                B3=EXP(-B*EK(K))        
                B4=(B**2+B1)*B3/B2

                QZ=B4*SIN(PIS*EK(K)) 
                Q(I,J,K)=QXY*QZ
             END DO
          END DO
       END DO

       !----------------------------------------------

    ENDIF
    !****************************************************8
    RETURN
  END  SUBROUTINE FORCE
  !-------------------------------------------------------------      

  !*******************************************************************
  SUBROUTINE CONDI(&
       IM      ,&! INTEGER      , INTENT(IN   ) :: IM
       JM      ,&! INTEGER      , INTENT(IN   ) :: JM
       KM      ,&! INTEGER      , INTENT(IN   ) :: KM
       KIM     ,&! INTEGER      , INTENT(IN   ) :: KIM
       IM1     ,&! INTEGER      , INTENT(IN   ) :: IM1
       JM1     ,&! INTEGER      , INTENT(IN   ) :: JM1
       DY      ,&! REAL(KIND=r8), INTENT(IN   ) :: DY
       N       ,&! REAL(KIND=r8), INTENT(IN   ) :: N    !N=1./RT
       OMEGA   ,&! REAL(KIND=r8), INTENT(IN   ) :: OMEGA!OMEGA=7.29E-5
       PT      ,&! REAL(KIND=r8), INTENT(IN   ) :: PT
       RAD     ,&! REAL(KIND=r8), INTENT(IN   ) :: RAD   !RAD=PIS/180
       RS      ,&! REAL(KIND=r8), INTENT(IN   ) :: RS    !    RS=287.
       EK      ,&! REAL(KIND=r8), INTENT(IN   ) :: EK(KM)
       EKI     ,&! REAL(KIND=r8), INTENT(IN   ) :: EKI(KM)
       FLAT    ,&! REAL(KIND=r8), INTENT(IN   ) :: FLAT(JM)
       NS      ,&! REAL(KIND=r8), INTENT(IN   ) :: NS(IM,JM)
       PRF     ,&! REAL(KIND=r8), INTENT(IN   ) :: PRF(IM,JM)
       PHI     ,&! REAL(KIND=r8), INTENT(INOUT) :: PHI(IM,JM,KIM)
       PS      ,&! REAL(KIND=r8), INTENT(INOUT) :: PS(IM,JM)
       PIF     ,&! REAL(KIND=r8), INTENT(INOUT) :: PIF(IM,JM)
       PHI1    ,&! REAL(KIND=r8), INTENT(INOUT) :: PHI1(IM,JM,KIM)
       TF      ,&! REAL(KIND=r8), INTENT(INOUT) :: TF(IM,JM,KM)
       UF      ,&! REAL(KIND=r8), INTENT(INOUT) :: UF(IM1,JM1,KM)
       KSE      &! INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
       ) 
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: IM
    INTEGER      , INTENT(IN   ) :: JM
    INTEGER      , INTENT(IN   ) :: KM
    INTEGER      , INTENT(IN   ) :: KIM
    INTEGER      , INTENT(IN   ) :: IM1
    INTEGER      , INTENT(IN   ) :: JM1
    REAL(KIND=r8), INTENT(IN   ) :: DY
    REAL(KIND=r8), INTENT(IN   ) :: N    !N=1./RT
    REAL(KIND=r8), INTENT(IN   ) :: OMEGA!OMEGA=7.29E-5
    REAL(KIND=r8), INTENT(IN   ) :: PT
    REAL(KIND=r8), INTENT(IN   ) :: RAD   !RAD=PIS/180
    REAL(KIND=r8), INTENT(IN   ) :: RS    !    RS=287.
    REAL(KIND=r8), INTENT(IN   ) :: EK(KM)
    REAL(KIND=r8), INTENT(IN   ) :: EKI(KM)
    REAL(KIND=r8), INTENT(IN   ) :: FLAT(JM)
    REAL(KIND=r8), INTENT(IN   ) :: NS(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: PRF(IM,JM)
    REAL(KIND=r8), INTENT(INOUT) :: PHI(IM,JM,KIM)
    REAL(KIND=r8), INTENT(INOUT) :: PS(IM,JM)
    REAL(KIND=r8), INTENT(INOUT) :: PIF(IM,JM)
    REAL(KIND=r8), INTENT(INOUT) :: PHI1(IM,JM,KIM)
    REAL(KIND=r8), INTENT(INOUT) :: TF(IM,JM,KM)
    REAL(KIND=r8), INTENT(INOUT) :: UF(IM1,JM1,KM)
    INTEGER      , INTENT(IN   ) :: KSE(IM,JM)

    REAL(KIND=r8) :: T0(13,8)
    REAL(KIND=r8) :: P0(13)
    REAL(KIND=r8) :: V0(13,7)
    REAL(KIND=r8) :: V1(49,7)
    REAL(KIND=r8) :: PS1(49)
    REAL(KIND=r8) :: T1(49,8)
    REAL(KIND=r8) :: AAA
    REAL(KIND=r8) :: ANG
    REAL(KIND=r8) :: BBB
    REAL(KIND=r8) :: BG1
    REAL(KIND=r8) :: BG2
    REAL(KIND=r8) :: CCC
    REAL(KIND=r8) :: DDD
    REAL(KIND=r8) :: FF1
    REAL(KIND=r8) :: FF2
    INTEGER       :: KTRO
    INTEGER       :: LLL
    INTEGER       :: i,j,k
    !2345678
    ! REAL  NS,N,M,MI
    ! DOUBLE PRECISION EP,WP,PHI,TPI,TU,TV,TT,PS                            
    ! DOUBLE PRECISION PIF,UF,VF,TF,PIA,UA,VA,TA,PIC,UC,VC,TC
    ! DOUBLE PRECISION GPX,GPY,GA,GF,GC
    ! DOUBLE PRECISION EKI,EK,PHIC
    ! 
    !PARAMETER(IM=146,JM=44,KM=7,IM1=IM,JM1=JM-1,KIM=KM+1,PT=20.)
    !PARAMETER(OMEGA=7.29E-5)
    !PARAMETER(RT=6.365E+6,PIS=3.141596,RAD=PIS/180.,RS=287.)
    !PARAMETER(N=1./RT,DLON=2.5,DLAT=2.5,DX=DLON*RAD,DY=DLAT*RAD) 

    !COMMON/F/PIF(IM,JM),UF(IM1,JM1,KM),VF(IM1,JM1,KM),TF(IM,JM,KM)
    !COMMON/A/PIA(IM,JM),UA(IM1,JM1,KM),VA(IM1,JM1,KM),TA(IM,JM,KM)
    !COMMON/Z/PIZ(IM,JM),UZ(IM1,JM1,KM),VZ(IM1,JM1,KM),TZ(IM,JM,KM)

    !COMMON/C/PIC(IM,JM),UC(IM1,JM1,KM),VC(IM1,JM1,KM),TC(IM,JM,KM)
    !COMMON/G/GC(IM1,JM1,KM),GA(IM1,JM1,KM),GF(IM1,JM1,KM)
    !COMMON/TEN/TPI(IM,JM),TU(IM1,JM1,KM),TV(IM1,JM1,KM),TT(IM,JM,KM)    
    !COMMON/COORD/NS(IM,JM),M(JM),MI(JM),FLAT(JM),FLON(IM)
    !COMMON/BLOCOS/KSE(IM,JM),KSU(IM1,JM1),KSV(IM1,JM1)
    !COMMON/ETAS/EK(KM),EKI(KIM)
    !COMMON/TOPOG/PHIS(IM,JM),PRF(IM,JM),PS(IM,JM)
    !COMMON/DIAG/EP(IM,JM,KIM),WP(IM,JM,KM)
    !COMMON/PHISS/PHI(IM,JM,KIM),PHIC(IM,JM,KIM)
    !COMMON/GRAD/GPX(IM1,JM1,KM),GPY(IM1,JM1,KM)                   
    !COMMON/AUX/WW(JM),Q(IM,JM,KM),ALFA(IM,JM,KM)
    !
    !DIMENSION PHI1(IM,JM,KIM),T0(13,8),P0(13),PS1(49)
    !DIMENSION T1(49,8)
    !DIMENSION V0(13,7),V1(49,7)
    !2345678            
    !---------------------------------------------------------------------
    !             PRESAO DE SUPERFICIE
    !             COM TOPO E SIN TOPO  e PIF
    !---------------------------------------------------------------------
    !                   lectur de p
    !                   -----------


    !---------------------------------------------------------------------
    !                   TEMPERATURA
    !--------------------------------------------------------------------
    !                  lectura de t y p
    !                 -----------------

    OPEN(90,FILE='tempe.dat',status='old')

    DO j=1,13

       READ(90,'(7x,8F5.1,1x,F7.1)')(T0(j,k),k=1,8),P0(j)

    END DO

    CLOSE(90)

    !-----------------------------------------------
    !                 presao
    !------------------------------------------------
    LLL=1
    DO I=1,12
       J=LLL

       PS1(J)=P0(I)

       PS1(J+1)=3.*P0(I)/4.+P0(I+1)/4.

       PS1(J+2)=(P0(I)+P0(I+1))*0.5

       PS1(J+3)=3.*P0(I+1)/4.+P0(I)/4.              
       LLL=LLL+4                 
    END DO

    !=============================================
    !             suavisacion de PS1
    !----------------------------------------------
    DO i=1,49-1
       PS1(i)=(PS1(i)+PS1(i+1))*0.5
    ENDDO

    DO J=1,JM          
       DO I=1,IM              
          PS(I,J)=PS1(J+2)
          PIF(I,J)=(PS(I,J)-PT)/NS(I,J)  

       END DO
    END DO

    !------------------------------------------------
    !           temperatura
    !------------------------------------------------

    DO K=1,8
       LLL=1
       DO I=1,12
          J=LLL

          T1(J,K)=T0(I,K)
          T1(J+1,K)=3.*T0(I,K)/4.+T0(I+1,K)/4.
          T1(J+2,K)=(T0(I,K)+T0(I+1,K))*0.5
          T1(J+3,K)=3.*T0(I+1,K)/4.+T0(I,K)/4.              
          LLL=LLL+4                 
       END DO
    END DO

    !-------------------------------------------
    !                       suavizando T1
    !------------------------------------------
    DO i=1,49-1
       DO k=1,8
          T1(i,K)=(T1(i,k)+T1(i+1,k))*0.5
       ENDDO
    ENDDO
    !2345678--------------------------------------------  

    DO J=1,JM          
       DO I=1,IM            
          TF(I,J,7)= T1(J+2,1)+273.3
          TF(I,J,6)= T1(J+2,2)+273.3
          TF(I,J,5)= T1(J+2,3)+273.3
          TF(I,J,4)= T1(J+2,4)+273.3
          TF(I,J,3)= T1(J+2,6)+273.3
          TF(I,J,2)= T1(J+2,7)+273.3
          TF(I,J,1)= T1(J+2,8)+273.3
          TF(i,j,7)=(TF(I,J,7)+TF(I,J,6))*0.5
          !         
       END DO
    END DO

    !----------------------------------------------------
    !       CALCLO DE ALTURA GEOPOTENCIALES PARA UUUUUUUU
    !-------------------------------------------------------

    DO J=1,JM       
       DO K=KM,1,-1
          FF1=PIF(1,J)*(EKI(K+1)-EKI(K))
          FF2=PIF(1,J)*(EKI(K+1)+EKI(K))+2.*PT
          PHI1(1,J,K)=PHI1(1,J,K+1)+2.*RS*TF(1,J,K)*FF1/FF2
       END DO
    END DO

    !----------------------------------------------------         
    DO J=1,JM                      
       DO K=KM+1,1,-1
          PHI(1,J,K)=PHI1(1,J,K)
       END DO
    END DO

    !2345678--------------------------------------------------
    DO J=2,JM1
       DO K=1,KM
          AAA=RS*TF(1,J,K)*EK(K)/(PIF(1,J)*EK(K)+PT)
          BBB=(PIF(1,J+1)-PIF(1,J-1))/(2.*DY)
          ANG=FLAT(J)*RAD
          IF(ANG.NE.0.)THEN
             CCC=N/(2.*SIN(ANG)*OMEGA)


             BG1=(PHI(1,J-1,K+1)+PHI(1,J-1,K))*0.5
             BG2=(PHI(1,J+1,K+1)+PHI(1,J+1,K))*0.5

             DDD=(BG2-BG1)/(2.*DY)
             UF(1,J,K)=-CCC*(DDD+AAA*BBB)

          ELSE

             UF(1,J,K)=UF(1,J-1,K)
          ENDIF
       END DO
    END DO

    !---------------------------------
    DO K=1,KM
       UF(1,1,K)=UF(1,2,K)
       UF(1,23,K)=(UF(1,22,K)+UF(1,24,K))*0.5
    END DO

    DO I=2,IM1
       DO J=1,JM1
          DO K=1,KM
             UF(I,J,K)=UF(1,J,K)
          END DO
       END DO
    END DO

    !-----------------------------------------------------------
    ! PRESAO DE SUPERFICIE CON TOPOGRAFIA
    !--------------------------------------------------
    DO I=1,IM
       DO J=1,JM

          !             PS(I,J)=(PS(I,J)-PT)*NS(I,J)+PT

          IF(KSE(I,J).NE.KIM)PS(I,J)=PRF(I,J)

          PIF(I,J)=(PS(I,J)-PT)/NS(I,J)
          DO K=1,KIM
             PHI1(I,J,K)=0. 
          END DO
       END DO
    END DO

    !-------------------------------------------------------
    !             lectura de vento para tropico
    !-------------------------------------------------------

    KTRO=2
    IF(KTRO.EQ.1)THEN      

       OPEN(95,FILE='vento.dat',status='old')

       DO j=1,13

          READ(95,98)(V0(j,k),k=1,7)

98        FORMAT(7x,7F5.1)

       END DO
       CLOSE(95)
       !------------------------------------------------------- 

       DO K=1,7
          LLL=1
          DO I=1,12
             J=LLL

             V1(J,K)=V0(I,K)
             V1(J+1,K)=3.*V0(I,K)/4.+V0(I+1,K)/4.
             V1(J+2,K)=(V0(I,K)+V0(I+1,K))*0.5
             V1(J+3,K)=3.*V0(I+1,K)/4.+V0(I,K)/4.              
             LLL=LLL+4                 
          END DO
       END DO
       !------------------------------------------
       !             suavizando u
       !----------------------------------------
       DO i=1,49-1
          DO k=1,7
             v1(i,K)=(v1(i,k)+v1(i+1,k))*0.5      !error grave
          ENDDO
       ENDDO
       !---------------------------------------------
       DO J=21,25
          DO I=1,IM1
             UF(I,J,7)= V1(J+2,1)
             UF(I,J,6)= V1(J+2,2)
             UF(I,J,5)= V1(J+2,3)
             UF(I,J,4)= V1(J+2,4)
             UF(I,J,3)= V1(J+2,5)
             UF(I,J,2)= V1(J+2,6)
             UF(I,J,1)= V1(J+2,7)       
          END DO
       END DO
    ENDIF
    !============================================                  
    !
    RETURN
  END SUBROUTINE CONDI

  !************************************************************************
  SUBROUTINE DIAG1(&
       NT              ,&!INTEGER      , INTENT(IN   ) :: NT
       KFRON           ,&!INTEGER      , INTENT(IN   ) :: KFRON 
       KCI             ,&!INTEGER      , INTENT(IN   ) :: KCI   
       DT              ,&!REAL(KIND=r8), INTENT(IN   ) :: DT    
       NT1             ,&!INTEGER      , INTENT(IN   ) :: NT1   
       IM              ,&!INTEGER      , INTENT(IN   ) :: IM
       JM              ,&!INTEGER      , INTENT(IN   ) :: JM
       KM              ,&!INTEGER      , INTENT(IN   ) :: KM
       KIM             ,&!INTEGER      , INTENT(IN   ) :: KIM
       IM1             ,&!INTEGER      , INTENT(IN   ) :: IM1
       JM1             ,&!INTEGER      , INTENT(IN   ) :: JM1
       DX              ,&!REAL(KIND=r8), INTENT(IN   ) :: DX
       DY              ,&!REAL(KIND=r8), INTENT(IN   ) :: DY
       PT              ,&!REAL(KIND=r8), INTENT(IN   ) :: PT
       RT              ,&!REAL(KIND=r8), INTENT(IN   ) :: RT  !RT=6.365E+6,
       N               ,&!REAL(KIND=r8), INTENT(IN   ) :: N   !N=1./RT
       WW              ,&!REAL(KIND=r8), INTENT(IN   ) :: WW(JM)
       M               ,&!REAL(KIND=r8), INTENT(IN   ) :: M(JM)
       MI              ,&!REAL(KIND=r8), INTENT(IN   ) :: MI(JM)
       EK              ,&!REAL(KIND=r8), INTENT(IN   ) :: EK(KIM)
       EKI             ,&!REAL(KIND=r8), INTENT(IN   ) :: EKI(KIM)
       PIA             ,&!REAL(KIND=r8), INTENT(IN   ) :: PIA(IM,JM)
       PIZ             ,&!REAL(KIND=r8), INTENT(IN   ) :: PIZ(IM,JM)
       PS              ,&!REAL(KIND=r8), INTENT(OUT  ) :: PS(IM,JM)
       DIV             ,&!REAL(KIND=r8), INTENT(OUT  ) :: DIV(IM,JM,KM)
       WP              ,&!REAL(KIND=r8), INTENT(OUT  ) :: WP(IM,JM,KM)
       NS              ,&!REAL(KIND=r8), INTENT(IN   ) :: NS(IM,JM)
       PIF             ,&!REAL(KIND=r8), INTENT(IN   ) :: PIF(IM,JM)
       UF              ,&!REAL(KIND=r8), INTENT(IN   ) :: UF(IM1,JM1,KM)
       VF              ,&!REAL(KIND=r8), INTENT(IN   ) :: VF(IM1,JM1,KM)
       TPI             ,&!REAL(KIND=r8), INTENT(INOUT) :: TPI(IM,JM)
       EP              ,&!REAL(KIND=r8), INTENT(INOUT) :: EP (IM,JM,KIM)
       EP1             ,&!REAL(KIND=r8), INTENT(INOUT) :: EP1(IM,JM,KIM)
       PIC             ,&!REAL(KIND=r8), INTENT(INOUT) :: PIC(IM,JM)
       KSE              &!INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
       )
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: NT
    INTEGER      , INTENT(IN   ) :: KFRON 
    INTEGER      , INTENT(IN   ) :: KCI   
    REAL(KIND=r8), INTENT(IN   ) :: DT    
    INTEGER      , INTENT(IN   ) :: NT1   
    INTEGER      , INTENT(IN   ) :: IM
    INTEGER      , INTENT(IN   ) :: JM
    INTEGER      , INTENT(IN   ) :: KM
    INTEGER      , INTENT(IN   ) :: KIM
    INTEGER      , INTENT(IN   ) :: IM1
    INTEGER      , INTENT(IN   ) :: JM1
    REAL(KIND=r8), INTENT(IN   ) :: DX
    REAL(KIND=r8), INTENT(IN   ) :: DY
    REAL(KIND=r8), INTENT(IN   ) :: PT
    REAL(KIND=r8), INTENT(IN   ) :: RT  !RT=6.365E+6,
    REAL(KIND=r8), INTENT(IN   ) :: N   !N=1./RT
    REAL(KIND=r8), INTENT(IN   ) :: WW(JM)
    REAL(KIND=r8), INTENT(IN   ) :: M(JM)
    REAL(KIND=r8), INTENT(IN   ) :: MI(JM)
    REAL(KIND=r8), INTENT(IN   ) :: EK(KIM)
    REAL(KIND=r8), INTENT(IN   ) :: EKI(KIM)
    REAL(KIND=r8), INTENT(IN   ) :: PIA(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: PIZ(IM,JM)
    REAL(KIND=r8), INTENT(OUT  ) :: PS (IM,JM)
    REAL(KIND=r8), INTENT(OUT  ) :: DIV(IM,JM,KM)
    REAL(KIND=r8), INTENT(OUT  ) :: WP (IM,JM,KM)
    REAL(KIND=r8), INTENT(IN   ) :: NS(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: PIF(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: UF(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: VF(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(INOUT) :: TPI(IM,JM)
    REAL(KIND=r8), INTENT(INOUT) :: EP (IM,JM,KIM)
    REAL(KIND=r8), INTENT(INOUT) :: EP1(IM,JM,KIM)
    REAL(KIND=r8), INTENT(INOUT) :: PIC(IM,JM)
    INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
    REAL(KIND=r8) :: A
    REAL(KIND=r8) :: A1
    REAL(KIND=r8) :: A2
    REAL(KIND=r8) :: AA
    REAL(KIND=r8) :: AX
    REAL(KIND=r8) :: AX1
    REAL(KIND=r8) :: AX2
    REAL(KIND=r8) :: AY
    REAL(KIND=r8) :: AY1
    REAL(KIND=r8) :: AY2
    REAL(KIND=r8) :: B
    REAL(KIND=r8) :: B1
    REAL(KIND=r8) :: B2
    REAL(KIND=r8) :: BB
    REAL(KIND=r8) :: CC
    REAL(KIND=r8) :: DD
    REAL(KIND=r8) :: EE
    REAL(KIND=r8) :: EPM
    REAL(KIND=r8) :: PIFX
    REAL(KIND=r8) :: PIFY
    REAL(KIND=r8) :: UMX
    REAL(KIND=r8) :: VMY

    INTEGER       :: KCI1
    INTEGER       :: KS
    INTEGER       :: NNT
    INTEGER       :: i,j,k

    !2345678
    !  REAL  NS,N,M,MI
    !  DOUBLE PRECISION EP,WP,PHI,TPI,TU,TV,TT,PS                            
    !  DOUBLE PRECISION PIF,UF,VF,TF,PIA,UA,VA,TA,PIC,UC,VC,TC
    !  DOUBLE PRECISION GPX,GPY,GA,GF,GC
    !  DOUBLE PRECISION EKI,EK,PHIC
    !   
    !  DOUBLE PRECISION UMX,VMY,PIFX,PIFY,EPM,EE,AX1,AX2,AX,AY1
    !  DOUBLE PRECISION DIV,AA,BB,CC,DD,A1,A2,A,B1,B2,B,AY2,AY,EP1
    !  PARAMETER(IM=146,JM=44,KM=7,IM1=IM,JM1=JM-1,KIM=KM+1,PT=20.) 
    !  PARAMETER(RT=6.365E+6,PIS=3.141596,RAD=PIS/180.,RS=287.)
    !  PARAMETER(N=1./RT,DLON=2.5,DLAT=2.5,DX=DLON*RAD,DY=DLAT*RAD)
    !
    ! COMMON/F/PIF(IM,JM),UF(IM1,JM1,KM),VF(IM1,JM1,KM),TF(IM,JM,KM)
    ! COMMON/A/PIA(IM,JM),UA(IM1,JM1,KM),VA(IM1,JM1,KM),TA(IM,JM,KM)
    ! COMMON/Z/PIZ(IM,JM),UZ(IM1,JM1,KM),VZ(IM1,JM1,KM),TZ(IM,JM,KM)

    ! COMMON/C/PIC(IM,JM),UC(IM1,JM1,KM),VC(IM1,JM1,KM),TC(IM,JM,KM)
    ! COMMON/G/GC(IM1,JM1,KM),GA(IM1,JM1,KM),GF(IM1,JM1,KM)
    ! COMMON/TEN/TPI(IM,JM),TU(IM1,JM1,KM),TV(IM1,JM1,KM),TT(IM,JM,KM)    
    ! COMMON/COORD/NS(IM,JM),M(JM),MI(JM),FLAT(JM),FLON(IM)
    ! COMMON/BLOCOS/KSE(IM,JM),KSU(IM1,JM1),KSV(IM1,JM1)
    ! COMMON/ETAS/EK(KM),EKI(KIM)
    ! COMMON/TOPOG/PHIS(IM,JM),PRF(IM,JM),PS(IM,JM)
    ! COMMON/DIAG/EP(IM,JM,KIM),WP(IM,JM,KM)
    ! COMMON/PHISS/PHI(IM,JM,KIM),PHIC(IM,JM,KIM)
    ! COMMON/GRAD/GPX(IM1,JM1,KM),GPY(IM1,JM1,KM)                   
    ! COMMON/AUX/WW(JM),Q(IM,JM,KM),ALFA(IM,JM,KM)
    !                                
    ! DIMENSION DIV(IM,JM,KM),EP1(IM,JM,KIM)
    !2345678-----------------------------------------------------
    A  =0.0
    A1  =0.0
    A2  =0.0
    AA  =0.0
    AX  =0.0
    AX1  =0.0
    AX2  =0.0
    AY  =0.0
    AY1  =0.0
    AY2  =0.0
    B  =0.0
    B1  =0.0
    B2  =0.0
    BB  =0.0
    CC  =0.0
    DD  =0.0
    EE  =0.0
    EPM  =0.0
    PIFX  =0.0
    PIFY  =0.0
    UMX  =0.0
    VMY  =0.0
    PS  =0.0
    DIV =0.0
    WP  =0.0
    DO I=1,IM
       DO J=1,JM
          PS(I,J)=PIF(I,J)*NS(I,J) + PT          
       END DO
    END DO

    !-------------------------------------------------------------
    !        PRIMEIRO  CALCULA-SE  A DIVERGENCIA
    !------------------------------------------------------------

    DO J=2,JM-1
       DO I=2,IM-1
          KS=KSE(I,J)
          DO K=1,KS-1
             AX1=(PIF(I,J)+PIF(I+1,J))*0.5
             AX2=(PIF(I,J)+PIF(I-1,J))*0.5
             AX=(AX1*UF(I,J,K)-AX2*UF(I-1,J,K))*RT/DX 
             AY1=(PIF(I,J)*MI(J)+PIF(I,J+1)*MI(J+1))*0.5
             AY2=(PIF(I,J)*MI(J)+PIF(I,J-1)*MI(J-1))*0.5
             AY=(AY1*VF(I,J,K)-AY2*VF(I,J-1,K))/DY
             DIV(I,J,K)=(AX+AY)*(EKI(K+1)-EKI(K))*M(J)*N

          END DO
       END DO
    END DO
    !--------------------------------------------------------------
    !                CALCULO DE LA TENDENCIA DE LA PRESAO P*
    !         P*(I,J)=(PS(I,J)-PT)/NS(I,J) ,SIRVE TAMBEM PARA (WP) 
    !-------------------------------------------------------------    
    DO I=2,IM-1
       DO J=2,JM-1
          KS=KSE(I,J)
          DO K=1,KS-1
             TPI(I,J) = TPI(I,J) - DIV(I,J,K)/NS(I,J)
          END DO
       END DO
    END DO

    !-----------------------------------------------------------
    !                CALCULO DE ETA PONTO EM K+1/2
    !-----------------------------------------------------------
    DO I=2,IM-1
       DO J=2,JM-1
          KS=KSE(I,J)
          DO  K=2,KS-1
             EP1(I,J,K)=EP1(I,J,K-1)-DIV(I,J,K-1)
          END DO
       END DO
    END DO

    DO I=2,IM-1
       DO J=2,JM-1
          KS=KSE(I,J)
          DO K=2,KS-1
             EE=-EKI(K)*TPI(I,J)
             EP(I,J,K) = (EE + EP1(I,J,K)) / PIF(I,J)
          END DO
       END DO
    END DO

    !-----------------------------------------------------------
    !           CALCULO DE VELOCIDADE VERTICAL OMEGA EM NIVEL K
    !----------------------------------------------------------
    !------------------------------------------------------------
    !               NIVELEIS PRINCIPAIS
    !               -------------------
    !PT=20.0 mb
    !
    !EK(1)=(100.-PT)/(1000.0-PT)
    !EK(2)=(200.-PT)/(1000.0-PT)
    !EK(3)=(300.-PT)/(1000.0-PT)
    !EK(4)=(500.-PT)/(1000.0-PT)
    !EK(5)=(700.-PT)/(1000.0-PT)
    !EK(6)=(850.-PT)/(1000.0-PT)
    !EK(7)=(960.-PT)/(1000.0-PT)
    !      PIF(I,J)=(PS(I,J)-PT)/NS(I,J)= 

    DO I=2,IM-1
       DO J=2,JM-1
          KS=KSE(I,J)
          DO K=1,KS-1
             UMX=(UF(I-1,J,K)+UF(I,J,K))*0.5
             VMY=(VF(I,J-1,K)+VF(I,J,K))*0.5
             PIFX=(PIF(I+1,J)-PIF(I-1,J))*EK(K)*UMX*RT*0.5/DX
             AA=PIF(I,J+1)*MI(J+1)-PIF(I,J-1)*MI(J-1)
             !            PIFY=AA*EK(K)*VMY*0.5/DY
             BB=PT*(MI(J+1)-MI(J-1))
             PIFY=(AA*EK(K)+BB)*VMY*0.5/DY              
             EPM=PIF(I,J)*(EP(I,J,K)+EP(I,J,K+1))*0.5
             WP(I,J,K)=EPM+EK(K)*TPI(I,J)+M(J)*N*(PIFX+PIFY)
          END DO
       END DO
    END DO

    !-----------------------------------------------------------
    !              completanmdo para ciclcioc  WP
    !-----------------------------------------------------------
    DO j=2,JM-1
       KS=KSE(1,J)
       DO  K=1,KS-1
          WP(1,J,K)= WP(IM-1,J,K)
       ENDDO
    ENDDO

    !--------------------------------------------------------
    !                    PROGNOSTICO DE PI
    !-----------------------------------------------------------

    NNT=NT
    IF(NT.LE.NT1)THEN                   
       !-----------------------------------------------------------

       DO J=2,JM-1
          DO  I=2,IM-1
             IF(MOD(NNT,2).EQ.0)THEN
                PIC(I,J)=PIF(I,J)-DT*TPI(I,J)
             ELSE
                PIC(I,J)=PIF(I,J)+DT*TPI(I,J)
             ENDIF
          END DO
       END DO

       !-----------------------------------------------------------
       !             CONTORNO
       !-----------------------------------------------------------
       DO J=1,JM
          DO I=1,IM
             PIC(1,J)=PIC(2,J)
             PIC(IM,J)=PIC(IM-1,J)
             !2345678----------------------------------------------------
             PIC(I,1)=PIC(I,2)
             PIC(I,JM)=PIC(I,JM1)
          END DO
       END DO

       !-----------------------------------------------------------
    ELSE     !lepafrog
       !------------------------------------------------------------
       DO J=2,JM1
          DO I=2,IM-1
             PIC(I,J)=PIA(I,J)+2.*DT*TPI(I,J)
          END DO
       END DO
       PRINT*,'fronteira este oeste'

       !***************************************************
       !           fronteira este oeste
       !****************************************************
       KCI1=KCI
       IF(KCI1.EQ.1)THEN

          ! ciclico    
          DO J=1,JM
             PIC(1,J)=PIC(IM-1,J)
             PIC(IM,J)=PIC(2,J)
          END DO
       PRINT*,'fronteira este oeste ciclico'

       ELSE  

       PRINT*,' R A D I A C I O N A L   L E S T E  O E S T E  ciclico'

          !-----------------------------------------------------------
          !             R A D I A C I O N A L   L E S T E  O E S T E
          !------------------------------------------------------------
          !               LESTE
          !---------------------------------------------------------------
          DO J=2,JM1
             B2=PIF(IM-1,J)+PIZ(IM-1,J)-2.*PIA(IM-1-1,J)
             IF(B2.NE.0.)THEN
                B1=-(PIF(IM-1,J)-PIZ(IM-1,J))
                B=B1/B2
                IF(B.LT.0.)B=0.
                IF(B.GT.1.)B=1.
             ELSE
                B=0.
             ENDIF
             AA=(1.-B)/(1.+B)
             BB=2.*B/(1.+B)
             PIC(IM,J)=AA*PIA(IM,J)+BB*PIF(IM-1,J)
          END DO

          !---------------------------------------------------------------
          !                OESTE
          !----------------------------------------------------------------
          DO J=2,JM1
             A2=PIF(2,J)+PIZ(2,J)-2.*PIA(3,J)
             IF(A2.NE.0.)THEN
                A1=PIF(2,J)-PIZ(2,J)
                A=A1/A2
                IF(A.GT.0.)A=0.
                IF(A.LT.-1.)A=-1.
             ELSE
                A=0.
             ENDIF
             CC=(1.+A)/(1.-A)
             DD=2.*A/(1.-A)
             PIC(1,J)=CC*PIA(1,J)-DD*PIF(2,J)
          END DO

          !------------------------------------------------------------
       ENDIF

   PRINT*,'FRONTEIRA NORTE E SUL',MAXVAL(PIC),MINVAL(PIC)

       !           FRONTEIRA NORTE E SUL
       !------------------------------------------------------------
       !             IF(KFRON-2)100,200,300 
       IF(KFRON-2<0)THEN
          !-----------------------------------------------------------
          !             NEWMAN
          !------------------------------------------------------------
          DO I=1,IM
             PIC(I,1)=PIC(I,2)
             PIC(I,JM)=PIC(I,JM1)
          END DO

       ELSE IF(KFRON-2==0)THEN
          !             GOTO 150
          !----------------------------------------------------------- 
          !              SPONJE
          !------------------------------------------------------------
          DO I=1,IM
             DO J=1,5
                PIC(I,J)=PIF(I,J)+WW(J)*TPI(I,J)
             END DO
          END DO

          DO I=1,IM
             DO J=JM-5,JM
                PIC(I,J)=PIF(I,J)+WW(J)*TPI(I,J)
             END DO
          END DO

       ELSE IF(KFRON-2>0)THEN
          !            GOTO 150
          !--------------------------------------------------------------
          !          RADIACIONAL NO SUL
          !------------------------------------------------------------                
          DO I=1,IM
             A2=PIF(I,2)+PIZ(I,2)-2.*PIA(I,3)
             IF(A2.NE.0.)THEN
                A1=PIF(I,2)-PIZ(I,2)
                A=A1/A2
                IF(A.GT.0.)A=0.
                IF(A.LT.-1.)A=-1.
             ELSE
                A=0.
             ENDIF
             CC=(1.+A)/(1.-A)
             DD=2.*A/(1.-A)
             PIC(I,1)=CC*PIA(I,1)-DD*PIF(I,2)
          END DO
          !-----------------------------------------------------------                 
          !          RADIACIONAL NORTE
          !----------------------------------------------------------
          DO I=1,IM
             B2=PIF(I,JM1)+PIZ(I,JM1)-2.*PIA(I,JM1-1)
             IF(B2.NE.0.)THEN
                B1=-(PIF(I,JM1)-PIZ(I,JM1))
                B=B1/B2
                IF(B.LT.0.)B=0.
                IF(B.GT.1.)B=1.
             ELSE
                B=0.
             ENDIF
             AA=(1.-B)/(1.+B)
             BB=2.*B/(1.+B)
             PIC(I,JM)=AA*PIA(I,JM)+BB*PIF(I,JM1)
          END DO

          !-------------------- ----------------------------------
       ENDIF
       !-------------------------------------------------------
       !         zeros para nao afectar as outras somas
       !-------------------------------------------------------
    ENDIF

    !150       CONTINUE

    DO I=1,IM
       DO J=1,JM
          TPI(I,J)=0.
          DO K=1,KIM                           
             EP1(I,J,K)=0.          
          END DO
       END DO
    END DO

    !-------------------------------------------

    RETURN
  END SUBROUTINE DIAG1

  !**********************************************************************
  SUBROUTINE DIAG2(&
       NT      ,&! INTEGER      , INTENT(IN   ) :: NT
       IM      ,&! INTEGER      , INTENT(IN   ) :: IM
       JM      ,&! INTEGER      , INTENT(IN   ) :: JM
       KM      ,&! INTEGER      , INTENT(IN   ) :: KM
       KIM     ,&! INTEGER      , INTENT(IN   ) :: KIM
       IM1     ,&! INTEGER      , INTENT(IN   ) :: IM1
       JM1     ,&! INTEGER      , INTENT(IN   ) :: JM1
       DX      ,&! REAL(KIND=r8), INTENT(IN   ) :: DX
       DY      ,&! REAL(KIND=r8), INTENT(IN   ) :: DY
       PT      ,&! REAL(KIND=r8), INTENT(IN   ) :: PT
       RS      ,&! REAL(KIND=r8), INTENT(IN   ) :: RS        !RS=287.
       RT      ,&! REAL(KIND=r8), INTENT(IN   ) :: RT        !RT=6.365E+6
       NS      ,&! REAL(KIND=r8), INTENT(IN   ) :: NS(IM,JM)
       PS      ,&! REAL(KIND=r8), INTENT(IN   ) :: PS(IM,JM)
       KSE     ,&! INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
       EK      ,&! REAL(KIND=r8), INTENT(IN   ) :: EK(KM)
       EKI     ,&! REAL(KIND=r8), INTENT(IN   ) :: EKI(KIM)
       TF      ,&! REAL(KIND=r8), INTENT(INOUT) :: TF(IM,JM,KM)
       PIF     ,&! REAL(KIND=r8), INTENT(IN   ) :: PIF(IM,JM)
       PHIS    ,&! REAL(KIND=r8), INTENT(IN   ) :: PHIS(IM,JM)
       ALFA    ,&! REAL(KIND=r8), INTENT(OUT  ) :: ALFA(IM,JM,KM)
       PHI     ,&! REAL(KIND=r8), INTENT(OUT  ) :: PHI(IM,JM,KIM)
       GPX     ,&! REAL(KIND=r8), INTENT(OUT  ) :: GPX(IM1,JM1,KM)
       GPY     ,&! REAL(KIND=r8), INTENT(OUT  ) :: GPY(IM1,JM1,KM)
       GF       &! REAL(KIND=r8), INTENT(OUT  ) :: GF(IM1,JM1,KM)
       )
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: NT
    INTEGER      , INTENT(IN   ) :: IM
    INTEGER      , INTENT(IN   ) :: JM
    INTEGER      , INTENT(IN   ) :: KM
    INTEGER      , INTENT(IN   ) :: KIM
    INTEGER      , INTENT(IN   ) :: IM1
    INTEGER      , INTENT(IN   ) :: JM1
    REAL(KIND=r8), INTENT(IN   ) :: DX
    REAL(KIND=r8), INTENT(IN   ) :: DY
    REAL(KIND=r8), INTENT(IN   ) :: PT
    REAL(KIND=r8), INTENT(IN   ) :: RS       !RS=287.
    REAL(KIND=r8), INTENT(IN   ) :: RT        !RT=6.365E+6
    REAL(KIND=r8), INTENT(IN   ) :: NS(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: PS(IM,JM)
    INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: EK(KM)
    REAL(KIND=r8), INTENT(IN   ) :: EKI(KIM)
    REAL(KIND=r8), INTENT(INOUT) :: TF(IM,JM,KM)
    REAL(KIND=r8), INTENT(IN   ) :: PIF(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: PHIS(IM,JM)
    REAL(KIND=r8), INTENT(OUT  ) :: ALFA(IM,JM,KM)
    REAL(KIND=r8), INTENT(OUT  ) :: PHI (IM,JM,KIM)
    REAL(KIND=r8), INTENT(OUT  ) :: GPX (IM1,JM1,KM)
    REAL(KIND=r8), INTENT(OUT  ) :: GPY (IM1,JM1,KM)
    REAL(KIND=r8), INTENT(OUT  ) :: GF  (IM1,JM1,KM)
    REAL(KIND=r8) :: ALFA1
    REAL(KIND=r8) :: DTX
    REAL(KIND=r8) :: EE1
    REAL(KIND=r8) :: EE2
    REAL(KIND=r8) :: EX1
    REAL(KIND=r8) :: EX2
    REAL(KIND=r8) :: EY1
    REAL(KIND=r8) :: EY2
    REAL(KIND=r8) :: FF1
    REAL(KIND=r8) :: FF2
    REAL(KIND=r8) :: FPM
    REAL(KIND=r8) :: FPMY
    REAL(KIND=r8) :: GPX1
    REAL(KIND=r8) :: GPY1
    REAL(KIND=r8) :: PHI10
    REAL(KIND=r8) :: PHI11
    REAL(KIND=r8) :: PHI2
    REAL(KIND=r8) :: PHI3
    REAL(KIND=r8) :: PHI4
    REAL(KIND=r8) :: PIF1
    REAL(KIND=r8) :: PHI1(IM,JM,KIM)
    INTEGER :: KSIG1
    INTEGER :: NT1
    INTEGER :: KS
    INTEGER :: i,j,k
    !--------------------------------------------------------------------
    !            PHI.......... ALTURA GEOPOTENCIAL
    !            GPX.......... GRADIENTE DE PRESAO EM X
    !            GPY.......... GRADIENTE DE PRESAO EM Y 
    !-------------------------------------------------------------------
    !2345678
    ! REAL  NS,N,M,MI
    ! DOUBLE PRECISION EP,WP,PHI,TPI,TU,TV,TT,PS                            
    ! DOUBLE PRECISION PIF,UF,VF,TF,PIA,UA,VA,TA,PIC,UC,VC,TC
    ! DOUBLE PRECISION GPX,GPY,GA,GF,GC
    ! DOUBLE PRECISION EKI,EK,PHIC

    ! DOUBLE PRECISION ALFA1,PIF1,PHI4,PHI3,PHI2,FF1,FF2,FPMY
    ! DOUBLE PRECISION GPY1,GPX1,FPM,EE1,EE2,EX1,EX2,EY1,EY2,PHI1
    ! PARAMETER(IM=146,JM=44,KM=7,IM1=IM,JM1=JM-1,KIM=KM+1,PT=20.) 
    ! PARAMETER(RT=6.365E+6,PIS=3.141596,RS=287.,DLAT=2.5,DLON=2.5)
    ! PARAMETER(RAD=PIS/180.,N=1./RT,DX=DLON*RAD,DY=DLAT*RAD)

    !COMMON/F/PIF(IM,JM),UF(IM1,JM1,KM),VF(IM1,JM1,KM),TF(IM,JM,KM)
    !COMMON/A/PIA(IM,JM),UA(IM1,JM1,KM),VA(IM1,JM1,KM),TA(IM,JM,KM)
    !COMMON/Z/PIZ(IM,JM),UZ(IM1,JM1,KM),VZ(IM1,JM1,KM),TZ(IM,JM,KM)

    !COMMON/C/PIC(IM,JM),UC(IM1,JM1,KM),VC(IM1,JM1,KM),TC(IM,JM,KM)
    !COMMON/G/GC(IM1,JM1,KM),GA(IM1,JM1,KM),GF(IM1,JM1,KM)
    !COMMON/TEN/TPI(IM,JM),TU(IM1,JM1,KM),TV(IM1,JM1,KM),TT(IM,JM,KM)    
    !COMMON/COORD/NS(IM,JM),M(JM),MI(JM),FLAT(JM),FLON(IM)
    !COMMON/BLOCOS/KSE(IM,JM),KSU(IM1,JM1),KSV(IM1,JM1)
    !COMMON/ETAS/EK(KM),EKI(KIM)
    !COMMON/TOPOG/PHIS(IM,JM),PRF(IM,JM),PS(IM,JM)
    !COMMON/DIAG/EP(IM,JM,KIM),WP(IM,JM,KM)
    !COMMON/PHISS/PHI(IM,JM,KIM),PHIC(IM,JM,KIM)
    !COMMON/GRAD/GPX(IM1,JM1,KM),GPY(IM1,JM1,KM)                   
    !COMMON/AUX/WW(JM),Q(IM,JM,KM),ALFA(IM,JM,KM)
    !                                                
    !DIMENSION PHI1(IM,JM,KIM)
    !2345678----------------------------------------------------
    !        CALCULO DE  GRADIENTE DE PRESAO EM I+1/2
    !-----------------------------------------------------------
    ALFA=0.0
    PHI =0.0
    GPX =0.0
    GPY =0.0
    GF=0.0  
    ALFA1=0.0
    DTX=0.0
    EE1=0.0
    EE2=0.0
    EX1=0.0
    EX2=0.0
    EY1=0.0
    EY2=0.0
    FF1=0.0
    FF2=0.0
    FPM=0.0
    FPMY=0.0
    GPX1=0.0
    GPY1=0.0
    PHI10=0.0
    PHI11=0.0
    PHI2=0.0
    PHI3=0.0
    PHI4=0.0
    PIF1=0.0
    PHI1=0.0
    DO I=1,IM-1
       DO J=1,JM-1
          KS=KSE(I,J)
          DO K=1,KS-1
             !----------------------------------------------
             !        GRADIENTE DE PRESAO EM X
             !----------------------------------------------
             EE1=1./NS(I+1,J)
             EE2=1./NS(I,J)    
             EX1=EE1+EE2
             EX2=EE1-EE2
             FPM=(PS(I,J)+PS(I+1,J)-2.*PT)*0.5
             GPX1=EX1*(PS(I+1,J)-PS(I,J))/(2.*DX)
             GPX(I,J,K)=(GPX1+FPM*EX2/DX)*EK(K)
             !----------------------------------------------
             !        GRADIENTE DE PRESAO EM Y
             !----------------------------------------------                                 
             EE1=1./NS(I,J+1)
             EE2=1./NS(I,J)    
             EY1=EE1+EE2
             EY2=EE1-EE2
             FPMY=(PS(I,J+1)+PS(I,J)-2.*PT)*0.5
             GPY1=EY1*(PS(I,J+1)-PS(I,J))/(2.*DY)
             GPY(I,J,K)=(GPY1+FPMY*EY2/DY)*EK(K)
          END DO
       END DO
    END DO
   PRINT*,'GRADIENTE DE PRESAO ',MAXVAL(GPY),MINVAL(GPY),MAXVAL(GPX),MINVAL(GPX)

    !----------------------------------------------------
    !       CALCLO DE ALTURA GEOPOTENCIALES
    !-------------------------------------------------------
    DO I=1,IM
       DO J=1,JM
          KS=KSE(I,J)
          DO K=KS-1,1,-1
             FF1=PIF(I,J)*(EKI(K+1)-EKI(K))
             FF2=PIF(I,J)*(EKI(K+1)+EKI(K))+2.*PT
             PHI1(I,J,K)=PHI1(I,J,K+1)+2.*RS*TF(I,J,K)*FF1/FF2
          END DO
       END DO
    END DO

    !----------------------------------------------------
    DO I=1,IM
       DO J=1,JM
          KS=KSE(I,J)
          DO K=KS,1,-1
             PHI(I,J,K)=PHI1(I,J,K)+PHIS(I,J)
          END DO
       END DO
    END DO

    !--------------------------------------------------------
    !        CALCULO DE ALFA
    !---------------------------------------------------------
    DO I=1,IM
       DO J=1,JM
          KS=KSE(I,J)
          DO K=1,KS-1
             ALFA(I,J,K)=RS*TF(I,J,K)/(PIF(I,J)*EK(K)+PT)
          END DO
       END DO
    END DO

    !---------------------------------------------------------
    !    CALCULO DE TEMPERATURA INICIAL  EM MREPOSO
    !                 PARA SIGMA CON TOPOGRAFIA
    !--------------------------------------------------------
    KSIG1=3
    NT1=NT
    IF(KSIG1.EQ.1.AND.NT1.EQ.1)THEN

       DO I=2,IM-1
          DO J=2,JM-1
             DO K=1,KM

                PHI10=(PHI(I,J,K+1)-PHI(I,J,K))

                PHI11=EKI(K+1)-EKI(K)
                DTX=(PHI10/PHI11)*(EK(K)*PIF(I,J)+PT)
                TF(I,J,K)=-DTX/PIF(I,J) 

             END DO
          END DO
       END DO

    ENDIF
    !----------------------------------------------------------------
    !                   CALCULO DE GRADIENTES PARA U
    !----------------------------------------------------------------
    DO I=2,IM1-1
       DO J=2,JM1-1
          KS=KSE(I,J)
          DO K=1,KS-1
             PHI4=(PHI(I,J,K)+PHI(I,J,K+1))*0.5
             PHI2=(PHI(I+1,J,K)+PHI(I+1,J,K+1))*0.5
             PHI3=(PHI2-PHI4)/DX
             !----------------------------------------------------------------
             ALFA1=(ALFA(I+1,J,K)+ALFA(I,J,K))*0.5
             PIF1=(PIF(I,J)+PIF(I+1,J))*0.5
             GF(I,J,K)=PIF1*(PHI3+ALFA1*GPX(I,J,K))*RT  
          END DO
       END DO
    END DO

    !------------------------------------------------------------------
    DO I=1,IM
       DO J=1,JM
          DO K=1,KIM
             PHI1(I,J,K)=0.
          END DO
       END DO
    END DO

    !------------------------------------------------------------ 

    RETURN
  END SUBROUTINE DIAG2

  !*****************************************************************
  !2345678

  SUBROUTINE TEMPE( &
       NT      ,&!INTEGER      , INTENT(IN ) :: NT
       IM      ,&!INTEGER      , INTENT(IN ) :: IM
       JM      ,&!INTEGER      , INTENT(IN ) :: JM
       KM      ,&!INTEGER      , INTENT(IN ) :: KM
       KIM     ,&!INTEGER      , INTENT(IN ) :: KIM
       IM1     ,&!INTEGER      , INTENT(IN ) :: IM1
       JM1     ,&!INTEGER      , INTENT(IN ) :: JM1
       KCI     ,&!INTEGER      , INTENT(IN ) :: KCI
       KFRON   ,&!INTEGER      , INTENT(IN ) :: KFRON
       KTOP    ,&!INTEGER      , INTENT(IN ) :: KTOP
       KQ      ,&!INTEGER      , INTENT(IN ) :: KQ
       DT      ,&!REAL(KIND=r8), INTENT(IN ) :: DT
       DX      ,&!REAL(KIND=r8), INTENT(IN ) :: DX
       DY      ,&!REAL(KIND=r8), INTENT(IN ) :: DY
       CP      ,&!REAL(KIND=r8), INTENT(IN ) :: CP
       TMAX    ,&!REAL(KIND=r8), INTENT(IN ) :: TMAX
       QMXY    ,&!REAL(KIND=r8), INTENT(IN ) :: QMXY  !QMXY=PIS/4
       QMZ     ,&!REAL(KIND=r8), INTENT(IN ) :: QMZ   !QMZ=2./PIS
       N       ,&!REAL(KIND=r8), INTENT(IN ) :: N     !N=1./RT
       RAD     ,&!REAL(KIND=r8), INTENT(IN ) :: RAD   !RAD=PIS/180.
       RT      ,&!REAL(KIND=r8), INTENT(IN ) :: RT    !RT=6.365E+6
       NT1     ,&!INTEGER      , INTENT(IN ) :: NT1
       NT2     ,&!INTEGER      , INTENT(IN ) :: NT2
       WW      ,&!REAL(KIND=r8), INTENT(IN ) :: WW(JM)
       FLAT    ,&!REAL(KIND=r8), INTENT(IN ) :: FLAT(JM)
       M       ,&!REAL(KIND=r8), INTENT(IN ) :: M(JM)
       MI      ,&!REAL(KIND=r8), INTENT(IN ) :: MI(JM)
       EKI     ,&!REAL(KIND=r8), INTENT(IN ) :: EKI(KIM)
       KSE     ,&!INTEGER      , INTENT(IN ) :: KSE(IM,JM)
       UF      ,&!REAL(KIND=r8), INTENT(IN ) :: UF(IM1,JM1,KM)
       VF      ,&!REAL(KIND=r8), INTENT(IN ) :: VF(IM1,JM1,KM)
       TF      ,&!REAL(KIND=r8), INTENT(IN ) :: TF(IM,JM,KM)
       Q       ,&!REAL(KIND=r8), INTENT(IN ) :: Q(IM,JM,KM)
       PIF     ,&!REAL(KIND=r8), INTENT(IN ) :: PIF(IM,JM)
       EP      ,&!REAL(KIND=r8), INTENT(IN ) :: EP(IM,JM,KIM)
       WP      ,&!REAL(KIND=r8), INTENT(IN ) :: WP(IM,JM,KM)
       ALFA    ,&!REAL(KIND=r8), INTENT(IN ) :: ALFA(IM,JM,KM) 
       PIC     ,&!REAL(KIND=r8), INTENT(IN ) :: PIC(IM,JM)
       PIA     ,&!REAL(KIND=r8), INTENT(IN ) :: PIA(IM,JM)
       TA      ,&!REAL(KIND=r8), INTENT(IN ) :: TA(IM,JM,KM)
       TZ      ,&!REAL(KIND=r8), INTENT(IN ) :: TZ(IM,JM,KM)
       TC      ,&!REAL(KIND=r8), INTENT(INOUT) :: TC(IM,JM,KM)
       TT       &!REAL(KIND=r8), INTENT(INOUT) :: TT(IM,JM,KM)  
       )
    INTEGER      , INTENT(IN   ) :: NT
    INTEGER      , INTENT(IN   ) :: IM
    INTEGER      , INTENT(IN   ) :: JM
    INTEGER      , INTENT(IN   ) :: KM
    INTEGER      , INTENT(IN   ) :: KIM
    INTEGER      , INTENT(IN   ) :: IM1
    INTEGER      , INTENT(IN   ) :: JM1
    INTEGER      , INTENT(IN   ) :: KCI
    INTEGER      , INTENT(IN   ) :: KFRON
    INTEGER      , INTENT(IN   ) :: KTOP
    INTEGER      , INTENT(IN   ) :: KQ
    REAL(KIND=r8), INTENT(IN   ) :: DT
    REAL(KIND=r8), INTENT(IN   ) :: DX
    REAL(KIND=r8), INTENT(IN   ) :: DY
    REAL(KIND=r8), INTENT(IN   ) :: CP
    REAL(KIND=r8), INTENT(IN   ) :: TMAX
    REAL(KIND=r8), INTENT(IN   ) :: QMXY  !QMXY=PIS/4
    REAL(KIND=r8), INTENT(IN   ) :: QMZ   !QMZ=2./PIS
    REAL(KIND=r8), INTENT(IN   ) :: N     !N=1./RT
    REAL(KIND=r8), INTENT(IN   ) :: RAD   !RAD=PIS/180.
    REAL(KIND=r8), INTENT(IN   ) :: RT    !RT=6.365E+6
    INTEGER      , INTENT(IN   ) :: NT1
    INTEGER      , INTENT(IN   ) :: NT2
    REAL(KIND=r8), INTENT(IN   ) :: WW(JM)
    REAL(KIND=r8), INTENT(IN   ) :: FLAT(JM)
    REAL(KIND=r8), INTENT(IN   ) :: M(JM)
    REAL(KIND=r8), INTENT(IN   ) :: MI(JM)
    REAL(KIND=r8), INTENT(IN   ) :: EKI(KIM)
    INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: UF(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: VF(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: TF(IM,JM,KM)
    REAL(KIND=r8), INTENT(IN   ) :: Q(IM,JM,KM)
    REAL(KIND=r8), INTENT(IN   ) :: PIF(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: EP(IM,JM,KIM)
    REAL(KIND=r8), INTENT(IN   ) :: WP(IM,JM,KM)
    REAL(KIND=r8), INTENT(IN   ) :: ALFA(IM,JM,KM) 
    REAL(KIND=r8), INTENT(IN   ) :: PIC(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: PIA(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: TA(IM,JM,KM)
    REAL(KIND=r8), INTENT(IN   ) :: TZ(IM,JM,KM)
    REAL(KIND=r8), INTENT(INOUT) :: TC(IM,JM,KM)
    REAL(KIND=r8), INTENT(INOUT) :: TT(IM,JM,KM)  


    REAL(KIND=r8) :: QQ  (IM,JM,KM)
    REAL(KIND=r8) :: A
    REAL(KIND=r8) :: A1
    REAL(KIND=r8) :: A2
    REAL(KIND=r8) :: A4
    REAL(KIND=r8) :: AA
    REAL(KIND=r8) :: B
    REAL(KIND=r8) :: B1
    REAL(KIND=r8) :: B2
    REAL(KIND=r8) :: BB
    REAL(KIND=r8) :: ADTX
    REAL(KIND=r8) :: ADTY
    REAL(KIND=r8) :: ADTZ
    REAL(KIND=r8) :: ANG
    REAL(KIND=r8) :: AT1
    REAL(KIND=r8) :: AT2
    REAL(KIND=r8) :: BT1
    REAL(KIND=r8) :: BT2
    REAL(KIND=r8) :: BT3
    REAL(KIND=r8) :: BTY1
    REAL(KIND=r8) :: BTY2
    REAL(KIND=r8) :: CC
    REAL(KIND=r8) :: CT1
    REAL(KIND=r8) :: CT2
    REAL(KIND=r8) :: DD
    REAL(KIND=r8) :: DD1
    REAL(KIND=r8) :: DENO
    REAL(KIND=r8) :: Q0
    REAL(KIND=r8) :: QMT
    REAL(KIND=r8) :: QT
    REAL(KIND=r8) :: TADI
    REAL(KIND=r8) :: TEMP
    REAL(KIND=r8) :: WQ
    INTEGER       :: KCI1
    INTEGER       :: KK1
    INTEGER       :: KK2
    INTEGER       :: KK3
    INTEGER       :: klate
    INTEGER       :: KS
    INTEGER       :: KSX
    INTEGER       :: KSY
    INTEGER       :: NNT
    INTEGER       :: NT8
    INTEGER       :: NT88
    INTEGER       :: NT99
    INTEGER       :: NTEE
    INTEGER       :: i,j,k
    ! REAL  NS,N,M,MI
    !DOUBLE PRECISION EP,WP,PHI,TPI,TU,TV,TT,PS
    !DOUBLE PRECISION PIF,UF,VF,TF,PIA,UA,VA,TA,PIC,UC,VC,TC
    !DOUBLE PRECISION GPX,GPY,GA,GF,GC
    !DOUBLE PRECISION EKI,EK,PHIC
    ! 
    ! DOUBLE PRECISION AA,BB,CC,DD,A1,A2,A,B1,B2,B,DENO,CT1,CT2
    ! DOUBLE PRECISION ADTZ,BT1,BT2,BTY1,BT3,ADTY,AT1,AT2,ADTX,WQ
    ! PARAMETER(IM=146,JM=44,KM=7,IM1=IM,JM1=JM-1,KIM=KM+1,PT=20.) 
    ! PARAMETER(IM2=IM-2,JM2=JM-2)
    ! PARAMETER(RT=6.365E+6,PIS=3.141596,RS=287.,RCP=0.286,CP=1004.)
    ! PARAMETER(DLAT=2.5,DLON=2.5)
    ! PARAMETER(RAD=PIS/180.,N=1./RT,DX=DLON*RAD,DY=DLAT*RAD)  

    !PARAMETER(QMXY=PIS/4.,QMZ=2./PIS)
    !COMMON/F/PIF(IM,JM),UF(IM1,JM1,KM),VF(IM1,JM1,KM),TF(IM,JM,KM)
    !COMMON/A/PIA(IM,JM),UA(IM1,JM1,KM),VA(IM1,JM1,KM),TA(IM,JM,KM)
    !COMMON/Z/PIZ(IM,JM),UZ(IM1,JM1,KM),VZ(IM1,JM1,KM),TZ(IM,JM,KM)

    !COMMON/C/PIC(IM,JM),UC(IM1,JM1,KM),VC(IM1,JM1,KM),TC(IM,JM,KM)
    !COMMON/G/GC(IM1,JM1,KM),GA(IM1,JM1,KM),GF(IM1,JM1,KM)
    !COMMON/TEN/TPI(IM,JM),TU(IM1,JM1,KM),TV(IM1,JM1,KM),TT(IM,JM,KM)    
    !COMMON/COORD/NS(IM,JM),M(JM),MI(JM),FLAT(JM),FLON(IM)
    !COMMON/BLOCOS/KSE(IM,JM),KSU(IM1,JM1),KSV(IM1,JM1)
    !COMMON/ETAS/EK(KM),EKI(KIM)
    !COMMON/TOPOG/PHIS(IM,JM),PRF(IM,JM),PS(IM,JM)
    !COMMON/DIAG/EP(IM,JM,KIM),WP(IM,JM,KM)
    !COMMON/PHISS/PHI(IM,JM,KIM),PHIC(IM,JM,KIM)
    !COMMON/GRAD/GPX(IM1,JM1,KM),GPY(IM1,JM1,KM)                   
    !COMMON/AUX/WW(JM),Q(IM,JM,KM),ALFA(IM,JM,KM) 
    !      
    !DIMENSION QQ(IM,JM,KM),QXYZT(IM2,JM2,KM)   
    !DIMENSION QTP(120)
    !2345678
    QQ=0.0
    A=0.0
    A1=0.0
    A2=0.0
    A4=0.0
    AA=0.0
    B=0.0
    B1=0.0
    B2=0.0
    BB=0.0
    ADTX=0.0
    ADTY=0.0
    ADTZ=0.0
    ANG=0.0
    AT1=0.0
    AT2=0.0
    BT1=0.0
    BT2=0.0
    BT3=0.0
    BTY1=0.0
    BTY2=0.0
    CC=0.0
    CT1=0.0
    CT2=0.0
    DD=0.0
    DD1=0.0
    DENO=0.0
    Q0=0.0
    QMT=0.0
    QT=0.0
    TADI=0.0
    TEMP=0.0
    WQ=0.0
    KCI1=0
    KK1=0
    KK2=0
    KK3=0
    klate=0
    KS=0
    KSX=0
    KSY=0
    NNT=0
    NT8=0
    NT88=0
    NT99=0
    NTEE=0
    IF(NT.GE.NT2)THEN
       !**************************************************************
       !                CUBICA INVERTIDA
       !                ----------------
       KK1=KQ
       IF(KK1.EQ.1)THEN

          NT88=288+NT2-1
          NT99=NT88+1
          IF(NT.LE.NT88)THEN

             TADI=3.-0.125
             TEMP=TADI-(NT-NT2)*DT/28800.
             A4=8./3.                      
             QT=0.5*(A4**3)*(TEMP**2)*EXP(-A4*TEMP)

             QMT=1.12/3                  !para 48h
             DD1=QMXY*QMZ*QMT*24.*3600
             Q0=TMAX*CP/DD1          
             IF(NT.GT.288-30)QT=0.2169      
             QT=Q0*QT


          ELSE         

             IF(NT.EQ.NT99)NT8=0
             NT8=NT8+1  

             TADI=3.-0.125

             TEMP=TADI-(NT8-1)*DT/28800.                                                                 
             A4=8./3.                     
             QT=0.5*(A4**3)*(TEMP**2)*EXP(-A4*TEMP)
             QMT=1.12/3         
             DD1=QMXY*QMZ*QMT*24.*3600
             Q0=TMAX*CP/DD1 
             IF(NT8.LE.96)QT=0.2169
             IF(NT8.GE.258)QT=0.2169
             QT=Q0*QT

             IF(NT8.EQ.288)NT8=0
          ENDIF

       ENDIF

       !************************************************************
       !                estacionaria con pico 17horas
       !-----------------------------------------------------------

       !                 KK2=KQ
       KK2=2
       IF(KK2.EQ.2)THEN

          NTEE=160+NT2-1
          IF(NT.LE.NTEE)THEN

             TADI=3.-0.125
             TEMP=TADI-(NT-NT2)*DT/28800.
             A4=8./3.                      
             QT=0.5*(A4**3)*(TEMP**2)*EXP(-A4*TEMP)

             QMT=1.12/3                  !para 48h
             DD1=QMXY*QMZ*QMT*24.*3600
             Q0=TMAX*CP/DD1 
             QT=Q0*QT

          ELSE

             QT=0.373

          ENDIF
       ENDIF
       !**************************************************************

       !            PRINT QT
       !            --------
       !
       !             IF(NT.LE.1441)THEN
       !                       
       !                  NT0=NT-1
       !                  IF(MOD(NT0,12).EQ.0.)THEN 
       !                  II=II+1
       !                  QTP(II)=QT*0.5*24*3600./1004.
       !                        ENDIF
       !              ENDIF
       !
       !2345678------------------------------------------------------
       !                   if(nt.eq.1442)then
       !        OPEN(39,FILE='QT.dat',STATUS='UNKNOWN')
       !
       !
       !          DO 129 I=1,120
       !          WRITE(39,111)I,QTP(I)
       !111       FORMAT(1X,I3,3X,F6.1)
       !129       CONTINUE
       !                   ENDIF

       !***************************************************************
       !         CUBICA  COM 12 HORAS DE PICO a4=4/3
       !----------------------------------------------------------- 
       KK3=KQ 
       IF(KK3.EQ.3)THEN
          TEMP=(NT-NT2)*DT/28800.
          A4=4/3.
          QT=0.5*(A4**3)*(TEMP**2)*EXP(-A4*TEMP)
          QMT=1./9.
          DD1=QMXY*QMZ*QMT*24.*3600
          Q0=TMAX*CP/DD1
          QT=Q0*QT 
       ENDIF
       !--------------------------------------------------------------------         
       DO I=1,IM
          DO J=1,JM
             DO K=1,KM
                QQ(I,J,K)=Q(I,J,K)*QT             
             END DO
          END DO
       END DO

       !*********************************************************************

    ENDIF
    !           PRINT  QQQQQQQQQQQQQQQQQQ
    !------------------------------------
    !           DO 391 I=2,IM-1
    !           DO 391 J=2,JM-1
    !           DO 391 K=KM,1,-1
    !          QXYZT(I-1,J-1,KIM-K)=Q0*24.*3600*Q(I,J,K)*QMT/CP                       
    !391       CONTINUE
    !2345678------------------------------------------------------
    !       OPEN(40,FILE='calor.dat',STATUS='UNKNOWN',
    !     * form='unformatted',access='direct',recl=IM2*JM2*KM)                               
    !       WRITE(40,rec=1)QXYZT
    !     
    !
    !************************************************************************
    !        AVECCION HORIXONTAL
    !------------------------------------------------------------------
    !2345678           EN X
    DO I=2,IM-1
       DO J=2,JM-1
          KS=KSE(I,J)
          DO K=1,KS-1
             AT1=(PIF(I,J)*TF(I,J,K)+PIF(I+1,J)*TF(I+1,J,K))*UF(I,J,K)*0.5
             AT2=(PIF(I,J)*TF(I,J,K)+PIF(I-1,J)*TF(I-1,J,K))*UF(I-1,J,K)*0.5
             ADTX=(AT1-AT2)*RT/DX
             !--------------------------------------------------------------------------
             !            EN Y
             BT1=PIF(I,J)*TF(I,J,K)*MI(J)
             BT2=PIF(I,J+1)*TF(I,J+1,K)*MI(J+1)
             BTY1=(BT1+BT2)*0.5
             BT3=PIF(I,J-1)*TF(I,J-1,K)*MI(J-1)
             BTY2=(BT1+BT3)*0.5
             ADTY=(BTY1*VF(I,J,K)-BTY2*VF(I,J-1,K))/DY
             !--------------------------------------------------------------------------
             !            EN Z
             DENO=2.*M(J)*N*(EKI(K+1)-EKI(K))
             IF(K.NE.KS-1)THEN
                CT1=EP(I,J,K+1)*(TF(I,J,K)+TF(I,J,K+1))
             ELSE
                CT1=0.
             ENDIF

             IF(K.NE.1)THEN
                CT2=EP(I,J,K)*(TF(I,J,K)+TF(I,J,K-1))
             ELSE
                CT2=0.
             ENDIF
             !------------------------------------------------------------------------
             ADTZ=PIF(I,J)*(CT1-CT2)/DENO
             !------------------------------------------------------------------------
             !             OMEGA EM T
             !           WT=ALFA(I,J,K)*WP(I,J,K)*PIF(I,J)*RT*MI(J)/CP
             !             CALOR LATENTE
             !               QQ(i,j,k)=0.
             !-------------------------------
             !      cuando uno qiere calor latente
             !---------------------------------
             klate=1440    !5dias depois debe progar ondas
             IF(NT .GE. klate) QQ(i,j,k)=0.

             WQ=PIF(I,J)*(ALFA(I,J,K)*WP(I,J,K)+QQ(I,J,K))*MI(J)*RT/CP
             !-------------------------------------------------------------------------
             ANG=FLAT(J)*RAD


             !------------------------------------------------------------------------

             !               TERMO DIFUSION
             !--------------------------------------------------------------------
             !
             !              DTX=(TA(I+1,J,K)-2.*TA(I,J,K)+TA(I-1,J,K))/(DX**2)
             !
             !              DTY=(TA(J+1,J,K)-2.*TA(I,J,K)+TA(J-1,J,K))/(DY**2)
             !
             !              DTY2=TAN(ANG)*(TA(J+1,J,K)-TA(J-1,J,K))/(2.*DY)
             !              
             !              FT=PIF(I,J)*CDF*(DTX*M(J)**2+DTY*N**2-DTY2*N**2)
             !-------------------------------------------------------------------

             !           TT(I,J,K)=(-ADTX-ADTY-ADTZ+WQ)*M(J)*N+FT
             TT(I,J,K)=(-ADTX-ADTY-ADTZ+WQ)*M(J)*N
          END DO
       END DO
    END DO

    !------------------------------------------------------------------------
    !             PRONOSTICO DA TEMPERATURA
    !----------------------------------------------------------------------

    NNT=NT
    IF(NT .LE. NT1)THEN

       DO J=2,JM-1
          DO I=2,IM-1
             KS=KSE(I,J)
             DO K=1,KS-1
                IF(MOD(NNT,2).EQ.0)THEN
                   TC(I,J,K)=(TF(I,J,K)*PIF(I,J)-DT*TT(I,J,K))/PIC(I,J) 
                ELSE
                   TC(I,J,K)=(TF(I,J,K)*PIF(I,J)+DT*TT(I,J,K))/PIC(I,J) 
                ENDIF
             END DO
          END DO
       END DO

       !---------------------------------------------------------------------
       !              FRONTERA NEWMAN NORTE/SUR E OESTE E ESTE
       DO K=1,KM
          DO J=1,JM
             DO I=1,IM
                TC(1,J,K)=TC(2,J,K)
                TC(IM,J,K)=TC(IM-1,J,K)
                TC(I,1,K)=TC(I,2,K)
                TC(I,JM,K)=TC(I,JM1,K)
             END DO
          END DO
       END DO

       !----------------------------------------------------------------------
    ELSE
       !--------------LEAPFROG------------------------------------------------
       DO J=2,JM-1
          DO I=2,IM-1
             KS=KSE(I,J)
             DO  K=1,KS-1
                TC(I,J,K)=(TA(I,J,K)*PIA(I,J)+2.*DT*TT(I,J,K))/PIC(I,J)
             END DO
          END DO
       END DO

       !****************************************************************
       !             frter seste oeste
       !********************************************************************
       KCI1=KCI
       IF(KCI.EQ.1)THEN

          !           CICLICO
          DO K=1,KM
             DO J=1,JM           
                TC(1,J,K)=TC(IM-1,J,K)
                TC(IM,J,K)=TC(2,J,K)
             END DO
          END DO

       ELSE
          !------------------------------------------------------------------
          !         R A D I A C I O N A L     L E S T E   O E S T E
          !-----------------------------------------------------------------
          !                  LESTE
          !-----------------------------------------------------------------
          DO K=1,KM
             DO J=2,JM1
                B2=TF(IM-1,J,K)+TZ(IM-1,J,K)-2.*TA(IM-1-1,J,K)
                IF(B2.NE.0.)THEN
                   B1=-(TF(IM-1,J,K)-TZ(IM-1,J,K))
                   B=B1/B2
                   IF(B.LT.0.)B=0.
                   IF(B.GT.1.)B=1.
                ELSE
                   B=0.
                ENDIF
                AA=(1.-B)/(1.+B)
                BB=2.*B/(1.+B)
                TC(IM,J,K)=AA*TA(IM,J,K)+BB*TF(IM-1,J,K)
             END DO
          END DO

          !-------------------------------------------------------------
          !               OESTE
          !------------------------------------------------------------
          DO K=1,KM
             DO J=2,JM1
                A2=TF(2,J,K)+TZ(2,J,K)-2.*TA(3,J,K)
                IF(A2.NE.0.)THEN
                   A1=TF(2,J,K)-TZ(2,J,K)
                   A=A1/A2
                   IF(A.GT.0.)A=0.
                   IF(A.LT.-1.)A=-1.
                ELSE
                   A=0.
                ENDIF
                CC=(1.+A)/(1.-A)
                DD=2.*A/(1.-A)
                TC(1,J,K)=CC*TA(1,J,K)-DD*TF(2,J,K)
             END DO
          END DO

       ENDIF


       !-------------------------------------------------------------
       !         NORTE E SUL
       !--------------------------------------------------------------
       !IF(KFRON-2)100,200,300
       IF(KFRON-2 < 0)THEN
          !------------------------------------------------------------
          !              NEWMAN
          !------------------------------------------------------------
          DO K=1,KM
             DO I=1,IM
                DO J=1,JM
                   TC(I,1,K)=TC(I,2,K)
                   TC(I,JM,K)=TC(I,JM1,K)
                END DO
             END DO
          END DO

       ELSE IF(KFRON-2 == 0)THEN
          !-----------------------------------------------------------
          !             SPONJA
          !------------------------------------------------------------
          DO K=1,KM
             DO I=1,IM
                DO J=1,5
                   TC(I,J,K)=TF(I,J,K)+WW(J)*TT(I,J,K)
                END DO
             END DO
          END DO

          DO K=1,KM
             DO I=1,IM
                DO J=JM-5,JM
                   TC(I,J,K)=TF(I,J,K)+WW(J)*TT(I,J,K)
                END DO
             END DO
          END DO

       ELSE IF(KFRON-2 > 0)THEN
          !------------------------------------------------------------
          !            FRONTERA SUR RADIACIONAL
          !-------------------------------------------------------------
          DO K=1,KM
             DO I=1,IM
                A2=TF(I,2,K)+TZ(I,2,K)-2.*TA(I,3,K)
                IF(A2.NE.0.)THEN
                   A1=TF(I,2,K)-TZ(I,2,K)
                   A=A1/A2
                   IF(A.GT.0.)A=0.
                   IF(A.LT.-1.)A=-1.
                ELSE
                   A=0.
                ENDIF
                CC=(1.+A)/(1.-A)
                DD=2.*A/(1.-A)
                TC(I,1,K)=CC*TA(I,1,K)-DD*TF(I,2,K)
             END DO
          END DO

          !------------------------------------------------------------
          !           FRONTERA NORTE RADIACIONAL
          !-------------------------------------------------------------
          DO K=1,KM
             DO I=1,IM

                B2=TF(I,JM1,K)+TZ(I,JM1,K)-2.*TA(I,JM1-1,K)

                IF(B2.NE.0.)THEN
                   B1=-(TF(I,JM1,K)-TZ(I,JM1,K))
                   B=B1/B2
                   IF(B.LT.0.)B=0.
                   IF(B.GT.1.)B=1.
                ELSE
                   B=0.
                ENDIF

                AA=(1.-B)/(1.+B)
                BB=2.*B/(1.+B)    
                TC(I,JM,K)=AA*TA(I,JM,K)+BB*TF(I,JM1,K)

             END DO
          END DO

          !-----------------------------------------------------------
       END IF
    ENDIF

    !150                              ENDIF
    !---------------------------------------------------------------
    !               variacao temporal da forcante termica
    !------------------ ---------------------------------------------
    !----------------------------------------------------------
    !---------------------------------------------------------
    !            TEMPERATURAS DENTRO DA TOPOGRAFIA
    !---------------------------------------------------------
    IF(KTOP.EQ.1)THEN

       DO I=2,IM-1
          DO J=2,JM1
             KS=KSE(I,J)
             KSX=KSE(I+1,J)
             KSY=KSE(I,J+1)
             IF((KSX-KS).LT.0) TC(I+1,J,KS-1)=TC(I,J,KS-1)
             IF((KSX-KS).GT.0) TC(I,J,KS)=TC(I+1,J,KS)
             IF((KSY-KS).LT.0) TC(I,J+1,KS-1)=TC(I,J,KS-1)
             IF((KSY-KS).GT.0) TC(I,J,KS)=TC(I,J+1,KS)
          END DO
       END DO

    ENDIF

    !------------------------------------------------------------

    RETURN
  END SUBROUTINE TEMPE

  !
  !*************************************************************
  SUBROUTINE CAMP(&
       NT        ,&!INTEGER     , INTENT(IN   ) :: NT
       IM        ,&!INTEGER     , INTENT(IN   ) :: IM
       JM        ,&!INTEGER     , INTENT(IN   ) :: JM
       KM        ,&!INTEGER     , INTENT(IN   ) :: KM
       KIM       ,&!INTEGER     , INTENT(IN   ) :: KIM
       IM1       ,&!INTEGER     , INTENT(IN   ) :: IM1
       JM1       ,&!INTEGER     , INTENT(IN   ) :: JM1
       PT        ,&!REAL(KIND=r8), INTENT(IN   ) :: PT
       RS        ,&!REAL(KIND=r8), INTENT(IN   ) :: RS
       RT        ,&!REAL(KIND=r8), INTENT(IN   ) :: RT
       DX        ,&!REAL(KIND=r8), INTENT(IN   ) :: DX
       PS        ,&!REAL(KIND=r8), INTENT(OUT  ) :: PS(IM,JM)
       NS        ,&!REAL(KIND=r8), INTENT(IN   ) :: NS(IM,JM)
       PIC       ,&!REAL(KIND=r8), INTENT(IN   ) :: PIC(IM,JM)
       GPXC      ,&!REAL(KIND=r8), INTENT(OUT  ) :: GPXC(IM1,JM1,KM)
       EK        ,&!REAL(KIND=r8), INTENT(IN   ) :: EK(KM)
       EKI       ,&!REAL(KIND=r8), INTENT(IN   ) :: EKI(KIM)
       PHIS      ,&!REAL(KIND=r8), INTENT(IN   ) :: PHIS(IM,JM)
       PHIC1     ,&!REAL(KIND=r8), INTENT(INOUT) :: PHIC1(IM,JM,KIM)
       PHIC      ,&!REAL(KIND=r8), INTENT(INOUT) :: PHIC(IM,JM,KIM)
       TC        ,&!REAL(KIND=r8), INTENT(IN   ) :: TC(IM,JM,KM)
       ALFAC     ,&!REAL(KIND=r8), INTENT(INOUT) :: ALFAC(IM,JM,KM)
       GA        ,&!REAL(KIND=r8), INTENT(IN   ) :: GA(IM1,JM1,KM)
       GC        ,&!REAL(KIND=r8), INTENT(INOUT) :: GC(IM1,JM1,KM)
       GF        ,&!REAL(KIND=r8), INTENT(INOUT) :: GF(IM1,JM1,KM)
       KSE        &!INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
       )
    INTEGER      , INTENT(IN   ) :: NT
    INTEGER      , INTENT(IN   ) :: IM
    INTEGER      , INTENT(IN   ) :: JM
    INTEGER      , INTENT(IN   ) :: KM
    INTEGER      , INTENT(IN   ) :: KIM
    INTEGER      , INTENT(IN   ) :: IM1
    INTEGER      , INTENT(IN   ) :: JM1
    REAL(KIND=r8), INTENT(IN   ) :: PT
    REAL(KIND=r8), INTENT(IN   ) :: RS
    REAL(KIND=r8), INTENT(IN   ) :: RT
    REAL(KIND=r8), INTENT(IN   ) :: DX
    REAL(KIND=r8), INTENT(OUT  ) :: PS(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: NS(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: PIC(IM,JM)
    REAL(KIND=r8), INTENT(OUT  ) :: GPXC(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: EK(KM)
    REAL(KIND=r8), INTENT(IN   ) :: EKI(KIM)
    REAL(KIND=r8), INTENT(IN   ) :: PHIS(IM,JM)
    REAL(KIND=r8), INTENT(INOUT) :: PHIC1(IM,JM,KIM)
    REAL(KIND=r8), INTENT(INOUT) :: PHIC(IM,JM,KIM)
    REAL(KIND=r8), INTENT(IN   ) :: TC(IM,JM,KM)
    REAL(KIND=r8), INTENT(INOUT) :: ALFAC(IM,JM,KM)
    REAL(KIND=r8), INTENT(IN   ) :: GA(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(INOUT) :: GC(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(INOUT) :: GF(IM1,JM1,KM)
    INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
    REAL(KIND=r8) :: ALFA1
    REAL(KIND=r8) :: EE1
    REAL(KIND=r8) :: EE2
    REAL(KIND=r8) :: EX1
    REAL(KIND=r8) :: EX2
    REAL(KIND=r8) :: FF1
    REAL(KIND=r8) :: FF2
    REAL(KIND=r8) :: FMI
    REAL(KIND=r8) :: FMI1
    REAL(KIND=r8) :: FPM
    REAL(KIND=r8) :: GPX1
    REAL(KIND=r8) :: PHI2
    REAL(KIND=r8) :: PHI3
    REAL(KIND=r8) :: PHI4
    REAL(KIND=r8) :: PIF1
    INTEGER :: KS
    INTEGER :: i,j,k
    !REAL  NS,N,M,MI
    !DOUBLE PRECISION EP,WP,PHI,TPI,TU,TV,TT,PS                           
    !DOUBLE PRECISION PIF,UF,VF,TF,PIA,UA,VA,TA,PIC,UC,VC,TC
    !DOUBLE PRECISION GPX,GPY,GA,GF,GC
    !DOUBLE PRECISION EKI,EK,PHIC

    !DOUBLE PRECISION ALFA1,PIF1,GPXC,PHIC1,ALFAC
    !DOUBLE PRECISION EE1,EE2,EX2,FPM,GPX1,FF1,FF2,PHI4,PHI2,PHI3
    !PARAMETER(IM=146,JM=44,KM=7,IM1=IM,JM1=JM-1,KIM=KM+1,PT=20.) 
    !PARAMETER(RT=6.365E+6,PIS=3.141596,RS=287.,DLON=2.5,DLAT=2.5)
    !PARAMETER(RAD=PIS/180.,N=1./RT,DX=DLON*RAD,DY=DLAT*RAD)      
    !COMMON/F/PIF(IM,JM),UF(IM1,JM1,KM),VF(IM1,JM1,KM),TF(IM,JM,KM)
    !COMMON/A/PIA(IM,JM),UA(IM1,JM1,KM),VA(IM1,JM1,KM),TA(IM,JM,KM)
    !COMMON/Z/PIZ(IM,JM),UZ(IM1,JM1,KM),VZ(IM1,JM1,KM),TZ(IM,JM,KM)

    !COMMON/C/PIC(IM,JM),UC(IM1,JM1,KM),VC(IM1,JM1,KM),TC(IM,JM,KM)
    !COMMON/G/GC(IM1,JM1,KM),GA(IM1,JM1,KM),GF(IM1,JM1,KM)
    !COMMON/TEN/TPI(IM,JM),TU(IM1,JM1,KM),TV(IM1,JM1,KM),TT(IM,JM,KM)    
    !COMMON/COORD/NS(IM,JM),M(JM),MI(JM),FLAT(JM),FLON(IM)
    !COMMON/BLOCOS/KSE(IM,JM),KSU(IM1,JM1),KSV(IM1,JM1)
    !COMMON/ETAS/EK(KM),EKI(KIM)
    !COMMON/TOPOG/PHIS(IM,JM),PRF(IM,JM),PS(IM,JM)
    !COMMON/DIAG/EP(IM,JM,KIM),WP(IM,JM,KM)
    !COMMON/PHISS/PHI(IM,JM,KIM),PHIC(IM,JM,KIM)
    !COMMON/GRAD/GPX(IM1,JM1,KM),GPY(IM1,JM1,KM)                   
    !COMMON/AUX/WW(JM),Q(IM,JM,KM),ALFA(IM,JM,KM)

    !DIMENSION GPXC(IM1,JM1,KM),PHIC1(IM,JM,KIM),ALFAC(IM,JM,KM)
    !DIMENSION GPYC(IM1,JM1,KM)     
    !2345678-----------------------------------------------------------------
    PS=0.0
    GPXC=0.0
    ALFA1=0.0
    EE1=0.0
    EE2=0.0
    EX1=0.0
    EX2=0.0
    FF1=0.0
    FF2=0.0
    FMI=0.0
    FMI1=0.0
    FPM=0.0
    GPX1=0.0
    PHI2=0.0
    PHI3=0.0
    PHI4=0.0
    PIF1=0.0
    DO I=1,IM
       DO J=1,JM
          PS(I,J)=PIC(I,J)*NS(I,J)+PT
       END DO
    END DO

    !              CALCULO DE  GRADIENTE DE PRESAO EM I+1/2
    DO I=1,IM-1
       DO J=1,JM-1
          KS=KSE(I,J)
          DO K=1,KS-1
             !----------------------------------------------
             !                      GRADIENTE DE PRSAO EM X
             !----------------------------------------------
             EE1=1./NS(I+1,J)
             EE2=1./NS(I,J)    
             EX1=EE1+EE2
             EX2=EE1-EE2
             FPM=(PS(I,J)+PS(I+1,J)-2.*PT)*0.5
             GPX1=EX1*(PS(I+1,J)-PS(I,J))/(2.*DX)
             GPXC(I,J,K)=(GPX1+FPM*EX2/DX)*EK(K)
          END DO
       END DO
    END DO

    !----------------------------------------------------

    !----------------------------------------------------
    !       CALCLO DE ALTURAS GEOPOTENCIALES
    !-------------------------------------------------------
    DO I=1,IM
       DO J=1,JM
          KS=KSE(I,J)
          DO K=KS-1,1,-1
             FF1=PIC(I,J)*(EKI(K+1)-EKI(K))
             FF2=PIC(I,J)*(EKI(K+1)+EKI(K))+2.*PT
             PHIC1(I,J,K)=PHIC1(I,J,K+1)+2.*RS*TC(I,J,K)*FF1/FF2
          END DO
       END DO
    END DO

    !----------------------------------------------------
    DO I=1,IM
       DO J=1,JM
          KS=KSE(I,J)
          DO K=KS,1,-1
             PHIC(I,J,K)=PHIC1(I,J,K)+PHIS(I,J)
          END DO
       END DO
    END DO

    !----------------------------------------------------
    DO I=1,IM
       DO J=1,JM
          DO K=1,KIM
             PHIC1(I,J,K)=0.                             
          END DO
       END DO
    END DO

    !----------------------------------------------------------------
    !                   CALCULO DE ALFA
    DO I=1,IM
       DO J=1,JM
          KS=KSE(I,J)
          DO K=1,KS-1
             ALFAC(I,J,K)=RS*TC(I,J,K)/(PIC(I,J)*EK(K)+PT)
          END DO
       END DO
    END DO

    !----------------------------------------------------------------
    !                   CALCULO DE GRADIENTES PARA U
    !----------------------------------------------------------------
    DO I=2,IM1-1
       DO J=2,JM1-1
          KS=KSE(I,J)
          DO K=1,KS-1
             PHI4=(PHIC(I,J,K)+PHIC(I,J,K+1))*0.5
             PHI2=(PHIC(I+1,J,K)+PHIC(I+1,J,K+1))*0.5
             PHI3=(PHI2-PHI4)/DX
             !----------------------------------------------------------------
             ALFA1=(ALFAC(I+1,J,K)+ALFAC(I,J,K))*0.5
             PIF1=(PIC(I,J)+PIC(I+1,J))*0.5
             GC(I,J,K)=PIF1*(PHI3+ALFA1*GPXC(I,J,K))*RT  
          END DO
       END DO
    END DO

    !----------------------------------------------------------------

    !---------------------------------------------------------------
    !                  MEDIAS PONDERADAS
    !----------------------------------------------------------------
    IF(NT.NE.1)THEN
       DO I=2,IM1-1
          DO J=2,JM1-1
             KS=KSE(I,J)
             DO K=1,KS-1

                FMI=0.2          !maximo e 0.25
                FMI1=1.-2.*FMI

                GF(I,J,K)=FMI1*GF(I,J,K)+FMI*(GC(I,J,K)+GA(I,J,K))
             END DO
          END DO
       END DO

    ENDIF
    !---------------------------------------------------------------
    RETURN
  END SUBROUTINE CAMP

  !2345678****************************************************************
  SUBROUTINE TVENTO(&
       IM       ,&! INTEGER      , INTENT(IN   ) :: IM
       JM       ,&! INTEGER      , INTENT(IN   ) :: JM
       KM       ,&! INTEGER      , INTENT(IN   ) :: KM
       KIM      ,&! INTEGER      , INTENT(IN   ) :: KIM
       IM1      ,&! INTEGER      , INTENT(IN   ) :: IM1
       JM1      ,&! INTEGER      , INTENT(IN   ) :: JM1
       DX       ,&! REAL(KIND=r8), INTENT(IN   ) :: DX
       DY       ,&! REAL(KIND=r8), INTENT(IN   ) :: DY
       GRA      ,&! REAL(KIND=r8), INTENT(IN   ) :: GRA
       OMEGA    ,&! REAL(KIND=r8), INTENT(IN   ) :: OMEGA
       RAD      ,&! REAL(KIND=r8), INTENT(IN   ) :: RAD
       RT       ,&! REAL(KIND=r8), INTENT(IN   ) :: RT
       N        ,&! REAL(KIND=r8), INTENT(IN   ) :: N
       KSE      ,&! INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
       UF       ,&! REAL(KIND=r8), INTENT(IN   ) :: UF(IM1,JM1,KM)
       VF       ,&! REAL(KIND=r8), INTENT(IN   ) :: VF(IM1,JM1,KM)
       PIF      ,&! REAL(KIND=r8), INTENT(IN   ) :: PIF(IM,JM)
       FLAT     ,&! REAL(KIND=r8), INTENT(IN   ) :: FLAT(JM)
       M        ,&! REAL(KIND=r8), INTENT(IN   ) :: M(JM)
       MI       ,&! REAL(KIND=r8), INTENT(IN   ) :: MI(JM)
       EK       ,&! REAL(KIND=r8), INTENT(IN   ) :: EK(KM)
       EKI      ,&! REAL(KIND=r8), INTENT(IN   ) :: EKI(KIM)
       EP       ,&! REAL(KIND=r8), INTENT(IN   ) :: EP(IM,JM,KIM)
       PHIS     ,&! REAL(KIND=r8), INTENT(IN   ) :: PHIS(IM,JM)
       ALFA     ,&! REAL(KIND=r8), INTENT(IN   ) :: ALFA(IM,JM,KM)
       GPY      ,&! REAL(KIND=r8), INTENT(IN   ) :: GPY(IM1,JM1,KM) 
       PHI      ,&! REAL(KIND=r8), INTENT(IN   ) :: PHI(IM,JM,KIM)
       TU       ,&! REAL(KIND=r8), INTENT(OUT  ) :: TU(IM1,JM1,KM)
       TV       ,&! REAL(KIND=r8), INTENT(OUT  ) :: TV(IM1,JM1,KM)
       GF        &! REAL(KIND=r8), INTENT(IN   ) :: GF(IM1,JM1,KM)
       )
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: IM
    INTEGER      , INTENT(IN   ) :: JM
    INTEGER      , INTENT(IN   ) :: KM
    INTEGER      , INTENT(IN   ) :: KIM
    INTEGER      , INTENT(IN   ) :: IM1
    INTEGER      , INTENT(IN   ) :: JM1
    REAL(KIND=r8), INTENT(IN   ) :: DX
    REAL(KIND=r8), INTENT(IN   ) :: DY
    REAL(KIND=r8), INTENT(IN   ) :: GRA
    REAL(KIND=r8), INTENT(IN   ) :: OMEGA
    REAL(KIND=r8), INTENT(IN   ) :: RAD
    REAL(KIND=r8), INTENT(IN   ) :: RT
    REAL(KIND=r8), INTENT(IN   ) :: N
    INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: UF(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: VF(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: PIF(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: FLAT(JM)
    REAL(KIND=r8), INTENT(IN   ) :: M(JM)
    REAL(KIND=r8), INTENT(IN   ) :: MI(JM)
    REAL(KIND=r8), INTENT(IN   ) :: EK(KM)
    REAL(KIND=r8), INTENT(IN   ) :: EKI(KIM)
    REAL(KIND=r8), INTENT(IN   ) :: EP(IM,JM,KIM)
    REAL(KIND=r8), INTENT(IN   ) :: PHIS(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: ALFA(IM,JM,KM)
    REAL(KIND=r8), INTENT(IN   ) :: GPY(IM1,JM1,KM) 
    REAL(KIND=r8), INTENT(IN   ) :: PHI(IM,JM,KIM)
    REAL(KIND=r8), INTENT(OUT  ) :: TU(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(OUT  ) :: TV(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: GF(IM1,JM1,KM)

    REAL(KIND=r8) :: ADUX
    REAL(KIND=r8) :: ADUY
    REAL(KIND=r8) :: ADUZ
    REAL(KIND=r8) :: ADVX
    REAL(KIND=r8) :: ADVY
    REAL(KIND=r8) :: ADVZ
    REAL(KIND=r8) :: ADZ1
    REAL(KIND=r8) :: ADZ2
    REAL(KIND=r8) :: ALF
    REAL(KIND=r8) :: ANG
    REAL(KIND=r8) :: AP1
    REAL(KIND=r8) :: AP2
    REAL(KIND=r8) :: AVZ1
    REAL(KIND=r8) :: AVZ2
    REAL(KIND=r8) :: BG1
    REAL(KIND=r8) :: BG2
    REAL(KIND=r8) :: BP1
    REAL(KIND=r8) :: BP2
    REAL(KIND=r8) :: BP3
    REAL(KIND=r8) :: BPY1
    REAL(KIND=r8) :: BPY2
    REAL(KIND=r8) :: CD
    REAL(KIND=r8) :: COU
    REAL(KIND=r8) :: COV
    REAL(KIND=r8) :: COV1
    REAL(KIND=r8) :: CV1
    REAL(KIND=r8) :: CV2
    REAL(KIND=r8) :: CV3
    REAL(KIND=r8) :: CVX1
    REAL(KIND=r8) :: CVX2
    REAL(KIND=r8) :: DEN1
    REAL(KIND=r8) :: DEN2
    REAL(KIND=r8) :: DEN3
    REAL(KIND=r8) :: DV1
    REAL(KIND=r8) :: DV2
    REAL(KIND=r8) :: FF
    REAL(KIND=r8) :: FU
    REAL(KIND=r8) :: FV
    REAL(KIND=r8) :: GGV
    REAL(KIND=r8) :: GPGU
    REAL(KIND=r8) :: GPGV
    REAL(KIND=r8) :: GPV
    REAL(KIND=r8) :: PHISU
    REAL(KIND=r8) :: PHISV
    REAL(KIND=r8) :: UFM1
    REAL(KIND=r8) :: UFM2
    REAL(KIND=r8) :: UK1
    REAL(KIND=r8) :: UK3
    REAL(KIND=r8) :: VK1
    REAL(KIND=r8) :: UMU
    REAL(KIND=r8) :: UUM
    REAL(KIND=r8) :: VFM1
    REAL(KIND=r8) :: VFM2
    INTEGER       :: KS
    INTEGER       :: i,j,k
    ! REAL  NS,N,M,MI
    ! DOUBLE PRECISION EP,WP,PHI,TPI,TU,TV,TT,PS                            
    ! DOUBLE PRECISION PIF,UF,VF,TF,PIA,UA,VA,TA,PIC,UC,VC,TC
    ! DOUBLE PRECISION GPX,GPY,GA,GF,GC
    ! DOUBLE PRECISION EKI,EK,PHIC

    ! 
    ! DOUBLE PRECISION ADUX,ADUY,ADUZ,COU,GPGU
    ! DOUBLE PRECISION ADVX,ADVY,ADVZ,COV,GPGV,GPV,BP1
    ! DOUBLE PRECISION CV1,CV2,CVX1,CV3,CVX2,UFM1,UFM2,DV1,DV2,VK1
    ! DOUBLE PRECISION AVZ1,AVZ2,UUM,COV1,BG1,BG2,GGV,AP1,AP2,BP2
    ! DOUBLE PRECISION BPY1,BP3,BPY2,VFM2,VFM1,ADZ2,ADZ1,FF
    ! DOUBLE PRECISION DEN1,DEN2,UK3,FU,FV

    !PARAMETER(OMEGA=7.292E-5,GRA=9.8)
    !PARAMETER(IM=146,JM=44,KM=7,IM1=IM,JM1=JM-1,KIM=KM+1,PT=20.) 
    !PARAMETER(RT=6.365E+6,PIS=3.141596,RS=287.,DLON=2.5,DLAT=2.5)
    !PARAMETER(RAD=PIS/180.,N=1./RT,DX=DLON*RAD,DY=DLAT*RAD)      
    !COMMON/F/PIF(IM,JM),UF(IM1,JM1,KM),VF(IM1,JM1,KM),TF(IM,JM,KM)
    !COMMON/A/PIA(IM,JM),UA(IM1,JM1,KM),VA(IM1,JM1,KM),TA(IM,JM,KM)
    !COMMON/Z/PIZ(IM,JM),UZ(IM1,JM1,KM),VZ(IM1,JM1,KM),TZ(IM,JM,KM)

    !COMMON/C/PIC(IM,JM),UC(IM1,JM1,KM),VC(IM1,JM1,KM),TC(IM,JM,KM)
    !COMMON/G/GC(IM1,JM1,KM),GA(IM1,JM1,KM),GF(IM1,JM1,KM)
    !COMMON/TEN/TPI(IM,JM),TU(IM1,JM1,KM),TV(IM1,JM1,KM),TT(IM,JM,KM)    
    !COMMON/COORD/NS(IM,JM),M(JM),MI(JM),FLAT(JM),FLON(IM)
    !COMMON/BLOCOS/KSE(IM,JM),KSU(IM1,JM1),KSV(IM1,JM1)
    !COMMON/ETAS/EK(KM),EKI(KIM)
    !COMMON/TOPOG/PHIS(IM,JM),PRF(IM,JM),PS(IM,JM)
    !COMMON/DIAG/EP(IM,JM,KIM),WP(IM,JM,KM)
    !COMMON/PHISS/PHI(IM,JM,KIM),PHIC(IM,JM,KIM)
    !COMMON/GRAD/GPX(IM1,JM1,KM),GPY(IM1,JM1,KM)                   
    !COMMON/AUX/WW(JM),Q(IM,JM,KM),ALFA(IM,JM,KM)


    !2345678----------------------------
    !----               ADVECCAO EM v
    !--------------------------------------    
    !          EM X
    TU=0.0
    TV=0.0
    ADUX=0.0
    ADUY=0.0
    ADUZ=0.0
    ADVX=0.0
    ADVY=0.0
    ADVZ=0.0
    ADZ1=0.0
    ADZ2=0.0
    ALF=0.0
    ANG=0.0
    AP1=0.0
    AP2=0.0
    AVZ1=0.0
    AVZ2=0.0
    BG1=0.0
    BG2=0.0
    BP1=0.0
    BP2=0.0
    BP3=0.0
    BPY1=0.0
    BPY2=0.0
    CD=0.0
    COU=0.0
    COV=0.0
    COV1=0.0
    CV1=0.0
    CV2=0.0
    CV3=0.0
    CVX1=0.0
    CVX2=0.0
    DEN1=0.0
    DEN2=0.0
    DEN3=0.0
    DV1=0.0
    DV2=0.0
    FF=0.0
    FU=0.0
    FV=0.0
    GGV=0.0
    GPGU=0.0
    GPGV=0.0
    GPV=0.0
    PHISU=0.0
    PHISV=0.0
    UFM1=0.0
    UFM2=0.0
    UK1=0.0
    UK3=0.0
    VK1=0.0
    UMU=0.0
    UUM=0.0
    VFM1=0.0
    VFM2=0.0
    DO I=2,IM1-1
       DO J=2,JM1-1
          KS=KSE(I,J)
          DO K=1,KS-1
             !                  DX=DLON*RAD
             !                  ADVECCAO EM X
             CV1  = VF(I  ,J,K)*(PIF(I  ,J) + PIF(I  ,J+1))*0.5
             CV2  = VF(I+1,J,K)*(PIF(I+1,J) + PIF(I+1,J+1))*0.5
             CVX1 = (CV1+CV2)*0.5
 
             CV3  = VF(I-1,J,K)*(PIF(I-1,J)+PIF(I-1,J+1))*0.5
             CVX2 = (CV1+CV3)*0.5
 
             UFM1=(UF(I-1,J,K) + UF(I-1,J+1,K))*0.5
             UFM2=(UF(I  ,J,K) + UF(I  ,J+1,K))*0.5
 
             ADVX=(CVX1*UFM2 - CVX2*UFM1)*RT/DX
             !-----------------------------------------------
             !                  ADVECCAO EM Y
             !                  DY=DLAT*RAD
             !2345678----------------------------------------------
             DV1=(VF(I,J,K)+VF(I,J+1,K))*0.5
             DV2=(VF(I,J,K)+VF(I,J-1,K))*0.5
             ADVY=(PIF(I,J+1)*MI(J+1)*DV1**2-PIF(I,J)*MI(J)*DV2**2)/DY
             !----------------------------------------------------------
             !                   ADVECCAO EM Z
             !----------------------------------------------------------

             !EKI             NIVEIS INTERMEDIOS(K+1/2)  DETERMINACAO DE ETAS
             !N=1./RT
             !M(J)= 1/(RT*COS(ANG))
	     !EP =  CALCULO DE ETA PONTO EM K+1/2

             DEN2=M(J)*N*(EKI(K+1)-EKI(K))
             VK1=(PIF(I,J)+PIF(I,J+1))*0.5
             IF(K.NE.KS-1)THEN
                AVZ1=EP(I,J,K+1)*(VF(I,J,K)+VF(I,J,K+1))*0.5
             ELSE
                AVZ1=0.
             ENDIF
             IF(K.NE.1)THEN
                AVZ2=EP(I,J,K)*(VF(I,J,K)+VF(I,J,K-1))*0.5
             ELSE
                AVZ2=0.
             ENDIF
             ADVZ=VK1*(AVZ1-AVZ2)/DEN2
             !-----------------------------------------------------------------

             !                       TERMO CORIOLIS
             !2345678-----------------------------------------------------------
             ANG=FLAT(J)*RAD
             UUM = (UF(I,J,K)+UF(I-1,J,K)+UF(I,J+1,K)+UF(I-1,J+1,K))*0.25 
             !
             COV1 = 2.*OMEGA*(MI(J+1)+MI(J))*0.5+UUM
             COV  = COV1*(PIF(I,J)+PIF(I,J+1))*UUM*SIN(ANG)*RT*0.5
             !-------------------------------------------------------------
             !                TERMO DIFUSION
             !----------------------------------------------------------------
             !              DVX=(VA(I+1,J,K)-2.*VA(I,J,K)+VA(I-1,J,K))/(DX**2)
             !
             !              DVY=(VA(J+1,J,K)-2.*VA(I,J,K)+VA(J-1,J,K))/(DY**2)
             !
             !              DVY2=TAN(ANG)*(VA(J+1,J,K)-VA(J-1,J,K))/(2.*DY)
             !              
             !                PXM=(PIF(I,J)+PIF(I,J+1))*0.5
             !              
             !              FY=PXM*CDF*(DVX*M(J)**2+DVY*N**2-DVY2*N**2)
             !*****************************************************************
             IF(K.EQ.KS-1)THEN
                !           FRICCION IN THE PLANETARY  BOUNDARY LAYER
                !           -----------------------------------------
                !            
                PHISV=(PHIS(I,J)+PHIS(I,J+1))*0.5
                CD=(1.+0.0001*PHISV)*1.0E-3

                ALF=(ALFA(I,J,KS-1)+ALFA(I,J+1,KS-1))*0.5
                DEN3=EK(KS-1)-EK(KS-2)
                FV=-2.*GRA*CD*VF(I,J,KS-1)*ABS(VF(I,J,KS-1))/(ALF*DEN3)

             ENDIF

             !****************************************************************

             !                      GRADIENTE DE PRESAO
             !-------------------------------------------------------------
             GPV = (ALFA(I,J,K) + ALFA(I,J+1,K))*GPY(I,J,K)*0.5
             !-----------------------------------------------------------
             !                      GRADIENTE DE GEOPOTENCIAL
             !2345678-----------------------------------------------------------
             BG1=(PHI(I,J,K+1)+PHI(I,J,K))*0.5
             BG2=(PHI(I,J+1,K+1)+PHI(I,J+1,K))*0.5
             GGV=(BG2-BG1)/DY
             GPGV=(PIF(I,J)*MI(J)+PIF(I,J+1)*MI(J+1))*(GPV+GGV)*0.5
             !           TV(I,J,K)=-(ADVX+ADVY+ADVZ+COV+GPGV)*M(J)*N+FY+FV
             TV(I,J,K)=-(ADVX+ADVY+ADVZ+COV+GPGV)*M(J)*N + FV

             !                      write(*,*)TV(i,j,7),FV
          END DO
       END DO
    END DO


    !----------------------------------------------------------------
    !                            TENDENDCIA DE U
    !----------------------------------------------------------------
    ! 2345678
    DO I=2,IM1-1
       DO J=2,JM1-1
          KS=KSE(I,J)
          DO K=1,KS-1
             !----------------------------------------------------------------
             !            ADVECAO DE U EM X
             !----------------------------------------------------------------
             AP1=(UF(I,J,K)+UF(I+1,J,K))*0.5
             AP2=(UF(I,J,K)+UF(I-1,J,K))*0.5
             !          
             ADUX=(PIF(I+1,J)*AP1**2-PIF(I,J)*AP2**2)/(N*DX)
             !----------------------------------------------------------------
             !           ADVECAO DE U EM Y
             !----------------------------------------------------------------
             BP1=UF(I,J,K)*MI(J)*(PIF(I,J)+PIF(I+1,J))*0.5
             BP2=UF(I,J+1,K)*MI(J+1)*(PIF(I,J+1)+PIF(I+1,J+1))*0.5
             BPY1=(BP1+BP2)*0.5
             BP3=UF(I,J-1,K)*MI(J-1)*(PIF(I,J-1)+PIF(I+1,J-1))*0.5
             BPY2=(BP1+BP3)*0.5
             VFM1=(VF(I,J,K)+VF(I+1,J,K))*0.5
             VFM2=(VF(I,J-1,K)+VF(I+1,J-1,K))*0.5
             ADUY=(BPY1*VFM1-BPY2*VFM2)/DY
             !---------------------------------------------------------------
             !          ADVECCAO DE U EM Z
             !----------------------------------------------------------------
             DEN1=M(J)*N*(EKI(K+1)-EKI(K))
             UK1=(PIF(I,J)+PIF(I+1,J))*0.5
             IF(K.NE.KS-1)THEN
                ADZ1=EP(I,J,K+1)*(UF(I,J,K)+UF(I,J,K+1))*0.5
             ELSE
                ADZ1=0.
             ENDIF
             IF(K.NE.1)THEN
                UK3=(PIF(I,J)+PIF(I+1,J))*0.5
                ADZ2=EP(I,J,K)*(UF(I,J,K)+UF(I,J,K-1))*0.5
             ELSE
                ADZ2=0.
             ENDIF
             ADUZ=UK1*(ADZ1-ADZ2)/DEN1
             !--------------------------------------------------------------------
             !           TERMO CORIOLIS
             !2345678-------------------------------------------------------------
             ANG=FLAT(J)*RAD
             UMU=(PIF(I,J)+PIF(I+1,J))*0.5
             FF=(2.*OMEGA*MI(J)+UF(I,J,K))*SIN(ANG)*UMU*RT
             COU=-FF*(VF(I,J,K)+VF(I,J-1,K)+VF(I+1,J,K)+VF(I+1,J-1,K))*0.25
             !--------------------------------------------------------------------
             !                TERMO DIFUSION
             !--------------------------------------------------------------------
             !              DUX=(UA(I+1,J,K)-2.*UA(I,J,K)+UA(I-1,J,K))/(DX**2)
             !
             !              DUY=(UA(J+1,J,K)-2.*UA(I,J,K)+UA(J-1,J,K))/(DY**2)
             !
             !              DUY2=TAN(ANG)*(UA(J+1,J,K)-UA(J-1,J,K))/(2.*DY)
             !              
             !              FX=UMU*CDF*(DUX*M(J)**2+DUY*N**2-DUY2*N**2)
             !--------------------------------------------------------------
             !********************************************************************
             IF(K.EQ.KS-1)THEN
                !           FRICCION IN THE PLANETARY  BOUNDARY LAYER
                !           -----------------------------------------
                !            
                PHISU=(PHIS(I,J)+PHIS(I+1,J))*0.5
                CD=(1.+0.0001*PHISU)*1.0E-3
                ALF=(ALFA(I,J,KS-1)+ALFA(I+1,J,KS-1))*0.5
                DEN3=EK(KS-1)-EK(KS-2)
                FU=-2.*GRA*CD*UF(I,J,KS-1)*ABS(UF(I,J,KS-1))/(ALF*DEN3)


             ENDIF

             !******************************************************************

             !           CAMPANA PARA GRADIENTE DE PRESAO E GEOPOTENCIAL
             !--------------------------------------------------------------------
             GPGU=GF(I,J,K)
             !        TU(I,J,K)=-(ADUX+ADUY+ADUZ+COU+GPGU)*M(J)*N+FX+FU    
             TU(I,J,K)=-(ADUX+ADUY+ADUZ+COU+GPGU)*M(J)*N+FU
          END DO
       END DO
    END DO

    !--------------------------------------------------------------------
    !                  RETURN
  END SUBROUTINE TVENTO


  !********************************************************************
  SUBROUTINE EULER(&
       IM       , &!INTEGER      , INTENT(IN   ) :: IM
       JM       , &!INTEGER      , INTENT(IN   ) :: JM
       KM       , &!INTEGER      , INTENT(IN   ) :: KM
       IM1      , &!INTEGER      , INTENT(IN   ) :: IM1
       JM1      , &!INTEGER      , INTENT(IN   ) :: JM1
       NT       , &!INTEGER      , INTENT(IN   ) :: NT
       DT       , &!REAL(KIND=r8), INTENT(IN   ) :: DT
       KSU      , &!INTEGER      , INTENT(IN   ) :: KSU(IM1,JM1)
       KSV      , &!INTEGER      , INTENT(IN   ) :: KSV(IM1,JM1)
       UC       , &!REAL(KIND=r8), INTENT(INOUT) :: UC(IM1,JM1,KM)
       VC       , &!REAL(KIND=r8), INTENT(INOUT) :: VC(IM1,JM1,KM)
       UF       , &!REAL(KIND=r8), INTENT(IN   ) :: UF(IM1,JM1,KM)
       VF       , &!REAL(KIND=r8), INTENT(IN   ) :: VF(IM1,JM1,KM)
       TU       , &!REAL(KIND=r8), INTENT(IN   ) :: TU(IM1,JM1,KM)
       TV       , &!REAL(KIND=r8), INTENT(IN   ) :: TV(IM1,JM1,KM)
       PIC      , &!REAL(KIND=r8), INTENT(IN   ) :: PIC(IM,JM)
       PIF        )!REAL(KIND=r8), INTENT(IN   ) :: PIF(IM,JM)

    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: IM
    INTEGER      , INTENT(IN   ) :: JM
    INTEGER      , INTENT(IN   ) :: KM
    INTEGER      , INTENT(IN   ) :: IM1
    INTEGER      , INTENT(IN   ) :: JM1
    INTEGER      , INTENT(IN   ) :: NT
    REAL(KIND=r8), INTENT(IN   ) :: DT
    INTEGER      , INTENT(IN   ) :: KSU(IM1,JM1)
    INTEGER      , INTENT(IN   ) :: KSV(IM1,JM1)
    REAL(KIND=r8), INTENT(INOUT) :: UC(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(INOUT) :: VC(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: UF(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: VF(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: TU(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: TV(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: PIC(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: PIF(IM,JM)


    REAL(KIND=r8) :: PIFM
    REAL(KIND=r8) :: PIFMJ
    REAL(KIND=r8) :: PTU
    REAL(KIND=r8) :: PTV
    INTEGER :: KS
    INTEGER :: NTT
    INTEGER :: i,j,k
    !PARAMETER(IM=146,JM=44,KM=7,IM1=IM,JM1=JM-1,KIM=KM+1,PT=20.)
    !!2345678 
    !REAL  NS,N,M,MI
    !DOUBLE PRECISION EP,WP,PHI,TPI,TU,!,TT,PS
    !DOUBLE PRECISION PIF,UF,VF,TF,PIA,UA,VA,TA,PIC,UC,VC,TC
    !DOUBLE PRECISION GPX,GPY,GA,GF,GC
    !DOUBLE PRECISION EKI,EK,PHIC
    !DOUBLE PRECISION PTU,PIFM,PTV,PIFMJ 
    !COMMON/F/PIF(IM,JM),UF(IM1,JM1,KM),VF(IM1,JM1,KM),TF(IM,JM,KM)
    !COMMON/A/PIA(IM,JM),UA(IM1,JM1,KM),VA(IM1,JM1,KM),TA(IM,JM,KM)
    !COMMON/Z/PIZ(IM,JM),UZ(IM1,JM1,KM),VZ(IM1,JM1,KM),TZ(IM,JM,KM)

    !COMMON/C/PIC(IM,JM),UC(IM1,JM1,KM),VC(IM1,JM1,KM),TC(IM,JM,KM)
    !COMMON/G/GC(IM1,JM1,KM),GA(IM1,JM1,KM),GF(IM1,JM1,KM)
    !COMMON/TEN/TPI(IM,JM),TU(IM1,JM1,KM),TV(IM1,JM1,KM),TT(IM,JM,KM)
    !COMMON/COORD/NS(IM,JM),M(JM),MI(JM),FLAT(JM),FLON(IM)
    !COMMON/BLOCOS/KSE(IM,JM),KSU(IM1,JM1),KSV(IM1,JM1)
    !COMMON/ETAS/EK(KM),EKI(KIM)
    !COMMON/TOPOG/PHIS(IM,JM),PRF(IM,JM),PS(IM,JM)
    !COMMON/DIAG/EP(IM,JM,KIM),WP(IM,JM,KM)
    !COMMON/PHISS/PHI(IM,JM,KIM),PHIC(IM,JM,KIM)
    !COMMON/GRAD/GPX(IM1,JM1,KM),GPY(IM1,JM1,KM)
    !COMMON/AUX/WW(JM),Q(IM,JM,KM),ALFA(IM,JM,KM)

    NTT=NT   
    !------------------------------------------------------------------------
    DO J=2,JM1-1
       DO I=2,IM1-1
          KS=KSU(I,J)
          DO K=1,KS-1
             PTU=(PIC(I,J)+PIC(I+1,J))*0.5
             PIFM=(PIF(I,J)+PIF(I+1,J))*0.5

             IF(MOD(NTT,2).EQ.0)THEN              
                UC(I,J,K)=(UF(I,J,K)*PIFM-DT*TU(I,J,K))/PTU
             ELSE
                UC(I,J,K)=(UF(I,J,K)*PIFM+DT*TU(I,J,K))/PTU
             ENDIF
          END DO
       END DO
    END DO
    !-----------------------------------------------------------------------
    DO J=2,JM1-1
       DO I=2,IM1-1
          KS=KSV(I,J)
          DO K=1,KS-1
             PTV=(PIC(I,J)+PIC(I,J+1))*0.5
             PIFMJ=(PIF(I,J)+PIF(I,J+1))*0.5
             IF(MOD(NTT,2).EQ.0)THEN 
                VC(I,J,K)=(VF(I,J,K)*PIFMJ-DT*TV(I,J,K))/PTV
             ELSE   

                !            VC(I,J,K)=(VF(I,J,K)*PIFMJ+DT*TV(I,J,K))/PTV
                VC(I,J,K)=0.        !SIN TOPOGRAFIA
             ENDIF
             !         
          END DO
       END DO
    END DO
    !-----------------------------------------------------------------------
    !               CONTORNO OESTE LESTE NORTE E SUL
    !---------------------------------------------------------------------
    DO K=1,KM
       DO J=1,JM1
          DO I=1,IM1
             UC(1,J,K)=UC(2,J,K)
             UC(IM1,J,K)=UC(IM1-1,J,K)
             UC(I,1,K)=UC(I,2,K)
             UC(I,JM1,K)=UC(I,JM1-1,K)
             !------------------------------------------------------------------- 
             VC(1,J,K)=VC(2,J,K)
             VC(IM1,J,K)=VC(IM1-1,J,K)
             VC(I,1,K)=0.
             VC(I,JM1,K)=0.
          END DO
       END DO
    END DO
    !---------------------------------------------------------------------
    RETURN
  END SUBROUTINE EULER

  !**********************************************************************
  !2345678
  SUBROUTINE LEAP( &
       IM          ,&!INTEGER      , INTENT(IN   ) :: IM
       JM          ,&!INTEGER      , INTENT(IN   ) :: JM
       KM          ,&!INTEGER      , INTENT(IN   ) :: KM
       IM1         ,&!INTEGER      , INTENT(IN   ) :: IM1
       JM1         ,&!INTEGER      , INTENT(IN   ) :: JM1
       KFRON       ,&!INTEGER      , INTENT(IN   ) :: KFRON
       KCI         ,&!INTEGER      , INTENT(IN   ) :: KCI
       DT          ,&!REAL(KIND=r8), INTENT(IN   ) :: DT
       KSU         ,&!INTEGER      , INTENT(IN   ) :: KSU(IM1,JM1)
       KSV         ,&!INTEGER      , INTENT(IN   ) :: KSV(IM1,JM1)
       PIC         ,&!REAL(KIND=r8), INTENT(IN   ) :: PIC(IM,JM)
       PIA         ,&!REAL(KIND=r8), INTENT(IN   ) :: PIA(IM,JM)
       UA          ,&!REAL(KIND=r8), INTENT(IN   ) :: UA(IM1,JM1,KM)
       VA          ,&!REAL(KIND=r8), INTENT(IN   ) :: VA(IM1,JM1,KM)
       UC          ,&!REAL(KIND=r8), INTENT(INOUT) :: UC(IM1,JM1,KM)
       VC          ,&!REAL(KIND=r8), INTENT(INOUT) :: VC(IM1,JM1,KM)
       UF          ,&!REAL(KIND=r8), INTENT(IN   ) :: UF(IM1,JM1,KM)
       VF          ,&!REAL(KIND=r8), INTENT(IN   ) :: VF(IM1,JM1,KM)
       UZ          ,&!REAL(KIND=r8), INTENT(IN   ) :: UZ(IM1,JM1,KM)
       VZ          ,&!REAL(KIND=r8), INTENT(IN   ) :: VZ(IM1,JM1,KM)
       TU          ,&!REAL(KIND=r8), INTENT(IN   ) :: TU(IM1,JM1,KM)
       TV          ,&!REAL(KIND=r8), INTENT(IN   ) :: TV(IM1,JM1,KM)
       WW           &!REAL(KIND=r8), INTENT(IN   ) :: WW(JM)
       )
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: IM
    INTEGER      , INTENT(IN   ) :: JM
    INTEGER      , INTENT(IN   ) :: KM
    INTEGER      , INTENT(IN   ) :: IM1
    INTEGER      , INTENT(IN   ) :: JM1
    INTEGER      , INTENT(IN   ) :: KFRON
    INTEGER      , INTENT(IN   ) :: KCI
    REAL(KIND=r8), INTENT(IN   ) :: DT
    INTEGER      , INTENT(IN   ) :: KSU(IM1,JM1)
    INTEGER      , INTENT(IN   ) :: KSV(IM1,JM1)
    REAL(KIND=r8), INTENT(IN   ) :: PIC(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: PIA(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: UA(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: VA(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(INOUT) :: UC(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(INOUT) :: VC(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: UF(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: VF(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: UZ(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: VZ(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: TU(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: TV(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: WW(JM)
    REAL(KIND=r8) :: A
    REAL(KIND=r8) :: A1
    REAL(KIND=r8) :: A2
    REAL(KIND=r8) :: AA
    REAL(KIND=r8) :: B
    REAL(KIND=r8) :: B1
    REAL(KIND=r8) :: B2
    REAL(KIND=r8) :: BB
    REAL(KIND=r8) :: CC
    REAL(KIND=r8) :: DD
    REAL(KIND=r8) :: PTU 
    REAL(KIND=r8) :: PIAM
    REAL(KIND=r8) :: PIAMJ
    REAL(KIND=r8) :: PTV
    INTEGER       :: KCI1
    INTEGER       :: KS

    INTEGER       :: i,j,k
    !REAL  NS,N,M,MI
    !DOUBLE PRECISION EP,WP,PHI,TPI,TU,TV,TT,PS                            
    !DOUBLE PRECISION PIF,UF,VF,TF,PIA,UA,VA,TA,PIC,UC,VC,TC
    !DOUBLE PRECISION GPX,GPY,GA,GF,GC
    !DOUBLE PRECISION EKI,EK,PHIC

    !DOUBLE PRECISION A,B
    !DOUBLE PRECISION PTU,PIAM,PTV,PIAMJ,AA,BB,CC,DD,A1,A2,B1,B2

    !PARAMETER(IM=146,JM=44,KM=7,IM1=IM,JM1=JM-1,KIM=KM+1,PT=20.)
    !COMMON/F/PIF(IM,JM),UF(IM1,JM1,KM),VF(IM1,JM1,KM),TF(IM,JM,KM)
    !COMMON/A/PIA(IM,JM),UA(IM1,JM1,KM),VA(IM1,JM1,KM),TA(IM,JM,KM)
    !COMMON/Z/PIZ(IM,JM),UZ(IM1,JM1,KM),VZ(IM1,JM1,KM),TZ(IM,JM,KM)

    !COMMON/C/PIC(IM,JM),UC(IM1,JM1,KM),VC(IM1,JM1,KM),TC(IM,JM,KM)
    !COMMON/G/GC(IM1,JM1,KM),GA(IM1,JM1,KM),GF(IM1,JM1,KM)
    !COMMON/TEN/TPI(IM,JM),TU(IM1,JM1,KM),TV(IM1,JM1,KM),TT(IM,JM,KM)    
    !COMMON/COORD/NS(IM,JM),M(JM),MI(JM),FLAT(JM),FLON(IM)
    !COMMON/BLOCOS/KSE(IM,JM),KSU(IM1,JM1),KSV(IM1,JM1)
    !COMMON/ETAS/EK(KM),EKI(KIM)
    !COMMON/TOPOG/PHIS(IM,JM),PRF(IM,JM),PS(IM,JM)
    !COMMON/DIAG/EP(IM,JM,KIM),WP(IM,JM,KM)
    !COMMON/PHISS/PHI(IM,JM,KIM),PHIC(IM,JM,KIM)
    !COMMON/GRAD/GPX(IM1,JM1,KM),GPY(IM1,JM1,KM)                   
    !COMMON/AUX/WW(JM),Q(IM,JM,KM),ALFA(IM,JM,KM)


    !2345678               
    !------------------------------------------------------------------------
    DO J=2,JM1-1
       DO I=2,IM1-1
          KS=KSU(I,J)           
          DO K=1,KS-1
             PTU =(PIC(I,J)+PIC(I+1,J))*0.5
             PIAM=(PIA(I,J)+PIA(I+1,J))*0.5
             UC(I,J,K)=(UA(I,J,K)*PIAM+2.*DT*TU(I,J,K))/PTU      
          END DO
       END DO
    END DO

    !------------------------------------------------------------------                      
    DO J=2,JM1-1
       DO I=2,IM1-1
          KS=KSV(I,J)
          DO K=1,KS-1
             PTV=(PIC(I,J)+PIC(I,J+1))*0.5
             PIAMJ=(PIA(I,J)+PIA(I,J+1))*0.5
             VC(I,J,K)=(VA(I,J,K)*PIAMJ+2.*DT*TV(I,J,K))/PTV           
          END DO
       END DO
    END DO

    !-------------------------------------------------------------------
    !**************************************************************************
    !                froneteira este oeste
    !********************************************************************

    KCI1=KCI
    IF(KCI1.EQ.1)THEN
       !          ciclico
       DO K=1,KM
          DO J=1,JM1           
             UC(1,J,K)=UC(IM1-1,J,K)
             UC(IM1,J,K)=UC(2,J,K)        
             VC(1,J,K)=VC(IM1-1,J,K)
             VC(IM1,J,K)=VC(2,J,K)
          END DO
       END DO

    ELSE

       !--------------------------------------------------------------
       !           R A D I A C I O N A L   L E S T E   O E S T E
       !--------------------------------------------------------------
       !                  LESTE UUUUUUUUUU
       !---------------------------------------------------------------
       DO K=1,KM
          DO J=2,JM1-1
             B2=UF(IM1-1,J,K)+UZ(IM1-1,J,K)-2.*UA(IM1-2,J,K)
             IF(B2.NE.0.)THEN
                B1=-(UF(IM1-1,J,K)-UZ(IM1-1,J,K))
                B=B1/B2
                IF(B.LT.0.)B=0.
                IF(B.GT.1.)B=1.
             ELSE
                B=0.
             ENDIF
             AA=(1.-B)/(1.+B)
             BB=2.*B/(1.+B)
             UC(IM1,J,K)=AA*UA(IM1,J,K)+BB*UF(IM1-1,J,K)
          END DO
       END DO
       !-----------------------------------------------------------------
       !                OESTE UUUUUUUUUUU
       !-----------------------------------------------------------------
       DO K=1,KM
          DO J=2,JM1-1
             A2=UF(2,J,K)+UZ(2,J,K)-2.*UA(3,J,K)
             IF(A2.NE.0.)THEN
                A1=UF(2,J,K)-UZ(2,J,K)
                A=A1/A2
                IF(A.GT.0.)A=0.
                IF(A.LT.-1.)A=-1.
             ELSE
                A=0.
             ENDIF
             CC=(1.+A)/(1.-A)
             DD=2.*A/(1.-A)
             UC(1,J,K)=CC*UA(1,J,K)-DD*UF(2,J,K)
          END DO
       END DO

       !------------------------------------------------------------
       !                  LESTE VVVVVVVVVV
       !-----------------------------------------------------------------
       DO K=1,KM
          DO J=2,JM1-1
             B2=VF(IM1-1,J,K)+VZ(IM1-1,J,K)-2.*VA(IM1-2,J,K)
             IF(B2.NE.0.)THEN
                B1=-(VF(IM1-1,J,K)-VZ(IM1-1,J,K))
                B=B1/B2
                IF(B.LT.0.)B=0.
                IF(B.GT.1.)B=1.
             ELSE
                B=0.
             ENDIF
             AA=(1.-B)/(1.+B)
             BB=2.*B/(1.+B)
             VC(IM1,J,K)=AA*VA(IM1,J,K)+BB*VF(IM1-1,J,K)
          END DO
       END DO

       !-----------------------------------------------------------------
       !                OESTE VVVVVVVV
       !-----------------------------------------------------------------
       DO K=1,KM
          DO J=2,JM1-1
             A2=VF(2,J,K)+VZ(2,J,K)-2.*VA(3,J,K)
             IF(A2.NE.0.)THEN
                A1=VF(2,J,K)-VZ(2,J,K)
                A=A1/A2
                IF(A.GT.0.)A=0.
                IF(A.LT.-1.)A=-1.
             ELSE
                A=0.
             ENDIF
             CC=(1.+A)/(1.-A)
             DD=2.*A/(1.-A)
             VC(1,J,K)=CC*VA(1,J,K)-DD*VF(2,J,K)
          END DO
       END DO

    ENDIF
    !      if ( b**2-4*a*c) 100,200,300
    !  100    print *, 'Roots are Complex'
    !      go to 400
    !  200    print *, 'Single Real Root' 
    !      go to 400
    !  300    print *, 'Two Real Roots'
    !  400 continue

    !--------------------------------------------------------------
    !IF(KFRON-2)100,200,300
    IF(KFRON-2 < 0)THEN
       !------------------------------------------------------------
       !                NEWMAN
       !------------------------------------------------------------
       DO K=1,KM
          DO I=1,IM1               
             UC(I,1,K)=UC(I,2,K)
             UC(I,JM1,K)=UC(I,JM1-1,K)
             VC(I,1,K)=0.
             VC(I,JM1,K)=0.
          END DO
       END DO

    ELSE IF(KFRON-2 == 0)THEN
       !-----------------------------------------------------------
       !               SPONJA
       !-------------------------------------------------------------
       DO K=1,KM
          DO I=1,IM1
             DO J=1,5
                UC(I,J,K)=UF(I,J,K)+WW(J)*TU(I,J,K)
                VC(I,J,K)=VF(I,J,K)+WW(J)*TV(I,J,K)
             END DO
          END DO
       END DO

       DO K=1,KM
          DO I=1,IM1
             DO J=JM1-5,JM1
                UC(I,J,K)=UF(I,J,K)+WW(J+1)*TU(I,J,K)
                VC(I,J,K)=VF(I,J,K)+WW(J+1)*TV(I,J,K)
             END DO
          END DO
       END DO

    ELSE IF(KFRON-2 > 0)THEN
       !-----------------------------------------------------------
       !           RADIACIONAL SUR UUUUUUUUUUUUUUU
       !--------------------------------------------------------------
       DO K=1,KM
          DO I=1,IM1
             A2=UF(I,2,K)+UZ(I,2,K)-2.*UA(I,3,K)
             IF(A2.NE.0.)THEN
                A1=UF(I,2,K)-UZ(I,2,K)
                A=A1/A2
                IF(A.GT.0.)A=0.
                IF(A.LT.-1.)A=-1.
             ELSE
                A=0.
             ENDIF
             CC=(1.+A)/(1.-A)
             DD=2.*A/(1.-A)
             UC(I,1,K)=CC*UA(I,1,K)-DD*UF(I,2,K)
          END DO
       END DO

       !-------------------------------------------------------------
       !            RADIACIONAL NORTE UUUUUUUUUUUU
       !------------------------------------------------------------
       DO K=1,KM
          DO I=1,IM1
             B2=UF(I,JM1-1,K)+UZ(I,JM1-1,K)-2.*UA(I,JM1-2,K)
             IF(B2.NE.0.)THEN
                B1=-(UF(I,JM1-1,K)-UZ(I,JM1-1,K))
                B=B1/B2
                IF(B.LT.0.)B=0.
                IF(B.GT.1.)B=1.
             ELSE
                B=0.
             ENDIF
             AA=(1.-B)/(1.+B)
             BB=2.*B/(1.+B)
             UC(I,JM1,K)=AA*UA(I,JM1,K)+BB*UF(I,JM1-1,K)     
          END DO
       END DO

       !-----------------------------------------------------------
       !             RADIACIONAL SUR VVVVVVVVVVV
       !-----------------------------------------------------------
       DO K=1,KM
          DO I=1,IM1
             A2=VF(I,2,K)+VZ(I,2,K)-2.*VA(I,3,K)
             IF(A2.NE.0.)THEN
                A1=VF(I,2,K)-VZ(I,2,K)
                A=A1/A2
                IF(A.GT.0.)A=0.
                IF(A.LT.-1.)A=-1.
             ELSE
                A=0.
             ENDIF
             CC=(1.+A)/(1.-A)
             DD=2.*A/(1.-A)
             VC(I,1,K)=CC*VA(I,1,K)-DD*VF(I,2,K)
          END DO
       END DO

       !---------------------------------------------------------------
       !            RADIACIONAL NORTE VVVVVVVVVV
       !---------------------------------------------------------------
       DO K=1,KM
          DO I=1,IM1
             B2=VF(I,JM1-1,K)+VZ(I,JM1-1,K)-2.*VA(I,JM1-2,K)
             IF(B2.NE.0.)THEN 
                B1=-(VF(I,JM1-1,K)-VZ(I,JM1-1,K))
                B=B1/B2
                IF(B.LT.0.)B=0.
                IF(B.GT.1.)B=1.
             ELSE
                B=0.
             ENDIF
             AA=(1.-B)/(1.+B)
             BB=2.*B/(1.+B)    
             VC(I,JM1,K)=AA*VA(I,JM1,K)+BB*VF(I,JM1-1,K)
          END DO
       END DO

       !----------------------------------------------------------------


    END IF
  END SUBROUTINE LEAP

  !*****************************************************************
  !        GERA DADOS PARA O GRAFICO
  SUBROUTINE SWRITE(&
       IREC1    ,&!INTEGER      , INTENT(INOUT) :: IREC1
       IM       ,&!INTEGER      , INTENT(IN   ) :: IM
       JM       ,&!INTEGER      , INTENT(IN   ) :: JM
       KM       ,&!INTEGER      , INTENT(IN   ) :: KM
       IM1      ,&!INTEGER      , INTENT(IN   ) :: IM1
       JM1      ,&!INTEGER      , INTENT(IN   ) :: JM1
       IM2      ,&!INTEGER      , INTENT(IN   ) :: IM2
       JM2      ,&!INTEGER      , INTENT(IN   ) :: JM2
       KIM      ,&!INTEGER      , INTENT(IN   ) :: KIM
       PT       ,&!REAL(KIND=r8), INTENT(IN   ) :: PT
       KSE      ,&!INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
       KSU      ,&!INTEGER      , INTENT(IN   ) :: KSU(IM1,JM1)
       KSV      ,&!INTEGER      , INTENT(IN   ) :: KSV(IM1,JM1)
       PIC      ,&!REAL(KIND=r8), INTENT(IN   ) :: PIC(IM,JM)
       EKI      ,&!REAL(KIND=r8), INTENT(IN   ) :: EKI(KIM)
       PHIC     ,&!REAL(KIND=r8), INTENT(IN   ) :: PHIC(IM,JM,KIM)
       WP       ,&!REAL(KIND=r8), INTENT(IN   ) :: WP(IM,JM,KM)
       UC       ,&!REAL(KIND=r8), INTENT(IN   ) :: UC(IM1,JM1,KM)
       VC       )!REAL(KIND=r8), INTENT(IN   ) :: VC(IM1,JM1,KM)
    IMPLICIT NONE
    INTEGER      , INTENT(INOUT) :: IREC1
    INTEGER      , INTENT(IN   ) :: IM
    INTEGER      , INTENT(IN   ) :: JM
    INTEGER      , INTENT(IN   ) :: KM
    INTEGER      , INTENT(IN   ) :: IM1
    INTEGER      , INTENT(IN   ) :: JM1
    INTEGER      , INTENT(IN   ) :: IM2
    INTEGER      , INTENT(IN   ) :: JM2
    INTEGER      , INTENT(IN   ) :: KIM
    REAL(KIND=r8), INTENT(IN   ) :: PT
    INTEGER      , INTENT(IN   ) :: KSE(IM,JM)
    INTEGER      , INTENT(IN   ) :: KSU(IM1,JM1)
    INTEGER      , INTENT(IN   ) :: KSV(IM1,JM1)
    REAL(KIND=r8), INTENT(IN   ) :: PIC(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: EKI(KIM)
    REAL(KIND=r8), INTENT(IN   ) :: PHIC(IM,JM,KIM)
    REAL(KIND=r8), INTENT(IN   ) :: WP(IM,JM,KM)
    REAL(KIND=r8), INTENT(IN   ) :: UC(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: VC(IM1,JM1,KM)
    REAL(KIND=r8) :: DD1
    REAL(KIND=r8) :: DD2
    REAL(KIND=r8) :: DD3
    REAL(KIND=r8) :: NS5
    REAL(KIND=r8) :: NS6
    REAL(KIND=r8) :: NS7
    INTEGER       :: lrec
    INTEGER       :: KS
    INTEGER       :: i,j,k

    !REAL  NS,N,M,MI

    !        DOUBLE PRECISION EP,WP,PHI,TPI,TU,TV,TT,PS                            
    !        DOUBLE PRECISION PIF,UF,VF,TF,PIA,UA,VA,TA,PIC,UC,VC,TC
    !        DOUBLE PRECISION GPX,GPY,GA,GF,GC
    !        DOUBLE PRECISION EKI,EK,PHIC
    !      PARAMETER(IM=146,JM=44,KM=7,IM1=IM,JM1=JM-1,KIM=KM+1,PT=20.)
    !      PARAMETER(IM2=IM-2,JM2=JM-2,RS=287)
    !------------------------------------------------------------------
    !      COMMON/F/PIF(IM,JM),UF(IM1,JM1,KM),VF(IM1,JM1,KM),TF(IM,JM,KM)
    !      COMMON/A/PIA(IM,JM),UA(IM1,JM1,KM),VA(IM1,JM1,KM),TA(IM,JM,KM)
    !      COMMON/Z/PIZ(IM,JM),UZ(IM1,JM1,KM),VZ(IM1,JM1,KM),TZ(IM,JM,KM)

    !      COMMON/C/PIC(IM,JM),UC(IM1,JM1,KM),VC(IM1,JM1,KM),TC(IM,JM,KM)
    !      COMMON/G/GC(IM1,JM1,KM),GA(IM1,JM1,KM),GF(IM1,JM1,KM)
    !      COMMON/TEN/TPI(IM,JM),TU(IM1,JM1,KM),TV(IM1,JM1,KM),TT(IM,JM,KM)    
    !      COMMON/COORD/NS(IM,JM),M(JM),MI(JM),FLAT(JM),FLON(IM)
    !      COMMON/BLOCOS/KSE(IM,JM),KSU(IM1,JM1),KSV(IM1,JM1)
    !      COMMON/ETAS/EK(KM),EKI(KIM)
    !      COMMON/TOPOG/PHIS(IM,JM),PRF(IM,JM),PS(IM,JM)
    !      COMMON/DIAG/EP(IM,JM,KIM),WP(IM,JM,KM)
    !      COMMON/PHISS/PHI(IM,JM,KIM),PHIC(IM,JM,KIM)
    !      COMMON/GRAD/GPX(IM1,JM1,KM),GPY(IM1,JM1,KM)                   
    !      COMMON/AUX/WW(JM),Q(IM,JM,KM),ALFA(IM,JM,KM)

    REAL(KIND=r8) ::  DPHI(IM2,JM2,KM)
    REAL(KIND=r8) ::  VX(IM2,JM2,KM)
    REAL(KIND=r8) ::  VY(IM2,JM2,KM)
    REAL(KIND=r8) ::  PP(IM,JM,KIM)
    REAL(KIND=r8) ::  SP(KM)
    REAL(KIND=r8) ::  PHIG(IM,JM,KM)
    REAL(KIND=r8) ::  SS(KM),SSM(KM)
    REAL(KIND=r8) ::  VZZ(IM2,JM2,KM)
    REAL(KIND=r8) ::  VZ1(IM2,JM2,KM)
    REAL(KIND=r8) ::  DPHI1(IM2,JM2,KM)
    REAL(KIND=r8) ::  VX1(IM2,JM2,KM)
    REAL(KIND=r8) ::  VY1(IM2,JM2,KM)
    REAL(KIND=r8) ::  UCC(IM1,JM1,KM)
    REAL(KIND=r8) ::  VCC(IM1,JM1,KM)
    REAL(KIND=r8) ::  WCC(IM,JM,KM)
    REAL(KIND=r4) ::  AUX_r4(IM2,JM2,KM)

    !-----------------------------------------------------------------------
    !          PP=PRESAO EM CADA NIVEL K+1/2
    !          PHI=ALT GEOP. EM CADA NIVEL K+1/2
    !          PHIG=ALT GEOPOTENCIAL EM NINEL "K"
    !          DPHI= DESVIO DE PHI EM NIVEL "K"
    !          VX  VY COMPO. DO VENTO EN X Y EM (I,J) DE PHI
    !          SP= NIVEL ISOBARICO(SUPERFICIE ISOBARICA)
    !          SS=SUMATORIA DE PHI EN CADA NIVEL
    !          SSM= MEDIA DE GEOPO EN CADA NIVEL "K"(ONDE SE ESTA INTER)
    !-----------------------------------------------------------------------
    DO I=1,IM
       DO J=1,JM
          DO K=1,KIM
             PP(I,J,K)=EKI(K)*PIC(I,J)+PT       
          END DO
       END DO
    END DO

    !-----------------------------------------------------------------------
    SP(1)=100.
    SP(2)=200.
    SP(3)=300
    SP(4)=500.
    SP(5)=700.
    SP(6)=850.
    SP(7)=960.
    !----------------------------------------------------------------------
    !         INTERPOLACAIO A NIVELIES ISOBARICOS SP de phi
    !----------------------------------------------------------------------          
    DO I=1,IM
       DO J=1,JM
          KS=KSE(I,J)
          DO K=1,KS-1
             DD1=SP(K)-PP(I,J,K)     
             DD2=PP(I,J,K+1)-PP(I,J,K)
             DD3=DD1/DD2
             PHIG(I,J,K)=PHIC(I,J,K)+(PHIC(I,J,K+1)-PHIC(I,J,K))*DD3        
          END DO
       END DO
    END DO

    !------------------------- ----------------------------------
    !             interpolacion de wP
    !-----------------------------------------------------------
    WRITE(*,*)'pasei www'
    DO I=1,IM
       DO J=1,JM
          KS=KSE(I,J)
          DO K=1,KS-1
             DD1=SP(K)-PP(I,J,K)     
             DD2=PP(I,J,K+1)-PP(I,J,K)
             DD3=DD1/DD2
             IF(K.LT.7)THEN
                WCC(I,J,K)=WP(I,J,K)+(WP(I,J,K+1)-WP(I,J,K))*DD3 
             ELSE
                WCC(I,J,K)=WP(I,J,K)
             ENDIF
          END DO
       END DO
    END DO

    WRITE(*,*)'pasei www2'
    !--------------------------------------------------------------------
    !      INTERPOLACAOA NIVEIS ISOVARICOS DE U 
    !---------------------------------------------------------------------
    DO I=1,IM1
       DO J=1,JM1
          KS=KSU(I,J)
          DO K=1,KS-1
             DD1=SP(K)-PP(I,J,K)     
             DD2=PP(I,J,K+1)-PP(I,J,K)
             DD3=DD1/DD2
             IF(K.LT.7)THEN
                UCC(I,J,K)=UC(I,J,K)+(UC(I,J,K+1)-UC(I,J,K))*DD3 
             ELSE
                UCC(I,J,K)=UC(I,J,K)
             ENDIF
          END DO
       END DO
    END DO

    WRITE(*,*)'pasei www3'
    !-------------------------------------------------------------------
    ! INTERPOLACAOA NIVEIS ISOVARICOS DE V 
    !---------------------------------------------------------------------
    DO I=1,IM1
       DO J=1,JM1
          KS=KSV(I,J)
          DO K=1,KS-1
             DD1=SP(K)-PP(I,J,K)     
             DD2=PP(I,J,K+1)-PP(I,J,K)
             DD3=DD1/DD2
             IF(K.LT.7)THEN
                VCC(I,J,K)=VC(I,J,K)+(VC(I,J,K+1)-VC(I,J,K))*DD3 
             ELSE
                VCC(I,J,K)=VC(I,J,K)
             ENDIF
          END DO
       END DO
    END DO

    WRITE(*,*)'pasei www4'
    !-----------------------------------------------------------------
    !          CALCULO DOS DESVIOS DE PHI
    !           SS(K)   SUMA DE TODAS AS GEOPOTECIASI EM NIVEL K
    !-----------------------------------------------------------------
    DO K=1,KM
       SS(K)=0.
    END DO

    NS5=0
    NS6=0  
    NS7=0          
    !---------------------------------------------------------------

    DO K=1,4
       DO J=1,JM
          DO I=1,IM-2        
             SS(K)=SS(K)+PHIG(I,J,K)
          END DO
       END DO
    END DO

    SSM(1)=SS(1)/(JM*IM2)
    SSM(2)=SS(2)/(JM*IM2)
    SSM(3)=SS(3)/(JM*IM2) 
    SSM(4)=SS(4)/(JM*IM2) 
    !------------------------------------------------------------------
    DO I=1,IM-2
       DO J=1,JM
          KS=KSE(I,J)
          IF(KS.GE.6)THEN
             NS5=NS5+1    
             SS(5)=SS(5)+PHIG(I,J,5)
          ENDIF
          !-----------------------------------------------------------------
          IF(KS.GE.7)THEN
             NS6=NS6+1
             SS(6)=SS(6)+PHIG(I,J,6)
          ENDIF
          !----------------------------------------------------------------
          IF(KS.GE.8)THEN
             NS7=NS7+1    
             SS(7)=SS(7)+PHIG(I,J,7)
          ENDIF
       END DO
    END DO

    !---------------------------------------------------------
    SSM(5)=SS(5)/NS5
    SSM(6)=SS(6)/NS6
    SSM(7)=SS(7)/NS7
    !-----------------------------------------------------
    !           U en puntos de PHI
    !-----------------------------------------------------
    DO I=2,IM-2
       DO J=2,JM-1
          DO K=1,KM        
             KS=KSE(I,J)
             IF(K.GT.KS-1)THEN
                VX(I,J-1,K)=-9.99E10
             ELSE
                VX(I,J-1,K)=(UCC(I-1,J-1,K)+UCC(I,J-1,K))*0.5
             ENDIF
          END DO
       END DO
    END DO

    !-----------------------------------------------
    ! completando el primer punto
    !----------------------------------------------                     
    DO j=2,Jm-1
       DO k=1,KM
          KS=KSE(1,J)
          IF(K.GT.KS-1)THEN
             VX(1,J-1,K)=-9.99E10
          ELSE
             VX(1,J-1,K)=(UCC(IM-2,J-1,K)+UCC(1,J-1,K))*0.5              
          ENDIF
       ENDDO
    ENDDO
    !2345678
    !-----------------------------------------------------
    !           v en puntos de PHI
    !-----------------------------------------------------
    DO I=1,IM-2
       DO J=2,JM-1
          DO K=1,KM        
             KS=KSE(I,J)
             IF(K.GT.KS-1)THEN         
                VY(I,J-1,K)=-9.99E10   
             ELSE           
                VY(I,J-1,K)=(VCC(I,J-1,K)+VCC(I,J,K))*0.5
             ENDIF
          END DO
       END DO
    END DO

    !2345678
    !-----------------------------------------
    !           calculo de desvios de phi y Wp
    !-----------------------------------------              
    DO I=1,IM-2
       DO J=2,JM1
          DO K=1,KM                          
             KS=KSE(I,J)
             IF(K.GT.KS-1)THEN
                DPHI(I,J-1,K)=-9.99E10
                VZZ(I,J-1,K)=-9.99E10
             ELSE
                !           DPHI(I,J-1,K)=PHIG(I,J,K)/9.8
                DPHI(I,J-1,K)=(PHIG(I,J,K)-SSM(K))/9.8                   
                VZZ(I,J-1,K)=WCC(I,J,K)
             ENDIF
          END DO
       END DO
    END DO

    !-------------------------------------------------      
    WRITE(*,*)'pasei www5'
    !---------------------------------------------------------
    !        CAMBIO PARA GRAFICO DE 960 PARA 100 (1-7)
    !----------------------------------------------------------


    DO  I=1,IM2
       DO  J=1,JM2
          DO  K=KM,1,-1     
             DPHI1(i,j,8-k)=DPHI(i,j,k)
             VX1(i,j,8-k)=VX(i,j,k)
             VY1(i,j,8-k)=VY(i,j,k)
             VZ1(i,j,8-k)=VZZ(i,j,k)
          END DO
       END DO
    END DO

    !--------------------------------------------------------
    !       GRAVACAO DOS DADOS NAO FORMATADO
    !          contopomedia.............con topo y media icinial
    !          sintopomedia.......... sin topog.y media inicial 
    !          contoporeposo..........con topor. inicial reposo
    !          sintoporeposo..........sin topog. inicial reposo
    !          00......................sin forcante
    !          01..................... 06 h. pico seno un ciclo
    !          02..................... 06 h. pico seno varios ciclos
    !          03..................... 12 h. pico seno un ciclo.
    !          04..................... estacionaria
    !          05..................... 06 horas de pico cub.
    !          06..................... 12horas de un solo pico cub.
    !          07..................... 12 h. pico seno ciclico.
    !-----------------------------------------------------------
    !          -04           cada 4 horas
    !          -12           cada 12 horas
    !          -24           cada 24 horas
    !2345678---------------------------------------------------
    !      OPEN(20,FILE='sintopomedia00-12h.dat',STATUS='UNKNOWN',
    !      OPEN(20,FILE='contopomedia00-12h.dat',STATUS='UNKNOWN',

    !      OPEN(20,FILE='sintopomedia04-12h.dat',STATUS='UNKNOWN',
    !      OPEN(20,FILE='contopomedia02-12h.dat',STATUS='UNKNOWN',

    !      OPEN(20,FILE='sintporeposo06-12h.dat',STATUS='UNKNOWN',
    !      OPEN(20,FILE='contopomedia06-12h.dat',STATUS='UNKNOWN',
    !--------------------------------------------------------------
    !2345678
    IF(IREC1==0)THEN
       INQUIRE(IOLENGTH=lrec)AUX_r4
       OPEN(20,FILE='eta.dat',STATUS='UNKNOWN',form='unformatted',access='direct',recl=lrec)
    END IF

    irec1=irec1+1
    AUX_r4=REAL(DPHI1,kind=r4)
    WRITE(20,rec=irec1)AUX_r4

    irec1=irec1+1
    AUX_r4=REAL(VX1,kind=r4)
    WRITE(20,rec=irec1)AUX_r4

    irec1=irec1+1
    AUX_r4=REAL(VY1,kind=r4)
    WRITE(20,rec=irec1)AUX_r4

    irec1=irec1+1
    AUX_r4=REAL(VZ1,kind=r4)
    WRITE(20,rec=irec1)AUX_r4

    !----------------------------------------------
    RETURN
  END SUBROUTINE SWRITE


  !2345678****************************************************************
  SUBROUTINE FILTER(&
       KTOP      ,&! INTEGER, INTENT(IN   ) ::  KTOP
       CF        ,&! REAL(KIND=r8), INTENT(IN   ) ::  CF  
       IM        ,&! INTEGER, INTENT(IN   ) ::  IM
       JM        ,&! INTEGER, INTENT(IN   ) ::  JM
       KM        ,&! INTEGER, INTENT(IN   ) ::  KM
       IM1       ,&! INTEGER, INTENT(IN   ) ::  IM1
       JM1       ,&! INTEGER, INTENT(IN   ) ::  JM1
       KSE       ,&! INTEGER      , INTENT(IN   ) :: KSE(IM,JM)  !COMMON/BLOCOS/KSE(IM,JM)
       KSU       ,&! INTEGER      , INTENT(IN   ) :: KSU(IM1,JM1)
       KSV       ,&! INTEGER      , INTENT(IN   ) :: KSV(IM1,JM1)
       TA        ,&! REAL(KIND=r8), INTENT(IN   ) :: TA(IM,JM,KM)
       TC        ,&! REAL(KIND=r8), INTENT(INOUT) :: TC(IM,JM,KM)
       TF        ,&! REAL(KIND=r8), INTENT(INOUT) :: TF(IM,JM,KM)
       UA        ,&! REAL(KIND=r8), INTENT(IN   ) :: UA(IM1,JM1,KM)
       VA        ,&! REAL(KIND=r8), INTENT(IN   ) :: VA(IM1,JM1,KM)
       UC        ,&! REAL(KIND=r8), INTENT(INOUT) :: UC(IM1,JM1,KM)
       VC        ,&! REAL(KIND=r8), INTENT(INOUT) :: VC(IM1,JM1,KM)
       UF        ,&! REAL(KIND=r8), INTENT(INOUT) :: UF(IM1,JM1,KM)
       VF        ,&! REAL(KIND=r8), INTENT(INOUT) :: VF(IM1,JM1,KM)
       PIF       ,&! REAL(KIND=r8), INTENT(INOUT) :: PIF(IM,JM)
       PIC       ,&! REAL(KIND=r8), INTENT(INOUT) :: PIC(IM,JM)
       PIA        )! REAL(KIND=r8), INTENT(IN   ) :: PIA(IM,JM)
    IMPLICIT NONE
    INTEGER, INTENT(IN   ) ::  KTOP
    REAL(KIND=r8), INTENT(IN   ) ::  CF  
    INTEGER, INTENT(IN   ) ::  IM
    INTEGER, INTENT(IN   ) ::  JM
    INTEGER, INTENT(IN   ) ::  KM
    INTEGER, INTENT(IN   ) ::  IM1
    INTEGER, INTENT(IN   ) ::  JM1
    INTEGER      , INTENT(IN   ) :: KSE(IM,JM)  !COMMON/BLOCOS/KSE(IM,JM)
    INTEGER      , INTENT(IN   ) :: KSU(IM1,JM1)
    INTEGER      , INTENT(IN   ) :: KSV(IM1,JM1)
    REAL(KIND=r8), INTENT(IN   ) :: TA(IM,JM,KM)
    REAL(KIND=r8), INTENT(INOUT) :: TC(IM,JM,KM)
    REAL(KIND=r8), INTENT(INOUT) :: TF(IM,JM,KM)
    REAL(KIND=r8), INTENT(IN   ) :: UA(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: VA(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(INOUT) :: UC(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(INOUT) :: VC(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(INOUT) :: UF(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(INOUT) :: VF(IM1,JM1,KM)
    REAL(KIND=r8), INTENT(INOUT) :: PIF(IM,JM)
    REAL(KIND=r8), INTENT(INOUT) :: PIC(IM,JM)
    REAL(KIND=r8), INTENT(IN   ) :: PIA(IM,JM)

    REAL(KIND=r8) :: AA
    REAL(KIND=r8) :: BB
    REAL(KIND=r8) :: CC
    REAL(KIND=r8) :: DD

    INTEGER :: J,I,K
    INTEGER :: KSX
    INTEGER :: KS
    INTEGER :: KSY


    !------------------------------------------------------------
    !           FILTRO DE SHUMAN
    !2345678------------------------------------------------------------

    DO  J=2,JM-1
       DO  I=2,IM-1
          KS=KSE(I,J)
          DO  K=1,KS-1          
             AA=PIC(I-1,J)+PIC(I+1,J)+PIC(I,J-1)+PIC(I,J+1)
             PIC(I,J)=PIC(I,J)+CF*(AA-4.*PIC(I,J))*0.25
             BB=TC(I-1,J,K)+TC(I+1,J,K)+TC(I,J-1,K)+TC(I,J+1,K)
             TC(I,J,K)=TC(I,J,K)+CF*(BB-4.*TC(I,J,K))*0.25
          END DO
       END DO
    END DO

    ! ----------------------------------------------------------
    IF(KTOP.EQ.1)THEN
       DO I=2,IM-1
          DO J=2,JM-1
             KS=KSE(I,J)
             KSX=KSE(I+1,J)
             KSY=KSE(I,J+1)

             IF((KSX-KS) .LT. 0 ) TC(I+1,J,KS-1) = TC(I,J,KS-1)
             IF((KSX-KS).GT.0)TC(I,J,KS)=TC(I+1,J,KS)
             IF((KSY-KS).LT.0)TC(I,J+1,KS-1)=TC(I,J,KS-1)
             IF((KSY-KS).GT.0)TC(I,J,KS)=TC(I,J+1,KS)
          END DO
       END DO

       !==============================
    ENDIF
    !------------------------------------------------------------
    DO J=2,JM1-1
       DO I=2,IM1-1
          KS=KSU(I,J)
          DO K=1,KS-1
             CC=UC(I-1,J,K)+UC(I+1,J,K)+UC(I,J-1,K)+UC(I,J+1,K)
             UC(I,J,K)=UC(I,J,K)+CF*(CC-4.*UC(I,J,K))*0.25
          END DO
       END DO
    END DO

    !-----------------------------------------------------------
    DO J=2,JM1-1
       DO I=2,IM1-1
          KS=KSV(I,J)
          DO K=1,KS-1
             DD=VC(I-1,J,K)+VC(I+1,J,K)+VC(I,J-1,K)+VC(I,J+1,K)
             VC(I,J,K)=VC(I,J,K)+CF*(DD-4.*VC(I,J,K))*0.25
          END DO
       END DO
    END DO

    !------------------------------------------------------------
    !            post filtro fronteira oeste
    !-----------------------------------------------------------

    !-----------------------------------------------------------                                   
    !              FILTRO DE ASELIN
    !-----------------------------------------------------------
    DO  J=1,JM
       DO  I=1,IM
          KS=KSE(I,J)
          DO  K=1,KS-1
             PIF(I,J)=0.9*PIF(I,J)+0.05*(PIC(I,J)+PIA(I,J))
             TF(I,J,K)=0.9*TF(I,J,K)+0.05*(TC(I,J,K)+TA(I,J,K))
          END DO
       END DO
    END DO

    !-------------------------------------------------------------
    DO  J=1,JM1
       DO  I=1,IM1
          KS=KSE(I,J)
          DO  K=1,KS-1
             VF(I,J,K)=0.9*VF(I,J,K)+0.05*(VC(I,J,K)+VA(I,J,K)) 
             UF(I,J,K)=0.9*UF(I,J,K)+0.05*(UC(I,J,K)+UA(I,J,K))
          END DO
       END DO
    END DO

    !-------------------------------------------------------------                       
    RETURN
  END SUBROUTINE FILTER


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !******************************************************************
  SUBROUTINE TROCA(&
       IM         ,&!!& INTEGER      , INTENT(IN   ) :: IM
       JM         ,&!!& INTEGER      , INTENT(IN   ) :: JM
       KM         ,&!!& INTEGER      , INTENT(IN   ) :: KM
       IM1        ,&!!& INTEGER      , INTENT(IN   ) :: IM1
       JM1        ,&!!& INTEGER      , INTENT(IN   ) :: JM1
       PIF        ,&!!& REAL(KIND=r8), INTENT(INOUT) :: PIF (IM,JM)        !       COMMON/F/
       TA         ,&!!& REAL(KIND=r8), INTENT(INOUT) :: TA  (IM,JM,KM)
       TC         ,&!!& REAL(KIND=r8), INTENT(IN   ) :: TC  (IM,JM,KM)
       UC         ,&!!& REAL(KIND=r8), INTENT(IN   ) :: UC  (IM1,JM1,KM)
       VC         ,&!!& REAL(KIND=r8), INTENT(IN   ) :: VC  (IM1,JM1,KM)
       GC         ,&!!& REAL(KIND=r8), INTENT(IN   ) :: GC  (IM1,JM1,KM)
       UF         ,&!!& REAL(KIND=r8), INTENT(INOUT) :: UF  (IM1,JM1,KM) ! COMMON/F
       VF         ,&!!& REAL(KIND=r8), INTENT(INOUT) :: VF  (IM1,JM1,KM) ! COMMON/F
       GF         ,&!!& REAL(KIND=r8), INTENT(INOUT) :: GF  (IM1,JM1,KM)
       TF         ,&!!& REAL(KIND=r8), INTENT(INOUT) :: TF  (IM ,JM ,KM) ! COMMON/F/
       UZ         ,&!!& REAL(KIND=r8), INTENT(OUT  ) :: UZ  (IM1,JM1,KM)
       VZ         ,&!!& REAL(KIND=r8), INTENT(OUT  ) :: VZ  (IM1,JM1,KM)
       PIA        ,&!!& REAL(KIND=r8), INTENT(INOUT) :: PIA (IM,JM)        !       COMMON/A/
       UA         ,&!!& REAL(KIND=r8), INTENT(INOUT) :: UA  (IM1,JM1,KM) ! COMMON/A/
       VA         ,&!!& REAL(KIND=r8), INTENT(INOUT) :: VA  (IM1,JM1,KM) ! COMMON/A/
       GA         ,&!!& REAL(KIND=r8), INTENT(OUT  ) :: GA  (IM1,JM1,KM)
       PIZ        ,&!!& REAL(KIND=r8), INTENT(OUT  ) :: PIZ (IM,JM)        !       COMMON/Z/
       TZ         ,&!!& REAL(KIND=r8), INTENT(OUT  ) :: TZ  (IM,JM,KM) ! COMMON/Z/
       PIC         )!  REAL(KIND=r8), INTENT(INOUT) :: PIC (IM,JM)        !       COMMON/C/
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: IM
    INTEGER      , INTENT(IN   ) :: JM
    INTEGER      , INTENT(IN   ) :: KM
    INTEGER      , INTENT(IN   ) :: IM1
    INTEGER      , INTENT(IN   ) :: JM1
    REAL(KIND=r8), INTENT(INOUT) :: PIF (IM,JM)        !       COMMON/F/
    REAL(KIND=r8), INTENT(INOUT) :: TA  (IM,JM,KM)
    REAL(KIND=r8), INTENT(IN   ) :: TC  (IM,JM,KM)
    REAL(KIND=r8), INTENT(IN   ) :: UC  (IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: VC  (IM1,JM1,KM)
    REAL(KIND=r8), INTENT(IN   ) :: GC  (IM1,JM1,KM)
    REAL(KIND=r8), INTENT(INOUT) :: UF  (IM1,JM1,KM)    !       COMMON/F
    REAL(KIND=r8), INTENT(INOUT) :: VF  (IM1,JM1,KM)    !       COMMON/F
    REAL(KIND=r8), INTENT(INOUT) :: GF  (IM1,JM1,KM)
    REAL(KIND=r8), INTENT(OUT  ) :: UZ  (IM1,JM1,KM)
    REAL(KIND=r8), INTENT(OUT  ) :: VZ  (IM1,JM1,KM)
    REAL(KIND=r8), INTENT(INOUT) :: TF  (IM,JM,KM)      !       COMMON/F/
    REAL(KIND=r8), INTENT(INOUT) :: PIA (IM,JM)        !       COMMON/A/
    REAL(KIND=r8), INTENT(INOUT) :: UA  (IM1,JM1,KM)    !       COMMON/A/
    REAL(KIND=r8), INTENT(INOUT) :: VA  (IM1,JM1,KM)    !       COMMON/A/
    REAL(KIND=r8), INTENT(OUT  ) :: GA  (IM1,JM1,KM)
    REAL(KIND=r8), INTENT(OUT  ) :: PIZ (IM,JM)        !       COMMON/Z/
    REAL(KIND=r8), INTENT(OUT  ) :: TZ  (IM,JM,KM)      !       COMMON/Z/
    REAL(KIND=r8), INTENT(INOUT) :: PIC (IM,JM)        !       COMMON/C/


    INTEGER :: K
    INTEGER :: J
    INTEGER :: I  
    !---------------------------------
    DO K=1,KM  
       DO J=1,JM
          DO I=1,IM
             PIZ (I,J)  =PIA(I,J)
             TZ  (I,J,K)=TA (I,J,K)
          END DO
       END DO
    END DO

    !-------------------------------
    DO  K=1,KM
       DO  J=1,JM
          DO  I=1,IM
             PIA(I,J)=PIF(I,J)
             TA(I,J,K)=TF(I,J,K)
          END DO
       END DO
    END DO

    !-------------------------------
    DO  K=1,KM
       DO  J=1,JM
          DO  I=1,IM
             PIF(I,J)=PIC(I,J)
             TF(I,J,K)=TC(I,J,K)
          END DO
       END DO
    END DO

    !--------------------------------

    DO  K=1,KM
       DO  J=1,JM1
          DO  I=1,IM1
             UZ(I,J,K)=UA(I,J,K)
             VZ(I,J,K)=VA(I,J,K)
          END DO
       END DO
    END DO
    !-----------------------------

    DO  K=1,KM
       DO  J=1,JM1
          DO  I=1,IM1
             UA(I,J,K)=UF(I,J,K)
             VA(I,J,K)=VF(I,J,K)
             GA(I,J,K)=GF(I,J,K) 
          END DO
       END DO
    END DO

    !--------------------------------------
    DO  K=1,KM
       DO  J=1,JM1
          DO  I=1,IM1
             UF(I,J,K)=UC(I,J,K)
             VF(I,J,K)=VC(I,J,K)
             GF(I,J,K)=GC(I,J,K)
          END DO
       END DO
    END DO

    !-------------------------------
    RETURN
  END SUBROUTINE TROCA

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

END MODULE MODULE_ETA_MODEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM Main
 USE MODULE_ETA_MODEL, ONLY : eta
  IMPLICIT NONE
  ! Selecting Kinds
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(15)! Kind for 128-bits Real Numbers

      INTEGER      ,  PARAMETER :: IM=146
      INTEGER      ,  PARAMETER :: JM=44
      INTEGER      ,  PARAMETER :: KM=7
      INTEGER      ,  PARAMETER :: IM1=IM
      INTEGER      ,  PARAMETER :: JM1=JM-1
      INTEGER      ,  PARAMETER :: KIM=KM+1
      REAL(KIND=r8),  PARAMETER :: PT=20.0

        
!     REAL(KIND=r8) :: PIF(IM,JM)
!     REAL(KIND=r8) :: UF(IM1,JM1,KM)
!     REAL(KIND=r8) :: VF(IM1,JM1,KM)
!     REAL(KIND=r8) :: TF(IM,JM,KM)
!     REAL(KIND=r8) :: PIA(IM,JM)
!     REAL(KIND=r8) :: UA(IM1,JM1,KM)
!     REAL(KIND=r8) :: VA(IM1,JM1,KM)
!     REAL(KIND=r8) :: TA(IM,JM,KM)
!     REAL(KIND=r8) :: PIZ(IM,JM)
!     REAL(KIND=r8) :: UZ(IM1,JM1,KM)
!     REAL(KIND=r8) :: VZ(IM1,JM1,KM)
!     REAL(KIND=r8) :: TZ(IM,JM,KM)

!     REAL(KIND=r8) :: PIC(IM,JM)
!     REAL(KIND=r8) :: UC(IM1,JM1,KM)
!     REAL(KIND=r8) :: VC(IM1,JM1,KM)
!     REAL(KIND=r8) :: TC(IM,JM,KM)
!     REAL(KIND=r8) :: GC(IM1,JM1,KM)
!     REAL(KIND=r8) :: GA(IM1,JM1,KM)
!     REAL(KIND=r8) :: GF(IM1,JM1,KM)
!     REAL(KIND=r8) :: TPI(IM,JM)
!     REAL(KIND=r8) :: TU(IM1,JM1,KM)
!     REAL(KIND=r8) :: TV(IM1,JM1,KM)
!     REAL(KIND=r8) :: TT(IM,JM,KM)

!     REAL(KIND=r8) :: NS(IM,JM)
!     REAL(KIND=r8) :: M(JM)
!     REAL(KIND=r8) :: MI(JM)

!     REAL(KIND=r8) :: FLAT(JM)
!     REAL(KIND=r8) :: FLON(IM)

!     REAL(KIND=r8) :: KSE(IM,JM)
!     REAL(KIND=r8) :: KSU(IM1,JM1)
!     REAL(KIND=r8) :: KSV(IM1,JM1)

!     REAL(KIND=r8) :: EK(KM)
!     REAL(KIND=r8) :: EKI(KIM)
!     REAL(KIND=r8) :: PHIS(IM,JM)
!     REAL(KIND=r8) :: PRF(IM,JM)
!     REAL(KIND=r8) :: PS(IM,JM)
!     REAL(KIND=r8) :: EP(IM,JM,KIM)
!     REAL(KIND=r8) :: WP(IM,JM,KM)
!     REAL(KIND=r8) :: PHI(IM,JM,KIM)
!     REAL(KIND=r8) :: PHIC(IM,JM,KIM)
!     REAL(KIND=r8) :: GPX(IM1,JM1,KM)
!     REAL(KIND=r8) :: GPY(IM1,JM1,KM)                   
!     REAL(KIND=r8) :: WW(JM)
!     REAL(KIND=r8) :: Q(IM,JM,KM)
!     REAL(KIND=r8) :: ALFA(IM,JM,KM)
 
 CALL Init() 
 CALL Run()  
 CALL Finalize()  

CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 SUBROUTINE Init()  
  IMPLICIT NONE

 END SUBROUTINE Init  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 SUBROUTINE Run()  
  IMPLICIT NONE
  CALL eta()
 END SUBROUTINE Run  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 SUBROUTINE Finalize()  
  IMPLICIT NONE

 END SUBROUTINE Finalize  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

END PROGRAM Main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
