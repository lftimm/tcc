!******************************************************************************C
!*                                                                            *C
!*      SIMULACAO NUMERICA DE ESCOAMENTOS INCOMPRESSIVEIS E PROBLEMAS DE      *C
!*    INTERAÇÃO FLUIDO-ESTRUTRA COM CORPOS IMERSOS INDESLOCÁVEIS ATRAVÉS DO   *C
!*                        METODO DOS ELEMENTOS FINITOS                        *C
!*                                                                            *C
!*    Características do modelo:                                              *C
!*                                                                            *C
!*                                                                            *C
!*  - Modelo numerico baseado no esquema de Taylor-Galerkin                   *C
!*    * Discretização temporal por Series de Taylor (até termos de segunda    *C
!*      ordem - Tensor de balanco difusivo) num processo explicito. A equação *C
!*      de conservacao de massa e obtida a partir da sua forma compressivel   *C
!*      com o auxilio do principio da pseudo-compressibilidade (Chorin,1967), *C
!*      produzindo uma equacao onde a pressao aparece explicitamente.         *C
!*    * Discretização espacial pelo Metodo dos Elementos Finitos usando       *C
!*      elementos hexaedricos de oito nos com integração reduzida e controle  *C
!*      de modos espurios ("hourglass").                                      *C
!*                                                                            *C
!******************************************************************************C
!*                                                                            *C
!*                             ALEXANDRE LUIS BRAUN                           *C
!*                          DOUTORANDO - PPGEC - UFRGS                        *C
!*                                                                            *C
!******************************************************************************C

!    ---------------------------------------------------------------------------
!    ------ BLOCO "A" ------ BLOCO "A" ------ BLOCO "A" ------ BLOCO "A" -------
!                        || DECLARACAO E DIMENSIONAMENTO ||
!    ---------------------------------------------------------------------------

USE DFPORT                
IMPLICIT real(8) (A-H,O-Z)
IMPLICIT INTEGER*4 (I-N)

! ***************************************************************************
!       PARAMETROS DE DIMENSIONAMENTO DO PROBLEMA A SER RESOLVIDO
! ***************************************************************************


 PARAMETER (NEMAX=16,NEMAP=16,NNOS=50,NNOP=50,&
            MNOBJ=0,NNMAX=0,NCMAX=0,NNCSX=0,NFCSX=0,&
            NBU=0,NBV=0,NBW=0,NBP=0)


!PARAMETER (NEMAX=6801,NEMAP=6800,NNOS=13981,NNOP=13980,&
!           MNOBJ=1,NNMAX=82,NCMAX=40,NNCSX=320,NFCSX=160,&
!           NBU=682,NBV=682,NBW=13980,NBP=82)


!   Problemas Tacoma
!PARAMETER (NEMAX=473001,NEMAP=473000,NNOS=495363,NNOP=495362,&
!           MNOBJ=1,NNMAX=0,NCMAX=0,NNCSX=10742,NFCSX=10480,&
!           NBU=18860,NBV=18860,NBW=42104,NBP=2296)


! Problema Multibuilding
!PARAMETER (NEMAX=1268273,NEMAP=1268272,NNOS=1318611,NNOP=1318610,&
!           MNOBJ=7,NNMAX=0,NCMAX=0,NNCSX=6270,NFCSX=6202,&
!           NBU=50110,NBV=67134,NBW=50110,NBP=9234)


! Problema Ventil
!PARAMETER (NEMAX=1037661,NEMAP=1037660,NNOS=1073919,NNOP=1073920,&
!           MNOBJ=0,NNMAX=0,NCMAX=0,NNCSX=0,NFCSX=0,&
!           NBU=36759,NBV=51009,NBW=36759,NBP=7068)

! Problema building2
!PARAMETER (NEMAX=752641,NEMAP=752640,NNOS=781519,NNOP=781520,&
!           MNOBJ=2,NNMAX=0,NCMAX=0,NNCSX=4247,NFCSX=4200,&
!           NBU=26179,NBV=37219,NBW=26179,NBP=7503)


!***********************************************************************************
!------>> OBSERVAÇÃO IMPORTANTE:
!         ----------------------
! Quando houver mais que um corpo imerso, as declarações NFCSX e NNCSX
! deverão tomar, respectivamente, o maior valor de faces de contorno
! sólido e o maior valor de nós de contorno sólido dentre os corpos existentes.
!***********************************************************************************  

!***********************************************************************************
!                             DEFINICAO DOS PARAMETROS
!
! NEMAX:   Numero de elementos da malha;
! NEMAP:   Numero (par) de elementos para dimensionamento, para que nao haja
!          arrays de dimensao multipla de 8;
! NNOS:    Numero de nos da malha;
! NNOP:    Nunero (par) de nos para dimensionamento, para que nao haja
!          arrays de dimensao multipla de 8;
! NNMAX:   Numero de nos com condicoes de contorno de 3a. especie (Simetria 
!          ou Outflow);
! NCMAX:   Numero de faces de elementos com condicoes de contorno de
!          3a.especie;
! NFCSX:   Numero de faces de contorno solido;
! NNCSX:   Número de nós de contorno sólido;
! NBU:     Numero de nos com condicao de contorno para a 1a componente de velocidade;
! NBV:     Numero de nos com condicao de contorno para a 2a componente de velocidade;
! NBW:     Numero de nos com condicao de contorno para a 3a componente de velocidade;
! NBP:     Numero de nos com condicao de contorno para a pressão;
! MNOBJ:   Número de objetos imersos no escoamento.
!************************************************************************************

real(8) ::  elapsed_time

real(8) ::  J11,J12,J13,J21,J22,J23,J31,J32,J33
      
real(8) ::  JCN11,JCN12,JCN13,JCN14,JCN15,JCN16,JCN17,JCN18,&
            JCN21,JCN22,JCN23,JCN24,JCN25,JCN26,JCN27,JCN28,&
            JCN31,JCN32,JCN33,JCN34,JCN35,JCN36,JCN37,JCN38

real(8) ::  JCN(8,3)

real(8) ::  JINV11,JINV12,JINV13,JINV21,JINV22,JINV23,JINV31,JINV32,JINV33
	            
real(8) ::  JB11,JB12,JB13,JB21,JB22,JB23,JB31,JB32,JB33,LRUNP,LRVNP,LRWNP,LRPNP,MRP

real(8) ::  DistB(12),uu(12),vv(12),ww(12),fs(12),s11(12),s12(12),s13(12),s22(12),s23(12),s33(12)


real(8) ::  MMM(8,8),MLUM(NNOS)

real(8) ::  JIN(9,NEMAX)

real(8) ::  Mx,My,Mz,Lchar(MNOBJ),Dchar(MNOBJ),Xobj(MNOBJ),Yobj(MNOBJ),Zobj(MNOBJ)

DIMENSION TAO(6,NFCSX),T11ARG(0:NEMAP),T22ARG(0:NEMAP),T33ARG(0:NEMAP),T12ARG(0:NEMAP),&
          T23ARG(0:NEMAP),T13ARG(0:NEMAP),VTV1(NNCSX),VTV2(NNCSX),VTV3(NNCSX),TTT1(NNCSX),TTT2(NNCSX),&
          TTT3(NNCSX),ARSuav(MNOBJ,NNCSX),ArAux(0:NEMAX),CF1(NNCSX),CF2(NNCSX),CF3(NNCSX),CPsuav(NNCSX),&
          CP(NNCSX),VTT1(NNCSX),VTT2(NNCSX),VTT3(NNCSX)
    
CHARACTER PNAME*6,PTNAME*14
     
DIMENSION BCU(NBU),BCV(NBV),BCW(NBW),BCP(NBP),IBCV(NBV),IBCU(NBU),IBCW(NBW),IBCP(NBP)

DIMENSION KONE(NEMAX*8+1),NEIBOR(0:16,NNOS)

! DIMENSION NEIBOL(0:16,NNOS)
    
DIMENSION BB1(8,8),BB2(8,8),BB3(8,8),DD11(8,8),DD12(8,8),DD13(8,8),&
          DD21(8,8),DD22(8,8),DD23(8,8),DD31(8,8),DD32(8,8),DD33(8,8),&
          EXP2(8,8),EXP1(8,8),EXP3(8,8),DIFMA(8,8),HOURG(8,8)

real(8) :: jbb(3,3),jbt(3,3),add(3,3),dd(3,3),uvwp(8,4),uvwt(8,4)

DIMENSION DZET(8),ETA(8),XI(8),SIGMA1(8),SIGMA2(8),SIGMA3(8),SIGMA4(8)

DIMENSION D11(8,8),D12(8,8),D13(8,8),D22(8,8),D23(8,8),D33(8,8),B123(24,8)

DIMENSION pp(9,NEMAX),uvwprom(3,NEMAX),PNV(NEMAX),DETERJ(NEMAX)

DIMENSION VDTurb(NEMAX),VtSuav(NNOS),Vtaux(0:NEMAX)

DIMENSION xyz(3,NNOS),PRSuav(NNOS),uvw(4,NNOS),crug(nnos),crvg(nnos),crwg(nnos),rrg(4,nnos)

real(8) :: xyzt(8,3),nk(8)
          
! *****************************************************************************************************
! Todos os vetores e matrizes que serao montados com ajuda das matrizes NEIBOR e NEIBOL devem 
! ser dimensionados levando-se em conta que a 1ra posicao corresponde ao subindice zero, uma 
! vez que os valores de NEIBOR e NEIBOL podem ser zero.
! Todos estes "arrays" devem  ser zerados na posicao correspondente ao subindice zero( VER BLOCO "C" ).  
! ***************************************************************************************************** 

dimension rru(8),rrv(8),rrw(8),rrp(8)

DIMENSION PPROME(0:NEMAP),VOLU(0:NEMAP)

DIMENSION CSNOR1(MNOBJ,NNCSX),CSNOR2(MNOBJ,NNCSX),CSNOR3(MNOBJ,NNCSX),CNORM(3,NNOS),cnort(3,8),&
          AREAC(NCMAX),IEL(NCMAX),FCONTOR(NCMAX*8+1)

DIMENSION KCONTS(MNOBJ,NFCSX*8+1),IELCS(MNOBJ,NFCSX),AREACS(MNOBJ,NFCSX)

DIMENSION IBNOBJ(MNOBJ,NNCSX),NFCN(MNOBJ),NNCN(MNOBJ)

DIMENSION S(6,NNOS),Smod(NNOS),FILTR2(NNOS),Cxt(NNOS),ND(8),NDSUP(0:12,NNOS),uvwf(3,NNOS),&
               SmodE(NEMAX),FILTR1(NEMAX),FILTS1(NNOS),dudx(NNOS),dudy(NNOS),&
               dudz(NNOS),dvdx(NNOS),dvdy(NNOS),dvdz(NNOS),dwdx(NNOS),dwdy(NNOS),dwdz(NNOS)

! DIMENSION dudxf(NNOS),dudyf(NNOS),dudzf(NNOS),dvdxf(NNOS),dvdyf(NNOS),&
!           dvdzf(NNOS),dwdxf(NNOS),dwdyf(NNOS),dwdzf(NNOS)


! DIMENSION SL11(NNOS),SL22(NNOS),SL33(NNOS),SL12(NNOS),SL13(NNOS),SL23(NNOS)


real(8) ::  SL11,SL22,SL33,SL12,SL13,SL23

DIMENSION JND(3)

real(8) ::  NDDIST(12,NNOS)

DIMENSION Umed(NNOS),Vmed(NNOS),Wmed(NNOS),PRmed(NNOS),PRSmed(NNOS),VTmed(NNOS)

DIMENSION NDTYPE(NNOS)
     
! DIMENSION FFFpG(NNOS),DP(NNOS)

! Arrays do novo modelo:

DIMENSION adv(8,8),TB1(8,8),TB2(8,8),TB3(8,8),ag1(8,8),ag2(8,8),ag3(8,8),&
          btd(8,8),bte(8,8),btn(8,8),dfp(8,8),CRU(8),CRV(8),CRW(8),&
          UVWa(4,NNOS),PRi(NNOS),PRf(NNOS),Gr1P(NEMAX),Gr2P(NEMAX),Gr3P(NEMAX)

real(8) :: MASS(8,8),masl(8,8)

real(8) :: MASD


DATA XI/-1.D0,1.D0,1.D0,-1.D0,-1.D0,1.D0,1.D0,-1.D0/, &
     ETA/-1.D0,-1.D0,1.D0,1.D0,-1.D0,-1.D0,1.D0,1.D0/, &
     DZET/-1.D0,-1.D0,-1.D0,-1.D0,1.D0,1.D0,1.D0,1.D0/ 
     
DATA SIGMA1/1.D0,1.D0,-1.D0,-1.D0,-1.D0,-1.D0,1.D0,1.D0/, &
     SIGMA2/1.D0,-1.D0,-1.D0,1.D0,-1.D0,1.D0,1.D0,-1.D0/, &
     SIGMA3/1.D0,-1.D0,1.D0,-1.D0,1.D0,-1.D0,1.D0,-1.D0/, &
     SIGMA4/-1.D0,1.D0,-1.D0,1.D0,1.D0,-1.D0,1.D0,-1.D0/
                             
! --------------------------------------------------------------------------- 
! ------ BLOCO "B" ------ BLOCO "B" ------ BLOCO "B" ------ BLOCO "B" -------
! ---------------------------------------------------------------------------
!                        || LEITURA DE DADOS ||
! ---------------------------------------------------------------------------
     
! DADOS DA MALHA - contornos solidos e contornos de Neumann 
! -----------------------------------------------------------

! NORCON    - número de nós com condição de contorno de Neumann  
! LCONTCARA -	número de faces de elementos com condição de contorno de Neumann 
! NFCS      -	número de faces de elementos de contorno sólido
! NNCS      -	número de nós de contorno sólido
! NOBJ      - número de objetos sólidos imersos

elapsed_time = TIMEF()

NORCON=NNMAX
LCONTCARA=NCMAX
NFCS=NFCSX
NNCS=NNCSX
NOBJ=MNOBJ


! ARQUIVO "inicia.par" - PARAMETROS DE inicialização 
! ----------------------------------------------------
      
OPEN(2,FILE='inicia.par',STATUS='OLD')
  REWIND (2)

  READ(2,*) PNAME
  READ(2,*) NCOEF 
  READ(2,*) NIR
  READ(2,*) NTR
  READ(2,*) NFILE
  READ(2,*) DtMAX
  READ(2,*) TPOAC
  READ(2,*) TOLTPO
  READ(2,*) NROTPO
  READ(2,*) CSEGUR
  READ(2,*) CONTROL1
  READ(2,*) NPASS
  READ(2,*) INDTURB
  READ(2,*) ELUMP
     
CLOSE(2)

! ***************************************************************************
!                          DEFINICAO DOS PARAMETROS

! PNAME:    Nome do projeto
! NCOEF:    Nro de intervalos entre registros de coefs. aerodinamicos;
! NIR:      Nro de intervalos entre registros;
! NTR:      Nro total de registros; 
! NFILE:    Nro do ultimo registro;
! DtMAX:    Estimativa inicial para o incremento de tempo;
! TPOAC:    Tempo atual;
! TOLTPO:   Tolerancia para o residuo ( Termino da Simulacao );
! NROTP:    Numero de passos de tempo apos atingir a tolerancia para o 
!           residuo para considerar-se o estado estacionario alcancado;
! CSEGUR:   Coeficiente de seguranca para determinar o incremento de tempo;
! CONTROL1  Parametro de controle para zerar ou nao alguma variavel. Pode
!           assumir os valores 0.0D+00 ou 1.0D+00.
! NPASS:    Número de passos de tempo a partir do qual serão calculados os
!           campos médios de velocidade e pressão. NPASSMED = NTR*NIR - NPASS
!           será o valor tomado para o cálculo da média. Para o caso em que
!           não se deseje este cálculo, usa-se NPASS = NTR*NIR
! INDTURB:  Indica se será ou não usado algum modelo de turbulência na
!           análise e, caso seja, faz a escolha entre os modelos CLÁSSICO e
!           DINÂMICO de SMAGORINSKY para as sub-escalas com simulação direta
!           para as grandes escalas, de acordo com os seguintes índices:
!   INDTURB = 0  --->  análise sem a presença da turbulência;
!   INDTURB = 1  --->  análise com turbulência e modelo sub-malha CLÁSSICO;
!   INDTURB = 2  --->  análise com turbulência e modelo sub-malha DINÂMICO; 
! ELUMP:    parametro seletivo de massa (normalmente adotado 0,9)
! ***************************************************************************

NNAUX = NTR*NIR
NPASSMED = NNAUX - NPASS
PASSMED = DFLOAT(NPASSMED)      

	
! ***************************************************************************
!   Atribui-se a variavel "PNAME" aos 6 1ros caracteres da variavel "PTNAME" 
! ***************************************************************************
     
PTNAME(1:6)=PNAME


!        ARQUIVO '.pro' - PROPRIEDADES DO FLUIDO, DO FLUXO E DIMENSÕES
!       ---------------------------------------------------------------
!                 CARACTERÍSTICAS DO(S) CORPO(S) IMERSO(S)
!                ------------------------------------------

! ***************************************************************************
!  Atribui-se a extensao ".pro" aos 4 ultimos caracteres da variavel "PTNAME"
! ***************************************************************************

 PTNAME(7:14)='.pro    '
      
OPEN(2,FILE=PTNAME,STATUS='OLD')
  
  REWIND (2)
  READ(2,*) VInf
  READ(2,*) VelSom
  READ(2,*) ViscCin
  READ(2,*) ViscVol
  READ(2,*) PInf
  READ(2,*) RHOInf

  IF(INDTURB.EQ.1) READ(2,*) Cs
    
  IF(MNOBJ.NE.0) THEN
      
    DO NOB=1,NOBJ

      READ(2,*) Lchar(NOB)
      READ(2,*) Dchar(NOB)
      READ(2,*) Xobj(NOB)
      READ(2,*) Yobj(NOB)
      READ(2,*) Zobj(NOB)
      READ(2,*) NFCN(NOB),NNCN(NOB)

    END DO

  END IF
            
CLOSE(2)

! Pdin é a pressão dinâmica baseada na velocidade de corrente não perturbada (Vinf):

Pdin = 0.5D0*RHOInf*Vinf*Vinf

! ***************************************************************************
!                          DEFINICAO DAS PROPRIEDADES
!
! Vinf:      Modulo do vetor velocidade da corrente nao perturbada (Infinito).
! VelSom:    Velocidade de propagação do som;
! ViscCin:   Viscosidade cinematica do fluido [m^2/s];
! ViscVol:   Viscosidade Volumetrica por unidade de massa especifica do fluido [m^2/s];
! Pinf:      Pressao da corrente nao perturbada (Infinito); 
! RHOInf:    Densidade do fluido;
! Cs:        Coeficiente de Smagorinsky (valores usuais de 0.1 a 0.22).
!            Válido somente para o modelo clássico de Smagorinsky (INDTURB=1).
! Lchar:     Comprimento característico do corpo imerso NOB;
! Dchar:     Dimensão característica do corpo imerso NOB
! (X,Y,Z)obj: Coordenadas do centro de massa do corpo imerso NOB
! NFCN,NNCN: Número de faces e de nós de contorno do corpo NOB   
! ***************************************************************************

!            ARQUIVO '.con' - CONETIVIDADES DOS ELEMENTOS
!            --------------------------------------------

! ***************************************************************************
! Atribui-se a extensao ".con" aos 4 ultimos caracteres da variavel "PTNAME"
! ***************************************************************************

PTNAME(7:14)='.con    '
OPEN (2,FILE=PTNAME,STATUS='OLD')

  REWIND (2)
  READ(2,*) NELEM      
  DO I=1,NELEM                                                        
    II=8*(I-1)
    READ(2,*) KONE(II+1),KONE(II+2),KONE(II+3),KONE(II+4),&
              KONE(II+5),KONE(II+6),KONE(II+7),KONE(II+8) 
  END DO
    
CLOSE(2)

!                   ARQUIVO '.cor'  - COORDENADAS DOS NOS 
!                   -------------------------------------

! ***************************************************************************
! Atribui-se  "00.cor" aos 6 ultimos caracteres da variavel "PTNAME"
! ***************************************************************************

PTNAME(7:9)='000'
PTNAME(10:14)='.cor '

OPEN(2,FILE=PTNAME,STATUS='OLD')

  REWIND (2)
  READ(2,*) NNM
  DO I=1,NNM
    READ(2,*) xyz(1,I),xyz(2,I),xyz(3,I)
  END DO
      
CLOSE(2)

!       ARQUIVO 'FF.sup' - AREAS E CONETIVIDADES DAS FACES DE CONTORNO  
!       --------------------------------------------------------------
!                            COM FORCAS DE SUPERFICIE
!                            ------------------------   

! ***************************************************************************
! Atribui-se "FF.sup" aos 6 ultimos caracteres da variavel "PTNAME"
! ***************************************************************************
     
PTNAME(7:14)='FF.sup  '

OPEN (2,FILE=PTNAME,STATUS='OLD')

  REWIND (2)
  DO I=1,LCONTCARA                                                        
    READ(2,*) IEL(I)
    II=8*(I-1)
    READ(2,*) AREAC(I),FCONTOR(II+1),FCONTOR(II+2),FCONTOR(II+3),&
              FCONTOR(II+4),FCONTOR(II+5),FCONTOR(II+6),FCONTOR(II+7),FCONTOR(II+8) 
  END DO
      
CLOSE(2)

!          ARQUIVO 'FF.nnn' - VETORES NORMAIS AOS NOS DE CONTORNO COM       
!          ---------------    ----------------------------------------
!                            CONDICOES DE NEUMANN
!                           ----------------------

! ***************************************************************************
! Atribui-se  "FF.nnn" aos 6 ultimos caracteres da variavel "PTNAME"
! ***************************************************************************

cnorm = 0.0D0
PTNAME(7:14)='FF.nnn  '

OPEN(2,FILE=PTNAME,STATUS='OLD')

  DO I=1,NORCON
    READ (2,*) NNOCON,cnorm(1,NNOCON),cnorm(2,NNOCON),cnorm(3,NNOCON)
  END DO

CLOSE(2)

IF(MNOBJ.NE.0) THEN

  DO NOB=1,NOBJ

    NB1 = NOB/10 + 48
    NB2 = NOB - (NOB/10)*10 + 48

!         ARQUIVO 'CS##.sup' - AREAS E CONETIVIDADES DAS FACES DE CONTORNO
!        ------------------------------------------------------------------ 
!                                      SÓLIDO
!                                     --------

!   ***************************************************************************
!   Atribui-se "CS##.sup" aos 8 últimos caracteres da variavel "PTNAME"
!   ***************************************************************************

    PTNAME(7:8)='CS'
    PTNAME(9:10)=CHAR(NB1)//CHAR(NB2)
    PTNAME(11:14)='.sup'

    OPEN (2,FILE=PTNAME,STATUS='OLD')

      REWIND (2)
      II=8*NFCS
      DO I=1,II
        KCONTS(NOB,I)=0
      END DO
      DO I=1,NFCS
        IELCS(NOB,I)=0
        AREACS(NOB,I)=0.0D0
      END DO
      DO I=1,NFCS
        IF(I.GT.NFCN(NOB)) EXIT
        READ(2,*) IELCS(NOB,I)
        II=8*(I-1)
        READ(2,*) AREACS(NOB,I),KCONTS(NOB,II+1),KCONTS(NOB,II+2),&
                  KCONTS(NOB,II+3),KCONTS(NOB,II+4),KCONTS(NOB,II+5),&
                  KCONTS(NOB,II+6),KCONTS(NOB,II+7),KCONTS(NOB,II+8)
      END DO

    CLOSE(2)

!         ARQUIVO 'CS##.nnn' -  VETORES NORMAIS AOS NOS DE CONTORNO SÓLIDO      
!        -------------------   --------------------------------------------

!   ***************************************************************************
!   Atribui-se  "CS##.nnn" aos 8 últimos caracteres da variavel "PTNAME"
!   ***************************************************************************

    PTNAME(11:14)='.nnn'

    CSNOR1(NOB,:)=0.0D0
    CSNOR2(NOB,:)=0.0D0
    CSNOR3(NOB,:)=0.0D0
    IBNOBJ(NOB,:)=0


    OPEN(2,FILE=PTNAME,STATUS='OLD')

      DO I=1,NNCS
        IF(I.GT.NNCN(NOB)) EXIT
        READ (2,*) NNOCON,CSNOR1(NOB,I),CSNOR2(NOB,I),CSNOR3(NOB,I)
        IBNOBJ(NOB,I)=NNOCON
      END DO

    CLOSE(2)

   END DO

END IF


! NUMERACAO DOS ARQUIVOS DE ACORDO COM O NUMERO DE REGISTRO CORRENTE 
! ------------------------------------------------------------------

NF1=NFILE/100+48
NF2=(NFILE-(NFILE/100)*100)/10+48
NF3=NFILE-(NFILE/10)*10+48
PTNAME(7:9)=CHAR(NF1)//CHAR(NF2)//CHAR(NF3)

! ****************************************************************************
! Esta e a forma de atribuir as setima e oitava posicoes dos caracteres
! do nome dos arquivos de um determinado registro, exatamente o numero desse
! registro. A dupla barra e a forma sintatica de indicar e separar cada um 
! dos caracteres. ( Ver no "Help" a funcao  CHAR.)     
! ***************************************************************************

!              ARQUIVOS COM OS VALORES NODAIS DAS VARIAVEIS FÍSICAS        
!           ----------------------------------------------------------

!                 ARQUIVO '.v': VALORES DE "U" "V" E "W" NOS NÓS 
!               --------------------------------------------------
!                 ARQUIVO '.pr': VALORES DE "PR" E PRSuav NOS NÓS
!               ---------------------------------------------------

! ***************************************************************************
! Atribui-se  aos 3 últimos caracteres da variavel "PTNAME" a extensao
! equivalente a variavel fisica armazenada em cada arquivo
! ***************************************************************************

PTNAME(10:14)='.v   '

OPEN(2,FILE=PTNAME,STATUS='OLD')

  REWIND (2)
  DO I=1,NNM
    READ(2,*)  uvw(1,I),uvw(2,I),uvw(3,I)
  END DO
      
CLOSE(2)

PTNAME(10:14)='.pr  '

OPEN(3,FILE=PTNAME,STATUS='OLD')

  REWIND (3)
  DO I=1,NNM
    READ(3,*)  uvw(4,I),PRSuav(I)
  END DO

CLOSE(3)

! ***************************************************************************



!            ARQUIVOS COM AS CONDICOES DE CONTORNO DAS VARIAVIES FISICAS 
!           -------------------------------------------------------------

!              ARQUIVO 'CC.bv': CONDICOES DE CONTORNO PARA "U" "V" E "W" 
!              ---------------------------------------------------------
!                  ARQUIVO 'CC.bp': CONDICOES DE CONTORNO PARA "PR"
!                  ------------------------------------------------     

! ***************************************************************************
! Atribui-se  aos  ultimos caracteres da variavel "PTNAME" os caracteres
! "CC" e a extensao equivalente a variavel fisica armazenada em cada arquivo
! ***************************************************************************

PTNAME(7:14) ='CC.bv   '

OPEN(2,FILE=PTNAME,STATUS='OLD')
 
  REWIND (2)


! ***************************************************************************
! NBCU e a quantidade de nos com condicoes de contorno (neste caso -> U).
! IBCU(I) e o Nro do no global respectivo ao no de contorno corrente "I".
! BCU(I)  e o valor da variavel prescrita (neste caso -> U).
! ***************************************************************************

  READ(2,*) NBCU
  IF (NBCU.NE.0) THEN
    DO I=1,NBCU
      READ(2,*) IBCU(I),BCU(I)
    END DO
  END IF

  READ(2,*) NBCV
  IF (NBCV.NE.0) THEN
    DO I=1,NBCV
      READ(2,*) IBCV(I),BCV(I)
    END DO
  END IF

  READ(2,*) NBCW
  IF (NBCW.NE.0) THEN
    DO I=1,NBCW
      READ(2,*) IBCW(I),BCW(I)
    END DO
  END IF

CLOSE(2)

PTNAME(7:14) ='CC.bp   '

OPEN(3,FILE=PTNAME,STATUS='OLD')

  REWIND (3)

  READ(3,*) NBCP
  IF (NBCP.NE.0) THEN
    DO I=1,NBCP
      READ(3,*) IBCP(I),BCP(I)
    END DO
  END IF

CLOSE(3)


! ***************************************************************************
! Formação do vetor NDTYPE(NNOS) que qualifica cada um dos nós da malha 
! conforme a existência ou não de condições de contormo a eles impostas:
!
! NDTYPE(I) =  0  ---> nó sem restrições para a velocidade
! NDTYPE(I) =  1  ---> nó com restrição na direção X da velocidade
! NDTYPE(I) =  2  ---> nó com restrição na direção Y da velocidade
! NDTYPE(I) =  3  ---> nó com restrição na direção Z da velocidade
! NDTYPE(I) =  12 ---> nó com restrições nas direções X e Y da velocidade
! NDTYPE(I) =  13 ---> nó com restrições nas direções X e Z da velocidade
! NDTYPE(I) =  23 ---> nó com restrições nas direções Y e Z da velocidade
! NDTYPE(I) = 123 ---> nó com restrições nas direções X, Y e Z da velocidade
!      
! ***************************************************************************

DO I=1,NNM

  JJ=0
  KK=0
  LL=0

  DO J=1,NBCU
    IF(I.EQ.IBCU(J)) THEN
      JJ=I
      EXIT
    END IF
  END DO

  DO KI=1,NBCV
    IF(I.EQ.IBCV(KI)) THEN
      KK=I
      EXIT
    END IF
  END DO

  DO L=1,NBCW
    IF(I.EQ.IBCW(L)) THEN
      LL=I
      EXIT
    END IF
  END DO

  IF(((I.EQ.JJ).AND.(I.EQ.KK).AND.(I.EQ.LL))) NDTYPE(I)=123
  IF(((I.EQ.JJ).AND.(I.EQ.KK).AND.(I.NE.LL))) NDTYPE(I)=12
  IF(((I.EQ.JJ).AND.(I.NE.KK).AND.(I.EQ.LL))) NDTYPE(I)=13
  IF(((I.NE.JJ).AND.(I.EQ.KK).AND.(I.EQ.LL))) NDTYPE(I)=23
  IF(((I.EQ.JJ).AND.(I.NE.KK).AND.(I.NE.LL))) NDTYPE(I)=1
  IF(((I.NE.JJ).AND.(I.EQ.KK).AND.(I.NE.LL))) NDTYPE(I)=2
  IF(((I.NE.JJ).AND.(I.NE.KK).AND.(I.EQ.LL))) NDTYPE(I)=3
  IF(((I.NE.JJ).AND.(I.NE.KK).AND.(I.NE.LL))) NDTYPE(I)=0

END DO



! CONDICOES DE CONTORNO INICIAIS PARA VELOCIDADE E PRESSÃO
! ----------------------------------------------------------

IF (NBCU.NE.0) THEN
  DO I= 1,NBCU
    NNO= IBCU(I)
    uvw(1,NNO)= BCU(I)
  END DO
END IF

IF (NBCV.NE.0) THEN
  DO I= 1,NBCV
    NNO= IBCV(I)
    uvw(2,NNO)= BCV(I)
  END DO
END IF

IF (NBCW.NE.0) THEN
  DO I= 1,NBCW
    NNO= IBCW(I)
    uvw(3,NNO)= BCW(I)
  END DO
END IF

IF (NBCP.NE.0) THEN
  DO I= 1,NBCP
    NNO= IBCP(I)
    uvw(4,NNO)= BCP(I)
  END DO
END IF


!***** ZERAMENTO DOS VETORES CONTENDO AS VARIÁVEIS DE CAMPO MÉDIO:

IF(NPASS.NE.NNAUX) THEN

  Umed = 0.0D0
  Vmed = 0.0D0
  Wmed = 0.0D0
  PRmed = 0.0D0
  VTmed = 0.0D0
  PRSmed = 0.0D0

END IF



! ---------------------------------------------------------------------------
! ------ BLOCO "C" ------ BLOCO "C" ------ BLOCO "C" ------ BLOCO "C" -------
! ---------------------------------------------------------------------------
!                    || MONTAJEM DE NEIBOL E NEIBOR ||
! ---------------------------------------------------------------------------

! As Matrizes NEIBOR e NEIBOL servem para o processo de montagem dos vetores
! Globais. Esta montagem deve dar-se em forma vetorizada para aumentar a
! eficiencia do programa portanto, por isso utilizam-se estes "Arrays".
!    
! IPOS e um indice para correlacionar as correspondentes posicoes das matrizes 
!     [ Por exemplo NEIBOR(3,NODE) se corresponde com NEIBOL(3,NODE)].
! NODE indica o Numero de No Global.
! 
! NEIBOR(IPOS,NODE): *na sua primeira posicao, NEIBOR(0,NODE),contem a
!                     quantidade de elementos concorrentes ao no global
!                     "NODE". SE CONSIDERA QUE NAO PODE HAVER MAIS DE 16 
!                     ELEMENTOS CONCORRENTES AO MESMO NO (Valido apenas  
!                     para malhas estruturadas).
!                    *nas posicoes restantes, [ NEIBOR(1,NODE) a 
!                     NEIBOR(16,NODE) ], contem o numero de elemento. 
!                     IPOS vai de 1 a 16 pois no maximo pode haver 
!                     16 elementos concorrentes a um determinado no. 
!
! NEIBOL(IPOS,NODE): *na sua primeira posicao, NEIBOL(0,NODE),contem a
!                     quantidade de elementos concorrentes ao no global
!                     "NODE".
!                    *nas posicoes restantes, [ NEIBOL(1,NODE) a 
!                     NEIBOL(16,NODE) ], contem o numero local de nodo do
!                     elemento dado pela matriz NEIBOR, correspondente 
!                     ao no global NODE.


DO I=1,NNOS
  NEIBOR(0,I)=0
END DO

  
DO IP=1,8
  DO I=1,NELEM
    I1 = 8*I - 8 + IP
    NODE = KONE(I1)
    NEIBOR(0,NODE) = NEIBOR(0,NODE) + 1
!   NEIBOL(0,NODE) = NEIBOR(0,NODE)
    IPOS = NEIBOR(0,NODE)
    NEIBOR(IPOS,NODE) = I
!   NEIBOL(IPOS,NODE) = IP
  END DO 
END DO 


! FORMACAO DE MATRIZES ATRAVES DO PRODUTO DAS FUNCOES DE INTERPOLACAO E SUAS  DERIVADAS
! *************************************************************************************
     
do j=1,8

  DO i=1,8   

    MMM(I,J)= (1.D0+ETA(I)*ETA(J)/3.D0)*(1.D0+DZET(I)*DZET(J)/3.D0)*(1.D0+XI(I)*XI(J)/3.D0) 
      
    DIFMA(I,J)= -MMM(I,J)
        
    BB1(I,J)= XI(J)*(1.D0+ETA(I)*ETA(J)/3.D0)*(1.D0+DZET(I)*DZET(J)/3.D0)
    BB2(I,J)= ETA(J)*(1.D0+DZET(I)*DZET(J)/3.D0)*(1.D0+XI(I)*XI(J)/3.D0)
    BB3(I,J)= DZET(J)*(1.D0+XI(I)*XI(J)/3.D0)*(1.D0+ETA(I)*ETA(J)/3.D0)

    TB1(I,J)= XI(I)*(1.D0+ETA(I)*ETA(J)/3.D0)*(1.D0+DZET(I)*DZET(J) /3.D0)
    TB2(I,J)= ETA(I)*(1.D0+DZET(I)*DZET(J)/3.D0)*(1.D0+XI(I)*XI(J)/3.D0)
    TB3(I,J)= DZET(I)*(1.D0+XI(I)*XI(J)/3.D0)*(1.D0+ETA(I)*ETA(J) /3.D0)
       
    EXP3(I,J)= (1.D0+DZET(I)*DZET(J)/3.D0)
    EXP2(I,J)= (1.D0+ETA(I)*ETA(J)/3.D0)
    EXP1(I,J)= (1.D0+XI(I)*XI(J)/3.D0)

    DD11(I,J)= XI(I)*BB1(I,J)
    DD22(I,J)= ETA(I)*BB2(I,J)
    DD33(I,J)= DZET(I)*BB3(I,J)
    DD12(I,J)= XI(I)*ETA(J)*EXP3(I,J)
    DD13(I,J)= XI(I)*DZET(J)*EXP2(I,J)
    DD21(I,J)= XI(J)*ETA(I)*EXP3(I,J)
    DD23(I,J)= ETA(I)*DZET(J)*EXP1(I,J)
    DD31(I,J)= XI(J)*DZET(I)*EXP2(I,J)
    DD32(I,J)= ETA(J)*DZET(I)*EXP1(I,J)
   
  END DO
       
  DIFMA(j,j)= 8.D0 - MMM(j,j)

end do
	 
	 
! CALCULO DA MATRIZ UTILIZADA PARA CONTROLE DE HOURGLASSING
! *********************************************************

DO I=1, 8
  DO J=1, 8 
    HOURG(I,J)=SIGMA1(I)*SIGMA1(J)+SIGMA2(I)*SIGMA2(J)+SIGMA3(I)*SIGMA3(J)+SIGMA4(I)*SIGMA4(J)
  END DO
END DO

	
! COEFICENTES UTILIZADOS PARA CALCULO DOS TERMOS DIFUSIVOS E FORCAS VISCOSAS
! **************************************************************************
! (MATRIZES "Dij" E VETORES "Sj")
! *******************************

COEF1= 2.0D0*ViscCin+ViscVol
COEF2= ViscCin  
COEF3= ViscVol

! ZERAMENTO DOS VETORES E MATRIZES QUE SERAO MONTADOS USANDO AS MATRIZES 
! **********************************************************************
! NEIBOR E NEIBOL, NA POSICAO CORRESPONDENTE AO SUBINDICE ZERO
! ************************************************************

VOLU(0)  =0.0D0
PPROME(0)=0.0D0
Vtaux(0) =0.0D0

! VOLUMES MAXIMO E MINIMO E DELTA T 
! *********************************

VOLMIN= 1.0D6
VOLMAX= 0.0D0


DO IELEM=1,NELEM

  II= (IELEM-1)*8
       
  do i=1,8
    xyzt(i,1) = xyz(1,kone(ii+i))
    xyzt(i,2) = xyz(2,kone(ii+i))
    xyzt(i,3) = xyz(3,kone(ii+i))
  end do	 
       
! MATRIZ JACOBIANA NO PONTO:(XI=0,ETA=0,DZET=0) 
! -----------------------------------------------

  j11 = dot_product(xi,xyzt(:,1))*0.125d0
  j12 = dot_product(xi,xyzt(:,2))*0.125d0
  j13 = dot_product(xi,xyzt(:,3))*0.125d0

  j21 = dot_product(eta,xyzt(:,1))*0.125d0
  j22 = dot_product(eta,xyzt(:,2))*0.125d0
  j23 = dot_product(eta,xyzt(:,3))*0.125d0

  j31 = dot_product(dzet,xyzt(:,1))*0.125d0
  j32 = dot_product(dzet,xyzt(:,2))*0.125d0
  j33 = dot_product(dzet,xyzt(:,3))*0.125d0



! DETERMINANTE JACOBIANO NO PONTO:(XI=0,ETA=0,DZET=0) 
! ------------------------------------------------------ 

  DETERJ(IELEM)= J11*J22*J33-J11*J23*J32-J21*J12*J33+J21*J13*J32+J31*J12*J23-J31*J13*J22


! INVERSA DA MATRIZ JACOBIANA NO PONTO:(XI=0,ETA=0,DZET=0)
! ----------------------------------------------------------- 

! ***************************************************************************
!  -1                                                           
! J  = JB / DET J ... sendo JB a trasposta da matriz adjunta de J. 
! ~    ~        ~            ~                                  ~
! ***************************************************************************

  JB11= (J22*J33-J23*J32)
  JB12=-(J12*J33-J13*J32)
  JB13= (J12*J23-J13*J22)
  JB21=-(J21*J33-J23*J31)
  JB22= (J11*J33-J13*J31)
  JB23=-(J11*J23-J13*J21)
  JB31= (J21*J32-J22*J31)
  JB32=-(J11*J32-J12*J31)
  JB33= (J11*J22-J12*J21)

  jin(1,IELEM)= JB11/DETERJ(IELEM)
  jin(2,IELEM)= JB12/DETERJ(IELEM)
  jin(3,IELEM)= JB13/DETERJ(IELEM)
  jin(4,IELEM)= JB21/DETERJ(IELEM)
  jin(5,IELEM)= JB22/DETERJ(IELEM)
  jin(6,IELEM)= JB23/DETERJ(IELEM)
  jin(7,IELEM)= JB31/DETERJ(IELEM)
  jin(8,IELEM)= JB32/DETERJ(IELEM)
  jin(9,IELEM)= JB33/DETERJ(IELEM)

! ***************************************************************************
!                           VOLUME DO ELEMENTO
!                          --------------------
!               (Vetor que contem o volume de cada elemento) 
!
! ***************************************************************************

  VOLU(IELEM)= 8.D0*DETERJ(IELEM)

  IF (VOLU(IELEM).GT.VOLMAX) VOLMAX= VOLU(IELEM)
  IF (VOLU(IELEM).LT.VOLMIN) VOLMIN= VOLU(IELEM)

! ***************************************************************************
! Necessita-se dos volumes maximo e minimo para as formulas de controle de
! modos espureos (HOURGLASS)
! ***************************************************************************

  DL12=DSQRT((xyzt(2,1)-xyzt(1,1))*(xyzt(2,1)-xyzt(1,1))+(xyzt(2,2)-&
              xyzt(1,2))*(xyzt(2,2)-xyzt(1,2))+(xyzt(2,3)-xyzt(1,3))&
              *(xyzt(2,3)-xyzt(1,3)))
  DL23=DSQRT((xyzt(3,1)-xyzt(2,1))*(xyzt(3,1)-xyzt(2,1))+(xyzt(3,2)-&
              xyzt(2,2))*(xyzt(3,2)-xyzt(2,2))+(xyzt(3,3)-xyzt(2,3))&
              *(xyzt(3,3)-xyzt(2,3)))
  DL34=DSQRT((xyzt(4,1)-xyzt(3,1))*(xyzt(4,1)-xyzt(3,1))+(xyzt(4,2)-&
              xyzt(3,2))*(xyzt(4,2)-xyzt(3,2))+(xyzt(4,3)-xyzt(3,3))&
              *(xyzt(4,3)-xyzt(3,3)))
  DL41=DSQRT((xyzt(1,1)-xyzt(4,1))*(xyzt(1,1)-xyzt(4,1))+(xyzt(1,2)-&
              xyzt(4,2))*(xyzt(1,2)-xyzt(4,2))+(xyzt(1,3)-xyzt(4,3))&
              *(xyzt(1,3)-xyzt(4,3)))
  DL15=DSQRT((xyzt(5,1)-xyzt(1,1))*(xyzt(5,1)-xyzt(1,1))+(xyzt(5,2)-&
              xyzt(1,2))*(xyzt(5,2)-xyzt(1,2))+(xyzt(5,3)-xyzt(1,3))&
              *(xyzt(5,3)-xyzt(1,3)))
  DL26=DSQRT((xyzt(6,1)-xyzt(2,1))*(xyzt(6,1)-xyzt(2,1))+(xyzt(6,2)-&
              xyzt(2,2))*(xyzt(6,2)-xyzt(2,2))+(xyzt(6,3)-xyzt(2,3))&
              *(xyzt(6,3)-xyzt(2,3)))
  DL37=DSQRT((xyzt(7,1)-xyzt(3,1))*(xyzt(7,1)-xyzt(3,1))+(xyzt(7,2)-&
              xyzt(3,2))*(xyzt(7,2)-xyzt(3,2))+(xyzt(7,3)-xyzt(3,3))&
              *(xyzt(7,3)-xyzt(3,3)))
  DL48=DSQRT((xyzt(8,1)-xyzt(4,1))*(xyzt(8,1)-xyzt(4,1))+(xyzt(8,2)-&
              xyzt(4,2))*(xyzt(8,2)-xyzt(4,2))+(xyzt(8,3)-xyzt(4,3))&
              *(xyzt(8,3)-xyzt(4,3)))
  DL56=DSQRT((xyzt(6,1)-xyzt(5,1))*(xyzt(6,1)-xyzt(5,1))+(xyzt(6,2)-&
              xyzt(5,2))*(xyzt(6,2)-xyzt(5,2))+(xyzt(6,3)-xyzt(5,3))&
              *(xyzt(6,3)-xyzt(5,3)))
  DL67=DSQRT((xyzt(7,1)-xyzt(6,1))*(xyzt(7,1)-xyzt(6,1))+(xyzt(7,2)-&
              xyzt(6,2))*(xyzt(7,2)-xyzt(6,2))+(xyzt(7,3)-xyzt(6,3))&
              *(xyzt(7,3)-xyzt(6,3)))
  DL78=DSQRT((xyzt(8,1)-xyzt(7,1))*(xyzt(8,1)-xyzt(7,1))+(xyzt(8,2)-&
              xyzt(7,2))*(xyzt(8,2)-xyzt(7,2))+(xyzt(8,3)-xyzt(7,3))&
              *(xyzt(8,3)-xyzt(7,3)))
  DL85=DSQRT((xyzt(5,1)-xyzt(8,1))*(xyzt(5,1)-xyzt(8,1))+(xyzt(5,2)-&
              xyzt(8,2))*(xyzt(5,2)-xyzt(8,2))+(xyzt(5,3)-xyzt(8,3))&
              *(xyzt(5,3)-xyzt(8,3)))

  DELTAL = DMIN1( DL12,DL23,DL34,DL41,DL15,DL26,DL37,DL48,DL56,DL67,DL78,DL85 )

! CONDIÇÃO DE ESTABILIDADE PARA ESCOAMENTOS COM CONVECÇÃO DOMINANTE     
      
  DeltE=CSEGUR*DELTAL/(VelSom+Vinf)

! CONDIÇÃO DE ESTABILIDADE PARA ESCOAMENTOS COM DIFUSÃO DOMINANTE
       
! DeltE=CSEGUR*DELTAL**2.0D0*RHOinf/(2.0D0*ViscCin)

  IF (DeltE.LT.DtMAX) DtMAX=DeltE


END DO
 


! ************************

IF(IndTurb.EQ.2) THEN

! ************************

! SELEÇÃO DOS NÓS PARA COMPOSIÇÃO DOS FILTROS INDIVIDUAIS DE CADA UM DOS NÓS
! **************************************************************************
! EXISTENTES NA MALHA PARA O PROCESSO DE SEGUNDA FILTRAGEM DO MODELO DINÂMICO
! ***************************************************************************
! DE TURBULÊNCIA
! **************

! ***** Zeramento da primeira posição de arrays usados	no modelo dinâmico:

  NDSUP = 0
  NDDIST = 0.0D0

  DO I=1,NNM

    II=NEIBOR(0,I)
    VOLtot=0.0D0
    JCOUNT=0

!   ***** Parte do algoritmo referente a nós com apenas 1 elemento em comum:

    IF(II.EQ.1) THEN

      NLM=NEIBOR(1,I)

      VOLtot = VOLtot + VOLU(NLM)

      do j=1,8
	    nd(j) = kone(8*nlm-8+j)
	  end do

!     ***** Tomando o nó diametralmente oposto a I no mesmo plano:

      DO J=1,8

        JJ=ND(J)

        IF(JJ.EQ.I) THEN

          IF(J.LE.2) THEN

            JK=J+2
            JK2=J+4
            JK3=8*(NLM-1)+JK2
            JK4=KONE(JK3)

            JCOUNT=JCOUNT +1
            NDSUP(0,I)=JCOUNT
            NDSUP(JCOUNT,I)=JK4

            Dx = xyz(1,JK4) - xyz(1,I)
            Dy = xyz(2,JK4) - xyz(2,I)
            Dz = xyz(3,JK4) - xyz(3,I)

            NDDIST(JCOUNT,I) = DSQRT(Dx*Dx + Dy*Dy + Dz*Dz)

          END IF

          IF((J.GT.2).AND.(J.LE.4)) THEN

            JK=J-2
            JK2=J+4
            JK3=8*(NLM-1)+JK2
            JK4=KONE(JK3)

            JCOUNT=JCOUNT +1
            NDSUP(0,I)=JCOUNT
            NDSUP(JCOUNT,I)=JK4

            Dx = xyz(1,JK4) - xyz(1,I)
            Dy = xyz(2,JK4) - xyz(2,I)
            Dz = xyz(3,JK4) - xyz(3,I)

            NDDIST(JCOUNT,I) = DSQRT(Dx*Dx + Dy*Dy +  Dz*Dz)

          END IF

          IF((J.GT.4).AND.(J.LE.6)) THEN

            JK=J+2
            JK2=J-4
            JK3=8*(NLM-1)+JK2
            JK4=KONE(JK3)

            JCOUNT=JCOUNT +1
            NDSUP(0,I)=JCOUNT
            NDSUP(JCOUNT,I)=JK4

            Dx = xyz(1,JK4) - xyz(1,I)
            Dy = xyz(2,JK4) - xyz(2,I)
            Dz = xyz(3,JK4) - xyz(3,I)

            NDDIST(JCOUNT,I) = DSQRT(Dx*Dx + Dy*Dy +  Dz*Dz)

          end if

          IF(J.GT.6) THEN

            JK=J-2
            JK2=J-4
            JK3=8*(NLM-1)+JK2
            JK4=KONE(JK3)

            JCOUNT=JCOUNT +1
            NDSUP(0,I)=JCOUNT
            NDSUP(JCOUNT,I)=JK4

            Dx = xyz(1,JK4) - xyz(1,I)
            Dy = xyz(2,JK4) - xyz(2,I)
            Dz = xyz(3,JK4) - xyz(3,I)

            NDDIST(JCOUNT,I) = DSQRT(Dx*Dx + Dy*Dy +   Dz*Dz)

          END IF

          NI=8*(NLM-1)+JK
          IN=KONE(NI)

          EXIT
               
        END IF

      END DO


      II=NEIBOR(0,IN)


      DO J=1,II

        NEL=NEIBOR(J,IN)

        IF(NEL.NE.NLM) THEN

          VOLtot = VOLtot + VOLU(NEL)

          DO IJ=1,8

            JI=8*(NEL-1)+IJ
            JJ=KONE(JI)

            IF(JJ.EQ.IN) THEN

              IF(IJ.LE.2) THEN

                JK=IJ+2
                JL1=8*(NEL-1)+JK
                JL2=KONE(JL1)

                JCOUNT=JCOUNT +1
                NDSUP(0,I)=JCOUNT
                NDSUP(JCOUNT,I)=JL2

                Dx = xyz(1,JL2) - xyz(1,I)
                Dy = xyz(2,JL2) - xyz(2,I)
                Dz = xyz(3,JL2) - xyz(3,I)

                NDDIST(JCOUNT,I) = DSQRT(Dx*Dx + Dy*Dy + Dz*Dz)

              END IF

              IF((IJ.GT.2).AND.(IJ.LE.4)) THEN

                JK=IJ-2
                JL1=8*(NEL-1)+JK
                JL2=KONE(JL1)

                JCOUNT=JCOUNT +1
                NDSUP(0,I)=JCOUNT
                NDSUP(JCOUNT,I)=JL2

                Dx = xyz(1,JL2) - xyz(1,I)
                Dy = xyz(2,JL2) - xyz(2,I)
                Dz = xyz(3,JL2) - xyz(3,I)

                NDDIST(JCOUNT,I) = DSQRT(Dx*Dx + Dy*Dy + Dz*Dz)

              END IF

              IF((IJ.GT.4).AND.(IJ.LE.6)) THEN

                JK=IJ+2
                JL1=8*(NEL-1)+JK
                JL2=KONE(JL1)

                JCOUNT=JCOUNT +1
                NDSUP(0,I)=JCOUNT
                NDSUP(JCOUNT,I)=JL2
 
                Dx = xyz(1,JL2) - xyz(1,I)
                Dy = xyz(2,JL2) - xyz(2,I)
                Dz = xyz(3,JL2) - xyz(3,I)

                NDDIST(JCOUNT,I) = DSQRT(Dx*Dx + Dy*Dy + Dz*Dz)

              END IF

              IF(IJ.GT.6) THEN

                JK=IJ-2
                JL1=8*(NEL-1)+JK
                JL2=KONE(JL1)

                JCOUNT=JCOUNT +1
                NDSUP(0,I)=JCOUNT
                NDSUP(JCOUNT,I)=JL2

                Dx = xyz(1,JL2) - xyz(1,I)
                Dy = xyz(2,JL2) - xyz(2,I)
                Dz = xyz(3,JL2) - xyz(3,I)

                NDDIST(JCOUNT,I) = DSQRT(Dx*Dx + Dy*Dy +  Dz*Dz)

              END IF
                        
              EXIT

            END IF

          END DO

        END IF

      END DO


     ELSE

       DO J=1,II

         NEL=NEIBOR(J,I)
         VOLtot=Voltot+VOLU(NEL)

         JI=8*(NEL-1)
         N1=KONE(JI+1)
         N2=KONE(JI+2)
         N3=KONE(JI+3)
         N4=KONE(JI+4)
         N5=KONE(JI+5)
         N6=KONE(JI+6)
         N7=KONE(JI+7)
         N8=KONE(JI+8)

         IF(I.EQ.N1) THEN
           JND(1)=N2
           JND(2)=N4
           JND(3)=N5
         END IF

         IF(I.EQ.N2) THEN
           JND(1)=N1
           JND(2)=N3
           JND(3)=N6
         END IF

         IF(I.EQ.N3) THEN
           JND(1)=N2
           JND(2)=N4
           JND(3)=N7
         END IF

         IF(I.EQ.N4) THEN
           JND(1)=N3
           JND(2)=N1
           JND(3)=N8
         END IF

         IF(I.EQ.N5) THEN
           JND(1)=N6
           JND(2)=N8
           JND(3)=N1
         END IF

         IF(I.EQ.N6) THEN
           JND(1)=N5
           JND(2)=N7
           JND(3)=N2
         END IF

         IF(I.EQ.N7) THEN
           JND(1)=N6
           JND(2)=N8
           JND(3)=N3
         END IF

         IF(I.EQ.N8) THEN
           JND(1)=N7
           JND(2)=N5
           JND(3)=N4
         END IF

         IF(J.EQ.1) THEN

           DO IJ=1,3

             JJ=JND(IJ)
             JCOUNT=JCOUNT+1
             NDSUP(0,I)=JCOUNT
             NDSUP(JCOUNT,I)=JJ

             Dx = xyz(1,JJ) - xyz(1,I)
             Dy = xyz(2,JJ) - xyz(2,I)
             Dz = xyz(3,JJ) - xyz(3,I)

             NDDIST(JCOUNT,I) = DSQRT(Dx*Dx + Dy*Dy +  Dz*Dz)

           END DO

         ELSE

           DO IJ=1,3

             JPAR=0
             JJ=JND(IJ)

             DO IK=1,JCOUNT

               JK=NDSUP(IK,I)

               IF(JJ.EQ.JK) THEN
                 JPAR=1
                 EXIT
               end if

             end do

            IF(JPAR.EQ.0) THEN

              JCOUNT=JCOUNT+1
              NDSUP(0,I)=JCOUNT
              NDSUP(JCOUNT,I)=JJ

              Dx = xyz(1,JJ) - xyz(1,I)
              Dy = xyz(2,JJ) - xyz(2,I)
              Dz = xyz(3,JJ) - xyz(3,I)

              NDDIST(JCOUNT,I) = DSQRT(Dx*Dx + Dy*Dy +  Dz*Dz)

            end if

          end do

        end if

      end do

    end if

    FILTR2(I)=(VOLtot)**(1.0D0/3.0D0)

  end do


! ***********************

end if

! ***********************








                                                               
! CALCULO DO INTERVALO DE TEMPO QUE SERA ADOTADO 
! **********************************************         
      
Delt=DtMAX

! CALCULO DO TEMPO TOTAL

TF= TPOAC+NIR*(NTR-NFILE)*Delt

! ***************************************************************************
! O Tempo total e igual a NIR*NTR*Delt; se todavia nao tem-se gravado 
! nenhum registro, NFILE = TPOAC = 0 e assim obtem-se o mesmo. Porem, se 
! inicia-se a simulacao a partir de um certo registro previamente obtido, 
! TPOAC contem o tempo ("variavel fisica") alcancado na rodada anterior. 
! ***************************************************************************

TEMPO= TPOAC

! *****************************************************************************
! VNVX: Relacao entre o menor e o maior volume usado na formula de Hourglassing. 
! *****************************************************************************

VNVX=VOLMIN/VOLMAX	    
 

! MONTAJEM DO VETOR QUE CONTEM OS ELEMENTOS DA DIAGONAL PRINCIPAL DA MATRIZ
! *************************************************************************
! DE MASSA CONCENTRADA     
! ********************

mlum = 0.0d0
do i=1,nelem
  do j=i*8-7,i*8
    mlum(kone(j)) = mlum(kone(j))+volu(i)
  end do
end do
mlum=mlum*.125d0


! **********************

IF(MNOBJ.NE.0) THEN

! **********************

! *******************************************************************************
! Suavização das áreas de contorno sólido utilizadas no cálculo dos 
! coeficientes aerodinâmicos:
! *******************************************************************************

  DO NOB=1,NOBJ

    ArAux=0.0D0

    DO I=1,NFCS

      IF(I.GT.NFCN(NOB)) EXIT
      IELM=IELCS(NOB,I)
      ArAux(IELM)=(AREACS(NOB,I)*0.25D0)*VOLU(IELM)*0.125D0

    end do

    DO J=1,NNCS

      IF(J.GT.NNCN(NOB)) EXIT

      I=IBNOBJ(NOB,J)

      K10 = NEIBOR(1,I)
      K20 = NEIBOR(2,I)
      K30 = NEIBOR(3,I)
      K40 = NEIBOR(4,I)
      K50 = NEIBOR(5,I)
      K60 = NEIBOR(6,I)
      K70 = NEIBOR(7,I)
      K80 = NEIBOR(8,I)
      K90 = NEIBOR(9,I)
      K100 = NEIBOR(10,I)
      K110 = NEIBOR(11,I)
      K120 = NEIBOR(12,I)
      K130 = NEIBOR(13,I)
      K140 = NEIBOR(14,I)
      K150 = NEIBOR(15,I)
      K160 = NEIBOR(16,I)

      ARSuav(NOB,J)= ( ArAux(K10)+ArAux(K20)+ArAux(K30)+ArAux(K40)+ArAux(K50)+ArAux(K60) +&
                       ArAux(K70)+ArAux(K80)+ArAux(K90)+ArAux(K100)+ArAux(K110)+ArAux(K120) +&
                       ArAux(K130)+ArAux(K140)+ArAux(K150)+ArAux(K160) ) / MLUM(I)

      Prod = DFLOAT(NEIBOR(0,I))

      ARSuav(NOB,J)=ARSuav(NOB,J)*Prod
            

    end do

  end do

! **********************

end if

! **********************
        


! ABERTURA DOS ARQUIVOS DE RESIDUO E ARMAZENAGEM DOS COEFICIENTES AERODINÂMICOS
! *****************************************************************************

! **********************

IF(MNOBJ.NE.0) THEN

! **********************

  NOPF=29

  DO NOB=1,NOBJ

    NB1 = NOB/10 + 48
    NB2 = NOB - (NOB/10)*10 + 48

    NOPF=NOPF+1
    PTNAME(7:9)='.cd'
    PTNAME(10:11)=CHAR(NB1)//CHAR(NB2)
    PTNAME(12:14)='   '
    OPEN(NOPF,FILE=PTNAME,STATUS='UNKNOWN')
    REWIND(NOPF)

    NOPF=NOPF+1
    PTNAME(7:9)='.cl'
    PTNAME(10:11)=CHAR(NB1)//CHAR(NB2)
    PTNAME(12:14)='   '
    OPEN(NOPF,FILE=PTNAME,STATUS='UNKNOWN')
    REWIND(NOPF)

    NOPF=NOPF+1
    PTNAME(7:9)='.cz'
    PTNAME(10:11)=CHAR(NB1)//CHAR(NB2)
    PTNAME(12:14)='   '
    OPEN(NOPF,FILE=PTNAME,STATUS='UNKNOWN')
    REWIND(NOPF)

    NOPF=NOPF+1
    PTNAME(7:10)='.cmx'
    PTNAME(11:12)=CHAR(NB1)//CHAR(NB2)
    PTNAME(13:14)='  '
    OPEN(NOPF,FILE=PTNAME,STATUS='UNKNOWN')
    REWIND(NOPF)

    NOPF=NOPF+1
    PTNAME(7:10)='.cmy'
    PTNAME(11:12)=CHAR(NB1)//CHAR(NB2)
    PTNAME(13:14)='  '
    OPEN(NOPF,FILE=PTNAME,STATUS='UNKNOWN')
    REWIND(NOPF)

    NOPF=NOPF+1
    PTNAME(7:10)='.cmz'
    PTNAME(11:12)=CHAR(NB1)//CHAR(NB2)
    PTNAME(13:14)='  '
    OPEN(NOPF,FILE=PTNAME,STATUS='UNKNOWN')
    REWIND(NOPF)

    NOPF=NOPF+1
    PTNAME(7:9)='.cp'
    PTNAME(10:11)=CHAR(NB1)//CHAR(NB2)
    PTNAME(12:14)='   '
    OPEN(NOPF,FILE=PTNAME,STATUS='UNKNOWN')
    REWIND(NOPF)

    NOPF=NOPF+1
    PTNAME(7:9)='.cf'
    PTNAME(10:11)=CHAR(NB1)//CHAR(NB2)
    PTNAME(12:14)='   '
    OPEN(NOPF,FILE=PTNAME,STATUS='UNKNOWN')
    REWIND(NOPF)

  end do

! **********************

end if

! **********************
 						   

OPEN(14,FILE='residuo.dat',STATUS='UNKNOWN')
REWIND(14)
      
INDTPO= 0
NNR= 0
NNNRRR= 0
NCOUNTCOEF=0
VDTURB=0.0D0


! ****************************************************
! Inicio do ciclo de tempo
! ****************************************************

DO WHILE ((TEMPO.LT.TF).AND.(INDTPO.LE.NROTPO))

! ****************************************************
! Inicio do ciclo de tempo
! ****************************************************

  TEMPO= TEMPO+Delt
  NNR= NNR+1
  NNNRRR= NNNRRR+1
  NCOUNTCOEF=NCOUNTCOEF+1
  write(*,910) NNNRRR,tempo

! Retencao da pressao inicial para o calculo do residuo:

  Pri = uvw(4,:)


! **********************
! ***** PRIMEIRO PASSO:
! **********************


! Derivadas de velocidade Pij e componentes de velocidade media no elemento
! *************************************************************************
! UPROM,VPROM e WPROM e gradientes de pressao GriP:
! *************************************************

  DO IELEM=1,NELEM

    JINV11= JIN(1,IELEM)
    JINV12= JIN(2,IELEM)
    JINV13= JIN(3,IELEM)
    JINV21= JIN(4,IELEM)
    JINV22= JIN(5,IELEM)
    JINV23= JIN(6,IELEM)
    JINV31= JIN(7,IELEM)
    JINV32= JIN(8,IELEM)
    JINV33= JIN(9,IELEM)

    II= (IELEM-1)*8
       
	do j=1,8
	  uvwp(j,1) = uvw(1,kone(ii+j))
	  uvwp(j,2) = uvw(2,kone(ii+j))
	  uvwp(j,3) = uvw(3,kone(ii+j))
	  uvwp(j,4) = uvw(4,kone(ii+j))
    end do





!   ***************************************************************************
!   Somatorio dos produtos dos elementos da matriz jacobiana pelas
!   coordenadas naturais; para serem usados em produtos formados abaixo. 
!   ***************************************************************************
      
    JCN(1,1) = -JINV11-JINV12-JINV13
    JCN(2,1) =  JINV11-JINV12-JINV13 
    JCN(3,1) =  JINV11+JINV12-JINV13 
    JCN(4,1) = -JINV11+JINV12-JINV13
    JCN(5,1) = -JINV11-JINV12+JINV13
    JCN(6,1) =  JINV11-JINV12+JINV13 
    JCN(7,1) =  JINV11+JINV12+JINV13 
    JCN(8,1) = -JINV11+JINV12+JINV13
    
	JCN(1,2) = -JINV21-JINV22-JINV23
    JCN(2,2) =  JINV21-JINV22-JINV23 
    JCN(3,2) =  JINV21+JINV22-JINV23 
    JCN(4,2) = -JINV21+JINV22-JINV23
    JCN(5,2) = -JINV21-JINV22+JINV23
    JCN(6,2) =  JINV21-JINV22+JINV23 
    JCN(7,2) =  JINV21+JINV22+JINV23 
    JCN(8,2) = -JINV21+JINV22+JINV23
     
    JCN(1,3) = -JINV31-JINV32-JINV33 
    JCN(2,3) =  JINV31-JINV32-JINV33
    JCN(3,3) =  JINV31+JINV32-JINV33       
    JCN(4,3) = -JINV31+JINV32-JINV33
    JCN(5,3) = -JINV31-JINV32+JINV33
    JCN(6,3) =  JINV31-JINV32+JINV33
    JCN(7,3) =  JINV31+JINV32+JINV33
    JCN(8,3) = -JINV31+JINV32+JINV33

!   **************************************************************************
!   Media dos valores nodais dos produtos das coord naturais pela
!   inversa de J(0) e pelas componentes de velocidade; CC Newman
!   **************************************************************************

    pp(1,IELEM)=dot_product(jcn(:,1),uvwp(:,1))*0.125D0
    pp(4,IELEM)=dot_product(jcn(:,2),uvwp(:,1))*0.125D0
    pp(7,IELEM)=dot_product(jcn(:,3),uvwp(:,1))*0.125D0
    pp(2,IELEM)=dot_product(jcn(:,1),uvwp(:,2))*0.125D0
    pp(5,IELEM)=dot_product(jcn(:,2),uvwp(:,2))*0.125D0
    pp(8,IELEM)=dot_product(jcn(:,3),uvwp(:,2))*0.125D0
    pp(3,IELEM)=dot_product(jcn(:,1),uvwp(:,3))*0.125D0
    pp(6,IELEM)=dot_product(jcn(:,2),uvwp(:,3))*0.125D0
    pp(9,IELEM)=dot_product(jcn(:,3),uvwp(:,3))*0.125D0




!   ***************************************************************************
!   Media dos valores nodais dos produtos das coordenadas naturais,
!   a inversa de J(0) e pelas tres componentes de velocidade; CC Newman
!   ***************************************************************************

!   divergência da velocidade (dVk/dXk):
       
	PNV(IELEM) = pp(1,IELEM) + pp(5,IELEM) + pp(9,IELEM)

!   gradientes de pressao (dp/dXk):

    Gr1P(IELEM) = dot_product(jcn(:,1),uvwp(:,4)) * 0.125D0
    Gr2P(IELEM) = dot_product(jcn(:,2),uvwp(:,4)) * 0.125D0
    Gr3P(IELEM) = dot_product(jcn(:,3),uvwp(:,4)) * 0.125D0

!   ***************************************************************************
!   Media dos valores nodais das componentes de velocidade.
!   ***************************************************************************

    uvwprom(1,ielem)= sum(uvwp(:,1))*0.125D0
    uvwprom(2,ielem)= sum(uvwp(:,2))*0.125D0
    uvwprom(3,ielem)= sum(uvwp(:,3))*0.125D0



  END DO




! Insercao do modelo de turbulencia e obtencao das matrizes difusivas Dij:
! ************************************************************************

  IF(IndTurb.NE.0) THEN

! *******************************************************************************
! --> Simulação Direta de Grandes Escalas (LES) com Modelo de Smagorinsky para 
!     as escalas inferiores a resolução da malha.
!
!     CASO 1 (IndTurb = 1):
!     ---------------------
! --> CÁLCULO DA VISCOSIDADE TURBULENTA SEGUNDO O MODELO CLÁSSICO DE SMAGORINSKY
!     E OBTENÇÃO DAS MATRIZES DE DIFUSIVIDADE CONSIDERANDO A PRESENÇA DA
!     VISCOSIDADE TURBULENTA ADICIONADA A VISCOSIDADE DINÂMICA.
!
!     CASO 2 (IndTurb = 2):
!     ---------------------
! --> CÁLCULO DA VISCOSIDADE TURBULENTA SEGUNDO O MODELO DINÂMICO DE SMAGORINSKY
!     E OBTENÇÃO DAS MATRIZES DE DIFUSIVIDADE CONSIDERANDO A PRESENÇA DA
!     VISCOSIDADE TURBULENTA ADICIONADA A VISCOSIDADE DINÂMICA.
! *******************************************************************************

    SELECT CASE(IndTurb)


!     ****        <<<<<<<<<<<<< MODELO CLÁSSICO DE SMAGORINSKY >>>>>>>>>>>>
      CASE(1)

        DO NELM=1,NELEM

          dux=pp(1,NELM)
          duy=pp(4,NELM)
          duz=pp(7,NELM)
          dvx=pp(2,NELM)
          dvy=pp(5,NELM)
          dvz=pp(8,NELM)
          dwx=pp(3,NELM)
          dwy=pp(6,NELM)
          dwz=pp(9,NELM)

!         DSijSij = 2 * Sij * Sij:


          DSijSij = 2.0D0*(dux*dux+dvy*dvy+dwz*dwz)+2.0D0*(duy*dvx+duz*dwx+dvz*dwy)+ &
                           duy*duy+dvx*dvx+duz*duz+dwx*dwx+dvz*dvz+dwy*dwy

          Filtro = (VOLU(NELM))**(1.0D0/3.0D0)

          VDTurb(NELM) = Cs*Cs*Filtro*Filtro * DSQRT(DSijSij)

        end do



!     *****        <<<<<<<<<<<<< MODELO DINÂMICO DE SMAGORINSKY >>>>>>>>>>>>
      CASE(2)


!       ***** --> Determinação das componentes de velocidade nodais de segundo filtro;
!       ***** --> Determinação das componentes do tensor de Leonard Lij nodais;
!       ***** --> Suaviazação do centro do elemento para os nós das derivadas das componentes
!       *****     de velocidade de primeiro e segundo filtro;

!       ***** Obtenção das componentes de velocidade de segundo filtro e do tensor de
!       ***** Leonard Lij, em nível nodal:

        DO I=1,NNM

          DistB = 0.0d0
	      uu    = 0.0d0
	      vv    = 0.0d0
	      ww    = 0.0d0

          do j = 1,12
            jj=ndsup(j,i)
            if(jj.ne.0) then
              uu(j)=uvw(1,jj)
              vv(j)=uvw(2,jj)
	          ww(j)=uvw(3,jj)
	          DistB(j)=1.0D0/NDDIST(j,I)
	        end if
	      end do

          Dtot = sum(DistB)

          Ufa = dot_product(uu,DistB)/Dtot
          Vfa = dot_product(vv,DistB)/Dtot
          Wfa = dot_product(ww,DistB)/Dtot
 
!         UUfa = dot_product(uu*uu,DistB)/Dtot
!         VVfa = dot_product(vv*vv,DistB)/Dtot
!         WWfa = dot_product(ww*ww,DistB)/Dtot
!         UVfa = dot_product(uu*vv,DistB)/Dtot
!         UWfa = dot_product(uu*ww,DistB)/Dtot
!         VWfa = dot_product(vv*ww,DistB)/Dtot


!         ***** Fator de ponderação segundo Silveira Neto (p):

          pond = 0.0D0

!         uvwf(1,I) = Ufa * (1.0D0 - pond) + pond * uvw(1,I)
!         uvwf(2,I) = Vfa * (1.0D0 - pond) + pond * uvw(2,I)
!         uvwf(3,I) = Wfa * (1.0D0 - pond) + pond * uvw(3,I)
!         UUf = UUfa * (1.0D0 - pond) + pond * (uvw(1,I) * uvw(1,I))
!         VVf = VVfa * (1.0D0 - pond) + pond * (uvw(2,I) * uvw(2,I))
!         WWf = WWfa * (1.0D0 - pond) + pond * (uvw(3,I) * uvw(3,I))
!         UVf = UVfa * (1.0D0 - pond) + pond * (uvw(1,I) * uvw(2,I))
!         UWf = UWfa * (1.0D0 - pond) + pond * (uvw(1,I) * uvw(3,I))
!         VWf = VWfa * (1.0D0 - pond) + pond * (uvw(2,I) * uvw(3,I))

          uvwf(1,I) = Ufa 
          uvwf(2,I) = Vfa
          uvwf(3,I) = Wfa


!         UUf = UUfa 
!         VVf = VVfa
!         WWf = WWfa
!         UVf = UVfa 
!         UWf = UWfa 
!         VWf = VWfa 

!         ***** Determinação das componentes do tensor de Leonard Lij em nível nodal: 

!         SL11(I) = UUf - (uvwf(1,I)*Uvwf(1,I))
!         SL22(I) = VVf - (uVwf(2,I)*uVwf(2,I))
!         SL33(I) = WWf - (uvWf(3,I)*uvWf(3,I))
!         SL12(I) = UVf - (Uvwf(1,I)*uVwf(2,I))
!         SL13(I) = UWf - (Uvwf(1,I)*uvWf(3,I))
!         SL23(I) = VWf - (uVwf(2,I)*uvWf(3,I))


!         *****  Aplicação das condições de contorno sobre os tensores de segundo filtro SLij:
!         --> Todas as componentes do tensor serão anuladas para nós pertencentes a 
!         *****    superfícies sólidas (paredes e/ou corpos imersos, ambos fixos).

          IF((uvw(1,I).EQ.0.0D0).AND.(uvw(2,I).EQ.0.0D0).AND.(uvw(3,I).EQ.0.0D0)) THEN

!           SL11(I) = 0.0D0
!           SL22(I) = 0.0D0
!           SL33(I) = 0.0D0
!           SL12(I) = 0.0D0
!           SL13(I) = 0.0D0
!           SL23(I) = 0.0D0
                     
	        Uvwf(1,I) = 0.0D0
            uVwf(2,I) = 0.0D0
            uvWf(3,I) = 0.0D0

          end if

        end do

!       ***** --> Suaviazação do centro do elemento para os nós das 
!       *****     derivadas das componentes de velocidade de primeiro filtro e da dimensão
!       *****     característica do primeiro filtro;
!       ***** --> Cálculo do módulo do tensor taxa de deformação para o primeiro filtro em
!       *****     nível de elemento

        filts1 = 0.0d0
	    dudx = 0.0d0
	    dudy = 0.0d0
        dudz = 0.0d0
	    dvdx = 0.0d0
        dvdy = 0.0d0
	    dvdz = 0.0d0
 	    dwdx = 0.0d0
	    dwdy = 0.0d0
	    dwdz = 0.0d0

        DO NELM = 1,NELEM

	      VOLU1 = VOLU(NELM)*.125D0

          II= (NELM-1)*8
       
          dux=pp(1,NELM)
          duy=pp(4,NELM)
          duz=pp(7,NELM)
          dvx=pp(2,NELM)
          dvy=pp(5,NELM)
          dvz=pp(8,NELM)
          dwx=pp(3,NELM)
          dwy=pp(6,NELM)
          dwz=pp(9,NELM)

          SmodE(NELM)=DSQRT(2.0D0*(dux*dux+dvy*dvy+dwz*dwz)+2.0D0*(duy*dvx+duz*dwx+dvz*dwy)+ &
                                   duy*duy+dvx*dvx+duz*duz+dwx*dwx+dvz*dvz+dwy*dwy)


          FILTR1(NELM) = VOLU(NELM)**(1.0D0/3.0D0)

          dux = dux*volu1
          duy = duy*volu1
          duz = duz*volu1
          dvx = dvx*volu1
          dvy = dvy*volu1
          dvz = dvz*volu1
          dwx = dwx*volu1
          dwy = dwy*volu1
          dwz = dwz*volu1

          filtr3 =FILTR1(NELM) * volu1


          do j=1,8

            Filts1(kone(ii+j))=Filts1(kone(ii+j)) + filtr3
			 
	        dudx(kone(ii+j))=dudx(kone(ii+j))+dux
	        dudy(kone(ii+j))=dudy(kone(ii+j))+duy
	        dudz(kone(ii+j))=dudz(kone(ii+j))+duz

            dvdx(kone(ii+j))=dvdx(kone(ii+j))+dvx
	        dvdy(kone(ii+j))=dvdy(kone(ii+j))+dvy
	        dvdz(kone(ii+j))=dvdz(kone(ii+j))+dvz
      
	        dwdx(kone(ii+j))=dwdx(kone(ii+j))+dwx
	        dwdy(kone(ii+j))=dwdy(kone(ii+j))+dwy
	        dwdz(kone(ii+j))=dwdz(kone(ii+j))+dwz

          end do

        end do

        Filts1 = Filts1/mlum
        dudx = dudx/mlum
  	    dudy = dudy/mlum
  	    dudz = dudz/mlum
 	    dvdx = dvdx/mlum
	    dvdy = dvdy/mlum
	    dvdz = dvdz/mlum
 	    dwdx = dwdx/mlum
	    dwdy = dwdy/mlum
	    dwdz = dwdz/mlum



!       ***** Fim da Suavização para o primeiro filtro!



!       ***** --> Determinação das componentes do tensor taxa de deformação Sij e seu módulo
!       *****     para o primeiro filtro em nível nodal.

        DO I=1,NNM

          S(1,I) = dudx(I)
          S(2,I) = 0.5D0 * (dudy(I) + dvdx(I))
          S(3,I) = 0.5D0 * (dudz(I) + dwdx(I))
          S(4,I) = dvdy(I)
          S(5,I) = 0.5D0 * (dvdz(I) + dwdy(I))
          S(6,I) = dwdz(I)

          Smod(I)=DSQRT(2.0D0*(dudx(I)*dudx(i)+dvdy(I)*dvdy(i)+dwdz(I)*dwdz(i))+ &
                        2.0D0*(dudy(I)*dvdx(I)+dudz(I)*dwdx(I)+dvdz(I)*dwdy(I))+ &
                               dudy(I)*dudy(i)+dvdx(I)*dvdx(i)+dudz(I)*dudz(i)+ &
                               dwdx(I)*dwdx(i)+dvdz(I)*dvdz(i)+dwdy(I)*dwdy(i))

        end do


!       ***** --> Suaviazação do centro do elemento para os nós das 
!       *****     derivadas das componentes de velocidade de segundo filtro;

	    dudx = 0.0d0
        dudy = 0.0d0
	    dudz = 0.0d0
        dvdx = 0.0d0
        dvdy = 0.0d0
        dvdz = 0.0d0
        dwdx = 0.0d0
	    dwdy = 0.0d0
	    dwdz = 0.0d0


        DO NELM= 1,NELEM

          volu1 = volu(nelm)*.015625d0

          JINV11= JIN(1,NELM)
          JINV12= JIN(2,NELM)
          JINV13= JIN(3,NELM)
          JINV21= JIN(4,NELM)
          JINV22= JIN(5,NELM)
          JINV23= JIN(6,NELM)
          JINV31= JIN(7,NELM)
          JINV32= JIN(8,NELM)
          JINV33= JIN(9,NELM)

          II= (NELM-1)*8

	      do j=1,8
	        uvwp(j,1) = uvwf(1,kone(ii+j))
	        uvwp(j,2) = uvwf(2,kone(ii+j))
	        uvwp(j,3) = uvwf(3,kone(ii+j))
          end do
       

!         ***************************************************************************
!         Somatorio dos produtos dos elementos da matriz jacobiana pelas
!         coordenadas naturais; para serem usados em produtos formados abaixo. 
!         ***************************************************************************
      
	      JCN(1,1) = -JINV11-JINV12-JINV13
          JCN(2,1) =  JINV11-JINV12-JINV13 
          JCN(3,1) =  JINV11+JINV12-JINV13 
          JCN(4,1) = -JINV11+JINV12-JINV13
          JCN(5,1) = -JINV11-JINV12+JINV13
          JCN(6,1) =  JINV11-JINV12+JINV13 
          JCN(7,1) =  JINV11+JINV12+JINV13 
          JCN(8,1) = -JINV11+JINV12+JINV13
     
          JCN(1,2) = -JINV21-JINV22-JINV23
          JCN(2,2) =  JINV21-JINV22-JINV23 
          JCN(3,2) =  JINV21+JINV22-JINV23 
          JCN(4,2) = -JINV21+JINV22-JINV23
          JCN(5,2) = -JINV21-JINV22+JINV23
          JCN(6,2) =  JINV21-JINV22+JINV23 
          JCN(7,2) =  JINV21+JINV22+JINV23 
          JCN(8,2) = -JINV21+JINV22+JINV23
     
          JCN(1,3) = -JINV31-JINV32-JINV33 
          JCN(2,3) =  JINV31-JINV32-JINV33
          JCN(3,3) =  JINV31+JINV32-JINV33       
          JCN(4,3) = -JINV31+JINV32-JINV33
          JCN(5,3) = -JINV31-JINV32+JINV33
          JCN(6,3) =  JINV31-JINV32+JINV33
          JCN(7,3) =  JINV31+JINV32+JINV33
          JCN(8,3) = -JINV31+JINV32+JINV33

!         **************************************************************************
!  	      Media dos valores nodais dos produtos das coord naturais pela
!         inversa de J(0) e pelas componentes de velocidade; CC Newman
!         **************************************************************************
      

          dux = dot_product(jcn(:,1),uvwp(:,1))*VOLU1
          duy = dot_product(jcn(:,2),uvwp(:,1))*VOLU1
          duz = dot_product(jcn(:,3),uvwp(:,1))*VOLU1

          dvx = dot_product(jcn(:,1),uvwp(:,2))*VOLU1
          dvy = dot_product(jcn(:,2),uvwp(:,2))*VOLU1
          dvz = dot_product(jcn(:,3),uvwp(:,2))*VOLU1

          dwx = dot_product(jcn(:,1),uvwp(:,3))*VOLU1
          dwy = dot_product(jcn(:,2),uvwp(:,3))*VOLU1
          dwz = dot_product(jcn(:,3),uvwp(:,3))*VOLU1

          do j=1,8

            dudx(kone(ii+j))=dudx(kone(ii+j))+dux
	        dudy(kone(ii+j))=dudy(kone(ii+j))+duy
	        dudz(kone(ii+j))=dudz(kone(ii+j))+duz

            dvdx(kone(ii+j))=dvdx(kone(ii+j))+dvx
	        dvdy(kone(ii+j))=dvdy(kone(ii+j))+dvy
	        dvdz(kone(ii+j))=dvdz(kone(ii+j))+dvz
      
	        dwdx(kone(ii+j))=dwdx(kone(ii+j))+dwx
	        dwdy(kone(ii+j))=dwdy(kone(ii+j))+dwy
	        dwdz(kone(ii+j))=dwdz(kone(ii+j))+dwz

          end do

        END DO

        dudx = dudx/mlum
        dudy = dudy/mlum
	    dudz = dudz/mlum
 	    dvdx = dvdx/mlum
	    dvdy = dvdy/mlum
	    dvdz = dvdz/mlum
 	    dwdx = dwdx/mlum
	    dwdy = dwdy/mlum
	    dwdz = dwdz/mlum

!       ***** Fim da Suavização para o segundo filtro!
 
!       ***** --> Determinação do produto <(Filtro1**2)*|Sij|*Sij> para o segundo filtro.
!       ***** --> Determinação do coeficiente dinâmico em nível nodal: C(x,t)

        DO I=1,NNM

          DistB = 0.0d0
          fs    = 0.0d0
	      s11   = 0.0d0
	      s22   = 0.0d0
          s33   = 0.0d0
	      s12   = 0.0d0
	      s13   = 0.0d0
	      s23   = 0.0d0

          do j = 1,12
            jj=ndsup(j,i)
            if(jj.ne.0) then
              s11(j)=s(1,jj)
              s12(j)=s(2,jj)
              s13(j)=s(3,jj)
              s22(j)=s(4,jj)
              s23(j)=s(5,jj)
              s33(j)=s(6,jj)
              DistB(j)=1.0D0/NDDIST(j,I)
              fs(j)=FILTS1(jj)*filts1(jj)*Smod(jj)*DistB(j)
	        end if
	      end do

          Dtot = sum(DistB)

          SS11fa = dot_product(fs,s11)/Dtot
          SS22fa = dot_product(fs,s22)/Dtot
          SS33fa = dot_product(fs,s33)/Dtot
          SS12fa = dot_product(fs,s12)/Dtot
          SS13fa = dot_product(fs,s13)/Dtot
          SS23fa = dot_product(fs,s23)/Dtot

!         ***** Fator de ponderação segundo Silveira Neto (p):

          pond = 0.0D0

!         fs11 = filts1(i) * filts1(i) * pond * Smod(i)

!         SS11f = SS11fa * (1.0D0 - pond) + fs11 * S11(I)
!         SS22f = SS22fa * (1.0D0 - pond) + fs11 * S22(I) 
!         SS33f = SS33fa * (1.0D0 - pond) + fs11 * S33(I) 
!         SS12f = SS12fa * (1.0D0 - pond) + fs11 * S12(I) 
!         SS13f = SS13fa * (1.0D0 - pond) + fs11 * S13(I) 
!         SS23f = SS23fa * (1.0D0 - pond) + fs11 * S23(I) 

          SS11f = SS11fa
          SS22f = SS22fa
          SS33f = SS33fa 
          SS12f = SS12fa
          SS13f = SS13fa
          SS23f = SS23fa

!         ***** --> Determinação das componentes do tensor taxa de deformação Sij e seu módulo
!         *****     para o segundo filtro:

          S11f = dudx(I)
          S12f = 0.5D0 * (dudy(I) + dvdx(I))
          S13f = 0.5D0 * (dudz(I) + dwdx(I))
          S22f = dvdy(I)
          S23f = 0.5D0 * (dvdz(I) + dwdy(I))
          S33f = dwdz(I)

          Sfmod=DSQRT(2.0D0*(dudx(I)*dudx(i)+dvdy(I)*dvdy(i)+dwdz(I)*dwdz(i))+ &
                      2.0D0*(dudy(I)*dvdx(I)+dudz(I)*dwdx(I)+dvdz(I)*dwdy(I))+ &
                             dudy(I)*dudy(i)+dvdx(I)*dvdx(i)+dudz(I)*dudz(i)+ &
                             dwdx(I)*dwdx(i)+dvdz(I)*dvdz(i)+dwdy(I)*dwdy(i))


!         ***** Cálculo do Tensor Mij:

          filtr4 = filtr2(i)*filtr2(i)*Sfmod  
          SM11 = filtr4 * S11f - SS11f
          SM22 = filtr4 * S22f - SS22f
          SM33 = filtr4 * S33f - SS33f
          SM12 = filtr4 * S12f - SS12f
          SM13 = filtr4 * S13f - SS13f
          SM23 = filtr4 * S23f - SS23f


!         *****        <<<<<<<<<<<<< MODELO DINÂMICO DE SMAGORINSKY >>>>>>>>>>>>

!         ***** Obtenção das componentes de velocidade de segundo filtro e do tensor de
!         ***** Leonard Lij, em nível nodal:

          DistB = 0.0d0
          uu    = 0.0d0
	      vv    = 0.0d0
	      ww    = 0.0d0

          do j = 1,12
            jj=ndsup(j,i)
            if(jj.ne.0) then
              uu(j)=uvw(1,jj)
              vv(j)=uvw(2,jj)
	          ww(j)=uvw(3,jj)
	          DistB(j)=1.0D0/NDDIST(j,I)
            end if
	      end do

          Dtot = sum(DistB)

          UUfa = dot_product(uu*uu,DistB)/Dtot
          VVfa = dot_product(vv*vv,DistB)/Dtot
          WWfa = dot_product(ww*ww,DistB)/Dtot
          UVfa = dot_product(uu*vv,DistB)/Dtot
          UWfa = dot_product(uu*ww,DistB)/Dtot
          VWfa = dot_product(vv*ww,DistB)/Dtot

!         ***** Fator de ponderação segundo Silveira Neto (p):

          pond = 0.0D0

!         UUf = UUfa * (1.0D0 - pond) + pond * (uvw(1,I) * uvw(1,I))
!         VVf = VVfa * (1.0D0 - pond) + pond * (uvw(2,I) * uvw(2,I))
!         WWf = WWfa * (1.0D0 - pond) + pond * (uvw(3,I) * uvw(3,I))
!         UVf = UVfa * (1.0D0 - pond) + pond * (uvw(1,I) * uvw(2,I))
!         UWf = UWfa * (1.0D0 - pond) + pond * (uvw(1,I) * uvw(3,I))
!         VWf = VWfa * (1.0D0 - pond) + pond * (uvw(2,I) * uvw(3,I))

          UUf = UUfa 
          VVf = VVfa
          WWf = WWfa
          UVf = UVfa 
          UWf = UWfa 
          VWf = VWfa 

!         ***** Determinação das componentes do tensor de Leonard Lij em nível nodal: 

          SL11 = UUf - (uvwf(1,I)*Uvwf(1,I))
          SL22 = VVf - (uVwf(2,I)*uVwf(2,I))
          SL33 = WWf - (uvWf(3,I)*uvWf(3,I))
          SL12 = UVf - (Uvwf(1,I)*uVwf(2,I))
          SL13 = UWf - (Uvwf(1,I)*uvWf(3,I))
          SL23 = VWf - (uVwf(2,I)*uvWf(3,I))


!         *****  Aplicação das condições de contorno sobre os tensores de segundo filtro SLij:
!         --> Todas as componentes do tensor serão anuladas para nós pertencentes a 
!         *****    superfícies sólidas (paredes e/ou corpos imersos, ambos fixos).

          IF(((uvw(1,I).EQ.0.0D0).AND.(uvw(2,I).EQ.0.0D0).AND. (uvw(3,I).EQ.0.0D0))) THEN
 
            SL11 = 0.0D0
            SL22 = 0.0D0
            SL33 = 0.0D0
            SL12 = 0.0D0
            SL13 = 0.0D0
            SL23 = 0.0D0

          end if

!         *****  Cálculo do coeficiente de Smagorinsky Dinâmico - C(x,t):
            
!         Snum = SL11(I)*SM11 + SL22(I)*SM22 + SL33(I)*SM33 + 2.0D0*SL12(I)*SM12 + 2.0D0*SL13(I)*SM13 + &
!                2.0D0*SL23(I)*SM23


          Snum = SL11*SM11+SL22*SM22+SL33*SM33+2.0D0*SL12*SM12+2.0D0*SL13*SM13+2.0D0*SL23*SM23

          Sden = SM11*sm11+SM22*sm22+SM33*sm33+2.0D0*(SM12*sm12+SM13*sm13+SM23*sm23) 

          IF(Sden.EQ.0.0D0) THEN
            Cxt(I) = 0.0D0
          ELSE
            Cxt(I)	= -0.5D0 * (Snum/Sden)
          end if

        end do


!       *****  Cálculo da viscosidade turbulenta segundo o modelo sub-malha
!       ***** dinâmico de Smagorinsky para o elemento IELEM:

        do IELEM=1,NELEM 

          IE = 8*(IELEM-1)
 
          NE1=KONE(IE+1)
          NE2=KONE(IE+2)
          NE3=KONE(IE+3)
          NE4=KONE(IE+4)
          NE5=KONE(IE+5)
          NE6=KONE(IE+6)
          NE7=KONE(IE+7)
          NE8=KONE(IE+8)

          VDTurb(IELEM) = 0.125D0*(Cxt(NE1)+Cxt(NE2)+Cxt(NE3)+Cxt(NE4)+Cxt(NE5)+Cxt(NE6) +&
                                   Cxt(NE7)+Cxt(NE8))*(FILTR1(IELEM)*FILTR1(IELEM))*SmodE(IELEM)

        end do


      CASE DEFAULT

    END SELECT



!   Matrizes difusivas Dij:
!   ***********************

!    IF(IndTurb.EQ.2) THEN

!     *****  Cálculo da viscosidade turbulenta segundo o modelo sub-malha
!     ***** dinâmico de Smagorinsky para o elemento IELEM:

!      DO IELEM=1,NELEM 

!        IE = 8*(IELEM-1)
 
!        NE1=KONE(IE+1)
!        NE2=KONE(IE+2)
!        NE3=KONE(IE+3)
!        NE4=KONE(IE+4)
!        NE5=KONE(IE+5)
!        NE6=KONE(IE+6)
!        NE7=KONE(IE+7)
!        NE8=KONE(IE+8)

!        VDTurb(IELEM) = 0.125D0*(Cxt(NE1)+Cxt(NE2)+Cxt(NE3)+Cxt(NE4)+Cxt(NE5)+Cxt(NE6) +&
!                                 Cxt(NE7)+Cxt(NE8))*(FILTR1(IELEM)*FILTR1(IELEM))*SmodE(IELEM)

!      END DO

!    end if  

  end if




! Formacao dos vetores locais RRn e RRP representando o lado direito das
! **********************************************************************
! equacoes de conservacao:
! ************************


  rrg = 0.0d0

  DO IELEM=1,NELEM 
  
    II= (IELEM-1)*8
    do j=1,8
      uvwp(j,1) = uvw(1,kone(ii+j))
      uvwp(j,2) = uvw(2,kone(ii+j))
      uvwp(j,3) = uvw(3,kone(ii+j))
      uvwp(j,4) = uvw(4,kone(ii+j))
    end do

!   ***************************
!   Matriz de massa consistente
!   ***************************

    VC= VOLU(IELEM)*0.015625D0

    mass=vc*transpose(mmm)
	masd=volu(ielem)/8.0d0
	masl=(1.0d0-elump)*mass
	do i=1,8
	  masl(i,i)=masl(i,i)+elump*masd
	end do


!   **** Controle de modos espurios (HOURGLASSING) - CHG:

!   ***************************************************************************
!                        Chg dado pelo criterio CHRISTON 
!   ***************************************************************************

    CHG= ( VOLU(IELEM)**(1.0D0/3.D0) ) * 1.0D0 * 0.5D0 * VOLU(IELEM)&
           * ( 1.0D0 - ( VNVX**(1.D0/3.D0) ) )*CONTROL1


    Aux = RHOInf*(VelSom**2)
    AuxADV = 0.0D0

    if(IndTurb.eq.0) then
      aux1 = ViscCin
    else 
      AUX1 = ViscCin + VDTurb(IELEM)
      IF(AUX1.LT.0.0D0) AUX1 = ViscCin
    end if
         
    COEF1= 2.0D0 * AUX1 + ViscVol
    COEF2= AUX1 
    COEF3= ViscVol         


    jbb(1,1) = JIN(1,IELEM) * DETERJ(IELEM)
    jbb(1,2) = JIN(2,IELEM) * DETERJ(IELEM)
    jbb(1,3) = JIN(3,IELEM) * DETERJ(IELEM)
    jbb(2,1) = JIN(4,IELEM) * DETERJ(IELEM)
    jbb(2,2) = JIN(5,IELEM) * DETERJ(IELEM)
    jbb(2,3) = JIN(6,IELEM) * DETERJ(IELEM)
    jbb(3,1) = JIN(7,IELEM) * DETERJ(IELEM)
    jbb(3,2) = JIN(8,IELEM) * DETERJ(IELEM)
    jbb(3,3) = JIN(9,IELEM) * DETERJ(IELEM)

    jbt = transpose(jbb)


    AUX1 = Delt*0.25D0*uvwprom(1,ielem)*uvwprom(1,ielem)
    AUX2 = Delt*0.25D0*uvwprom(2,ielem)*uvwprom(2,ielem)
    AUX3 = Delt*0.25D0*uvwprom(3,ielem)*uvwprom(3,ielem)
    AUX12 = Delt*0.25D0*uvwprom(1,ielem)*uvwprom(2,ielem)
    AUX13 = Delt*0.25D0*uvwprom(1,ielem)*uvwprom(3,ielem)
    AUX23 = Delt*0.25D0*uvwprom(2,ielem)*uvwprom(3,ielem)


!   bt = 1.0D6
!   AUX4 = VDTurb(IELEM)/bt
     
	AUX4 = 0.0D0


!   **************************************************************************************
!   Matrizes de termos do tensor de balanco difusivo e de turbulencia na equacao de massa:
!
!   Eq. de quantidade de mov.:
!     BTD - adicionada aos termos volumetricos das matrizes difusivas (Dii):
!
!   Eq. da pressao:
!     BTE - termos volumetricos do tensor difusivo
!     BTN - termos cruzados do tensor de balanco difusivo
!     DFP - termos volumetricos de turbulencia - laplaciano da pressao
!   **************************************************************************************


!   Matriz do termo de pressao integrado por partes AGi:

    do j=1,8
	  do jj=1,8
	    AG1(jj,J)=(jbb(1,1)*TB1(jj,J)+jbb(1,2)*TB2(jj,J)+ jbb(1,3)*TB3(jj,J))*.125D0
        AG2(jj,J)=(jbb(2,1)*TB1(jj,J)+jbb(2,2)*TB2(jj,J)+jbb(2,3)*TB3(jj,J))*.125D0
        AG3(jj,J)=(jbb(3,1)*TB1(jj,J)+jbb(3,2)*TB2(jj,J)+jbb(3,3)*TB3(jj,J))*.125D0
      end do
	end do

!   Matrizes auxiliares Bi utilizadas no termo de adveccao e de contorno de difusao:

    do J=1,8
      do jj=1,8
        jjj=jj*3
        b123(jjj-2,j)=(jbb(1,1)*BB1(jj,J)+jbb(1,2)*BB2(jj,J)+ jbb(1,3)*BB3(jj,J))*.125D0
        b123(jjj-1,j)=(jbb(2,1)*BB1(jj,J)+jbb(2,2)*BB2(jj,J)+ jbb(2,3)*BB3(jj,J))*.125D0
        b123(jjj,j)  =(jbb(3,1)*BB1(jj,J)+jbb(3,2)*BB2(jj,J)+ jbb(3,3)*BB3(jj,J))*.125D0

!       Matriz Advectiva ADV:

        adv(j,jj)=b123(jjj-2,J)*uvwprom(1,ielem)+b123(jjj-1,J)*uvwprom(2,ielem)+b123(jjj,J)*uvwprom(3,ielem)

	  end do
    end do


    do j=1,8
      do jj=1,8
     
  	    dd(1,1) = dd11(jj,j)
   	    dd(1,2) = dd12(jj,j)
	    dd(1,3) = dd13(jj,j)
	    dd(2,1) = dd21(jj,j)
 	    dd(2,2) = dd22(jj,j)
	    dd(2,3) = dd23(jj,j)
	    dd(3,1) = dd31(jj,j)
 	    dd(3,2) = dd32(jj,j)
	    dd(3,3) = dd33(jj,j)
 
        add = matmul(matmul(jbb,dd),jbt)

        add(1,1) = add(1,1) +chg*hourg(jj,j)
	    add(2,2) = add(2,2) +chg*hourg(jj,j)
	    add(3,3) = add(3,3) +chg*hourg(jj,j)

        D11(j,jj)= (COEF1*add(1,1)+ COEF2*(add(2,2)+add(3,3)))/VOLU(IELEM)
        D22(j,jj)= (COEF1*add(2,2)+ COEF2*(add(3,3)+add(1,1)))/VOLU(IELEM)
        D33(j,jj)= (COEF1*add(3,3)+ COEF2*(add(1,1)+add(2,2)))/VOLU(IELEM)
        D12(j,jj)= (COEF2*add(1,2)+ COEF3*add(2,1))/VOLU(IELEM)
        D13(j,jj)= (COEF2*add(1,3)+ COEF3*add(3,1))/VOLU(IELEM)
        D23(j,jj)= (COEF2*add(2,3)+ COEF3*add(3,2))/VOLU(IELEM)
 
        BTD(j,jj)=(AUX1*add(1,1)+AUX2*add(2,2)+ AUX3*add(3,3)+AUX12*(add(1,2)+&
                   add(2,1))+AUX13*(add(1,3)+ add(3,1))+AUX23*(add(2,3)+add(3,2)))/VOLU(IELEM)

        BTE(j,jj)=(AUX1*add(1,1)+AUX2*add(2,2)+AUX3*add(3,3))/VOLU(IELEM)

!       DFP(j,jj)= AUX4*(add(1,1)+add(2,2)+ add(3,3))/VOLU(IELEM)
        DFP(j,jj) = 0.0d0

        BTN(j,jj)=(AUX12*(add(1,2)+add(2,1))+AUX13*(add(1,3)+add(3,1))+&
                   AUX23*(add(2,3)+add(3,2)))/VOLU(IELEM)


      
      end do
    end do      


   

    do i=1,8

      LRUNP=-dot_product(adv(:,i)+d11(:,i)+btd(:,i),uvwp(:,1))-dot_product(d12(:,i),uvwp(:,2)) &
            -dot_product(d13(:,i),uvwp(:,3))+(1.0D0/RHOInf)*dot_product(ag1(i,:),uvwp(:,4))

      LRVNP=-dot_product(adv(:,i)+d22(:,i)+btd(:,i),uvwp(:,2))-dot_product(d12(i,:),uvwp(:,1)) &
            -dot_product(d23(:,i),uvwp(:,3))+(1.0D0/RHOInf)*dot_product(ag2(i,:),uvwp(:,4)) 

      LRWNP=-dot_product(adv(:,i)+d33(:,i)+btd(:,i),uvwp(:,3))-dot_product(d13(i,:),uvwp(:,1)) &
            -dot_product(d23(i,:),uvwp(:,2))+(1.0D0/RHOInf)*dot_product(ag3(i,:),uvwp(:,4))

      LRPNP=-AuxADV*dot_product(adv(:,i),uvwp(:,4))-Aux*(dot_product(ag1(:,i),uvwp(:,1))+&
            dot_product(ag2(:,i),uvwp(:,2))+dot_product(ag3(:,i),uvwp(:,3)))-&
            dot_product(bte(:,i)+dfp(:,i)+btn(:,i),uvwp(:,4))

                   	   
!     ***** Termos de massa - 

!     MRU=MASD(1,IELEM)*uvwp(1,1)+MASD(2,IELEM)*uvwp(1,2)+&           
!     MASD(3,IELEM)*uvwp(1,3)+MASD(4,IELEM)*uvwp(1,4)+&
!     MASD(5,IELEM)*uvwp(1,5)+MASD(6,IELEM)*uvwp(1,6)+&
!     MASD(7,IELEM)*uvwp(1,7)+MASD(8,IELEM)*uvwp(1,8)

!     MRV=MASD(1,IELEM)*uvwp(2,1)+MASD(2,IELEM)*uvwp(2,2)+&
!     MASD(3,IELEM)*uvwp(2,3)+MASD(4,IELEM)*uvwp(2,4)+&
!     MASD(5,IELEM)*uvwp(2,5)+MASD(6,IELEM)*uvwp(2,6)+&
!     MASD(7,IELEM)*uvwp(2,7)+MASD(8,IELEM)*uvwp(2,8)

!     MRW=MASD(1,IELEM)*uvwp(3,1)+MASD(2,IELEM)*uvwp(3,2)+&
!     MASD(3,IELEM)*uvwp(3,3)+MASD(4,IELEM)*uvwp(3,4)+&
!     MASD(5,IELEM)*uvwp(3,5)+MASD(6,IELEM)*uvwp(3,6)+&
!     MASD(7,IELEM)*uvwp(3,7)+MASD(8,IELEM)*uvwp(3,8)

!     MRP=MASS(I,1)*uvwp(1,4)+MASS(I,2)*uvwp(2,4)+&      
!     MASS(I,3)*uvwp(3,4)+MASS(I,4)*uvwp(4,4)+&
!     MASS(I,5)*uvwp(5,4)+MASS(I,6)*uvwp(6,4)+&
!     MASS(I,7)*uvwp(7,4)+MASS(I,8)*uvwp(8,4)


      MRP=dot_product(masl(:,i),uvwp(:,4))
 
      rru(i) = Delt * 0.5D0 * LRUNP
      rrv(i) = Delt * 0.5D0 * LRVNP
      rrw(i) = Delt * 0.5D0 * LRWNP
      rrp(i) = MRP + Delt * 0.5D0 * LRPNP

   
!     Formacao dos vetores globais RRnG e RRPG representando o lado direito das
!     *************************************************************************
!     equacoes de conservacao:
!     ************************

    end do 


    do j=1,8
      rrg(1,kone(ii+j)) = rrg(1,kone(ii+j))+rru(j)
      rrg(2,kone(ii+j)) = rrg(2,kone(ii+j))+rrv(j)
      rrg(3,kone(ii+j)) = rrg(3,kone(ii+j))+rrw(j)
      rrg(4,kone(ii+j)) = rrg(4,kone(ii+j))+rrp(j)
    end do

  END DO



! FORMACAO DOS VETORES DE CARGAS EQUIVALENTES AS ACOES DE SUPERFICIE:
! *******************************************************************


! FFFpG = 0.0d0

  DO J=1,LCONTCARA

    IELEM=IEL(J)
    III= (J-1)*8 
 
    NK(1)= AREAC(J)*FCONTOR(III+1)
    NK(2)= AREAC(J)*FCONTOR(III+2)
    NK(3)= AREAC(J)*FCONTOR(III+3)
    NK(4)= AREAC(J)*FCONTOR(III+4)
    NK(5)= AREAC(J)*FCONTOR(III+5)
    NK(6)= AREAC(J)*FCONTOR(III+6)
    NK(7)= AREAC(J)*FCONTOR(III+7)
    NK(8)= AREAC(J)*FCONTOR(III+8)
       
    II= (IELEM-1)*8

    do jj=1,8
      cnort(1,jj)=cnorm(1,kone(ii+jj))
      cnort(2,jj)=cnorm(2,kone(ii+jj))
      cnort(3,jj)=cnorm(3,kone(ii+jj))
    end do
       
    N1= KONE(II+1)
    N2= KONE(II+2)
    N3= KONE(II+3)
    N4= KONE(II+4)
    N5= KONE(II+5)
    N6= KONE(II+6)
    N7= KONE(II+7)
    N8= KONE(II+8)

    do jj=1,8
      uvwp(jj,4)=uvw(4,kone(ii+jj))
    end do


!   Pressao media no elemento


    PRmd = sum(uvwp(:,4))*0.125D0


!   Constantes do tensor de balanco difusivo

    AUX1 = Delt*0.25D0*uvwprom(1,ielem)*uvwprom(1,ielem)
    AUX2 = Delt*0.25D0*uvwprom(2,ielem)*uvwprom(2,ielem)
    AUX3 = Delt*0.25D0*uvwprom(3,ielem)*uvwprom(3,ielem)
    AUX12 = Delt*0.25D0*uvwprom(1,ielem)*uvwprom(2,ielem)
    AUX13 = Delt*0.25D0*uvwprom(1,ielem)*uvwprom(3,ielem)
    AUX23 = Delt*0.25D0*uvwprom(2,ielem)*uvwprom(3,ielem)
    AUX4 = 0.0D0

!   Constantes de difusao

!    IF(IndTurb.NE.0) THEN

      AUXv = ViscCin + VDTurb(IELEM)   !viscosidade total

      IF(AUXv.LT.0.0D0) THEN
        AUXv = ViscCin
      end if

      COEF1= 2.0D0 * AUXv + ViscVol
      COEF2= AUXv 
      COEF3= ViscVol

!     bt = 1.0D6
!     AUX4 = VDTurb(IELEM)/bt

!    end if


    do jj=1,8
 
!     Termos de contorno da eq. de momentum

      IF (FCONTOR(III+jj)/=0.0D0)	THEN

        D11c = ( COEF1 + AUX1 ) * cnort(1,jj) * pp(1,IELEM) + &
               ( COEF2 + AUX2 ) * cnort(2,jj) * pp(4,IELEM) +&
               ( COEF2 + AUX3 ) * cnort(3,jj) * pp(7,IELEM) +&
        AUX12 * ( pp(1,IELEM) * cnort(2,jj) + pp(4,IELEM) *cnort(1,jj))+&
        AUX13 * ( pp(1,IELEM) * cnort(3,jj) + pp(7,IELEM) *cnort(1,jj))+ & 
        AUX23 * ( pp(4,IELEM) * cnort(3,jj) + pp(7,IELEM) *cnort(2,jj))

        D22c = ( COEF1 + AUX2 ) * cnort(2,jj) * pp(5,IELEM) + &
               ( COEF2 + AUX1 ) * cnort(1,jj) * pp(2,IELEM) +&
               ( COEF2 + AUX3 ) * cnort(3,jj) * pp(8,IELEM) + &   
        AUX12 * ( pp(2,IELEM) * cnort(2,jj) + pp(5,IELEM) *cnort(1,jj))+&
        AUX13 * ( pp(2,IELEM) * cnort(3,jj) + pp(8,IELEM) *cnort(1,jj))+ & 
        AUX23 * ( pp(5,IELEM) * cnort(3,jj) + pp(8,IELEM) *cnort(2,jj))

        D33c = ( COEF1 + AUX3 ) * cnort(3,jj) * pp(9,IELEM) + &
               ( COEF2 + AUX1 ) * cnort(1,jj) * pp(3,IELEM) +&
               ( COEF2 + AUX2 ) * cnort(2,jj) * pp(6,IELEM) + &   
        AUX12 * ( pp(3,IELEM) * cnort(2,jj) + pp(6,IELEM) *cnort(1,jj))+&
        AUX13 * ( pp(3,IELEM) * cnort(3,jj) + pp(9,IELEM) *cnort(1,jj))+ & 
        AUX23 * ( pp(6,IELEM) * cnort(3,jj) + pp(9,IELEM) *cnort(2,jj))

        D12c = COEF2*pp(2,IELEM)*cnort(2,jj) + COEF3*pp(5,IELEM)*cnort(1,jj)
        D21c = COEF2*pp(4,IELEM)*cnort(1,jj) + COEF3*pp(1,IELEM)*cnort(2,jj)
        D13c = COEF2*pp(3,IELEM)*cnort(3,jj) + COEF3*pp(9,IELEM)*cnort(1,jj)
        D31c = COEF2*pp(7,IELEM)*cnort(1,jj) + COEF3*pp(1,IELEM)*cnort(3,jj)
        D23c = COEF2*pp(6,IELEM)*cnort(3,jj) + COEF3*pp(9,IELEM)*cnort(2,jj)
        D32c = COEF2*pp(8,IELEM)*cnort(2,jj) + COEF3*pp(5,IELEM)*cnort(3,jj)
   
        FF1= D11c + D12c + D13c - PRmd*(1.0D0/RHOInf)*cnort(1,jj)
        FF2= D21c + D22c + D23c - PRmd*(1.0D0/RHOInf)*cnort(2,jj)
        FF3= D31c + D32c + D33c - PRmd*(1.0D0/RHOInf)*cnort(3,jj)
	 
        rrg(1,KONE(II+jj))=rrg(1,KONE(II+jj))+FF1*NK(jj)*Delt*0.5D0
        rrg(2,KONE(II+jj))=rrg(2,KONE(II+jj))+FF2*NK(jj)*Delt*0.5D0 
        rrg(3,KONE(II+jj))=rrg(3,KONE(II+jj))+FF3*NK(jj)*Delt*0.5D0

!       Termos de contorno da eq. de massa

!       FFpE = ( AUX4 + AUX1 ) * Gr1P(IELEM) * cnort(1,jj) + &
!              ( AUX4 + AUX2 ) * Gr2P(IELEM) * cnort(2,jj) + &
!              ( AUX4 + AUX3 ) * Gr3P(IELEM) * cnort(3,jj)

!       FFpN = AUX12*(Gr1P(IELEM)*cnort(2,jj) + Gr2P(IELEM)*cnort(1,jj))+
!              AUX13*(Gr1P(IELEM)*cnort(3,jj) + Gr3P(IELEM)*cnort(1,jj))+
!              AUX23*(Gr2P(IELEM)*cnort(3,jj) + Gr3P(IELEM)*cnort(2,jj))

!       FFFpG(KONE(II+jj))=FFFpG(KONE(II+jj))+(FFpE+FFpN)*NK1*Delt*0.5D0
	  
      END IF

	end do

  END DO



                   

! CALCULO DOS VETORES GLOBAIS Ua Va Wa e PRa
! ******************************************

  DO I=1,NNM

    Aux0 = 0.0D0

!   UVWa(1,I)= uvw(1,I) + ( RRUG(I)+FFFf1G(I) ) / MLUM(I)
!   UVWa(2,I)= uvw(2,I) + ( RRVG(I)+FFFf2G(I) ) / MLUM(I)
!   UVWa(3,I)= uvw(3,I) + ( RRWG(I)+FFFf3G(I) ) / MLUM(I)
!   UVWa(4,I) = ( RRPG(I)+FFFpG(I)*Aux0 ) / MLUM(I)

    UVWa(1,I)= uvw(1,I) + RRG(1,I)/MLUM(I)
    UVWa(2,I)= uvw(2,I) + RRG(2,I)/MLUM(I)
    UVWa(3,I)= uvw(3,I) + RRG(3,I)/MLUM(I)
    UVWa(4,I) = RRG(4,I)/MLUM(I)

  END DO



! Correcao das componentes de velocidade 
! **************************************

  crug = 0.0d0
  crvg = 0.0d0
  crwg = 0.0d0
	
  DO IELEM=1,NELEM 
  
    II= (IELEM-1)*8
      
    N1= KONE(II+1)
    N2= KONE(II+2)
    N3= KONE(II+3)
    N4= KONE(II+4)
    N5= KONE(II+5)
    N6= KONE(II+6)
    N7= KONE(II+7)
    N8= KONE(II+8)

    JB11= JIN(1,IELEM) * DETERJ(IELEM)
    JB12= JIN(2,IELEM) * DETERJ(IELEM)
    JB13= JIN(3,IELEM) * DETERJ(IELEM)
    JB21= JIN(4,IELEM) * DETERJ(IELEM)
    JB22= JIN(5,IELEM) * DETERJ(IELEM)
    JB23= JIN(6,IELEM) * DETERJ(IELEM)
    JB31= JIN(7,IELEM) * DETERJ(IELEM)
    JB32= JIN(8,IELEM) * DETERJ(IELEM)
    JB33= JIN(9,IELEM) * DETERJ(IELEM)

    do j=1,8
      uvwp(j,4) = uvwa(4,kone(ii+j))-uvw(4,kone(ii+j))
    end do

!   Matrizes auxiliares Bi utilizadas no termo de adveccao e de contorno de difusao:

    do j=1,8
      do jj=1,8
        jjj=jj*3

        b123(jjj-2,j)=(JB11*BB1(jj,J)+JB12*BB2(jj,J)+JB13*BB3(jj,J))*.125D0
        b123(jjj-1,j)=(JB21*BB1(jj,J)+JB22*BB2(jj,J)+JB23*BB3(jj,J))*.125D0
        b123(jjj,j)  =(JB31*BB1(jj,J)+JB32*BB2(jj,J)+JB33*BB3(jj,J))*.125D0

      end do
	end do

    do i=1,8

      jj=i*3-2
      cru(i)=dot_product(b123(jj,:),uvwp(:,4))
      jj=jj+1
      crv(i)=dot_product(b123(jj,:),uvwp(:,4))
      jj=jj+1
      crw(i)=dot_product(b123(jj,:),uvwp(:,4))

    end do

    crug(n1) = crug(n1) + cru(1) 
    crug(n2) = crug(n2) + cru(2) 
    crug(n3) = crug(n3) + cru(3) 
    crug(n4) = crug(n4) + cru(4) 
    crug(n5) = crug(n5) + cru(5) 
    crug(n6) = crug(n6) + cru(6) 
    crug(n7) = crug(n7) + cru(7) 
    crug(n8) = crug(n8) + cru(8) 

    crvg(n1) = crvg(n1) + crv(1) 
    crvg(n2) = crvg(n2) + crv(2) 
    crvg(n3) = crvg(n3) + crv(3) 
    crvg(n4) = crvg(n4) + crv(4) 
    crvg(n5) = crvg(n5) + crv(5) 
    crvg(n6) = crvg(n6) + crv(6) 
    crvg(n7) = crvg(n7) + crv(7) 
    crvg(n8) = crvg(n8) + crv(8) 

    crwg(n1) = crwg(n1) + crw(1) 
    crwg(n2) = crwg(n2) + crw(2) 
    crwg(n3) = crwg(n3) + crw(3) 
    crwg(n4) = crwg(n4) + crw(4) 
    crwg(n5) = crwg(n5) + crw(5) 
    crwg(n6) = crwg(n6) + crw(6) 
    crwg(n7) = crwg(n7) + crw(7) 
    crwg(n8) = crwg(n8) + crw(8) 

  end do

  do I=1,NNM
    UVWa(1,I)=UVWa(1,I) - (1.0D0/RHOInf)*((Delt)/4.0D0)*CRUG(i)
    UVWa(2,I)=UVWa(2,I) - (1.0D0/RHOInf)*((Delt)/4.0D0)*CRVG(i)
    UVWa(3,I)=UVWa(3,I) - (1.0D0/RHOInf)*((Delt)/4.0D0)*CRWG(i)
  end do

! Aplicacao das condicoes de contorno sobre as variaveis de campo
! ***************************************************************

  IF (NBCU.NE.0) THEN
    DO I= 1,NBCU
      NNO= IBCU(I)
      UVWa(1,NNO)= BCU(I) 
    END DO
  END IF

  IF (NBCV.NE.0) THEN
    DO I= 1,NBCV
      NNO= IBCV(I)
      UVWa(2,NNO)= BCV(I)      
    END DO
  END IF

  IF (NBCW.NE.0) THEN
    DO I= 1,NBCW
      NNO= IBCW(I)
      UVWa(3,NNO)= BCW(I) 
    END DO
  END IF

  IF (NBCP.NE.0) THEN
    DO I= 1,NBCP
      NNO= IBCP(I)
      UVWa(4,NNO)= BCP(I) 
    END DO
  END IF



! ******************************************************************************
! ***** FIM DO PRIMEIRO PASSO***************************************************
! ******************************************************************************






! ******************************************************************************
! ***** SEGUNDO PASSO **********************************************************
! ******************************************************************************


! Derivadas de velocidade Pij e componentes de velocidade media no elemento
! *************************************************************************
! UPROM,VPROM e WPROM e gradientes de pressao GriP:
! *************************************************

  DO IELEM=1,NELEM

    JINV11= JIN(1,IELEM)
    JINV12= JIN(2,IELEM)
    JINV13= JIN(3,IELEM)
    JINV21= JIN(4,IELEM)
    JINV22= JIN(5,IELEM)
    JINV23= JIN(6,IELEM)
    JINV31= JIN(7,IELEM)
    JINV32= JIN(8,IELEM)
    JINV33= JIN(9,IELEM)

    II= (IELEM-1)*8

    do j=1,8
	  uvwp(j,1) = uvwa(1,kone(ii+j))
	  uvwp(j,2) = uvwa(2,kone(ii+j))
	  uvwp(j,3) = uvwa(3,kone(ii+j))
	  uvwp(j,4) = uvwa(4,kone(ii+j))
    end do
       
    N1= KONE(II+1)
    N2= KONE(II+2)
    N3= KONE(II+3)
    N4= KONE(II+4)
    N5= KONE(II+5)
    N6= KONE(II+6)
    N7= KONE(II+7)
    N8= KONE(II+8)

!   ***************************************************************************
!   Somatorio dos produtos dos elementos da matriz jacobiana pelas
!   coordenadas naturais; para serem usados em produtos formados abaixo. 
!   ***************************************************************************

	JCN(1,1) = -JINV11-JINV12-JINV13
    JCN(2,1) =  JINV11-JINV12-JINV13 
    JCN(3,1) =  JINV11+JINV12-JINV13 
    JCN(4,1) = -JINV11+JINV12-JINV13
    JCN(5,1) = -JINV11-JINV12+JINV13
    JCN(6,1) =  JINV11-JINV12+JINV13 
    JCN(7,1) =  JINV11+JINV12+JINV13 
    JCN(8,1) = -JINV11+JINV12+JINV13
     
    JCN(1,2) = -JINV21-JINV22-JINV23
    JCN(2,2) =  JINV21-JINV22-JINV23 
    JCN(3,2) =  JINV21+JINV22-JINV23 
    JCN(4,2) = -JINV21+JINV22-JINV23
    JCN(5,2) = -JINV21-JINV22+JINV23
    JCN(6,2) =  JINV21-JINV22+JINV23 
    JCN(7,2) =  JINV21+JINV22+JINV23 
    JCN(8,2) = -JINV21+JINV22+JINV23
     
    JCN(1,3) = -JINV31-JINV32-JINV33 
    JCN(2,3) =  JINV31-JINV32-JINV33
    JCN(3,3) =  JINV31+JINV32-JINV33       
    JCN(4,3) = -JINV31+JINV32-JINV33
    JCN(5,3) = -JINV31-JINV32+JINV33
    JCN(6,3) =  JINV31-JINV32+JINV33
    JCN(7,3) =  JINV31+JINV32+JINV33
    JCN(8,3) = -JINV31+JINV32+JINV33


!   **************************************************************************
!   Media dos valores nodais dos produtos das coord naturais pela
!   inversa de J(0) e pelas componentes de velocidade; CC Newman
!   **************************************************************************
      
    pp(1,IELEM)=dot_product(jcn(:,1),uvwp(:,1))*0.125D0
    pp(4,IELEM)=dot_product(jcn(:,2),uvwp(:,1))*0.125D0
    pp(7,IELEM)=dot_product(jcn(:,3),uvwp(:,1))*0.125D0
    pp(2,IELEM)=dot_product(jcn(:,1),uvwp(:,2))*0.125D0
    pp(5,IELEM)=dot_product(jcn(:,2),uvwp(:,2))*0.125D0
    pp(8,IELEM)=dot_product(jcn(:,3),uvwp(:,2))*0.125D0
    pp(3,IELEM)=dot_product(jcn(:,1),uvwp(:,3))*0.125D0
    pp(6,IELEM)=dot_product(jcn(:,2),uvwp(:,3))*0.125D0
    pp(9,IELEM)=dot_product(jcn(:,3),uvwp(:,3))*0.125D0

!   ***************************************************************************
!   Media dos valores nodais dos produtos das coordenadas naturais,
!   a inversa de J(0) e pelas tres componentes de velocidade; CC Newman
!   ***************************************************************************

!   divergência da velocidade (dVk/dXk):
       
    PNV(IELEM) = pp(1,IELEM) + pp(5,IELEM) + pp(9,IELEM)

!   gradientes de pressao (dp/dXk):

    Gr1P(IELEM) = dot_product(jcn(:,1),uvwp(:,4))*0.125D0
    Gr2P(IELEM) = dot_product(jcn(:,2),uvwp(:,4))*0.125D0
    Gr3P(IELEM) = dot_product(jcn(:,3),uvwp(:,4))*0.125D0

!   ***************************************************************************
!   Media dos valores nodais das componentes de velocidade.
!   ***************************************************************************

    uvwprom(1,ielem)= sum(uvwp(:,1))*0.125D0
    uvwprom(2,ielem)= sum(uvwp(:,2))*0.125D0
    uvwprom(3,ielem)= sum(uvwp(:,3))*0.125D0

  END DO


! Formacao dos vetores locais RRn e RRP representando o lado direito das
! **********************************************************************
! equacoes de conservacao:
! ************************


  rrg = 0.0d0

  DO IELEM=1,NELEM 
  
    II= (IELEM-1)*8
       
    do j=1,8
	  uvwp(j,1) = uvwa(1,kone(ii+j))
	  uvwp(j,2) = uvwa(2,kone(ii+j))
	  uvwp(j,3) = uvwa(3,kone(ii+j))
	  uvwp(j,4) = uvwa(4,kone(ii+j))
    end do

	do j=1,8
	  uvwt(j,4) = uvw(4,kone(ii+j))
    end do


!   ***************************
!   Matriz de massa consistente
!   ***************************

    VC= VOLU(IELEM)*0.015625D0

    mass=vc*transpose(mmm)
	masd=volu(ielem)/8.0d0
	masl=(1.0d0-elump)*mass
	do i=1,8
	  masl(i,i)=masl(i,i)+elump*masd
	end do


    Aux = RHOInf*(VelSom**2)
    AuxADV = 0.0D0


!   ***** Controle de modos espurios (HOURGLASSING) - CHG:

!   ***************************************************************************
!                         Chg dado pelo criterio CHRISTON 
!   ***************************************************************************

    CHG= (VOLU(IELEM)**(1.0D0/3.D0))*1.0D0*0.5D0*VOLU(IELEM)*(1.0D0-(VNVX**(1.D0/3.D0)))*CONTROL1

    jbb(1,1) = JIN(1,IELEM) * DETERJ(IELEM)
    jbb(1,2) = JIN(2,IELEM) * DETERJ(IELEM)
    jbb(1,3) = JIN(3,IELEM) * DETERJ(IELEM)
    jbb(2,1) = JIN(4,IELEM) * DETERJ(IELEM)
    jbb(2,2) = JIN(5,IELEM) * DETERJ(IELEM)
    jbb(2,3) = JIN(6,IELEM) * DETERJ(IELEM)
    jbb(3,1) = JIN(7,IELEM) * DETERJ(IELEM)
    jbb(3,2) = JIN(8,IELEM) * DETERJ(IELEM)
    jbb(3,3) = JIN(9,IELEM) * DETERJ(IELEM)

    jbt=transpose(jbb)

!   *************************************************************************************
!   Matrizes de termos do tensor de balanco difusivo e de turbulencia na equacao de massa 
!
!   Eq. de quantidade de mov.:
!     BTD - adicionada aos termos volumetricos das matrizes difusivas (Dii):
!
!   Eq. da pressao:
!     BTE - termos volumetricos do tensor difusivo
!     BTN - termos cruzados do tensor de balanco difusivo
!     DFP - termos volumetricos de turbulencia - laplaciano da pressao
!   *************************************************************************************


!   Matriz do termo de pressao integrado por partes AGi:

    do j=1,8
      do jj=1,8
	
	    AG1(jj,J)=(jbb(1,1)*TB1(jj,J)+jbb(1,2)*TB2(jj,J)+jbb(1,3)*TB3(jj,J))*.125D0
        AG2(jj,J)=(jbb(2,1)*TB1(jj,J)+jbb(2,2)*TB2(jj,J)+ jbb(2,3)*TB3(jj,J))*.125D0
        AG3(jj,J)=(jbb(3,1)*TB1(jj,J)+jbb(3,2)*TB2(jj,J)+ jbb(3,3)*TB3(jj,J))*.125D0

	  end do
	end do


!   Matrizes auxiliares Bi utilizadas no termo de adveccao e de contorno de difusao:

    do J=1,8
      do jj=1,8
        jjj=jj*3
        b123(jjj-2,j)=(jbb(1,1)*BB1(jj,J)+jbb(1,2)*BB2(jj,J)+ jbb(1,3)*BB3(jj,J))*.125D0
        b123(jjj-1,j)=(jbb(2,1)*BB1(jj,J)+jbb(2,2)*BB2(jj,J)+ jbb(2,3)*BB3(jj,J))*.125D0
        b123(jjj,j)  =(jbb(3,1)*BB1(jj,J)+jbb(3,2)*BB2(jj,J)+ jbb(3,3)*BB3(jj,J))*.125D0

!       Matriz Advectiva ADV:

        adv(j,jj)=b123(jjj-2,J)*uvwprom(1,ielem)+ b123(jjj-1,J)*uvwprom(2,ielem)+&
                  b123(jjj,J)*uvwprom(3,ielem)

	  end do

    end do

    if(IndTurb.eq.0) then
	  aux1 = ViscCin
	else 
      AUX1 = ViscCin + VDTurb(IELEM)
      IF(AUX1.LT.0.0D0) AUX1 = ViscCin
    end if

    COEF1= 2.0D0 * AUX1 + ViscVol
    COEF2= AUX1 
    COEF3= ViscVol         

    do j=1,8
      do jj=1,8

        dd(1,1) = dd11(jj,j)
 	    dd(1,2) = dd12(jj,j)
	    dd(1,3) = dd13(jj,j)
	    dd(2,1) = dd21(jj,j)
 	    dd(2,2) = dd22(jj,j)
	    dd(2,3) = dd23(jj,j)
	    dd(3,1) = dd31(jj,j)
 	    dd(3,2) = dd32(jj,j)
	    dd(3,3) = dd33(jj,j)
 
        add = matmul(matmul(jbb,dd),jbt)

        add(1,1) = add(1,1) +chg*hourg(jj,j)
	    add(2,2) = add(2,2) +chg*hourg(jj,j)
	    add(3,3) = add(3,3) +chg*hourg(jj,j)

        D11(j,jj)= (COEF1*add(1,1)+ COEF2*(add(2,2)+add(3,3)))/VOLU(IELEM)
        D22(j,jj)= (COEF1*add(2,2)+ COEF2*(add(3,3)+add(1,1)))/VOLU(IELEM)
        D33(j,jj)= (COEF1*add(3,3)+ COEF2*(add(1,1)+add(2,2)))/VOLU(IELEM)
        D12(j,jj)= (COEF2*add(1,2)+ COEF3*add(2,1))/VOLU(IELEM)
        D13(j,jj)= (COEF2*add(1,3)+ COEF3*add(3,1))/VOLU(IELEM)
        D23(j,jj)= (COEF2*add(2,3)+ COEF3*add(3,2))/VOLU(IELEM)
        DFP(j,jj)=(AUX4*(add(1,1)+add(2,2)+ add(3,3)))/VOLU(IELEM)

      end do
    end do

    do i=1,8

      LRUNP=-dot_product(adv(:,i)+d11(:,i),uvwp(:,1))-dot_product(d12(:,i),uvwp(:,2)) &
            -dot_product(d13(:,i),uvwp(:,3)) +(1.0D0/RHOInf)*dot_product(ag1(i,:),uvwp(:,4))
 
      LRVNP=-dot_product(adv(:,i)+d22(:,i),uvwp(:,2))-dot_product(d12(i,:),uvwp(:,1))&
            -dot_product(d23(:,i),uvwp(:,3))+(1.0D0/RHOInf)*dot_product(ag2(i,:),uvwp(:,4))

      LRWNP=-dot_product(adv(:,i)+d33(:,i),uvwp(:,3)) -dot_product(d13(i,:),uvwp(:,1))&
            -dot_product(d23(i,:),uvwp(:,2))+(1.0D0/RHOInf)*dot_product(ag3(i,:),uvwp(:,4))
 
      LRPNP=-AuxADV*dot_product(adv(:,i),uvwp(:,4)) -Aux*(dot_product(ag1(:,i),uvwp(:,1))+&
                    dot_product(ag2(:,i),uvwp(:,2))+ dot_product(ag3(:,i),uvwp(:,3)))-&
                    dot_product(dfp(:,i),uvwp(:,4))
                   	   
!     ***** Termos de massa - 

!     MRU=MASD(1,IELEM)*uvw(1,n1)+MASD(2,IELEM)*uvw(1,n2)+ &
!         MASD(3,IELEM)*uvw(1,n3)+MASD(4,IELEM)*uvw(1,n4)+ &
!         MASD(5,IELEM)*uvw(1,n5)+MASD(6,IELEM)*uvw(1,n6)+ &
!         MASD(7,IELEM)*uvw(1,n7)+MASD(8,IELEM)*uvw(1,n8)

!     MRV=MASD(1,IELEM)*uvw(2,n1)+MASD(2,IELEM)*uvw(2,n2)+&           
!         MASD(3,IELEM)*uvw(2,n3)+MASD(4,IELEM)*uvw(2,n4)+&
!         MASD(5,IELEM)*uvw(2,n5)+MASD(6,IELEM)*uvw(2,n6)+&
!         MASD(7,IELEM)*uvw(2,n7)+MASD(8,IELEM)*uvw(2,n8)

!     MRW=MASD(1,IELEM)*uvw(3,n1)+MASD(2,IELEM)*uvw(3,n2)+ &          
!         MASD(3,IELEM)*uvw(3,n3)+MASD(4,IELEM)*uvw(3,n4)+ &
!         MASD(5,IELEM)*uvw(3,n5)+MASD(6,IELEM)*uvw(3,n6)+ &
!         MASD(7,IELEM)*uvw(3,n7)+MASD(8,IELEM)*uvw(3,n8)

!     MRP=MASS(1,I)*uvw(4,n1)+MASS(2,I)*uvw(4,n2)+ &    
!         MASS(3,I)*uvw(4,n3)+MASS(4,I)*uvw(4,n4)+ &
!         MASS(5,I)*uvw(4,n5)+MASS(6,I)*uvw(4,n6)+ &
!         MASS(7,I)*uvw(4,n7)+MASS(8,I)*uvw(4,n8)


      MRP=dot_product(masl(:,i),uvwt(:,4))

!     RRU(i,IELEM) = MRU*0.0D0 + Delt * LRUNP
!     RRV(i,IELEM) = MRV*0.0D0 + Delt * LRVNP
!     RRW(i,IELEM) = MRW*0.0D0 + Delt * LRWNP

      rru(i)=Delt * LRUNP
      rrv(i)=Delt * LRVNP
      rrw(i)=Delt * LRWNP
      rrp(i)=MRP + Delt * LRPNP

    end do

    do j=1,8
      rrg(1,kone(ii+j)) = rrg(1,kone(ii+j)) + rru(j)
      rrg(2,kone(ii+j)) = rrg(2,kone(ii+j)) + rrv(j)
      rrg(3,kone(ii+j)) = rrg(3,kone(ii+j)) + rrw(j)
      rrg(4,kone(ii+j)) = rrg(4,kone(ii+j)) + rrp(j)
    end do

  END DO


! FORMACAO DOS VETORES DE CARGAS EQUIVALENTES AS ACOES DE SUPERFICIE:
! *******************************************************************

! FFFpG = 0.0d0

  DO J=1,LCONTCARA

    IELEM=IEL(J)
    III= (J-1)*8 
 
    NK(1)= AREAC(J)*FCONTOR(III+1)
    NK(2)= AREAC(J)*FCONTOR(III+2)
    NK(3)= AREAC(J)*FCONTOR(III+3)
    NK(4)= AREAC(J)*FCONTOR(III+4)
    NK(5)= AREAC(J)*FCONTOR(III+5)
    NK(6)= AREAC(J)*FCONTOR(III+6)
    NK(7)= AREAC(J)*FCONTOR(III+7)
    NK(8)= AREAC(J)*FCONTOR(III+8)
       
    II= (IELEM-1)*8

    N1= KONE(II+1)
    N2= KONE(II+2)
    N3= KONE(II+3)
    N4= KONE(II+4)
    N5= KONE(II+5)
    N6= KONE(II+6)
    N7= KONE(II+7)
    N8= KONE(II+8)

!   Pressao media no elemento

    PRmd = (UVWa(4,N1)+UVWa(4,N2)+UVWa(4,N3)+UVWa(4,N4)+UVWa(4,N5)+&
            UVWa(4,N6)+UVWa(4,N7)+UVWa(4,N8))*0.125D0

    AUX4 = 0.0D0

    do jj=1,8
      cnort(1,jj)=cnorm(1,kone(ii+jj))
      cnort(2,jj)=cnorm(2,kone(ii+jj))
      cnort(3,jj)=cnorm(3,kone(ii+jj))
    end do


!   Constantes de difusao

!    IF(IndTurb.NE.0) THEN

      AUXv = ViscCin + VDTurb(IELEM)   !viscosidade total

      IF(AUXv.LT.0.0D0) THEN
        AUXv = ViscCin
       end if

      COEF1= 2.0D0 * AUXv + ViscVol
      COEF2= AUXv 
      COEF3= ViscVol

!     bt = 1.0D6
!     AUX4 = VDTurb(IELEM)/bt

!    end if

    do jj=1,8

      IF (FCONTOR(III+jj)/=0.0D0)	THEN

!       Termos de contorno da eq. de momentum

        D11c = COEF1*cnort(1,jj)*pp(1,IELEM)+COEF2*cnort(2,jj)*pp(4,IELEM)+ &
               COEF2*cnort(3,jj)*pp(7,IELEM)

        D22c = COEF1*cnort(2,jj)*pp(5,IELEM)+COEF2*cnort(1,jj)*pp(2,IELEM)+ &
               COEF2*cnort(3,jj)*pp(8,IELEM)

        D33c = COEF1*cnort(3,jj)*pp(9,IELEM)+COEF2*cnort(1,jj)*pp(3,IELEM)+ &
               COEF2*cnort(2,jj)*pp(6,IELEM)

        D12c = COEF2*pp(2,IELEM)*cnort(2,jj) + COEF3*pp(5,IELEM)*cnort(1,jj)
        D21c = COEF2*pp(4,IELEM)*cnort(1,jj) + COEF3*pp(1,IELEM)*cnort(2,jj)
        D13c = COEF2*pp(3,IELEM)*cnort(3,jj) + COEF3*pp(9,IELEM)*cnort(1,jj)
        D31c = COEF2*pp(7,IELEM)*cnort(1,jj) + COEF3*pp(1,IELEM)*cnort(3,jj)
        D23c = COEF2*pp(6,IELEM)*cnort(3,jj) + COEF3*pp(9,IELEM)*cnort(2,jj)
        D32c = COEF2*pp(8,IELEM)*cnort(2,jj) + COEF3*pp(5,IELEM)*cnort(3,jj)
	   
        FF1= D11c + D12c + D13c - PRmd*(1.0D0/RHOInf)*cnort(1,jj)
        FF2= D21c + D22c + D23c - PRmd*(1.0D0/RHOInf)*cnort(2,jj)
        FF3= D31c + D32c + D33c - PRmd*(1.0D0/RHOInf)*cnort(3,jj)
 
        rrg(1,kone(ii+jj))=rrg(1,kone(ii+jj)) + FF1*NK(jj)*Delt
        rrg(2,kone(ii+jj))=rrg(2,kone(ii+jj)) + FF2*NK(jj)*Delt 
        rrg(3,kone(ii+jj))=rrg(3,kone(ii+jj)) + FF3*NK(jj)*Delt

!       Termos de contorno da eq. de massa

!       FFpE = AUX4 * Gr1P(IELEM) * cnort(1,jj) + AUX4 * Gr2P(IELEM) * cnort(2,jj) + &
!              AUX4 * Gr3P(IELEM) * cnort(3,jj)

!       FFFpg(kone(ii+jj))=FFFpG(kone(ii+jj)) + FFpE*NK1*Delt
	  
      END IF

    end do

  END DO


                    

! CALCULO DOS VETORES GLOBAIS U V W e PR: 
! ***************************************

  DO I=1,NNM

    Aux0 = 0.0D0

!   uvw(1,I)= uvw(1,I) + ( RRUG(I)+FFFf1G(I) ) / MLUM(I)
!   uvw(2,I)= uvw(2,I) + ( RRVG(I)+FFFf2G(I) ) / MLUM(I)
!   uvw(3,I)= uvw(3,I) + ( RRWG(I)+FFFf3G(I) ) / MLUM(I)
!   uvw(4,I) = ( RRPG(I)+FFFpG(I)*Aux0 ) / MLUM(I)


    uvw(1,I)= uvw(1,I) + ( RRG(1,I) ) / MLUM(I)
    uvw(2,I)= uvw(2,I) + ( RRG(2,I) ) / MLUM(I)
    uvw(3,I)= uvw(3,I) + ( RRG(3,I) ) / MLUM(I)
    uvw(4,I) = RRG(4,I) / MLUM(I)

  END DO


! Aplicacao das condicoes de contorno sobre as variaveis de campo
! ***************************************************************

  IF (NBCU.NE.0) THEN
    DO I= 1,NBCU
      NNO= IBCU(I)
      uvw(1,NNO)= BCU(I) 
    END DO
  END IF

  IF (NBCV.NE.0) THEN
    DO I= 1,NBCV
      NNO= IBCV(I)
      uvw(2,NNO)= BCV(I)      
    END DO
  END IF

  IF (NBCW.NE.0) THEN
    DO I= 1,NBCW
      NNO= IBCW(I)
      uvw(3,NNO)= BCW(I) 
    END DO
  END IF

  IF (NBCP.NE.0) THEN
    DO I= 1,NBCP
      NNO= IBCP(I)
      uvw(4,NNO)= BCP(I) 
    END DO
  END IF

! Retencao da pressao final para o calculo do residuo:

  Prf = uvw(4,:)


! SUAVIZACAO DE PRESSOES E DA VISCOSIDADE TURBULENTA
! **************************************************

  DO IELEM= 1,NELEM 

    II= (IELEM-1)*8
      
    N1= KONE(II+1)
    N2= KONE(II+2)
    N3= KONE(II+3)
    N4= KONE(II+4)
    N5= KONE(II+5)
    N6= KONE(II+6)
    N7= KONE(II+7)
    N8= KONE(II+8)

    PPROME(IELEM)=(uvw(4,n1)+uvw(4,n2)+uvw(4,n3)+uvw(4,n4)+uvw(4,n5)+&
                   uvw(4,n6)+uvw(4,n7)+uvw(4,n8))*VOLU(IELEM)* 0.015625D0

    Vtaux(IELEM)= VDTurb(IELEM) * VOLU(IELEM) * 0.125D0

  END DO

  DO I=1,NNM

    K10 = NEIBOR(1,I)
    K20 = NEIBOR(2,I)
    K30 = NEIBOR(3,I)
    K40 = NEIBOR(4,I)
    K50 = NEIBOR(5,I)
    K60 = NEIBOR(6,I)
    K70 = NEIBOR(7,I)
    K80 = NEIBOR(8,I)
    K90 = NEIBOR(9,I)
    K100 = NEIBOR(10,I)
    K110 = NEIBOR(11,I)
    K120 = NEIBOR(12,I)
    K130 = NEIBOR(13,I)
    K140 = NEIBOR(14,I)
    K150 = NEIBOR(15,I)
    K160 = NEIBOR(16,I)

    PRSuav(I)=( PPROME(K10)+PPROME(K20)+PPROME(K30)+PPROME(K40)+PPROME(K50)+PPROME(K60)+&
                PPROME(K70)+PPROME(K80)+PPROME(K90)+PPROME(K100)+PPROME(K110)+PPROME(K120)+&
                PPROME(K130)+PPROME(K140)+PPROME(K150)+PPROME(K160) ) / MLUM(I)

    VtSuav(I)=( Vtaux(K10)+Vtaux(K20)+Vtaux(K30)+Vtaux(K40)+Vtaux(K50)+Vtaux(K60)+&
                Vtaux(K70)+Vtaux(K80)+Vtaux(K90)+Vtaux(K100)+Vtaux(K110)+Vtaux(K120)+&
                Vtaux(K130)+Vtaux(K140)+Vtaux(K150)+Vtaux(K160) ) / MLUM(I)

  END DO
 

! CALCULO DOS CAMPOS MEDIOS DE VELOCIDADE, PRESSAO E VISCOSIDADE PARA UM 
! **********************************************************************
! DETERMINADO PERIODO DE TEMPO (PASSMED)
! **************************************         


  IF(NPASS.NE.NNAUX.and.NNNRRR.GT.NPASS) THEN
    DO I=1,NNM
      Umed(I) = Umed(I) + uvw(1,I)/PASSMED
      Vmed(I) = Vmed(I) + uvw(2,I)/PASSMED
      Wmed(I) = Wmed(I) + uvw(3,I)/PASSMED
      PRmed(I) = PRmed(I) + uvw(4,I)/PASSMED
      VTmed(I) = VTmed(I) + VtSuav(I)/PASSMED
      PRSmed(I) = PRSmed(I) + PRSuav(I)/PASSMED
    end do
  end if

! ************************************
! CONDICAO PARA GRAVACAO DE RESULTADOS
! ************************************

  IF(NNR.EQ.NIR) THEN

!   ***************************************************************************
!   NNR e um contador dentro do laco de tempo que armazena a quantidade de
!   vezes que entrou-se neste ciclo e NIR e a quantidade de passos de tempo
!   entre registros. Logo se NNR nao e igual a NIR, nao se efetua nenhuma
!   gravacao. Se for igual, e feita a gravacao, soma-se um a NFILE, e zera-se 
!   NNR, para o proximo registro.        
!   ***************************************************************************

!   GRAVACAO DE RESULTADOS
!   **********************

    NFILE=NFILE+1
    NNR=0
    TPOAC=TEMPO

    NF1=NFILE/100+48
    NF2=(NFILE-(NFILE/100)*100)/10+48
    NF3=NFILE-(NFILE/10)*10+48
    PTNAME(7:9)=CHAR(NF1)//CHAR(NF2)//CHAR(NF3)

!   ARQUIVO '.v': VALORES DE "U" "V" E "W" NOS NOS 
!   --------------------------------------------------
!   ARQUIVO '.pr': VALORES DE "PR" E "PRSuav" NOS NOS
!   -----------------------------------------------------
!   ARQUIVO '.vts': VALORES DE "VtSuav" NOS NOS
!   -----------------------------------------------

    PTNAME(10:14)='.v   '
    OPEN(2,FILE=PTNAME,STATUS='UNKNOWN')
    REWIND (2)

    PTNAME(10:14)='.pr  '
    OPEN(3,FILE=PTNAME,STATUS='UNKNOWN')
    REWIND (3)

    PTNAME(10:14)='.vts '
    OPEN(4,FILE=PTNAME,STATUS='UNKNOWN')
    REWIND (4)

    DO I=1,NNM
      WRITE(2,*) uvw(1,I),uvw(2,I),uvw(3,I)
      WRITE(3,*) uvw(4,I),PRSuav(I)
      WRITE(4,*) VtSuav(I)

    END DO

    CLOSE(2)
    CLOSE(3)
    CLOSE(4)


  end if

! ***********************************************************************
! FIM DA CONDICAO PARA CALCULO DE COEFICENTES E GRAVACAO DE RESULTADOS 
! ***********************************************************************



  IF(MNOBJ.NE.0) THEN

    IF(NCOUNTCOEF.EQ.NCOEF) THEN

      NCOUNTCOEF=0
      NARQ=29

      DO NOB = 1,NOBJ

!       ****************************************
!       CALCULO DE COEFICIENTES AERODINAMICOS
!       ****************************************

        Fx = 0.0D0
        Fy = 0.0D0
        Fz = 0.0D0
        Mx = 0.0D0
        My = 0.0D0
        Mz = 0.0D0

!       *******************************************************************************
!       Procedimentos referentes à obtenção dos valores nodais das componentes
!       do tensor de tensões viscosas para elementos de fluido em contato com
!       contornos sólidos. Suavização	do centro para os nós pertencentes à face
!       de contorno correspondente.
!       *******************************************************************************

!       Cálculo das tensões nos elementos de contorno sólido: 

        DO J=1,NFCS

          IF(J.GT.NFCN(NOB)) EXIT

          IELEM=IELCS(NOB,J)
  
          JINV11= JIN(1,IELEM)
          JINV12= JIN(2,IELEM)
          JINV13= JIN(3,IELEM)
          JINV21= JIN(4,IELEM)
          JINV22= JIN(5,IELEM)
          JINV23= JIN(6,IELEM)
          JINV31= JIN(7,IELEM)
          JINV32= JIN(8,IELEM)
          JINV33= JIN(9,IELEM)

          II= (IELEM-1)*8
        
          N1= KONE(II+1)
          N2= KONE(II+2)
          N3= KONE(II+3)
          N4= KONE(II+4)
          N5= KONE(II+5)
          N6= KONE(II+6)
          N7= KONE(II+7)
          N8= KONE(II+8)

!         ***************************************************************************
!   	  Somatorio dos produtos dos elementos da inversa da matriz jacobiana 
!         pelas coordenadas naturais; para serem usados em produtos formados abaixo. 
!         ***************************************************************************
      
          JCN11 = -JINV11-JINV12-JINV13
          JCN12 =  JINV11-JINV12-JINV13 
          JCN13 =  JINV11+JINV12-JINV13 
          JCN14 = -JINV11+JINV12-JINV13
          JCN15 = -JINV11-JINV12+JINV13
          JCN16 =  JINV11-JINV12+JINV13 
          JCN17 =  JINV11+JINV12+JINV13 
          JCN18 = -JINV11+JINV12+JINV13
     
          JCN21 = -JINV21-JINV22-JINV23
          JCN22 =  JINV21-JINV22-JINV23 
          JCN23 =  JINV21+JINV22-JINV23 
          JCN24 = -JINV21+JINV22-JINV23
          JCN25 = -JINV21-JINV22+JINV23
          JCN26 =  JINV21-JINV22+JINV23 
          JCN27 =  JINV21+JINV22+JINV23 
          JCN28 = -JINV21+JINV22+JINV23
     
          JCN31 = -JINV31-JINV32-JINV33 
          JCN32 =  JINV31-JINV32-JINV33
          JCN33 =  JINV31+JINV32-JINV33          
          JCN34 = -JINV31+JINV32-JINV33
          JCN35 = -JINV31-JINV32+JINV33
          JCN36 =  JINV31-JINV32+JINV33
          JCN37 =  JINV31+JINV32+JINV33
          JCN38 = -JINV31+JINV32+JINV33

!         **************************************************************************
!	      Media dos valores nodais dos produtos das coord naturais pela
!         inversa de J(0) e pelas componentes de velocidade.
!         **************************************************************************
      
          pp(1,IELEM)=(JCN11*uvw(1,n1) + JCN12*uvw(1,n2) + JCN13*uvw(1,n3) +&
                       JCN14*uvw(1,n4)+  JCN15*uvw(1,n5) + JCN16*uvw(1,n6) +&
                       JCN17*uvw(1,n7) + JCN18*uvw(1,n8))*0.125D0
     
          pp(4,IELEM)=(JCN21*uvw(1,n1) + JCN22*uvw(1,n2) + JCN23*uvw(1,n3) +&
                       JCN24*uvw(1,n4)+  JCN25*uvw(1,n5) + JCN26*uvw(1,n6) +&
                       JCN27*uvw(1,n7) + JCN28*uvw(1,n8))*0.125D0
     
          pp(7,IELEM)=(JCN31*uvw(1,n1) + JCN32*uvw(1,n2) + JCN33*uvw(1,n3) +&
                       JCN34*uvw(1,n4)+  JCN35*uvw(1,n5) + JCN36*uvw(1,n6) +&
                       JCN37*uvw(1,n7) + JCN38*uvw(1,n8))*0.125D0


          pp(2,IELEM)=(JCN11*uvw(2,n1) + JCN12*uvw(2,n2) + JCN13*uvw(2,n3) +&
                       JCN14*uvw(2,n4) + JCN15*uvw(2,n5) + JCN16*uvw(2,n6) +&
                       JCN17*uvw(2,n7) + JCN18*uvw(2,n8))*0.125D0
     
          pp(5,IELEM)=(JCN21*uvw(2,n1) + JCN22*uvw(2,n2) + JCN23*uvw(2,n3) +&
                       JCN24*uvw(2,n4) + JCN25*uvw(2,n5) + JCN26*uvw(2,n6) +&
                       JCN27*uvw(2,n7) + JCN28*uvw(2,n8))*0.125D0
     
          pp(8,IELEM)=(JCN31*uvw(2,n1) + JCN32*uvw(2,n2) + JCN33*uvw(2,n3) +&
                       JCN34*uvw(2,n4) + JCN35*uvw(2,n5) + JCN36*uvw(2,n6) +&
                       JCN37*uvw(2,n7) + JCN38*uvw(2,n8))*0.125D0

          pp(3,IELEM)=(JCN11*uvw(3,n1) + JCN12*uvw(3,n2) + JCN13*uvw(3,n3) +&
                       JCN14*uvw(3,n4) + JCN15*uvw(3,n5) + JCN16*uvw(3,n6) +&
                       JCN17*uvw(3,n7) + JCN18*uvw(3,n8))*0.125D0
     
          pp(6,IELEM)=(JCN21*uvw(3,n1) + JCN22*uvw(3,n2) + JCN23*uvw(3,n3) +&
                       JCN24*uvw(3,n4) + JCN25*uvw(3,n5) + JCN26*uvw(3,n6) +&
                       JCN27*uvw(3,n7) + JCN28*uvw(3,n8))*0.125D0
     
          pp(9,IELEM)=(JCN31*uvw(3,n1) + JCN32*uvw(3,n2) + JCN33*uvw(3,n3) +&
                       JCN34*uvw(3,n4) + JCN35*uvw(3,n5) + JCN36*uvw(3,n6) +&
                       JCN37*uvw(3,n7) + JCN38*uvw(3,n8))*0.125D0

!         divergência da velocidade (dVk/dXk):
       
          PNV(IELEM) = pp(1,IELEM) + pp(5,IELEM) + pp(9,IELEM) 

!         ***************************************************************************
!	      Componentes do tensor de tensões viscosas do elemento IELEM tomados no
!         centro do memsmo. Considera-se que estas tensões são as tensões atuantes
!         na face de contorno!!
!         ***************************************************************************

!         Constantes de difusao

!          IF(IndTurb.NE.0) THEN
            AUXv = ViscCin + VDTurb(IELEM)   !viscosidade total
            IF(AUXv.LT.0.0D0) THEN
              AUXv = ViscCin
            end if
            COEF1= (2.0D0 * RHOInf * AUXv) + ViscVol
            COEF2= AUXv * RHOInf
            COEF3= ViscVol
!          end if

          TAO(1,J)= 2.0D0*COEF2*pp(1,IELEM) + COEF3*PNV(IELEM)
          TAO(2,J)= 2.0D0*COEF2*pp(5,IELEM) + COEF3*PNV(IELEM)
          TAO(3,J)= 2.0D0*COEF2*pp(9,IELEM) + COEF3*PNV(IELEM)
          TAO(4,J)= COEF2*( pp(2,IELEM) + pp(4,IELEM) )
          TAO(5,J)= COEF2*( pp(7,IELEM) + pp(3,IELEM) )
          TAO(6,J)= COEF2*( pp(8,IELEM) + pp(6,IELEM) )

        END DO

!       Processo de suavização das tensões viscosas:

        DO J=0,NELEM
          T11ARG(J)=0.0D0
          T22ARG(J)=0.0D0
          T33ARG(J)=0.0D0
          T12ARG(J)=0.0D0
          T13ARG(J)=0.0D0
          T23ARG(J)=0.0D0
        end do


        DO J=1,NFCS
          IF(J.GT.NFCN(NOB)) EXIT
          IELEM=IELCS(NOB,J)
          T11ARG(IELEM)= (TAO(1,J)*VOLU(IELEM))/8.0D0
          T22ARG(IELEM)= (TAO(2,J)*VOLU(IELEM))/8.0D0
          T33ARG(IELEM)= (TAO(3,J)*VOLU(IELEM))/8.0D0
          T12ARG(IELEM)= (TAO(4,J)*VOLU(IELEM))/8.0D0
          T13ARG(IELEM)= (TAO(5,J)*VOLU(IELEM))/8.0D0
          T23ARG(IELEM)= (TAO(6,J)*VOLU(IELEM))/8.0D0
        end do

!
!       Obtenção das componentes de tensões viscosas suavizadas nos nós de 
!       contorno sólido;
!       Obtenção do vetor de tensões totais (tensões viscosas + pressão
!       hidrostática) nos nós de contorno sólido;
!       Obtenção do vetor de tensões viscosas nos nós de contorno sólido. 

        DO I=1,NNCS

          IF(I.GT.NNCN(NOB)) EXIT
          J=IBNOBJ(NOB,I)

          K10 = NEIBOR(1,J)
          K20 = NEIBOR(2,J)
          K30 = NEIBOR(3,J)
          K40 = NEIBOR(4,J)
          K50 = NEIBOR(5,J)
          K60 = NEIBOR(6,J)
          K70 = NEIBOR(7,J)
          K80 = NEIBOR(8,J)
          K90 = NEIBOR(9,J)
          K100 = NEIBOR(10,J)
          K110 = NEIBOR(11,J)
          K120 = NEIBOR(12,J)
          K130 = NEIBOR(13,J)
          K140 = NEIBOR(14,J)
          K150 = NEIBOR(15,J)
          K160 = NEIBOR(16,J)

          T11S = ( T11ARG(K10)+T11ARG(K20)+T11ARG(K30)+T11ARG(K40)+T11ARG(K50)+T11ARG(K60)+ &
                   T11ARG(K70)+T11ARG(K80)+T11ARG(K90)+T11ARG(K100)+T11ARG(K110)+T11ARG(K120)+ &
                   T11ARG(K130)+T11ARG(K140)+T11ARG(K150)+T11ARG(K160) ) / MLUM(J)

          T22S = ( T22ARG(K10)+T22ARG(K20)+T22ARG(K30)+T22ARG(K40)+T22ARG(K50)+T22ARG(K60)+ &
                   T22ARG(K70)+T22ARG(K80)+T22ARG(K90)+T22ARG(K100)+T22ARG(K110)+T22ARG(K120)+&
                   T22ARG(K130)+T22ARG(K140)+T22ARG(K150)+T22ARG(K160) ) / MLUM(J)

          T33S = ( T33ARG(K10)+T33ARG(K20)+T33ARG(K30)+T33ARG(K40)+T33ARG(K50)+T33ARG(K60)+ &
                   T33ARG(K70)+T33ARG(K80)+T33ARG(K90)+T33ARG(K100)+T33ARG(K110)+T33ARG(K120)+ &
                   T33ARG(K130)+T33ARG(K140)+T33ARG(K150)+T33ARG(K160) ) / MLUM(J)

          T12S = ( T12ARG(K10)+T12ARG(K20)+T12ARG(K30)+T12ARG(K40)+T12ARG(K50)+T12ARG(K60)+ &
                   T12ARG(K70)+T12ARG(K80)+T12ARG(K90)+T12ARG(K100)+T12ARG(K110)+T12ARG(K120)+ &
                   T12ARG(K130)+T12ARG(K140)+T12ARG(K150)+T12ARG(K160) ) / MLUM(J)

          T13S = ( T13ARG(K10)+T13ARG(K20)+T13ARG(K30)+T13ARG(K40)+T13ARG(K50)+T13ARG(K60)+ &
                   T13ARG(K70)+T13ARG(K80)+T13ARG(K90)+T13ARG(K100)+T13ARG(K110)+T13ARG(K120)+ &
                   T13ARG(K130)+T13ARG(K140)+T13ARG(K150)+T13ARG(K160) ) / MLUM(J)

          T23S = ( T23ARG(K10)+T23ARG(K20)+T23ARG(K30)+T23ARG(K40)+T23ARG(K50)+T23ARG(K60)+ &
                   T23ARG(K70)+T23ARG(K80)+T23ARG(K90)+T23ARG(K100)+T23ARG(K110)+T23ARG(K120)+ &
                   T23ARG(K130)+T23ARG(K140)+T23ARG(K150)+T23ARG(K160) ) / MLUM(J)
  

          VTV1(I)= -T11S*CSNOR1(NOB,I) - T12S*CSNOR2(NOB,I) - T13S*CSNOR3(NOB,I)
          VTV2(I)= -T12S*CSNOR1(NOB,I) - T22S*CSNOR2(NOB,I) - T23S*CSNOR3(NOB,I)
          VTV3(I)= -T13S*CSNOR1(NOB,I) - T23S*CSNOR2(NOB,I) - T33S*CSNOR3(NOB,I)

          VTT1(I)= -(T11S - uvw(4,J))*CSNOR1(NOB,I) - T12S*CSNOR2(NOB,I) - T13S*CSNOR3(NOB,I)
          VTT2(I)= -T12S*CSNOR1(NOB,I) -  (T22S - uvw(4,J))*CSNOR2(NOB,I) - T23S*CSNOR3(NOB,I)
          VTT3(I)= -T13S*CSNOR1(NOB,I) - T23S*CSNOR2(NOB,I)- (T33S - uvw(4,J))*CSNOR3(NOB,I)

        END DO

!       Obtenção dos coeficientes aerodinâmicos:

        DO I=1,NNCS

          IF(I.GT.NNCN(NOB)) EXIT
          J=IBNOBJ(NOB,I)
 
          DeltX = xyz(1,J) - Xobj(NOB)
          DeltY = xyz(2,J) - Yobj(NOB)
          DeltZ = xyz(3,J) - Zobj(NOB)
 
          F1 = VTT1(I)*ARSuav(NOB,I)
          F2 = VTT2(I)*ARSuav(NOB,I)
          F3 = VTT3(I)*ARSuav(NOB,I)
          
          Fx = Fx + F1
          Fy = Fy + F2 
          Fz = Fz + F3
 
!         Momento aerodinâmico: sinal positivo pela regra da mão direita no sentido
!         dos eixos coordenados

          Mx = Mx + (F3*DeltY - F2*DeltZ)
          My = My + (F1*DeltZ - F3*DeltX)
          Mz = Mz + (F2*DeltX - F1*DeltY)

        end do
 
        CD = Fx	/ (Pdin*Dchar(NOB)*Lchar(NOB))
        CL = Fy	/ (Pdin*Dchar(NOB)*Lchar(NOB))
        CZ = Fz / (Pdin*Dchar(NOB)*Lchar(NOB))
        CMx = Mx	/ (Pdin*(Dchar(NOB)**2)*Lchar(NOB))
        CMy = My	/ (Pdin*(Dchar(NOB)**2)*Lchar(NOB))
        CMz = Mz	/ (Pdin*(Dchar(NOB)**2)*Lchar(NOB))

!       Obtenção da componente tangencial do tensor de tensões viscosas que atua 
!       sobre a superfície de contorno sólido. 

        DO I=1,NNCS

          IF(I.GT.NNCN(NOB)) EXIT
          
          CSN21= VTV2(I)*CSNOR3(NOB,I) - VTV3(I)*CSNOR2(NOB,I)     
          CSN22= VTV3(I)*CSNOR1(NOB,I) - VTV1(I)*CSNOR3(NOB,I)     
          CSN23= VTV1(I)*CSNOR2(NOB,I) - VTV2(I)*CSNOR1(NOB,I)     
          CSN31= CSNOR2(NOB,I)*CSN23 - CSNOR3(NOB,I)*CSN22     
          CSN32= CSNOR3(NOB,I)*CSN21 - CSNOR1(NOB,I)*CSN23     
          CSN33= CSNOR1(NOB,I)*CSN22 - CSNOR2(NOB,I)*CSN21 

!         NÃO É VTVi(I) e sim CSNORi(I)????!!!!

          T_CSN3= VTV1(I)*CSN31 + VTV2(I)*CSN32 + VTV3(I)*CSN33    
          CSN3_CSN3= CSN31*CSN31 + CSN32*CSN32 + CSN33*CSN33

          FFACT= 0.0D0
          IF (DABS(CSN3_CSN3).GE.1.0D-7)  FFACT= T_CSN3 / CSN3_CSN3

          TTT1(I)= FFACT * CSN31
          TTT2(I)= FFACT * CSN32
          TTT3(I)= FFACT * CSN33

        END DO

!       Obtenção dos coeficientes de pressão e de fricção ("skin friction")  

        DO I=1,NNCS

          IF(I.GT.NNCN(NOB)) EXIT
          NOCOSO= IBNOBJ(NOB,I)

          CP(I)=(uvw(4,NOCOSO)-PINF)/(Pdin)
          CPsuav(I)=(PRSuav(NOCOSO)-PINF)/(Pdin)

          CF1(I)=TTT1(I)/(Pdin)  
          CF2(I)=TTT2(I)/(Pdin)
          CF3(I)=TTT3(I)/(Pdin)

        END DO


        NARQ=NARQ+1
        WRITE(NARQ,*) CD

        NARQ=NARQ+1
        WRITE(NARQ,*) CL

        NARQ=NARQ+1
        WRITE(NARQ,*) CZ

        NARQ=NARQ+1
        WRITE(NARQ,*) CMx

        NARQ=NARQ+1
        WRITE(NARQ,*) CMy

        NARQ=NARQ+1
        WRITE(NARQ,*) CMz

!       Coeficiente de pressão: Cp(i) = P(i) - Pref / Pdin

        NARQ=NARQ+1
        WRITE(NARQ,*) TEMPO

        DO J=1,NNCS
          IF(J.GT.NNCN(NOB)) EXIT
          L=IBNOBJ(NOB,J)
          WRITE(NARQ,*) L,CP(J),CPsuav(J)
        END DO

        WRITE(NARQ,*) ' '

!       Coeficiente "skin friction": CF(i) = TVisc(i) / Pdin

        NARQ=NARQ+1
        WRITE(NARQ,*) TEMPO

        DO J=1,NNCS
          IF(J.GT.NNCN(NOB)) EXIT
          L=IBNOBJ(NOB,J)
          WRITE(NARQ,*) L,CF1(J),CF2(J),CF3(J)
        END DO

        WRITE(NARQ,*) ' '

      end do

    end if

  END IF									

! ***********************************************************************
! FIM DA CONDICAO PARA CALCULO DE COEFICENTES E GRAVACAO DE RESULTADOS
! ***********************************************************************


! CALCULO DA CONVERGENCIA TEMPORAL E TOMADA DE DECISIAO SOBRE A CONTINUIDADE
! **************************************************************************
! OU FIM DO PROCESSAMENTO
! ************************

  SDETP=0.0D0
  SDTP=0.0D0
      
  IF (NNNRRR.GT.1) THEN

    DO I=1,NNM
      RES1=PRf(I)-PRi(I)
      RESCUA1=RES1*RES1
      SDETP=SDETP+RESCUA1
   
      PRES=PRi(I)
      PRESCUA=PRES*PRES
      SDTP=SDTP+PRESCUA
    END DO

    SDETP=DSQRT(SDETP/SDTP)

!   GRAVACAO DO ARQUIVO DE RESIDUOS
!   ---------------------------------
                    
    ICANT=NIR/4
    RNNN=FLOAT(NNNRRR)
    CANT=FLOAT(ICANT)
    RESS=RNNN/CANT
    IRESS=AINT(RESS)
    RESSE=FLOAT(IRESS)
    VARIND=RESS-RESSE

!   float - converte um número inteiro em um real
!   aint - trunca o argumento

    IF (VARIND.EQ.0.0D0) WRITE(14,202) TEMPO,NNNRRR,SDETP

    INDTPO=0

    IF (SDETP.LE.TOLTPO) INDTPO=INDTPO+1

  END IF


!     ****************************************************
!     Final do ciclo de tempo
!     ****************************************************

end do

!     ****************************************************
!     Final do ciclo de tempo
!     ****************************************************

CLOSE (14)
IF(MNOBJ.NE.0) THEN
  DO I=30,NOPF
    CLOSE(I)
  end do
end if

! **************************************************************************
! REGISTRO DAS VARIAVIES APOS O FINAL DO CICLO DE TEMPO
! **************************************************************************

NFILE=NFILE+1
TPOAC=TEMPO

! ARQUIVO COM TEMPO E PROXIMO NFILE
! ---------------------------------

PTNAME(7:14)='.steady '
OPEN(2,FILE=PTNAME,STATUS='UNKNOWN')
REWIND (2)
 
NF1=NFILE/100+48
NF2=(NFILE-(NFILE/100)*100)/10+48
NF3=NFILE-(NFILE/10)*10+48
PTNAME(7:9)=CHAR(NF1)//CHAR(NF2)//CHAR(NF3)

PTNAME(10:14)='.v   '
OPEN(3,FILE=PTNAME,STATUS='UNKNOWN')
REWIND (3)

PTNAME(10:14)='.pr  '
OPEN(4,FILE=PTNAME,STATUS='UNKNOWN')
REWIND (4)

PTNAME(10:14)='.vts '
OPEN(5,FILE=PTNAME,STATUS='UNKNOWN')
REWIND (5)

WRITE(2,909)NFILE,TEMPO

DO I=1,NNM
  WRITE(3,*) uvw(1,I),uvw(2,I),uvw(3,I)
  WRITE(4,*) uvw(4,I),PRSUAV(I)
  WRITE(5,*) VtSuav(I)
END DO

CLOSE(2)
CLOSE(3)
CLOSE(4)
CLOSE(5)


! **************************************************************************
! REGISTRO DAS VARIAVEIS DE CAMPO MÉDIO
! *************************************************************************

IF(NPASS.NE.NNAUX) THEN

  PTNAME(7:14)='.vm     '
  OPEN(3,FILE=PTNAME,STATUS='UNKNOWN')
  REWIND (3)

  PTNAME(7:14)='.prm    '
  OPEN(4,FILE=PTNAME,STATUS='UNKNOWN')
  REWIND (4)

  PTNAME(7:14)='.vtsm   '
  OPEN(5,FILE=PTNAME,STATUS='UNKNOWN')
  REWIND (5)

  DO I=1,NNM
    WRITE(3,*) Umed(I),Vmed(I),Wmed(I)
    WRITE(4,*) PRmed(I),PRSmed(I)
    WRITE(5,*) VTmed(I)
  END DO

  CLOSE(3)
  CLOSE(4)
  CLOSE(5)
  
end if

elapsed_time = TIMEF()
WRITE(*,*)
WRITE(*,*) '**********************************************'
WRITE(*,*) 'Tempo total de execucao:',elapsed_time
WRITE(*,*) '**********************************************'

! **************************************************************************
! FORMATS
! *********************************

  101 FORMAT (3(D19.12))
  201 FORMAT (2(D19.12))
  301 FORMAT (D19.12)
  202 FORMAT (E15.7,I10,E15.7,E15.7)
  303 FORMAT (I8)
  404 FORMAT (3(E12.5))
  505 FORMAT (8(I8))
  606 FORMAT (6(E12.5))
  707 FORMAT (3(D13.5,1X))
  808 FORMAT (10(E12.5))
  909 FORMAT (I3,E12.5)
  910 FORMAT ('Passo de Tempo:',I8,' Tempo: ',D10.3)
      

END
