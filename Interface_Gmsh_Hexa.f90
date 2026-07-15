module functions
    contains
!***********************************************************
!			FUNÇÕES
!***********************************************************
    
	function PowerLaw(z, Vref, zd, zref, p) result(V)

		implicit none
    
		real(8), intent(in) :: z, Vref, zd, zref, p
		real(8) :: V
    
		if (z .lt. zd) then
			V = 0.0d0
		else
			V = Vref * ( (z - zd) / zref )**p
		end if
    
	end function PowerLaw
	
	function LogLaw(z, Ustar, zd, z0) result(V)

		implicit none
    
		real(8), intent(in) :: z, Ustar, zd, z0
		real(8) :: V
    
		if ((z-zd) .lt. z0) then
			V = 0.0d0
		else
			V = (Ustar / 0.4d0) * log( (z - zd) / z0 )
		end if
    
	end function LogLaw

end module
    
program Interface_Gmsh_Hexa

    !**********************************************************************
    ! Programa de interface entre os dados da malha gerados pelo GMSH e o
    ! programa de CFD em MEF com elementos hexaedricos de 8 nós denominado
    ! fluexp2p_Hermann_nuevo.f90 ou fluexp2p_dissert_guilherme.f90.
    !
    ! Autor: Ígor Marini Peter - PPGEC - UFRGS
    ! Baseado no código de interface para elementos tetrahedricos última
    ! vez modificado por Gabriela Bianchin - PPGEC - UFRGS
    !
    ! Última modificacao: 13/04/2023
    !**********************************************************************
    ! Obs.: ______ refere-se a um nome de 6 dígitos dado ao projeto e aparece na maioria dos arquivos
    !
    ! Arquivos de Entrada:
    !   101 - ______.phys -> Inclui os dados das condições de contorno através dos grupos físicos
    !   102 - ______.msh -> Arquivo de exportacao da malha do GMSH em formato .msh versao 2.2
    !
    ! Arquivos de Saída:
    !   201 - inicia.par        -> Parametros de inicializacao da simulacao numerica
    !   202 - ______.pro        -> Propriedades do fluido, do fluxo e dimensões características dos corpos imersos
    !   203 - ______.con        -> Conectividade dos elementos
    !   204 - ______000.cor  	-> Coordenadas iniciais dos nós da malha
    !   205 - ______FF.sup      -> Areas e conectividades das faces de contorno com forcas de suprfície (nao estritamente necessario)
    !   206 - ______FF.nnn      -> Vetores normais aos nós de contorno com condicões de Neumann
    !               ** Arquivos FF sao necessarios devido às condicões de contorno da integracao por partes
    !               ** podendo ser desprezados se o final do domínimo de calculo for longe o suficiente
    !   301 - ______CS##.sup    -> Áreas e conectividades das faces de contorno sólido
    !   302 - ______CS##.nnn    -> Vetores normais aos nós de controno sólido
    !               ** Os arquivos CS sao necessarios quando ha um ou mais corpos submersos, sendo utilizados
    !               ** para calcular a pressao no corpo devido ao escoamento. Deve-se gerar arquivos tanto quanto
    !               ** necessarios para representar todos os corpos submersos, com o número ## caracterizando o
    !               ** número do corpo submerso
    !   207 - ______000.v       -> Campo de velocidades inicial do fluido
    !   208 - ______000.pr      -> Campo de pressões inicial do fluido
    !   209 - ______000.t       -> Campo de temperatura inicial do fluido
    !   210 - ______000.rha     -> Campo de densidade de especia no fluido
    !   211 - ______CC.bv  		-> Condicões de contorno de velocidade do escoamento
    !   212 - ______CC.bp  		-> Condicões de contorno de pressao do escoamento
    !   213 - ______CC.bt       -> Condicões de contorno de temperatura do escoamento
    !   214 - ______CC.brha     -> Condicões de controno de especie do escoamento
    !   215 - ______CE.fan      -> Condicões especiais ???
    !   216 - ______CE.faq      -> Condicões especiais ???
    !   217 - ______CE.raq      -> Condicões especiais ???
    !**********************************************************************
	use functions
    
    !**********************************************************************
    ! FORMATOS
    !**********************************************************************
    
    !*** CONSOLE ***!
1   FORMAT (A,T80)		! formato para escrever no cmd
    
    !*** INICIA.PAR E PNAME.PRO ***!
2   FORMAT (A,T10,A,T6)		! formato para nome do projeto OpenMP
3   FORMAT (A,T10,I19)		! formato para inteiro do projeto OpenMP
4   FORMAT (A,T10,D20.12)	! formato para real do projeto OpenMP
5   FORMAT (A,T10,3D20.12)  ! formato para três reais do projeto OpenMP
6   FORMAT (A,T6)			! formato para nome do projeto CUDA
7   FORMAT (I19)			! formato para inteiro do projeto CUDA
8   FORMAT (D20.12)			! formato para real do projeto CUDA
9   FORMAT (3D20.12)		! formato para três reais do projeto CUDA
10  FORMAT (2I19)			! formato para dois inteiros do projeto CUDA
    
    !*** ARQUIVOS DE SAÍDA ***!
11  FORMAT (I10)			! formato para inteiro, usado para armazenar contadores
12  FORMAT (3D20.12)		! formato para 3 números decimais, usado nas coordenadas dos nós e no campo de velocidades
13  FORMAT (2D20.12)		! formato para 2 números decimais, usado nos campos de pressão, temperatura e espécie
14  FORMAT (8I10)			! formato para armazenamento das conectividades dos elementos, 8 números inteiros com até 9 dígitos
15  FORMAT (I10, D20.12)	! formato para inteiro e real, usado nos arquivos CC
16  FORMAT (D20.12, 4I10)	! formato para real e 4 inteiros, usado nos arquivos CS
17	FORMAT (I10, 3D20.12)	! formato para inteiro e 3 reais, usado nos arquivos CS    
18  FORMAT (D20.12)			! formato para 1 número decimal, usado para arquivos FF, CC e CE quando zerados
    
    !*** PROGRESSO ***!
100 FORMAT (I10, ' de ', I10, ' quads | ' I10, ' de ', I10, ' hexas | ', F10.2, '%')

	!*** NEIBOR ***!
200 FORMAT ('NEIBOR - ', I10)
    
    !**********************************************************************
    ! DECLARAÇÃO DAS VARIÁVEIS DO PROGRAMA
    !**********************************************************************
    implicit none
    
    logical :: bool, T, RHA           		! boolean
    
    character (6) :: pname                  ! nome do projeto
    character (14) :: ptname                ! nome do arquivo
    character (200) :: header               ! header para leitura de arquivos
    character (200) :: garbc                ! caracter desnecessario
    
	integer :: CodeType, INDTURB
    integer :: i, j, k1, k2, n, n1, n2, n3, n4, n5, n6, n7, n8, face, group, counter, iele   ! contadores e auxiliares inteiros
	integer :: count1, count2, count3, count4
    integer :: ntn, nquad, nhexa, ntele     ! número de nós, quadrilateros, hexahedros e elementos
    integer :: nV1potnodes, nV3potnodes		! número de nós com CC por lei de potência
	integer :: nV1lognodes, nV3lognodes		! número de nós com CC por lei logaritmica
	integer :: nV1nodes, nV2nodes, nV3nodes, nPnodes	! número de nós com condicões de contorno
	integer :: MNOBJ						! número de objetos submersos no escoamento
	integer :: MAXcountfIFE					! número máximo de faces em um objeto submerso
    integer :: garbi                        ! inteiro desnecessario
	integer :: NB1, NB2
    

	integer :: quadnodes(4)					! vetor auxiliar com os nós de um elementos quadrilatero
	
    integer, allocatable :: nccV1potnodes(:), nccV3potnodes(:)	! vetores com os nós com CC por lei de potência
	integer, allocatable :: nccV1lognodes(:), nccV3lognodes(:)	! vetores com os nós com CC por lei logaritmica
	integer, allocatable :: nccV1nodes(:), nccV2nodes(:), nccV3nodes(:), nccPnodes(:)			! vetores com os nós com condicao de contorno
	integer, allocatable :: CSsupELEM(:,:)		! matriz com o número dos elementos a serem armazenados no arquivo CS##.sup
												! onde CSsupELEM(2,3) é o terceiro elemento do CS02.sup
	integer, allocatable :: CSsupKONE(:,:,:)	! matriz com a conectividade dos elementos quadriláteros a serem armazenados no arquivo CS##.sup
												! onde CSsupKONE(2,3,4) é o número do quarto nó, do terceiro elemento do CS02.sup
	integer, allocatable :: CSnnnNODE(:,:)		! matriz com o número dos nós a serem armazenados no arquivo CS##.nnn
												! onde CSnnnNODE(2,3) é o terceiro nó do CS02.nnn
    integer, allocatable :: CSsupCOUNT(:)		! vetor que armazena o número de elementos quadriláteros em cada corpo submerso
	integer, allocatable :: NEIBOR(:) 	! Verificação de quantos elementos tem em volta de cada nó
	
	real(8) :: garbr                           ! real desnecessario
	real(8) :: r, r2, r3, NNN1, NNN2, NNN3, ccvalue, V1, V2, V3, P
    
    real(8), allocatable :: vccV1potnodes(:), vccV3potnodes(:)		! vetores com os valores da CC por lei de potência em cada nó
	real(8), allocatable :: vccV1lognodes(:), vccV3lognodes(:)		! vetores com os valores da CC por lei logaritmica em cada nó
	real(8), allocatable :: vccV1nodes(:), vccV2nodes(:), vccV3nodes(:), vccPnodes(:)	! vetores com os valores da condicao de contorno em cada nó
    real(8), allocatable :: xx(:), yy(:), zz(:)    ! vetores de coordenadas dos nós
	real(8), allocatable :: CSsupAREA(:,:)		! matriz com a área dos elementos quadriláteros a serem armazenados no arquivo CS##.sup
											! onde CSsupAREA(2,3) é a área do terceiro elemento do CS02.sup
	real(8), allocatable :: CSnnnNORM(:,:,:)	! matriz com as componentes normais nos nós do objeto submerso a serem armazenados no arquivo CS##.nnn
											! onde CSnnn(1,2,3) é o terceiro componente do segundo nó do CS01.nnn		
    
	! Variaveis da leitura dos grupos físicos e das CC
    integer 				:: countPhysGroups	! número de grupos físicos
    integer 				:: countCCTypes(4) 	! vetor com o número de CC nulas, constante, potência e log
	integer					:: initialFieldGroup	! número do grupo com as CC que representam o campo inicial de velocidade e pressão
	integer, allocatable	:: indexCC(:,:)			! indexagem das CC
	integer, allocatable	:: indexNullCC(:,:)		! indexagem das CC nulas
	integer, allocatable	:: indexConstCC(:,:)	! indexagem das CC constantes
	integer, allocatable	:: indexPowerCC(:,:)	! indexagem das CC por Lei de Potência
	integer, allocatable	:: indexLogCC(:,:)		! indexagem das CC por Lei Logaritmica
													! os quatro são número do grupo seguido por número da variavel (0=V1, 1=V2, 2=V3, 3=P)
	integer, allocatable	:: OBJPhysGroup(:)		! vetor que armazena o núemro do grupo físico de cada objeto submerso
	real(8), allocatable	:: valueConstCC(:) 		! valores das CC constantes
	real(8), allocatable	:: valuePowerCC(:,:)	! valores das CC por Lei de Potência (Vref, zd, zref, p)
	real(8), allocatable	:: valueLogCC(:,:)		! valroes das CC por Lei Logaritmica (u*, zd, z0)
	
	! Variaveis dos objetos submersos
	real(8), allocatable	:: objCenter(:,:)		! coordenadas x, y e z do centro de cada objeto, usado para definir a direção normal
	
												
	!**********************************************************************
    ! DADOS INICIAIS
    !**********************************************************************
    write (*, 1) 'INFORMATIVO SOBRE O PROJETO' 
    write (*, 1) 'Digite o nome do projeto (pname) - em seis caracteres:'
    read (*, *) pname
	ptname (1:6) = pname
    
    !write (*, 1) 'Qual e o perfil da condicao de contorno da entrada?'
    !write (*, 1) 'Uniforme = digite 0'
    !write (*, 1) 'Lei Potencial = digite 1'
    !write (*, 1) 'Lei Logaritmica = digite 2'
    !read(*,*) inle3
    
    !write (*, *) 'Qual e a condicao de contorno da saida?'
    !write (*, *) 'Pressao nula = digite 0'
    !write (*, *) 'Gradiente nulo (Outflow) = digite 1'
    !write (*, *) 'Velocidade uniforme = digite 2'
    !write (*, *) 'Velocidade paralela e pressao nulas = digite 3'
    !write (*, *) 'pressao nula nos pontos do extremo dominio inferior = digite 4'
    !read(*,*) outl
    
    ! TODO - O original implementa ainda questao de IFE, Temperatura, Especie e Cascas
    
    write(*,1) 'Informe a quantidade de elementos quadrilateros do problema:'
    read(*,*) nquad
    write(*,1) 'Informe a quantidade de elementos hexaedricos do problema:'
    read(*,*) nhexa
    write(*,1) 'Informe o tipo de algoritmo que sera utilizado:'
    write(*,1) '	OpenMP	= 0'
    write(*,1) '	CUDA 	= 1'
	write(*,1) '	CUDA_HEXA_IGOR	= 2'
    read(*,*) CodeType

    write(*,1) 'O problema envolve campo de temperaturas?'
    write(*,1) '	Nao = 0'
    write(*,1) '	Sim = 1'
    read(*,*) n1
    if (n1.eq.0) then
        T = .FALSE.
    else
        T = .TRUE.
    end if
    write(*,1) 'O problema envolve transporte de particulas?'
    write(*,1) '	Nao = 0'
    write(*,1) '	Sim = 1'
    read(*,*) n1
    if (n1.eq.0) then
        write(*,1) '[DEBUG] ZERO transporte'
        RHA = .FALSE.
    else
        RHA = .TRUE.
    end if
    
    
    !*********************************************************************
    ! DADOS DE ENTRADA
    !*********************************************************************
    
    ! TODO - Seria melhor se os dados de entrada fossem identificado atraves de headers ao inves de utilizar garbc
    !open (101, file="dadosV2.dat", status='old')
    !read (101,*) garbc
    !read (101,*) garbc
    !read (101,*) nV1pot, nV3pot   ! número de faces com condições de controno V1 e V3 por lei de potência
    
    !allocate( nccV1pot(nV1pot), vccV1pot(nV1pot,4), nccV3pot(nV3pot), vccV3pot(nV3pot,4) )
    
    !read (101,*) garbc
    !do i = 1,nV1pot
    !    read (101,*) nccV1pot(i), vccV1pot(i,1), vccV1pot(i,2), vccV1pot(i,3), vccV1pot(i,4)
    !end do
    
    !read (101,*) garbc
    !do i = 1,nV3pot
    !    read (101,*) nccV3pot(i), vccV3pot(i,1), vccV3pot(i,2), vccV3pot(i,3), vccV3pot(i,4)
    !end do
    
    !read (101,*) garbc
    !read (101,*) nV1, nV2, nV3, nP   ! número de faces com condicões de contorno V1, V2, V3 e P
    
    !allocate( nccV1(nV1), vccV1(nV1), nccV2(nV2), vccV2(nV2), nccV3(nV3), vccV3(nV3), nccP(nP), vccP(nP) )
    
    ! Ler o número das faces correspondentes a cada tipo
    ! O número da face pode ser obtida atraves da geometria pelo GMSH
    !read (101,*) garbc
    !do i = 1,nV1
    !    read (101,*) nccV1(i), vccV1(i)
    !end do
    !read (101,*) garbc
    !do i = 1,nV2
    !    read (101,*) nccV2(i), vccV2(i)
    !end do
    !read (101,*) garbc
    !do i = 1,nV3
    !    read (101,*) nccV3(i), vccV3(i)
    !end do
    !read (101,*) garbc
    !do i = 1,nP
    !    read (101,*) nccP(i), vccP(i)
    !end do
	
	! Ler os dados dos objetos submersos
    !read (101,*) garbc
	!read (101,*) MNOBJ
	
	!allocate( countfIFE(MNOBJ) )
	
	! Para cada objeto submerso, registrar o número de faces que compoem o mesmo
	!MAXcountfIFE = 0 ! contador para registrar o número de faces do objeto submerso com maior número de faces
	!read (101,*) garbc
	!do i=1,MNOBJ
	!	read (101,*) countfIFE(i)
	!	if (countfIFE(i) .gt. MAXcountfIFE) then
	!		MAXcountfIFE = countfIFE(i)
	!	end if
	!end do
	
	!allocate( nfIFE(MNOBJ,MAXcountfIFE) )
	!nfIFE = 0
	
	! Para cada objeto submerso, registrar o número das faces que compoem o mesmo
	!do i=1,MNOBJ
	!	read (101,*) garbc
	!	do j=1,countfIFE(i)
	!		read (101,*) nfIFE(i,j)
	!	end do
	!end do
	
    !close (101)
    
    !****************************************************************
    ! LEITURA DO ARQUIVO DE CONDIÇÕES DE CONTORNO .PHYS
    !****************************************************************
	
    ptname (7:14) = '.phys   '
    open (101, file=ptname, status='old')
    
    bool = .TRUE.
    do while ( bool )
        read (101,*) header
        
        select case (header)
        case ('$GroupCount')
            ! Armazenar número de grupos físicos e o número de cada tipo de CC
            read (101,*) countPhysGroups, countCCTypes(1), countCCTypes(2), countCCTypes(3), countCCTypes(4)
			allocate ( indexNullCC(countCCTypes(1),2), indexConstCC(countCCTypes(2),2), &
						indexPowerCC(countCCTypes(3),2), indexLogCC(countCCTypes(4),2), &
						valueConstCC(countCCTypes(2)), valuePowerCC(countCCTypes(3),4), &
						valueLogCC(countCCTypes(4),3), indexCC(CountPhysGroups,5) )
            read (101,*) garbc	! $EndGroupCount
        case ('$GroupType')
            ! Armazenar os índices das condições de contorno
			! Contadores
			count1 = 0
			count2 = 0
			count3 = 0
			count4 = 0
			do i=1,countPhysGroups
				read (101,*) n1, n2, n3, n4, n5
				
				indexCC(i,1) = n1
				indexCC(i,2) = n2
				indexCC(i,3) = n3
				indexCC(i,4) = n4
				indexCC(i,5) = n5
				
				! V1
				select case (n2)
				case (1)
					count1 = count1 + 1
					indexNullCC(count1, 1) = n1
					indexNullCC(count1, 2) = 0
				case (2)
					count2 = count2 + 1
					indexConstCC(count2, 1) = n1
					indexConstCC(count2, 2) = 0
				case (3)
					count3 = count3 + 1
					indexPowerCC(count3, 1) = n1
					indexPowerCC(count3, 2) = 0
				case (4)
					count4 = count4 + 1
					indexLogCC(count4, 1) = n1
					indexLogCC(count4, 2) = 0
				end select
				
				! V2
				select case (n3)
				case (1)
					count1 = count1 + 1
					indexNullCC(count1, 1) = n1
					indexNullCC(count1, 2) = 1
				case (2)
					count2 = count2 + 1
					indexConstCC(count2, 1) = n1
					indexConstCC(count2, 2) = 1
				case (3)
					count3 = count3 + 1
					indexPowerCC(count3, 1) = n1
					indexPowerCC(count3, 2) = 1
				case (4)
					count4 = count4 + 1
					indexLogCC(count4, 1) = n1
					indexLogCC(count4, 2) = 1
				end select
				
				! V3
				select case (n4)
				case (1)
					count1 = count1 + 1
					indexNullCC(count1, 1) = n1
					indexNullCC(count1, 2) = 2
				case (2)
					count2 = count2 + 1
					indexConstCC(count2, 1) = n1
					indexConstCC(count2, 2) = 2
				case (3)
					count3 = count3 + 1
					indexPowerCC(count3, 1) = n1
					indexPowerCC(count3, 2) = 2
				case (4)
					count4 = count4 + 1
					indexLogCC(count4, 1) = n1
					indexLogCC(count4, 2) = 2
				end select
				
				! P
				select case (n5)
				case (1)
					count1 = count1 + 1
					indexNullCC(count1, 1) = n1
					indexNullCC(count1, 2) = 3
				case (2)
					count2 = count2 + 1
					indexConstCC(count2, 1) = n1
					indexConstCC(count2, 2) = 3
				case (3)
					count3 = count3 + 1
					indexPowerCC(count3, 1) = n1
					indexPowerCC(count3, 2) = 3
				case (4)
					count4 = count4 + 1
					indexLogCC(count4, 1) = n1
					indexLogCC(count4, 2) = 3
				end select
				
			end do
			
			read (101,*) garbc	! $EndGroupType
        case ('$ConstantValues')
            ! Armazenar os valores das CC constantes
			do i=1,countCCTypes(2)
				read (101,*) valueConstCC(i)
			end do
			
			read (101,*) garbc	! $EndConstantValues
        case ('$PowerLawValues')
            ! Armazenar os valores das CC por Lei de Potência
			do i=1,countCCTypes(3)
				read (101,*) valuePowerCC(i,1), valuePowerCC(i,2), valuePowerCC(i,3), valuePowerCC(i,4)
			end do
			
			read (101,*) garbc	! $EndPowerLawValues
        case ('$LogLawValues')
            ! Armazenar os valores das CC por Lei Logaritmica
			do i=1,countCCTypes(4)
				read (101,*) valueLogCC(i,1), valueLogCC(i,2), valueLogCC(i,3)
			end do
			
			read (101,*) garbc	! $EndLogLawValues
		case ('$ObjectCount')
			! Armazenar o número de objetos submersos MNOBJ
			read (101,*) MNOBJ
			allocate ( OBJPhysGroup(MNOBJ) )
			read (101,*) garbc	! $EndObjectCount
		case ('$ObjectProp')
			! Armazenar o número dos grupos físicos de cada objeto submerso
			do i=1,MNOBJ
				read (101,*) OBJPhysGroup(i)
			end do
			
			read (101,*) garbc	! $EndObjectProp
		case ('$InitialFieldGroup')
			! Armazenar o número do grupo físico com as CC correspondentes ao campo inicial
			read (101,*) initialFieldGroup
			
			read (101,*) garbc	! $EndInitialFieldGroup
		case ('$End')
			bool = .FALSE.
        end select
    end do
    
	!****************************************************************
	! CRIAÇÃO DO ARQUIVO DE PARÂMETROS DE INICIALIZAÇÃO INICIA.PAR
	!****************************************************************
	
	write (*, 1) 'Gerar arquivo inicia.par?        0 = Nao   1 = Sim'
	read (*, *) i
	
	if (i .eq. 1) then
		open (201, file='inicia.par', status='unknown')
        if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (201, 6) pname
        else
			write (201, 2) 'PNAME: ', pname
        end if
		
		write (*, 1) 'Os dados de tempo a seguir sao utilizados para determinar &
		a duracao da simulacao e a frequencia da gravacao dos resultados em disco.'
		write (*, 1) 'Tenha em mente que o tempo total da simulacao e dado por:'
		write (*, 1) '     Ttot = NIR * NTR * DtMAX'
		write (*, 1) 'Entao escolha os valores de acordo com o tempo de simulacao desejado.'
		write (*, 1) ''
	
		write (*, 1) 'Indique o numero de intervalos entre registros de coeficientes aerodinamicos NCOEF:'
		write (*, 1) 'Esse valor e utilizado para estabelecer o periodo de &
		armazenamento dos coeficientes aerodinamicos para corpos imersos no &
		escoamento. E usual que os mesmo sejam armazenados com uma frequencia &
		maior que os dados dos campos de velocidade e pressao do escoamento'
		read (*, *) n1
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (201, 7) n1
        else
            write (201, 3) 'NCOEF:', n1
        end if
		write (*, 1) ''
		
		write (*, 1) 'Indique o numero de intervalos entre registros NIR:'
		write (*, 1) 'Este e o numero de passos de tempo decorridos para que &
		seja feito o registro dos campos de velocidade e pressao do escoamento.'
		read (*, *) n1
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (201, 7) n1
        else
            write (201, 3) 'NIR:', n1
        end if
		write (*, 1) ''
	
		write (*, 1) 'Indique o numero total de registros NTR:'
		write (*, 1) 'Registros mais frequentes permitem analisar resultados &
		preliminares antes do fim da simulacao, e permitem retomar analises no &
		caso de problemas, como quedas de energia, a partir do registro mais recente.'
		read (*, *) n1
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (201, 7) n1
        else
            write (201, 3) 'NTR:', n1
        end if
		write (*, 1) ''
	
		write (*, 1) 'Indique o numero do registro mais recente NFILE:'
		write (*, 1) 'Este valor deve ser utilizado quando se deseja continuar &
		uma simulacao anterior a partir de um registro, sendo NFILE o numero do &
		registro. Caso a simulacao seja nova, NFILE = 0.'
		read (*, *) n1
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (201, 7) n1
        else
            write (201, 3) 'NFILE:', n1
        end if
		write (*, 1) ''
	
		write (*, 1) 'Indique a estimativa inicial para o incremento de tempo DtMAX:'
		write (*, 1) 'Como regra geral, DtMAX = DeltaX / (c + Vinf), sendo &
		DeltaX o tamanho característico do menor elemento da malha, c a &
		velocidade do som no escoamento e Vinf a velocidade no infinito, ou &
		a velocidade de referencia do problema. Em certos casos, c pode ser &
		manipulado para obter valores maiores para DtMAX. Valores de DtMAX que &
		sejam superiores a estimativa acima podem levar a problemas de estabilidade.'
		read (*, *) r
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (201, 8) r
        else
            write (201, 4) 'DtMAX:', r
        end if
		write (*, 1) ''
		
		! TPOAC não é utilizado no Hermann_nuevo
		!write (*, *) 'Indique o tempo atual TPOAC:'
		!read (*, *) r
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (201, 8) 0.0
        else
            write (201, 4) 'TPOAC:', 0.0
        end if
		
		!write (*, 1) 'Indique a tolerancia para o residuo TOLTPO:'
		!write (*, 1) 'O calculo do residuo e uma forma de interromper a &
		!simulacao ao se identificar que o escoamento esta estacionario. &
		!Recomenda-se um residuo maximo de 10^-6.'
		!read (*, *) r
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (201, 8) 0.000001
        else
            write (201, 4) 'TOLTPO:', 0.000001
        end if
		!write (*, 1) ''
		
		!write (*, 1) 'Indique o numero de passos de tempo apos atingir a tolerancia NROTPO:'
		!write (*, 1) 'Esta variavel determina quantos passos de tempo percorrer &
		!para garantir que o escoamento seja de fato estacionario. Recomenda-se &
		!adotar um valor minimo de 100.'
		!read (*, *) n1
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (201, 7) 100
        else
            write (201, 3) 'NROTPO:', 100
        end if
		!write (*, 1) ''
		
		write (*, 1) 'Indique o coeficiente de seguranca para determinar o incremento de tempo CSEGUR:'
		write (*, 1) 'Este trata de coeficiente de seguranca que leva em consideracao &
		as diferentes direcoes que uma particula pode passar pelo elemento. &
		Para escoamentos laminares admite-se valores proximos a 0,7, mas para &
		escoamento turbulentos, recomenda-se valores entre 0,1 e 0,3.'
		read (*, *) r
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (201, 8) r
        else
            write (201, 4) 'CSEGUR:', r
        end if
		write (*, 1) ''
		
		!write (*, 1) 'Indique o parametro de controle de modos expurios CONTROL1:'
		!write (*, 1) 'Este parametro define se sera utilizada integracao &
		!reduzida para controle dos modos expurios. Deve ser sempre igual a 0 &
		!quando nao se deseja integracao reduzida e 1 quando se deseja integracao &
		!reduzida.'
		!read (*, *) r
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (201, 8) 1.0
        else
            write (201, 4) 'CONTROL1:', 1.0
        end if
		!write (*, 1) ''
		
		write (*, 1) 'Indique o numero de passos de tempo a partir do qual &
		serao calculados os campos medios de velocidade e pressao NPASS:'
		read (*, *) n1
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (201, 7) n1
        else
            write (201, 3) 'NPASS:', n1
        end if
		
		write (*, 1) 'Indique se sera ou nao usado algum modelo de turbulencia na analise INDTURB:'
		write (*, 1) '     INDTURB = 0  --->  Analise sem a presenca de turbulencia'
		write (*, 1) '     INDTURB = 1  --->  Analise com turbulencia e modelo sub-malha classico'
		write (*, 1) '     INDTURB = 2  --->  Analise com turbulencia e modelo sub-malha dinamico, com 1 ponto de integracao'
		write (*, 1) '     INDTURB = 3  --->  Analise com turbulencia e modelo sub-malha dinamico, com 8 pontos de integracao'
		read (*, *) INDTURB
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (201, 7) INDTURB
        else
            write (201, 3) 'INDTURB:', INDTURB
        end if
		
		!write (*, *) 'Indique o parametro seletivo de massa na equacao de continuidade (padrao 0,9) ELUMP1:'
		!read (*, *) r
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (201, 8) 0.9
        else
            write (201, 4) 'ELUMP1:', 0.9
        end if
		
		!write (*, *) 'Indique o parametro seletivo de massa na equacao de energia (padrao 0,9) ELUMP2:'
		!read (*, *) r
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (201, 8) 0.9
        else
            write (201, 4) 'ELUMP2:', 0.9
        end if
		
		!write (*, *) 'Indique o parametro seletivo de massa na equacao de especie (padrao 0,9) ELUMP3:'
		!read (*, *) r
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (201, 8) 0.9
        else
            write (201, 4) 'ELUMP3:', 0.9
        end if
		
		
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            !write (*, 1) 'Indique o numero de threads por bloco de processamento nThreads:'
            !write (*, 1) '	Recomenda-se utilizar uma potencia de 2'
            !read (*, *) n1
            write (201, 7) 128
        else
            ! Para o registro das componentes flutuantes do campo de velocidades e pressoes
			! Nao implementado em Hermann_nuevo
			!write (*, *) 'Indique o valor NRMS:'
			!read (*, *) n1
            write (201, 3) 'NRMS:', 0
        end if
		
		close (201)
	end if
    
	!****************************************************************
	! CRIAÇÃO DO ARQUIVO DE PROPRIEDADES DO FLUIDO .PRO
	!****************************************************************
	
	write (*, 1) 'Gerar arquivo de propriedades do fluido?        0 = Nao   1 = Sim'
	read (*, *) i
	
	if (i .eq. 1) then
		ptname (7:14) = '.pro    '
		open (202, file=ptname, status='unknown')
		
		write (*, 1) 'Indique o modulo do vetor velocidade da corrente nao perturbada (Infinito) Vinf:'
		write (*, 1) 'Este valor da referencia para a velocidade do escoamento.'
		read (*, *) r
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (202, 8) r
        else
            write (202, 4) 'Vinf:', r
        end if
		write (*, 1) ''
		
		write (*, 1) 'Indique a velocidade de propagacao do som VelSom:'
		write (*, 1) 'E comum a utilizacao de valores baixos para aumentar o &
		passo de tempo e melhorar a estabilidade do metodo. De forma geral, &
		procura-se adotar velocidades de forma que a velocidade de referencia &
		atinja Mach = 0,2, ou seja, VelSom = 5*Vinf.'
		write (*, 1) '     Valores reais:'
		write (*, 1) '        Ar (20C) -> 343 m/s'
		write (*, 1) '        Agua     -> 1480 m/s'
		read (*, *) r
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (202, 8) r
        else
            write (202, 4) 'VelSom:', r
        end if
		write (*, 1) ''
		
		write (*, 1) 'Indique a viscosidade cinematica do fluido [m^2/s] ViscCin:'
		write (*, 1) 'Este valor pode ser manipulado para se obter o numero de &
		Reynolds desejado para uma simulacao.'
		write (*, 1) '     Valores reais:'
		write (*, 1) '        Ar (20C) -> 15,06 * 10^-6 m^2/s'
		write (*, 1) '        Agua     -> 1,004 * 10^-6 m^2/s'
		read (*, *) r
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (202, 8) r
        else
            write (202, 4) 'ViscCin:', r
        end if
		write (*, 1) ''
		
		! Este valor só é utilizado para escoamentos compressiveis
		!write (*, *) 'Indique a viscosidade volumetrica por unidade de massa especifica do fluido [m²/s] ViscVol:'
		!read (*, *) r
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (202, 8) 0.0
        else
            write (202, 4) 'ViscVol:', 0.0
        end if
		
		if (CodeType .ne. 2) then
			write (*, 1) 'Indique a pressao da corrente nao perturbada (Infinito) Pinf:'
			read (*, *) r
			if (CodeType .eq. 1) then 
				write (202, 8) r
			else
				write (202, 4) 'Pinf:', r
			end if
		end if
		write (*, 1) ''
		
		write (*, 1) 'Indique a densidade do fluido RHOInf:'
		write (*, 1) '     Valores reais:'
		write (*, 1) '        Ar (0C e 100 kPa)      -> 1,2754 kg/m^3'
		write (*, 1) '        Ar (20C e 101,325 kPa) -> 1,2041 kg/m^3'
		write (*, 1) '        Agua                    -> 997 kg/m^3'
		read (*, *) r
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (202, 8) r
        else
            write (202, 4) 'RHOInf:', r
        end if
		write (*, 1) ''
		
		! Não serão feitas analises termodinamicas
		!write (*, *) 'Indique a temperatura do fluido TInf:'
		!read (*, *) r
		if (CodeType .ne. 2) then
			if (CodeType .eq. 1) then 
				write (202, 8) 0.0
			else
				write (202, 4) 'Tinf:', 0.0
			end if
		end if
		
		!write (*, *) 'Indique o calor específica a volume constante Cv:'
		!write (*, *) '     Valores tipicos:'
		!write (*, *) '        Ar (20C) -> 1012 J/kg'
		!write (*, *) '        Agua      -> 4186 J/kg'
		!read (*, *) r
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (202, 8) 0.00000001
        else
            write (202, 4) 'Cv:', 0.00000001
        end if
		
		!write (*, *) 'Indique o coeficiente de condutividade termica Kdif:'
		!write (*, *) '     Valores tipicos:'
		!write (*, *) '        Ar (27C)   -> 0,03'
		!write (*, *) '        Agua (27C) -> 0,61'
		!read (*, *) r
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (202, 8) 0.00000001
        else
            write (202, 4) 'Kdif:', 0.00000001
        end if
		
		!write (*, *) 'Indique o coeficiente de expansao volumetrica do fluido Beta:'
		!write (*, *) '     Valores tipicos:'
		!write (*, *) '        Ar   -> 3,39 * 10^-4'
		!write (*, *) '        Agua -> 1,3 * 10^-4'
		!read (*, *) r
		if (CodeType .ne. 2) then
			if (CodeType .eq. 1) then 
				write (202, 8) 0.0
			else
				write (202, 4) 'Beta:', 0.0
			end if
		end if
		
		!write (*, 1) 'Indique as componentes do vetor de aceleracao da gravidade segundo os eixos coordenados gr1, gr2 e gr3:'
		!write (*, 1) 'A nao ser que senha uma aplicacao especifica em mente, se &
		!utiliza os tres iguais a zero.'
		!read (*, *) r, r2, r3
		if (CodeType .ne. 2) then
			if (CodeType .eq. 1) then 
				write (202, 9) 0.0, 0.0, 0.0
			else
				write (202, 5) 'gri:', 0.0, 0.0, 0.0
			end if
		end if
		!write (*, 1) ''
		
		if (INDTURB .eq. 1) then
			write (*, 1) 'Indique o coeficiente de Smagorinsky Cs:'
			read (*, *) r
		else
			r = 0.1
		end if
		
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
			write (202, 8) r
		else
			write (202, 4) 'Cs:', r
		end if
		write (*, 1) ''
		
		! TODO - Verificar se Kdift e Dabt ainda sao validos, Hermann_nuevo conta com Prtlt e Schit no lugar
		!write (*, *) 'Indique o coeficiente de condutividade termica turbulenta Kdift:'
		!read (*, *) r
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (202, 8) 0.00000001
        else
            write (202, 4) 'Kdift:', 0.00000001
        end if
		
		!write (*, *) 'Indique o coeficiente de difusividade massica turbulenta Dabt:'
		!read (*, *) r
		if (CodeType .ne. 2) then
			if (CodeType .eq. 1) then 
				write (202, 8) 0.00000001
			else
				write (202, 4) 'Dabt:', 0.00000001
			end if
		end if
		
		!write (*, *) 'Indique o coeficiente de difusidade massica Dab:'
		!read (*, *) r
		if (CodeType .eq. 1 .or. CodeType .eq. 2) then
            write (202, 8) 0.00000001
        else
            write (202, 4) 'Dab:', 0.00000001
        end if
		
		!write (*, *) 'Indique a densidade da especie A RHOinfA:'
		!read (*, *) r
		if (CodeType .ne. 2) then
			if (CodeType .eq. 1) then 
				write (202, 8) 0.0
			else
				write (202, 4) 'RHOinfA:', 0.0
			end if
		end if
		
		write (*, 1) 'Indique o numero de objetos submersos:'
		read (*, *) MNOBJ
        if (CodeType .eq. 0) then
			write (202, 3) 'MNOBJ:', MNOBJ
        end if
		write (*, 1) ''
		
		allocate ( objCenter(MNOBJ,3) )
		
		do i=1,MNOBJ
			write (*, 3) 'OBJETO', i
			write (*, 1) ''
			
			write (*, 1) 'Informe o comprimento caracteristico do corpo imerso Lchar:'
			read (*, *) r
			if (CodeType .eq. 1 .or. CodeType .eq. 2) then
				write (202, 8) r
			else
				write (202, 4) 'Lchar:', r
			end if
			write (*, 1) ''
			
			write (*, 1) 'Informe a dimensao caracteristica do corpo imerso Dchar:'
			read (*, *) r
			if (CodeType .eq. 1 .or. CodeType .eq. 2) then
				write (202, 8) r
			else
				write (202, 4) 'Dchar:', r
			end if
			write (*, 1) ''
			
			write (*, 1) 'Informe a coordenada X do centro de gravidade do corpo imerso Xobj:'
			read (*, *) objCenter(i,1)
			if (CodeType .eq. 1 .or. CodeType .eq. 2) then
				write (202, 8) objCenter(i,1)
			else
				write (202, 4) 'Xobj:', objCenter(i,1)
			end if
			write (*, 1) ''
			
			write (*, 1) 'Informe a coordenada Y do centro de gravidade do corpo imerso Yobj:'
			read (*, *) objCenter(i,2)
			if (CodeType .eq. 1 .or. CodeType .eq. 2) then
				write (202, 8) objCenter(i,2)
			else
				write (202, 4) 'Yobj:', objCenter(i,2)
			end if
			write (*, 1) ''
			
			write (*, 1) 'Informe a coordenada Z do centro de gravidade do corpo imerso Zobj:'
			read (*, *) objCenter(i,3)
			if (CodeType .eq. 1 .or. CodeType .eq. 2) then
				write (202, 8) objCenter(i,3)
			else
				write (202, 4) 'Zobj:', objCenter(i,3)
			end if
			write (*, 1) ''
			
			write (*, 1) 'Informe o numero de faces de contorno do corpo imerso NFCN:'
			read (*, *) n1
            if (CodeType .eq. 0) then
				write (202, 3) 'NFCN:', n1
            end if
			write (*, 1) ''
			
			write (*, 1) 'Informe o numero de nos de contorno do corpo imerso NNCN:'
			read (*, *) n2
			if (CodeType .eq. 1 .or. CodeType .eq. 2) then
				write (202, 10) n1, n2
			else
				write (202, 3) 'NNCN:', n2
			end if
			write (*, 1) ''
			
		end do
		
		close (202)
	end if
	
    !****************************************************************
    ! LEITURA DO ARQUIVO DA MALHA .MSH
    !****************************************************************
    
    ptname (7:14) = '.msh    '
    open (102, file=ptname, status='old')
    
    bool = .TRUE.
    do while ( bool )
        read (102,*) header
		select case (header)
        case ('$MeshFormat')
            ! TODO - Na realidade seria interessante armazenar essas informacões
            read (102,*) garbr, garbi, garbi
            read (102,*) garbc
		
		case ('$PhysicalNames')
			read (102,*) n1
			do i=1,n1
				read (102,*) garbi, garbi, garbc
			end do
			read (102,*) garbc
            
        case ('$Nodes')
            
            !*****************************************************
            ! LEITURA E ARMAZENAMENTO DAS COORDENADAS DOS NÓS
            !*****************************************************
			write(*,*) header
            
            ptname (7:14) = '000.cor '
            open (204, file=ptname, status='unknown')
			
			ptname (10:14) = '.v   '
            open (207, file=ptname, status='unknown')
			
			ptname (10:14) = '.pr  '
            open (208, file=ptname, status='unknown')
			
            if (T) then
				ptname (10:14) = '.t   '
				open (209, file=ptname, status='unknown')
            end if
	
            if (RHA) then
				ptname (10:14) = '.rha '
				open (210, file=ptname, status='unknown')
            end if
                
            ! Lógica para ler os dados dos nós da malha
            read (102, *) ntn
			write (204, 11) ntn
            
            ! Alocando o tamanho dos vetores com as coordenadas dos nós
            allocate( xx(ntn), yy(ntn), zz(ntn) )
            do i=1,ntn
                read (102, *) garbi, xx(i), yy(i), zz(i)   ! o número do nó nao e relevante de armazenar
                write (204, 12) xx(i), yy(i), zz(i)
				
				! Condições iniciais
				select case (indexCC(initialFieldGroup,2))
				case (2)
					! Campo Constante
					do j=1,countCCTypes(2)
						if (indexConstCC(j,1) .eq. initialFieldGroup .and. indexConstCC(j,2) .eq. 0) then
							exit
						end if
					end do
					
					V1 = valueConstCC(j)
				case (3)
					! Campo em Lei de Potência
					do j=1,countCCTypes(3)
						if (indexPowerCC(j,1) .eq. initialFieldGroup .and. indexPowerCC(j,2) .eq. 0) then
							exit
						end if
					end do
					
					V1 = PowerLaw(yy(i), valuePowerCC(j,1), valuePowerCC(j,2), valuePowerCC(j,3), valuePowerCC(j,4))
				case (4)
					! Campo em Lei Logaritmica
					! TODO - Implementar
					V1 = -10000.0
				case default
					! Campo nulo
					V1 = 0.0
				end select
				
				select case (indexCC(initialFieldGroup,3))
				case (2)
					! Campo Constante
					do j=1,countCCTypes(2)
						if (indexConstCC(j,1) .eq. initialFieldGroup .and. indexConstCC(j,2) .eq. 1) then
							exit
						end if
					end do
					
					V2 = valueConstCC(j)
				case default
					! Campo nulo
					V2 = 0.0
				end select
				
				select case (indexCC(initialFieldGroup,4))
				case (2)
					! Campo Constante
					do j=1,countCCTypes(2)
						if (indexConstCC(j,1) .eq. initialFieldGroup .and. indexConstCC(j,2) .eq. 2) then
							exit
						end if
					end do
					
					V3 = valueConstCC(j)
				case (3)
					! Campo em Lei de Potência
					do j=1,countCCTypes(3)
						if (indexPowerCC(j,1) .eq. initialFieldGroup .and. indexPowerCC(j,2) .eq. 2) then
							exit
						end if
					end do
					
					V3 = PowerLaw(yy(i), valuePowerCC(j,1), valuePowerCC(j,2), valuePowerCC(j,3), valuePowerCC(j,4))
				case (4)
					! Campo em Lei Logaritmica
					! TODO - Implementar
					V3 = -10000.0
				case default
					! Campo nulo
					V3 = 0.0
				end select
				
				select case (indexCC(initialFieldGroup,5))
				case (2)
					! Campo Constante
					do j=1,countCCTypes(2)
						if (indexConstCC(j,1) .eq. initialFieldGroup .and. indexConstCC(j,2) .eq. 3) then
							exit
						end if
					end do
					
					P = valueConstCC(j)
				case default
					! Campo nulo
					P = 0.0
				end select
				
				write (207, 12) V1, V2, V3
				write (208, 13) P, P
                
                if (T) write (209, 13) 0.0, 0.0
				if (RHA) write (210, 13) 0.0, 0.0
            end do
            ! TODO - É possível implementar uma verificacao para ver se essa linha e de fato '$EndNodes'
            read (102, *) garbc
            
            close (204)
			close (207)
			close (208)
            if (T) close (209)
			if (RHA) close (210)
            
        case ('$Elements')
                        
            !*****************************************************
            ! LEITURA E ARMAZENAMENTO DOS DADOS DOS ELEMENTOS
            !*****************************************************
			write(*,*) header
            
            ! TODO - Adicionar lógica para os elementos
            read (102,*) ntele
			write(*,*) ntele
            
            ! TODO - Idealmente seria conveniente identificar o tipo de elemento pelo índice do arquivo .msh
            ! Calculando o número de elementos pontuais e unidimensionais
            n = ntele - nquad - nhexa
			write(*,*) n
            
            ! Elementos pontuais e unidimensionais nao sao necessarios
            do i=1,n
                read (102,*) garbi
            end do
            
			!************************************
            ! ELEMENTOS QUADRILÁTEROS
			!************************************
			
			write(*,*) 'Quads'
			
			! Alocando e inicializando o vetor auxiliar
			allocate( nccV1potnodes(ntn), nccV3potnodes(ntn), &
					vccV1potnodes(ntn), vccV3potnodes(ntn), &
					nccV1lognodes(ntn), nccV3lognodes(ntn), &
					vccV1lognodes(ntn), vccV3lognodes(ntn), &
					nccV1nodes(ntn), nccV2nodes(ntn), nccV3nodes(ntn), nccPnodes(ntn), &
					vccV1nodes(ntn), vccV2nodes(ntn), vccV3nodes(ntn), vccPnodes(ntn), &
					CSsupELEM(MNOBJ,nquad), CSsupKONE(MNOBJ,nquad,4), &
					CSsupAREA(MNOBJ,nquad), CSsupCOUNT(MNOBJ), &
					CSnnnNODE(MNOBJ,ntn), CSnnnNORM(MNOBJ,ntn,3) )
            nccV1potnodes = 0
            nccV3potnodes = 0
            vccV1potnodes = 0
            vccV3potnodes = 0
			nccV1nodes = 0
			nccV2nodes = 0
			nccV3nodes = 0
			nccPnodes = 0
			vccV1nodes = 0
			vccV2nodes = 0
			vccV3nodes = 0
			vccPnodes = 0
            nV1potnodes = 0
            nV3potnodes = 0
			nV1lognodes = 0
			nV3lognodes = 0
			nV1nodes = 0
			nV2nodes = 0
			nV3nodes = 0
			nPnodes = 0
			CSsupELEM = 0
			CSsupKONE = 0
			CSsupAREA = 0
			CSsupCOUNT = 0
			CSnnnNODE = 0
			CSnnnNORM = 0
			
            do i=1,nquad
                read (102,*) garbi, garbi, garbi, group, face, n1, n2, n3, n4 ! n1, n2, n3 e n4 sao os nós do elemento
				
				write(*,100) i, nquad, 0, nhexa, 100 * real(i) / real(nquad + nhexa)
                
				! V1
				select case (indexCC(group,2))
				case (1)
					! CC Nula
					
					! Verificar se os nós ja estao na lista
					quadnodes(1) = n1
					quadnodes(2) = n2
					quadnodes(3) = n3
					quadnodes(4) = n4
						
					call AddNewCCNodes(ntn, nV1nodes, 4, nccV1nodes, quadnodes)
                    call UpdateUniformCCValues(ntn, nccV1nodes, quadnodes, 0.0d0, vccV1nodes)
				case (2)
					! CC Constante
					
					! Verificar se os nós ja estao na lista
					quadnodes(1) = n1
					quadnodes(2) = n2
					quadnodes(3) = n3
					quadnodes(4) = n4
					
					do j=1,countCCTypes(2)
						if (indexConstCC(j,1) .eq. group .and. indexConstCC(j,2) .eq. 0) then
							exit
						end if
					end do
						
					call AddNewCCNodes(ntn, nV1nodes, 4, nccV1nodes, quadnodes)
                    call UpdateUniformCCValues(ntn, nccV1nodes, quadnodes, valueConstCC(j), vccV1nodes)
				case (3)
					! CC Potência
					
					! Verificar se os nós ja estao na lista
					quadnodes(1) = n1
					quadnodes(2) = n2
					quadnodes(3) = n3
					quadnodes(4) = n4
					
					do j=1,countCCTypes(3)
						if (indexPowerCC(j,1) .eq. group .and. indexPowerCC(j,2) .eq. 0) then
							exit
						end if
					end do
                    
					call AddNewCCNodes(ntn, nV1potnodes, 4, nccV1potnodes, quadnodes)
                    call UpdatePowerLawCCValues(ntn, nccV1potnodes, n1, yy(n1), valuePowerCC(j,:), vccV1potnodes)
                    call UpdatePowerLawCCValues(ntn, nccV1potnodes, n2, yy(n2), valuePowerCC(j,:), vccV1potnodes)
                    call UpdatePowerLawCCValues(ntn, nccV1potnodes, n3, yy(n3), valuePowerCC(j,:), vccV1potnodes)
                    call UpdatePowerLawCCValues(ntn, nccV1potnodes, n4, yy(n4), valuePowerCC(j,:), vccV1potnodes)
				case (4)
					! CC Log
					
					! Verificar se os nós ja estao na lista
					quadnodes(1) = n1
					quadnodes(2) = n2
					quadnodes(3) = n3
					quadnodes(4) = n4
					
					do j=1,countCCTypes(4)
						if (indexLogCC(j,1) .eq. group .and. indexLogCC(j,2) .eq. 0) then
							exit
						end if
					end do
                    
					call AddNewCCNodes(ntn, nV1lognodes, 4, nccV1lognodes, quadnodes)
                    call UpdateLogLawCCValues(ntn, nccV1lognodes, n1, yy(n1), valueLogCC(j,:), vccV1lognodes)
                    call UpdateLogLawCCValues(ntn, nccV1lognodes, n2, yy(n2), valueLogCC(j,:), vccV1lognodes)
                    call UpdateLogLawCCValues(ntn, nccV1lognodes, n3, yy(n3), valueLogCC(j,:), vccV1lognodes)
                    call UpdateLogLawCCValues(ntn, nccV1lognodes, n4, yy(n4), valueLogCC(j,:), vccV1lognodes)
				end select
                
				! V2
				select case (indexCC(group,3))
				case (1)
					! CC Nula
					
					! Verificar se os nós ja estao na lista
					quadnodes(1) = n1
					quadnodes(2) = n2
					quadnodes(3) = n3
					quadnodes(4) = n4
						
					call AddNewCCNodes(ntn, nV2nodes, 4, nccV2nodes, quadnodes)
                    call UpdateUniformCCValues(ntn, nccV2nodes, quadnodes, 0.0d0, vccV2nodes)
				case (2)
					! CC Constante
					
					! Verificar se os nós ja estao na lista
					quadnodes(1) = n1
					quadnodes(2) = n2
					quadnodes(3) = n3
					quadnodes(4) = n4
					
					do j=1,countCCTypes(2)
						if (indexConstCC(j,1) .eq. group .and. indexConstCC(j,2) .eq. 1) then
							exit
						end if
					end do
						
					call AddNewCCNodes(ntn, nV2nodes, 4, nccV2nodes, quadnodes)
                    call UpdateUniformCCValues(ntn, nccV2nodes, quadnodes, valueConstCC(j), vccV2nodes)
				end select
				
				! V3
				select case (indexCC(group,4))
				case (1)
					! CC Nula
					
					! Verificar se os nós ja estao na lista
					quadnodes(1) = n1
					quadnodes(2) = n2
					quadnodes(3) = n3
					quadnodes(4) = n4
						
					call AddNewCCNodes(ntn, nV3nodes, 4, nccV3nodes, quadnodes)
                    call UpdateUniformCCValues(ntn, nccV3nodes, quadnodes, 0.0d0, vccV3nodes)
				case (2)
					! CC Constante
					
					! Verificar se os nós ja estao na lista
					quadnodes(1) = n1
					quadnodes(2) = n2
					quadnodes(3) = n3
					quadnodes(4) = n4
					
					do j=1,countCCTypes(2)
						if (indexConstCC(j,1) .eq. group .and. indexConstCC(j,2) .eq. 2) then
							exit
						end if
					end do
						
					call AddNewCCNodes(ntn, nV3nodes, 4, nccV3nodes, quadnodes)
                    call UpdateUniformCCValues(ntn, nccV3nodes, quadnodes, valueConstCC(j), vccV3nodes)
				case (3)
					! CC Potência
					
					! Verificar se os nós ja estao na lista
					quadnodes(1) = n1
					quadnodes(2) = n2
					quadnodes(3) = n3
					quadnodes(4) = n4
					
					do j=1,countCCTypes(3)
						if (indexPowerCC(j,1) .eq. group .and. indexPowerCC(j,2) .eq. 2) then
							exit
						end if
					end do
                    
					call AddNewCCNodes(ntn, nV3potnodes, 4, nccV3potnodes, quadnodes)
                    call UpdatePowerLawCCValues(ntn, nccV3potnodes, n1, yy(n1), valuePowerCC(j,:), vccV3potnodes)
                    call UpdatePowerLawCCValues(ntn, nccV3potnodes, n2, yy(n2), valuePowerCC(j,:), vccV3potnodes)
                    call UpdatePowerLawCCValues(ntn, nccV3potnodes, n3, yy(n3), valuePowerCC(j,:), vccV3potnodes)
                    call UpdatePowerLawCCValues(ntn, nccV3potnodes, n4, yy(n4), valuePowerCC(j,:), vccV3potnodes)
				case (4)
					! CC Log
					
					! Verificar se os nós ja estao na lista
					quadnodes(1) = n1
					quadnodes(2) = n2
					quadnodes(3) = n3
					quadnodes(4) = n4
					
					do j=1,countCCTypes(4)
						if (indexLogCC(j,1) .eq. group .and. indexLogCC(j,2) .eq. 2) then
							exit
						end if
					end do
                    
					call AddNewCCNodes(ntn, nV3lognodes, 4, nccV3lognodes, quadnodes)
                    call UpdateLogLawCCValues(ntn, nccV3lognodes, n1, yy(n1), valueLogCC(j,:), vccV3lognodes)
                    call UpdateLogLawCCValues(ntn, nccV3lognodes, n2, yy(n2), valueLogCC(j,:), vccV3lognodes)
                    call UpdateLogLawCCValues(ntn, nccV3lognodes, n3, yy(n3), valueLogCC(j,:), vccV3lognodes)
                    call UpdateLogLawCCValues(ntn, nccV3lognodes, n4, yy(n4), valueLogCC(j,:), vccV3lognodes)
				end select
                
				! P
				select case (indexCC(group,5))
				case (1)
					! CC Nula
					
					! Verificar se os nós ja estao na lista
					quadnodes(1) = n1
					quadnodes(2) = n2
					quadnodes(3) = n3
					quadnodes(4) = n4
						
					call AddNewCCNodes(ntn, nPnodes, 4, nccPnodes, quadnodes)
                    call UpdateUniformCCValues(ntn, nccPnodes, quadnodes, 0.0d0, vccPnodes)
				case (2)
					! CC Constante
					
					! Verificar se os nós ja estao na lista
					quadnodes(1) = n1
					quadnodes(2) = n2
					quadnodes(3) = n3
					quadnodes(4) = n4
					
					do j=1,countCCTypes(2)
						if (indexConstCC(j,1) .eq. group .and. indexConstCC(j,2) .eq. 3) then
							exit
						end if
					end do
						
					call AddNewCCNodes(ntn, nPnodes, 4, nccPnodes, quadnodes)
                    call UpdateUniformCCValues(ntn, nccPnodes, quadnodes, valueConstCC(j), vccPnodes)
				end select

				! Objetos Submersos
				do j=1,MNOBJ
					if (group .eq. OBJPhysGroup(j)) then
						! Aumentar o contador que armazena o número de elementos do corpo submerso
						CSsupCOUNT(j) = CSsupCOUNT(j) + 1
						
						! Calcular a área do elemento
						call GetQuadArea(CSsupAREA(j,CSsupCOUNT(j)), &
											xx(n1), xx(n2), xx(n3), xx(n4), &
											yy(n1), yy(n2), yy(n3), yy(n4), &
											zz(n1), zz(n2), zz(n3), zz(n4) )
											
						! Conectividade
						CSsupKONE(j,CSsupCOUNT(j),1) = n1
						CSsupKONE(j,CSsupCOUNT(j),2) = n2
						CSsupKONE(j,CSsupCOUNT(j),3) = n3
						CSsupKONE(j,CSsupCOUNT(j),4) = n4
						
						! Verificar quais nós já tem algum vetor normal armazenado
						quadnodes = 0
						counter = 0
						do k2=1,ntn
							if (CSnnnNODE(j,k2) .eq. n1) then
								quadnodes(1) = 1
							else if (CSnnnNODE(j,k2) .eq. n2) then
								quadnodes(2) = 1
							else if (CSnnnNODE(j,k2) .eq. n3) then
								quadnodes(3) = 1
							else if (CSnnnNODE(j,k2) .eq. n4) then
								quadnodes(4) = 1
							else if (CSnnnNODE(j,k2) .eq. 0) then
								counter = k2
								exit
							end if
						end do
						
						! Para os nós que não tiver vetor normal já armazenado, incluir na lista
						if (quadnodes(1) .eq. 0) then
							CSnnnNODE(j,counter) = n1
							counter = counter + 1
						end if
						if (quadnodes(2) .eq. 0) then
							CSnnnNODE(j,counter) = n2
							counter = counter + 1
						end if
						if (quadnodes(3) .eq. 0) then
							CSnnnNODE(j,counter) = n3
							counter = counter + 1
						end if
						if (quadnodes(4) .eq. 0) then
							CSnnnNODE(j,counter) = n4
							counter = counter + 1
						end if
						
						! Calcular as componentes do vetor normal de cada nó
						! Todos os nós usam a mesma subrotina, mas com permutação ciclica
						! As componentes são adicionadas a componentes já armazenadas e normalizadas antes de armazenar no arquivo
						do k2=1,ntn
							if (CSnnnNODE(j,k2) .eq. n1) then
								call GetNodeNormal(NNN1, NNN2, NNN3, &
											xx(n1), xx(n2), xx(n3), xx(n4), &
											yy(n1), yy(n2), yy(n3), yy(n4), &
											zz(n1), zz(n2), zz(n3), zz(n4) )
								! Corrigir a direção da normal com base nas coordenadas do centroide do objeto
								call CorrectNormalDirection(NNN1, NNN2, NNN3, xx(n1), yy(n1), zz(n1), objCenter(j,1), objCenter(j,2), objCenter(j,3))
								CSnnnNORM(j,k2,1) = CSnnnNORM(j,k2,1) + NNN1
								CSnnnNORM(j,k2,2) = CSnnnNORM(j,k2,2) + NNN2
								CSnnnNORM(j,k2,3) = CSnnnNORM(j,k2,3) + NNN3
							else if (CSnnnNODE(j,k2) .eq. n2) then
								call GetNodeNormal(NNN1, NNN2, NNN3, &
											xx(n2), xx(n3), xx(n4), xx(n1), &
											yy(n2), yy(n3), yy(n4), yy(n1), &
											zz(n2), zz(n3), zz(n4), zz(n1) )
								! Corrigir a direção da normal com base nas coordenadas do centroide do objeto
								call CorrectNormalDirection(NNN1, NNN2, NNN3, xx(n2), yy(n2), zz(n2), objCenter(j,1), objCenter(j,2), objCenter(j,3))
								CSnnnNORM(j,k2,1) = CSnnnNORM(j,k2,1) + NNN1
								CSnnnNORM(j,k2,2) = CSnnnNORM(j,k2,2) + NNN2
								CSnnnNORM(j,k2,3) = CSnnnNORM(j,k2,3) + NNN3
							else if (CSnnnNODE(j,k2) .eq. n3) then
								call GetNodeNormal(NNN1, NNN2, NNN3, &
											xx(n3), xx(n4), xx(n1), xx(n2), &
											yy(n3), yy(n4), yy(n1), yy(n2), &
											zz(n3), zz(n4), zz(n1), zz(n2) )
								! Corrigir a direção da normal com base nas coordenadas do centroide do objeto
								call CorrectNormalDirection(NNN1, NNN2, NNN3, xx(n3), yy(n3), zz(n3), objCenter(j,1), objCenter(j,2), objCenter(j,3))
								CSnnnNORM(j,k2,1) = CSnnnNORM(j,k2,1) + NNN1
								CSnnnNORM(j,k2,2) = CSnnnNORM(j,k2,2) + NNN2
								CSnnnNORM(j,k2,3) = CSnnnNORM(j,k2,3) + NNN3
							else if (CSnnnNODE(j,k2) .eq. n4) then
								call GetNodeNormal(NNN1, NNN2, NNN3, &
											xx(n4), xx(n1), xx(n2), xx(n3), &
											yy(n4), yy(n1), yy(n2), yy(n3), &
											zz(n4), zz(n1), zz(n2), zz(n3) )
								! Corrigir a direção da normal com base nas coordenadas do centroide do objeto
								call CorrectNormalDirection(NNN1, NNN2, NNN3, xx(n4), yy(n4), zz(n4), objCenter(j,1), objCenter(j,2), objCenter(j,3))
								CSnnnNORM(j,k2,1) = CSnnnNORM(j,k2,1) + NNN1
								CSnnnNORM(j,k2,2) = CSnnnNORM(j,k2,2) + NNN2
								CSnnnNORM(j,k2,3) = CSnnnNORM(j,k2,3) + NNN3
							else if (CSnnnNODE(j,k2) .eq. 0) then
								exit
							end if
						end do
					end if
				end do
				
			end do
			
			do i=1,MNOBJ
				! Normalizar os vetores normais
				do j=1,ntn
					if (CSnnnNODE(i,j) .ne. 0) then
						r = dsqrt( CSnnnNORM(i,j,1) * CSnnnNORM(i,j,1) + &
								CSnnnNORM(i,j,2) * CSnnnNORM(i,j,2) + &
								CSnnnNORM(i,j,3) * CSnnnNORM(i,j,3) )
						CSnnnNORM(i,j,1) = CSnnnNORM(i,j,1) / r
						CSnnnNORM(i,j,2) = CSnnnNORM(i,j,2) / r
						CSnnnNORM(i,j,3) = CSnnnNORM(i,j,3) / r
					else
						exit
					end if
				end do
			end do
			
			! Armazenamento dos dados de contorno V1, V2 e V3
			ptname (7:14) = 'CC.bv   '
            open (211, file=ptname, status='unknown')
			
			! Primeiramente armazena o número de nós com a respectiva condicao de contorno.
			! Em seguida, para cada nó, se armazena o número global do nó e o valor da condicao.
			write (211, 11) nV1potnodes + nV1lognodes + nV1nodes
            do i=1,nV1potnodes
				write (211, 15) nccV1potnodes(i), vccV1potnodes(i)
			end do
			do i=1,nV1lognodes
				write (211, 15) nccV1lognodes(i), vccV1lognodes(i)
			end do
			do i=1,nV1nodes
				write (211, 15) nccV1nodes(i), vccV1nodes(i)
			end do
			
			write (211, *) nV2nodes
			do i=1,nV2nodes
				write (211, 15) nccV2nodes(i), vccV2nodes(i)
			end do
			
			write (211, *) nV3potnodes + nV3lognodes + nV3nodes
            do i=1,nV3potnodes
				write (211, 15) nccV3potnodes(i), vccV3potnodes(i)
			end do
			do i=1,nV3lognodes
				write (211, 15) nccV3lognodes(i), vccV3lognodes(i)
			end do
			do i=1,nV3nodes
				write (211, 15) nccV3nodes(i), vccV3nodes(i)
			end do
			
			close (211)
			
			! Armazenamento dos dados de contorno P
			ptname (7:14) = 'CC.bp   '
            open (212, file=ptname, status='unknown')
			
			write (212, 11) nPnodes
			do i=1,nPnodes
				write (212, 15) nccPnodes(i), vccPnodes(i)
			end do

			close (212)
			
            !*****************************************
            ! CONECTIVIDADE DOS ELEMENTOS HEXAÉDRICOS
			!*****************************************
			
			ptname (7:14) = '.con    '
			open (203, file=ptname, status='unknown')
			write (203, 11) nhexa
			
			allocate ( NEIBOR(ntn) )
			NEIBOR = 0
			
			do i=1,nhexa
                write(*,100) nquad, nquad, i, nhexa, 100 * real(nquad + i) / real(nquad + nhexa)
                
				! TODO - Armazenar a conectividade dos elementos
				read (102, *) garbi, garbi, garbi, garbi, garbi, n1, n2, n3, n4, n5, n6, n7, n8 ! n1, n2, n3, n4, n5, n6, n7 e n8 sao os nós do elemento
				write (203, 14) n1, n2, n3, n4, n5, n6, n7, n8
				
				NEIBOR(n1) = NEIBOR(n1) + 1
				NEIBOR(n2) = NEIBOR(n2) + 1
				NEIBOR(n3) = NEIBOR(n3) + 1
				NEIBOR(n4) = NEIBOR(n4) + 1
				NEIBOR(n5) = NEIBOR(n5) + 1
				NEIBOR(n6) = NEIBOR(n6) + 1
				NEIBOR(n7) = NEIBOR(n7) + 1
				NEIBOR(n8) = NEIBOR(n8) + 1
				
				! Se houver corpo submersos:
				! Encontrar o número do elemento hexa que cada face do contorno sólido faz parte
				! 	IELCS
				do j=1,MNOBJ
					! Para cada face de controno sólido k1
					do k1=1,CSsupCOUNT(j)
						counter = 0
						! Para cada nó k2 da face k1
						do k2=1,4
							if (CSsupKONE(j,k1,k2) .eq. n1) then
								counter = counter + 1
							else if (CSsupKONE(j,k1,k2) .eq. n2) then
								counter = counter + 1
							else if (CSsupKONE(j,k1,k2) .eq. n3) then
								counter = counter + 1
							else if (CSsupKONE(j,k1,k2) .eq. n4) then
								counter = counter + 1
							else if (CSsupKONE(j,k1,k2) .eq. n5) then
								counter = counter + 1
							else if (CSsupKONE(j,k1,k2) .eq. n6) then
								counter = counter + 1
							else if (CSsupKONE(j,k1,k2) .eq. n7) then
								counter = counter + 1
							else if (CSsupKONE(j,k1,k2) .eq. n8) then
								counter = counter + 1
							end if

						end do
						! Se os 4 nós fizerem parte do elemento hexa
						! Significa que esse é o elemento de contorno e o index é i
						if (counter .eq. 4) then
							CSsupELEM(j,k1) = i
						end if
					end do
				end do
			end do
			close (203)
			
			ptname (7:14) = 'CS      '
			do j=1,MNOBJ
				! Escrever arquivos do objeto submerso
				! Código retirado e adaptado de PatransferIFE3D.for
				NB1 = j/10 + 48
				NB2 = j - (j/10)*10 + 48
				
				
				ptname (9:10) = char(NB1)//char(NB2)
				ptname (11:14) = '.sup'
				open (301, file=ptname, status='unknown')
				ptname (11:14) = '.nnn'
				open (302, file=ptname, status='unknown')
				
				! Armazenamento do aqruivo .sup
				! Com o número do elemento hexa do contorno sólido do objeto j
				! E a conectividade da face k1
				do k1=1,CSsupCOUNT(j)
					write (301, 11) CSsupELEM(j,k1)
					write (301, 16) CSsupAREA(j,k1), CSsupKONE(j,k1,1), &
							CSsupKONE(j,k1,2), CSsupKONE(j,k1,3), CSsupKONE(j,k1,4)
				end do
					
				! Amrazenamento do arquivo .nnn
				! Com o número do nó k1 do contorno sólido do objeto j
				! E as componentes do vetor normal ao nós
				do k1=1,ntn
					if (CSnnnNODE(j,k1) .ne. 0) then
						write (302, 17) CSnnnNODE(j,k1), CSnnnNORM(j,k1,1), &
										CSnnnNORM(j,k1,2), CSnnnNORM(j,k1,3)
					else
						exit
					end if
				end do
				
				close (301)
				close (302)
            end do
            bool = .FALSE.
        case default
            ! TODO - Adicionar error handling caso o programa encontre algo que nao espera
            bool = .FALSE.
        end select
    end do
    
	close (102)
	
	j = 0
	k1 = 0
	do i=1,ntn
		if (NEIBOR(i) .gt. j) then
			j = NEIBOR(i)
			k1 = i
		end if
	end do
	write(*, 200) j
	write(*, 200) k1
	
	!****************************************************
	! GERAR OS ARQUIVOS QUE NÃO SÃO UTILIZADOS
	!****************************************************
	
	ptname(7:14) = 'FF.sup  '
	open (205, file=ptname, status='unknown')
    write (205, 18) 0.0
	close (205)
	
	ptname(7:14) = 'FF.nnn  '
	open (206, file=ptname, status='unknown')
    write (206, 18) 0.0
	close (206)
	
	ptname(7:14) = 'CC.bt   '
	open (213, file=ptname, status='unknown')
	write (213, 18) 0.0
	close (213)
	
	ptname(7:14) = 'CC.brha '
	open (214, file=ptname, status='unknown')
	write (214, 18) 0.0
	close (214)
	
	ptname(7:14) = 'CE.fan  '
	open (215, file=ptname, status='unknown')
	write (215, 18) 0.0
	close (215)
	
	ptname(7:14) = 'CE.faq  '
	open (216, file=ptname, status='unknown')
	write (216, 18) 0.0
	close (216)
	
	ptname(7:14) = 'CE.raq  '
	open (217, file=ptname, status='unknown')
	write (217, 18) 0.0
	close (217)
    
end program

!***********************************************************
!			SUBROTINAS
!***********************************************************

subroutine AddNewCCNodes(ntn, ntot, nnew, master, new)

	implicit none
	
	integer :: i, j
	integer :: ntn      ! Número total de nós
	integer :: ntot		! Número total de nós ja existentes
	integer :: nnew		! Número de novos nós a serem analisados
	
	integer, dimension(ntn) :: master	! Vetor de nós ja existentes
	integer, dimension(nnew) :: new		! Vetor de novos nós
							
	! EXPLICAÇÃO: Para cada novo nó, a subrotina esta buscando nos nós existentes
	! 				no vetor 'master' se ha um nó com o mesmo número. Se houver,
	!				ele muda o número do nó para 0 no vetor 'new' e segue adiante.
	!				No final, para cada elemento em 'new' que nao for zero o algoritmo
	!				adiciona um novo elemento ao vetor 'master'. Isso garante que todos
	!				os nós sejam adicionados sem duplicatas.
	do i=1,nnew
		do j=1,ntot
			if (new(i) .eq. master(j)) then
				new(i) = 0
				exit
			end if
		end do
	end do
			
	do i=1,nnew
		if (new(i) .ne. 0) then
			ntot = ntot + 1
			master(ntot) = new(i)
		end if
	end do
	return
end subroutine AddNewCCNodes

subroutine UpdateUniformCCValues(ntn, master, quadnodes, ccvalue, values)
	
	implicit none
    
    !******************************
    ! Variáveis externas
    !******************************
    integer, intent(in) :: ntn			! Número total de nós
    integer, intent(in) :: master(ntn)  ! Vetor de nós com CC
    integer, intent(in) :: quadnodes(4)	! Vetor com os nós a serem atualizados
    
    real(8), intent(in) :: ccvalue		! Valor a ser inserido
    
    real(8) :: values(ntn)	! Vetor com os valores de CC a ser modificado
    
    !******************************
    ! Variáveis internas
    !******************************
    integer :: i,j
    
    do i=1,4
        do j=1,ntn
            if (quadnodes(i) .eq. master(j)) then
                values(j) = ccvalue
                exit
            end if
        end do
    end do

end subroutine UpdateUniformCCValues
    
subroutine UpdatePowerLawCCValues(ntn, master, node, z, prop, values)
	
	use functions
	implicit none
    
    !******************************
    ! Variáveis externas
    !******************************
    integer :: ntn			! Número total de nós
    integer, dimension(ntn) :: master  ! Vetor de nós com CC
    integer :: node			! Nó a ser atualizado
    
    real(8) :: z		! Vetor com a altura dos nós
    real(8), dimension(4) :: prop	! Propriedades Vref, zd, zref e p dos nós
    
    real(8), dimension(ntn) :: values	! Vetor com os valores de CC a ser modificado
        
    !******************************
    ! Variáveis internas
    !******************************
    integer :: i
    real(8) :: Vref, zd, zref, p
        
    do i=1,ntn
        if (node .eq. master(i)) then
            Vref = prop(1)
            zd = prop(2)
            zref = prop(3)
            p = prop(4)
            
            values(i) = PowerLaw(z, Vref, zd, zref, p)
            exit
        end if
    end do

end subroutine UpdatePowerLawCCValues

subroutine UpdateLogLawCCValues(ntn, master, node, z, prop, values)
	
	use functions
	implicit none
    
    !******************************
    ! Variáveis externas
    !******************************
    integer :: ntn			! Número total de nós
    integer, dimension(ntn) :: master  ! Vetor de nós com CC
    integer :: node			! Nó a ser atualizado
    
    real(8) :: z		! Vetor com a altura dos nós
    real(8), dimension(3) :: prop	! Propriedades Ustar, zd e z0 dos nós
    
    real(8), dimension(ntn) :: values	! Vetor com os valores de CC a ser modificado
        
    !******************************
    ! Variáveis internas
    !******************************
    integer :: i
    real(8) :: Ustar, zd, z0
        
    do i=1,ntn
        if (node .eq. master(i)) then
            Ustar = prop(1)
            zd = prop(2)
            z0 = prop(3)
            
            values(i) = LogLaw(z, Ustar, zd, z0)
            exit
        end if
    end do

end subroutine UpdateLogLawCCValues

subroutine GetQuadArea(area, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4)
	real(8) :: area, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
	real(8) :: v1a1x, v1a1y, v1a1z, v1a2x, v1a2y, v1a2z
	real(8) :: v2a1x, v2a1y, v2a1z, v2a2x, v2a2y, v2a2z
	real(8) :: AR1X, AR1Y, AR1Z, AR2X, AR2Y, AR2Z
	real(8) :: AREA1, AREA2
	
	! Código retirado do algoritmo PatransferIFE3D.for e adaptado
	v1a1x= x2 - x1
	v1a1y= y2 - y1
    v1a1z= z2 - z1
      
	v2a1x= x4 - x1
	v2a1y= y4 - y1
	v2a1z= z4 - z1

	v1a2x= x4 - x3
	v1a2y= y4 - y3
	v1a2z= z4 - z3

	v2a2x= x2 - x3
	v2a2y= y2 - y3
	v2a2z= z2 - z3
	
	AR1X= v1a1y * v2a1z - v1a1z * v2a1y
	AR1Y= v1a1z * v2a1x - v1a1x * v2a1z
	AR1Z= v1a1x * v2a1y - v1a1y * v2a1x

	AR2X= v1a2y * v2a2z - v1a2z * v2a2y
	AR2Y= v1a2z * v2a2x - v1a2x * v2a2z
	AR2Z= v1a2x * v2a2y - v1a2y * v2a2x
	
	AREA1= DSQRT ( AR1X*AR1X + AR1Y*AR1Y + AR1Z*AR1Z )
    AREA2= DSQRT ( AR2X*AR2X + AR2Y*AR2Y + AR2Z*AR2Z )
	 
    area = (AREA1/2.0D0) + (AREA2/2.0D0)
	return
end subroutine

subroutine GetNodeNormal(NNN1, NNN2, NNN3, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4)
	implicit none
	real(8) :: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
	real(8) :: V1_C1, V1_C2, V1_C3, V2_C1, V2_C2, V2_C3, V3_C1, V3_C2, V3_C3
	real(8) :: N1_C1, N1_C2, N1_C3, N2_C1, N2_C2, N2_C3
	real(8) :: NNN1, NNN2, NNN3, MODULL
	
	! Código retirado do algoritmo PatransferIFE3D.for e adaptado
	!write (*, *) 'COORD', x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
			! 8.5, -0.5, 0 ; 8.5, -0.5, 0.75 ; 9, -0.5, 0.75 ; 9, -0.5, 0
	
	V1_C1= x2 - x1	! 0
	V1_C2= y2 - y1	! 0
	V1_C3= z2 - z1	! 0.75

	V2_C1= x3 - x1	! 0.5
	V2_C2= y3 - y1	! 0
	V2_C3= z3 - z1	! 0.75

	V3_C1= x4 - x1	! 0.5
	V3_C2= y4 - y1  ! 0
	V3_C3= z4 - z1  ! 0

	N1_C1= V1_C2 * V2_C3 - V2_C2 * V1_C3	! 0
	N1_C2= V1_C3 * V2_C1 - V2_C3 * V1_C1	! 0.375
	N1_C3= V1_C1 * V2_C2 - V2_C1 * V1_C2    ! 0            

	N2_C1= V2_C2 * V3_C3 - V3_C2 * V2_C3	! 0
	N2_C2= V2_C3 * V3_C1 - V3_C3 * V2_C1	! 0.375
	N2_C3= V2_C1 * V3_C2 - V3_C1 * V2_C2	! 0

	NNN1= N1_C1 + N2_C1		! 0
	NNN2= N1_C2 + N2_C2		! 0.75
	NNN3= N1_C3 + N2_C3		! 0
	
	!write (*, *) 'NNN', NNN1, NNN2, NNN3

	MODULL=DSQRT( NNN1*NNN1 + NNN2*NNN2 + NNN3*NNN3 )
	
	!write (*, *) 'MODULL', MODULL

	NNN1= -NNN1/MODULL
	NNN2= -NNN2/MODULL
	NNN3= -NNN3/MODULL
	return
end subroutine

! Corrigir a direção da normal com base nas coordenadas do centroide do objeto
subroutine CorrectNormalDirection(NNN1, NNN2, NNN3, x, y, z, xobj, yobj, zobj)
	implicit none
	real(8) NNN1, NNN2, NNN3, x, y, z, xobj, yobj, zobj
	
	if ((xobj-x) * NNN1 .lt. 0.0) then
		NNN1 = -1 * NNN1
	end if
	
	if ((yobj-y) * NNN2 .lt. 0.0) then
		NNN2 = -1 * NNN2
	end if
	
	if ((zobj-z) * NNN3 .lt. 0.0) then
		NNN3 = -1 * NNN3
	end if
	
end subroutine CorrectNormalDirection
