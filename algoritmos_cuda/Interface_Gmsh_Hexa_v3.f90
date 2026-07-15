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
	
	function LogLaw(z, Ustar, z0) result(V)

		implicit none
    
		real(8), intent(in) :: z, Ustar, z0
		real(8) :: V
    
		if (z .lt. 1.0d-6) then
			V = 0.0d0
		else
			V = (Ustar / 0.4d0) * log( (z + z0) / z0 )
		end if
    
	end function LogLaw

end module
    
program Interface_Gmsh_Hexa

    !**********************************************************************
    ! Programa de interface entre os dados da malha gerados pelo GMSH e o
    ! programa de CFD em MEF com elementos hexaedricos de 8 nós.
	! Adaptado para uso com o algoritmo CUDA_HEXA_IGOR_v1.
	! Inclui condições de entrada para turbulência.
    !
    ! Última modificacao: 11/12/2023
    !**********************************************************************
    ! Obs.: ______ refere-se a um nome de 6 dígitos dado ao projeto e aparece na maioria dos arquivos
    !
    ! Arquivos de Entrada:
    !   101 - ______.phys -> Inclui os dados das condições de contorno através dos grupos físicos
    !   102 - ______.msh -> Arquivo de exportacao da malha do GMSH em formato .msh versao 2.2
	!	103 - INTERFACE.dat -> Arquivo com as respostas das perguntas deste código, para quando não se deseja responder manualmente
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
	!	218 - ______.inflow		-> Condições da turbulência na entrada
	!	219 - ______.bvflut		-> Lista dos nós com turbulência
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
    
    logical :: bool, T, RHA, readFromFile	! boolean
    
    character (6) :: pname                  ! nome do projeto
    character (14) :: ptname                ! nome do arquivo
    character (200) :: header               ! header para leitura de arquivos
    character (200) :: garbc                ! caracter desnecessario
	
	! Variáveis da INTERFACE.dat
	integer 				:: NQUAD, NHEXA, NCOEF, NIR, NTR, NFILE, NPASS, INDTURB, inflowTurb, MNOBJ
	real(8) 				:: DtMAX, CSEGUR, VInf, VelSom, ViscCin, RHOInf
	integer, allocatable	:: NFCN(:), NNCN(:)
	real(8), allocatable	:: Lchar(:), Dchar(:), Xobj(:), Yobj(:), Zobj(:)
	
	integer	:: Ninflow, Nfseg, Nfran
	real(8) :: fmin, fmax, Uref, hrefU, alphaU, Dgama, hrefI, hrefL
	real(8) :: Cxyz(3), Iref(3), alphaI(3), Lref(3), alphaL(3)
    
    integer :: i, j, k1, k2, n, n1, n2, n3, n4, n5, n6, n7, n8, face, group, counter, iele   ! contadores e auxiliares inteiros
	integer :: count1, count2, count3, count4
    integer :: ntn, ntele     ! número de nós, quadrilateros, hexahedros e elementos
    integer :: nV1potnodes, nV2potnodes		! número de nós com CC por lei de potência
	integer :: nV1lognodes, nV2lognodes		! número de nós com CC por lei logaritmica
	integer :: nV1nodes, nV2nodes, nV3nodes, nPnodes	! número de nós com condicões de contorno
	integer :: MAXcountfIFE					! número máximo de faces em um objeto submerso
    integer :: garbi                        ! inteiro desnecessario
	integer :: NB1, NB2
	integer :: inflowTurbGroup
    

	integer :: quadnodes(4)					! vetor auxiliar com os nós de um elementos quadrilatero
	
    integer, allocatable :: nccV1potnodes(:), nccV2potnodes(:)	! vetores com os nós com CC por lei de potência
	integer, allocatable :: nccV1lognodes(:), nccV2lognodes(:)	! vetores com os nós com CC por lei logaritmica
	integer, allocatable :: nccV1nodes(:), nccV2nodes(:), nccV3nodes(:), nccPnodes(:)			! vetores com os nós com condicao de contorno
	integer, allocatable :: CSsupELEM(:,:)		! matriz com o número dos elementos a serem armazenados no arquivo CS##.sup
												! onde CSsupELEM(2,3) é o terceiro elemento do CS02.sup
	integer, allocatable :: CSsupKONE(:,:,:)	! matriz com a conectividade dos elementos quadriláteros a serem armazenados no arquivo CS##.sup
												! onde CSsupKONE(2,3,4) é o número do quarto nó, do terceiro elemento do CS02.sup
	integer, allocatable :: CSnnnNODE(:,:)		! matriz com o número dos nós a serem armazenados no arquivo CS##.nnn
												! onde CSnnnNODE(2,3) é o terceiro nó do CS02.nnn
    integer, allocatable :: CSsupCOUNT(:)		! vetor que armazena o número de elementos quadriláteros em cada corpo submerso
	integer, allocatable :: NEIBOR(:) 	! Verificação de quantos elementos tem em volta de cada nó
	integer, allocatable :: inflowNodes(:)
	
	real(8) :: garbr                           ! real desnecessario
	real(8) :: r, r2, r3, NNN1, NNN2, NNN3, ccvalue, V1, V2, V3, P
    
    real(8), allocatable :: vccV1potnodes(:), vccV2potnodes(:)		! vetores com os valores da CC por lei de potência em cada nó
	real(8), allocatable :: vccV1lognodes(:), vccV2lognodes(:)		! vetores com os valores da CC por lei logaritmica em cada nó
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
	
												
	!**********************************************************************
    ! DADOS INICIAIS
    !**********************************************************************
    write (*, 1) 'INFORMATIVO SOBRE O PROJETO' 
    write (*, 1) 'Digite o nome do projeto (pname) - em seis caracteres:'
    read (*, *) pname
	ptname (1:6) = pname
	
	open (103, file='INTERFACE.dat', status='old')

	read(103,*) NQUAD
	read(103,*) NHEXA
    
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
						valueLogCC(countCCTypes(4),2), indexCC(CountPhysGroups,5) )
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
				read (101,*) valueLogCC(i,1), valueLogCC(i,2)
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
		case ('$InflowTurbGroup')
			! Armazenar o número do grupo físico com turbulência na entrada
			read (101,*) inflowTurbGroup
			
			read (101,*) garbc	! $EndInflowTurbGroup
		case ('$End')
			bool = .FALSE.
        end select
    end do
	
	close(101)
    
	!****************************************************************
	! CRIAÇÃO DO ARQUIVO DE PARÂMETROS DE INICIALIZAÇÃO INICIA.PAR
	!****************************************************************
	
	read (103, *) NCOEF
	read (103, *) NIR
	read (103, *) NTR
	read (103, *) NFILE
	read (103, *) DtMAX
	read (103, *) CSEGUR
	read (103, *) NPASS
	read (103, *) INDTURB
	read (103, *) inflowTurb
	read (103, *) VInf
	read (103, *) VelSom
	read (103, *) ViscCin
	read (103, *) RHOInf
	read (103, *) MNOBJ
	
	allocate( Lchar(MNOBJ), Dchar(MNOBJ), Xobj(MNOBJ), Yobj(MNOBJ), Zobj(MNOBJ), NFCN(MNOBJ), NNCN(MNOBJ) )
	
	do i=1,MNOBJ
		read (103, *) Lchar(i)
		read (103, *) Dchar(i)
		read (103, *) Xobj(i)
		read (103, *) Yobj(i)
		read (103, *) Zobj(i)
		read (103, *) NFCN(i)
		read (103, *) NNCN(i)
	end do
	
	open (201, file='inicia.par', status='unknown')
	
	write (201, 6) pname
	write (201, 7) NCOEF
	write (201, 7) NIR
	write (201, 7) NTR
	write (201, 7) NFILE
	write (201, 8) DtMAX
	write (201, 8) 0.0 !TPOAC
	write (201, 8) 0.000001 !TOLTPO
	write (201, 7) 100 !NROTPO
	write (201, 8) CSEGUR
	write (201, 8) 1.0 !CONTROL1
	write (201, 7) NPASS
	write (201, 7) INDTURB
	write (201, 7) inflowTurb
	write (201, 8) 0.9 !ELUMP1
	write (201, 8) 0.9 !ELUMP2
	write (201, 8) 0.9 !ELUMP3
	write (201, 7) 128 !tBlock
	
	close (201)
    
	!****************************************************************
	! CRIAÇÃO DO ARQUIVO DE PROPRIEDADES DO FLUIDO .PRO
	!****************************************************************

	ptname (7:14) = '.pro    '
	open (202, file=ptname, status='unknown')
	
	write (202, 8) VInf
	write (202, 8) VelSom
	write (202, 8) ViscCin
	write (202, 8) 0.0 !ViscVol
	write (202, 8) RHOInf
	write (202, 8) 0.00000001 !Cv
	write (202, 8) 0.00000001 !Kdif
	write (202, 8) 0.15 !Cs
	write (202, 8) 0.00000001 !Prtlt
	write (202, 8) 0.00000001 !Dab
		
	do i=1,MNOBJ
		write (202, 8) Lchar(i)
		write (202, 8) Dchar(i)
		write (202, 8) Xobj(i)
		write (202, 8) Yobj(i)
		write (202, 8) Zobj(i)
		write (202, 7) NFCN(i)
		write (202, 7) NNCN(i)
	end do
	
	close (202)
	
	!****************************************************************
	! CRIAÇÃO DO ARQUIVO DE TURBULÊNCIA NA ENTRADA .inflow
	!****************************************************************
	
	if (inflowTurb .eq. 1) then
		ptname (7:14) = '.inflow '
		open (218, file=ptname, status='unknown')
		
		read (103, *) Ninflow, Nfseg, Nfran
		read (103, *) fmin, fmax
		read (103, *) Uref, hrefU, alphaU
		read (103, *) Dgama
		read (103, *) Cxyz(1), Cxyz(2), Cxyz(3)
		read (103, *) Iref(1), Iref(2), Iref(3)
		read (103, *) hrefI
		read (103, *) alphaI(1), alphaI(2), alphaI(3)
		read (103, *) Lref(1), Lref(2), Lref(3)
		read (103, *) hrefL
		read (103, *) alphaL(1), alphaL(2), alphaL(3)
		
		write (218, *) Ninflow, Nfseg, Nfran
		write (218, *) fmin, fmax
		write (218, *) Uref, hrefU, alphaU
		write (218, *) Dgama
		write (218, *) Cxyz(1), Cxyz(2), Cxyz(3)
		write (218, *) Iref(1), Iref(2), Iref(3)
		write (218, *) hrefI
		write (218, *) alphaI(1), alphaI(2), alphaI(3)
		write (218, *) Lref(1), Lref(2), Lref(3)
		write (218, *) hrefL
		write (218, *) alphaL(1), alphaL(2), alphaL(3)
		
		close (218)
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
					
					V1 = PowerLaw(zz(i), valuePowerCC(j,1), valuePowerCC(j,2), valuePowerCC(j,3), valuePowerCC(j,4))
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
				case (3)
					! Campo em Lei de Potência
					do j=1,countCCTypes(3)
						if (indexPowerCC(j,1) .eq. initialFieldGroup .and. indexPowerCC(j,2) .eq. 1) then
							exit
						end if
					end do
					
					V3 = PowerLaw(zz(i), valuePowerCC(j,1), valuePowerCC(j,2), valuePowerCC(j,3), valuePowerCC(j,4))
				case (4)
					! Campo em Lei Logaritmica
					! TODO - Implementar
					V3 = -10000.0
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
                
            end do
            ! TODO - É possível implementar uma verificacao para ver se essa linha e de fato '$EndNodes'
            read (102, *) garbc
            
            close (204)
			close (207)
			close (208)
            
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
			allocate( nccV1potnodes(ntn), nccV2potnodes(ntn), &
					vccV1potnodes(ntn), vccV2potnodes(ntn), &
					nccV1lognodes(ntn), nccV2lognodes(ntn), &
					vccV1lognodes(ntn), vccV2lognodes(ntn), &
					nccV1nodes(ntn), nccV2nodes(ntn), nccV3nodes(ntn), nccPnodes(ntn), &
					vccV1nodes(ntn), vccV2nodes(ntn), vccV3nodes(ntn), vccPnodes(ntn), &
					CSsupELEM(MNOBJ,nquad), CSsupKONE(MNOBJ,nquad,4), &
					CSsupAREA(MNOBJ,nquad), CSsupCOUNT(MNOBJ), &
					CSnnnNODE(MNOBJ,ntn), CSnnnNORM(MNOBJ,ntn,3), inflowNodes(ntn) )
            nccV1potnodes = 0
            nccV2potnodes = 0
            vccV1potnodes = 0
            vccV2potnodes = 0
			nccV1nodes = 0
			nccV2nodes = 0
			nccV3nodes = 0
			nccPnodes = 0
			vccV1nodes = 0
			vccV2nodes = 0
			vccV3nodes = 0
			vccPnodes = 0
            nV1potnodes = 0
            nV2potnodes = 0
			nV1lognodes = 0
			nV2lognodes = 0
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
			inflowNodes = 0
			
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
					
					! TODO - Refatorar
                    call UpdatePowerLawCCValues(ntn, nccV1potnodes, n1, zz(n1), valuePowerCC(j,:), vccV1potnodes)
                    call UpdatePowerLawCCValues(ntn, nccV1potnodes, n2, zz(n2), valuePowerCC(j,:), vccV1potnodes)
                    call UpdatePowerLawCCValues(ntn, nccV1potnodes, n3, zz(n3), valuePowerCC(j,:), vccV1potnodes)
                    call UpdatePowerLawCCValues(ntn, nccV1potnodes, n4, zz(n4), valuePowerCC(j,:), vccV1potnodes)
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
                    call UpdateLogLawCCValues(ntn, nccV1lognodes, n1, zz(n1), valueLogCC(j,:), vccV1lognodes)
                    call UpdateLogLawCCValues(ntn, nccV1lognodes, n2, zz(n2), valueLogCC(j,:), vccV1lognodes)
                    call UpdateLogLawCCValues(ntn, nccV1lognodes, n3, zz(n3), valueLogCC(j,:), vccV1lognodes)
                    call UpdateLogLawCCValues(ntn, nccV1lognodes, n4, zz(n4), valueLogCC(j,:), vccV1lognodes)
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
                    
					call AddNewCCNodes(ntn, nV2potnodes, 4, nccV2potnodes, quadnodes)
					
					! TODO - Refatorar
                    call UpdatePowerLawCCValues(ntn, nccV2potnodes, n1, zz(n1), valuePowerCC(j,:), vccV2potnodes)
                    call UpdatePowerLawCCValues(ntn, nccV2potnodes, n2, zz(n2), valuePowerCC(j,:), vccV2potnodes)
                    call UpdatePowerLawCCValues(ntn, nccV2potnodes, n3, zz(n3), valuePowerCC(j,:), vccV2potnodes)
                    call UpdatePowerLawCCValues(ntn, nccV2potnodes, n4, zz(n4), valuePowerCC(j,:), vccV2potnodes)
				case (4)
					! CC Log
					
					! Verificar se os nós ja estao na lista
					quadnodes(1) = n1
					quadnodes(2) = n2
					quadnodes(3) = n3
					quadnodes(4) = n4
					
					do j=1,countCCTypes(4)
						if (indexLogCC(j,1) .eq. group .and. indexLogCC(j,2) .eq. 1) then
							exit
						end if
					end do
                    
					call AddNewCCNodes(ntn, nV2lognodes, 4, nccV2lognodes, quadnodes)
                    call UpdateLogLawCCValues(ntn, nccV2lognodes, n1, zz(n1), valueLogCC(j,:), vccV2lognodes)
                    call UpdateLogLawCCValues(ntn, nccV2lognodes, n2, zz(n2), valueLogCC(j,:), vccV2lognodes)
                    call UpdateLogLawCCValues(ntn, nccV2lognodes, n3, zz(n3), valueLogCC(j,:), vccV2lognodes)
                    call UpdateLogLawCCValues(ntn, nccV2lognodes, n4, zz(n4), valueLogCC(j,:), vccV2lognodes)
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
								call CorrectNormalDirection(NNN1, NNN2, NNN3, xx(n1), yy(n1), zz(n1), Xobj(j), Yobj(j), Zobj(j))
								CSnnnNORM(j,k2,1) = CSnnnNORM(j,k2,1) + NNN1
								CSnnnNORM(j,k2,2) = CSnnnNORM(j,k2,2) + NNN2
								CSnnnNORM(j,k2,3) = CSnnnNORM(j,k2,3) + NNN3
							else if (CSnnnNODE(j,k2) .eq. n2) then
								call GetNodeNormal(NNN1, NNN2, NNN3, &
											xx(n2), xx(n3), xx(n4), xx(n1), &
											yy(n2), yy(n3), yy(n4), yy(n1), &
											zz(n2), zz(n3), zz(n4), zz(n1) )
								! Corrigir a direção da normal com base nas coordenadas do centroide do objeto
								call CorrectNormalDirection(NNN1, NNN2, NNN3, xx(n2), yy(n2), zz(n2), Xobj(j), Yobj(j), Zobj(j))
								CSnnnNORM(j,k2,1) = CSnnnNORM(j,k2,1) + NNN1
								CSnnnNORM(j,k2,2) = CSnnnNORM(j,k2,2) + NNN2
								CSnnnNORM(j,k2,3) = CSnnnNORM(j,k2,3) + NNN3
							else if (CSnnnNODE(j,k2) .eq. n3) then
								call GetNodeNormal(NNN1, NNN2, NNN3, &
											xx(n3), xx(n4), xx(n1), xx(n2), &
											yy(n3), yy(n4), yy(n1), yy(n2), &
											zz(n3), zz(n4), zz(n1), zz(n2) )
								! Corrigir a direção da normal com base nas coordenadas do centroide do objeto
								call CorrectNormalDirection(NNN1, NNN2, NNN3, xx(n3), yy(n3), zz(n3), Xobj(j), Yobj(j), Zobj(j))
								CSnnnNORM(j,k2,1) = CSnnnNORM(j,k2,1) + NNN1
								CSnnnNORM(j,k2,2) = CSnnnNORM(j,k2,2) + NNN2
								CSnnnNORM(j,k2,3) = CSnnnNORM(j,k2,3) + NNN3
							else if (CSnnnNODE(j,k2) .eq. n4) then
								call GetNodeNormal(NNN1, NNN2, NNN3, &
											xx(n4), xx(n1), xx(n2), xx(n3), &
											yy(n4), yy(n1), yy(n2), yy(n3), &
											zz(n4), zz(n1), zz(n2), zz(n3) )
								! Corrigir a direção da normal com base nas coordenadas do centroide do objeto
								call CorrectNormalDirection(NNN1, NNN2, NNN3, xx(n4), yy(n4), zz(n4), Xobj(j), Yobj(j), Zobj(j))
								CSnnnNORM(j,k2,1) = CSnnnNORM(j,k2,1) + NNN1
								CSnnnNORM(j,k2,2) = CSnnnNORM(j,k2,2) + NNN2
								CSnnnNORM(j,k2,3) = CSnnnNORM(j,k2,3) + NNN3
							else if (CSnnnNODE(j,k2) .eq. 0) then
								exit
							end if
						end do
					end if
				end do
				
				! Turbulência na Entrada
				if (inflowTurb .eq. 1) then
					if (group.eq.inflowTurbGroup) then
						! Verificar se os nós já estão na lista
						quadnodes(1) = n1
						quadnodes(2) = n2
						quadnodes(3) = n3
						quadnodes(4) = n4
						count1 = 0
						do j=1,Ninflow
							if (inflowNodes(j).eq.0) then
								count1 = j
								exit
							end if
							
							do k1=1,4
								if (inflowNodes(j).eq.quadnodes(k1)) then
									quadnodes(k1) = 0
								end if
							end do
							
						end do
						
						! Adicionar nós restantes
						do j=1,4
							if (quadnodes(j).ne.0) then
								inflowNodes(count1) = quadnodes(j)
								count1 = count1 + 1
							end if
						end do
					
					end if
				end if
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
			
			write (211, *) nV2nodes + nV2potnodes + nV2lognodes
			do i=1,nV2nodes
				write (211, 15) nccV2nodes(i), vccV2nodes(i)
			end do
			do i=1,nV2potnodes
				write (211, 15) nccV2potnodes(i), vccV2potnodes(i)
			end do
			do i=1,nV2lognodes
				write (211, 15) nccV2lognodes(i), vccV2lognodes(i)
			end do
			
			write (211, *) nV3nodes
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
			
			! Armazenamento dos nós com turbulência
			ptname (7:14) = '.bvflut '
			open (219, file=ptname, status='unknown')
			
			do i=1,Ninflow
				write (219, 11) inflowNodes(i)
			end do
			
			close (219)
			
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

	do i=1, ntn
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
	
	! ptname(7:14) = 'FF.sup  '
	! open (205, file=ptname, status='unknown')
    ! write (205, 18) 0.0
	! close (205)
	
	! ptname(7:14) = 'FF.nnn  '
	! open (206, file=ptname, status='unknown')
    ! write (206, 18) 0.0
	! close (206)
	
	! ptname(7:14) = 'CC.bt   '
	! open (213, file=ptname, status='unknown')
	! write (213, 18) 0.0
	! close (213)
	
	! ptname(7:14) = 'CC.brha '
	! open (214, file=ptname, status='unknown')
	! write (214, 18) 0.0
	! close (214)
	
	! ptname(7:14) = 'CE.fan  '
	! open (215, file=ptname, status='unknown')
	! write (215, 18) 0.0
	! close (215)
	
	! ptname(7:14) = 'CE.faq  '
	! open (216, file=ptname, status='unknown')
	! write (216, 18) 0.0
	! close (216)
	
	! ptname(7:14) = 'CE.raq  '
	! open (217, file=ptname, status='unknown')
	! write (217, 18) 0.0
	! close (217)
    
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
    
    real(8) :: z		! Altura do nó
    real(8), dimension(2) :: prop	! Propriedades Ustar e z0 dos nós
    
    real(8), dimension(ntn) :: values	! Vetor com os valores de CC a ser modificado
        
    !******************************
    ! Variáveis internas
    !******************************
    integer :: i
    real(8) :: Ustar, z0
        
    do i=1,ntn
        if (node .eq. master(i)) then
            Ustar = prop(1)
            z0 = prop(2)
            
            values(i) = LogLaw(z, Ustar, z0)
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
	real(8) :: NNN1, NNN2, NNN3, x, y, z, xobj, yobj, zobj
	real(8) :: norm_vector(3), centroid_vector(3)
	
	norm_vector = (/ NNN1, NNN2, NNN3 /)
	centroid_vector = (/ xobj-x, yobj-y, zobj-z /)
	
	if (DOT_PRODUCT(norm_vector, centroid_vector) .lt. 0.0) then
		NNN1 = -1 * NNN1
		NNN2 = -1 * NNN2
		NNN3 = -1 * NNN3
	end if
	
end subroutine CorrectNormalDirection
