program Post_VTK
    implicit none

    19 FORMAT(9I14)
    21 FORMAT(E14.5)
    23 FORMAT(3E14.5)

    INTEGER :: I
    INTEGER :: nnodes
    INTEGER :: nhexa

    INTEGER, ALLOCATABLE :: con(:,:)

    REAL(8), ALLOCATABLE :: xyz(:,:)
    REAL(8), ALLOCATABLE :: P(:)
    REAL(8), ALLOCATABLE :: uvw(:,:)

    CHARACTER(6) :: pname
    CHARACTER(3) :: iter
    CHARACTER(14) :: filename
    
    WRITE(*,*) "Nome do projeto em 6 caracteres"
    READ(*,*) pname

    WRITE(*,*) "Número do registro em 3 caracteres"
    READ(*,*) iter

    filename(1:6) = pname

    ! Leitura das coordenadas dos nós
    filename(7:14) = "000.cor "
    OPEN(1, file=filename, status="old")

    READ(1,*) nnodes

    ALLOCATE ( xyz(nnodes,3), P(nnodes), uvw(nnodes,3))

    do I = 1, nnodes
        READ(1,*) xyz(I,1), xyz(I,2), xyz(I,3)
    end do

    CLOSE(1)

    ! Leitura da conectividade da malha
    filename(7:14) = ".con    "
    OPEN(1, file=filename, status="old")

    READ(1,*) nhexa

    ALLOCATE ( con(nhexa,8) )

    do I = 1, nhexa
        READ(1,*) con(I,1), con(I,2), con(I,3), con(I,4), &
        con(I,5), con(I,6), con(I,7), con(I,8)
    end do

    CLOSE(1)

    ! Leitura da pressão
    filename(7:9) = iter
    filename(10:14) = ".pr  "
    OPEN(1, file=filename, status="old")

    do I = 1, nnodes
        READ(1,*) P(I)
    end do

    CLOSE(1)

    ! Leitura das velocidades
    filename(10:14) = ".v   "
    OPEN(1, file=filename, status="old")

    do I = 1, nnodes
        READ(1,*) uvw(I,1), uvw(I,2), uvw(I,3)
    end do

    CLOSE(1)

    ! Criação do arquivo .vtk
    filename(10:14) = ".vtk "
    OPEN(1, file=filename, status="unknown")

    WRITE(1,*) "# vtk DataFile Version 2.0"
    WRITE(1,*) "Unstructured Grid ", pname
    WRITE(1,*) "ASCII"
    WRITE(1,*) "DATASET UNSTRUCTURED_GRID"
    WRITE(1,*) ""

    WRITE(1,*) "POINTS", nnodes, "float"
    do I = 1, nnodes
        WRITE(1,23) xyz(I,1), xyz(I,2), xyz(I,3)
    end do
    WRITE(1,*) ""

    WRITE(1,*) "CELLS", nhexa, 9*nhexa
    do I = 1, nhexa
        WRITE(1,19) 8, con(I,1)-1, con(I,2)-1, con(I,3)-1, con(I,4)-1, &
        con(I,5)-1, con(I,6)-1, con(I,7)-1, con(I,8)-1
    end do
    WRITE(1,*) ""

    WRITE(1,*) "CELL_TYPES", nhexa
    do I = 1, nhexa
        WRITE(1,*) "12"
    end do
    WRITE(1,*) ""

    WRITE(1,*) "POINT_DATA", nnodes
    WRITE(1,*) "SCALARS Pressão float 1"
    WRITE(1,*) "LOOKUP_TABLE default"
    do I = 1, nnodes
        WRITE(1,21) p(I)
    end do
    WRITE(1,*) ""

    WRITE(1,*) "VECTORS Velocidade float"
    do I = 1, nnodes
        WRITE(1,23) uvw(I,1), uvw(I,2), uvw(I,3)
    end do
    WRITE(1,*) ""

    CLOSE(1)
    
end program Post_VTK