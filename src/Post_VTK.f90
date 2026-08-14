program Post_VTK
    implicit none

    19 FORMAT(9I14)
    21 FORMAT(E14.5)
    23 FORMAT(3E14.5)

    integer :: I
    integer :: nnodes
    integer :: nhexa

    integer, allocatable :: con(:,:)

    REAL(8), allocatable :: xyz(:,:)
    REAL(8), allocatable :: P(:)
    REAL(8), allocatable :: uvw(:,:)

    character(6) :: pname
    character(3) :: iter
    character(3) :: file_num
    character(14) :: filename
    
    write(*,*) "Nome do projeto em 6 caracteres"
    read(*,*) pname

    write(*,*) "Número do registro inicial em 3 caracteres"
    read(*,*) iter

    write(*,*) "Número de registros, não preencher = 1"
    read(int,*) file_num

    filename(1:6) = pname

    ! Leitura das coordenadas dos nós
    filename(7:14) = "000.cor "
    open(1, file=filename, status="old")

    read(1,*) nnodes

    allocate ( xyz(nnodes,3), P(nnodes), uvw(nnodes,3))

    do I = 1, nnodes
        read(1,*) xyz(I,1), xyz(I,2), xyz(I,3)
    end do

        close(1)

        ! Leitura da conectividade da malha
        filename(7:14) = ".con    "
        open(1, file=filename, status="old")

        read(1,*) nhexa

            allocate ( con(nhexa,8) )

            do I = 1, nhexa
                read(1,*) con(I,1), con(I,2), con(I,3), con(I,4), &
                con(I,5), con(I,6), con(I,7), con(I,8)
            end do

            close(1)

            ! Leitura da pressão
            filename(7:9) = iter
            filename(10:14) = ".pr  "
            open(1, file=filename, status="old")

            do I = 1, nnodes
                read(1,*) P(I)
            end do

            close(1)

            ! Leitura das velocidades
            filename(10:14) = ".v   "
            open(1, file=filename, status="old")

            do I = 1, nnodes
                read(1,*) uvw(I,1), uvw(I,2), uvw(I,3)
            end do

            close(1)

            ! Criação do arquivo .vtk
            filename(10:14) = ".vtk "
            open(1, file=filename, status="unknown")

        write(1,'(A)') "# vtk DataFile Version 2.0"
        write(1,*) "Unstructured Grid ", pname
        write(1,*) "ASCII"
        write(1,*) "DATASET UNSTRUCTURED_GRID"
        write(1,*) ""

        write(1,*) "POINTS", nnodes, "float"
        do I = 1, nnodes
            write(1,23) xyz(I,1), xyz(I,2), xyz(I,3)
        end do
        write(1,*) ""

        write(1,*) "CELLS", nhexa, 9*nhexa
        do I = 1, nhexa
            write(1,19) 8, con(I,1)-1, con(I,2)-1, con(I,3)-1, con(I,4)-1, &
            con(I,5)-1, con(I,6)-1, con(I,7)-1, con(I,8)-1
        end do
        write(1,*) ""

        write(1,*) "CELL_TYPES", nhexa
        do I = 1, nhexa
            write(1,*) "12"
        end do
        write(1,*) ""

        write(1,*) "POINT_DATA", nnodes
        write(1,*) "SCALARS Pressão float 1"
        write(1,*) "LOOKUP_TABLE default"
        do I = 1, nnodes
            write(1,21) p(I)
        end do
        write(1,*) ""

        write(1,*) "VECTORS Velocidade float"
        do I = 1, nnodes
            write(1,23) uvw(I,1), uvw(I,2), uvw(I,3)
        end do
        write(1,*) ""

        close(1)
    
end program Post_VTK
