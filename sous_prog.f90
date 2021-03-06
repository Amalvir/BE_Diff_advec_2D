subroutine lecture_donnee(p,m)
        ! Lis les données d'entrée du problème dans le fichier donnee.dat
        ! Pour les mettre dans p (struc phys) et m (struc maillage).
        use m_type
        implicit none

        type(maillage), intent(out) :: m
        type(phys), intent(out) :: p

        write(*,*) "[I] Lecture des données."
        open(10,file='donnee.dat')
        
        read(10,*) p%L
        read(10,*) p%D
        read(10,*) p%C0
        read(10,*) p%C1
        read(10,*) p%alph
        read(10,*) m%Nx
        read(10,*) m%Ny
        close(10)

end subroutine lecture_donnee

subroutine def_maillage(p,m)
        ! Calcul les coordonnées des noeuds du maillage m grâce aux données de p (stuc phys) et m. Maillage régulier selon x mais irrégulier selon y.
        use m_type
        implicit none

        type(maillage), intent(inout) :: m
        type(phys), intent(in) :: p
        integer :: i
        real :: y_reg, PI

        PI = acos(-1.)

        allocate(m%xn(m%Nx, m%Ny))
        allocate(m%yn(m%Nx, m%Ny))
        !allocate(m%dyn(m%Ny-1))

        m%dx = p%L/real(m%Nx-1)
        write(*,*) "[I] Maillage régulier selon x"
        do i=1,m%Nx
                m%xn(i,:) = (i-1)*m%dx
        end do

        write(*,*) "[I] Maillage irrégulier selon y"
        do i=1,m%Ny
                y_reg = real(i-1)*p%L/real(m%Ny-1)
                m%yn(:,i) = y_reg - p%L/(3.*PI)*sin((2.*PI*y_reg)/p%L)
                !m%dyn(i-1) = m%yn(1,i) - m%yn(1,i-1)
                !m%dyv(i) = 0.5*(m%dyn(i-1) + m%dyn(i))
        end do

end subroutine def_maillage

subroutine vitesse(p,m)
        use m_type
        implicit none

        type(maillage), intent(inout) :: m
        type(phys), intent(in) :: p
        integer :: i,j
        real :: dyn,PI

        PI = acos(-1.)

        allocate(m%u(m%nx,m%ny-1), m%v(m%nx-1,m%ny))

        do i = 1,m%nx-1
                do j = 1,m%ny-1
                dyn = m%yn(i,j+1) - m%yn(i,j)
                m%u(i,j)=p%alph*sin(PI*m%xn(i,j)/p%L)*cos(PI*(m%yn(i,j) + dyn/2.)/p%L) ! beta = 0
                m%v(i,j)=-p%alph*cos(PI*(m%xn(i,j)+m%dx/2.)/p%L)*sin(PI*m%yn(i,j)/p%L)
                end do
        end do
        do j = 1,m%ny-1
                dyn = m%yn(m%nx,j+1) - m%yn(m%nx,j)
                m%u(m%nx,j)=p%alph*sin(PI*m%xn(m%nx,j)/p%L)*cos(PI*(m%yn(m%nx,j) + dyn/2.)/p%L)
        end do
        do i=1,m%nx-1
                m%v(i,m%ny)=-p%alph*cos(PI*(m%xn(i,m%ny)+m%dx/2.)/p%L)*sin(PI*m%yn(i,m%ny)/p%L)
        end do

end subroutine vitesse

