subroutine lecture_donnee(p,m)
        ! Lis les données d'entrée du problème dans le fichier donnee.dat
        ! Pour les mettre dans p (struc phys) et m (struc maillage).
        use m_type
        implicit none

        type(maillage), intent(out) :: m
        type(phys), intent(out) :: p

        open(10,file='donnee.dat')
        
        read(10,*) p%L
        read(10,*) p%t_tot
        read(10,*) p%D
        read(10,*) p%C0
        read(10,*) p%C1
        read(10,*) p%alph
        read(10,*) m%nx
        read(10,*) m%ny
        read(10,*) p%cfl
        read(10,*) p%r
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

        m%dx = p%L/real(m%Nx-1)
        do i=1,m%Nx
                m%xn(i,:) = real(i-1)*m%dx
        end do

        do i=1,m%Ny
                y_reg = real(i-1)*p%L/real(m%Ny-1)
                m%yn(:,i) = y_reg - p%L/(3.*PI)*sin((2.*PI*y_reg)/p%L)
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
                dyn = m%yn(i,j+1) - m%yn(i,j)
                m%u(m%nx,j)=p%alph*sin(PI*m%xn(m%nx,j)/p%L)*cos(PI*(m%yn(m%nx,j) + dyn/2.)/p%L)
        end do
        do i=1,m%nx-1
                m%v(i,m%ny)=-p%alph*cos(PI*(m%xn(i,m%ny)+m%dx/2.)/p%L)*sin(PI*m%yn(i,m%ny)/p%L)
        end do
end subroutine vitesse

subroutine concentration(p, m, c)
        use m_type
        implicit none

        type(maillage), intent(in) :: m
        type(phys), intent(in) :: p
        type(conc), intent(inout) :: c
        integer :: i,j
        real :: vol, adv, diff
        
        call cal_Fe_adv(c, m)
        call cal_Fo_adv(c, m)
        call cal_Fn_adv(c, m, p)
        call cal_Fs_adv(c, m, p)
        
        call cal_Fe_diff(c, m)
        call cal_Fo_diff(c, m)
        call cal_Fn_diff(c, m, p)
        call cal_Fs_diff(c, m, p)
        
        do i = 1,m%nx-1
                do j = 1,m%ny-1
                        vol = m%dx*(m%yn(i,j+1) - m%yn(i,j))
                        adv = c%Fo_adv(i,j) + c%Fe_adv(i,j) + c%Fs_adv(i,j) + c%Fn_adv(i,j)
                        diff = c%Fo_diff(i,j) + c%Fe_diff(i,j) + c%Fs_diff(i,j) + c%Fn_diff(i,j)
                        c%mat_c(i,j) = c%mat_c(i,j) - m%dt/vol*(adv - p%D*diff)
                end do
        end do
        deallocate(c%Fe_adv,c%Fo_adv,c%Fn_adv,c%Fs_adv)
        deallocate(c%Fe_diff,c%Fo_diff,c%Fn_diff,c%Fs_diff)
end subroutine concentration

subroutine pdt(p,m)
	use m_type
	implicit none 

	type(phys), intent(in) :: p
	type(maillage), intent(inout) :: m
	real, dimension(m%nx-1,m%ny-1) :: T
	real :: dyn
	integer :: i,j
	
	do i = 1,m%nx-1
	    do j = 1,m%ny-1
		dyn = m%yn(i,j+1) - m%yn(i,j)
		T(i,j) = abs(m%u(i,j))/(p%CFL*m%dx) + abs(m%v(i,j))/(p%CFL*dyn) &
		+ p%D/(p%R*(m%dx**2.))+p%D/(p%R*(dyn**2.))
	    end do
	end do
	m%dt = minval(1./T)
	m%nt = int(p%t_tot/m%dt)
	write(*,*)
	write(*,*) "Temps total (s)", p%t_tot
	write(*,*) "dt (s)", m%dt
	write(*,*) "Nombre points", m%nt
	write(*,*)
end subroutine pdt
