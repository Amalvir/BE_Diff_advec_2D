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
        read(10,*) p%t_tot
        read(10,*) p%D
        read(10,*) p%C0
        read(10,*) p%C1
        read(10,*) p%alph
        read(10,*) m%nx
        read(10,*) m%ny
        read(10,*) m%nt
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
                dyn = m%yn(i,j+1) - m%yn(i,j)
                m%u(m%nx,j)=p%alph*sin(PI*m%xn(m%nx,j)/p%L)*cos(PI*(m%yn(m%nx,j) + dyn/2.)/p%L)
        end do
        do i=1,m%nx-1
                m%v(i,m%ny)=-p%alph*cos(PI*(m%xn(i,m%ny)+m%dx/2.)/p%L)*sin(PI*m%yn(i,m%ny)/p%L)
        end do

end subroutine vitesse

subroutine concentration(p, m, c, conc_n)
        use m_type
        implicit none

        type(maillage), intent(in) :: m
        type(phys), intent(in) :: p
        type(conc), intent(inout) :: c
        real, dimension(m%nx-1,M%ny-1), intent(in) :: conc_n
        integer :: i,j
        
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
			c%mat_c(i,j) = conc_n(i,j) - m%dt/m%v(i,j)*(c%Fo_adv(i,j) + c%Fe_adv(i,j) + c%Fs_adv(i,j) + c%Fn_adv(i,j) &
			- p%D*(c%Fo_diff(i,j) + c%Fe_diff(i,j) + c%Fs_diff(i,j) + c%Fn_diff(i,j)))
		end do
	end do
	
	deallocate(c%Fe_adv,c%Fo_adv,c%Fn_adv,c%Fs_adv)
	deallocate(c%Fe_diff,c%Fo_diff,c%Fn_diff,c%Fs_diff)
end subroutine concentration

! Subroutine pour le flux Est advectif, faut faire celle pour le reste...

subroutine cal_Fe_adv(c, m)
        use m_type
        implicit none
        
        type(conc), intent(inout) :: c
        type(maillage), intent(in) :: m
        integer :: i,j
        real :: Se

        allocate(c%Fe_adv(m%nx-1,m%ny-1))

        do i = 1,m%nx-2
                do j = 1,m%ny-2
                        Se = m%yn(i,j+1) - m%yn(i,j)
                        if (m%u(i,j) >= 0.) then
                                c%Fe_adv(i,j) = c%mat_c(i,j)*m%u(i,j)*Se
                        else
                                c%Fe_adv(i,j) = c%mat_c(i+1,j)*m%u(i,j)*Se
                        end if
                end do
        end do 

        ! Conditions aux limites
        
        do j = 1,m%ny-1
                Se = m%yn(m%nx-1,j+1) - m%yn(m%nx-1,j) ! Je sais pas par quoi remplacer i
                if (m%u(m%nx-1,j) >= 0.) then
                        c%Fe_adv(m%nx-1,j) = c%mat_c(m%nx-1,j)*m%u(m%nx-1,j)*Se
                else
                        c%Fe_adv(m%nx-1,j) = 0.
                end if
        end do
end subroutine cal_Fe_adv

subroutine cal_Fo_adv(c,m)
    use m_type 
    implicit none

    type(conc), intent(inout) :: c
        type(maillage), intent(in) :: m
        integer :: i,j
        real :: So

        allocate(c%Fo_adv(m%nx-1,m%ny-1))

        do i = 2,m%nx-1
                do j = 2,m%ny-1
                        So = m%yn(i,j+1) - m%yn(i,j)
                        if (m%u(i,j) >= 0.) then
                                c%Fo_adv(i,j) = -c%mat_c(i-1,j)*m%u(i,j)*So
                        else
                                c%Fo_adv(i,j) = -c%mat_c(i,j)*m%u(i,j)*So
                        end if
                end do
        end do 

        ! Conditions aux limites
        
        do j = 1,m%ny-1
                So = m%yn(1,j+1) - m%yn(1,j)
                if (m%u(i,j) >= 0.) then 
                        c%Fo_adv(1,j) = 0.
                else
                        c%Fo_adv(1,j) = -c%mat_c(1,j)*m%u(1,j)*So
                end if
        end do
end subroutine cal_Fo_adv

subroutine cal_Fn_adv(c,m,p)
    use m_type 
    implicit none

        type(conc), intent(inout) :: c
        type(maillage), intent(in) :: m
        type(phys), intent(in) :: p
        integer :: i,j
        real :: Sn

        allocate(c%Fn_adv(m%nx-1,m%ny-1))
        Sn = m%dx
        do i = 1,m%nx-2
                do j = 1,m%ny-2
                        if (m%v(i,j) >= 0.) then
                                c%Fn_adv(i,j) = c%mat_c(i,j)*m%v(i,j)*Sn
                        else
                                c%Fn_adv(i,j) = c%mat_c(i,j+1)*m%v(i,j)*Sn
                        end if
                end do
        end do 

        ! Conditions aux limites
        
        do i = 1,m%nx-1
                if (m%v(i,j) >= 0.) then 
                        c%Fn_adv(i,m%ny-1) = c%mat_c(i,m%ny-1)*m%v(i,m%ny-1)*Sn
                else
                        c%Fn_adv(i,m%ny-1) = p%C0*m%v(i,m%ny-1)*Sn
                end if
        end do
end subroutine cal_Fn_adv

subroutine cal_Fs_adv(c,m,p)
    use m_type 
    implicit none

        type(conc), intent(inout) :: c
        type(maillage), intent(in) :: m
        type(phys), intent(in) :: p
        integer :: i,j
        real :: Ss

        allocate(c%Fs_adv(m%nx-1,m%ny-1))
        Ss = m%dx
        do i = 2,m%nx-1
                do j = 2,m%ny-1
                        if (m%v(i,j) >= 0.) then
                                c%Fs_adv(i,j) = -c%mat_c(i,j-1)*m%v(i,j)*Ss
                        else
                                c%Fs_adv(i,j) = -c%mat_c(i,j)*m%v(i,j)*Ss
                        end if
                end do
        end do 

        ! Conditions aux limites
        
        do i = 1,m%nx-1
                if (m%v(i,j) >= 0.) then 
                        c%Fs_adv(i,1) = -p%C1*m%v(i,1)*Ss
                else
                        c%Fs_adv(i,1) = -c%mat_c(i,1)*m%v(i,1)*Ss
                end if
        end do
end subroutine cal_Fs_adv

! Flux diffusif

subroutine cal_Fe_diff(c,m)
    use m_type 
    implicit none

    type(conc), intent(inout) :: c
        type(maillage), intent(in) :: m
        integer :: i,j
        real :: Se

        allocate(c%Fe_diff(m%nx-1,m%ny-1))

        do i = 1,m%nx-2
                do j = 1,m%ny-1
                        Se = m%yn(i+1,j+1) - m%yn(i+1,j)                        
                        c%Fe_diff(i,j) = (c%mat_c(i+1,j) - c%mat_c(i,j))/m%dx * Se
                end do
        end do 

        ! Conditions aux limites
        c%Fe_diff(m%nx-1,:) = 0.

end subroutine cal_Fe_diff

subroutine cal_Fo_diff(c,m)
    use m_type 
    implicit none

    type(conc), intent(inout) :: c
        type(maillage), intent(in) :: m
        integer :: i,j
        real :: So

        allocate(c%Fo_diff(m%nx-1,m%ny-1))

        do i = 2,m%nx-1
                do j = 1,m%ny-1
                        So = m%yn(i,j+1) - m%yn(i,j)                        
                        c%Fe_diff(i,j) = (c%mat_c(i,j) - c%mat_c(i-1,j))/m%dx * So
                end do
        end do 

        ! Conditions aux limites
        c%Fo_diff(1,:) = 0.

end subroutine cal_Fo_diff

subroutine cal_Fn_diff(c,m,p)
    use m_type 
    implicit none

    type(conc), intent(inout) :: c
        type(maillage), intent(in) :: m
        type(phys), intent(in) :: p
        integer :: i,j
        real :: Sn, dyv
	
	Sn = m%dx
        allocate(c%Fn_diff(m%nx-1,m%ny-1))

        do i = 1,m%nx-1
                do j = 1,m%ny-2
                	dyv = 1./2.*(m%yn(i,j+2) - m%yn(i,j))
                        c%Fn_diff(i,j) = (c%mat_c(i,j+1) - c%mat_c(i,j))/dyv * Sn
                end do
        end do 

        ! Conditions aux limites
      
        do i = 1,m%nx-1
        	dyv = 1./2.*(m%yn(i,m%ny) - m%yn(i,m%ny-1))
        	c%Fn_diff(i,m%ny-1) = (p%c0 - c%mat_c(i,m%ny-1))/dyv * Sn
        end do
end subroutine cal_Fn_diff

subroutine cal_Fs_diff(c,m,p)
    use m_type 
    implicit none

    type(conc), intent(inout) :: c
        type(maillage), intent(in) :: m
        type(phys), intent(in) :: p
        integer :: i,j
        real :: Ss, dyv
	
	Ss = m%dx
        allocate(c%Fs_diff(m%nx-1,m%ny-1))

        do i = 1,m%nx-1
                do j = 2,m%ny-1
                	dyv = 1./2.*(m%yn(i,j+1) - m%yn(i,j-1))
                        c%Fn_diff(i,j) = (c%mat_c(i,j) - c%mat_c(i,j-1))/dyv * Ss
                end do
        end do 

        ! Conditions aux limites
      
        do i = 1,m%nx-1
        	dyv = 1./2.*(m%yn(i,2) - m%yn(i,1))
        	c%Fn_diff(i,1) = (c%mat_c(i,1) - p%c1)/dyv * Ss
        end do
end subroutine cal_Fs_diff

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
                T(i,j) = abs(m%u(i,j))/(p%CFL*m%dx)+abs(m%v(i,j))/(p%CFL*dyn)&
                +p%D/(p%R*(m%dx**2))+p%D/(p%R*(dyn**2))
            end do
        end do
        m%dt = 1./maxval(T)
end subroutine pdt


