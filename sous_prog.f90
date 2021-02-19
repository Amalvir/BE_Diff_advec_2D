subroutine lecture_donnee(p,m)
        use m_type
        implicit none

        type(maillage), intent(out) :: m
        type(phys), intent(out) :: p

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
        use m_type
        implicit none

        type(maillage), intent(inout) :: m
        type(phys), intent(in) :: p
        integer :: i
        real :: y_reg, PI

        PI = acos(-1.)

        allocate(m%xn(m%Nx, m%Ny))
        allocate(m%yn(m%Nx, m%Ny))
        allocate(m%dyn(m%Ny-1))

        do i=1,m%Nx
        m%xn(i,:) = (i-1)*p%L/real(m%Nx-1)
        end do

        m%yn(:,1) = 0.
        do i=2,m%Ny
        y_reg = real(i-1)*p%L/real(m%Ny-1)
        m%yn(:,i) = y_reg - p%L/(3.*PI)*sin((2.*PI*y_reg)/p%L)
        m%dyn(i-1) = m%yn(1,i) - m%yn(1,i-1)
        !m%dyv(i) = 0.5*(m%dyn(i-1) + m%dyn(i))
        end do

end subroutine def_maillage

