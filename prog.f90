program main
        ! Programme principal. 
        use m_type
        implicit none

        type(maillage) :: m
        type(phys) :: p
        type(conc) :: c
        real, dimension(:,:), allocatable :: T
        integer :: i
        real :: dt

        call lecture_donnee(p,m)
        call def_maillage(p,m)
        call vitesse(p,m)

        allocate(c%mat_c(m%nx-1,m%ny-1))
        ! Conditions initiales
        c%mat_c(:,:) = 0.
        c%mat_c(:,1) = p%c1
        c%mat_c(:,m%ny-1) = p%c0

        write(*,*) "[I] Exportation des données."
        call VTSWriter(0.,1,m%Nx,m%Ny,m%xn,m%yn,c%mat_c,m%u,m%v,'ini')
        dt = p%t_tot/real(m%nt-1)
        do i=2,m%nt
                conc_n = c%mat_c
                call concentration(c, m, conc_n, dt)
                if i == m%nt then
                        call VTSWriter(p%t_tot,m%nt,m%nx,m%ny,m%xn,m%yn,c%mat_c,m%u,m%v,'end')
                else
                        call VTSWriter(real(i-1)*dt,i,m%nx,m%ny,m%xn,m%yn,c%mat_c,m%u,m%v,'int')
                end if
        end do
        deallocate(m%xn,m%yn,c%mat_c,m%u,m%v)
        write(*,*) "[I] Fini."
end program main
