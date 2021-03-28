program main
        ! Programme principal. 
        use m_type
        implicit none

        type(maillage) :: m
        type(phys) :: p
        type(conc) :: c
        integer :: i
        
        write(*,*) "[I] Lecture des données."
        call lecture_donnee(p,m)
        write(*,*) "[I] Définition du maillage."
        call def_maillage(p,m)
        write(*,*) "[I] Calcul des vitesses."
        call vitesse(p,m)

        allocate(c%mat_c(m%nx-1,m%ny-1))
        ! Conditions initiales
        c%mat_c(:,:) = 0.
        write(*,*) "[I] Calcul du pas de temps dt."
        call pdt(p, m)
        write(*,*) "[I] Exportation des données."
        call VTSWriter(0.,0,m%Nx,m%Ny,m%xn,m%yn,c%mat_c,m%u,m%v,'ini')

        do i=1,m%nt
                call concentration(p, m, c) ! La matrice C(i,j)^n+1
                if (MOD(i,m%nt/100) == 0) then
                	call VTSWriter(real(i)*m%dt,i,m%nx,m%ny,m%xn,m%yn,c%mat_c,m%u,m%v,'int')
                end if
        end do
        call VTSWriter(m%nt*m%dt,m%nt,m%nx,m%ny,m%xn,m%yn,c%mat_c,m%u,m%v,'end')
        deallocate(m%xn,m%yn,c%mat_c,m%u,m%v)
        write(*,*) "[I] Fini."
end program main
