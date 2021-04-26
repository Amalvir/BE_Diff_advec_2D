program main
        ! Programme principal. 
        use m_type
        implicit none

        type(maillage) :: m
        type(phys) :: p
        type(conc) :: c
        integer :: i,j
        logical :: cond
        real :: theorique,dyn
        
        write(*,*) "[I] Lecture des données."
        call lecture_donnee(p,m)
        write(*,*) "[I] Définition du maillage."
        call def_maillage(p,m)
        write(*,*) "[I] Calcul des vitesses."
        call vitesse(p,m)

        allocate(c%mat_c(m%nx-1,m%ny-1))
        ! Conditions initiales
        c%mat_c(:,1:(m%ny-1)/2) = p%c1
        c%mat_c(:,(m%ny-1)/2+1:m%ny-1) = p%c0
        write(*,*) "[I] Calcul du pas de temps dt."
        call pdt(p, m)
        write(*,*) "[I] Exportation des données."
        call VTSWriter(0.,0,m%Nx,m%Ny,m%xn,m%yn,c%mat_c,m%u,m%v,'ini')

        open(10,file="test_diff.csv")
        do i=1,m%nt
                call concentration(p, m, c) ! La matrice C(i,j)^n+1

                write(10,*) real(i)*m%dt
                write(10,'(A18,A18,A18)') "Hauteur","Numérique","Éxacte"
                do j=1,m%ny-1
                    dyn = m%yn(m%nx/2,j+1) - m%yn(m%nx/2,j)
                    write(10,*) m%yn(m%nx/2,j)+dyn/2., c%mat_c(m%nx/2,j), theorique(m%yn(m%nx/2,j)+dyn/2.,real(i)*m%dt,p)
                end do
                if (m%nt < 99) then
                        cond = .True.
                else
                        cond = MOD(i,m%nt/99) == 0
                end if

                if (cond) then
                        call VTSWriter(real(i)*m%dt,i,m%nx,m%ny,m%xn,m%yn,c%mat_c,m%u,m%v,'int')
                end if
        end do

        close(10)

        call VTSWriter(m%nt*m%dt,m%nt,m%nx,m%ny,m%xn,m%yn,c%mat_c,m%u,m%v,'end')
        deallocate(m%xn,m%yn,c%mat_c,m%u,m%v)
        write(*,*) "[I] Fini."
end program main
